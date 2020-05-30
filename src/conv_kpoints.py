################################################################
#    Processing PWscf input file
#
#
#
#
################################################################


import numpy
from units import units
from globals import globals
from plot_output import plot_output

###########################################
#  CLASS conv_ecut2d
###########################################
class conv_kpoints:

  @staticmethod
  def run():
  
    # Make and save file
    kf = pwscf_input()
    kf.load("input_template.in", globals.dirs['td'])
    kf.set_dirs()
    kf.load_config(globals.kconv['type'], globals.kconv['size'], globals.kconv['alat_expanded_pw'])
    kf.set_cp_arr(globals.kconv['cp'])
    kf.set_ecutwfc(globals.ecut['ecutwfc_converged_pw'])
    kf.set_ecutrho(globals.ecut['ecutrho_converged_pw']) 
    kf.rand_vary(globals.ecut['rand_variance'], globals.ecut['rand_seed'])
    kf.save("kpoints.in", globals.dirs['td'])  
    
    
    globals.log_fh.write('CONVERGE KPOINTS \n')
    globals.log_fh.write('================ \n')
    
    k = globals.kconv['k_min']
    i = 0
    data_w = 0
    while(k<=globals.kconv['k_max']):
      k_str = str(k) + ' ' + str(k) + ' ' + str(k) + ' ' + str(globals.kconv['k_offset'][0]) + ' ' + str(globals.kconv['k_offset'][1]) + ' ' + str(globals.kconv['k_offset'][2])
      
      data_h = 0
      for j in range(len(globals.kconv['smearing_pw'])):
 
        # Make and save file
        file_name = globals.file_name() + '.in'
        kfile = pwscf_input()
        kfile.load(globals.dirs['td'] + "/" + 'kpoints.in')
        kfile.set_k_points(globals.kconv['k_type'], k_str)
        kfile.set_degauss(globals.kconv['smearing_pw'][j])
        kfile.save(file_name, globals.dirs['kpoints'])  
    
        # Run pwscf
        
        runfile = []
        runfile.append(globals.dirs['kpoints'] + '/' + file_name)
        log, files_out, run_list = pwscf_exec.execute(runfile)
                
        fpo = pwscf_output(files_out[0]['file'])  
        cpu, wall = fpo.get_times()
        globals.set_times(cpu, wall)
        globals.kconv['energy_ry'][data_h,data_w] = fpo.get_energy_per_atom()
        globals.kconv['force_rybohr'][data_h,data_w] = fpo.get_force_per_atom()   
        globals.kconv['energy_ev'][data_h,data_w] = units.convert(globals.pw_units['energy'], 'EV', fpo.get_energy_per_atom())
        globals.kconv['force_evang'][data_h,data_w] = units.convert(globals.pw_units['force'], 'EV/ANG', fpo.get_force_per_atom()) 
            
        globals.log_fh.write(str(k_str) + ' ' + str(globals.kconv['smearing_pw'][j]) + ' ' + str(globals.kconv['energy_ry'][data_h,data_w]) + ' ' + str(globals.kconv['force_rybohr'][data_h,data_w])  + '\n')
                         
        data_h = data_h + 1
      k = k + globals.kconv['k_inc']
      data_w = data_w + 1
  

    
    globals.log_fh.write('\n')
    

    # Save csv
    numpy.savetxt(globals.dirs['csv'] + '/kpoints_energy_ry.csv', globals.kconv['energy_ry'][0:data_h,0:data_w], fmt='%18.10f', delimiter=',')
    numpy.savetxt(globals.dirs['csv'] + '/kpoints_energy_ev.csv', globals.kconv['energy_ev'][0:data_h,0:data_w], fmt='%18.10f', delimiter=',')
    numpy.savetxt(globals.dirs['csv'] + '/kpoints_force_ry.csv', globals.kconv['force_rybohr'][0:data_h,0:data_w], fmt='%18.10f', delimiter=',')
    numpy.savetxt(globals.dirs['csv'] + '/kpoints_force_ev.csv', globals.kconv['force_evang'][0:data_h,0:data_w], fmt='%18.10f', delimiter=',')

    fh = open(globals.dirs['csv'] + '/kpoints_xy_lims.csv', 'w')
    fh.write(str(globals.kconv['k_min']) + '\n')
    fh.write(str(globals.kconv['k_max']) + '\n')
    for j in range(len(globals.kconv['smearing_pw'])):
      fh.write(str(str(globals.kconv['smearing_pw'][j])) + '\n')
    fh.write('\n')
    fh.close()
    



    globals.kconv['data_w'] = data_w
    globals.kconv['data_h'] = data_h
    
    
    

  @staticmethod
  def conv_w(data_in, h, w):
    out = numpy.zeros((h, w-1),)
    for i in range(w-1):
      out[:,i] = abs(data_in[:,i+1] - data_in[:,i])
    return out
    
  @staticmethod
  def conv_h(data_in, h, w):
    out = numpy.zeros((h-1, w),)
    for i in range(h-1):
      out[i,:] = abs(data_in[i+1,:] - data_in[i,:])
    return out
    
  @staticmethod
  def half_scale(arr_in):
    arr_out = numpy.zeros((len(arr_in) - 1,),)
    for i in range(len(arr_out)):
      arr_out[i] = 0.5 * (arr_in[i+1] + arr_in[i])
    return arr_out

