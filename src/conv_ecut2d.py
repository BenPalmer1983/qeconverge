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
from plots import plots
from matplotlib.colors import LinearSegmentedColormap

###########################################
#  CLASS conv_ecut2d
###########################################
class conv_ecut2d:

  @staticmethod
  def run():
  
    globals.log_fh.write('CONVERGE ECUT 2D \n')
    globals.log_fh.write('================ \n')
    print(globals.ecut2d['wfc_min_pw'], globals.ecut2d['wfc_max_pw'])
    print(globals.ecut2d['rho_min_pw'], globals.ecut2d['rho_max_pw'])
    ecutwfc = globals.ecut2d['wfc_min_pw']
    data_h = 0
    while(ecutwfc<=globals.ecut2d['wfc_max_pw']):
      ecutrho_m = globals.ecut2d['rho_min_pw']
      data_w = 0
      while(ecutrho_m<=globals.ecut2d['rho_max_pw']):
        # Get file name
        file_name = globals.file_name() + '.in'
        
        ecutrho = ecutwfc * ecutrho_m
        
        # Make and save file
        ecutfile = pwscf_input()
        ecutfile.load(globals.dirs['td'] + "/" + 'ecut.in')
        ecutfile.set_ecutwfc(ecutwfc)
        ecutfile.set_ecutrho(ecutrho)
        ecutfile.save(file_name, globals.dirs['ecut'])  
        
                
                  
        runfile = []
        runfile.append(globals.dirs['ecut'] + '/' + file_name)
        log, files_out, run_list = pwscf_exec.execute(runfile)
        
          
        fpo = pwscf_output(files_out[0]['file'])  
        cpu, wall = fpo.get_times()
        globals.set_times(cpu, wall)
        globals.ecut2d['energy_ry'][data_h,data_w] = fpo.get_energy_per_atom()
        globals.ecut2d['force_rybohr'][data_h,data_w] = fpo.get_force_per_atom()   
        globals.ecut2d['energy_ev'][data_h,data_w] = units.convert(globals.pw_units['energy'], 'EV', fpo.get_energy_per_atom())
        globals.ecut2d['force_evang'][data_h,data_w] = units.convert(globals.pw_units['force'], 'EV/ANG', fpo.get_force_per_atom()) 

        ecutrho_m = ecutrho_m + globals.ecut2d['rho_inc_pw'] 
        data_w = data_w + 1
        
      ecutwfc = ecutwfc + globals.ecut2d['wfc_inc_pw']        
      data_h = data_h + 1
      
      
    x = numpy.linspace(globals.ecut2d['rho_min_pw'], globals.ecut2d['rho_max_pw'], data_h)
    y = numpy.linspace(globals.ecut2d['wfc_min_pw'], globals.ecut2d['wfc_max_pw'], data_w)
      
    
    globals.log_fh.write('\n')
    

    # Save csv
    numpy.savetxt(globals.dirs['csv'] + '/ecut2d_energy_ry.csv', globals.ecut2d['energy_ry'][0:data_h,0:data_w], fmt='%18.10f', delimiter=',')
    numpy.savetxt(globals.dirs['csv'] + '/ecut2d_energy_ev.csv', globals.ecut2d['energy_ev'][0:data_h,0:data_w], fmt='%18.10f', delimiter=',')
    numpy.savetxt(globals.dirs['csv'] + '/ecut2d_force_ry.csv', globals.ecut2d['force_rybohr'][0:data_h,0:data_w], fmt='%18.10f', delimiter=',')
    numpy.savetxt(globals.dirs['csv'] + '/ecut2d_force_ev.csv', globals.ecut2d['force_evang'][0:data_h,0:data_w], fmt='%18.10f', delimiter=',')
    
    
    fh = open(globals.dirs['csv'] + '/ecut2d_xy_lims.csv', 'w')
    fh.write(str(globals.ecut2d['wfc_min_pw']) + '\n')
    fh.write(str(globals.ecut2d['wfc_max_pw']) + '\n')
    fh.write(str(globals.ecut2d['rho_min_pw']) + '\n')
    fh.write(str(globals.ecut2d['rho_max_pw']) + '\n')
    fh.close()
    
    
    
    globals.ecut2d['data_w'] = data_w
    globals.ecut2d['data_h'] = data_h
    
    
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
    
