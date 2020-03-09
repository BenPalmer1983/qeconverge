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

###########################################
#  CLASS conv_ecut2d
###########################################
class conv_ecut2d:

  @staticmethod
  def run():
  
    globals.log_fh.write('CONVERGE ECUT 2D \n')
    globals.log_fh.write('================ \n')
    

    ecutrho = globals.ecut2d['rho_min_pw']
    data_h = 0
    while(ecutrho<=globals.ecut2d['rho_max_pw']):
      ecutwfc = globals.ecut2d['wfc_min_pw']
      data_w = 0
      while(ecutwfc<=globals.ecut2d['wfc_max_pw']):
        # Get file name
        file_name = globals.file_name() + '.in'
        
        # Make and save file
        ecutfile = pwscf_input()
        ecutfile.load(globals.dirs['td'] + "/" + 'ecut.in')
        ecutfile.set_ecutwfc(ecutwfc)
        ecutfile.set_ecutrho(ecutrho)
        ecutfile.save(file_name, globals.dirs['ecut'])  
                
        # Run pwscf
        runfile = []
        runfile.append(globals.dirs['ecut'] + '/' + file_name)
        log, files_out, run_list = pwscf_exec.execute(runfile)
                
        fpo = pwscf_output(files_out[0]['file'])  
        globals.ecut2d['energy_ry'][data_h,data_w] = fpo.get_energy_per_atom()
        globals.ecut2d['force_rybohr'][data_h,data_w] = fpo.get_force_per_atom()   
        globals.ecut2d['energy_ev'][data_h,data_w] = units.convert(globals.pw_units['energy'], 'EV', fpo.get_energy_per_atom())
        globals.ecut2d['force_evang'][data_h,data_w] = units.convert(globals.pw_units['force'], 'EV/ANG', fpo.get_force_per_atom())    
        
        globals.log_fh.write(str(ecutwfc) + ' ' + str(ecutrho) + ' ' + str(globals.ecut2d['energy_ry'][data_h,data_w]) + ' ' + str(globals.ecut2d['force_rybohr'][data_h,data_w])  + '\n')
                         
        ecutwfc = ecutwfc + globals.ecut2d['wfc_inc_pw']
        data_w = data_w + 1
        
      ecutrho = ecutrho + globals.ecut2d['rho_inc_pw'] 
      data_h = data_h + 1

    
    globals.log_fh.write('\n')
    

    # Save csv
    numpy.savetxt(globals.dirs['csv'] + '/ecut2d_energy_ry.csv', globals.ecut2d['energy_ry'][0:data_h,0:data_w], fmt='%10.5f', delimiter=',')
    numpy.savetxt(globals.dirs['csv'] + '/ecut2d_energy_ev.csv', globals.ecut2d['energy_ev'][0:data_h,0:data_w], fmt='%10.5f', delimiter=',')
    numpy.savetxt(globals.dirs['csv'] + '/ecut2d_force_ry.csv', globals.ecut2d['force_rybohr'][0:data_h,0:data_w], fmt='%10.5f', delimiter=',')
    numpy.savetxt(globals.dirs['csv'] + '/ecut2d_force_ev.csv', globals.ecut2d['force_evang'][0:data_h,0:data_w], fmt='%10.5f', delimiter=',')
    

    # PLOT
    
    x = numpy.linspace(globals.ecut2d['wfc_min_pw'], globals.ecut2d['wfc_max_pw'], data_w)
    y = numpy.linspace(globals.ecut2d['rho_min_pw'], globals.ecut2d['rho_max_pw'], data_h)
    
    plt.clf()    
    plt.contourf(x, y, globals.ecut2d['energy_ry'][0:data_h,0:data_w], 20, cmap='Greys')
    plt.colorbar()
    plt.savefig(globals.dirs['plots'] + '/' + 'ecut2d_energy_ry.svg')

    plt.clf()
    plt.contourf(x, y, globals.ecut2d['force_rybohr'][0:data_h,0:data_w], 20, cmap='Greys')
    plt.colorbar()
    plt.savefig(globals.dirs['plots'] + '/' + 'ecut2d_force_rybohr.svg')

    a = units.convert(globals.pw_units['energy'], 'EV', globals.ecut2d['wfc_min_pw'])
    b = units.convert(globals.pw_units['energy'], 'EV', globals.ecut2d['wfc_max_pw']) 
    c = units.convert(globals.pw_units['charge_density'], 'ANG-3', globals.ecut2d['rho_min_pw'])
    d = units.convert(globals.pw_units['charge_density'], 'ANG-3', globals.ecut2d['rho_max_pw']) 
    x = numpy.linspace(a, b, data_w)
    y = numpy.linspace(c, d, data_h)
    
    plt.clf()    
    plt.contourf(x, y, globals.ecut2d['energy_ev'][0:data_h,0:data_w], 20, cmap='Greys')
    plt.colorbar()
    plt.savefig(globals.dirs['plots'] + '/' + 'ecut2d_energy_ev.svg')

    plt.clf()
    plt.contourf(x, y, globals.ecut2d['force_evang'][0:data_h,0:data_w], 20, cmap='Greys')
    plt.colorbar()
    plt.savefig(globals.dirs['plots'] + '/' + 'ecut2d_force_evang.svg')
    