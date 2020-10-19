################################################################
#    Processing PWscf input file
#
#
#
#
################################################################


import numpy
from units import units
from plot_output import plot_output

###########################################
#  CLASS pwscf_converg
###########################################
class conv_ecut:

  @staticmethod
  def run():
  
    print("ECUT")
    globals.log_fh.write('CONVERGE ECUT \n')
    globals.log_fh.write('============= \n')
    
    
    ef = pwscf_input()
    ef.load("input_template.in", globals.dirs['td'])
    ef.set_dirs()
    #ef.load_config(globals.ecut['type'], globals.ecut['size'], globals.ecut['alat_expanded_pw'])
    
    s = {
         'type': globals.ecut['type'],
         'labels': None,
         'size_x': globals.ecut['size'],
         'size_y': globals.ecut['size'],
         'size_z': globals.ecut['size'],
        }       
    ef.set_alat(globals.ecut['alat_pw'])
    ef.set_config(s)
    
    ef.set_cp_arr(globals.ecut['cp'])
    ef.set_k_points(globals.ecut['kpoints_type'], globals.ecut['kpoints_val'])
    ef.rand_vary(globals.ecut['rand_variance'], globals.ecut['rand_seed'])
    ef.save("ecut.in", globals.dirs['td'])  
    
    
    #################################
    # 1. INCREASE ECUTWFC 
    #################################
    
    globals.log_fh.write('1. INCREASE ECUTWFC \n')
    
    converged = 0
    ecutwfc = globals.ecut['wfc_min_pw']

    # CONVERGENCE FOR EACH CHANGE BY 1 RY 
    econv_threshold = globals.ecut['c_energy_pw'] * globals.ecut['wfc_inc']
    fconv_threshold = globals.ecut['c_force_pw'] * globals.ecut['wfc_inc']
    
    
    globals.log_fh.write('econv_threshold: ' + str(econv_threshold) + ' \n')
    globals.log_fh.write('fconv_threshold: ' + str(fconv_threshold) + ' \n')
    

    ecutwfc_v = None
    ecutrho_v = None
    counter = 0   
    globals.log_fh.write('Data: \n')
    while((ecutwfc <= globals.ecut['wfc_max_pw'] and converged < 2) or counter < 5):
      counter = counter + 1

      # Get file name
      file_name = globals.file_name() + '.in'

      # Make and save file
      ecutfile = pwscf_input()
      ecutfile.load(globals.dirs['td'] + "/" + 'ecut.in')
      ecutfile.set_ecutwfc(ecutwfc)
      ecutfile.set_ecutrho(4 * ecutwfc)
      ecutfile.save(file_name, globals.dirs['ecut'])  

      # Run pwscf
      runfile = []
      runfile.append(globals.dirs['ecut'] + '/' + file_name)
      log, files_out, run_list = pwscf_exec.execute(runfile)

      fpo = pwscf_output(files_out[0]['file'])  
      cpu, wall = fpo.get_times()
      globals.set_times(cpu, wall)
      
      n = globals.ecut['data_len'][0]
      globals.ecut['data_len'][0] = globals.ecut['data_len'][0] + 1

      globals.ecut['data'][0,n] = ecutwfc
      globals.ecut['data'][1,n] = 4 * ecutwfc
      globals.ecut['data'][2,n] = fpo.get_energy_per_atom()
      globals.ecut['data'][3,n] = fpo.get_force_per_atom()
      
      econv = None
      fconv = None
      if(counter > 1):
        econv = round(abs(globals.ecut['data'][2,n] - globals.ecut['data'][2,n-1]),8)
        fconv = round(abs(globals.ecut['data'][3,n] - globals.ecut['data'][3,n-1]),8)
        if(converged < 2):
          ecutwfc_v = ecutwfc - globals.ecut['wfc_inc_pw']
        if(abs(globals.ecut['data'][2,n] - globals.ecut['data'][2,n-1]) <= econv_threshold and abs(globals.ecut['data'][3,n] - globals.ecut['data'][3,n-1]) <= fconv_threshold):
          converged = converged + 1
        else:
          converged = 0

      globals.log_fh.write('INC WFC ' + 
                           std.str_padded(counter, 5) + 
                           std.str_padded(converged, 5) + 
                           std.str_padded(globals.ecut['data'][0,n], 10) + 
                           std.str_padded(globals.ecut['data'][1,n], 10) + 
                           std.str_padded(globals.ecut['data'][2,n], 17) + 
                           std.str_padded(globals.ecut['data'][3,n], 17) + 
                           std.str_padded(econv, 17) + 
                           std.str_padded(fconv, 17) + 
                           str(ecutwfc_v) + ' ' +
                           str(ecutrho_v) + '\n')
      
      # Increment
      ecutwfc = ecutwfc + globals.ecut['wfc_inc_pw']   
      
    # Record value  
    globals.log_fh.write('Ecutwfc: ' + str(ecutwfc_v) + ' \n')
    globals.log_fh.write('Ecutrho: ' + str(4 * ecutwfc_v) + ' \n')
    
    
    
    #################################
    # 2. DECREASE ECUTRHO 
    #################################
    
    globals.log_fh.write('2. DECREASE ECUTRHO \n')
    
    notconverged = 0
    ecutwfc = ecutwfc_v
    ecutrho = 4.0 * ecutwfc_v
    ecutrho_min = globals.ecut['wr_min_ratio'] * ecutwfc
    
    # CONVERGENCE FOR EACH CHANGE BY 1 RY 
    econv_threshold = globals.ecut['c_energy_pw'] * globals.ecut['rho_dec']
    fconv_threshold = globals.ecut['c_force_pw'] * globals.ecut['rho_dec']
    
    
    globals.log_fh.write('econv_threshold: ' + str(econv_threshold) + ' \n')
    globals.log_fh.write('fconv_threshold: ' + str(fconv_threshold) + ' \n')
    

    ecutrho_v = ecutrho
    counter = 0   
    globals.log_fh.write('Data: \n')
    while((ecutrho >= ecutrho_min and notconverged < 1) or counter < 5):
      counter = counter + 1

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
      cpu, wall = fpo.get_times()
      globals.set_times(cpu, wall)
      
      n = globals.ecut['data_len'][1]
      globals.ecut['data_len'][1] = globals.ecut['data_len'][1] + 1

      globals.ecut['data'][5,n] = ecutwfc
      globals.ecut['data'][6,n] = ecutrho
      globals.ecut['data'][7,n] = fpo.get_energy_per_atom()
      globals.ecut['data'][8,n] = fpo.get_force_per_atom()
     
      #globals.log_fh.write(str(globals.ecut['data'][5,n]) + ' ' + str(globals.ecut['data'][6,n]) + ' ' + str(globals.ecut['data'][7,n]) + ' ' + str(globals.ecut['data'][8,n]) + ' ' + '\n')
      econv = None
      fconv = None
      if(counter > 1):
        econv = round(abs(globals.ecut['data'][7,n] - globals.ecut['data'][7,n-1]),8)
        fconv = round(abs(globals.ecut['data'][8,n] - globals.ecut['data'][8,n-1]),8)
        if(abs(globals.ecut['data'][7,n] - globals.ecut['data'][7,n-1]) > econv_threshold or abs(globals.ecut['data'][8,n] - globals.ecut['data'][8,n-1]) > fconv_threshold):
          notconverged = notconverged + 1       
        ecutrho_v = ecutrho + globals.ecut['rho_dec_pw']  # BEST IS ALWAYS LAST ECUTRHO (i.e. this ecutrho + the difference)
          

      globals.log_fh.write('DEC RHO ' + 
                           std.str_padded(counter, 5) + 
                           std.str_padded(converged, 5) + 
                           std.str_padded(globals.ecut['data'][5,n], 10) + 
                           std.str_padded(globals.ecut['data'][6,n], 10) + 
                           std.str_padded(globals.ecut['data'][7,n], 17) + 
                           std.str_padded(globals.ecut['data'][8,n], 17) + 
                           std.str_padded(econv, 17) + 
                           std.str_padded(fconv, 17) + 
                           str(ecutwfc_v) + ' ' +
                           str(ecutrho_v) + '\n')

      
      # Decrease
      ecutrho = ecutrho - globals.ecut['rho_dec_pw']   
    
    # Record value  
    globals.log_fh.write('Ecutwfc: ' + str(ecutwfc_v) + ' \n')
    globals.log_fh.write('Ecutrho: ' + str(ecutrho_v) + ' \n')
    
    
    #################################
    # 3. DECREASE ECUTWFC 
    #################################
    
    globals.log_fh.write('3. DECREASE ECUTWFC \n')
    
    unconverged = 0
    ecutwfc = ecutwfc_v
    ecutrho = ecutrho_v

    # CONVERGENCE FOR EACH CHANGE BY 1 RY 
    econv_threshold = globals.ecut['c_energy_pw'] * globals.ecut['wfc_dec']
    fconv_threshold = globals.ecut['c_force_pw'] * globals.ecut['wfc_dec']
    
    
    globals.log_fh.write('econv_threshold: ' + str(econv_threshold) + ' \n')
    globals.log_fh.write('fconv_threshold: ' + str(fconv_threshold) + ' \n')
    

    counter = 0   
    globals.log_fh.write('Data: \n')
    while((ecutwfc >= globals.ecut['wfc_min_pw'] and unconverged < 1) or counter < 5):
      counter = counter + 1

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
      cpu, wall = fpo.get_times()
      globals.set_times(cpu, wall)    
      
      n = globals.ecut['data_len'][2]
      globals.ecut['data_len'][2] = globals.ecut['data_len'][2] + 1

      globals.ecut['data'][10,n] = ecutwfc
      globals.ecut['data'][11,n] = ecutrho
      globals.ecut['data'][12,n] = fpo.get_energy_per_atom()
      globals.ecut['data'][13,n] = fpo.get_force_per_atom()
     
      #globals.log_fh.write(str(globals.ecut['data'][10,n]) + ' ' + str(globals.ecut['data'][11,n]) + ' ' + str(globals.ecut['data'][12,n]) + ' ' + str(globals.ecut['data'][13,n]) + ' ' + '\n')

      econv = None
      fconv = None
      if(counter > 1):
        econv = round(abs(globals.ecut['data'][12,n] - globals.ecut['data'][12,n-1]),8)
        fconv = round(abs(globals.ecut['data'][13,n] - globals.ecut['data'][13,n-1]),8)
        if(abs(globals.ecut['data'][12,n] - globals.ecut['data'][12,n-1]) <= econv_threshold and abs(globals.ecut['data'][13,n] - globals.ecut['data'][13,n-1]) <= fconv_threshold):
          notconverged = notconverged + 1  
        ecutwfc_v = ecutwfc + globals.ecut['wfc_dec_pw']  # BEST IS ALWAYS LAST ECUTRHO (i.e. this ecutrho + the difference)  

      globals.log_fh.write('DEC WFC ' + 
                           std.str_padded(counter, 5) + 
                           std.str_padded(converged, 5) + 
                           std.str_padded(globals.ecut['data'][10,n], 10) + 
                           std.str_padded(globals.ecut['data'][11,n], 10) + 
                           std.str_padded(globals.ecut['data'][12,n], 17) + 
                           std.str_padded(globals.ecut['data'][13,n], 17) + 
                           std.str_padded(econv, 17) + 
                           std.str_padded(fconv, 17) + 
                           str(ecutwfc_v) + ' ' +
                           str(ecutrho_v) + '\n')
      
      # Increment
      ecutwfc = ecutwfc - globals.ecut['wfc_dec_pw']   
      
    # Record value 
    globals.log_fh.write('Ecutwfc: ' + str(ecutwfc_v) + ' \n') 
    globals.log_fh.write('Ecutrho: ' + str(ecutrho_v) + ' \n')    
    
    globals.ecut['ecutwfc_converged_pw'] = ecutwfc_v
    globals.ecut['ecutrho_converged_pw'] = ecutrho_v
    
    
    globals.log_fh.write('\n')
    
    fh = open(globals.dirs['wd'] + '/' + 'ecut_converge.txt', 'w')
    fh.write('Ecutwfc: ' + str(ecutwfc_v) + ' \n') 
    fh.write('Ecutrho: ' + str(ecutrho_v) + ' \n') 
    fh.close()
    
    # CONVERSIONS
    
    n = globals.ecut['data_len'][0]    
    globals.ecut['data_relative'][0,0:n] = globals.ecut['data'][0,0:n]
    globals.ecut['data_relative'][1,0:n] = globals.ecut['data'][1,0:n]
    globals.ecut['data_relative'][2,0:n] = globals.ecut['data'][2,0:n] - numpy.amin(globals.ecut['data'][2,0:n])
    globals.ecut['data_relative'][3,0:n] = globals.ecut['data'][3,0:n] - numpy.amin(globals.ecut['data'][3,0:n])    
    for i in range(n):
      globals.ecut['data_ev'][0,i] = units.convert(globals.pw_units['energy'], 'EV', globals.ecut['data'][0,i])
      globals.ecut['data_ev'][1,i] = units.convert(globals.pw_units['charge_density'], 'ANG-3', globals.ecut['data'][1,i])
      globals.ecut['data_ev'][2,i] = units.convert(globals.pw_units['energy'], 'EV', globals.ecut['data'][2,i])
      globals.ecut['data_ev'][3,i] = units.convert(globals.pw_units['force'], 'EV/ANG', globals.ecut['data'][3,i])
        
    n = globals.ecut['data_len'][1]    
    globals.ecut['data_relative'][5,0:n] = globals.ecut['data'][5,0:n]
    globals.ecut['data_relative'][6,0:n] = globals.ecut['data'][6,0:n]
    globals.ecut['data_relative'][7,0:n] = globals.ecut['data'][7,0:n] - numpy.amin(globals.ecut['data'][7,0:n])
    globals.ecut['data_relative'][8,0:n] = globals.ecut['data'][8,0:n] - numpy.amin(globals.ecut['data'][8,0:n])    
    for i in range(n):
      globals.ecut['data_ev'][5,i] = units.convert(globals.pw_units['energy'], 'EV', globals.ecut['data'][5,i])
      globals.ecut['data_ev'][6,i] = units.convert(globals.pw_units['charge_density'], 'ANG-3', globals.ecut['data'][6,i])
      globals.ecut['data_ev'][7,i] = units.convert(globals.pw_units['energy'], 'EV', globals.ecut['data'][7,i])
      globals.ecut['data_ev'][8,i] = units.convert(globals.pw_units['force'], 'EV/ANG', globals.ecut['data'][8,i])
        
    n = globals.ecut['data_len'][2]    
    globals.ecut['data_relative'][10,0:n] = globals.ecut['data'][10,0:n]
    globals.ecut['data_relative'][11,0:n] = globals.ecut['data'][11,0:n]
    globals.ecut['data_relative'][12,0:n] = globals.ecut['data'][12,0:n] - numpy.amin(globals.ecut['data'][12,0:n])
    globals.ecut['data_relative'][13,0:n] = globals.ecut['data'][13,0:n] - numpy.amin(globals.ecut['data'][13,0:n])   
    for i in range(n):
      globals.ecut['data_ev'][10,i] = units.convert(globals.pw_units['energy'], 'EV', globals.ecut['data'][10,i])
      globals.ecut['data_ev'][11,i] = units.convert(globals.pw_units['charge_density'], 'ANG-3', globals.ecut['data'][11,i])
      globals.ecut['data_ev'][12,i] = units.convert(globals.pw_units['energy'], 'EV', globals.ecut['data'][12,i])
      globals.ecut['data_ev'][13,i] = units.convert(globals.pw_units['force'], 'EV/ANG', globals.ecut['data'][13,i])
    
    
    
    # PLOT CONVERGENCE RY

    plt.clf()    
    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')
    fig, axs = plt.subplots(3, 2, figsize=(12,9))
    fig.tight_layout(pad=5.0)    
    fig.suptitle('Convergence - Relative Values')
    
    n = globals.ecut['data_len'][0]    
    axs[0, 0].plot(globals.ecut['data'][0,0:n], globals.ecut['data'][2,0:n], color='k', ls='solid')
    axs[0, 0].set_title('Increasing Ecutwfc')
    axs[0, 0].set_xlabel('Ecutwfc (RY)')
    axs[0, 0].set_ylabel('Energy (RY)')
    axs[0, 1].plot(globals.ecut['data'][0,0:n], globals.ecut['data'][3,0:n], color='k', ls='solid')
    axs[0, 1].set_title('Increasing Ecutwfc')
    axs[0, 1].set_xlabel('Ecutwfc (RY)')
    axs[0, 1].set_ylabel('Force (RY/BOHR)')
    
    n = globals.ecut['data_len'][1]    
    axs[1, 0].plot(globals.ecut['data'][6,0:n], globals.ecut['data'][7,0:n], color='k', ls='solid')
    axs[1, 0].set_title('Decreasing Ecutrho')
    axs[1, 0].set_xlabel('Ecutrho (BOHR-3)')
    axs[1, 0].set_ylabel('Energy (RY)')    
    axs[1, 1].plot(globals.ecut['data'][6,0:n], globals.ecut['data'][8,0:n], color='k', ls='solid')
    axs[1, 1].set_title('Decreasing Ecutrho')
    axs[1, 1].set_xlabel('Ecutrho (BOHR-3)')
    axs[1, 1].set_ylabel('Force (RY/BOHR)')
    
    n = globals.ecut['data_len'][2]    
    axs[2, 0].plot(globals.ecut['data'][10,0:n], globals.ecut['data'][12,0:n], color='k', ls='solid')
    axs[2, 0].set_title('Increasing Ecutwfc')
    axs[2, 0].set_xlabel('Ecutwfc (RY)')
    axs[2, 0].set_ylabel('Energy (RY)')    
    axs[2, 1].plot(globals.ecut['data'][10,0:n], globals.ecut['data'][13,0:n], color='k', ls='solid')
    axs[2, 1].set_title('Increasing Ecutwfc')
    axs[2, 1].set_xlabel('Ecutwfc (RY)')
    axs[2, 1].set_ylabel('Force (RY/BOHR)')

    #plt.savefig(globals.dirs['plots'] + '/' + 'convergence_ry.svg')
    plot_output.plot(plt, 'ecut_convergence_ry')
  
    
    # PLOT CONVERGENCE REL RY
    
    plt.clf()    
    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')
    fig, axs = plt.subplots(3, 2, figsize=(12,9))
    fig.tight_layout(pad=5.0)    
    fig.suptitle('Convergence - Adjusted Values')
    
    n = globals.ecut['data_len'][0]    
    axs[0, 0].plot(globals.ecut['data_relative'][0,0:n], globals.ecut['data_relative'][2,0:n], color='k', ls='solid')
    axs[0, 0].set_title('Increasing Ecutwfc')
    axs[0, 0].set_xlabel('Ecutwfc (RY)')
    axs[0, 0].set_ylabel('Energy (RY)')
    axs[0, 1].plot(globals.ecut['data_relative'][0,0:n], globals.ecut['data_relative'][3,0:n], color='k', ls='solid')
    axs[0, 1].set_title('Increasing Ecutwfc')
    axs[0, 1].set_xlabel('Ecutwfc (RY)')
    axs[0, 1].set_ylabel('Force (RY/BOHR)')
    
    n = globals.ecut['data_len'][1]    
    axs[1, 0].plot(globals.ecut['data_relative'][6,0:n], globals.ecut['data_relative'][7,0:n], color='k', ls='solid')
    axs[1, 0].set_title('Decreasing Ecutrho')
    axs[1, 0].set_xlabel('Ecutrho (BOHR-3)')
    axs[1, 0].set_ylabel('Energy (RY)')    
    axs[1, 1].plot(globals.ecut['data_relative'][6,0:n], globals.ecut['data_relative'][8,0:n], color='k', ls='solid')
    axs[1, 1].set_title('Decreasing Ecutrho')
    axs[1, 1].set_xlabel('Ecutrho (BOHR-3)')
    axs[1, 1].set_ylabel('Force (RY/BOHR)')
    
    n = globals.ecut['data_len'][2]    
    axs[2, 0].plot(globals.ecut['data_relative'][10,0:n], globals.ecut['data_relative'][12,0:n], color='k', ls='solid')
    axs[2, 0].set_title('Increasing Ecutwfc')
    axs[2, 0].set_xlabel('Ecutwfc (RY)')
    axs[2, 0].set_ylabel('Energy (RY)')    
    axs[2, 1].plot(globals.ecut['data_relative'][10,0:n], globals.ecut['data_relative'][13,0:n], color='k', ls='solid')
    axs[2, 1].set_title('Increasing Ecutwfc')
    axs[2, 1].set_xlabel('Ecutwfc (RY)')
    axs[2, 1].set_ylabel('Force (RY/BOHR)')

    #plt.savefig(globals.dirs['plots'] + '/' + 'convergence_adjusted_ry.svg')
    plot_output.plot(plt, 'ecut_convergence_adjusted_ry')
  
    
    # PLOT CONVERGENCE EV
    
    plt.clf()    
    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')
    fig, axs = plt.subplots(3, 2, figsize=(12,9))
    fig.tight_layout(pad=5.0)    
    fig.suptitle('Convergence - Adjusted Values')
    
    n = globals.ecut['data_len'][0]    
    axs[0, 0].plot(globals.ecut['data_ev'][0,0:n], globals.ecut['data_ev'][2,0:n], color='k', ls='solid')
    axs[0, 0].set_title('Increasing Ecutwfc')
    axs[0, 0].set_xlabel('Ecutwfc (eV)')
    axs[0, 0].set_ylabel('Energy (eV)')
    axs[0, 1].plot(globals.ecut['data_ev'][0,0:n], globals.ecut['data_ev'][3,0:n], color='k', ls='solid')
    axs[0, 1].set_title('Increasing Ecutwfc')
    axs[0, 1].set_xlabel('Ecutwfc (eV)')
    axs[0, 1].set_ylabel('Force (eV/ang)')
    
    n = globals.ecut['data_len'][1]    
    axs[1, 0].plot(globals.ecut['data_ev'][6,0:n], globals.ecut['data_ev'][7,0:n], color='k', ls='solid')
    axs[1, 0].set_title('Decreasing Ecutrho')
    axs[1, 0].set_xlabel('Ecutrho (ang-3)')
    axs[1, 0].set_ylabel('Energy (eV)')    
    axs[1, 1].plot(globals.ecut['data_ev'][6,0:n], globals.ecut['data_ev'][8,0:n], color='k', ls='solid')
    axs[1, 1].set_title('Decreasing Ecutrho')
    axs[1, 1].set_xlabel('Ecutrho (ang-3)')
    axs[1, 1].set_ylabel('Force (eV/ang)')
    
    n = globals.ecut['data_len'][2]    
    axs[2, 0].plot(globals.ecut['data_ev'][10,0:n], globals.ecut['data_ev'][12,0:n], color='k', ls='solid')
    axs[2, 0].set_title('Increasing Ecutwfc')
    axs[2, 0].set_xlabel('Ecutwfc (eV)')
    axs[2, 0].set_ylabel('Energy (eV)')    
    axs[2, 1].plot(globals.ecut['data_ev'][10,0:n], globals.ecut['data_ev'][13,0:n], color='k', ls='solid')
    axs[2, 1].set_title('Increasing Ecutwfc')
    axs[2, 1].set_xlabel('Ecutwfc (eV)')
    axs[2, 1].set_ylabel('Force (eV/ang)')

    #plt.savefig(globals.dirs['plots'] + '/' + 'convergence_ev.svg')
    plot_output.plot(plt, 'ecut_convergence_ev')
    
    

