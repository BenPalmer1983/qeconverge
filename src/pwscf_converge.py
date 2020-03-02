################################################################
#    Converging
#
#
#
#
################################################################



import time
import os
import sys
import numpy
import datetime
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits import mplot3d
from pwscf_input import pwscf_input
from pwscf_output import pwscf_output
from pwscf_exec import pwscf_exec
from lma import lma
from eos_fitting import eos_fitting
from fitting import fitting
from log import log
from read_config import read_config
from misc import misc
from units import units


######################
# 
######################

class pwscf_converge:

  def __init__(self):
    self.log = log()
    self.config = {}
    self.results = {}
    

  def run(self):
    
    # Start
    self.log.add("Running")

    # Load config file and defaults (or exit)
    result = self.load_file()
    if(result == False):
      self.exit()

    print(self.config)



    # Make dir and template
    self.make_dirs()

    # Prepare template
    self.prep_template()

    # Set file name counter
    self.fname_count = 10000

    # Run ecut
    self.ecut_conv()
    self.ecut_2d()

    # Run k-points
    self.kpoints_2d()

    # Exit
    self.exit()





  def make_dirs(self):
    self.wd = "wd"
    self.td = "wd/template"
    self.ecutd = "wd/ecut"
    self.kpointsd = "wd/kpoints"
    self.runpwscfd = "wd/runpwscf"
    self.plotsd = "wd/plots"
    self.csvd = "wd/csv"

    misc.makedir(self.wd) 
    misc.makedir(self.td) 
    misc.makedir(self.ecutd) 
    misc.makedir(self.kpointsd) 
    misc.makedir(self.runpwscfd) 
    misc.makedir(self.plotsd) 
    misc.makedir(self.csvd) 


  def prep_template(self):

    # TEMPLATE FOR ECUT

    template_file = "template.in"
    if("template" in self.config.keys()):
      template_file = read_config.get(self.config, 'template','file')


    if(not os.path.exists(template_file)):
      self.log.add("No template file")
      self.exit()

    config_type = read_config.get(self.config, 'configecut','type')
    config_size = read_config.get(self.config, 'configecut','size')
    config_alat = read_config.get(self.config, 'configecut','alat')
    config_alat_units = read_config.get(self.config, 'configecut','alat_units')
    config_cp = read_config.get_list(self.config, 'configecut','cp')
    config_rand = read_config.get_list(self.config, 'configecut','rand')


    # UNIT CONVERSIONS
    if(config_alat is not None):
      config_alat = units.convert(config_alat_units,'BOHR',float(config_alat))   # CONVERT INTO BOHR
      config_alat = float(config_size) * float(config_alat)


    # LOAD TEMPLATE FILE
    template = pwscf_input()
    template.load(template_file)
    template.set_dirs()
 
    # Set config
    if(config_type != None and config_size != None and config_alat != None):
      template.load_config(config_type.upper(), int(config_size), float(config_alat))

    # Set cell parameters
    if(config_cp != None):
      template.set_cp_arr(config_cp)

    if(config_rand != None):
      print(config_rand)
      template.rand_vary(config_rand[0], config_rand[1])
      #template.set_cp_arr(config_cp)


    template.save("input_template.in", self.td, 7)  
    self.tfile = self.td + "/" + "input_template.in"
    self.log.add("Template file loaded") 
    self.log.add(self.tfile)


    # TEMPLATE FOR KPOINTS
    config_type = read_config.get(self.config, 'configkpoint','type')
    config_size = read_config.get(self.config, 'configkpoint','size')
    config_alat = read_config.get(self.config, 'configkpoint','alat')
    config_alat_units = read_config.get(self.config, 'configkpoint','alat_units')
    config_cp = read_config.get_list(self.config, 'configkpoint','cp')
    config_rand = read_config.get_list(self.config, 'configkpoint','rand')


    # UNIT CONVERSIONS
    if(config_alat is not None):
      config_alat = units.convert(config_alat_units,'BOHR',float(config_alat))   # CONVERT INTO BOHR
      config_alat = float(config_size) * float(config_alat)

    # LOAD TEMPLATE FILE
    template_k = pwscf_input()
    template_k.load(template_file)
    template_k.set_dirs()
 
    # Set config
    if(config_type != None and config_size != None and config_alat != None):
      template_k.load_config(config_type.upper(), int(config_size), float(config_alat))

    # Set cell parameters
    if(config_cp != None):
      template_k.set_cp_arr(config_cp)

    if(config_rand != None):
      template_k.rand_vary(config_rand[0], config_rand[1])

    template_k.save("input_template_k.in", self.td, 7)  
    self.tkfile = self.td + "/" + "input_template_k.in"
    self.log.add("Template k-points file loaded") 
    self.log.add(self.tkfile)




  def ecut_conv(self):    

    # LOAD SETTINGS FROM INPUT FILE
    # ecutconv wfc_min=25 wfc_inc=5
    settings = {'run': True, 'econv': 0.1, 'econv_units': 'RY', 'fconv': 0.1, 'fconv_units': 'RY/BOHR', 'ecutwfcmin': 30, 'ecutwfcmax': '250', 'minwfcrhoratio': 2.0,}

    settings = self.add_to_dict(settings, 'econv', read_config.get(self.config, 'converge', 'energy'))
    settings = self.add_to_dict(settings, 'econv_units', read_config.get(self.config, 'converge', 'energy_units'))
    settings = self.add_to_dict(settings, 'fconv', read_config.get(self.config, 'converge', 'force'))
    settings = self.add_to_dict(settings, 'fconv_units', read_config.get(self.config, 'converge', 'force_units'))
    settings = self.add_to_dict(settings, 'run', read_config.get(self.config, 'ecutconv', 'run'))
    settings = self.add_to_dict(settings, 'wfc_min', read_config.get(self.config, 'ecutconv', 'wfc_min'))
    settings = self.add_to_dict(settings, 'wfc_inc', read_config.get(self.config, 'ecutconv', 'wfc_inc'))
    settings = self.add_to_dict(settings, 'wfc_dec', read_config.get(self.config, 'ecutconv', 'wfc_dec'))
    settings = self.add_to_dict(settings, 'wfc_max', read_config.get(self.config, 'ecutconv', 'wfc_max'))
    settings = self.add_to_dict(settings, 'rho_dec', read_config.get(self.config, 'ecutconv', 'rho_dec'))
    settings = self.add_to_dict(settings, 'kpointstype', read_config.get(self.config, 'ecutconv', 'kpointstype'))
    settings = self.add_to_dict(settings, 'kpointsval', read_config.get(self.config, 'ecutconv', 'kpointsval'))
    settings = self.add_to_dict(settings, 'minwfcrhoratio', read_config.get(self.config, 'ecutconv', 'minwfcrhoratio'))


    self.results = {'ecutwfc': 50, 'ecutrho': 100,}
    if(settings['run'] == False or settings['run'].lower()[0:1] == "f"):
      return None


    print("ECUT_CONV")

    plot_data = [[[],[]],[[],[]],[[],[]],[[],[]],[[],[]],[[],[]]]
    plot_data_conv = [[[],[]],[[],[]],[[],[]],[[],[]],[[],[]],[[],[]]]
    
    pdata = numpy.zeros((100, 1000),)
    pdata_count = numpy.zeros((100,),dtype=numpy.int32)
    

    # INCREASE ECUTWFC AND ECUTRHO     
    wfc_min = float(settings['wfc_min'])
    wfc_inc = float(settings['wfc_inc'])
    wfc_dec = float(settings['wfc_dec'])
    wfc_max = float(settings['wfc_max'])
    rho_dec = float(settings['rho_dec'])
    minwfcrhoratio = float(settings['minwfcrhoratio'])

    # CONVERT TO RY
    wfc_min = units.convert(self.units_energy,'RY',float(wfc_min))   # CONVERT INTO RY
    wfc_inc = units.convert(self.units_energy,'RY',float(wfc_inc))   # CONVERT INTO RY
    wfc_dec = units.convert(self.units_energy,'RY',float(wfc_dec))   # CONVERT INTO RY
    wfc_max = units.convert(self.units_energy,'RY',float(wfc_max))   # CONVERT INTO RY



    #################################
    # 1. INCREASE ECUTWFC 
    #################################
    converged = 0
    elast = None
    flast = None
    ecutwfc = wfc_min

    # CONVERGENCE FOR EACH CHANGE BY 1 RY 
    econv_threshold = self.econv_threshold_ry * wfc_inc
    fconv_threshold = self.fconv_threshold * wfc_inc

    ecutwfc_v = None
    counter = 0
    while((ecutwfc <= wfc_max and converged < 2) or counter < 5):
      counter = counter + 1

      # Get file name
      file_name = self.fname()

      # Make and save file
      ecutfile = pwscf_input()
      ecutfile.load(self.tfile)
      ecutfile.set_ecutwfc(ecutwfc)
      ecutfile.set_ecutrho(4 * ecutwfc)
      ecutfile.set_k_points(settings['kpointstype'], settings['kpointsval'])
      ecutfile.save(file_name, self.ecutd, 7)  

      # Run pwscf
      results = pwscf_exec.execute(file_name, self.ecutd, self.ecutd, self.runpwscfd, self.log)

      # Get energy and force
      energy = results[0]['energy_per_atom']
      force = results[0]['total_force']

      # Save plot data
      pdata[0,counter-1] = ecutwfc
      pdata[1,counter-1] = energy
      pdata[2,counter-1] = force

      if(elast != None):
        if(converged < 2):
          ecutwfc_v = ecutwfc - wfc_inc
        econv = abs(elast - energy)
        fconv = abs(flast - force)
        if(econv <= econv_threshold and fconv <= fconv_threshold):
          converged = converged + 1
        else:
          converged = 0

      # Save last energy
      elast = energy
      flast = force

      # Increment
      ecutwfc = ecutwfc + wfc_inc
      
    pdata_count[0:2] = counter
    pdata_count[3:5] = counter - 1
    for i in range(counter-1):
      pdata[3,i] = 0.5 * (pdata[0,i+1] + pdata[0,i])
      pdata[4,i] = abs(pdata[1,i+1] - pdata[1,i]) / wfc_inc
      pdata[5,i] = abs(pdata[2,i+1] - pdata[2,i]) / wfc_inc 
      
      

    print(ecutwfc_v, 4.0 * ecutwfc_v)
    self.log.add("################################") 
    self.log.add("Ecut Converge Stage 1") 
    self.log.add("steps = " + str(counter)) 
    self.log.add("ecutwfc = " + str(ecutwfc_v)) 
    self.log.add("ecutrho = " + str(4.0 * ecutwfc_v)) 
    self.log.add("################################") 
    self.log.add("") 
    


    #################################
    # 2. DECREASE ECUTRHO
    #################################
    ecutrho = 4.0 * ecutwfc_v
    converged = 0
    elast = None
    flast = None
    ecutrhomin = minwfcrhoratio * ecutwfc_v

    # CONVERGENCE FOR EACH CHANGE BY 1 RY 
    econv_threshold = self.econv_threshold_ry * rho_dec
    fconv_threshold = self.fconv_threshold * rho_dec

    counter = 0
    ecutrho_v = None
    while((ecutrho >= ecutwfc_v and converged < 2) or counter < 5):
      counter = counter + 1

      # Get file name
      file_name = self.fname()

      # Make and save file
      ecutfile = pwscf_input()
      ecutfile.load(self.tfile)
      ecutfile.set_ecutwfc(ecutwfc_v)
      ecutfile.set_ecutrho(ecutrho)
      ecutfile.set_k_points(settings['kpointstype'], settings['kpointsval'])
      ecutfile.save(file_name, self.ecutd, 7)  


      # Run pwscf
      results = pwscf_exec.execute(file_name, self.ecutd, self.ecutd, self.runpwscfd, self.log)

      # Get energy and force
      energy = results[0]['energy_per_atom']
      force = results[0]['total_force']

      # Save plot data
      pdata[6,counter-1] = ecutrho
      pdata[7,counter-1] = energy
      pdata[8,counter-1] = force

      # If this pwscf failed
      if(energy == None):
        if(ecutrho_v == None):
          ecutrho_v = ecutrho + rho_dec
        converged = 2
      # If less than equal to min, set min
      elif(ecutrho <= ecutrhomin):
        if(ecutrho_v == None):
          ecutrho_v = ecutrhomin
        converged = 2
      # Threshold met
      elif(elast != None):
        econv = abs(elast - energy)
        fconv = abs(flast - force)
        if(econv > econv_threshold or fconv > fconv_threshold):
          converged = 2
          if(ecutrho_v == None):
            ecutrho_v = ecutrho + rho_dec
      

      # Save last energy
      elast = energy
      flast = force

      # Decrease
      ecutrho = ecutrho - rho_dec
      
      
    pdata_count[6:8] = counter
    pdata_count[9:11] = counter - 1
    for i in range(counter-1):
      pdata[9,i] = 0.5 * (pdata[6,i+1] + pdata[6,i])
      pdata[10,i] = abs(pdata[7,i+1] - pdata[7,i]) / rho_dec
      pdata[11,i] = abs(pdata[8,i+1] - pdata[8,i]) / rho_dec
      

    print(ecutwfc_v, ecutrho_v)
    self.log.add("################################") 
    self.log.add("Ecut Converge Stage 2") 
    self.log.add("steps = " + str(counter)) 
    self.log.add("ecutwfc = " + str(ecutwfc_v)) 
    self.log.add("ecutrho = " + str(ecutrho_v)) 
    self.log.add("################################") 
    self.log.add("") 

    # CONVERGENCE FOR EACH CHANGE BY 1 RY 
    econv_threshold = self.econv_threshold_ry * wfc_dec
    fconv_threshold = self.fconv_threshold * wfc_dec



    #################################
    # 3. DECREASE ECUTWFC
    #################################
    ecutwfc = ecutwfc_v
    ecutrho = ecutrho_v
    converged = 0
    elast = None
    flast = None

    ecutwfc_v = None
    counter = 0
    while((ecutwfc >= 0 and converged < 2) or counter < 5):
      counter = counter + 1
      # Get file name
      file_name = self.fname()

      # Make and save file
      ecutfile = pwscf_input()
      ecutfile.load(self.tfile)
      ecutfile.set_ecutwfc(ecutwfc)
      ecutfile.set_ecutrho(ecutrho)
      ecutfile.set_k_points(settings['kpointstype'], settings['kpointsval'])
      ecutfile.save(file_name, self.ecutd, 7)  

      # Run pwscf
      results = pwscf_exec.execute(file_name, self.ecutd, self.ecutd, self.runpwscfd, self.log)

      # Get energy and force
      energy = results[0]['energy_per_atom']
      force = results[0]['total_force']

      # Save plot data
      pdata[12,counter-1] = ecutwfc
      pdata[13,counter-1] = energy
      pdata[14,counter-1] = force
      
      if(energy == None):
        converged = 2
        if(ecutwfc_v == None):
          ecutwfc_v = ecutwfc + wfc_dec

      elif(elast != None):
        econv = abs(elast - energy)
        fconv = abs(flast - force)
        if(econv > econv_threshold or fconv > fconv_threshold):
          converged = 2 
          if(ecutwfc_v == None):
            ecutwfc_v = ecutwfc + wfc_dec  

      # Save last energy
      elast = energy
      flast = force

      # Decrease
      ecutwfc = ecutwfc - wfc_dec

      
    pdata_count[12:14] = counter
    pdata_count[15:17] = counter - 1
    for i in range(counter-1):
      pdata[15,i] = 0.5 * (pdata[12,i+1] + pdata[12,i])
      pdata[16,i] = abs(pdata[13,i+1] - pdata[13,i]) / wfc_dec
      pdata[17,i] = abs(pdata[14,i+1] - pdata[14,i]) / wfc_dec

    print(ecutwfc_v, ecutrho_v)
    self.log.add("################################") 
    self.log.add("Ecut Converge Stage 3") 
    self.log.add("steps = " + str(counter)) 
    self.log.add("ecutwfc = " + str(ecutwfc_v)) 
    self.log.add("ecutrho = " + str(ecutrho_v)) 
    self.log.add("################################") 
    self.log.add("") 

    p_list = [1,2,4,5,7,8,10,11,13,14,16,17]
    for p in p_list:
      p_min = numpy.amin(pdata[p,:])
      pdata[p,:] = pdata[p,:] - p_min
    
    
    # PLOT 1

    plt.clf()
    
    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')

    fig, axs = plt.subplots(3, 2, figsize=(12,9))
    fig.tight_layout(pad=5.0)
    
    fig.suptitle('Convergence - Relative Values')
    
    c = pdata_count[0]
    axs[0, 0].plot(pdata[0,1:c], pdata[1,1:c], color='k', ls='solid')
    axs[0, 0].set_title('Increasing Ecutwfc')
    axs[0, 0].set_xlabel('Ecutwfc / eV')
    axs[0, 0].set_ylabel('Energy / eV')

    axs[0, 1].plot(pdata[0,1:c], pdata[2,1:c], color='k', ls='solid')
    axs[0, 1].set_title('Increasing Ecutwfc')
    axs[0, 1].set_xlabel('Ecutwfc / eV')
    axs[0, 1].set_ylabel('Force / RY/BOHR')
    
    c = pdata_count[6]
    axs[1, 0].plot(pdata[6,1:c], pdata[7,1:c], color='k', ls='solid')
    axs[1, 0].set_title('Decreasing Ecutrho')
    axs[1, 0].set_xlabel('Ecutrho / bohr-1')
    axs[1, 0].set_ylabel('Energy / eV')
    
    axs[1, 1].plot(pdata[6,1:c], pdata[8,1:c], color='k', ls='solid')
    axs[1, 1].set_title('Decreasing Ecutrho')
    axs[1, 1].set_xlabel('Ecutrho / bohr-1')
    axs[1, 1].set_ylabel('Force / RY/BOHR')
    
    c = pdata_count[12]
    axs[2, 0].plot(pdata[12,1:c], pdata[13,1:c], color='k', ls='solid')
    axs[2, 0].set_title('Increasing Ecutwfc')
    axs[2, 0].set_xlabel('Ecutwfc / eV')
    axs[2, 0].set_ylabel('Energy / eV')
    
    axs[2, 1].plot(pdata[12,1:c], pdata[14,1:c], color='k', ls='solid')
    axs[2, 1].set_title('Increasing Ecutwfc')
    axs[2, 1].set_xlabel('Ecutwfc / eV')
    axs[2, 1].set_ylabel('Force / RY/BOHR')

    plt.savefig(self.plotsd + '/' + 'convergence.svg')
    
       

    # PLOT 2

    plt.clf()
    
    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')

    fig, axs = plt.subplots(3, 2, figsize=(12,9))
    fig.tight_layout(pad=5.0)
    
    fig.suptitle('Convergence - Differences')
    
    c = pdata_count[3]
    axs[0, 0].plot(pdata[3,1:c], pdata[4,1:c], color='k', ls='solid')
    axs[0, 0].set_title('Increasing Ecutwfc')
    axs[0, 0].set_xlabel('Ecutwfc / eV')
    axs[0, 0].set_ylabel('|Energy| / eV')

    axs[0, 1].plot(pdata[3,1:c], pdata[5,1:c], color='k', ls='solid')
    axs[0, 1].set_title('Increasing Ecutwfc')
    axs[0, 1].set_xlabel('Ecutwfc / eV')
    axs[0, 1].set_ylabel('|Force| / RY/BOHR')
    
    c = pdata_count[9]
    axs[1, 0].plot(pdata[9,1:c], pdata[10,1:c], color='k', ls='solid')
    axs[1, 0].set_title('Decreasing Ecutrho')
    axs[1, 0].set_xlabel('Ecutrho / bohr-1')
    axs[1, 0].set_ylabel('|Energy| / eV')
    
    axs[1, 1].plot(pdata[9,1:c], pdata[11,1:c], color='k', ls='solid')
    axs[1, 1].set_title('Decreasing Ecutrho')
    axs[1, 1].set_xlabel('Ecutrho / bohr-1')
    axs[1, 1].set_ylabel('|Force| / RY/BOHR')
    
    c = pdata_count[15]
    axs[2, 0].plot(pdata[15,1:c], pdata[16,1:c], color='k', ls='solid')
    axs[2, 0].set_title('Increasing Ecutwfc')
    axs[2, 0].set_xlabel('Ecutwfc / eV')
    axs[2, 0].set_ylabel('|Energy| / eV')
    
    axs[2, 1].plot(pdata[15,1:c], pdata[17,1:c], color='k', ls='solid')
    axs[2, 1].set_title('Increasing Ecutwfc')
    axs[2, 1].set_xlabel('Ecutwfc / eV')
    axs[2, 1].set_ylabel('|Force| / RY/BOHR')

    plt.savefig(self.plotsd + '/' + 'convergence_differences.svg')


    # Store results
    self.results = {'ecutwfc': ecutwfc_v, 'ecutrho': ecutrho_v,}


    fh = open(self.wd + "/ECUT_WFC_RHO.txt", "w")
    fh.write('ecutwfc: ' + str(ecutwfc_v) + '\n')
    fh.write('ecutrho: ' + str(ecutrho_v) + '\n')
    fh.close()



  def ecut_2d(self):  


    # LOAD SETTINGS FROM INPUT FILE

    settings = {'run': True, 'ecutwfcmin': 20, 'ecutwfcmax': 100, 'ecutwfcinc': 5, 'ecutrhomin': 80, 'ecutrhomax': 400, 'ecutrhoinc': 20}

    settings = self.add_to_dict(settings, 'run', read_config.get(self.config, 'ecut2d', 'run'))
    settings = self.add_to_dict(settings, 'ecutwfcmin', read_config.get(self.config, 'ecut2d', 'wfc1'))
    settings = self.add_to_dict(settings, 'ecutwfcmax', read_config.get(self.config, 'ecut2d', 'wfc2'))
    settings = self.add_to_dict(settings, 'ecutwfcinc', read_config.get(self.config, 'ecut2d', 'wfc3'))
    settings = self.add_to_dict(settings, 'ecutrhomin', read_config.get(self.config, 'ecut2d', 'rho1'))
    settings = self.add_to_dict(settings, 'ecutrhomax', read_config.get(self.config, 'ecut2d', 'rho2'))
    settings = self.add_to_dict(settings, 'ecutrhoinc', read_config.get(self.config, 'ecut2d', 'rho3'))

    if(settings['run'] == False or settings['run'].lower()[0:1] == "f"):
      return None

    print("ECUT_2D")  


    # INCREASE ECUTWFC AND ECUTRHO 
    
    ecutwfc = float(settings['ecutwfcmin'])
    wfcinc = float(settings['ecutwfcinc'])
    wfcend = float(settings['ecutwfcmax'])

    # CONVERT TO RY
    
    
    ecutrho = float(settings['ecutrhomin'])
    rhoinc = float(settings['ecutrhoinc'])
    rhoend = float(settings['ecutrhomax'])


    e = []   
    f = []  

    # Loop through energies
    converged = False

    i = 0
    ecut_rhos = []
    while(ecutrho <= rhoend):
      e.append([])
      f.append([])

      ecut_rhos.append(ecutrho)
      if(i == 0):        
        ecut_wfcs = []

      ecutwfc = float(settings['ecutwfcmin'])
      while(ecutwfc <= wfcend):
        if(i == 0):        
          ecut_wfcs.append(ecutwfc)
 
        # Get file name
        file_name = self.fname()

        # Make and save file
        ecutfile = pwscf_input()
        ecutfile.load(self.tfile)
        ecutfile.set_ecutwfc(ecutwfc)
        ecutfile.set_ecutrho(ecutrho)
        ecutfile.save(file_name, self.ecutd, 7)  

        # Run pwscf
        results = pwscf_exec.execute(file_name, self.ecutd, self.ecutd, self.runpwscfd, self.log)

        energy = results[0]['energy_per_atom']
        force = results[0]['total_force']

        e[i].append(energy)
        f[i].append(force)

        # Incrementc
        ecutwfc = ecutwfc + wfcinc

      # Increment
      i = i + 1
      ecutrho = ecutrho + rhoinc

    # Write to csv
    fh = open(self.csvd + "/ecut_energy.csv", "w")
    fd = ""
    for j in range(len(ecut_wfcs)):
      fd = fd + "," + str(ecut_wfcs[j])
    fd = fd + "\n"
    for i in range(len(ecut_rhos)):
      fd = fd + str(ecut_rhos[i])
      for j in range(len(ecut_wfcs)):
        fd = fd + "," + str(e[i][j])
      fd = fd + "\n"
    fh.write(fd)
    fh.close()  

    # Write to csv
    fh = open(self.csvd + "/ecut_force.csv", "w")
    fd = ""
    for j in range(len(ecut_wfcs)):
      fd = fd + "," + str(ecut_wfcs[j])
    fd = fd + "\n"
    for i in range(len(ecut_rhos)):
      fd = fd + str(ecut_rhos[i])
      for j in range(len(ecut_wfcs)):
        fd = fd + "," + str(f[i][j])
      fd = fd + "\n"
    fh.write(fd)
    fh.close()  
    

    cE = []
    for i in range(len(e)):
      cE.append([])
      for j in range(len(e[0])):
        if(e[i][j] == None):
          cE[i].append(None)
        else:
          cE[i].append(self.avg(e[i][j], e[i][:]))

    cF = []
    for i in range(len(f)):
      cF.append([])
      for j in range(len(f[0])):
        if(f[i][j] == None):
          cF[i].append(None)
        else:
          cF[i].append(self.avg(f[i][j], f[i][:]))

    
    X = []
    ecutwfc = float(settings['ecutwfcmin'])
    while(ecutwfc <= wfcend):
      X.append(ecutwfc)
      ecutwfc = ecutwfc + wfcinc
    
    Y = []
    ecutrho = float(settings['ecutrhomin'])
    while(ecutrho <= rhoend):
      Y.append(ecutrho)
      ecutrho = ecutrho + rhoinc
 

    # Set zero as minimum value
    e_min = None
    for i in range(len(cE)):
      for j in range(len(cE[i])):
        if(cE[i][j] != None and (e_min == None or e_min > cE[i][j])):
          e_min = cE[i][j]
    for i in range(len(cE)):
      for j in range(len(cE[i])):
        if(cE[i][j] != None):
          cE[i][j] = cE[i][j] - e_min
 

    # Set zero as minimum value
    f_min = None
    for i in range(len(cF)):
      for j in range(len(cF[i])):
        if(cF[i][j] != None and (f_min == None or f_min > cF[i][j])):
          f_min = cE[i][j]
    for i in range(len(cF)):
      for j in range(len(cF[i])):
        if(cF[i][j] != None):
          cF[i][j] = cF[i][j] - f_min


    # Plot
    plt.clf()    
    plt.contourf(X, Y, cE, 20, cmap='Greys')
    plt.colorbar()
    plt.savefig(self.plotsd + '/' + 'ecut_energy.eps', format='eps')

    plt.clf()
    plt.contourf(X, Y, cF, 20, cmap='Greys')
    plt.colorbar()
    plt.savefig(self.plotsd + '/' + 'ecut_force.eps', format='eps')



    



  def kpoints_2d(self):    

    settings = {'run': True, }

    settings = self.add_to_dict(settings, 'run', read_config.get(self.config, 'kconv', 'run'))
    settings = self.add_to_dict(settings, 'k_min', read_config.get(self.config, 'kconv', 'k_min'))
    settings = self.add_to_dict(settings, 'k_inc', read_config.get(self.config, 'kconv', 'k_inc'))
    settings = self.add_to_dict(settings, 'k_max', read_config.get(self.config, 'kconv', 'k_max'))
    settings = self.add_to_dict(settings, 'k_type', read_config.get(self.config, 'kconv', 'k_type'))
    smear_widths = self.config['kconv'][0]['smearing']

    if(settings['run'] == False or settings['run'].lower()[0:1] == "f"):
      return None

    print("KPOINTS_2D")  

    e = []
    f = []

    k_start = int(settings['k_min'])
    k_inc = int(settings['k_inc'])
    k_end = int(settings['k_max'])
    k_type = settings['k_type']

    i = 0
    k_points = []
    k = k_start
    while(k <= k_end):
      e.append([])
      f.append([])

      k_points.append(k)
      if(i == 0):        
        smear_points = []

      for s in smear_widths:
        if(i == 0):        
          smear_points.append(s)

        # Get file name
        file_name = self.fname()

        # Make and save file
        ecutfile = pwscf_input()
        ecutfile.load(self.tkfile)
        ecutfile.set_ecutwfc(self.results['ecutwfc'])
        ecutfile.set_ecutrho(self.results['ecutrho'])
        ecutfile.set_degauss(s)
        ecutfile.set_k_points(k_type, [k,k,k,1,1,1])
        ecutfile.save(file_name, self.kpointsd, 7)  

        # Run pwscf
        results = pwscf_exec.execute(file_name, self.kpointsd, self.kpointsd, self.runpwscfd, self.log)

        if(results != None):
          energy = results[0]['energy_per_atom']
          force = results[0]['total_force']
        else:
          energy = None
          force = None

        e[i].append(energy)
        f[i].append(force)

      # Increment k
      i = i + 1
      k = k + k_inc


    # Write to csv
    fh = open(self.csvd + "/kpoints_energy.csv", "w")
    fd = ""
    for j in range(len(smear_points)):
      fd = fd + "," + str(smear_points[j])
    fd = fd + "\n"
    for i in range(len(k_points)):
      fd = fd + str(smear_points[i])
      for j in range(len(smear_points)):
        fd = fd + "," + str(e[i][j])
      fd = fd + "\n"
    fh.write(fd)
    fh.close() 


    # Write to csv
    fh = open(self.csvd + "/kpoints_force.csv", "w")
    fd = ""
    for j in range(len(smear_points)):
      fd = fd + "," + str(smear_points[j])
    fd = fd + "\n"
    for i in range(len(k_points)):
      fd = fd + str(smear_points[i])
      for j in range(len(smear_points)):
        fd = fd + "," + str(f[i][j])
      fd = fd + "\n"
    fh.write(fd)
    fh.close() 


    cE = []
    for i in range(len(e)):
      cE.append([])
      for j in range(len(e[0])):
        if(e[i][j] == None):
          cE[i].append(None)
        else:
          cE[i].append(self.avg(e[i][j], e[i][:]))

    cF = []
    for i in range(len(f)):
      cF.append([])
      for j in range(len(f[0])):
        if(f[i][j] == None):
          cF[i].append(None)
        else:
          cF[i].append(self.avg(f[i][j], f[i][:]))
    

    X = []
    for s in smear_widths:
      X.append(float(s))
    
    Y = []
    k = k_start
    while(k <= k_end):
      Y.append(k)
      k = k + k_inc
    

    # Set zero as minimum value
    e_min = None
    for i in range(len(cE)):
      for j in range(len(cE[i])):
        if(cE[i][j] != None and (e_min == None or e_min > cE[i][j])):
          e_min = cE[i][j]
    for i in range(len(cE)):
      for j in range(len(cE[i])):
        if(cE[i][j] != None):
          cE[i][j] = cE[i][j] - e_min
 

    # Set zero as minimum value
    f_min = None
    for i in range(len(cF)):
      for j in range(len(cF[i])):
        if(cF[i][j] != None and (f_min == None or f_min > cF[i][j])):
          f_min = cE[i][j]
    for i in range(len(cF)):
      for j in range(len(cF[i])):
        if(cF[i][j] != None):
          cF[i][j] = cF[i][j] - f_min



    plt.clf()
    
    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')

    plt.figure(figsize=(6,4))
 
    #plt.contourf(X, Y, e)
    #plt.legend(fancybox=True, framealpha=1, shadow=True, borderpad=1)
    plt.contourf(X, Y, cE, 20, cmap='Greys')
    plt.colorbar()
    plt.savefig(self.plotsd + '/' + 'k_points_energy.eps', format='eps')

    plt.clf()
    #plt.contourf(X, Y, f)
    #plt.legend(fancybox=True, framealpha=1, shadow=True, borderpad=1)
    plt.contourf(X, Y, cF, 20, cmap='Greys')
    plt.colorbar()
    plt.savefig(self.plotsd + '/' + 'k_points_force.eps', format='eps')
    

  def avg(self, item, itemlist):
    n = 0
    nsum = 0.0
    for nitem in itemlist:
      if(nitem != None):
        n = n + 1
        nsum = nsum + nitem
    return item - nsum / n

  def load_file(self):

    if(len(sys.argv) < 2):
      return False
    if(not os.path.exists(sys.argv[1])):
      return False
    self.log.add("Config file: " + sys.argv[1])
    self.config = read_config.read_file(sys.argv[1])
    
    self.units_len = self.set_value(read_config.get(self.config, 'units', 'length'), "BOHR")
    self.units_energy = self.set_value(read_config.get(self.config, 'units', 'energy'), "RY")
    
    self.econv_threshold = float(self.set_value(read_config.get(self.config, 'converge', 'energy'), 1.0E-6))
    self.fconv_threshold = float(self.set_value(read_config.get(self.config, 'converge', 'force'), 1.0E-5))


    # THRESHOLDS IN RY
    self.econv_threshold_ry = units.convert(self.units_energy,'RY',float(self.econv_threshold))   # CONVERT INTO RY

  
   

    self.log.add("Length units: " + str(self.units_len))
    self.log.add("Energy units: " + str(self.units_energy))
    self.log.add("Energy convergence: " + str(self.econv_threshold) + " " + str(self.units_energy))
    self.log.add("Force convergence: " + str(self.fconv_threshold))
    self.log.add("Energy convergence: " + str(self.econv_threshold_ry) + " RY")









  def exit(self):
    # Output log
    self.log.output()
    print("END")
    # Exit Program
    sys.exit()


  def fname(self):
    self.fname_count = self.fname_count + 1    
    return "file_" + str(self.fname_count) + ".in"

  def add_to_dict(self, dict_in, key, value, default_value=None):
    if(value != None):
      dict_in[key] = value
    elif(default_value != None):
      dict_in[key] = default_value
    return dict_in

  def set_value(self, inp, inp_d=None):
    if(inp==None):
      inp = inp_d
    if(inp!=None):
      inp = inp.upper()
    return inp






