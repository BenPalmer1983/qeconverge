################################################################
#    Processing PWscf input file
#
#
#
#
################################################################


#!/bin/python3
########################################################################
import os
import datetime
import re
import sys
import time
import numpy
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits import mplot3d
import hashlib
import random
import math
from shutil import copyfile

###########################################
#  CLASS pwscf_converg
###########################################
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
    axs[0, 0].set_xlabel('Ecutwfc /eV')
    axs[0, 0].set_ylabel('|Energy| / eV')

    axs[0, 1].plot(pdata[0,1:c], pdata[2,1:c], color='k', ls='solid')
    axs[0, 1].set_title('Increasing Ecutwfc')
    axs[0, 1].set_xlabel('Ecutwfc /eV')
    axs[0, 1].set_ylabel('|Energy| / eV')
    
    c = pdata_count[6]
    axs[1, 0].plot(pdata[6,1:c], pdata[7,1:c], color='k', ls='solid')
    axs[1, 0].set_title('Decreasing Ecutrho')
    axs[1, 0].set_xlabel('Ecutrho / bohr-1')
    axs[1, 0].set_ylabel('|Energy| / eV')
    
    axs[1, 1].plot(pdata[6,1:c], pdata[8,1:c], color='k', ls='solid')
    axs[1, 1].set_title('Decreasing Ecutrho')
    axs[1, 1].set_xlabel('Ecutrho / bohr-1')
    axs[1, 1].set_ylabel('|Energy| / eV')
    
    c = pdata_count[12]
    axs[2, 0].plot(pdata[12,1:c], pdata[13,1:c], color='k', ls='solid')
    axs[2, 0].set_title('Increasing Ecutwfc')
    axs[2, 0].set_xlabel('Ecutwfc /eV')
    axs[2, 0].set_ylabel('|Energy| / eV')
    
    axs[2, 1].plot(pdata[12,1:c], pdata[14,1:c], color='k', ls='solid')
    axs[2, 1].set_title('Increasing Ecutwfc')
    axs[2, 1].set_xlabel('Ecutwfc /eV')
    axs[2, 1].set_ylabel('|Energy| / eV')

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
    axs[0, 0].set_xlabel('Ecutwfc /eV')
    axs[0, 0].set_ylabel('|Energy| / eV')

    axs[0, 1].plot(pdata[3,1:c], pdata[5,1:c], color='k', ls='solid')
    axs[0, 1].set_title('Increasing Ecutwfc')
    axs[0, 1].set_xlabel('Ecutwfc /eV')
    axs[0, 1].set_ylabel('|Energy| / eV')
    
    c = pdata_count[9]
    axs[1, 0].plot(pdata[9,1:c], pdata[10,1:c], color='k', ls='solid')
    axs[1, 0].set_title('Decreasing Ecutrho')
    axs[1, 0].set_xlabel('Ecutrho / bohr-1')
    axs[1, 0].set_ylabel('|Energy| / eV')
    
    axs[1, 1].plot(pdata[9,1:c], pdata[11,1:c], color='k', ls='solid')
    axs[1, 1].set_title('Decreasing Ecutrho')
    axs[1, 1].set_xlabel('Ecutrho / bohr-1')
    axs[1, 1].set_ylabel('|Energy| / eV')
    
    c = pdata_count[15]
    axs[2, 0].plot(pdata[15,1:c], pdata[16,1:c], color='k', ls='solid')
    axs[2, 0].set_title('Increasing Ecutwfc')
    axs[2, 0].set_xlabel('Ecutwfc /eV')
    axs[2, 0].set_ylabel('|Energy| / eV')
    
    axs[2, 1].plot(pdata[15,1:c], pdata[17,1:c], color='k', ls='solid')
    axs[2, 1].set_title('Increasing Ecutwfc')
    axs[2, 1].set_xlabel('Ecutwfc /eV')
    axs[2, 1].set_ylabel('|Energy| / eV')

    plt.savefig(self.plotsd + '/' + 'convergence_differences.svg')

# Store results
    self.results = {'ecutwfc': ecutwfc_v, 'ecutrho': ecutrho_v,}

    fh = open(self.wd + "/ECUT_WFC_RHO.txt", "w")
    fh.write('ecutwfc: ' + str(ecutwfc_v))
    fh.write('ecutrho: ' + str(ecutrho_v))
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

###########################################
#  CLASS pwscf_inpu
###########################################
class pwscf_input:

  def __init__(self, file_name=None, file_dir=None):
    self.file_data = []
    self.file_name = None
    self.dir_name = None
    self.reset()
    self.defaults()
    self.rand = rand_dist()
    self.rand.gheat()
    self.rand_seed_set = 0
    if(file_name != None):
      self.load(file_name, file_dir)

  def reset(self):

# Control
    self.control = {
      "calculation": None,
      "title": None,
      "verbosity": None,
      "restart_mode": None,
      "wf_collect": None,
      "nstep": None,
      "iprint": None,
      "tstress": None,
      "tprnfor": None,
      "dt": None,
      "outdir": None,
      "wfcdir": None,
      "prefix": None,
      "lkpoint_dir": None,
      "max_seconds": None,
      "etot_conv_thr": None,
      "forc_conv_thr": None,
      "disk_io": None,
      "pseudo_dir": None,
      "tefield": None,
      "dipfield": None,
      "lelfield": None,
      "nberrycyc": None,
      "lorbm": None,
      "lberry": None,
      "gdir": None,
      "nppstr": None,
      "lfcpopt": None,
      "gate": None
    }

# SYSTEM
    self.system = {
      "ibrav": None,
      "celldm": None,
      "A": None,
      "B": None,
      "C": None,
      "cosAB": None,
      "cosAC": None,
      "cosBC": None,
      "nat": None,
      "ntyp": None,
      "nbnd": None,
      "tot_charge": None,
      "starting_charge": None,
      "tot_magnetization": None,
      "starting_magnetization": None,
      "ecutwfc": None,
      "ecutrho": None,
      "ecutfock": None,
      "nr1": None,
      "nr2": None,
      "nr3": None,
      "nr1s": None,
      "nr2s": None,
      "nr3s": None,
      "nosym": None,
      "nosym_evc": None,
      "noinv": None,
      "no_t_rev": None,
      "force_symmorphic": None,
      "use_all_frac": None,
      "occupations": None,
      "one_atom_occupations": None,
      "starting_spin_angle": None,
      "degauss": None,
      "smearing": None,
      "nspin": None,
      "noncolin": None,
      "ecfixed": None,
      "qcutz": None,
      "q2sigma": None,
      "input_dft": None,
      "exx_fraction": None,
      "screening_parameter": None,
      "exxdiv_treatment": None,
      "x_gamma_extrapolation": None,
      "ecutvcut": None,
      "nqx1": None,
      "nqx2": None,
      "nqx3": None,
      "lda_plus_u": None,
      "lda_plus_u_kind": None,
      
      "xdm": None,
      "xdm_a1": None,
      "xdm_a2": None,
      "space_group": None,
      "uniqueb": None,
      "origin_choice": None,
      "rhombohedral": None,
      "zgate": None,
      "relaxz": None,
      "block": None,
      "block_1": None,
      "block_2": None,
      "block_height": None
    }
    
# ELECTRONS
    self.electrons = {
      "electron_maxstep": None,
      "scf_must_converge": None,
      "conv_thr": None,
      "adaptive_thr": None,
      "conv_thr_init": None,
      "conv_thr_multi": None,
      "mixing_mode": None,
      "mixing_beta": None,
      "mixing_ndim": None,
      "mixing_fixed_ns": None,
      "diagonalization": None,
      "ortho_para": None,
      "diago_thr_init": None,
      "diago_cg_maxiter": None,
      "diago_david_ndim": None,
      "diago_full_acc": None,
      "efield": None,
      "efield_cart": None,
      "efield_phase": None,
      "startingpot": None,
      "startingwfc": None,
      "tqr": None
    }

# IONS
    self.ions = {
      "ion_dynamics": None,
      "ion_positions": None,
      "pot_extrapolation": None,
      "wfc_extrapolation": None,
      "remove_rigid_rot": None,
      "ion_temperature": None,
      "tempw": None,
      "tolp": None,
      "delta_t": None,
      "nraise": None,
      "refold_pos": None,
      "upscale": None,
      "bfgs_ndim": None,
      "trust_radius_max": None,
      "trust_radius_min": None,
      "trust_radius_ini": None,
      "w_1": None,
      "w_2": None
    }

# CELL
    self.cell = {
      "cell_dynamics": None,
      "press": None,
      "wmass": None,
      "cell_factor": None,
      "press_conv_thr": None,
      "cell_dofree": None
    }

# Other lists
    self.atomic_species = []
    self.atomic_positions = []
    self.k_points = []
    self.cell_parameters = []

# File
    self.file = ""

  def defaults(self):
      
    try:
      self.scratch_dir = os.environ['PWSCF_SCRATCH']
    except:
      self.scratch_dir = '/opt/scratch'

    try:
      self.pp_dir = os.environ['PWSCF_PP']
    except:
      self.pp_dir = '/opt/pp'

#  Load data from file
  def load(self, file_name, file_dir=None): 
    self.file_name = file_name
    self.dir_name = file_dir
    if(file_dir != None):
      self.file_path = file_dir + "/" + file_name
    else:  
      self.file_path = file_name
    data = self.load_from_file(self.file_path)
    self.load_data(data)
    
# Load from a block of data (text, file etc)
  def load_data(self, data):  
# Store data into file_data list
    self.file_data.append(data)
    
# Reset data store
    self.reset()
    
# Clean
    data = pwscf_input.clean(data)

# split
    data = data.split("\n")    
    
# Load keywords
###################################
    keywords = []
# Load Keywords
    for line in data:
      line = line.strip()
      if(len(line)>0):
# Remove trailing comma
        if(line[-1] == ","):
          line = line[0:-1]
        fields = line.split("=")
        if(len(fields) == 2):
          field_lc = fields[0].lower()
          keyword, id = pwscf_input.process_keyword(field_lc)          
          pwscf_input.add_keyword(keywords, keyword, id, fields[1])
      
    for pair in keywords:
      if(pair[0] in self.control):
        self.control[pair[0]] = pair[1]
      elif(pair[0] in self.system):
        self.system[pair[0]] = pair[1]
      elif(pair[0] in self.electrons):
        self.electrons[pair[0]] = pair[1]
      elif(pair[0] in self.ions):
        self.ions[pair[0]] = pair[1]
      elif(pair[0] in self.cell):
        self.cell[pair[0]] = pair[1]        

# Load atomic species
###################################
    n_species = 0    
    if(self.system['ntyp'] != None):
      try:
        n_species = int(self.system['ntyp'])
      except:  
        n_species = 0        
    if(n_species > 0):
      counter = 0
      for line in data:
        line = line.strip()
        if(line.upper()[0:14] == "ATOMIC_SPECIES"):
          counter = counter + 1
        elif(counter > 0 and counter <= n_species and line != ""):
          counter = counter + 1
          self.atomic_species.append(pwscf_input.fields(line))

# Load atomic positions
###################################
    n_atoms = 0    
    if(self.system['nat'] != None):
      try:
        n_atoms = int(self.system['nat'])
      except:  
        n_atoms = 0        
    if(n_atoms > 0):
      counter = 0
      for line in data:
        line = line.strip()
        if(line.upper()[0:16] == "ATOMIC_POSITIONS"):
          fields = pwscf_input.fields(line)
          if(len(fields) == 2):
            self.atomic_positions.append(fields[1])
          counter = counter + 1
        elif(counter > 0 and counter <= n_atoms and line != ""):
          counter = counter + 1
          self.atomic_positions.append(pwscf_input.fields(line))

# k_points
###################################
    flag = 0
    for line in data:
      line = line.strip()
      if(line.upper()[0:8] == "K_POINTS"):           
        fields = pwscf_input.fields(line)
        k_points_type = fields[1]
        self.k_points.append(k_points_type)
        if(k_points_type.upper() == "AUTOMATIC"):
          flag = 1
      elif(flag > 0):
        flag = flag - 1        
        fields = pwscf_input.fields(line)
        self.k_points.append(fields)
        
# cell parameters
###################################
    flag = 0
    for line in data:
      line = line.strip()
      if(line.upper()[0:15] == "CELL_PARAMETERS"): 
        fields = pwscf_input.fields(line)
        self.cell_parameters.append(fields[1])

        flag = 3
      elif(flag>0):        
        fields = pwscf_input.fields(line)
        self.cell_parameters.append(fields)

    self.make()

#  Run as it's own program
  def run(self):
    self.reset()

    option = ""
    file_name = ""

    if(len(sys.argv) > 1 and sys.argv[1] is not None):
      option = sys.argv[1]

    if(len(sys.argv) > 2 and sys.argv[2] is not None):
      file_name = sys.argv[2]

    if(option.lower().strip() == "" or option.lower().strip() == "interactive"):
      self.menu()
      exit()
    elif(option.lower().strip() == "quiet"):
      print("Quiet")
    else:
      return 0

#################################
# READ/LOAD input file
#################################

  def load_from_file(self, file_name):
# Init variable
    file_data = ""

# Read it in line by line
    fh = open(file_name, "r")
    for file_row in fh:
      file_data = file_data + file_row.strip() + '\n'

    return file_data

#################################
# MAKE input file
#################################

  def make(self, rounding=None):
  
    now = datetime.datetime.now()
    time_now = str(now.hour) + ":" + str(now.minute) + "   " + str(now.day) + "/" + str(now.month) + "/" + str(now.year)
    file = "! Edited " + time_now + "\n"
    
# CONTROL
    file += "&CONTROL \n"
    for key in sorted(self.control.keys()):
      file += pwscf_input.make_line(key, self.control[key])
    file += "/ \n"
    
# SYSTEM
    file += "&SYSTEM \n"
    for key in sorted(self.system.keys()):
      file += pwscf_input.make_line(key, self.system[key]) 
    file += "/ \n"

# ELECTRONS
    file += "&ELECTRONS \n"
    for key in sorted(self.electrons.keys()):
      value = self.electrons[key]
      if(value != None):
        file += key + " = " + value + ", \n"      
    file += "/ \n"

# IONS
    file += "&IONS \n"
    for key in sorted(self.ions.keys()):
      value = self.ions[key]
      if(value != None):
        file += key + " = " + value + ", \n"      
    file += "/ \n"

# CELL
    file += "&CELL \n"
    for key in sorted(self.cell.keys()):
      value = self.cell[key]
      if(value != None):
        file += key + " = " + value + ", \n"      
    file += "/ \n"

# ATOMIC_SPECIES
    file += "ATOMIC_SPECIES \n"
    for species in self.atomic_species:      
      for field in species:
        file += field + " "
      file += "\n"

# ATOMIC_POSITIONS
    header = 0
    for position in self.atomic_positions:      
      if(header == 0):
        file += "ATOMIC_POSITIONS "
        file += position + "\n"
        header = 1
#elif(header == 1):
#  file += position[1] + "\n"
#  header = 2
      elif(header == 1):  
        j = 0
        for field in position:
          j = j + 1
          if(rounding == None or j == 1):
            file += field + "   "            
          else:
            file += str(round(float(field), int(rounding))) + "   "
        file += "\n"

# K_POINTS
    file += "K_POINTS " + self.k_points[0]
    file += "\n"
    for i in range(1,len(self.k_points)):
      for point in self.k_points[i]:
        file += point + " "
      file += "\n"
        
# K_POINTS
    file += "CELL_PARAMETERS " + self.cell_parameters[0]
    file += "\n"
    for i in range(1,len(self.cell_parameters)):
      for point in self.cell_parameters[i]:
        file += point + " "
      file += "\n"
      
# Process
    file = file.strip()
    
# Store data into file_data list
    self.file_data.append(file)

  def print_out(self):
    self.make()
    print(self.file_data[-1])
    
  def print_history(self):
    for file in self.file_data:
      print(file)
      print()
    
############################
#  Save, Save/Load Original
############################
    
  def save(self, file=None, dir=None, rounding=None):
# Build latest version of file
    self.make(rounding)
    
    if(file == None):
      file = self.file_name
    if(dir == None):  
      dir = self.dir_name
      
    self.file_name = file
    self.dir_name = dir
    
    if(dir == None):
      path = file
    else:  
      if (not os.path.exists(dir)):
        os.makedirs(dir)
      path = dir + "/" + file
      
# Write latest data
    fh = open(path, "w")
    fh.write(self.file_data[-1])
    fh.close()
    
  def save_original(self, file=None, dir=None):
    if(file == None):
      file = self.file_name
    if(dir == None):  
      dir = self.dir_name
      
    self.file_name = file
    self.dir_name = dir
    
    if(dir == None):
      path = file
    else:  
      if (not os.path.exists(dir)):
        os.makedirs(dir)
      path = dir + "/" + file
      
# Write latest data
    fh = open(path, "w")
    fh.write(self.file_data[0])
    fh.close()
    
  def load_original(self, file=None, dir=None):
    self.load_data(self.file_data[0])
    
############################
#  Set
############################
  
# Control
####################
    
  def set_calculation(self, type=None):
    if(type == None):
      type = "SCF"    
    self.control['calculation'] = '"' + type.lower() + '"'
    
  def set_var(self, var, value):
    for key in sorted(self.control.keys()):
      if(key.upper().strip() == var.upper().strip()):
        self.control['key'] = value
        return True
    return False

  def set_scratch_dir(self, scratch_dir = None):
    if(scratch_dir != None):
      self.scratch_dir = scratch_dir
    self.control['outdir'] = '"' + self.scratch_dir + '"'

  def set_pp_dir(self, pp_dir = None):
    if(pp_dir != None):
      self.pp_dir = pp_dir
    self.control['pseudo_dir'] = '"' + self.pp_dir + '"'

  def set_dirs(self, scratch_dir = None, pp_dir = None):
    self.set_scratch_dir(scratch_dir)
    self.set_pp_dir(pp_dir)

  def set_prefix(self, pin = None):    
    if(pin == None):
      pin = pwscf_input.rand_string()
    self.control['prefix'] = '"' + pin + '"'

# System
####################
    
  def set_alat(self, alat, variance=0.0, rand_seed=0):
    if(variance != 0.0):
      if(rand_seed == 1 and self.rand_seed_set == 0):
        self.rand_seed_set = 1
        self.rand.randomSeed()
      alat = alat + variance * self.rand.rng()  
      
    self.system['celldm'][0] = str(alat)
    
  def set_degauss(self, degauss):  
    self.system['degauss'] = str(degauss)

  def set_ecutrho(self, ecutrho):  
    self.system['ecutrho'] = str(ecutrho)

  def set_ecutwfc(self, ecutwfc):  
    self.system['ecutwfc'] = str(ecutwfc)

  def set_nosym(self, nosym=False):  
    if(nosym == True or nosym.lower() == ".true."):
      self.system['nosym'] = ".TRUE."
    else:
      self.system['nosym'] = ".FALSE."

  def set_nspin(self, value):
    if(value == 2 or value == "2"):
      self.system['nspin'] = "2"
    elif(value == 4 or value == "4"):
      self.system['nspin'] = "4"
    else:  
      self.system['nspin'] = "1"

  def set_tot_magnetization(self, tot_magnetization):  
    self.system['tot_magnetization'] = str(tot_magnetization)
      
  def set_as_isolated(self, alat=10.0):
    self.set_alat(alat)
    self.set_nosym(True)
    self.set_nspin(2)
    self.set_tot_magnetization(0)
    self.set_cp_identity()
    self.load_config("ISOLATED")

# Cell Parameters
####################

# Just the 3x3 array
  def set_cell_parameters(self, cell_in):
    type = self.cell_parameters[0]
    self.cell_parameters = [type]
    for row in cell_in:
      new_row = []
      for cell in row:
        new_row.append(str(cell))
      self.cell_parameters.append(new_row) 
      
# Just the 3x3 array
  def set_cp_arr(self, cp):
    if(len(cp) == 6):
      cp = pwscf_standard.unvoight(cp)
    self.set_cell_parameters(cp)
    
  def set_cp_strain(self, cp, strain):
    cp = numpy.matmul(strain, cp)
    self.set_cell_parameters(cp)
      
# Copy entire list [,[,,],[,,],[,,]]
  def set_cp(self, cp):
    if(isinstance(cp, numpy.ndarray)):
      cp_list = []
      for i in range(3):
        cp_list.append([cp[i,0],cp[i,1], cp[i,2]])
      self.set_cell_parameters(cp_list)
    elif(isinstance(cp, list)):
      self.cell_parameters = cp  
    
  def set_cp_identity(self):
    self.set_cell_parameters([[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]])
    
  def set_cp_zeros(self):
    self.set_cell_parameters([[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]])
    
  def nomalise_cell_parameters(self):
    self.system['celldm'][0] = str(float(self.system['celldm'][0]) * float(self.cell_parameters[1][0]))
    d = float(self.cell_parameters[1][0])
    for i in range(1,4):
      for j in range(0,3): 
        self.cell_parameters[i][j] = str(float(self.cell_parameters[i][j]) / d)
        
# K-Points
####################
  def set_k_points(self, k_points_type, points_list):
    if(type(points_list) == list):
      for i in range(len(points_list)):
        points_list[i] = str(points_list[i])
    if(type(points_list) == str):
      points_list = pwscf_input.clean(points_list)
      points_list = points_list.split(" ")
    self.k_points = []
    self.k_points.append(k_points_type)
    self.k_points.append(points_list)
    
############################
#  Load From Output
############################
  
  def load_from_output(self, file_path):
# Need to add a fail safe, and load atom coords
  
# Load output file
    output_file = pwscf_output(file_path)
    alat, relaxed_cell = output_file.get_norm_relaxed()
    
# Set alat and cell parameters
    self.set_alat(alat)
    self.set_cell_parameters(relaxed_cell)
    
############################
#  Load Atom Config
############################

  def load_config(self, type="FCC", size=1, alat=None, cp=None): 
    type = type.upper()
    result = {"CrystalAtoms": 0, "TotalAtoms": 0}
    
# FCC
    if(type == "FCC"):
      labels = []
      for row in self.atomic_species:
        labels.append(row[0])
      atoms, c_atoms, n_atoms = pwscf_standard.fcc(labels,size) 
      self.atomic_positions = atoms
      self.system['nat'] = str(len(self.atomic_positions) - 1)
    
# BCC
    if(type == "BCC"):
      labels = []
      for row in self.atomic_species:
        labels.append(row[0])
      atoms, c_atoms, n_atoms = pwscf_standard.bcc(labels,size) 
      self.atomic_positions = atoms
      self.system['nat'] = str(len(self.atomic_positions) - 1)
    
# ISOLATED
    if(type == "ISOLATED"):
      labels = []
      for row in self.atomic_species:
        labels.append(row[0])
      atoms, c_atoms, n_atoms = pwscf_standard.isolated(labels)  
      self.atomic_positions = atoms
      self.system['nat'] = str(len(self.atomic_positions) - 1)      

# ALAT
    if(alat != None):
      self.set_alat(alat)

# CP
    if(cp != None):
      self.set_cp_arr(cp)
      
# Misc
    if(type == "FCC_VARIED_DOUBLE"):
      labels = []
      for row in self.atomic_species:
        labels.append(row[0])
      atoms, c_atoms, n_atoms = pwscf_standard.fcc_varied_double(labels) 
      self.atomic_positions = atoms
      self.system['nat'] = str(len(self.atomic_positions) - 1)    
    
# Make
    self.make()
    
# Return
    return c_atoms, n_atoms
    
  def load_custom_config(self, coords, atoms=None):  
    atom_count = len(coords)
    self.system['nat'] = str(atom_count)   
    
#for coord in coords:
  
############################
#  Vary Config
############################
  
  def vary_atom(self, n, dx, dy, dz):
    n = 1 + n % (len(self.atomic_positions) - 1)
    self.atomic_positions[n][1] = str(float(self.atomic_positions[n][1]) + float(dx))
    self.atomic_positions[n][2] = str(float(self.atomic_positions[n][2]) + float(dy))
    self.atomic_positions[n][3] = str(float(self.atomic_positions[n][3]) + float(dz))
  
  def rand_vary(self, amount= 0.01, rand_seed=0):
    amount = float(amount)
    rand_seed = int(rand_seed)
    if(rand_seed == 1 and self.rand_seed_set == 0):
      self.rand_seed_set = 1
      self.rand.randomSeed()
  
    c = numpy.zeros((3,3))
    for i in range(1,4):
      for j in range(0,3): 
        c[i-1,j] = self.cell_parameters[i][j]
    c_inv = numpy.linalg.inv(c)

    for n in range(1, len(self.atomic_positions)):
      r = numpy.zeros((3))
      r[0] = (amount * self.rand.rng()) * (1.0 / float(self.system['celldm'][0]))
      r[1] = (amount * self.rand.rng()) * (1.0 / float(self.system['celldm'][0]))
      r[2] = (amount * self.rand.rng()) * (1.0 / float(self.system['celldm'][0]))
    
      r = numpy.matmul(c_inv, r)
      self.atomic_positions[n][1] = str(float(self.atomic_positions[n][1]) + r[0])
      self.atomic_positions[n][2] = str(float(self.atomic_positions[n][2]) + r[1])
      self.atomic_positions[n][3] = str(float(self.atomic_positions[n][3]) + r[2])

############################
#  Get
############################

  def get_path(self):
    if(self.file_name == None):
      file = "pwscf.in"
    else:
      file = self.file_name
    if(self.dir_name == None):  
      path = file
    else:
      path = self.dir_name + "/" + file
    return path

  def get_file_name(self):
    return self.file_name

  def get_cp_array(self):
    cp = numpy.zeros((3,3))
    for i in range(3):
      for j in range(3):
        cp[i,j] = float(self.cell_parameters[i+1][j])
    return cp    

  def get_data(self, make=False):
    if(make):
      self.make()
    return self.file_data[-1]
    
  def get_nat(self):  
    return int(self.system['nat'])

  def get_alat(self):
    return self.system['celldm'][0]

#################################
# Signature
#################################
    
  def signature(self):
    skip = ['outdir','prefix','wfcdir','pseudo_dir','max_seconds','disk_io']
  
# CONTROL
    file = "&CONTROL \n"
    for key in sorted(self.control.keys()):
      read = True
      if(key in skip):
        read = False
        
      if(read):
        file += pwscf_input.make_line(key, self.control[key])
    file += "/ \n"    
# SYSTEM
    file += "&SYSTEM \n"
    for key in sorted(self.system.keys()):
      file += pwscf_input.make_line(key, self.system[key]) 
    file += "/ \n"
# ELECTRONS
    file += "&ELECTRONS \n"
    for key in sorted(self.electrons.keys()):
      value = self.electrons[key]
      if(value != None):
        file += key + " = " + value + ", \n"      
    file += "/ \n"
# IONS
    file += "&IONS \n"
    for key in sorted(self.ions.keys()):
      value = self.ions[key]
      if(value != None):
        file += key + " = " + value + ", \n"      
    file += "/ \n"
# CELL
    file += "&CELL \n"
    for key in sorted(self.cell.keys()):
      value = self.cell[key]
      if(value != None):
        file += key + " = " + value + ", \n"      
    file += "/ \n"
# ATOMIC_SPECIES
    file += "ATOMIC_SPECIES \n"
    for species in self.atomic_species:      
      for field in species:
        file += str(field) + " "
      file += "\n"
# ATOMIC_POSITIONS
    header = 0
    for position in self.atomic_positions:      
      if(header == 0):
        file += "ATOMIC_POSITIONS "
        file += position + "\n"
        header = 1
      elif(header == 1):  
        for field in position:
          file += str(field) + "   "
        file += "\n"
# K_POINTS
    file += "K_POINTS " + self.k_points[0]
    file += "\n"
    for i in range(1,len(self.k_points)):
      for point in self.k_points[i]:
        file += point + " "
      file += "\n"    
# K_POINTS
    file += "CELL_PARAMETERS " + self.cell_parameters[0]
    file += "\n"
    for i in range(1,len(self.cell_parameters)):
      for point in self.cell_parameters[i]:
        file += point + " "
      file += "\n"
    
# String being hashed must be converted to utf-8 encoding
    input_string = file.encode('utf-8')
# Make hash object
    my_hash = hashlib.sha512()
# Update
    my_hash.update(input_string)
# Return hash
    return my_hash.hexdigest()

#################################
# Interactive
#################################

  def menu(self):
    while(True):
      choice = self.print_menu().upper()
      print(choice)
      if(choice == "X"):
        exit()
      elif(choice == "1"):
        self.i_load()
      elif(choice == "2"):
        self.i_display()

  def print_menu(self):
    pwscf_input.header("Menu")
    print("1. Load File")
    print("2. Display File")
    print("X. Exit")
    return input("Choice: ")

  def i_load(self):
    pwscf_input.header("Load Input File")
    file_name = input("Enter file name: ")
    self.load(file_name)
    input()

  def i_display(self):
    self.make()
    pwscf_input.header("Display File")
    print(self.file)
    input()

#################################
# Help
#################################

  def help(self):
    print("HELP")

#################################
# Static Methods
#################################

  @staticmethod
  def is_pwscf(file_path): 
    file_size = os.path.getsize(file_path)
    if(file_size > 500000):
      return False

    fh = open(file_path, 'r')
    counter = 0
    for line in fh:
      line = line.strip().lower()
      if(line[0:8] == "&control"):
        counter = counter + 1
      if(line[0:7] == "&system"):
        counter = counter + 1
      if(line[0:10] == "&electrons"):
        counter = counter + 1
      if(line[0:5] == "&ions"):
        counter = counter + 1
      if(line[0:5] == "&cell"):
        counter = counter + 1
      if(line[0:14] == "&atomic_species"):
        counter = counter + 1
      if(line[0:16] == "&atomic_positions"):
        counter = counter + 1
      if(line[0:6] == "outdir"):
        counter = counter + 1
      if(line[0:10] == "pseudo_dir"):
        counter = counter + 1
    if(counter > 5):
      return True
    return False

  @staticmethod
  def remove_spaces(input_string):
    return input_string.replace(" ", "")
 
  @staticmethod
  def fields(input_string):
    input_string = input_string.strip()
    output_string = ""
    last = None
    for character in input_string:
      if(character != " " or (character == " " and last != " ")):
        output_string += character
    return output_string.split(" ")
    
  @staticmethod
  def check_keyword(line, keyword):
    if(line.upper()[0:len(keyword)] == keyword.upper()):
      return True
    return False

  @staticmethod
  def clear_screen():
    os.system('cls' if os.name == 'nt' else 'clear')

  @staticmethod
  def header(sub_title=""):
    pwscf_input.clear_screen()
    print("==========================================================")
    print("                    PWscf Input Editor                    ")
    print("==========================================================")
    print()
    print(sub_title)
    print()
    print()
    print()

  @staticmethod
  def process_keyword(str_in):
    str_in = str_in.lower().strip()
    str_in = pwscf_input.remove_spaces(str_in)
    id = None
    keyword = ""
    flag = 0
    for character in str_in:
      if(character == "("):
        id = ""
        flag = 1
      elif(character == ")"):
        flag = 2
      elif(flag == 0):
        keyword += character
      elif(flag == 1):
        id = id + character
    if(id != None):
      try:
        id = int(id)
      except:
        id = None
    return keyword, id  

  @staticmethod
  def add_keyword(keywords, keyword, id, value):
    if(id == None):
      added = False
      for i in range(len(keywords)):
        if(keywords[i][0] == keyword):
          added = True
          keywords[i][1] = keyword
      if(added == False):
        keywords.append([keyword, value])
    else:   
      n = None
      for i in range(len(keywords)):
        if(keywords[i][0] == keyword):
          n = i
          break
      if(n == None):    
        keywords.append([keyword,[None]])
        n = len(keywords) - 1
        
      while(len(keywords[n][1]) < id):
        keywords[n][1].append(None)

      keywords[n][1][id-1] = value  

  @staticmethod
  def make_line(key, value):
    output = ""
    if(value != None):
       if(isinstance(value, (list,))):
         for i in range(len(value)):
           if(value[i] != None):
             output += key + "(" + str(i+1) + ") = " + value[i] + ", \n"                
       else:
         output += key + " = " + value + ", \n"   
    return output    

  @staticmethod
  def coord_format(float_in):
    pad = "              "
    value = str(round(float_in, 6)).strip()
    return value
    
  @staticmethod
  def label_format(label):  
    pad = "              "
    label = label.strip()
    return label
    
  @staticmethod
  def clean(str_in):  
    str_out = ""
    l = len(str_in)
    for i in range(l):
# Last, Next, This
      if(i == 0):
        last = None
      else:
        last = str_in[i-1]
      if(i < (l-1)):
        next = str_in[i+1]
      else:  
        next = None
      char = str_in[i]
    
# Check
      ok = True
      if(last == " " and char == " "):
        ok = False
      elif(last == "\n" and char == "\n"):
        ok = False
      elif(last == "\n" and char == " "):
        ok = False
      elif(char == " " and next == "\n"):
        ok = False
      elif(last == "=" and char == " "):
        ok = False
      elif(char == " " and next == "="):
        ok = False
        
# Add to string
      if(ok):
        str_out += char
    return str_out    
    
  @staticmethod
  def rand_string(len_in=16):      
    output = ""
    char_set = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOQRSTUVWXYZ"
    for i in range(len_in):
      r = random.randint(0,len(char_set)-1)
      output += char_set[r]
    return output
    
#data = re.sub(r'\s\s+', ' ', data)
#data = re.sub(r'\s=\s', '=', data)

###########################################
#  CLASS rand_dis
###########################################
class rand_dist:

  def __init__(self, seed=None):
# RNG vars
    if seed is None:
      self.seed = 12791244
    else:
      self.seed = seed
    self.xn = self.seed
    self.m = 4294967296
    self.a = 1103515245
    self.c = 12345
#Dist Type vars
    self.lower = 0.0e0
    self.upper = 1.0e0
    self.distType = "flat"
# make a new distribution table
    self.table = RandDistTable()

  def randomSeed(self):
    currentTime = time.time()
    myhash = hashlib.md5()
    message = "randnum"+str(currentTime)
    myhash.update(message.encode())
    seed = myhash.hexdigest()
    seedInt = int("0x"+seed[0:8],0)
    self.xn = (self.xn + seedInt) % self.m

  def setRange(self,lower,upper):
    self.lower = lower
    self.upper = upper

  def flat(self):
    self.distType = "flat"
    self.table.set = False

  def sqrt(self):
    self.distType = "sqrt"
    self.table.set = False

  def gheat(self, lower=-1.0e0, upper=1.0e0, p1=1.0, p2=0.0, p3=1.0, p4=4.0):
    self.distType = "gheat"
    self.table.set = False
    self.lower = lower
    self.upper = upper
    self.p1 = p1  # sigma
    self.p2 = p2  # mu
    self.p3 = p3  # factor
    self.p4 = p4  # factor

  def doubleGaussian(self, lower=-1.0e0, upper=1.0e0, p1=1.0, p2=0.0, p3=1.0, p4=4.0, p5=1.5, p6=0.8, p7=1.0, p8=4.0):
    self.distType = "doubleGaussian"
    self.table.set = False
    self.lower = lower
    self.upper = upper
    self.p1 = p1
    self.p2 = p2
    self.p3 = p3
    self.p4 = p4
    self.p5 = p5
    self.p6 = p6
    self.p7 = p7
    self.p8 = p8

  def getFloat(self):
# Get a random float between 0 and 1
    self.xn = (self.a * self.xn + self.c) % self.m
    randFloat = self.xn / self.m
    return randFloat

  def rng(self):
# Choose function
    if(self.distType=="flat"):
      self.randFloat = self.getFloat()
      randOut = self.lower + self.randFloat * (self.upper - self.lower)
# Square root
    if(self.distType=="sqrt"):
      self.randFloat = self.getFloat()
      randOut = self.lower + self.randFloat * (self.upper - self.lower)
      randOut = math.sqrt(abs(randOut))
# gheat
    if(self.distType=="gheat"):
      randOut = self.randDistF()
# doubleGaussian
    if(self.distType=="doubleGaussian"):
      randOut = self.randDistF()
# store and return output
    self.randOut = randOut
    return randOut

  def randDistF(self):
# If a table hasn't been built for the function, build it now
    if(self.table.set==False):
      self.rand_MakeTable()
#self.table.display()

# Loop attempts
    loopTrials = True
    i = 0
    while(loopTrials and i<10000):
      i = i + 1
      randA = self.getFloat()  # Block
      randB = self.getFloat()  # xVal
      randC = self.getFloat()  # yVal
      if(randA==0):
        randBlock = 0
      else:
        randBlock = math.ceil(randA * self.table.size)-1
      xA = self.table.points_xA[randBlock]
      x_m = self.table.points_x_m[randBlock]
      xB = self.table.points_xB[randBlock]
      yA = self.table.points_yA[randBlock]
      y_m = self.table.points_y_m[randBlock]
      yB = self.table.points_yB[randBlock]
      yMax = self.table.points_yMax[randBlock]
# set x,y coords
      x = xA + randB * (xB - xA)
      y = randC * yMax
# Set up list
      interpPoints = [[0 for x in range(2)] for y in range(3)]
      interpPoints[0][0] = xA
      interpPoints[0][1] = yA
      interpPoints[1][0] = x_m
      interpPoints[1][1] = y_m
      interpPoints[2][0] = xB
      interpPoints[2][1] = yB
# Interp
      fx = self.interp(x,interpPoints)
      if(y<=fx):
        loopTrials = False
      randNumber = x

    return randNumber

  def rand_MakeTable(self):
# Calculate total area under function
    self.table.initPoints(256)
    randRange = self.upper - self.lower
    xInc = randRange / (self.table.targetSize - 1)
    xA = self.lower
    xB = xA + xInc
    area = 0.0
    for i in range(0,self.table.targetSize):
### Choose function
###---------------------------------------------
      if(self.distType=="gheat"):
        yA = self.rand_Gaussian(xA)
        yB = self.rand_Gaussian(xB)
      if(self.distType=="doubleGaussian"):
        yA = self.rand_DoubleGaussian(xA)
        yB = self.rand_DoubleGaussian(xB)
###--------------------------------------------
# Calc area
      area = area + 0.5 * xInc * (yA + yB)
# Increment
      xA = xB
      xB = xA + xInc
# Area when split into blocks
    blockArea = area / self.table.targetSize
    xA = self.lower
    runLoop = True
    i = 0
    while(runLoop):
      i = i + 1
### Choose function
###---------------------------------------------
      if(self.distType=="gheat"):
        yA = self.rand_Gaussian(xA)
      if(self.distType=="doubleGaussian"):
        yA = self.rand_DoubleGaussian(xA)
###--------------------------------------------
# Find xB
      xInc = randRange * 0.0001
      xB = xA
      area = 0.0
      yMax = yA
      findArea = True
      while(findArea):
        xB = xB + xInc
        yB = self.rand_Gaussian(xB)
        if(yB>yA):
          yMax = yB
        area = area + 0.5 * xInc * yMax
        if(area>blockArea):
          findArea = False # Area has been found
      xB = xA + blockArea / yMax
      if(xB>self.upper):
        xB = self.upper
        yMax = area/(xB - xA)
        runLoop = False
# set table size
        self.table.size = i
#Run calculations
### Choose function
###--------------------------------------------
      x_m = 0.5 * (xA + xB)
      if(self.distType=="gheat"):
        yA = self.rand_Gaussian(xA)
        y_m = self.rand_Gaussian(x_m)
        yB = self.rand_Gaussian(xB)
      if(self.distType=="doubleGaussian"):
        yA = self.rand_DoubleGaussian(xA)
        y_m = self.rand_DoubleGaussian(x_m)
        yB = self.rand_DoubleGaussian(xB)
###--------------------------------------------
# Store data points
      self.table.points_yMax.append(yMax)
      self.table.points_xA.append(xA)
      self.table.points_x_m.append(x_m)
      self.table.points_xB.append(xB)
      self.table.points_yA.append(yA)
      self.table.points_y_m.append(y_m)
      self.table.points_yB.append(yB)
# Set xA for next loop
      xA = xB
# Set table to true
    self.table.set = True

  def rand_Gaussian(self, x):
    fx = self.p3 * (0.398942 / self.p1)
    expV = (-0.5/(self.p1*self.p1)) * (self.p4 * self.p4 * (x - self.p2) * (x - self.p2))
    fx = fx * math.exp(expV)
    return fx

  def rand_DoubleGaussian(self, x):
    fxA = self.p3 * (0.398942 / self.p1)
    expV = (-0.5/(self.p1*self.p1)) * (self.p4 * self.p4 * (x - self.p2) * (x - self.p2))
    fxA = fxA * math.exp(expV)
    fxB = self.p7 * (0.398942 / self.p5)
    expV = (-0.5/(self.p5*self.p5)) * (self.p8 * self.p8 * (x - self.p6) * (x - self.p6))
    fxB = fxB * math.exp(expV)
    fx = fxA + fxB
    return fx

  def makeTally(self, tallySize=50, sampleSize=100000):
# Declare list
    self.tally = []
    self.tallyX = []
    self.tallySize = tallySize
    sampleSize = 100000
    self.sampleSize = sampleSize
    halfIncrement = (self.upper - self.lower) / (2 * tallySize)
# Loop through
    for i in range(0,tallySize):
      self.tally.append(0)
      xVal = self.lower + (2 * i + 1) * halfIncrement
      self.tallyX.append(xVal)
    for i in range(0,self.sampleSize):
      randNum = self.rng()
      binNum = randNum - self.lower
      binNum = binNum / (self.upper - self.lower)
      binNum = math.floor(binNum * tallySize)
# Force to be in range
      if(binNum>=tallySize):
        binNum = tallySize - 1
      if(binNum<0):
        binNum = 0
      self.tally[binNum] = self.tally[binNum] + 1

  def displayTally(self):
    print("==============================")
    print("     Tally: ",self.distType)
    print("==============================")
    for i in range(0,self.tallySize):
#print(i,"  ",self.tally[i],"  ",self.tallyX[i],"  ",self.tally[i])
      print('{0:3G},{1:8G},{2:5f},{3:8G}'.format(i,self.tally[i],self.tallyX[i],self.tally[i]))
    print("==============================")
    print()

  @staticmethod
  def interp(x, points):
# Lagrange interpolation
    output = 0.0
    count = len(points)
    coefficients = []
# Make coefficients
    for i in range(0,count):
      numerator = 1.0
      denominator = 1.0
      for j in range(0,count):
        if(i!=j):
          numerator = numerator * (x - points[j][0])
          denominator = denominator * (points[i][0] - points[j][0])
      coefficients.append(numerator/denominator)
# Calculate y
    y = 0
    for i in range(0,count):
      y = y + points[i][1] * coefficients[i]
    return y

###########################################
#  CLASS RandDistTabl
###########################################
class RandDistTable:
  def __init__(self):
    self.set = False

  def initPoints(self,targetSize=256):
    self.size = 0
    self.targetSize = targetSize
    self.points_yMax = []
    self.points_xA = []
    self.points_x_m = []
    self.points_xB = []
    self.points_yA = []
    self.points_y_m = []
    self.points_yB = []

  def display(self):
    print("======================")
    print("        Table         ")
    print("======================")
    print("Size:   ", self.size )
    print("")
    for i in range(0,self.size):
      print(i,"     ",end="")
      print(self.points_yMax[i],"     ",end="")
      print(self.points_xA[i],self.points_x_m[i],self.points_xB[i],"",sep="    ",end="")
      print(self.points_yA[i],self.points_y_m[i],self.points_yB[i],"",sep="    ",end="\n")
    print("")
    print("")
    print("======================")

############################################################
# Random Number Generator (testing)
############################################################

###########################################
#  CLASS Test_RandDis
###########################################
class Test_RandDist:

  @staticmethod
  def testRandFloatSingle():
    rd = rand_dist()
    print(rd.rng())
    print(rd.rng())
    print(rd.rng())
    rd.randomSeed()
    print(rd.rng())
    print(rd.rng())
    print(rd.rng())

  @staticmethod
  def testRandFloat():
    rd = rand_dist()
    rd.makeTally()
    rd.displayTally()

  @staticmethod
  def testRandSqrt():
    rd = rand_dist()
    rd.sqrt()
    rd.makeTally()
    rd.displayTally()

  @staticmethod
  def testRandGheat():
    rd = rand_dist()
    rd.gheat(-1.0,1.0,1.0,0.0,1.0,4.0)
    rd.makeTally()
    rd.displayTally()

  @staticmethod
  def testRandDoublegaussian():
    rd = rand_dist()
    rd.doubleGaussian(-1.0,1.0,1.0,0.0,1.0,4.0,1.5,0.8,1.0,4.0)
    rd.makeTally()
    rd.displayTally()

###########################################
#  CLASS pwscf_outpu
###########################################
class pwscf_output:

  def __init__(self, file_in=None):
    self.reset()
    if(file_in != None):
      self.load(file_in)

  def reset(self):
  
    self.z = numpy.zeros((3,3))

# Control
    self.data = {
      "ok": False,
      "job_done": False,
      "error": False,
      "type": None,
      "summary": None,
      "mpi_processes": None,
      "bravais_lattice_index": None,
      "alat": None,
      "volume": None,
      "electrons": None, 
      "electrons_up": None, 
      "electrons_down": None, 
      
      "crystal_in": numpy.zeros((3,3)),
      "crystal_calc": numpy.zeros((3,3)),
      
      "total_energy": None,
      "density_full": None,
      "density": None,
      "stress": numpy.zeros((3,3)),
      "stress_sum": None,
      
      "cpu_time": None,
      "wall_time": None
      
    }

#  Load, and use in another program
  def load(self, file_name): 
# Load data from file
    data = self.load_from_file(file_name)
    
# Read through data
    self.load_data(data)

# Load from a block of data (text, file etc)
  def load_data(self, data):  
    
# Reset data store
    self.reset()
    
# split
    data = data.split("\n")    
    
# OK
###################################
    self.data['ok'] = False
    counter = 0
    for line in data:
      line = line.strip()
      if(pwscf_output.compare(line, "JOB DONE.")):
        counter = counter + 1
      if(pwscf_output.compare(line, "Exit code:")):
        counter = counter - 1
      if(pwscf_output.compare(line, "convergence NOT achieved")):
        counter = counter - 1
    if(counter == 1):
      self.data['ok'] = True
    
# Calc Type
###################################
    self.data['type'] = "SCF"
    for line in data:
      line = line.strip()
      if(line[0:23] == "A final scf calculation"):
        self.data['type'] = "VC-RELAX"
        
# Load
###################################
    n = 0
    counter = 0
    while(n < len(data)):
      n, line, line_uc = self.next_line(n, data)
      if(line != ""):
        counter += 1
        if(counter == 1):
          self.data['summary'] = line
        else:
          if(pwscf_output.compare(line, "Number of MPI processes:")):
            self.data['mpi_processes'] = pwscf_output.extract(line, ":", "", "i") 
        
          if(pwscf_output.compare(line, "bravais-lattice index     =")):
            self.data['bravais_lattice_index'] = pwscf_output.extract(line, "=", "", "i") 
            
          if(pwscf_output.compare(line, "lattice parameter (alat)  =")):
            self.data['alat'] = pwscf_output.extract(line, "=", "a.u.", "f")  
            
          if(pwscf_output.compare(line, "unit-cell volume          =")):
            self.data['volume'] = pwscf_output.extract(line, "=", "(a.u.)^3", "f")     
            
          if(pwscf_output.compare(line, "number of atoms/cell      =")):
            self.data['nat'] = pwscf_output.extract(line, "=", "", "i")  
            
          if(pwscf_output.compare(line, "number of atomic types    =")):
            self.data['types'] = pwscf_output.extract(line, "=", "", "i")  
            
          if(pwscf_output.compare(line, "number of electrons       =")):
            str_e = pwscf_output.extract(line, "=", "", "s")
            e, eu, ed = pwscf_output.electron_string(str_e)
            self.data['electrons'] = e
            self.data['electrons_up'] = eu
            self.data['electrons_down'] = ed
          
          if(pwscf_output.compare(line, "number of Kohn-Sham states=")):
            self.data['ks_states'] = pwscf_output.extract(line, "=", "", "i")   
            
          if(pwscf_output.compare(line, "kinetic-energy cutoff     =")):
            self.data['ecutwfc'] = pwscf_output.extract(line, "=", "Ry", "f")  
            
          if(pwscf_output.compare(line, "charge density cutoff     =")):
            self.data['ecutrho'] = pwscf_output.extract(line, "=", "Ry", "f")   
        
          if(pwscf_output.compare(line.strip(), "crystal axes:") and pwscf_output.is_zero(self.data['crystal_in'])):            
            for j in range(3):              
              n, line, line_uc = self.next_line(n, data)
              fields = pwscf_output.extract(line, "= (", ")", "s", " ")
              self.data['crystal_in'][j,:] = fields  
              self.data['crystal_calc'][j,:] = fields  
          
          if(pwscf_output.compare(line.strip(), "crystal axes:")):            
            for j in range(3):              
              n, line, line_uc = self.next_line(n, data)
              fields = pwscf_output.extract(line, "= (", ")", "s", " ")
              self.data['crystal_calc'][j,:] = fields            
        
          if(pwscf_output.compare(line, "!    total energy")):
            self.data['total_energy'] = pwscf_output.extract(line, "=", "Ry", "f")
            
          if(pwscf_output.compare(line, "Total force =")):
            self.data['total_force'] = pwscf_output.extract(line, "=", "T", "f")
            
          if(pwscf_output.compare(line, "total   stress  (Ry/bohr**3)")):        
            self.data['stress_sum'] = 0.0
            for j in range(3):              
              n, line, line_uc = self.next_line(n, data)   
              fields = pwscf_output.extract(line, "", "", "f", " ", True)  
              self.data['stress'][j,0] = fields[0] 
              self.data['stress'][j,1] = fields[1] 
              self.data['stress'][j,2] = fields[2]
              self.data['stress_sum'] = self.data['stress_sum'] + abs(fields[0]) + abs(fields[1]) + abs(fields[2])
            
#                  "stress": numpy.zeros((3,3)),
#      "stress_sum": None,
            
          if(pwscf_output.compare(line, "density = ")):
            self.data['density_full'] = pwscf_output.extract(line, "=", "", "s")
            self.data['density'] = pwscf_output.extract(line, "=", "g/cm^3", "f")
          
          if(pwscf_output.compare(line, "PWSCF        :")):
            self.data['cpu_time'] = pwscf_output.extract(line, ":", "CPU", "s")
            
          if(pwscf_output.compare(line, "PWSCF        :")):
            self.data['wall_time'] = pwscf_output.extract(line, "CPU", "WALL", "s")
          
          if(pwscf_output.compare(line, "JOB DONE.")):
            self.data['job_done'] = True
            
          if(pwscf_output.compare(line, "Exit code:")):
            self.data['error'] = True  
             
  def next_line(self, n, data):
    if(n < len(data)):
      line = data[n].strip()
      line_uc = line.upper()
      n = n + 1
      return n, line, line_uc
    else:
      n = n + 1
      return n, None, None
    
  def store(self, store, line, field, n=0):
    l, f = pwscf_output.read_line(line, field)  
    if(l != False):
      self.data[store] = f[n]

#  Run as it's own program
  def run(self):
    self.reset()

    option = ""
    file_name = ""

    if(len(sys.argv) > 1 and sys.argv[1] is not None):
      option = sys.argv[1]

    if(len(sys.argv) > 2 and sys.argv[2] is not None):
      file_name = sys.argv[2]

    if(option.lower().strip() == "" or option.lower().strip() == "interactive"):
      self.menu()
      exit()
    elif(option.lower().strip() == "quiet"):
      print("Quiet")
    else:
      return 0

#################################
# READ/LOAD input file
#################################

  def load_from_file(self, file_name):
# Init variable
    file_data = ""

# Read it in line by line
    fh = open(file_name, "r")
    for file_row in fh:
      file_data = file_data + file_row.strip() + '\n'

    return file_data

#################################
# Get
#################################

  def get_alat(self):
    return self.data['alat']
    
  def get_volume(self):
    return self.data['volume']  
  
  def get_total_energy(self):
    return self.data['total_energy']  
    
  def get_energy_per_atom(self):
    if(self.data['total_energy'] is None):
      return None
    elif(self.data['nat'] is None):
      return None
    else:
      return (float(self.data['total_energy']) / float(self.data['nat']))    
  
  def get_total_force(self):
    return self.data['total_force']  
    
  def get_force_per_atom(self):
    return (float(self.data['total_force']) / float(self.data['nat']))  
  
  def get_density(self):
    return self.data['density']  
    
  def get_cell_parameters(self):
    cp = ['alat', 
          [str(self.data['crystal_calc'][0,0]), str(self.data['crystal_calc'][0,1]), str(self.data['crystal_calc'][0,2])], 
          [str(self.data['crystal_calc'][1,0]), str(self.data['crystal_calc'][1,1]), str(self.data['crystal_calc'][1,2])], 
          [str(self.data['crystal_calc'][2,0]), str(self.data['crystal_calc'][2,1]), str(self.data['crystal_calc'][2,2])]]
    return cp

# Return relaxed unit vector
  def get_cell_array(self):
    return self.data['crystal_calc']

# return alat and normalised unit vector
  def get_norm_relaxed(self):
    alat = self.data['alat']    
    cp = numpy.copy(self.data['crystal_calc'])
    f = cp[0,0]
    alat = alat * f
    cp = cp / f
    return alat, cp

# Get stress
  def get_stress(self):
    return self.data['stress']
    
  def get_stress_sum(self):
    return self.data['stress_sum']
  
  def get_job_done(self):
    return self.data['job_done']

  def get_ok(self):
    return self.data['ok']

#################################
# Interactive
#################################

  def menu(self):
    while(True):
      choice = self.print_menu().upper()
      print(choice)
      if(choice == "X"):
        exit()
      elif(choice == "1"):
        self.i_load()
      elif(choice == "2"):
        self.i_display()

  def print_menu(self):
    pwscf_output.header("Menu")
    print("1. Load File")
    print("2. Display File")
    print("X. Exit")
    return input("Choice: ")

  def i_load(self):
    pwscf_output.header("Load Output File")
    file_name = input("Enter file name: ")
    self.load(file_name)
    print("File loaded.")
    input()

  def i_display(self):
    pwscf_output.header("Display File")
    self.output_details()
    input()

  def output_details(self):
    print("Output")
    for key in sorted(self.data.keys()):
      value = self.data[key]
      print(key, ":  ", value)
    
#for key, value in self.data.items():
#  print(key, ":  ", value)

#################################
# Static Methods
#################################

  @staticmethod
  def remove_spaces(input_string):
    return input_string.replace(" ", "")
    
  @staticmethod
  def extract(input_string, start=None, end=None, type=None, split=None, trim=False):
    if(start == ""):
      start = None
    if(end == ""):
      end = None
    if(trim):
      input_string = input_string.strip()
    
# Start/End
    start_n = None
    end_n = None
      
    if(start == None and end == None):   
      start_n = 0
      end_n = len(input_string)
    elif(start == None and end != None):  
      end_l = len(end)   
      start_n = 0
      for n in range(len(input_string)):
        if(input_string[n:n+end_l] == end[0:end_l]):
          end_n = n
          break
    elif(start != None and end == None):  
      start_l = len(start)
      end_n = len(input_string)
      for n in range(len(input_string)):
        if(input_string[n:n+start_l] == start[0:start_l]):
          start_n = n + start_l
    else:  
      start_l = len(start)
      end_l = len(end)  
    
      for n in range(len(input_string)):
        if(input_string[n:n+start_l] == start[0:start_l]):
          start_n = n + start_l
        if(start_n != None and input_string[n:n+end_l] == end[0:end_l]):
          end_n = n
          break
        
# Read
    result = input_string[start_n:end_n].strip()       

# Split
    if(split != None):
      if(split == " "):
#result = re.sub(r'\s\s+', ' ', result)
        result = pwscf_output.single_spaces(result)
      result = result.split(split)
      for i in range(len(result)):
        if(type.lower() == "f"):
          result[i] = float(result[i])
        elif(type.lower() == "i"):
          result[i] = int(result[i])
        
    else:  
      if(type.lower() == "f"):
        result = float(result)
      elif(type.lower() == "i"):
        result = int(result)
        
# Return
    return result
      
  @staticmethod
  def compare(line, field):
    line = line.strip()
    line = line.upper() 
    
    field = field.strip()
    field = field.upper()
    
    f_len = len(field)
    if(len(line) >= f_len and line[0:f_len] == field[0:f_len]):
      return True
    return False
    
  @staticmethod
  def read_line(line, field):
    line = line.strip()
#line = re.sub(r'\s\s+', ' ', line)
#line = re.sub(r'\s=\s', '=', line)
    line = pwscf_output.clean(line)
    line_uc = line.upper() 
    
    field = field.strip()
#field = re.sub(r'\s\s+', ' ', field)
#field = re.sub(r'\s=\s', '=', field)
    field = pwscf_output.clean(field)
    field = field.upper()
    
    f_len = len(field)
    if(len(line_uc) >= f_len and line_uc[0:f_len] == field[0:f_len]):
      output = line[f_len:].strip()
      fields = output.split(" ")
      return output, fields      
    return False, False
    
  @staticmethod
  def fields(input_string):
    input_string = input_string.strip()
    output_string = ""
    last = None
    for character in input_string:
      if(character != " " or (character == " " and last != " ")):
        output_string += character
    return output_string.split(" ")
    
  @staticmethod
  def check_keyword(line, keyword):
    if(line.upper()[0:len(keyword)] == keyword.upper()):
      return True
    return False

  @staticmethod
  def clear_screen():
    os.system('cls' if os.name == 'nt' else 'clear')

  @staticmethod
  def header(sub_title=""):
    pwscf_output.clear_screen()
    print("==========================================================")
    print("                    PWscf Input Editor                    ")
    print("==========================================================")
    print()
    print(sub_title)
    print()
    print()
    print()

  @staticmethod
  def process_keyword(str_in):
    str_in = str_in.lower().strip()
    str_in = pwscf_output.remove_spaces(str_in)
    id = None
    keyword = ""
    flag = 0
    for character in str_in:
      if(character == "("):
        id = ""
        flag = 1
      elif(character == ")"):
        flag = 2
      elif(flag == 0):
        keyword += character
      elif(flag == 1):
        id = id + character
    if(id != None):
      try:
        id = int(id)
      except:
        id = None
    return keyword, id  

  @staticmethod
  def add_keyword(keywords, keyword, id, value):
    if(id == None):
      added = False
      for i in range(len(keywords)):
        if(keywords[i][0] == keyword):
          added = True
          keywords[i][1] = keyword
      if(added == False):
        keywords.append([keyword, value])
    else:   
      n = None
      for i in range(len(keywords)):
        if(keywords[i][0] == keyword):
          n = i
          break
      if(n == None):    
        keywords.append([keyword,[None]])
        n = len(keywords) - 1
        
      while(len(keywords[n][1]) < id):
        keywords[n][1].append(None)

      keywords[n][1][id-1] = value  

  @staticmethod
  def make_line(key, value):
    output = ""
    if(value != None):
       if(isinstance(value, (list,))):
         for i in range(len(value)):
           if(value[i] != None):
             output += key + "(" + str(i+1) + ") = " + value[i] + ", \n"                
       else:
         output += key + " = " + value + ", \n"   
    return output    

  @staticmethod
  def coord_format(float_in):
    pad = "              "
    value = str(round(float_in, 6)).strip()
    return value
    
  @staticmethod
  def label_format(label):  
    pad = "              "
    label = label.strip()
    return label
    
  @staticmethod
  def is_zero(arr):
    for i in range(arr.shape[0]):
      for j in range(arr.shape[1]):
        if(arr[i, j] != 0.0):
          return False
    return True
    
  @staticmethod
  def clean(str_in):  
    str_out = ""
    l = len(str_in)
    for i in range(l):
# Last, Next, This
      if(i == 0):
        last = None
      else:
        last = str_in[i-1]
      if(i < (l-1)):
        next = str_in[i+1]
      else:  
        next = None
      char = str_in[i]
    
# Check
      ok = True
      if(last == " " and char == " "):
        ok = False
      elif(last == "\n" and char == "\n"):
        ok = False
      elif(last == "\n" and char == " "):
        ok = False
      elif(char == " " and next == "\n"):
        ok = False
      elif(last == "=" and char == " "):
        ok = False
      elif(char == " " and next == "="):
        ok = False
        
# Add to string
      if(ok):
        str_out += char
    return str_out    
    
  @staticmethod
  def electron_string(str_in):
    arr = str_in.split("(up:")
    e = arr[0]
    if(len(arr) == 1):
      return e.strip(), None, None
    if(len(arr)==2):
      arr_b = arr[1].split(", down:")
      eu = arr_b[0]
      arr_c = arr_b[1].split(")")
      ed = arr_c[0]
      return e.strip(), eu.strip(), ed.strip()
  
    print("TEST")
    return "","",""
    
  @staticmethod
  def single_spaces(str_in):
    str_out = ""
    last = None
    for char in str_in:
      if(char != " " or (char == " " and last != " ")):
        str_out = str_out + char
      last = char
    return str_out
    
###########################################
#  CLASS pwscf_standar
###########################################
class pwscf_standard:

  @staticmethod
  def fcc(label, size=1):
    c_atoms = 4
    n_atoms = c_atoms * size**3
    if(not isinstance(label, (list,))):
      label = [label]
    atoms = ['crystal']
    for x in range(size):
      for y in range(size):
        for z in range(size):
          coords = [[str((x+0.0)/size), str((y+0.0)/size), str((z+0.0)/size)],
                    [str((x+0.5)/size), str((y+0.5)/size), str((z+0.0)/size)],
                    [str((x+0.5)/size), str((y+0.0)/size), str((z+0.5)/size)],
                    [str((x+0.0)/size), str((y+0.5)/size), str((z+0.5)/size)]]
          for i in range(len(coords)):
            atoms.append([label[i % len(label)],coords[i][0],coords[i][1],coords[i][2]])
    return atoms, c_atoms, n_atoms     

  @staticmethod
  def bcc(label, size=1):
    c_atoms = 2
    n_atoms = c_atoms * size**3
    if(not isinstance(label, (list,))):
      label = [label]
    atoms = ['crystal']
    for x in range(size):
      for y in range(size):
        for z in range(size):
          coords = [[str((x+0.0)/size), str((y+0.0)/size), str((z+0.0)/size)],
                    [str((x+0.5)/size), str((y+0.5)/size), str((z+0.5)/size)]]
          for i in range(len(coords)):
            atoms.append([label[i % len(label)],coords[i][0],coords[i][1],coords[i][2]])
    return atoms, c_atoms, n_atoms    

  @staticmethod
  def isolated(label, size=1):   
    c_atoms = 1
    n_atoms = 1
    if(not isinstance(label, (list,))):
      label = [label]
    atoms = ['crystal']
    atoms.append([label[0], "0.5", "0.5", "0.5"])
    return atoms, c_atoms, n_atoms    
    
  @staticmethod
  def unvoight(cp_in):
    cp_out = [[cp_in[0], cp_in[5], cp_in[4]] , [cp_in[3], cp_in[1], cp_in[3]] , [cp_in[4], cp_in[5], cp_in[2]]]
    return cp_out

###########################################
#  CLASS pwscf_exe
###########################################
class pwscf_exec:

  @staticmethod
  def execute(files_in, input_dir=None, output_dir=None, working_dir=None, log=None, set_dirs=True, allow_cache=True):
    try:
# Make list of files
      data = []
      if type(files_in) is str:
        data = pwscf_exec.add_files(data, files_in, input_dir, output_dir, working_dir, set_dirs, allow_cache)
      elif type(files_in) is list:
        for file_in in files_in:
          data = pwscf_exec.add_files(data, file_in, input_dir, output_dir, working_dir, set_dirs, allow_cache)

# Run
      pwscf_exec.execute_or_cache(data, log)

# Process Results
      data = pwscf_exec.results(data, log)

# Return
      return data
    except:
      return None

  @staticmethod
  def add_files(data, file_in, input_dir=None, output_dir=None, working_dir=None, set_dirs=True, allow_cache=True):   

# Get environment settings
    s = pwscf_settings.load()    

# INPUT FILE/PATH
    input_path = file_in
    if(input_dir != None):
      input_path = pwscf_exec.join_file_dir(file_in, input_dir)
    input_file = pwscf_exec.file_name_only(input_path)
    input_dir = pwscf_exec.dir_only(input_path)

# OUTPUT FILE/PATH
    output_file = pwscf_exec.file_name_only_no_ext(input_path) + ".out"
    if(output_dir != None):
      output_path = pwscf_exec.join_file_dir(output_file, output_dir)
    else:
      output_path = pwscf_exec.join_file_dir(output_file, input_dir)
    output_dir = pwscf_exec.dir_only(output_path)

# Absolute paths
    cwd = os.getcwd()
    input_path = pwscf_exec.join_file_dir(input_path, cwd)     
    output_path = pwscf_exec.join_file_dir(output_path, cwd)     

    if(working_dir == "None"):
      working_dir = cwd
    else:
      if(working_dir[0] != "/"):
        working_dir = pwscf_exec.join_file_dir(working_dir, cwd)       
    wd_ifp = pwscf_exec.join_file_dir(input_file, working_dir)    
    wd_ifn = pwscf_exec.file_name_only(wd_ifp)  
    wd_ifd = pwscf_exec.dir_only(wd_ifp)
 
# CREATE TEMPORARY INPUT FILE
    pw_in = pwscf_input()
    pw_in.load(input_path)

# SET DIRS
    if(set_dirs):
      if('pwscf_scratch' in s.keys() and 'pwscf_pp' in s.keys()):
        pw_in.set_dirs(s['pwscf_scratch'], s['pwscf_pp'])

# Get Signature
    sig = pw_in.signature()

# SAVE WORKING INPUT FILE
    pw_in.save(wd_ifn, wd_ifd)

# Make command
    proc_count = 1
    if('proc_count' in s.keys()):
      proc_count = s['proc_count']

    pwscf_bin = 1
    if('pwscf_bin' in s.keys()):
      pwscf_bin = s['pwscf_bin']
    
# MAKE COMMAND
    cmd = 'mpirun -n ' + proc_count + ' ' + pwscf_bin + ' -i ' + wd_ifp + ' > ' + output_path

# CACHE FILE
    pwscf_cache = None
    if('pwscf_cache' in s.keys()):
      pwscf_cache = pwscf_exec.abs_dir(s['pwscf_cache'])
      cache_file_in = pwscf_cache + "/" + sig + ".in"
      cache_file_out = pwscf_cache + "/" + sig + ".out"
      cache_exists = False

      if(os.path.exists(cache_file_in) and os.path.exists(cache_file_out)):
        cache_exists = True

# Set up dictionary
    entry = {
    'input_file': input_file, 
    'input_dir': input_dir, 
    'input_path': input_path, 
    'output_file': output_file, 
    'output_dir': output_dir,
    'output_path': output_path, 
    'wd': working_dir,
    'wd_ifn': wd_ifn,
    'wd_ifd': wd_ifd,
    'wd_ifp': wd_ifp,
    'signature': sig,
    'set_dirs': set_dirs,
    'allow_cache': allow_cache,
    'cache_dir': pwscf_cache,
    'cache_exists': cache_exists,
    'cache_file_in': cache_file_in,
    'cache_file_out': cache_file_out,
    'cmd': cmd,
    'calculation_successful': False,
    'total_energy': None,
    'energy_per_atom': None,
    'total_force': None,

    }
 
    data.append(entry)
    return data

  @staticmethod
  def execute_or_cache(data, log):
# Loop through entries
    for entry in data:
      if(entry['allow_cache'] and entry['cache_exists']):
        if(log != None):
          log.add("Cache used " + entry['cache_file_out'])
        copyfile(entry['cache_file_out'], entry['output_path'])
      else:
# Log
        if(log != None):
          log.add(entry['cmd'])

# Run
        os.system(entry['cmd'])

# Load file
        pwo = pwscf_output(entry['output_path'])
        if(pwo.get_job_done()):
          entry['calculation_successful'] = True

# Copy
        copyfile(entry['wd_ifp'], entry['cache_file_in'])
        copyfile(entry['output_path'], entry['cache_file_out'])
        if(log != None):
          log.add("Cached " + entry['cache_file_in'])
          log.add("Cached " + entry['cache_file_out'])
          if(not pwo.get_job_done()):
            log.add("PWscf Failed")

  @staticmethod
  def results(data, log):
# Loop through entries
    for entry in data:
      print(entry['output_path'])
      pwo = pwscf_output(entry['output_path'])
      if(pwo.get_job_done()):
        entry['total_energy'] = pwo.get_total_energy()
        entry['energy_per_atom'] = pwo.get_energy_per_atom()
        entry['total_force'] = pwo.get_total_force()
    return data

  @staticmethod
  def file_name_only(file_path):
    l1 = file_path.split("/")    
    return l1[-1]

  @staticmethod
  def file_name_only_no_ext(file_path):
    l1 = file_path.split("/")    
    l2 = l1[-1].split(".")
    file_name = l1[-1][:-(len(l2)+1)]
    return file_name
    
  @staticmethod
  def dir_only(file_path):
    l1 = file_path.split("/")   
    return file_path[:-(len(l1[-1])+1)]

  @staticmethod
  def join_file_dir(file_in, dir_in):
    a = dir_in
    if(dir_in[-1] == "/"):
      a = dir_in[:-2]
    b = file_in
    if(file_in[0] == "/"):
      b = file_in[1:]
    return a + "/" + b

  @staticmethod
  def abs_dir(dir_in):
    if(dir_in[0] != "/"):
      dir_in = join_file_dir(dir_in, os.getcwd())
    return dir_in

###########################################
#  CLASS lo
###########################################
class log:

  def __init__(self, file_name="log.txt"):
    self.data = []
    self.start_time = time.time()
    self.now = datetime.datetime.now()
    self.time_now = self.get_now()
    self.file_name = file_name

    fh = open(self.file_name, 'w')
    fh.write("####################################################################\n")  
    fh.write("                              LOG FILE                              \n")    
    fh.write("####################################################################\n")   
    fh.write("\n")    
    fh.write(self.time_now + "\n")     
    fh.write("\n")    
    fh.write("\n")   
    fh.close()

  def set_time(self, time_in):
    self.start_time = time_in

  def add(self, log_line):
    if(type(log_line) == list):
      for line in log_line:
        self.add(line)
    else:
      time_str = str(round(time.time() - self.start_time, 4))
      p = 11 - len(time_str)
      pad = "            "
      line = time_str + ":" + pad[0:p] + log_line
      fh = open(self.file_name, 'a')
      fh.write(line + "\n")    
      fh.close()

  def log(self, line):
    self.add(line)

  def output(self, file_name="log.txt"):
    pass

  def get_now(self):
    now = datetime.datetime.now()
    hour = str(now.hour)
    if(now.hour < 10):
      hour = "0" + hour
    minute = str(now.minute)
    if(now.minute < 10):
      minute = "0" + minute
    day = str(now.day)
    if(now.day < 10):
      day = "0" + day
    month = str(now.month)
    if(now.month < 10):
      month = "0" + month
    year = str(now.year)
    return hour + ":" + minute + "   " + day + "/" + month + "/" + year
    
###########################################
#  CLASS pwscf_setting
###########################################
class pwscf_settings:
 
  @staticmethod
  def load():
# Default Settings
        
    s = {
    "omp_num_threads":  1,
    "proc_count":  1,
    "pwscf_scratch": '',
    "pwscf_pp": '',
    "pwscf_cache": '',
    "pwscf_bin": '',
    }

    s['omp_num_threads'] = os.environ['OMP_NUM_THREADS']
    s['proc_count'] = os.environ['PROC_COUNT']
    s['pwscf_scratch'] = os.environ['PWSCF_SCRATCH']
    s['pwscf_pp'] = os.environ['PWSCF_PP']
    s['pwscf_cache'] = os.environ['PWSCF_CACHE']
    s['pwscf_bin'] = os.environ['PWSCF_BIN']

    return s

###########################################
#  CLASS lm
###########################################
class lma:

  @staticmethod
  def run(f, p, settings_in={}, data = None):
# Load Data
    data = lma.load(data)
  
# Default Settings
    settings = {
    'conv': 1.0E-7,
    'h': 0.00001,
    'lcut': 0.1,
    'outer_cycles': 50,
    'inner_cycles': 100,
    'covariance': False,
    }
  
# Load new settings
    for key in settings_in:
      if(key in settings):
        settings[key] = settings_in[key]
  
# calculate
    return lma.calc(f, p, settings, data)
    
  @staticmethod
  def calc(f, p, settings, data):
    return lma.outer_cycle(f, p, settings, data)
    
  @staticmethod
  def outer_cycle(f, p, settings, data):
# (JTJ+Lambda*diag(JTJ)) P = (-1*JTR)
# (H+Lambda*diag(H)) P = (-1*JTR)
    i = 0
    outer_cycles = settings['outer_cycles']
    is_converged = False
    lam = 0
    lam_cut = 0
    
    while(i < outer_cycles):
      i = i + 1
      
      R = lma.residual(f, p, data)        # R
      J = lma.jacobian(f, p, data, settings)        # J
      W = lma.covariance(R, settings['covariance'])
      JT = numpy.transpose(J)
      JTW = numpy.matmul(JT, W)
      JTWJ = numpy.matmul(JTW, J)
      JTWR = numpy.matmul(JTW, R)
      rss = lma.calc_rss_from_r(R)
      
      if(i == 1):
        lam = lma.calc_lambda(J)
        lam_cut = lam
# Run calc
      p_new, lam, new_rss = lma.inner_cycle(f, p, settings, data, J, JT, R, W, JTW, JTWJ, JTWR, lam, lam_cut, rss)
      
      if(new_rss < rss):
        p = p_new        
        is_converged = lma.converged(new_rss, rss, settings)
        rss = new_rss      
        if(is_converged):
          i = outer_cycles      
    return p, rss         
       
  @staticmethod
  def inner_cycle(f, p, settings, data, J, JT,  R, W, JTW, JTWJ, JTWR, lam, lam_cut, best_rss):
    inner_cycles = settings['inner_cycles']
    for n in range(inner_cycles):  
      try:
        p_new = lma.lma_calc(p, J, JT, R, W, JTW, JTWJ, JTWR, lam, lam_cut) 
        rss = lma.calc_rss(f, p_new, data)
        if(rss < best_rss):
          p = p_new
          lam = lam * 0.2
          best_rss = rss
        else:
          if(lam < lam_cut):
            lam = 2.0 * lam_cut
          lam = lam * 1.5
      except:    
        if(lam < lam_cut):
          lam = 2.0 * lam_cut
          lam = lam * 1.5 
    return p, lam, best_rss
  
  @staticmethod
  def converged(new_rss, last_rss, settings):    
    if(new_rss == 0.0 or abs(new_rss - last_rss) <= settings['conv']):
      return True
    return False
      
  @staticmethod
  def lma_calc(p, J, JT, R, W, JTW, JTWJ, JTWR, lam, lam_cut):
    JTWJ_diag = lam * numpy.diag(JTWJ)
    A = JTWJ + JTWJ_diag
    B = -1 * JTWR
    p_change = numpy.linalg.solve(A, B)
    return (p + p_change)
    
  @staticmethod
  def jacobian(f, p, data, settings):
    dl = len(data)
    pl = len(p)
    J = numpy.zeros((dl, pl), dtype=numpy.float64)    
    h = lma.step(p, settings['h'])
    for i in range(0, dl):
      for j in range(0, pl):
# Reset parameters
        p_new = numpy.copy(p)
# Vary jth parameter
        p_new[j] = p_new[j] + h[j]        
# Calc J matrix
        J[i,j] = (f(p_new, data[i, 0]) - data[i, 1]) / h[j]    
    return J
  
  @staticmethod
  def jacobian_transpose(J):
    return numpy.transpose(self.J)
    
  @staticmethod
  def hessian(J, JT):
    return numpy.matmul(JT, J)

  @staticmethod
  def dampening(l, H):
    damp = numpy.identity(len(H), dtype=numpy.float64)
    for i in range(0,len(H)):
      damp[i,i] = l * H[i,i]    
  
  @staticmethod
  def step(p, dh):
    h = numpy.zeros((len(p)), dtype=numpy.float64)  
    for i in range(len(h)):
      if(p[i] == 0.0):
        h[i] = dh
      else:
        h[i] = dh * p[i]
    return h
    
  @staticmethod
  def residual(f, p, data):
# Calculate residual
    return f(p, data[:, 0]) - data[:, 1] 
    
  def covariance(R, on=True):
    W = numpy.identity((len(R)), dtype=numpy.float64)
    if(on):
      build = True
      for i in range(len(R)):
        if(R[i] == 0.0):
          build = False
          break
      if(build):
        for i in range(len(R)):
          W[i,i] = 1.0 / R[i]**2
    return W
    
  @staticmethod
  def calc_lambda(J):
    JT = numpy.transpose(J)
    JTJ = numpy.matmul(JT, J)
    JTJ_inv = numpy.linalg.inv(JTJ)
    trace = numpy.trace(JTJ_inv)
    l = trace**(-1)
    return l

  @staticmethod
  def calc_rss(f, p, data):
    return sum((f(p, data[:, 0]) - data[:, 1])**2)

  @staticmethod
  def calc_rss_from_r(R):
    return sum(R**2)
  
###################################
#  Load Data
###################################

  @staticmethod
  def load(file_in=None):
    try:
      if(not isinstance(file_in, (list,)) and not isinstance(file_in,numpy.ndarray) and file_in == None):
        if(len(sys.argv) > 1 and sys.argv[1] is not None):
          file_in = sys.argv[1]
          if(not os.path.isfile(file_in)):
            return None
        else:
          return None 
    except:
      pass

# List
    if(isinstance(file_in, (list,))):
      if(len(file_in) > 2 or (len(file_in) == 2 and len(file_in[0]) == 2)):
        data = numpy.zeros((len(file_in), 2), dtype=numpy.float64) 
        for row in range(len(file_in)):
          data[row,0] = float(file_in[row][0])
          data[row,1] = float(file_in[row][1])
        return data
      if(len(file_in) == 2 and len(file_in[0]) > 2):
        data = numpy.zeros((len(file_in[0]), 2), dtype=numpy.float64) 
        for row in range(len(file_in[0])):
          data[row,0] = float(file_in[0][row])
          data[row,1] = float(file_in[1][row])
        return data
# If array
    elif(isinstance(file_in,numpy.ndarray)): 
      data = numpy.zeros((len(file_in), 2), dtype=numpy.float64) 
      for row in range(len(file_in)):
        data[row,0] = float(file_in[row,0])
        data[row,1] = float(file_in[row,1])
      return data    
# If file
    elif(os.path.isfile(file_in)):
# Init variable
      file_data = ""

# Read it in line by line
      fh = open(file_in, "r")
      for file_row in fh:
        if(file_row.strip() != ""):
          file_data = file_data + file_row.strip() + '\n'
# Clean
      file_data = lma.clean(file_data)
      
# Read in data
      lines = file_data.split('\n')
  
# Load to list first
      data_list = []
      for line in lines:
        fields = line.split(',')
        if(len(fields) == 2):
          try: 
            data_list.append([float(fields[0]), float(fields[1])])
          except:
            pass
            
# Define and load into array
      data = numpy.zeros((len(data_list), 2), dtype=numpy.float64) 
      for row in range(len(data_list)):
        data[row,0] = float(data_list[row][0])
        data[row,1] = float(data_list[row][1])
      return data

  @staticmethod
  def clean(str_in):  
    str_out = ""
    l = len(str_in)
    for i in range(l):
# Last, Next, This
      if(i == 0):
        last = None
      else:
        last = str_in[i-1]
      if(i < (l-1)):
        next = str_in[i+1]
      else:  
        next = None
      char = str_in[i]
      
# Tabs and pipes to commas
      if(char == "\t"):
        char = ","
      if(char == "|"):
        char = ","
    
# Check
      ok = True
      if(last == " " and char == " "):
        ok = False
      elif(last == "\n" and char == "\n"):
        ok = False
      elif(last == "\n" and char == " "):
        ok = False
      elif(char == " " and next == "\n"):
        ok = False
      elif(char == "," and next == ","):
        ok = False

# Add to string
      if(ok):
        str_out += char
    return str_out      
  
###########################################
#  CLASS eos_fittin
###########################################
class eos_fitting:

  @staticmethod
  def bm(data=None):
    return eos_fitting.eos_fit(eos_fitting.bm_calc, data)
    
  def pt(data=None):
    return eos_fitting.eos_fit(eos_fitting.pt_calc, data)
    
  def vinet(data=None):
    return eos_fitting.eos_fit(eos_fitting.vinet_calc, data)
    
  def all(data=None):
    data = eos_fitting.load(data)
    x, y = eos_fitting.eos_fit(eos_fitting.bm_calc, data) 
    plt.clf()
    plt.plot(data[:,0], data[:,1], 'r+', label="E-V")
    plt.plot(x, y, 'b-', label="E-V fit")
    plt.savefig('1.eps', format='eps')  
    x, y = eos_fitting.eos_fit(eos_fitting.bm_calc, data) 
    plt.clf()
    plt.plot(data[:,0], data[:,1], 'r+', label="E-V")
    plt.plot(x, y, 'b-', label="E-V fit")
    plt.savefig('2.eps', format='eps')  
    x, y = eos_fitting.eos_fit(eos_fitting.bm_calc, data) 
    plt.clf()
    plt.plot(data[:,0], data[:,1], 'r+', label="E-V")
    plt.plot(x, y, 'b-', label="E-V fit")
    plt.savefig('3.eps', format='eps')  

  @staticmethod
  def eos_fit(f, data=None, verbose=True):
    data = eos_fitting.load(data)
    V = data[:,0]
    E = data[:,1]
  
    eos = {"E0": 0.0, "V0": 0.0, "B0": 0.0, "B0P": 0.0}

# 2nd order polynomial fit
    poly = numpy.polyfit(V, E, 2)
# Starting points
    eos['V0'] = (-1 * poly[1]) / (2 * poly[0])
    eos['E0'] = (poly[0] * eos['V0'] * eos['V0']) + (poly[1] * eos['V0']) + poly[2]
    eos['B0'] = 2.0 * poly[0] * eos['V0']
    eos['B0P'] = 2.0
    
    p = numpy.zeros((4))
    p[0] = eos['V0']
    p[1] = eos['E0']
    p[2] = eos['B0']
    p[3] = eos['B0P']
    
    print(p)

    rss = fopt.calc_rss(f, p, data)
    print(rss)
#p, rss = fopt.run(f, p, {}, data)

    datafit = numpy.zeros((101,2))
    datafit[:,0] = numpy.linspace(data[0,0],data[-1,0], 101)
    datafit[:,1] = f(p, datafit[:,0])   
    return p, rss, datafit
   
  @staticmethod
  def bm_calc(p, V):
    V0 = p[0]
    E0 = p[1]
    B0 = p[2]
    B0P = p[3]
    eta = (V/V0)**(1/3.0)
    return E0 + (9/16.0) * (B0 * V0) * ((eta*eta - 1)*(eta*eta - 1)) * (6.0 + B0P * (eta * eta - 1) - 4 * eta * eta ) 
    
  @staticmethod
  def pt_calc(p, V):
    V0 = p[0]
    E0 = p[1]
    B0 = p[2]
    B0P = p[3]
    eta = (V/V0)
    return E0 + ((B0 * V0) / 2) * (numpy.ln(eta))**2 + ((B0 * V0) / 6) * (numpy.ln(eta))**3 * (B0P - 2)
    
  @staticmethod
  def vinet_calc(p, V):
    V0 = p[0]
    E0 = p[1]
    B0 = p[2]
    B0P = p[3]
    eta = (V/V0)**(1/3.0)
    return E0 + ((2 * B0 * V0) / (B0P - 1)**2) * (2 - (5 + 3 * eta * (B0P - 1) - 3 * B0P) * numpy.exp(-(3.0/2.0) * (B0P - 1) * (eta - 1)))

###################################
#  Load Data
###################################

#   Loads data into numpy array
#   Tries to read from file if file_in is set
#   If file_in is a list or numpy array it reads it in directly
#   If file_in is none, it checks if a file has been specified as an argument

  @staticmethod
  def load(file_in=None):
    if(not isinstance(file_in, (list,)) and not isinstance(file_in,numpy.ndarray) and file_in == None):
      if(len(sys.argv) > 1 and sys.argv[1] is not None):
        file_in = sys.argv[1]
        if(not os.path.isfile(file_in)):
          return None
      else:
        return None 

# List
    if(isinstance(file_in, (list,))):
      if(len(file_in) > 2 or (len(file_in) == 2 and len(file_in[0]) == 2)):
        data = numpy.zeros((len(file_in), 2)) 
        for row in range(len(file_in)):
          data[row,0] = float(file_in[row][0])
          data[row,1] = float(file_in[row][1])
        return data
      if(len(file_in) == 2 and len(file_in[0]) > 2):
        data = numpy.zeros((len(file_in[0]), 2)) 
        for row in range(len(file_in[0])):
          data[row,0] = float(file_in[0][row])
          data[row,1] = float(file_in[1][row])
        return data
# If array
    elif(isinstance(file_in,numpy.ndarray)): 
      data = numpy.zeros((len(file_in), 2)) 
      for row in range(len(file_in)):
        data[row,0] = float(file_in[row,0])
        data[row,1] = float(file_in[row,1])
      return data    
# If file
    elif(os.path.isfile(file_in)):
# Init variable
      file_data = ""

# Read it in line by line
      fh = open(file_in, "r")
      for file_row in fh:
        if(file_row.strip() != ""):
          file_data = file_data + file_row.strip() + '\n'
# Clean
      file_data = eos_fitting.clean(file_data)
      
# Read in data
      lines = file_data.split('\n')
  
# Load to list first
      data_list = []
      for line in lines:
        fields = line.split(',')
        if(len(fields) == 2):
          try: 
            data_list.append([float(fields[0]), float(fields[1])])
          except:
            pass
            
# Define and load into array
      data = numpy.zeros((len(data_list), 2)) 
      for row in range(len(data_list)):
        data[row,0] = float(data_list[row][0])
        data[row,1] = float(data_list[row][1])
      return data      

  @staticmethod
  def clean(str_in):  
    str_out = ""
    l = len(str_in)
    for i in range(l):
# Last, Next, This
      if(i == 0):
        last = None
      else:
        last = str_in[i-1]
      if(i < (l-1)):
        next = str_in[i+1]
      else:  
        next = None
      char = str_in[i]
      
# Tabs and pipes to commas
      if(char == "\t"):
        char = ","
      if(char == "|"):
        char = ","
    
# Check
      ok = True
      if(last == " " and char == " "):
        ok = False
      elif(last == "\n" and char == "\n"):
        ok = False
      elif(last == "\n" and char == " "):
        ok = False
      elif(char == " " and next == "\n"):
        ok = False
      elif(char == "," and next == ","):
        ok = False

# Add to string
      if(ok):
        str_out += char
    return str_out  

###########################################
#  CLASS fop
###########################################
class fopt:

  @staticmethod
  def run(f, p, settings = {}, data=None):
    data = fopt.load()
    print(data)
    n = 0
    while(n<5):
      try: 
        return fopt.action(f, p, settings, data)
      except:
        n = n + 1
  
  @staticmethod
  def action(f, p, settings = {}, data=None):
    default = {
    'ga': True,
    'sa': True,
    'lma': True,
    'ga_bounds': 25.0,
    'pl': None,
    'pu': None,
    'ga_poolsize': 100,
    'ga_poolsamples': 250,
    'ga_loops': 1000,
    'sa_tempmin': 0.1,
    'sa_tempmax': 10.0,
    'sa_factor': 0.9,
    'sa_loops': 100,
    'lma_inner': 50,
    'lma_outer': 100
    }
  
    for key in settings:
      if(key in default):
        default[key] = settings[key]
    print("GA start")
# Step 1 - ga
    fit_ga = genetic(data)
    fit_ga.set_f(f)
    fit_ga.set_p(p)
    if(default['pl'] != None and default['pu'] != None):
      fit_ga.set_bounds_lu(default['pl'], default['pu'])
    else:    
      fit_ga.set_bounds(default['ga_bounds'])
    fit_ga.set_pool(default['ga_poolsize'], default['ga_poolsamples'])
    fit_ga.set_loops(default['ga_loops'])
    p, rss = fit_ga.run()
    pl, pu = fit_ga.get_lu()
    print("GA: ", rss)
    
# Step 2 - sa
    fit_sa = simulated_annealing(data)
    fit_sa.set_f(f)
    fit_sa.set_settings({"temp_max": default['sa_tempmax'], 
                         "temp_min": default['sa_tempmin'], 
                         "factor": default['sa_factor'], 
                         "iterations": default['sa_loops']})
    fit_sa.set_p(p)
    fit_sa.set_pl(pl)
    fit_sa.set_pu(pu)
    p, rss = fit_sa.run()
    print("SA: ", rss)

# Step 3 - lma
#pb, rss = ng.run(f, p, {}, data)
#print("ng: ", rss)
    """
    fit_lma = lma(data)
    fit_lma.set_fit(f, p)
    fit_lma.set_cycles(default['lma_inner'], default['lma_outer'])
    p, rss = fit_lma.calc()
    print("LMA: ", rss)
    """
    
    p, rss = lma.run(f, p, {}, data)
    print("LMA: ", rss)
    
    return p, rss
    
###################################
#  Calc RSS
###################################

  @staticmethod
  def calc_rss(f, p, data):
# If there's an error, just set as a very high rss to avoid program crashing
    try: 
      rss = sum((data[:,1] - f(p, data[:,0]))**2) 
    except:
      rss = 1.0e99
    return rss     
    
###################################
#  Load Data
###################################

  @staticmethod
  def load(file_in=None):
    try:
      if(file_in == None):
        if(len(sys.argv) > 1 and sys.argv[1] is not None):
          file_in = sys.argv[1]
          if(not os.path.isfile(file_in)):
            return None
        else:
          return None 
    except:
      pass

# List
    if(isinstance(file_in, (list,))):
      if(len(file_in) > 2 or (len(file_in) == 2 and len(file_in[0]) == 2)):
        data = numpy.zeros((len(file_in), 2)) 
        for row in range(len(file_in)):
          data[row,0] = float(file_in[row][0])
          data[row,1] = float(file_in[row][1])
        return data
      if(len(file_in) == 2 and len(file_in[0]) > 2):
        data = numpy.zeros((len(file_in[0]), 2)) 
        for row in range(len(file_in[0])):
          data[row,0] = float(file_in[0][row])
          data[row,1] = float(file_in[1][row])
        return data
# If array
    elif(isinstance(file_in,numpy.ndarray)): 
      data = numpy.zeros((len(file_in), 2)) 
      for row in range(len(file_in)):
        data[row,0] = float(file_in[row,0])
        data[row,1] = float(file_in[row,1])
      return data    
# If file
    elif(os.path.isfile(file_in)):
# Init variable
      file_data = ""

# Read it in line by line
      fh = open(file_in, "r")
      for file_row in fh:
        if(file_row.strip() != ""):
          file_data = file_data + file_row.strip() + '\n'
# Clean
      file_data = fopt.clean(file_data)
      
# Read in data
      lines = file_data.split('\n')
  
# Load to list first
      data_list = []
      for line in lines:
        fields = line.split(',')
        if(len(fields) == 2):
          try: 
            data_list.append([float(fields[0]), float(fields[1])])
          except:
            pass
            
# Define and load into array
      data = numpy.zeros((len(data_list), 2)) 
      for row in range(len(data_list)):
        data[row,0] = float(data_list[row][0])
        data[row,1] = float(data_list[row][1])
      return data

  @staticmethod
  def clean(str_in):  
    str_out = ""
    l = len(str_in)
    for i in range(l):
# Last, Next, This
      if(i == 0):
        last = None
      else:
        last = str_in[i-1]
      if(i < (l-1)):
        next = str_in[i+1]
      else:  
        next = None
      char = str_in[i]
      
# Tabs and pipes to commas
      if(char == "\t"):
        char = ","
      if(char == "|"):
        char = ","
    
# Check
      ok = True
      if(last == " " and char == " "):
        ok = False
      elif(last == "\n" and char == "\n"):
        ok = False
      elif(last == "\n" and char == " "):
        ok = False
      elif(char == " " and next == "\n"):
        ok = False
      elif(char == "," and next == ","):
        ok = False

# Add to string
      if(ok):
        str_out += char
    return str_out      
  
###########################################
#  CLASS geneti
###########################################
class genetic:

  @staticmethod
  def run(f, p_in, settings_in={}, data = None):
# Load Data
    data = genetic.load(data)
  
# Default Settings
    settings = {
    'p_diff': None,
    'p_percent': 50,
    'pool_size': 200,
    'pool_seed_factor': 3,
    'breed_events': 2,
    'generations': 10,
    'mutations': 10,
    'reseeds': 10
    }
    
# Load settings
    settings = genetic.load_settings(settings_in, settings)    
    
    p = numpy.asarray(p_in)
    pl, pu = genetic.set_bounds(p, settings)

    pool = genetic.init_pool(f, pl, pu, settings, data)
    
    genetic.evolve(f, pl, pu, settings, data, pool)
    return pool['p'][0,:], pool['rss'][0]
    
  @staticmethod
  def load_settings(settings_in, settings): 
# Load new settings
    for key in settings_in:
      if(key in settings):
        settings[key] = settings_in[key]
    return settings
    
  @staticmethod
  def set_bounds(p, settings):  
    if(settings['p_diff'] != None):
      p_diff = numpy.asarray(settings['p_diff'])
      pl = numpy.copy(p) - p_diff
      pu = numpy.copy(p) + p_diff  
    else:
      pl = numpy.copy(p) * (1 - settings['p_percent']/100)
      pu = numpy.copy(p) * (1 + settings['p_percent']/100)
    return pl, pu
    
  @staticmethod
  def make_pool(pool_size, param_size):
    pool = {}
    pool['p'] = numpy.zeros((pool_size, param_size))
    pool['rss'] = numpy.zeros((pool_size))
    pool['rss'][:] = -1.0    
    return pool
    
  @staticmethod
  def copy_pool(pool_in):
    pool_out = {}
    pool_out['p'] = numpy.copy(pool_in['p']) 
    pool_out['rss'] = numpy.copy(pool_in['rss']) 
    return pool_out  

  @staticmethod
  def init_pool(f, pl, pu, settings, data):  
    pool = genetic.make_pool(settings['pool_size'], len(pl)) 
    n = 0 
    while(pool['rss'][-1] == -1.0):
      p = genetic.seed(pl, pu)
      rss = genetic.calc_rss(f, p, data)
      genetic.add_to_pool(pool, p, rss)
      n = n + 1  
  
    total_loops = int(settings['pool_seed_factor'] * n)  
    while(n <= total_loops):
      p = genetic.seed(pl, pu)
      rss = genetic.calc_rss(f, p, data)
      genetic.add_to_pool(pool, p, rss)
      n = n + 1   
    return pool

  @staticmethod
  def add_to_pool(pool, p, rss):
    if(genetic.in_pool(pool, p)):
      return 0
    pool_len = len(pool['rss'])
# If worse than last entry, just leave
    if(rss > pool['rss'][-1] and pool['rss'][-1] != -1.0):
      return 0
    
    if(pool['rss'][0] == -1.0):
      pool['rss'][0] = rss
      pool['p'][0] = p
      return 0     
    
# Find where to insert
    for n in range(pool_len):
      if(rss < pool['rss'][n] or pool['rss'][n] == -1.0):
        pool['rss'] = numpy.insert(pool['rss'], n, rss)[0:pool_len]
        pool['p'] = numpy.insert(pool['p'], n, p, 0)[0:pool_len,:]
        break
    
  @staticmethod
  def in_pool(pool, p):
    for i in range(len(pool['p'])):
      for j in range(len(p)):
        same = True
        if(p[j] != pool['p'][i, j]):
          same = False
          break
      if(same == True):
        return True
    return False
  
  @staticmethod
  def seed(pl, pu):
    return pl[:] + (pu[:] - pl[:]) * numpy.random.rand(len(pl))
    
  @staticmethod
  def evolve(f, pl, pu, settings, data, pool):
    pool_len = len(pool['p'])  
    for g in range(settings['generations']):
      pool = genetic.next_generation(f, pl, pu, settings, data, pool)
      pool = genetic.mutate_pool(f, pl, pu, settings, data, pool)
      pool = genetic.reseed_pool(f, pl, pu, settings, data, pool)
    return pool   
    
  def next_generation(f, pl, pu, settings, data, pool):
      new_pool = genetic.copy_pool(pool) 
      breed_events = int(settings['breed_events'] * len(pool['p']))
      for n in range(breed_events):
        p = genetic.breed(pool, pl, pu)
        rss = genetic.calc_rss(f, p, data)
        genetic.add_to_pool(new_pool, p, rss)  
      pool = genetic.copy_pool(new_pool)
      return pool
      
  def mutate_pool(f, pl, pu, settings, data, pool):
    dp = pu - pl 
    pool_len = len(pool['p'])
    for m in range(settings['mutations']):
      n = numpy.random.randint(0, pool_len) 
      p = genetic.mutate(pool['p'][n,:], dp)
      if(not genetic.in_pool(pool, p)):
        rss = genetic.calc_rss(f, p, data)
        genetic.add_to_pool(pool, p, rss)
    return pool  
      
  def mutate(p, dp):
    return p[:] + (dp) * (0.5 - numpy.random.rand(len(p)))
    
  def reseed_pool(f, pl, pu, settings, data, pool):
    n = 0
    while(n < settings['reseeds']):
      p = genetic.seed(pl, pu)
      if(not genetic.in_pool(pool, p)):
        rss = genetic.calc_rss(f, p, data)
        genetic.add_to_pool(pool, p, rss)
        n = n + 1  
    return pool    
    
  @staticmethod
  def breed(pool, pl, pu, mutate=0.05): 
    pool_len = len(pool['p'])
    p_len = len(pool['p'][0])    
    
    if(numpy.random.uniform(0.0,1.0) <= mutate):
      return genetic.seed(pl, pu)
    else:
      p_x = numpy.random.randint(0, pool_len)
      p_y = p_x
      while(p_y == p_x):
        p_y = numpy.random.randint(0, pool_len)      
      p = numpy.copy(pool['p'][p_x])
      for i in range(0, p_len):
        if(numpy.random.randint(0, 99) < 50):
          p[i] =  pool['p'][p_y, i]
      if(genetic.in_pool(pool, p)):
# No clones
        return genetic.seed(pl, pu)
      return p
   
  @staticmethod
  def calc_rss(f, p, data):
    try: 
      rss = sum((data[:,1] - f(p, data[:,0]))**2) 
      return rss
    except:
      return 1.0e99

###################################
#  Load Data
###################################

  @staticmethod
  def load(file_in=None):
    try:
      if(file_in == None):
        if(len(sys.argv) > 1 and sys.argv[1] is not None):
          file_in = sys.argv[1]
          if(not os.path.isfile(file_in)):
            return None
        else:
          return None 
    except:
      pass

# List
    if(isinstance(file_in, (list,))):
      if(len(file_in) > 2 or (len(file_in) == 2 and len(file_in[0]) == 2)):
        data = numpy.zeros((len(file_in), 2)) 
        for row in range(len(file_in)):
          data[row,0] = float(file_in[row][0])
          data[row,1] = float(file_in[row][1])
        return data
      if(len(file_in) == 2 and len(file_in[0]) > 2):
        data = numpy.zeros((len(file_in[0]), 2)) 
        for row in range(len(file_in[0])):
          data[row,0] = float(file_in[0][row])
          data[row,1] = float(file_in[1][row])
        return data
# If array
    elif(isinstance(file_in,numpy.ndarray)): 
      data = numpy.zeros((len(file_in), 2)) 
      for row in range(len(file_in)):
        data[row,0] = float(file_in[row,0])
        data[row,1] = float(file_in[row,1])
      return data    
# If file
    elif(os.path.isfile(file_in)):
# Init variable
      file_data = ""

# Read it in line by line
      fh = open(file_in, "r")
      for file_row in fh:
        if(file_row.strip() != ""):
          file_data = file_data + file_row.strip() + '\n'
# Clean
      file_data = genetic.clean(file_data)
      
# Read in data
      lines = file_data.split('\n')
  
# Load to list first
      data_list = []
      for line in lines:
        fields = line.split(',')
        if(len(fields) == 2):
          try: 
            data_list.append([float(fields[0]), float(fields[1])])
          except:
            pass
            
# Define and load into array
      data = numpy.zeros((len(data_list), 2)) 
      for row in range(len(data_list)):
        data[row,0] = float(data_list[row][0])
        data[row,1] = float(data_list[row][1])
      return data

  @staticmethod
  def clean(str_in):  
    str_out = ""
    l = len(str_in)
    for i in range(l):
# Last, Next, This
      if(i == 0):
        last = None
      else:
        last = str_in[i-1]
      if(i < (l-1)):
        next = str_in[i+1]
      else:  
        next = None
      char = str_in[i]
      
# Tabs and pipes to commas
      if(char == "\t"):
        char = ","
      if(char == "|"):
        char = ","
    
# Check
      ok = True
      if(last == " " and char == " "):
        ok = False
      elif(last == "\n" and char == "\n"):
        ok = False
      elif(last == "\n" and char == " "):
        ok = False
      elif(char == " " and next == "\n"):
        ok = False
      elif(char == "," and next == ","):
        ok = False

# Add to string
      if(ok):
        str_out += char
    return str_out      
  
###########################################
#  CLASS simulated_annealin
###########################################
class simulated_annealing:

  @staticmethod
  def run(f, p_in, settings_in={}, data = None):
# Load Data
    data = simulated_annealing.load(data)
  
# Default Settings
    settings = {
    'p_diff': None,
    'p_percent': 50,
    "temp_max": 10, 
    "temp_min": 0.1,
    "temp_min_loops": 20,
    "factor": 0.90,
    "iterations": 1000,
    "p_factor": 0.90
    }
    
# Load settings
    settings = simulated_annealing.load_settings(settings_in, settings)    
    
    p = numpy.asarray(p_in)
    pl, pu = simulated_annealing.set_bounds(p, settings)

    p, rss = simulated_annealing.solve(f, pl, pu, settings, data) 
     
    return p, rss
#return pool['p'][0,:], pool['rss'][0]
    
  @staticmethod
  def load_settings(settings_in, settings): 
# Load new settings
    for key in settings_in:
      if(key in settings):
        settings[key] = settings_in[key]
    return settings
        
  @staticmethod
  def set_bounds(p, settings):  
    if(settings['p_diff'] != None):
      p_diff = numpy.asarray(settings['p_diff'])
      pl = numpy.copy(p) - p_diff
      pu = numpy.copy(p) + p_diff  
    else:
      pl = numpy.copy(p) * (1 - settings['p_percent']/100)
      pu = numpy.copy(p) * (1 + settings['p_percent']/100)
    return pl, pu
        
  @staticmethod
  def solve(f, pl, pu, settings, data): 
    t = settings["temp_max"]
    rss_best = None
    p_best = None
    
# Initial parameters
    dp = pu[:] - pl[:]
    p = pl[:] + 0.5 * dp[:] 
    p_best = numpy.copy(p)
    rss_best = simulated_annealing.calc_rss(f, p_best, data)
    
    t_loop = True
    n = 0
    while(t_loop):
      n = n + 1
      rss_best, p_best = simulated_annealing.solve_loop(f, p, dp, t, settings, data, p_best, rss_best)
      t = t * settings["factor"]
      dp = dp * settings["p_factor"]
      
      if(t < settings["temp_min"]):
        t_loop = False
      if(settings["temp_min_loops"] != None and n < settings["temp_min_loops"]):
        t_loop = True
      
# Return
    return p_best, rss_best  
    
  @staticmethod 
  def solve_loop(f, p, dp, t, settings, data, p_best, rss_best):   
# Loop through temperature reductions
    n = 0
    p = p_best
    rss = rss_best
    bad = 0
    while(n<settings["iterations"]):
      p_trial = simulated_annealing.trial(p, dp)
      rss_trial = simulated_annealing.calc_rss(f, p_trial, data)
      if(rss_trial < rss):
        rss = rss_trial
        p = numpy.copy(p_trial)
        if(rss_best == None or rss < rss_best): 
          rss_best = rss
          p_best = numpy.copy(p)
      elif(numpy.exp((rss_best - rss_trial) / t) > random.random()):  
        rss = rss_trial
        p = numpy.copy(p_trial)
        bad = bad + 1
      n = n + 1
    return rss_best, p_best

  @staticmethod
  def trial(p, dp):
    r = numpy.random.rand(len(p))
    return p + (r - 0.5) * dp[:]
    
  @staticmethod
  def calc_rss(f, p, data):
    try: 
      rss = sum((data[:,1] - f(p, data[:,0]))**2) 
      return rss
    except:
      return 1.0e99

###################################
#  Load Data
###################################

  @staticmethod
  def load(file_in=None):
    try:
      if(file_in == None):
        if(len(sys.argv) > 1 and sys.argv[1] is not None):
          file_in = sys.argv[1]
          if(not os.path.isfile(file_in)):
            return None
        else:
          return None 
    except:
      pass

# List
    if(isinstance(file_in, (list,))):
      if(len(file_in) > 2 or (len(file_in) == 2 and len(file_in[0]) == 2)):
        data = numpy.zeros((len(file_in), 2)) 
        for row in range(len(file_in)):
          data[row,0] = float(file_in[row][0])
          data[row,1] = float(file_in[row][1])
        return data
      if(len(file_in) == 2 and len(file_in[0]) > 2):
        data = numpy.zeros((len(file_in[0]), 2)) 
        for row in range(len(file_in[0])):
          data[row,0] = float(file_in[0][row])
          data[row,1] = float(file_in[1][row])
        return data
# If array
    elif(isinstance(file_in,numpy.ndarray)): 
      data = numpy.zeros((len(file_in), 2)) 
      for row in range(len(file_in)):
        data[row,0] = float(file_in[row,0])
        data[row,1] = float(file_in[row,1])
      return data    
# If file
    elif(os.path.isfile(file_in)):
# Init variable
      file_data = ""

# Read it in line by line
      fh = open(file_in, "r")
      for file_row in fh:
        if(file_row.strip() != ""):
          file_data = file_data + file_row.strip() + '\n'
# Clean
      file_data = simulated_annealing.clean(file_data)
      
# Read in data
      lines = file_data.split('\n')
  
# Load to list first
      data_list = []
      for line in lines:
        fields = line.split(',')
        if(len(fields) == 2):
          try: 
            data_list.append([float(fields[0]), float(fields[1])])
          except:
            pass
            
# Define and load into array
      data = numpy.zeros((len(data_list), 2)) 
      for row in range(len(data_list)):
        data[row,0] = float(data_list[row][0])
        data[row,1] = float(data_list[row][1])
      return data

  @staticmethod
  def clean(str_in):  
    str_out = ""
    l = len(str_in)
    for i in range(l):
# Last, Next, This
      if(i == 0):
        last = None
      else:
        last = str_in[i-1]
      if(i < (l-1)):
        next = str_in[i+1]
      else:  
        next = None
      char = str_in[i]
      
# Tabs and pipes to commas
      if(char == "\t"):
        char = ","
      if(char == "|"):
        char = ","
    
# Check
      ok = True
      if(last == " " and char == " "):
        ok = False
      elif(last == "\n" and char == "\n"):
        ok = False
      elif(last == "\n" and char == " "):
        ok = False
      elif(char == " " and next == "\n"):
        ok = False
      elif(char == "," and next == ","):
        ok = False

# Add to string
      if(ok):
        str_out += char
    return str_out      
  
###########################################
#  CLASS n
###########################################
class ng:

  @staticmethod
  def run(f, p, settings_in={}, data = None, fp = None):
# Load Data
    data = ng.load(data)
  
# Default Settings
    settings = {
    'conv': 1.0e-9,
    'h': 0.00001,
    'lcut': 0.1,
    'cycles': 100,
    'differentiator': 'cd'   
    }
  
# Load new settings
    for key in settings_in:
      if(key in settings):
        settings[key] = settings_in[key]
  
# calculate
    return ng.calc(f, p, settings, data, fp)
    
  @staticmethod
  def calc(f, p, settings, data, fp):
    i = 0
    cycles = settings['cycles']
    lam = 0
    lam_cut = 0
    
    while(i < cycles):
      i = i + 1
      R = ng.residual(f, p, data)        # R
      J = ng.jacobian(f, p, data, settings, fp)        # J
      JT = numpy.transpose(J)
      JTJ = numpy.matmul(JT, J)
      JTR = numpy.matmul(JT, R)      
      rss = ng.calc_rss_from_r(R)
      print(rss)
      if(i > 1 and (rss > rss_last or abs(rss - rss_last) <= settings['conv'])):
        return p_last, rss_last  
      p_change = numpy.linalg.solve(JTJ, -JTR)
      p = p + p_change
      p_last = numpy.copy(p)
      rss_last = rss
      
    return p, rss         

  """
    
  @staticmethod
  def outer_cycle(f, p, settings, data):
# (JTJ+Lambda*diag(JTJ)) P = (-1*JTR)
# (H+Lambda*diag(H)) P = (-1*JTR)
    i = 0
    outer_cycles = settings['outer_cycles']
    is_converged = False
    lam = 0
    lam_cut = 0
    
    while(i < outer_cycles and is_converged == False):
      i = i + 1
      R = lma.residual(f, p, data)        # R
      J = lma.jacobian(f, p, data, settings)        # J
      W = lma.covariance(R, settings['covariance'])
      JT = numpy.transpose(J)
      JTW = numpy.matmul(JT, W)
      JTWJ = numpy.matmul(JTW, J)
      JTWR = numpy.matmul(JTW, R)
      rss = lma.calc_rss_from_r(R)
      
      if(i == 1):
        lam = lma.calc_lambda(J)
        lam_cut = lam
      p_new, lam, new_rss = lma.inner_cycle(f, p, settings, data, J, JT, R, W, JTW, JTWJ, JTWR, lam, lam_cut, rss)
      
      if(new_rss < rss):
        p = p_new        
        is_converged = lma.converged(new_rss, rss, settings)
        rss = new_rss
    return p, rss
  
  @staticmethod
  def lma_calc(p, J, JT, R, W, JTW, JTWJ, JTWR, lam, lam_cut):
    JTWJ_diag = lam * numpy.diag(JTWJ)
    A = JTWJ + JTWJ_diag
    B = -1 * JTWR
    print(A)
    print(B)
    print()
    p_change = numpy.linalg.solve(A, B)
    return (p + p_change)
  """
    
  @staticmethod
  def jacobian(f, p, data, settings, fp):
    if(fp == None and settings['differentiator'] == 'fd'):
      dl = len(data)
      pl = len(p)
      J = numpy.zeros((dl, pl), dtype=numpy.float64)    
      h = ng.step(p, settings['h'])
      for i in range(0, dl):
        for j in range(0, pl):
# Reset parameters
          p_new = numpy.copy(p)
# Vary jth parameter
          p_new[j] = p_new[j] + h[j]        
# Calc J matrix
          J[i,j] = (f(p_new, data[i, 0]) - data[i, 1]) / h[j]    
    elif(fp == None and settings['differentiator'] == 'cd'):
      dl = len(data)
      pl = len(p)
      J = numpy.zeros((dl, pl), dtype=numpy.float64)    
      h = ng.step(p, settings['h'])
      for i in range(0, dl):
        for j in range(0, pl):
# Reset parameters
          p_b = numpy.copy(p)
          p_f = numpy.copy(p)
# Vary jth parameter
          p_b[j] = p_b[j] - 0.5 * h[j]  
          p_f[j] = p_f[j] + 0.5 * h[j]        
# Calc J matrix
          J[i,j] = (f(p_f, data[i, 0]) - f(p_b, data[i, 0])) / h[j]    
    else:
      dl = len(data)
      pl = len(p)
      J = numpy.zeros((dl, pl), dtype=numpy.float64)
      for i in range(0, dl):
# Calc J matrix
        J[i,:] = fp(p, data[i, 0])    
    return J
  
  @staticmethod
  def jacobian_transpose(J):
    return numpy.transpose(self.J)
    
  @staticmethod
  def hessian(J, JT):
    return numpy.matmul(JT, J)

  @staticmethod
  def dampening(l, H):
    damp = numpy.identity(len(H), dtype=numpy.float64)
    for i in range(0,len(H)):
      damp[i,i] = l * H[i,i]    
  
  @staticmethod
  def step(p, dh):
    h = numpy.zeros((len(p)), dtype=numpy.float64)  
    for i in range(len(h)):
      if(p[i] == 0.0):
        h[i] = dh
      else:
        h[i] = dh * p[i]
    return h
    
  @staticmethod
  def residual(f, p, data):
# Calculate residual
    return f(p, data[:, 0]) - data[:, 1] 
    
  def covariance(R, on=True):
    W = numpy.identity((len(R)), dtype=numpy.float64)
    if(on):
      build = True
      for i in range(len(R)):
        if(R[i] == 0.0):
          build = False
          break
      if(build):
        for i in range(len(R)):
          W[i,i] = 1.0 / R[i]**2
    return W
    
  @staticmethod
  def calc_lambda(J):
    JT = numpy.transpose(J)
    JTJ = numpy.matmul(JT, J)
    JTJ_inv = numpy.linalg.inv(JTJ)
    trace = numpy.trace(JTJ_inv)
    l = trace**(-1)
    return l

  @staticmethod
  def calc_rss(f, p, data):
    return sum((f(p, data[:, 0]) - data[:, 1])**2)

  @staticmethod
  def calc_rss_from_r(R):
    return sum(R**2)
  
###################################
#  Load Data
###################################

  @staticmethod
  def load(file_in=None):
    try:
      if(not isinstance(file_in, (list,)) and not isinstance(file_in,numpy.ndarray) and file_in == None):
        if(len(sys.argv) > 1 and sys.argv[1] is not None):
          file_in = sys.argv[1]
          if(not os.path.isfile(file_in)):
            return None
        else:
          return None 
    except:
      pass

# List
    if(isinstance(file_in, (list,))):
      if(len(file_in) > 2 or (len(file_in) == 2 and len(file_in[0]) == 2)):
        data = numpy.zeros((len(file_in), 2), dtype=numpy.float64) 
        for row in range(len(file_in)):
          data[row,0] = float(file_in[row][0])
          data[row,1] = float(file_in[row][1])
        return data
      if(len(file_in) == 2 and len(file_in[0]) > 2):
        data = numpy.zeros((len(file_in[0]), 2), dtype=numpy.float64) 
        for row in range(len(file_in[0])):
          data[row,0] = float(file_in[0][row])
          data[row,1] = float(file_in[1][row])
        return data
# If array
    elif(isinstance(file_in,numpy.ndarray)): 
      data = numpy.zeros((len(file_in), 2), dtype=numpy.float64) 
      for row in range(len(file_in)):
        data[row,0] = float(file_in[row,0])
        data[row,1] = float(file_in[row,1])
      return data    
# If file
    elif(os.path.isfile(file_in)):
# Init variable
      file_data = ""

# Read it in line by line
      fh = open(file_in, "r")
      for file_row in fh:
        if(file_row.strip() != ""):
          file_data = file_data + file_row.strip() + '\n'
# Clean
      file_data = ng.clean(file_data)
      
# Read in data
      lines = file_data.split('\n')
  
# Load to list first
      data_list = []
      for line in lines:
        fields = line.split(',')
        if(len(fields) == 2):
          try: 
            data_list.append([float(fields[0]), float(fields[1])])
          except:
            pass
            
# Define and load into array
      data = numpy.zeros((len(data_list), 2), dtype=numpy.float64) 
      for row in range(len(data_list)):
        data[row,0] = float(data_list[row][0])
        data[row,1] = float(data_list[row][1])
      return data

  @staticmethod
  def clean(str_in):  
    str_out = ""
    l = len(str_in)
    for i in range(l):
# Last, Next, This
      if(i == 0):
        last = None
      else:
        last = str_in[i-1]
      if(i < (l-1)):
        next = str_in[i+1]
      else:  
        next = None
      char = str_in[i]
      
# Tabs and pipes to commas
      if(char == "\t"):
        char = ","
      if(char == "|"):
        char = ","
    
# Check
      ok = True
      if(last == " " and char == " "):
        ok = False
      elif(last == "\n" and char == "\n"):
        ok = False
      elif(last == "\n" and char == " "):
        ok = False
      elif(char == " " and next == "\n"):
        ok = False
      elif(char == "," and next == ","):
        ok = False

# Add to string
      if(ok):
        str_out += char
    return str_out      
  
###########################################
#  CLASS fittin
###########################################
class fitting:

#################################
#   Exponential Fit
#################################

  @staticmethod
  def single_exp(p, x):
    return p[0] + p[1] * numpy.exp(p[2] * x + p[3])  

  @staticmethod
  def fit_exp(data = None, points = 101):
    result = {
    'p': numpy.zeros((1)),
    'rss': 1e99,
    'data': numpy.zeros((1))
    }
  
# Load data
    data = fitting.load(data)
    y_min = numpy.min(data[:,1])
    y_max = numpy.max(data[:,1])
    dy = y_max - y_min
    y_mid = (y_min + y_max) / 2

    p = [y_mid,0,0,0]
    pl = [y_mid - abs(dy),-2 * abs(dy),-5,-1]
    pu = [y_mid + abs(dy),2 * abs(dy),5,1]
    pd = [abs(dy), abs(dy), 5, 1]
    p, rss = hybrid_sa.run(fitting.single_exp, p, {'p_diff': [abs(dy), abs(dy), 5, 1]})

# LMA
    p, rss = lma.run(fitting.single_exp, p, {}, data)
    result['p'] = p    
    result['rss'] = rss    
    result['data'] = data
    result['data_fit'] = fitting.data_points(fitting.single_exp, p, points, data)
    return result
    
###################################
#  Load Data
###################################

#   Loads data into numpy array
#   Tries to read from file if file_in is set
#   If file_in is a list or numpy array it reads it in directly
#   If file_in is none, it checks if a file has been specified as an argument

  @staticmethod
  def data_points(f, p, points, data):
    data_out = numpy.zeros((points,2))
    data_out[:,0] = numpy.linspace(data[0,0],data[-1,0],points)
    data_out[:,1] = f(p, data_out[:,0])
    return data_out
    
  @staticmethod
  def load(file_in=None):
    if(file_in == None):
      if(len(sys.argv) > 1 and sys.argv[1] is not None):
        file_in = sys.argv[1]
        if(not os.path.isfile(file_in)):
          return None
      else:
        return None 

# List
    if(isinstance(file_in, (list,))):
      if(len(file_in) > 2 or (len(file_in) == 2 and len(file_in[0]) == 2)):
        data = numpy.zeros((len(file_in), 2)) 
        for row in range(len(file_in)):
          data[row,0] = float(file_in[row][0])
          data[row,1] = float(file_in[row][1])
        return data
      if(len(file_in) == 2 and len(file_in[0]) > 2):
        data = numpy.zeros((len(file_in[0]), 2)) 
        for row in range(len(file_in[0])):
          data[row,0] = float(file_in[0][row])
          data[row,1] = float(file_in[1][row])
        return data
# If array
    elif(isinstance(file_in,numpy.ndarray)): 
      data = numpy.zeros((len(file_in), 2)) 
      for row in range(len(file_in)):
        data[row,0] = float(file_in[row,0])
        data[row,1] = float(file_in[row,1])
      return data    
# If file
    elif(os.path.isfile(file_in)):
# Init variable
      file_data = ""

# Read it in line by line
      fh = open(file_in, "r")
      for file_row in fh:
        if(file_row.strip() != ""):
          file_data = file_data + file_row.strip() + '\n'
# Clean
      file_data = fitting.clean(file_data)
      
# Read in data
      lines = file_data.split('\n')
  
# Load to list first
      data_list = []
      for line in lines:
        fields = line.split(',')
        if(len(fields) == 2):
          try: 
            data_list.append([float(fields[0]), float(fields[1])])
          except:
            pass
            
# Define and load into array
      data = numpy.zeros((len(data_list), 2)) 
      for row in range(len(data_list)):
        data[row,0] = float(data_list[row][0])
        data[row,1] = float(data_list[row][1])
      return data      

  @staticmethod
  def clean(str_in):  
    str_out = ""
    l = len(str_in)
    for i in range(l):
# Last, Next, This
      if(i == 0):
        last = None
      else:
        last = str_in[i-1]
      if(i < (l-1)):
        next = str_in[i+1]
      else:  
        next = None
      char = str_in[i]
      
# Tabs and pipes to commas
      if(char == "\t"):
        char = ","
      if(char == "|"):
        char = ","
    
# Check
      ok = True
      if(last == " " and char == " "):
        ok = False
      elif(last == "\n" and char == "\n"):
        ok = False
      elif(last == "\n" and char == " "):
        ok = False
      elif(char == " " and next == "\n"):
        ok = False
      elif(char == "," and next == ","):
        ok = False

# Add to string
      if(ok):
        str_out += char
    return str_out    

###########################################
#  CLASS hybrid_s
###########################################
class hybrid_sa:

  @staticmethod
  def run(f, p_in, settings_in={}, data = None):
# Load Data
    data = hybrid_sa.load(data)
  
# Default Settings
    settings = {
    'p_diff': None,
    'p_percent': 50,
    "temp_max": 10, 
    "temp_min": 0.1,
    "temp_min_loops": 10,
    "factor": 0.90,
    "iterations": 2000,
    "p_factor": 0.90,
    'pool_size': 50,
    'pool_seed_factor': 3,
    'breed_events': 3
    }
    
# Load settings
    settings = hybrid_sa.load_settings(settings_in, settings)    
    
    p = numpy.asarray(p_in)
    dp = hybrid_sa.set_bounds(p, settings)

    p, rss = hybrid_sa.solve(f, p, dp, settings, data) 
     
    return p, rss
#return pool['p'][0,:], pool['rss'][0]
    
  @staticmethod
  def load_settings(settings_in, settings): 
# Load new settings
    for key in settings_in:
      if(key in settings):
        settings[key] = settings_in[key]
    return settings
        
  @staticmethod
  def set_bounds(p, settings):  
    if(settings['p_diff'] != None):
      dp = numpy.asarray(settings['p_diff'])
    else:
      pl = numpy.copy(p) * (1 - settings['p_percent']/100)
      pu = numpy.copy(p) * (1 + settings['p_percent']/100)
      dp = pu - pl
    return dp
    
#########################################
#
#   Main loop with decreasing temperature
#
#########################################
  
  @staticmethod
  def solve(f, p, dp, settings, data): 
    pool_size = settings['pool_size']
    pool = hybrid_sa.make_pool(pool_size, len(p))
    hybrid_sa.init_pool(f, p, dp, settings, data, pool)
    
    t = settings["temp_max"]
    rss_best = None
    p_best, rss_best = hybrid_sa.best_in_pool(pool)
    
    t_loop = True
    n = 0
    while(t_loop):
      n = n + 1    
# Find solution at temperature
      rss_best, p_best = hybrid_sa.solve_loop_at_temp(f, pool, dp, t, settings, data, p_best, rss_best)
      
# Breed  solutions
      pool = hybrid_sa.next_generation(f, dp, settings, data, pool)

      t = t * settings["factor"]
      dp = dp * settings["p_factor"]    
      
      if(t < settings["temp_min"]):
        t_loop = False
      if(settings["temp_min_loops"] != None and n < settings["temp_min_loops"]):
        t_loop = True

    return p_best, rss_best
  
#########################################
#
#   SA Loop
#
#########################################
  
  @staticmethod
  def solve_loop_at_temp(f, pool, dp, t, settings, data, p_best, rss_best):   
# Loop through temperature reductions
    n = 0
    while(n<settings["iterations"]):
      pn = n%settings['pool_size']
      rss = pool['rss'][pn]
      p_trial = hybrid_sa.next_p(pool['p'][pn,:], dp)
      rss_trial = hybrid_sa.calc_rss(f, p_trial, data)
      if(rss_trial < rss):
        rss = rss_trial
        pool['p'][pn,:] = numpy.copy(p_trial)
        if(rss_best == None or rss < rss_best): 
          rss_best = rss
          p_best = numpy.copy(p_trial)
      elif(numpy.exp((rss_best - rss_trial) / t) > random.random()):  
        rss = rss_trial
        pool['p'][pn,:] = numpy.copy(p_trial)      
      n = n + 1
    return rss_best, p_best

#########################################
#
#   GA next gen and breeding
#
#########################################
  
  @staticmethod
  def next_generation(f, dp, settings, data, pool):
    new_pool = hybrid_sa.copy_pool(pool) 
    breed_events = int(settings['breed_events'] * len(pool['p']))
    for n in range(breed_events):
      p = hybrid_sa.breed(pool, dp)
      rss = hybrid_sa.calc_rss(f, p, data)
      hybrid_sa.add_to_pool(new_pool, p, rss)  
    pool = hybrid_sa.copy_pool(new_pool)
    return pool  
  
  @staticmethod
  def breed(pool, dp, mutate=0.05): 
    pool_len = len(pool['p'])
    p_len = len(pool['p'][0])    
    
    if(numpy.random.uniform(0.0,1.0) <= mutate):
      p_x = numpy.random.randint(0, pool_len)
      return hybrid_sa.mutant(pool['p'][p_x], dp)
    else:
      p_x = numpy.random.randint(0, pool_len)
      p_y = p_x
      while(p_y == p_x):
        p_y = numpy.random.randint(0, pool_len)      
      p = numpy.copy(pool['p'][p_x])
      for i in range(0, p_len):
        if(numpy.random.randint(0, 99) < 50):
          p[i] =  pool['p'][p_y, i]
      if(hybrid_sa.in_pool(pool, p)):
# No clones
        return hybrid_sa.mutant(pool['p'][p_x], dp)
      return p  

#########################################
#
#   Pool
#
#########################################
      
  @staticmethod
  def make_pool(pool_size, param_size):
    pool = {}
    pool['p'] = numpy.zeros((pool_size, param_size))
    pool['rss'] = numpy.zeros((pool_size))
    pool['rss'][:] = -1.0    
    return pool
    
  @staticmethod
  def copy_pool(pool_in):
    pool_out = {}
    pool_out['p'] = numpy.copy(pool_in['p']) 
    pool_out['rss'] = numpy.copy(pool_in['rss']) 
    return pool_out
    
  @staticmethod  
  def init_pool(f, p, dp, settings, data, pool):
    for n in range(settings['pool_size'] * settings['pool_seed_factor']):
      p_trial = hybrid_sa.mutant(p, dp)
      rss = hybrid_sa.calc_rss(f, p_trial, data)
      hybrid_sa.add_to_pool(pool, p_trial, rss)
  
  @staticmethod
  def add_to_pool(pool, p, rss):
    if(hybrid_sa.in_pool(pool, p)):
      return 0
    pool_len = len(pool['rss'])
# If worse than last entry, just leave
    if(rss > pool['rss'][-1] and pool['rss'][-1] != -1.0):
      return 0
    
    if(pool['rss'][0] == -1.0):
      pool['rss'][0] = rss
      pool['p'][0] = p
      return 0     
    
# Find where to insert
    for n in range(pool_len):
      if(rss < pool['rss'][n] or pool['rss'][n] == -1.0):
        pool['rss'] = numpy.insert(pool['rss'], n, rss)[0:pool_len]
        pool['p'] = numpy.insert(pool['p'], n, p, 0)[0:pool_len,:]
        break 

  @staticmethod
  def in_pool(pool, p):
    for i in range(len(pool['p'])):
      for j in range(len(p)):
        same = True
        if(p[j] != pool['p'][i, j]):
          same = False
          break
      if(same == True):
        return True
    return False       
        
  @staticmethod
  def best_in_pool(pool):      
    return numpy.copy(pool['p'][0,:]), pool['rss'][0]
        
  @staticmethod
  def mutant(p, dp):
    r = numpy.random.rand(len(p))
    return p[:] + (r[:] - 0.5) * dp[:]
  
  @staticmethod
  def next_p(p, dp):
    p_out = numpy.copy(p)
    r = numpy.random.rand()
    pn = numpy.random.randint(0, len(p_out))    
    p_out[pn] = p_out[pn] + (r - 0.5) * dp[pn]
    return p_out
    
  @staticmethod
  def calc_rss(f, p, data):
    try: 
      rss = sum((data[:,1] - f(p, data[:,0]))**2) 
      return rss
    except:
      return 1.0e99

###################################
#  Load Data
###################################

  @staticmethod
  def load(file_in=None):
    try:
      if(file_in == None):
        if(len(sys.argv) > 1 and sys.argv[1] is not None):
          file_in = sys.argv[1]
          if(not os.path.isfile(file_in)):
            return None
        else:
          return None 
    except:
      pass

# List
    if(isinstance(file_in, (list,))):
      if(len(file_in) > 2 or (len(file_in) == 2 and len(file_in[0]) == 2)):
        data = numpy.zeros((len(file_in), 2)) 
        for row in range(len(file_in)):
          data[row,0] = float(file_in[row][0])
          data[row,1] = float(file_in[row][1])
        return data
      if(len(file_in) == 2 and len(file_in[0]) > 2):
        data = numpy.zeros((len(file_in[0]), 2)) 
        for row in range(len(file_in[0])):
          data[row,0] = float(file_in[0][row])
          data[row,1] = float(file_in[1][row])
        return data
# If array
    elif(isinstance(file_in,numpy.ndarray)): 
      data = numpy.zeros((len(file_in), 2)) 
      for row in range(len(file_in)):
        data[row,0] = float(file_in[row,0])
        data[row,1] = float(file_in[row,1])
      return data    
# If file
    elif(os.path.isfile(file_in)):
# Init variable
      file_data = ""

# Read it in line by line
      fh = open(file_in, "r")
      for file_row in fh:
        if(file_row.strip() != ""):
          file_data = file_data + file_row.strip() + '\n'
# Clean
      file_data = hybrid_sa.clean(file_data)
      
# Read in data
      lines = file_data.split('\n')
  
# Load to list first
      data_list = []
      for line in lines:
        fields = line.split(',')
        if(len(fields) == 2):
          try: 
            data_list.append([float(fields[0]), float(fields[1])])
          except:
            pass
            
# Define and load into array
      data = numpy.zeros((len(data_list), 2)) 
      for row in range(len(data_list)):
        data[row,0] = float(data_list[row][0])
        data[row,1] = float(data_list[row][1])
      return data

  @staticmethod
  def clean(str_in):  
    str_out = ""
    l = len(str_in)
    for i in range(l):
# Last, Next, This
      if(i == 0):
        last = None
      else:
        last = str_in[i-1]
      if(i < (l-1)):
        next = str_in[i+1]
      else:  
        next = None
      char = str_in[i]
      
# Tabs and pipes to commas
      if(char == "\t"):
        char = ","
      if(char == "|"):
        char = ","
    
# Check
      ok = True
      if(last == " " and char == " "):
        ok = False
      elif(last == "\n" and char == "\n"):
        ok = False
      elif(last == "\n" and char == " "):
        ok = False
      elif(char == " " and next == "\n"):
        ok = False
      elif(char == "," and next == ","):
        ok = False

# Add to string
      if(ok):
        str_out += char
    return str_out      
  
###########################################
#  CLASS read_confi
###########################################
class read_config:
  
  @staticmethod
  def read_file(file_path):
    input = {}
    fh = open(file_path, 'r')
    for line in fh:
      input = read_config.read_line(input, line)
    fh.close()
    return input
    
  @staticmethod
  def read_line(input, line):
    line = line.strip()
    if(len(line) > 0 and line[0] != "#"):
      cmd, data = read_config.get_fields_dict(line)
      if(cmd not in input.keys()):
        input[cmd] = []
      input[cmd].append(data)
    return input
        
  @staticmethod  
  def get_fields_dict(line): 
    fields = read_config.split_by(line, ' ')    
    dict = {}    
    dict['COMMAND'] = fields[0]
    
    for field in fields:
      f = field.split("=")
      if(len(f) == 2):
        fb = read_config.split_by(f[1], ',')
        if(len(fb) == 1):
          dict[f[0].lower()] = [f[1]]
        elif(len(fb) > 1):
          dict[f[0].lower()] = fb
    
    return fields[0], dict
    
  @staticmethod        
  def get(data, c1, c2, n1=0, n2=0):
    if(c1 not in data.keys()):
      return None
    c1_list = data[c1]
    if(n1>=len(c1_list)):
      return ''
    if(c2 not in data[c1][n1].keys()):
      return None
    c2_list = data[c1][n1][c2]
    if(n2>=len(c2_list)):
      return ''
    return data[c1][n1][c2][n2]

  @staticmethod        
  def get_list(data, c1, c2, n1=0):
    if(c1 not in data.keys()):
      return None
    c1_list = data[c1]
    if(n1>=len(c1_list)):
      return ''
    if(c2 not in data[c1][n1].keys()):
      return None
    return data[c1][n1][c2]
        
  @staticmethod  
  def split_by(line, sep=' ', ignore_double_sep=True):
    last_char = None
    in_quotes = 0
    fields = []
    temp_line = ""
    
    for char in line:
      if(char == "'" and in_quotes == 0 and last_char != "\\"):
        in_quotes = 1
      elif(char == "'" and in_quotes == 1 and last_char != "\\"):
        in_quotes = 0
      elif(char == '"' and in_quotes == 0 and last_char != "\\"):
        in_quotes = 2
      elif(char == '"' and in_quotes == 2 and last_char != "\\"):
        in_quotes = 0
      elif(in_quotes > 0):
        temp_line = temp_line + char
      elif(in_quotes == 0 and char != sep):
        temp_line = temp_line + char
      elif(char == sep and last_char == sep and ignore_double_sep):
        pass
      elif(char == sep):
        fields.append(temp_line)
        temp_line = "" 
    if(temp_line != ""):
      fields.append(temp_line)
    
    return fields

###########################################
#  CLASS mis
###########################################
class misc:

  @staticmethod
  def makedir(dir):
    try:
      os.mkdir(dir)
    except:
      pass

###########################################
#  CLASS unit
###########################################
class units:
  
  @staticmethod
  def convert(conv_from, conv_to, value_in):
    try:
      value_in = float(value_in)
    except:
      return None
    conv_from = conv_from.upper()
    conv_to = conv_to.upper()

# LENGTH METERS
    length = {
    'M': 1.0,
    'CM': 100,
    'MM': 1E3,
    'UM': 1E6,
    'NM': 1E9,
    'ANG': 1E10,
    'BOHR': 1.89E10,
    }

# ENERGY J
    energy = {
    'J': 1.0,
    'EV': 6.2415E18,
    'RY': 4.5874E17,
    }

# FORCE N
    force = {
    'N': 1.0,
    'RY/BOHR': 2.4276e7,
    'EV/ANG':6.2414E8,    
    }
    
# VELOCITY
    velocity = {
    'M/S': 1.0,
    'MPH': 2.25,    
    }
    
# PRESSURE
    pressure = {
    'PA': 1.0,
    'GPA': 1.0E-9,    
    'BAR': 1.0E-5,    
    'ATMOSPHERE': 9.8692E-6,    
    'PSI': 1.45038E-4, 
    'KBAR': 1.0E-8,   
    'RY/BOHR3': 6.857E-14,   
    'EV/ANG3': 6.241E-12
    }
    
    unit_list = [length, energy, force, velocity, pressure]
    
    for l in unit_list:
      if(conv_from in l.keys() and conv_to in l.keys()):
        return round((l[conv_to] / l[conv_from]) * float(value_in),9)
  
"""  
  @staticmethod
  def convert(conv_from, conv_to, value_in):

    conv_from = conv_from.upper()
    conv_to = conv_to.upper()

# METERS
    length = {
    'M': 1.0,
    'CM': 100,
    'MM': 1E3,
    'UM': 1E6,
    'NM': 1E9,
    'ANG': 1E10,
    'BOHR': 1.89E10,
    }

# J
    energy = {
    'J': 1.0,
    'EV': 6.242E18,
    'RY': 4.5874E17,
    }

    if(conv_from in length.keys() and conv_to in length.keys()):
      return round((length[conv_to] / length[conv_from]) * float(value_in),9)

    if(conv_from in energy.keys() and conv_to in energy.keys()):
      return round((energy[conv_to] / energy[conv_from]) * float(value_in),9)
"""

###########################################
###########################################
#  MAIN
###########################################
###########################################

new_run = pwscf_converge()

new_run.run()

