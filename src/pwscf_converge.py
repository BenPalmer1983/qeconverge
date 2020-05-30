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
from std import std
from matplotlib import cm
from mpl_toolkits import mplot3d
from pwscf_input import pwscf_input
from pwscf_output import pwscf_output
from pwscf_exec import pwscf_exec
from read_config import read_config
from units import units
from globals import globals
from read_input import read_input
from conv_ecut import conv_ecut
from conv_ecut2d import conv_ecut2d
from conv_kpoints import conv_kpoints
from make_plots import make_plots
from post_plot import post_plot

######################
# 
######################

class pwscf_converge:

  def __init__(self):
    # Make dir and template
    self.make_dirs()
    
    

  def run(self):
    print("QE Converge")
    print(len(sys.argv))
    print(sys.argv)
  
    if(len(sys.argv) == 1):
      pwscf_converge.run_calc()
    elif(len(sys.argv)>2):
      if(sys.argv[2] == "calc"):
        pwscf_converge.run_calc()
      elif(sys.argv[2] == "plot"):
        pwscf_converge.run_plot()
  
  
  def run_calc():
    print("Running Calc")
    
    # Open log file
    globals.log_fh = open(globals.dirs['logs'] + '/qeeos.log', 'w')
    globals.log_fh.write('QE CONVERGE \n')
  
    # Start time
    globals.times['start'] = time.time()

    # Read Input
    read_input.run()
    
    # LOAD TEMPLATE FILE
    template = pwscf_input()
    template.load(globals.d['pwscf_template'])
    template.set_dirs()
    template.save("input_template.in", globals.dirs['td'])  
    
    # ECUT CONV
    conv_ecut.run()
    
    # ECUT 2D CONV
    conv_ecut2d.run()
    
    # KPOINTS CONV
    conv_kpoints.run()
    
    
    globals.log_fh.write('\nCalculations complete.\n\n')
    globals.log_fh.write('CPU min:     ' + str(globals.pwscf_times['cpu_min']) + '\n')
    globals.log_fh.write('CPU max:     ' + str(globals.pwscf_times['cpu_max']) + '\n')
    globals.log_fh.write('CPU total:   ' + str(globals.pwscf_times['cpu_total']) + '\n')
    globals.log_fh.write('\n')
    globals.log_fh.write('Wall min:    ' + str(globals.pwscf_times['wall_min']) + '\n')
    globals.log_fh.write('Wall max:    ' + str(globals.pwscf_times['wall_max']) + '\n')
    globals.log_fh.write('Wall total:  ' + str(globals.pwscf_times['wall_total']) + '\n')
    globals.log_fh.write('\n')
    
    
    
    # End time
    globals.times['end'] = time.time()
    globals.times['duration'] = globals.times['end'] - globals.times['start']
    globals.log_fh.write('Time:  ' + str(globals.times['duration']) + '\n')
    globals.log_fh.close()
    
    
    
    # Plots
    make_plots.run()
    
    
    # Exit Program
    sys.exit()
    

     
  def run_plot():
    print("Running Plot")
    post_plot.run()
  
  
  
  
  
  
  
  
  
  
  
  

  def make_dirs(self):  
    for d in globals.dirs.keys():
      dir = globals.dirs[d]
      std.make_dir(dir)
  

