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

######################
# 
######################

class pwscf_converge:

  def __init__(self):
    # Make dir and template
    self.make_dirs()
    
    

  def run(self):
    
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
    
    
    # End time
    globals.times['end'] = time.time()
    globals.times['duration'] = globals.times['end'] - globals.times['start']
    globals.log_fh.write('Time:  ' + str(globals.times['duration']) + '\n')
    globals.log_fh.close()
    
    # Exit Program
    sys.exit()
    

     

  def make_dirs(self):
  
    for d in globals.dirs.keys():
      dir = globals.dirs[d]
      std.make_dir(dir)
  

