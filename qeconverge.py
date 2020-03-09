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
  
###########################################
#  CLASS st
###########################################
class std:

  @staticmethod
  def file_to_list(file_name, clean=False):
# Init variable
    file_data = []
# Read it in line by line
    fh = open(file_name, "r")
    for line in fh:
      if(clean):
        line = line.strip()
        if(line != ""):
          file_data.append(line)          
      else:
        file_data.append(line[0:-1])
# Return
    return file_data
    
  @staticmethod
  def split_fields(line, sep=" "):
    out = line.split(sep)
    key = out[0]
    value = out[1]
    value_out = ''    
    indata = False
    for char in value:
      if(indata and char != '"'):
        value_out = value_out + char
      elif(indata and char == '"'):
        indata = False
      elif(not indata and char == '"'):
        indata = True
    return key, value_out
    
  @staticmethod
  def one_space(line, sep=" "):
    out = ''   
    indata = 0
    last_char = None
    for char in line:
      if(indata == 1 and char != "'" and last_char != "\\"):
        out = out + char
      elif(indata == 1 and char == "'" and last_char != "\\"):
        out = out + char
        indata = 0
      elif(indata == 2 and char != '"' and last_char != "\\"):
        out = out + char
      elif(indata == 2 and char == '"' and last_char != "\\"):
        out = out + char
        indata = 0
      elif(indata == 0 and not (char == " " and last_char == " ")):
        out = out + char
    return out   
    
  @staticmethod
  def to_fields(line, sep=" "):
    out = []
    temp = ''
    indata = 0
    last_char = None
    for char in line:
      if(indata == 1 and char != "'" and last_char != "\\"):
        temp = temp + char
      elif(indata == 1 and char == "'" and last_char != "\\"):
        temp = temp + char
        indata = 0
      elif(indata == 2 and char != '"' and last_char != "\\"):
        temp = temp + char
      elif(indata == 2 and char == '"' and last_char != "\\"):
        temp = temp + char
        indata = 0
      elif(indata == 0 and not (char == sep and last_char == sep)):
        if(char == sep):
          temp = temp.strip()
          if(temp != ""):
            out.append(temp)
            temp = ''
        else:
          temp = temp + char
    
    temp = temp.strip()
    if(temp != ""):
      out.append(temp)      
    return out    
    
  @staticmethod
  def make_dir(dir):
    try:
      if(not os.path.exists(dir) and dir.strip() != ''):
        os.mkdir(dir) 
        return True
      return False
    except:
      return False
    
  @staticmethod
  def remove_comments(content):
    data = ''
    i = 0
    for line in content:
      if(i > 0):
        data += '\n'
      data += line
      i = i + 1
    out = ''
    indata = 0
    incomment = 0
    for i in range(len(data)):
# Get char and next char
      char = data[i]
      next = None
      prev = None
      if(i < len(data)-1):
        next = data[i + 1]
      if(i > 0):
        prev = data[i - 1]
# If in '  '
      if(indata == 1 and char != "'" and last_char != "\\"):
        out = out + char
      elif(indata == 1 and char == "'" and last_char != "\\"):
        out = out + char
        indata = 0
# If in "  "
      elif(indata == 2 and char != '"' and last_char != "\\"):
        out = out + char
      elif(indata == 2 and char == '"' and last_char != "\\"):
        out = out + char
        indata = 0
      elif(indata == 0):
        if(incomment == 0 and char == "/" and next == "/"):
          incomment = 1
        elif(incomment == 1 and char == "\n"):
          incomment = 0
        if(incomment == 0 and char == "!"):
          incomment = 2
        elif(incomment == 2 and char == "\n"):
          incomment = 0
        if(incomment == 0 and char == "/" and next == "*"):
          incomment = 3
        elif(incomment == 3 and prev == "*" and char == "/"):
          incomment = 0
        elif(incomment == 0):
          out = out + char  
    return out.split("\n")    
    
# Remove comments from a block of data/text
  @staticmethod
  def remove_comments_data(data):
    out = ""
    n = 0
    inquotes = 0
    incomment = 0
    while n < len(data):
# Get char and next char
      char = data[n]
      next = None
      prev = None
      if(n < len(data)-1):
        next = data[n + 1]
      if(n > 0):
        prev = data[n - 1]
        
# If in '  '
      if(inquotes == 1 and char != "'" and last_char != "\\"):
        out = out + char
      elif(inquotes == 1 and char == "'" and last_char != "\\"):
        out = out + char
        inquotes = 0
# If in "  "
      elif(inquotes == 2 and char != '"' and last_char != "\\"):
        out = out + char
      elif(inquotes == 2 and char == '"' and last_char != "\\"):
        out = out + char
        inquotes = 0
# If not inside quotes
      elif(inquotes == 0):
# Comment on a line
        if(incomment == 0 and char == "/" and next == "/"):
          incomment = 1
        elif(incomment == 0 and char == "!"):
          incomment = 1
        elif(incomment == 0 and char == "#"):
          incomment = 1    
# Comment on line close
        elif(incomment == 1 and char == "\n"):
          incomment = 0
# Comment block
        elif(incomment == 0 and char == "/" and next == "*"):
          incomment = 3
        elif(incomment == 3 and prev == "*" and char == "/"):
          incomment = 0
        elif(incomment == 0):
          out = out + char  
# Increment counter
      n = n + 1
    return out        

# Single spaces, tabs to spaces
  @staticmethod
  def prep_data(content):
    out = []
    for line in content:
      line_new = std.prep_data_line(line)
      if(line_new != ''):
        out.append(line_new)
    return out  
      
  @staticmethod
  def prep_data_line(line): 
    temp = ''
    indata = 0
    last_char = None
    for char in line:
      if(char == '\t'):
        char = ' '
      if(indata == 1 and char != "'" and last_char != "\\"):
        temp = temp + char
      elif(indata == 1 and char == "'" and last_char != "\\"):
        temp = temp + char
        indata = 0
      elif(indata == 2 and char != '"' and last_char != "\\"):
        temp = temp + char
      elif(indata == 2 and char == '"' and last_char != "\\"):
        temp = temp + char
        indata = 0
      elif(indata == 0 and not (char == ' ' and last_char == ' ')):
        temp = temp + char       
      last_char = char  
    return temp.strip()    
    
  @staticmethod
  def remove_quotes(inp): 
    if(isinstance(inp, list)):    
      for i in range(len(inp)):
        inp[i] = std.remove_quotes(inp[i])        
      return inp
    else:
      inp = inp.strip()
      if(inp[0] == '"' and inp[-1] == '"'):
        return inp[1:-1]
      if(inp[0] == "'" and inp[-1] == "'"):
        return inp[1:-1]
      return inp
      
  @staticmethod
  def config_file_to_list(file_name):
# Init variable
    file_data = []
# Read it in line by line
    fh = open(file_name, "r")
    for line in fh:
      if(line.strip() != ""):
        line = line.strip()
        line = std.remove_comments(line)
        line = std.prep_data_line(line)
        fields = std.to_fields(line)
        file_data.append(fields)         
# Return
    file_data = std.remove_quotes(file_data)
    return file_data
    
  @staticmethod
  def get_dir(file_path):
    directory = ''
    read = False
    for i in range(len(file_path)):
      if(read):
        directory = file_path[-1-i] + directory
      if(file_path[-1-i] == "/"):
        read = True
    return directory
  
  @staticmethod
  def read_csv(filename):
    data = []
    if(os.path.isfile(filename)):
# Read from file into memory
      fh = open(filename, 'r')
      file_data = ""
      for line in fh:
        file_data = file_data + line
      fh.close()
# Remove comments
      file_data = std.remove_comments_data(file_data)
# Read Data
      lines = file_data.split("\n")
      for line in lines:
        line = line.strip()
        if(line != ""):
          data.append(line.split(","))  
    return data
    
  @staticmethod
  def write_csv(filename, arr):  
    fh = open(filename, 'w')
    for i in range(len(arr)):
      for j in range(len(arr[i])):
        fh.write(std.float_padded(arr[i,j],8) + " ")
      fh.write("\n")
    fh.close()
  
  @staticmethod
  def option(input):
    input = input.strip().upper()
    if(input[0:1] == "Y"):
      return True
    elif(input[0:2] == "ON"):
      return True
    elif(input[0:1] == "T"):
      return True
    else:
      return False
    
  @staticmethod
  def float_padded(inp, pad=7):
    out = float(inp)
    out = round(out, pad-3)
    out = str(out)  
    while(len(out)<pad):
      out = out + " "      
    return out[0:pad]
    
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

  def make(self):
  
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
        for field in position:
          file += field + "   "
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
    
  def save(self, file=None, dir=None):
# Build latest version of file
    self.make()
    
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

  def set_dirs(self, scratch_dir = None, pp_dir = None):
    if(scratch_dir != None):
      self.scratch_dir = scratch_dir
    if(pp_dir != None):
      self.scratch_dir = pp_dir
   
    self.control['outdir'] = '"' + self.scratch_dir + '"'
    self.control['pseudo_dir'] = '"' + self.pp_dir + '"'

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
    
# SC
    if(type == "SC"):
      labels = []
      for row in self.atomic_species:
        labels.append(row[0])
      atoms, c_atoms, n_atoms = pwscf_standard.sc(labels,size) 
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
    
# FCC
    if(type == "FCC"):
      labels = []
      for row in self.atomic_species:
        labels.append(row[0])
      atoms, c_atoms, n_atoms = pwscf_standard.fcc(labels,size) 
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
      self.calcs()

  def reset(self):
  
    self.z = numpy.zeros((3,3))
    
# Important so store in it's own variable
    self.atom_count = 1   

# Control
    self.data = {
      "ok": False,
      "job_done": False,
      "error": False,
      "error_code": None,
      "converged": False,
      "converged_in": None,
      
      "type": None,
      "summary": None,
      "mpi_processes": None,
      "threads_per_mpi_process": None,
      
      "species": {},
      "scf_settings": None,      
      "crystals": [],
      "results": [],
      
      "initial_positions": None,   
      "total_energy": None,
      "density_full": None,
      "stress": numpy.zeros((3,3)),
      "stress_sum": None,      
      "cpu_time": None,
      "wall_time": None,   
      "xyz": [],
      
      "mass_per_crystal": 0.0,
      "density": [],
    }
    
# Defaults
    self.xyz_units = 'evang'
    self.stress_units = 'gpa'
    
# Constants
    self.avogadro = 6.02214086E23

  def scf_settings(self):
    return {
    "bravais_lattice_index": None,
    "alat": None,
    "volume": None,
    "electrons": None, 
    "electrons_up": None, 
    "electrons_down": None, 
    "ecut_wfc": None, 
    "ecut_rho": None, 
    "convergence_threshold": None, 
    "mixing_beta": None, 
    "atomic_species": {},
    }
    
  def scf_crystal(self):
    return {
    "alat": 0.0,
    "cell_parameters": numpy.zeros((3,3)),
    "position_units": None,
    "atomic_labels": [],
    "atomic_positions": numpy.zeros((self.atom_count,3)),    
    "alat_adj": 0.0,
    "cell_parameters_adj": numpy.zeros((3,3)),
    "crystal_positions": numpy.zeros((self.atom_count,3)),
    "cell_volume": 0.0,
    "cell_density": 0.0,
    }

  def scf_results(self):
    return {
    "energy": 0.0,
    "total_force": 0.0,
    "stress": numpy.zeros((3,3)),
    "forces": numpy.zeros((self.atom_count,3)),
    "f_on": False,
    "s_on": False,
    }

#  Load, and use in another program
  def load(self, file_name): 
  
# Load data from file
    data_file = self.load_from_file(file_name)
    self.d = data_file.split("\n") 
  
#print(self.d)
  
# Reset data store
    self.reset()
       
# Load
    self.load_status()
    self.load_type()
    self.load_count()
    self.load_cpuinfo()
    self.load_crystal()
    self.load_scf_settings()
    self.load_results()
    self.load_species()
    
#print(len(self.data['crystals']))
#print(len(self.data['results']))
#self.load_scf('final_scf')
#self.load_results('initial_scf')
#self.load_results('final_scf')
    
  def load_status(self):    
# OK
###################################
    self.data['ok'] = False
    counter = 0

    for line in self.d:
      line = line.strip()
      if(pwscf_output.compare(line, "JOB DONE.")):
        self.data['job_done'] = True
      if(pwscf_output.compare(line, "Exit code:")):
        self.data['error'] = True              
        try:
          self.data['error_code'] = int(line[13:-1])
        except:
          pass
      if(pwscf_output.compare(line, "convergence has been achieved")):
        self.data['converged'] = True
        try:
          self.data['converged_in'] = int(line[39:42])
        except:
          pass
    if(self.data['job_done']):
      self.data['ok'] = True
  
  def load_type(self):
# Calc Type
###################################
    self.data['type'] = "SCF"
    for line in self.d:
      line = line.strip()
      if(line[0:23] == "A final scf calculation"):
        self.data['type'] = "VC-RELAX"
    
  def load_count(self):
    n = 0
    while(n < len(self.d)):
      n, line, line_uc = self.next_line(n, self.d)
      if(pwscf_output.compare(line, "number of atoms/cell      =")):
        count = pwscf_output.extract(line, "=", "", "i")  
        try:
          self.atom_count = int(count)
          return self.atom_count 
        except:
          return 0
          
  def load_cpuinfo(self):
    counter = 0
    n = 0
    while(n < len(self.d)):
      n, line, line_uc = self.next_line(n, self.d)
      if(line != ""):
        counter += 1
        if(counter == 1):
          self.data['summary'] = line
        else:
          if(pwscf_output.compare(line, "Number of MPI processes:")):
            self.data['mpi_processes'] = pwscf_output.extract(line, ":", "", "i") 
          elif(pwscf_output.compare(line, "Threads/MPI process:")):
            self.data['threads_per_mpi_process'] = pwscf_output.extract(line, ":", "", "i") 
   
  def load_species(self):
    n = 0
    while(n < len(self.d)):
      n, line, line_uc = self.next_line(n, self.d)
      if(line[0:14] == "atomic species"):
        n, line, line_uc = self.next_line(n, self.d)
        while(line.strip() != ""):      
          line = pwscf_output.single_spaces(line).strip()
          f = line.split(" ")
          self.data['scf_settings'][f[0]] = [float(f[1]), float(f[2])]   # Valence electrons, atomic mass
          n, line, line_uc = self.next_line(n, self.d)
        break
      
#self.data['scf_settings']
#
#  atomic species   valence    mass     pseudopotential
#  Al             3.00    26.98200     Al( 1.00)
 
  def load_scf_settings(self):
    self.data['scf_settings'] = self.scf_settings()  
    
    n = 0
    while(n < len(self.d)):
      n, line, line_uc = self.next_line(n, self.d)
      if(pwscf_output.compare(line, "bravais-lattice index     =")):
        self.data['scf_settings']['bravais_lattice_index'] = pwscf_output.extract(line, "=", "", "i") 
      elif(pwscf_output.compare(line, "lattice parameter (alat)  =")):
        self.data['scf_settings']['alat'] = pwscf_output.extract(line, "=", "a.u.", "f")  
      elif(pwscf_output.compare(line, "unit-cell volume          =")):
        self.data['scf_settings']['volume'] = pwscf_output.extract(line, "=", "(a.u.)^3", "f")     
      elif(pwscf_output.compare(line, "number of atoms/cell      =")):
        self.data['scf_settings']['nat'] = pwscf_output.extract(line, "=", "", "i")  
      elif(pwscf_output.compare(line, "number of atomic types    =")):
        self.data['scf_settings']['types'] = pwscf_output.extract(line, "=", "", "i")  
      elif(pwscf_output.compare(line, "number of electrons       =")):
        str_e = pwscf_output.extract(line, "=", "", "s")
        e, eu, ed = pwscf_output.electron_string(str_e)
        self.data['scf_settings']['electrons'] = e
        self.data['scf_settings']['electrons_up'] = eu
        self.data['scf_settings']['electrons_down'] = ed
      elif(pwscf_output.compare(line, "kinetic-energy cutoff     =")):
        self.data['scf_settings']['ecut_wfc'] = pwscf_output.extract(line, "=", "Ry", "f")  
      elif(pwscf_output.compare(line, "charge density cutoff     =")):
        self.data['scf_settings']['ecut_rho'] = pwscf_output.extract(line, "=", "Ry", "f")  
      elif(pwscf_output.compare(line, "convergence threshold     =")):
        self.data['scf_settings']['convergence_threshold'] = pwscf_output.extract(line, "=", "", "f")  
      elif(pwscf_output.compare(line, "mixing beta               =")):
        self.data['scf_settings']['mixing_beta'] = pwscf_output.extract(line, "=", "", "f")  
      elif(("atomic species" in line) and ("valence" in line) and ("mass" in line) and ("pseudopotential" in line)):
        loop = True
        while(loop):
          n, line, line_uc = self.next_line(n, self.d)
          if(line.strip() == ""):
            loop = False
          else:  
            line = pwscf_output.single_spaces(line).strip()
            line_arr = line.split(" ")
            self.data['scf_settings']['atomic_species'][line_arr[0]] = {}
            self.data['scf_settings']['atomic_species'][line_arr[0]]['valence'] = line_arr[1]
            self.data['scf_settings']['atomic_species'][line_arr[0]]['mass'] = line_arr[2]

# End of file/loop
        n = len(self.d)

###################################
# LOAD CRYSTALS FROM OUTPUT FILE
###################################

  def load_crystal(self):
# Make new list for crystals
    self.data['crystals'] = []
    
# FIRST
    n = 0
    crystal = self.scf_crystal()
    while(n < len(self.d)):
      n, line, line_uc = self.next_line(n, self.d)  
      if(line[0:10] == "celldm(1)="):
        crystal['alat'] = float(line[10:21].strip())
      elif(pwscf_output.compare(line.strip(), "crystal axes:")):
        for j in range(3):              
          n, line, line_uc = self.next_line(n, self.d)
          fields = pwscf_output.extract(line, "= (", ")", "s", " ")
          for i in range(len(fields)):
            crystal['cell_parameters'][j, i] = float(fields[i])
      elif(pwscf_output.compare(line.strip(), "Cartesian axes")):
        n, line, line_uc = self.next_line(n, self.d)
        n, line, line_uc = self.next_line(n, self.d)
        
# Unit
        crystal['position_units'] = "alat"
        
        loop = True 
        k = 0
        while(loop):
          n, line, line_uc = self.next_line(n, self.d)   
          if(line.strip() == ""):
            loop = False
          else:
            line_arr = line.split("tau(")
            label = line_arr[0][-15:].strip()
            crystal['atomic_labels'].append(label)
            
            coords = line_arr[1]
            x = float(coords[9:21])
            y = float(coords[22:33])
            z = float(coords[34:44])
            crystal['atomic_positions'][k, 0] = x
            crystal['atomic_positions'][k, 1] = y
            crystal['atomic_positions'][k, 2] = z 
            
# Increment
            k = k + 1
        
# Cell Volume
        crystal['cell_volume'] = pwscf_output.cell_volume(crystal['alat'], crystal['cell_parameters'][:,:])    
        
# Add/Save
        self.data['crystals'].append(crystal)
        n = len(self.d)
    
# MIDDLE
    n = 0
    while(n < len(self.d)):
      n, line, line_uc = self.next_line(n, self.d)   
      if(line[0:15] == "CELL_PARAMETERS"):
# Create
        crystal = self.scf_crystal()
        
# Get alat
        line_arr = line.split("=")
        line_arr = line_arr[1].split(")")
        crystal['alat'] = float(line_arr[0].strip())
        
#Cell Parameters
        for j in range(3):              
          n, line, line_uc = self.next_line(n, self.d)
          line = pwscf_output.single_spaces(line)
          fields = line.split(" ")
          for i in range(len(fields)):
            crystal['cell_parameters'][j, i] = float(fields[i])
      elif(line[0:9] == "density ="):  
        crystal['density'] = float(line[10:23])      
            
      elif(line[0:16] == "ATOMIC_POSITIONS"):     
        
# Unit
        crystal['position_units'] = "crystal"
        
# Read Coords
        loop = True 
        k = 0
        while(loop):
          n, line, line_uc = self.next_line(n, self.d)   
          if(line.strip() == ""):
            loop = False
          elif(line.strip() == "End final coordinates"):
            loop = False
          else:
            line = pwscf_output.single_spaces(line)
            line_arr = line.split(" ")
            crystal['atomic_labels'].append(line_arr[0])
            
            crystal['atomic_positions'][k, 0] = float(line_arr[1])
            crystal['atomic_positions'][k, 1] = float(line_arr[2])
            crystal['atomic_positions'][k, 2] = float(line_arr[3]) 
            
# Increment
            k = k + 1
            
# Cell Volume
        crystal['cell_volume'] = pwscf_output.cell_volume(crystal['alat'], crystal['cell_parameters'][:,:])    
            
# Add/Save
        self.data['crystals'].append(crystal)
    
# END
    n = 0
    crystal = self.scf_crystal()
    d = 0
    while(n < len(self.d)):
      n, line, line_uc = self.next_line(n, self.d)   
      if(line[0:10] == "celldm(1)="):
        d = d + 1
        if(d == 2):
          crystal['alat'] = float(line[10:21].strip())
      elif(d == 2 and pwscf_output.compare(line.strip(), "crystal axes:")):
        for j in range(3):              
          n, line, line_uc = self.next_line(n, self.d)
          fields = pwscf_output.extract(line, "= (", ")", "s", " ")
          for i in range(len(fields)):
            crystal['cell_parameters'][j, i] = float(fields[i])
     
      elif(d == 2 and pwscf_output.compare(line.strip(), "Cartesian axes")):      
        
# Unit
        crystal['position_units'] = "alat"
        
# Read coords
        n, line, line_uc = self.next_line(n, self.d)
        n, line, line_uc = self.next_line(n, self.d)
        
        loop = True 
        k = 0
        while(loop):
          n, line, line_uc = self.next_line(n, self.d)   
          if(line.strip() == ""):
            loop = False
          else:
            line_arr = line.split("tau(")
            label = line_arr[0][-15:].strip()
            crystal['atomic_labels'].append(label)
            
            coords = line_arr[1]
            x = float(coords[9:21])
            y = float(coords[22:33])
            z = float(coords[34:44])
            crystal['atomic_positions'][k, 0] = x
            crystal['atomic_positions'][k, 1] = y
            crystal['atomic_positions'][k, 2] = z 
            
# Increment
            k = k + 1
            
# Cell Volume
        crystal['cell_volume'] = pwscf_output.cell_volume(crystal['alat'], crystal['cell_parameters'][:,:])    
            
# Add/Save
        self.data['crystals'].append(crystal)
        n = len(self.d)
    
# Loop through crystals
    for i in range(len(self.data['crystals'])):
# Adjust alat and cell_parameters so celldm(1)=1.0
      factor = 1.0 / self.data['crystals'][i]['cell_parameters'][0, 0]
      
      self.data['crystals'][i]['alat_adj'] = self.data['crystals'][i]['alat'] * self.data['crystals'][i]['cell_parameters'][0, 0]
      self.data['crystals'][i]['cell_parameters_adj'][:, :] = factor * self.data['crystals'][i]['cell_parameters'][:, :]
    
# Make crystal_positions
      if(self.data['crystals'][i]['position_units'] == 'crystal'):
        self.data['crystals'][i]['crystal_positions'][:,:] = self.data['crystals'][i]['atomic_positions'][:,:] 
      elif(self.data['crystals'][i]['position_units'] == 'alat'):  
        minv = numpy.linalg.inv(self.data['crystals'][i]['cell_parameters'][:, :])
        for j in range(len(self.data['crystals'][i]['atomic_positions'])):
          self.data['crystals'][i]['crystal_positions'][j, :] = numpy.matmul(minv[:,:], self.data['crystals'][i]['atomic_positions'][j, :])

  def load_results(self):
  
# Make new list for results
    self.data['results'] = []
    
    n = 0
    while(n < len(self.d)):
      n, line, line_uc = self.next_line(n, self.d)
      
# READ ENERGY
      if(pwscf_output.compare(line, "!    total energy")):
# Create dictionary
        results = self.scf_results()  
        results['energy'] = pwscf_output.extract(line, "=", "Ry", "f") 
        
# READ FORCES
      elif(pwscf_output.compare(line, "Forces acting on atoms")):
        n, line, line_uc = self.next_line(n, self.d)  
        loop = True
        f = 0
        while(loop):
          n, line, line_uc = self.next_line(n, self.d)  
          if(line.strip() == ""):
            loop = False
          else:
            line_arr = line.split("force =")
            fields = pwscf_output.single_spaces(line_arr[1].strip()).split(" ")
            results['forces'][f,0] = float(fields[0])
            results['forces'][f,1] = float(fields[1])
            results['forces'][f,2] = float(fields[2])
            f = f + 1
        if(f>0):
          results['f_on'] = True
        
# READ TOTAL FORCE
      elif(pwscf_output.compare(line, "Total force =")):
        results['total_force'] = pwscf_output.extract(line, "=", "T", "f")
        
# READ STRESS
      elif(pwscf_output.compare(line, "total   stress  (Ry/bohr**3)")):  
        results['s_on'] = True      
        for j in range(3):              
          n, line, line_uc = self.next_line(n, self.d)  
          fields = pwscf_output.extract(line, "", "", "f", " ", True)  
          results['stress'][j,0] = fields[3] 
          results['stress'][j,1] = fields[4] 
          results['stress'][j,2] = fields[5]

#SAVE
        self.data['results'].append(results)

  def calcs(self):
  
# MASS PER CELL
    self.data['mass_per_crystal'] = 0.0    
    for l in self.data['crystals'][0]['atomic_labels']:
      self.data['mass_per_crystal'] = self.data['mass_per_crystal'] + self.data['scf_settings'][l][1]    
      
# CELL VOLUMES
    for i in range(len(self.data['crystals'])):
      self.data['crystals'][i]['cell_volume'] = pwscf_output.cell_volume(self.data['crystals'][i]['alat'], self.data['crystals'][i]['cell_parameters'][:,:])  
      self.data['crystals'][i]['cell_density'] = pwscf_output.cell_density(self.data['crystals'][i]['cell_volume'], self.data['mass_per_crystal'])
      
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

  def get_nat(self): 
    return self.atom_count

  def get_alat(self):
    return self.data['alat']
    
  def get_volume(self):
    return self.data['scf_settings']['volume'] 
    
  def get_volume_per_atom(self):
    return self.data['scf_settings']['volume'] / self.atom_count
  
  def get_total_energy(self, n=None):
    if(n == None):
      return self.data['results'][-1]['energy']  
    else:
      return self.data['results'][n]['energy']  
    
  def get_energy_per_atom(self, n=None):
    if(n == None):
      return self.data['results'][-1]['energy'] / self.atom_count  
    else:
      return self.data['results'][n]['energy'] / self.atom_count 
    
  def get_total_force(self):
    if(n == None):
      return self.data['results'][-1]['total_force']  
    else:
      return self.data['results'][n]['total_force']
    
  def get_force_per_atom(self, n=None):
    if(n == None):
      return self.data['results'][-1]['total_force'] / self.atom_count  
    else:
      return self.data['results'][n]['total_force'] / self.atom_count  
  
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
    alat = self.data['crystals'][-2]['alat_adj']
    cp = self.data['crystals'][-2]['cell_parameters_adj']
    return alat, cp

# Get stress
  def get_stress(self):
    return self.data['stress']
    
  def get_stress_sum(self):
    return self.data['stress_sum']
  
  def get_job_done(self):
    return self.data['job_done']
    
  def get_job_error(self):
    return self.data['error']
    
  def get_job_converged(self):
    return self.data['converged']

  def get_ok(self):
    return self.data['ok']

# GET - VC-RELAXED

  def get_vc_relax_alat(self):
    if(self.data['type']=="VC-RELAX"):
      return self.data['crystals'][-1]['alat_adj']
    else:
      return None

  def get_vc_relax_cp(self):
    if(self.data['type']=="VC-RELAX"):
      return self.data['crystals'][-1]['cell_parameters_adj']
    else:
      return None

  def get_vc_relax_volume(self):
    if(self.data['type']=="VC-RELAX"):
      return self.data['crystals'][-1]['cell_volume']
    else:
      return None

  def get_vc_relax_density(self):
    if(self.data['type']=="VC-RELAX"):
      return self.data['crystals'][-1]['cell_density']
    else:
      return None
      
  def get_mass_per_crystal(self):
    return self.data['mass_per_crystal']
      
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
    print("=======================================================================")
    for key in sorted(self.data.keys()):
      value = self.data[key]
      print(key, ":  ", value)
    print("=======================================================================")
    print()
    
  def xyz_evang(self):
    self.xyz_units = 'ev/ang'
    self.stress_units = 'gpa'

  def xyz_stress_gpa(self):
    self.stress_units = 'gpa'

  def make_xyz(self, option=None):
    self.xyz = []
    if(option == None):
      for rn in range(len(self.data['results'])):
        option = rn + 1
    elif(option == -1):
      option = len(self.data['results'])
    else:
      option = (option - 1) % len(self.data['results']) + 1
    self.make_xyz_inner(option)
    return self.xyz
    
  def make_xyz_inner(self, option):
    if(len(self.data['results'])==0):
      return False
    if(len(self.data['crystals'])==0):
      return False  
   
    rn = (option - 1) % len(self.data['results'])
    cn = rn
    if(rn == len(self.data['results']) - 1):
      cn = len(self.data['crystals']) - 1
      
    crystal = self.data['crystals'][cn]
    result = self.data['results'][rn]
    settings = self.data['scf_settings']
    species = settings['atomic_species']
    
# Add new list and set counter n
    self.xyz.append([])
    n = len(self.xyz) - 1
    
# Add data
    self.xyz[n].append("#ALAT " + str(crystal['alat_adj']))
    self.xyz[n].append("#X " + str(crystal['cell_parameters_adj'][0][0]) + " " + str(crystal['cell_parameters_adj'][0][1]) + " " + str(crystal['cell_parameters_adj'][0][2]))
    self.xyz[n].append("#Y " + str(crystal['cell_parameters_adj'][1][0]) + " " + str(crystal['cell_parameters_adj'][1][1]) + " " + str(crystal['cell_parameters_adj'][1][2]))
    self.xyz[n].append("#Z " + str(crystal['cell_parameters_adj'][2][0]) + " " + str(crystal['cell_parameters_adj'][2][1]) + " " + str(crystal['cell_parameters_adj'][2][2]))
    
# Just use 1 1 1
    self.xyz[n].append("#C 2 2 2")
    self.xyz[n].append("#RCUT 6.5")
    self.xyz[n].append("#L_UNITS bohr")
    self.xyz[n].append("#E_UNITS ry")
    self.xyz[n].append("#F_UNITS ry/bohr")
    self.xyz[n].append("#S_UNITS kbar")
    
    self.xyz[n].append("#E " + str(result['energy']))
    
    if(result['s_on']):
      self.xyz[n].append("#SX " + str(result['stress'][0,0]) + " " +  str(result['stress'][0,1]) + " " + str(result['stress'][0,2]))
      self.xyz[n].append("#SY " + str(result['stress'][1,0]) + " " +  str(result['stress'][1,1]) + " " + str(result['stress'][1,2]))
      self.xyz[n].append("#SZ " + str(result['stress'][2,0]) + " " +  str(result['stress'][2,1]) + " " + str(result['stress'][2,2]))
      
    for label in species.keys():
      self.xyz[n].append("#M " + label + " " + str(species[label]['mass']))
    
    for i in range(self.atom_count):
      line = crystal['atomic_labels'][i]
      line = line + " " + str(crystal['crystal_positions'][i,0]) + " " + str(crystal['crystal_positions'][i,1]) + " " + str(crystal['crystal_positions'][i,2])      
      if(result['f_on']):   
        line = line + " " + str(result['forces'][i,0]) + " " + str(result['forces'][i,1]) + " " + str(result['forces'][i,2])
#"s_on": False,
      self.xyz[n].append(line)
      
# Return
    return self.xyz[n]
    
  def get_data(self, file=None): 
  
    """
# Important so store in it's own variable
    self.atom_count = 1    

# Control
    self.data = {
      "ok": False,
      "job_done": False,
      "error": False,
      "type": None,
      "summary": None,
      "mpi_processes": None,
      "threads_per_mpi_process": None,
      
      "scf_settings": None,      
      "crystals": [],
      "results": [],
      
      "initial_positions": None,   
      "total_energy": None,
      "density_full": None,
      "density": None,
      "stress": numpy.zeros((3,3)),
      "stress_sum": None,      
      "cpu_time": None,
      "wall_time": None,   
      "xyz": [],
    }
    """
    
    out = "##############################################################################################################\n"
    out = out + "atom count:                  " + str(self.atom_count) + "\n"
    out = out + "ok:                          " + str(self.data['ok']) + "\n"
    out = out + "job_done:                    " + str(self.data['job_done']) + "\n"
    out = out + "error:                       " + str(self.data['error']) + "\n"
    out = out + "type:                        " + str(self.data['type']) + "\n"
    out = out + "summary:                     " + str(self.data['summary']) + "\n"
    out = out + "mpi_processes:               " + str(self.data['mpi_processes']) + "\n"
    out = out + "threads_per_mpi_process:     " + str(self.data['threads_per_mpi_process']) + "\n"
    out = out + "scf_settings:                " + str(len(self.data['scf_settings'])) + "\n"
    
    for k in self.data['scf_settings'].keys():    
      out = out + "                             " + k + '  '+ str(self.data['scf_settings'][k]) + "\n"  
    
    out = out + "crystals (count):            " + str(len(self.data['crystals'])) + "\n"
    out = out + "results (count):             " + str(len(self.data['results'])) + "\n"
    i = 0
    for i in range(len(self.data['crystals'])-1):
      out = out + "##############################################################################################################\n"
      c = self.data['crystals'][i]
      n = i + 1
      if(n == 1):
        out = out + "crystal " + str(n) + " (input):\n"         
      elif(n == len(self.data['crystals']) - 1):
        out = out + "crystal " + str(n) + " (relaxed):\n" 
      else:
        out = out + "crystal " + str(n) + ":\n"  
      out = out + "##############################################################################################################\n"
      out = out + "alat:                        " + str(c['alat']) + "\n"  
      out = out + "cp:                          " + str(c['cell_parameters'][0,0]) + " " + str(c['cell_parameters'][0,1]) + " " + str(c['cell_parameters'][0,2]) + "\n" 
      out = out + "                             " + str(c['cell_parameters'][1,0]) + " " + str(c['cell_parameters'][1,1]) + " " + str(c['cell_parameters'][1,2]) + "\n"
      out = out + "                             " + str(c['cell_parameters'][2,0]) + " " + str(c['cell_parameters'][2,1]) + " " + str(c['cell_parameters'][2,2]) + "\n" 
      out = out + "alat adjusted:               " + str(c['alat_adj']) + "\n"  
      out = out + "cp adjusted:                 " + str(c['cell_parameters_adj'][0,0]) + " " + str(c['cell_parameters_adj'][0,1]) + " " + str(c['cell_parameters_adj'][0,2]) + "\n" 
      out = out + "                             " + str(c['cell_parameters_adj'][1,0]) + " " + str(c['cell_parameters_adj'][1,1]) + " " + str(c['cell_parameters_adj'][1,2]) + "\n"
      out = out + "                             " + str(c['cell_parameters_adj'][2,0]) + " " + str(c['cell_parameters_adj'][2,1]) + " " + str(c['cell_parameters_adj'][2,2]) + "\n" 
      if(i<len(self.data['crystals'])-2):
        s = self.data['results'][i]
        out = out + "energy:                      " + str(s['energy']) + "\n"  
        out = out + "total_force:                 " + str(s['total_force']) + "\n"  
        out = out + "stress:                      " + str(s['stress'][0,0]) + " " + str(s['stress'][0,1]) + " " + str(s['stress'][0,2]) + "\n"  
        out = out + "stress:                      " + str(s['stress'][1,0]) + " " + str(s['stress'][1,1]) + " " + str(s['stress'][1,2]) + "\n"  
        out = out + "stress:                      " + str(s['stress'][2,0]) + " " + str(s['stress'][2,1]) + " " + str(s['stress'][2,2]) + "\n"   
        out = out + "alat positions,  crystal positions,  forces:   " + "\n" 
        for m in range(self.atom_count):
          out = out + "                             " + str(c['atomic_positions'][m,0]) + " " + str(c['atomic_positions'][m,1]) + " " + str(c['atomic_positions'][m,2]) + "    " + str(c['crystal_positions'][m,0]) + " " + str(c['crystal_positions'][m,1]) + " " + str(c['crystal_positions'][m,2]) + "    " + str(s['forces'][m,0]) + " " + str(s['forces'][m,1]) + " " + str(s['forces'][m,2]) + "\n"
      
      out = out +     "cell volume:                 " + str(c['cell_volume']) + "\n"  
      out = out +     "density:                     " + str(c['cell_density']) + "\n" 
       
    c = self.data['crystals'][-1]
    s = self.data['results'][-1]
    n = len(self.data['crystals'])
    
    out = out + "##############################################################################################################\n"
    out = out + "crystal " + str(n) + " (final): \n"     
    out = out + "##############################################################################################################\n"
    out = out + "alat:                        " + str(c['alat']) + "\n" 
    out = out + "cp:                          " + str(c['cell_parameters'][0,0]) + " " + str(c['cell_parameters'][0,1]) + " " + str(c['cell_parameters'][0,2]) + "\n" 
    out = out + "                             " + str(c['cell_parameters'][1,0]) + " " + str(c['cell_parameters'][1,1]) + " " + str(c['cell_parameters'][1,2]) + "\n"
    out = out + "                             " + str(c['cell_parameters'][2,0]) + " " + str(c['cell_parameters'][2,1]) + " " + str(c['cell_parameters'][2,2]) + "\n" 
    out = out + "alat adjusted:               " + str(c['alat_adj']) + "\n"  
    out = out + "cp adjusted:                 " + str(c['cell_parameters_adj'][0,0]) + " " + str(c['cell_parameters_adj'][0,1]) + " " + str(c['cell_parameters_adj'][0,2]) + "\n" 
    out = out + "                             " + str(c['cell_parameters_adj'][1,0]) + " " + str(c['cell_parameters_adj'][1,1]) + " " + str(c['cell_parameters_adj'][1,2]) + "\n"
    out = out + "                             " + str(c['cell_parameters_adj'][2,0]) + " " + str(c['cell_parameters_adj'][2,1]) + " " + str(c['cell_parameters_adj'][2,2]) + "\n" 

    out = out + "energy:                      " + str(s['energy']) + "\n"  
    out = out + "total_force:                 " + str(s['total_force']) + "\n"  
    out = out + "stress:                      " + str(s['stress'][0,0]) + " " + str(s['stress'][0,1]) + " " + str(s['stress'][0,2]) + "\n"  
    out = out + "stress:                      " + str(s['stress'][1,0]) + " " + str(s['stress'][1,1]) + " " + str(s['stress'][1,2]) + "\n"  
    out = out + "stress:                      " + str(s['stress'][2,0]) + " " + str(s['stress'][2,1]) + " " + str(s['stress'][2,2]) + "\n"   
    out = out + "alat positions,  crystal positions,  forces:   " + "\n" 
    for m in range(self.atom_count):
      out = out + "                             " + str(c['atomic_positions'][m,0]) + " " + str(c['atomic_positions'][m,1]) + " " + str(c['atomic_positions'][m,2]) + "    " + str(c['crystal_positions'][m,0]) + " " + str(c['crystal_positions'][m,1]) + " " + str(c['crystal_positions'][m,2]) + "    " + str(s['forces'][m,0]) + " " + str(s['forces'][m,1]) + " " + str(s['forces'][m,2]) + "\n"
         
    out = out +     "cell volume:                 " + str(c['cell_volume']) + "\n"  
    out = out +     "density:                     " + str(c['cell_density']) + "\n" 
    
    fh = open(file,'w')
    fh.write(out)
    fh.close()
    
    return out
    
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
    
  @staticmethod
  def cell_volume(alat, cp):
    cp_alat = numpy.zeros((3,3,),)
    cp_alat[:,:] = alat * cp[:,:]
    v = numpy.dot(cp_alat[0,:],numpy.cross(cp_alat[1,:], cp_alat[2,:]))
    return v
    
  @staticmethod
  def cell_density(v, mpc):
    v_m3 = v * 1.48036e-31
    m = mpc * 1.66054E-027
    rho = m / v_m3
    return rho
    
"""
        
  def aaa():
    
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
              self.data['stress'][j,0] = fields[3] 
              self.data['stress'][j,1] = fields[4] 
              self.data['stress'][j,2] = fields[5]
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
"""

###########################################
#  CLASS pwscf_standar
###########################################
class pwscf_standard:

  @staticmethod
  def sc(label, size=1):
    c_atoms = 4
    n_atoms = c_atoms * size**3
    if(not isinstance(label, (list,))):
      label = [label]
    atoms = ['crystal']
    for x in range(size):
      for y in range(size):
        for z in range(size):
          coords = [[str((x+0.0)/size), str((y+0.0)/size), str((z+0.0)/size)]]
          for i in range(len(coords)):
            atoms.append([label[i % len(label)],coords[i][0],coords[i][1],coords[i][2]])
    return atoms, c_atoms, n_atoms   

  @staticmethod
  def bcc(label, size=1):
    c_atoms = 4
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
  def run():
    print("Running")

# Get current working directory
    cwd = os.getcwd()  

# Files list
    files = []

# Read file list
    files = pwscf_exec.make_file_list(cwd, files)

# Run
    out = pwscf_exec.execute(files)
   
  @staticmethod
  def make_file_list(path, files):
    for file_name in os.listdir(path):
      full_path = path + "/" + file_name
      if(os.path.isfile(full_path)):
        if(pwscf_input.is_pwscf(full_path)):
          print(full_path)
          files.append(full_path)
      if(os.path.isdir(full_path)):
        files = pwscf_exec.make_file_list(full_path, files)
    return files

  @staticmethod
  def execute(files_in, log=None, set_dirs=True, allow_cache=True):
    files = []
    if type(files_in) is str:
      files.append(files_in)
    elif type(files_in) is list:
      files = files_in
    else:
      return log
      
# Prepare run list
    run_list = pwscf_exec.pre_run(files_in, set_dirs)
  
# Run
    log, list_out, run_list = pwscf_exec.run_pwscf_list(run_list, log, allow_cache, True)

# Save
    pwscf_exec.save_run_list(run_list)
    
# Return log
    return log, list_out, run_list

  @staticmethod
  def pre_run(files_in, set_dirs):
# Create run list
    run_list = []
    
# Get environment settings
    s = pwscf_settings.load()    
    cache = None
    if('pwscf_cache' in s):
      cache = s['pwscf_cache']
    pwscf_scratch = None
    if('pwscf_scratch' in s):
      pwscf_scratch = s['pwscf_scratch']
    pwscf_pp = None
    if('pwscf_cache' in s):
      pwscf_pp = s['pwscf_pp']
     
    for file_path in files_in:  
      file_name = pwscf_exec.file_name_only(file_path)
      path = pwscf_exec.file_path_only(file_path)
    
      pw_in = pwscf_input()
      pw_in.load(file_path)
      pw_in.set_dirs(pwscf_scratch, pwscf_pp)
      pw_in.save()

      sig = pw_in.signature()
      cache_file = None
      
      use_cache = False
      if(cache is not None):
        cache_file_in = cache + "/" + sig + ".in"
        cache_file_out = cache + "/" + sig + ".out"
        
        use_cache = True
        if (not os.path.exists(cache_file_in)):
          use_cache = False
        if (not os.path.exists(cache_file_out)):
          use_cache = False
      
      cmd = 'mpirun -n ' + s['proc_count'] + ' ' + s['pwscf_bin'] + ' -i ' + file_path + ' > ' + path + '/' + file_name + '.out'
      
      run_list.append({
        'file_name': file_name, 
        'path': path, 
        'file_path_in': file_path, 
        'file_path_out': path + '/' + file_name + '.out', 
        'hash': sig,
        'use_cache': use_cache,
        'cache_in': cache_file_in,
        'cache_out': cache_file_out,
        'cmd': cmd,
        'cache_used': False,
        })
      
# return
    return run_list
    
  @staticmethod
  def run_pwscf_list(run_list, log, allow_cache, run=True):
# Get environment settings
    s = pwscf_settings.load()    
    cache = None
    if('pwscf_cache' in s):
      cache = s['pwscf_cache']
      
    list_out = []
    for i in range(len(run_list)):
      run_file = run_list[i]
    
# Use cache?
      use_cache = run_file['use_cache']
    
# Recheck for cache
      if(cache is not None):
        if(os.path.exists(run_file['cache_in']) and os.path.exists(run_file['cache_out'])):
          use_cache = True        
      
      if(use_cache and allow_cache):
        run_list[i]['cache_used'] = True
        copyfile(run_file['cache_out'], run_file['file_path_out'])
        list_out.append({"file": run_file['file_path_out'], "status": "complete",})
        if(log is not None):
          log.add("##START PWSCF##")
          log.add("PWscf Run")
          log.add(run_file['cmd'])
          log.add("Cache used")
          log.add("##END PWSCF##")
          log.add("")
      else:  
      
# Run
        if(run):
          os.system(run_file['cmd'])
          pwo = pwscf_output(run_file['file_path_out'])
          if(pwo.get_job_done()):
            copyfile(run_file['file_path_in'], run_file['cache_in'])
            copyfile(run_file['file_path_out'], run_file['cache_out'])
            list_out.append({"file": run_file['file_path_out'], "status": "complete",})
          else:
            list_out.append({"file": run_file['file_path_out'], "status": "failed",})          
        
          if(log is not None):
            log.add("##START PWSCF##")
            log.add("PWscf Run")
            log.add(run_file['cmd'])
            if(pwo.get_job_done()):
              log.add("Run successful")
            log.add("##END PWSCF##")
            log.add("")

    return log, list_out, run_list
      
  @staticmethod
  def file_name_only(file_path):
    l1 = file_path.split("/")    
    l2 = l1[-1].split(".")
    file_name = l1[-1][:-(len(l2)+1)]
    return file_name
    
  @staticmethod
  def file_path_only(file_path):
    l1 = file_path.split("/")   
    path = file_path[:-(len(l1[-1])+1)]
    return path

  @staticmethod
  def save_run_list(run_list):  
# Save Run List
    cwd = os.getcwd() 

    fh = open(cwd + "/run_list.txt", "w")

    for r in run_list:
      field = "File Name:"
      while(len(field)<32):
        field = field + " "
      fh.write(field + r['file_name'] + "\n")
      field = "Path:"
      while(len(field)<32):
        field = field + " "
      fh.write(field + r['path'] + "\n")
      field = "File In:"
      while(len(field)<32):
        field = field + " "
      fh.write(field + r['file_path_in'] + "\n")
      field = "File Out:"
      while(len(field)<32):
        field = field + " "
      fh.write(field + r['file_path_out'] + "\n")
      field = "Hash:"
      while(len(field)<32):
        field = field + " "
      fh.write(field + r['hash'] + "\n")

      field = "Use Cache?:"
      while(len(field)<32):
        field = field + " "
      if(r['use_cache']):
        fh.write(field + "True\n")
      else:
        fh.write(field + "False\n")

      field = "Cache Used:"
      while(len(field)<32):
        field = field + " "
      if(r['cache_used']):
        fh.write(field + "True\n")
      else:
        fh.write(field + "False\n")

      field = "Cache In:"
      while(len(field)<32):
        field = field + " "
      fh.write(field + r['cache_in'] + "\n")

      field = "Cache Out:"
      while(len(field)<32):
        field = field + " "
      fh.write(field + r['cache_out'] + "\n")

      field = "CMD:"
      while(len(field)<32):
        field = field + " "
      fh.write(field + r['cmd'] + "\n")
      fh.write("\n")
      fh.write("\n")

    fh.close()

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
#  CLASS read_confi
###########################################
class read_config:
  
  @staticmethod
  def read_file(file_path):
  
# Input dictionary
    input = {}
  
# READ DATA
    d = []
    fh = open(file_path, 'r')
    for line in fh:
      line = line.strip()
      if(len(line) > 0 and line[0] != "#"):
        d.append(line)      
    fh.close()
    
# Count commands
    commands = {}
    for line in d:
      fields = read_config.split_by(line, ' ')
      c = fields[0].lower()
      if(c in commands.keys()):
        commands[c] = commands[c] + 1
      else:
        commands[c] = 1
        
# Prepare input dictionary
    for k in commands.keys():
      if(commands[k] == 1):
        input[k] = None
      else:
        input[k] = []
    
# Read Data into input
    for line in d:
      fields = read_config.split_by(line, ' ')
      fkey = fields[0].lower()
      
      fd_size = {}
      for i in range(1, len(fields)):
        f = fields[i]
        fs = f.split("=")
        fc = fs[0].lower()
        if(fc in fd_size.keys()):
          fd_size[fc] = fd_size[fc] + 1
        else:
          fd_size[fc] = 1
          
# Prepare dictionary
      fd = {} 
      for k in fd_size.keys():
        if(fd_size[k] == 1):
          fd[k] = None
        else:
          fd[k] = []        
        
      for i in range(1, len(fields)):
        f = fields[i]
        fs = f.split("=")     
        fc = fs[0].lower()        
        fs = read_config.split_by(fs[1], ',')         
        fs = read_config.store(fs)
        
        if(fd_size[fc] == 1):
          if(len(fs) == 1):
            fd[fc] = read_config.store(fs[0])
          else:
            fd[fc] = read_config.store(fs)
        else:
          if(len(fs) == 1):
            fd[fc].append(read_config.store(fs[0]))
          else:
            fd[fc].append(read_config.store(fs))
            
      if(commands[fkey] == 1):
        input[fkey] = fd
      else:
        input[fkey].append(fd)  

    return input
        
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
    
  @staticmethod
  def store(inp):  
    if(isinstance(inp, list)):
      for i in range(len(inp)):
        try:
          if('.' in inp[i]  or 'e' in inp[i]):
            inp[i] = float(inp[i])
          else:
            inp[i] = int(inp[i])
        except:
          pass
    else:
      try:
        if('.' in inp or 'e' in inp):
          inp = float(inp)
        else:
          inp = int(inp)
      except:
        pass
    return inp
      
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
    
# CHARGE DENSITY (UNIT CHARGE PER VOLUME - ANG^3)
    charge_density = {
    'ANG-3': 1.0,
    'BOHR-3': 0.14812,    
    }
    
    unit_list = [length, energy, force, velocity, pressure, charge_density]
    
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
#  CLASS global
###########################################
class globals: 
  
  dirs = {
         'wd': 'wd',
         'td': 'wd/template',
         'ecut': 'wd/ecut',
         'kpoints': 'wd/kpoints',
         'runpwscf': 'wd/runpwscf',
         'plots': 'wd/plots',
         'csv': 'wd/csv',  
         'logs': 'wd/logs',  
         }
  
  times = {
          'start' : 0.0,
          'end' : 0.0,
          'duration' : 0.0,
          }
          
  inp = {}    
  
  pw_units = {
             'energy': 'RY',
             'length': 'BOHR',
             'force': 'RY/BOHR',
             'charge_density': 'BOHR-3',
             } 
          
  d = {
      'input_file': 'input.in',
      'pwscf_template': '',       
         
      }
      
  ecut = {
         'c_force_units': 'RY/BOHR',
         'c_force': 1.0E-5,
         'c_force_pw': None,
         'c_energy_units': 'RY/BOHR',
         'c_energy': 1.0E-5,
         'c_energy_pw': None,
         
         'run': True,
         'wfc_units': 'RY',
         'rho_units': 'BOHR-3',
         'wfc_min': 50,
         'wfc_max': 200,         
         'wfc_inc': 5, 
         'wfc_dec': 1,
         'rho_dec': 5,         
         'wfc_min_pw': None,
         'wfc_max_pw': None,         
         'wfc_inc_pw': None, 
         'wfc_dec_pw': None,
         'rho_dec_pw': None,
         'wr_min_ratio': 2.0,
         
         'kpoints_type': 'automatic',
         'kpoints_val': '5 5 5 1 1 1',
         
         'type': 'sc',
         'size': 1,
         'alat': 1.0,
         'alat_units': 'bohr',
         'alat_pw': None,
         'alat_expanded_pw': None,
         'cp': numpy.zeros((6,),dtype=numpy.float64,), 
         'rand_variance': 0.01,
         'rand_seed': 0,
         
         'data': numpy.zeros((20,1000),dtype=numpy.float64,), 
         'data_relative': numpy.zeros((20,1000),dtype=numpy.float64,), 
         'data_ev': numpy.zeros((20,1000),dtype=numpy.float64,), 
         'data_len': numpy.zeros((3),dtype=numpy.int32,), 
         
         'ecutwfc_converged_pw' : None,
         'ecutrho_converged_pw' : None,
         
         }
         
  ecut2d = {       
           'wfc_units': 'RY',
           'rho_units': 'BOHR-3',
           'wfc_min': 40,
           'wfc_max': 70,
           'wfc_inc': 5,
           'rho_min': 100,
           'rho_max': 300,
           'rho_inc': 20,
           'wfc_min_pw': None,
           'wfc_max_pw': None,
           'wfc_inc_pw': None,
           'rho_min_pw': None,
           'rho_max_pw': None,
           'rho_inc_pw': None,
           'energy_ry': numpy.zeros((200,200),dtype=numpy.float64,), 
           'force_rybohr': numpy.zeros((200,200),dtype=numpy.float64,), 
           'energy_ev': numpy.zeros((200,200),dtype=numpy.float64,), 
           'force_evang': numpy.zeros((200,200),dtype=numpy.float64,), 
           'data_w': 0,
           'data_h': 0,
           }
      
  kconv = {
          'ecutwfc': None,
          'ecutrho': None,
          
          'run': True,
          'k_min': 3,
          'k_max': 7,
          'k_inc': 2,
          'k_type': 'automatic',
          'k_offset': [1,1,1],
          'smearing': [0.001,0.002,0.003,0.004,0.005],
          'smearing_pw': [],
          'smearing_units': 'RY',
          'type': 'sc',
          'size': 1,
          'alat': 1.0,
          'alat_units': 'bohr',
          'alat_pw': None,
          'alat_expanded_pw': None,
          'cp': numpy.zeros((6,),dtype=numpy.float64,), 
          'rand_variance': 0.01,
          'rand_seed': 0,
          
          'energy_ry': numpy.zeros((200,200),dtype=numpy.float64,), 
          'force_rybohr': numpy.zeros((200,200),dtype=numpy.float64,), 
          'energy_ev': numpy.zeros((200,200),dtype=numpy.float64,), 
          'force_evang': numpy.zeros((200,200),dtype=numpy.float64,), 
          'data_w': 0,
          'data_h': 0,
          }
         
  file_counter = 0   
  
  log_fh = None
         
  def file_name():
    globals.file_counter = globals.file_counter + 1
    name = "file_"
    file_counter_str = str(globals.file_counter)
    while(len(file_counter_str) < 6):
      file_counter_str = '0' + file_counter_str
    name = name + file_counter_str    
    return name
         
###########################################
#  CLASS read_inpu
###########################################
class read_input:

  @staticmethod
  def run():
    globals.log_fh.write('\n')
    globals.log_fh.write('READ INPUT \n')
    globals.log_fh.write('========== \n')
  
# Read input
    if(not os.path.isfile(globals.d['input_file'])):
      globals.log_fh.write('No input file \n')
      exit()
    globals.inp = read_config.read_file(globals.d['input_file'])
    
    try:
      globals.d['pwscf_template'] = globals.inp['template']['file']
    except:
      pass      
    
    if(not os.path.isfile(globals.d['pwscf_template'])):
      globals.log_fh.write('No template file \n')
      exit()

# CONVERGE
    try:
      globals.ecut['c_force_units'] = globals.inp['converge']['force_units']
    except:
      pass
    try:
      globals.ecut['c_force'] = globals.inp['converge']['force']
    except:
      pass
    try:
      globals.ecut['c_energy_units'] = globals.inp['converge']['energy_units']
    except:
      pass
    try:
      globals.ecut['c_energy'] = globals.inp['converge']['energy']
    except:
      pass 
      
# ECUTCONV
    try:
      globals.ecut['run'] = std.option(globals.inp['ecutconv']['run'])
    except:
      pass
    try:
      globals.ecut['wfc_units'] = globals.inp['ecutconv']['wfc_units']
    except:
      pass
    try:
      globals.ecut['rho_units'] = globals.inp['ecutconv']['rho_units']
    except:
      pass
    try:
      globals.ecut['wfc_min'] = globals.inp['ecutconv']['wfc_min']
    except:
      pass
    try:
      globals.ecut['wfc_max'] = globals.inp['ecutconv']['wfc_max']
    except:
      pass
    try:
      globals.ecut['wfc_inc'] = globals.inp['ecutconv']['wfc_inc']
    except:
      pass
    try:
      globals.ecut['wfc_dec'] = globals.inp['ecutconv']['wfc_dec']
    except:
      pass
    try:
      globals.ecut['rho_dec'] = globals.inp['ecutconv']['rho_dec']
    except:
      pass
    try:
      globals.ecut['wr_min_ratio'] = globals.inp['ecutconv']['wr_min_ratio']
    except:
      pass
    try:
      globals.ecut['kpoints_type'] = globals.inp['ecutconv']['kpoints_type']
    except:
      pass
    try:
      globals.ecut['kpoints_val'] = globals.inp['ecutconv']['kpoints_val']
    except:
      pass
      
# ECUT CRYSTAL
    try:
      globals.ecut['type'] = globals.inp['configecut']['type']
    except:
      pass
    try:
      globals.ecut['size'] = globals.inp['configecut']['size']
    except:
      pass
    try:
      globals.ecut['alat'] = globals.inp['configecut']['alat']
    except:
      pass
    try:
      globals.ecut['alat_units'] = globals.inp['configecut']['alat_units']
    except:
      pass
    try:
      for i in range(6):
        globals.ecut['cp'][i] = globals.inp['configecut']['cp'][i]
    except:
      pass
    try:
      globals.ecut['rand_variance'] = globals.inp['configecut']['rand_variance']
    except:
      pass
    try:
      globals.ecut['rand_seed'] = globals.inp['configecut']['rand_seed']
    except:
      pass
    
# KPOINTSCONV:configkpoint
    try:
      globals.kconv['run'] = std.option(globals.inp['kconv']['run'])
    except:
      pass
    try:
      globals.kconv['k_min'] = globals.inp['kconv']['k_min']
    except:
      pass   
    try:
      globals.kconv['k_inc'] = globals.inp['kconv']['k_inc']
    except:
      pass 
    try:
      globals.kconv['k_max'] = globals.inp['kconv']['k_max']
    except:
      pass 
    try:
      globals.kconv['k_type'] = globals.inp['kconv']['k_type']
    except:
      pass 
    try:
      globals.kconv['smearing'] = globals.inp['kconv']['smearing']
    except:
      pass 
    try:
      globals.kconv['smearing_units'] = globals.inp['kconv']['smearing_units']
    except:
      pass 
      
# KPOINTS CRYSTAL
    try:
      globals.kconv['type'] = globals.inp['configkpoint']['type']
    except:
      pass
    try:
      globals.kconv['size'] = globals.inp['configkpoint']['size']
    except:
      pass
    try:
      globals.kconv['alat'] = globals.inp['configkpoint']['alat']
    except:
      pass
    try:
      globals.kconv['alat_units'] = globals.inp['configkpoint']['alat_units']
    except:
      pass
    try:
      for i in range(6):
        globals.kconv['cp'][i] = globals.inp['configkpoint']['cp'][i]
    except:
      pass
    try:
      globals.kconv['rand_variance'] = globals.inp['configkpoint']['rand_variance']
    except:
      pass
    try:
      globals.kconv['rand_seed'] = globals.inp['configkpoint']['rand_seed']
    except:
      pass
    
    """
          'run': True,
          'k_min': 3,
          'k_max': 7,
          'k_inc': 2,
          'k_type': 'automatic',
          'smearing': [0.001,0.002,0.003,0.004,0.005]
          
          'type': 'sc',
          'size': 1,
          'alat': 1.0,
          'alat_units': 'bohr',
          'cp': numpy.zeros((6,),dtype=numpy.float64,), 
          'rand_variance': 0.01,
          'rand_seed': 0,
         }
    """
      
# UNIT CONVERSIONS
    globals.ecut['c_force_pw'] = units.convert(globals.ecut['c_force_units'], globals.pw_units['force'], globals.ecut['c_force'])
    globals.ecut['c_energy_pw'] = units.convert(globals.ecut['c_energy_units'], globals.pw_units['energy'], globals.ecut['c_energy'])
    globals.ecut['wfc_min_pw'] = units.convert(globals.ecut['wfc_units'], globals.pw_units['energy'], globals.ecut['wfc_min'])
    globals.ecut['wfc_max_pw'] = units.convert(globals.ecut['wfc_units'], globals.pw_units['energy'], globals.ecut['wfc_max'])
    globals.ecut['wfc_inc_pw'] = units.convert(globals.ecut['wfc_units'], globals.pw_units['energy'], globals.ecut['wfc_inc'])
    globals.ecut['wfc_dec_pw'] = units.convert(globals.ecut['wfc_units'], globals.pw_units['energy'], globals.ecut['wfc_dec'])
    globals.ecut['rho_dec_pw'] = units.convert(globals.ecut['rho_units'], globals.pw_units['charge_density'], globals.ecut['rho_dec'])
    globals.ecut['alat_pw'] = units.convert(globals.ecut['alat_units'], globals.pw_units['length'], globals.ecut['alat'])    
    globals.ecut['alat_expanded_pw'] = globals.ecut['size'] * globals.ecut['alat_pw']
            
    globals.ecut2d['wfc_min_pw'] = units.convert(globals.ecut2d['wfc_units'], globals.pw_units['energy'], globals.ecut2d['wfc_min'])
    globals.ecut2d['wfc_max_pw'] = units.convert(globals.ecut2d['wfc_units'], globals.pw_units['energy'], globals.ecut2d['wfc_max'])
    globals.ecut2d['wfc_inc_pw'] = units.convert(globals.ecut2d['wfc_units'], globals.pw_units['energy'], globals.ecut2d['wfc_inc'])
    globals.ecut2d['rho_min_pw'] = units.convert(globals.ecut2d['rho_units'], globals.pw_units['charge_density'], globals.ecut2d['rho_min'])
    globals.ecut2d['rho_max_pw'] = units.convert(globals.ecut2d['rho_units'], globals.pw_units['charge_density'], globals.ecut2d['rho_max'])
    globals.ecut2d['rho_inc_pw'] = units.convert(globals.ecut2d['rho_units'], globals.pw_units['charge_density'], globals.ecut2d['rho_inc'])
    
    for i in range(len(globals.kconv['smearing'])):
      globals.kconv['smearing_pw'].append(units.convert(globals.kconv['smearing_units'], globals.pw_units['energy'], globals.kconv['smearing'][i]))
    globals.kconv['alat_pw'] = units.convert(globals.kconv['alat_units'], globals.pw_units['length'], globals.kconv['alat'])  
    globals.kconv['alat_expanded_pw'] = globals.kconv['size'] * globals.kconv['alat_pw']    
    
    globals.log_fh.write('Input files \n')  
    globals.log_fh.write('Input:      '+ globals.d['input_file'] + '\n')
    globals.log_fh.write('Template:   '+ globals.d['pwscf_template'] + '\n')
    globals.log_fh.write('\n')  
      
    globals.log_fh.write('Ecut \n')
    for k in globals.ecut.keys():
      globals.log_fh.write(str(k) + ':   ' + str(globals.ecut[k]) + '\n')
    globals.log_fh.write('\n')
      
    globals.log_fh.write('Ecut 2D \n')
    for k in globals.ecut2d.keys():
      globals.log_fh.write(str(k) + ':   ' + str(globals.ecut2d[k]) + '\n')
    globals.log_fh.write('\n')  
      
    globals.log_fh.write('Kpoint\n')
    for k in globals.kconv.keys():
      globals.log_fh.write(str(k) + ':   ' + str(globals.kconv[k]) + '\n')
    globals.log_fh.write('\n')
       
###########################################
#  CLASS conv_ecu
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
    ef.load_config(globals.ecut['type'], globals.ecut['size'], globals.ecut['alat_expanded_pw'])
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
      
      n = globals.ecut['data_len'][0]
      globals.ecut['data_len'][0] = globals.ecut['data_len'][0] + 1

      globals.ecut['data'][0,n] = ecutwfc
      globals.ecut['data'][1,n] = 4 * ecutwfc
      globals.ecut['data'][2,n] = fpo.get_energy_per_atom()
      globals.ecut['data'][3,n] = fpo.get_force_per_atom()
     
      globals.log_fh.write(str(globals.ecut['data'][0,n]) + ' ' + str(globals.ecut['data'][1,n]) + ' ' + str(globals.ecut['data'][2,n]) + ' ' + str(globals.ecut['data'][3,n]) + ' ' + '\n')

      if(counter > 1):
        if(converged < 2):
          ecutwfc_v = ecutwfc - globals.ecut['wfc_inc_pw']
        if(abs(globals.ecut['data'][2,n] - globals.ecut['data'][2,n-1]) <= econv_threshold and abs(globals.ecut['data'][3,n] - globals.ecut['data'][3,n-1]) <= fconv_threshold):
          converged = converged + 1
        else:
          converged = 0

# Increment
      ecutwfc = ecutwfc + globals.ecut['wfc_inc_pw']   
      
# Record value
    globals.log_fh.write('Ecutwfc: ' + str(ecutwfc_v) + ' \n')
    globals.log_fh.write('Ecutrho: ' + str(4 * ecutwfc_v) + ' \n')
    
#################################
# 2. DECREASE ECUTRHO
#################################
    
    globals.log_fh.write('2. DECREASE ECUTRHO \n')
    
    converged = 0
    ecutwfc = ecutwfc_v
    ecutrho = 4.0 * ecutwfc_v
    ecutrho_min = globals.ecut['wr_min_ratio'] * ecutwfc
    
# CONVERGENCE FOR EACH CHANGE BY 1 RY
    econv_threshold = globals.ecut['c_energy_pw'] * globals.ecut['rho_dec']
    fconv_threshold = globals.ecut['c_force_pw'] * globals.ecut['rho_dec']
    
    globals.log_fh.write('econv_threshold: ' + str(econv_threshold) + ' \n')
    globals.log_fh.write('fconv_threshold: ' + str(fconv_threshold) + ' \n')
    
    ecutrho_v = None
    counter = 0   
    globals.log_fh.write('Data: \n')
    while((ecutrho >= ecutrho_min and converged < 1) or counter < 5):
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
      
      n = globals.ecut['data_len'][1]
      globals.ecut['data_len'][1] = globals.ecut['data_len'][1] + 1

      globals.ecut['data'][5,n] = ecutwfc
      globals.ecut['data'][6,n] = ecutrho
      globals.ecut['data'][7,n] = fpo.get_energy_per_atom()
      globals.ecut['data'][8,n] = fpo.get_force_per_atom()
     
      globals.log_fh.write(str(globals.ecut['data'][5,n]) + ' ' + str(globals.ecut['data'][6,n]) + ' ' + str(globals.ecut['data'][7,n]) + ' ' + str(globals.ecut['data'][8,n]) + ' ' + '\n')

      if(counter > 1):
        if(converged < 1):
          ecutrho_v = ecutrho + globals.ecut['rho_dec_pw']
        if(abs(globals.ecut['data'][7,n] - globals.ecut['data'][7,n-1]) <= econv_threshold and abs(globals.ecut['data'][8,n] - globals.ecut['data'][8,n-1]) <= fconv_threshold):
          converged = converged + 1
        else:
          converged = 0

# Decrease
      ecutrho = ecutrho - globals.ecut['rho_dec_pw']   
    
# Record value
    globals.log_fh.write('Ecutwfc: ' + str(ecutwfc_v) + ' \n')
    globals.log_fh.write('Ecutrho: ' + str(ecutrho_v) + ' \n')
    
#################################
# 3. DECREASE ECUTWFC
#################################
    
    globals.log_fh.write('3. DECREASE ECUTWFC \n')
    
    converged = 0
    ecutwfc = ecutwfc_v
    ecutrho = ecutrho_v

# CONVERGENCE FOR EACH CHANGE BY 1 RY
    econv_threshold = globals.ecut['c_energy_pw'] * globals.ecut['wfc_dec']
    fconv_threshold = globals.ecut['c_force_pw'] * globals.ecut['wfc_dec']
    
    globals.log_fh.write('econv_threshold: ' + str(econv_threshold) + ' \n')
    globals.log_fh.write('fconv_threshold: ' + str(fconv_threshold) + ' \n')
    
    counter = 0   
    globals.log_fh.write('Data: \n')
    while((ecutwfc >= globals.ecut['wfc_min_pw'] and converged < 1) or counter < 5):
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
      
      n = globals.ecut['data_len'][2]
      globals.ecut['data_len'][2] = globals.ecut['data_len'][2] + 1

      globals.ecut['data'][10,n] = ecutwfc
      globals.ecut['data'][11,n] = ecutrho
      globals.ecut['data'][12,n] = fpo.get_energy_per_atom()
      globals.ecut['data'][13,n] = fpo.get_force_per_atom()
     
      globals.log_fh.write(str(globals.ecut['data'][10,n]) + ' ' + str(globals.ecut['data'][11,n]) + ' ' + str(globals.ecut['data'][12,n]) + ' ' + str(globals.ecut['data'][13,n]) + ' ' + '\n')

      if(counter > 1):
        if(converged < 1):
          ecutwfc_v = ecutwfc + globals.ecut['wfc_dec_pw']
        if(abs(globals.ecut['data'][12,n] - globals.ecut['data'][12,n-1]) <= econv_threshold and abs(globals.ecut['data'][13,n] - globals.ecut['data'][13,n-1]) <= fconv_threshold):
          converged = converged + 1
        else:
          converged = 0

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

    plt.savefig(globals.dirs['plots'] + '/' + 'convergence_ry.svg')
  
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

    plt.savefig(globals.dirs['plots'] + '/' + 'convergence_adjusted_ry.svg')
  
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

    plt.savefig(globals.dirs['plots'] + '/' + 'convergence_ev.svg')
    
###########################################
#  CLASS conv_ecut2
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
    
###########################################
#  CLASS conv_kpoint
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
    numpy.savetxt(globals.dirs['csv'] + '/kpoints_energy_ry.csv', globals.kconv['energy_ry'][0:data_h,0:data_w], fmt='%10.5f', delimiter=',')
    numpy.savetxt(globals.dirs['csv'] + '/kpoints_energy_ev.csv', globals.kconv['energy_ev'][0:data_h,0:data_w], fmt='%10.5f', delimiter=',')
    numpy.savetxt(globals.dirs['csv'] + '/kpoints_force_ry.csv', globals.kconv['force_rybohr'][0:data_h,0:data_w], fmt='%10.5f', delimiter=',')
    numpy.savetxt(globals.dirs['csv'] + '/kpoints_force_ev.csv', globals.kconv['force_evang'][0:data_h,0:data_w], fmt='%10.5f', delimiter=',')

# PLOT
    
    x = numpy.linspace(globals.kconv['k_min'], globals.kconv['k_max'], data_w)
    y = numpy.zeros((len(globals.kconv['smearing_pw']),),)
    y[:] = globals.kconv['smearing_pw'][:]
    
    plt.clf()    
    plt.contourf(x, y, globals.kconv['energy_ry'][0:data_h,0:data_w], 20, cmap='Greys')
    plt.colorbar()
    plt.savefig(globals.dirs['plots'] + '/' + 'kpoints_energy_ry.svg')

    plt.clf()
    plt.contourf(x, y, globals.kconv['force_rybohr'][0:data_h,0:data_w], 20, cmap='Greys')
    plt.colorbar()
    plt.savefig(globals.dirs['plots'] + '/' + 'kpoints_force_rybohr.svg')

    x = numpy.linspace(globals.kconv['k_min'], globals.kconv['k_max'], data_w)
    for i in range(len(y)):
      y[i] = units.convert(globals.pw_units['energy'], 'EV', y[i])
    
    plt.clf()    
    plt.contourf(x, y, globals.kconv['energy_ev'][0:data_h,0:data_w], 20, cmap='Greys')
    plt.colorbar()
    plt.savefig(globals.dirs['plots'] + '/' + 'kpoints_energy_ev.svg')

    plt.clf()
    plt.contourf(x, y, globals.kconv['force_evang'][0:data_h,0:data_w], 20, cmap='Greys')
    plt.colorbar()
    plt.savefig(globals.dirs['plots'] + '/' + 'kpoints_force_evang.svg')

###########################################
###########################################
#  MAIN
###########################################
###########################################

new_run = pwscf_converge()

new_run.run()

