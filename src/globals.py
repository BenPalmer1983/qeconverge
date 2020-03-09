##########
# GLOBALS

import numpy

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
         
         