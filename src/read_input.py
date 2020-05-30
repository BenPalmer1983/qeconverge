################################################################
#    Processing PWscf input file
#
#
#
#
################################################################


import numpy
import os
from units import units

###########################################
#  CLASS pwscf_converg
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
      
      
      
      
    try:
      globals.ecut2d['wfc_min'] = globals.inp['ecut2d']['wfc1']
    except:
      pass
    try:
      globals.ecut2d['wfc_max'] = globals.inp['ecut2d']['wfc2']
    except:
      pass
    try:
      globals.ecut2d['wfc_inc'] = globals.inp['ecut2d']['wcf3']
    except:
      pass
    try:
      globals.ecut2d['rho_min'] = globals.inp['ecut2d']['rho1']
    except:
      pass
    try:
      globals.ecut2d['rho_max'] = globals.inp['ecut2d']['rho2']
    except:
      pass
    try:
      globals.ecut2d['rho_inc'] = globals.inp['ecut2d']['rho3']
    except:
      pass    
    try:
      globals.ecut2d['wr_min_ratio'] = globals.inp['ecut2d']['wr_min_ratio']
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
       

