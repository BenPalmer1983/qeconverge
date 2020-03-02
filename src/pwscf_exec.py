################################################################
#    Processing PWscf input file
#
#
#
#
################################################################

import os
import datetime
import re
import sys
import hashlib 
from shutil import copyfile


from log import log
from pwscf_input import pwscf_input
from pwscf_output import pwscf_output
from pwscf_settings import pwscf_settings

############################
#  pwscf_input
############################


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




      
    
