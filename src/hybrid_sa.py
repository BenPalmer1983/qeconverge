#!/usr/bin/python3
################################################################################
#
#  Simulated Annealing
#
################################################################################

# 
# Example useage
#
# See End
#

import sys
import os
import random
import numpy

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
  
  