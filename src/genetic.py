#!/usr/bin/python3
################################################################################
#
#  Genetic
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
  
  






