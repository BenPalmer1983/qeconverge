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
  
  