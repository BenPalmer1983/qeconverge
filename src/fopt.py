import numpy
import sys
import os
include = os.environ['PYCLASSES']
sys.path.append(include)
from genetic import genetic
from simulated_annealing import simulated_annealing
from lma import lma
from ng import ng



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
  
  
      