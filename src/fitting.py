

import numpy
import random
import sys
import os
from hybrid_sa import hybrid_sa
from lma import lma

class fitting:

#################################
#   Exponential Fit
#################################

  @staticmethod
  def single_exp(p, x):
    return p[0] + p[1] * numpy.exp(p[2] * x + p[3])  

  @staticmethod
  def fit_exp(data = None, points = 101):
    result = {
    'p': numpy.zeros((1)),
    'rss': 1e99,
    'data': numpy.zeros((1))
    }
  
    # Load data
    data = fitting.load(data)
    y_min = numpy.min(data[:,1])
    y_max = numpy.max(data[:,1])
    dy = y_max - y_min
    y_mid = (y_min + y_max) / 2

    p = [y_mid,0,0,0]
    pl = [y_mid - abs(dy),-2 * abs(dy),-5,-1]
    pu = [y_mid + abs(dy),2 * abs(dy),5,1]
    pd = [abs(dy), abs(dy), 5, 1]
    p, rss = hybrid_sa.run(fitting.single_exp, p, {'p_diff': [abs(dy), abs(dy), 5, 1]})

    # LMA
    p, rss = lma.run(fitting.single_exp, p, {}, data)
    result['p'] = p    
    result['rss'] = rss    
    result['data'] = data
    result['data_fit'] = fitting.data_points(fitting.single_exp, p, points, data)
    return result
    







###################################
#  Load Data
###################################

#   Loads data into numpy array
#   Tries to read from file if file_in is set
#   If file_in is a list or numpy array it reads it in directly
#   If file_in is none, it checks if a file has been specified as an argument

  @staticmethod
  def data_points(f, p, points, data):
    data_out = numpy.zeros((points,2))
    data_out[:,0] = numpy.linspace(data[0,0],data[-1,0],points)
    data_out[:,1] = f(p, data_out[:,0])
    return data_out
    

  @staticmethod
  def load(file_in=None):
    if(file_in == None):
      if(len(sys.argv) > 1 and sys.argv[1] is not None):
        file_in = sys.argv[1]
        if(not os.path.isfile(file_in)):
          return None
      else:
        return None 

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
      file_data = fitting.clean(file_data)
      
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