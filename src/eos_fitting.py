

import numpy
import random
import sys
import os
import matplotlib.pyplot as plt
from fopt import fopt
from lma import lma


class eos_fitting:

  
  @staticmethod
  def bm(data=None):
    return eos_fitting.eos_fit(eos_fitting.bm_calc, data)
    
  def pt(data=None):
    return eos_fitting.eos_fit(eos_fitting.pt_calc, data)
    
  def vinet(data=None):
    return eos_fitting.eos_fit(eos_fitting.vinet_calc, data)
    
  def all(data=None):
    data = eos_fitting.load(data)
    x, y = eos_fitting.eos_fit(eos_fitting.bm_calc, data) 
    plt.clf()
    plt.plot(data[:,0], data[:,1], 'r+', label="E-V")
    plt.plot(x, y, 'b-', label="E-V fit")
    plt.savefig('1.eps', format='eps')  
    x, y = eos_fitting.eos_fit(eos_fitting.bm_calc, data) 
    plt.clf()
    plt.plot(data[:,0], data[:,1], 'r+', label="E-V")
    plt.plot(x, y, 'b-', label="E-V fit")
    plt.savefig('2.eps', format='eps')  
    x, y = eos_fitting.eos_fit(eos_fitting.bm_calc, data) 
    plt.clf()
    plt.plot(data[:,0], data[:,1], 'r+', label="E-V")
    plt.plot(x, y, 'b-', label="E-V fit")
    plt.savefig('3.eps', format='eps')  

  @staticmethod
  def eos_fit(f, data=None, verbose=True):
    data = eos_fitting.load(data)
    V = data[:,0]
    E = data[:,1]
  
    eos = {"E0": 0.0, "V0": 0.0, "B0": 0.0, "B0P": 0.0}

    # 2nd order polynomial fit
    poly = numpy.polyfit(V, E, 2)
    # Starting points
    eos['V0'] = (-1 * poly[1]) / (2 * poly[0])
    eos['E0'] = (poly[0] * eos['V0'] * eos['V0']) + (poly[1] * eos['V0']) + poly[2]
    eos['B0'] = 2.0 * poly[0] * eos['V0']
    eos['B0P'] = 2.0
    
    p = numpy.zeros((4))
    p[0] = eos['V0']
    p[1] = eos['E0']
    p[2] = eos['B0']
    p[3] = eos['B0P']
    
    print(p)

    rss = fopt.calc_rss(f, p, data)
    print(rss)
    #p, rss = fopt.run(f, p, {}, data) 

    datafit = numpy.zeros((101,2))
    datafit[:,0] = numpy.linspace(data[0,0],data[-1,0], 101)
    datafit[:,1] = f(p, datafit[:,0])   
    return p, rss, datafit
   
  @staticmethod
  def bm_calc(p, V):
    V0 = p[0]
    E0 = p[1]
    B0 = p[2]
    B0P = p[3]
    eta = (V/V0)**(1/3.0)
    return E0 + (9/16.0) * (B0 * V0) * ((eta*eta - 1)*(eta*eta - 1)) * (6.0 + B0P * (eta * eta - 1) - 4 * eta * eta ) 
    
    
  @staticmethod
  def pt_calc(p, V):
    V0 = p[0]
    E0 = p[1]
    B0 = p[2]
    B0P = p[3]
    eta = (V/V0)
    return E0 + ((B0 * V0) / 2) * (numpy.ln(eta))**2 + ((B0 * V0) / 6) * (numpy.ln(eta))**3 * (B0P - 2)
    
  @staticmethod
  def vinet_calc(p, V):
    V0 = p[0]
    E0 = p[1]
    B0 = p[2]
    B0P = p[3]
    eta = (V/V0)**(1/3.0)
    return E0 + ((2 * B0 * V0) / (B0P - 1)**2) * (2 - (5 + 3 * eta * (B0P - 1) - 3 * B0P) * numpy.exp(-(3.0/2.0) * (B0P - 1) * (eta - 1)))

###################################
#  Load Data
###################################

#   Loads data into numpy array
#   Tries to read from file if file_in is set
#   If file_in is a list or numpy array it reads it in directly
#   If file_in is none, it checks if a file has been specified as an argument


  @staticmethod
  def load(file_in=None):
    if(not isinstance(file_in, (list,)) and not isinstance(file_in,numpy.ndarray) and file_in == None):
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
      file_data = eos_fitting.clean(file_data)
      
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
