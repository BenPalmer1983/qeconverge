#!/usr/bin/python3
################################################################################
#
#  LMA Optimisation
#
################################################################################

# Calculates the variables of a function so that it minimises the sum
# of squares between the function and the data points
#
# LMA algorithm
# Levenberg, K.: 1944, Quarterly of Applied Mathematics 2, 164:
#    Dampened Newton Gauss Algorithm
# Marquardt, D.: 1963, Journal of the Society for Industrial and Applied Mathematics 11(2), 431:
#    Introduction of Lambda (Marquardt Parameter)
# R. Fletcher 1971, A modified Marquardt Subroutine for Non-linear Least Squares
#    Starting Lambda and lower lambda cutoff (where it is replaced by 0)
# M. Transtrum, J. Sethna 2012, Improvements to the Levenberg-Marquardt algorithm for nonlinear
# least-squares minimization:
#    Delayed gratification scheme for varying lambda
# Jens Jessen-Hansen 2011, Levenberg-Marquardts algorithm Project in Numerical Methods:
#    Independent/Diagonal covariance matrix for weighting
# Parameter limits implemented
# Finite difference method used to calculate J



# 
# Example useage
#
# See End
#

import sys
import os
import random
import numpy
import time


class lma:

  @staticmethod
  def run(f, p, settings_in={}, data = None):
    # Load Data
    data = lma.load(data)
  
    # Default Settings
    settings = {
    'conv': 1.0E-7,
    'h': 0.00001,
    'lcut': 0.1,
    'outer_cycles': 50,
    'inner_cycles': 100,
    'covariance': False,
    }
  
    # Load new settings
    for key in settings_in:
      if(key in settings):
        settings[key] = settings_in[key]
  
    # calculate
    return lma.calc(f, p, settings, data)
    
  
  @staticmethod
  def calc(f, p, settings, data):
    return lma.outer_cycle(f, p, settings, data)
    
    
  @staticmethod
  def outer_cycle(f, p, settings, data):
    # (JTJ+Lambda*diag(JTJ)) P = (-1*JTR)
    # (H+Lambda*diag(H)) P = (-1*JTR)
    i = 0
    outer_cycles = settings['outer_cycles']
    is_converged = False
    lam = 0
    lam_cut = 0
    
    while(i < outer_cycles):
      i = i + 1
      
      R = lma.residual(f, p, data)        # R
      J = lma.jacobian(f, p, data, settings)        # J
      W = lma.covariance(R, settings['covariance'])
      JT = numpy.transpose(J)
      JTW = numpy.matmul(JT, W)
      JTWJ = numpy.matmul(JTW, J)
      JTWR = numpy.matmul(JTW, R)
      rss = lma.calc_rss_from_r(R)
      
      if(i == 1):
        lam = lma.calc_lambda(J)
        lam_cut = lam
      # Run calc
      p_new, lam, new_rss = lma.inner_cycle(f, p, settings, data, J, JT, R, W, JTW, JTWJ, JTWR, lam, lam_cut, rss)
      
      if(new_rss < rss):
        p = p_new        
        is_converged = lma.converged(new_rss, rss, settings)
        rss = new_rss      
        if(is_converged):
          i = outer_cycles      
    return p, rss         
       
    
   
  @staticmethod
  def inner_cycle(f, p, settings, data, J, JT,  R, W, JTW, JTWJ, JTWR, lam, lam_cut, best_rss):
    inner_cycles = settings['inner_cycles']
    for n in range(inner_cycles):  
      try:
        p_new = lma.lma_calc(p, J, JT, R, W, JTW, JTWJ, JTWR, lam, lam_cut) 
        rss = lma.calc_rss(f, p_new, data)
        if(rss < best_rss):
          p = p_new
          lam = lam * 0.2
          best_rss = rss
        else:
          if(lam < lam_cut):
            lam = 2.0 * lam_cut
          lam = lam * 1.5
      except:    
        if(lam < lam_cut):
          lam = 2.0 * lam_cut
          lam = lam * 1.5 
    return p, lam, best_rss
  
  @staticmethod
  def converged(new_rss, last_rss, settings):    
    if(new_rss == 0.0 or abs(new_rss - last_rss) <= settings['conv']):
      return True
    return False
      
  @staticmethod
  def lma_calc(p, J, JT, R, W, JTW, JTWJ, JTWR, lam, lam_cut):
    JTWJ_diag = lam * numpy.diag(JTWJ)
    A = JTWJ + JTWJ_diag
    B = -1 * JTWR
    p_change = numpy.linalg.solve(A, B)
    return (p + p_change)
    
  @staticmethod
  def jacobian(f, p, data, settings):
    dl = len(data)
    pl = len(p)
    J = numpy.zeros((dl, pl), dtype=numpy.float64)    
    h = lma.step(p, settings['h'])
    for i in range(0, dl):
      for j in range(0, pl):
        # Reset parameters
        p_new = numpy.copy(p)
        # Vary jth parameter
        p_new[j] = p_new[j] + h[j]        
        # Calc J matrix
        J[i,j] = (f(p_new, data[i, 0]) - data[i, 1]) / h[j]    
    return J
  
  @staticmethod
  def jacobian_transpose(J):
    return numpy.transpose(self.J)
    
  @staticmethod
  def hessian(J, JT):
    return numpy.matmul(JT, J)

  @staticmethod
  def dampening(l, H):
    damp = numpy.identity(len(H), dtype=numpy.float64)
    for i in range(0,len(H)):
      damp[i,i] = l * H[i,i]    
  
  @staticmethod
  def step(p, dh):
    h = numpy.zeros((len(p)), dtype=numpy.float64)  
    for i in range(len(h)):
      if(p[i] == 0.0):
        h[i] = dh
      else:
        h[i] = dh * p[i]
    return h
    
  @staticmethod
  def residual(f, p, data):
    # Calculate residual
    return f(p, data[:, 0]) - data[:, 1] 
    
  def covariance(R, on=True):
    W = numpy.identity((len(R)), dtype=numpy.float64)
    if(on):
      build = True
      for i in range(len(R)):
        if(R[i] == 0.0):
          build = False
          break
      if(build):
        for i in range(len(R)):
          W[i,i] = 1.0 / R[i]**2
    return W
    
  @staticmethod
  def calc_lambda(J):
    JT = numpy.transpose(J)
    JTJ = numpy.matmul(JT, J)
    JTJ_inv = numpy.linalg.inv(JTJ)
    trace = numpy.trace(JTJ_inv)
    l = trace**(-1)
    return l

  @staticmethod
  def calc_rss(f, p, data):
    return sum((f(p, data[:, 0]) - data[:, 1])**2)

  @staticmethod
  def calc_rss_from_r(R):
    return sum(R**2)
  
###################################
 #  Load Data
 ###################################

  @staticmethod
  def load(file_in=None):
    try:
      if(not isinstance(file_in, (list,)) and not isinstance(file_in,numpy.ndarray) and file_in == None):
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
        data = numpy.zeros((len(file_in), 2), dtype=numpy.float64) 
        for row in range(len(file_in)):
          data[row,0] = float(file_in[row][0])
          data[row,1] = float(file_in[row][1])
        return data
      if(len(file_in) == 2 and len(file_in[0]) > 2):
        data = numpy.zeros((len(file_in[0]), 2), dtype=numpy.float64) 
        for row in range(len(file_in[0])):
          data[row,0] = float(file_in[0][row])
          data[row,1] = float(file_in[1][row])
        return data
    # If array   
    elif(isinstance(file_in,numpy.ndarray)): 
      data = numpy.zeros((len(file_in), 2), dtype=numpy.float64) 
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
      file_data = lma.clean(file_data)
      
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
      data = numpy.zeros((len(data_list), 2), dtype=numpy.float64) 
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
  