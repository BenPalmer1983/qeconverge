################################################################
#    Processing PWscf output file
#
#
#
#
################################################################



import os
import datetime
import re
import sys
import numpy

############################
#  pwscf_output
############################

class pwscf_output:

  def __init__(self, file_in=None):
    self.reset()
    if(file_in != None):
      self.load(file_in)

  def reset(self):
  
    self.z = numpy.zeros((3,3))

    # Control
    self.data = {
      "ok": False,
      "job_done": False,
      "error": False,
      "type": None,
      "summary": None,
      "mpi_processes": None,
      "bravais_lattice_index": None,
      "alat": None,
      "volume": None,
      "electrons": None, 
      "electrons_up": None, 
      "electrons_down": None, 
      
      
      "crystal_in": numpy.zeros((3,3)),
      "crystal_calc": numpy.zeros((3,3)),
      
      
      "total_energy": None,
      "density_full": None,
      "density": None,
      "stress": numpy.zeros((3,3)),
      "stress_sum": None,
      
      "cpu_time": None,
      "wall_time": None
      
    }



  #  Load, and use in another program    
  def load(self, file_name): 
    # Load data from file
    data = self.load_from_file(file_name)
    
    # Read through data
    self.load_data(data)



  # Load from a block of data (text, file etc)
  def load_data(self, data):  
    
    # Reset data store
    self.reset()
    
    # split
    data = data.split("\n")    
    
    
    
    # OK
    ###################################
    self.data['ok'] = False
    counter = 0
    for line in data:
      line = line.strip()
      if(pwscf_output.compare(line, "JOB DONE.")):
        counter = counter + 1
      if(pwscf_output.compare(line, "Exit code:")):
        counter = counter - 1
      if(pwscf_output.compare(line, "convergence NOT achieved")):
        counter = counter - 1
    if(counter == 1):
      self.data['ok'] = True
    
    # Calc Type
    ###################################
    self.data['type'] = "SCF"
    for line in data:
      line = line.strip()
      if(line[0:23] == "A final scf calculation"):
        self.data['type'] = "VC-RELAX"
        
    
    # Load
    ###################################
    n = 0
    counter = 0
    while(n < len(data)):
      n, line, line_uc = self.next_line(n, data)
      if(line != ""):
        counter += 1
        if(counter == 1):
          self.data['summary'] = line
        else:
          if(pwscf_output.compare(line, "Number of MPI processes:")):
            self.data['mpi_processes'] = pwscf_output.extract(line, ":", "", "i") 
        
          if(pwscf_output.compare(line, "bravais-lattice index     =")):
            self.data['bravais_lattice_index'] = pwscf_output.extract(line, "=", "", "i") 
            
          if(pwscf_output.compare(line, "lattice parameter (alat)  =")):
            self.data['alat'] = pwscf_output.extract(line, "=", "a.u.", "f")  
            
          if(pwscf_output.compare(line, "unit-cell volume          =")):
            self.data['volume'] = pwscf_output.extract(line, "=", "(a.u.)^3", "f")     
            
          if(pwscf_output.compare(line, "number of atoms/cell      =")):
            self.data['nat'] = pwscf_output.extract(line, "=", "", "i")  
            
          if(pwscf_output.compare(line, "number of atomic types    =")):
            self.data['types'] = pwscf_output.extract(line, "=", "", "i")  
            
          if(pwscf_output.compare(line, "number of electrons       =")):
            str_e = pwscf_output.extract(line, "=", "", "s")
            e, eu, ed = pwscf_output.electron_string(str_e)
            self.data['electrons'] = e
            self.data['electrons_up'] = eu
            self.data['electrons_down'] = ed
          
            
          if(pwscf_output.compare(line, "number of Kohn-Sham states=")):
            self.data['ks_states'] = pwscf_output.extract(line, "=", "", "i")   
            
          if(pwscf_output.compare(line, "kinetic-energy cutoff     =")):
            self.data['ecutwfc'] = pwscf_output.extract(line, "=", "Ry", "f")  
            
          if(pwscf_output.compare(line, "charge density cutoff     =")):
            self.data['ecutrho'] = pwscf_output.extract(line, "=", "Ry", "f")   
        
          if(pwscf_output.compare(line.strip(), "crystal axes:") and pwscf_output.is_zero(self.data['crystal_in'])):            
            for j in range(3):              
              n, line, line_uc = self.next_line(n, data)
              fields = pwscf_output.extract(line, "= (", ")", "s", " ")
              self.data['crystal_in'][j,:] = fields  
              self.data['crystal_calc'][j,:] = fields  
          
          if(pwscf_output.compare(line.strip(), "crystal axes:")):            
            for j in range(3):              
              n, line, line_uc = self.next_line(n, data)
              fields = pwscf_output.extract(line, "= (", ")", "s", " ")
              self.data['crystal_calc'][j,:] = fields            
        
          
        
        
          if(pwscf_output.compare(line, "!    total energy")):
            self.data['total_energy'] = pwscf_output.extract(line, "=", "Ry", "f")
            
          if(pwscf_output.compare(line, "Total force =")):
            self.data['total_force'] = pwscf_output.extract(line, "=", "T", "f")
            
            
            
          if(pwscf_output.compare(line, "total   stress  (Ry/bohr**3)")):        
            self.data['stress_sum'] = 0.0
            for j in range(3):              
              n, line, line_uc = self.next_line(n, data)   
              fields = pwscf_output.extract(line, "", "", "f", " ", True)  
              self.data['stress'][j,0] = fields[0] 
              self.data['stress'][j,1] = fields[1] 
              self.data['stress'][j,2] = fields[2]
              self.data['stress_sum'] = self.data['stress_sum'] + abs(fields[0]) + abs(fields[1]) + abs(fields[2])
            
            
#                  "stress": numpy.zeros((3,3)),
#      "stress_sum": None,
            
          if(pwscf_output.compare(line, "density = ")):
            self.data['density_full'] = pwscf_output.extract(line, "=", "", "s")
            self.data['density'] = pwscf_output.extract(line, "=", "g/cm^3", "f")
          
          
            
          if(pwscf_output.compare(line, "PWSCF        :")):
            self.data['cpu_time'] = pwscf_output.extract(line, ":", "CPU", "s")
            
          if(pwscf_output.compare(line, "PWSCF        :")):
            self.data['wall_time'] = pwscf_output.extract(line, "CPU", "WALL", "s")
          
          if(pwscf_output.compare(line, "JOB DONE.")):
            self.data['job_done'] = True
            
                    
          if(pwscf_output.compare(line, "Exit code:")):
            self.data['error'] = True  
             

  def next_line(self, n, data):
    if(n < len(data)):
      line = data[n].strip()
      line_uc = line.upper()
      n = n + 1
      return n, line, line_uc
    else:
      n = n + 1
      return n, None, None
    
        
  def store(self, store, line, field, n=0):
    l, f = pwscf_output.read_line(line, field)  
    if(l != False):
      self.data[store] = f[n]

  #  Run as it's own program
  def run(self):
    self.reset()

    option = ""
    file_name = ""

    if(len(sys.argv) > 1 and sys.argv[1] is not None):
      option = sys.argv[1]

    if(len(sys.argv) > 2 and sys.argv[2] is not None):
      file_name = sys.argv[2]

    if(option.lower().strip() == "" or option.lower().strip() == "interactive"):
      self.menu()
      exit()
    elif(option.lower().strip() == "quiet"):
      print("Quiet")
    else:
      return 0






#################################
# READ/LOAD input file
#################################

  def load_from_file(self, file_name):
    # Init variable
    file_data = ""

    # Read it in line by line
    fh = open(file_name, "r")
    for file_row in fh:
      file_data = file_data + file_row.strip() + '\n'

    return file_data




#################################
# Get
#################################

  def get_alat(self):
    return self.data['alat']
    
  def get_volume(self):
    return self.data['volume']  
  
  def get_total_energy(self):
    return self.data['total_energy']  
    
  def get_energy_per_atom(self):
    if(self.data['total_energy'] is None):
      return None
    elif(self.data['nat'] is None):
      return None
    else:
      return (float(self.data['total_energy']) / float(self.data['nat']))    
  
  def get_total_force(self):
    return self.data['total_force']  
    
  def get_force_per_atom(self):
    return (float(self.data['total_force']) / float(self.data['nat']))  
  
  def get_density(self):
    return self.data['density']  
    
  def get_cell_parameters(self):
    cp = ['alat', 
          [str(self.data['crystal_calc'][0,0]), str(self.data['crystal_calc'][0,1]), str(self.data['crystal_calc'][0,2])], 
          [str(self.data['crystal_calc'][1,0]), str(self.data['crystal_calc'][1,1]), str(self.data['crystal_calc'][1,2])], 
          [str(self.data['crystal_calc'][2,0]), str(self.data['crystal_calc'][2,1]), str(self.data['crystal_calc'][2,2])]]
    return cp

  # Return relaxed unit vector
  def get_cell_array(self):
    return self.data['crystal_calc']

  # return alat and normalised unit vector
  def get_norm_relaxed(self):
    alat = self.data['alat']    
    cp = numpy.copy(self.data['crystal_calc'])
    f = cp[0,0]
    alat = alat * f
    cp = cp / f
    return alat, cp

  # Get stress
  def get_stress(self):
    return self.data['stress']
    
  def get_stress_sum(self):
    return self.data['stress_sum']
  

  def get_job_done(self):
    return self.data['job_done']

  def get_ok(self):
    return self.data['ok']


#################################
# Interactive
#################################

  def menu(self):
    while(True):
      choice = self.print_menu().upper()
      print(choice)
      if(choice == "X"):
        exit()
      elif(choice == "1"):
        self.i_load()
      elif(choice == "2"):
        self.i_display()


  def print_menu(self):
    pwscf_output.header("Menu")
    print("1. Load File")
    print("2. Display File")
    print("X. Exit")
    return input("Choice: ")


  def i_load(self):
    pwscf_output.header("Load Output File")
    file_name = input("Enter file name: ")
    self.load(file_name)
    print("File loaded.")
    input()


  def i_display(self):
    pwscf_output.header("Display File")
    self.output_details()
    input()


  def output_details(self):
    print("Output")
    for key in sorted(self.data.keys()):
      value = self.data[key]
      print(key, ":  ", value)
    
    
    #for key, value in self.data.items():
    #  print(key, ":  ", value)



#################################
# Static Methods
#################################


  @staticmethod
  def remove_spaces(input_string):
    return input_string.replace(" ", "")
    
    

  @staticmethod
  def extract(input_string, start=None, end=None, type=None, split=None, trim=False):
    if(start == ""):
      start = None
    if(end == ""):
      end = None
    if(trim):
      input_string = input_string.strip()
    
    # Start/End
    start_n = None
    end_n = None
      
    if(start == None and end == None):   
      start_n = 0
      end_n = len(input_string)
    elif(start == None and end != None):  
      end_l = len(end)   
      start_n = 0
      for n in range(len(input_string)):
        if(input_string[n:n+end_l] == end[0:end_l]):
          end_n = n
          break
    elif(start != None and end == None):  
      start_l = len(start)
      end_n = len(input_string)
      for n in range(len(input_string)):
        if(input_string[n:n+start_l] == start[0:start_l]):
          start_n = n + start_l
    else:  
      start_l = len(start)
      end_l = len(end)  
    
      for n in range(len(input_string)):
        if(input_string[n:n+start_l] == start[0:start_l]):
          start_n = n + start_l
        if(start_n != None and input_string[n:n+end_l] == end[0:end_l]):
          end_n = n
          break
        
    # Read
    result = input_string[start_n:end_n].strip()       

      
    # Split
    if(split != None):
      if(split == " "):
        #result = re.sub(r'\s\s+', ' ', result)
        result = pwscf_output.single_spaces(result)
      result = result.split(split)
      for i in range(len(result)):
        if(type.lower() == "f"):
          result[i] = float(result[i])
        elif(type.lower() == "i"):
          result[i] = int(result[i])
        
    else:  
      if(type.lower() == "f"):
        result = float(result)
      elif(type.lower() == "i"):
        result = int(result)
        
    # Return
    return result
      
    
  @staticmethod
  def compare(line, field):
    line = line.strip()
    line = line.upper() 
    
    field = field.strip()
    field = field.upper()
    
    f_len = len(field)
    if(len(line) >= f_len and line[0:f_len] == field[0:f_len]):
      return True
    return False
    
    
  @staticmethod
  def read_line(line, field):
    line = line.strip()
    #line = re.sub(r'\s\s+', ' ', line)    
    #line = re.sub(r'\s=\s', '=', line)
    line = pwscf_output.clean(line)
    line_uc = line.upper() 
    
    field = field.strip()
    #field = re.sub(r'\s\s+', ' ', field)
    #field = re.sub(r'\s=\s', '=', field)
    field = pwscf_output.clean(field)
    field = field.upper()
    
    f_len = len(field)
    if(len(line_uc) >= f_len and line_uc[0:f_len] == field[0:f_len]):
      output = line[f_len:].strip()
      fields = output.split(" ")
      return output, fields      
    return False, False
    
 
  @staticmethod
  def fields(input_string):
    input_string = input_string.strip()
    output_string = ""
    last = None
    for character in input_string:
      if(character != " " or (character == " " and last != " ")):
        output_string += character
    return output_string.split(" ")
    
  @staticmethod
  def check_keyword(line, keyword):
    if(line.upper()[0:len(keyword)] == keyword.upper()):
      return True
    return False


  @staticmethod
  def clear_screen():
    os.system('cls' if os.name == 'nt' else 'clear')


  @staticmethod
  def header(sub_title=""):
    pwscf_output.clear_screen()
    print("==========================================================")
    print("                    PWscf Input Editor                    ")
    print("==========================================================")
    print()
    print(sub_title)
    print()
    print()
    print()

  @staticmethod
  def process_keyword(str_in):
    str_in = str_in.lower().strip()
    str_in = pwscf_output.remove_spaces(str_in)
    id = None
    keyword = ""
    flag = 0
    for character in str_in:
      if(character == "("):
        id = ""
        flag = 1
      elif(character == ")"):
        flag = 2
      elif(flag == 0):
        keyword += character
      elif(flag == 1):
        id = id + character
    if(id != None):
      try:
        id = int(id)
      except:
        id = None
    return keyword, id  

  @staticmethod
  def add_keyword(keywords, keyword, id, value):
    if(id == None):
      added = False
      for i in range(len(keywords)):
        if(keywords[i][0] == keyword):
          added = True
          keywords[i][1] = keyword
      if(added == False):
        keywords.append([keyword, value])
    else:   
      n = None
      for i in range(len(keywords)):
        if(keywords[i][0] == keyword):
          n = i
          break
      if(n == None):    
        keywords.append([keyword,[None]])
        n = len(keywords) - 1
        
      while(len(keywords[n][1]) < id):
        keywords[n][1].append(None)

      keywords[n][1][id-1] = value  



  @staticmethod
  def make_line(key, value):
    output = ""
    if(value != None):
       if(isinstance(value, (list,))):
         for i in range(len(value)):
           if(value[i] != None):
             output += key + "(" + str(i+1) + ") = " + value[i] + ", \n"                
       else:
         output += key + " = " + value + ", \n"   
    return output    

  @staticmethod
  def coord_format(float_in):
    pad = "              "
    value = str(round(float_in, 6)).strip()
    return value
    
  @staticmethod
  def label_format(label):  
    pad = "              "
    label = label.strip()
    return label
    
    

  @staticmethod
  def is_zero(arr):
    for i in range(arr.shape[0]):
      for j in range(arr.shape[1]):
        if(arr[i, j] != 0.0):
          return False
    return True
    
    
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
      elif(last == "=" and char == " "):
        ok = False
      elif(char == " " and next == "="):
        ok = False
        
      # Add to string  
      if(ok):
        str_out += char
    return str_out    
    
    
    
  @staticmethod
  def electron_string(str_in):
    arr = str_in.split("(up:")
    e = arr[0]
    if(len(arr) == 1):
      return e.strip(), None, None
    if(len(arr)==2):
      arr_b = arr[1].split(", down:")
      eu = arr_b[0]
      arr_c = arr_b[1].split(")")
      ed = arr_c[0]
      return e.strip(), eu.strip(), ed.strip()
  
    print("TEST")
    return "","",""
    
    
    
  @staticmethod
  def single_spaces(str_in):
    str_out = ""
    last = None
    for char in str_in:
      if(char != " " or (char == " " and last != " ")):
        str_out = str_out + char
      last = char
    return str_out
    
    
