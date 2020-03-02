
class read_config:
  
  @staticmethod
  def read_file(file_path):
    input = {}
    fh = open(file_path, 'r')
    for line in fh:
      input = read_config.read_line(input, line)
    fh.close()
    return input
    
  @staticmethod
  def read_line(input, line):
    line = line.strip()
    if(len(line) > 0 and line[0] != "#"):
      cmd, data = read_config.get_fields_dict(line)
      if(cmd not in input.keys()):
        input[cmd] = []
      input[cmd].append(data)
    return input
        
  @staticmethod  
  def get_fields_dict(line): 
    fields = read_config.split_by(line, ' ')    
    dict = {}    
    dict['COMMAND'] = fields[0]
    
    for field in fields:
      f = field.split("=")
      if(len(f) == 2):
        fb = read_config.split_by(f[1], ',')
        if(len(fb) == 1):
          dict[f[0].lower()] = [f[1]]
        elif(len(fb) > 1):
          dict[f[0].lower()] = fb
    
    return fields[0], dict
    
  @staticmethod        
  def get(data, c1, c2, n1=0, n2=0):
    if(c1 not in data.keys()):
      return None
    c1_list = data[c1]
    if(n1>=len(c1_list)):
      return ''
    if(c2 not in data[c1][n1].keys()):
      return None
    c2_list = data[c1][n1][c2]
    if(n2>=len(c2_list)):
      return ''
    return data[c1][n1][c2][n2]

  @staticmethod        
  def get_list(data, c1, c2, n1=0):
    if(c1 not in data.keys()):
      return None
    c1_list = data[c1]
    if(n1>=len(c1_list)):
      return ''
    if(c2 not in data[c1][n1].keys()):
      return None
    return data[c1][n1][c2]
        
        
    
  @staticmethod  
  def split_by(line, sep=' ', ignore_double_sep=True):
    last_char = None
    in_quotes = 0
    fields = []
    temp_line = ""
    
    for char in line:
      if(char == "'" and in_quotes == 0 and last_char != "\\"):
        in_quotes = 1
      elif(char == "'" and in_quotes == 1 and last_char != "\\"):
        in_quotes = 0
      elif(char == '"' and in_quotes == 0 and last_char != "\\"):
        in_quotes = 2
      elif(char == '"' and in_quotes == 2 and last_char != "\\"):
        in_quotes = 0
      elif(in_quotes > 0):
        temp_line = temp_line + char
      elif(in_quotes == 0 and char != sep):
        temp_line = temp_line + char
      elif(char == sep and last_char == sep and ignore_double_sep):
        pass
      elif(char == sep):
        fields.append(temp_line)
        temp_line = "" 
    if(temp_line != ""):
      fields.append(temp_line)
    
    return fields
