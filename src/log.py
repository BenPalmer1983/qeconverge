################################################################
#    Log
#
#
#
#
################################################################



import time
import datetime


class log:

  def __init__(self, file_name="log.txt"):
    self.data = []
    self.start_time = time.time()
    self.now = datetime.datetime.now()
    self.time_now = self.get_now()
    self.file_name = file_name

    fh = open(self.file_name, 'w')
    fh.write("####################################################################\n")  
    fh.write("                              LOG FILE                              \n")    
    fh.write("####################################################################\n")   
    fh.write("\n")    
    fh.write(self.time_now + "\n")     
    fh.write("\n")    
    fh.write("\n")   
    fh.close()

  
  def set_time(self, time_in):
    self.start_time = time_in

  def add(self, log_line):
    if(type(log_line) == list):
      for line in log_line:
        self.add(line)
    else:
      time_str = str(round(time.time() - self.start_time, 4))
      p = 11 - len(time_str)
      pad = "            "
      line = time_str + ":" + pad[0:p] + log_line
      fh = open(self.file_name, 'a')
      fh.write(line + "\n")    
      fh.close()

  def log(self, line):
    self.add(line)

  def output(self, file_name="log.txt"):
    pass



  def get_now(self):
    now = datetime.datetime.now()
    hour = str(now.hour)
    if(now.hour < 10):
      hour = "0" + hour
    minute = str(now.minute)
    if(now.minute < 10):
      minute = "0" + minute
    day = str(now.day)
    if(now.day < 10):
      day = "0" + day
    month = str(now.month)
    if(now.month < 10):
      month = "0" + month
    year = str(now.year)
    return hour + ":" + minute + "   " + day + "/" + month + "/" + year
    
