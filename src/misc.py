import os

#################################
# Misc functions
#################################

class misc:

  @staticmethod
  def makedir(dir):
    try:
      os.mkdir(dir)
    except:
      pass