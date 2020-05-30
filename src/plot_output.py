import numpy
from units import units
from globals import globals

###########################################
#  CLASS conv_ecut2d
###########################################
class plot_output:

  @staticmethod
  def plot(plt, file_name):    
    plt.savefig(globals.dirs['plots_svg'] + '/' + file_name + '.svg', format='svg')
    plt.savefig(globals.dirs['plots_eps'] + '/' + file_name + '.eps', format='eps')
  
    plt.rcParams["figure.figsize"] = [4,3]
    plt.savefig(globals.dirs['plots_eps_small'] + '/' + file_name + '.eps', format='eps')
  
  