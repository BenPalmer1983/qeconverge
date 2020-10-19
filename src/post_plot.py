import numpy
from units import units
from globals import globals
from plot_output import plot_output
from plots import plots
from matplotlib.colors import LinearSegmentedColormap
from square import square

###########################################
#  CLASS conv_ecut2d
###########################################
class post_plot:

  @staticmethod
  def run():
  
    # ECUT LIMITS
    
    fh = open(globals.dirs['csv'] + "/ecut2d_xy_lims.csv", 'r')
    data = []
    for line in fh:
      line = line.strip()
      if(line != ""):
        data.append(float(line.strip()))
  
    wfc_min = data[0]
    wfc_max = data[1]
    rho_min = data[2]
    rho_max = data[3]
  
  
  
    # ECUT LIMITS
    
    fh = open(globals.dirs['csv'] + "/kpoints_xy_lims.csv", 'r')
    data = []
    for line in fh:
      line = line.strip()
      if(line != ""):
        data.append(float(line.strip()))
  
    kpoints_min = data[0]
    kpoints_max = data[1]
    smear_min = data[2]
    smear_max = data[-1]
  
    Z = std.csv_to_array(globals.dirs['csv'] + "/ecut2d_energy_ry.csv")
    Z_ecut2d_e_in = numpy.copy(Z)
    Z_ecut2d_e = square.square_data(Z, 200, 200, 2)
    Z_ecut2d_adj_e = Z_ecut2d_e - numpy.amin(Z_ecut2d_e)
    Z_ecut2d_rhoconv_e = square.square_data(post_plot.conv_w(Z - numpy.amin(Z)), 200, 200, 2)
    Z_ecut2d_wfcconv_e = square.square_data(post_plot.conv_h(Z - numpy.amin(Z)), 200, 200, 2)

    
    Z = std.csv_to_array(globals.dirs['csv'] + "/ecut2d_force_ry.csv")
    Z_ecut2d_f_in = numpy.copy(Z)
    Z_ecut2d_f = square.square_data(Z, 200, 200, 2)
    Z_ecut2d_adj_f = Z_ecut2d_f - numpy.amin(Z_ecut2d_f)
    Z_ecut2d_rhoconv_f = square.square_data(post_plot.conv_w(Z - numpy.amin(Z)), 200, 200, 2)
    Z_ecut2d_wfcconv_f = square.square_data(post_plot.conv_h(Z - numpy.amin(Z)), 200, 200, 2)
  
    Z = std.csv_to_array(globals.dirs['csv'] + "/kpoints_energy_ry.csv")
    Z_kpoints2d_e_in = numpy.copy(Z)
    #Z = numpy.flip(Z, 1)
    Z_kpoints2d_e = square.square_data(Z, 200, 200, 2)
    Z_kpoints2d_adj_e = Z_kpoints2d_e - numpy.amin(Z_kpoints2d_e)
    Z_kpoints2d_kpointconv_e = square.square_data(post_plot.conv_w(Z - numpy.amin(Z)), 200, 200, 2)
    Z_kpoints2d_degaussconv_e = square.square_data(post_plot.conv_h(Z - numpy.amin(Z)), 200, 200, 2)
    
    Z = std.csv_to_array(globals.dirs['csv'] + "/kpoints_force_ry.csv")
    #Z = numpy.flip(Z, 1)
    Z_kpoints2d_f_in = numpy.copy(Z)
    Z_kpoints2d_f = square.square_data(Z, 200, 200, 2)
    Z_kpoints2d_adj_f = Z_kpoints2d_f - numpy.amin(Z_kpoints2d_f)
    Z_kpoints2d_kpointconv_f = square.square_data(post_plot.conv_w(Z - numpy.amin(Z)), 200, 200, 2)
    Z_kpoints2d_degaussconv_f = square.square_data(post_plot.conv_h(Z - numpy.amin(Z)), 200, 200, 2)
  
    
    # Plot original uninterpolated data
    x_rho = numpy.linspace(rho_min, rho_max, len(Z_ecut2d_e_in[0,:]))
    y_wfc = numpy.linspace(wfc_min, wfc_max, len(Z_ecut2d_e_in))
    
    x_kpoints = numpy.linspace(kpoints_min, kpoints_max, len(Z_kpoints2d_e_in[0,:]))
    y_smear = numpy.linspace(smear_min, smear_max, len(Z_kpoints2d_e_in))
    
    print(len(Z_ecut2d_e_in[0,:]), len(Z_ecut2d_e_in))
    print(len(Z_kpoints2d_e_in[0,:]), len(Z_kpoints2d_e_in))
    
 
    sa = {'title': 'Energy (Ry) Ecutwfc and Ecutrho', 'x_axis': 'Ecutrho/Ecutwfc', 'y_axis': 'Ecutwfc/Ry', 
         'x': x_rho, 'y': y_wfc, 'Z': Z_ecut2d_e_in, 'cmap': 'viridis'}     
    sb = {'title': 'Force (Ry/Bohr) Ecutwfc and Ecutrho', 'x_axis': 'Ecutrho/Ecutwfc', 'y_axis': 'Ecutwfc/Ry', 
         'x': x_rho, 'y': y_wfc, 'Z': Z_ecut2d_f_in, 'cmap': 'viridis'}     
    sc = {'title': 'Energy (Ry) Kpoint and Degauss', 'x_axis': 'K-points', 'y_axis': 'Degauss/Ry', 
         'x': x_kpoints, 'y': y_smear, 'Z': Z_kpoints2d_e_in, 'cmap': 'viridis'}          
    sd = {'title': 'Force (Ry/Bohr) Kpoint and Degauss', 'x_axis': 'K-points', 'y_axis': 'Degauss/Ry', 
         'x': x_kpoints, 'y': y_smear, 'Z': Z_kpoints2d_f_in, 'cmap': 'viridis'}         
    plots.clear('original_summary_colour', 'Ecutwfc and Ecutrho, Kpoint and Degauss', False)
    plots.add_colour_plot(sa)
    plots.add_colour_plot(sb)
    plots.add_colour_plot(sc)
    plots.add_colour_plot(sd)
    plots.make(globals.dirs['plots_eps'], globals.dirs['plots_svg'])
  
  
  
 
    # AXES
 
    x_rho = numpy.linspace(rho_min, rho_max, 200)
    y_wfc = numpy.linspace(wfc_min, wfc_max, 200)
    
    x_kpoints = numpy.linspace(kpoints_min, kpoints_max, 200)
    y_smear = numpy.linspace(smear_min, smear_max, 200)
    

    sa = {'title': 'Energy (Ry) Ecutwfc and Ecutrho', 'x_axis': 'Ecutrho/Ecutwfc', 'y_axis': 'Ecutwfc/Ry', 
         'x': x_rho, 'y': y_wfc, 'Z': Z_ecut2d_adj_e, 'cmap': 'viridis'}     
    sb = {'title': 'Force (Ry/Bohr) Ecutwfc and Ecutrho', 'x_axis': 'Ecutrho/Ecutwfc', 'y_axis': 'Ecutwfc/Ry', 
         'x': x_rho, 'y': y_wfc, 'Z': Z_ecut2d_adj_f, 'cmap': 'viridis'}     
    sc = {'title': 'Energy (Ry) Kpoint and Degauss', 'x_axis': 'K-points', 'y_axis': 'Degauss/Ry', 
         'x': x_kpoints, 'y': y_smear, 'Z': Z_kpoints2d_adj_e, 'cmap': 'viridis'}          
    sd = {'title': 'Force (Ry/Bohr) Kpoint and Degauss', 'x_axis': 'K-points', 'y_axis': 'Degauss/Ry', 
         'x': x_kpoints, 'y': y_smear, 'Z': Z_kpoints2d_adj_f, 'cmap': 'viridis'}         
    plots.clear('interpolated_summary_colour', 'Ecutwfc and Ecutrho, Kpoint and Degauss', False)
    plots.add_colour_plot(sa)
    plots.add_colour_plot(sb)
    plots.add_colour_plot(sc)
    plots.add_colour_plot(sd)
    plots.make(globals.dirs['plots_eps'], globals.dirs['plots_svg'])
 
    # ECUT ENERGY
    
    s = {'title': 'Energy (Ry) vs Ecutwfc and Ecutrho', 'x_axis': 'Ecutrho/Ecutwfc', 'y_axis': 'Ecutwfc/Ry', 
         'x': x_rho, 'y': y_wfc, 'Z': Z_ecut2d_e, 'cmap': 'Greys'}    
    plots.clear('energy_grey', '', True)
    plots.add_colour_plot(s)
    plots.make(globals.dirs['plots_eps'], globals.dirs['plots_svg'])
    
    s = {'title': 'Energy (Ry) vs Ecutwfc and Ecutrho', 'x_axis': 'Ecutrho/E, h, wcutwfc', 'y_axis': 'Ecutwfc/Ry', 
         'x': x_rho, 'y': y_wfc, 'Z': Z_ecut2d_e, 'cmap': 'viridis'}    
    plots.clear('energy_colour', '', True)
    plots.add_colour_plot(s)
    plots.make(globals.dirs['plots_eps'], globals.dirs['plots_svg'])
    
    s = {'title': 'Energy (Ry) (Adjusted) vs Ecutwfc and Ecutrho', 'x_axis': 'Ecutrho/Ecutwfc', 'y_axis': 'Ecutwfc/Ry', 
         'x': x_rho, 'y': y_wfc, 'Z': Z_ecut2d_adj_e, 'cmap': 'Greys'}    
    plots.clear('energy_adjusted_grey', '', True)
    plots.add_colour_plot(s)
    plots.make(globals.dirs['plots_eps'], globals.dirs['plots_svg'])
    
    s = {'title': 'Energy (Ry) (Adjusted) vs Ecutwfc and Ecutrho', 'x_axis': 'Ecutrho/Ecutwfc', 'y_axis': 'Ecutwfc/Ry', 
         'x': x_rho, 'y': y_wfc, 'Z': Z_ecut2d_adj_e, 'cmap': 'viridis'}    
    plots.clear('energy_adjusted_colour', '', True)
    plots.add_colour_plot(s)
    plots.make(globals.dirs['plots_eps'], globals.dirs['plots_svg'])
    
    s = {'title': 'Energy (Ry) Converge wrt Rho', 'x_axis': 'Ecutrho/Ecutwfc', 'y_axis': 'Ecutwfc/Ry', 
         'x': x_rho, 'y': y_wfc, 'Z': Z_ecut2d_rhoconv_e, 'cmap': 'Greys'}    
    plots.clear('energy_conv_wrt_rho_grey', '', True)
    plots.add_colour_plot(s)
    plots.make(globals.dirs['plots_eps'], globals.dirs['plots_svg'])
    
    s = {'title': 'Energy (Ry) Converge wrt Rho', 'x_axis': 'Ecutrho/Ecutwfc', 'y_axis': 'Ecutwfc/Ry', 
         'x': x_rho, 'y': y_wfc, 'Z': Z_ecut2d_rhoconv_e, 'cmap': 'viridis'}    
    plots.clear('energy_conv_wrt_rho_colour', '', True)
    plots.add_colour_plot(s)
    plots.make(globals.dirs['plots_eps'], globals.dirs['plots_svg'])
  
    
    s = {'title': 'Energy (Ry) Converge wrt Wfc', 'x_axis': 'Ecutrho/Ecutwfc', 'y_axis': 'Ecutwfc/Ry', 
         'x': x_rho, 'y': y_wfc, 'Z': Z_ecut2d_wfcconv_e, 'cmap': 'Greys'}    
    plots.clear('energy_conv_wrt_wfc_grey', '', True)
    plots.add_colour_plot(s)
    plots.make(globals.dirs['plots_eps'], globals.dirs['plots_svg'])
    
    s = {'title': 'Energy (Ry) Converge wrt Wfc', 'x_axis': 'Ecutrho/Ecutwfc', 'y_axis': 'Ecutwfc/Ry', 
         'x': x_rho, 'y': y_wfc, 'Z': Z_ecut2d_wfcconv_e, 'cmap': 'viridis'}    
    plots.clear('energy_conv_wrt_wfc_colour', '', True)
    plots.add_colour_plot(s)
    plots.make(globals.dirs['plots_eps'], globals.dirs['plots_svg'])
  
  
    
  
    # ECUT FORCE
    
    s = {'title': 'Force (Ry/Bohr) vs Ecutwfc and Ecutrho', 'x_axis': 'Ecutrho/Ecutwfc', 'y_axis': 'Ecutwfc/Ry', 
         'x': x_rho, 'y': y_wfc, 'Z': Z_ecut2d_f, 'cmap': 'Greys'}    
    plots.clear('force_grey', '', True)
    plots.add_colour_plot(s)
    plots.make(globals.dirs['plots_eps'], globals.dirs['plots_svg'])
    
    s = {'title': 'Force (Ry/Bohr) vs Ecutwfc and Ecutrho', 'x_axis': 'Ecutrho/E, h, wcutwfc', 'y_axis': 'Ecutwfc/Ry', 
         'x': x_rho, 'y': y_wfc, 'Z': Z_ecut2d_f, 'cmap': 'viridis'}    
    plots.clear('force_colour', '', True)
    plots.add_colour_plot(s)
    plots.make(globals.dirs['plots_eps'], globals.dirs['plots_svg'])
    
    s = {'title': 'Force (Ry/Bohr) (Adjusted) vs Ecutwfc and Ecutrho', 'x_axis': 'Ecutrho/Ecutwfc', 'y_axis': 'Ecutwfc/Ry', 
         'x': x_rho, 'y': y_wfc, 'Z': Z_ecut2d_adj_f, 'cmap': 'Greys'}    
    plots.clear('force_adjusted_grey', '', True)
    plots.add_colour_plot(s)
    plots.make(globals.dirs['plots_eps'], globals.dirs['plots_svg'])
    
    s = {'title': 'Force (Ry/Bohr) (Adjusted) vs Ecutwfc and Ecutrho', 'x_axis': 'Ecutrho/Ecutwfc', 'y_axis': 'Ecutwfc/Ry', 
         'x': x_rho, 'y': y_wfc, 'Z': Z_ecut2d_adj_f, 'cmap': 'viridis'}    
    plots.clear('force_adjusted_colour', '', True)
    plots.add_colour_plot(s)
    plots.make(globals.dirs['plots_eps'], globals.dirs['plots_svg'])
    
    s = {'title': 'Force (Ry/Bohr) Converge wrt Rho', 'x_axis': 'Ecutrho/Ecutwfc', 'y_axis': 'Ecutwfc/Ry', 
         'x': x_rho, 'y': y_wfc, 'Z': Z_ecut2d_rhoconv_f, 'cmap': 'Greys'}    
    plots.clear('force_conv_wrt_rho_grey', '', True)
    plots.add_colour_plot(s)
    plots.make(globals.dirs['plots_eps'], globals.dirs['plots_svg'])
    
    s = {'title': 'Force (Ry/Bohr) Converge wrt Rho', 'x_axis': 'Ecutrho/Ecutwfc', 'y_axis': 'Ecutwfc/Ry', 
         'x': x_rho, 'y': y_wfc, 'Z': Z_ecut2d_rhoconv_f, 'cmap': 'viridis'}    
    plots.clear('force_conv_wrt_rho_colour', '', True)
    plots.add_colour_plot(s)
    plots.make(globals.dirs['plots_eps'], globals.dirs['plots_svg'])
  
    
    s = {'title': 'Force (Ry/Bohr) Converge wrt Wfc', 'x_axis': 'Ecutrho/Ecutwfc', 'y_axis': 'Ecutwfc/Ry', 
         'x': x_rho, 'y': y_wfc, 'Z': Z_ecut2d_wfcconv_f, 'cmap': 'Greys'}    
    plots.clear('force_conv_wrt_wfc_grey', '', True)
    plots.add_colour_plot(s)
    plots.make(globals.dirs['plots_eps'], globals.dirs['plots_svg'])
    
    s = {'title': 'Force (Ry/Bohr) Converge wrt Wfc', 'x_axis': 'Ecutrho/Ecutwfc', 'y_axis': 'Ecutwfc/Ry', 
         'x': x_rho, 'y': y_wfc, 'Z': Z_ecut2d_wfcconv_f, 'cmap': 'viridis'}    
    plots.clear('force_conv_wrt_wfc_colour', '', True)
    plots.add_colour_plot(s)
    plots.make(globals.dirs['plots_eps'], globals.dirs['plots_svg'])
  
  
    
    sa = {'title': 'Energy (Ry) Converge wrt Wfc', 'x_axis': 'Ecutrho/Ecutwfc', 'y_axis': 'Ecutwfc/Ry', 
         'x': x_rho, 'y': y_wfc, 'Z': Z_ecut2d_wfcconv_e, 'cmap': 'Greys'} 
    sb = {'title': 'Energy (Ry) Converge wrt Rho', 'x_axis': 'Ecutrho/Ecutwfc', 'y_axis': 'Ecutwfc/Ry', 
         'x': x_rho, 'y': y_wfc, 'Z': Z_ecut2d_rhoconv_e, 'cmap': 'Greys'}  
    sc = {'title': 'Force (Ry/Bohr) Converge wrt Wfc', 'x_axis': 'Ecutrho/Ecutwfc', 'y_axis': 'Ecutwfc/Ry', 
         'x': x_rho, 'y': y_wfc, 'Z': Z_ecut2d_wfcconv_f, 'cmap': 'Greys'}  
    sd = {'title': 'Force (Ry/Bohr) Converge wrt Rho', 'x_axis': 'Ecutrho/Ecutwfc', 'y_axis': 'Ecutwfc/Ry', 
         'x': x_rho, 'y': y_wfc, 'Z': Z_ecut2d_rhoconv_f, 'cmap': 'Greys'} 
        
    plots.clear('summary_grey', 'Convergence vs Ecutwfc and Ecutrho', True)
    plots.add_colour_plot(sa)
    plots.add_colour_plot(sb)
    plots.add_colour_plot(sc)
    plots.add_colour_plot(sd)
    plots.make(globals.dirs['plots_eps'], globals.dirs['plots_svg'])
  
  
    
    sa = {'title': 'Energy (Ry) Converge wrt Wfc', 'x_axis': 'Ecutrho/Ecutwfc', 'y_axis': 'Ecutwfc/Ry', 
         'x': x_rho, 'y': y_wfc, 'Z': Z_ecut2d_wfcconv_e, 'cmap': 'viridis'}     
    sb = {'title': 'Energy (Ry) Converge wrt Rho', 'x_axis': 'Ecutrho/Ecutwfc', 'y_axis': 'Ecutwfc/Ry', 
         'x': x_rho, 'y': y_wfc, 'Z': Z_ecut2d_rhoconv_e, 'cmap': 'viridis'}     
    sc = {'title': 'Force (Ry/Bohr) Converge wrt Wfc', 'x_axis': 'Ecutrho/Ecutwfc', 'y_axis': 'Ecutwfc/Ry', 
         'x': x_rho, 'y': y_wfc, 'Z': Z_ecut2d_wfcconv_f, 'cmap': 'viridis'}          
    sd = {'title': 'Force (Ry/Bohr) Converge wrt Rho', 'x_axis': 'Ecutrho/Ecutwfc', 'y_axis': 'Ecutwfc/Ry', 
         'x': x_rho, 'y': y_wfc, 'Z': Z_ecut2d_rhoconv_f, 'cmap': 'viridis'}         
    plots.clear('summary_colour', 'Convergence vs Ecutwfc and Ecutrho', True)
    plots.add_colour_plot(sa)
    plots.add_colour_plot(sb)
    plots.add_colour_plot(sc)
    plots.add_colour_plot(sd)
    plots.make(globals.dirs['plots_eps'], globals.dirs['plots_svg'])
    
  
  
    
  
    # KPOINTS ENERGY
    
    s = {'title': 'Energy (Ry) vs Ecutwfc and Ecutrho', 'x_axis': 'K-Points', 'y_axis': 'Degauss/Ry', 
         'x': x_kpoints, 'y': y_smear, 'Z': Z_kpoints2d_e, 'cmap': 'Greys'}    
    plots.clear('kpoints_energy_grey', '', True)
    plots.add_colour_plot(s)
    plots.make(globals.dirs['plots_eps'], globals.dirs['plots_svg'])
     
    s = {'title': 'Energy (Ry) vs Ecutwfc and Ecutrho', 'x_axis': 'K-Points', 'y_axis': 'Degauss/Ry', 
         'x': x_kpoints, 'y': y_smear, 'Z': Z_kpoints2d_e, 'cmap': 'viridis'}    
    plots.clear('kpoints_energy_colour', '', True)
    plots.add_colour_plot(s)
    plots.make(globals.dirs['plots_eps'], globals.dirs['plots_svg'])
     
    s = {'title': 'Energy (Ry) vs Ecutwfc and Ecutrho', 'x_axis': 'K-Points', 'y_axis': 'Degauss/Ry', 
         'x': x_kpoints, 'y': y_smear, 'Z': Z_kpoints2d_adj_e, 'cmap': 'Greys'}    
    plots.clear('kpoints_energy_adjusted_grey', '', True)
    plots.add_colour_plot(s)
    plots.make(globals.dirs['plots_eps'], globals.dirs['plots_svg'])
     
    s = {'title': 'Energy (Ry) vs Ecutwfc and Ecutrho', 'x_axis': 'K-Points', 'y_axis': 'Degauss/Ry', 
         'x': x_kpoints, 'y': y_smear, 'Z': Z_kpoints2d_adj_e, 'cmap': 'viridis'}    
    plots.clear('kpoints_energy_adjusted_colour', '', True)
    plots.add_colour_plot(s)
    plots.make(globals.dirs['plots_eps'], globals.dirs['plots_svg'])
    
    
    
    
    s = {'title': 'Energy (Ry) Converge wrt Kpoint', 'x_axis': 'K-Points', 'y_axis': 'Degauss/Ry', 
         'x': x_kpoints, 'y': y_smear, 'Z': Z_kpoints2d_kpointconv_e, 'cmap': 'Greys'}    
    plots.clear('kpoints_energy_conv_wrt_kpoint_grey', '', True)
    plots.add_colour_plot(s)
    plots.make(globals.dirs['plots_eps'], globals.dirs['plots_svg'])
    
    s = {'title': 'Energy (Ry) Converge wrt Kpoint', 'x_axis': 'K-Points', 'y_axis': 'Degauss/Ry', 
         'x': x_kpoints, 'y': y_smear, 'Z': Z_kpoints2d_kpointconv_e, 'cmap': 'viridis'}    
    plots.clear('kpoints_energy_conv_wrt_kpoint_colour', '', True)
    plots.add_colour_plot(s)
    plots.make(globals.dirs['plots_eps'], globals.dirs['plots_svg'])
  
    
    s = {'title': 'Energy (Ry) Converge wrt Degauss', 'x_axis': 'K-Points', 'y_axis': 'Degauss/Ry', 
         'x': x_kpoints, 'y': y_smear, 'Z': Z_kpoints2d_degaussconv_e, 'cmap': 'Greys'}    
    plots.clear('kpoints_energy_conv_wrt_degauss_grey', '', True)
    plots.add_colour_plot(s)
    plots.make(globals.dirs['plots_eps'], globals.dirs['plots_svg'])
    
    s = {'title': 'Energy (Ry) Converge wrt Degauss', 'x_axis': 'K-Points', 'y_axis': 'Degauss/Ry', 
         'x': x_kpoints, 'y': y_smear, 'Z': Z_kpoints2d_degaussconv_e, 'cmap': 'viridis'}    
    plots.clear('kpoints_energy_conv_wrt_degauss_colour', '', True)
    plots.add_colour_plot(s)
    plots.make(globals.dirs['plots_eps'], globals.dirs['plots_svg'])
    
    
    
    
    
    
  
    # KPOINTS FORCES
    
    s = {'title': 'Force (Ry/Bohr) vs Ecutwfc and Ecutrho', 'x_axis': 'K-Points', 'y_axis': 'Degauss/Ry', 
         'x': x_kpoints, 'y': y_smear, 'Z': Z_kpoints2d_e, 'cmap': 'Greys'}    
    plots.clear('kpoints_force_grey', '', True)
    plots.add_colour_plot(s)
    plots.make(globals.dirs['plots_eps'], globals.dirs['plots_svg'])
     
    s = {'title': 'Force (Ry/Bohr) vs Ecutwfc and Ecutrho', 'x_axis': 'K-Points', 'y_axis': 'Degauss/Ry', 
         'x': x_kpoints, 'y': y_smear, 'Z': Z_kpoints2d_e, 'cmap': 'viridis'}    
    plots.clear('kpoints_force_colour', '', True)
    plots.add_colour_plot(s)
    plots.make(globals.dirs['plots_eps'], globals.dirs['plots_svg'])
     
    s = {'title': 'Force (Ry/Bohr) vs Ecutwfc and Ecutrho', 'x_axis': 'K-Points', 'y_axis': 'Degauss/Ry', 
         'x': x_kpoints, 'y': y_smear, 'Z': Z_kpoints2d_adj_e, 'cmap': 'Greys'}    
    plots.clear('kpoints_force_adjusted_grey', '', True)
    plots.add_colour_plot(s)
    plots.make(globals.dirs['plots_eps'], globals.dirs['plots_svg'])
     
    s = {'title': 'Force (Ry/Bohr) vs Ecutwfc and Ecutrho', 'x_axis': 'K-Points', 'y_axis': 'Degauss/Ry', 
         'x': x_kpoints, 'y': y_smear, 'Z': Z_kpoints2d_adj_e, 'cmap': 'viridis'}    
    plots.clear('kpoints_force_adjusted_colour', '', True)
    plots.add_colour_plot(s)
    plots.make(globals.dirs['plots_eps'], globals.dirs['plots_svg'])
    
    
    
    s = {'title': 'Force (Ry/Bohr) Converge wrt Kpoint', 'x_axis': 'K-Points', 'y_axis': 'Degauss/Ry', 
         'x': x_kpoints, 'y': y_smear, 'Z': Z_kpoints2d_kpointconv_f, 'cmap': 'Greys'}    
    plots.clear('kpoints_force_conv_wrt_kpoint_grey', '', True)
    plots.add_colour_plot(s)
    plots.make(globals.dirs['plots_eps'], globals.dirs['plots_svg'])
    
    s = {'title': 'Force (Ry/Bohr) Converge wrt Kpoint', 'x_axis': 'K-Points', 'y_axis': 'Degauss/Ry', 
         'x': x_kpoints, 'y': y_smear, 'Z': Z_kpoints2d_kpointconv_f, 'cmap': 'viridis'}    
    plots.clear('kpoints_force_conv_wrt_kpoint_colour', '', True)
    plots.add_colour_plot(s)
    plots.make(globals.dirs['plots_eps'], globals.dirs['plots_svg'])
      
    s = {'title': 'Force (Ry/Bohr) Converge wrt Degauss', 'x_axis': 'K-Points', 'y_axis': 'Degauss/Ry', 
         'x': x_kpoints, 'y': y_smear, 'Z': Z_kpoints2d_degaussconv_f, 'cmap': 'Greys'}    
    plots.clear('kpoints_force_conv_wrt_degauss_grey', '', True)
    plots.add_colour_plot(s)
    plots.make(globals.dirs['plots_eps'], globals.dirs['plots_svg'])
    
    s = {'title': 'Force (Ry/Bohr) Converge wrt Degauss', 'x_axis': 'K-Points', 'y_axis': 'Degauss/Ry', 
         'x': x_kpoints, 'y': y_smear, 'Z': Z_kpoints2d_degaussconv_f, 'cmap': 'viridis'}    
    plots.clear('kpoints_force_conv_wrt_degauss_colour', '', True)
    plots.add_colour_plot(s)
    plots.make(globals.dirs['plots_eps'], globals.dirs['plots_svg'])
    
    """
    s = {'title': 'Force (Ry/Bohr) vs Ecutwfc and Ecutrho', 'x_axis': 'Ecutrho/E, h, wcutwfc', 'y_axis': 'Ecutwfc/Ry', 
         'x': x_kpoints, 'y': y_smear, 'Z': Z_ecut2d_f, 'cmap': 'viridis'}    
    plots.clear('energy_colour', '', True)
    plots.add_colour_plot(s)
    plots.make(globals.dirs['plots_eps'], globals.dirs['plots_svg'])
    """
    
     
    
    sa = {'title': 'Energy (Ry) Converge wrt Kpoint', 'x_axis': 'Kpoint', 'y_axis': 'Degauss/Ry', 
         'x': x_kpoints, 'y': y_smear, 'Z': Z_kpoints2d_kpointconv_e, 'cmap': 'Greys'}     
    sb = {'title': 'Energy (Ry) Converge wrt Degauss', 'x_axis': 'Kpoint', 'y_axis': 'Degauss/Ry', 
         'x': x_kpoints, 'y': y_smear, 'Z': Z_kpoints2d_degaussconv_e, 'cmap': 'Greys'}     
    sc = {'title': 'Force (Ry/Bohr) Converge wrt Kpoint', 'x_axis': 'Kpoint', 'y_axis': 'Degauss/Ry', 
         'x': x_kpoints, 'y': y_smear, 'Z': Z_kpoints2d_kpointconv_f, 'cmap': 'Greys'}          
    sd = {'title': 'Force (Ry/Bohr) Converge wrt Degauss', 'x_axis': 'Kpoint', 'y_axis': 'Degauss/Ry', 
         'x': x_kpoints, 'y': y_smear, 'Z': Z_kpoints2d_degaussconv_f, 'cmap': 'Greys'}         
    plots.clear('kpoint_summary_grey', 'Convergence vs Kpoint and Degauss', True)
    plots.add_colour_plot(sa)
    plots.add_colour_plot(sb)
    plots.add_colour_plot(sc)
    plots.add_colour_plot(sd)
    plots.make(globals.dirs['plots_eps'], globals.dirs['plots_svg'])       
    
    sa = {'title': 'Energy (Ry) Converge wrt Kpoint', 'x_axis': 'Kpoint', 'y_axis': 'Degauss/Ry', 
         'x': x_kpoints, 'y': y_smear, 'Z': Z_kpoints2d_kpointconv_e, 'cmap': 'viridis'}     
    sb = {'title': 'Energy (Ry) Converge wrt Degauss', 'x_axis': 'Kpoint', 'y_axis': 'Degauss/Ry', 
         'x': x_kpoints, 'y': y_smear, 'Z': Z_kpoints2d_degaussconv_e, 'cmap': 'viridis'}     
    sc = {'title': 'Force (Ry/Bohr) Converge wrt Kpoint', 'x_axis': 'Kpoint', 'y_axis': 'Degauss/Ry', 
         'x': x_kpoints, 'y': y_smear, 'Z': Z_kpoints2d_kpointconv_f, 'cmap': 'viridis'}          
    sd = {'title': 'Force (Ry/Bohr) Converge wrt Degauss', 'x_axis': 'Kpoint', 'y_axis': 'Degauss/Ry', 
         'x': x_kpoints, 'y': y_smear, 'Z': Z_kpoints2d_degaussconv_f, 'cmap': 'viridis'}  
    plots.clear('kpoint_summary_colour', 'Convergence vs Kpoint and Degauss', True)
    plots.add_colour_plot(sa)
    plots.add_colour_plot(sb)
    plots.add_colour_plot(sc)
    plots.add_colour_plot(sd)
    plots.make(globals.dirs['plots_eps'], globals.dirs['plots_svg'])
    
    
    
    
  @staticmethod
  def conv_w(data_in):
    h = len(data_in)
    w = len(data_in[0,:])
    out = numpy.zeros((h, w-1),)
    for i in range(w-1):
      out[:,i] = abs(data_in[:,i+1] - data_in[:,i])
    return out
    
  @staticmethod
  def conv_h(data_in):
    h = len(data_in)
    w = len(data_in[0,:])
    out = numpy.zeros((h-1, w),)
    for i in range(h-1):
      out[i,:] = abs(data_in[i+1,:] - data_in[i,:])
    return out
  
  
  
  
  
  
  
  
  