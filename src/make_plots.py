
import numpy
from units import units
from globals import globals
from plot_output import plot_output
from plots import plots
from matplotlib.colors import LinearSegmentedColormap
from scipy import interpolate

###########################################
#  CLASS conv_ecut2d
###########################################
class make_plots:

  @staticmethod
  def run():
  
    
    
    #############################################################################################
    #############################################################################################    
    # ECUT
    #############################################################################################
    #############################################################################################

    
    
    #############################################################################################
    #############################################################################################    
    # ECUT 2D
    #############################################################################################
    #############################################################################################


    # ENERGY
    #####################
    
    data_h = globals.ecut2d['data_h']
    data_w = globals.ecut2d['data_w']
    
    x = numpy.linspace(globals.ecut2d['rho_min_pw'], globals.ecut2d['rho_max_pw'], data_h)
    y = numpy.linspace(globals.ecut2d['wfc_min_pw'], globals.ecut2d['wfc_max_pw'], data_w)
    Z = globals.ecut2d['energy_ry'][0:data_h,0:data_w]    
    
    s = {'title': 'Energy (Ry) vs Ecutwfc and Ecutrho', 'x_axis': 'Ecutrho/Ecutwfc', 'y_axis': 'Ecutwfc/Ry', 'x': x, 'y': y, 'Z': Z, 'cmap': 'Greys'}    
    plots.clear('energy_grey', '', False)
    plots.add_colour_plot(s)
    plots.make(globals.dirs['plots_svg'], globals.dirs['plots_eps'])
    
    s = {'title': 'Energy (Ry) vs Ecutwfc and Ecutrho', 'x_axis': 'Ecutrho/Ecutwfc', 'y_axis': 'Ecutwfc/Ry', 'x': x, 'y': y, 'Z': Z, 'cmap': 'viridis'}    
    plots.clear('energy_colour', '', False)
    plots.add_colour_plot(s)
    plots.make(globals.dirs['plots_svg'], globals.dirs['plots_eps'])
    
    Z = globals.ecut2d['energy_ry'][0:data_h,0:data_w] - numpy.amin(globals.ecut2d['energy_ry'][0:data_h,0:data_w])
     
    s = {'title': 'Energy (Ry) (Adjusted) vs Ecutwfc and Ecutrho', 'x_axis': 'Ecutrho/Ecutwfc', 'y_axis': 'Ecutwfc/Ry', 'x': x, 'y': y, 'Z': Z, 'cmap': 'Greys'}    
    plots.clear('energy_adjusted_grey', '', False)
    plots.add_colour_plot(s)
    plots.make(globals.dirs['plots_svg'], globals.dirs['plots_eps'])
    
    s = {'title': 'Energy (Ry) (Adjusted) vs Ecutwfc and Ecutrho', 'x_axis': 'Ecutrho/Ecutwfc', 'y_axis': 'Ecutwfc/Ry', 'x': x, 'y': y, 'Z': Z, 'cmap': 'viridis'}    
    plots.clear('energy_adjusted_colour', '', False)
    plots.add_colour_plot(s)
    plots.make(globals.dirs['plots_svg'], globals.dirs['plots_eps'])
    
    # ENERGY CONVERGENCE
    #####################
    
    x_wfc = conv_ecut2d.half_scale(x)
    Z = conv_ecut2d.conv_h(globals.ecut2d['energy_ry'][0:data_h,0:data_w] - numpy.amin(globals.ecut2d['energy_ry'][0:data_h,0:data_w]), data_h, data_w) 
    s = {'title': 'Energy (Ry) Convergence wrt Ecutwfc vs Ecutwfc and Ecutrho', 'x_axis': 'Ecutrho/Ecutwfc', 'y_axis': 'Ecutwfc/Ry', 'x': x_wfc, 'y': y, 'Z': Z, 'cmap': 'viridis'}  
    plots.clear('energy_ecutwfc_convergence_colour', '', False)
    plots.add_colour_plot(s)
    plots.make(globals.dirs['plots_svg'], globals.dirs['plots_eps'])
    
    
    y_rho = conv_ecut2d.half_scale(y)
    Z = conv_ecut2d.conv_w(globals.ecut2d['energy_ry'][0:data_h,0:data_w] - numpy.amin(globals.ecut2d['energy_ry'][0:data_h,0:data_w]), data_h, data_w) 
    s = {'title': 'Energy (Ry) Convergence wrt Ecutrho vs Ecutwfc and Ecutrho', 'x_axis': 'Ecutrho/Ecutwfc', 'y_axis': 'Ecutwfc/Ry', 'x': x, 'y': y_rho, 'Z': Z, 'cmap': 'viridis'}  
    plots.clear('energy_ecutrho_convergence_colour', '', False)
    plots.add_colour_plot(s)
    plots.make(globals.dirs['plots_svg'], globals.dirs['plots_eps'])
    
    
    
    # FORCES
    #####################
    
    x = numpy.linspace(globals.ecut2d['rho_min_pw'], globals.ecut2d['rho_max_pw'], data_h)
    y = numpy.linspace(globals.ecut2d['wfc_min_pw'], globals.ecut2d['wfc_max_pw'], data_w)
    Z = globals.ecut2d['force_rybohr'][0:data_h,0:data_w]
    
    
    s = {'title': 'Force (Ry/Bohr) vs Ecutwfc and Ecutrho', 'x_axis': 'Ecutrho/Ecutwfc', 'y_axis': 'Ecutwfc/Ry', 'x': x, 'y': y, 'Z': Z, 'cmap': 'Greys'}    
    plots.clear('force_grey', '', False)
    plots.add_colour_plot(s)
    plots.make(globals.dirs['plots_svg'], globals.dirs['plots_eps'])
    
    s = {'title': 'Force (Ry/Bohr) vs Ecutwfc and Ecutrho', 'x_axis': 'Ecutrho/Ecutwfc', 'y_axis': 'Ecutwfc/Ry', 'x': x, 'y': y, 'Z': Z, 'cmap': 'viridis'}    
    plots.clear('force_colour', '', False)
    plots.add_colour_plot(s)
    plots.make(globals.dirs['plots_svg'], globals.dirs['plots_eps'])
    
    
    Z = globals.ecut2d['force_rybohr'][0:data_h,0:data_w] - numpy.amin(globals.ecut2d['force_rybohr'][0:data_h,0:data_w])
    
    s = {'title': 'Force (Ry/Bohr) (Adjusted) vs Ecutwfc and Ecutrho', 'x_axis': 'Ecutrho/Ecutwfc', 'y_axis': 'Ecutwfc/Ry', 'x': x, 'y': y, 'Z': Z, 'cmap': 'Greys'}    
    plots.clear('force_adjusted_colour', '', False)
    plots.add_colour_plot(s)
    plots.make(globals.dirs['plots_svg'], globals.dirs['plots_eps'])
        
    s = {'title': 'Force (Ry/Bohr) (Adjusted) vs Ecutwfc and Ecutrho', 'x_axis': 'Ecutrho/Ecutwfc', 'y_axis': 'Ecutwfc/Ry', 'x': x, 'y': y, 'Z': Z, 'cmap': 'viridis'}    
    plots.clear('force_adjusted_colour', '', False)
    plots.add_colour_plot(s)
    plots.make(globals.dirs['plots_svg'], globals.dirs['plots_eps'])
    
    
    
    # FORCE CONVERGENCE
    #####################
    
    x_wfc = conv_ecut2d.half_scale(x)
    Z = conv_ecut2d.conv_h(globals.ecut2d['force_rybohr'][0:data_h,0:data_w] - numpy.amin(globals.ecut2d['force_rybohr'][0:data_h,0:data_w]), data_h, data_w) 
    s = {'title': 'Force (Ry/Bohr) Convergence wrt Ecutwfc vs Ecutwfc and Ecutrho', 'x_axis': 'Ecutrho/Ecutwfc', 'y_axis': 'Ecutwfc/Ry', 'x': x_wfc, 'y': y, 'Z': Z, 'cmap': 'viridis'}  
    plots.clear('force_ecutwfc_convergence_colour', '', False)
    plots.add_colour_plot(s)
    plots.make(globals.dirs['plots_svg'], globals.dirs['plots_eps'])
    
    
    y_rho = conv_ecut2d.half_scale(y)
    Z = conv_ecut2d.conv_w(globals.ecut2d['force_rybohr'][0:data_h,0:data_w] - numpy.amin(globals.ecut2d['force_rybohr'][0:data_h,0:data_w]), data_h, data_w) 
    s = {'title': 'Force (Ry/Bohr) Convergence wrt Ecutrho vs Ecutwfc and Ecutrho', 'x_axis': 'Ecutrho/Ecutwfc', 'y_axis': 'Ecutwfc/Ry', 'x': x, 'y': y_rho, 'Z': Z, 'cmap': 'viridis'}  
    plots.clear('force_ecutrho_convergence_colour', '', False)
    plots.add_colour_plot(s)
    plots.make(globals.dirs['plots_svg'], globals.dirs['plots_eps'])
    
    
    
    
    # FORCE CONVERGENCE
    #####################
    
    
    x = numpy.linspace(globals.ecut2d['rho_min_pw'], globals.ecut2d['rho_max_pw'], data_h)
    y = numpy.linspace(globals.ecut2d['wfc_min_pw'], globals.ecut2d['wfc_max_pw'], data_w)
    
    x_wfc = conv_ecut2d.half_scale(x)
    y_rho = conv_ecut2d.half_scale(y)
    
    Za = conv_ecut2d.conv_h(globals.ecut2d['energy_ry'][0:data_h,0:data_w] - numpy.amin(globals.ecut2d['energy_ry'][0:data_h,0:data_w]), data_h, data_w)
    Zb = conv_ecut2d.conv_w(globals.ecut2d['energy_ry'][0:data_h,0:data_w] - numpy.amin(globals.ecut2d['energy_ry'][0:data_h,0:data_w]), data_h, data_w) 
    Zc = conv_ecut2d.conv_h(globals.ecut2d['force_rybohr'][0:data_h,0:data_w] - numpy.amin(globals.ecut2d['force_rybohr'][0:data_h,0:data_w]), data_h, data_w) 
    Zd = conv_ecut2d.conv_w(globals.ecut2d['force_rybohr'][0:data_h,0:data_w] - numpy.amin(globals.ecut2d['force_rybohr'][0:data_h,0:data_w]), data_h, data_w) 
    
    sa = {'title': 'Energy (Ry) Convergence wrt Ecutwfc', 'x_axis': 'Ecutrho/Ecutwfc', 'y_axis': 'Ecutwfc/Ry', 'x': x_wfc, 'y': y, 'Z': Za, 'cmap': 'Greys'} 
    sb = {'title': 'Energy (Ry) Convergence wrt Ecutrho', 'x_axis': 'Ecutrho/Ecutwfc', 'y_axis': 'Ecutwfc/Ry', 'x': x, 'y': y_rho, 'Z': Zb, 'cmap': 'Greys'} 
    sc = {'title': 'Force (Ry/Bohr) Convergence wrt Ecutwfc', 'x_axis': 'Ecutrho/Ecutwfc', 'y_axis': 'Ecutwfc/Ry', 'x': x_wfc, 'y': y, 'Z': Zc, 'cmap': 'Greys'} 
    sd = {'title': 'Force (Ry/Bohr) Convergence wrt Ecutrho', 'x_axis': 'Ecutrho/Ecutwfc', 'y_axis': 'Ecutwfc/Ry', 'x': x, 'y': y_rho, 'Z': Zd, 'cmap': 'Greys'}    
    plots.clear('summary_grey', 'Convergence vs Ecutwfc and Ecutrho', False)
    plots.add_colour_plot(sa)
    plots.add_colour_plot(sb)
    plots.add_colour_plot(sc)
    plots.add_colour_plot(sd)
    plots.make(globals.dirs['plots_svg'], globals.dirs['plots_eps'])
    
    
    sa = {'title': 'Energy (Ry) Convergence wrt Ecutwfc', 'x_axis': 'Ecutrho/Ecutwfc', 'y_axis': 'Ecutwfc/Ry', 'x': x_wfc, 'y': y, 'Z': Za, 'cmap': 'viridis'} 
    sb = {'title': 'Energy (Ry) Convergence wrt Ecutrho', 'x_axis': 'Ecutrho/Ecutwfc', 'y_axis': 'Ecutwfc/Ry', 'x': x, 'y': y_rho, 'Z': Zb, 'cmap': 'viridis'} 
    sc = {'title': 'Force (Ry/Bohr) Convergence wrt Ecutwfc', 'x_axis': 'Ecutrho/Ecutwfc', 'y_axis': 'Ecutwfc/Ry', 'x': x_wfc, 'y': y, 'Z': Zc, 'cmap': 'viridis'} 
    sd = {'title': 'Force (Ry/Bohr) Convergence wrt Ecutrho', 'x_axis': 'Ecutrho/Ecutwfc', 'y_axis': 'Ecutwfc/Ry', 'x': x, 'y': y_rho, 'Z': Zd, 'cmap': 'viridis'}    
    plots.clear('summary_colour', 'Convergence vs Ecutwfc and Ecutrho', False)
    plots.add_colour_plot(sa)
    plots.add_colour_plot(sb)
    plots.add_colour_plot(sc)
    plots.add_colour_plot(sd)
    plots.make(globals.dirs['plots_svg'], globals.dirs['plots_eps'])
    
    
    
    
    
    
    #############################################################################################
    #############################################################################################    
    # K-POINTS
    #############################################################################################
    #############################################################################################
    
    
    
    data_h = globals.kconv['data_h']
    data_w = globals.kconv['data_w']


    y = numpy.linspace(globals.kconv['k_min'], globals.kconv['k_max'], data_w)   
    x = numpy.zeros((len(globals.kconv['smearing_pw']),),)
    x[:] = globals.kconv['smearing_pw'][:]
    Z = globals.kconv['energy_ry'][0:data_h,0:data_w]    
    
    s = {'title': 'Energy (Ry) vs K-points and Degauss', 'x_axis': 'Degauss', 'y_axis': 'K-points', 
         'x': x, 'y': y, 'Z': Z, 'cmap': 'Greys'}  
    plots.clear('kpoints_energy_grey', '', False)
    plots.add_colour_plot(s)
    plots.make(globals.dirs['plots_svg'], globals.dirs['plots_eps'])

    s = {'title': 'Energy (Ry) vs K-points and Degauss', 'x_axis': 'Degauss', 'y_axis': 'K-points',  
         'x': x, 'y': y, 'Z': Z, 'cmap': 'viridis'}  
    plots.clear('kpoints_energy_colour', '', False)
    plots.add_colour_plot(s)
    plots.make(globals.dirs['plots_svg'], globals.dirs['plots_eps'])
    
    
    Z = globals.kconv['energy_ry'][0:data_h,0:data_w] - numpy.amin(globals.kconv['energy_ry'][0:data_h,0:data_w])
    
    s = {'title': 'Energy (Ry) (Adjusted) vs K-points and Degauss', 'x_axis': 'Degauss', 'y_axis': 'K-points',  
         'x': x, 'y': y, 'Z': Z, 'cmap': 'Greys'}  
    plots.clear('kpoints_energy_adjusted_grey', '', False)
    plots.add_colour_plot(s)
    plots.make(globals.dirs['plots_svg'], globals.dirs['plots_eps'])

    s = {'title': 'Energy (Ry) (Adjusted) vs K-points and Degauss', 'x_axis': 'Degauss', 'y_axis': 'K-points', 
         'x': x, 'y': y, 'Z': Z, 'cmap': 'viridis'}  
    plots.clear('kpoints_energy_adjusted_colour', '', False)
    plots.add_colour_plot(s)
    plots.make(globals.dirs['plots_svg'], globals.dirs['plots_eps'])

    
    
    x_degauss = conv_kpoints.half_scale(x)
    Z = conv_kpoints.conv_h(globals.kconv['energy_ry'][0:data_h,0:data_w] - numpy.amin(globals.kconv['energy_ry'][0:data_h,0:data_w]), data_h, data_w) 
        
    s = {'title': 'Energy (Ry) Convergence wrt Degauss vs K-points and Degauss', 'x_axis': 'Degauss', 'y_axis': 'K-points', 
         'x': x_degauss, 'y': y, 'Z': Z, 'cmap': 'Greys'}  
    plots.clear('kpoints_degauss_energy_convergence_grey', '', False)
    plots.add_colour_plot(s)
    plots.make(globals.dirs['plots_svg'], globals.dirs['plots_eps'])

    s = {'title': 'Energy (Ry) Convergence wrt Degauss vs K-points and Degauss', 'x_axis': 'Degauss', 'y_axis': 'K-points', 
         'x': x_degauss, 'y': y, 'Z': Z, 'cmap': 'viridis'} 
    plots.clear('kpoints_degauss_energy_convergence_colour', '', False)
    plots.add_colour_plot(s)
    plots.make(globals.dirs['plots_svg'], globals.dirs['plots_eps'])

    
    
    y_kpoints = conv_kpoints.half_scale(y)
    Z = conv_kpoints.conv_w(globals.kconv['energy_ry'][0:data_h,0:data_w] - numpy.amin(globals.kconv['energy_ry'][0:data_h,0:data_w]), data_h, data_w) 
        
    s = {'title': 'Energy (Ry) Convergence wrt Degauss vs K-points and Degauss', 'x_axis': 'Degauss', 'y_axis': 'K-points', 
         'x': x, 'y': y_kpoints, 'Z': Z, 'cmap': 'Greys'}  
    plots.clear('kpoints_kpoint_energy_convergence_grey', '', False)
    plots.add_colour_plot(s)
    plots.make(globals.dirs['plots_svg'], globals.dirs['plots_eps'])

    s = {'title': 'Energy (Ry) Convergence wrt Degauss vs K-points and Degauss', 'x_axis': 'Degauss', 'y_axis': 'K-points', 
         'x': x, 'y': y_kpoints, 'Z': Z, 'cmap': 'viridis'} 
    plots.clear('kpoints_kpoint_energy_convergence_colour', '', False)
    plots.add_colour_plot(s)
    plots.make(globals.dirs['plots_svg'], globals.dirs['plots_eps'])






    y = numpy.linspace(globals.kconv['k_min'], globals.kconv['k_max'], data_w)   
    x = numpy.zeros((len(globals.kconv['smearing_pw']),),)
    x[:] = globals.kconv['smearing_pw'][:]
    Z = globals.kconv['force_rybohr'][0:data_h,0:data_w]    
    
    s = {'title': 'Force (Ry/Bohr) vs K-points and Degauss', 'x_axis': 'Degauss', 'y_axis': 'K-points', 
         'x': x, 'y': y, 'Z': Z, 'cmap': 'Greys'}  
    plots.clear('kpoints_force_grey', '', False)
    plots.add_colour_plot(s)
    plots.make(globals.dirs['plots_svg'], globals.dirs['plots_eps'])
    
    s = {'title': 'Force (Ry/Bohr) vs K-points and Degauss', 'x_axis': 'Degauss', 'y_axis': 'K-points', 
         'x': x, 'y': y, 'Z': Z, 'cmap': 'viridis'}  
    plots.clear('kpoints_force_colour', '', False)
    plots.add_colour_plot(s)
    plots.make(globals.dirs['plots_svg'], globals.dirs['plots_eps'])
    

    Z = globals.kconv['force_rybohr'][0:data_h,0:data_w] - numpy.amin(globals.kconv['force_rybohr'][0:data_h,0:data_w])
    
    s = {'title': 'Force (Ry/Bohr) (Adjusted) vs K-points and Degauss', 'x_axis': 'Degauss', 'y_axis': 'K-points', 
         'x': x, 'y': y, 'Z': Z, 'cmap': 'Greys'}  
    plots.clear('kpoints_force_grey', '', False)
    plots.add_colour_plot(s)
    plots.make(globals.dirs['plots_svg'], globals.dirs['plots_eps'])
    
    s = {'title': 'Force (Ry/Bohr) (Adjusted) vs K-points and Degauss', 'x_axis': 'Degauss', 'y_axis': 'K-points', 
         'x': x, 'y': y, 'Z': Z, 'cmap': 'viridis'}  
    plots.clear('kpoints_force_colour', '', False)
    plots.add_colour_plot(s)
    plots.make(globals.dirs['plots_svg'], globals.dirs['plots_eps'])


    
    
    x_kpoints = conv_kpoints.half_scale(x)
    Z = conv_kpoints.conv_h(globals.kconv['energy_ry'][0:data_h,0:data_w] - numpy.amin(globals.kconv['energy_ry'][0:data_h,0:data_w]), data_h, data_w) 
        
    s = {'title': 'Force (Ry/Bohr) Convergence wrt Degauss vs K-points and Degauss', 'x_axis': 'Degauss', 'y_axis': 'K-points', 
         'x': x_kpoints, 'y': y, 'Z': Z, 'cmap': 'Greys'}  
    plots.clear('kpoints_degauss_force_convergence_grey', '', False)
    plots.add_colour_plot(s)
    plots.make(globals.dirs['plots_svg'], globals.dirs['plots_eps'])

    s = {'title': 'Force (Ry/Bohr) Convergence wrt Degauss vs K-points and Degauss', 'x_axis': 'Degauss', 'y_axis': 'K-points', 
         'x': x_kpoints, 'y': y, 'Z': Z, 'cmap': 'viridis'} 
    plots.clear('kpoints_degauss_force_convergence_colour', '', False)
    plots.add_colour_plot(s)
    plots.make(globals.dirs['plots_svg'], globals.dirs['plots_eps'])

    y_kpoints = conv_kpoints.half_scale(y)
    Z = conv_kpoints.conv_w(globals.kconv['energy_ry'][0:data_h,0:data_w] - numpy.amin(globals.kconv['energy_ry'][0:data_h,0:data_w]), data_h, data_w) 
        
    s = {'title': 'Force (Ry/Bohr) Convergence wrt Kpoint vs K-points and Degauss', 'x_axis': 'Degauss', 'y_axis': 'K-points', 
         'x': x, 'y': y_kpoints, 'Z': Z, 'cmap': 'Greys'}  
    plots.clear('kpoints_kpoint_force_convergence_grey', '', False)
    plots.add_colour_plot(s)
    plots.make(globals.dirs['plots_svg'], globals.dirs['plots_eps'])

    s = {'title': 'Force (Ry/Bohr) Convergence wrt Kpoint vs K-points and Degauss', 'x_axis': 'Degauss', 'y_axis': 'K-points', 
         'x': x, 'y': y_kpoints, 'Z': Z, 'cmap': 'viridis'} 
    plots.clear('kpoints_kpoint_force_convergence_colour', '', False)
    plots.add_colour_plot(s)
    plots.make(globals.dirs['plots_svg'], globals.dirs['plots_eps'])
    
    
    
    
    
    
    
    
    