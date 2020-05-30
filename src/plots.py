###########################################
#  CLASS plot
###########################################
class plots:
  
  d = None
  
  def clear(filename=None, title=None, threed=None):
    plots.d = {
              'filename': 'plot',
              'pw': 1,
              'ph': 1,
              'title': '',
              'subplots': [],
              '3d': False,
              }
    if(filename is not None):
      plots.d['filename'] = filename 
    if(title is not None):
      plots.d['title'] = title 
    if(threed is not None):
      plots.d['3d'] = threed 

# Custom colour maps?
    colors = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
    cmap_name = 'rgb_list'
    cm = LinearSegmentedColormap.from_list(cmap_name, colors, N=128)

  def add_colour_plot(s_in):

    s = {
        'x': None,
        'y': None,
        'Z': None,
        'x_axis': '',
        'y_axis': '',
        'title': '',
        'cmap': 'viridis',
        }

    plots.d['subplots'].append(s)

# Read in data
    for k in s_in.keys():
      if(k in s.keys()):
        s[k] = s_in[k]

  def make(out_dir_eps=None, out_dir_svg=None):
    if(out_dir_eps is not None):
      out_dir_eps = out_dir_eps + '/'
    else:
      out_dir_eps = ''

    if(out_dir_svg is not None):
      out_dir_svg = out_dir_svg + '/'
    else:
      out_dir_svg = ''

    title = plots.d['title']
    size = len(plots.d['subplots'])
    pw = plots.d['pw']
    ph = plots.d['ph']
    p_size = pw * ph

    if(size > pw * ph):
      pw = int(numpy.ceil(numpy.sqrt(size)))
      ph = int(numpy.ceil(size / pw))
      p_size = pw * ph

####################
#  2D colour map
####################
  
# plot settings
    plt.clf()    
    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')
    plt.xticks(fontsize=9)
    fig, axs = plt.subplots(ph, pw, figsize=(12,9))
    fig.tight_layout(pad=5.0)
    fig.suptitle(title)  
    for n in range(size):      
      sub_title = plots.d['subplots'][n]['title']
      x_axis = plots.d['subplots'][n]['x_axis']
      y_axis = plots.d['subplots'][n]['y_axis']
      x = plots.d['subplots'][n]['x']
      y = plots.d['subplots'][n]['y']
      Z = plots.d['subplots'][n]['Z']
      if(p_size == 1):
        ax_n = axs
      elif(p_size > 1 and ph == 1):
        ax_n = axs[n]
      else:
        ax_n = axs[int(numpy.floor(n / pw)), n % pw]

      cmap = plots.d['subplots'][n]['cmap']
      ax_n.set_title(sub_title)

      Z_plot = numpy.copy(Z)
      cf = ax_n.contourf(x, y, Z_plot, 128, cmap=cmap)
      ax_n.set_xlabel(x_axis)
      ax_n.set_ylabel(y_axis)
      plt.colorbar(cf, ax=ax_n)
      if(n == size):
        break

    plt.savefig(out_dir_svg + plots.d['filename'] + '.svg', format='svg')
    plt.savefig(out_dir_eps + plots.d['filename'] + '.eps', format='eps')
    plt.close('all')

####################
#  3D colour map
####################

    if(plots.d['3d']):

# plot settings
      plt.clf()    
      plt.rc('font', family='serif')
      plt.rc('xtick', labelsize='x-small')
      plt.rc('ytick', labelsize='x-small')

      fig = plt.figure(figsize=(12,9))
      fig.suptitle(title)      #plt.axes()
      
      for n in range(size):
        x = plots.d['subplots'][n]['x']
        y = plots.d['subplots'][n]['y']

        if(len(Z) != len(Z[0,:])):
          x_new, y_new, Z_new = plots.square_data(x, y, Z)
        else:
          x_new, y_new, Z_new = x, y, Z

        xm = numpy.outer(x_new, numpy.ones(len(x_new)))
        ym = numpy.outer(y_new, numpy.ones(len(y_new))).transpose()

        sub_title = plots.d['subplots'][n]['title']
        x_axis = plots.d['subplots'][n]['x_axis']
        y_axis = plots.d['subplots'][n]['y_axis']

        ax = fig.add_subplot(ph, pw, n+1, projection='3d')
        ax.set_title(sub_title)
        ax.set_xlabel(x_axis)
        ax.set_ylabel(y_axis)
        cmap = plots.d['subplots'][n]['cmap']

        cs = ax.plot_surface(xm, ym, Z_new, rstride=1, cstride=1, cmap=cmap,
                         linewidth=0, antialiased=False)
        fig.colorbar(cs, shrink=0.5, aspect=10)

      plt.savefig(out_dir_svg + plots.d['filename'] + '_3d.svg', format='svg')
      plt.savefig(out_dir_eps + plots.d['filename'] + '_3d.eps', format='eps')
      plt.close('all')

  def interp(xi, x, y):
    yi = 0.0
    for i in range(4):
      li = 1.0
      for j in range(4):
        if(i != j):
          li = li * (xi - x[j]) / (x[i] - x[j])
      yi = yi + li * y[i]
    return yi

  def stretch(x, y, l):
    x_out = numpy.linspace(x[0], x[-1], l)
    y_out = numpy.zeros((l,),)
    xn = 0
    for n in range(l):
      while(not (x_out[n] <= x[xn+1] and x_out[n] >= x[xn])):
        xn = xn + 1
      nn = xn
      if(nn < 0):
        nn = 0
      elif(nn + 4 > len(x)):
        nn = len(x) - 4
      y_out[n] = plots.interp(x_out[n], x[nn:nn+4], y[nn:nn+4])      
    return x_out, y_out

  def square_data(x, y, Z):
    if(len(x) == len(y)):
      return Z
    if(len(x) > len(y)):
      y_out = numpy.linspace(y[0], y[-1], len(x))
      Z_new = numpy.zeros((len(x), len(x),),)
      for i in range(len(x)):
        y_stretch, Z_new[i, :] = plots.stretch(y, Z[i, :], len(x))
      return x, y_out, Z_new
    if(len(x) < len(y)):
      x_out = numpy.linspace(x[0], x[-1], len(y))
      Z_new = numpy.zeros((len(y), len(y),),)
      for i in range(len(y)):
        x_stretch, Z_new[:, i] = plots.stretch(x, Z[:, i], len(y))
      return x_out, y, Z_new

"""

Example use

x = numpy.outer(numpy.linspace(-2, 2, 30), numpy.ones(30))
y = x.copy().T # transpose
Z = numpy.cos(x ** 2 + y ** 2)

sa = {
    'title': 'A',
    'x': x,
    'y': y,
    'Z': Z,
    'x_axis': 'X',
    'y_axis': 'Y',
    }
sb = {
    'title': 'B',
    'x': x,
    'y': y,
    'Z': Z,
    }
sc = {
    'title': 'C',
    'x': x,
    'y': y,
    'Z': Z,
    }
sd = {
    'title': 'D',
    'x': x,
    'y': y,
    'Z': Z,
    }
se = {
    'title': 'E',
    'x': x,
    'y': y,
    'Z': Z,
    }

plots.clear()
plots.add_colour_plot(sa)
plots.add_colour_plot(sb)
plots.add_colour_plot(sc)
plots.add_colour_plot(sd)
plots.add_colour_plot(se)
plots.make()

plots.clear('new')
plots.add_colour_plot(sa)
plots.make()

"""
