class square:

  def interp(xi, x, y, np):
    yi = 0.0
    for i in range(np):
      li = 1.0
      for j in range(np):
        if(i != j):
          li = li * (xi - x[j]) / (x[i] - x[j])
      
      yi = yi + li * y[i]
    return yi
    

  def stretch(x, y, l, np):
    x_out = numpy.linspace(x[0], x[-1], l)
    y_out = numpy.zeros((l,),)
    xn = 0
    for n in range(l):
      while(not (x_out[n] <= x[xn+1] and x_out[n] >= x[xn])):
        xn = xn + 1
      nn = xn
      if(nn < 0):
        nn = 0
      elif(nn + np > len(x)):
        nn = len(x) - np
      y_out[n] = square.interp(x_out[n], x[nn:nn+np], y[nn:nn+np], np)      
    return x_out, y_out

  def square_data(Z, l, w, np=4):
    lz = len(Z)
    wz = len(Z[0,:])
    Z_at = numpy.zeros((lz, w,),)
    Z_a = numpy.zeros((l, w,),)
    Z_bt = numpy.zeros((l, wz,),)
    Z_b = numpy.zeros((l, w,),)

    x = numpy.linspace(0.0, 1.0, wz)
    y = numpy.linspace(0.0, 1.0, lz)
    
    for li in range(lz): 
      x_out, Z_at[li,:] = square.stretch(x, Z[li,:], w, np)  
    for wi in range(w): 
      y_out, Z_a[:,wi] = square.stretch(y, Z_at[:,wi], l, np)


    x = numpy.linspace(0.0, 1.0, wz)
    y = numpy.linspace(0.0, 1.0, lz)
    
    for wi in range(wz): 
      y_out, Z_bt[:,wi] = square.stretch(y, Z[:,wi], l, np)
    for li in range(l): 
      x_out, Z_b[li,:] = square.stretch(x, Z_bt[li,:], w, np)

    Z_out = 0.5 * (Z_b + Z_a)
    
    return Z_out