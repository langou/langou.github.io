def polyval_nested_w_base_points( c, b, x ):
  d = np.size(c)
  px = c[ d-1 ] * np.ones( np.shape(x) )
  for i in range( d-2, -1, -1 ):
    px = px * ( x - b[i] ) + c[i]
  return px
