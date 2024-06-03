### This work is licensed under a Creative Commons Attribution 4.0 International License
### https://creativecommons.org/licenses/by/4.0/
### Copyright (c) 2021, Julien Langou. All rights reserved.

### def x_to_the_power_of_n( x, n ):
###   return y

### purpose
# `x_to_the_power_of_n` computes `y = x**n` (`x` to the power of `n`)
# `n` must be an integer

### Uses binary expansion of `n` and the powers of the powers of 2 of `x` (x^1, x^2, x^4, etc.) 
### Computational complexity is in between `log2(n)` and `2*log2(n)` (depending on how many 1's are in the base 2 of `n`)

def x_to_the_power_of_n( x, n ):
  
  m = n
  y = 1
  z = x
  
  if ( m%2 ):
    y = y * z
  m = m//2

  while( m >= 1):
    z = z * z
    if ( m%2 ):
      y = y * z
    m = m//2

  return y
