#!/usr/bin/env python
# This file is part of pyAstroLib.
#
# pyAstroLib is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# pyAstroLib is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU Lesser Public License
# along with pyAstroLib.  If not, see <http://www.gnu.org/licenses/>.

import numpy as n

def adstring(ra_dec, dec="", precision="", truncate=""):
  

  # Testing Input Parameters

  if n.size(ra_dec) == 1:
    if not dec:
      dec = ra_dec
    else:
      ra = ra_dec
  elif n.size(ra_dec) == 2:
    if n.size(dec) < 2:
      ra = n.mod(ra_dec[0], 360)
      if (not precision):
        precision = dec
      dec = ra_dec[1]
  else: ra = ra_dec
  
  if n.size(ra) != n.size(dec):
    raise TypeError, 'ERROR - RA and Declination do not have equal number of elements'

  if n.size(ra) == n.size(dec):
      badrange = n.where( (dec > 90 ) or (dec < -90)  )
      if n.size(badrange) > 0:
        print "Warning:  Some declination values are out of valid range (-90 < dec <90)"
                    



  return 0

def radec(ra, dec, hours=""):
  """radec(ra, dec, hours="")
  Converts RA and Dec from decimal to sexigesimal units
  Returns a tuple (ihr, imin, xsec, imn, wsc)
  
  INPUTS:
    ra - right ascension, float or array, in degrees unless 
         "hours" is set
    dec - declination in decimal degrees, float or array, same
          number of elements as ra     
  
  OPTIONAL INPUT:
    hours - if set to true, then the right ascension input should
            be set to hours instead of degrees
            
  OUTPUTS:
    ihr - right ascension hours (float or array)
    imin - right ascension minutes (float or array)
    xsec - right ascension seconds (float or array)
    ideg - declination degrees (float or array)
    imn - declination minutes (float or array)
    xsc - declination seconds (float or array)         
  
  >>> radec(0,0) 
  array(0,0,0,0,0)
  """


  # Compute RA
  if(hours):
    ra =  n.mod(ra, 24)
    ra = ra + 24*(n.less(ra, 0) )
    ihr = n.fix(ra)
    xmin = n.abs(ra*60.0 - ihr*60.0)
  else:
    ra = n.mod(ra, 360)
    ra = ra + 360*(n.less(ra, 0))
    ihr = n.fix(ra/15.0)
    xmin = n.abs(ra*4.0 - ihr*60.0)

  imin = n.fix(xmin)
  xsec = (xmin - imin)*60.0

  # Compute Dec
  ideg = n.fix(dec)
  xmn = n.abs(dec - ideg)*60.0
  imn = n.fix(xmn)
  xsc = (xmn - imn)*60.0

  # Testing for Special Case of Zero Degrees

  zero_deg = n.equal(ideg, 0)  & n.less(dec, 0) 
  imn = imn - 2*imn*n.fix( zero_deg * ( n.not_equal(imn,0) ) )
  xsc = xsc - 2 *xsc*zero_deg*(n.equal(imn, 0) )  

  return ihr, imin, xsec, ideg, imn, xsc

def sixty(scalar):
    """ sixty(scalar)
    Converts a decimal number to sexigesimal
    returns a three element array
    
    INPUTS: 
      scalar - Decimal Quantity
      
    OUTPUTS: 
      result - three element array of sexigesimal equivalent.
               If the scalar is negative, the first non-zero
               element in the result will have a negative sign
               
    COMPATIBILITY NOTES:
      1. Since python does not support negative zeros, the trail
         sign keyword from the IDL procedure has been removed
    
    >>>sixty(53)
    [53.0, 0.0, 0.0]
    """
    if n.size(scalar) != 1:
            raise TypeError, 'ERROR - Parameter must contain 1 element'
    
    ss = n.abs(3600.0 * scalar)
    mm = n.abs(60.0 * scalar)
    dd = n.abs(scalar)
    result = n.zeros(3, float)  
    
    result[0] = n.fix(dd)
    result[1] = n.fix(mm - 60.0*result[0] )
    result[2] = ss - 3600.0*result[0] - 60.0*result[1]
    
    if scalar < 0.0:
        if result[0] != 0:
            result[0] = -result[0]
        elif result[1] != 0:
            result[1] = -result[1]
        else: 
            result[2] = -result[2]
    
    return result

def ten(dd, mm="", ss=""):
    
    return result

if __name__ == '__main__':
  # The following two lines will test all functions in the module.
  # Run `python astro.py -v` to see verbose output.  It's probably best to rely
  # on this kind of testing once you believe the code to already be functional
  # since it is a great method for building but difficult to work with durring
  # development.
  import doctest
  doctest.testmod()
