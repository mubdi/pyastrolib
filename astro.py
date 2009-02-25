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
    print "ERROR - RA and Declination do not have equal number of elements"
    return false




  return true


if __name__ == '__main__':
  # The following two lines will test all functions in the module.
  # Run `python astro.py -v` to see verbose output.  It's probably best to rely
  # on this kind of testing once you believe the code to already be functional
  # since it is a great method for building but difficult to work with durring
  # development.
  import doctest
  doctest.testmod()
