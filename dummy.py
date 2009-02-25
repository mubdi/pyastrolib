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

def dummy1():
  """This function does nothing
  >>> dummy1() == "hello"
  True
  >>> type(dummy1())
  <type 'str'>"""
  return "hello"
  
def dummy2(a,b):
  """This function needs more explination... talk about what comes in and what comes out.  Next comes the doctest:
  
  >>> dummy2(1,2) == None
  True"""
  return None
  
def dummy3(input):
  """Example of complicated docstring:
  >>> import numpy as n
  >>> x = n.arange(.1,1,.4)
  >>> dummy3(x) == x
  array([ True,  True,  True], dtype=bool)"""
  return input

if __name__ == '__main__':
  # The following two lines will test all functions in the module.
  # Run `python astro.py -v` to see verbose output.  It's probably best to rely
  # on this kind of testing once you believe the code to already be functional
  # since it is a great method for building but difficult to work with durring
  # development.
  import doctest
  doctest.testmod()
  