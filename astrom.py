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
import pyfits as f

#Defining Angle Constants
pi = n.pi
pi2 = n.pi/2.0
radeg = 180.0/pi

def wcs_getpole(crval, lonpole, theta0, latpole=""):
    """ wcs_getpole(crval, lonpole, theta0, latpole="")
    Compute the coordinates of the native pole for a non-polar projection
    Returns at tuple (alpha_p, delta_p)
    
    INPUTS: 
      crval - 2 element vector containing standard system coordinates
              (the longitude and latitude) of the reference point in
              degrees
      lonpole - native longitude of the celestial North Pole (degrees)
      theta0 - native latitude of the fiducial point
      
    OPTIONAL INPUTS:
      latpole - native latitude of the celestial North Pole (degrees)
      
    OUTPUTS:
      alpha_p - celestial longitude of the native pole (degrees)
          
      delta_p - celestial latitude of the native pole (degrees)
      
    NOTES: 
      1. Converted from the IDL astrolib procedure, last updated
         February 2004.
      2. (From IDL:) For non-polar (cylindrical or conic) projections, the 
         native pole is not at the reference point, and WCS_GETPOLE is used 
         to determine the position of the native pole.
         See section 2.4 of the paper 
         "Representation of Celestial Coordinates in FITS" by Calabretta 
         Greisen (2002, A&A, 395, 1077, also available at  
         http://fits.gsfc.nasa.gov/fits_wcs.html    Called by WCS_ROTATE

    """
    
    alpha_0 = crval[0]/radeg
    delta_0 = crval[1]/radeg
    
    if theta0 == 90:
        return alpha_0, delta_0

    phi_p = lonpole/radeg
    sp = n.sin(phi_p)
    cp = n.cos(phi_p)
    sd = n.sin(delta_0)
    cd = n.cos(delta_0)
    tand = tan(delta_0)
    
    if theta0 == 0:
        if (delta_0 == 0) and (lonpole == 90.0):
            delta_p = latpole
        else: 
            delta_p = n.arccos(sd/cp) 
        
        if latpole != 90.0:
            if n.abs(latpole + delta_p) < n.abs(latpole - delta_p):
                delta_p = -delta_p 
        
        if (lonpole == 180.0) or (cd == 0.0):
            alpha_p = alpha0
        else:
            alpha_p = alpha_0 - n.arctan2(sp/cd, -n.tan(delta_p)*tand)
    else:
        
        ctheta = n.cos(theta0/radeg)
        stheta = n.sin(theta0/radeg)
        
        term1 = n.arctan2(stheta, ctheta*cp)
        term2 = n.arccos(sd/( n.sqrt(1.0 - (ctheta**2 * sp**2) ) ) )
        
        if term2 == 0.0: 
            delta_p = term1
        else:
            delta_p1 = n.abs( (term1 + term2) *radeg)
            delta_p2 = n.abs( (term1 - term2) *radeg)

            if (delta_p1 > 90.0) and (delta_p2 > 90.0): 
                raise RuntimeError, 'No Valid Solution'  
            elif (delta_p1 <= 90.0) and (delta_p2 > 90.0):
                delta_p = term1 + term2
            elif (delta_p1 > 90.0) and (delta_p2 <= 90.0):
                delta_p = term1 - term2
            else:
                delta_p1 = (term1 + term2)*radeg
                delta_p2 = (term1 - term2)*radeg
                if n.abs(latpole - delta_p1) < n.abs(latpole - delta_p2):
                    delta_p = term1 + term2
                else:
                    delta_p = term1 - term2
            
            if (cd == 0.0):
                alpha_p = alpha_0
            else:
                sdelt = n.sin(delta_p)
                if sdelt == 1: 
                    alpha_p = alpha_0 - phi_p - pi
                elif sdelt == -1:
                    alpha_p = alpha_0 - phi_p
                else:
                    alpha_p = alpha_0 - n.arctan2( (stheta-n.sin(delta_p)*sd)/(n.cos(delta_p)*cd), sp*ctheta/cd )
            
    return alpha_p, delta_p

def wcs_rotate(x, y, crval, longpole=180.0, latpole="", reverse="", theta0=""):
    """ wcs_rotate(x, y, crval, lonpole=180.0, latpole="", reverse="", theta0="")
    Rotate between standard (e.g. celestial) and native coordinates
    Returns at tuple (a, b)
    
    INPUTS: 
      x - longitude of data in the native system (degrees) or in the
          standard celestial coordinate system (degrees) if reverse
          is set, scalar or vector
      y - latitude of data in the native system (degrees) or in the
          standard celestial coordinate system (degrees) if reverse
          is set, scalar or vector
      crval - 2 element vector containing standard system coordinates
              (the longitude and latitude) of the reference point in
              degrees

      
    OPTIONAL INPUTS:
      lonpole - native longitude of the celestial North Pole (degrees)
      latpole - native latitude of the celestial North Pole (degrees)
      reverse - if set to non-zero, then the input are in the standard
                system and the output are in the native system
      theta0 - native latitude of the fiducial point
      
    OUTPUTS:
      a - longitude of data in the standard celestial coordinate system 
          (degrees) or in the native system (degrees) if reverse is set,
          scalar or vector
      b - latitude of data in the standard celestial coordinate system 
          (degrees) or in the native system (degrees) if reverse is set,
          scalar or vector
          
    NOTES: 
      1. Converted from the IDL astrolib procedure, last updated
         May 2004.
      2. (From IDL:) Computes a spherical coordinate rotation between 
         native coordinates and  standard celestial coordinate system 
         (celestial, Galactic, or ecliptic).   Applies the equations 
         in Appendix B of the paper "Representation of Celestial 
         Coordinates in FITS" by Calabretta Greisen (2002, A&A, 
         395, 1077). Also see http://fits.gsfc.nasa.gov/fits_wcs.html

    """

    if reverse: 
        phi, theta = n.array(x), n.array(y)
    else: 
        longitude, latitude = n.array(x), n.array(y)

    
    # Longpole is the longitude in the native system of the North Pole
    # in the standard system
    
    phi_p = longpole/radeg
    sp = n.sin(phi_p)
    cp = n.cos(phi_p)
    
    if theta0 == 90.0:
        alpha_p = crval[0]/radeg
        delta_p = crval[1]/radeg
    else:
        alpha_p, delta_p = wcs_getpole(crval, longpole, theta0, latpole=latpole)
        
    # compute useful quantities relating to reference angles
    sa = n.sin(alpha_p)
    ca = n.cos(alpha_p)
    sd = n.sin(delta_p)
    cd = n.cos(delta_p)  
    
    # calculate rotation matrix
    
    r = n.array([[-sa*sp - ca*cp*sd,  ca*sp - sa*cp*sd, cp*cd ] ,
                 [ sa*cp - ca*sp*sd, -ca*cp - sa*sp*sd, sp*cd ] ,
                 [ ca*cd           ,  sa*cd           , sd    ] ])

    # solve the set of equations for each datum point
    
    if reverse:
        latitude = phi
        longitude = theta
        g = n.where( n.isfinite(phi) & n.isfinite(theta) )
        Ng = n.size(g)
        if Ng == 0: return
        phi1 = phi[g]/radeg
        theta1 = theta[g]/radeg        
        r = r.transpose()
    else:
        phi = longitude
        phi1 = longitude/radeg 
        theta1 = latitude/radeg
    
    l = n.cos(theta1)*n.cos(phi1)
    m = n.cos(theta1)*n.sin(phi1)
    n = n.sin(theta1)
    
    # Find Solution to the system of equations and put it in b
    # Can't use matrix notation in case l,m,n are vectors
    
    b0 = n.array(r[0,0]*l + r[1,0]*m + r[2,0]*n)
    b1 = n.array(r[0,1]*l + r[1,1]*m + r[2,1]*n)
    b2 = n.array(r[0,2]*l + r[1,2]*m + r[2,2]*n)
    
    b2[n.where(b2 < -1.0)] = -1.0
    b2[n.where(b2 > 1.0)] = 1.0
    
    if reverse: 
        latitude[g] = n.arcsin(b2)*radeg
        longitude[g] = n.arctan2(b1, b0)*radeg
        a, b =  longitude, latitude
    else:
        theta = n.arcsin(b2)*radeg
        phi = n.arctan2(b1, b0)*radeg     
        a, b = phi, theta
    
    return a, b 

def wcssph2xy(longitude, latitude, map_type=0, ctype="", face="", pv2="", crval="", crxy="", longpole="", latpole="", north_offset="", south_offset=""):
    
    
    return x, y, badindex
    
if __name__ == '__main__':
  # The following two lines will test all functions in the module.
  # Run `python astrom.py -v` to see verbose output.  It's probably best to rely
  # on this kind of testing once you believe the code to already be functional
  # since it is a great method for building but difficult to work with during
  # development.
  import doctest
  doctest.testmod()
    
    