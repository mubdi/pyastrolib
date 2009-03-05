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

import numpy
import scipy.ndimage
import pyfits



def aper(image, xc, yc, mags, errap, sky, skyerr, phpadu, apr, skyrad, badpix=None, SETSKYVAL=None, PRINT=True, SILENT=True, FLUX=None, EXACT=None, Nan=numpy.nan, READNOISE=None, MEANBACK=None, CLIPSIG=None, MAXITER=None, CONVERGE_NUM=None):
    """
    """
    dimx,dimy = image.shape # Number of columns and rows in image array
    naper = apr.size # Number of apertures
    nstars = xc.size # Number of stars
    area = PI*apr**2 # Area of each aperture

    bigrad = apr + 0.5 # Big radius for subpixel approximation
    smallrad = apr/numpy.sqrt(2.) - 0.5 # Small radius for subpixel approximation

    rinsq =  skyrad[0]**2 # inner sky radius squared
    routsq = skyrad[1]**2 # outer sky radius squared
    
    ###
    # Compute the limits of the submatrix.   Do all stars in vector notation.
    ###
    lx = (xc-skyrad[1]).astype(int) # Lower limit X direction
    ux = (xc+skyrad[1]).astype(int) # Upper limit X direction
    nx = ux-lx+1 # Number of pixels X direction
    ly = (yc-skyrad[1]).astype(int) # Lower limit Y direction
    uy = (yc+skyrad[1]).astype(int) # Upper limit Y direction
    ny = uy-ly +1 # Number of pixels Y direction
    dx = xc-lx # X coordinate of star's centroid in subarray
    dy = yc-ly # Y coordinate of star's centroid in subarray

    badstar = (lx < 0) * (ux > dimx) * (ly < 0) * (uy > dimy) # Stars too close to the edge
    
    nbox = numpy.ceil(skyrad)
    
    temp = numpy.arange(nbox)-middle # make 1D indices of the kernel box
    temp = numpy.resize(temp, (nbox,nbox)) # make 2D indices of the kernel box (for x or y)
    offsetsx = numpy.resize(temp, (nfound,nbox,nbox)) # make 2D indices of the kernel box for x, repeated nfound times
    offsetsy = numpy.resize(temp.T, (nfound,nbox,nbox)) # make 2D indices of the kernel box for y, repeated nfound times
    offsetsx = (indx + offsetsx.swapaxes(0,-1)).swapaxes(0,-1) # make it relative to image coordinate
    offsetsy = (indy + offsetsy.swapaxes(0,-1)).swapaxes(0,-1) # make it relative to image coordinate
    offsets_vals = image[offsetsx,offsetsy] # a (nfound, nbox, nbox) array of values (i.e. the kernel box values for each nfound candidate)

    return None


def find(image, hmin, fwhm, roundlim=[-1.,1.], sharplim=[0.2,1.]):
    """find(image, hmin, fwhm, roundlim=[-1.,1.], sharplim=[0.2,1.])
    Identifies stars in an image.
    Returns a tupple (x, y, flux, sharpness, roundness).
    
    image: 2D array containing the image
    hmin: Minimum threshold for detection. Should be 3-4 sigma above background RMS.
    fwhm: FWHM to be used for the convolution filter. Should be the same as the PSF FWHM.
    roundlim: Threshold for the roundness criteria.
    sharplim: Threshold for the sharpness criteria.
    
    Note: Pyfits imports images with x and y inverted with respect to IDL's convention.
    
    >>> import pyfits, pyastrolib
    >>> image = pyfits.getdata('test.fits')
    >>> (x, y, flux, sharpness, roundness) = pyastrolib.find(image, 15, 5.)
    """
    ###
    # Setting the convolution kernel
    ###
    maxbox = 13 # Maximum size of convolution box in pixels
   
    n_x = image.shape[0] # x dimension
    n_y = image.shape[1] # y dimension
    
    sigmatofwhm = 2*numpy.sqrt(2*numpy.log(2))
    radius = 1.5 * fwhm / sigmatofwhm # Radius is 1.5 sigma
    sigsq = (fwhm/sigmatofwhm)**2 # sigma squared
    nhalf = int(radius) # Center of the kernel
    nbox = 2*nhalf+1 # Number of pixels inside of convolution box
    middle = nhalf # Index of central pixel
    
    lastro = n_x - nhalf
    lastcl = n_y - nhalf
    
    kern_x, kern_y = numpy.ix_(numpy.arange(nbox),numpy.arange(nbox)) # x,y coordinates of the kernel
    g = (kern_x-nhalf)**2+(kern_y-nhalf)**2 # Compute the square of the distance to the center
    mask = g <= radius**2 # We make a mask to select the inner circle of radius "radius"
    nmask = mask.sum()
    g = numpy.exp(-0.5*g/sigsq) # We make the 2D gaussian profile
    
    ###
    # Convolving the image with a kernel representing a gaussian (which is assumed to be the psf)
    ###
    c = g*mask # For the kernel, values further than "radius" are equal to zero
    c[mask] = (c[mask] - c[mask].mean())/(c[mask].var() * nmask) # We normalize the gaussian kernel
    
    c1 = g[nhalf] # c1 will be used to the test the roundness
    sumc1 = c1.mean()
    sumc1sq = (c1**2).sum() - sumc1
    c1 = (c1-c1.mean())/((c1**2).sum() - c1.mean())
    
    h = scipy.ndimage.convolve(image,c,mode='constant',cval=0.0) # Convolve image with kernel "c"
    h[:nhalf,:] = 0 # Set the sides to zero in order to avoid border effects
    h[-nhalf:,:] = 0
    h[:,:nhalf] = 0
    h[:,-nhalf:] = 0
    
    mask[middle,middle] = False # From now on we exclude the central pixel
    nmask = mask.sum() # so the number of valid pixels is reduced by 1
    goodx,goody = mask.nonzero() # "good" identifies position of valid pixels
    
    ###
    # Identifying the point source candidates that stand above the background
    ###
    indx,indy = (h >= hmin).nonzero() # we identify point that are above the threshold, image coordinate
    nfound = indx.size # nfound is the number of candidates
    
    offsetsx = numpy.resize(goodx-middle,(nfound,nmask)) # a (nfound, nmask) array of good positions in the mask, mask coordinate
    offsetsx = indx + offsetsx.T # a (nmask, nfound) array of positions in the mask for each candidate, image coordinate
    offsetsy = numpy.resize(goody-middle,(nfound,nmask)) # a (nfound, nmask) array of good positions in the mask, mask coordinate
    offsetsy = indy + offsetsy.T # a (nmask, nfound) array of positions in the mask for each candidate, image coordinate
    offsets_vals = h[offsetsx,offsetsy] # a (nmask, nfound) array of mask values around each candidate
    vals = h[indx,indy] # a (nfound) array of the intensity of each candidate
    
    ###
    # Identifying the candidates are local maxima
    ###
    ind_goodcandidates = ((vals - offsets_vals) > 0).all(axis=0) # a (nfound) array identifying the candidates whose values are above the mask (i.e. neighboring) pixels, candidate coordinate
    nfound = ind_goodcandidates.sum() # update the number of candidates
    indx = indx[ind_goodcandidates] # a (nfound) array of x indices of good candidates, image coordinate
    indy = indy[ind_goodcandidates] # a (nfound) array of y indices of good candidates, image coordinate
    
    ###
    # Identifying the candidates that meet the sharpness criteria
    ###
    d = h[indx,indy] # a (nfound) array of the intensity of good candidates
    d_image = image[indx,indy] # a (nfound) array of the intensity of good candidates in the original image (before convolution)
    offsetsx = offsetsx[:,ind_goodcandidates] # a (nmask, nfound) array of positions in the mask for each candidate, image coordinate
    offsetsy = offsetsy[:,ind_goodcandidates] # a (nmask, nfound) array of positions in the mask for each candidate, image coordinate
    offsets_vals = image[offsetsx,offsetsy]
    sharp1 = (d_image - offsets_vals.sum(0)/nmask) / d
    ind_goodcandidates = (sharp1 > sharplim[0]) * (sharp1 < sharplim[1]) # a (nfound) array of candidates that meet the sharpness criteria
    nfound = ind_goodcandidates.sum() # update the number of candidates
    indx = indx[ind_goodcandidates] # a (nfound) array of x indices of good candidates, image coordinate
    indy = indy[ind_goodcandidates] # a (nfound) array of y indices of good candidates, image coordinate
    sharp1 = sharp1[ind_goodcandidates] # update sharp1 with the good candidates
    
    ###
    # Identifying the candidates that meet the roundness criteria
    ###
    temp = numpy.arange(nbox)-middle # make 1D indices of the kernel box
    temp = numpy.resize(temp, (nbox,nbox)) # make 2D indices of the kernel box (for x or y)
    offsetsx = numpy.resize(temp, (nfound,nbox,nbox)) # make 2D indices of the kernel box for x, repeated nfound times
    offsetsy = numpy.resize(temp.T, (nfound,nbox,nbox)) # make 2D indices of the kernel box for y, repeated nfound times
    offsetsx = (indx + offsetsx.swapaxes(0,-1)).swapaxes(0,-1) # make it relative to image coordinate
    offsetsy = (indy + offsetsy.swapaxes(0,-1)).swapaxes(0,-1) # make it relative to image coordinate
    offsets_vals = image[offsetsx,offsetsy] # a (nfound, nbox, nbox) array of values (i.e. the kernel box values for each nfound candidate)
    dx = (offsets_vals.sum(1)*c1).sum(1)
    dy = (offsets_vals.sum(2)*c1).sum(1)
    around = 2*(dx-dy)/(dx+dy)
    ind_goodcandidates = (around > roundlim[0]) * (around < roundlim[1]) * (dx >= 0.) * (dy >= 0.) # a (nfound) array of candidates that meet the roundness criteria
    nfound = ind_goodcandidates.sum() # update the number of candidates
    indx = indx[ind_goodcandidates] # a (nfound) array of x indices of good candidates, image coordinate
    indy = indy[ind_goodcandidates] # a (nfound) array of y indices of good candidates, image coordinate
    sharp1 = sharp1[ind_goodcandidates] # update sharp1 with the good candidates
    around = around[ind_goodcandidates] # update around with the good candidates
    offsets_vals = offsets_vals[ind_goodcandidates] # update offsets_vals with good candidates
    
    ###
    # Compute the approximate flux
    ###
    wt = nhalf - abs(numpy.arange(nbox, dtype=float)-nhalf) + 1
    xwt = numpy.resize(wt, (nbox,nbox)).T
    ywt = xwt.T
    sgx = (g*xwt).sum(0)
    sgy = (g*ywt).sum(1)
    p = wt.sum()
    sumgx = (wt*sgy).sum()
    sumgy = (wt*sgx).sum()
    sumgsqx = (wt*sgx*sgx).sum()
    sumgsqy = (wt*sgy*sgy).sum()
    vec = nhalf - numpy.arange(nbox)
    dgdx = sgy*vec
    dgdy = sgx*vec
    sdgdxs = (wt*dgdx**2).sum()
    sdgdys = (wt*dgdy**2).sum()
    sdgdx = (wt*dgdx).sum()
    sdgdy = (wt*dgdy).sum()
    sgdgdx = (wt*sgy*dgdx).sum()
    sgdgdy = (wt*sgx*dgdy).sum()
    
    sd = (offsets_vals*ywt).sum(2)
    sumgd = (wt*sgy*sd).sum(1)
    sumd = (wt*sd).sum(1)
    sddgdx = (wt*sd*dgdx).sum(1)
    hx = (sumgd - sumgx*sumd/p) / (sumgsqy - sumgx**2/p)
    skylvl = (sumd - hx*sumgx)/p
    dx = (sgdgdx - (sddgdx-sdgdx*(hx*sumgx + skylvl*p)))/(hx*sdgdxs/sigsq)

    sd = (offsets_vals*xwt).sum(1)
    sumgd = (wt*sgx*sd).sum(1)
    sumd = (wt*sd).sum(1)
    sddgdy = (wt*sd*dgdy).sum(1)
    hy = (sumgd - sumgy*sumd/p) / (sumgsqx - sumgy**2/p)
    skylvl = (sumd - hy*sumgy)/p
    dy = (sgdgdy - (sddgdy-sdgdy*(hy*sumgy + skylvl*p)))/(hy*sdgdys/sigsq)

    x = indx + dx
    y = indy + dy
    flux = h[indx,indy]
    
    return x, y, flux, sharp1, around





if __name__ == '__main__':
  # The following two lines will test all functions in the module.
  # Run `python astro.py -v` to see verbose output.  It's probably best to rely
  # on this kind of testing once you believe the code to already be functional
  # since it is a great method for building but difficult to work with durring
  # development.
  import doctest
  doctest.testmod()


