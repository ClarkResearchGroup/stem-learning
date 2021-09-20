#
#   ftiff2py.py
#
#  read floatTIFF images using PIL 
#   (tested with python 2.7)
#
#------------------------------------------------------------------------
#Copyright 2018 Earl J. Kirkland
#
#This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#---------------------- NO WARRANTY ------------------
#THIS PROGRAM IS PROVIDED AS-IS WITH ABSOLUTELY NO WARRANTY
#OR GUARANTEE OF ANY KIND, EITHER EXPRESSED OR IMPLIED,
#INCLUDING BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
#MERCHANABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
#IN NO EVENT SHALL THE AUTHOR BE LIABLE
#FOR DAMAGES RESULTING FROM THE USE OR INABILITY TO USE THIS
#PROGRAM (INCLUDING BUT NOT LIMITED TO LOSS OF DATA OR DATA
#BEING RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR
#THIRD PARTIES OR A FAILURE OF THE PROGRAM TO OPERATE WITH
#ANY OTHER PROGRAM). 
#-----------------------------------------------------------------------------
#
#  started 26-may-2018 E. Kirkland
#  last updated 26-may-2018 ejk
#
# uses numpy, scipy, matplotlib packages
import numpy as np
from PIL import Image
from pylab import *

#  a few parameter offsets from slicelib.hpp
pNPIX    =   0   # number of pix 1 for real and 2 for complex
pXBTILT  =   5   # x beam tilt in rad
pYBTILT  =   6   # y beam tilt in rad
pDEFOCUS =   11  # defocus in Angstroms
pASTIG   =   12  # astigmatism in Angstroms
pTHETA   =   13  # angle of astigmatism in radians
pDX      =   14  # dimension of pixel in x direction in Angstroms
pDY      =   15  # dimension of pixel in y direction in Angstroms
pENERGY  =   16  # beam energy in keV
pOAPERT  =   17  # objective aperture semi-angle in radians
                 #    also see pOAPMIN below
pCS      =   18  # spherical aberration in Angstroms
pWAVEL   =   19  # electron wavelength in Angstroms
pNX      =   104 # (int) main image size (transmission function)
pNY      =   105

#---------- read_fTIFF() -------------------------------
#  read a floatTIFF 32 bit floating point image
#  remember that python reads in matrix order (row,col) not image (x,y)
def read_fTIFF( filename ):
    img = Image.open( filename, mode='r' )
    img.seek( 1 )    #  32 bit floating point image
    pix = np.array(img)
    img.seek( 2 )    #  32 bit floating point parameters
    param = transpose(np.array(img))  #  convert to 1D array
    if( 2 == param[pNPIX] ):   # a complex image
        ny = pix.shape[0]
        nx = pix.shape[1]/2
        repix = pix[:,0:nx] # re, im parts are side by side
        impix = pix[:,nx:2*nx]
        cpix = repix + 1j * impix
        return cpix, param
    else:
        return pix, param

infile = "si3n4.tif"   # a real image
#infile = "siauto.tif"  # a complex image
pix, param = read_fTIFF( infile )

# ------- show the image ----------
nx = pix.shape[1]    #  = param[pNX] = image size in pixels
ny = pix.shape[0]    #  = param[pNY]
dx = param[pDX]   #  pixel size (in Ang.)
dy = param[pDY]
print "image= ",infile, "nx, ny = ", nx, ", ", ny
print "pixel size = ", dx, ", ", dy

img = imshow( real(pix),  extent=(0,nx*dx,0,ny*dy) )
#img = imshow( imag(pix),  extent=(0,nx*dx,0,ny*dy) )
img.set_cmap( 'gray' )
colorbar()
savefig('test2.png')
show()
