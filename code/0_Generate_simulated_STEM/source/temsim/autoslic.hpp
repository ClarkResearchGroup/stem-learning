/*
    *** autoslic.hpp ***        (normal C++)
    *** autoslic_cuda.hpp ***   (CUDA code version)

  NOTE: the CUDA (GPU) version does the bandwith limit in a different
        order and will get a slightly different numerical result.

------------------------------------------------------------------------
Copyright 1998-2019 Earl J. Kirkland

This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

---------------------- NO WARRANTY ------------------
THIS PROGRAM IS PROVIDED AS-IS WITH ABSOLUTELY NO WARRANTY
OR GUARANTEE OF ANY KIND, EITHER EXPRESSED OR IMPLIED,
INCLUDING BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
IN NO EVENT SHALL THE AUTHOR BE LIABLE
FOR DAMAGES RESULTING FROM THE USE OR INABILITY TO USE THIS
PROGRAM (INCLUDING BUT NOT LIMITED TO LOSS OF DATA OR DATA
BEING RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR
THIRD PARTIES OR A FAILURE OF THE PROGRAM TO OPERATE WITH
ANY OTHER PROGRAM). 
------------------------------------------------------------------------

  ANSI C and TIFF version
  this version uses FFTW 3 (net about a factor of 2X faster)

  FFTW choses an optimum form of the FFT at run time so there
  is some variation in execution speed depending on what else 
  the CPU is doing during this planning stage

  see:   www.fftw.org

  on Windows file libfftw3f-3.dll must be in the PATH

  on Linux build as:
  g++ -O -fopenmp -o autoslic autoslic.cpp autosliccmd.cpp slicelib.o
                       tiffsubs.o  cfpix.o -lfftw3f

  Transmit an electron wave through a specimen using the
  multislce method with automatic slicing.  Read in the (x,y,z)
  coordinates of the whole specimen and break into slices
  on-the-fly.

  started 24-july-1996 E. Kirkland
  working 19feb-1997 ejk
  added look-up-table vzatomLUT() for 3X-4X increase 
        in speed 23-may-1997 ejk
  put bandwith limit inside trlayer() 1-oct-1997 ejk
  added Gaussian thermal displacements 1-oct-1997 ejk
  removed /sqrt(3) from Thermal rms displacements 
    to be consistent with Int'l X-ray tables 22-dec-1997 ejk
  corrected zmin/max error with thermal displac. 24-dec-1997 ejk
  fixed small aliasing problem 5-jan-1998 ejk
  added unit cell replication option and moved ReadXYZcoord()
    into slicelib.c  11-jan-1998 ejk
  added astigmatism and modify to use different set of
    random offsets on each illum. angle with partial coherence
         5-feb-1998 ejk
  fix typo in z range message with partial coherence and
    thermal vibrations 9-july-1998 ejk
  update memory allocation routines 19-nov-1999 ejk
  change void main() to int main() for better portability
         22-jan-2000 ejk
  fixed bug in zmin/zmax calculation in coherent mode
     (move to after sortByZ() - it was before ) 8-jan-2002 ejk
  add cross section option (in non-partial coherence mode only)
        27-may-2005 ejk
  convet to faster sortByZ() 8-feb-2006 ejk
  move sortbyz() to slicelib.c 5-sep-2006 ejk
  add echo on y position in pixels for xz mode 4-may-2007 ejk
  update data type of nxl,nyl to be consistent with new tiffsubs
     17-jul-2007 ejk
  move xz depthpix save to be after transmit+propagate to get a
     full slice and proper anti-aliasing and also be consisten
     with what you get doing it by hand  and increase possible
     slices output (nz was off by one) 24-jan-2008 ejk
  change propagation range to be whole unit cell not just
     range of atoms to treat sparsely populated spec.
     better (consistent with autostem) 23-mar-2008 ejk
  take small things out of loop in trlayer() 14-may-2008 ejk
  parameterize vzatomLUT() vs r^2 instead of r to avoid a lot of sqrt()
      calls (a little faster)  6-jun-2008 ejk
  move vzatomLUT() to slicelib.c  11-jun-2008 ejk
  convert to GPL 3-jul-2008 ejk
  add Cs5 (and Cs->Cs3) 15-dec-2009 ejk
  get return value of scanf() to remove warnings from gcc 4.4
     and convert to 4 char TAB size formatting 21-feb-2010 ejk
  add parallel computing of a few parts 21-feb-2010 ejk
  start conversion to faster FFTW 24-feb-2010 ejk
  move some things into slicelibW.c to share 6-mar-2010 ejk
  fix sign convention in FFTW 21-mar-2010 ejk
  update comments 4-apr-2010 ejk
  add option to average over many frozen phonon
      configurations 3-aug-2010 ejk
  add multipole aberrations to probe 12-may-2011 ejk
  start conversion to floatTIFF.cpp and C++ 28-may-2012 ejk
  working 3-jun-2012 ejk
  convert to cfpix/fftw class from raw fftw 13-nov-2012 to 21-nov-2012 ejk
  move calculation into a class with separate command line front end
      29-may-2013 ejk
  convert to string message 9-sep-2013 ejk
  change RNG seed argument to referenece so it get updated for 
      successive calls 21-sep-2013 ejk
  move toString() to slicelib from here 28-nov-2013 ejk
  add abbPhase2D() to calculate the 2D phase abb function 21-aug-2014 ejk 
  convert malloc1D() etc to vector<> 28-jun-2016 ejk
  split into separate simple calculation with openMP multithreading
     over phonon configurations 10-aug-2016 ejk
  add calculateCBED_TDS()  16-aug-2016 ejk
  add verbose()+echo to turn off unneccesary
        terminal echo 13-jan-2017 ejk
  convert max angle to obj. apert. (instead of 2/3 max)
      in abbPhase2D() 20-oct-2017 ejk
  normalize CBED to nwobble for easier comparison
      2-sep-2018 ejk
  delete unused variable and start GPU/cuda mode 8-oct-2018 ejk
  GPU/cuda code working 20-dec-2018 ejk
  add record of z coord. for each beam values 16-mar-2019 ejk
  add abbError() (copy from autostem) 11-aug-2019 ejk
  add last few changes back into cuda code 29-aug-2019 ekj

  ax,by,cz  = unit cell size in x,y
  BW     = Antialiasing bandwidth limit factor
  acmin  = minimum illumination angle
  acmax  = maximum illumination angle
  Cs     = spherical aberration coefficient
  df0    = defocus (mean value)
  sgmaf = defocus spread (standard deviation)
  dfdelt = sampling interval for defocus integration
  
  this file is formatted for a TAB size of 4 characters 
  
*/

#ifndef AUTOSLIC_HPP   // only include this file if its not already
#define AUTOSLIC_HPP   // remember that this has been included


#include <cstdio>  /* ANSI C libraries */
#include <cstdlib>
#include <cstring>
#include <cmath>

#include <string>   // STD string class
#include <vector>

using namespace std;

#include "cfpix.hpp"        // complex image handler with FFT
#include "slicelib.hpp"     // misc. routines for multislice
#include "probe.hpp"        //  for CBED

//#define ASL_USE_CUDA    // define to use nvidia cuda

#ifdef ASL_USE_CUDA
#include <cuda.h>
#include <cuda_device_runtime_api.h>
#include <cuda_runtime.h>
#include <cufft.h>
#else
#define USE_OPENMP  /* define to use openMP */
#endif


//------------------------------------------------------------------
class autoslic{

public:
    
    autoslic( );         // constructor functions
    
    ~autoslic();        //  destructor function

    // (input) control flags- should be set externally
    //   things that don't fit in param[]
    int lbeams, lcross, lpartl, lstart, lwobble, lcbed;

    int nillum;   //  (output) number of illumination angles used

    //  add random aberration tuning pi/4 errors for 2nd through 5th order
    void abbError( vector<float> &p1, int np, int NPARAM, 
        unsigned long &iseed, int echo, double scale );

    //  calculate aberation phase errors 
    void abbPhase2D( cfpix &pix, vectorf &param, int multiMode );

    //  main calculation
    int calculate(cfpix &pix, cfpix &wave0, cfpix &depthpix,
        vectorf &param, int multiMode, int natom,
        vectori &Znum, vectorf &x, vectorf &y, vectorf &z, vectorf &occ,
        cfpix &beams, vectori &hb, vectori &kb, int nbout, float ycross, int verbose   );

    //  Partial Coherent image calculation
    int calculatePartial(cfpix &pix,
        vectorf &param, int multiMode, int natom, unsigned long *iseed,
        vectori &Znum, vectorf &x, vectorf &y, vectorf &z, vectorf &occ, vectorf &wobble,
        float dfdelt   );

    //  CBED calculation with frozen phonons
    int calculateCBED_TDS(cfpix &pix,
        vectorf &param, int multiMode, int natom, unsigned long *iseed,
        vectori &Znum, vectorf &x, vectorf &y, vectorf &z, vectorf &occ,
        vectorf &wobble );

    //  separate initialize misc. variables 
    //  must be called before calculate()
    void initAS( vectorf &param, vectori Znum, int natom );

    void quiet() { echo = 0 ; };          // turn off unneeded terminal output
    void verbose() { echo = 1; };         // turn on unneeded terminal output

private:

        int NZMAX, echo;
        float BW, pi, k2max;
        double twopi, ABERR;

        vectorf kx, ky, xpos, ypos, kx2, ky2;

        cfpix cprop;           // complex propagator in Fourier space
        
        void trlayer(   const vectorf &x, const vectorf &y, const vectorf &occ,
            const vectori &Znum, const int natom, const int istart,
            const float ax, const float by,
            const float kev, cfpix &trans, const int nx, const int ny,
            const vectorf &kx2, const vectorf &ky2,
            double *phirms, int *nbeams, const float k2max   );

        std::string sbuffer;

        void messageAS( std::string &smsg,  int level = 0 );  // common error message handler

#ifdef ASL_USE_CUDA
        //  D prefix = on device, and H prefix = on Host
        float *Hcbed, *Dcbed, *Dphimin, *Dphimax;
        float *Dkx, *Dky, *Dkx2, *Dky2;
        float *Dspec, *Hspec, *DpotnR, *HpotnR;
        double *Hsums, *Dsums, *Hfparams, *Dfparams;
        cufftComplex *Dwave, *Hwave, *Dprop, *Hprop, *Htrans, *Dtrans, *DpotnC;
        cufftHandle cuplan, cuplanTc2r;
        int checkCudaErr( const char msg[] );
        double k2maxt;  //  really only used for cuda version
#endif


}; // end autoslic::

#endif
