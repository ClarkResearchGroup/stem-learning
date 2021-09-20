/*
    *** autoslic.cpp ***       (normal C++)
    *** autoslic_cuda.cu ***   (CUDA code version)

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

  ANSI C++ and TIFF version
  this version uses FFTW 3 (net about a factor of 2X faster)

  There are two multithreaded (openMP) modes.  One is coarse grain
  (multithread-1) and multithreaded over frozen phonon configurations
  which is much faster for frozen phonon calculation. The other mode
  is fine grain (multithread-2) and multithreaded over atoms in the
  specimen transmission function. This mode is a little faster for
  every mode. Only one can be uncommented (enabled) as both can't
  be used at the same time.

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
  fix typo in starting wave function 11-jun-2013 ejk
  fix bug to restore orginal df value, and add code to handle
     sigmf=0 problem 29-jun-2013 ejk
  fix minor format issue in nbeams message (%ld to %d) 19-jul-2013 ejk
  convert to string message 9-sep-2013 ejk
  change RNG seed argument to referenece so it get updated for 
      successive calls 21-sep-2013 ejk
  move toString() to slicelib from here 28-nov-2013 ejk
  add abbPhase2D() to calculate the 2D phase abb function 21-aug-2014 ejk 
  convert malloc1D() etc to vector<> 28-jun-2016 ejk
  split into separate simple calculation with openMP multithreading
     over phonon configurations 10-aug-2016 ejk
  add calculateCBED_TDS()  16,28-aug-2016 ejk
  add verbose()+echo to turn off unneccesary
        terminal echo 13-jan-2017 ejk
  add normal diffraction mode to calculate_CBED_TDS()
    using lcbed flag 12-jul-2017 ejk
  convert max angle to obj. apert. (instead of 2/3 max)
      in abbPhase2D() 20-oct-2017 ejk
  normalize CBED to nwobble for easier comparison
      2-sep-2018 ejk
  delete unused variable and start GPU/cuda mode 8-oct-2018 ejk
  add option to trlayer() to return just real potential if k2max<0
      (to match autostem.cpp) 25-oct-2018 ejk
  GPU/cuda code working 20-dec-2018 ejk
  add record of z coord. for each beam values 16-mar-2019 ejk
  fix output size in xz cross section (remove few extra lines)
     - has been here since C++ conversion (cuda not tested) 7-apr-2019 ejk
  add total intensity output inautoslic::calculateCBED_TDS 26-jun-2019 
  add abbError() (copy from autostem) 11-aug-2019 ejk
  add last few changes back into cuda code 29-aug-2019 ekj
  fix small bug in beams output in cuda mode 12-sep-2019 ejk
  fix normalization in diff mode in calculateCBED_TDS() 13-sep-2019

  ax,by,cz  = unit cell size in x,y()
  BW     = Antialiasing bandwidth limit factor
  acmin  = minimum illumination angle
  acmax  = maximum illumination angle
  Cs     = spherical aberration coefficient
  df0    = defocus (mean value)
  sgmaf = defocus spread (standard deviation)
  dfdelt = sampling interval for defocus integration
  
  this file is formatted for a TAB size of 4 characters 
  
*/

#include "slicelib.hpp"    // misc. routines for multislice
#include "cfpix.hpp"       // complex image handler with FFT

//
//  to select cuda/GPU code include autoslic_cuda.hpp
//  and the filename extention for this file must end in .cu
//  (.cpp for normal C++) - its better put other selections in .hpp so
//  the calling program has more information - use only one .hpp
//
//#include "autoslic_cuda.hpp"    // header for this class
#include "autoslic.hpp"    // header for this class

#include <sstream>  // string streams

//   set up cuda things
#ifdef ASL_USE_CUDA
#include "cudaSlice.cu"
//
#define AST_USE_CUATOMPOT  // to use cuAtompot() - usually do NOT do this
#ifdef AST_USE_CUATOMPOT
extern vector< vector<double> > fparams;  // scattering factor table
#endif
#else
#define USE_OPENMP      // define to use openMP 
#endif


//=============================================================
//---------------  creator and destructor --------------

autoslic::autoslic()
{
        BW= (2.0F/3.0F);  // bandwidth limit
        ABERR= 1.0e-4;    // max error for a,b

        NZMAX= 103;   // max atomic number Z

        twopi = 2.0 * (4.0 * atan( 1.0 ));

        // init control flags
        lcross = 0;
        lpartl = 0;
        lstart = 0;
        lwobble = 0;
        lbeams = 0;
        lcbed = 1;  // default to original CBED

        echo = 1;   // >0 to echo status 

        pi = (float) (4.0 * atan( 1.0 ));

        return;

};   //  end autoslic::autoslic()

autoslic::~autoslic()
{
}

//=============================================================
//   same subroutine as in autostem.cpp
//
/*--------------------- abbError() ----------------------------*/
/*
  add random aberration tuning pi/4 errors for 2nd through 5th order

  p1 = input parameter array (modified)
  np = number of params to change (np=9 for order 2,3 and np=22 for order 2,3,4,5)
  NPARAM = total number of parameters
  iseed = the current randim number seed
  echo = echo diagnostics if == 1
  scale = scale error as scale * pi/4 (usually should be 1.0)
          (for testing tuning accuracy sensitivity)

  modified from abbError() in mcprobea.cpp

*/
void autoslic::abbError( vector<float> &p1, int np, int NPARAM, 
    unsigned long &iseed, int echo, double scale )
{
    int i;
    double dchiMax[6], aTune, at3, wavl, x;

    const int paramToVary[]={ pC21a, pC21b, pC23a, pC23b,
                        pCS, pC32a, pC32b, pC34a, pC34b,
                        pC41a, pC41b, pC43a, pC43b, pC45a, pC45b,
                        pCS5, pC52a, pC52b, pC54a, pC54b, pC56a, pC56b};

    const int orderToVary[]={2,2,2,2, 3,3,3,3,3, 4,4,4,4,4,4, 5,5,5,5,5,5,5};

    if( (np < 1 ) || (np>22) || (np>NPARAM) ) return;

    wavl = p1[pWAVEL];
    aTune = p1[pOAPERT];

    //------ calculate max aberr. error for scale*pi/4 phase at max obj. angle
    at3 = aTune * aTune * aTune;
    dchiMax[2] = scale * (3.0*wavl)/( 8.0 * at3 );          // 2nd order
    dchiMax[3] = scale * (4.0*wavl)/( 8.0 * at3 * aTune);   // 3rd order
    dchiMax[4] = scale * (5.0*wavl)/( 8.0 * at3 * aTune * aTune);   // 4th order
    dchiMax[5] = scale * (6.0*wavl)/( 8.0 * at3 * at3 );    // 5th order
    if( 1 == echo ) {
        sbuffer = "max. aberr. error for order 2,3,4,5 = \n"+toString(dchiMax[2])
            + ", " +toString(dchiMax[3]) + ", " + toString(dchiMax[4]) 
            + ", " +toString(dchiMax[5])+ " Ang.";
        //messageAST( sbuffer, 0 );
        messageAS( sbuffer, 0 );
    }
 
    for( i=0; i<np; i++) {
        x = (float) ( 2.0*ranflat(&iseed) - 1.0);  //  range = -1 to +1
        p1[paramToVary[i]] += (float) ( x  * dchiMax[orderToVary[i]] );
    }

    return;

}  //  end abbError()
//=============================================================
/*  abbPhase2D()

  calculate the phase of the abberation function in 2D
  in the objective aperture plane 
  (but out to the max angle allowed by sampling)
  - mainly just to look at

  added 21-aug-2014 ejk
  convert max angle to obj. apert. (instead of 2/3 max)
      20-oct-2017 ejk

  ab2D() = will get the 2D image of the abb. phase
  param[] = holds image parameters
  multiMode = flag, if not 0 then include all multipole abberations

*/
void autoslic::abbPhase2D( cfpix &ab2D,vectorf &param, int multiMode )
{
    int ix, iy, ixmid, iymid, nx, ny;
    float k2, k2max, v0, wavlen, ax, by, t, aobj;
    double chi0, alx, aly;

    vectorf kx, ky, xpos, ypos, kx2, ky2;

    // ---- get setup parameters from param[]
    ax = param[ pAX ];
    by = param[ pBY ];
    nx = ToInt( param[ pNX ] );
    ny = ToInt( param[ pNY ] );
    v0 = param[pENERGY];                // electron beam energy in keV
    aobj = param[ pOAPERT ];

    wavlen = (float) wavelength( v0 );

    //----- calculate spatial frequencies and positions for future use 

    kx.resize( nx );
    kx2.resize( nx );
    xpos.resize( nx );
    freqn( kx, kx2, xpos, nx, ax );

    ky.resize( ny );
    ky2.resize( ny );
    ypos.resize( ny );
    freqn( ky, ky2, ypos, ny, by );

    //  rearrange frequencies before calculation so we don't
    //  have to rearrange the whole image
    ixmid = nx/2;
    iymid = ny/2;
    for( ix=0; ix<ixmid; ix++) {
        t=kx[ix];  kx[ix] = kx[ix+ixmid];  kx[ix+ixmid] =t;
        t=kx2[ix]; kx2[ix]= kx2[ix+ixmid]; kx2[ix+ixmid]=t;
    }
    for( iy=0; iy<iymid; iy++) {
        t=ky[iy];  ky[iy] = ky[iy+iymid];  ky[iy+iymid] =t;
        t=ky2[iy]; ky2[iy]= ky2[iy+iymid]; ky2[iy+iymid]=t;
    }

    //-----  make array right size if needed
    ab2D.resize( nx, ny );
    //ab2D = 0.0F;   //  grey background
    ab2D = pi;    // white background

    //-----  calculate max angles
    k2max = aobj/ wavlen; 
    k2max = k2max*k2max;
    
    //------- calcualte the phase and file the array
    for( ix=0; ix<nx; ix++) {
        alx = wavlen * kx[ix];  // x component of angle alpha
        for( iy=0; iy<ny; iy++) {
            aly = wavlen * ky[iy];  // y component of angle alpha
            k2 = kx2[ix] + ky2[iy];
            if( k2 <= k2max ) {
                chi0 = (2.0*pi/wavlen) * chi( param, 
                               alx, aly, multiMode );
                //  make phase modulo 2pi - should be a better way to do this (?)
                //   remember that % operator only works with int
                while( chi0 < -pi) chi0 += twopi;
                while( chi0 >  pi) chi0 -= twopi;

                ab2D.re(ix, iy ) = (float) chi0;
            }
        } // end for( iy=...
    }  // end for(ix=....

    //----------- end --------------------

    return;
 
}  // end abbPhase2D()

//=============================================================
/*  calculate()

calculate one coherent pass of incident wave through the specimen

keep this as simple as possible and thread safe (run in parallel if possible)
loop over this for more complicated things

NOTE: FFTW measure is NOT thread safe so cfpix.init() which calls fftw-measure
cannot be in this routine (init cfpix before calling this subroutine)
- calling init/measure on same size pix in every thread is a bad idea anyway

initAS() must be called once before calling this and pix
must have been .init()'d for FFT before calling this subroutine

  input:
        wave0 = complex image with starting image, left unchanged
                (ignored if lstart = 0 or partial coherence calculated)
        param[] = image parameters (most will not be changed here)
        multimode = flag controlling multipole aberrations
        natom = number of atoms
        x[],y[],z[] = atomic coord (will be sorted)
        Znum[] atomic number of each atom (will be sorted)
        occ[] = occupancy of each atomic site (will be sorted)
        wobble[] = thermal oscillation amplitude  (will be sorted)
        hb[],kb[] = (h,k) indexes of beams to monitor during
                          propagation
                (ignore if lbeams = 0)
        nbout  = number of beams to record
        ycross = y position to save xz cross section
        verbose = if >0 then print a lot otherwise keep mostly silent

  mode flags:  lbeams, lcross, lpartl, lstart, lwobble

  output:
        pix = complex image to get results, orig. data lost and may be resized
                (complex for coherent calc. and real for partial coherence)
                must be of right size anf have fft init()
        beams = track specified beams as wave propagates
                (ignored if lbeams = 0)
*       depthpix = will be resized and get xz cross section image
                save intensity in real part (imag part not used)
                (ignored if lcross = 0)

    return value:  >=0 for success and <0 for failure
*/
#ifndef ASL_USE_CUDA  //  plain non-cuda code
int autoslic::calculate(cfpix &pix, cfpix &wave0, cfpix &depthpix,
        vectorf &param, int multiMode, int natom,
        vectori &Znum, vectorf &x, vectorf &y, vectorf &z, vectorf &occ,
        cfpix &beams, vectori &hb, vectori &kb, int nbout, float ycross, int verbose )
{
    int i, ix, iy, iz, nx, ny, nz, iycross, istart, nbeams,
        ib, na, islice, nzbeams, nzout;

    float xmin,xmax, ymin, ymax, zmin, zmax;
    float scale, v0, wavlen, ax, by, tctx;

    double sum, zslice, deltaz, phirms;

    string sbuf;   //  need local copy to run in parallel

    vectori Znum2, hbeam, kbeam;

    cfpix wave;            // complex probe wave functions
    cfpix trans;           // complex transmission functions
    cfpix depthpix0;       // temp depth pix to get size right

    // ---- get setup parameters from param[]
    ax = param[ pAX ];
    by = param[ pBY ];
    nx = ToInt( param[ pNX ] );
    ny = ToInt( param[ pNY ] );
    v0 = param[pENERGY];                // electron beam energy in keV
    deltaz = param[ pDELTAZ ];          // slice thickness

    if( (nx < 1) || (ny < 1) || (ax<0.0) || (by<0.0) ){
        sbuf="bad size parameters in autoslic::calculate()";
        messageAS( sbuf );
        return( -1 );    //  bad situation ????
    }
            
    /*  calculate relativistic factor and electron wavelength */
    wavlen = (float) wavelength( v0 );
    //printf("electron wavelength = %g Angstroms\n", wavlen);

    /*  calculate the total specimen volume and echo */
    xmin = xmax = x[0];
    ymin = ymax = y[0];
    zmin = zmax = z[0];

    for( i=0; i<natom; i++) {
        if( x[i] < xmin ) xmin = x[i];
        if( x[i] > xmax ) xmax = x[i];
        if( y[i] < ymin ) ymin = y[i];
        if( y[i] > ymax ) ymax = y[i];
        if( z[i] < zmin ) zmin = z[i];
        if( z[i] > zmax ) zmax = z[i];
    }
    // --- leave this in main calling program
    //sprintf(stemp, "Total specimen range is\n %g to %g in x\n"
    //       " %g to %g in y\n %g to %g in z",
    //         xmin, xmax, ymin, ymax, zmin, zmax );
    //messageAS( stemp );

    /*---- allocate some more arrays and initialize wavefunction ----*/

    if( (nbout > 0) && (lbeams ==1) ) {
        hbeam.resize( nbout );
        kbeam.resize( nbout );
    }

    trans.resize( nx, ny );
    wave.resize( nx, ny );

    if( (lstart == 0) || (nx!= wave0.nx()) || (ny!=wave0.ny()) ) {
             wave = 1.0F;  
    } else { wave = wave0; }

    trans.copyInit( pix );  // remember that init() is NOT thread safe
    wave.copyInit( pix );   //  must be after "wave = wave0"

    if( lcross == 1 ) {
        /* nz may be too small with thermal vibrations so add a few extra */
        nz = (int) ( (zmax-zmin)/ deltaz + 3.5);
        depthpix0.resize( nx, nz );
        for( ix=0; ix<nx; ix++)
        for( iz=0; iz<nz; iz++)  depthpix0.re(ix,iz) = depthpix0.im(ix,iz) = 0.0F;
        iycross = (int) ( 0.5 + (ny * ycross / by));
        while( iycross < 0 ) iycross += ny;
        iycross = iycross%ny;  /* make periodic in ny */
        sbuf = "save xz cross section at iy= "+toString(iycross)+" pixels";
        messageAS( sbuf );
    }

    /*  ------  */
 
    k2max = nx/(2.0F*ax);
    tctx = ny/(2.0F*by);
    if( tctx < k2max ) k2max = tctx;
    k2max = BW * k2max;
    if( verbose > 0 ) {
        sbuf= "Bandwidth limited to a real space resolution of "+toString(1.0F/k2max)
                +" Angstroms";
        messageAS( sbuf );
        sbuf = "   (= " + toString(wavlen*k2max*1000.0F)
                + " mrad)  for symmetrical anti-aliasing.";
        messageAS( sbuf );
    }
    k2max = k2max*k2max;
 
/*  iterate the multislice algorithm proper

   NOTE: zero freg is in the bottom left corner and
     expands into all other corners - not in the center
     this is required for the FFT - don't waste time rearranging

*/

/* ---- start coherent method below ----------------
        (remember that wave was initialize above) */

    if( lbeams ==1 ) {
        nzbeams = (int) ( (zmax-zmin)/ deltaz + 3.5);
        beams.resize(nbout+1, nzbeams );   // to save values
        beams = 0;
        for(ib=0; ib<nbout; ib++) {
                hbeam[ib] = hb[ib];
                kbeam[ib] = kb[ib];
        }

        //  make them all positive just in case
        for( ib=0; ib<nbout; ib++) {
            if( hbeam[ib] < 0 ) hbeam[ib] = nx + hbeam[ib];
            if( kbeam[ib] < 0 ) kbeam[ib] = ny + kbeam[ib];
            if( hbeam[ib] < 0 ) hbeam[ib] = 0;
            if( kbeam[ib] < 0 ) kbeam[ib] = 0;
            if( hbeam[ib] > nx-1 ) hbeam[ib] = nx-1;
            if( kbeam[ib] > ny-1 ) kbeam[ib] = ny-1;
        }
    }  // end if( lbeams....

    if( verbose > 0 ) {
        sbuf = "Sorting atoms by depth...";
        messageAS( sbuf );
    }
    sortByZ( x, y, z, occ, Znum, natom );

    if( lwobble == 1 ){
        zmin = z[0];        /* reset zmin/max after wobble */
        zmax = z[natom-1];
        if( verbose > 0 ) {
            sbuf ="Thickness range "
                " is "+toString(zmin)+" to "+toString(zmax)+" (in z)";
            messageAS( sbuf );
        }
    }

    scale = 1.0F / ( ((float)nx) * ((float)ny) );

    zslice = 0.75*deltaz;  /*  start a little before top of unit cell */
    istart = 0;
    islice = 1;

    while( (istart < natom) && ( zslice < (zmax+deltaz) ) ) {

        /* find range of atoms for current slice */
        na = 0;
        for(i=istart; i<natom; i++) 
        if( z[i] < zslice ) na++; else break;

        /* calculate transmission function, skip if layer empty */
        if( na > 0 ) {
            trlayer( x, y, occ,
                Znum, na, istart, ax, by, v0, trans,
                nx, ny, kx2, ky2, &phirms, &nbeams, k2max );
   
            wave *= trans;    //  transmit
        }

        /*  bandwidth limit */
        wave.fft();

        if( (lbeams== 1) && (islice<nzbeams) && (islice>0) )  {
            for( ib=0; ib<nbout; ib++) {
                beams.re(ib,islice-1) = scale*wave.re(hbeam[ib],kbeam[ib] );   // real
                beams.im(ib,islice-1) = scale*wave.im(hbeam[ib],kbeam[ib] );   // imag
            }
            beams.re(nbout, islice-1) = (float) zslice;  //  save actual z coord.
            beams.im(nbout, islice-1) = (float) +1;  //  track actual beams- in case last is missed 
        }

        /* remember: prop needed here to get anti-aliasing right */
        wave *= cprop;
        wave.ifft();

        /* save depth cross section if requested */
        if( (lcross == 1) && (islice<=nz) ) {
            for( ix=0; ix<nx; ix++) {
                depthpix0.re(ix, islice-1) = 
                    wave.re(ix,iycross)*wave.re(ix,iycross)
                        + wave.im(ix,iycross)*wave.im(ix,iycross);
            }
            nzout = islice;
        }

        if( verbose > 0 ) {
            sum = 0.0;
            for( ix=0; ix<nx; ix++) {
                for( iy=0; iy<ny; iy++)
                    sum += wave.re(ix,iy)*wave.re(ix,iy) +
                        wave.im(ix,iy)*wave.im(ix,iy);
            }
            sum = sum * scale;

            sbuf = "z= " + toString(zslice)+" A, " + toString(nbeams) + " beams, "
                    + toString(na)+" coord., \n"
                    + "     aver. phase= "+toString(phirms)
                    +", total intensity = "+toString(sum) ;
            messageAS( sbuf );
        }   //  if( verbose > 0....

        zslice += deltaz;
        istart += na;
        islice++;

    } /* end while(istart<natom..) */

    pix.resize(nx,ny);
    pix = wave;

    //----------- end and exit --------------------

    // save depth cross section if requested and get rid of possible extra lines
    if( 1 == lcross ) {
        depthpix.resize( nx, nzout );  //  return to calling program
        for( ix=0; ix<nx; ix++)        //  depthpix0 will get deallocated
        for( iz=0; iz<nzout; iz++)  {
            depthpix.re(ix,iz) = depthpix0.re(ix,iz);
            depthpix.im(ix,iz) = depthpix0.im(ix,iz);
        }
    }

    return( +1 );

};   //  end autoslic::calculate() -plain C/C++ version

//---------------------------------------------------------------------------------------

#elif defined(ASL_USE_CUDA) //  cuda code
//  check last CUDA error
int autoslic::checkCudaErr( const char msg[] )
{
    if( cudaGetLastError() != cudaSuccess ) {
        std::string sbuffer = "cuda errror: " + std::string(msg);
            messageAS( sbuffer, 0 );
        exit(0 );  //  should do something better here (?)
    }
    return +1;
};

// NOTE:  cuda version is not thread-safe (nvcc doesn't support openMP either)
//        - so put all of init inside calculate() for now because its easier 
//           to manage - maybe reorg. in future

int autoslic::calculate(cfpix &pix, cfpix &wave0, cfpix &depthpix,
        vectorf &param, int multiMode, int natom,
        vectori &Znum, vectorf &x, vectorf &y, vectorf &z, vectorf &occ,
        cfpix &beams, vectori &hb, vectori &kb, int nbout, float ycross, int verbose )
{
    int i, ix, iy, iz, nx, ny, nz, iycross, istart, nbeams,
        ib, na, islice, nzbeams, nzout;

    float xmin,xmax, ymin, ymax, zmin, zmax;
    float scale, v0, wavlen, ax, by, tctx;

    double sum, zslice, deltaz, phirms;

    string sbuf;   //  need local copy to run in parallel

    vectori Znum2, hbeam, kbeam;

    cfpix wave;            // complex probe wave functions
    cfpix trans;           // complex transmission functions
    cfpix depthpix0;       // temp septh pix to get size right

    // ---- get setup parameters from param[]
    ax = param[ pAX ];
    by = param[ pBY ];
    nx = ToInt( param[ pNX ] );
    ny = ToInt( param[ pNY ] );
    v0 = param[pENERGY];          // electron beam energy in keV
    deltaz = param[ pDELTAZ ];          // slice thickness

    if( (nx < 1) || (ny < 1) || (ax<0.0) || (by<0.0) ){
        sbuf="bad size parameters in autoslic::calculate()";
        messageAS( sbuf );
        return( -1 );    //  bad situation ????
    }
            
    /*  calculate relativistic factor and electron wavelength */
    wavlen = (float) wavelength( v0 );
    //printf("electron wavelength = %g Angstroms\n", wavlen);

    /*  calculate the total specimen volume and echo */
    xmin = xmax = x[0];
    ymin = ymax = y[0];
    zmin = zmax = z[0];

    for( i=0; i<natom; i++) {
        if( x[i] < xmin ) xmin = x[i];
        if( x[i] > xmax ) xmax = x[i];
        if( y[i] < ymin ) ymin = y[i];
        if( y[i] > ymax ) ymax = y[i];
        if( z[i] < zmin ) zmin = z[i];
        if( z[i] > zmax ) zmax = z[i];
    }
    // --- leave this in main calling program
    //sprintf(stemp, "Total specimen range is\n %g to %g in x\n"
    //       " %g to %g in y\n %g to %g in z",
    //         xmin, xmax, ymin, ymax, zmin, zmax );
    //messageAS( stemp );

    /*---- allocate some more arrays and initialize wavefunction ----*/

    if( (nbout > 0) && (lbeams ==1) ) {
        hbeam.resize( nbout );
        kbeam.resize( nbout );
    }

    trans.resize( nx, ny );
    wave.resize( nx, ny );

    trans.copyInit( pix );  // remember that init() is NOT thread safe
    wave.copyInit( pix );   //  must be after "wave = wave0"

    if( (lstart == 0) || (nx!= wave0.nx()) || (ny!=wave0.ny()) ) {
             wave = 1.0F;  
    } else { wave = wave0; }

    //  make cuda plan on device for probe size and transmission size
    if( cufftPlan2d( &cuplan, nx, ny, CUFFT_C2C) != CUFFT_SUCCESS) {
    sbuffer = "cuda error; unable to create plan P\n" ; 
        messageAS( sbuffer, 0 );
        return( -1 );
    }

    //   allocate trans memory on the device 
    cudaMalloc( (void**)&Dtrans, sizeof(cufftComplex)*nx*ny);
    checkCudaErr( "failed to allocate cuda array Dtrans 1" );
    cudaMalloc( (void**)&DpotnR, sizeof(float)*nx*ny );
    checkCudaErr( "failed to allocate cuda array DpotnR" );
    HpotnR = new float[ nx*ny ];

    //  allocate host memory for transmission function
    Htrans = new cufftComplex[ nx*ny ];
    if( NULL == Htrans ) {
        sbuffer = "cannot allocate cuda host array Htrans\n";
        messageAS( sbuffer, 0 );
        return( -1 );
    }

    //  - have to do this here to know when to delete Dprop on device
    //  sometime do this in Hprop[] without cprop in initAS() (?)
    //  - really should not do this everytime but for now leave it here
    //  allocate propagator on host and device
    Hprop = new cufftComplex [nx*ny];
    cudaMalloc( (void**)&Dprop, sizeof(cufftComplex)*nx*ny );
    checkCudaErr( "failed to allocate device memory for Dprop");

    //  calculate 2D propagator function
    float scale2 = 1.0/(nx*ny);   // for FFT scaling
    for( ix=0; ix<nx; ix++) {
        for( iy=0; iy<ny; iy++) {
            Hprop[iy + ix*ny].x = cprop.re(ix,iy) * scale2;
            Hprop[iy + ix*ny].y = cprop.im(ix,iy) * scale2;

    } }  //  end for( iy=0... for( ix=0...
    //--- copy to device ---------
    cudaMemcpy( Dprop, Hprop, nx*ny * sizeof(cufftComplex),
             cudaMemcpyHostToDevice );
    checkCudaErr( "cannot copy prop from host to device");
    delete [] Hprop;

    //  allocate host memory for detector sums
    Hsums= new double[ ny ];
    if( NULL == Hsums ) {
        sbuffer = "cannot allocate cuda host array Hsums\n";
        messageAS( sbuffer, 0 );
        return( -1 );
    }
    cudaMalloc( (void**)&Dsums, sizeof(double)*ny );
    checkCudaErr( "failed to allocate device memory for Dsums");

    //   allocate CBED memory on host and device 
    Hcbed = new float[ nx*ny ];
    cudaMalloc( (void**)&Dcbed, sizeof(float)*nx*ny);
    checkCudaErr( "failed to allocate cuda array Dcbed" );

    //  allocate kx[],ky[] kx2[], ky2[] on device
    cudaMalloc( (void**)&Dkx, sizeof(float)*nx );
    checkCudaErr( "failed to allocate device memory for Dkx");
    cudaMalloc( (void**)&Dky, sizeof(float)*ny );
    checkCudaErr( "failed to allocate device memory for Dky");
    cudaMalloc( (void**)&Dkx2, sizeof(float)*nx );
    checkCudaErr( "failed to allocate device memory for Dkx2");
    cudaMalloc( (void**)&Dky2, sizeof(float)*ny );
    checkCudaErr( "failed to allocate device memory for Dky2");

    //  ------ calculate spatial freq on GPU
    int threadsPerBlock, blocksPerGrid;
    threadsPerBlock = 512;  // max thread/block usually 1024 =dev specific
    blocksPerGrid = (nx + threadsPerBlock -1 )/threadsPerBlock;  //  round up
    cuFreq<<<blocksPerGrid, threadsPerBlock>>>( Dkx, Dkx2, nx, ax );
    checkCudaErr( "cuFreq failed for Dkx");
    blocksPerGrid = (ny + threadsPerBlock -1 )/threadsPerBlock;  //  round up
    cuFreq<<<blocksPerGrid, threadsPerBlock>>>( Dky, Dky2, ny, by );
    checkCudaErr( "cuFreq failed for Dky");
 
   //   allocate wave memory on the device 
    cudaMalloc( (void**)&Dwave, sizeof(cufftComplex)*nx*ny);
    checkCudaErr( "failed to allocate cuda array Dwave" );

    //  allocate host memory for wave function
    Hwave = new cufftComplex[ nx*ny ];
    if( NULL == Hwave ) {
        sbuffer = "cannot allocate cuda host array Hwave\n";
        messageAS( sbuffer, 0 );
        return( -1 );
    }

    for( ix=0; ix<nx; ix++) for( iy=0; iy<ny; iy++) {
        Hwave[ iy + ix*ny ].x = wave.re(ix,iy);
        Hwave[ iy + ix*ny ].y = wave.im(ix,iy);
    }
    //  copy this complex wave to device
    cudaMemcpy( Dwave, Hwave, nx*ny*sizeof(cufftComplex),
                        cudaMemcpyHostToDevice );
    checkCudaErr( "cannot copy host to device for Dwave");


    if( lcross == 1 ) {
        /* nz may be too small with thermal vibrations so add a few extra */
        nz = (int) ( (zmax-zmin)/ deltaz + 3.5);
        depthpix0.resize( nx, nz );
        for( ix=0; ix<nx; ix++)
        for( iz=0; iz<nz; iz++)  depthpix0.re(ix,iz) = depthpix0.im(ix,iz) = 0.0F;
        iycross = (int) ( 0.5 + (ny * ycross / by));
        while( iycross < 0 ) iycross += ny;
        iycross = iycross%ny;  /* make periodic in ny */
        sbuf = "save xz cross section at iy= "+toString(iycross)+" pixels";
        messageAS( sbuf );
    }

    /*  ------  */
 
    k2max = nx/(2.0F*ax);
    tctx = ny/(2.0F*by);
    if( tctx < k2max ) k2max = tctx;
    k2maxt = k2max;
    k2max = BW * k2max;
    if( verbose > 0 ) {
        sbuf= "Bandwidth limited to a real space resolution of "+toString(1.0F/k2max)
                +" Angstroms";
        messageAS( sbuf );
        sbuf = "   (= " + toString(wavlen*k2max*1000.0F)
                + " mrad)  for symmetrical anti-aliasing.";
        messageAS( sbuf );
    }
    k2max = k2max*k2max;
    k2maxt = k2maxt * k2maxt;  // for trans function potential

#ifdef AST_USE_CUATOMPOT 
    if( cufftPlan2d( &cuplanTc2r, nx, ny, CUFFT_C2R) != CUFFT_SUCCESS) {
    sbuffer = "cuda error; unable to create plan Tc2r\n" ; 
        messageAS( sbuffer, 0 );
        return( -1 );
    }

    //  allocate fparams on on host and device
    float fex = featom( 12, 0.1 );  //  force init of fparams[][]
    const int NPMAX=   12;  // number of parameters for each Z
    int natomin = natom;  //  copt from autostem_cuda.cpp
    Hfparams = new double [(NZMAX+1)*NPMAX];
    cudaMalloc( (void**)&Dfparams, sizeof(double)*(NZMAX+1)*NPMAX );
    checkCudaErr( "failed to allocate device memory for fparams");
    for( ix=0; ix< NPMAX; ix++) for(iy=1; iy<=NZMAX; iy++)
        Hfparams [ix + iy*NPMAX] = fparams[iy][ix];

    //--- copy to device ---------
    cudaMemcpy( Dfparams, Hfparams, (NZMAX+1)*NPMAX * sizeof(double),
             cudaMemcpyHostToDevice );
    checkCudaErr( "cannot copy Hfparams from host to device");
    delete [] Hfparams;

    //  allocate host memory for atom positions (plain array not vector to copy)
    Hspec = new float[ 4*natomin ];
    if( NULL == Hspec ) {  //  should throw an exception but just in case
        sbuffer = "cannot allocate cuda host array for Hspec\n";
        messageAS( sbuffer, 0 );
        return( -1 );
    }
    //  make associated device memory
    cudaMalloc( (void**)&Dspec, 4*natomin*sizeof(float) );
    checkCudaErr( "failed to allocate device memory for Dspec");

    for( ix=0; ix< natom; ix++) {  // make simple arrays to copy
        Hspec[0 + 4*ix] = x[ix]; // pack all 4 into one to min
        Hspec[1 + 4*ix] = y[ix]; // number of transfers
        Hspec[2 + 4*ix] = occ[ix];
        Hspec[3 + 4*ix] = Znum[ix];
    }
    cudaMemcpy(Dspec, Hspec, 4*natom*sizeof(float), cudaMemcpyHostToDevice );
    checkCudaErr( "cannot copy Dspec from host to device");

    cudaMalloc( (void**)&DpotnC, sizeof(cufftComplex)*nx*(ny/2 + 1) );
    checkCudaErr( "failed to allocate cuda array DpotnC" );

#endif
 
/*  iterate the multislice algorithm proper

   NOTE: zero freg is in the bottom left corner and
     expands into all other corners - not in the center
     this is required for the FFT - don't waste time rearranging

*/

/* ---- start coherent method below ----------------
        (remember that wave was initialize above) */

    if( lbeams ==1 ) {
        nzbeams = (int) ( (zmax-zmin)/ deltaz + 3.5);
        beams.resize(nbout+1, nzbeams );   // to save values
        for(ib=0; ib<nbout; ib++) {
                hbeam[ib] = hb[ib];
                kbeam[ib] = kb[ib];
        }

        //  make them all positive just in case
        for( ib=0; ib<nbout; ib++) {
            if( hbeam[ib] < 0 ) hbeam[ib] = nx + hbeam[ib];
            if( kbeam[ib] < 0 ) kbeam[ib] = ny + kbeam[ib];
            if( hbeam[ib] < 0 ) hbeam[ib] = 0;
            if( kbeam[ib] < 0 ) kbeam[ib] = 0;
            if( hbeam[ib] > nx-1 ) hbeam[ib] = nx-1;
            if( kbeam[ib] > ny-1 ) kbeam[ib] = ny-1;
        }
    }

    if( verbose > 0 ) {
        sbuf = "Sorting atoms by depth...";
        messageAS( sbuf );
    }
    sortByZ( x, y, z, occ, Znum, natom );

    if( lwobble == 1 ){
        zmin = z[0];        /* reset zmin/max after wobble */
        zmax = z[natom-1];
        if( verbose > 0 ) {
            sbuf ="Thickness range "
                " is "+toString(zmin)+" to "+toString(zmax)+" (in z)";
            messageAS( sbuf );
        }
    }

    scale = 1.0F / ( ((float)nx) * ((float)ny) );

    zslice = 0.75*deltaz;  /*  start a little before top of unit cell */
    istart = 0;
    islice = 1;

    // ------ CUDA thread index parameters - may need to change for differenet GPUs (?)
    //  for probe as 1D
    int N = nx * ny;
    threadsPerBlock = 512;  // max thread/block usually 1024 =dev specific
    blocksPerGrid = (N + threadsPerBlock -1 )/threadsPerBlock;  //  round up

    //  for transmission function as 2D
    int thr3 = 16;              
    dim3 threads3(thr3,thr3);
    int blk3xt =  (nx + thr3 -1 )/thr3;  //  round up
    int blk3yt =  (ny + thr3 -1 )/thr3;  //  round up
    dim3 blocks3t(blk3xt, blk3yt);

    while( (istart < natom) && ( zslice < (zmax+deltaz) ) ) {

        /* find range of atoms for current slice */
        na = 0;
        for(i=istart; i<natom; i++) 
        if( z[i] < zslice ) na++; else break;

        /* calculate transmission function, skip if layer empty */
        if( na > 0 ) {

#ifdef AST_USE_CUATOMPOT
            //  -- its about the same speed doing the reciprocal space
            //     potential on GPU and real space potential on host
            //  -- so leave both here for future reference

            //  do atomic potnetial on GPU in reciprocal space not real space
            //   to work with fine grain parallel cores

            float mm0 = 1.0F + v0/511.0F;
            float scale2 = wavlen * mm0/(ax*by);

            cuAtompot<<<blocks3t, threads3>>>(  DpotnC, Dspec, na, istart, 
                (float)ax, (float)by, (float)v0, 
                nx, ny, Dkx, Dky, Dkx2, Dky2, (float) k2maxt, Dfparams, scale2 );
            phirms = 0;
            nbeams = 0;  // ???? need this from somewhere else
            checkCudaErr( "cuAtompot() failed");

            cufftExecC2R( cuplanTc2r, DpotnC, DpotnR );
            checkCudaErr( "cuFFTc2r forward failed");
#else
            //  calculate just real valued sigma*Vz on host - usually faster
            //???? if( istart == 0 ) trlayer( x, y, occ,
            trlayer( x, y, occ,
                Znum, na, istart, ax, by, v0, trans,
                nx, ny, kx2, ky2, &phirms, &nbeams, -1.0F );

            nbeams = 0;  // ???? need this from somewhere else
            for( ix=0; ix<nx; ix++) for( iy=0; iy<ny; iy++)
                HpotnR[ iy + ix*ny ] = trans.re(ix,iy);

            //  copy this real potential to device
            cudaMemcpy( DpotnR, HpotnR, nx*ny*sizeof(cufftReal),
                        cudaMemcpyHostToDevice );
            checkCudaErr( "cannot copy host to device for DpotnR");
#endif

            cuPhasegrating<<<blocks3t, threads3>>>( DpotnR, Dtrans, nx, ny );
            checkCudaErr( "cuPhasegrating() failed");

            cufftExecC2C( cuplan, Dtrans, Dtrans, CUFFT_FORWARD);
            checkCudaErr( "cuFFT forward failed");

            cuBWlimit<<<blocks3t, threads3>>>( Dtrans,
                Dkx2, Dky2, (float) k2max, nx, ny );
            checkCudaErr( "cuBWlimit() failed");

            cufftExecC2C( cuplan, Dtrans, Dtrans, CUFFT_INVERSE);
            checkCudaErr( "cuFFT inverse failed");

            //   wave *= trans;    //  transmit
            cmplVecMul<<<blocksPerGrid, threadsPerBlock>>>( Dwave,
                    Dtrans, Dwave, N );
            checkCudaErr( "cmplVecMul() failed");
      

        }

        /*  bandwidth limit */
        //wave.fft();
        cufftExecC2C( cuplan, Dwave, Dwave, CUFFT_FORWARD);
        checkCudaErr( "cuFFT forward failed");

        if( (lbeams== 1) && (islice<nzbeams) && (islice>0) )  {
            //  copy back to host - will be very slow - should do this better
            cudaMemcpy( Hwave, Dwave, nx*ny*sizeof(cufftComplex),
                        cudaMemcpyDeviceToHost );
            checkCudaErr( "cannot copy host to device for Dwave"); 
            for( ib=0; ib<nbout; ib++) {
                //beams.re(ib,islice-1) = scale*wave.re(hbeam[ib],kbeam[ib] );   // real
                //beams.im(ib,islice-1) = scale*wave.im(hbeam[ib],kbeam[ib] );   // imag
                beams.re(ib,islice-1) = scale*Hwave[ kbeam[ib] + hbeam[ib]*ny].x;   // real
                beams.im(ib,islice-1) = scale*Hwave[ kbeam[ib] + hbeam[ib]*ny].y;   // imag
            }
            beams.re(nbout, islice-1) = (float) zslice;  //  save actual z coord.
        }

        /* remember: prop needed here to get anti-aliasing right */
        //wave *= cprop;
        //wave.ifft();

        //  multiplied by the propagator function 
        cmplVecMul<<<blocksPerGrid, threadsPerBlock>>>( Dwave,
                Dprop, Dwave, N );
        checkCudaErr( "cmplVecMul() failed");

        cufftExecC2C( cuplan, Dwave, Dwave, CUFFT_INVERSE);
        checkCudaErr( "cuFFT inverse failed");

        /* save depth cross section if requested */
        if( (lcross == 1) && (islice<=nz) ) {
            //  copy back to host - will be very slow - should do this better
            cudaMemcpy( Hwave, Dwave, nx*ny*sizeof(cufftComplex),
                        cudaMemcpyDeviceToHost );
            checkCudaErr( "cannot copy host to device for Dwave"); 
            for( ix=0; ix<nx; ix++) {
                depthpix0.re(ix, islice-1) = 
                    //wave.re(ix,iycross)*wave.re(ix,iycross)
                    //    + wave.im(ix,iycross)*wave.im(ix,iycross);
                    Hwave[iycross+ix*ny].x*Hwave[iycross+ix*ny].x
                        + Hwave[iycross+ix*ny].y*Hwave[iycross+ix*ny].y;
            }
            nzout = islice;
        }

        if( verbose > 0 ) {

            //  do mag. sq on GPU device (use STEM detector functions)
            magSqPix<<<blocks3t, threads3>>>( Dcbed, Dwave, nx, ny );
            checkCudaErr( "magSqPix failed");

            zeroDbleArray<<<ny, 1>>>( Dsums, ny );
            checkCudaErr( "zeroDbleArray failed");

            int TOTAL = 3;  //  from autostem detector modes
            //  integrate in real space NOT Fourier space - but still works
            //  Dkx,Dky,Dkx2,Dky2 not used for total so not initialized
            //  k2min,max and phiMin,Max = not used either but need variable
            integCBED<<<blocks3t, threads3>>>( Dsums, Dcbed, nx, ny,
                        TOTAL, Dkx, Dky, Dkx2, Dky2, 0.0, 0.0, 0.0, 0.0 );
            checkCudaErr( "integCBED failed");

            cudaMemcpy( Hsums, Dsums, ny*sizeof(double),
                        cudaMemcpyDeviceToHost );
            checkCudaErr( "cannot copy sums device to host");
            sum = 0.0;
            for( iy=0; iy<ny; iy++) sum += Hsums[iy];
            sum = sum * scale;

            sbuf = "z= " + toString(zslice)+" A, " + toString(nbeams) + " beams, "
                    + toString(na)+" coord., \n"
                    + "     aver. phase= "+toString(phirms)
                    +", total intensity = "+toString(sum) ;
            messageAS( sbuf );
        }

        zslice += deltaz;
        istart += na;
        islice++;

    } /* end while(istart<natom..) */

    pix.resize(nx,ny);
    //pix = wave;
    cudaMemcpy( Hwave, Dwave, nx*ny*sizeof(cufftComplex), cudaMemcpyDeviceToHost );
    checkCudaErr( "cannot copy probe device to host");
    for( ix=0; ix<nx; ix++) 
        for( iy=0; iy<ny; iy++) {
            pix.re(ix,iy) = Hwave[ iy + ix*ny].x;
            pix.im(ix,iy) = Hwave[ iy + ix*ny].y;
    }


    //-------------  free GPU memory etc.

    cufftDestroy( cuplan );
    cufftDestroy( cuplanTc2r );
    cudaFree( Dkx );
    cudaFree( Dky );
    cudaFree( Dkx2 );
    cudaFree( Dky2 );
    cudaFree( Dprop );
    cudaFree( Dtrans );
    cudaFree( DpotnR );
    cudaFree( Dcbed );
    cudaFree( Dwave );
    cudaFree( Dsums );
    delete [] Hcbed;
    delete [] Htrans;
    delete [] Hwave;
    delete [] Hsums;
#ifdef AST_USE_CUATOMPOT 
    cudaFree( Dfparams );
    cudaFree( Dspec );
    cudaFree( DpotnC );
    delete [] Hspec;
#endif
 
    //----------- end and exit --------------------
    // save depth cross section if requested and get rid of possible extra lines
    if( 1 == lcross ) {
        depthpix.resize( nx, nzout );  //  return to calling program
        for( ix=0; ix<nx; ix++)        //  depthpix0 will get deallocated
        for( iz=0; iz<nzout; iz++)  {
            depthpix.re(ix,iz) = depthpix0.re(ix,iz);
            depthpix.im(ix,iz) = depthpix0.im(ix,iz);
        }
    }

    return( +1 );

};   //  end autoslic::calculate() - cuda version
#endif

//=============================================================
/*  calculatePartial()

loop over calculate() to sum a partially coherent image
multithread over phonon configurations (nwobble)

  input:
        param[] = image parameters (most will not be changed here)
        multimode = flag controlling multipole aberrations
        natom = number of atoms
        x[],y[],z[] = atomic coord
        Znum[] atomic number of each atom
        occ[] = occupancy of each atomic site
        wobble[] = thermal oscillation amplitude 
        hb[],kb[] = (h,k) indexes of beams to monitor during
                          propagation
                (ignore if lbeams = 0)
        nbout  = number of beams to record
        ycross = y position to save xz cross section

  mode flags:  lbeams, lcross, lpartl, lstart, lwobble

  output:
        pix = complex image to get results, orig. data lost and may be resized
                (complex for coherent calc. and real for partial coherence)
        beams = track specified beams as wave propagates
                (ignored if lbeams = 0)
*       depthpix = will be resized and get xz cross section image
                save intensity in real part (imag part not used)
                (ignored if lcross = 0)

    return value 
        >= 0 for success
        < 0 for failure
*/
int autoslic::calculatePartial(cfpix &pix, vectorf &param, int multiMode, 
        int natom, unsigned long *iseed, vectori &Znum, 
        vectorf &x, vectorf &y, vectorf &z, vectorf &occ, vectorf &wobble,
        float dfdelt )
{
    int i, ix, iy, nx, ny, nwobble, np, iverbose,
        nacx,nacy, iqx, iqy, iwobble, ndf, idf, n1, n2, nbout;

    float wmin, wmax, xmin,xmax, ymin, ymax, zmin, zmax;
    float k2, k2max, scale, v0, wavlen, rx, ry,
        ax, by, ycross, tctx, acmin, acmax, df, df0, sigmaf,
        aobj, qx, qy, qy2, q2, q2min, q2max, sumdf, pdf, k2maxo,
        temperature;
    float tr, ti, wr, wi;

    double sum, xdf, chi0, t, alx, aly;

    vectori hbeam, kbeam;

    cfpix *wave;            // complex probe wave functions
    cfpix *temp, *pixw;     // complex scratch wave function

    cfpix depthpix, beams;  //  dummy argument needed for calculate()

    // ---- get setup parameters from param[]
    ax = param[ pAX ];
    by = param[ pBY ];
    nx = ToInt( param[ pNX ] );
    ny = ToInt( param[ pNY ] );
    v0 = param[pENERGY];                // electron beam energy in keV
    df0 = param[pDEFOCUS];              // defocus
    sigmaf = param[pDDF];
    acmax = param[pCAPERT];             // condencer angles
    acmin = param[pCAPERTMIN];
    aobj = param[ pOAPERT ];            // objective aperture
    temperature = param[ pTEMPER ];     // temperature
    nwobble = ToInt( param[ pNWOBBLE ] );       //  number config. to average

    lbeams = 0;  //  can't record beams in this mode
    nbout = 0;
    ycross = 10.0F;  //  dummy variable for calculate()

    if( nwobble < 1 ) nwobble = 1;      //  has to be at least one configuration
    if( 0 == lwobble) nwobble = 1;      //  don't use frozen phonon offsets

    if( (nx < 1) || (ny < 1) || (ax<0.0) || (by<0.0) ){
        sbuffer="bad size parameters in autoslic::calculate()";
        messageAS( sbuffer );
        return( -1 );
    }
            
    /*  calculate relativistic factor and electron wavelength */
    wavlen = (float) wavelength( v0 );

    /*  calculate the total specimen volume and echo - may not be needed ????? */
    xmin = xmax = x[0];
    ymin = ymax = y[0];
    zmin = zmax = z[0];
    wmin = wmax = wobble[0];

    for( i=0; i<natom; i++) {
        if( x[i] < xmin ) xmin = x[i];
        if( x[i] > xmax ) xmax = x[i];
        if( y[i] < ymin ) ymin = y[i];
        if( y[i] > ymax ) ymax = y[i];
        if( z[i] < zmin ) zmin = z[i];
        if( z[i] > zmax ) zmax = z[i];
        if( wobble[i] < wmin ) wmin = wobble[i];
        if( wobble[i] > wmax ) wmax = wobble[i];
    }
    // --- leave this in main calling program
    //sprintf(stemp, "Total specimen range is\n %g to %g in x\n"
    //       " %g to %g in y\n %g to %g in z",
    //         xmin, xmax, ymin, ymax, zmin, zmax );
    //messageAS( stemp );

    // --- leave this in main calling program
    //if( lwobble == 1 ) {
    //    sprintf(stemp, "Range of thermal rms displacements (300K) = %g to %g\n",
    //        wmin, wmax );
    //    messageAS( stemp );t
    //}

/*  calculate spatial frequencies and positions for future use */

    rx = 1.0F/ax;       //???  may not be needed anymore
    ry = 1.0F/by;

    initAS( param, Znum, natom );  //  init for calculate()

    /*---- allocate some more arrays and initialize wavefunction ----*/

    wave = new cfpix[ nwobble ];
    temp = new cfpix[ nwobble ];
    pixw = new cfpix[ nwobble ];
    if( (NULL == wave) || (NULL == temp) || (NULL == pixw) ) {
        sbuffer = "Cannot allocate wave,temp,pix array";
        messageAS( sbuffer, 2 );
        return( -1 );
    }
    wave[0].resize(nx, ny );
    wave[0].init();
    temp[0].resize( nx, ny );
    temp[0].copyInit( wave[0] );
    pixw[0].resize( nx, ny );
    pixw[0].copyInit( wave[0] );

    for( iwobble=1; iwobble<nwobble; iwobble++){
        wave[iwobble].resize( nx, ny );
        temp[iwobble].resize( nx, ny );
        pixw[iwobble].resize( nx, ny );

        wave[iwobble].copyInit( wave[0] );
        temp[iwobble].copyInit( wave[0] );
        pixw[iwobble].copyInit( wave[0] );
    }
 
    k2max = nx/(2.0F*ax);
    tctx = ny/(2.0F*by);
    if( tctx < k2max ) k2max = tctx;
    k2max = BW * k2max;
    if( echo > 0 ) {
        sbuffer= "Bandwidth limited to a real space resolution of "+toString(1.0F/k2max)
                +" Angstroms";
        messageAS( sbuffer );
        sbuffer= "   (= " + toString(wavlen*k2max*1000.0F)
                + " mrad)  for symmetrical anti-aliasing.";
        messageAS( sbuffer );
    }
    k2max = k2max*k2max;

/*  iterate the multislice algorithm proper

   NOTE: zero freg is in the bottom left corner and
     expands into all other corners - not in the center
     this is required for the FFT - don't waste time rearranging

  partial coherence method
   force the integrals to include the origin and to be symmetric
   about the origin and to have the same periodic boundary
   conditions as the sampling grid
*/
        
    //sprintf(stemp,"Illumination angle sampling (in mrad) = %f, %f\n",
    //    1000.*rx*wavlen, 1000.*ry*wavlen);
    //messageAS( stemp );

    pix.resize(nx,ny);
    pix = 0.0F;     // start with zero and sum into this pix

    if( fabs( (double) dfdelt ) < 1.0 ) ndf = 1; 
    else ndf = (int) ( ( 2.5F * sigmaf ) / dfdelt );

    nacx = (int) ( ( acmax / ( wavlen * rx ) ) + 1.5F );
    nacy = (int) ( ( acmax / ( wavlen * ry ) ) + 1.5F );

    q2max = acmax / wavlen;
    q2max = q2max*q2max;

    q2min = acmin / wavlen;
    q2min = q2min*q2min;

    k2maxo = aobj / wavlen;
    k2maxo = k2maxo*k2maxo;

    // for Monte Carlo stuff
    //  make nwobble set of coord serially so RNG works
    vector< vector<float> > x2( nwobble, x);
    vector< vector<float> > y2( nwobble, y);
    vector< vector<float> > z2( nwobble, z);
    vector< vector<float> > occ2( nwobble, occ);
    vector< vector<int> > Znum2( nwobble, Znum);
    vector< vector<float> > param2( nwobble, param );

    /*  add random thermal displacements scaled by temperature
            if requested 
        remember that initial wobble is at 300K for each direction */

    np = (int) param.size();
    scale = (float) sqrt(temperature/300.0) ;

    for( iwobble=0; iwobble<nwobble; iwobble++)
            for( i=0; i<np; i++) param2[iwobble][i] = param[i];

    if( (lwobble == 1) ) {          //  add frozen phonon random displacements
        for( iwobble=0; iwobble<nwobble; iwobble++) {
            for( i=0; i<natom; i++) {
                x2[iwobble][i] = x[i] + (float)(wobble[i]*rangauss(iseed)*scale);
                y2[iwobble][i] = y[i] + (float)(wobble[i]*rangauss(iseed)*scale);
                z2[iwobble][i] = z[i] + (float)(wobble[i]*rangauss(iseed)*scale);
                occ2[iwobble][i] = occ[i];
                Znum2[iwobble][i] = Znum[i];
            }
        }
    } else {                        //  just copy the original
        for( i=0; i<natom; i++) {
            x2[0][i] = x[i];
            y2[0][i] = y[i];
            z2[0][i] = z[i];
            occ2[0][i] = occ[i];
            Znum2[0][i] = Znum[i];
        }
    }

    vector< string > str( nwobble );    //  must have separate string for each thread

    //  same constant for all threads so take outside of loop
    if( fabs( (double) sigmaf ) < 1.0 ) n1 = n2 = 0;
    else {
            n1 = -ndf;
            n2 = ndf;
    }

    nillum = 0;

    //---  make separate thread for each TDS configuration
    //---  multithread-1
#pragma omp parallel for private(iqx,iqy,qx,qy,qy2,q2,t,ix,iy,tr,ti,wr,wi,alx,aly,idf,xdf,pdf,sum,k2)
    for( iwobble=0; iwobble<nwobble; iwobble++) {
        if( (lwobble == 1) && ( echo > 0 ) ) {
            str[iwobble] = "configuration # " + toString( iwobble+1 );
           messageAS( str[iwobble] );
        }
        pixw[iwobble] = 0.0F;
        //  integrate over the illumination angles
        for( iqy= -nacy; iqy<=nacy; iqy++) {
            qy = iqy * ry;
            qy2 = qy * qy;
        
            for( iqx= -nacx; iqx<=nacx; iqx++) {
                qx = iqx * rx;
                q2 = qx*qx + qy2;
        
                if( (q2 <= q2max) && (q2 >= q2min) ) {
                    if( 0 == iwobble) nillum += 1;
                    for( ix=0; ix<nx; ix++) {
                        for( iy=0; iy<ny; iy++) {
                            t = 2.0*pi*( qx*xpos[ix] + qy*ypos[iy] );
                            wave[iwobble].re(ix,iy) = (float) cos(t);  // real
                            wave[iwobble].im(ix,iy) = (float) sin(t);  // imag
                        }
                    }

                    //-----  transmit thru the specimen with this configuration
                    iverbose = 0;   //  turn off echo in calculate()
                    lstart = 1;     //  must start calculate() from this wave
                    calculate( temp[iwobble], wave[iwobble], depthpix, param2[iwobble], 
                        multiMode, natom, Znum2[iwobble],
                        x2[iwobble], y2[iwobble], z2[iwobble],
                        occ2[iwobble], beams, hbeam, kbeam, nbout,
                        ycross, iverbose );
          
                    wave[iwobble] = temp[iwobble];  //  copy back results (not efficient?)

                    sum = 0.0;
                    for( ix=0; ix<nx; ix++) {
                        for( iy=0; iy<ny; iy++)
                            sum += wave[iwobble].re(ix,iy)*wave[iwobble].re(ix,iy)
                                + wave[iwobble].im(ix,iy)*wave[iwobble].im(ix,iy);
                    }
                    sum = sum / ( ((float)nx) * ((float)ny) );

                    if( (0 == iwobble) && ( echo > 0 ) ) {
                        str[iwobble]=  "Illum. angle = " + toString(1000.*qx*wavlen) +
                            ", "+toString(1000.*qy*wavlen) +
                            " mrad, integ. intensity= "+toString(sum);
                        messageAS( str[iwobble] );
                    }
            
                    //-------- integrate over +/- 2.5 sigma of defocus ------------ 
                    //   should convert to Gauss-Hermite quadrature sometime
                    wave[iwobble].fft();
                    if( iwobble == 0 ) sumdf = 0.0F;

                    for( idf= n1; idf<=n2; idf++) {
                        param2[iwobble][pDEFOCUS] = df = df0 + idf*dfdelt;
            
                        for( ix=0; ix<nx; ix++) {
                            alx = wavlen * kx[ix];  // x component of angle alpha
                            for( iy=0; iy<ny; iy++) {
                                aly = wavlen * ky[iy];  // y component of angle alpha
                                k2 = kx2[ix] + ky2[iy];
                                if( k2 <= k2maxo ) {
                                    chi0 = (2.0*pi/wavlen) * chi( param2[iwobble], 
                                            alx, aly, multiMode );
                                    tr = (float)  cos(chi0);
                                    ti = (float) -sin(chi0);
                                    wr = wave[iwobble].re(ix,iy);
                                    wi = wave[iwobble].im(ix,iy);
                                    temp[iwobble].re(ix,iy) = wr*tr - wi*ti;
                                    temp[iwobble].im(ix,iy) = wr*ti + wi*tr;
                                } else {
                                    temp[iwobble].re(ix,iy) = 0.0F;  // real
                                    temp[iwobble].im(ix,iy) = 0.0F;  // imag
                                }
                            }  /*  end for( iy=0... ) */
                        }   /*  end for( ix=0... ) */

                        temp[iwobble].ifft();
            
                        if( (0==n1) && (0==n2) ) pdf = 1;
                        else {
                            xdf = (double) ( (df - df0) /sigmaf );
                            pdf = (float) exp( -0.5 * xdf*xdf );
                        }
                        if( iwobble == 0 ) sumdf += pdf;
            
                        for( ix=0; ix<nx; ix++) {
                            for( iy=0; iy<ny; iy++) {
                                wr = temp[iwobble].re(ix,iy);
                                wi = temp[iwobble].im(ix,iy);
                                pixw[iwobble].re(ix,iy) += pdf* ( wr*wr + wi*wi );
                            }
                        }
            
                    }/* end for(idf..) */

                    param2[iwobble][ pDEFOCUS ] = df0;  // return to original value

                }/* end if( q2...) */
        
            } /* end for( iqx..) */
        } /* end for( iqy..) */
    } /* end for( iwobble...) */

    //----  put these in main calling program if neede
    //sprintf(stemp, "Total number of illumination angle = %ld",
    //        nillum);
    //message ( stemp );
    //sprintf(stemp, "Total number of defocus values = %d", 2*ndf+1);
    //messageAS( stemp );

    //  sum results from each thread
    pix = pixw[0];
    if( nwobble > 1 ) for( iwobble=1; iwobble<nwobble; iwobble++) {
        pix += pixw[iwobble];
    }

    // scale the whole sum
    scale = 1.0F / (sumdf *(float)(nillum*nwobble)); 
    for( ix=0; ix<nx; ix++)
    for( iy=0; iy<ny; iy++) {
        pix.re(ix,iy) *= scale;
    }

    //----------- end --------------------

    return( +1 );

};   //  end autoslic::calculatePartial()


//=============================================================
/*  calculateCBED_TDS()

loop over calculate() to sum a CBED with frozen phonons
multithread over phonon configurations (nwobble)

  input:
        param[] = image parameters (most will not be changed here)
        multimode = flag controlling multipole aberrations
        natom = number of atoms
        x[],y[],z[] = atomic coord
        Znum[] atomic number of each atom
        occ[] = occupancy of each atomic site
        wobble[] = thermal oscillation amplitude 

  mode flags:  lbeams, lcross, lpartl, lstart, lwobble

    lcbed = 1 for CBED and 0 for normal e-diffraction

  output:
        pix = complex image to get results, orig. data lost and may be resized
                (complex for coherent calc. and real for partial coherence)

    return value 
        >= 0 for success
        < 0 for failure
*/
int autoslic::calculateCBED_TDS(cfpix &pix, vectorf &param, 
        int multiMode, int natom, unsigned long *iseed, vectori &Znum, 
        vectorf &x, vectorf &y, vectorf &z, vectorf &occ, vectorf &wobble )
{
    int i, ix, iy, nx, ny, nwobble, iverbose, ismoth,
         iwobble, npixels, nbout;

    float wmin, wmax, xmin,xmax, ymin, ymax, zmin, zmax;
    float  scale, v0, wavlen, rx, ry, rx2,ry2,
        ax, by, ycross, xp, yp,
        aobj, temperature, k2maxo;
    float tr, ti, ds;

    double pixel, sum;

    vectori hbeam, kbeam;

    cfpix wave0;        // complex probe wave functions
    cfpix *temp;        // complex scratch wave function

    cfpix depthpix, beams;  //  dummy argument needed for calculate()

    probe prb;

    // ---- get setup parameters from param[]
    ax = param[ pAX ];
    by = param[ pBY ];
    nx = ToInt( param[ pNX ] );
    ny = ToInt( param[ pNY ] );
    v0 = param[pENERGY];                // electron beam energy in keV
    aobj = param[ pOAPERT ];            // objective aperture
    temperature = param[ pTEMPER ];     // temperature
    nwobble = ToInt( param[ pNWOBBLE ] );       //  number config. to average

    xp = param[ pPPOSX ];   //   probe position in Angst.
    yp = param[ pPPOSY ];

    lbeams = 0;  //  can't record beams in this mode
    nbout = 0;
    ycross = 10.0F;  //  dummy variable for calculate()

    if( nwobble < 1 ) nwobble = 1;      //  has to be at least one configuration

    if( (nx < 1) || (ny < 1) || (ax<0.0) || (by<0.0) ){
        sbuffer="bad size parameters in autoslic.calculateCBED_TDS()";
        messageAS( sbuffer );
        return( -1 );
    }
            
    /*  calculate relativistic factor and electron wavelength */
    wavlen = (float) wavelength( v0 );
    //printf("electron wavelength = %g Angstroms\n", wavlen);

    /*  calculate the total specimen volume and echo */
    xmin = xmax = x[0];
    ymin = ymax = y[0];
    zmin = zmax = z[0];
    wmin = wmax = wobble[0];

    for( i=0; i<natom; i++) {
        if( x[i] < xmin ) xmin = x[i];
        if( x[i] > xmax ) xmax = x[i];
        if( y[i] < ymin ) ymin = y[i];
        if( y[i] > ymax ) ymax = y[i];
        if( z[i] < zmin ) zmin = z[i];
        if( z[i] > zmax ) zmax = z[i];
        if( wobble[i] < wmin ) wmin = wobble[i];
        if( wobble[i] > wmax ) wmax = wobble[i];
    }

    //  calculate spatial frequencies and positions for future use
    rx = 1.0F/ax;
    rx2= rx*rx;
    ry = 1.0F/by;
    ry2= ry*ry;
    pixel = ( rx2 + ry2 );
    ismoth = 0;

    initAS( param, Znum, natom );  //  init for calculate()

    //---- allocate some more arrays and initialize wavefunction ----

    temp = new cfpix[ nwobble ];
    if( (NULL == temp) ) {
        sbuffer = "Cannot allocate temp array in autoslic.calculateCBED_TDS()";
        messageAS( sbuffer, 2 );
        return( -1 );
    }
    wave0.resize( nx, ny );
    wave0.init();

    for( iwobble=0; iwobble<nwobble; iwobble++){
        temp[iwobble].resize( nx, ny );
        temp[iwobble].copyInit( wave0 );
    }
 
    k2maxo = aobj / wavlen;     //  max in obj. aperture
    k2maxo = k2maxo * k2maxo;

/*  iterate the multislice algorithm proper

   NOTE: zero freg is in the bottom left corner and
     expands into all other corners - not in the center
     this is required for the FFT - don't waste time rearranging

  partial coherence method
   force the integrals to include the origin and to be symmetric
   about the origin and to have the same periodic boundary
   conditions as the sampling grid
*/
        
    //sprintf(stemp,"Illumination angle sampling (in mrad) = %f, %f\n",
    //    1000.*rx*wavlen, 1000.*ry*wavlen);
    //messageAS( stemp );

    pix.resize( nx, ny );
    pix.copyInit( wave0 );

    // for Monte Carlo stuff
    //  make nwobble set of coord serially so RNG works
    vector< vector<float> > x2( nwobble, x);
    vector< vector<float> > y2( nwobble, y);
    vector< vector<float> > z2( nwobble, z);
    vector< vector<float> > occ2( nwobble, occ);
    vector< vector<int> > Znum2( nwobble, Znum);

    /*  add random thermal displacements scaled by temperature
            if requested 
        remember that initial wobble is at 300K for each direction */

    scale = (float) sqrt(temperature/300.0) ;

    for( iwobble=0; iwobble<nwobble; iwobble++) {
        for( i=0; i<natom; i++) {
            x2[iwobble][i] = x[i] + (float)(wobble[i]*rangauss(iseed)*scale);
            y2[iwobble][i] = y[i] + (float)(wobble[i]*rangauss(iseed)*scale);
            z2[iwobble][i] = z[i] + (float)(wobble[i]*rangauss(iseed)*scale);
            occ2[iwobble][i] = occ[i];
            Znum2[iwobble][i] = Znum[i];
        }
    }

    if( 0 == lcbed ) {      // normal elect. diffraction
        wave0 = (float) ( 1.0/sqrt( ( (double)(nx)*((double)ny) ) ) );
    }  else {               //  CBED
        //---- make the incident probe wave function
        //   remember that probe is normalize in real space
        npixels = prb.makeProbe( wave0, nx, ny, xp, yp, 
                param, wavlen, k2maxo, pixel, multiMode, ismoth, kx, kx2, ky, ky2);
        sbuffer = toString(npixels)+" pixels in obj. apert.";
        messageAS( sbuffer );
    }

    vector< string > str( nwobble );    //  must have separate string for each thread

    //---  make separate thread for each TDS configuration
    //---  multithread-1
    iverbose = 0;       //  turn off echo in calculate()
    lstart = 1;         //  must start calculate() from this wave
#pragma omp parallel for private(ix,iy,tr,ti,sum,ds)
    for( iwobble=0; iwobble<nwobble; iwobble++) {
        if( (lwobble == 1) ) {
            str[iwobble] = "configuration # " + toString( iwobble+1 );
            messageAS( str[iwobble] );
        }
        
        //-----  transmit thru the specimen with this configuration
        calculate( temp[iwobble], wave0, depthpix, param, 
            multiMode, natom, Znum2[iwobble],
            x2[iwobble], y2[iwobble], z2[iwobble]
            ,occ2[iwobble], beams, hbeam, kbeam, nbout,
            ycross, iverbose );

            // convert to diffraction pattern intensity
            temp[iwobble].fft(); 
            sum = 0.0;
            for( ix=0; ix<nx; ix++) for(iy=0; iy<ny; iy++) {
                tr = temp[iwobble].re(ix,iy);
                ti = temp[iwobble].im(ix,iy);
                temp[iwobble].re(ix,iy) = ds = tr*tr + ti*ti;  // intensity in CBED
                //temp[iwobble].im(ix,iy) = 0.0F;  // should not be needed
                sum += ds;
            }
            str[iwobble] = "# " + toString( iwobble+1 )+ " total intensity = " + toString( sum/(nx*ny) );
            messageAS( str[iwobble] );
                   
     } /* end for( iwobble...) */

    //  add intensity of each phonon config to the total in pix
    pix = 0.0F;
    for( ix=0; ix<nx; ix++) for(iy=0; iy<ny; iy++) {
        for( iwobble=0; iwobble<nwobble; iwobble++) {
            pix.re(ix,iy) += temp[iwobble].re(ix,iy);
        }
        pix.re(ix,iy) = pix.re(ix,iy) / ((float)nwobble);
    }

    pix.invert2D();  // put zero in the center
    nillum = 1;

    //----------- end --------------------

    return( +1 );

};   //  end autoslic::calculateCBED_TDS()

//=============================================================
/* -------------------  initAS() -------------------

   separate variable initializaiton to avoid repeated init
   in multithreaded mode

   call once before calculate()
*/
void autoslic::initAS( vectorf &param, vectori Znum, int natom )
{
    int nx, ny, ix, iy;

    float ax, by, ctiltx, ctilty, tctx, tcty, scale, v0, wavlen, tt;
    double deltaz, tx, t;

    // ---- get setup parameters from param[]
    ax = param[ pAX ];
    by = param[ pBY ];
    nx = ToInt( param[ pNX ] );
    ny = ToInt( param[ pNY ] );
    ctiltx = param[ pXCTILT ];          // crystal tilt
    ctilty = param[ pYCTILT ];
    deltaz = param[ pDELTAZ ];          // slice thickness
    v0 = param[pENERGY];                // electron beam energy in keV
    
    // -------- spatial frequency
    kx.resize( nx );
    kx2.resize( nx );
    xpos.resize( nx );
    freqn( kx, kx2, xpos, nx, ax );

    ky.resize( ny );
    ky2.resize( ny );
    ypos.resize( ny );
    freqn( ky, ky2, ypos, ny, by );
  
    // ------- calculate propagator ----------
    tctx = (float) (2.0 * tan(ctiltx));
    tcty = (float) (2.0 * tan(ctilty));

    scale = pi * ((float)deltaz);
    wavlen = (float) wavelength( v0 );

    k2max = nx/(2.0F*ax);  //  must be same as calculate() - should make global
    tt = ny/(2.0F*by);
    if( tt < k2max ) k2max = tt;
    k2max = BW * k2max;
    k2max = k2max*k2max;

    cprop.resize( nx, ny );
    for( ix=0; ix<nx; ix++) {
        tx = ( kx2[ix]*wavlen - kx[ix]*tctx );
        for( iy=0; iy<ny; iy++) {
            if( (kx2[ix] + ky2[iy]) < k2max ) {
                t = scale * ( tx + ky2[iy]*wavlen - ky[iy]*tcty );
                cprop.re(ix,iy) = (float)  cos(t);
                cprop.im(ix,iy) = (float) -sin(t);
            } else {
                cprop.re(ix,iy) = 0.0F;
                cprop.im(ix,iy) = 0.0F;
            }  //  end if( kx2[ix]... 
        } // end for(iy..) 
    } // end for(ix..) 

#ifdef USE_OPENMP
    /*  force LUT init. to avoid redundant init in parallel form */ 
    int i;
    double rsq, vz;
    rsq = 0.5;  /* arbitrary position */   
    for( i=0; i<natom; i++) vz =  vzatomLUT( Znum[i], rsq );
#endif

}  // end autoslic::initAS()


//=============================================================
/* -------------------  messageAS() -------------------
   message output
   direct all output message here to redirect to the command line
   or a GUI status line or message box when appropriate

   msg[] = character string with message to disply
   level = level of seriousness
        0 = simple status message
        1 = significant warning
        2 = possibly fatal error
*/
void autoslic::messageAS( std::string &smsg,  int level )
{
        messageSL( smsg.c_str(), level );  //  just call slicelib version for now

}  // end autoslic::messageAS()

//=============================================================
/*--------------------- trlayer() -----------------------*/
/*   same subroutine in autoslic.cpp and autostem.cpp

  Calculate complex specimen transmission function
  for one layer using real space projected atomic potentials

  x[],y[] = real array of atomic coordinates
  occ[]   = real array of occupancies
  Znum[]  = array of atomic numbers
  istart  = starting index of atom coord.
  natom   = number of atoms
  ax, by  = size of transmission function in Angstroms
  kev     = beam energy in keV
  trans   = 2D array to get complex specimen
        transmission function
  nx, ny  = dimensions of transmission functions
  *phirms = average phase shift of projected atomic potential
  *nbeams = will get number of Fourier coefficients
  k2max   = square of max k = bandwidth limit

  convert to cfpix class for trans 10-nov-2012 ejk
  convert arrays to vector<> 25-dec-2017 ejk
  add option to return just potential if k2max<0 28-sep-2018 ejk

*/
void autoslic::trlayer(  const vectorf &x, const vectorf &y, const vectorf &occ,
            const vectori &Znum, const int natom, const int istart,
            const float ax, const float by,
            const float kev, cfpix &trans, const int nx, const int ny,
            const vectorf &kx2, const vectorf &ky2,
            double *phirms, int *nbeams, const float k2max  )
{
    int idx, idy, i, ixo, iyo, ix, iy, ixw, iyw, nx1, nx2, ny1, ny2;
    float k2;
    double r, rx2, rsq, vz, rmin, rminsq, sum, scale, scalex, scaley;

    const double rmax=3.0, rmax2=rmax*rmax; /* max atomic radius in Angstroms */

    scale = sigma( kev ) / 1000.0;  /* in 1/(volt-Angstroms) */

    scalex = ax/nx;
    scaley = by/ny;

    /* min radius to avoid  singularity */
    rmin = ax/((double)nx);
    r = by/((double)ny);
    rmin =  0.25 * sqrt( 0.5*(rmin*rmin + r*r) );
    rminsq = rmin*rmin;

    idx = (int) ( nx*rmax/ax ) + 1;
    idy = (int) ( ny*rmax/by ) + 1;

    for( ix=0; ix<nx; ix++) {
        for( iy=0; iy<ny; iy++)
           trans.re(ix,iy) = 0.0F;    // real part
    }

/*  run this in parallel  */
//  there is only a small improvement (maybe 20%) here but prevents multithreading 
//  at a higher level so don't
/*#pragma omp parallel for private(ix,iy,ixo,iyo,nx1,nx2,ny1,ny2,rx2,ixw,iyw,vz,rsq) */
//----- multithread-2
//#pragma omp parallel for private(ix,iy,ixo,iyo,nx1,nx2,ny1,ny2,rx2,ixw,iyw,vz,rsq)
    for( i=istart; i<(istart+natom); i++) {
       ixo = (int) ( x[i]/scalex );
       iyo = (int) ( y[i]/scaley );
       nx1 = ixo - idx;
       nx2 = ixo + idx;
       ny1 = iyo - idy;
       ny2 = iyo + idy;

    /* add proj. atomic potential at a local region near its center
       taking advantage of small range of atomic potential */

       for( ix=nx1; ix<=nx2; ix++) {
          rx2 = x[i] - ((double)ix)*scalex;
          rx2 = rx2 * rx2;
          ixw = ix;
          while( ixw < 0 ) ixw = ixw + nx;
          ixw = ixw % nx;
          for( iy=ny1; iy<=ny2; iy++) {
             rsq = y[i] - ((double)iy)*scaley;
             rsq = rx2 + rsq*rsq;
             if( rsq <= rmax2 ) {
                iyw = iy;
                while( iyw < 0 ) iyw = iyw + ny;
                iyw = iyw % ny;
                if( rsq < rminsq ) rsq = rminsq;
                /* r = sqrt( r );
                vz = occ[i] * scale * vzatom( Znum[i], r ); slow */
                vz = occ[i] * vzatomLUT( Znum[i], rsq );
                trans.re(ixw,iyw) += (float) vz;
             }
          } /* end for(iy... */
       }  /* end for(ix... */

    } /* end for(i=0... */

    if( k2max > 0.0 ){
        /* convert phase to a complex transmission function */
        sum = 0;
        for( ix=0; ix<nx; ix++) {
            for( iy=0; iy<ny; iy++) {
                vz = scale * trans.re(ix,iy);
                sum += vz;
                trans.re(ix,iy) = (float) cos( vz );
                trans.im(ix,iy) = (float) sin( vz );
            }
        }

        /* bandwidth limit the transmission function */
        *nbeams = 0;
        trans.fft();
        for( ix=0; ix<nx; ix++) {
            for( iy=0; iy<ny; iy++) {
                k2 = ky2[iy] + kx2[ix];
                if (k2 < k2max) *nbeams += 1;
                else trans.re(ix,iy) = trans.im(ix,iy) = 0.0F;
            }
        }
        trans.ifft();

    } else {

        //  just scale to sigma * Vz and return real valued potential
        sum = 0;
        for( ix=0; ix<nx; ix++) 
            for( iy=0; iy<ny; iy++) {
                trans.re(ix,iy) *= (float) scale;
                sum += trans.re(ix,iy);
            }
    }

    *phirms = sum / ( ((double)nx)*((double)ny) );
    
    return;

 }  /* end autoslic::trlayer() */
 
