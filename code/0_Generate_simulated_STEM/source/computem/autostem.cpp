/*      
    *** autostem.cpp ***       (normal C++)
    *** autostem_cuda.cu ***   (CUDA code version)

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

-----------------------------------------------------------------------------

  This code has both the openMP code and CUDA GPU code which is selected
  with symbol AST_USE_CUDA (can use only one at a time!)

  Remember that file name must end in .cu for CUDA version

  The CUDA code has only been tested under Ubuntu 16.04 with a GTX 1070 (8GB memory0
  and my not work in other configurations

  C++ and TIFF version
  this version uses FFTW 3 (net about a factor of 2X faster)
  see:   www.fftw.org

  to output the position averaged CBED pattern include the filename
  after the program name on the command line (in 2D image mode only)

  on Windows file libfftw3f-3.dll must be in the PATH

  on Linux build as:
  g++ -O -fopenmp -o autostem autostem.cpp slicelob.o
                      floatTIFF.o  cfpix.o -lfftw3f

  Calculate STEM images and line scans from a non-periodic
  distribution of (x,y,z) atomic coordinates using repetitive multislice
  
  The transmission function for each slice can take a lot of
  computer time.  This program propagates the STEM probe for many
  beam position through the specimen at the same time to avoid recalculating
  the specimen transmission function for each probe position. This
  requires a huge amount of memory.  In 1D line scan mode a 2D probe
  wave function is stored for all positions.  In 2D image mode the 2D
  probe wave functions for a whole line are stored simulataneously.

  this file is formatted for a tab size of 4 characters

  multithreaded code using openMP may be turned on/off with
  symbol USE_OPENMP (ignore undefined pragma warning if off)

  Ref:
  [1] Zhiheng Yu, Philip Batson, John Silcox, "Artifacts in Aberration
      Corrected ADF-STEM Imaging", Ultramicroscopy 96 (2003) p.275-284
  [2] J. M. LeBeau et al, "Position averaged convergent beam electron
      diff.: Theory and Applications", Ultramic. 110 (2010) p.118-125

  started from stemslic.c and autoslic.c 19-jul-1998 E. Kirkland
  convert to simultaneous transmission of many probes at the same
    time to reuse transmission functions  28-jul-1998 ejk
  finished 2D image mode (1D mode works OK) 29-jul-1998 ejk
  last modified 29-jul-1998 ejk
  add Mote-Carlo integration of source size 24-aug-1998 ejk
  fixed typo in random dither in y direction of probe position
    1-feb-1999 ejk
  fixed error in nprobes in image mode (should be nyout
      but it was nxout- at top question) 3-july-2001 ejk
  update memory allocation routines,
     change void main() to int main() for better portability,
     and add 5th order spherical aberration  3-mar-2004 ejk
  start modification to keep multiple thickness's 21-jul-2005 ejk
      - in working form 26-jul-2005 ejk
  put in faster sorting routine sortByZ() 6-feb-2006 ejk
  add some openMP multithreading 23-may-2006 ejk
  move openMP stuff into conditional compilation
     so I only need one version of the code,
     and move sortbyz() to sliclib 23-aug-2006 ejk
  add periodic() to put probe position back into the supercell
     and recal. x,y position without source size wobble
     for output in 1D mode 31-aug-2006 ejk
  change propagation range to be whole unit cell not just
     range of atoms to treat sparsely populated spec.
     better 14-sep-2006 ejk
  change range to start at -0.25*deltaz 26-sep-2006 ejk
  minor cosmetic changes 11-jul-2007 ejk
  minor fixes for openMP stuff 27-feb-2008 ejk
  move vzatomLUT() to slicelib.c and reformat for TAB size of 4 char
               11-jun-2008 ejk
  convert to GPL 3-jul-2008 ejk
  fix bug in multithreading (add more private var.) 
      5-nov-2008 ejk
  add confocal mode 31-jul-2009 E. Kirkland
  offset confocal collector by zslice (same ref as obj) 4-aug-2009 ejk
  fix param listing in confocal mode 3,6-dec-2009 ejk
  get return value of scanf() to remove warnings from gcc 4.4 
         10-feb-2010 ejk
  start conversion to FFTW  10-feb-2010 ejk
     in working form 16-feb-2010 ejk
  move vz LUT init (for openMP) to top 22-feb-2010 ejk
  fix sign convention in FFTW 21-mar-2010 ejk
  update comments 4-apr-2010 ejk
  add multipole aberrations to probe (but not collector yet)
       9-may-2011 ejk
  change detector max test to < from <= so adding many
     ADF detectors together will work 6-jul-2011 ejk
  add trap for undersampling the probe+aperture 3-mar-2012 ejk
  start conversion to floatTIFF and C++ 24-may-2012 ejk
  add option to save position averaged CBED 30-jun-2012 ejk
  convert to cfpix/fftw class from raw fftw 07-nov-2012 to 10-nov-2012 ejk
  fix problem with different probe size (in pixels) that was
     created when multipole aberr added 7-apr-2013 ejk
  start conversion to autostem class with separate cmd line front end
       27-jul-2013 ejk
  add CountBeams() 17-aug-2013 ejk
  change obj. apert. to radians in param[] to be consistent with other programs
      -was in mrad  22-aug-2013 ejk
  add string message 23-aug-2013 ejk
  add confocal parameters 24-aug-2013 ejk
  change probe to new and add delete and return status value in calculate
      14-sep-2013 ejk
  change RNG seed argument to referenece so it get updated for successive calls
       an lverbose 21-sep-2013 ejk
  fix confocal detection 5-oct-2013 ejk
  convert detector angles to radians 19-oct-2013 ejk
  remove tString() from here and use slicelib 5-dec-2013 ejk
  convert malloc1D() to vector<> 5-jul-2016 ejk
  convert malloc2D(), malloc3D() to new2D(), new3D() 24-sep-2017 ejk
  change arguments in trlayer() from float* to vector<>,
        and add istart 25-dec-2017 ejk
  add abbError() 28-feb-2018, 1-mar-2018 ejk
  fix typo in error message in STEMsignals() 7-mar-2018 ejk
  change calculate() arg to reference and add scale arg
     option to abbError() 29-mar-2018 ejk
  add error trap/message for delta z too small 31-mar-2018 ejk
  add lverbose test on more mesages for silent running
       22-apr-2018 ejk
  add segmented ADF detector 28-apr-2018 ejk
  fix one more calculate() arg reference 13-may-2018 ejk
  start GPU/cuda mode 3-jul-2018 ejk
  cuda code working 31-jul-2018, updated 3-16-aug-2018 ejk
  move cuda code to cudasubs.cu (can be used in other prog,)
      20-aug-2018  ejk
  add integCBED() GPU/cuda  30-aug-2018 ejk
  add fparams to GPU/cuda/cuAtompot() 10-sep-2018 ejk
  cuAtompot() working (but not faster) 21-sep-2018 ejk
  add option to trlayer() to return just real potential if k2max<0
      28-sep-2018 ejk
  convert propagator from 2 x 1D to a 2D array with less
    roundoff error 2-oct-2018
  fix small bug in phi=atan2() for segmented ADF detector
          3-oct-2018 ejk
  fix a few small typo's (comments and unused code) 11-dec-2018 ejk
  comments out sourcesize lines so nvcc doesn't complain 30-jan-2019 ejk

    this file is formatted for a TAB size of 4 characters 
*/


#include <ctime>  /* ANSI C libraries used */

#include <sstream>      // string streams

//
//  to select cuda/GPU code include autostem_cuda.hpp
//  and the filename extention for this file must end in .cu
//  (.cpp for normal C++) - its better put other selections in .hpp so
//  the calling program has more information - use only one .hpp
//
//#include "autostem_cuda.hpp"   // header for cuda nvcc version of this class
#include "autostem.hpp"        // header for C++ version of this class

//   set up cuda things
#ifdef AST_USE_CUDA
#include "cudaSlice.cu"
//#define AST_USE_CUATOMPOT  // to use cuAtompot() - usually do NOT do this
#ifdef AST_USE_CUATOMPOT
extern vector< vector<double> > fparams;  // scattering factor table
#endif
#else
#define USE_OPENMP      // define to use openMP 
#endif


#define MANY_ABERR      //  define to include many aberrations 


//=============================================================
//---------------  creator and destructor --------------

autostem::autostem()
{

        BW= (2.0F/3.0F);  // bandwidth limit
        ABERR= 1.0e-4;    // max error for a,b

        NZMAX= 103;    // max atomic number Z

        NRMAX=   100;    // number in look-up-table in vzatomLUT
        RMIN=    0.01;   // r (in Ang) range of LUT for vzatomLUT()
        RMAX=    5.0;

        pi = 4.0 * atan( 1.0 );

        // init control flags
        lwobble = l1d = lxzimage = lpacbed = 0;
        doConfocal = xFALSE;

        return;

}   //  end autostem::autostem()

autostem::~autostem()
{
}

//=============================================================
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
void autostem::abbError( vector<float> &p1, int np, int NPARAM, 
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
        messageAST( sbuffer, 0 );
    }
 
    for( i=0; i<np; i++) {
        x = (float) ( 2.0*ranflat(&iseed) - 1.0);  //  range = -1 to +1
        p1[paramToVary[i]] += (float) ( x  * dchiMax[orderToVary[i]] );
    }

    return;

}  //  end abbError()

//=============================================================
/*  calculate()

  input:
        param[] = image parameters (most will not be changed here)
        multimode = flag controlling multipole aberrations
        natomin = number of atoms
        x[],y[],z[] = atomic coord
        Znum[] atomic number of each atom
        occ[] = occupancy of each atomic site
        wobble[] = thermal oscillation amplitude
        iseed = random number seed - will be changed here

        xi,xf,yi,yf = range of output
        nxout, nyout = size of output in pixels 
                (nxout must be =1 for 1D mode)

        ThickSave[] = thickness levels to save
        nThick = number of thickness levels to save

        almin[], almax[] = detector angles in radians or Angst.
        collectorMode[] = type of detectors (ADF vs. confocal)
        ndetect = number of detector geometries

  mode flags:  lwobble, l1d, lxzimage, lpacbed

  output:
        pacbedPix[][] = (real) 2D image with position averaged CBED (nxprobe x nyprobe)
                                if lpacbed == 1

        pixr[][][] = (real) array of 2D output images
                        size; (ndetect*nThick) x nxout x nyout
                        indexed as pixr[id+it*ndetect][ix][iy]
                    for 1D line scan nxout=1 and nyout is length of line in pixels

        rmin[][], rmax[][] = range of output pixr indexed as [it][id]
                                size; nThick x ndetect

  return value is +1 for successs and negative for failure 
*/
int autostem::calculate( vectorf &param, int multiMode, int natomin, unsigned long *iseed,
        vectori &Znum, vectorf &xa, vectorf &ya, vectorf &za, vectorf &occ, vectorf &wobble,
        double xi, double xf, double yi, double yf, int nxout, int nyout,
        vectord &ThickSave, int nThick,
        vectord &almin, vectord &almax, vectori &collectorMode, int ndetect,
        vectord &phiMin, vectord &phiMax,
        float ***pixr, float  **rmin, float **rmax,
        float **pacbedPix )
{
    int ix, iy, i, idetect, iwobble, nwobble,
        nprobes, ip, it, nbeamp, nbeampo, ix2, iy2;

    float prr, pri, temp, temperature;

    double scale, sum, wx, w, ***detect,
       tctx, tcty, dx, dy, ctiltx, ctilty, k2maxa, k2maxb, k2;

    //double sourcesize, sourceFWHM;  //  MC source size is not practical

    vectord x, y, sums;

    // ---- get setup parameters from param[]
    ax = param[ pAX ];
    by = param[ pBY ];
    cz = param[ pCZ ];
    nx = ToInt( param[ pNX ] );         // size of transmission function (in pixels)
    ny = ToInt( param[ pNY ] );
    keV = param[ pENERGY ];             // electron beam energy in keV
    df = param[pDEFOCUS];               // defocus
    ctiltx = param[ pXCTILT ];          // crystal tilt
    ctilty = param[ pYCTILT ];
    apert1 = param[ pOAPMIN ];          // obj. apert size in radians (min, max for annular apert)
    apert2 = param[ pOAPERT ];
    temperature = param[ pTEMPER ];             // temperature
    nwobble = ToInt( param[ pNWOBBLE ] );       // number config. to average
    deltaz = param[ pDELTAZ ];                  // slice thickness
    //sourceFWHM = param[ pSOURCE ];              // source size  - doesn't work so not really used

    nxprobe = ToInt( param[ pNXPRB ] );  // probe size in pixels
    nyprobe = ToInt( param[ pNYPRB ] );
    
    //  confocal collector lens parameters
    dfC = param[ pCDF ];                // collector defocus
    dfa2C = param[ pCDFA2 ];            // collector astig 2nd order
    dfa2phiC = param[ pCDFA2PHI ];
    dfa3C = param[ pCDFA3 ];            // collector astig 2nd order
    dfa3phiC = param[ pCDFA3PHI ];
    Cs3C = param[ pCCS3 ];              // collector spherical aberr.
    Cs5C =  param[ pCCS5 ];
    apert1C = param[ pCCAPMIN ];        //  collector apert. in radians
    apert2C = param[ pCCAPMAX ];

    natom = natomin;
    wavlen = wavelength( keV );

    //  this code now works more consistently for both 1D and 2D
    nprobes = nyout;
    if( nxout < 1 ) {
        sbuffer = "nxout must be > 1 in autostem but it is "+toString(nxout);
        messageAST( sbuffer, 2 );
        return( -1 );
    }
    if( nyout < 1 ) {
        sbuffer = "nyout must be > 1 in autostem but it is "+toString(nyout);
        messageAST( sbuffer, 2 );
        return( -2 );
    }

    if( deltaz < 0.1F ) { 
        sbuffer = "delta z is too small; it is "+toString(deltaz);
        messageAST( sbuffer, 2 );
        return( -3 );
    }

    //  remember that nwobble must be at least one even if there is no phonon wobble
    if( nwobble < 1 ) nwobble = 1;

    /* convert FWHM to standard deviation 
            by dividing by 2*sqrt(2*ln(2)) 
        --  monte carlo of position does NOT converge very well for source size so
            do NOT use this, but leave the code here in case I ever figure out a better way */
    //sourcesize = sourceFWHM / 2.354820045;

    //  see if confocal needed
    doConfocal = xFALSE;
    for( i=0; i<ndetect; i++)
            if( CONFOCAL == collectorMode[i]  ) doConfocal = xTRUE;

#ifdef USE_OPENMP
    /*  force LUT init. to avoid redundant init in parallel form */ 
    double rsq = 0.5;  /* arbitrary position */ 
    double vz;  
    for( i=0; i<natom; i++) vz =  vzatomLUT( Znum[i], rsq );
#endif

    if( lwobble == 0 ) {
        if( 0 != lverbose ) {
            sbuffer = "Sorting atoms by depth...";
            messageAST( sbuffer, 0 );
        }
        sortByZ( xa, ya, za, occ, Znum, natom );
    }
    /* to add random offsets */
    xa2.resize( natom );
    ya2.resize( natom );
    za2.resize( natom );
    Znum2.resize( natom );
    occ2.resize( natom );

    /*  check that requested probe size is not bigger 
        than transmission function size (or too small)
    */
    if( (nxprobe > nx) || (nxprobe < 2) ) {
        nxprobe = nx;
        sbuffer = "Probe size reset to nx= " + toString( nxprobe);
        messageAST( sbuffer, 0 );
    }

    if( (nyprobe > ny) || (nyprobe < 2) ) {
        nyprobe = ny;
        sbuffer = "probe size reset to ny= " + toString( nyprobe );
        messageAST( sbuffer, 0 );
    }

    /*  calculate spatial frequencies for future use
        (one set for transmission function and one for probe
        wavefunction)
    NOTE: zero freg is in the bottom left corner and
        expands into all other corners - not in the center
        this is required for FFT - don't waste time rearranging

    remember : the x values are the same for both sets
    
    x2, y2 are used for confocal
    
    */

    kx.resize( nx );
    ky.resize( ny );
    kx2.resize( nx );
    ky2.resize( ny );
    xp.resize( nx );
    yp.resize( ny );

    freqn( kx, kx2, xp, nx, ax );
    freqn( ky, ky2, yp, ny, by );

    kxp.resize( nxprobe );
    kyp.resize( nyprobe );
    kxp2.resize( nxprobe );
    kyp2.resize( nyprobe );
    
    freqn( kxp, kxp2, xp, nxprobe, ax*((double)nxprobe)/nx );
    freqn( kyp, kyp2, yp, nyprobe, by*((double)nyprobe)/ny );

    // impose anti-aliasing bandwidth limit on transmission functions

    sum = ((double)nx)/(2.0*ax);
    k2maxp = ((double)ny)/(2.0*by);
    if( sum < k2maxp ) k2maxp = sum;
    k2maxt = k2maxp;  //  for trans. function potential in cuAtompot()
    k2maxp= BW * k2maxp;
    //printf("Bandwidth limited to a real space resolution of %f Angstroms\n",   //???
    //                 1.0F/k2maxp);
    //printf("   (= %.2f mrad)  for symmetrical anti-aliasing.\n",
    //     wavlen*k2maxp*1000.0F);
    k2maxp = k2maxp * k2maxp;  // for probe 
    k2maxt = k2maxt * k2maxt;  // for trans function potential

    /*  initialize propagator */
    cprop.resize( nxprobe, nyprobe );
    tctx = 2.0 * tan(ctiltx);
    tcty = 2.0 * tan(ctilty);
    scale = pi * deltaz;
    for( ix=0; ix<nxprobe; ix++) {
        wx = ( kxp2[ix]*wavlen - kxp[ix]*tctx );
        for( iy=0; iy<nyprobe; iy++) {
            if( (kxp2[ix] + kyp2[iy]) < k2maxp ) {
                w = scale * ( wx + kyp2[iy]*wavlen - kyp[iy]*tcty );
                cprop.re(ix,iy) = (float)  cos(w);
                cprop.im(ix,iy) = (float) -sin(w);
            } else {
                cprop.re(ix,iy) = 0.0F;
                cprop.im(ix,iy) = 0.0F;
            }  //  end if( kx2[ix]... 
        } // end for(iy..) 
    } // end for(ix..) 

    /*   calculate number of pixels in the probe and obj. apert. */
    k2maxa = apert1 /wavlen;
    k2maxa = k2maxa *k2maxa;
    k2maxb = apert2 /wavlen;
    k2maxb = k2maxb * k2maxb;
    
    nbeamp = nbeampo = 0;
    for( iy=0; iy<nyprobe; iy++)
    for( ix=0; ix<nxprobe; ix++) {
        k2 = kyp2[iy] + kxp2[ix];
        if( k2 < k2maxp ) nbeamp++;
        if( (k2 >= k2maxa) && (k2 <= k2maxb) ) nbeampo++;
    }

    //  output this in outer calling program if wanted - not every time here
    //sbuffer = "Number of symmetrical anti-aliasing beams in probe= "
    //       + toString(nbeamp);
    //messageAST( sbuffer, 0 );
    //sbuffer = "Number of beams in probe aperture= " + toString(nbeampo);
    //messageAST( sbuffer, 0 );

    if( nbeamp < 200 ) {
        sbuffer = "WARNING: the probe is under sampled, this is a bad idea...";
        messageAST( sbuffer, 1 );
    }
    if( nbeampo < 100 ) {
        sbuffer = "WARNING: the probe aperture is under sampled, this is a bad idea...";
        messageAST( sbuffer, 2 );
        exit( EXIT_FAILURE );
    }

    /*  convert aperture dimensions */

    k2min.resize( ndetect );
    k2max.resize( ndetect );

    for( idetect=0; idetect<ndetect; idetect++) {
        if( (ADF == collectorMode[idetect]) || (ADF_SEG == collectorMode[idetect]) ) {
            k2max[idetect] = almax[idetect]/wavlen;
            k2max[idetect] = k2max[idetect] * k2max[idetect];
            k2min[idetect] = almin[idetect]/wavlen;
            k2min[idetect] = k2min[idetect] * k2min[idetect];
        } else if( CONFOCAL == collectorMode[idetect] ) {
            k2max[idetect] = almax[idetect] * almax[idetect];
            k2min[idetect] = almin[idetect] * almin[idetect];
        }
    }

    /*  init the min/max record of total integrated intensity */

    totmin =  10.0;
    totmax = -10.0;
    detect  = new3D<double>( nThick, ndetect, nprobes,
        "detect" );
    sums.resize( nprobes ); 

    // allocate probe wave function and transmission function arrays

#ifndef AST_USE_CUDA
    probe = new cfpix[ nprobes ];
    if( NULL == probe ) {
        sbuffer = "Cannot allocate probe array";
        messageAST( sbuffer, 2 );
        exit( EXIT_FAILURE );
    }
    ix = probe[0].resize(nxprobe, nyprobe );
    if( ix < 0 ) {
        sbuffer = "Cannot allocate probe array storage";
        messageAST( sbuffer, 2 );
        exit( EXIT_FAILURE );
    }
    probe[0].init();
    if( nprobes > 1 ) for( ip=1; ip<nprobes; ip++){
        ix = probe[ip].resize(nxprobe, nyprobe );
        if( ix < 0 ) {
            sbuffer = "Cannot allocate probe array storage";
            messageAST( sbuffer, 2 );
            exit( EXIT_FAILURE );    //  should do something better here ??
        }
        probe[ip].copyInit( probe[0] );
    }

#elif defined(AST_USE_CUDA)

    // ------- cuda setup -------------------

    cudaSetDevice( 0 );  // select GPU number 0
    cudaDeviceReset();   // reset

    //  make cuda plan on device for probe size and transmission size
    if( cufftPlan2d( &cuplanP, nxprobe, nyprobe, CUFFT_C2C) != CUFFT_SUCCESS) {
    sbuffer = "cuda error; unable to create plan P\n" ; 
        messageAST( sbuffer, 0 );
        return( -1 );
    }
    if( cufftPlan2d( &cuplanT, nx, ny, CUFFT_C2C) != CUFFT_SUCCESS) {
    sbuffer = "cuda error; unable to create plan T\n" ; 
        messageAST( sbuffer, 0 );
        return( -1 );
    }

   //   allocate trans memory on the device 
    cudaMalloc( (void**)&Dtrans, sizeof(cufftComplex)*nx*ny);
    checkCudaErr( "failed to allocate cuda array Dtrans" );
    cudaMalloc( (void**)&DpotnR, sizeof(float)*nx*ny );
    checkCudaErr( "failed to allocate cuda array Dtrans" );
    HpotnR = new float[ nx*ny ];

    //   allocate CBED memory on host and device 
    Hcbed = new float[ nxprobe*nyprobe ];
    cudaMalloc( (void**)&Dcbed, sizeof(float)*nxprobe*nyprobe);
    checkCudaErr( "failed to allocate cuda array Dcbed" );

    //  allocate nprobe * probes on GPU (very big!)
    cudaMalloc( (void**)&Dprobe, sizeof(cufftComplex)*nxprobe*nyprobe * nprobes );
    checkCudaErr( "failed to allocate device memory for Dprobe");
    //  allocate scratch array on device to build probes
    cudaMalloc( (void**)&Dprobe1, sizeof(cufftComplex)*nxprobe*nyprobe );
    checkCudaErr( "failed to allocate device memory for Dprobe");

    //  allocate kxp[],kyp[], kxp2[], kyp2[]  on device
    cudaMalloc( (void**)&Dkxp, sizeof(float)*nxprobe );
    checkCudaErr( "failed to allocate device memory for Dkxp");
    cudaMalloc( (void**)&Dkyp, sizeof(float)*nyprobe );
    checkCudaErr( "failed to allocate device memory for Dkyp");
    cudaMalloc( (void**)&Dkxp2, sizeof(float)*nxprobe );
    checkCudaErr( "failed to allocate device memory for Dkxp2");
    cudaMalloc( (void**)&Dkyp2, sizeof(float)*nyprobe );
    checkCudaErr( "failed to allocate device memory for Dkyp2");

    // ------ caclulate spatial freq on GPU
    int threadsPerBlock, blocksPerGrid;
    threadsPerBlock = 512;  // max thread/block usually 1024 =dev specific
    blocksPerGrid = (nxprobe + threadsPerBlock -1 )/threadsPerBlock;  //  round up
    cuFreq<<<blocksPerGrid, threadsPerBlock>>>( Dkxp, Dkxp2, nxprobe, ax*((double)nxprobe)/nx );
    checkCudaErr( "cuFreq failed for Dkxp");
    blocksPerGrid = (nyprobe + threadsPerBlock -1 )/threadsPerBlock;  //  round up
    cuFreq<<<blocksPerGrid, threadsPerBlock>>>( Dkyp, Dkyp2, nyprobe, by*((double)nyprobe)/ny );
    checkCudaErr( "cuFreq failed for Dkyp");

    //  allocate kx[],ky[] kx2[], ky2[] on device
    cudaMalloc( (void**)&Dkx, sizeof(float)*nx );
    checkCudaErr( "failed to allocate device memory for Dkx");
    cudaMalloc( (void**)&Dky, sizeof(float)*ny );
    checkCudaErr( "failed to allocate device memory for Dky");
    cudaMalloc( (void**)&Dkx2, sizeof(float)*nx );
    checkCudaErr( "failed to allocate device memory for Dkx2");
    cudaMalloc( (void**)&Dky2, sizeof(float)*ny );
    checkCudaErr( "failed to allocate device memory for Dky2");

    //  ------ caclulate spatial freq on GPU
    blocksPerGrid = (nx + threadsPerBlock -1 )/threadsPerBlock;  //  round up
    cuFreq<<<blocksPerGrid, threadsPerBlock>>>( Dkx, Dkx2, nx, ax );
    checkCudaErr( "cuFreq failed for Dkx");
    blocksPerGrid = (ny + threadsPerBlock -1 )/threadsPerBlock;  //  round up
    cuFreq<<<blocksPerGrid, threadsPerBlock>>>( Dky, Dky2, ny, by );
    checkCudaErr( "cuFreq failed for Dky");

#ifdef AST_USE_CUATOMPOT 
    if( cufftPlan2d( &cuplanTc2r, nx, ny, CUFFT_C2R) != CUFFT_SUCCESS) {
    sbuffer = "cuda error; unable to create plan Tc2r\n" ; 
        messageAST( sbuffer, 0 );
        return( -1 );
    }

    //  allocate fparams on on host and device
    float fex = featom( 12, 0.1 );  //  force init of fparams[][]
    const int NPMAX=   12;  // number of parameters for each Z
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
        messageAST( sbuffer, 0 );
        return( -1 );
    }
    //  make associated device memory
    cudaMalloc( (void**)&Dspec, 4*natomin*sizeof(float) );
    checkCudaErr( "failed to allocate device memory for Dspec");

    cudaMalloc( (void**)&DpotnC, sizeof(cufftComplex)*nx*(ny/2 + 1) );
    checkCudaErr( "failed to allocate cuda array DpotnC" );
#endif

    //  allocate propagator on host and device
    Hprop = new cufftComplex [nxprobe*nyprobe];
    cudaMalloc( (void**)&Dprop, sizeof(cufftComplex)*nxprobe*nyprobe );
    checkCudaErr( "failed to allocate device memory for Dprop");

    //  calculate 2D propagator function
    scale = pi * deltaz;
    double scale2 = 1.0/(nxprobe*nyprobe);   // for FFT scaling

    for( ix=0; ix<nxprobe; ix++) {
        for( iy=0; iy<nyprobe; iy++) {
            Hprop[iy + ix*nyprobe].x = cprop.re(ix,iy) * scale2;
            Hprop[iy + ix*nyprobe].y = cprop.im(ix,iy) * scale2;

    } }  //  end for( iy=0... for( ix=0...
    //--- copy to device ---------
    cudaMemcpy( Dprop, Hprop, nxprobe*nyprobe * sizeof(cufftComplex),
             cudaMemcpyHostToDevice );
    checkCudaErr( "cannot copy prop from host to device");
    delete [] Hprop;

    //  allocate host memory for one probe
    Hprobe = new cufftComplex[ nxprobe*nyprobe ];
    if( NULL == Hprobe ) {
        sbuffer = "cannot allocate cuda host array Hprobe\n";
        messageAST( sbuffer, 0 );
        return( -1 );
    }

    //  allocate host memory for transmission function
    Htrans = new cufftComplex[ nx*ny ];
    if( NULL == Htrans ) {
        sbuffer = "cannot allocate cuda host array Htrans\n";
        messageAST( sbuffer, 0 );
        return( -1 );
    }

    //  allocate host memory for detector sums
    Hsums= new double[ nyprobe ];
    if( NULL == Hsums ) {
        sbuffer = "cannot allocate cuda host array Hsums\n";
        messageAST( sbuffer, 0 );
        return( -1 );
    }
    cudaMalloc( (void**)&Dsums, sizeof(double)*nyprobe );
    checkCudaErr( "failed to allocate device memory for Dsums");

    //  scratch pix for confocal mode
    probe0.resize( nxprobe, nyprobe );
    probe0.init( );

#endif

    trans.resize( nx, ny );
    trans.init();

    if( lpacbed == xTRUE ) {
        for( ix=0; ix<nxprobe; ix++) for( iy=0; iy<nyprobe; iy++)
                pacbedPix[ix][iy] = 0;
    }

/* ------------- start here for a full image output -------------- */
/*
  do one whole line at once NOT the whole image (which may be huge)
*/
    if( l1d == 0 ) {
       sbuffer = "output file size in pixels is " + toString(nxout) +" x " 
               + toString(nyout);
       messageAST( sbuffer, 0 );
       if( nprobes != nyout ) {
           sbuffer = "Error, nprobes= "+ toString(nprobes)+" must be the same as nyout= "
               +toString(nyout)+", in image mode.";
           messageAST( sbuffer, 2 );
           exit( 0 );
       }
       /* double up first index to mimic a 4D array */
       for( i=0; i<(nThick*ndetect); i++) {
           for( ix=0; ix<nxout; ix++)
           for( iy=0; iy<nyout; iy++)
            pixr[i][ix][iy] = 0.0F;
       }

       /*  iterate the multislice algorithm proper for each
           position of the focused probe */

       if( nxout > 1 ) dx = (xf-xi)/((double)(nxout-1));
       else dx = 1.0;
       if( nyout > 1 ) dy = (yf-yi)/((double)(nyout-1));
       else dy = 1.0;
       x.resize( nprobes );
       y.resize( nprobes );

        /*  add random thermal displacements 
               scaled by temperature if requested 
            remember that initial wobble is at 300K for
               each direction */

        for( iwobble=0; iwobble<nwobble; iwobble++) {
            if( lwobble == 1 ){
                scale = (float) sqrt(temperature/300.0) ;
                for( i=0; i<natom; i++) {
                    xa2[i] = xa[i] + 
                        (float)(wobble[i]*rangauss(iseed)*scale);
                    ya2[i] = ya[i] + 
                        (float)(wobble[i]*rangauss(iseed)*scale);
                    za2[i] = za[i] + 
                            (float)(wobble[i]*rangauss(iseed)*scale);
                    occ2[i] = occ[i];
                    Znum2[i] = Znum[i];
                }
                sortByZ( xa2, ya2, za2, occ2, Znum2, natom );
                sbuffer = "configuration # " + toString( iwobble+1 );
                messageAST( sbuffer, 0 );
                sbuffer = "The new range of z is "
                    + toString(za2[0]) + " to " + toString( za2[natom-1] );
                messageAST( sbuffer, 0 );
            } else for( i=0; i<natom; i++) {
                xa2[i] = xa[i];
                ya2[i] = ya[i];
                za2[i] = za[i];
                occ2[i] = occ[i];
                Znum2[i] = Znum[i];
            }
            zmin = za2[0];  /* reset zmin/max after wobble */
            zmax = za2[natom-1];
    
            for( ix=0; ix<nxout; ix++) {
    
                for( iy=0; iy<nyout; iy++) {
                    x[iy] = xi + dx * ((double) ix);
                        //  + sourcesize * rangauss(iseed);  - does not converge well
                    y[iy] = yi + dy * ((double) iy);
                        //  + sourcesize * rangauss(iseed);  - does not converge well
                    x[iy] = periodic( x[iy], ax );   /* put back in supercell */
                    y[iy] = periodic( y[iy], by );   /* if necessary */
                }

                sbuffer =  "calculate line " + toString(ix);  //  new 3-jun-2015 ejk
                messageAST( sbuffer, 0 );

                STEMsignals( x, y, nyout, param, multiMode, detect, ndetect, 
                    ThickSave, nThick, sums, collectorMode, phiMin, phiMax );
                for( iy=0; iy<nyout; iy++) {
                    if( sums[iy] < totmin ) totmin = sums[iy];
                    if( sums[iy] > totmax ) totmax = sums[iy];
                    for( it=0; it<nThick; it++){
                        for( idetect=0; idetect<ndetect; idetect++)
                        pixr[idetect + it*ndetect][ix][iy] += (float)
                            (detect[it][idetect][iy]/((double)nwobble));
                    }
                    if( sums[iy] < 0.9) {
                        sbuffer = "Warning integrated intensity too small, = "
                           + toString(sums[iy])+" at "+toString(x[iy])+", "+toString(y[iy]);
                        messageAST( sbuffer, 0 );
                    }
                    if( sums[iy] > 1.1) {
                        sbuffer =  "Warning integrated intensity too large, = "
                           + toString(sums[iy])+" at "+toString(x[iy])+", "+toString(y[iy]);
                        messageAST( sbuffer, 0 );
                    }
                }

                /*   sum position averaged CBED if requested 
                     - assume probe still left from stemsignal()  */
#ifndef AST_USE_CUDA
                if( lpacbed == xTRUE ) {
                    for( iy=0; iy<nyout; iy++) {
                        for( ix2=0; ix2<nxprobe; ix2++)
                        for( iy2=0; iy2<nyprobe; iy2++) {
                            prr = probe[iy].re(ix2,iy2);
                            pri = probe[iy].im(ix2,iy2);
                            pacbedPix[ix2][iy2] += (prr*prr + pri*pri);
                       }
                    }
                }   /*  end if( lpacbed.... */
#elif defined(AST_USE_CUDA)
                if( lpacbed == xTRUE ) {
                   for( iy=0; iy<nyout; iy++) {
                        //--- copy to device ---------
                    cudaMemcpy( Hprobe, &Dprobe[iy*nx*ny], nxprobe*nyprobe*sizeof(cufftComplex),
                    cudaMemcpyDeviceToHost );
                    checkCudaErr( "cannot copy probe device to host");
                        for( ix2=0; ix2<nxprobe; ix2++)
                        for( iy2=0; iy2<nyprobe; iy2++) {
                            prr = Hprobe[iy2 + ix2*nyprobe].x;   //  probe[iy].re(ix2,iy2);
                            pri = Hprobe[iy2 + ix2*nyprobe].y;   //  probe[iy].im(ix2,iy2);
                            pacbedPix[ix2][iy2] += (prr*prr + pri*pri);
                       }
                    }
                 }   /*  end if( lpacbed.... */

#endif            
            } /* end for(ix...) */
    
        } /* end for(iwobble... ) */

        /*  find range to output data files  */
        for( it=0; it<nThick; it++)
        for( i=0; i<ndetect; i++) {
            rmin[it][i] = rmax[it][i] = pixr[i+it*ndetect][0][0];
            for( ix=0; ix<nxout; ix++)
            for( iy=0; iy<nyout; iy++) {
                temp = pixr[i+it*ndetect][ix][iy];
                if( temp < rmin[it][i] )rmin[it][i] = (float) temp;
                if( temp > rmax[it][i] )rmax[it][i] = (float) temp;
            }
        }
        if( lpacbed == xTRUE ) {
            invert2D( pacbedPix, nxprobe, nyprobe );  /*  put zero in middle */
         }

    /* ------------- start here for 1d line scan ---------------- */

    } else { // the only other posibility is; if ( l1d == 1 ) {

       if( lpacbed == xTRUE ) {
          sbuffer = "warning: cannot do pos. aver. CBED in 1d";
          messageAST( sbuffer, 0 );
       }
       if( nxout > 1 ) {
               sbuffer="nxout must be 1 in 1D mode but is "+toString(nxout);
               messageAST( sbuffer, 0 );
       }
       if( nyout > 1 ) dx = (xf-xi)/((double)(nyout-1));
       else dx = 1.0;
       if( nyout > 1 ) dy = (yf-yi)/((double)(nyout-1));
       else dy = 1.0;
       x.resize( nprobes );
       y.resize( nprobes );
       for( ip=0; ip<nyout; ip++) {
            for( it=0; it<nThick; it++)
            for( idetect=0; idetect<ndetect; idetect++)
                    pixr[idetect+it*ndetect][0][ip] = 0.0F;
       }

       /*  add random thermal displacements scaled by temperature
            if requested 
        remember that initial wobble is at 300K for each direction */

       for( iwobble=0; iwobble<nwobble; iwobble++) {
    
            if( lwobble == 1 ){
                scale = (float) sqrt(temperature/300.0) ;
                for( i=0; i<natom; i++) {
                    xa2[i] = xa[i] + 
                        (float)(wobble[i]*rangauss(iseed)*scale);
                    ya2[i] = ya[i] + 
                        (float)(wobble[i]*rangauss(iseed)*scale);
                    za2[i] = za[i] + 
                        (float)(wobble[i]*rangauss(iseed)*scale);
                    occ2[i] = occ[i];
                    Znum2[i] = Znum[i];
                }
                sbuffer = "configuration # " + toString( iwobble+1 );
                messageAST( sbuffer, 0 );
                sortByZ( xa2, ya2, za2, occ2, Znum2, natom );
                sbuffer = "The new range of z is "
                    + toString(za2[0]) + " to " + toString( za2[natom-1] );
                messageAST( sbuffer, 0 );
            } else for( i=0; i<natom; i++) {
                xa2[i] = xa[i];
                ya2[i] = ya[i];
                za2[i] = za[i];
                occ2[i] = occ[i];
                Znum2[i] = Znum[i];
            }
            zmin = za2[0];      /* reset zmin/max after wobble */
            zmax = za2[natom-1];
            for( ip=0; ip<nyout; ip++) {
                x[ip] = xi + dx * ((double)ip);
                            //  + sourcesize * rangauss(iseed);  - does not converge well
                y[ip] = yi + dy * ((double)ip);
                            //  + sourcesize * rangauss(iseed);  - does not converge well
                x[ip] = periodic( x[ip], ax );   // put back in supercell
                y[ip] = periodic( y[ip], by );   // if necessary
            }
         
            STEMsignals( x, y, nprobes, param, multiMode, detect, ndetect, 
                ThickSave, nThick, sums, collectorMode, phiMin, phiMax );
            for( ip=0; ip<nprobes; ip++) {
                if( sums[ip] < totmin ) totmin = sums[ip];
                if( sums[ip] > totmax ) totmax = sums[ip];
                for( it=0; it<nThick; it++){
                for( idetect=0; idetect<ndetect; idetect++)
                   //  nwobble should be small so its prob. safe to sum into single prec. var.
                   pixr[idetect+it*ndetect][0][ip] += 
                        (float) ( detect[it][idetect][ip]/((double)nwobble) );
                }
                if( sums[ip] < 0.9) {
                    sbuffer = "Warning integrated intensity too small, = "
                           + toString(sums[ip])+" at "+toString(x[ip])+", "+toString(y[ip]);
                    messageAST( sbuffer, 0 );
                }
                if( sums[ip] > 1.1) {
                    sbuffer =  "Warning integrated intensity too large, = "
                           + toString(sums[ip])+" at "+toString(x[ip])+", "+toString(y[ip]);
                    messageAST( sbuffer, 0 );
                }
            }

       }  /* end for(iwobble... */


    } /* end if( l1d.. ) */

    //----------- end:  free scratch arrays and exit --------------------

    delete3D<double>( detect, nThick, ndetect );

#ifndef AST_USE_CUDA
    delete [] probe;
#elif defined(AST_USE_CUDA)
    cufftDestroy( cuplanP );
    cufftDestroy( cuplanT );
    cufftDestroy( cuplanTc2r );
    cudaFree( Dprobe );
    cudaFree( Dprobe1 );
    cudaFree( Dkxp );
    cudaFree( Dkyp );
    cudaFree( Dkxp2 );
    cudaFree( Dkyp2 );
    cudaFree( Dkx );
    cudaFree( Dky );
    cudaFree( Dkx2 );
    cudaFree( Dky2 );
    cudaFree( Dprop );
    cudaFree( Dtrans );
    cudaFree( Dcbed );
    cudaFree( Dsums );
#ifdef AST_USE_CUATOMPOT 
    cudaFree( Dfparams );
    cudaFree( Dspec );
    cudaFree( DpotnC );
    delete [] Hspec;
#endif
    cudaFree( DpotnR );
    delete [] Hprobe;
    delete [] Htrans;
    delete [] Hcbed;
    delete [] Hsums;
    delete [] HpotnR;
#endif

    return( + 1 );
         
}  // end autostem::calculate()

/* -------------------  CountBeams() -------------------

   count the number of beams (Fourier coefficients) in the
   transmission function and the probe for informational purposes
   - to communicate to the main calling program
   -  must be recalculated here (should be identical to that
      used in calculate()
  input:
        param[] = array of parameters

  output:
        nbeamp  = number of symmetrical beams in the probe
        nbeampo = number of beams in the probe aperture
        res  = bandwidth limited resolution in Angstroms
        almax = alpha max in radians
*/
void autostem::CountBeams( vectorf &param, int &nbeamp, int &nbeampo, float &res, float &almax )
{
    //  use all local array variables so I don't accidentally
    //     disturb the main calculation
    int ix, iy, nx1, ny1, nxprobe1, nyprobe1;
    float ax, by, wavl, apert1, apert2;
    double sum, k2maxp, k2maxa, k2maxb, k2;

    vectorf xp1, yp1, kxp1, kyp1, kxp21, kyp21;

    ax = param[ pAX ];                  // supercell size in Ang.
    by = param[ pBY ];

    nx1 = ToInt( param[ pNX ] );        // size of transmission function (in pixels)
    ny1 = ToInt( param[ pNY ] );

    nxprobe1 = ToInt( param[ pNXPRB ] );        // probe size in pixels
    nyprobe1 = ToInt( param[ pNYPRB ] );

    wavl = (float) wavelength( param[ pENERGY ] );      //  electron wavelength in Ang.

    apert1 = param[ pOAPMIN ];  // obj. apert size in radians (min, max for annular apert)
    apert2 = param[ pOAPERT ];

    //  transmission function sampling
    // impose anti-aliasing bandwidth limit on transmission functions
    sum = ((double)nx1)/(2.0*ax);
    k2maxp = ((double)ny1)/(2.0*by);
    if( sum < k2maxp ) k2maxp = sum;
    k2maxp= BW * k2maxp;
    res = (float) ( 1.0/k2maxp);
    almax = (float) (wavl*k2maxp);
    k2maxp = k2maxp * k2maxp;

    //  probe sampling
    //  the only way to do this is to do part of the calculation and
    //    throw away the results (kx,ky etc.) - but only a small bit wasted
    xp1.resize( nxprobe1 );
    yp1.resize( nyprobe1 );
    kxp1.resize( nxprobe1 );
    kyp1.resize( nyprobe1 );
    kxp21.resize( nxprobe1 );
    kyp21.resize( nyprobe1 );
    
    freqn( kxp1, kxp21, xp1, nxprobe1, ax*((double)nxprobe1)/nx1 );
    freqn( kyp1, kyp21, yp1, nyprobe1, by*((double)nyprobe1)/ny1 );

    /*   calculate number of pixels in the probe and obj. apert. */
    k2maxa = apert1 /wavl;
    k2maxa = k2maxa *k2maxa;
    k2maxb = apert2 /wavl;
    k2maxb = k2maxb * k2maxb;
    
    nbeamp = nbeampo = 0;
    for( iy=0; iy<nyprobe1; iy++)
    for( ix=0; ix<nxprobe1; ix++) {
        k2 = kyp21[iy] + kxp21[ix];
        if( k2 < k2maxp ) nbeamp++;
        if( (k2 >= k2maxa) && (k2 <= k2maxb) ) nbeampo++;
    }

};  // end autostem::CountBeams()


/* -------------------  messageAST() -------------------
   message output
   direct all output message here to redirect to the command line
   or a GUI status line or message box when appropriate

   msg[] = character string with message to disply
   level = level of seriousness
        0 = simple status message
        1 = significant warning
        2 = possibly fatal error
*/
void autostem::messageAST( std::string &smsg,  int level )
{
        messageSL( smsg.c_str(), level );  //  just call slicelib version for now

}  // end autostem::messageAST()

/*------------------------ periodic() ---------------------*/
/*
     make probe positions periodic in the supercell
     in case some wobble off the edge with source size of user excess

    pos = input position (x or y);
    size = supercell size ( 0 to size)

    return positive value  0 <= x < size
*/
double autostem::periodic( double pos, double size )
{
    double x=pos;
    while( x < 0 ) x += size;
    x = fmod( x, size );
    return( x );
}

/*------------------------ STEMsignals() ---------------------*/
/*

  NOTE: this is NOT the same as STEMsignal() in stemslice.c

  subroutine to calculate the stem signal arising from a given
  probe position

  iterate the multislice algorithm proper for each position of
  the focused probe

  This version uses massive amounts of memory to avoid
  recalculating the transmission functions more than necessary

   note zero freg is in the bottom left corner and
     expands into all other corners - not in the center
     this is required for fft - don't waste time rearranging

     change to propagate thru whole unit cell not just
     range of atoms on 14-sep-2006 ejk
     add multipole aberrations 9-may-2011 ejk
     change to cfpix for probe and trans 10-nov-2012 ejk

  x[],y[]     = real positions of the incident probe
  npos        = int number of positions
  param[]     = parameters of probe
  multiMode   = flag to add multipole aberrations
  detect[][][]= real array to get signal into each detector
            for each probe position and thickness
  ndetect     = number of detector geometries
  ThickSave[] = thicknesses at which to save data (other than the last)
  nThick      = number of thickness levels (including the last)
  sum         = real total integrated intensity
  
  the assumed global variables are:
  
  nxprobe,nyprobe   = int size of probe wavefunction in pixels
  nx,ny         = int size of transmission function in pixels
  trans                  = float complex transmission function
  propxr[][], propxi[][] = float real,imag propagator vs x
  propyr[][], propyi[][] = float real,imag propagator vs y
  ax,by,cz      = float unit cell size in Angs
  kxp[], kyp[]      = float spatial frequencies vs x, y
  kxp2[], kyp2[]    = float square of kxp[], kyp[]
  xp[], yp[]        = float real space positions in probe (confocal)
  apert1, apert2    = double min,max objective aperture (in radians)
  k2maxp            = double max spatial freq of probe squared
  pi                = double constant PI
  wavlen            = double electron wavelength in Angs
  df                = double defocus (in Ang)
  Cs3,Cs5           = double spherical aberration (in Ang)

  xa[],ya[],za[]    = atom coordinates
  occ[]         = atomic occupancy
  Znum[]        = atomic numbers
  natom         = number of atoms
  deltaz        = slice thickness
  v0            = beam energy
  nbeamt        = number of beams in transmission function
  zmin, zmax    = range of z coord. of the atoms
  nslice        = number of slices
  doConfocal    = flag indicating confocal is needed
  k2min[], k2max[] = detector angles in 1/Ang
  
    NOTE:  too many thing come in as globals, but...

*/
#ifndef AST_USE_CUDA
void autostem::STEMsignals( vectord &x, vectord &y, int npos, vectorf &p,
         int multiMode, double ***detect, int ndetect,
         vectord &ThickSave, int nThick, vectord &sum, vectori &collectorMode,
         vectord &phiMin, vectord &phiMax )
{
    int ix, iy, ixt, iyt, idetect,  ixmid, iymid;
    int istart, na, ip, i, it;

    long nxl, nyl;

    float scale, prr, pri, tr, ti;

    double  chi0, chi1, k2maxa, k2maxb,
        w, k2, phi, phirms, alx, aly;
    double sum0, sum1, delta, zslice, totalz;

    vectori ixoff, iyoff;
    vectord xoff, yoff;
    
    /* extra for confocal */
    float hr, hi;
    double chi2C, chi3C, k2maxaC, k2maxbC, r2, rx2;
    cfpix cpix;            /* complex confocal image */

    /* ------ make sure x,y are ok ------ */

    for( ip=0; ip<npos; ip++) {
        if( (x[ip] < 0.0) || (x[ip] > ax) ||
            (y[ip] < 0.0) || (y[ip] > by) ) {
            sum[ip] = -1.2345;
            sbuffer = "bad x,y in STEMsignals() = \n" + toString(x[ip]) + toString(y[ip]);
            messageAST( sbuffer, 0 );
            return;
        }
    }

    /*  generate a probe at position x,y which must be inside the
        0->ax and 0->by region
        normalize to have unit integrated intensity
        probe position (x,y) = (0,0) = lower left corner

    NOTE: The probe wavefunction is nxprobe x nyprobe which may
        be smaller than the transmission function.
        Position the probe near the center of the nxprobe x nyprobe
        region because it does not wrap around properly
        (i.e. manually wrap around in the nx x ny region if
         the probe is near a boundary of the nx x ny region))
    */

    ixmid = nxprobe/2;
    iymid = nyprobe/2;
    chi1 = pi * wavlen;
    k2maxa = apert1 /wavlen;
    k2maxa = k2maxa *k2maxa;
    k2maxb = apert2 /wavlen;
    k2maxb = k2maxb * k2maxb;
    
    /* extra for confocal */
    chi2C = 0.5 * Cs3C *wavlen*wavlen;
    chi3C = Cs5C * wavlen*wavlen*wavlen*wavlen /3.0;
    k2maxaC = apert1C /wavlen;
    k2maxaC = k2maxaC *k2maxaC;
    k2maxbC = apert2C /wavlen;
    k2maxbC = k2maxbC * k2maxbC;

    ixoff.resize( npos );
    iyoff.resize( npos );
    xoff.resize( npos );
    yoff.resize( npos );

    /* ------- calculate all of the probe wave functions at once ------
        to reuse the transmission functions which takes a long
        time to calculate*/

/*  paralleling this loop has little effect */
#pragma omp parallel for private(ix,iy,sum0,k2,w,chi0,scale,tr,ti,alx,aly) 
    for( ip=0; ip<npos; ip++) {
        ixoff[ip] = (int) floor( x[ip]*((double)nx) / ax ) - ixmid;
        xoff[ip]  = x[ip] - ax*((double)ixoff[ip])/((double)nx);

        iyoff[ip] = (int) floor( y[ip]*((double)ny) / by ) - iymid;
        yoff[ip]  = y[ip] - by*((double)iyoff[ip])/((double)ny);

        sum0 = 0.0;
        for( ix=0; ix<nxprobe; ix++) {
            alx = wavlen * kxp[ix];  /* x component of angle alpha */
            for( iy=0; iy<nyprobe; iy++) {
                aly = wavlen * kyp[iy];  /* y component of angle alpha */
                k2 = kxp2[ix] + kyp2[iy];
                if( (k2 >= k2maxa) && (k2 <= k2maxb) ) {
                    w = 2.*pi* ( xoff[ip]*kxp[ix] + yoff[ip]*kyp[iy] );
                    chi0 = (2.0*pi/wavlen) * chi( p, alx, aly, multiMode );
                    chi0= - chi0 + w;
                    probe[ip].re(ix,iy) = tr = (float) cos( chi0 );
                    probe[ip].im(ix,iy) = ti = (float) sin( chi0 );
                    sum0 += (double) (tr*tr + ti*ti);
                } else {
                    probe[ip].re(ix,iy) = 0.0F;
                    probe[ip].im(ix,iy) = 0.0F;
                }
            }
        }  /* end for( ix... */

        scale = (float) ( 1.0/sqrt(sum0) );  
        probe[ip] *= scale;

    }  /* end for( ip...) */

    /* -------- transmit thru nslice layers ------------------------
        beware ixoff,iyoff must be with one nx,ny
            of the array bounds
    */          


    nxl = (long) nx;
    nyl = (long) ny;
    
    scale = 1.0F / ( ((float)nx) * ((float)ny) );

    zslice = 0.75*deltaz;  /*  start a little before top of unit cell */
    istart = 0;
    nslice = 0;
    it = 0;     /* thickness level index */          

    if( zmax > cz ) totalz = zmax;
        else totalz = cz;
    if( 0 != lverbose ) {
        sbuffer= "specimen range is 0 to "+ toString(totalz) + " Ang.";
        messageAST( sbuffer, 0 );
    }
    
    /* range of unit cell */
    while(  (zslice < (totalz+0.25*deltaz)) || (istart<natom) ) {

       /* find range of atoms for current slice */
       na = 0;
       for(i=istart; i<natom; i++)
            if( za2[i] < zslice ) na++; else break;

       if( 0 != lverbose ) {
           sbuffer= "slice ending at z= "+toString(zslice)+" Ang. with "+toString(na)+" atoms";
           messageAST( sbuffer, 0 );
           //ss.str("");
           //ss << "slice ending at z= " << zslice << " Ang. with " << na << " atoms";
           //messageAST( ss.str(), 0 );
       }

       /* calculate transmission function and bandwidth limit */
       if( na > 0 ) trlayer( xa2, ya2, occ2,
            Znum2, na, istart, (float)ax, (float)by, (float)keV,
            trans, nxl, nyl, &phirms, &nbeamt, (float) k2maxp );

       /*----- one multislice trans/prop cycle for all probes ---- */
#pragma omp parallel for private(ix,iy,ixt,iyt,prr,pri)
       for( ip=0; ip<npos; ip++) {
           /* apply transmission function if there are atoms in this slice */
           if( na > 0 ) {
                probe[ip].ifft();
                for( ix=0; ix<nxprobe; ix++) {
                    ixt = ix + ixoff[ip];
                    if( ixt >= nx ) ixt = ixt - nx;
                    else if( ixt < 0 ) ixt = ixt + nx;
                    for( iy=0; iy<nyprobe; iy++) {
                        iyt = iy + iyoff[ip];
                        if( iyt >= ny ) iyt = iyt - ny;
                        else if( iyt < 0 ) iyt = iyt + ny;
                        prr = probe[ip].re(ix,iy);
                        pri = probe[ip].im(ix,iy);
                        probe[ip].re(ix,iy) =  prr*trans.re(ixt,iyt)
                                              - pri*trans.im(ixt,iyt);  // real
                        probe[ip].im(ix,iy) =  prr*trans.im(ixt,iyt)
                                              + pri*trans.re(ixt,iyt);  // imag
                    } /* end for(iy...) */
                }  /* end for(ix...) */
                probe[ip].fft();
           }
    
            /*  multiplied by the propagator function */
            probe[ip] *= cprop;

        }  /* end  for( ip=... */

        /*  if this is a good thickness level then save the ADF or confocal signals
           - remember that the last level may be off by one layer with
              thermal displacements so special case it
        */
       
        /*  look at all values because they may not be in order */
        for( it = 0; it<nThick; it++ ) 
        if( fabs(ThickSave[it]-zslice)<fabs(0.5*deltaz)) {
            
            if( 0 != lverbose ) {
                sbuffer= "save ADF/confocal signals, thickness level " + toString(it);  // diagnostic
                messageAST( sbuffer, 0 );
            }
    
            /*  loop over all probes again */
#pragma omp parallel for private(ix,iy,idetect,prr,pri,delta,k2,cpix,phi,chi0,hr,hi,sum0,sum1,rx2,r2)
            for( ip=0; ip<npos; ip++) {

                /*  zero sum count */
                sum[ip] = 0.0;
                for(ix=0; ix<ndetect; ix++) detect[it][ix][ip] = 0.0;
   
                /*  sum intensity incident on the ADF detector
                        and calculate total integrated intensity
                    - changed detector limits to >= min and < max
                       so many concentric ADF detectors sum correctly
                       7-jul-2011
                */
                for( ix=0; ix<nxprobe; ix++) {
                    for( iy=0; iy<nyprobe; iy++) {

                        prr = probe[ip].re(ix,iy);
                        pri = probe[ip].im(ix,iy);
                        delta = prr*prr + pri*pri;
                        sum[ip] += delta;

                        k2 = kxp2[ix] + kyp2[iy];
                        phi = atan2( kyp[iy], kxp[ix] );  //  for ADF_SEG detector

                        for( idetect=0; idetect<ndetect; idetect++) {
                            if( ADF == collectorMode[idetect] ) {
                                if( (k2 >= k2min[idetect] ) &&
                                (k2 < k2max[idetect] ) )
                                detect[it][idetect][ip] += delta;
                            } else if( ADF_SEG == collectorMode[idetect] ) {
                                if( (k2 >= k2min[idetect] ) && (k2 < k2max[idetect] )
                                    && (phi >= phiMin[idetect]) && (phi< phiMax[idetect]) )
                                detect[it][idetect][ip] += delta;
                            }
                        }
                    } /* end for(iy..) */
                }  /* end for(ix...) */

                /*  transform back if confocal needed 
                    - use copy of probe so original can continue in use  */
                if( doConfocal == xTRUE ) {
                    /*  allocate/deallocate here so openMP will work 
                        otherwise have to allocate nxprobe cpix arrays 
                        - a littel slow but a lot less memory */
                    cpix.resize(nxprobe,nyprobe);
                    cpix.copyInit(probe[0]);
                    sum0 = 0;
                    for( ix=0; ix<nxprobe; ix++) {
                        for( iy=0; iy<nyprobe; iy++) {
                            k2 = kxp2[ix] + kyp2[iy];
                            if( (k2 >= k2maxaC) && (k2 <= k2maxbC) ) {
                                phi = atan2( ky[iy], kx[ix] );
                                /*  offset defocus by zslice so both lens referenced to 
                                   entrance surface of specimen  */
                                chi0 = chi1*k2* ( (chi2C + chi3C*k2)*k2 - dfC + zslice
                                    + dfa2C*sin( 2.0*(phi-dfa2phiC) ) 
                                    + 2.0F*dfa3C*wavlen*sqrt(k2)*
                                    sin( 3.0*(phi-dfa3phiC) )/3.0 );
                                chi0= - chi0;
                                hr = (float) cos( chi0 );
                                hi = (float) sin( chi0 );
                                prr = probe[ip].re(ix,iy);  // real
                                pri = probe[ip].im(ix,iy);  // imag
                                cpix.re(ix,iy) = prr*hr -pri*hi;
                                cpix.im(ix,iy) = prr*hi +pri*hr;
                                sum0 += prr*prr + pri*pri;
                            } else {
                                cpix.re(ix,iy) = 0.0F;
                                cpix.im(ix,iy) = 0.0F;
                            }
                        }  /*  end for( iy... )  */
                    }  /*  end for( ix... )  */
            
                    cpix.ifft();
                    
                    /* find normalization constant 
                      i.e. correct for constants in the FFT */
                    sum1 = 0.0;
                    for( ix=0; ix<nxprobe; ix++) {
                         for( iy=0; iy<nyprobe; iy++) {
                            prr = cpix.re(ix,iy);
                            pri = cpix.im(ix,iy);
                            sum1 += prr*prr + pri*pri;
                        }
                   }
            
                   /* integrate over real space detector and normalize */
                   for( ix=0; ix<nxprobe; ix++) {
                           rx2 = xoff[ip] - xp[ix];
                           rx2 = rx2*rx2;
                           for( iy=0; iy<nyprobe; iy++) {
                               r2 = yoff[ip] - yp[iy];
                               r2 = rx2 + r2*r2;
                               prr = cpix.re(ix,iy);
                               pri = cpix.im(ix,iy);
                               delta = prr*prr + pri*pri;
                               for( idetect=0; idetect<ndetect; idetect++) {
                                   if( CONFOCAL == collectorMode[idetect] ) {
                                     if( (r2 >= k2min[idetect] ) &&
                                               (r2 < k2max[idetect] ) )
                                         detect[it][idetect][ip] += delta*(sum0/sum1);
                                   }
                               }
                           }  /* end for( iy... )*/
                   }  /*  end for( ix....) */

                }  /* end if( doConfocal==TRUE) */
            
            }  /* end for( ip.. */

        }  /* end if( ((it...*/

        sbuffer = "xxx bad";
        nslice++;
        zslice += deltaz;
        istart += na;

    }  /* end while( istart...) */

    return;

}// end autostem::STEMsignals() - openMP version

#elif defined(AST_USE_CUDA)

//-------  checkCudaErr() -------------
//  check last CUDA error
int autostem::checkCudaErr( const char msg[] )
{
    if( cudaGetLastError() != cudaSuccess ) {
        std::string sbuffer = "cuda errror: " + std::string(msg);
            messageAST( sbuffer, 0 );
        exit(0 );  //  should do something better here (?)
    }
    return +1;
}

void autostem::STEMsignals( vectord &x, vectord &y, int npos, vectorf &p,
         int multiMode, double ***detect, int ndetect,
         vectord &ThickSave, int nThick, vectord &sum, vectori &collectorMode,
         vectord &phiMin, vectord &phiMax )
{
    int ix, iy, idetect,  ixmid, iymid;
    int istart, na, ip, i, it;

    long nxl, nyl;

    float scale, prr, pri, tr, ti;

    double  chi0, chi1, k2maxa, k2maxb,
        k2, phi, phirms, alx, aly;
    double sum0, sum1, delta, zslice, totalz;

    vectori ixoff, iyoff;
    vectord xoff, yoff;
     
    /* extra for confocal */
    float hr, hi;
    double chi2C, chi3C, k2maxaC, k2maxbC, r2, rx2;
 
    /* ------ make sure x,y are ok ------ */

    for( ip=0; ip<npos; ip++) {
        if( (x[ip] < 0.0) || (x[ip] > ax) ||
            (y[ip] < 0.0) || (y[ip] > by) ) {
            sum[ip] = -1.2345;
            sbuffer = "bad x,y in STEMsignals() = \n" + toString(x[ip]) + toString(y[ip]);
            messageAST( sbuffer, 0 );
            return;
        }
    }

    /*  generate a probe at position x,y which must be inside the
        0->ax and 0->by region
        normalize to have unit integrated intensity

        probe position (x,y) = (0,0) = lower left corner

    NOTE: The probe wavefunction is nxprobe x nyprobe which may
        be smaller than the transmission function.
        Position the probe near the center of the nxprobe x nyprobe
        region because it does not wrap around properly
        (i.e. manually wrap around in the nx x ny region if
         the probe is near a boundary of the nx x ny region))
    */

    ixmid = nxprobe/2;
    iymid = nyprobe/2;
    chi1 = pi * wavlen;
    k2maxa = apert1 /wavlen;
    k2maxa = k2maxa *k2maxa;
    k2maxb = apert2 /wavlen;
    k2maxb = k2maxb * k2maxb;
    
    /* extra for confocal */
    chi2C = 0.5 * Cs3C *wavlen*wavlen;
    chi3C = Cs5C * wavlen*wavlen*wavlen*wavlen /3.0;
    k2maxaC = apert1C /wavlen;
    k2maxaC = k2maxaC *k2maxaC;
    k2maxbC = apert2C /wavlen;
    k2maxbC = k2maxbC * k2maxbC;

    ixoff.resize( npos );
    iyoff.resize( npos );
    xoff.resize( npos );
    yoff.resize( npos );

    //  calcualate a single unshifted probe
    sum0 = 0.0;
    for( ix=0; ix<nxprobe; ix++) {
        alx = wavlen * kxp[ix];  /* x component of angle alpha */
        for( iy=0; iy<nyprobe; iy++) {
            aly = wavlen * kyp[iy];  /* y component of angle alpha */
            k2 = kxp2[ix] + kyp2[iy];
            if( (k2 >= k2maxa) && (k2 <= k2maxb) ) {
                chi0 = (2.0*pi/wavlen) * chi( p, alx, aly, multiMode );
                chi0= - chi0;
                Hprobe[iy + ix*nyprobe].x = tr = (float) cos( chi0 );  // real
                Hprobe[iy + ix*nyprobe].y = ti = (float) sin( chi0 );  // imag
                sum0 += (double) (tr*tr + ti*ti);
            } else {
                Hprobe[iy + ix*nyprobe].x = 0.0F;
                Hprobe[iy + ix*nyprobe].y = 0.0F;
            }
        }
    }  /* end for( ix... */

    scale = (float) ( 1.0/sqrt(sum0) );  
    for( ix=0; ix<nxprobe; ix++) {
       for( iy=0; iy<nyprobe; iy++) {
            Hprobe[iy + ix*nyprobe].x *= scale;
            Hprobe[iy + ix*nyprobe].y *= scale;
    }  }

    //  copy unshifted probe to device
    cudaMemcpy( Dprobe1, Hprobe, nxprobe*nyprobe*sizeof(cufftComplex),
             cudaMemcpyHostToDevice );
    checkCudaErr( "cannot copy host to device for Dprobe1");

    // ------ CUDA thread index parameters - may need to change for differenet GPUs (?)
    //  for probe as 1D
    int N = nxprobe * nyprobe;
    int threadsPerBlock = 512;  // max thread/block usually 1024 =dev specific
    int blocksPerGrid = (N + threadsPerBlock -1 )/threadsPerBlock;  //  round up

    int thr3 = 16;              //  for probe as 2D
    dim3 threads3(thr3,thr3);
    int blk3x =  (nxprobe + thr3 -1 )/thr3;  //  round up
    int blk3y =  (nyprobe + thr3 -1 )/thr3;  //  round up
    dim3 blocks3(blk3x, blk3y);

    //  for transmission function as 2D
    int blk3xt =  (nx + thr3 -1 )/thr3;  //  round up
    int blk3yt =  (ny + thr3 -1 )/thr3;  //  round up
    dim3 blocks3t(blk3xt, blk3yt);

    //  for CBED sums
    int blocksp =  (nyprobe + thr3 -1 )/thr3;  //  round up

#ifdef AST_USE_CUATOMPOT 
    for( ix=0; ix< natom; ix++) {  // make simple arrays to copy
        Hspec[0 + 4*ix] = xa2[ix]; // pack all 4 into one to min
        Hspec[1 + 4*ix] = ya2[ix]; // number of transfers
        Hspec[2 + 4*ix] = occ2[ix];
        Hspec[3 + 4*ix] = Znum2[ix];
    }
    cudaMemcpy(Dspec, Hspec, 4*natom*sizeof(float), cudaMemcpyHostToDevice );
    checkCudaErr( "cannot copy Dspec from host to device");
#endif
 
    /* ------- calculate all of the probe wave functions at once ------
        to reuse the transmission functions which takes a long
        time to calculate*/
    //  only need one scratch probe in host and nprobe probes on GPU device
    for( ip=0; ip<npos; ip++) {
        ixoff[ip] = (int) floor( x[ip]*((double)nx) / ax ) - ixmid;
        xoff[ip]  = x[ip] - ax*((double)ixoff[ip])/((double)nx);

        iyoff[ip] = (int) floor( y[ip]*((double)ny) / by ) - iymid;
        yoff[ip]  = y[ip] - by*((double)iyoff[ip])/((double)ny);

        //  shift initial probe device - avoid host calc. and big data transfer
        probeShift<<<blocks3, threads3>>>( &Dprobe[ip*nxprobe*nyprobe], Dprobe1,
                nxprobe, nyprobe, xoff[ip], yoff[ip], Dkxp, Dkyp );
        checkCudaErr( "cannot shift probe on device");

    }  /* end for( ip...) */

    /* -------- transmit thru nslice layers ------------------------
        beware ixoff,iyoff must be with one nx,ny
            of the array bounds
    */

    nxl = (long) nx;
    nyl = (long) ny;
    
    scale = 1.0F / ( ((float)nx) * ((float)ny) );

    zslice = 0.75*deltaz;  /*  start a little before top of unit cell */
    istart = 0;
    nslice = 0;
    it = 0;     /* thickness level index */          

    if( zmax > cz ) totalz = zmax;
        else totalz = cz;
    if( 0 != lverbose ) {
        sbuffer= "specimen range is 0 to "+ toString(totalz) + " Ang.";
        messageAST( sbuffer, 0 );
    }
    
    /* range of unit cell */
    while(  (zslice < (totalz+0.25*deltaz)) || (istart<natom) ) {

        // find range of atoms for current slice
        na = 0;
        for(i=istart; i<natom; i++)
            if( za2[i] < zslice ) na++; else break;

        if( 0 != lverbose ) {
           sbuffer= "slice ending at z= "+toString(zslice)+" Ang. with "+toString(na)+" atoms";
           messageAST( sbuffer, 0 );
        }

        // calculate transmission function and bandwidth limit 
        //  trans is used for many probe so doing it on the host is not a huge losss
        //    and its complicated to do on the GPU
        if( na > 0 ) {

#ifdef AST_USE_CUATOMPOT
            //  -- its about the same speed doing the reciprocal space
            //     potential on GPU and real space potential on host
            //  -- so leave both here for future reference

            //  do atomic potnetial on GPU in reciprocal space not real space
            //   to work with fine grain parallel cores

            float mm0 = 1.0F + keV/511.0F;
            float scale2 = wavlen * mm0/(ax*by);

            cuAtompot<<<blocks3t, threads3>>>(  DpotnC, Dspec, na, istart, 
                (float)ax, (float)by, (float)keV, 
                nx, ny, Dkx, Dky, Dkx2, Dky2, (float) k2maxt, Dfparams, scale2 );
            phirms = 0;
            nbeamt = 0;  // ???? need this from somewhere else
            checkCudaErr( "cuAtompot() failed");

            cufftExecC2R( cuplanTc2r, DpotnC, DpotnR );
            checkCudaErr( "cuFFTc2r forward failed");
#else
            //  calculate just real valued sigma*Vz on host - usually faster
            trlayer( xa2, ya2, occ2,
                Znum2, na, istart, (float)ax, (float)by, (float)keV, 
                trans, nxl, nyl, &phirms, &nbeamt, -1.0F );

            nbeamt = 0;  // ???? need this from somewhere else
            for( ix=0; ix<nx; ix++) for( iy=0; iy<ny; iy++)
                HpotnR[ iy + ix*ny ] = trans.re(ix,iy);

            //  copy this real potential to device
            cudaMemcpy( DpotnR, HpotnR, nx*ny*sizeof(cufftReal),
                        cudaMemcpyHostToDevice );
            checkCudaErr( "cannot copy host to device for DpotnR");
#endif
            cuPhasegrating<<<blocks3t, threads3>>>( DpotnR, Dtrans, nx, ny );
            checkCudaErr( "cuPhasegrating() failed");

            cufftExecC2C( cuplanT, Dtrans, Dtrans, CUFFT_FORWARD);
            checkCudaErr( "cuFFT forward failed");

            cuBWlimit<<<blocks3t, threads3>>>( Dtrans,
                Dkx2, Dky2, (float) k2maxp, nx, ny );
            checkCudaErr( "cuBWlimit() failed");

            cufftExecC2C( cuplanT, Dtrans, Dtrans, CUFFT_INVERSE);
            checkCudaErr( "cuFFT inverse failed");

        }

        //----- one multislice trans/prop cycle for all probes ---- 
        for( ip=0; ip<npos; ip++) {

            /* apply transmission function if there are atoms in this slice */
            if( na > 0 ) {

                cufftExecC2C( cuplanP, &Dprobe[ip*nxprobe*nyprobe], 
                    &Dprobe[ip*nxprobe*nyprobe], CUFFT_INVERSE);
                checkCudaErr( "cuFFT inverse failed");

                cmplPixMul<<<blocks3, threads3>>>(  Dtrans, &Dprobe[ip*nxprobe*nyprobe],
                    nx, ny, nxprobe, nyprobe, ixoff[ip], iyoff[ip] );
                checkCudaErr( "cmplPixMul() failed");

                cufftExecC2C( cuplanP, &Dprobe[ip*nxprobe*nyprobe], 
                    &Dprobe[ip*nxprobe*nyprobe], CUFFT_FORWARD);
                checkCudaErr( "cuFFT forward failed");
            } // end if( na>....
    
            //  multiplied by the propagator function 
            cmplVecMul<<<blocksPerGrid, threadsPerBlock>>>( &Dprobe[ip*nxprobe*nyprobe],
                    Dprop, &Dprobe[ip*nxprobe*nyprobe], N );
            checkCudaErr( "cmplVecMul() failed");

        }  // end  for( ip=... 

        /*  if this is a good thickness level then save the ADF or confocal signals
           - remember that the last level may be off by one layer with
              thermal displacements so special case it
        */
    
        //  look at all values because they may not be in order 
        for( it = 0; it<nThick; it++ ) 
        if( fabs(ThickSave[it]-zslice)<fabs(0.5*deltaz)) {
            
            if( 0 != lverbose ) {
                sbuffer= "save ADF/confocal signals, thickness level " + toString(it);  // diagnostic
                messageAST( sbuffer, 0 );
            }

            /*  loop over all probes again */
            for( ip=0; ip<npos; ip++) {

                //  do mag. sq on device 
                magSqPix<<<blocks3, threads3>>>( Dcbed, &Dprobe[ip*nxprobe*nyprobe],
                     nxprobe, nyprobe );
                checkCudaErr( "magSqPix failed");

                for( idetect=0; idetect<ndetect; idetect++) {

                    zeroDbleArray<<<blocksp, thr3>>>( Dsums, nyprobe );
                    checkCudaErr( "zeroDbleArray failed");

                    integCBED<<<blocks3, threads3>>>( Dsums, Dcbed, nxprobe, nyprobe,
                        collectorMode[idetect], Dkxp, Dkyp, Dkxp2, Dkyp2, 
                        k2min[idetect], k2max[idetect], phiMin[idetect], phiMax[idetect] );
                    checkCudaErr( "integCBED failed");

                    cudaMemcpy( Hsums, Dsums, nyprobe*sizeof(double),
                        cudaMemcpyDeviceToHost );
                    checkCudaErr( "cannot copy sums device to host");
                    for( iy=1; iy<nyprobe; iy++) Hsums[0] += Hsums[iy];

                    detect[it][idetect][ip] = Hsums[0];

                }  //  end for( idetect=...

                //  still need to convert this to cuda/gpu but its not used very often
                //      so leave old code for now  30-aug-2018 ejk
                /*  transfer back if confocal needed 
                    - use copy of probe so original can continue in use  */
                if( doConfocal == xTRUE ) {
                    cudaMemcpy( Hprobe, &Dprobe[ip*nxprobe*nyprobe],
                        nxprobe*nyprobe*sizeof(cufftComplex), cudaMemcpyDeviceToHost );
                    checkCudaErr( "cannot copy probe device to host");
                    sum0 = 0;
                    for( ix=0; ix<nxprobe; ix++) {
                        for( iy=0; iy<nyprobe; iy++) {
                            k2 = kxp2[ix] + kyp2[iy];
                            if( (k2 >= k2maxaC) && (k2 <= k2maxbC) ) {
                                phi = atan2( ky[iy], kx[ix] );
                                /*  offset defocus by zslice so both lens referenced to 
                                   entrance surface of specimen  */
                                chi0 = chi1*k2* ( (chi2C + chi3C*k2)*k2 - dfC + zslice
                                    + dfa2C*sin( 2.0*(phi-dfa2phiC) ) 
                                    + 2.0F*dfa3C*wavlen*sqrt(k2)*
                                    sin( 3.0*(phi-dfa3phiC) )/3.0 );
                                chi0= - chi0;
                                hr = (float) cos( chi0 );
                                hi = (float) sin( chi0 );
                                prr = Hprobe[iy + ix*nyprobe].x;   // real
                                pri = Hprobe[iy + ix*nyprobe].y;   // imag
                                probe0.re(ix,iy) = prr*hr -pri*hi;
                                probe0.im(ix,iy) = prr*hi +pri*hr;
                                sum0 += prr*prr + pri*pri;
                            } else {
                                probe0.re(ix,iy) = 0.0F;
                                probe0.im(ix,iy) = 0.0F;
                            }
                        }  /*  end for( iy... )  */
                    }  /*  end for( ix... )  */
            
                    probe0.ifft();
                     
                    /* find normalization constant 
                      i.e. correct for constants in the FFT */
                    sum1 = 0.0;
                    for( ix=0; ix<nxprobe; ix++) {
                         for( iy=0; iy<nyprobe; iy++) {
                            prr = probe0.re(ix,iy);
                            pri = probe0.im(ix,iy);
                            sum1 += prr*prr + pri*pri;
                        }
                   }
            
                   /* integrate over real space detector and normalize */
                   for( ix=0; ix<nxprobe; ix++) {
                           rx2 = xoff[ip] - xp[ix];
                           rx2 = rx2*rx2;
                           for( iy=0; iy<nyprobe; iy++) {
                               r2 = yoff[ip] - yp[iy];
                               r2 = rx2 + r2*r2;
                               prr = probe0.re(ix,iy);
                               pri = probe0.im(ix,iy);
                               delta = prr*prr + pri*pri;
                               for( idetect=0; idetect<ndetect; idetect++) {
                                   if( CONFOCAL == collectorMode[idetect] ) {
                                     if( (r2 >= k2min[idetect] ) &&
                                               (r2 < k2max[idetect] ) )
                                         detect[it][idetect][ip] += delta*(sum0/sum1);
                                   }
                               }
                           }  /* end for( iy... )*/
                   }  /*  end for( ix....) */

                }  /* end if( doConfocal==TRUE) */
            
            }  /* end for( ip.. */

        }  /* end if( ((it...*/

        sbuffer = "xxx bad";
        nslice++;
        zslice += deltaz;
        istart += na;

    }  /* end while( istart...) */

    //  final integrated sum
    for( ip=0; ip<npos; ip++) {

        //  do mag. sq on device 
        magSqPix<<<blocks3, threads3>>>( Dcbed, &Dprobe[ip*nxprobe*nyprobe],
             nxprobe, nyprobe );
        checkCudaErr( "magSqPix failed");

        sum[ip] = 0.0;
        zeroDbleArray<<<nyprobe, 1>>>( Dsums, nyprobe );
        checkCudaErr( "zeroDbleArray failed");

        integCBED<<<blocks3, threads3>>>( Dsums, Dcbed, nxprobe, nyprobe,
                        TOTAL, Dkxp, Dkyp, Dkxp2, Dkyp2, 
                        k2min[0], k2max[0], phiMin[0], phiMax[0] );
        checkCudaErr( "integCBED failed");

        cudaMemcpy( Hsums, Dsums, nyprobe*sizeof(double),
                        cudaMemcpyDeviceToHost );
        checkCudaErr( "cannot copy sums device to host");
        for( iy=0; iy<nyprobe; iy++) sum[ip] += Hsums[iy];
    }

    return;

}/* end autostem::STEMsignals() */
#endif

/*--------------------- trlayer() -----------------------*/
/*  same subroutine in autoslic.c and autostem.c
    (arguments have drifted a little over time)

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
void autostem::trlayer( const vectorf &x, const vectorf &y, const vectorf &occ,
    const vectori &Znum, const int natom, const int istart,
    const float ax, const float by, const float kev,
    cfpix &trans, const long nx, const long ny,
    double *phirms, long *nbeams, const float k2max )
{
    int idx, idy, i, ixo, iyo, ix, iy, ixw, iyw, nx1, nx2, ny1, ny2;
    float k2;
    double r, rx2, rsq, vz, rmin, rmin2, sum, scale, scalex, scaley;

    const double rmax=3.0, rmax2=rmax*rmax; /* max atomic radius in Angstroms */

    scale = sigma( kev ) / 1000.0;  /* in 1/(volt-Angstroms) */

    scalex = ax/nx;
    scaley = by/ny;

    /* min radius to avoid  singularity */
    rmin = ax/((double)nx);
    r = by/((double)ny);
    rmin =  0.25 * sqrt( 0.5*(rmin*rmin + r*r) );
    rmin2 = rmin*rmin;

    idx = (int) ( nx*rmax/ax ) + 1;
    idy = (int) ( ny*rmax/by ) + 1;

    for( ix=0; ix<nx; ix++) {
        for( iy=0; iy<ny; iy++)
           trans.re(ix,iy) = 0.0F;
    }
    
/*  paralleling this loop has little effect   */
/*#pragma omp parallel for private(ix,iy,ixo,iyo,nx1,nx2,ny1,ny2,rx2,r,ixw,iyw,vz,rsq)*/
#pragma omp parallel for private(ix,iy,ixo,iyo,nx1,nx2,ny1,ny2,rx2,ixw,iyw,vz,rsq)
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
                  if( rsq < rmin2 ) rsq = rmin2;
                  /*r = sqrt( rx2 + r*r );
                  vz = occ[i] * vzatom( Znum[i], r ); slow */
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

 };  /* end autostem::trlayer() */

/*------------------------- invert2D() ----------------------*/
/*
        rearrange pix with corners moved to center (a la FFT's)

         pix[ix][iy] = real array with image
         nx,ny = range of image 0<ix<(nx-1) and 0<iy<(ny-1)

*/
void autostem::invert2D( float** pix, long nx, long ny )
{
#define SWAP(a,b)       {t=a; a=b; b=t;}

        long ix, iy, ixmid, iymid;
        float t;

        ixmid = nx/2;
        iymid = ny/2;

        for( ix=0; ix<nx; ix++) 
        for( iy=0; iy<iymid; iy++)
                SWAP( pix[ix][iy], pix[ix][iy+iymid] );

        for( ix=0; ix<ixmid; ix++) 
        for( iy=0; iy<ny; iy++)
                SWAP( pix[ix][iy], pix[ix+ixmid][iy] );

#undef SWAP
};  // end autostem::invert2D()
