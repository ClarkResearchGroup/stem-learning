/*      *** autostemcmd.cpp ***

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
   ANSI-C and TIFF version
  command line front end for autostem.cpp

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
  remove nprobes as a separate variable (will be nyout) 15-aug-2013 ejk
  fix minor format issues 5-oct-2013 ejk
  save output size in param[] 6-oct-2013 ejk
  add param[] for mode 14-oct-2013 ejk
  convert detector angles to radians 19-oct-2013 ejk
  convert to streams and strings 28-apr-2014 ejk
  fix small bug in info .txt files 4-may-2014 ejk
  change nbeamt to long to be consistent with autostem.cpp
     and mode to mAUTOSTEM 25-oct-2015 ejk
   convert malloc1D() etc to vector<> 5-jul-2016 ejk
  convert malloc2D(), malloc3D() to new2D(), new3D() 24-sep-2017 ejk
  add abbError() 28-feb-2018, 1-mar-2018  ejk
  start segmented detector 28-apr-2018 ejk
  add pPMAXDET,pPMINDET params for azim. STEM detect. angles (seg. detector)
       12-may-2018 ejk
  remove sourceFWHM question because it is not practical 20-jan-2019 ejk
  force multiMode=1 with abbError() (which makes more sense)
      and add more info output in 1D mode 20-apr-2019 ejk
  fix bug with init Rng in abbError() mode 21-apr-2019 ejk
  NOTE: confocal values changed by small amount after fixing small error
     in xpos calculation in slicelib/freqn()  21-aug-2019 ejk
  
*/

#include <cstdio>  /* ANSI C libraries used */
#include <cstdlib>
#include <cstring>
#include <cmath> 
#include <ctime>

#include <string>
#include <iostream>  //  C++ stream IO
#include <fstream>
#include <iomanip>   //  to format the output
#include <vector>

using namespace std;

#include "cfpix.hpp"      // complex image handler with FFT 
#include "slicelib.hpp"   // misc. routines for multislice
#include "floatTIFF.hpp"  // file I/O routines in TIFF format

//  use only one;  the calculation part
//#include "autostem_cuda.hpp"   // header for cuda nvcc version of this class
#include "autostem.hpp"        // header for C++ version of this class

#define MANY_ABERR      //  define to include many aberrations

//
#define USE_OPENMP      /* define to use openMP */

#ifdef USE_OPENMP
#include <omp.h>
/*  get wall time for benchmarking openMP */
#define walltim() ( omp_get_wtime() )
double walltimer;
#endif

enum{ TRUE=1, FALSE=0};

int main( int argc, char *argv[ ] )
{  
    string filein, fileout, fileoutpre, description, cmode, cline,
        pacbedFile;
    const string version = "21-aug-2019 (ejk)";

    int ix, iy, i, idetect, nxout, nyout,
        ncellx, ncelly, ncellz, nwobble, ndetect, ip, nThick, it,
        done, status, multiMode;
    int nx, ny, nxprobe, nyprobe, nslice, natom;

    int l1d=0, lwobble=0, lxzimage=0, labErr=0, NPARAM, np, echo;
    int lpacbed, ixo,iyo, ix2,iy2, nx1,nx2, ny1,ny2, nxout2,nyout2;
    int doConfocal, doSegment;  // for confocal, segmented detector mode

    int nbeamp, nbeampo;
    long  ltime, nbeamt;
    unsigned long iseed;

    float ***pixr, temp, pmin, pmax, aimin, aimax;
    float wmin, wmax, xmin,xmax, ymin, ymax, temperature;
    float rmin0, rmax0, scalef, res, thetamax;
    float ax, by, cz;                         //  specimen dimensions
    float dfa2C, dfa2phiC, dfa3C, dfa3phiC;   // astigmatism parameters
    float  **rmin, **rmax;

    float zmin, zmax;
    float **pacbedPix;        //  to save position averaged CBED 

    double dxp, dyp;

    double xi,xf, yi,yf, dx, dy, totmin, totmax,
       ctiltx, ctilty, timer, sourceFWHM,  vz;

    double wavlen, Cs3,Cs5, df,apert1, apert2, pi, keV;
    double deltaz;
    double Cs3C, Cs5C, dfC, apert1C, apert2C; // aberrations of collector lens

    vectori collectorMode;  // for ADF, confocal or segmented
    vectori Znum;
    vectorf param;
    vectorf xa, ya, za, occ, wobble;

    //  almin/max = detector polar angles, thetaMin/Max= detector azimuthal angles
    vectord ThickSave, almin, almax, x, y, phiMin, phiMax;

    ofstream fp;
    
    floatTIFF myFile;

    autostem ast;    // "stack smashing detected" error is from here !! ????

/* start by announcing version etc */

    cout << "autostem (ADF,confocal,segmented) version dated " << version  << endl;
    cout << "Copyright (C) 1998-2019 Earl J. Kirkland"  << endl;
    cout <<  "This program is provided AS-IS with ABSOLUTELY NO WARRANTY\n "
            << " under the GNU general public license\n"  << endl;

    cout <<  "Calculate STEM images using FFTW" << endl;
#ifdef USE_OPENMP
    cout <<  "and multithreaded using openMP" << endl;
#endif
    cout << "\n";

/*    get option to save position averaged CBED pattern */
    if( argc >= 2 ) {
        pacbedFile =  argv[1]; 
        cout << "calculate 2D position averaged CBED pattern (arb. units) in file "
                << pacbedFile << endl;
        lpacbed = TRUE;
    } else {
        cout << "To calculate 2D pos. aver. CBED, include a file name on command line." << endl;
        lpacbed = FALSE;
    }

/*  get simulation options */

    pi = 4.0 * atan( 1.0 );
    NPARAM = myFile.maxParam();
    param.resize( NPARAM );
    for( i=0; i<NPARAM; i++) param[i] = 0.0F;

    cout << "Name of file with input atomic "
          << "potential in x,y,z format:" << endl;
    cin >> filein;

    cout << "Replicate unit cell by NCELLX,NCELLY,NCELLZ :" << endl;
    cin >> ncellx >> ncelly >> ncellz;
    if( ncellx < 1 ) ncellx = 1;
    if( ncelly < 1 ) ncelly = 1;
    if( ncellz < 1 ) ncellz = 1;

/*  get more parameter etc */

    cout << "STEM probe parameters, V0(kv), Cs3(mm), Cs5(mm),"
           << " df(Angstroms), apert1,2(mrad) :" << endl;
    cin >> keV >> Cs3 >> Cs5 >> df >> apert1 >> apert2;
    param[pDEFOCUS] = (float) df;
    param[pCS]  = (float) ( Cs3*1.0e7 );
    param[pCS5] = (float) ( Cs5*1.0e7 );
    apert1 *= 0.001;   // convert to radians
    apert2 *= 0.001;

#ifdef MANY_ABERR
    /*   get higher order aberrations if necessary */
    cout << "type higher order aber. name (as C32a, etc.) followed\n"
        << " by a value in mm. (END to end)" << endl;
    done = multiMode = 0;
    do{
        cin >> cline;
        if( cline.compare( "END" ) == 0 ) {
            done = 1;
        } else {
            cin >> vz;
            status = readCnm( cline, param, vz );        
            if( status < 0 ) {
                cout << "unrecognized aberration, exit..." << endl;
                exit( EXIT_SUCCESS );
            } else multiMode = 1;
        }
    } while( !done );

#endif

    labErr = askYN("Do you want to add random pi/4 tunning errors for orders 2 to 5");

    wavlen = wavelength( keV );
    cout << "wavelength = " << wavlen << " Angstroms" << endl;
    if( apert1 > apert2 ) {
       cout << "Bad probe aperture specification." << endl;
       cout << "apert1 must be less than apert2." << endl;
       cout << "apert1= " << apert1 << ", apert2 = " << apert2 << endl;
       exit( 0 );
    }

    /*  FFTW is not restricted to 2^n sizes so don't need test anymore */
    cout << "Size of specimen transmission function"
          << " Nx,Ny in pixels : \n" << endl;
    cin >> nx >> ny;

    cout << "Size of probe wave function"
          << " Nx,Ny in pixels : " << endl;
    cin >> nxprobe >> nyprobe;

    cout << "Crystal tilt x,y in mrad. :" << endl;
    cin >> ctiltx >> ctilty;
    ctiltx = ctiltx * 0.001;
    ctilty = ctilty * 0.001;

    l1d = askYN("Do you want to calculate a 1D line scan");

    if( l1d == 1 ) {
        lxzimage = askYN("Do you want to save all depth information as xz image");
        nThick = 1;
    } else {
        do { cout << "Number of thickness levels to save, including"
                " the end(>=1):" << endl;
             cin >> nThick;
        } while (nThick <= 0);
        ThickSave.resize( nThick );
        if( nThick > 1 ) {
           cout << "type thickness (in Ang.) of " << (nThick-1)
               << " intermediate layers :" << endl; 
           for( it=0; it<(nThick-1); it++) cin >> ThickSave[it];
        }
    }

    cout << "File name prefix to get output of STEM multislice result "
        << "(no extension):" << endl;
    cin >> fileoutpre;

    do {    cout << "Number of detector geometries (>=1):" << endl;
            cin >> ndetect;
    } while (ndetect <= 0);

    almin.resize( ndetect );
    almax.resize( ndetect );
    collectorMode.resize( ndetect );
    phiMin.resize( ndetect );
    phiMax.resize( ndetect );

    doConfocal = FALSE;  // flag confocal mode
    doSegment = FALSE;  // flag segmented (i.e. quadrant) detector mode

    for( idetect=0; idetect<ndetect; idetect++) {
        cout << "Detector " << idetect+1 <<", type: min max polar angles(mrad)"
               << " or radius(Ang.) " << endl;
        cout << "followed by m, A, or seg (ADF, confocal, segmented-mrad)" 
               << endl;
        cout << "add phimin phimax (degrees) if segmented" << endl;
        cin >> almin[idetect] >> almax[idetect] >> cmode;
        if( (cmode.compare("m") == 0) || (cmode.compare("M")==0 ) ) {
                collectorMode[idetect] = ADF;
                cout << "normal ADF detector" << endl;
                almin[idetect] *= 0.001;  // convert to radians
                almax[idetect] *= 0.001;
                if( almin[idetect] > almax[idetect] ) {
                    cout << "bad angles, exit..." << endl;
                    return( 0 );
                }
        } else if( (cmode.compare("a")==0) || (cmode.compare("A")==0) ) {
                collectorMode[idetect] = CONFOCAL;
                cout << "confocal detector" << endl;
                doConfocal = TRUE;
        } else if( (cmode.compare("seg")==0) || (cmode.compare("SEG")==0) ) {
                collectorMode[idetect] = ADF_SEG;
                cin >> phiMin[idetect] >> phiMax[idetect];
                cout << "ADF segment detector phi= " << phiMin[idetect] 
                    << " to " << phiMax[idetect] << endl;
                //  remember atan2() goes from -pi(-180) to +pi(+180)
                if( (phiMin[idetect] < -180 ) || (phiMin[idetect] > 180 )
                    || (phiMax[idetect] < -180 ) || (phiMax[idetect] > 180 ) ){
                        cout << "phi must be between -180 and +180 deg., exit..." << endl;
                        return( 0 );
                }
                phiMin[idetect] *= (pi/180.0);
                phiMax[idetect] *= (pi/180.0);
                almin[idetect] *= 0.001;  // convert to radians
                almax[idetect] *= 0.001;
                doSegment = TRUE;
                if( (phiMin[idetect] > phiMax[idetect]) ||
                    (almin[idetect] > almax[idetect]) ) {
                    cout << "bad angles, exit..." << endl;
                    return( 0 );
                }
        } else {
                cout << "unrecognized collector mode = " << cmode << endl;
                exit( 0 );
        }
    }  /* end for(idetect=.. */

    Cs3C = Cs5C = dfC = apert1C = apert2C = 0;
    dfa2C = dfa2phiC = dfa3C = dfa3phiC = 0;
    if( doConfocal == TRUE ) {
        cout << "Collector lens parameters, Cs3(mm), Cs5(mm),"
                << " df(Angstroms), apert1,2(mrad) :" << endl;
        cin >> Cs3C >> Cs5C >> dfC >> apert1C >> apert2C;
        cout << "Magnitude and angle of 2-fold astigmatism"
                      <<  " (in Ang. and degrees):" << endl;
        cin >> dfa2C >> dfa2phiC;
        dfa2phiC = (float) (dfa2phiC * pi /180.0F);
        cout << "Magnitude and angle of 3-fold astigmatism"
              <<  " (in Ang. and degrees):" << endl;
        cin >> dfa3C >> dfa3phiC;
        dfa3phiC = (float) (dfa3phiC * pi /180.0F);

        if( apert1C > apert2C ) {
                cout << "Bad collector aperture specification." << endl;
                cout << "apert1 must be less than apert2." << endl;
                cout << "apert1= " << apert1C << ", apert2 = " << apert2C << endl;
                exit( 0 );
        }
    }  /* end if( doConfocal...  */

    if( l1d == 1 ) {
        cout << "xi, xf, yi, yf, nout :" << endl;
        cin >> xi >> xf >> yi >> yf >> nyout;
        nxout = 1;
    }else {
        cout << "xi,xf,yi,yf, nxout,nyout :" << endl;
        cin >> xi >> xf >> yi >> yf >> nxout >> nyout;
    }

    /*  remember that the slice thickness must be > atom size
        to use projected atomic potential */
    cout << "Slice thickness (in Angstroms):" << endl;
    cin >> deltaz;
    if( deltaz < 1.0 ) {
        cout << "WARNING: this slice thickness is probably too thin"
            << " for autostem to work properly." << endl;
    }

    lwobble = askYN("Do you want to include thermal vibrations");
    if( lwobble == 1 ) {
        cout << "Type the temperature in degrees K:" << endl;
        cin >> temperature;
        /* get random number seed from time if available 
            otherwise ask for a seed */
        cout << "Type number of configurations to average over:" << endl;
        cin >> nwobble;
        if( nwobble < 1 ) nwobble = 1;
        //  source size doesn't work here - requires too many MC samples
        //  but leave part here so I don't forget 30-jan-2019 ejk
        //cout << "Type source size (FWHM in Ang.):" << endl;
        //cin >> sourceFWHM;
        sourceFWHM = 0.0;
    } else {
        temperature = 0.0F;
        nwobble = 1;
        sourceFWHM = 0.0;
        iseed = 0;
    }

    if( (TRUE ==lwobble ) || (TRUE  == labErr) ) {  // for random numbers
        ltime = (long) time( NULL );
        iseed = (unsigned) ltime;
        if( ltime == -1 ) {
            cout << "Type initial seed for random number generator:" << endl;
            cin >> iseed;
        } else {
            cout << "Random number seed initialized to " << iseed << endl;
        }
    }

    timer = cputim();   /* get initial CPU time */
#ifdef USE_OPENMP
    walltimer = walltim();  /* wall time for opneMP */
#endif

    /*  calculate relativistic factor and electron wavelength */
    wavlen = (float) wavelength( keV );
    cout << "electron wavelength = " << wavlen << " Angstroms" << endl;

    /*---- read in specimen coordinates and scattering factors----- */

    natom = ReadXYZcoord( filein.c_str(), ncellx, ncelly, ncellz,
        &ax, &by, &cz, Znum, xa, ya, za, occ, wobble,
        description );

    cout << natom <<" atomic coordinates read in" << endl;
    cout << description << endl;

    cout << "Lattice constant a,b,c = " << ax << ", " << by << ", " << cz << endl;

    /* calculate thickness levels to save (1D mode) or check range (2D mode) */
    if( l1d == 1 ) {
        if( lxzimage == 1 ) {
            /* save all thickness levels  */
            nThick = (int) ( cz/deltaz + 0.5 );
            ThickSave.resize( nThick );
            for( it=0; it<nThick; it++) {
                ThickSave[it] = deltaz*(it+1);
            }
        } else {
            nThick = 1;
            ThickSave.resize( nThick );
            ThickSave[0] = cz;
        }
        cout << "save up to " << nThick << " thickness levels" << endl;  // diagnostic
    } else {
        ThickSave[nThick-1] = cz;  /*  always save the last level */
        for( it=0; it<(nThick-1); it++) 
        if( (ThickSave[it] < 0.0) || (ThickSave[it] > cz) ) {
            cout << "Bad thickness level = "<< ThickSave[it] << " A,"
                << "allowed range= 0.0 to " << cz << " A" << endl;
            exit( 0 );
        }
    }  /* end if( l1d == ... */

    /*  calculate the total specimen volume and echo */
    xmin = xmax = xa[0];
    ymin = ymax = ya[0];
    zmin = zmax = za[0];
    wmin = wmax = wobble[0];

    for( i=0; i<natom; i++) {
        if( xa[i] < xmin ) xmin = xa[i];
        if( xa[i] > xmax ) xmax = xa[i];
        if( ya[i] < ymin ) ymin = ya[i];
        if( ya[i] > ymax ) ymax = ya[i];
        if( za[i] < zmin ) zmin = za[i];
        if( za[i] > zmax ) zmax = za[i];
        if( wobble[i] < wmin ) wmin = wobble[i];
        if( wobble[i] > wmax ) wmax = wobble[i];
    }
    cout << "Total specimen range is\n " << xmin << " to " << xmax << " in x\n"
           << ymin << " to " << ymax << " in y\n"
           << zmin << " to " << zmax << " in z" << endl;
    if( lwobble == 1 )
        cout << "Range of thermal rms displacements (300K) = " 
             << wmin << " to " << wmax << endl;

    /*  check for valid scan coordinates  */

    if( (xi < 0.0) || (xi > ax) ||
        (xf < 0.0) || (xf > ax) ||
        (yi < 0.0) || (yi > by) ||
        (yf < 0.0) || (yf > by) ) {
            cout << "WARNING: Coordinates out of range; will be made periodic."<< endl;
            cout << "xi,xf,yi,yf= "<< xi<< ", "<< xf<< ", "<< yi << ", "<< yf << endl;
    }

    /*  check that requested probe size is not bigger 
        than transmission function size (or too small)
    */
    if( (nxprobe > nx) || (nxprobe < 2) ) {
        nxprobe = nx;
        cout << "Probe size reset to nx = " << nxprobe << endl;
    }

    if( (nyprobe > ny) || (nyprobe < 2) ) {
        nyprobe = ny;
        cout << "probe size reset to ny = " << nyprobe  << endl;
    }

    param[ pAX ] = ax;
    param[ pBY ] = by;
    param[ pCZ ] = cz;
    param[ pNX ] = (float) nx;      // size of transmission function (in pixels)
    param[ pNY ] = (float) ny;
    param[ pNXPRB ] = (float) nxprobe;    // probe size in pixels
    param[ pNYPRB ] = (float) nyprobe;
    param[ pENERGY ] = (float) keV;
    param[ pWAVEL ]   = (float) wavlen;
    param[ pOAPMIN ] = (float) apert1;      // objective aperture in radians
    param[ pOAPERT ] = (float) apert2;

    param[ pCDF ] = (float) dfC;      // collector defocus
    param[ pCDFA2 ] = (float) dfa2C;      // collector astig 2nd order
    param[ pCDFA2PHI ] = (float) dfa2phiC;
    param[ pCDFA3 ] = (float) dfa3C;      // collector astig 2nd order
    param[ pCDFA3PHI ] = dfa3phiC;
    param[ pCCS3 ] = (float) ( Cs3C * 1.0e7 );    // collector spherical aberr.
    param[ pCCS5 ] =  (float) ( Cs5C * 1.0e7 );
    param[ pCCAPMIN ] = (float) ( apert1C * 0.001 );  //  collector apert. in radians
    param[ pCCAPMAX ] = (float) ( apert2C * 0.001 );

    np = 22;  // number of abberations to add tuning errors (22 for 2-5th order)
    echo = 1;  // print some information
    if( TRUE == labErr ) {
        ast.abbError( param, np, NPARAM, iseed, echo );
        multiMode = 1;
        cout << "add random pi/4 aberration tuning errors" << endl;
    }

    ast.CountBeams( param, nbeamp, nbeampo, res, thetamax );

    //  transmission function sampling
    cout << "Bandwidth limited to a real space resolution of "
        << res << " Angstroms" << endl;
    cout << "   (= " << thetamax*1000.0F << 
        " mrad)  for symmetrical anti-aliasing." << endl;

    //  probe sampling
    cout << "Number of symmetrical anti-aliasing "
          << "beams in probe = " << nbeamp << endl;
    cout << "Number of beams in probe aperture = " << nbeampo << endl;
    if( nbeamp < 200 ) {
        cout << "WARNING: the probe is under sampled, this is a bad idea..." << endl;
    }
    if( nbeampo < 100 ) {
        cout << "WARNING: the probe aperture is under sampled, this is a bad idea..." << endl;
        exit( 0 );
    }

    /* init the min/max record of total integrated intensity */
    rmin = new2D<float>( nThick, ndetect, "rmin" );
    rmax = new2D<float>( nThick, ndetect, "rmax" );

    if( lpacbed == TRUE ) {
        if(  l1d == 1 ) {
           cout << "Cannot do pos. aver. CBED in 1d, exit...." << endl;
           exit( 0 );
        }
        pacbedPix = new2D<float>( nxprobe, nyprobe, "pacbedPix" );
        for( ix=0; ix<nxprobe; ix++) for( iy=0; iy<nyprobe; iy++)
                pacbedPix[ix][iy] = 0;
    } else pacbedPix = NULL;

    //   set autostem modes
    ast.l1d = l1d;
    ast.lpacbed = lpacbed;
    ast.lwobble = lwobble;
    ast.lxzimage = lxzimage;
    //????? ast.lverbose = 1;
    ast.lverbose = 0;
   
    // ---- setup parameters in param[] - some already set above
    param[ pAX ] = ax;
    param[ pBY ] = by;
    param[ pNX ] = (float) nx;    // size of transmission function (in pixels)
    param[ pNY ] = (float) ny;
    param[ pXCTILT ] = (float) ctiltx;      // crystal tilt
    param[ pYCTILT ] = (float) ctilty;
    param[ pOAPERT ] = (float) apert2;      // objective aperture
    param[ pTEMPER ] = temperature;     // temperature
    param[ pNWOBBLE ] = (float) nwobble;    //  number config. to average
    param[ pDELTAZ ] = (float) deltaz;      // slice thickness
    param[ pSOURCE ] = (float) sourceFWHM;  // source size  -- should get rid of this ??? doesn't work

    param[ pNXPRB ] = (float) nxprobe;      // probe size in pixels
    param[ pNYPRB ] = (float) nyprobe;

    cout << "output file size in pixels is "
          << nxout << " x " << nyout  << endl;
    // double up first index to mimic a 4D array
    pixr = new3D<float>( ndetect*nThick, nxout, nyout, "pixr" );

   //  do the autostem calculation
   ast.calculate( param, multiMode, natom, &iseed, 
        Znum, xa,ya,za, occ, wobble,
        xi,xf,yi,yf, nxout, nyout,
        ThickSave, nThick,
        almin, almax, collectorMode, ndetect,
        phiMin, phiMax,
        pixr, rmin, rmax, pacbedPix );

    nslice = (int) ((zmax-zmin)/deltaz + 0.5);   // may be off by 1 or 2 with wobble
    nbeamt = ast.nbeamt;   //  ??? get beam count - should do this better
    totmin = ast.totmin;
    totmax = ast.totmax;

    param[pMODE] = mAUTOSTEM;  // save mode = autostem

    // ------------- start here for a full image output --------------
    if( l1d == 0 ) {

       dx = (xf-xi)/((double)(nxout-1));  // pixels size for image output
       dy = (yf-yi)/((double)(nyout-1));

        /*  directory file listing parameters for each image file */
        fileout = fileoutpre + ".txt";
        fp.open( fileout.c_str() );
        if( fp.bad() ) {
            cout << "Cannot open output file " << fileout << endl;
            exit( 0 );
        }
     
       fp << "C" << endl;
       fp << "C   output of autostem version " << version << endl;
       fp << "C" << endl;
       fp << "C   nslice= " << nslice << endl;
       fp << "C deltaz= " << deltaz << ", file in= " << filein  << endl;
       fp << "C V0= "<< keV <<", Cs3= "<< Cs3 <<", Cs5= "<< Cs5 <<", df= "<< df << endl;
       fp << "C Apert= " << apert1*1000.0 << " mrad to " << apert2*1000.0 << " mrad" << endl;
       if( doConfocal == TRUE ) {
          fp << "C confocal lens Cs3= "<< Cs3C << ", Cs5= " << Cs5C << ", df= " << dfC << endl;
          fp << "C confocal lens apert= " << apert1C << " mrad to " << apert2C << " mrad" << endl;
          fp << "C confocal dfa2C= "<< dfa2C << ", dfa2phiC= " << dfa2phiC  
              << ", dfa3C= " << dfa3C << ", dfa3phiC= " << dfa3phiC  << endl;
       }
       fp << "C Crystal tilt x,y= " << ctiltx << ", " << ctilty << endl;

       for(  idetect=0; idetect<ndetect; idetect++) {
            if( ADF == collectorMode[idetect]  )
                fp << "C Detector " << idetect << ", Almin= " << almin[idetect]*1000.0
                    << " mrad, Almax= " << almax[idetect]*1000.0  << " mrad" << endl;
            else if( CONFOCAL == collectorMode[idetect] )
                fp << "C Detector " << idetect << ", cmin= " << almin[idetect]
                    << " Angst, cmax= " << almax[idetect] << " Angst." << endl;
            else if( ADF_SEG == collectorMode[idetect] )
                fp << "C Detector " << idetect << ", Almin= " << almin[idetect]*1000.0
                    << " mrad, Almax= " << almax[idetect]*1000.0  << " mrad, " 
                    << "phimin= " << phiMin[idetect]*180.0/pi << ", phimax= " 
                    << phiMax[idetect]*180.0/pi  << " deg." << endl;
       }

       fp << "C ax= " << ax << " A, by= " << by << " A, cz= " << cz << " A" << endl;
       fp << "C Number of symmetrical anti-aliasing "
            << "beams in probe wave function= " << nbeamp  << endl;
       fp << "C with a resolution (in Angstroms) = " << res  << endl;
       if( lwobble == 1 ) {
          fp << "C Number of thermal configurations = " << nwobble << endl;
          fp << "C Source size = " << sourceFWHM << " Ang. (FWHM)" << endl;
       }
       if( TRUE == labErr ) {
          fp << "C  add pi/4 aberr. tuning errors" << endl;
	   }
       fp << endl;

        /*  store params plus min and max */
        param[pIMAX]    = 0.0F;
        param[pIMIN]    = 0.0F;
        param[pXCTILT]  = (float) ctiltx;
        param[pYCTILT]  = (float) ctilty;
        param[pDEFOCUS] = (float) df;
        param[pDX]      = (float) dx;
        param[pDY]      = (float) dy;
        param[pENERGY]  = (float) keV;
        param[pOAPERT]  = (float) apert2;
        param[pWAVEL]   = (float) wavlen;
        param[pNSLICES] = (float) nslice;
        param[ pNXOUT ] = (float) nxout;  // size of output (in pixels)
        param[ pNYOUT ] = (float) nyout;

        for( i=0; i<NPARAM; i++) myFile.setParam( i, param[i] );
        myFile.setnpix( 1 );
        myFile.resize( nxout, nyout );
        aimin = aimax = 0.0F;

        for( it=0; it<nThick; it++)
        for( i=0; i<ndetect; i++) {
            //sprintf( fileout, "%s%03d%03d.tif", fileoutpre, i, it );   // plain C
            fileout = fileoutpre + toString(i) + "_" + toString(it) + ".tif";
            cout << fileout << ": output pix range : " << rmin[it][i] << " to " 
                << rmax[it][i] << endl;
            myFile.setParam( pRMAX, rmax[it][i] );
            myFile.setParam( pRMIN, rmin[it][i] );
            myFile.setParam(pMINDET, 0 );
            myFile.setParam(pMAXDET, 0 );
            myFile.setParam(pPMINDET, 0 );
            myFile.setParam(pPMAXDET, 0 );
            myFile.setParam(pCMINDET, 0 );
            myFile.setParam(pCMAXDET, 0 );
            if( ADF == collectorMode[i]  ) {
                myFile.setParam( pMINDET, (float) ( almin[i] ) );
                myFile.setParam( pMAXDET, (float) ( almax[i] ) );
            } else if( CONFOCAL  == collectorMode[i] ) {
                myFile.setParam(pCMINDET, (float) almin[i] );
                myFile.setParam(pCMAXDET, (float) almax[i] );
            } else if( ADF_SEG == collectorMode[i] ) {
                myFile.setParam( pMINDET, (float) ( almin[i] ) );
                myFile.setParam( pMAXDET, (float) ( almax[i] ) );
                myFile.setParam( pPMINDET, (float) ( phiMin[i] ) );
                myFile.setParam( pPMAXDET, (float) ( phiMax[i] ) );
            }
            for( ix=0; ix<nxout; ix++) for( iy=0; iy<nyout; iy++)  
                myFile( ix, iy ) = pixr[i+it*ndetect][ix][iy];
            if( myFile.write( fileout.c_str(), rmin[it][i], rmax[it][i], aimin, aimax,
                (float) dx, (float) dy ) != 1 ) {
                    cout << "Cannot write output file " << fileout << endl;
            }

            if( ADF == collectorMode[i]  )
                fp << "file: " << fileout 
                    << ", detector= " << almin[i]*1000.0 << " to " << almax[i]*1000.0 << " mrad, "
                    << "thicknes= " << ThickSave[it] << " A, range= " << rmin[it][i] 
                    << " to " << rmax[it][i] << endl;
            else if( CONFOCAL == collectorMode[i] )
                fp <<"file: " << fileout
                    << ", detector= " << almin[i] << " to " << almax[i] << " Angst., "
                    << "thicknes= " << ThickSave[it] << " A, range= " << rmin[it][i]
                    << " to " << rmax[it][i] << endl;
            else if( ADF_SEG == collectorMode[i]  )
                fp << "file: " << fileout 
                    << ", detector= " << almin[i]*1000.0 << " to " << almax[i]*1000.0 << " mrad, "
                    << "thicknes= " << ThickSave[it] << " A, range= " << rmin[it][i] 
                    << " to " << rmax[it][i] << endl;
        }  /*  end for(i=... */

        fp.close();

        /*   save pos. aver. CBED if needed */
        if( lpacbed == TRUE ) {
            myFile.setnpix( 1 );
            nx1 =    nxprobe / 6;
            nx2 = (5*nxprobe) / 6;  /*  cut out center portion without anti-aliasing zeros */
            ny1 =    nyprobe / 6;
            ny2 = (5*nyprobe) / 6;
            nxout2 = nx2 - nx1 + 1;
            nyout2 = ny2 - ny1 + 1;
            if( (nxout2<1) || (nyout2<1) ) {
                    nx1 = ny1 = 0;
                    nx2 = nxout2 = nxprobe;
                    ny2 = nyout2 = nyprobe;
            }
            myFile.resize( nxout2, nyout2 );
            aimin = aimax = 0.0F;
            scalef = (float) (1.0/(nxout2*nyout2) );
            ixo = 0;
            for( ix2=nx1; ix2<=nx2; ix2++) {
                iyo = 0;
                for( iy2=ny1; iy2<=ny2; iy2++) {
                    myFile(ixo,iyo++) = scalef * pacbedPix[ix2][iy2];
                }  ixo++;
            }
            rmin0 = myFile.min(0);
            rmax0 = myFile.max(0);
            dxp = ax*((double)nxprobe)/nx;
            dyp = by*((double)nyprobe)/ny;
            dxp = 1.0/dxp;
            dyp = 1.0/dyp;
            myFile.setParam( pRMAX, rmax0 );
            myFile.setParam( pRMIN, rmin0 );
            myFile.setParam( pDX, (float) dxp);
            myFile.setParam( pDY, (float) dyp );
            cout << "pos. averg. CBED (unaliased) size " << nxout2 << " x " << nyout2 << " pixels\n"
                << " and range (arb. units): " << rmin0 << " to " << rmax0 << endl;
            if( myFile.write( pacbedFile.c_str(), rmin0, rmax0, aimin, aimax,
                (float) dxp, (float) dyp ) != 1 ) {
                cout << "Cannot write output file " << pacbedFile << endl;
            }
        }   /*  end if( lpacbed.... */

    /* ------------- start here for 1d line scan output ---------------- */

    } else if ( l1d == 1 ) {

       dx = (xf-xi)/((double)(nyout-1));
       dy = (yf-yi)/((double)(nyout-1));
       x.resize( nyout );
       y.resize( nyout );
       for( ip=0; ip<nyout; ip++) {
           /* recalculate mean x,y without source size wobble */
           x[ip] = xi + dx * ((double)ip);
           y[ip] = yi + dy * ((double)ip);
       }

       /* ------ first output text data ---------------- */
       fileout = fileoutpre + ".dat";
       cout << "output file= " << fileout << endl;

       fp.open( fileout.c_str() );
       if( fp.bad() ) {
           cout << "Cannot open output file " << fileout << endl;
             exit( 0 );
       }

       fp << "C" << endl;
       fp << "C   output of autostem version " << version << endl;
       fp << "C" << endl;
       fp << "C   nslice= " << nslice << endl;
       fp << "C deltaz= " << deltaz << ", file in= " << filein  << endl;
       fp << "C V0= "<< keV << ", Cs3= "<< Cs3<< ", Cs5= "<< Cs5<< ", df= "<< df  << endl;
       fp << "C Apert= " << apert1*1000.0 << " mrad to " << apert2*1000.0 << " mrad" << endl;
       if( doConfocal == TRUE ) {
          fp << "C confocal lens Cs3= "<< Cs3C<< ", Cs5= "<< Cs5C << ", df= "<< dfC << endl;
          fp << "C confocal lens apert= " << apert1C << " mrad to " <<apert2C << " mrad" << endl;
          fp << "C confocal dfa2C= " << dfa2C << ", dfa2phiC= " << dfa2phiC << ", dfa3C= "
              << dfa3C << ", dfa3phiC= " << dfa3phiC << endl;
       }
       fp << "C Crystal tilt x,y= " << ctiltx << ", " << ctilty << endl;

       for(  idetect=0; idetect<ndetect; idetect++) {
            if( ADF == collectorMode[idetect]  )
                fp << "C Detector " << idetect << ", Almin= " << almin[idetect]*1000.0 
                   << " mrad, Almax= " << almax[idetect]*1000.0 << " mrad" << endl;
            else if( CONFOCAL == collectorMode[idetect] )
                fp << "C Detector " <<idetect << ", cmin= " << almin[idetect] 
                    << " Angst, cmax= " << almax[idetect] << " Angst." << endl;
            else if( ADF_SEG == collectorMode[idetect] )
                fp << "C Detector " << idetect << ", Almin= " << almin[idetect]*1000.0
                    << " mrad, Almax= " << almax[idetect]*1000.0  << " mrad, " 
                    << "phimin= " << phiMin[idetect]*180.0/pi << ", phimax= " 
                    << phiMax[idetect]*180.0/pi  << " deg." << endl;
       }

       fp << "C ax= " << ax << " A, by= " << by << " A, cz= " << cz << " A" << endl;
       fp << "C Number of symmetrical anti-aliasing "
               << "beams in probe wave function= " << nbeamp  << endl;
       fp << "C with a resolution (in Angstroms) = " << res  << endl;
       if( lwobble == 1 ) {
              fp << "C Number of thermal configurations = " << nwobble << endl;
              fp << "C Source size = " << sourceFWHM <<" Ang. (FWHM) " << endl;
       }
       if( TRUE == labErr ) {
          fp << "C  add pi/4 aberr. tuning errors" << endl;
	   }
       fp << "C     x      y     signal" << endl;

       //  remember there is only one thickness for 1D mode
       for( ip=0; ip<nyout; ip++) {
           /* recalculate mean x,y without source size wobble */
           x[ip] = xi + dx * ((double)ip);
           y[ip] = yi + dy * ((double)ip);
           fp << setw(14) << x[ip] << " " << setw(14) << y[ip];
           for(i=0; i<ndetect; i++)
               fp << " " << setw(14) << pixr[i+(nThick-1)*ndetect][0][ip];
           fp <<  endl;
       }

       fp.close();

       /* ------ next output xz image data ---------------- */

       if( lxzimage == 1 ) {
    
           /*  directory file listing parameters for each image file */
            fileout = fileoutpre + ".txt";
            fp.open( fileout.c_str() );
            if( fp.bad() ) {
                cout << "Cannot open output file " << fileout << endl;
                exit( 0 );
            }
   
            /*  store params plus min and max */
            param[pIMAX]    = 0.0F;
            param[pIMIN]    = 0.0F;
            param[pXCTILT]  = (float) ctiltx;
            param[pYCTILT]  = (float) ctilty;
            param[pDEFOCUS] = (float) df;
            param[pDX]  = (float) dx;
            param[pDY]  = (float) dy;
            param[pENERGY]  = (float) keV;
            param[pOAPERT]  = (float) apert2;
            param[pWAVEL]   = (float) wavlen;
            param[pNSLICES] = (float) nslice;
    
            myFile.setnpix( 1 );
            myFile.resize( nyout, nThick );
            aimin = aimax = 0.0F;
            for( i=0; i<NPARAM; i++) myFile.setParam( i, param[i] );
  
            for( idetect=0; idetect<ndetect; idetect++){
                // sprintf( fileout, "%s%03d.tif", fileoutpre, idetect );  // plain C
                fileout = fileoutpre + toString(idetect) + ".tif";
                cout << "output file= " << fileout << endl;
    
                /* convert to float and fix pixel order */
                pmin = pmax = (float) pixr[idetect][0][0];
                for( ix=0; ix<nyout; ix++)
                for( iy=0; iy<nThick; iy++) {
                    temp = myFile( ix, iy ) = pixr[idetect+iy*ndetect][0][ix];
                    if( temp < pmin )pmin = temp;
                    if( temp > pmax )pmax = temp;
                }
    
                cout << fileout << ": output pix range : " 
                           << pmin << " to " << pmax << endl;
                myFile.setParam(pRMAX, pmax );
                myFile.setParam(pRMIN, pmin );
                myFile.setParam(pMINDET, 0 );
                myFile.setParam(pMAXDET, 0 );
                myFile.setParam(pCMINDET, 0 );
                myFile.setParam(pCMAXDET, 0 );
                if( collectorMode[idetect] == ADF ) {
                    myFile.setParam(pMINDET, (float) almin[idetect] );
                    myFile.setParam(pMAXDET, (float) almax[idetect] );
                } else if( collectorMode[idetect] == CONFOCAL ) {
                    myFile.setParam(pCMINDET, (float) almin[idetect] );
                    myFile.setParam(pCMAXDET, (float) almax[idetect] );
                }
                if(  myFile.write( fileout.c_str(), pmin, pmax,
                    aimin, aimax, (float) dx, (float) dy ) != 1 ) {
                        cout << "Cannot write output file " << fileout << endl;
                }
                if( ADF == collectorMode[idetect]  )
                    fp << "file: " << fileout << ", detector= "
                        << almin[idetect]*1000.0 << " to " << almax[idetect]*1000.0 
                        << " mrad, range= " << pmin 
                        << " to " << pmax << endl;
                else if( CONFOCAL == collectorMode[idetect] )
                    fp <<"file: " << fileout
                        << ", detector= " << almin[idetect] << " to " << almax[idetect]  
                        << " Angst., range= " << pmin
                        << " to " << pmax << endl;
           }
            fp.close();
 
       }  /* end if( lxzimage==1... */

    } /* end if( l1d.. ) */

    if( nbeamt > 0 ) {
        cout << "Number of symmetrical anti-aliasing "
            "beams in trans. function = " << nbeamt << endl;
    }

    /*  echo min/max of total integrated intensity */
    cout << "The total integrated intensity range was:\n" << endl;
    cout << "   " << totmin << " to " << totmax << "\n" << endl;

    cout << "CPU time = " << cputim()-timer << " sec." << endl;
#ifdef USE_OPENMP
    cout << "wall time = " << walltim() - walltimer << " sec." << endl;
#endif

    return( EXIT_SUCCESS );

}  /* end main() */

