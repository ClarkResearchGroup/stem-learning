/*
      *** autosliccmd.cpp ***        (normal C++)
      *** autosliccmd_cuda.cpp ***   (CUDA code version)

Front end for autoslic.cpp containing a command line user interface.

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
  gcc -O -fopenmp -o autoslic autoslic.c slicelib.o
                       tiffsubs.o  cfpix.o -lfftw3f

  Transmit an electron wave through a specimen using the
  multislce method with automatic slicing.  Read in the (x,y,z)
  coordinates of the whole specimen and break into slices
  on-the-fly.

  started 24-july-1996 E. Kirkland
  working 19feb-1997 ejk
  last revised 19-feb-1997 ejk
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
  change RNG seed argument to referenece so it get updated for 
      successive calls 21-sep-2013 ejk
  fix %ld to %d format in printout of aslice.nillum  12-oct-2013 ejk
  add param[] for mode 15-oct-2013 ejk
  convert to streams and strings 22-mar-2014 ejk
  fix minor formatting issue 2-jun-2014 ejk
  add optional aberration phase output 24,31-aug-2014 ejk
  add 1 to slice for beams output because it start after first slice
       21-feb-2015 ejk
  convert malloc1D() etc to vector<> 2-jul-2016 ejk
  add trap for bad xyz file read failure 30-jul-2016 ejk
  convert to TDS multi-threading model 6-aug-2016 ejk
  add CBED with frozen phonon TDS 18,28-aug-2016 ejk
  add normal diffraction with TDS 15-jul-2017 ejk
  fix dx,dy error for CBED+diffract (must invert) 
      and add trap for ax,by,cz <= 0 on 23-aug-2018 ejk
  add option to output CBED on log scale and
       azimuthal average 26-aug-2018 ejk
  start GPU/cuda mode 8-oct-2018, working 21-dec-2019 ejk
  fix dx,dy bug on preview pix output 14-feb-2019 ejk
  add record of z coord. for each beam values 16-mar-2019 ejk
  put back verbose mode in exit wave mode 28-mar-2019 ejk

  ax,by,cz  = unit cell size in x,y
  acmin  = minimum illumination angle
  acmax  = maximum illumination angle
  Cs     = spherical aberration coefficient
  df0    = defocus (mean value)
  sgmaf = defocus spread (standard deviation)
  dfdelt = sampling interval for defocus integration
  
  this file is formatted for a TAB size of 4 characters 
  
*/

#include <cstdio>  /* ANSI C libraries */
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <ctime>

#include <string>
#include <iostream>  //  C++ stream IO
#include <fstream>
#include <iomanip>   //  to format the output
#include <vector>   // STD vector class

using namespace std;

#include "cfpix.hpp"       // complex image handler with FFT
#include "slicelib.hpp"    // misc. routines for multislice
#include "floatTIFF.hpp"   // file I/O routines in TIFF format

//  use only one;  the calculation part
//#include "autoslic_cuda.hpp"    //  header for cuda nvcc version of this program
#include "autoslic.hpp"    //  header for C++ version of this class

#define MANY_ABERR      //  define to include many aberrations

#ifdef USE_OPENMP
#include <omp.h>
/*  get wall time for benchmarking openMP */
#define walltim() ( omp_get_wtime() )
double walltimer;
#endif

const int NSMAX= 1000;   // max number of slices
const int NCMAX= 1024;   // max characters in file names
const int NZMAX= 103;    // max atomic number Z

enum{ TRUE=1, FALSE=0};

int main()
{
    string filein, fileout, filestart, filebeam, filecross, cline, description;
  
    const string version = "28-mar-2019 (ejk)";

    //  change labPhase to calculate the aberration phase image
    //int lstart=0, lpartl=0, lbeams=0, lwobble=0, lcross=0, nwobble=1, labPhase=1;
    int lstart=0, lpartl=0, lbeams=0, lwobble=0, lcross=0, nwobble=1, labPhase=0,
        ldiffract=0, logpix=0, lCBED=0;
    int ix, iy, iz, nx, ny, nzout, i, nslic0, islice, nsum, nh, verbose,
        ndf, nbout, ib, ncellx, ncelly, ncellz, NPARAM, ixmid, iymid;
    int nillum, nzbeams;
    int natom, done, status, multiMode;
    long  ltime;

    unsigned long iseed;

    float v0, mm0, wavlen, ax, by, cz, pi,
        rmin, rmax, aimin, aimax, ctiltx, ctilty,
        acmin, acmax, Cs3, Cs5, df0, sigmaf, dfdelt, aobj,
        temperature, ycross, dx, dy, xp, yp, pixc;

    const float tiny=1.0e-4F;
   
    float wmin, wmax, xmin,xmax, ymin, ymax, zmin, zmax;

    double timer, deltaz, vz;
    double sum, rx, ry, ry2;

    vector<int> nhist;

    vectori hbeam, kbeam, Znum;
    vectorf x, y, z, occ, wobble, param, sparam;

    vector<double> hist;
        
    ofstream fp1, fp;

    cfpix pix;      // to get results of calculation
    cfpix wave0;    // initial wavefunction (if used)
    cfpix depthpix; // to get xz cross section results
    cfpix beams;    // to get valuse of requested beams during propagation

    floatTIFF myFile;   //  file input/output
    autoslic aslice;    // has the calculation engine

/*  echo version date and get input file name */

    cout << "autoslic(e) version dated " << version << endl;
    cout << "Copyright (C) 1998-2019 Earl J. Kirkland"  << endl;
    cout << "This program is provided AS-IS with ABSOLUTELY NO WARRANTY"  << endl;
    cout <<     " under the GNU general public license" << endl << endl;

    cout << "perform CTEM multislice with automatic slicing and FFTW" << endl;
#ifdef USE_OPENMP
    cout << "and multithreaded using openMP" << endl;
#endif
    cout <<  " "  << endl;

    pi = (float) (4.0 * atan( 1.0 ));
    NPARAM = myFile.maxParam();
    param.resize( NPARAM );
    sparam.resize( NPARAM );
    for( ix=0; ix<NPARAM; ix++ ) param[ix] = 0.0F;

    cout << "Name of file with input atomic "
           << "coord. in x,y,z format:"  << endl;
    cin >> filein;

/*  get simulation options */

    cout << "Replicate unit cell by NCELLX,NCELLY,NCELLZ :" << endl;
    cin >> ncellx >> ncelly >> ncellz;
    if( ncellx < 1 ) ncellx = 1;
    if( ncelly < 1 ) ncelly = 1;
    if( ncellz < 1 ) ncellz = 1;

    cout << "Name of file to get binary output of multislice result:" << endl;
    cin >> fileout ;

    lstart = 0;
    //  this if() is to convoluted; should fix sometime (15-jul-2017 ejk)
    lpartl = askYN("Do you want to include partial coherence");
    if( 0 == lpartl ) {
        lCBED = askYN("Do you want to calculate CBED with TDS");
        if( 0 == lCBED ) {
            ldiffract = askYN("Do you want to calculate diffraction with TDS");
        }
    }

    //  this if() is too convoluted - should work on this sometime
    acmin = acmax = 0;
    xp = yp = 0;
    if( (1 == lpartl) || ( 1 == lCBED) ) {
        if( 0 == lCBED ) {
            cout << "Illumination angle min, max in mrad:" << endl;
            cin >> acmin >> acmax;
            acmin  = acmin  * 0.001F;
            acmax  = acmax  * 0.001F;
        }
        cout << "Spherical aberration Cs3, Cs5(in mm.):" << endl;
        cin >>  Cs3 >> Cs5;
        param[pCS]  = (float) ( Cs3*1.0e7 );
        param[pCS5] = (float) ( Cs5*1.0e7 );
        if( 0 == lCBED ) {
            cout << "Defocus, mean, standard deviation, and"
                 " sampling size (in Angstroms) =" << endl;
            cin >> df0 >> sigmaf >> dfdelt;
        } else {
            cout << "Defocus (in Angstroms) =" << endl;
            cin >> df0 ;
            sigmaf = 0;
            dfdelt = 0;
        }
        param[pDEFOCUS] = (float) df0;
        param[pDDF] = (float) sigmaf;
        cout << "Objective aperture (in mrad) =" << endl;
        cin >> aobj;
        aobj = (float) fabs( aobj * 0.001F);
#ifdef MANY_ABERR
        /*   get higher order aberrations if necessary */
        cout << "type higher order aber. name (as C32a, etc.) followed\n"
            << " by a value in mm. (END to end)" << endl;
        done = multiMode = 0;
        do{
            cin >> cline ;
            if( cline.compare( "END" ) == 0  ) {
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
        lstart = 0;
        if( 1 == lCBED ) {
            cout << "CBED probe postion (in Ang.):" << endl;
            cin >> xp >> yp;
        }
    } else if( (0 == lCBED) && (0 == ldiffract) ) {
        cout << "NOTE, the program image must also be run." << endl;
        lstart = askYN("Do you want to start from previous result");
    }

    if ( lstart == 1 ) {
        cout << "Name of file to start from:" << endl;
        cin >> filestart;
    } else {
        cout << "Incident beam energy in kev:" << endl;
        cin >> v0;
        cout << "Wavefunction size in pixels, Nx,Ny:" << endl;
        cin >> nx >> ny;
    }

    cout << "Crystal tilt x,y in mrad.:" << endl;
    cin >> ctiltx >> ctilty;
    ctiltx = ctiltx /1000;
    ctilty = ctilty /1000;

    /*  remember that the slice thickness must be > atom size
        to use projected atomic potential */
    cout << "Slice thickness (in Angstroms):" << endl;
    cin >> deltaz;
    if( deltaz < 1.0 ) {
        cout << "WARNING: this slice thickness is probably too thin"
            << " for autoslice to work properly." << endl;
    }

    if( (lpartl == 0) && (0 == lCBED) && ( 0 == ldiffract)) {
        lbeams = askYN("Do you want to record the (real,imag) value\n"
        " of selected beams vs. thickness");
        if( lbeams == 1 ) {
            cout << "Name of file for beams info:" << endl;
            cin >> filebeam;
            cout << "Number of beams:" << endl;
            cin >> nbout;
            if( nbout<1 ) nbout = 1;
            hbeam.resize( nbout );
            kbeam.resize( nbout );
            for( ib=0; ib<nbout; ib++) {
                cout << "Beam " << ib+1 << ", h,k=" << endl;
                cin >> hbeam[ib] >> kbeam[ib];
            }
        }
    }

    lwobble = askYN("Do you want to include thermal vibrations");
    if( lwobble == 1 ) {
        cout << "Type the temperature in degrees K:" << endl;
        cin >> temperature ;
        cout << "Type number of configurations to average over:" << endl;
        cin >>  nwobble;
        if( nwobble < 1 ) nwobble = 1;
        /* get random number seed from time if available 
            otherwise ask for a seed */
        ltime = (long) time( NULL );
        iseed = (unsigned) ltime;
        if( ltime == -1 ) {
            cout << "Type initial seed for random number generator:" << endl;
            cin >>  iseed;
        } else
            cout << "Random number seed initialized to " << iseed << endl;
    } else temperature = 0.0F;

    if( (0 == lpartl)  && (0 == lCBED) && (0 == ldiffract) ) {
        lcross = askYN("Do you want to output intensity vs. depth cross section");
        if( lcross == 1 ){
            cout << "Type name of file to get depth profile image:" << endl;
            cin >> filecross;
            cout << "Type y position of depth cross section (in Ang.):" << endl;
            cin >> ycross;
        }
    }

    if( (1 == lCBED) || (1 == ldiffract) ) {
        logpix = askYN("Do you want to output intensity on a log scale");
    }

/* start timing the actual computation just for fun */

    timer = cputim();
#ifdef USE_OPENMP
    walltimer = walltim();  /* wall time for openMP */
#endif

/*  get starting value of transmitted wavefunction if required
   (this can only be used in coherent mode)
    remember to save params for final output pix  */

    if ( lstart == 1 ) {
        if( myFile.read( filestart.c_str() ) != 1 ) {
            cout << "Cannot open input file: " << filestart << endl; 
            exit( 0 );
        }

        if( myFile.getnpix() != 2 ) {
           cout << "Input file " << filestart << " must be complex, can't continue." << endl;
           exit( 0 );
        }

        nx =  myFile.nx();
        ny =  myFile.ny();

        nx = nx/2;
        wave0.resize( nx, ny );

        //  save starting pix for later
        for( ix=0; ix<nx; ix++) for( iy=0; iy<ny; iy++) {
                wave0.re(ix,iy) = myFile(ix,iy);
                wave0.im(ix,iy) = myFile(ix+nx,iy);
        }

        //  save parameters to verify successive images are same size etc.
        for( i=0; i<NPARAM; i++) sparam[i] = myFile.getParam( i );

        ax = sparam[pDX] * nx;
        by = sparam[pDY] * ny;
        v0     = sparam[pENERGY];
        nslic0 = (int) sparam[pNSLICES];
        cout << "Starting pix range " << sparam[pRMIN] << " to " << sparam[pRMAX] 
              << " real\n" << "           " << sparam[pIMIN] << " to " 
              << sparam[pIMAX] << " imag" << endl;
        cout << "Beam voltage = " << v0 << " kV" << endl;
        cout << "Old crystal tilt x,y = " << 1000.*sparam[pXCTILT] <<  ", " 
            << 1000.*sparam[pYCTILT] << " mrad" << endl;

    } else nslic0 = 0;     /* end if( lstart...) */

/*  calculate relativistic factor and electron wavelength */

    mm0 = 1.0F + v0/511.0F;
    wavlen = (float) wavelength( v0 );
    cout << "electron wavelength = " << wavlen << " Angstroms" << endl;

/*  read in specimen coordinates and scattering factors */

    natom = ReadXYZcoord( filein.c_str(), ncellx, ncelly, ncellz,
        &ax, &by, &cz, Znum, x, y, z, occ, wobble,
        description );

    if( natom < 1 ) {
        cout << "no atom coord. read in from file " << filein << endl;
        cout << "can not continue...." << endl;
        return EXIT_SUCCESS;
    }

    if( (ax<tiny) || (by<tiny) || (cz<tiny) ) {
        cout << "Bad unit cell size ax,by,cz = " << ax << ", " << by << ", " << cz << endl;
        exit( 0 );
    }

    cout << natom << " atomic coordinates read in"  << endl;
    cout << description << endl;

    cout <<"Size in pixels Nx, Ny= " << nx << " x " << ny << " = " << nx*ny 
        << " beams" << endl;
    cout <<"Lattice constant a,b = " << ax << ", " << by << endl;

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
    cout << "Total specimen range is\n" 
        << xmin << " to " << xmax << " in x\n"
            << ymin << " to " << ymax << " in y\n"
        << zmin << " to " << zmax << " in z" << endl;
    if( lwobble == 1 )
        cout << "Range of thermal rms displacements (300K) = "
            << wmin << " to " << wmax << endl;

    // ---------  setup calculation -----
    //   set calculation flags
    aslice.lbeams = lbeams;
    aslice.lcross = lcross;
    aslice.lpartl = lpartl;
    aslice.lstart = lstart;
    aslice.lwobble = lwobble;

    //   set calculation parameters (some already set above)
    param[ pAX ] = ax;          // supercell size
    param[ pBY ] = by;
    param[ pNX ] = (float) nx;
    param[ pNY ] = (float) ny;
    param[pENERGY]   =  v0;
    param[ pDELTAZ ] = (float) deltaz;  // slice thickness
    param[ pOAPERT ] = aobj;
    param[ pXCTILT ] = ctiltx;      // crystal tilt
    param[ pYCTILT ] = ctilty;
    param[pCAPERT] = acmax;     // condencer angles
    param[pCAPERTMIN] = acmin;
    param[ pTEMPER ] = (float) fabs( temperature );
    param[ pNWOBBLE ] = (float) nwobble;    //  number config. to average
    param[pWAVEL] = wavlen;         //  probably recal. autoslice::calculate()
    param[pPPOSX] = xp;     //  CBED probe position
    param[pPPOSY] = yp;

    param[pMODE] = 6;  // save mode = autoslic

    if ( lpartl == 1 ) {
            param[pDEFOCUS] = df0;
            param[pOAPERT] = aobj;
            param[pDDF] = sigmaf;
            param[pCAPERT] = acmax;

            rx = 1.0F/ax;
            ry = 1.0F/by;
            cout << "Illumination angle sampling (in mrad) = "
                << 1000.*rx*wavlen << ", " << 1000.*ry*wavlen << "\n" << endl;
    }

    // ------- iterate the multislice algorithm proper -----------

    if( (0 == lpartl) && (0 == lCBED) && (0 == ldiffract)) {
        cout << "calculate exit wave function" << endl;
        aslice.initAS( param, Znum, natom );
        pix.resize( nx, ny );
        pix.init( );
        verbose = 1;
        aslice.calculate( pix, wave0, depthpix, param, multiMode, natom,
                Znum, x,y,z,occ, beams, hbeam, kbeam, nbout, ycross, verbose);
    } else if( (1 == lpartl) && (0 == lCBED) ) {
        cout << "calculate pix with partial coherence" << endl;
        nbout = 0;
        aslice.calculatePartial( pix, param, multiMode, natom, &iseed,
                Znum, x,y,z,occ,wobble, dfdelt );
    } else if( (1 == lCBED) && (0 == ldiffract) ) {
        cout << "calculate CBED" << endl;
        nbout = 0;
        aslice.lcbed = 1;
        aslice.calculateCBED_TDS( pix, param, multiMode, natom, &iseed,
                Znum, x,y,z,occ,wobble );
    }else if( (0 == lCBED) && (1 == ldiffract) ) {
        cout << "calcuate diffraction" << endl;
        nbout = 0;
        aslice.lcbed = 0;
        aslice.calculateCBED_TDS( pix, param, multiMode, natom, &iseed,
                Znum, x,y,z,occ,wobble );
    }
 
 
    if( lpartl == 1 ) {         //    with partial coherence
        nillum = aslice.nillum;
        cout << "Total number of illumination angle = "
                << nillum << endl;
        ndf = (int) ( ( 2.5F * sigmaf ) / dfdelt );  // recal same value in class
        cout << "Total number of defocus values = " << 2*ndf+1 << endl;
    }
   
    else if( lbeams ==1 ) {
            fp1.open( filebeam.c_str() );
            if( fp1.bad() ) {
                cout << "can't open file " << filebeam << endl;
                exit(0);
            }
            fp1 << " (h,k) = " ;
            for(ib=0; ib<nbout; ib++)
                fp1 << " (" << hbeam[ib] << "," << kbeam[ib] << ")";
            fp1 << endl;
            nzbeams = beams.ny();
            fp1 << "zslice(Ang.), (real,imag) (real,imag) ...\n" << endl;
            for( islice=0; islice<nzbeams; islice++) {
                if( beams.im(nbout, islice) > 0 ) {
                    fp1 << setw(10) << beams.re(nbout, islice);  // actual z coord.
                    for( ib=0; ib<nbout; ib++)  //  setprecision(4)
                        fp1 << "  " << setw(10) << beams.re(ib,islice)
                            << "  " << setw(10) << beams.im(ib,islice); 
                    fp1 << endl;
                }
            }
        fp1.close();

    } // end else 

/*  ------------------------------------------------------
    output results and find min and max to echo
*/

    //if( ldiffract > 0 ) {
    //  pix.re(nx/2,ny/2) = 0;  //  zero center beam ; doesn't help much (??)
    //  cout << "zero out the center beam (diffraction mode)" << endl;
    //}

    pix.findRange( rmin, rmax, aimin, aimax );

    // calculate azimuthal average, convert to log scale if requested
    if ( ( 1 == lCBED) || ( 1 == ldiffract) ) {

        // histogram the azimutal average
        hist.resize( (nx+ny) );
        nhist.resize( (nx+ny) );
        for( ix=0; ix<(nx+ny); ix++) {
            hist[ix] = 0.0;
            nhist[ix] = 0;
        }

        sum = 0.0;
        nsum = 0;
        nh = 0;
        ixmid = nx/2;
        iymid = ny/2;

        for( iy=0; iy<ny; iy++) {
            ry2 = (double) ( iy-iymid);
            ry2 = ry2*(ax/by);
            ry2 = ry2*ry2;
            for( ix=0; ix<nx; ix++) {
                pixc = pix.re( ix, iy );
                rx = (double) (ix-ixmid);
                i = (int) ( sqrt( rx*rx + ry2 ) + 0.5);
                hist[i] += pixc;
                nhist[i]++;
                if( i > nh ) nh = i;
                if( TRUE == logpix ) {
                    if( pixc > 1.e-10F)  pixc = 
                        (float) log( (double) fabs(pixc) );
                    else pixc = -23.0F;
                    pix.re( ix, iy ) = pixc;
                }
                if( (ix == 0) && (iy == 0) ) {
                    rmin = pixc;
                    rmax = rmin;
                } else if( (ix != ixmid) && (iy != iymid) ) {
                    if( pixc < rmin ) rmin = pixc;
                    if( pixc > rmax ) rmax = pixc;
                }
                if( (ix>(3*nx)/8) && (ix<(5*nx)/8) &&
                    (iy>(3*ny)/8) && (iy<(5*ny)/8) ) {
                    sum = sum + pixc;
                    nsum += 1;
                }

            }  /* end for ix... */
        } /* end for iy... */

//        rmin= (float) (0.05*rmin + 0.95*sum/nsum);   //  sometimes better
        rmin= (float) (0.10*rmin + 0.90*sum/nsum);   //  sometimes better

        cout << "write azimuthal averaged intensity vs. \n"
            "  spatial frequency k, into file azimuth.dat" << endl;
        fp.open( "azimuth.dat" );
        if( fp.bad() ) {
            cout << "cannot open file azimuthal.dat" << endl;
            exit( 0 );
        }
        for( i=0; i<=nh; i++) {
            hist[i] = hist[i] / nhist[i];
            fp << setw(16) << ((double)i)/ax << setw(16) << hist[i] << endl;
        }
        fp.close( );

    }  //  end if (  (1 == logpix).....

    param[pRMAX]  = rmax;
    param[pIMAX]  = aimax;
    param[pRMIN]  = rmin;
    param[pIMIN]  = aimin;

    if (  ( 1 == lCBED) || ( 1 == ldiffract) ) {
        param[pDX] = dx = (float) ( 1.0/ax );  //  for FFT space !!
        param[pDY] = dy = (float) ( 1.0/by );
    } else {
        param[pDX] = dx = (float) ( ax/((float)nx) );
        param[pDY] = dy = (float) ( by/((float)ny) );
    }

    //param[pNSLICES] = 0.0F;  /* ??? */

    for( ix=0; ix<NPARAM; ix++ ) myFile.setParam( ix, param[ix] );

    if ( (lpartl == 1) || ( 1 == lCBED) || ( 1 == ldiffract) ) {
        myFile.resize( nx, ny );
        myFile.setnpix( 1 );
        for( ix=0; ix<nx; ix++) for( iy=0; iy<ny; iy++)
            myFile(ix,iy) = pix.re(ix,iy);
    } else {
        myFile.resize( 2*nx, ny );
        myFile.setnpix( 2 );
        for( ix=0; ix<nx; ix++) for( iy=0; iy<ny; iy++) {
            myFile(ix,iy)    = pix.re(ix,iy);
            myFile(ix+nx,iy) = pix.im(ix,iy);
        }
    }

    i = myFile.write( fileout.c_str(), rmin, rmax, aimin, aimax, dx, dy );
    if( i != 1 ) cout << "autoslice cannot write TIF file " << fileout << endl;
    cout << "pix range " << rmin << " to " << rmax << " real,\n" <<
            "          " << aimin << " to " << aimax << " imag" << endl;

    /* ----- output depth cross section if requested ------- */
    if( lcross == 1 ){
        depthpix.findRange( rmin, rmax, aimin, aimax );
        myFile.setParam( pRMAX, rmax );
        myFile.setParam( pIMAX, 0.0F );
        myFile.setParam( pRMIN, rmin );
        myFile.setParam( pIMIN, 0.0F );
        myFile.setParam( pDY, dy = (float) ( deltaz ) );

        nzout = depthpix.ny();
        myFile.resize( nx, nzout );
        myFile.setnpix( 1 );
        for( ix=0; ix<nx; ix++) for( iz=0; iz<nzout; iz++) {
            myFile(ix,iz) = depthpix.re(ix,iz);
        }
        i = myFile.write( filecross.c_str(), rmin, rmax, aimin, aimax, dx, dy );

        if( i != 1 ) cout << "autoslice cannot write TIF file "
                << filecross << endl;
        cout << "depth pix range " << rmin << " to " << rmax << " real" << endl;
    }

    /* ----- output aberration function phase if requested ------- */
    if( (lpartl == 1) && (labPhase == 1) ){

        aslice.abbPhase2D( pix, param, multiMode );
        pix.findRange( rmin, rmax, aimin, aimax );
        myFile.setParam( pRMAX, rmax );
        myFile.setParam( pIMAX, 0.0F );
        myFile.setParam( pRMIN, rmin );
        myFile.setParam( pIMIN, 0.0F );
        dx = (float) ( 1.0/ax );  //  for FFT space !!
        dy = (float) ( 1.0/by );
        myFile.setParam( pDX, dx );
        myFile.setParam( pDY, dy );

        myFile.resize( nx, ny );
        myFile.setnpix( 1 );
        for( ix=0; ix<nx; ix++) for( iy=0; iy<ny; iy++) {
            myFile(ix,iy) = pix.re(ix,iy);
        }
        i = myFile.write( "abbPhase.tif", rmin, rmax, aimin, aimax, dx, dy );

        if( i != 1 ) cout << "autoslice cannot write TIF file "
                << filecross << endl;
        cout << "ab. phase pix range " << rmin << " to " << rmax << endl;
    }


    cout << "Total CPU time = " << cputim()-timer << " sec." << endl;
#ifdef USE_OPENMP
    cout << "wall time = " << walltim() - walltimer << " sec." << endl;
#endif

    return EXIT_SUCCESS;

} /* end main() */



