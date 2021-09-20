/*      incostemCMD.cpp

------------------------------------------------------------------------
Copyright 1998-2016 Earl J. Kirkland

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

  command line based front end for incostem calulation
  calculate images in the incoherent STEM approximation

  put the partial cross section (integrated over the ADF detector
  angles) at a single pixels (position of the corresponding atom)
  and convolve with the point spread function (focused probe intensity)

  reference:

  [1] E. Kirkland, "Advanced Computing in Electron Microscopy",
        Plenum 1998, 2nd edit. Springer 2010
  
  this file is formatted for a tab size of 4 char

  started 16-mar-1998 E. Kirkland
  change float1D() etc to malloc1D() etc. and get rid of
     memory.h 18-jun-2007 ejk
  add cputim() echo 27-jun-2007 ejk
  put in adaptive quadrature for r and k integration 
        - faster and more accurate 28-jun-07 ejk
  add source size 13-jul-2008 ejk
  start adding defocus spread (Cc effects) 1-mar-2009 ejk
  add Cs5 on 7-mar-2009 ejk
  add Gauss-Hermite integration of defocus spread 15-mar-2009 ejk
  change wording on initial dialog 14-nov-2009 ejk
  fix bug in FWHM conversion (factor of two error) 13-nov-2010 ejk
  start convert to FFTW  8-feb-2011 ejk
  adding many aberration mode (C32a etc)  13-feb-2011
  get return value of scanf() to remove warnings from gcc 4.4
      13-feb-2011 ejk
  convert many aberration mode to use chi() from slicelib
      and fix a few small typos 21-may-2011 ejk
  convert &lf to &lg reading multipole aberration and fix
      small typo in comments 20-jun-2011 ejk
  fix small typo's in comments and formatting 27-jul-2011 ejk
  convert to C++ and floatTIFF.cpp  14-aug-2012 ejk
  convert to cfpix/fftw class from raw fftw 26-feb-2013 to 30-oct-2012 ejk
  fix typo in final param[] saving loop  17-mar-2013 ejk 
  fix typo in final param[] updates - keep existing aberr.
     and get rid of redundant variable wavelen (use wavl) 30-mar-2013 ejk
  move main calculation into subroutine incostemCal() so I can put this 
    whole thing somewhere else if needed 13-apr-2013 ejk
  separate calculation into a separate class for easy use
     in other places 20-apr-2013 ejk
  convert calculate() to calculate2D() 25-jul-2013 ejk
  add param[] for mode 17-oct-2013 ejk
  convert to streams and strings 8-apr-2014 ejk
  add scaling for probe current and dwell time to calculate
      noise in image, and fix one leftover printf  26-oct-2014 ejk
  add param[] elements for pPROBET and pPROBEDT 11-oct-2015 ejk
  relinked with new incostem.cpp 30-jan-2016 ejk
  convert malloc1D() to vector<> 5-jul-2016 ejk
*/

#include <stdio.h>  /* standard ANSI libraries */
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <string>
#include <iostream>  //  C++ stream IO
#include <fstream>
#include <iomanip>   //  to format the output

using namespace std;

#include "cfpix.hpp"       /* complex image handler with FFT */
#include "slicelib.hpp"    /* misc. routines for multislice */
#include "floatTIFF.hpp"   /* file I/O routines in TIFF format */

#include "incostem.hpp"    //  the calculations 

#define MANY_ABERR      /*  define to include many aberrations */

int main()
{
    string filein, fileout, description, cline;

    int ix, iy, nx, ny, i, ncellx, ncelly, ScaleProbe, np,
        ncellz, done, natom, multiMode=0, status, NPARAM;

	vectori  Znum;
    vectorf x, y, z, occ, wobble, param;
   
    float wmin, wmax, xmin,xmax, ymin, ymax, zmin, zmax;
    float fdx, fdy, ax, by, cz, rmin, rmax, aimin, aimax;
    
    double vz, Cs3, Cs5, df0, ddf, aobj;
    double  keV, thetamin, thetamax, timer, dsource, pi;
    double twopi, sigmae, wavl, probeI, dwellTime, scalep;

    cfpix pix;  //  complex floating point images

    floatTIFF myFile;  // file handler

    incostem inc;

/* ------ echo version date and get input file name ---------- */

    cout << "incostem version dated 6-jul-2016 ejk" << endl
        << "calculate ADF-STEM images for very thin specimens" << endl
        << "using the incoherent image model" << endl;
    cout <<"Copyright (C) 1998-2016 Earl J. Kirkland" << endl;
    cout << "This program is provided AS-IS with ABSOLUTELY NO WARRANTY "
        << " under the GNU general public license" << endl << endl;

    NPARAM = myFile.maxParam();
    pi = (float) (4.0 * atan( 1.0 ));

    param.resize( NPARAM );
    for( i=0; i<NPARAM; i++) param[i] = 0.0f;

    cout << "Name of file with input atomic "
           << "potential in x,y,z format:" << endl;
    cin >> filein ;

/* ------ get simulation options ------ */

    cout << "Replicate unit cell by NCELLX,NCELLY,NCELLZ :" << endl;
    cin >> ncellx >> ncelly >> ncellz;
    if( ncellx < 1 ) ncellx = 1;
    if( ncelly < 1 ) ncelly = 1;
    if( ncellz < 1 ) ncellz = 1;

    cout << "Name of file to get binary output of adf-stem result:" << endl;
    cin >> fileout ;

    /*  any size works with FFTW so don't need to force power of 2 */ 
    cout << "Image size, Nx,Ny in pixels :" << endl;
    cin >> nx >> ny;

    cout << "STEM probe parameters, V0(kv), Cs3(mm), Cs5(mm),"
          << " df(Angstroms), apert(mrad):" << endl;
    cin >> keV >> Cs3 >> Cs5 >> df0 >> aobj;
    aobj *= 0.001;
    param[pDEFOCUS] = (float) df0;
    param[pCS]  = (float) ( Cs3*1.0e7 );
    param[pCS5] = (float) ( Cs5*1.0e7 );

    cout << "ADF detector angles thetamin, thetamax (in mrad) =" << endl;
    cin >> thetamin >> thetamax;
    thetamin *= 0.001;
    thetamax *= 0.001;
    param[pMINDET] = (float) thetamin;
    param[pMAXDET] = (float) thetamax;

    /*   get higher order aberrations if necessary */
    done = multiMode = 0;
#ifdef MANY_ABERR
    /*   get higher order aberrations if necessary */
    cout << "type higher order aber. name (as C32a, etc.) followed\n"
        << " by a value in mm. (END to end)" << endl;
    do{
        cin >> cline ;
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

    cout <<  "Source size at specimen (FWHM in Ang.):"<< endl;
    cin >> dsource;
    param[pSOURCE] = (float) dsource;

    cout <<  "Defocus spread (FWHM in Ang.):"<< endl;
    cin >> ddf;

    probeI = dwellTime = 0.0;
    ScaleProbe = askYN( "Do you want to add electron counting noise" );
    if( ScaleProbe > 0 ){
        cout << "Type total probe current (in pAmp) and dwell time (in microSec):" << endl;
        cin >> probeI >> dwellTime;
    }
    
    /* start timing the actual computation just for fun */
    timer = cputim();

/* ------ calculate relativistic electron wavelength ------ */

    twopi = 2.0 * 4.0 * atan( 1.0 );
    wavl = wavelength( keV );
    sigmae = sigma( keV )/ 1000.0;

    cout << "electron wavelength = " << wavl <<" Angstroms" << endl;

/* ------ read in specimen coordinates and scattering factors ------ */

    natom = ReadXYZcoord( filein.c_str(), ncellx, ncelly, ncellz,
        &ax, &by, &cz, Znum, x, y, z, occ, wobble,
        description );

    cout << natom << " atomic coordinates read in" << endl;
    cout << description  << endl;

    if( natom < 1 ) {
        cout << "Need more atoms to continue, exiting...." << endl;
        return EXIT_SUCCESS;
    }

    cout << "Lattice constant a= " << ax << ", b= " << by << endl;

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
    cout << "Total specimen range is " << xmin << " to " << xmax << " in x\n"
           << ymin << " to " << ymax << " in y\n"
           << zmin << " to " << zmax << " in z" << endl;

    //-------------------------------------------------
    // ---  main calculation ---------------
    param[pDEFOCUS] = (float) df0;
    param[pDDF] = (float) ddf;
    param[pWAVEL] = (float) wavl;
    param[pOAPERT] = (float) aobj;
    param[pENERGY] = (float) keV;
    param[pAX] = ax;
    param[pBY] = by;
    param[pDX] = fdx = (float) ( ax/((float)nx) );
    param[pDY] = fdy = (float) ( by/((float)ny) );
    param[pNX] = (float) nx;
    param[pNY] = (float) ny;

    param[pMODE] = 8;  // save mode = incostem

    inc.calculate2D( pix,  param, multiMode, natom, Znum, x, y, occ ); 

    //-------------------------------------------------
    // --- output section ---------------

    //   find pix range
    pix.findRange( rmin, rmax, aimin, aimax );

    timer = cputim() - timer;
    cout << "CPU time = " << timer << " sec." << endl;

    cout << "final pix range: " << rmin << " to " << rmax << " (real)\n"
        "          " << aimin << " to " << aimax << " (imag)" << endl;

    if( ScaleProbe > 0 ){
        // remember; 1 pAmp = 6.24146 MHz of electrons
        // incident current in electrons
        scalep = fabs(probeI * dwellTime * 6.24146);
        cout << "scale by " << scalep << endl;
        cout << "scaled pix range: " << scalep*rmin << " to " 
            << scalep*rmax << " (counts)" << endl;

        np = inc.addNoise( pix, nx, ny, probeI, dwellTime );
        cout << "Poisson noise added to " << np << " out of " 
            << nx*ny << " pixels" << endl;
        pix.findRange( rmin, rmax, aimin, aimax );
        cout << "final pix range: " << rmin << " to " 
            << rmax << " (counts)" << endl;
        param[pPROBEI]  = (float) (probeI * 1.0e-12);           //  convert to Amp.
        param[pPROBEDT]  = (float) (dwellTime * 1.0e-6);        //  convert to Sec.
    }

    //  update a extra parameters - remember that most already set
    param[pRMAX]  = rmax;
    param[pIMAX]  = aimax;
    param[pRMIN]  = rmin;
    param[pIMIN]  = aimin;

/*------  write the results into a file ------ */

    for( i=0; i<NPARAM; i++ ) myFile.setParam( i, param[i]);

    cout << "output only the real part" << endl;

    myFile.resize( nx, ny );
    myFile.setnpix( 1 );
    for( ix=0; ix<nx; ix++) for( iy=0; iy<ny; iy++)
        myFile(ix,iy) = pix.re(ix,iy); // save real part for output
    i = myFile.write( fileout.c_str(), rmin, rmax, aimin, aimax, fdx, fdy );

    if( i != 1 ) cout << "incostem cannot write TIF file "
            << fileout  << endl;


#ifdef footest

    double  rfe, ife, k, k2max;
    
    FILE *fp;

/* ------  quick test ------ */

    k = 0.0;
    za = 6;
    keV = 200.0;
    feMoliere( k, za, keV, &rfe, &ife );
    cout << "fe = %g, %g\n", rfe, ife );

/*------  output single atom cross section vs atomic number Z  ------ */

    cout << "write atomcs.dat...\n" );
    fp = fopen( "atomcs.dat", "w+" );

    keV = 200.0;
    wavl = wavelength( keV );
    thetamin = 0.050;
    thetamax = 0.200;
    for( za=1; za<=NZMAX; za++) {
        signal = atomsignal( za, keV, thetamin, thetamax );
        cout << "%5d    %g\n", za, signal );
        fprintf( fp, "%5d    %g\n", za, signal );
    }
    fclose( fp );
#endif

    cout << "all done" << endl;

}  /* end main() */


