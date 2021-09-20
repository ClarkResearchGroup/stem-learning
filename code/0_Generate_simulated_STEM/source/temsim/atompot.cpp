/*
      *** atompot.cpp ***

    WARNING: this program has been superceeded and should not be used
    for future application. It is minimally maintainened, not significantly
    updated and may be discontinued at any time (26-apr-2014).

------------------------------------------------------------------------
Copyright 1998-2017 Earl J. Kirkland

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

  ANSI-C and TIFF version
  this version uses FFTW 3

  FFTW choses an optimum form of the FFT at run time so there
  is some variation in execution speed depending on what else 
  the CPU is doing during this planning stage

  see:   www.fftw.org

  on Windows file libfftw3f-3.dll must be in the PATH

  on Linux build as:
  g++ -O -fopenmp -o atompot atompot.cpp slicelib.o
                       floatTIFF.o -lfftw3f -lm

  Calculate the atomic potentials using tabulated electron
  scattering factors

  The units of the output are such that when multiplied by
  Lambda*(M/M0) the result is the phase shift of the transmitted
  electron wave (i.e. a simple sum of the unadulterated Born
  approx. scattering amplitudes).

  The calculation is performed in a 2D planar slice as a
  prelude to the multislice calculation of electron micrographs.

  Gaussian random displacements may be added to atom coordinate
  to simulate thermal vibrations (i.e. one snap shoot).

    rewritten in C  13-June-1995 ejk
    C version working 28-jun-1995 ejk
    fix bug in natom calculation 7-july-1995
    minor rearrangments 13-July-1995 ejk
    fix small bug in seval on end of range 16-july-1995 ejk
    converted to TIFF file format 8-apr-1996 ejk
    put fetab, featom() in slicelib.c 26-july-1996 ejk
    switch to new fe params 16-jan-1997 ejk
    switched to rangauss() random number generator 22-may-1997 ejk
    removed comma delimiters on input 5-july-1997 ejk
    remove option to set Fourier resolution and fix at max
        25-sep-1997 ejk
    switch to temperature input on thermal vibrations
        3-oct-1997 ejk
    cosmetic changed to compile on MAC/codeWarrior 5-mar-1998 ejk
    changed memory allocator functions 6-nov-1999 ejk
    change void main() to int main() for better portability
         21-jan-2000 ejk
    fix small input read error under new version 
     Linux (set iz=-1 before sscanf()) 22-jan-2000 ejk
    polished code in scamp() a little 8-feb-2006 ejk
    convert to GPL 3-jul-2008 ejk
    get return value of scanf() to remove warnings from gcc 4.4
      and convert to 4 char TAB size formatting 29-mar-2010 ejk
    convert to faster FFTW 4-apr-2010 ejk
    convert to floatTIFF.cpp 18-mar-2012 to 15-apr-2012 ejk
    convert to cfpix/fftw class from raw fftw 14-oct-2012 ejk
    move ReadLine() to here (from sliclib) because its only
        used here now 26-apr-2014 ejk
    convert malloc1D() to vector<> 5-jul-2016 ejk
    convert to streams and strings 2,5-aug-2017 ejk

  This program is ANSI standard C++ and should be transportable
  although this is NOT guaranteed

  This source code is formatted for a tab size of 4.

*/

#include <cstdio>  /* ANSI C libraries */
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <ctime>

#include <string>
#include <iostream>  //  C++ stream IO
#include <sstream>   // string streams
#include <fstream>
#include <iomanip>   //  to format the output
#include <vector>    // STD vector class

using namespace std;

#include "cfpix.hpp"       /* image handler with FFTW 3*/
#include "floatTIFF.hpp"   /* file I/O routines in TIFF format */
#include "slicelib.hpp"    /* misc. routines for multislice */

const int NAMAX  = 2000;    /* max number of atoms of one type */
const int  NSMAX =  20;  /* max number of symmetry operations */
const int  NZMAX =  103; /* max atomic number */

int natom, nsx, nsy;  /* global vars for scamp() */
vectorf x, y, occ;

// define subroutines defined at end of file 
int ReadLine( ifstream &fp, istringstream &ss, string &s );
void scamp( float, float, double*, double* );


int main() 
{
    int ix,iy, nx, ny, ixmid, iymid,
        i, j, nsym, jj, iz, is, lwobble;
    int ncellx, ncelly, ncellz;
    long ncoeff, ltime;
    unsigned long iseed;

    float rmin, rmax, imin, imax, dx, dy;
    double scampr, scampi, sum, runtime, fe, ky2, k2, k2max,
        scale, ax,by,cz,total1, total2, rx2, ry2,
        temperature, scalet;

    string filein, fileot, cline;

    vectorf symx1, symx2, symy1, symy2, kx, ky, wobble;

    cfpix cpix;     //  complex floating point image with FFTW

    ifstream fp;
    istringstream isbuf;

    floatTIFF myFile;

/*  the following are the chemical symbols for the periodic table */

    const string symbol = 
        " HHeLiBe B C N O FNeNaMgAlSi P SCl"
        "Ar KCaScTi VCrMnFeCoNiCuZnGaGeAsSeBr"
        "KrRbSr YZrNbMoTcRuRhPdAgCdInSnSbTe"
        " IXeCsBaLaCePrNdPmSmEuGdTbDyHoErTm"
        "YbLuHfTa WReOsIrPtAuHgTlPbBiPoAtRn"
        "FrRaAcThPa UNpPuAmCmBkCfEsFmMdNoLr"
    ;

    //  Echo version date 
    cout << "atompot version dated 5-aug-2017 EJK"<< endl;
    cout << "Copyright (C) 1998-2017 Earl J. Kirkland" << endl;
    cout <<  "This program is provided AS-IS with ABSOLUTELY NO WARRANTY\n "
        " under the GNU general public license\n" << endl;

    cout << "Warning: this program is obsolete and may be discontinued soon\n"<< endl;

    cout <<  "calculate projected atomic potentials (to use in multislice)"<< endl;
    cout <<  "    using FFTW\n"<< endl;

    //  Get input file name etc. 
    cout << "Name of file with input crystal data :"<< endl;
    cin >> filein;

    cout << "Name of file to get binary output"
          " of atomic potential :"<< endl;
    cin >> fileot;

    cout << "Real space dimensions in pixels Nx, Ny :"<< endl;
    cin >> nx >> ny;

    cout << "Replicate unit cell by NCELLX,NCELLY,NCELLZ :"<< endl;
    cin >> ncellx >> ncelly >> ncellz;
    if( ncellx < 1 ) ncellx = 1;
    if( ncelly < 1 ) ncelly = 1;
    if( ncellz < 1 ) ncellz = 1;

/*  read in parameters for Monte Carlo displacements
      to simulate thermal motion 
    get the initial seed from the time counter if it exists
    otherwise ask user for a seed 
   - this should get a different random number sequence each time the
      program is run 
*/
    lwobble = askYN( "Do you want to add thermal displacements to atomic coord.?" );
    if( lwobble == 1 ) {
        cout << "Temperature in degrees K:" << endl;
        cin >> temperature;
        ltime = (long) time( NULL );
        iseed = (unsigned) ltime;
        if( ltime == -1 ) {
            cout << "Type initial seed for random number generator:"<< endl;
            cin >> iseed;
        } else
            cout <<  "Random number seed initialized to " << iseed << endl;
    }

/*  start timer */

    runtime = cputim();

/*  Start to read in specimen parameters
    remember to use fgets() in ReadLine() and NOT fscanf() 
    because fscanf() ignores newlines and 
        we need to sync on whole lines
*/

    fp.open( filein.c_str() );
    if( fp.bad() ) {
        cout << " can't open crystal data input file" << endl;
        exit(0 );
    }

    ReadLine( fp, isbuf, cline );
    isbuf >> ax >> by >> cz;

    cout << "2D lattice constants= " << ax << " x " << by << " Angstroms\n"
        " and propagation constant= " << cz << " Angstroms\n" << endl ;
    if( (ncellx > 1) || (ncelly > 1) || (ncellz > 1 ) ) {
        ax = ax * ncellx;
        by = by * ncelly;
        cz = cz * ncellz;
        cout << "Unit cell replicated to a= " << ax << ", b= " << by 
            << ", c= " << cz << "  Angstroms" << endl;
    }

    /* read in symmetry operations */

    ReadLine( fp, isbuf, cline );
    isbuf >> nsym;

    if( nsym > 0 ) {
        symx1.resize( nsym );
        symx2.resize( nsym );
        symy1.resize( nsym );
        symy2.resize( nsym );
        for( i=0; i<nsym; i++)  {
            ReadLine( fp, isbuf, cline );
            isbuf >> symx1[i] >> symx2[i] >> symy1[i] >> symy2[i];
       }
    }

/*  Calculate misc constants (Nx*Ny added because FFT2D
    performs scaling) also adjust k2max for circular symmetry
    remember that k2max is still in units of Angstroms
*/
    k2max = 2.0*ax/nx;
    scale = 2.0*by/ny;
    if( scale > k2max ) k2max =scale;

    cout << "Maximum symmetrical resolution set to " << k2max << " Angstroms" << endl;
    k2max = 1.0 / (k2max * k2max);

    rx2 = (1.0/ax);  rx2 = rx2*rx2;
    ry2 = (1.0/by);  ry2 = ry2*ry2;
    scale = ( ((double)nx) * ((double)ny) ) /(ax*by);
    ixmid = nx/2;
    iymid = ny/2;

    ky.resize( ny );
    kx.resize( nx );

    /*  will only need half of these but doesn't take much CPU */
    for( iy=0; iy<ny; iy++) {
        if( iy < iymid )  ky[iy] = (float) iy;
        else ky[iy] = (float)(iy-ny);
    }

    for( ix=0; ix<nx; ix++) {
        if( ix < ixmid )  kx[ix] = (float) ix;
        else kx[ix] = (float)(ix-nx);
    }

    //  allocate arrays
    cpix.resize( nx, ny/2+1 );

    cpix.init( 2 );   //  complex to real transform estimate mode - fast init, slow execution

    for( ix=0; ix<nx; ix++) 
    for( iy=0; iy<=iymid;iy++)  {
        cpix.re(ix,iy) = 0.0F;  // real 
        cpix.im(ix,iy) = 0.0F;  // imag
    }

    /* ---- alloc misc arrays ----- */
    x.resize( NAMAX );
    y.resize( NAMAX );
    occ.resize( NAMAX );
    wobble.resize( NAMAX );

    /* ------ read in actual coordinates ------ */

    total2 = 0.0;
    cout << " " << endl;

    // -------- start loop of atom types -----------------
    do {
        iz = -1;  // just in case
        ReadLine( fp, isbuf, cline );
        isbuf >> iz ;

        if( (cline.length() > 0) && (iz >= 1) && (iz <= NZMAX) ) {
            j = 0;
            total1 = 0.0;
            while( ReadLine( fp, isbuf, cline ) > 2 ) { 

                isbuf >> occ[j] >> x[j] >> y[j] >> wobble[j];
                total1 = total1 +  occ[j] * (nsym+1) * ncellx * ncelly;
                if( nsym > 0) {
                    if( (j+nsym+1) > NAMAX) {
                        cout << "Too many atoms" << endl;
                        cout << "  Maximum allowed = " << NAMAX << endl;
                        fp.close( );
                        exit(0);

                    }
                    for( jj=0; jj<nsym; jj++) {
                        occ[j+jj+1] = occ[j];
                        wobble[j+jj+1] = wobble[j];
                        x[j+jj+1] = (symx1[jj]*x[j] +symx2[jj] )/ncellx;
                        y[j+jj+1] = (symy1[jj]*y[j] +symy2[jj] )/ncelly;
                    }
                }

                x[j] = x[j] /ncellx;
                y[j] = y[j] /ncelly;
                j = j + nsym + 1;
                if( j > NAMAX ) {
                    cout << "Too many atoms" << endl;
                    cout << "  Maximum allowed = " << NAMAX << endl;
                    fp.close( );
                    exit( 0 );
                }
            }   // end while( ReadLine()>2 */

            natom = j ;
            is = 2*(iz-1);
            cout << setw(8) << total1 << " atoms with Z= " << setw(3) << iz 
              << " (" << symbol[is] << symbol[is+1] << ")" << endl;

            /*   if thermal vibrations are requested then we must expand
                ncellx,y here otherwise factor it inside scamp() 
                (its faster)
                 add random displacements to all positions
           */
            if( lwobble == 1 ) {
                if( (j-1)*ncellx*ncelly > NAMAX) {
                    cout << "Too many atoms\n" << endl;
                    cout << "  maximum allowed = " << NAMAX << endl;
                    exit( 0 );
                }

                if( ncellx > 1) {
                    for( ix=1; ix<ncellx; ix++)
                    for( i=0; i<natom; i++) {
                        x[j] = x[i] + ((float)ix)/((float)ncellx);
                        y[j] = y[i];
                        wobble[j] = wobble[i];
                        occ[j++] = occ[i];
                    }
                    natom = j - 1;
                }

                if( ncelly > 1) {
                    for( iy=1; iy<ncelly; iy++)
                    for( i=0; i<natom; i++) {
                        x[j] = x[i];
                        y[j] = y[i] + ((float)iy)/((float)ncelly);
                        wobble[j] = wobble[i];
                        occ[j++] = occ[i];
                    }
                    natom = j - 1;
                }

                /* scale thermal displacements to 300 deg. K 
                     and normalize to 3D - integrating over the z direction */
                scalet = sqrt( temperature / 300.0 );
                for( i=0; i<natom; i++) {
                    x[i] = x[i] + 
                        (float) (wobble[i] * scalet * rangauss( &iseed ) /ax);
                    y[i] = y[i] + 
                        (float) (wobble[i] * scalet * rangauss( &iseed ) /by);
                }

                nsx = 1;
                nsy = 1;

            } else {
                nsx = ncellx;
                nsy = ncelly;

            }  /* end if( lwobble == TRUE ) */

/*  calculate scattering amplitude in upper half plane (FFTW)
   NOTE zero freg is in the bottom left corner and
     expands into all other corners - not in the center
     this is required for FFT, high freq is in the center

    - add extra complex conjugation to correct for the opposite
    sign convention used by FFTW

*/
            ncoeff = 0;
            for( iy=0; iy<=iymid; iy++) {
                ky2 = ky[iy]*ky[iy] * ry2;
                for( ix=0; ix<nx; ix++) {
                    k2 = kx[ix]*kx[ix] * rx2 + ky2;
                    if( k2 <= k2max) {
                        fe = scale * featom( iz, k2 );
                        scamp( kx[ix], ky[iy], &scampr, &scampi ) ;
                        cpix.re(ix,iy) += (float) (scampr * fe);  // real  === bad here ???
                        cpix.im(ix,iy) += (float) (-scampi * fe); // imag 
                        ncoeff++;
                    }

                }  /* end for(ix */
            } /* end for(iy */

            total2 = total2 + total1;
            //goto More;

        }  /* end top if( iz>= 1) */
        
    /*  end loop over types */
    } while( iz>= 1 );  // end do{  Readline....

    cout << "\n   for a grand total of " << total2 << " atoms" << endl;
    fp.close( );

    cpix.ifft();   /*    inverse fft */

    //  output results and find min and max to echo 

    //  copy to floatTIFF pix to output with 
    myFile.resize( nx, ny );
    myFile.setnpix( 1 );
    rmin = rmax = cpix.rre(0,0);
    sum = 0.0;

    for( iy=0; iy<ny; iy++) {
        for( ix=0; ix<nx; ix++) {
            myFile(ix,iy) = cpix.rre(ix,iy);
            sum += myFile(ix,iy);
            if( myFile(ix,iy) < rmin ) rmin = myFile(ix,iy);
            if( myFile(ix,iy) > rmax ) rmax = myFile(ix,iy);
        }
    }

    myFile.zeroParam();
    myFile.setParam( pRMAX, rmax);
    myFile.setParam( pRMIN,  rmin);
    myFile.setParam( pC, (float) cz);
    myFile.setParam( pRES, (float) (1.0/sqrt(k2max)) );
    myFile.setParam( pDX,  dx = (float) ax / ((float)nx) );
    myFile.setParam( pDY,  dy = (float) by / ((float)ny) );

    imin = imax = 0;
    if( myFile.write( fileot.c_str(), rmin,rmax, imin,imax, dx, dy ) != 1 )
        cout <<  "atompot cannot write an output file." << endl;

    cout << " pix range " << rmin << " to " << rmax <<  endl;
    cout << ncoeff << " fourier coeff. calculated in right half plane"
            << endl;
    cout << "The average real space value was "
           << sum /(((float)nx)*((float)ny)) << endl;
          
    cout << "CPU time (excluding set-up) = " << cputim()-runtime 
        << " sec."  << endl;

    return EXIT_SUCCESS;

}  /* end main() */

/*--------------------- ReadLine() -----------------------*/
/*
    read a full line from a file and 
    return length of line read
        
    fp = input file stream
    ss = strinstream to get result
    cline = buffer string to use

    C++/streams is even more awkward than plain C 
*/
int ReadLine( ifstream &fp, istringstream &ss, string &s )
{
    getline( fp, s );   //  read a full line from the file

    // parsing a blank line fails, so don't try
    if( s.length() <= 0 ) {
        return( 0 );
    }

    ss.clear();         //  start fresh each time, just in case
    ss.str( s );
    ss.seekg( 0 );      // MUST rewind for next read !

    //  there is an error here - abort
    if( ss.fail() != 0 ) {
        cout << "--bad data in ReadLine(), cannot continue...." << endl;
        cout << "--ReadLine: s, length= " << s << ", " << s.length() << endl;
        cout << "--ReadLine: ss.str()= " << ss.str()  << endl;
        cout << "--ReadLine: pos, eof = " << ss.tellg() << ", " << ss.eof() << endl;
        cout << "--ReadLine: bad, fail = " << ss.bad() << ", " << ss.fail() << endl;
        exit( 0 );
    }

    return( int(s.length() ) );

}  //  end ReadLine() 

/*--------------------- scamp() -----------------------*/
/*
  Complex specimen scattering amplitude

  kx,ky = real scattering vectors
  scampr, scampi = real and imag. parts of scattering amplitude
           returned by this routine

  global variables used:
  x[natom], y[natom] = real array with atomic coord.
  occ[natom] = real array with occupations at each x,y
  natom      = integer number of atoms (all with same Z)
  nsx, nsy = integer number of cells to replicate in x,y

*/
void scamp( float kx, float ky, double *scampr, double *scampi )
{
    int i, ixc, iyc;
    double pi2, w1, w;
    double  scalex, scaley;

    double phasexr, phasexi, phaseyr, phaseyi;

    pi2 = 8.0 * atan( 1.0 );
    
    scalex = pi2*kx / ((double)nsx);
    scaley = pi2*ky / ((double)nsy);

    *scampr = *scampi = 0.0;

/*  use a trick to sum nsx + nsy terms
    instead of nsx * nsy
    - cos()/sin() could be done recursive for speed
*/
    for( i=0; i<natom; i++) {

        phasexr = phasexi = 0.0;
        w1 = pi2 * kx * x[i];
        for( ixc=0; ixc<nsx; ixc++) {
            w = w1 + ixc*scalex;
            phasexr += cos(w);
            phasexi += sin(w);
        }

        phaseyr = phaseyi = 0.0;
        w1 = pi2 * ky * y[i];
        for( iyc=0; iyc<nsy; iyc++) {
            w = w1 + iyc*scaley;
            phaseyr += cos(w);
            phaseyi += sin(w);
        }

        *scampr += ( phasexr*phaseyr - phasexi*phaseyi ) * occ[i];
        *scampi += ( phasexr*phaseyi + phasexi*phaseyr ) * occ[i];

    } /* end for(i=0... */

    return;

 }  // end scamp() 
