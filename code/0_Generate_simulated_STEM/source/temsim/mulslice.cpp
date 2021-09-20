/*
      *** mulslice.cpp ***

    WARNING: this program has been superceeded and should not be used
    for future application. It is minimally maintainened, not significantly
    updated and may be discontinued at any time (30-oct-2012).

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

  mulslic takes one argument that is the number of threads to use in
  FFTW.  If no argument is prsent then it defaults to NTHREADS_FFTW=4.

  ANSI C++ and TIFF version
  this version uses FFTW 3 (usually 2x to 3x faster)

  FFTW choses an optimum form of the FFT at run time so there
  is some variation in execution speed depending on what else 
  the CPU is doing during this planning stage

  see:   www.fftw.org

  on Windows, file libfftw3f-3.dll must be in the PATH

  on Linux build as:
  g++ -O -fopenmp -o mulslice mulslice.cpp slicelib.o 
                       floatTIFF.o cfpix.o -lfftw3f_threads -lfftw3f

  rewritten in C 26-july-1995 E. Kirkland
  in working form 5-Oct-1995 ejk
  switch to TIFF file format 5-may-1996 ejk
  remove slice expansion (i.e. put in atompot) 4-aug-1996 ejk
  add dx,dy to output parameters 19-feb-1997 ejk
  remove commas from input formats 11-july-1997 ejk
  fixed small problem with anti-aliasing 5-jan-1998 ejk
  added astigmatism in pc mode and inc. beam tilt 28-jan-1998 ejk
  fixed format of error message 16-feb-1998 ejk
  update memory allocation routines 13-nov-1999 ejk
  change void main() to int main() for better portability
         22-jan-2000 ejk
  change data type of nxl,nyl to long32 for compatibility with new
      tiffsubs.c  17-jul-2007 ejk
  convert to GPL 3-jul-2008 ejk
  convert to FFTW 20-mar-2010 ejk
  get return value of scanf() to remove warnings from gcc 4.4
     and convert to 4 char TAB size formatting 21-mar-2010 ejk
  add multithread option for FFTW 28-mar-2010 ejk
  add eleapsed time clock (for multithreading) 16-apr-2010 ejk
  convert to C++ and floatTIFF.cpp  18-mar-2012 to 15-apr-2012 ejk
  convert to cfpix/fftw class from raw fftw 14-oct-2012 to 30-oct-2012 ejk
  convert malloc1D() etc to vector<> 28-jun-2016 ejk
  convert to streams and strings 2,5-aug-2017 ejk
  remove malloc2D() 6-aug-2017 ejk
  fix format of one info output 4-may-2018 ejk
  move propagate() to here so slicelib does not depend on
       cfpix+fftw 30-jul-2019 ejk
  fix small typo's in user questions, and 
       add Cs5 to partial coherence mode 3-aug-2019 ejk

    PIX   = final pix for partial coherence mode
    WAVE  = current specimen transmitted wavefunction
    TRANS = single slice transmission function
    PROPX,PROPY  = propagator function factored as two 1D arrays
    TEMP  = scratch array

  This program calls  subroutines from slicelib.cpp,
    floatTIFF.cpp

  ax,by   = unit cell size in x,y
  BW      = Antialiasing bandwidth limit factor
  acmin   = minimum illumination angle
  acmax   = maximum illumination angle
  Cs3,Cs5 = 3rd, 5th order spherical aberration coefficient
  DF0     = defocus (mean value)
  SIGMAF  = defocus spread (standard deviation)
  DFDELT  = sampling interval for defocus integration
  
  this file is formatted for a TAB size of 4 characters 
  
*/

#include <stdio.h>  /* ANSI C libraries */
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include <string>
#include <iostream>  //  C++ stream IO
#include <sstream>   // string streams
#include <fstream>
#include <iomanip>   //  to format the output
#include <vector>    // STD vector class

using namespace std;

#include "cfpix.hpp"       /* complex image handler with FFT */
#include "slicelib.hpp"    /* misc. routines for multislice */
#include "floatTIFF.hpp"   /* file I/O routines in TIFF format */

const float BW= (2.0F/3.0F);   /* bandwidth limit */
const double ABERR= 1.0e-4;    /* max error for a,b */

const int NSMAX=  1000;  /* max number of slices */
const int NLMAX=    52;  /* maximum number of layers */

/* default number of threads for FFTW to use - >1 thread doesn't work with 32 bit MSVS 2010*/
const int NTHREADS_FFTW=  4;

//  declare subrountines at end
void propagate( cfpix &wave,
    vectorf &propxr, vectorf &propxi, vectorf &propyr, vectorf &propyi,
    vectorf &kx2, vectorf &ky2, float k2max, int nx, int ny );

//------------------------------------------------------------

int main( int argc, char *argv[ ] )
{
    int lstart=0, lpartl=0, lbeams=0;
    int ix,iy, nx,ny, nx2, ny2, ixmid,iymid, i, j, islice, nslice,
        nacx,nacy, iqx, iqy, ilayer, nlayer, npix, nthreads,
        nslic0, nslic1, ndf, idf, nbout, ib, NPARAM, is;
    long nbeams, nillum;
 
    float k2, k2max, scale, v0, vz, mm0, wavlen, rx, ry, dx, dy,
        ax, by, pi, rmin, rmax, aimin, aimax, x,y,
        ax2,by2, rx2,ry2, cztot, ctiltx, ctilty, tctx, tcty,
        acmin, acmax, Cs3, Cs5, df, df0, sigmaf, dfdelt, aobj,
        qx, qy, qy2, q2, q2min, q2max, sumdf, pdf, k2maxo,
        dfa2, dfa2phi, dfa3, dfa3phi, btiltx, btilty;

    float tr, ti, wr, wi;

    string fileout, filestart, filebeam;
    vector<string> filein;

    string ccin, ccin2;
    double sum, timer, etime, xdf, chi, chi1, chi2, chi3, phi, t;

    time_t InitTime, EndTime, pt;

    ofstream fp1;

    vectori layer, hbeam, kbeam;
    vectorf kx, ky, kx2, ky2, xpos, ypos, cz, param, sparam;

    cfpix wave, wave0, ctemp, pix, *trans;     //  complex floating point images
    floatTIFF myFile;

    /*  set up symbolic mapping
        this must be the same as in parlay */
    const string cname = 
        "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";

    /*  find number of threads option */
    if( argc < 2 ) nthreads = NTHREADS_FFTW;
    else  {
        istringstream iss( argv[1] );
        iss >> nthreads;
    }

    // ----- echo version date ----- 

    cout << "mulslice version dated 3-aug-2019 ejk" << endl;
    cout << "Copyright (C) 1998-2019 Earl J. Kirkland"  << endl;
    cout <<  "This program is provided AS-IS with ABSOLUTELY NO WARRANTY\n "
        " under the GNU general public license\n"  << endl;

    cout << "Warning: this program is obsolete and may be discontinued soon\n" << endl;

    cout << "perform traditional CTEM multislice calculation" << endl;
    cout << "  using FFTW with " << nthreads << " thread(s)\n" << endl;

    NPARAM = myFile.maxParam();
    pi = (float) (4.0 * atan( 1.0 ));
    param.resize( NPARAM );
    sparam.resize( NPARAM );

    /*  read in the layer stacking sequence and parse it
        multiple line continuation is signified with a  trailing '-'
        a trailing '/echo' displays the results of parsing
    */

    //-----  don't have multiline mode yet - forget it for now
    cout << "Type in the stacking sequence :" << endl;
    cin >> ccin;

    layer.resize( NSMAX );

    if( parlay( ccin.c_str(), layer, NSMAX, NLMAX, &nslice, 1)
       < 0 ) exit( 0 );

    /*  Find total number of layers  */

    nlayer = 0;
    for( i=0; i<nslice; i++) 
        if( layer[i] > nlayer ) nlayer = layer[i];

    nlayer += 1;

    /*  Get input file name etc. */

    cout << "\nType in the name of " << nlayer << " atomic potential layers :\n"
         << endl;

    filein.resize( nlayer );
    for( i=0; i<nlayer; i++) {
       cout << "Name of file with input atomic potential  :"
             << cname[i] << endl;
       cin >> filein[i];
    }

    /*  get more file names etc. */

    cout << "Name of file to get binary output of multislice result:" << endl;
    cin >> fileout;

    lpartl = askYN("Do you want to include partial coherence");

    if( lpartl == 1 ) {
        cout << "Illumination angle min, max in mrad:"<< endl;
        cin >> acmin >> acmax;
        acmin  = acmin  * 0.001F;
        acmax  = acmax  * 0.001F;
        cout << "Spherical aberration Cs3, Cs5 (in mm.):" << endl;
        cin >> Cs3 >> Cs5;
        Cs3 = Cs3 * 1.0e7F;   // convert mm. to Ang.
        Cs5 = Cs5 * 1.0e7F;
        cout << "Defocus, mean, standard deviation, and"
               " sampling size (in Ang.):"<< endl;
        cin >> df0 >> sigmaf >> dfdelt;
        cout << "Objective aperture (in mrad) =" << endl;
        cin >> aobj;
        aobj = aobj * 0.001F;
        cout <<  "Magnitude and angle of 2-fold astig."
            " (in Ang. and degrees):" << endl;
        cin >> dfa2 >> dfa2phi;
            dfa2phi = dfa2phi * pi /180.0F;
        cout <<  "Magnitude and angle of 3-fold astig."
            " (in Ang. and degrees):" << endl;
        cin >> dfa3 >> dfa3phi;
        dfa3phi = dfa3phi * pi /180.0F;
        lstart = 0;
    } else {
        cout << "NOTE, the program image must also be run." << endl;
        lstart = askYN("Do you want to start from previous result");
    }

    if ( lstart == 1 ) {
        cout << "Name of file to start from:" << endl;
        cin >> filestart;
    } else {
        cout << "Incident beam energy in kev:" << endl;
        cin >> v0;
    }

    cout << "Crystal tilt x,y in mrad.:" << endl;
    cin >> ctiltx  >> ctilty ;
    ctiltx = ctiltx * 0.001F;
    ctilty = ctilty * 0.001F;

    if( lpartl == 0 ) {
        cout << "Incident beam tilt x,y in mrad.:" << endl;
        cin >> btiltx >> btilty;
        btiltx = btiltx * 0.001F;
        btilty = btilty * 0.001F;

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
    timer = cputim();          // get initial CPU time
    InitTime = time( &pt );    // get initial wall time (test multithreading)

/*  get starting value of transmitted wavefunction if required
   (this can only be used in coherent mode)
    remember to save params for final output pix  */

    if ( lstart == 1 ) {
        if( myFile.read( filestart.c_str() ) != 1 ) {
            cout << "Cannot open input file: " << filestart << endl; 
            exit( 0 );
        }

        if( myFile.getnpix() != 2 ) {
           cout << "Input file " << filestart << " must be complex, can't continue."
               << endl;
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
        cout << "Starting pix range " << sparam[pRMIN] << " to " << sparam[pRMAX] << " real\n"
               "           " << sparam[pIMIN] << " to " << sparam[pIMAX] << " imag"
               << endl;
        cout << "Beam voltage = " << v0 << " kV" << endl;
        cout << "Old crystal tilt x,y = " << 1000.*sparam[pXCTILT]
                << ", " << 1000.*sparam[pYCTILT] << " mrad" << endl;

    } else {
      nslic0 = 0;
    }

    /*  calculate relativistic factor and electron wavelength */
    mm0 = 1.0F + v0/511.0F;
    wavlen = (float) wavelength( v0 );
    cout << "Wavelength = " << wavlen << " Angstroms" << endl;

/*  read in atomic potential and specimen parameters
    and calculate specimen transmission function
    for a single slice in transr,i */

    cz.resize( nlayer );
    trans = new cfpix[nlayer]; 
    if( NULL == trans ) {
        cout << "Cannot allocate transmission function memory" << endl;
        exit( 0 );
    }

    for( ilayer=0;  ilayer<nlayer; ilayer++ ) {

        myFile.zeroParam();
        if( (is= myFile.read( filein[ilayer].c_str() )) != 1){
            cout << "Cannot read input file " << filein[ilayer] 
                    << " (status= " << is << ")." << endl;
            exit(0);
        }

        nx2 = myFile.nx();
        ny2 = myFile.ny();
        npix = myFile.getnpix();
        for( i=0; i<NPARAM; i++) param[i] = myFile.getParam( i );

        trans[ilayer].resize(nx2, ny2 );

        if( 0 == ilayer ) trans[ilayer].init( 0, nthreads);
        else trans[ilayer].copyInit( trans[0] );

        if( npix != 1 ) {
            cout << "Input potential file " << filein[ilayer] << " is not real."
                 << endl;
            exit( 0 );
        }
 
        cz[ilayer] = param[pC];
        cout << "layer " << cname[ilayer] << ", cz = " << cz[ilayer] << endl;

        if ( ( lstart==1 ) || ( ilayer != 0) ) {
            if ( (nx!=nx2) || (ny!=ny2) ) {
                cout << "pix size incompatible." << endl;
                cout << "old size = " << nx << ", " << ny << endl;
                cout << "new size = " << nx2 << ", " << ny2 << endl;
                cout << "layer = " << cname[ilayer] << endl;
                exit( 0 );
            }
            ax2 = param[pDX] * nx;
            by2 = param[pDY] * ny;
            if( ( fabs( ax-ax2 ) > fabs(ABERR*ax) ) ||
                ( fabs( by-by2 ) > fabs(ABERR*by) ) ) {
                cout << "incompatible lattice constants" << endl;
                cout << "potential    a,b,c = " << ax2 << ", " << by2 
                    << ", " << cz[ilayer] <<  endl ;
                cout << "starting pix a,b,c = " << ax << ", " << by << endl;
                cout << "   layer = " << cname[ilayer] << endl;
                exit( 0 );
            }
        } else {
            nx = nx2;
            ny = ny2;
            ax = param[pDX] * nx;
            by = param[pDY] * ny;
        }

        scale = wavlen * mm0;
        for( iy=0; iy<ny; iy++) {
            for( ix=0; ix<nx; ix++) {
                vz= scale * myFile(ix, iy);
                trans[ilayer].re(ix,iy) = (float) cos(vz);
                trans[ilayer].im(ix,iy) = (float) sin(vz);
            }
        }

    }  /* end for(ilayer=... */

    cout << "Size in pixels Nx x Ny= " << nx << " x " << ny << " = " 
        << nx*ny << " beams" << endl;
    cout << "Lattice constant a = " << ax << ", b = " << by << endl;

/*  calculate the total specimen thickness and echo */

    cztot = 0.0F;
    for( islice=0; islice<nslice; islice++) 
        cztot += cz[ layer[islice] ];
    cout << "Total specimen thickness = " << cztot << " Angstroms" << endl;
    
/*  calculate spatial frequencies and positions for future use */

    rx = 1.0F/ax;
    rx2= rx*rx;
    ry = 1.0F/by;
    ry2= ry*ry;
    ixmid = nx/2;
    iymid = ny/2;

    kx.resize( nx );
    kx2.resize( nx );
    xpos.resize( nx );
    freqn( kx, kx2, xpos, nx, ax );

    ky.resize( ny );
    ky2.resize( ny );
    ypos.resize( ny );
    freqn( ky, ky2, ypos, ny, by );

/*  allocate some more arrays and initialize wavefunction */

    wave.resize( nx, ny );
    wave.copyInit( trans[0] );

    if( lstart == 0 ) {
        qx = btiltx / wavlen;   // add incident beam tilt
        qy = btilty / wavlen;
        for( ix=0; ix<nx; ix++)
        for( iy=0; iy<ny; iy++) {
            t = 2.0*pi*( qx*xpos[ix] + qy*ypos[iy] );
            wave.re(ix,iy) = (float) cos( t );  // real
            wave.im(ix,iy) = (float) sin( t );  // imag
        }
    } else {    // get incident wavefunctions
        for( ix=0; ix<nx; ix++)
        for( iy=0; iy<ny; iy++) {
            wave.re(ix,iy) = wave0.re(ix,iy);  //  real part
            wave.im(ix,iy) = wave0.im(ix,iy);  //  imag part
        }
    }  // end if( lstart...)

/*  calculate propagator function, and bandwidth limit the transmission
          function for anti-aliasing.
    set to zero outside sampling circle
*/
    k2max = nx/(2.0F*ax);
    tctx = ny/(2.0F*by);
    if( tctx < k2max ) k2max = tctx;
    k2max = BW * k2max;
    cout << "Bandwidth limited to a real space resolution of " << 1.0F/k2max
            << " Angstroms" << endl;
    cout << "   (= " << wavlen*k2max*1000.0F << " mrad)  for symmetrical anti-aliasing." 
        << endl;
    k2max = k2max*k2max;

    tctx = (float) (2.0 * tan(ctiltx));
    tcty = (float) (2.0 * tan(ctilty));

    // make array of propagator arrays
    vector< vector<float> > propxr( nlayer, kx );
    vector< vector<float> > propxi( nlayer, kx );
    vector< vector<float> > propyr( nlayer, ky );
    vector< vector<float> > propyi( nlayer, ky );

    for( ilayer=0; ilayer<nlayer; ilayer++) {

        scale = pi * cz[ilayer];
        nbeams = 0;

        for( ix=0; ix<nx; ix++) {
            t = scale * ( kx2[ix]*wavlen - kx[ix]*tctx );
            propxr[ilayer][ix] = (float)  cos(t);
            propxi[ilayer][ix] = (float) -sin(t);
        }
        for( iy=0; iy<ny; iy++) {
            t = scale * ( ky2[iy]*wavlen - ky[iy]*tcty );
            propyr[ilayer][iy] = (float)  cos(t);
            propyi[ilayer][iy] = (float) -sin(t);
        }

        trans[ilayer].fft();
        for( ix=0; ix<nx; ix++) for( iy=0; iy<ny; iy++) {
            k2= ky2[iy] + kx2[ix];
            if (k2 < k2max) nbeams++;
            else trans[ilayer].re(ix,iy) =
                     trans[ilayer].im(ix,iy) = 0.0F;
        }
        trans[ilayer].ifft();

    } /* end for(ilayer... */

    cout << "Number of symmetrical non-aliasing beams = " << nbeams << endl;

/*  iterate the multislice algorithm proper

   NOTE: zero freq. is in the bottom left corner and
     expande into all other corners - not in the center
     this is required for the FFT - don't waste time rearranging

-------- partial coherence method ----------------------------
   force the integrals to include the origin and to be symmetric
   about the origin and to have the same periodic boundary
   conditions as the sampling grid
*/
    if( lpartl == 1 ) {

        /*  for defocus integration */
        ctemp.resize(nx,ny);
        ctemp.copyInit( trans[0] );

        cout << "Illumination angle sampling (in mrad) = " << 1000.*rx*wavlen << 
            ", " << 1000.*ry*wavlen << "\n" << endl;

        pix.resize( nx, ny );
        pix = 0;
        
        scale = 1.0F / ( ((float)nx) * ((float)ny) );
        ndf = (int) ( ( 2.5F * sigmaf ) / dfdelt );

        nacx = (int) ( ( acmax / ( wavlen * rx ) ) + 1.5F );
        nacy = (int) ( ( acmax / ( wavlen * ry ) ) + 1.5F );

        q2max = acmax / wavlen;
        q2max = q2max*q2max;

        q2min = acmin / wavlen;
        q2min = q2min*q2min;

        k2maxo = aobj / wavlen;
        k2maxo = k2maxo*k2maxo;

        chi1 = pi * wavlen;
        chi2 = 0.5 * Cs3 * wavlen *wavlen;
        chi3 = Cs5 * wavlen*wavlen*wavlen*wavlen /3.0;
        nillum = 0;

        /*  integrate over the illumination angles */

        for( iqy= -nacy; iqy<=nacy; iqy++) {
            qy = iqy * ry;
            qy2 = qy * qy;

            for( iqx= -nacx; iqx<=nacx; iqx++) {
                qx = iqx * rx;
                q2 = qx*qx + qy2;

                if( (q2 <= q2max) && (q2 >= q2min) ) {
                    nillum += 1;
                    for( ix=0; ix<nx; ix++) for( iy=0; iy<ny; iy++) {
                        t = 2.0*pi*( qx*xpos[ix] + qy*ypos[iy] );
                        wave.re(ix,iy) = (float) cos(t);    // real
                        wave.im(ix,iy) = (float) sin(t);    // imag
                    }
                    for( islice=0; islice<nslice; islice++) {
                        ilayer = layer[islice];
                        wave *= trans[ilayer];   //  transmit
                        wave.fft();
                        propagate( wave, 
                            propxr[ilayer], propxi[ilayer],
                            propyr[ilayer], propyi[ilayer],
                            kx2,  ky2,  k2max, nx, ny );
                        wave.ifft();
                    }

                   sum = 0.0;
                   for( ix=0; ix<nx; ix++) for( iy=0; iy<ny; iy++)
                       sum += wave.re(ix,iy)*wave.re(ix,iy)
                            + wave.im(ix,iy)*wave.im(ix,iy);
                   sum = sum * scale;

                   cout << "Illum. angle = " << 1000.*qx*wavlen <<
                       ", " << 1000.*qy*wavlen << " mrad" 
                       << ", total intensity= " << sum << endl;

                   /*  integrate over +/- 2.5 sigma of defocus */

                   wave.fft();
                   sumdf = 0.0F;

                   for( idf= -ndf; idf<=ndf; idf++) {
                       df = df0 + idf*dfdelt;

                       for( ix=0; ix<nx; ix++)
                       for( iy=0; iy<ny; iy++) {
                           k2 = kx2[ix] + ky2[iy];
                           if( k2 <= k2maxo ) {
                               phi = atan2( ky[iy], kx[ix] ); 
                               chi = chi1*k2* ( (chi2 + chi3*k2)*k2 - df 
                                   + dfa2*sin( 2.0*(phi-dfa2phi) ) 
                                   + 2.0F*dfa3*wavlen*sqrt(k2)*
                                     sin( 3.0*(phi-dfa3phi) )/3.0 );
                               tr = (float)  cos(chi);
                               ti = (float) -sin(chi);
                               wr = wave.re(ix,iy);   // real
                               wi = wave.im(ix,iy);   // imag
                               ctemp.re(ix,iy) = wr*tr - wi*ti;
                               ctemp.im(ix,iy) = wr*ti + wi*tr;
                           } else {
                               ctemp.re(ix,iy) = 0.0F;
                               ctemp.im(ix,iy) = 0.0F;
                           }
                       }  /* end for( iy=...) */

                       ctemp.ifft();
 
                       xdf = (double) ( (df - df0) /sigmaf );
                       pdf = (float) exp( -0.5F * xdf*xdf );
                       sumdf += pdf;

                       for( ix=0; ix<nx; ix++) {
                           j = ix*ny;
                           for( iy=0; iy<ny; iy++) {
                               x = ctemp.re(ix,iy);  // real
                               y = ctemp.im(ix,iy);  // imag
                               pix.re(ix,iy) += pdf* ( x*x + y*y );
                           }
                       }

                   }/* end for(idf..) */
                }/* end if( q2...) */

            } /* end for( iqx..) */
        } /* end for( iqy..) */

        cout << "Total number of illumination angle = " << nillum << endl;
        cout << "Total number of defocus values = " << 2*ndf+1 << endl;
        scale = 1.0F / ( ((float)nillum) * sumdf );
        rmin  = pix.re(0,0) * scale;
        rmax  = rmin;
        aimin = 0.0F;
        aimax = 0.0F;

        for( ix=0; ix<nx; ix++)
        for( iy=0; iy<ny; iy++) {
            pix.re(ix,iy) *= scale;
            if( pix.re(ix,iy) < rmin ) rmin = pix.re(ix,iy);
            if( pix.re(ix,iy) > rmax ) rmax = pix.re(ix,iy);
        }

/* ------------------- coherent method -------------------- 
    remember that wave was initialized above */

    } else {
        if( lbeams == 1 ) {
            fp1.open( filebeam.c_str() );
            if( fp1.fail() ) {
                cout << "can't open file " << filebeam << endl;
                exit(0);
            }
            fp1 << " (h,k) = " ;
            for(ib=0; ib<nbout; ib++)
                fp1 << " (" << hbeam[ib] << "," << kbeam[ib] << ")";
            fp1 <<  endl;
            fp1 <<  "nslice, (real,imag) (real,imag) ..." << endl;
            for( ib=0; ib<nbout; ib++) {
                if( hbeam[ib] < 0 ) hbeam[ib] = nx + hbeam[ib];
                if( kbeam[ib] < 0 ) kbeam[ib] = ny + kbeam[ib];
                if( hbeam[ib] < 0 ) hbeam[ib] = 0;
                if( kbeam[ib] < 0 ) kbeam[ib] = 0;
                if( hbeam[ib] > nx-1 ) hbeam[ib] = nx-1;
                if( kbeam[ib] > ny-1 ) kbeam[ib] = ny-1;
            }
        }  /* end if( lbeams..) */
    
        nslic1 = nslic0;
        scale = 1.0F / ( ((float)nx) * ((float)ny) );

        for( islice=0; islice<nslice; islice++ ) {

            ilayer =layer[islice];
            wave *= trans[ilayer];   // transmit
            wave.fft();

            // remember: prop must be here to anti-alias
            propagate( wave,
                propxr[ilayer], propxi[ilayer],
                propyr[ilayer], propyi[ilayer],
                kx2,  ky2,  k2max, nx, ny );
            if( lbeams == 1 )  {
                fp1 << setw(5) <<  nslic0+1;
                for( ib=0; ib<nbout; ib++) 
                    fp1 << setw(10) << scale*wave.re(hbeam[ib],kbeam[ib])
                       << setw(10) <<scale*wave.im(hbeam[ib],kbeam[ib]);
                fp1 <<  endl;
            }
            wave.ifft();
 
            sum = 0.0;
            for( ix=0; ix<nx; ix++) {
                for( iy=0; iy<ny; iy++)
                    sum += wave.re(ix,iy)*wave.re(ix,iy)
                       + wave.im(ix,iy)*wave.im(ix,iy);
            }
            sum = sum * scale;

            nslic0 +=  1;
            cout << "slice " << setw(4) << nslic0 << ", layer = "
                << cname[ilayer] << ", integrated intensity = " << sum << endl;

        } /* end for(islice...) */

    wave.findRange( rmin, rmax, aimin, aimax );

    } /* end else .. coherent section */

/*  output results and find min and max to echo */
    if( lstart == 1 )
        for( ix=0; ix<NPARAM; ix++ ) myFile.setParam( ix, sparam[ix]);
    else
        for( ix=0; ix<NPARAM; ix++ ) myFile.setParam( ix, 0.0F);

    myFile.setParam( pMODE, mMULSLICE);
    myFile.setParam( pRMAX, rmax);
    myFile.setParam( pIMAX, aimax);
    myFile.setParam( pRMIN, rmin);
    myFile.setParam( pIMIN, aimin);
    myFile.setParam( pXCTILT, ctiltx);
    myFile.setParam( pYCTILT, ctilty);
    myFile.setParam( pENERGY, v0);
    myFile.setParam( pDX, dx = (float) ( ax/((float)nx) ) );
    myFile.setParam( pDY, dy = (float) ( by/((float)ny) ) );
    myFile.setParam( pWAVEL, wavlen );
    myFile.setParam( pNSLICES, (float) nslic0 );
    if ( lpartl == 1 ) {
        myFile.setParam( pDEFOCUS, df0 );
        myFile.setParam( pOAPERT, aobj );
        myFile.setParam( pCS, Cs3 );
        myFile.setParam( pCS5, Cs5 );
        myFile.setParam( pCAPERT, acmax );
        myFile.setParam( pDDF, sigmaf );
    }

    if ( lpartl == 1 ) {
        myFile.resize( nx, ny );
        myFile.setnpix( 1 );
        for( ix=0; ix<nx; ix++) for( iy=0; iy<ny; iy++)
            myFile(ix,iy) = pix.re(ix,iy);
        i = myFile.write( fileout.c_str(), rmin, rmax, aimin, aimax, dx, dy );
    } else {
        cout <<  "make output pix " << nx << " x " << ny << endl;
        myFile.resize( 2*nx, ny );
        myFile.setnpix( 2 );
        for( ix=0; ix<nx; ix++) for( iy=0; iy<ny; iy++) {
            myFile(ix,iy) = wave.re(ix,iy);
            myFile(ix+nx,iy) = wave.im(ix,iy);
        }
        i = myFile.write( fileout.c_str(), rmin, rmax, aimin, aimax, dx, dy );
    }
    if( i != 1 ) cout <<  "mulslice cannot write TIF file " << fileout << endl;

    cout <<  "pix range " <<  rmin << " to " << rmax  << " real,\n"
        << "      " <<  aimin << " to " << aimax  << " imag" << endl;

    /* include all CPU's */
    cout << "Total CPU time = " << cputim()-timer  << " sec." << endl;

    EndTime = time( &pt );
    etime = difftime( EndTime, InitTime );
    cout << "elapsed time = " << etime  << " sec." << endl;

    return EXIT_SUCCESS;

} /* end main() */

//
//  propagate should be the same in mulslice.cpp and stemslic.cpp
//
/*------------------------ propagate() ------------------------*/
/*
    propagate the wavefunction thru one layer

    wave  = complex wavefunction
    propxr,i[ix]     = real and imag. parts of x comp of propagator
    propyr,i[iy]     = real and imag. parts of y comp of propagator

    kx2[], ky2[]     = spatial frequency components
    k2max        = square of maximum k value 
    nx, ny       = size of array
    
    on entrance waver,i and 
         propxr,i/propyr,i are in reciprocal space
    
    only wave will be changed by this routine
*/
void propagate( cfpix &wave,
    vectorf &propxr, vectorf &propxi, vectorf &propyr, vectorf &propyi,
    vectorf &kx2, vectorf &ky2, float k2max, int nx, int ny )
{
    int ix, iy;
    float pxr, pxi, pyr, pyi, wr, wi, tr, ti;

    /*  multiplied by the propagator function */

/*  parallizing this loop usually runs slower (!)     
#pragma omp parallel for private(j,iy,pxr,pxi,pyr,pyi,wr,wi,tr,ti) */
    for( ix=0; ix<nx; ix++) {
        if( kx2[ix] < k2max ) {
            pxr = propxr[ix];
            pxi = propxi[ix];
            for( iy=0; iy<ny; iy++) {
                if( (kx2[ix] + ky2[iy]) < k2max ) {
                    pyr = propyr[iy];
                    pyi = propyi[iy];
                    wr = wave.re(ix,iy);   // real 
                    wi = wave.im(ix,iy);   // imag
                    tr = wr*pyr - wi*pyi;
                    ti = wr*pyi + wi*pyr;
                    wave.re(ix,iy) = tr*pxr - ti*pxi;
                    wave.im(ix,iy) = tr*pxi + ti*pxr;
                } else {
                    wave.re(ix,iy) = 0.0F;
                    wave.im(ix,iy) = 0.0F;
                }
            } /* end for(iy..) */

        } else for( iy=0; iy<ny; iy++) {
                    wave.re(ix,iy) = 0.0F;
                    wave.im(ix,iy) = 0.0F;
        }  /*  end if( kx2[ix]... */
            
    } /* end for(ix..) */

} /* end propagate() */
