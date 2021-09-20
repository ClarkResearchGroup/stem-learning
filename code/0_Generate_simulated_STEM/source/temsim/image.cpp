/*      *** image.cpp ***

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

  ANSI C and TIFF version
  this version uses FFTW 3

  see:   www.fftw.org

  on Windows file libfftw3f-3.dll must be in the PATH

  on Linux build as:
  g++ -O -o image image.cpp slicelib.o
                       floatTIFF.cpp cfpix.o -lfftw3f


  rewritten in ANSI-C 29-sep-1995 E. Kirkland
  convert to TIFF file format 20-may-1996 ejk
  remove commas in input format 11-july-1997 ejk
  fixed sign error 23-jan-1998 ejk
  add 2-fold and 3-fold astigmatism to coherent mode 24-jan-1998 ejk
  add "\n" to one printf() format 28-jan-1998 ejk
  changed scaling in diffraction pattern mode 2-feb-1998 ejk
  fixed sqrt() in diff patt scaling 19-feb-1998 ejk
  changed variable PI in tcross() to pi because Linux gcc
      predefines PI  9-feb-1999 ejk
  update memory allocation routines 13-nov-1999 ejk
  change void main() to int main() for better portability
         22-jan-2000 ejk
  change data type of nxl,nyl to long32 for compatibility with new
      tiffsubs.c  17-jul-2007 ejk
  convert to GPL 4-jul-2008 ejk
  get return value of scanf() to remove warnings from gcc 4.4
      and convert to 4 char TAB size formatting 10-apr-2010 ejk
  convert to faster FFTW 10,24-apr-2010, 8-may-2010 ejk
  parameterize extra param[] offset so they can be moved easily
      and not conflict with slicelib.h offsets 21-dec-2010 ejk
  change Cs to Cs3 and add Cs5 on 21-dec-2010 ejk
  fix sign error that caused PC images to be flipped in tcross
       17-jan-2011 ejk
  convert to C++ and floatTIFF 21-mar-2012 ejk
  convert to cfpix/fftw class from raw fftw 14-oct-2012 to 30-oct-2012 ejk
  convert to streams and strings 26-apr-2014 ejk
  convert tcross() extra param offsets to const int 21-jun-2015 ejk
  remove extra factor pi in p[p52] in tcross() 21-jun-2015 ejk
  convert malloc1D() to vector<> 11-jul-2016 ejk
  convert malloc2D() to cfpix and remove invert2D()
        (now in cfpix) 30-jul-2017 ejk

  This file is formatted for a tab size of 4 characters

*/
#include <cstdio>  /* ANSI C libraries */
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <ctime>

#include <string>
#include <iostream>  //  C++ stream IO
#include <vector>

using namespace std;

#include "cfpix.hpp"       /* complex image handler with FFT */
#include "slicelib.hpp"    /* misc. routines for multislice */
#include "floatTIFF.hpp"   /* file I/O routines in TIFF format */

const int MCPIX=   0;   /* coherent image mode */
const int MPCPIX=  1;   /* partial coherent mode */
const int MDIFFP=  2;   /* diffraction mode */

// define more parameter keys not in slicelib.h
// make symbolic offsets to merge with slicelib.h def's 
const int
    pB = 260,
    p51 =236,
    p52 =237,

    p31 =251,
    p32 =252,
    p33 =253,
    p34 =254,
    p35 =255,
    p36 =256,
    p37 =257,
    p38 =258,
    p39 =259;

/*  define subroutines at end of file */
void tcross( cfpix &pixin, cfpix &pixo,
        int nx, int ny, vectorf &p );

int main()
{
    string filein, fileout;
    const char version[] = "30-jul-2017 (ejk)";

    int ix, iy, nx, ny, ixmid, iymid, npix, ns,
        itens, mode, nsum, NPARAM;
    int lcenter=0, lapert=0;

    float kx,ky, kx2,ky2, k2, k2max, v0, wavlen, scale, pixc,
        ax, by, cz, rx, ry, pi, rmin, rmax, aimin, aimax,
        Cs3, Cs5, df,  alpha0, ddf, objlx, objly,
        tr, ti, wr, wi, dfa2, dfa2phi, dfa3, dfa3phi;

    double sum, time, chi, chi1, chi2, chi3, phi, clog;

    vectorf param;

    cfpix cpix, cpix2;     //  complex floating point image
    floatTIFF myFile;

/*  echo version date */

    cout << "image version dated " << version << endl;
    cout << "Copyright (C) 1998-2017 Earl J. Kirkland" << endl;
    cout <<  "This program is provided AS-IS with ABSOLUTELY NO WARRANTY\n "
        << " under the GNU general public license\n" << endl;

    cout <<  "calculate TEM images with defocus, using FFTW\n"<< endl;

    pi = (float) (4.0 * atan( 1.0 ));

/*  get input file name */

    cout << "Name of file with input multislice result:"<< endl;
    cin >> filein;

/*  open input file and check sizes etc. */

    if( myFile.tFloatTest( filein.c_str() ) != 1 ) {
        cout << "Cannot open input file " << filein << endl;
        exit( 0 );
    }

/*  ask mode */

    cout << "Type " << MCPIX << " for coherent real space image,\n"
        "  or " << MPCPIX  << " for partially coherent real space image,\n"
        "  or " << MDIFFP  << " for diffraction pattern output:\n" << endl;
    cin >> mode;

/*  get imaging parameters if required */

    Cs3 = 0.0F;
    Cs5 = 0.0F;
    df = 0.0F;
    k2max = 0.0F;
    objlx = 0.0F;
    objly = 0.0F;

    if( (mode == MCPIX) || (mode == MPCPIX) ) {
        cout << "Name of file to get defocused output:" << endl;
        cin >> fileout;

        cout << "Spherical aberration Cs3, Cs5 (in mm.):" << endl;
        cin >>  Cs3 >> Cs5;
        Cs3 *= 1.0e7F;   /*  convert to Angstroms */
        Cs5 *= 1.0e7F;

        cout << "Defocus in Angstroms:" << endl;
        cin >> df;

        cout << "Objective aperture size in mrad:" << endl;
        cin >> k2max;

        if( mode == MPCPIX ) {
            cout << "Illumination semi-angle in mrad:" << endl;
            cin >> alpha0;
            alpha0 = alpha0 * 0.001F;
            cout << "Defocus spread in Angstroms:" << endl;
            cin >> ddf;
        } else if( mode == MCPIX ) {
            cout <<  "Magnitude and angle of two-fold astigmatism"
                << " (in Angst. and degrees):" << endl;
            cin >>  dfa2 >> dfa2phi;
            dfa2phi = dfa2phi * pi /180.0F;
            cout <<  "Magnitude and angle of three-fold astigmatism"
                << " (in Angst. and degrees):" << endl;
            cin >> dfa3 >> dfa3phi;
            dfa3phi = dfa3phi * pi /180.0F;
            cout << "Objective lens and aperture"
                 << " center x,y in mrad\n"
                 << " (i.e. non-zero for dark field):" << endl;
            cin >> objlx >> objly;
        }

/*  get diffraction pattern parameters if required  */

    } else {
        cout << "Name of file to get diffraction pattern:" << endl;
        cin >>  fileout;

        lcenter = askYN("Do you want to include central beam");
        lapert = askYN("Do you want to impose the aperture");
        if ( lapert == 1 ) {
            cout << "Aperture size in mrad:" << endl;
            cin >> k2max;
            cout << "Objective lens and aperture"
                 << " center x,y in mrad\n"
                 << "  (i.e. non-zero for dark field):" << endl;
            cin >> objlx >> objly;
        }

        cout <<  "Type 0 for linear scale,\n"
             << "  or 1 to do logarithmic intensity scale:\n"
             << "  or 2 to do log(1+c*pixel) scale:" << endl;
        cin >> itens;
        if( itens == 2 ) {
            cout <<  "Type scaling constant c:" << endl;
            cin >>  clog;
        }
    }  /* end else... */

    /* ---- read in specimen parameters and multislice result ---- */

    time = cputim();    /* get CPU time for comparison */

    NPARAM = myFile.maxParam();
    param.resize( NPARAM );
    for( ix=0; ix<NPARAM; ix++) param[ix] = 0.0F;

    if( myFile.read( filein.c_str() ) != 1 ) {
        cout << "Cannot open input file " << filein  << endl;
        exit( 0 );
    }
    nx = myFile.nx();
    nx = nx/2;      /* should be complex */
    ny = myFile.ny();
    npix = myFile.getnpix();
    ixmid = nx/2;
    iymid = ny/2;

    if( npix != 2 ) {
        cout << "Input file " << filein << " must be complex, can't continue."
              << endl;
        exit( 0 );
    }

    //  crude port of old code using param[]
    for( ix=0; ix<NPARAM; ix++) param[ix] = myFile.getParam(ix);

    ax = param[pDX] * ((float) nx);
    by = param[pDY] * ((float) ny);
    if( param[pDY] <= 0.0F ) by = param[pDX] * ((float) ny);
    cz = param[pC];
    rx  = 1.0F / ax;
    ry  = 1.0F / by;

    v0 = param[pENERGY];
    wavlen = (float) wavelength( v0 );
    cout << "Starting pix energy = " << v0 << " keV"  << endl;

    rmax  = param[pRMAX];
    aimax = param[pIMAX];
    rmin  = param[pRMIN];
    aimin = param[pIMIN];
    cout << "Starting pix range " << rmin << "  " << rmax << " real\n"
            "                   " << aimin << "  " << aimax << " imag." << endl;

    param[pDEFOCUS] = df;
    param[pASTIG] = 0.0F;
    param[pTHETA] = 0.0F;
    param[pOAPERT] = k2max * 0.001F;
    param[pCS] = Cs3;
    param[pCS5] = Cs5;
    param[pWAVEL] = wavlen;

    k2max = k2max*0.001F/wavlen;
    k2max = k2max * k2max;

    objlx = objlx * 0.001F / wavlen;   /* convert to spatial freq. */
    objly = objly * 0.001F / wavlen;

    cpix.resize( nx, ny );
    cpix.init( 1 );   //  fast init but slow execution is better here

    //------  copy to FFTW style array -------
    for( ix=0; ix<nx; ix++)
        for( iy=0; iy<ny; iy++) {
            cpix.re(ix,iy) = myFile(ix,iy);    // real
            cpix.im(ix,iy) = myFile(ix+nx,iy); // imag
    }

    cpix.fft();   //  every option uses FT so do it here

/*------  coherent real space image mode ---------------------

    Convolute multislice result with objective lens
    aberration function (coherent case)
    with shifted objective lens for dark field

   NOTE zero freq. is in the bottom left corner and
     expands into all other corners - not in the center
     this is required for FFT
*/

    if( mode == MCPIX ) {
        cout << "calculate coherent image...." << endl;
        chi1  = pi * wavlen;
        chi2  = 0.5 * Cs3 * wavlen * wavlen;
        chi3  = Cs5 * wavlen*wavlen*wavlen*wavlen /3.0;

        ns = 0;
        for( ix=0; ix<nx; ix++) {
            kx = (float) ix;
            if( ix > ixmid ) kx = (float) (ix-nx);
            kx = kx*rx - objlx;
            kx2 = kx*kx;
            for( iy=0; iy<ny; iy++) {
                ky = (float) iy;
                if( iy > iymid ) ky = (float)(iy-ny);
                ky = ky*ry - objly;
                k2 = ky*ky + kx2;
                if ( k2 < k2max ) {
                    ns++;
                    phi = atan2( ky, kx );
                    chi = chi1*k2* ( (chi2 + chi3*k2)*k2 - df 
                        + dfa2*sin( 2.0*(phi-dfa2phi) ) 
                        + 2.0F*dfa3*wavlen*sqrt(k2)*
                        sin( 3.0*(phi-dfa3phi) )/3.0 );
                    tr = (float)  cos( chi );
                    ti = (float) -sin( chi );
                    wr = cpix.re(ix,iy);
                    wi = cpix.im(ix,iy);
                    cpix.re(ix,iy) = tr*wr - ti*wi;  /* real */
                    cpix.im(ix,iy) = tr*wi + ti*wr;  /* imag */
                 } else {
                    cpix.re(ix,iy) = 0.0F;  /* real */
                    cpix.im(ix,iy) = 0.0F;  /* imag */
                }
            }  /* end for( ix.... */
        }  /* end for(iy... */

        cpix.ifft();
        cout << "There were " << ns << " pixels inside the obj. apert." << endl;

        //  all options leave results in cpix

/* -------- calculate partially coherent image here  ----------------
    NOTE this assumes that the parameter offsets in txcoef()
    match those in slicelib.h !!!
*/
    } else if( mode == MPCPIX ) {
        cout << "calculate partially coherent image...." << endl;

        param[pCAPERT] = alpha0;
        param[pDDF] = ddf;
        param[ pB ] = 0.0F;

        cpix2.resize( nx, ny );
        cpix2.init( 1 );

        cpix.invert2D();
        tcross( cpix, cpix2, nx, ny, param );
        cpix2.invert2D();

        cpix2.ifft();

        //  all options leave results in cpix
        cpix = cpix2;

    //------------  do diffraction pattern here  ------------------

    } else {

         cout << "calculate diffraction pattern...." << endl;

        //  all options leave results in cpix
        cpix.invert2D();

        if( lcenter == 0 ) 
            cpix.re(ixmid,iymid) = cpix.im(ixmid,iymid) = 0.0F;

        if ( lapert ) {
            for( iy=0; iy<ny; iy++) {
                ky = (iy-iymid)*ry - objly;
                ky2 = ky*ky;
                for( ix=0; ix<nx; ix++) {
                    kx = (ix-ixmid)*rx - objlx;
                    k2 = kx*kx + ky2;
                    if ( k2 > k2max )
                        cpix.re(ix,iy) = cpix.im(ix,iy) = 0.0F;
                } /* end for ix */
            } /* end for iy... */
        }  /* end if( lapert )... */

    }  /* end else... = last mode */

/*  Output results and find min and max to echo */

    cout << "output image" << endl;
    scale = 1.0F / ( ((float)nx) * ((float)ny) );
    
    myFile.resize( nx, ny );
    myFile.setnpix( 1 );
    sum = 0.0;
    nsum = 0;
    for( iy=0; iy<ny; iy++) {
        for( ix=0; ix<nx; ix++) {
            tr = cpix.re(ix,iy);
            ti = cpix.im(ix,iy);
            if( mode == MPCPIX ) pixc = tr;
            else  pixc = tr*tr + ti*ti;
            if ( (mode == MDIFFP) &&  (itens == 1)) {
                if( pixc > 1.e-30F)  pixc = (float) log( (double) pixc );
                else pixc = -30.0F;
            } else if ( (mode == MDIFFP) && (itens == 2)) {
                pixc = (float) log( 1.0 + clog*sqrt((double)pixc) );
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
            myFile(ix,iy) = pixc;   //pixr[ix][iy] = pixc;
        }  /* end for ix... */
    } /* end for iy... */

    param[pRMAX] = rmax;
    param[pIMAX] = 0.0F;
    param[pRMIN] = rmin;
    param[pIMIN] = 0.0F;
    if ( mode == MDIFFP ) {
        param[pRMIN] = (float) (0.05*rmin + 0.95*sum/nsum);
        param[pDX] = 1.0F / (nx*param[pDX]);
        param[pDY] = 1.0F / (ny*param[pDY]);
    }

    //  crude port of old code using param[]
    for(ix=0; ix<NPARAM; ix++) myFile.setParam(ix, param[ix] );

    myFile.write( fileout.c_str(), rmin, rmax, 0.0F, 0.0F, param[pDX], param[pDY] );

    cout << "Pix range " << param[pRMIN] << " to " << rmax << endl;

    time = cputim() - time;
    cout << "Elapsed time = " << time << " sec." << endl;

    return EXIT_SUCCESS;

} /* end main() */

/*------------------- tcross() ------------------------------*/
/*  this version uses cfpix style arrays

 Subroutines to do exact nonlinear partial coherence transfer as in
      [1]  M.A. O'Keefe, 37th EMSA (1979) p.556
      [2]  K. Ishizuka, Ultramicroscopy 5 (1980) p.55-65

  both input and output are in Fourier space with zero frequency
  in the center (NX/2 , NY/2)

 perform 2D weighted convolution with transmission cross coefficient
  TXCOEF to calculate partially coherent CTEM images
    NOTE: this version uses Friedels law

  pixin  : complex input pix array
             range used = P(17) = aperture
  pixo   : complex output pix array
             range output= 2.0*P(17) = 2 x aperture
  nx, ny : actual image size to work with
  p[43]  : parameter array ( dx,dy, defocus etc.) explained in TXCOEF

  started 14-NOV-85 on CIPRES2 Earl J. Kirkland
      started by modifying XFERFN.FTN (from RSX CIPRES)
  changed TXCOEF to be more general (i.e. include alignment and leading
      factors so that it is complete - little or no speed lose)
     28-NOV-85 EJK
  fixed typo (DT3 to D3T) and sign (-D4) in TXCOEF 17-jun-1986 ejk
  changed sign convention to be consistent with forward propagation
    of electrons (see Self et.al. UM 11 (1983) p.35-52 
      - only TXCOEF effected   EJK 2-JUL-86
  converted to C 8-nov-1995 ejk
  fixed sign error 23-jan-1998 ejk
  convert to FFTW style arrays 8-may-2010 ejk
  fix order of psi/psi* indexes to fix flipped image error
      17-jan-2011 ejk
  convert to cfpix style arrays 3-nov-2012 ejk
  remove extra factor pi in p[p52] in tcross() 21-jun-2015 ejk

*/

void tcross( cfpix &pixin, cfpix &pixo,
        int nx, int ny, vectorf &p )
{
    int j, ix, iy, ixmid, iymid, ixmin, ixmax, iymin, iymax,
        ix2min, ix2max, iy0, iy2min, iy2max, ix1, ix2, iy1, iy2,
        ixi0, ixf0;
    float scale, k2p, k2pp, rx, ry, txcoefr, txcoefi,
        xr, xi, yr, yi, zr, zi, PI;

    void txcoef( float kxp, float kyp, float kxpp, float kypp,
         vectorf &p, float *txcoefr, float *txcoefi );

    vectori ixi, ixf;
    vectorf kxp, kyp, kxp2, kyp2;

/*  get scratch arrays and init params */

    ixi.resize( nx );
    ixf.resize( nx );
    kxp.resize( nx );
    kyp.resize( ny );
    kxp2.resize( nx );
    kyp2.resize( ny );

    PI = (float) (4.0 * atan( 1.0 ));
    xr = p[pOAPERT]/p[pWAVEL];
    p[p31] = xr*xr;
    xi = p[pWAVEL];
    p[p32] = 0.5F*PI* p[pCS]* xi*xi*xi;
    p[p33] = PI* p[pWAVEL] * p[pDEFOCUS];

    p[p51] = PI * p[pCS5]* xi*xi*xi*xi*xi/3.0F;
    p[p52] = p[pCS5]* xi*xi*xi*xi;

    xr = PI*p[pCAPERT]* p[pDDF];
    p[p34] = xr*xr;
    p[p35] = p[pCS]* p[pWAVEL]* p[pWAVEL];
    xr = PI*p[pCAPERT];
    p[p36] = xr*xr;
    xr =  PI* p[pWAVEL]* p[pDDF];
    p[p37] = xr*xr/4.0F;
    xr = PI*PI*p[pCAPERT]*p[pCAPERT]*p[pDDF] ;
    p[p38] = xr*xr;
    xr =  PI*p[pCAPERT]*p[pDDF] ;
    p[p39] = PI*( xr*xr )*p[pWAVEL];

    /* initialize  */

    rx = p[pDX]* nx;
    ry = p[pDY]* ny;
    if( ry <= 0.0F) ry = rx;
    rx = 1.0F/rx;
    ry = 1.0F/ry;

    scale=1.0F/( ((float)nx) * ((float)ny) );

    ixmid = nx/2;
    iymid = ny/2;

    /* find range of convolution */

    j = (int) ( p[pOAPERT] / (rx*p[pWAVEL]) +1.0F );
    ixmin  = ixmid - j;
    ixmax  = ixmid + j;
    ix2min = ixmid - 2*j;
    ix2max = ixmid + 2*j;
    j = (int) ( p[pOAPERT] / (ry*p[pWAVEL]) +1.0F );
    iymin  = iymid - j;
    iymax  = iymid + j;
    iy2min = iymid - 2*j;
    iy2max = iymid + 2*j;
    
    if( ixmin < 0 ) ixmin = 0;
    if( ixmax > (nx-1) ) ixmax = (nx-1);
    if( ix2min < 0 ) ix2min = 0;
    if( ix2max > (nx-1) ) ix2max = (nx-1);
    if( iymin < 0 ) iymin = 0;
    if( iymax > (ny-1) ) iymax = (ny-1);
    if( iy2min < 0 ) iy2min = 0;
    if( iy2max > (ny-1) ) iy2max = (ny-1);

    for ( ix=ix2min; ix<=ix2max; ix++) {
       ix1= ixmin - ix + ixmid;
       if( ixmin > ix1 ) ixi[ix] = ixmin;  else  ixi[ix] = ix1;
       ix1= ixmax - ix + ixmid;
       if( ixmax < ix1 ) ixf[ix] = ixmax;   else  ixf[ix] = ix1;
    }

    for( ix=ix2min; ix<=ix2max; ix++) {
        kxp[ix] = (ix-ixmid) * rx;
        kxp2[ix]= kxp[ix] * kxp[ix];
    }
    for( iy=iy2min; iy<=iy2max; iy++) {
        kyp[iy] = (iy-iymid) * ry;
        kyp2[iy]= kyp[iy] * kyp[iy];
    }

    for( ix=0; ix<nx; ix++) for( iy=0; iy<ny; iy++)
        pixo.re(ix,iy) = pixo.im(ix,iy) = 0.0F;

/*  do convolution in bottom half plane 
    the origin is at (ixmid,iymid)
    (ix,iy)   = output point = k
    (ix1,iy1) = input point from psi* = kp
    (ix2,iy2) = input point from psi  = kp + k
*/

    for( iy=iy2min; iy<=iymid; iy++ )
        for( iy1=iymin; iy1<=iymax; iy1++ ) {
            iy2 = iy1 + (iy - iymid);
            if( (iy2>=iymin) && (iy2<=iymax)) {
                for( ix=ix2min; ix<=ix2max; ix++ ) {
                    ixi0 = ixi[ix];
                    ixf0 = ixf[ix];
                    for( ix1=ixi0; ix1<=ixf0; ix1++ ) {
                        ix2= ix1 + (ix - ixmid);
                        k2p = kxp2[ix1] + kyp2[iy1];
                        k2pp= kxp2[ix2] + kyp2[iy2];
                        if ( (k2p<=p[p31]) && (k2pp<=p[p31]) ) {
                            txcoef( kxp[ix1],kyp[iy1],kxp[ix2],kyp[iy2],
                                p, &txcoefr, &txcoefi);
                            xr = pixin.re(ix1,iy1);  /* rpixin[ix1][iy1]; */
                            xi = -pixin.im(ix1,iy1); /* -ipixin[ix1][iy1]; */
                            yr = pixin.re(ix2,iy2);  /* rpixin[ix2][iy2]; */
                            yi = pixin.im(ix2,iy2);  /* ipixin[ix2][iy2]; */
                            zr = xr*yr - xi*yi;
                            zi = xr*yi + xi*yr;
                            pixo.re(ix,iy) += zr*txcoefr - zi*txcoefi;
                            pixo.im(ix,iy) += zr*txcoefi + zi*txcoefr;
                        }  /* end if( ( k2p... */
                    }  /* end for( ix1=... */
                }  /* end for( ix=... */
            }  /* end if( (iy2... */
        }  /* end for iy1... */

/*  scale result */

    for( ix= ix2min; ix<=ix2max; ix++) 
    for( iy= iy2min; iy<=iymid; iy++) {
        pixo.re(ix,iy) *= scale;
        pixo.im(ix,iy) *= scale;
    }

/*  Invoke Friedel's law to get top half plane  */

    iy0 = iymid+1;
    for( iy=iy0; iy<=iy2max; iy++) {
        iy1= iymid - (iy-iymid);
        for( ix=ix2min; ix<=ix2max; ix++) {
            ix1= ixmid - (ix-ixmid);
            pixo.re(ix,iy) =  pixo.re(ix1,iy1); /* rpixo[ix1][iy1]; */
            pixo.im(ix,iy) = -pixo.im(ix1,iy1); /* -ipixo[ix1][iy1]; */
        }
       pixo.re(0,iy) = pixo.im(0,iy) = 0.0F;  /* ??? */
    }

    return;
}     /*   end tcross() */

/*------------------------ txcoef -------------------------*/
/*
   The cross spectral transfer function as in:
    [1] M.A. O'Keefe, 37th EMSA (1979) p.556 EJK
    [2] K. Ishizuka, Ultramic. 5 (1980) p.55-65.

   switched to my derivation 5-DEC-83 EJK
   changed sign convention to be consistent with forward propagation
     of electrons (see Self et.al. UM 11 (1983) p.35-52 
       EJK 2-JUL-86
   converted to C  6-nov-1995 ejk
   fixed sign error 23-jan-1998 ejk
   fix sign to fix flipped image error
      17-jan-2011 ejk

  the following are the array index definitions
  (the order is historical - don't ask why! )

   P(11) = defocus
   P(12) = defocus astigmatism
   P(13) = angle of astigmatism
   P(14) = dx in A pixel dimension
   P(15) = dy in A pixel dimension
   P(17) = aperture in radians
   P(18) = Cs3 in A
   P(19) = wavelength in A
   P(21) = illumination semiangle in rad
   P(22) = defocus spread in A
   P(24) = Debye Waller temp factor/4 in A
   P(31) = P(17)/P(19) **2 maximum k^2 value
 
   P(32) = PI*P(18)* (P(19)**3) /2    ; coherent part
   P(33) = PI* P(19) * P(11)

   p(pCS5) = Cs5 in A         ; add 5th order 21-dec-2010 ejk
   p(p51) = PI*P(pCS5)* (P(19)**5) /3.0 
 
   P(34) = (PI* P(21)* P(22))**2
   P(35) = P(18)* P(19)* P(19)
   P(36) = (PI* P(21)) **2
   P(37) = (PI* P(19)* P(22))**2 /4.
   P(38) = PI^4 * P(21)^4 *P(22)^2
   P(39) = PI^3 * P(22)^2 * P(21)^2 *P(19)
 
   P(p52) = p(pCS5) * p(pWAVEL)^4 
     
   kp  = kprime      = argument of psi*
   kpp = kprime + k  = argument of psi
 
   txcoefr, txcoefi = returned real and imag. parts of function
        (i.e. C can't return a complex value )
*/

void txcoef( float kxp, float kyp, float kxpp, float kypp,
         vectorf &p, float *txcoefr, float *txcoefi )
{
    double kx, ky, k2 ,kp2, kpp2, chip, chipp,
        d1, d2, d3, d3t, d4, d4t, vx, vy, w, u, xr;

    kp2  = kxp*kxp + kyp*kyp;
    kpp2 = kxpp*kxpp + kypp*kypp;

    if( ( kp2 > p[p31] ) || ( kpp2 > p[p31] ) ) {
      *txcoefr = *txcoefi = 0.0F;
      return;
    }

    kx = kxpp - kxp;   /* output freq. k */
    ky = kypp - kyp;
    k2 = kx*kx + ky*ky;

    /*  v = Wc1/(2*pi*lambda)
        w = Wc2/(-pi*lambda)
        u = 1 + pi^2 * beta^2 * delta0^2 k^2
    */
/* -- add Cs5 below 21-dec-2010 ejk */
    vx = p[p52]*( kp2*kp2*kxp - kpp2*kpp2*kxpp )
        + p[p35]*( kp2*kxp - kpp2*kxpp ) + p[pDEFOCUS]*kx;
    vy = p[p52]*( kp2*kp2*kyp - kpp2*kpp2*kypp )
        + p[p35]*( kp2*kyp - kpp2*kypp ) + p[pDEFOCUS]*ky;
    w  = kp2 - kpp2;
    u  = 1.0 + p[p34]*k2;

    d1  = p[p36] * (vx*vx+vy*vy);
    d2  = p[p37] *w*w /u;
    d3t = vx*kx + vy*ky;
    d3  = p[p38]*d3t*d3t/u;
    d4t = p[p39]*w/u;
    d4  = d4t*d3t;

    /* ---- coherent chi  --------------------
         add 5th order Cs5 21-dec-2010 ejk */
    chip  = ( (p[p51]*kp2  + p[p32] )*kp2  - p[p33] ) *kp2;
    chipp = ( (p[p51]*kpp2 + p[p32] )*kpp2 - p[p33] ) *kpp2;

    //chip = chip - chipp - d4;  fix sign 21-jun-2015 ejk
    chip = chip - chipp + d4;
    xr = exp( -d1 -d2 +d3 -p[pB]*(kp2+kpp2) ) / sqrt(u);
    *txcoefr = (float) ( cos(chip) * xr );
    *txcoefi = (float) ( sin(chip) * xr );

    return;

}  /* end txcoef() */


