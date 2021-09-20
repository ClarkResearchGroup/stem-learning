/*      *** probe.cpp ***

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

    ANSI-C++ version

    Calculate a focused probe wavefunction in real space

    this file is formatted for a tab size of 4 characters

    rewritten in C 6-dec-1995 ejk
    fixed sign error in aberration function 1-mar-1997 ejk
    removed commas from keyboard input 3-oct-1997 ejk
    updated memory allocation routines 20-nov-1999 ejk
    change void main() to int main() for better portability
         22-jan-2000 ejk
    add Cs5=5th order aberration and astigmatism  19-jul-2005 ejk
    small cosmetic changes 18-jul-2007 ejk
    convert to GPL 3-jul-2008 ejk
        convert to large list aber. with coma option 23-nov-2008 ejk
    get return value of scanf() to remove warnings from gcc 4.4
      and convert to 4 char TAB size formatting 10-apr-2010 ejk
    convert to FFTW 9-may-2010 ejk
    fix C34a,b terms 10-may-2010 ejk
    change astig. parameters to a,b from mag.+angle 30-may-2010 ejk
    add more aberations to 5th order 30-jun-2010 to 4-jul-2010 ejk
    add probe size calculation 5-jul-2010 ejk
    split up into subroutine 16-jul-2010 ejk
    fix a few things in prbSize() 5-sep-2010 ejk
    switch to storing multipole aberr. in param[] to save in
       output file 2-may-2011 ejk
    add multiMode in chi() to avoid extra calculations if there
       are no multipole aberrations 8-may-2011 ejk
    convert to C++ and floatTIFF.cpp  22-mar-2012 ejk
    convert to cfpix/fftw class from raw fftw 5-nov-2012 ejk
    fix small bug in aspect ration of display image (actual
       floating point image was fine) and add probe position
       parameter symbolic constants  6-apr-2013 ejk
   start conversion to separate class file 1-jul-2013 ejk
   move invert2D() to here (was in incostem) 7-jul-2013 ejk
   move makeProbeIntensity() from mcprobe.cpp to here 24-jan-2016 ejk
   convert malloc1D() etc to vector<> 28-jun-2016 ejk
   convert last malloc2D to cfpix 29-jul-2017 ejk
   add non-zero probe position to makeProbeIntensity() 31-jan-2018 ejk
   remove probe size calc. from makeProbeIntensity() because it
      only works if probe near center - do externally when appropriate
      - change to return = int   2,3-feb-2018 ejk
   copy abbPhase2D( ) from autoslic to here (easer to use) 25-aug-2019 ejk
*/

#include "probe.hpp"   //  header for this class

//=============================================================
//---------------  creator and destructor --------------

probe::probe()
{
}

probe::~probe()
{
}

//=============================================================
/*  abbPhase2D()

  calculate the phase of the abberation function in 2D
  in the objective aperture plane 
  (but out to the max angle allowed by sampling)
  - mainly just to look at

  added 21-aug-2014 ejk
  convert max angle to obj. apert. (instead of 2/3 max)
      20-oct-2017 ejk
  copy from autoslic to here (and add var pi, twopi) 25-aug-2019 ejk

  ab2D() = will get the 2D image of the abb. phase
  param[] = holds image parameters
  multiMode = flag, if not 0 then include all multipole abberations

*/
void probe::abbPhase2D( cfpix &ab2D,vectorf &param, int multiMode )
{
    int ix, iy, ixmid, iymid, nx, ny;
    float k2, k2max, v0, wavlen, ax, by, t, aobj;
    double chi0, alx, aly, pi, twopi;

    vectorf kx, ky, xpos, ypos, kx2, ky2;

    // ---- get setup parameters from param[]
    ax = param[ pAX ];
    by = param[ pBY ];
    nx = ToInt( param[ pNX ] );
    ny = ToInt( param[ pNY ] );
    v0 = param[pENERGY];                // electron beam energy in keV
    aobj = param[ pOAPERT ];

    wavlen = (float) wavelength( v0 );
    pi = (4.0 * atan( 1.0 ));
    twopi = 2.0 * pi;

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
    ab2D = (float) pi;    // white background

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

/*   ------------------ makeProbe() ----------------------
    calculate aberration limited probe wavefunction in real space
    normalized to an integrated intensity of unity
  input:
    cpix[] = fftw array to get wave function
    nx, ny = size of wavefunciton in pixels
    xp, yp = probe position in Ang.
    p[]    = aberration values
    k2max  = max k^2
    pixel  = smoothing range
    multiMode = if not 0 then include multipole aberrations
    ismoth = flag to control smoothing at the edge

  returned value = number of pixels in the aperture

  add arguments  kx[], kx2[], ky[], ky2[] on 14-apr-2013 ejk
*/

int probe::makeProbe( cfpix &cpix, int nx, int ny, double xp, double yp,
    vectorf &p, double wavlen, double k2max, double pixel, int multiMode,
    int ismoth, vectorf &kx, vectorf &kx2, vectorf &ky, vectorf &ky2 )
{ 
    int ix, iy, npixels;
    float scale;
    double alx, aly, k2, chi0, pi, dx2p, dy2p, sum, tr, ti;

    /*    PIXEL = diagonal width of pixel squared
        if a pixel is on the aperture boundary give it a weight
        of 1/2 otherwise 1 or 0
    */
    npixels = 0;
    pi = 4.0 * atan( 1.0 );

    dx2p = 2.0*pi*xp;
    dy2p = 2.0*pi*yp;

    for( iy=0; iy<ny; iy++) {
        aly = wavlen * ky[iy];  /* y component of angle alpha */
        for( ix=0; ix<nx; ix++) {
            k2 = kx2[ix] + ky2[iy];
            alx = wavlen * kx[ix];  /* x component of angle alpha */
            if ( ( ismoth != 0) && 
                ( fabs(k2-k2max) <= pixel) ) {
                   chi0 = (2.0*pi/wavlen) * chi( p, alx, aly, multiMode )
                           - ( (dx2p*kx[ix]) + (dy2p*ky[iy]) );
                cpix.re(ix,iy) = (float) ( 0.5 * cos(chi0));  /* real */
                cpix.im(ix,iy) = (float) (-0.5 * sin(chi0));  /* imag */
                /* printf("smooth by 0.5 at ix=%d, iy=%d\n", ix, iy ); */
            } else if ( k2 <= k2max ) {
                   chi0 = (2.0*pi/wavlen) * chi( p, alx, aly, multiMode )
                           - ( (dx2p*kx[ix]) + (dy2p*ky[iy]) );
                cpix.re(ix,iy) = (float)  cos(chi0);  /* real */
                cpix.im(ix,iy) = (float) -sin(chi0);  /* imag */
                npixels++;
            } else {
                cpix.re(ix,iy) = cpix.im(ix,iy) = 0.0F;
            }
        }  /* end for( ix=0... )  */
    }  /* end for( iy=0... )  */

    cpix.ifft();  //  inverse transform back to real space

    // -----  normalize probe intensity -----
    sum = 0.0;
    for( ix=0; ix<nx; ix++) for( iy=0; iy<ny; iy++) {
            tr = cpix.re(ix,iy); 
            ti = cpix.im(ix,iy);
            sum += tr*tr + ti*ti;
    }
        
    // ----- Normalize probe intensity to unity and add weight ------------ 
    scale = (float)sqrt( 1.0 / sum ); 
    cpix *= scale;

    return( npixels );

}  // end probe::makeProbe()  


/*---------------------------- makeProbeIntensity() ----------------------------*/
/*
   calculate the probe intensity
   modified from incostem (no specimen)

   includes chromatic aberration defocus spread and source size

    started 13-apr-2013 ejk
    change name to just calculate() 21-apr-2013 ejk
    change to calculate2D() 25-jul-2013 ejk
    fix bug in reusing xpos[] if ny>nx 4-dec-2013 ejk
    convert to makeProbeIntensity() 23-jun-2015 ejk
    move to probe.cpp to be more generally available 24-jan-2016 ejk
    add optional probe size return 28-jan-2016 ejk
    add non-zero probe position to makeProbeIntensity() 31-jan-2018 ejk
    remove prb size calculation and return=void 3-feb-2018 ejk

    pix = 2D image to get results
    param[]  = array of probe aberrations and parameters
    multiMode = if 1 include multipole aberrations (0 to ignore)
    initFFT = if 1 then init the FFT otherwise copy from pix
    return < 0 for failure

*/
int probe::makeProbeIntensity( cfpix &pix, vectorf &param, int multiMode,
    int initFFT )
{
    int ix, iy, nx, ny, ixmid, iymid, idf, ndf, ismoth, npixels;

    float ax, by, tr, ti, scale, weight, k2;

    vectorf kx, ky, xpos, kx2, ky2;

    double sum, df, df0, ddf, ddf2, aobj, k2max, dx, dy, keV,
            xp, yp, pixel, pi, dsource, ds, ds2;

    double sigmae, wavl;

    /*  absiccas and weights for Gauss-Hermite Quadrature 
       with exp(-x*x) weighting of integrand 
       from Abramowitz and Stegun, and Numerical Recipes */
    const int NGH=9;   /* number of Gauss-Hermete coeff. to use  */
    double xGH[]={ 3.190993201781528, 2.266580584531843, 1.468553289216668,
        0.723551018752838, 0.000000000000000, -0.723551018752838,
        -1.468553289216668,-2.266580584531843,-3.190993201781528};
    double wGH[]={3.960697726326e-005, 4.943624275537e-003 ,8.847452739438e-002,
        4.326515590026e-001, 7.202352156061e-001, 4.326515590026e-001,
        8.847452739438e-002, 4.943624275537e-003, 3.960697726326e-005};

    cfpix trans, pixsq, ctf;
    std::string sbuff;

    // ---- get setup params from param[]
    ax = param[ pAX ];
    by = param[ pBY ];
    nx = ToInt( param[ pNX ] );
    ny = ToInt( param[ pNY ] );

    if( (nx < 1) || (ny < 1) || (ax<0.0) || (by<0.0) ){
        // gcc requires this step
        sbuff = "bad parameters (nx,ny,ax,by) in probe::makeProbeIntensity()";
        messagePR( sbuff );
        return( -1 );
    }

    df0 = param[pDEFOCUS];
    ddf = param[pDDF];
    wavl = param[pWAVEL];
    aobj = param[pOAPERT];
    keV = param[pENERGY];
    dsource = param[pSOURCE];
    xp = param[pPPOSX];  // probe position 
    yp = param[pPPOSY];

    pi = (float) (4.0 * atan( 1.0 ));

    wavl = wavelength( keV );
    sigmae = sigma( keV )/ 1000.0;

   // ------  calculate spatial frequencies and positions for future use ------ 

    ixmid = nx/2;
    iymid = ny/2;

    //  xpos[] and ypos[]
    kx.resize( nx );
    kx2.resize( nx );
    xpos.resize( nx );
    freqn( kx, kx2, xpos, nx, ax );

    //  reuse xpos[] if possible - its not used anywhere else in the function
    //  but remember to make it bigger if ny>nx
    if( ny > nx ) xpos.resize( ny );

    ky.resize( ny );
    ky2.resize( ny );
    freqn( ky, ky2, xpos, ny, by );

    /*------  make FFTW arrays and plans ------- */
    if( initFFT != 0 ) {
        trans.resize( nx, ny );
        trans.init( 1 );
        pix.resize( nx, ny );
        pix.copyInit( trans );
    } else {
        trans.resize( nx, ny );
        trans.copyInit( pix );
    }

    /* ------  Calculate the probe ------  */
    ctf.resize( nx, ny );
    pixsq.resize( nx, ny );

    k2max = aobj / wavl;
    k2max = k2max * k2max;
    dx = ( ax/((float)nx) );
    dy = ( by/((float)ny) );

    if( ddf > 1.0 ) {       /*  df integration parameters */
        ndf = NGH;
        /* ddf2 = sqrt(log(2.0)/(ddf*ddf));   convert from HWHM */
        ddf2 = sqrt(log(2.0)/(ddf*ddf/4.0));  /* convert from FWHM */
    } else {
        ndf = 1;
        ddf2 = 0.0;
    }

    ctf = 0.0;

    pixel = 0.0;  /* smoothing size = not used */
    ismoth = 0;

    /*---- integrate over defocus spread if ddf is large enough ----*/
    /* use Gauss-Hermite quadrature and convert exp(-a^2df^2) to exp(-u^2)  */
    for( idf=0; idf<ndf; idf++) {

        if( ndf > 1 ){
                df = df0 + xGH[idf]/ddf2;
                weight = (float) wGH[idf];
        }  else {
                df = df0;
                weight = 1.0;
        }
        /*  diagnostic */
        //sbuff = "df step "+ toString(idf)+ ", df= " +toString(df) + 
        //        ", weight= " + toString(weight);
        //messagePR( sbuff, 0 );
            
        param[pDEFOCUS] = (float) df;
            
        /* --------- calculate probe wavefunction -------- */
        npixels = makeProbe( trans, nx, ny, xp, yp, 
            param, wavl, k2max, pixel, multiMode, ismoth,
            kx, kx2, ky, ky2 );

        /* -----  normalize and save probe intensity ----- */
        //  just use real part of pixsq - wastes memory for imag part but is easier
        sum = 0.0;
        for( ix=0; ix<nx; ix++)
        for( iy=0; iy<ny; iy++) {
            tr = trans.re(ix,iy); 
            ti = trans.im(ix,iy);
            pixsq.re(ix,iy) = tr*tr + ti*ti;
            sum +=  pixsq.re(ix,iy);
        }
        
        /* ----- Normalize probe intensity to unity and add weight ------------ */
        scale = weight * ((float)sqrt( 1.0 / sum ) ); 
        
        for( ix=0; ix<nx; ix++) 
           for( iy=0; iy<ny; iy++) {
            ctf.re(ix,iy) +=  scale * pixsq.re(ix,iy);
        }

    }  //  end for(idf...)  

    param[pDEFOCUS] = (float) df0;   // put back the original defocus

    //---- convolve with source size -------
    if( dsource > 0.001 ) {   //  arb. small min. source size (may need to change in future)
        for( ix=0; ix<nx; ix++)
        for( iy=0; iy<ny; iy++) {
            trans.re(ix,iy) = ctf.re(ix,iy);
            trans.im(ix,iy) = 0;
        }
        trans.fft();

        dsource = 0.5* dsource;   // convert diameter to radius
        ds = pi*pi * dsource*dsource/log(2.0);  // source size factor- convert to FWHM
        for( ix=0; ix<nx; ix++)
        for( iy=0; iy<ny; iy++) {
            k2 = kx2[ix] + ky2[iy];
            ds2 = exp( -ds*k2 );
            trans.re(ix,iy) *= (float) ds2;
            trans.im(ix,iy) *= (float) ds2;
        }

        trans.ifft();
        for( ix=0; ix<nx; ix++)
        for( iy=0; iy<ny; iy++) {
             ctf.re(ix,iy) = trans.re(ix,iy);
        }
    }  //  end if( dsource >////

    //----  normalize ------------
    sum = 0.0;
    for( ix=0; ix<nx; ix++)     /* find integrated intensity */
    for( iy=0; iy<ny; iy++)  sum += ctf.re(ix,iy);

    scale = (float) (1.0 / (dx*dy*sum) );
    for( ix=0; ix<nx; ix++)     /* normalize integrated intensity */
    for( iy=0; iy<ny; iy++) {
        pix.re(ix,iy) = scale * ctf.re(ix,iy);
        pix.im(ix,iy) = 0.0F;
    }

    return( 0 );

}  // end probe::makeProbeIntensity()

/* -------------------  messagePR() -------------------
   message output
   direct all output message here to redirect to the command line
   or a GUI status line or message box when appropriate

   msg[] = character string with message to disply
   level = level of seriousness
        0 = simple status message
        1 = significant warning
        2 = possibly fatal error
*/
void probe::messagePR( std::string &smsg,  int level )
{
        messageSL( smsg.c_str(), level );  //  just call slicelib version for now

}  // end probe::messagePR()

/* -------------------  prbsSize() -------------------
   calculate probe size 

   may be wrong if probe is not near the center and wraps around

  input:
     pixsq = intensity in probe (in real part)
     nx, ny = size of probe in pixels
     xp, yp = original probe position (in Ang.)
     ax, by = size of cell in Ang.

  return value:  size of probe in Ang.

*/
double probe::prbSize( cfpix &pixsq, int nx, int ny,
    double xp, double yp, double ax, double by )
{
    int ix, iy, ir, nr;
    double dr, scale, rx, ry, x, y, y2;
    vectord ivsr;

    dr = ax/nx;
    scale = by/ny;
    if( scale < dr ) dr = scale;  //  smallest radial spacing
    dr = 0.5*dr;                  //  use sub pixel sampling
    if( ny > nx ) nr = ny;
    else nr = nx;
    nr = 2*nr;
    ivsr.resize( nr );
    for( ir=0; ir<nr; ir++) ivsr[ir] = 0.0;

    /*  find curent vs. radius = azimuthal average */
    /*    orig. probe position = (dx,dy)  */
    rx = ax/nx;
    ry = by/ny;
    for( iy=0; iy<ny; iy++) {
        y = (iy*ry) - yp;
        y2 = y*y;
        for( ix=0; ix<nx; ix++) {
            scale = pixsq.re(ix,iy);
            x = (ix*rx) - xp;
            ir = (int)( sqrt( x*x + y2 )/dr + 0.5);
            if( ir < 0 ) ir = 0;
            if( ir >= nr ) ir = nr-1;
            ivsr[ir] += scale;
       }
    }

    /*  integrate current */
    for( ir=1; ir<nr; ir++) ivsr[ir] += ivsr[ir-1];

    /*  find half current radius */
    scale = 0.5 * ivsr[nr-1];
    ix = 0;
    for( ir=0; ir<nr; ir++) {
        ix = ir;
        if( ivsr[ir] > scale ) break;
    }

    y = ivsr[ix] - ivsr[ix-1];   /* interpolate */
    if( fabs(y) > 0.0 ) y = (scale-ivsr[ix-1])/y;
    else y = 0.0;
    x = 2.0 * ( (ix-1)*dr + y*dr );

    return( x );

}  //  end probe::prbSize()

