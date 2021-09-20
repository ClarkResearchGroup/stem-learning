/*      *** slicelib.cpp ***

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

    function library for use with multislice programs
    
    should put each subroutine in a separate file and make
    a real obj library but its easier to manage this way

    askYN()       : ask a yes/no question and return true/false
    bessi0()      : modified Bessel function I0(x)
    bessk0()      : modified Bessel function K0(x)
    chi()         : return aberration function (needs 2pi/wave factor)
    cputim()      : return current CPU time in sec
    featom()      : return scattering factor for a given atom
    freqn()       : calculate spatial frequencies
    messageSL()   : message handler
    parlay()      : parse the layer structure
    propagate()   : propagate a 2D wavefunction
    propagateW()   : propagate a 2D wavefunction
    ranflat()     : return a random number with uniform distribution
    rangauss()    : return a random number with Gaussian distribution
    ranPoisson()  : return a random number with Poisson distribution
    readCnm()     : decypher aberr. of the form C34a etc.
    ReadfeTable() : read fe scattering factor table
    ReadXYZcoord(): read a set of (x,y,z) coordinates from a file
    scaleW()       : scale a 2D FFTW
    seval()       : Interpolate from cubic spline coefficients
    sigma()       : return the interaction parameter
    sortByZ()     : sort atomic x,y,z coord. by z
    splinh()      : fit a quasi-Hermite  cubic spline
    transmit()    : transmit a 2D wavefunction
    transmitW()   : transmit a 2D wavefunction with openMP
    vatom()       : return real space atomic potential (NOT projected)
    vzatom()      : return real space projected atomic potential
    vzatomLUT()   : same as vzatom() but with a look-up-table (faster)
    wavelength()  : return electron wavelength in A for given keV

    this file is formatted for a tab size of 4 characters

    started may 1996 E. Kirkland
    added askYN() 23-jun-1996 ejk
    add scattering factor routines featom, ReadfeTable, ReadLine
        seval, spinh    26-july-1996 ejk
    move random number generator here 3-aug-1996 ejk
    move freqn(), transmitt() and propagate() to here 4-aug-1996 ejk
    converted to new 12 parameter fe(k) from mcdf 16-jan-1997 ejk
    add bessk0(), bessi0(), vzatom() and wavelength() 18-feb-1997 ejk
    update random number generators 20-may-1997 ejk
    added more directories to look for "fparams.dat" 5-july-1997 ejk
    added vatom() and changed constants in vxatom() in about
            5th sig. fig.  24-nov-1997 ejk
    added sigma() 30-nov-1997 ejk
    change min Z to 1 for hydrogen 15-dec-1997 ejk
    fixed possible log(0) problem in rangaus() 10-jan-1998 ejk
    moved seval(), splinh() and ReadXYZcoord()
            into this library 11-jan-1998 ejk
    switched scattering factors to hardcoded table
       so you don't have to worry about an external data
       file (and also so it can't be changed accidentily)
          21-jan-199 ejk
    add nowarranty() 27-jan-1999 ejk
    put in new scattering factor parameters 9-oct-1999 ejk
    move malloc1D() and malloc2D() to here 6-nov-1999 ejk
    rearrange nowarranty 22-jan-2000 ejk
    add malloc3D(), free2D() and free3D()  20-jul-2005 ejk
    move sortbyz() to here, and add short nowaranty 23-aug-2006 ejk
    remove nowarranty() subroutine from this file 11-jul-2007 ejk
    fix a few small issues (scanf field width etc.) 16-mar-2008 ejk
    convert format to 4 char TAB size and add vzatomLUT() 11-jun-2008 ejk
    convert freqn() for use with n not a power of 2 on 21-feb-2010
    add subroutines for use with FFTW 18-mar-2010 ejk
    put conditional compilation around FFTW portion and add parameter
       offsets for multipole aberrations 12-apr-2011 ejk
    add readCnm() and chi() and ns to askYN() 8-may-2011 ejk
    fix error in C12a,b= astigmatism calculation 4-jul-2012 ejk
    add messageSL() and convert all printf output to go thru it
       for redirections under a GUI  14-apr-2013
   convert char[] and sprintf() to std::string  and
       add 3 versions of toString() 24-nov-2013 ejk
   change printfOK to USE_TERMINAL and
        start conversion to streams 10-mar-2014 ejk
   remove ReadLine() and use C++ getline() 4-apr-2014 ejk
   add ranPoisson() 27-sep-2015 ejk
   convert malloc1D() etc to vector<> 28-jun-2016 ejk
   remove malloc2D(), malloc3D(), free2D(), free3D() and
      make separate newnd.hpp with new/del template 22-sep-2017 ejk
   fix typo in comments for ranPoisson() 29-nov-2017 ejk
   update comments in ReadfeTable( ) and remove redundant declarations
        of ReadfeTable( )  24-apr-2018 ejk
   fix  freqn/xo small error 3-jul-2019 ejk
   remove propagate() so slicelib is not dependent on cfpix+fftw 29-jul-2019 ejk
*/


#include <stdio.h>  /* ANSI C libraries */
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <ctype.h>
#include <cmath>
#include <time.h>

#include <string>
#include <iostream>
#include <sstream>  // string streams
#include <fstream>  // STD file IO streams
#include <vector>   // STD vector class

using namespace std;

/*#include <omp.h>   for openMP testing */

#include "slicelib.hpp"  /* verify consistency and use some routines */

//#define wxGUI   //    set for wxWidgets graphical user interface

#ifdef wxGUI
#include "wx/wx.h"  // regular headers
//#include "wx/wxprec.h"  // for precompiled headers
#endif 

// for the atomic scattering factor tables
const int NPMAX=   12;  // number of parameters for each Z
const int NZMIN=   1;   // min Z 
const int NZMAX=   103; // max Z 

int feTableRead=0;      // flag to remember if the param file has been read
const int nl=3, ng=3;   // number of Lorenzians and Gaussians

vector<double> tmp( NPMAX );       // to get fe(k) parameters
vector< vector<double> > fparams( NZMAX+1, tmp );

/*--------------------- askYN() -----------------------------------*/
/*
    ask a yes/no question and return 1 or 0 for TRUE/FALSE

    message[] = question to ask
*/
#ifdef USE_TERMINAL
int askYN( const char message[] )
{
    string cline;

    cout << message << " (y/n) :" << endl;
    cin >> cline;

    if( (cline.compare("Y")==0) || (cline.compare("y")==0) ||
        (cline.compare("1")==0) ) return( 1 );  else return( 0 );

} /* end askYN() */
#endif

/*-------------------- bessi0() ---------------*/
/*
    modified Bessel function I0(x)
    see Abramowitz and Stegun page 379

    x = (double) real arguments

    12-feb-1997 E. Kirkland
 */
 double bessi0( double x )
 {
    int i;
    double ax, sum, t;
    
    static double i0a[] = { 1.0, 3.5156229, 3.0899424, 1.2067492,
        0.2659732, 0.0360768, 0.0045813 };

    static double i0b[] = { 0.39894228, 0.01328592, 0.00225319,
        -0.00157565, 0.00916281, -0.02057706, 0.02635537,
        -0.01647633, 0.00392377};

    ax = fabs( x );
    if( ax <= 3.75 ) {
        t = x / 3.75;
        t = t * t;
        sum = i0a[6];
        for( i=5; i>=0; i--) sum = sum*t + i0a[i]; 
    } else {
        t = 3.75 / ax;
        sum = i0b[8];
        for( i=7; i>=0; i--) sum = sum*t + i0b[i];
        sum = exp( ax ) * sum / sqrt( ax );
    }
    return( sum );

}  /* end bessi0() */

/*-------------------- bessk0() ---------------*/
/*
    modified Bessel function K0(x)
    see Abramowitz and Stegun page 380
    
    Note: K0(0) is not define and this function
    returns 1E20
 
    x = (double) real arguments
    
    this routine calls bessi0() = Bessel function I0(x)
    
    12-feb-1997 E. Kirkland
 */
 double bessk0( double x )
 {
    double bessi0(double);
 
    int i;
    double ax, x2, sum;
    static double k0a[] = { -0.57721566, 0.42278420, 0.23069756,
         0.03488590, 0.00262698, 0.00010750, 0.00000740};
        
    static double k0b[] = { 1.25331414, -0.07832358, 0.02189568,
         -0.01062446, 0.00587872, -0.00251540, 0.00053208};

    ax = fabs( x );
    if( (ax > 0.0)  && ( ax <=  2.0 ) ){
        x2 = ax/2.0;
        x2 = x2 * x2;
        sum = k0a[6];
        for( i=5; i>=0; i--) sum = sum*x2 + k0a[i];
        sum = -log(ax/2.0) * bessi0(x) + sum;
    } else if( ax > 2.0 ) {
        x2 = 2.0/ax;
        sum = k0b[6];
        for( i=5; i>=0; i--) sum = sum*x2 + k0b[i];
        sum = exp( -ax ) * sum / sqrt( ax );
    } else sum = 1.0e20;
    return ( sum );

}  /* end bessk0() */

/*------------------------ chi() ---------------------*/
/*  return the aberration function
    must multiply this result by 2pi/wavelength before using
    C10= -df is defocus with the opposite sign

    put in a subroutine so I only have to get it right once
    
    now includes thru 5th order 4-jul-2010 ejk
    convert to param array p[] 2-may-2011 ejk
    add multiMode 8-may-2011 ejk
    fix bug in astigmatism calculation 4-jul-2012 ejk

input:
    p = param[] array with aberration coeff. of chi (in Angstroms)
    alx, aly = x,y components of alpha in radians
    multiMode = if not 0 then use multipole aberrations
           set to 0 to speed up the calculation (no C34a etc.)

*/
double chi( vectorf &p, double alx, double aly, int multiMode )
{
    double theta, al, al2, aln;
    double w, ct, st, c2t, s2t, c3t, s3t, c4t, s4t, c5t, s5t,
        c6t, s6t;

    aln = al2 = alx*alx + aly*aly;    /*  alpha squared */

    /*  just rotationally symm. aberrations (faster) */
    w = ( ( al2*p[pCS5]/6.0 +  0.25*p[pCS] )*al2 - 0.5*p[pDEFOCUS] )*al2;

    if( multiMode != 0 ) {
            /* ---- first order ----- */
            theta = atan2( aly, alx );
            ct = cos( theta );
            st = sin( theta );
            c2t = ct*ct - st*st;    /*  cos/sin of 2*theta */
            s2t = 2.0*ct*st;
            w += al2*(  p[pC12a]*c2t + p[pC12b]*s2t )/2.0;

            al = sqrt( al2 );

            /* ---- second order ----- */
            /*   generate the other theta's recursively to reduce CPU time */
            aln = al2*al;  /* alpha^3 */
            c3t = ct*c2t - st*s2t;    /*  cos/sin of 3*theta */
            s3t = ct*s2t + st*c2t;
            w += aln*( p[pC21a]*ct + p[pC21b]*st + p[pC23a]*c3t + p[pC23b]*s3t )/3.0;

            /* ---- third order ----- */
            aln = al2*al2;  /* alpha^4 */
            c4t = ct*c3t - st*s3t;    /*  cos/sin of 4*theta */
            s4t = ct*s3t + st*c3t;
            w += aln*( p[pC32a]*c2t + p[pC32b]*s2t + p[pC34a]*c4t 
                     + p[pC34b]*s4t )/4.0;

             /* ---- fourth order ----- */
            aln = aln*al;  /* alpha^5 */
            c5t = ct*c4t - st*s4t;    /*  cos/sin of 5*theta */
            s5t = ct*s4t + st*c4t;
            w += aln*( p[pC41a]*ct + p[pC41b]*st +  p[pC43a]*c3t + p[pC43b]*s3t
                     + p[pC45a]*c5t + p[pC45b]*s5t )/5.0;

              /* ---- fifth order ----- */
            aln = aln*al;  /* alpha^6 */
            c6t = ct*c5t - st*s5t;    /*  cos/sin of 5*theta */
            s6t = ct*s5t + st*c5t;
            w += aln*(  p[pC52a]*c2t + p[pC52b]*s2t +  p[pC54a]*c4t + p[pC54b]*s4t
                      + p[pC56a]*c6t + p[pC56b]*s6t )/6.0;
           
#ifdef TEST_COS
            string stemp;
            ct = cos( 6*theta );    /*  quick test of recursion algebra */
            st = sin( 6*theta );
            if( fabs(ct-c6t)+fabs(st-s6t) > 1.0e-10 ) {
                stemp = "cos/sin= " +toString(ct)+", "+ toString(st),", c6t/s6t= "+
                        toString(c6t) +", "+ toString(s6t );
                messageSL( stemp.c_str() );
            }
#endif
    }

    return( w );
}

/*--------------------- cputim() -----------------------------------*/
/*
   retrieve current CPU time in seconds
*/

double cputim()
{
    return ( ( (double)clock() ) / ( (double)CLOCKS_PER_SEC) );

}  /* end cputim() */

/*--------------------- featom() -----------------------------------*/
/*
    return the electron scattering factor for atomic
    number Z at scattering angle k

    Z = atomic number 1 <= Z <= 103
    k2  = k*k where k =1/d = scattering angle (in 1/A)

  assumed global vars:

#define NZMIN   1    min Z in featom.tab 
#define NZMAX   103  max Z in featom.tab 

int feTableRead=0; = flag to remember if the param file has been read 
int nl=3, ng=3; = number of Lorenzians and Gaussians 
double fparams[][] = fe parameters

*/

double featom( int Z, double k2 )
{
   int i, nfe;
   double sum;

   if( (Z<NZMIN) || (Z>NZMAX) ) return( 0.0 );

    /* read in the table from a file if this is the
        first time this is called */
    if( feTableRead == 0 ) nfe = ReadfeTable();

    sum = 0.0;

    /* Lorenztians */
    for( i=0; i<2*nl; i+=2 )
        sum += fparams[Z][i]/( k2 + fparams[Z][i+1] );

    /* Gaussians */
    for( i=2*nl; i<2*(nl+ng); i+=2 )
        sum += fparams[Z][i]*exp( - k2 * fparams[Z][i+1] );

    return( sum );

}  /* end featom() */


/*------------------------ freqn() ------------------------*/
/*
    Calculate spatial frequencies for use with fft's
    NOTE: zero freg is in the bottom left corner and
        expands into all other corners - not in the center
        this is required for fft - don't waste time rearranging

    This routine must be called once for each direction

    ko[n]  = real array to get spatial frequencies
    ko2[n] = real array to get k[i]*k[i]
    xo[n]  = real array to get positions 
    nk     = integer number of pixels
    ak     = real full scale size of image in pixels
*/
void freqn( vectorf &ko, vectorf &ko2, vectorf &xo, int nk, double ak )
{
    int i, imid;

    /*   imid = nk/2;  only for nk=2^m - very small error otherwise */
    imid = (int) ( nk/2.0 + 0.5);   /* when nk may not be 2^m */

    for( i=0; i<nk; i++) {
        xo[i] = ((float) (i * ak) ) / ((float)(nk));
        if ( i > imid ) {
            ko[i]  = ((float)(i-nk)) / ((float)ak);
        } else {
            ko[i]  = ((float)i) / ((float)ak);
        }
        ko2[i] = ko[i] * ko[i];
    }

}  /*  end freqn() */


/*--------------------- messageSL() -----------------------------------*/
/*
   slicelib message output
   direct all output message here to redirect to the command line
   or a GUI status line or message box when appropriate

   msg[] = character string with message to disply
   level = level of seriousness
            0 = simple status message
        1 = significant warning
        2 = possibly fatal error
*/
void messageSL( const char msg[],  int level )
{
#ifdef USE_TERMINAL
    cout << msg << endl;  //simplest possible version
#endif

#ifdef wxGUI
    if( 0 == level ){
        //wxLogStatus( wxT(msg)  );
        wxLogStatus( msg );

    } else if( 1 == level ) {
        wxMessageBox( msg, wxT("slicelib"), wxOK | wxICON_INFORMATION ); 
    }
#endif

}

/*--------------------- parlay() -----------------------------------*/
/*
  subroutine to parse the atomic layer stacking sequence 
  for use with multislice programs.

  This converts layer structure definition of the form:
      2(abc)d
  into a sequence of numerical indices where a=1, b=2, c=3, ...
  The main attraction is the repeat operator n(...) where
  the structure inside the parenthesis (...) is repeated n times.
  for instance the above structure is translated into:
      1 2 3 1 2 3 4
  The parenthesis may be nested up to 100 levels (determined by
  nlmax). For instance  5( 2(abc) 3(cba) ) is also a valid structure
  definition. This is a compact way of specifying the layer structure
  for a multislice calculation. tab's may be present in the structure
  anywhere (they are ignored).

  This is done by pushing the position of each '(' and its repeat
  count onto a stack. Each ')' pops the last entry from the stack
  and invokes a duplication process.

  Layers refer to the distinquishable subset of different
  types of layers, whereas slices refer to the way these
  layers are sequentially stacked.

  fortran version started 22-dec-1989 earl j. kirkland
  added nested parenthesis by stacking entry points
     31-jan-1990 ejk
  added tab handling (in ascii and ebcdic) 1-feb-1990 ejk
  minor changes to error messages to make rs6000's happy
      7-july-1992 ejk
  converted to ANSI-C 26-jun-1995 ejk
  added include of stdlib.h 6-july-1995 ejk

   c         = input character string
   islice[nsmax] = integer array to get stacking sequence indicies
   nsmax     = (integer) size of layer
   lmax      = (integer) maximum allowed layer index
   nslice    = (integer) number of layers
   returned value= (integer) success/error code
               0 : success
              -1 : layer out of range
              -2 : missing left parenthesis
              -3 : parenthesis nested too deep
              -4 : bad repeat code
              -5 : unmatched right parenthesis
              -6 : invalid character
              -7 : too many layers
              -8 : incomplete stacking sequence

   fperr       = (int) if this is not a NULL then write
           error messages 

  NOTE: islice and nslice may be modified by this routine

  cname determines the mapping of characters into numbers.

*/
#ifdef USE_TERMINAL

#define NCMAX   132 // characters per line to read

#define NSTKMAX 100 /* maximum stack depth (i.e. maximum level
                 of nesting the parenthesis) */
#define NCHARS 52   /* maximum number of character symbols */

int parlay( const char c[], vectori &islice, int nsmax, int lmax,
            int *nslice, int fperr )
{
    /* define our own symbol sequence */
    const char cname[] =
       "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";

    int ic, lenc, name, i, j, istack, n,
        ipoint[NSTKMAX], repeat[NSTKMAX];

    /*  initialize misc constants  */

    *nslice = 0;
    istack = 0;
    lenc = (int) strlen( c );  /* length of c */

    for( ic=0; ic<lenc; ic++ ) {    /* loop over all characters in c */

       /*  skip over embedded spaces, tabs and weird characters */
       while( isspace( c[ic] )  ) ic++;

       /*  if its a character do this  */       
       if( isalpha( c[ic] ) ) {
            for(name=0; name<NCHARS; name++) 
               if( c[ic] == cname[name] ) break;
            if( name > lmax ) {
               if( fperr != 0 ) 
                  printf("Layer /%c/ out of range in the stacking sequence"
                " in subroutine parlay.\n", c[ic] );
                return( -1 );
            } else {
               if( *nslice >= nsmax ) {
                if( fperr != 0 ) 
                printf("Too many layers generated in the stacking"
                    " sequence in parlay.\n");
                    return( -7 );
               }
                islice[ (*nslice)++ ] = name;
            }
        
    /*  if its a number then extract repeat count up to the '('
    and save it until a ')' appears (i.e. push on the stack)
    */
       } else if ( isdigit( c[ic] ) ) {
        if( istack >= NSTKMAX ) {
                if( fperr != 0 ) 
                printf("Parenthesis nested too deep in "
                   "the stacking sequence in subroutine parlay.\n");
                return( -3 );
            }
            repeat[istack] = atoi( &c[ic] ) - 1; 
            while( isdigit(c[ic]) || isspace(c[ic]) ) ic++;
            if( c[ic] != '(' ) {
                if( fperr != 0 ) 
                   printf("Missing left parenthesis in "
                   "the stacking sequence in subroutine parlay.'\n");
                   for(i=0; i<=ic; i++) printf("%c", c[i]);
                   printf("?\n");
                   return( -2 );
            }
            ipoint[istack++] = *nslice;
    
        /*  if its a ')' then repeat back to the last '('
                    (i.e. pop the stack) */

        } else if( c[ic] == ')' ) {
           if( istack <= 0 ) {
                if( fperr != 0 ) 
                printf("Unmatched right parenthesis in "
                    "the stacking sequence in subroutine parlay.\n");
                return( -5 );
           }
           n = *nslice;
           istack--;
           for( j=0; j<repeat[istack]; j++)
            for(i=ipoint[istack]; i<n; i++){
               if( *nslice >= nsmax ) {
                    if( fperr != 0 ) 
                    printf("Too many layers generated in the stacking"
                        " sequence in parlay.\n");
                        return( -7 );
                }
                islice[ (*nslice)++ ] = islice[ i ];
            }
        } else {
           if( fperr != 0 ) 
           printf("Invalid character /%c/ encountered in "
              "the stacking sequence in subroutine parlay.\n",
              c[ic]);
           return( -6 );
        }

    } /* end for( ic... */
    
    if( istack != 0 ) {
        if( fperr != 0 ) 
        printf("incomplete stacking sequence in parlay.\n");
        return( -8 ); 
    } else return( 0 );

#undef NCMAX
#undef NSTKMAX 
#undef NCHARS

}  /* end parlay() */
#endif

/*---------------------------- ranflat -------------------------------*/
/*
    return a random number in the range 0.0->1.0
    with uniform distribution

    the 'Magic Numbers' are from 
        Numerical Recipes 2nd edit pg. 285
*/
double ranflat( unsigned long *iseed )
{
    static unsigned long a=1366, c=150889L, m=714025L;
    
    *iseed = ( a * (*iseed) + c ) % m;
    
    return( ((double) *iseed)/m);

}  /* end ranflat() */

/*-------------------- rangauss() -------------------------- */
/*
    Return a normally distributed random number with 
    zero mean and unit variance using Box-Muller method
    
    ranflat() is the source of uniform deviates

    ref.  Numerical Recipes, 2nd edit. page 289

    added log(0) test 10-jan-1998 E. Kirkland
*/
double rangauss( unsigned long *iseed )
{
    double ranflat( unsigned long* );
    double x1, x2, y;
    static double tpi=2.0*3.141592654;

    /* be careful to avoid taking log(0) */
    do{
        x1 = ranflat( iseed );
        x2 = ranflat( iseed );

    } while ( (x1 < 1.0e-30) || (x2 < 1.0e-30) );

    y = sqrt( - 2.0 * log(x1) ) * cos( tpi * x2 );

    return( y );

}  /* end rangauss() */


/*---------------------------- ranPoisson -------------------------------*/
/*
    return a random number with Poisson distribution
    and mean value mean
    
    There is a nice poisson RNG in c++11 but some compilers still do not
    have it, so make one here. The function lgamma()= log of gamma function
    is in the C99 standard and seems to have more widespread support
    (in mac osx 10.8 but NOT MSVS2010!).  Leave lgamma() code for future
    and approximate large means with Stirling's formula.
    
    A. C, Atkinson, "The Computer Generation of Poisson Random Variable"
         J. Royal Statistical Society, Series C (Applied Statistics),
         Vol. 28, No. 1, (1979) p. 29-35.
         
    D. E. Knuth, "The Art of Computer Programming, vol.2, Seminumerical
        Algorithms"  Addison-Wesley 1981, 1969, p. 132.

    calls ranflat()
    
    input: 
        mean = desired mean (can be fractional)
        iseed = random number seed
    
    started 23-sep-2015 ejk
    add large mean portion 26-sep-2015 ejk
*/
int ranPoisson( double mean, unsigned long *iseed )
{
    int n;
    static int lut=-1;
    static double oldMean=-100, oldEmean=0;
    static double alpha, beta, c, k, PI, lnf[256];
    
    //  negative mean is not allowed
    if( mean <= 0 ) return 0;
    
    if( lut < 0 ) {   // init log( n! ) look up table
        lnf[0] = lnf[1] = 0.0;
        for( int i=2; i<256; i++) lnf[i] = lnf[i-1] + log((double)i);
        lut = +1;
    }
    
    //---- use Atkinson method PM for small means
    //     (also Knuth vol.2)
    if( mean < 30 ) {
        
        double s;
        if( oldMean != mean ) {  //  save if same mean is repeated
            oldMean = mean;
            oldEmean = exp( -mean );
        }
        n= -1;
        s = 1.0;
        do{
           n = n + 1;
           s = s * ranflat( iseed );
        } while( s >= oldEmean );
        return(  n );
        
    //---  use Atkinson method PA for large means 
    } else {
        
        double u1, u2, x, y, lhs, rhs, temp, lnf0;
    
        if( oldMean != mean ) {  //  save if same mean is repeated
            oldMean = mean;
            PI=4.0*atan(1.0);
            beta = PI/sqrt(3.0*mean);
            alpha = beta*mean;
            c = 0.767 - 3.36/mean;
            k = log(c) - mean - log(beta);
        }
        
        while(1) {      //  hopefully this will end !!
            
            do{     u1 = ranflat( iseed );
                    x = (alpha - log((1.0 - u1)/u1))/beta;
            } while( x < -0.5 );
            n = (int) (x + 0.5);
            u2 = ranflat( iseed );
            y = alpha - beta*x;
            temp = 1.0 + exp(y);
            lhs = y + log(u2/(temp*temp));
            //  lgamma() is missing in some compilers so.....
            //rhs = k + n*log(mean) - lgamma(n+1); // gamma(n+1) = n!
            if( n < 0 ) continue;
            if( n < 255 ) lnf0 = lnf[n]; // log( n! ) from lookup table
            else {
                x = n;   //  Stirling Formula
                lnf0 = 0.5*log(2*PI) + (x+0.5)*log(x) - x + 1.0/(12.0*x);
            }
            rhs = k + n*log(mean) - lnf0;
            if (lhs <= rhs) return n;
        }
    }

}  // end ranPoisson()


/*--------------------- readCnm() -----------------------*/
/*
    convert aberration line to a number 

    cline = character line in style "C32a  2.4"
    param[] = parameter array to get parameter value

    return >= 0 for success and <0  if parameter not recognized

    convert to C++ string 30-mar-2014 ejk
    convert param to vector<float> 2-jul-2016 ejk
*/
int readCnm( string &cline, vectorf &param, double x )
{
    /*  so we can loop through them in some ordered manner 
       - next two array MUST be in identical order 
       - leave out C10 because its opp. sign of df=defocus */
    int aberIndex[]={pC12a, pC12b, pC21a, pC21b, pC23a, pC23b, pCS,
        pC32a, pC32b, pC34a, pC34b, pC41a, pC41b, pC43a, pC43b, pC45a, pC45b, 
        pCS5, pC52a, pC52b, pC54a, pC54b, pC56a, pC56b };

    string aberSymb[]= {"C12a", "C12b", "C21a", "C21b", "C23a",
        "C23b", "C30", "C32a", "C32b", "C34a", "C34b", "C41a", "C41b", 
        "C43a", "C43b", "C45a", "C45b", "C50", "C52a", "C52b", "C54a", 
        "C54b", "C56a", "C56b" };
    int i, status=-1, Naber=24;

    for(i=0; i<Naber; i++) {
        if( cline.compare( aberSymb[i] ) == 0 ) {
            param[aberIndex[i]] = (float) (x * 1.0e7);   /* mm. to Ang. */
            status = aberIndex[i];
        }
        if( status >= 0 ) break;
    }

    return( status );

}  /* end readCnm() */

/*--------------------- ReadfeTable() -----------------------*/
/*
    initialize electron scattering factors parameters

  the old version read electron scattering factors
  parameters from file fparam.dat hence the name
  
  the constants that must be defined above are

const int NPMAX=   12;  // number of parameters for each Z
const int NZMIN=   1;   // min Z 
const int NZMAX=   103; // max Z 

  assumed global vars:

int feTableRead=0; = flag to remember if the param file has been read 
int nl=3, ng=3; = number of Lorenzians and Gaussians 
vector< vector<double> > fparams( NZMAX+1, tmp );

*/
int ReadfeTable( )
{
    /*  hard code the scattering factor parameters so you
     don't have to keep track of the extra file and also
     to guarantee that the data won't be changed*/

   /* if the file has been initialized already then just return */
   if( feTableRead == 1 ) return(0);

    /*  scattering factor parameter data from fparams  */

    fparams[1][0] =  4.20298324e-003 ;
    fparams[1][1] =  2.25350888e-001 ;
    fparams[1][2] =  6.27762505e-002 ;
    fparams[1][3] =  2.25366950e-001 ;
    fparams[1][4] =  3.00907347e-002 ;
    fparams[1][5] =  2.25331756e-001 ;
    fparams[1][6] =  6.77756695e-002 ;
    fparams[1][7] =  4.38854001e+000 ;
    fparams[1][8] =  3.56609237e-003 ;
    fparams[1][9] =  4.03884823e-001 ;
    fparams[1][10] =  2.76135815e-002 ;
    fparams[1][11] =  1.44490166e+000 ;
    fparams[2][0] =  1.87543704e-005 ;
    fparams[2][1] =  2.12427997e-001 ;
    fparams[2][2] =  4.10595800e-004 ;
    fparams[2][3] =  3.32212279e-001 ;
    fparams[2][4] =  1.96300059e-001 ;
    fparams[2][5] =  5.17325152e-001 ;
    fparams[2][6] =  8.36015738e-003 ;
    fparams[2][7] =  3.66668239e-001 ;
    fparams[2][8] =  2.95102022e-002 ;
    fparams[2][9] =  1.37171827e+000 ;
    fparams[2][10] =  4.65928982e-007 ;
    fparams[2][11] =  3.75768025e+004 ;
    fparams[3][0] =  7.45843816e-002 ;
    fparams[3][1] =  8.81151424e-001 ;
    fparams[3][2] =  7.15382250e-002 ;
    fparams[3][3] =  4.59142904e-002 ;
    fparams[3][4] =  1.45315229e-001 ;
    fparams[3][5] =  8.81301714e-001 ;
    fparams[3][6] =  1.12125769e+000 ;
    fparams[3][7] =  1.88483665e+001 ;
    fparams[3][8] =  2.51736525e-003 ;
    fparams[3][9] =  1.59189995e-001 ;
    fparams[3][10] =  3.58434971e-001 ;
    fparams[3][11] =  6.12371000e+000 ;
    fparams[4][0] =  6.11642897e-002 ;
    fparams[4][1] =  9.90182132e-002 ;
    fparams[4][2] =  1.25755034e-001 ;
    fparams[4][3] =  9.90272412e-002 ;
    fparams[4][4] =  2.00831548e-001 ;
    fparams[4][5] =  1.87392509e+000 ;
    fparams[4][6] =  7.87242876e-001 ;
    fparams[4][7] =  9.32794929e+000 ;
    fparams[4][8] =  1.58847850e-003 ;
    fparams[4][9] =  8.91900236e-002 ;
    fparams[4][10] =  2.73962031e-001 ;
    fparams[4][11] =  3.20687658e+000 ;
    fparams[5][0] =  1.25716066e-001 ;
    fparams[5][1] =  1.48258830e-001 ;
    fparams[5][2] =  1.73314452e-001 ;
    fparams[5][3] =  1.48257216e-001 ;
    fparams[5][4] =  1.84774811e-001 ;
    fparams[5][5] =  3.34227311e+000 ;
    fparams[5][6] =  1.95250221e-001 ;
    fparams[5][7] =  1.97339463e+000 ;
    fparams[5][8] =  5.29642075e-001 ;
    fparams[5][9] =  5.70035553e+000 ;
    fparams[5][10] =  1.08230500e-003 ;
    fparams[5][11] =  5.64857237e-002 ;
    fparams[6][0] =  2.12080767e-001 ;
    fparams[6][1] =  2.08605417e-001 ;
    fparams[6][2] =  1.99811865e-001 ;
    fparams[6][3] =  2.08610186e-001 ;
    fparams[6][4] =  1.68254385e-001 ;
    fparams[6][5] =  5.57870773e+000 ;
    fparams[6][6] =  1.42048360e-001 ;
    fparams[6][7] =  1.33311887e+000 ;
    fparams[6][8] =  3.63830672e-001 ;
    fparams[6][9] =  3.80800263e+000 ;
    fparams[6][10] =  8.35012044e-004 ;
    fparams[6][11] =  4.03982620e-002 ;
    fparams[7][0] =  5.33015554e-001 ;
    fparams[7][1] =  2.90952515e-001 ;
    fparams[7][2] =  5.29008883e-002 ;
    fparams[7][3] =  1.03547896e+001 ;
    fparams[7][4] =  9.24159648e-002 ;
    fparams[7][5] =  1.03540028e+001 ;
    fparams[7][6] =  2.61799101e-001 ;
    fparams[7][7] =  2.76252723e+000 ;
    fparams[7][8] =  8.80262108e-004 ;
    fparams[7][9] =  3.47681236e-002 ;
    fparams[7][10] =  1.10166555e-001 ;
    fparams[7][11] =  9.93421736e-001 ;
    fparams[8][0] =  3.39969204e-001 ;
    fparams[8][1] =  3.81570280e-001 ;
    fparams[8][2] =  3.07570172e-001 ;
    fparams[8][3] =  3.81571436e-001 ;
    fparams[8][4] =  1.30369072e-001 ;
    fparams[8][5] =  1.91919745e+001 ;
    fparams[8][6] =  8.83326058e-002 ;
    fparams[8][7] =  7.60635525e-001 ;
    fparams[8][8] =  1.96586700e-001 ;
    fparams[8][9] =  2.07401094e+000 ;
    fparams[8][10] =  9.96220028e-004 ;
    fparams[8][11] =  3.03266869e-002 ;
    fparams[9][0] =  2.30560593e-001 ;
    fparams[9][1] =  4.80754213e-001 ;
    fparams[9][2] =  5.26889648e-001 ;
    fparams[9][3] =  4.80763895e-001 ;
    fparams[9][4] =  1.24346755e-001 ;
    fparams[9][5] =  3.95306720e+001 ;
    fparams[9][6] =  1.24616894e-003 ;
    fparams[9][7] =  2.62181803e-002 ;
    fparams[9][8] =  7.20452555e-002 ;
    fparams[9][9] =  5.92495593e-001 ;
    fparams[9][10] =  1.53075777e-001 ;
    fparams[9][11] =  1.59127671e+000 ;
    fparams[10][0] =  4.08371771e-001 ;
    fparams[10][1] =  5.88228627e-001 ;
    fparams[10][2] =  4.54418858e-001 ;
    fparams[10][3] =  5.88288655e-001 ;
    fparams[10][4] =  1.44564923e-001 ;
    fparams[10][5] =  1.21246013e+002 ;
    fparams[10][6] =  5.91531395e-002 ;
    fparams[10][7] =  4.63963540e-001 ;
    fparams[10][8] =  1.24003718e-001 ;
    fparams[10][9] =  1.23413025e+000 ;
    fparams[10][10] =  1.64986037e-003 ;
    fparams[10][11] =  2.05869217e-002 ;
    fparams[11][0] =  1.36471662e-001 ;
    fparams[11][1] =  4.99965301e-002 ;
    fparams[11][2] =  7.70677865e-001 ;
    fparams[11][3] =  8.81899664e-001 ;
    fparams[11][4] =  1.56862014e-001 ;
    fparams[11][5] =  1.61768579e+001 ;
    fparams[11][6] =  9.96821513e-001 ;
    fparams[11][7] =  2.00132610e+001 ;
    fparams[11][8] =  3.80304670e-002 ;
    fparams[11][9] =  2.60516254e-001 ;
    fparams[11][10] =  1.27685089e-001 ;
    fparams[11][11] =  6.99559329e-001 ;
    fparams[12][0] =  3.04384121e-001 ;
    fparams[12][1] =  8.42014377e-002 ;
    fparams[12][2] =  7.56270563e-001 ;
    fparams[12][3] =  1.64065598e+000 ;
    fparams[12][4] =  1.01164809e-001 ;
    fparams[12][5] =  2.97142975e+001 ;
    fparams[12][6] =  3.45203403e-002 ;
    fparams[12][7] =  2.16596094e-001 ;
    fparams[12][8] =  9.71751327e-001 ;
    fparams[12][9] =  1.21236852e+001 ;
    fparams[12][10] =  1.20593012e-001 ;
    fparams[12][11] =  5.60865838e-001 ;
    fparams[13][0] =  7.77419424e-001 ;
    fparams[13][1] =  2.71058227e+000 ;
    fparams[13][2] =  5.78312036e-002 ;
    fparams[13][3] =  7.17532098e+001 ;
    fparams[13][4] =  4.26386499e-001 ;
    fparams[13][5] =  9.13331555e-002 ;
    fparams[13][6] =  1.13407220e-001 ;
    fparams[13][7] =  4.48867451e-001 ;
    fparams[13][8] =  7.90114035e-001 ;
    fparams[13][9] =  8.66366718e+000 ;
    fparams[13][10] =  3.23293496e-002 ;
    fparams[13][11] =  1.78503463e-001 ;
    fparams[14][0] =  1.06543892e+000 ;
    fparams[14][1] =  1.04118455e+000 ;
    fparams[14][2] =  1.20143691e-001 ;
    fparams[14][3] =  6.87113368e+001 ;
    fparams[14][4] =  1.80915263e-001 ;
    fparams[14][5] =  8.87533926e-002 ;
    fparams[14][6] =  1.12065620e+000 ;
    fparams[14][7] =  3.70062619e+000 ;
    fparams[14][8] =  3.05452816e-002 ;
    fparams[14][9] =  2.14097897e-001 ;
    fparams[14][10] =  1.59963502e+000 ;
    fparams[14][11] =  9.99096638e+000 ;
    fparams[15][0] =  1.05284447e+000 ;
    fparams[15][1] =  1.31962590e+000 ;
    fparams[15][2] =  2.99440284e-001 ;
    fparams[15][3] =  1.28460520e-001 ;
    fparams[15][4] =  1.17460748e-001 ;
    fparams[15][5] =  1.02190163e+002 ;
    fparams[15][6] =  9.60643452e-001 ;
    fparams[15][7] =  2.87477555e+000 ;
    fparams[15][8] =  2.63555748e-002 ;
    fparams[15][9] =  1.82076844e-001 ;
    fparams[15][10] =  1.38059330e+000 ;
    fparams[15][11] =  7.49165526e+000 ;
    fparams[16][0] =  1.01646916e+000 ;
    fparams[16][1] =  1.69181965e+000 ;
    fparams[16][2] =  4.41766748e-001 ;
    fparams[16][3] =  1.74180288e-001 ;
    fparams[16][4] =  1.21503863e-001 ;
    fparams[16][5] =  1.67011091e+002 ;
    fparams[16][6] =  8.27966670e-001 ;
    fparams[16][7] =  2.30342810e+000 ;
    fparams[16][8] =  2.33022533e-002 ;
    fparams[16][9] =  1.56954150e-001 ;
    fparams[16][10] =  1.18302846e+000 ;
    fparams[16][11] =  5.85782891e+000 ;
    fparams[17][0] =  9.44221116e-001 ;
    fparams[17][1] =  2.40052374e-001 ;
    fparams[17][2] =  4.37322049e-001 ;
    fparams[17][3] =  9.30510439e+000 ;
    fparams[17][4] =  2.54547926e-001 ;
    fparams[17][5] =  9.30486346e+000 ;
    fparams[17][6] =  5.47763323e-002 ;
    fparams[17][7] =  1.68655688e-001 ;
    fparams[17][8] =  8.00087488e-001 ;
    fparams[17][9] =  2.97849774e+000 ;
    fparams[17][10] =  1.07488641e-002 ;
    fparams[17][11] =  6.84240646e-002 ;
    fparams[18][0] =  1.06983288e+000 ;
    fparams[18][1] =  2.87791022e-001 ;
    fparams[18][2] =  4.24631786e-001 ;
    fparams[18][3] =  1.24156957e+001 ;
    fparams[18][4] =  2.43897949e-001 ;
    fparams[18][5] =  1.24158868e+001 ;
    fparams[18][6] =  4.79446296e-002 ;
    fparams[18][7] =  1.36979796e-001 ;
    fparams[18][8] =  7.64958952e-001 ;
    fparams[18][9] =  2.43940729e+000 ;
    fparams[18][10] =  8.23128431e-003 ;
    fparams[18][11] =  5.27258749e-002 ;
    fparams[19][0] =  6.92717865e-001 ;
    fparams[19][1] =  7.10849990e+000 ;
    fparams[19][2] =  9.65161085e-001 ;
    fparams[19][3] =  3.57532901e-001 ;
    fparams[19][4] =  1.48466588e-001 ;
    fparams[19][5] =  3.93763275e-002 ;
    fparams[19][6] =  2.64645027e-002 ;
    fparams[19][7] =  1.03591321e-001 ;
    fparams[19][8] =  1.80883768e+000 ;
    fparams[19][9] =  3.22845199e+001 ;
    fparams[19][10] =  5.43900018e-001 ;
    fparams[19][11] =  1.67791374e+000 ;
    fparams[20][0] =  3.66902871e-001 ;
    fparams[20][1] =  6.14274129e-002 ;
    fparams[20][2] =  8.66378999e-001 ;
    fparams[20][3] =  5.70881727e-001 ;
    fparams[20][4] =  6.67203300e-001 ;
    fparams[20][5] =  7.82965639e+000 ;
    fparams[20][6] =  4.87743636e-001 ;
    fparams[20][7] =  1.32531318e+000 ;
    fparams[20][8] =  1.82406314e+000 ;
    fparams[20][9] =  2.10056032e+001 ;
    fparams[20][10] =  2.20248453e-002 ;
    fparams[20][11] =  9.11853450e-002 ;
    fparams[21][0] =  3.78871777e-001 ;
    fparams[21][1] =  6.98910162e-002 ;
    fparams[21][2] =  9.00022505e-001 ;
    fparams[21][3] =  5.21061541e-001 ;
    fparams[21][4] =  7.15288914e-001 ;
    fparams[21][5] =  7.87707920e+000 ;
    fparams[21][6] =  1.88640973e-002 ;
    fparams[21][7] =  8.17512708e-002 ;
    fparams[21][8] =  4.07945949e-001 ;
    fparams[21][9] =  1.11141388e+000 ;
    fparams[21][10] =  1.61786540e+000 ;
    fparams[21][11] =  1.80840759e+001 ;
    fparams[22][0] =  3.62383267e-001 ;
    fparams[22][1] =  7.54707114e-002 ;
    fparams[22][2] =  9.84232966e-001 ;
    fparams[22][3] =  4.97757309e-001 ;
    fparams[22][4] =  7.41715642e-001 ;
    fparams[22][5] =  8.17659391e+000 ;
    fparams[22][6] =  3.62555269e-001 ;
    fparams[22][7] =  9.55524906e-001 ;
    fparams[22][8] =  1.49159390e+000 ;
    fparams[22][9] =  1.62221677e+001 ;
    fparams[22][10] =  1.61659509e-002 ;
    fparams[22][11] =  7.33140839e-002 ;
    fparams[23][0] =  3.52961378e-001 ;
    fparams[23][1] =  8.19204103e-002 ;
    fparams[23][2] =  7.46791014e-001 ;
    fparams[23][3] =  8.81189511e+000 ;
    fparams[23][4] =  1.08364068e+000 ;
    fparams[23][5] =  5.10646075e-001 ;
    fparams[23][6] =  1.39013610e+000 ;
    fparams[23][7] =  1.48901841e+001 ;
    fparams[23][8] =  3.31273356e-001 ;
    fparams[23][9] =  8.38543079e-001 ;
    fparams[23][10] =  1.40422612e-002 ;
    fparams[23][11] =  6.57432678e-002 ;
    fparams[24][0] =  1.34348379e+000 ;
    fparams[24][1] =  1.25814353e+000 ;
    fparams[24][2] =  5.07040328e-001 ;
    fparams[24][3] =  1.15042811e+001 ;
    fparams[24][4] =  4.26358955e-001 ;
    fparams[24][5] =  8.53660389e-002 ;
    fparams[24][6] =  1.17241826e-002 ;
    fparams[24][7] =  6.00177061e-002 ;
    fparams[24][8] =  5.11966516e-001 ;
    fparams[24][9] =  1.53772451e+000 ;
    fparams[24][10] =  3.38285828e-001 ;
    fparams[24][11] =  6.62418319e-001 ;
    fparams[25][0] =  3.26697613e-001 ;
    fparams[25][1] =  8.88813083e-002 ;
    fparams[25][2] =  7.17297000e-001 ;
    fparams[25][3] =  1.11300198e+001 ;
    fparams[25][4] =  1.33212464e+000 ;
    fparams[25][5] =  5.82141104e-001 ;
    fparams[25][6] =  2.80801702e-001 ;
    fparams[25][7] =  6.71583145e-001 ;
    fparams[25][8] =  1.15499241e+000 ;
    fparams[25][9] =  1.26825395e+001 ;
    fparams[25][10] =  1.11984488e-002 ;
    fparams[25][11] =  5.32334467e-002 ;
    fparams[26][0] =  3.13454847e-001 ;
    fparams[26][1] =  8.99325756e-002 ;
    fparams[26][2] =  6.89290016e-001 ;
    fparams[26][3] =  1.30366038e+001 ;
    fparams[26][4] =  1.47141531e+000 ;
    fparams[26][5] =  6.33345291e-001 ;
    fparams[26][6] =  1.03298688e+000 ;
    fparams[26][7] =  1.16783425e+001 ;
    fparams[26][8] =  2.58280285e-001 ;
    fparams[26][9] =  6.09116446e-001 ;
    fparams[26][10] =  1.03460690e-002 ;
    fparams[26][11] =  4.81610627e-002 ;
    fparams[27][0] =  3.15878278e-001 ;
    fparams[27][1] =  9.46683246e-002 ;
    fparams[27][2] =  1.60139005e+000 ;
    fparams[27][3] =  6.99436449e-001 ;
    fparams[27][4] =  6.56394338e-001 ;
    fparams[27][5] =  1.56954403e+001 ;
    fparams[27][6] =  9.36746624e-001 ;
    fparams[27][7] =  1.09392410e+001 ;
    fparams[27][8] =  9.77562646e-003 ;
    fparams[27][9] =  4.37446816e-002 ;
    fparams[27][10] =  2.38378578e-001 ;
    fparams[27][11] =  5.56286483e-001 ;
    fparams[28][0] =  1.72254630e+000 ;
    fparams[28][1] =  7.76606908e-001 ;
    fparams[28][2] =  3.29543044e-001 ;
    fparams[28][3] =  1.02262360e-001 ;
    fparams[28][4] =  6.23007200e-001 ;
    fparams[28][5] =  1.94156207e+001 ;
    fparams[28][6] =  9.43496513e-003 ;
    fparams[28][7] =  3.98684596e-002 ;
    fparams[28][8] =  8.54063515e-001 ;
    fparams[28][9] =  1.04078166e+001 ;
    fparams[28][10] =  2.21073515e-001 ;
    fparams[28][11] =  5.10869330e-001 ;
    fparams[29][0] =  3.58774531e-001 ;
    fparams[29][1] =  1.06153463e-001 ;
    fparams[29][2] =  1.76181348e+000 ;
    fparams[29][3] =  1.01640995e+000 ;
    fparams[29][4] =  6.36905053e-001 ;
    fparams[29][5] =  1.53659093e+001 ;
    fparams[29][6] =  7.44930667e-003 ;
    fparams[29][7] =  3.85345989e-002 ;
    fparams[29][8] =  1.89002347e-001 ;
    fparams[29][9] =  3.98427790e-001 ;
    fparams[29][10] =  2.29619589e-001 ;
    fparams[29][11] =  9.01419843e-001 ;
    fparams[30][0] =  5.70893973e-001 ;
    fparams[30][1] =  1.26534614e-001 ;
    fparams[30][2] =  1.98908856e+000 ;
    fparams[30][3] =  2.17781965e+000 ;
    fparams[30][4] =  3.06060585e-001 ;
    fparams[30][5] =  3.78619003e+001 ;
    fparams[30][6] =  2.35600223e-001 ;
    fparams[30][7] =  3.67019041e-001 ;
    fparams[30][8] =  3.97061102e-001 ;
    fparams[30][9] =  8.66419596e-001 ;
    fparams[30][10] =  6.85657228e-003 ;
    fparams[30][11] =  3.35778823e-002 ;
    fparams[31][0] =  6.25528464e-001 ;
    fparams[31][1] =  1.10005650e-001 ;
    fparams[31][2] =  2.05302901e+000 ;
    fparams[31][3] =  2.41095786e+000 ;
    fparams[31][4] =  2.89608120e-001 ;
    fparams[31][5] =  4.78685736e+001 ;
    fparams[31][6] =  2.07910594e-001 ;
    fparams[31][7] =  3.27807224e-001 ;
    fparams[31][8] =  3.45079617e-001 ;
    fparams[31][9] =  7.43139061e-001 ;
    fparams[31][10] =  6.55634298e-003 ;
    fparams[31][11] =  3.09411369e-002 ;
    fparams[32][0] =  5.90952690e-001 ;
    fparams[32][1] =  1.18375976e-001 ;
    fparams[32][2] =  5.39980660e-001 ;
    fparams[32][3] =  7.18937433e+001 ;
    fparams[32][4] =  2.00626188e+000 ;
    fparams[32][5] =  1.39304889e+000 ;
    fparams[32][6] =  7.49705041e-001 ;
    fparams[32][7] =  6.89943350e+000 ;
    fparams[32][8] =  1.83581347e-001 ;
    fparams[32][9] =  3.64667232e-001 ;
    fparams[32][10] =  9.52190743e-003 ;
    fparams[32][11] =  2.69888650e-002 ;
    fparams[33][0] =  7.77875218e-001 ;
    fparams[33][1] =  1.50733157e-001 ;
    fparams[33][2] =  5.93848150e-001 ;
    fparams[33][3] =  1.42882209e+002 ;
    fparams[33][4] =  1.95918751e+000 ;
    fparams[33][5] =  1.74750339e+000 ;
    fparams[33][6] =  1.79880226e-001 ;
    fparams[33][7] =  3.31800852e-001 ;
    fparams[33][8] =  8.63267222e-001 ;
    fparams[33][9] =  5.85490274e+000 ;
    fparams[33][10] =  9.59053427e-003 ;
    fparams[33][11] =  2.33777569e-002 ;
    fparams[34][0] =  9.58390681e-001 ;
    fparams[34][1] =  1.83775557e-001 ;
    fparams[34][2] =  6.03851342e-001 ;
    fparams[34][3] =  1.96819224e+002 ;
    fparams[34][4] =  1.90828931e+000 ;
    fparams[34][5] =  2.15082053e+000 ;
    fparams[34][6] =  1.73885956e-001 ;
    fparams[34][7] =  3.00006024e-001 ;
    fparams[34][8] =  9.35265145e-001 ;
    fparams[34][9] =  4.92471215e+000 ;
    fparams[34][10] =  8.62254658e-003 ;
    fparams[34][11] =  2.12308108e-002 ;
    fparams[35][0] =  1.14136170e+000 ;
    fparams[35][1] =  2.18708710e-001 ;
    fparams[35][2] =  5.18118737e-001 ;
    fparams[35][3] =  1.93916682e+002 ;
    fparams[35][4] =  1.85731975e+000 ;
    fparams[35][5] =  2.65755396e+000 ;
    fparams[35][6] =  1.68217399e-001 ;
    fparams[35][7] =  2.71719918e-001 ;
    fparams[35][8] =  9.75705606e-001 ;
    fparams[35][9] =  4.19482500e+000 ;
    fparams[35][10] =  7.24187871e-003 ;
    fparams[35][11] =  1.99325718e-002 ;
    fparams[36][0] =  3.24386970e-001 ;
    fparams[36][1] =  6.31317973e+001 ;
    fparams[36][2] =  1.31732163e+000 ;
    fparams[36][3] =  2.54706036e-001 ;
    fparams[36][4] =  1.79912614e+000 ;
    fparams[36][5] =  3.23668394e+000 ;
    fparams[36][6] =  4.29961425e-003 ;
    fparams[36][7] =  1.98965610e-002 ;
    fparams[36][8] =  1.00429433e+000 ;
    fparams[36][9] =  3.61094513e+000 ;
    fparams[36][10] =  1.62188197e-001 ;
    fparams[36][11] =  2.45583672e-001 ;
    fparams[37][0] =  2.90445351e-001 ;
    fparams[37][1] =  3.68420227e-002 ;
    fparams[37][2] =  2.44201329e+000 ;
    fparams[37][3] =  1.16013332e+000 ;
    fparams[37][4] =  7.69435449e-001 ;
    fparams[37][5] =  1.69591472e+001 ;
    fparams[37][6] =  1.58687000e+000 ;
    fparams[37][7] =  2.53082574e+000 ;
    fparams[37][8] =  2.81617593e-003 ;
    fparams[37][9] =  1.88577417e-002 ;
    fparams[37][10] =  1.28663830e-001 ;
    fparams[37][11] =  2.10753969e-001 ;
    fparams[38][0] =  1.37373086e-002 ;
    fparams[38][1] =  1.87469061e-002 ;
    fparams[38][2] =  1.97548672e+000 ;
    fparams[38][3] =  6.36079230e+000 ;
    fparams[38][4] =  1.59261029e+000 ;
    fparams[38][5] =  2.21992482e-001 ;
    fparams[38][6] =  1.73263882e-001 ;
    fparams[38][7] =  2.01624958e-001 ;
    fparams[38][8] =  4.66280378e+000 ;
    fparams[38][9] =  2.53027803e+001 ;
    fparams[38][10] =  1.61265063e-003 ;
    fparams[38][11] =  1.53610568e-002 ;
    fparams[39][0] =  6.75302747e-001 ;
    fparams[39][1] =  6.54331847e-002 ;
    fparams[39][2] =  4.70286720e-001 ;
    fparams[39][3] =  1.06108709e+002 ;
    fparams[39][4] =  2.63497677e+000 ;
    fparams[39][5] =  2.06643540e+000 ;
    fparams[39][6] =  1.09621746e-001 ;
    fparams[39][7] =  1.93131925e-001 ;
    fparams[39][8] =  9.60348773e-001 ;
    fparams[39][9] =  1.63310938e+000 ;
    fparams[39][10] =  5.28921555e-003 ;
    fparams[39][11] =  1.66083821e-002 ;
    fparams[40][0] =  2.64365505e+000 ;
    fparams[40][1] =  2.20202699e+000 ;
    fparams[40][2] =  5.54225147e-001 ;
    fparams[40][3] =  1.78260107e+002 ;
    fparams[40][4] =  7.61376625e-001 ;
    fparams[40][5] =  7.67218745e-002 ;
    fparams[40][6] =  6.02946891e-003 ;
    fparams[40][7] =  1.55143296e-002 ;
    fparams[40][8] =  9.91630530e-002 ;
    fparams[40][9] =  1.76175995e-001 ;
    fparams[40][10] =  9.56782020e-001 ;
    fparams[40][11] =  1.54330682e+000 ;
    fparams[41][0] =  6.59532875e-001 ;
    fparams[41][1] =  8.66145490e-002 ;
    fparams[41][2] =  1.84545854e+000 ;
    fparams[41][3] =  5.94774398e+000 ;
    fparams[41][4] =  1.25584405e+000 ;
    fparams[41][5] =  6.40851475e-001 ;
    fparams[41][6] =  1.22253422e-001 ;
    fparams[41][7] =  1.66646050e-001 ;
    fparams[41][8] =  7.06638328e-001 ;
    fparams[41][9] =  1.62853268e+000 ;
    fparams[41][10] =  2.62381591e-003 ;
    fparams[41][11] =  8.26257859e-003 ;
    fparams[42][0] =  6.10160120e-001 ;
    fparams[42][1] =  9.11628054e-002 ;
    fparams[42][2] =  1.26544000e+000 ;
    fparams[42][3] =  5.06776025e-001 ;
    fparams[42][4] =  1.97428762e+000 ;
    fparams[42][5] =  5.89590381e+000 ;
    fparams[42][6] =  6.48028962e-001 ;
    fparams[42][7] =  1.46634108e+000 ;
    fparams[42][8] =  2.60380817e-003 ;
    fparams[42][9] =  7.84336311e-003 ;
    fparams[42][10] =  1.13887493e-001 ;
    fparams[42][11] =  1.55114340e-001 ;
    fparams[43][0] =  8.55189183e-001 ;
    fparams[43][1] =  1.02962151e-001 ;
    fparams[43][2] =  1.66219641e+000 ;
    fparams[43][3] =  7.64907000e+000 ;
    fparams[43][4] =  1.45575475e+000 ;
    fparams[43][5] =  1.01639987e+000 ;
    fparams[43][6] =  1.05445664e-001 ;
    fparams[43][7] =  1.42303338e-001 ;
    fparams[43][8] =  7.71657112e-001 ;
    fparams[43][9] =  1.34659349e+000 ;
    fparams[43][10] =  2.20992635e-003 ;
    fparams[43][11] =  7.90358976e-003 ;
    fparams[44][0] =  4.70847093e-001 ;
    fparams[44][1] =  9.33029874e-002 ;
    fparams[44][2] =  1.58180781e+000 ;
    fparams[44][3] =  4.52831347e-001 ;
    fparams[44][4] =  2.02419818e+000 ;
    fparams[44][5] =  7.11489023e+000 ;
    fparams[44][6] =  1.97036257e-003 ;
    fparams[44][7] =  7.56181595e-003 ;
    fparams[44][8] =  6.26912639e-001 ;
    fparams[44][9] =  1.25399858e+000 ;
    fparams[44][10] =  1.02641320e-001 ;
    fparams[44][11] =  1.33786087e-001 ;
    fparams[45][0] =  4.20051553e-001 ;
    fparams[45][1] =  9.38882628e-002 ;
    fparams[45][2] =  1.76266507e+000 ;
    fparams[45][3] =  4.64441687e-001 ;
    fparams[45][4] =  2.02735641e+000 ;
    fparams[45][5] =  8.19346046e+000 ;
    fparams[45][6] =  1.45487176e-003 ;
    fparams[45][7] =  7.82704517e-003 ;
    fparams[45][8] =  6.22809600e-001 ;
    fparams[45][9] =  1.17194153e+000 ;
    fparams[45][10] =  9.91529915e-002 ;
    fparams[45][11] =  1.24532839e-001 ;
    fparams[46][0] =  2.10475155e+000 ;
    fparams[46][1] =  8.68606470e+000 ;
    fparams[46][2] =  2.03884487e+000 ;
    fparams[46][3] =  3.78924449e-001 ;
    fparams[46][4] =  1.82067264e-001 ;
    fparams[46][5] =  1.42921634e-001 ;
    fparams[46][6] =  9.52040948e-002 ;
    fparams[46][7] =  1.17125900e-001 ;
    fparams[46][8] =  5.91445248e-001 ;
    fparams[46][9] =  1.07843808e+000 ;
    fparams[46][10] =  1.13328676e-003 ;
    fparams[46][11] =  7.80252092e-003 ;
    fparams[47][0] =  2.07981390e+000 ;
    fparams[47][1] =  9.92540297e+000 ;
    fparams[47][2] =  4.43170726e-001 ;
    fparams[47][3] =  1.04920104e-001 ;
    fparams[47][4] =  1.96515215e+000 ;
    fparams[47][5] =  6.40103839e-001 ;
    fparams[47][6] =  5.96130591e-001 ;
    fparams[47][7] =  8.89594790e-001 ;
    fparams[47][8] =  4.78016333e-001 ;
    fparams[47][9] =  1.98509407e+000 ;
    fparams[47][10] =  9.46458470e-002 ;
    fparams[47][11] =  1.12744464e-001 ;
    fparams[48][0] =  1.63657549e+000 ;
    fparams[48][1] =  1.24540381e+001 ;
    fparams[48][2] =  2.17927989e+000 ;
    fparams[48][3] =  1.45134660e+000 ;
    fparams[48][4] =  7.71300690e-001 ;
    fparams[48][5] =  1.26695757e-001 ;
    fparams[48][6] =  6.64193880e-001 ;
    fparams[48][7] =  7.77659202e-001 ;
    fparams[48][8] =  7.64563285e-001 ;
    fparams[48][9] =  1.66075210e+000 ;
    fparams[48][10] =  8.61126689e-002 ;
    fparams[48][11] =  1.05728357e-001 ;
    fparams[49][0] =  2.24820632e+000 ;
    fparams[49][1] =  1.51913507e+000 ;
    fparams[49][2] =  1.64706864e+000 ;
    fparams[49][3] =  1.30113424e+001 ;
    fparams[49][4] =  7.88679265e-001 ;
    fparams[49][5] =  1.06128184e-001 ;
    fparams[49][6] =  8.12579069e-002 ;
    fparams[49][7] =  9.94045620e-002 ;
    fparams[49][8] =  6.68280346e-001 ;
    fparams[49][9] =  1.49742063e+000 ;
    fparams[49][10] =  6.38467475e-001 ;
    fparams[49][11] =  7.18422635e-001 ;
    fparams[50][0] =  2.16644620e+000 ;
    fparams[50][1] =  1.13174909e+001 ;
    fparams[50][2] =  6.88691021e-001 ;
    fparams[50][3] =  1.10131285e-001 ;
    fparams[50][4] =  1.92431751e+000 ;
    fparams[50][5] =  6.74464853e-001 ;
    fparams[50][6] =  5.65359888e-001 ;
    fparams[50][7] =  7.33564610e-001 ;
    fparams[50][8] =  9.18683861e-001 ;
    fparams[50][9] =  1.02310312e+001 ;
    fparams[50][10] =  7.80542213e-002 ;
    fparams[50][11] =  9.31104308e-002 ;
    fparams[51][0] =  1.73662114e+000 ;
    fparams[51][1] =  8.84334719e-001 ;
    fparams[51][2] =  9.99871380e-001 ;
    fparams[51][3] =  1.38462121e-001 ;
    fparams[51][4] =  2.13972409e+000 ;
    fparams[51][5] =  1.19666432e+001 ;
    fparams[51][6] =  5.60566526e-001 ;
    fparams[51][7] =  6.72672880e-001 ;
    fparams[51][8] =  9.93772747e-001 ;
    fparams[51][9] =  8.72330411e+000 ;
    fparams[51][10] =  7.37374982e-002 ;
    fparams[51][11] =  8.78577715e-002 ;
    fparams[52][0] =  2.09383882e+000 ;
    fparams[52][1] =  1.26856869e+001 ;
    fparams[52][2] =  1.56940519e+000 ;
    fparams[52][3] =  1.21236537e+000 ;
    fparams[52][4] =  1.30941993e+000 ;
    fparams[52][5] =  1.66633292e-001 ;
    fparams[52][6] =  6.98067804e-002 ;
    fparams[52][7] =  8.30817576e-002 ;
    fparams[52][8] =  1.04969537e+000 ;
    fparams[52][9] =  7.43147857e+000 ;
    fparams[52][10] =  5.55594354e-001 ;
    fparams[52][11] =  6.17487676e-001 ;
    fparams[53][0] =  1.60186925e+000 ;
    fparams[53][1] =  1.95031538e-001 ;
    fparams[53][2] =  1.98510264e+000 ;
    fparams[53][3] =  1.36976183e+001 ;
    fparams[53][4] =  1.48226200e+000 ;
    fparams[53][5] =  1.80304795e+000 ;
    fparams[53][6] =  5.53807199e-001 ;
    fparams[53][7] =  5.67912340e-001 ;
    fparams[53][8] =  1.11728722e+000 ;
    fparams[53][9] =  6.40879878e+000 ;
    fparams[53][10] =  6.60720847e-002 ;
    fparams[53][11] =  7.86615429e-002 ;
    fparams[54][0] =  1.60015487e+000 ;
    fparams[54][1] =  2.92913354e+000 ;
    fparams[54][2] =  1.71644581e+000 ;
    fparams[54][3] =  1.55882990e+001 ;
    fparams[54][4] =  1.84968351e+000 ;
    fparams[54][5] =  2.22525983e-001 ;
    fparams[54][6] =  6.23813648e-002 ;
    fparams[54][7] =  7.45581223e-002 ;
    fparams[54][8] =  1.21387555e+000 ;
    fparams[54][9] =  5.56013271e+000 ;
    fparams[54][10] =  5.54051946e-001 ;
    fparams[54][11] =  5.21994521e-001 ;
    fparams[55][0] =  2.95236854e+000 ;
    fparams[55][1] =  6.01461952e+000 ;
    fparams[55][2] =  4.28105721e-001 ;
    fparams[55][3] =  4.64151246e+001 ;
    fparams[55][4] =  1.89599233e+000 ;
    fparams[55][5] =  1.80109756e-001 ;
    fparams[55][6] =  5.48012938e-002 ;
    fparams[55][7] =  7.12799633e-002 ;
    fparams[55][8] =  4.70838600e+000 ;
    fparams[55][9] =  4.56702799e+001 ;
    fparams[55][10] =  5.90356719e-001 ;
    fparams[55][11] =  4.70236310e-001 ;
    fparams[56][0] =  3.19434243e+000 ;
    fparams[56][1] =  9.27352241e+000 ;
    fparams[56][2] =  1.98289586e+000 ;
    fparams[56][3] =  2.28741632e-001 ;
    fparams[56][4] =  1.55121052e-001 ;
    fparams[56][5] =  3.82000231e-002 ;
    fparams[56][6] =  6.73222354e-002 ;
    fparams[56][7] =  7.30961745e-002 ;
    fparams[56][8] =  4.48474211e+000 ;
    fparams[56][9] =  2.95703565e+001 ;
    fparams[56][10] =  5.42674414e-001 ;
    fparams[56][11] =  4.08647015e-001 ;
    fparams[57][0] =  2.05036425e+000 ;
    fparams[57][1] =  2.20348417e-001 ;
    fparams[57][2] =  1.42114311e-001 ;
    fparams[57][3] =  3.96438056e-002 ;
    fparams[57][4] =  3.23538151e+000 ;
    fparams[57][5] =  9.56979169e+000 ;
    fparams[57][6] =  6.34683429e-002 ;
    fparams[57][7] =  6.92443091e-002 ;
    fparams[57][8] =  3.97960586e+000 ;
    fparams[57][9] =  2.53178406e+001 ;
    fparams[57][10] =  5.20116711e-001 ;
    fparams[57][11] =  3.83614098e-001 ;
    fparams[58][0] =  3.22990759e+000 ;
    fparams[58][1] =  9.94660135e+000 ;
    fparams[58][2] =  1.57618307e-001 ;
    fparams[58][3] =  4.15378676e-002 ;
    fparams[58][4] =  2.13477838e+000 ;
    fparams[58][5] =  2.40480572e-001 ;
    fparams[58][6] =  5.01907609e-001 ;
    fparams[58][7] =  3.66252019e-001 ;
    fparams[58][8] =  3.80889010e+000 ;
    fparams[58][9] =  2.43275968e+001 ;
    fparams[58][10] =  5.96625028e-002 ;
    fparams[58][11] =  6.59653503e-002 ;
    fparams[59][0] =  1.58189324e-001 ;
    fparams[59][1] =  3.91309056e-002 ;
    fparams[59][2] =  3.18141995e+000 ;
    fparams[59][3] =  1.04139545e+001 ;
    fparams[59][4] =  2.27622140e+000 ;
    fparams[59][5] =  2.81671757e-001 ;
    fparams[59][6] =  3.97705472e+000 ;
    fparams[59][7] =  2.61872978e+001 ;
    fparams[59][8] =  5.58448277e-002 ;
    fparams[59][9] =  6.30921695e-002 ;
    fparams[59][10] =  4.85207954e-001 ;
    fparams[59][11] =  3.54234369e-001 ;
    fparams[60][0] =  1.81379417e-001 ;
    fparams[60][1] =  4.37324793e-002 ;
    fparams[60][2] =  3.17616396e+000 ;
    fparams[60][3] =  1.07842572e+001 ;
    fparams[60][4] =  2.35221519e+000 ;
    fparams[60][5] =  3.05571833e-001 ;
    fparams[60][6] =  3.83125763e+000 ;
    fparams[60][7] =  2.54745408e+001 ;
    fparams[60][8] =  5.25889976e-002 ;
    fparams[60][9] =  6.02676073e-002 ;
    fparams[60][10] =  4.70090742e-001 ;
    fparams[60][11] =  3.39017003e-001 ;
    fparams[61][0] =  1.92986811e-001 ;
    fparams[61][1] =  4.37785970e-002 ;
    fparams[61][2] =  2.43756023e+000 ;
    fparams[61][3] =  3.29336996e-001 ;
    fparams[61][4] =  3.17248504e+000 ;
    fparams[61][5] =  1.11259996e+001 ;
    fparams[61][6] =  3.58105414e+000 ;
    fparams[61][7] =  2.46709586e+001 ;
    fparams[61][8] =  4.56529394e-001 ;
    fparams[61][9] =  3.24990282e-001 ;
    fparams[61][10] =  4.94812177e-002 ;
    fparams[61][11] =  5.76553100e-002 ;
    fparams[62][0] =  2.12002595e-001 ;
    fparams[62][1] =  4.57703608e-002 ;
    fparams[62][2] =  3.16891754e+000 ;
    fparams[62][3] =  1.14536599e+001 ;
    fparams[62][4] =  2.51503494e+000 ;
    fparams[62][5] =  3.55561054e-001 ;
    fparams[62][6] =  4.44080845e-001 ;
    fparams[62][7] =  3.11953363e-001 ;
    fparams[62][8] =  3.36742101e+000 ;
    fparams[62][9] =  2.40291435e+001 ;
    fparams[62][10] =  4.65652543e-002 ;
    fparams[62][11] =  5.52266819e-002 ;
    fparams[63][0] =  2.59355002e+000 ;
    fparams[63][1] =  3.82452612e-001 ;
    fparams[63][2] =  3.16557522e+000 ;
    fparams[63][3] =  1.17675155e+001 ;
    fparams[63][4] =  2.29402652e-001 ;
    fparams[63][5] =  4.76642249e-002 ;
    fparams[63][6] =  4.32257780e-001 ;
    fparams[63][7] =  2.99719833e-001 ;
    fparams[63][8] =  3.17261920e+000 ;
    fparams[63][9] =  2.34462738e+001 ;
    fparams[63][10] =  4.37958317e-002 ;
    fparams[63][11] =  5.29440680e-002 ;
    fparams[64][0] =  3.19144939e+000 ;
    fparams[64][1] =  1.20224655e+001 ;
    fparams[64][2] =  2.55766431e+000 ;
    fparams[64][3] =  4.08338876e-001 ;
    fparams[64][4] =  3.32681934e-001 ;
    fparams[64][5] =  5.85819814e-002 ;
    fparams[64][6] =  4.14243130e-002 ;
    fparams[64][7] =  5.06771477e-002 ;
    fparams[64][8] =  2.61036728e+000 ;
    fparams[64][9] =  1.99344244e+001 ;
    fparams[64][10] =  4.20526863e-001 ;
    fparams[64][11] =  2.85686240e-001 ;
    fparams[65][0] =  2.59407462e-001 ;
    fparams[65][1] =  5.04689354e-002 ;
    fparams[65][2] =  3.16177855e+000 ;
    fparams[65][3] =  1.23140183e+001 ;
    fparams[65][4] =  2.75095751e+000 ;
    fparams[65][5] =  4.38337626e-001 ;
    fparams[65][6] =  2.79247686e+000 ;
    fparams[65][7] =  2.23797309e+001 ;
    fparams[65][8] =  3.85931001e-002 ;
    fparams[65][9] =  4.87920992e-002 ;
    fparams[65][10] =  4.10881708e-001 ;
    fparams[65][11] =  2.77622892e-001 ;
    fparams[66][0] =  3.16055396e+000 ;
    fparams[66][1] =  1.25470414e+001 ;
    fparams[66][2] =  2.82751709e+000 ;
    fparams[66][3] =  4.67899094e-001 ;
    fparams[66][4] =  2.75140255e-001 ;
    fparams[66][5] =  5.23226982e-002 ;
    fparams[66][6] =  4.00967160e-001 ;
    fparams[66][7] =  2.67614884e-001 ;
    fparams[66][8] =  2.63110834e+000 ;
    fparams[66][9] =  2.19498166e+001 ;
    fparams[66][10] =  3.61333817e-002 ;
    fparams[66][11] =  4.68871497e-002 ;
    fparams[67][0] =  2.88642467e-001 ;
    fparams[67][1] =  5.40507687e-002 ;
    fparams[67][2] =  2.90567296e+000 ;
    fparams[67][3] =  4.97581077e-001 ;
    fparams[67][4] =  3.15960159e+000 ;
    fparams[67][5] =  1.27599505e+001 ;
    fparams[67][6] =  3.91280259e-001 ;
    fparams[67][7] =  2.58151831e-001 ;
    fparams[67][8] =  2.48596038e+000 ;
    fparams[67][9] =  2.15400972e+001 ;
    fparams[67][10] =  3.37664478e-002 ;
    fparams[67][11] =  4.50664323e-002 ;
    fparams[68][0] =  3.15573213e+000 ;
    fparams[68][1] =  1.29729009e+001 ;
    fparams[68][2] =  3.11519560e-001 ;
    fparams[68][3] =  5.81399387e-002 ;
    fparams[68][4] =  2.97722406e+000 ;
    fparams[68][5] =  5.31213394e-001 ;
    fparams[68][6] =  3.81563854e-001 ;
    fparams[68][7] =  2.49195776e-001 ;
    fparams[68][8] =  2.40247532e+000 ;
    fparams[68][9] =  2.13627616e+001 ;
    fparams[68][10] =  3.15224214e-002 ;
    fparams[68][11] =  4.33253257e-002 ;
    fparams[69][0] =  3.15591970e+000 ;
    fparams[69][1] =  1.31232407e+001 ;
    fparams[69][2] =  3.22544710e-001 ;
    fparams[69][3] =  5.97223323e-002 ;
    fparams[69][4] =  3.05569053e+000 ;
    fparams[69][5] =  5.61876773e-001 ;
    fparams[69][6] =  2.92845100e-002 ;
    fparams[69][7] =  4.16534255e-002 ;
    fparams[69][8] =  3.72487205e-001 ;
    fparams[69][9] =  2.40821967e-001 ;
    fparams[69][10] =  2.27833695e+000 ;
    fparams[69][11] =  2.10034185e+001 ;
    fparams[70][0] =  3.10794704e+000 ;
    fparams[70][1] =  6.06347847e-001 ;
    fparams[70][2] =  3.14091221e+000 ;
    fparams[70][3] =  1.33705269e+001 ;
    fparams[70][4] =  3.75660454e-001 ;
    fparams[70][5] =  7.29814740e-002 ;
    fparams[70][6] =  3.61901097e-001 ;
    fparams[70][7] =  2.32652051e-001 ;
    fparams[70][8] =  2.45409082e+000 ;
    fparams[70][9] =  2.12695209e+001 ;
    fparams[70][10] =  2.72383990e-002 ;
    fparams[70][11] =  3.99969597e-002 ;
    fparams[71][0] =  3.11446863e+000 ;
    fparams[71][1] =  1.38968881e+001 ;
    fparams[71][2] =  5.39634353e-001 ;
    fparams[71][3] =  8.91708508e-002 ;
    fparams[71][4] =  3.06460915e+000 ;
    fparams[71][5] =  6.79919563e-001 ;
    fparams[71][6] =  2.58563745e-002 ;
    fparams[71][7] =  3.82808522e-002 ;
    fparams[71][8] =  2.13983556e+000 ;
    fparams[71][9] =  1.80078788e+001 ;
    fparams[71][10] =  3.47788231e-001 ;
    fparams[71][11] =  2.22706591e-001 ;
    fparams[72][0] =  3.01166899e+000 ;
    fparams[72][1] =  7.10401889e-001 ;
    fparams[72][2] =  3.16284788e+000 ;
    fparams[72][3] =  1.38262192e+001 ;
    fparams[72][4] =  6.33421771e-001 ;
    fparams[72][5] =  9.48486572e-002 ;
    fparams[72][6] =  3.41417198e-001 ;
    fparams[72][7] =  2.14129678e-001 ;
    fparams[72][8] =  1.53566013e+000 ;
    fparams[72][9] =  1.55298698e+001 ;
    fparams[72][10] =  2.40723773e-002 ;
    fparams[72][11] =  3.67833690e-002 ;
    fparams[73][0] =  3.20236821e+000 ;
    fparams[73][1] =  1.38446369e+001 ;
    fparams[73][2] =  8.30098413e-001 ;
    fparams[73][3] =  1.18381581e-001 ;
    fparams[73][4] =  2.86552297e+000 ;
    fparams[73][5] =  7.66369118e-001 ;
    fparams[73][6] =  2.24813887e-002 ;
    fparams[73][7] =  3.52934622e-002 ;
    fparams[73][8] =  1.40165263e+000 ;
    fparams[73][9] =  1.46148877e+001 ;
    fparams[73][10] =  3.33740596e-001 ;
    fparams[73][11] =  2.05704486e-001 ;
    fparams[74][0] =  9.24906855e-001 ;
    fparams[74][1] =  1.28663377e-001 ;
    fparams[74][2] =  2.75554557e+000 ;
    fparams[74][3] =  7.65826479e-001 ;
    fparams[74][4] =  3.30440060e+000 ;
    fparams[74][5] =  1.34471170e+001 ;
    fparams[74][6] =  3.29973862e-001 ;
    fparams[74][7] =  1.98218895e-001 ;
    fparams[74][8] =  1.09916444e+000 ;
    fparams[74][9] =  1.35087534e+001 ;
    fparams[74][10] =  2.06498883e-002 ;
    fparams[74][11] =  3.38918459e-002 ;
    fparams[75][0] =  1.96952105e+000 ;
    fparams[75][1] =  4.98830620e+001 ;
    fparams[75][2] =  1.21726619e+000 ;
    fparams[75][3] =  1.33243809e-001 ;
    fparams[75][4] =  4.10391685e+000 ;
    fparams[75][5] =  1.84396916e+000 ;
    fparams[75][6] =  2.90791978e-002 ;
    fparams[75][7] =  2.84192813e-002 ;
    fparams[75][8] =  2.30696669e-001 ;
    fparams[75][9] =  1.90968784e-001 ;
    fparams[75][10] =  6.08840299e-001 ;
    fparams[75][11] =  1.37090356e+000 ;
    fparams[76][0] =  2.06385867e+000 ;
    fparams[76][1] =  4.05671697e+001 ;
    fparams[76][2] =  1.29603406e+000 ;
    fparams[76][3] =  1.46559047e-001 ;
    fparams[76][4] =  3.96920673e+000 ;
    fparams[76][5] =  1.82561596e+000 ;
    fparams[76][6] =  2.69835487e-002 ;
    fparams[76][7] =  2.84172045e-002 ;
    fparams[76][8] =  2.31083999e-001 ;
    fparams[76][9] =  1.79765184e-001 ;
    fparams[76][10] =  6.30466774e-001 ;
    fparams[76][11] =  1.38911543e+000 ;
    fparams[77][0] =  2.21522726e+000 ;
    fparams[77][1] =  3.24464090e+001 ;
    fparams[77][2] =  1.37573155e+000 ;
    fparams[77][3] =  1.60920048e-001 ;
    fparams[77][4] =  3.78244405e+000 ;
    fparams[77][5] =  1.78756553e+000 ;
    fparams[77][6] =  2.44643240e-002 ;
    fparams[77][7] =  2.82909938e-002 ;
    fparams[77][8] =  2.36932016e-001 ;
    fparams[77][9] =  1.70692368e-001 ;
    fparams[77][10] =  6.48471412e-001 ;
    fparams[77][11] =  1.37928390e+000 ;
    fparams[78][0] =  9.84697940e-001 ;
    fparams[78][1] =  1.60910839e-001 ;
    fparams[78][2] =  2.73987079e+000 ;
    fparams[78][3] =  7.18971667e-001 ;
    fparams[78][4] =  3.61696715e+000 ;
    fparams[78][5] =  1.29281016e+001 ;
    fparams[78][6] =  3.02885602e-001 ;
    fparams[78][7] =  1.70134854e-001 ;
    fparams[78][8] =  2.78370726e-001 ;
    fparams[78][9] =  1.49862703e+000 ;
    fparams[78][10] =  1.52124129e-002 ;
    fparams[78][11] =  2.83510822e-002 ;
    fparams[79][0] =  9.61263398e-001 ;
    fparams[79][1] =  1.70932277e-001 ;
    fparams[79][2] =  3.69581030e+000 ;
    fparams[79][3] =  1.29335319e+001 ;
    fparams[79][4] =  2.77567491e+000 ;
    fparams[79][5] =  6.89997070e-001 ;
    fparams[79][6] =  2.95414176e-001 ;
    fparams[79][7] =  1.63525510e-001 ;
    fparams[79][8] =  3.11475743e-001 ;
    fparams[79][9] =  1.39200901e+000 ;
    fparams[79][10] =  1.43237267e-002 ;
    fparams[79][11] =  2.71265337e-002 ;
    fparams[80][0] =  1.29200491e+000 ;
    fparams[80][1] =  1.83432865e-001 ;
    fparams[80][2] =  2.75161478e+000 ;
    fparams[80][3] =  9.42368371e-001 ;
    fparams[80][4] =  3.49387949e+000 ;
    fparams[80][5] =  1.46235654e+001 ;
    fparams[80][6] =  2.77304636e-001 ;
    fparams[80][7] =  1.55110144e-001 ;
    fparams[80][8] =  4.30232810e-001 ;
    fparams[80][9] =  1.28871670e+000 ;
    fparams[80][10] =  1.48294351e-002 ;
    fparams[80][11] =  2.61903834e-002 ;
    fparams[81][0] =  3.75964730e+000 ;
    fparams[81][1] =  1.35041513e+001 ;
    fparams[81][2] =  3.21195904e+000 ;
    fparams[81][3] =  6.66330993e-001 ;
    fparams[81][4] =  6.47767825e-001 ;
    fparams[81][5] =  9.22518234e-002 ;
    fparams[81][6] =  2.76123274e-001 ;
    fparams[81][7] =  1.50312897e-001 ;
    fparams[81][8] =  3.18838810e-001 ;
    fparams[81][9] =  1.12565588e+000 ;
    fparams[81][10] =  1.31668419e-002 ;
    fparams[81][11] =  2.48879842e-002 ;
    fparams[82][0] =  1.00795975e+000 ;
    fparams[82][1] =  1.17268427e-001 ;
    fparams[82][2] =  3.09796153e+000 ;
    fparams[82][3] =  8.80453235e-001 ;
    fparams[82][4] =  3.61296864e+000 ;
    fparams[82][5] =  1.47325812e+001 ;
    fparams[82][6] =  2.62401476e-001 ;
    fparams[82][7] =  1.43491014e-001 ;
    fparams[82][8] =  4.05621995e-001 ;
    fparams[82][9] =  1.04103506e+000 ;
    fparams[82][10] =  1.31812509e-002 ;
    fparams[82][11] =  2.39575415e-002 ;
    fparams[83][0] =  1.59826875e+000 ;
    fparams[83][1] =  1.56897471e-001 ;
    fparams[83][2] =  4.38233925e+000 ;
    fparams[83][3] =  2.47094692e+000 ;
    fparams[83][4] =  2.06074719e+000 ;
    fparams[83][5] =  5.72438972e+001 ;
    fparams[83][6] =  1.94426023e-001 ;
    fparams[83][7] =  1.32979109e-001 ;
    fparams[83][8] =  8.22704978e-001 ;
    fparams[83][9] =  9.56532528e-001 ;
    fparams[83][10] =  2.33226953e-002 ;
    fparams[83][11] =  2.23038435e-002 ;
    fparams[84][0] =  1.71463223e+000 ;
    fparams[84][1] =  9.79262841e+001 ;
    fparams[84][2] =  2.14115960e+000 ;
    fparams[84][3] =  2.10193717e-001 ;
    fparams[84][4] =  4.37512413e+000 ;
    fparams[84][5] =  3.66948812e+000 ;
    fparams[84][6] =  2.16216680e-002 ;
    fparams[84][7] =  1.98456144e-002 ;
    fparams[84][8] =  1.97843837e-001 ;
    fparams[84][9] =  1.33758807e-001 ;
    fparams[84][10] =  6.52047920e-001 ;
    fparams[84][11] =  7.80432104e-001 ;
    fparams[85][0] =  1.48047794e+000 ;
    fparams[85][1] =  1.25943919e+002 ;
    fparams[85][2] =  2.09174630e+000 ;
    fparams[85][3] =  1.83803008e-001 ;
    fparams[85][4] =  4.75246033e+000 ;
    fparams[85][5] =  4.19890596e+000 ;
    fparams[85][6] =  1.85643958e-002 ;
    fparams[85][7] =  1.81383503e-002 ;
    fparams[85][8] =  2.05859375e-001 ;
    fparams[85][9] =  1.33035404e-001 ;
    fparams[85][10] =  7.13540948e-001 ;
    fparams[85][11] =  7.03031938e-001 ;
    fparams[86][0] =  6.30022295e-001 ;
    fparams[86][1] =  1.40909762e-001 ;
    fparams[86][2] =  3.80962881e+000 ;
    fparams[86][3] =  3.08515540e+001 ;
    fparams[86][4] =  3.89756067e+000 ;
    fparams[86][5] =  6.51559763e-001 ;
    fparams[86][6] =  2.40755100e-001 ;
    fparams[86][7] =  1.08899672e-001 ;
    fparams[86][8] =  2.62868577e+000 ;
    fparams[86][9] =  6.42383261e+000 ;
    fparams[86][10] =  3.14285931e-002 ;
    fparams[86][11] =  2.42346699e-002 ;
    fparams[87][0] =  5.23288135e+000 ;
    fparams[87][1] =  8.60599536e+000 ;
    fparams[87][2] =  2.48604205e+000 ;
    fparams[87][3] =  3.04543982e-001 ;
    fparams[87][4] =  3.23431354e-001 ;
    fparams[87][5] =  3.87759096e-002 ;
    fparams[87][6] =  2.55403596e-001 ;
    fparams[87][7] =  1.28717724e-001 ;
    fparams[87][8] =  5.53607228e-001 ;
    fparams[87][9] =  5.36977452e-001 ;
    fparams[87][10] =  5.75278889e-003 ;
    fparams[87][11] =  1.29417790e-002 ;
    fparams[88][0] =  1.44192685e+000 ;
    fparams[88][1] =  1.18740873e-001 ;
    fparams[88][2] =  3.55291725e+000 ;
    fparams[88][3] =  1.01739750e+000 ;
    fparams[88][4] =  3.91259586e+000 ;
    fparams[88][5] =  6.31814783e+001 ;
    fparams[88][6] =  2.16173519e-001 ;
    fparams[88][7] =  9.55806441e-002 ;
    fparams[88][8] =  3.94191605e+000 ;
    fparams[88][9] =  3.50602732e+001 ;
    fparams[88][10] =  4.60422605e-002 ;
    fparams[88][11] =  2.20850385e-002 ;
    fparams[89][0] =  1.45864127e+000 ;
    fparams[89][1] =  1.07760494e-001 ;
    fparams[89][2] =  4.18945405e+000 ;
    fparams[89][3] =  8.89090649e+001 ;
    fparams[89][4] =  3.65866182e+000 ;
    fparams[89][5] =  1.05088931e+000 ;
    fparams[89][6] =  2.08479229e-001 ;
    fparams[89][7] =  9.09335557e-002 ;
    fparams[89][8] =  3.16528117e+000 ;
    fparams[89][9] =  3.13297788e+001 ;
    fparams[89][10] =  5.23892556e-002 ;
    fparams[89][11] =  2.08807697e-002 ;
    fparams[90][0] =  1.19014064e+000 ;
    fparams[90][1] =  7.73468729e-002 ;
    fparams[90][2] =  2.55380607e+000 ;
    fparams[90][3] =  6.59693681e-001 ;
    fparams[90][4] =  4.68110181e+000 ;
    fparams[90][5] =  1.28013896e+001 ;
    fparams[90][6] =  2.26121303e-001 ;
    fparams[90][7] =  1.08632194e-001 ;
    fparams[90][8] =  3.58250545e-001 ;
    fparams[90][9] =  4.56765664e-001 ;
    fparams[90][10] =  7.82263950e-003 ;
    fparams[90][11] =  1.62623474e-002 ;
    fparams[91][0] =  4.68537504e+000 ;
    fparams[91][1] =  1.44503632e+001 ;
    fparams[91][2] =  2.98413708e+000 ;
    fparams[91][3] =  5.56438592e-001 ;
    fparams[91][4] =  8.91988061e-001 ;
    fparams[91][5] =  6.69512914e-002 ;
    fparams[91][6] =  2.24825384e-001 ;
    fparams[91][7] =  1.03235396e-001 ;
    fparams[91][8] =  3.04444846e-001 ;
    fparams[91][9] =  4.27255647e-001 ;
    fparams[91][10] =  9.48162708e-003 ;
    fparams[91][11] =  1.77730611e-002 ;
    fparams[92][0] =  4.63343606e+000 ;
    fparams[92][1] =  1.63377267e+001 ;
    fparams[92][2] =  3.18157056e+000 ;
    fparams[92][3] =  5.69517868e-001 ;
    fparams[92][4] =  8.76455075e-001 ;
    fparams[92][5] =  6.88860012e-002 ;
    fparams[92][6] =  2.21685477e-001 ;
    fparams[92][7] =  9.84254550e-002 ;
    fparams[92][8] =  2.72917100e-001 ;
    fparams[92][9] =  4.09470917e-001 ;
    fparams[92][10] =  1.11737298e-002 ;
    fparams[92][11] =  1.86215410e-002 ;
    fparams[93][0] =  4.56773888e+000 ;
    fparams[93][1] =  1.90992795e+001 ;
    fparams[93][2] =  3.40325179e+000 ;
    fparams[93][3] =  5.90099634e-001 ;
    fparams[93][4] =  8.61841923e-001 ;
    fparams[93][5] =  7.03204851e-002 ;
    fparams[93][6] =  2.19728870e-001 ;
    fparams[93][7] =  9.36334280e-002 ;
    fparams[93][8] =  2.38176903e-001 ;
    fparams[93][9] =  3.93554882e-001 ;
    fparams[93][10] =  1.38306499e-002 ;
    fparams[93][11] =  1.94437286e-002 ;
    fparams[94][0] =  5.45671123e+000 ;
    fparams[94][1] =  1.01892720e+001 ;
    fparams[94][2] =  1.11687906e-001 ;
    fparams[94][3] =  3.98131313e-002 ;
    fparams[94][4] =  3.30260343e+000 ;
    fparams[94][5] =  3.14622212e-001 ;
    fparams[94][6] =  1.84568319e-001 ;
    fparams[94][7] =  1.04220860e-001 ;
    fparams[94][8] =  4.93644263e-001 ;
    fparams[94][9] =  4.63080540e-001 ;
    fparams[94][10] =  3.57484743e+000 ;
    fparams[94][11] =  2.19369542e+001 ;
    fparams[95][0] =  5.38321999e+000 ;
    fparams[95][1] =  1.07289857e+001 ;
    fparams[95][2] =  1.23343236e-001 ;
    fparams[95][3] =  4.15137806e-002 ;
    fparams[95][4] =  3.46469090e+000 ;
    fparams[95][5] =  3.39326208e-001 ;
    fparams[95][6] =  1.75437132e-001 ;
    fparams[95][7] =  9.98932346e-002 ;
    fparams[95][8] =  3.39800073e+000 ;
    fparams[95][9] =  2.11601535e+001 ;
    fparams[95][10] =  4.69459519e-001 ;
    fparams[95][11] =  4.51996970e-001 ;
    fparams[96][0] =  5.38402377e+000 ;
    fparams[96][1] =  1.11211419e+001 ;
    fparams[96][2] =  3.49861264e+000 ;
    fparams[96][3] =  3.56750210e-001 ;
    fparams[96][4] =  1.88039547e-001 ;
    fparams[96][5] =  5.39853583e-002 ;
    fparams[96][6] =  1.69143137e-001 ;
    fparams[96][7] =  9.60082633e-002 ;
    fparams[96][8] =  3.19595016e+000 ;
    fparams[96][9] =  1.80694389e+001 ;
    fparams[96][10] =  4.64393059e-001 ;
    fparams[96][11] =  4.36318197e-001 ;
    fparams[97][0] =  3.66090688e+000 ;
    fparams[97][1] =  3.84420906e-001 ;
    fparams[97][2] =  2.03054678e-001 ;
    fparams[97][3] =  5.48547131e-002 ;
    fparams[97][4] =  5.30697515e+000 ;
    fparams[97][5] =  1.17150262e+001 ;
    fparams[97][6] =  1.60934046e-001 ;
    fparams[97][7] =  9.21020329e-002 ;
    fparams[97][8] =  3.04808401e+000 ;
    fparams[97][9] =  1.73525367e+001 ;
    fparams[97][10] =  4.43610295e-001 ;
    fparams[97][11] =  4.27132359e-001 ;
    fparams[98][0] =  3.94150390e+000 ;
    fparams[98][1] =  4.18246722e-001 ;
    fparams[98][2] =  5.16915345e+000 ;
    fparams[98][3] =  1.25201788e+001 ;
    fparams[98][4] =  1.61941074e-001 ;
    fparams[98][5] =  4.81540117e-002 ;
    fparams[98][6] =  4.15299561e-001 ;
    fparams[98][7] =  4.24913856e-001 ;
    fparams[98][8] =  2.91761325e+000 ;
    fparams[98][9] =  1.90899693e+001 ;
    fparams[98][10] =  1.51474927e-001 ;
    fparams[98][11] =  8.81568925e-002 ;
    fparams[99][0] =  4.09780623e+000 ;
    fparams[99][1] =  4.46021145e-001 ;
    fparams[99][2] =  5.10079393e+000 ;
    fparams[99][3] =  1.31768613e+001 ;
    fparams[99][4] =  1.74617289e-001 ;
    fparams[99][5] =  5.02742829e-002 ;
    fparams[99][6] =  2.76774658e+000 ;
    fparams[99][7] =  1.84815393e+001 ;
    fparams[99][8] =  1.44496639e-001 ;
    fparams[99][9] =  8.46232592e-002 ;
    fparams[99][10] =  4.02772109e-001 ;
    fparams[99][11] =  4.17640100e-001 ;
    fparams[100][0] =  4.24934820e+000 ;
    fparams[100][1] =  4.75263933e-001 ;
    fparams[100][2] =  5.03556594e+000 ;
    fparams[100][3] =  1.38570834e+001 ;
    fparams[100][4] =  1.88920613e-001 ;
    fparams[100][5] =  5.26975158e-002 ;
    fparams[100][6] =  3.94356058e-001 ;
    fparams[100][7] =  4.11193751e-001 ;
    fparams[100][8] =  2.61213100e+000 ;
    fparams[100][9] =  1.78537905e+001 ;
    fparams[100][10] =  1.38001927e-001 ;
    fparams[100][11] =  8.12774434e-002 ;
    fparams[101][0] =  2.00942931e-001 ;
    fparams[101][1] =  5.48366518e-002 ;
    fparams[101][2] =  4.40119869e+000 ;
    fparams[101][3] =  5.04248434e-001 ;
    fparams[101][4] =  4.97250102e+000 ;
    fparams[101][5] =  1.45721366e+001 ;
    fparams[101][6] =  2.47530599e+000 ;
    fparams[101][7] =  1.72978308e+001 ;
    fparams[101][8] =  3.86883197e-001 ;
    fparams[101][9] =  4.05043898e-001 ;
    fparams[101][10] =  1.31936095e-001 ;
    fparams[101][11] =  7.80821071e-002 ;
    fparams[102][0] =  2.16052899e-001 ;
    fparams[102][1] =  5.83584058e-002 ;
    fparams[102][2] =  4.91106799e+000 ;
    fparams[102][3] =  1.53264212e+001 ;
    fparams[102][4] =  4.54862870e+000 ;
    fparams[102][5] =  5.34434760e-001 ;
    fparams[102][6] =  2.36114249e+000 ;
    fparams[102][7] =  1.68164803e+001 ;
    fparams[102][8] =  1.26277292e-001 ;
    fparams[102][9] =  7.50304633e-002 ;
    fparams[102][10] =  3.81364501e-001 ;
    fparams[102][11] =  3.99305852e-001 ;
    fparams[103][0] =  4.86738014e+000 ;
    fparams[103][1] =  1.60320520e+001 ;
    fparams[103][2] =  3.19974401e-001 ;
    fparams[103][3] =  6.70871138e-002 ;
    fparams[103][4] =  4.58872425e+000 ;
    fparams[103][5] =  5.77039373e-001 ;
    fparams[103][6] =  1.21482448e-001 ;
    fparams[103][7] =  7.22275899e-002 ;
    fparams[103][8] =  2.31639872e+000 ;
    fparams[103][9] =  1.41279737e+001 ;
    fparams[103][10] =  3.79258137e-001 ;
    fparams[103][11] =  3.89973484e-001 ;

   feTableRead = 1; /* remember that table has been read */
   return( (NZMAX-NZMIN+1)*NPMAX );

} /* end ReadfeTable() */

/*--------------------- ReadXYZcoord() -----------------------*/
/*
    read a set of (x,y,z) coordinates from a file
    and return number of coord read
    
    infile = name of input file to read from
    ncellx,y,z = number of unit cells in x,y,z to replicate
            (must be >= 1 )
    ax, by, cz = will get unit cell size from file
    x,y,z  = pointer to a pointer to get array
            of (x,y,z) coord
    occ    = pointer to a pointer to get array
            of occupancy
    wobble = rms thermal displacement (in Angstroms)
            at Temperature = 300 degrees K
    Znum  = pointer to a pointer to get array
            atomic numbers Z
    line1 = string to get 1st line of file
        with description

   NOTE: x,y,z,occ and Znum array are allocated by this routine
    because it is not known ahead of time home many points
    will be read in (i.e. this routine figures it out)

    only infile is unchanged by this routine

    convert to string and streams mar-2014 ejk
*/
int ReadXYZcoord( const char* infile, const int ncellx, const int ncelly,
    const int ncellz, float *ax, float *by, float *cz, vectori &Znum, 
    vectorf &x, vectorf &y, vectorf &z, vectorf &occ, vectorf &wobble,
    string &line1 )
{
   int i, j, iz, ncoord, na, ntotal;
   string buf, stemp;
   ifstream fp;

   vector<int> Znumv;
   vector<float> xv, yv, zv, occv, wobblev;

    //  abort if ncellx,y,z not valid
    if( (ncellx<1) || (ncelly<1) || (ncellz<1) ) {
        stemp= "invalid ncellx,y,z = "+toString(ncellx)+", "+toString(ncelly)+
            ", " + toString(ncellz) + " in ReadXYZcoord()",
        messageSL( stemp.c_str(), 0 );
        return( 0 );
    }

    fp.open( infile );
    if( fp.fail() ) {
        stemp= "ReadXYZcoord() cannot open file " + std::string( infile );
        messageSL( stemp.c_str(), 0 );
        return( 0 );
    }

    getline( fp, line1 );   // read file description

    // read unit cell size
    *cz = 0.0F;  // init to 0 because cz may be blank in file
    getline( fp, buf );
    istringstream isbuf( buf );
    isbuf >> *ax >> *by >> *cz;

    // read into a dynamically sized vector to avoid an extra pass through the
    //    file just to count number of atomic coordinates
    ncoord = 0;
    do {
        getline( fp, buf );
        float xx, yy, zz, oo, ww;
        istringstream isbuf( buf );
        isbuf >> iz ;
        if( (iz > 0) && (iz <= NZMAX ) ) {
                ncoord += 1;
                Znumv.push_back( iz );
                ww = 0.0F;
                isbuf >>  xx >> yy >> zz >> oo >> ww;
                xv.push_back( xx );
                yv.push_back( yy );
                zv.push_back( zz );
                occv.push_back( oo );
                wobblev.push_back( ww );
                
        } else if( iz > NZMAX ) {
                stemp= "Warning bad atomic number " + toString( iz ) + " in ReadXYZcoord()";
                messageSL( stemp.c_str() );
        } 
    } while ( iz > 0 );
    fp.close( );
    
    /* now that we know how many coordinates there are
        allocate the arrays */

    ntotal  = ncoord * ncellx * ncelly * ncellz;
    x.resize( ntotal );
    y.resize( ntotal );
    z.resize( ntotal );
    occ.resize( ntotal );
    wobble.resize( ntotal );
    Znum.resize( ntotal );

    //  save what was just read into final arrays
    for( j=0; j< ncoord; j++ ) {
        x[j] = xv[j];
        y[j] = yv[j];
        z[j] = zv[j];
        occ[j] = occv[j];
        wobble[j] = wobblev[j];
        Znum[j] = Znumv[j];
    }
    na = ncoord;

    //   replicate unit cell if needed
    if( (ncellx > 1) || (ncelly > 1) || (ncellz > 1 ) ) {

        if( ncellx > 1 ) {
            for( i=1; i<ncellx; i++)
            for( j=0; j<na; j++) {
                x[j + na*i]  = x[j] + i*(*ax);
                y[j + na*i]  = y[j];
                z[j + na*i]  = z[j];
                occ[j + na*i]     = occ[j];
                wobble[j + na*i]  = wobble[j];
                Znum[j + na*i]    = Znum[j];
            }
            na = na * ncellx;
            *ax = (*ax) * ncellx;
        }

        if( ncelly > 1 ) {
            for( i=1; i<ncelly; i++)
            for( j=0; j<na; j++) {
                x[j + na*i]  = x[j];
                y[j + na*i]  = y[j] + i*(*by);
                z[j + na*i]  = z[j];
                occ[j + na*i]     = occ[j];
                wobble[j + na*i]  = wobble[j];
                Znum[j + na*i]    = Znum[j];
            }
            na = na * ncelly;
            *by = (*by) * ncelly;
        }

        if( ncellz > 1 ) {
            for( i=1; i<ncellz; i++)
            for( j=0; j<na; j++) {
                 x[j + na*i]  = x[j];
                y[j + na*i]  = y[j];
                z[j + na*i]  = z[j] + i*(*cz);
                occ[j + na*i]     = occ[j];
                wobble[j + na*i]  = wobble[j];
                Znum[j + na*i]    = Znum[j];
            }
            na = na * ncellz;
            *cz = (*cz) * ncellz;
        }

    }  /* end if( ncellx,y,z > 1 ) */

    return( na );

}  /* end ReadXYZcoord() */

/*----------------------- seval() ----------------------*/
/*
    Interpolate from cubic spline coefficients

    E. Kirkland 4-JUL-85
    modified to do a binary search for efficiency 13-Oct-1994 ejk
    converted to C 26-jun-1995 ejk
    fixed problem on end-of-range 16-July-1995 ejk

    The inputs are:
        x[n] = array of x values in ascending order, each x[i] must
            be unique
        y[n] = array of y values corresponding to x[n]
        b[n] = array of spline coefficients for (x-x[i])
        c[n] = array of spline coefficients for (x-x[i])**2
        d[n] = array of spline coefficients for (x-x[i])**3
        n  = number of data points
        x0  = the x value to interpolate at
        (x[i] <= x <= x[i+1]) and all inputs remain unchanged

    The value returned is the interpolated y value.

    The coefficients b[i], c[i], d[i] refer to the x[i] to x[i+1]
    interval. NOTE that the last set of coefficients,
    b[n-1], c[n-1], d[n-1] are meaningless.
*/
double seval( vectord &x, vectord &y, vectord &b, vectord &c,
         vectord &d, int n, double x0 )
{
    int i, j, k;
    double z, seval1;

    /*  exit if x0 is outside the spline range */
    if( x0 <= x[0] ) i = 0;
    else if( x0 >= x[n-2] ) i = n-2;
    else { 
        i = 0;
        j = n;
        do{ k = ( i + j ) / 2 ;
            if( x0 < x[k] )  j = k;
            else if( x0 >= x[k] ) i = k;
        } while ( (j-i) > 1 );
    }
    
    z = x0 - x[i];
    seval1 = y[i] + ( b[i] + ( c[i] + d[i] *z ) *z) *z;

    return( seval1 );

} /* end seval() */

/*--------------------- sigma() -----------------------------------*/
/*
    return the interaction parameter sigma in radians/(kv-Angstroms)
    keep this is one place so I don't have to keep typing in these
    constants (that I can never remember anyhow)

    ref: Physics Vade Mecum, 2nd edit, edit. H. L. Anderson
        (The American Institute of Physics, New York) 1989
        page 4.

    kev = electron energy in keV

*/

double sigma( double kev )
{
    double s, pi, wavl, x;
    const double emass=510.99906; /* electron rest mass in keV */
    double wavelength( double kev );  /*  get electron wavelength */

    x = ( emass + kev ) / ( 2.0*emass + kev);
    wavl = wavelength( kev );
    pi = 4.0 * atan( 1.0 );
    
    s = 2.0 * pi * x / (wavl*kev);
    
    return( s );

}  /* end sigma() */

/*----------------- sortByZ() ------------------------------

    improved Shell sort modeled after prog. 6.5 (pg. 274) of
    R. Sedgewick, "Algorithms in C", 3rd edit. Addison-Wesley 1998
    
    x[], y[], z[]   = atom coordinates 
    occ[]           = occupancy of each atom
    Znum[]          = atomic number of each atom
    natom           = number of atoms
*/
void sortByZ( vectorf &x, vectorf &y, vectorf &z, vectorf &occ,
    vectori &Znum, int natom )
{
    int i, j, h, Znum2;
    float x2, y2, z2, occ2;

    for( h=1; h<=(natom-1)/9; h=3*h+1);

    for( ; h>0; h/=3 )
        for( i=h; i<natom; i++) {
            j = i;
            x2 = x[i];
            y2 = y[i];
            z2 = z[i];
            occ2 = occ[i];
            Znum2 = Znum[i];
            while( (j >= h) && ( z2 < z[j-h]) ) {
                x[j] = x[j-h];
                y[j] = y[j-h];
                z[j] = z[j-h];
                occ[j] = occ[j-h];
                Znum[j] = Znum[j-h];
                j -= h;
            }
            x[j] = x2;
            y[j] = y2;
            z[j] = z2;
            occ[j] = occ2;
            Znum[j] = Znum2;
        }

/*  Test sort routine -- DELETE this after awhile  */
    char mbad[]="Bad sort !";   //  gcc requires this step
    for( i=1; i<natom; i++) {
        if( z[i-1] > z[i] ) messageSL( mbad );
    }
 
}  /* end sortByZ() */

/*------------------ splinh() -----------------------------*/
/*
    fit a quasi-Hermite  cubic spline
    
    [1] Spline fit as in H.Akima, J. ACM 17(1970)p.589-602
        'A New Method of Interpolation and Smooth
        Curve Fitting Based on Local Procedures'

    [2] H.Akima, Comm. ACM, 15(1972)p.914-918

    E. Kirkland 4-JUL-85
    changed zero test to be a small nonzero number 8-jul-85 ejk
    converted to C 24-jun-1995 ejk
    convert to vector<double> 28-jun-2016 ejk

    The inputs are:
        x[n] = array of x values in ascending order, each X(I) must
            be unique
        y[n] = array of y values corresponding to X(N)
        n  = number of data points must be 2 or greater

    The outputs are (with z=x-x(i)):
        b[n] = array of spline coefficients for (x-x[i])
        c[n] = array of spline coefficients for (x-x[i])**2
        d[n] = array of spline coefficients for (x-x[i])**3
        ( x[i] <= x <= x[i+1] )
    To interpolate y(x) = yi + bi*z + c*z*z + d*z*z*z

    The coefficients b[i], c[i], d[i] refer to the x[i] to x[i+1]
    interval. NOTE that the last set of coefficients,
    b[n-1], c[n-1], d[n-1] are meaningless.
*/
void splinh( vectord &x, vectord &y,
         vectord &b, vectord &c, vectord &d, int n)
{
#define SMALL 1.0e-25

    int i, nm1, nm4;
    double m1, m2, m3, m4, m5, t1, t2, m54, m43, m32, m21, x43;

    if( n < 4) return;

    /* Do the first end point (special case),
       and get starting values */

    m5 = ( y[3] - y[2] ) / ( x[3] - x[2] ); /* mx = slope at pt x */
    m4 = ( y[2] - y[1] ) / ( x[2] - x[1] );
    m3 = ( y[1] - y[0] ) / ( x[1] - x[0] );

    m2 = m3 + m3 - m4;  /* eq. (9) of reference [1] */
    m1 = m2 + m2 - m3;

    m54 = fabs( m5 - m4);
    m43 = fabs( m4 - m3);
    m32 = fabs( m3 - m2);
    m21 = fabs( m2 - m1);

    if ( (m43+m21) > SMALL )
        t1 = ( m43*m2 + m21*m3 ) / ( m43 + m21 );
    else
        t1 = 0.5 * ( m2 + m3 );

    /*  Do everything up to the last end points */

    nm1 = n-1;
    nm4 = n-4;

    for( i=0; i<nm1; i++) {

        if( (m54+m32) > SMALL )
            t2= (m54*m3 + m32*m4) / (m54 + m32);
        else 
            t2 = 0.5* ( m3 + m4 );
      
        x43 = x[i+1] - x[i];
        b[i] = t1;
        c[i] = ( 3.0*m3 - t1 - t1 - t2 ) /x43;
        d[i] = ( t1 + t2 - m3 - m3 ) / ( x43*x43 );

        m1 = m2;
        m2 = m3;
        m3 = m4;
        m4 = m5;
        if( i < nm4 ) {
            m5 = ( y[i+4] - y[i+3] ) / ( x[i+4] - x[i+3] );
        } else {
            m5 = m4 + m4 - m3;
        }

        m21 = m32;
        m32 = m43;
        m43 = m54;
        m54 = fabs( m5 - m4 );
        t1 = t2;
    }

    return;

} /* end splinh() */

/*------------------------- toString( int ) ----------------------*/
/*
    convert a number into a string
*/
std::string toString( int i )
{
        std::stringstream ss;
        ss << i;
        return ss.str();

}  // end toString( int )

/*------------------------- toString( float ) ----------------------*/
/*
    convert a number into a string
*/
std::string toString( float x )
{
        std::stringstream ss;
        ss << x;
        return ss.str();

}  // end toString( float)

/*------------------------- toString( double ) ----------------------*/
/*
    convert a number into a string
*/
std::string toString( double x )
{
        std::stringstream ss;
        ss << x;
        return ss.str();

}  // end toString( double )


/*------------------------ transmit() ------------------------*/
/*
    transmit the wavefunction thru one layer

    waver,i[ix][iy]  = real and imaginary parts of wavefunction
    transr,i[ix][iy] = real and imag parts of transmission functions

    nx, ny = size of array
    
    on entrance waver,i and transr,i are in real space
    
    only waver,i will be changed by this routine
*/
void transmit( float** waver, float** wavei,
               float** transr, float** transi,
                int nx, int ny )
{
    int ix, iy;
    float wr, wi, tr, ti;
    
       for( ix=0; ix<nx; ix++) {
            for( iy=0; iy<ny; iy++) {
               wr = waver[ix][iy];
               wi = wavei[ix][iy];
               tr = transr[ix][iy];
               ti = transi[ix][iy];
               waver[ix][iy] = wr*tr - wi*ti;
               wavei[ix][iy] = wr*ti + wi*tr;
            } /* end for(iy...) */
       }  /* end for(ix...) */

} /* end transmit() */

/*--------------------- vatom() -----------------------------------*/
/*
    return the real space atomic potential (NOT projected)
    in volts for atomic number Z at radius r

    Z = atomic number 1 <= Z <= 103
    radius  = radius in Angstroms (MUST be > 0)

  assumed global vars:

#define NZMIN   1   = min Z in featom.tab 
#define NZMAX   103 = max Z in featom.tab 

int feTableRead=0; = flag to remember if the param file has been read 
int nl=3, ng=3; = number of Lorenzians and Gaussians 
double fparams[][] = fe parameters

  al and ag calculated using physical constants from:
    H. L. Anderson, editor "A Physicist's Desk Reference",
        2nd edition, Amer. Instit. Physics, 1989

  started from vzatom() 24-nov-1997 ejk
*/

double vatom( int Z, double radius )
{
   int i, nfe;
   double suml, sumg, x,t, r;

   /* Lorenzian, Gaussian constants */
   const double al=150.4121417, ag=266.5985798;
   const double pi=3.141592654;

   if( (Z<NZMIN) || (Z>NZMAX) ) return( 0.0 );

   /* read in the table from a file if this is the
    first time this is called */
   if( feTableRead == 0 ) nfe = ReadfeTable();

   r = fabs( radius );
   if( r < 1.0e-10 ) r = 1.0e-10;  /* avoid singularity at r=0 */
   suml = sumg = 0.0;

   /* Lorenztians */
   x = 2.0*pi*r;
   for( i=0; i<2*nl; i+=2 )
    suml += fparams[Z][i]* exp( -x*sqrt(fparams[Z][i+1]) );

   /* Gaussians */
   x = pi*r;
   x = x*x;
   for( i=2*nl; i<2*(nl+ng); i+=2 ) {
        t = sqrt( fparams[Z][i+1] );
        t = t*t*t;
        sumg += fparams[Z][i]*exp(-x/fparams[Z][i+1]) / t;
   }

   return( al*suml/r + ag*sumg );

}  /* end vatom() */

/*--------------------- vzatom() -----------------------------------*/
/*
    return the real space projected atomic potential
    in volt-Angstroms for atomic number Z at radius r

    Z = atomic number 1 <= Z <= 103
    radius  = radius in Angstroms (MUST be > 0)

  assumed global vars:

#define NZMIN   1   = min Z in featom.tab 
#define NZMAX   103 = max Z in featom.tab

int feTableRead=0; = flag to remember if the param file has been read 
int nl=3, ng=3; = number of Lorenzians and Gaussians 
double fparams[][] = fe parameters

  al and ag calculated using physical constants from:
    H. L. Anderson, editor "A Physicist's Desk Reference",
        2nd edition, Amer. Instit. Physics, 1989
*/

double vzatom( int Z, double radius )
{
   int i, nfe;
   double suml, sumg, x, r;

   /* Lorenzian, Gaussian constants */
   const double al=300.8242834, ag=150.4121417;
   const double pi=3.141592654;

   if( (Z<NZMIN) || (Z>NZMAX) ) return( 0.0 );

   /* read in the table from a file if this is the
    first time this is called */
   if( feTableRead == 0 ) nfe = ReadfeTable();

   r = fabs( radius );
   if( r < 1.0e-10 ) r = 1.0e-10;  /* avoid singularity at r=0 */
   suml = sumg = 0.0;

   /* Lorenztians */
   x = 2.0*pi*r;
   for( i=0; i<2*nl; i+=2 )
        suml += fparams[Z][i]* bessk0( x*sqrt(fparams[Z][i+1]) );

   /* Gaussians */
   x = pi*r;
   x = x*x;
   for( i=2*nl; i<2*(nl+ng); i+=2 )
        sumg += fparams[Z][i]*exp(-x/fparams[Z][i+1]) / fparams[Z][i+1];

   return( al*suml + ag*sumg );

}  /* end vzatom() */


/*--------------------- vzatomLUT() -----------------------------------*/
/*
    return the (real space) projected atomic potential for atomic
    number Z at radius r (in Angstroms)

    this mimics vzatom() in slicelib.c but uses a look-up-table
    with cubic spline interpolation to make it run about 2X-4X faster

    started 23-may-1997 E. Kirkland
    fix Z range to allow Hydrogen 1-jan-1998 ejk
    switch to r^2 to avoid a lot of sqrt() and reduce CPU time
             6-may-2008 ejk
    convert to vector<double> 28-jun-2016 ejk

    Z = atomic number 1 <= Z <= 98
    rsq = square of (radius in Angstroms)
*/

double vzatomLUT( int Z, double rsq )
{
    std::string stemp;
    int i, iz;
    double dlnr, vz, r;
    const static double RMIN= 0.01;    /* r (in Ang) range of LUT for vzatomLUT() */
    const static double RMAX= 5.0;
    const static int NRMAX= 100;       /* number of in look-up-table in vzatomLUT */

    /* spline interpolation coeff. */
    static int splineInit=0;
    static vector<int> nspline( NZMAX );

    static vectord splinx( NRMAX, 0 );
    static vector<vectord > spliny( NZMAX, splinx ); 
    static vector<vectord > splinb( NZMAX, splinx ); 
    static vector<vectord > splinc( NZMAX, splinx ); 
    static vector<vectord > splind( NZMAX, splinx ); 

    if( splineInit == 0 ) {
 
        /*  generate a set of logarithmic r values */
        dlnr = log(RMAX/RMIN)/(NRMAX-1);
        for( i=0; i<NRMAX; i++)
            splinx[i] = RMIN * exp( i * dlnr );
        stemp= "fit from r= " + toString(splinx[0]) + " to r= " + toString(splinx[NRMAX-1]);
        messageSL( stemp.c_str(), 0 );
        for( i=0; i<NRMAX; i++)
        splinx[i] = splinx[i] * splinx[i];   /* use r^2 not r for speed */

        for( i=0; i<NZMAX; i++) nspline[i] = 0;
        splineInit = 1;     /* remember that this has been done */
    }

    iz = Z - 1;  /* convert atomic number to array index */
    if( (Z < 1) || ( Z > NZMAX) ) {
        stemp= "Bad atomic number " + toString(Z) + " in vzatomLUT()";
        messageSL( stemp.c_str(), 2 );
        exit( 0 );
    }

    /* if this atomic number has not been called before
                generate the spline coefficients */
    if( nspline[iz] == 0 ) {

    for( i=0; i<NRMAX; i++) {
        r = sqrt( splinx[i] );
        spliny[iz][i] = vzatom( Z, r );
    }
        nspline[iz] = NRMAX;
        splinh( splinx, spliny[iz], splinb[iz],
                  splinc[iz], splind[iz], NRMAX);
    }

    /* now that everything is set up find the
                scattering factor by interpolation in the table 
    */

    vz = seval( splinx, spliny[iz], splinb[iz],
                splinc[iz], splind[iz], nspline[iz], rsq );

    return( vz );

}  /* end vzatomLUT() */

/*--------------------- wavelength() -----------------------------------*/
/*
    return the electron wavelength (in Angstroms)
    keep this is one place so I don't have to keep typing in these
    constants (that I can never remember anyhow)

    ref: Physics Vade Mecum, 2nd edit, edit. H. L. Anderson
        (The American Institute of Physics, New York) 1989
        page 4.

    kev = electron energy in keV

*/

double wavelength( double kev )
{
    double w;
    const double emass=510.99906; /* electron rest mass in keV */
    const double hc=12.3984244; /* Planck's const x speed of light*/

    /* electron wavelength in Angstroms */
    w = hc/sqrt( kev * ( 2*emass + kev ) );

    return( w );

}  /* end wavelength() */
