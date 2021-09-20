/*      *** slicelib.hpp ***

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
    ranflat()     : return a random number with uniform distribution
    rangauss()    : return a random number with Gaussian distribution
    ranPoisson()  : return a random number with Poisson distribution
    readCnm()     : decypher aberr. of the form C34a etc.
    ReadfeTable() : read fe scattering factor table
    ReadXYZcoord(): read a set of (x,y,z) coordinates from a file
    seval()       : Interpolate from cubic spline coefficients
    sigma()       : return the interaction parameter
    sortByZ()     : sort atomic x,y,z coord. by z
    splinh()      : fit a quasi-Hermite  cubic spline
    toString()    : convert numbers (int,float,double)to string
    ToInt()       : convert float to int (round properly)
    transmit()    : transmit a 2D wavefunction
    vatom()       : return real space atomic potential (NOT projected)
    vzatom()      : return real space projected atomic potential
    vzatomLUT()   : same as vzatom() but with a look-up-table (faster)
    wavelength()  : return electron wavelength in A for given keV


    this file is formatted for a tab size of 8 characters

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
    add nowarranty() 27-jan-1999 ejk
    put in new scattering factor parameters 9-oct-1999 ejk
    move malloc1D() and malloc2D() to here 6-nov-1999 ejk
    rearrange nowarranty 22-jan-2000 ejk
    add malloc3D(), free2D() and free3D()  20-jul-2005 ejk
    move sortbyz() to here 23-aug-2006 ejk
    remove nowarranty() subroutine from this file 11-jul-2007 ejk
    fix a few small issues (scanf field width etc.) 16-mar-2008 ejk
    convert format to 4 char TAB size and add vzatomLUT() 11-jun-2008 ejk
    add subroutines for use with FFTW 18-mar-2010 ejk
    put conditional compilation around FFTW portion and add parameter
       offsets for multipole aberrations 12-apr-2011 ejk
    add readCnm() and chi() 8-may-2011 ejk
    remove explicit FFTW reference, scaleW(), transportW(), propagateW(),
        add cfpix class in propagate() 4-nov-2012 ejk
    add ToInt() and start putting in the calculation parameters
       such as nx, ny etc. to capture the whole calc. 
       and change #define to const int, and add pSOURCE
         13-apr-2013 ejk
    add printfOK conditional to exclude code that won't work in a GUI
         14-apr-2013 ejk
    add pTEMPER, pNWOBBLE 2-jun-2013,  pDELTAZ 5-jun-2013 ejk
    add pOAPMIN 17-aug-2013 ejk
    add pCZ+ confocal24-aug-2013, pMODE 14-oct-2013 ejk
   convert char[] and sprintf() to std::string  and
       add 3 versions of toString() 24-nov-2013 ejk
   change printfOK to USE_TERMINAL and
        start conversion to streams 10-mar-2014 ejk
   remove ReadLine() and use C++ getline() 4-apr-2014 ejk
   add enum for MODE's 27-sep-2014 ejk
   add ranPoisson() 27-sep-2015 ejk
   add pPROBEI, pPROBET, pDOSE  4-oct-2015 ejk
   convert malloc1D() etc to vector<> 28-jun-2016 ejk
   remove malloc2D(), malloc3D(), free2D(), free3D() and
      make separate newnd.hpp with new/del template 22-sep-2017 ejk
   fix typo in comments for ranPoisson() 29-nov-2017 ejk
   update comments in ReadfeTable( ) and remove redundant declarations
        of ReadfeTable( )  24-apr-2018 ejk
   add pPMAXDET,pPMINDET params for azim. STEM detect. angles (seg. detector)
       12-may-2018 ejk
   fix definition of seval() 17-feb-2019 ejk
   remove propagate() so slicelib is not dependent on cfpix+fftw 29-jul-2019 ejk
*/

#ifndef SLICELIB_HPP   // only include this file if its not already

#define SLICELIB_HPP   // remember that this has been included

#include <string>       // STD string class
#include <vector>   // STD vector class

using namespace std;

#define USE_TERMINAL   //  define to allow code to use printf (not allowed inside a GUI)

// shorthand data types to save some typing
typedef vector<double> vectord;
typedef vector<float> vectorf;
typedef vector<int> vectori;

// define symbols for the parameter offsets
//  (could use an enum here too)
// some of these are redundant because they can be calculated from other parameters
// this has evolved over a long time with params added as needed, 
//             so the order is not ideal - sorry
const int
        pNPIX    =   0,   // number of pix 1 for real and 2 for complex
        pRMAX    =   1,   // maximum value of the real part of the image
        pIMAX    =   2,   // maximum value of the imaginary part of the image
        pRMIN    =   3,   // minimum value of the real part of the image
        pIMIN    =   4,   // minimum value of the imaginary part of the image
        pXBTILT  =   5,   // x beam tilt in rad
        pYBTILT  =   6,   // y beam tilt in rad
        pC       =   7,   // c unit cell dimension in Angstroms
        pRES     =   8,   // real space resolution in atompot
        pXCTILT  =   9,   // x crystal tilt in rad
        pYCTILT  =   10,  // y crystal tilt in rad
        pDEFOCUS =   11,  // defocus in Angstroms
        pASTIG   =   12,  // astigmatism in Angstroms
        pTHETA   =   13,  // angle of astigmatism in radians
        pDX      =   14,  // dimension of pixel in x direction in Angstroms
        pDY      =   15,  // dimension of pixel in y direction in Angstroms
        pENERGY  =   16,  // beam energy in keV
        pOAPERT  =   17,  // objective aperture semi-angle in radians
                          //     also see pOAPMIN below
        pCS      =   18,  // spherical aberration in Angstroms
        pWAVEL   =   19,  // electron wavelength in Angstroms
        pCAPERT  =   21,  // condenser (CTEM) illumination angle in radians
                          //    see also pCAPERTMIN below
        pDDF     =   22,  // defocus spread in Angstroms
        pNSLICES =   29,  // number of slices
        pDELTAZ  =   30,  // slice thickness in Angstroms
        pMINDET =    31,  // minimum detector angle (STEM) in radians
        pMAXDET =    32,  // maximum detector angle (STEM) in radians
        pPMAXDET =   33,  // max. detector azimuthal angle (phi, STEM) in radians
        pPMINDET =   34,  // min. detector azimuthal angle (phi, STEM) in radians

//      pC10   =     35,  //C10 = -df so don't use 
//      pC30   =     18,  //3rd order spherical - reuse from above
//      pC50   =     50,  //5th order spherical 

        pCS5   =     35,  // 5th order Cs5

        //  multipole aberrations (in Ang) 
        pC12a  =     12,  // astig - switch to these from above
        pC12b  =     13,  // because running short of parameters

        pC21a =      36,  // coma
        pC21b =      37,
        pC23a =      38,  // three-fold astig
        pC23b =      39,

        pC32a =      40,
        pC32b =      41,
        pC34a =      42,
        pC34b =      43,

        pC41a =      44,
        pC41b =      45,
        pC43a =      46,
        pC43b =      47,
        pC45a =      48,
        pC45b =      49,

        pC52a =      51,
        pC52b =      52,
        pC54a =      53,
        pC54b =      54,
        pC56a =      55,
        pC56b =      56,
        
        // sampling parameters
        pAX   =      101,  //  supercell size x, y
        pBY   =      102,
        pCZ   =      103,
        pNX   =      104,  // (int) main image size (transmission function)
        pNY   =      105,
        pNXPRB =     106,  // (int) probe size in pixels
        pNYPRB =     107,
        pNXOUT =     108,  // (int) output size (in pixels)
        pNYOUT =     109,
        pNWOBBLE=    111,  // (int) number of thermal configurations

        //  misc
        pTEMPER=     120,  // temperature (deg. K) for thermal vibrations
        pCAPERTMIN=  121,  // min. condencer angle in radians (CTEM)
        pSOURCE=     122,  // source size (in Angstroms) - for STEM
        pPPOSX=      123,  // probe position in Ang.
        pPPOSY=      124,  // probe position in Ang.
        pOAPMIN=     125,  // min. obj aperture (in radians) for annular ap.
        pMODE=       126,  // mode or program used
                           // 0=unknown, 1=atompot, 2=mulslice, 3=image,
                           // 4=probe, 5=stemslice, 6=autoslic, 7=autostem,
                           // 8=incostem
        pPROBEI=     127,  // probe current in Amp (STEM)
        pPROBEDT=    128,  // probe dwell time in Sec (STEM)
        pDOSE=       129,  // dose in electron/sqAng. (CTEM)
        //  confocal aberration - should convert to multipole form in future
        pCDF =       150,  // defocus
        pCDFA2 =     151,  // two-fold astig
        pCDFA2PHI =  152,
        pCDFA3 =     153,  // three-fold astig
        pCDFA3PHI =  154, 
        pCCS3 =      155,  //  3rd order spherical
        pCCS5 =      156,  //  5th order spherical
        pCCAPMIN =   157,  //  min aperture (in radians)
        pCCAPMAX =   158,  //  max aperture  (in radians)
        pCMINDET =   159,  // minimum detector angle (CONFOCAL) in Angstroms
        pCMAXDET =   160;  // maximum detector angle (CONFOCAL) in Angstroms

//  define calculation modes
enum{
    mUNKNOWN = 0,
    mATOMPOT = 1,
    mMULSLICE = 2,
    mIMAGE = 3,
    mPROBE = 4,
    mSTEMSLICE = 5,
    mAUTOSLICE = 6,
    mAUTOSTEM = 7,
    mINCOSTEM = 8,
    mMODEL3D = 9,
    mCBED = 10,
    mABBPHASE = 11
};

//-----  convert float back to integer and allow for small round off errror
//   neeeded to store integer parames in a float array
inline int ToInt( const float x ) { return( int( x + 0.1F ) ); }

/*--------------------- askYN() -----------------------------------*/
/*
    ask a yes/no question and return 1 or 0 for TRUE/FALSE

    message[] = question to ask
*/
#ifdef USE_TERMINAL
int askYN( const char message[] );
#endif

/*-------------------- bessi0() ---------------*/
/*
    modified Bessel function I0(x)
    see Abramowitz and Stegun page 379

    x = (double) real arguments

    12-feb-1997 E. Kirkland
 */
 double bessi0( double x );

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
 double bessk0( double x );

/*------------------------ chi() ---------------------*/
/*  return the aberration function
    must multiply this result by 2pi/wavelength before using
    C10= -df is defocus with the opposite sign

    put in a subroutine so I only have to get it right once
    
    now includes thru 5th order 4-jul-2010 ejk
    convert to param array p[] 2-may-2011 ejk
    add multiMode 8-may-2011 ejk

input:
    p = param[] array with aberration coeff. of chi (in Angstroms)
    alx, aly = x,y components of alpha in radians
    multiMode = if not 0 then use multipole aberrations
           set to 0 to speed up the calculation (no C34a etc.)

*/
double chi( vectorf &p, double alx, double aly, int multiMode );


/*--------------------- cputim() -----------------------------------*/
/*
   retrieve current CPU time in seconds
*/

double cputim();

/*--------------------- featom() -----------------------------------*/
/*
    return the electron scattering factor for atomic
    number Z at scattering angle k

    Z = atomic number 2 <= Z <= 103
    k2  = k*k where k =1/d = scattering angle (in 1/A)

  assumed global vars:

int feTableRead=0; = flag to remember if the param file has been read 
int nl=3, ng=3; = number of Lorenzians and Gaussians 
double fparams[][] = fe parameters

*/

double featom( int iz, double s );

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
void freqn( vectorf &ko, vectorf &ko2, vectorf &xo, int nk, double ak );


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
void messageSL( const char msg[],  int level = 0 );

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

   c         = input character string
*  islice(nsmax) = integer array to get stacking sequence indicies
   nsmax     = (integer) size of layer
   lmax      = (integer) maximum allowed layer index
*  nslice    = (integer) number of layers
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

   a * means that these variables may be modified by this routinec

  cname determines the mapping of characters into numbers.

*/
#ifdef USE_TERMINAL
int parlay( const char c[], vectori &islice, int nsmax, int lmax,
            int *nslice, int fperr );
#endif

/*---------------------------- ranflat -------------------------------*/
/*
    return a random number in the range 0.0->1.0
    with uniform distribution

    the 'Magic Numbers' are from 
        Numerical Recipes 2nd edit pg. 285
*/
double ranflat( unsigned long *iseed );

/*-------------------- rangauss() -------------------------- */
/*
    Return a normally distributed random number with 
    zero mean and unit variance using Box-Muller method
    
    ranflat() is the source of uniform deviates

    ref.  Numerical Recipes, 2nd edit. page 289
*/
double rangauss( unsigned long *iseed );


/*---------------------------- ranPoisson -------------------------------*/
/*
    return a random number with Poisson distribution
    and mean value mean
    
    There is a nice poisson RNG in c++11 but some compilers still do not
    have it, so make one here. The function lgamma()= log of gamma function
    is in the C99 standard and seems to have more widespread support
    (in mac osx 10.8 but NOT MSVS2010!).  Leave lgamma() code for future
    and approximate large means as Gaussian (which is really good above
    about 100).
    
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
int ranPoisson( double mean, unsigned long *iseed );


/*--------------------- readCnm() -----------------------*/
/*
    convert aberration line to a number 

    cline = character line in style "C32a  2.4"
    param[] = parameter array to get parameter value

    return 1 for success and 0 if parameter not recognized
*/
int readCnm( string &cline, vectorf &param, double x );


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
int ReadfeTable( );
 
 /*--------------------- ReadXYZcoord() -----------------------*/
/*
    read a set of (x,y,z) coordinates from a file
    and return number of coord. read
    
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
    line1 = char array to get 1st line of file
        with description
    nline1 = number of char in line1[]

   NOTE: x,y,z,occ and Znum array are allocated by this routine
    because it is not known ahead of time home many points
    will be read in (i.e. this routine figures it out)

    only infile is unchanged by this routine

    convert to string and streams mar-2014 ejk
*/
int ReadXYZcoord( const char* infile, const int ncellx, const int ncelly,
    const int ncellz, float *ax, float *by, float *cz, vectori &Znum, 
    vectorf &x, vectorf &y, vectorf &z, vectorf &occ, vectorf &wobble,
    string &line1 );

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
         vectord &d, int n, double x0 );

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

double sigma( double kev );

/*----------------- sortByZ() ------------------------------

    improved Shell sort modeled after prog. 6.5 (pg. 274) of
    R. Sedgewick, "Algorithms in C", 3rd edit. Addison-Wesley 1998
    
    x[], y[], z[]   = atom coordinates 
    occ[]           = occupancy of each atom
    Znum[]          = atomic number of each atom
    natom           = number of atoms
*/
void sortByZ( vectorf &x, vectorf &y, vectorf &z, vectorf &occ,
    vectori &Znum, int natom );

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
         vectord &b, vectord &c, vectord &d, int n);

/*------------------------- toString( int ) ----------------------*/
/*
    convert a number into a string
*/
std::string toString( int i );

/*------------------------- toString( float ) ----------------------*/
/*
    convert a number into a string
*/
std::string toString( float x );

/*------------------------- toString( double ) ----------------------*/
/*
    convert a number into a string
*/
std::string toString( double x );

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
                int nx, int ny );

/*--------------------- vatom() -----------------------------------*/
/*
    return the real space atomic potential (NOT projected)
    in volts for atomic number Z at radius r

    Z = atomic number 2 <= Z <= 103
    radius  = radius in Angstroms (MUST be > 0)

  assumed global vars:

int feTableRead=0; = flag to remember if the param file has been read 
int nl=3, ng=3; = number of Lorenzians and Gaussians 
double fparams[][] = fe parameters

  al and ag calculated using physical constants from:
    H. L. Anderson, editor "A Physicist's Desk Reference",
        2nd edition, Amer. Instit. Physics, 1989

  started from vzatom() 24-nov-1997 ejk
*/

double vatom( int Z, double radius );

/*--------------------- vzatom() -----------------------------------*/
/*
    return the real space projected atomic potential
    in volt-Angstroms for atomic number Z at radius r

    Z = atomic number 2 <= Z <= 103
    radius  = radius in Angstroms (MUST be > 0)

  assumed global vars:

int feTableRead=0; = flag to remember if the param file has been read 
int nl=3, ng=3; = number of Lorenzians and Gaussians 
double fparams[][] = fe parameters

  al and ag calculated using physical constants from:
    H. L. Anderson, editor "A Physicist's Desk Reference",
        2nd edition, Amer. Instit. Physics, 1989
*/

double vzatom( int Z, double radius );

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

    Z = atomic number 1 <= Z <= 98
    rsq = square of (radius in Angstroms)
*/

double vzatomLUT( int Z, double rsq );


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

double wavelength( double kev );

#endif