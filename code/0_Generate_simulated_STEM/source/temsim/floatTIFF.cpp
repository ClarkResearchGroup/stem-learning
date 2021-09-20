/*              *** floatTIFF.cpp ***

 ------------------------------------------------------------------------
Copyright 2012-2014 Earl J. Kirkland

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

   C++ class to read/write floating point images as TIFF format files (in ANSI-C++)
   and manage the memory storage.  There are three (3) images in each file.
   The first image is a simple 8 bit image with square pixels for display only.
   The second image is a 32 bit floating point image with possibly rectangular
   pixels for calculations.  The third image is a 2048 by 1 "image" of
   32 bit floating point parameters related t the calculation.

   This is mainly an easy to use high level front end derived from the old
   tiffsubs.c non-class library.

   The computer is assumed to have either big-endian or little-endian
   byte ordering and to be a 32/64 bit machine (i.e. a long must be 32 bits).

   These extended TIFF routines assume that the computer hardware uses 
   IEEE floating point (these may still work if this is not true but the 
   image files will not be transporable).

   These routines are thought to conform to TIFF version 6.0 
   (dated June 3, 1992) from:

      Adobe Developers Association
      Adobe Systems Incorporated
      1585 Charleston Road
      P.O. Box 7900
      Mountain View, CA 94039-7900

      http://www.adobe.com/Support/TechNotes.html
      ftp://ftp.adobe.com/pub/adobe/DeveloperSupport/TechNotes/PDFfiles

   NOTE-1: This software is thought to be correct but absolutely no guarantee
        whether written or implied is given.  This software is supplied on a
        user-beware "as is" basis.

   NOTE-2: These routines read in either byte order (Intel=little endian or
        Motorola=big endian) but write only the byte order of the computer
        they are running on.  It is assumed that the computer supports
        byte addressing and ASCII.
 
   NOTE-3:  Various symbol definitions (at top of file) can be changed
       for different computers and modes (the symbol CompByteOrder
       will be automatically determined at run time):

The source code is formatted for a tab size of 4.

----------------------------------------------------------
The public member functions are:

getParam()      : return value of parameter
getnpix()       : return value of npix
getDateTime()   :  get date and time of image

max()           : return maximum in image
min()           : return minimum value in image

nx()            : return horz. size of current image (in pixels)
ny()            : return vert. size of current image (in pixels)

maxParam()      : return maximum nu,ber of parameters

operator()(ix,iy) : return reference to pixel at (ix,iy)
                  :  a complex image is stored side-by-side (real-image) 
                  :  as one image with npix=2

read( file )    : read file with name 'file' into memory buffer 
                    (allocate memory buffer if necessary)

resize()        : resize current in-memory image (current data may be lost)

setParam()      : change value of parameter
setnpix()       : change value of npix (must be 1 for real or 2 for complex)

tFloatTest      : check if a file has a valid TIFF header
                        call before topenFloat() to avoid problems

write(file)     : write current image in memory to file with name 'file'

zeroParam()     : zero all parameter

----------------------------------------------------------

The private routines in this file are:

tclose:  close currently open tiff file and
             deallocate internal storage

tcreateFloatPixFile: create a floating point image file with an
    8 bit image at the beginning for viewing (can be read by
    treadFloatPix)

tifferr:  common error handler - print messages
        (internal use only)

topenFloat: open an extended TIFF floating point image
         for reading

tread: read bytes from file with byte reversal if necessary
        (internal use only)

treadIFD:  read the IFD or directory of currently open file
        (internal use only)

treadPix: read a 'standard' TIFF image (8 or 16 bit integer)
        file into memory

tsetByteOrder:  determine the byte order of the computer that
       this is running on and verify data type sizes

----------------------------------------------------------

The floating point image routines use an extended TIFF format
(i.e. they conform to the official TIFF standard but are not a
common permutation).  A floating point image (real or 
complex=real+imaginary) and a parameter array are stored in
a single file.  The first image in the file is a standard
8 bit greyscale image which most TIFF readers should be able
to handle.  This 8 bit image may be expanded in one direction
using bilinear interpolation to get square pixels.  The second
image in the file is stored as 32 bit IEEE floating point (for
simulation data) and the third image is one line with a 64
element floating point parameter array.

The order of the data in the file is:

<TIFF header>
unsigned 8 bit image data with square pixels
<IFD-1>
32 bit floating point image data (pixels may be rectangular)
<IFD-2>
32 bit parameter data
<IFD-3>
----------------------------------------------------------
----------------------------------------------------------

   started class version from tiffsubs.c 5-nov-2003 E. Kirkland
   switch from long to long32 with typdef to int/long so it
      will work on 64 bit machines 20-may-2008 ejk
   work on write(), min()/max() subroutines 11-mar-2012 ejk
   change to PMAX parameter length,
      add getDateTime() and zeroParam() 18-mar-2012 ejk
   add maxParam() member 21-mar-2012
   modified for wxWidgets (remove printf() for error messages )
      --   only change tifferr() and ifdef symbols 5-feb-2013 ejk
   convert sprintf()+char[] to strings 1-feb-2014 ejk
   convert to streams IO 1-feb-2014 ejk
   convert remaining buffers to use new/delete 2-feb-2014 ejk
   update wxGUI error message in tifferr() 10-feb-2014 ejk
   last modified 10-feb-2014 ejk

*/

#include <cstdlib>
#include <cstdio>
#include <cerrno>
#include <cstring>
#include <ctime>

#include "floatTIFF.hpp"    // class definition + inline functions here

//#define wxGUI   //    set for wxWidgets graphical user interface

#ifdef wxGUI
#include "wx/wx.h"  // regular headers
//#include "wx/wxprec.h"  // for precompiled headers
#endif  


//-------------- constants ----------------------------------------
// define the byte order types
const int tBigEndian=0x4d4d;        // = 'MM' = Motorola
const int tLittleEndian=0x4949;     // = 'II' = Intel

const int sizes[] = { 0, 1, 1, 2, 4, 8 };  // size in bytes of each TIFF data type

//---------------------------------------------------------------
//------------- begin functions here ----------------------
//---------------------------------------------------------------

//------------------ constructor --------------------------------

floatTIFF::floatTIFF()
{
    nIFD=0;
    nstrips=0;
    nstripcounts=0;
    stempo=NULL;
    stempc=NULL;
    ifd1 = NULL;

    StripByteCounts = NULL;
    StripOffsets = NULL;
    
    DateTime = "No-Date-Given";
    
    nxl = nyl = 0;      // initialize to "no-data" condition

    npix = 1;  // default to real

    PMAX = 512;   // maximum number of parameters
    param = new float[PMAX];
    if( NULL == param ) {   //  may exit() here
        tifferr( "floatTIFF cannot open allocate param array." );
    }

    tsetByteOrder();
}

//------------------ destructor ---------------------------------
floatTIFF::~floatTIFF()
{
    if( ifd1 != NULL ) {
        delete [] ifd1;  nIFD = 0;   ifd1 = NULL;
    }

    if( StripByteCounts != NULL) {  nstrips = 0;
        delete [] StripByteCounts; StripByteCounts = NULL; }

    if( StripOffsets != NULL) { nstripcounts = 0;
        delete [] StripOffsets;  StripOffsets = NULL; }

    if( stempo != NULL ) { delete [] stempo; stempo = NULL; }
    if( stempc != NULL ) { delete [] stempc; stempc = NULL; }

    if( (nxl>0) && (nyl>0) ) delete [] data;

    delete [] param;
}


//--------------------- getParam() ----------------------------------

float floatTIFF::getParam( const int i )
{

    if( (i>=0) && (i<=PMAX) ) return( param[i] );
    else return( 0 );

}  // end floatTIFF::getParam()

//-------------------- max() ----------------------------------
//  ipix = pix number ( 0 for real part and 1 for imag part )
//   ignored if there is only one pix
float floatTIFF::max( int ipix )
{   uslong ix, iy, nx, xo=0;
    float val, t;

    if( npix > 1 ){
        nx = nxl /npix ;
        if( ipix > 0 ) xo = nxl/2;
    }else {
        nx = nxl;
    }

    val = operator()(xo,0);
    for( iy=0; iy<nyl; iy++)
        for( ix=0; ix<nx; ix++){
            t = operator()(ix+xo,iy);
            if( t > val ) val = t;
        }

    return( val );
}

//-------------------- min() ----------------------------------
//  ipix = pix number ( 0 for real part and 1 for imag part )
//   ignored if there is only one pix
float floatTIFF::min( int ipix )
{   uslong ix, iy, nx, xo=0;
    float val, t;

    if( npix > 1 ){
        nx = nxl /npix ;
        if( ipix > 0 ) xo = nxl/2;
    }else {
        nx = nxl;
    }

    val = operator()(xo,0);
    for( iy=0; iy<nyl; iy++)
        for( ix=0; ix<nx; ix++){
            t = operator()(ix+xo,iy);
            if( t < val ) val = t;
        }

    return( val );
}


//--------------------- read() ----------------------------------
int floatTIFF::read( const char *filname )
{
    int is=-102;
    if( (is= topenFloat( filname )) != 1 ) {
        tifferr( "floatTIFF::read() cannot open input file." );
        return( is );
    }

    resize( ImageWidth, ImageLength );

    //------ read the image data --------------------------
    if( (is= treadFloatPix()) != 1 ) {
        tifferr( "floatTIFF::read() cannot read input file." );
        tclose();
        return( is );
    }

    //----- read() is succesful ---------------------------
    tclose();

    return ( 1 );

}  // end floatTIFF::read()


//--------------------- resize() ----------------------------------
// note: existing data (if any) may be destroyed

int floatTIFF::resize( const int nx, const int ny )
{

    //--- resize buffer ------------------------
    if( (nx != int(nxl)) || (ny != int(nyl)) ){
        if( (nxl != 0) && (nyl != 0) ) delete [] data;
        ImageWidth = nx;
        ImageLength = ny;
        data = new float [ ImageWidth*ImageLength ];
        if( NULL == data ) {
            tifferr( "Cannot allocate image storage in floatTIFF::resize()\n");
            tclose();
            return( -3 );
        }
        nxl = ImageWidth;
        nyl = ImageLength;
    }

    return( +1 );

}  // end floatTIFF::resize()

//--------------------- setParam() ----------------------------------

void floatTIFF::setParam( const int i, const float val )
{

    if( (i>=0) && (i<=PMAX) ) param[i] = val;

}  // end floatTIFF::setParam()


//--------------------- setnpix() ----------------------------------

void floatTIFF::setnpix( const int i )
{

    if( (i>=1) && (i<=2) )  npix = i;

}  // end floatTIFF::setParam()


/*--------------------- tclose -------------------*/
/*
    close a TIFF file opened by topen
    and free allocated data storage

    return 1 for success and </=0 for failure
*/

int floatTIFF::tclose( )
{
    tfp.close();
    if( tfp.fail() ) {
        tifferr( "tclose() cannot close the file.");
        return( -1 ); }

    if( ifd1 != NULL ) {
        delete [] ifd1;  nIFD = 0;   ifd1 = NULL;
    }

    if( StripByteCounts != NULL) {  nstrips = 0;
        delete [] StripByteCounts; StripByteCounts = NULL; }

    if( StripOffsets != NULL) { nstripcounts = 0;
        delete [] StripOffsets;  StripOffsets = NULL; }

    if( stempo != NULL ) { delete [] stempo; stempo = NULL; }
    if( stempc != NULL ) { delete [] stempc; stempc = NULL; }

    return(1);

}  /* end tclose */


/*--------------------- tFloatTest() -------------------*/
/*
   Test that a file is an extended TIFF file and
   has a floating point image (not an exhaustive test)

   filename  = pointer to string with filename

   return 1 for success and </=0 for failure
*/

int floatTIFF::tFloatTest( const char *filename )
{
   int iread;

    tsetByteOrder();
    if( tifftest( filename ) != 1 ) return( -1 );

    tfp.open( filename, ios_base::binary );
    if( tfp.fail() ) {
        tifferr( "tFloatTest() bad filename." );
        return( -2 ); }

    tfp.read( (char*)(&FileByteOrder), sizeof(short) );

    if( (tfp.fail()) || 
       !( (FileByteOrder==tLittleEndian) || 
                      (FileByteOrder==tBigEndian) ) ) {
        tfp.close( );
        return( -3 ); }

    tread( &Version, 2, 1, 2L );
    if( Version != 42 ) {       // the TIFF magic number
        tfp.close( );
        return( -3 ); }

    // get location of IFD  and read it

    tread( &CurrentIFD, 4, 1, 4L );
    iread = treadIFD( CurrentIFD );

    if( iread != 1 ) {
        tclose();
        return( -4 ); }

    // find the first floating point image
    while( ( SampleFormat != 3 ) &&
           ( BitsPerSample[0] != 8*sizeof(float) ) &&
           ( NextIFD > 0 ) ) {
        if( treadIFD( NextIFD ) < 1 ) {
            tclose( );;
            return( -5 );
        }
    }

    tclose();
    if( ( SampleFormat != 3 ) ||
           ( BitsPerSample[0] != 8*sizeof(float) ) ) return( -6 );

    return( 1 );

}  // end tFloatTest()


/*--------------------- tifferr() ------------------*/
/*
    common error handler - print messages
    only for internal use 
    DO NOT CALL FROM MAIN PROGRAM
*/

void floatTIFF::tifferr( const string &error_text )
{
#ifdef wxGUI
    string stemp;
    
    stemp = "floatTIFF Error: " + error_text;
    wxMessageBox( stemp, wxT("floatTIFF"),
                wxOK | wxICON_INFORMATION );  //??? , wxGetApp().GetTopFrame() );
#else

    //if( errno != 0 )
    //    fprintf( stderr, "System Error: %s\n", strerror( errno ) );    // ????
    cout << "floatTIFF Error: " <<  error_text.c_str() << endl;

#endif

    return;

}  /* end tifferr() */

/*--------------------- tifftest() -------------------*/
/*
   check if a file has a valid TIFF header

   return a 1 if it is a valid TIFF header, </=0 otherwise
   filename  = pointer to string with filename
*/
int floatTIFF::tifftest( const char *filename )
{
    ifstream tfp3;
    short byte01, byte23;

    tfp3.open( filename, ios_base::binary );
    if( tfp3.fail() ) {
       tifferr( "Bad filename in tifftest");
       return( -1 ); }

    tfp3.read( (char*)(&byte01), sizeof(short) );
    tfp3.read( (char*)(&byte23), sizeof(short) );

    if( tfp3.fail() ) {
       tifferr( "tifftest cannot read header.");
       return(-3); }

    tfp3.close();
    if( tfp3.fail() ) {
       tifferr( "tifftest cannot close file.");
       return( -2 ); }

   /*   tLittleEndian = "II" = least significant byte first 
    tBigEndian    = "MM" = most significant byte first
    42 = 0x2A = TIFF magic number  */

    if( (byte01!=tLittleEndian) && 
        (byte01!=tBigEndian) ) return(-4);
    
    if( (byte01==CompByteOrder) && (byte23!=0x002A) ) return(-5);
    if( (byte01!=CompByteOrder) && (byte23!=0x2A00) ) return(-5);

    return( +1 );

}  /* end tifftest */


/*--------------------- tinter() ----------------------------*/
/*
  Bilinear interpolation from data array

  nx,ny    = integer dimension of the image
  x,y      = real coord of image point to interpolate at
        0->(nx-1) and 0->(ny-1)
*/
float floatTIFF::tinter( long32 nx, long32 ny, double x, double y )
{
    long32 ix, iy, ix2, iy2;
    float  a, b, c, d, x1, y1, ans;

    ix = (int) x;
    iy = (int) y;
    if( ix > (nx-2) ) ix = (nx-2);
    if( ix <    0 )   ix = 0;
    if( iy > (ny-2) ) iy = (ny-2);
    if( iy <    0 )   iy = 0;
    ix2 = ix + 1;
    iy2 = iy + 1;
    x1 = (float) ix;
    y1 = (float) iy;

    d = operator()(ix,iy)  - operator()(ix,iy2) - operator()(ix2,iy)
                 + operator()(ix2,iy2);
    c = operator()(ix,iy2) - operator()(ix,iy)  - d*x1;
    b = operator()(ix2,iy) - operator()(ix,iy)  - d*y1;
    a = operator()(ix,iy)  - b*x1 - c*y1 - d*x1*y1;

    ans = (float) ( a + b*x + c*y + d*x*y );
    
    return( ans );

}  /* end tinter() */


/*--------------------- topenFloat() -------------------*/
/*
   open a extended TIFF file for reading as a floating
   point image and check that the header is a valid TIFF header

   filename  = pointer to string with filename

   return 1 for success and </=0 for failure
*/

int floatTIFF::topenFloat( const char *filename )
{
    int iread;

    tsetByteOrder();
    tfp.open( filename, ios_base::binary );
    if( tfp.fail() ) {
        tifferr( "topenFloat() bad filename." );
        return( -1 ); }

    tfp.read( (char*)(&FileByteOrder), sizeof(short) );

    if( (tfp.fail()) || 
       !( (FileByteOrder==tLittleEndian) || 
                      (FileByteOrder==tBigEndian) ) ) {
        tifferr( "Not a TIFF file in topenFloat().");
        return( -2 ); }

    tread( &Version, 2, 1, 2L );
    if( Version != 42 ) {           // the TIFF magic number
        tifferr( "Not a TIFF file in topenFloat().");
        return( -3 ); }

    /* get location of IFD  and read it */

    tread( &CurrentIFD, 4, 1, 4L );
    iread = treadIFD( CurrentIFD );

    if( iread != 1 ) {
        tifferr( "topenFloat() can't read first IFD.");
        return( iread ); }

    /* find the first floating point image   */
    while( ( SampleFormat != 3 ) &&
           ( BitsPerSample[0] != 8*sizeof(float) ) &&
           ( NextIFD > 0 ) ) {
        if( treadIFD( NextIFD ) < 1 ) {
            tifferr("Cannot read IFD in treadFloatPix.");
            return( -4 );
        }
    }

    if( ( SampleFormat != 3 ) ||
           ( BitsPerSample[0] != 8*sizeof(float) ) ) return( -5 );

    return( 1 );

}  /* end topenFloat() */


/*--------------------- tread() -------------------*/
/*
    this routine mimics the ANSI fread function
    however it fixes the byte order if required

    read n elements of length size into the buffer
    pointed to by bufptr from the current TIFF file
    at offset position (in bytes)

    if offset < 0 then do sequential reads from
    current position

    return 1 for success and </=0 for failure
    
    for internal use only
    DO NOT CALL FROM MAIN PROGRAM
*/

int floatTIFF::tread( void *bufptr, int size, int n, long32 offset )
{
   int i, j;
   unsigned char *temp, *ctmp;

    if( offset > 0 ) {
        tfp.seekg( offset, ios_base::beg );
        if( tfp.fail() ) {
            tifferr("Bad file seek in tread().");
            return( -1 ); } }
    
    tfp.read( (char*)bufptr, size*n );
    if( tfp.fail() ) {
        tifferr("Bad file read in tread().");
        return( -2 ); }

    if ( ( FileByteOrder != CompByteOrder ) && ( size > 1 ) )
      { ctmp = (unsigned char*) new char[ size ];
        temp = (unsigned char*) bufptr;
        for( i=0; i<n; i++) {
            for( j=0; j<size; j++) ctmp[j] = temp[i*size + j];
            for( j=0; j<size; j++) temp[i*size + j] = 
                    ctmp[size-1-j];
        }
        delete [] ctmp;
      }

    return( 1 );
        
}  /* end tread */

/*--------------------- treadFloatPix() -------------------*/
/*
    read in the current image file opend by topen()
    (i.e. topen() must be called before this routine)
    --- assumed to be a floating point image made by
        tcreateFloatPixFile() 

    this routine returns +1 for success 
        and </= 0 for failure

    NOTE: this routine makes some non-standard assumptions
    about the image file (i.e. floating point pixel data
    in a specific order - as created by tcreatFloatPix)

    return +1 for success, 0 or negative for fatal errors
        and >+1 for non-fatal errors

   started 25-mar-1996 ejk
   converted to class 17-nov-2003 ejk
   change param length to PMAX (=2048 in .hpp) 18-mar-2012 ejk
   fix maximum # param test 5-apr-2012 ejk

*/
int floatTIFF::treadFloatPix( )
{   
   float *FStripBuf;
   uslong ix, iy, iyi, iyf, nbytes;
   int istrip, BytesPerPixel, status, np;

    /* check that this is a floating point image */

    if( (BitsPerSample[0] != 8*sizeof(float)) || (SampleFormat!=3) ) {
       tifferr("floatTIFF::treadFloatPix cannot find floating point image.");
       return( -2 );
    }

    if( nstrips != nstripcounts ) {
       tifferr( "nstrips and nstripcounts do not agree in floatTIFF::treadFloatPix.");
       /* return( -3 );  leave return for rigorous TIFF */
    }

    if( SamplesPerPixel != 1 ) {
       tifferr("floatTIFF::treadFloatPix cannot handle more than 1 sample/pixel.");
       return( -4 ); }

    if( (nxl < (signed) ImageWidth ) || ( nyl < (signed) ImageLength ) ) {
       tifferr("Pix array size too small in treadFloatPix.");
       return( -5 ); }
       
    BytesPerPixel = sizeof(float);

    /* if the image is properly broken into strips then
        use them, otherwise read it in one row at a time
        to minimize the size of the intermediate storage
        buffer so that this will work on DOS type machines too */

    if( nstrips > 1 )
       nbytes = RowsPerStrip * ImageWidth;
    else
       nbytes = ImageWidth;

    FStripBuf = new float[ nbytes ];
    if( NULL == FStripBuf ) {
       tifferr( "floatTIFF::treadFloatPix unable to allocate pix temporary memory.");
       return( -6 ); }
    nbytes = nbytes * BytesPerPixel;

    if( nstrips > 1 ) {    /* use stipes if present  */
       iyi = 0;
       for( istrip=0; istrip<nstrips; istrip++ ) {
           if( tread( FStripBuf, BytesPerPixel, 
               (int) (StripByteCounts[istrip]/BytesPerPixel),
               (long32) StripOffsets[istrip]  ) != 1) {
                   delete [] FStripBuf;
                   tifferr("Bad file read in treadFloatPix.");
                   return( -7 ); }
          iyf = iyi + RowsPerStrip;
          if( iyf > ImageLength ) iyf = ImageLength;
          for( iy=iyi; iy<iyf; iy++) {
              for( ix=0; ix<ImageWidth; ix++)
              operator()(ix,iy) = FStripBuf[ix + (iy-iyi)*ImageWidth ];
          }
          iyi = iyf;
       } /* end for istrip loop */

    } else {        /* otherwise break it up into rows */
       tfp.seekg( StripOffsets[0], ios_base::beg );
       if( tfp.fail() ) {
            delete [] FStripBuf;
            tifferr("Bad file seek in floatTIFF::treadFloatPix.");
            return( -8 ); }
       for( iy=0; iy<ImageLength; iy++ ) {
          if( tread( FStripBuf, BytesPerPixel,   /* -1 = don't seek */
                 (int)(nbytes/BytesPerPixel), -1L ) != 1 ) {
                 delete [] FStripBuf;
              tifferr("Bad file read in floatTIFF::treadFloatPix");
              return( -9 ); }
         for( ix=0; ix<ImageWidth; ix++) operator()(ix,iy) = FStripBuf[ix];
       }
    } /* end if nstrip */

    delete [] FStripBuf;

   /*  read the parameter array from the 3rd IFD */

    if( NextIFD > 0 ) {
       if( status=treadIFD( NextIFD ) < 1 ) {
        tifferr("Cannot read parameter IFD in floatTIFF::treadFloatPix.");
        return( +2 );  /* wrong but not fatal */
       }
    } else {
        tifferr("No parameters available in floatTIFF::treadFloatPix.");
        return( +3 ); /* wrong but non-fatal */
    }

    if( (BitsPerSample[0] != 8*sizeof(float)) || (SampleFormat!=3) ) {
       tifferr("floatTIFF::treadFloatPix cannot read parameters.");
       return( +4 ); /* wrong but not fatal */
    }

    np = (int) StripByteCounts[0]/BytesPerPixel;   // number of parameters in file
    if( PMAX < np ) np = PMAX;  //   read up to PMAX parameters into param[] array
    if( tread( param, BytesPerPixel, np,
        (long32) StripOffsets[0]  ) != 1) {
            tifferr("Bad file read in floatTIFF::treadFloatPix.");
            return( -10 ); }
    npix = (int) ( param[0] + 0.5 );

   /*  all done */

   return( 1 );  /* success */

}  /* end treadFloatPix */

/*--------------------- treadIFD() -------------------*/
/*
    read the tiff directory structure (IFD)
    this is called by topen()

    return 1 for success and </=0 for failure

    for internal use only
    DO NOT CALL FROM MAIN PROGRAM
*/
int floatTIFF::treadIFD( uslong IFDoffset )
{
        usshort n1, *stemp;
        int i, j, n, ierr;
        uslong ltemp;

        string s, s2;   //  temp storage for error text
        stringstream ss;

        tread( &n1, 2, 1, IFDoffset);
        n = (int) n1;

        if( nIFD <= 0 ) ifd1  = (IFD*) new IFD[ n ];
        else if( n > nIFD ) {
                delete [] ifd1;
                ifd1  = (IFD*) new IFD[ n ];
        }

        if( ifd1 == NULL ){
                tifferr( "Unable to allocate memory in treadIFD(1).");
                return( -1 ); }

        nIFD = short(n);

        for( i=0; i<nIFD; i++ ){
                tread( &ifd1[i].tag   , 2, 1, -1L );
                tread( &ifd1[i].type  , 2, 1, -1L );
                tread( &ifd1[i].length, 4, 1, -1L );
                /* to get word order right */
                if( (ifd1[i].type == 3) && (ifd1[i].length <= 2 ) )
                        tread( &ifd1[i].value, 2, 2, -1L );
                else 
                        tread( &ifd1[i].value, 4, 1, -1L ); 
        }

        stemp = (usshort*) &ltemp;  /* for data type 3 separations */

        tread( &NextIFD, 4, 1, -1L);  /* point to next IFD */

        /* set defaults */

        BitsPerSample[0] = 1;
        BitsPerSample[1] = 0;
        BitsPerSample[2] = 0;
        Compression = 1;    /* no compression */
        DateTime= "No-Date-Given";
        FillOrder = 1;
        MaxSampleValue = 1; /* tag no longer recommended */
        MinSampleValue = 0; /* tag no longer recommended */
        NewSubFileType = 0;
        nstrips = 0;
        nstripcounts = 0;
        Orientation = 1; /* tags not recommended */
        PhotometricInterpretation = 1;  /* 0=black, 1+=white */
        PlanarConfiguration = 1;
        Predictor = 1;
        ResolutionUnit = 2;
        RowsPerStrip = 0; /* default is really to infinite- see below */
        SamplesPerPixel = 1;
        SampleFormat = 1;  /* default to unsigned integer */
        SubfileType = 0; /* tag no longer recommended */
        XPosition[0] = 0;
        XPosition[1] = 0;
        YPosition[0] = 0;
        YPosition[1] = 0;
        XResolution[0] = 0;
        XResolution[1] = 0;
        YResolution[0] = 0;
        YResolution[1] = 0;

        for( i=0; i<nIFD; i++ ) { ierr = 0;
           switch ( ifd1[i].tag ) {
           case 254:                       /* NewSubfileType tag */
                if( (ifd1[i].type != 4) || 
                    (ifd1[i].length != 1 ) ) {
                        ierr = -254; break; }
                NewSubfileType = ifd1[i].value;
                break;
           case 255:                       /* SubfileType tag */
                if( (ifd1[i].type != 3) || (ifd1[i].length != 1 ) )
                        { ierr = -255; break; }
                ltemp = ifd1[i].value;
                SubfileType =  stemp[0];
                break;
           case 256:                       /* ImageWidth tag */
                if((ifd1[i].length != 1 ) ) {
                        ierr = -256; break; }
                if( ifd1[i].type == 3 ) {  /* type = short */
                   ltemp = ifd1[i].value;
                   ImageWidth =  stemp[0]; }
                else if( ifd1[i].type == 4 )   /* type = long */
                   ImageWidth = ifd1[i].value;
                else { ierr = -256; break; }
                break;
           case 257:                       /* ImageLength tag */
                if((ifd1[i].length != 1 ) ) {
                        ierr = -257; break; }
                if( ifd1[i].type == 3 ) {  /* type = short */
                   ltemp = ifd1[i].value;
                   ImageLength =  stemp[0]; }
                else if( ifd1[i].type == 4 )   /* type = long */
                   ImageLength = ifd1[i].value;
                else { ierr = -257; break; }
                break;
           case 258:                       /* BitsPerSample tag */
                if( (ifd1[i].type != 3) || 
                    (ifd1[i].length < 1 ) || (ifd1[i].length > 3) ) {
                        ierr = -258; 
                        ss << "BitsPerSample type= " << ifd1[i].type << ", length= " 
                                << ifd1[i].length << endl;
                        s = ss.str();
                        break; }
                if( ifd1[i].length == 1 ) {
                   ltemp = ifd1[i].value;
                   BitsPerSample[0] = stemp[0];
                   if( (BitsPerSample[0] < 1)
                     || (BitsPerSample[0] > 32 ) )  {
                      ss << "Bad BitsPerSample = " << BitsPerSample[0] 
                                << ", will try 8." << endl;
                      s2 = ss.str();
                      tifferr( s2 );
                      BitsPerSample[0] = 8;  }
                   }
                else if( ifd1[i].length == 2 ) {
                   ltemp = ifd1[i].value;
                   BitsPerSample[0] = stemp[0];
                   BitsPerSample[1] = stemp[1];  }
                else if( ifd1[i].length == 3 )
                   tread( BitsPerSample, 2, 3, ifd1[i].value );
                break;
           case 259:                       /* Compression tag */
                if( (ifd1[i].type != 3) || 
                    (ifd1[i].length != 1 ) ) {
                        ierr = -259;
                        ss << "Compression type= " << ifd1[i].type << 
                                ", length= " << ifd1[i].length << endl;
                        s = ss.str();
                        break; }
                ltemp = ifd1[i].value;
                Compression = stemp[0];
                break;
           case 262:             /* PhotometricInterpretation tag */
                if( (ifd1[i].type != 3) || 
                    (ifd1[i].length != 1 ) ) {
                        ierr = -262; break; }
                ltemp = ifd1[i].value;
                PhotometricInterpretation = stemp[0];
                break;
           case 266:                    /* FillOrder tag */
                if( (ifd1[i].type != 3) || 
                    (ifd1[i].length != 1 ) ) {
                        ierr = -266; break; }
                ltemp = ifd1[i].value;
                FillOrder = stemp[0];
                break;
           case 273:                       /* StripOffsets tag */
                n = (int) nstrips;
                nstrips = (int) ifd1[i].length;
                if( StripOffsets == NULL ) {
                        StripOffsets = (uslong*) new long32[ nstrips ];
                        stempo = (usshort*) new short[ nstrips ];
                } else if( nstrips > n ) {
                        delete [] StripOffsets;
                        delete [] stempo;
                        StripOffsets = (uslong*) new long32[ nstrips ];
                        stempo = (usshort*) new short[ nstrips ]; }
                if( (StripOffsets == NULL) || ( stempo == NULL) ) {
                   tifferr( "treadIFD is unable to allocate Offsets"
                           " memory(2).");
                   return( -2 ); }
                   
                /*  check if IFDvalues is a value or an offset */
                if( (ifd1[i].length * sizes[ ifd1[i].type ]) > 4 ) {
                   if( ifd1[i].type == 3 ) {  /* type = short */
                      tread( stempo, 2, nstrips, ifd1[i].value );
                      for( j=nstrips-1; j>=0; j--)
                       StripOffsets[j] = (uslong) stempo[j]; }
                   else if( ifd1[i].type == 4 )   /* type = long */
                      tread( StripOffsets, 4, nstrips, ifd1[i].value );
                   else { ierr = -273; break; } }
                else {
                   if( ifd1[i].type == 3 ) {  /* type = short */
                    ltemp = ifd1[i].value;
                    StripOffsets[0] = stemp[0]; }
                   else if( ifd1[i].type == 4 )   /* type = long */
                    StripOffsets[0] = ifd1[i].value;
                   else { ierr = -273; break; } }
                break;
          case 274:                       /* Orientation tag */
                if( (ifd1[i].type != 3) || (ifd1[i].length != 1 ) )
                        { ierr = -274; break; }
                ltemp = ifd1[i].value;
                Orientation =  stemp[0];
                break;
           case 277:                       /* SamplesPerPixel tag */
                if( (ifd1[i].type != 3) || 
                    (ifd1[i].length != 1 ) ) {
                        ierr = -277; break; }
                ltemp = ifd1[i].value;
                SamplesPerPixel = stemp[0];
                break;
           case 278:                       /* RowsPerStrip tag */
                if( ifd1[i].length != 1 ) {
                        ierr = -278;
                        ss << "RowsPerStrip type= " << ifd1[i].type << 
                                ", length= " << ifd1[i].length << endl;
                        s = ss.str();
                        break; }
                if( ifd1[i].type == 3 ) {  /* type = short */
                   ltemp = ifd1[i].value;
                   RowsPerStrip =  stemp[0]; }
                else if( ifd1[i].type == 4 )   /* type = long */
                   RowsPerStrip = ifd1[i].value;
                else { ierr = -278; break; }
                break;
           case 279:                       /* StripByteCounts tag */
                n = nstripcounts;
                nstripcounts = (int) ifd1[i].length;
                if( StripByteCounts == NULL ) {
                        StripByteCounts = (uslong*) new long32[ nstrips ];
                        stempc = (usshort*) new long32[ nstrips ];
                } else if( nstripcounts > n ) {
                        delete [] StripByteCounts;
                        delete [] stempc;
                        StripByteCounts = (uslong*) new long32[ nstrips ];
                        stempc = (usshort*) new long32[ nstrips ]; }
                if( (StripByteCounts == NULL) || ( stempc == NULL) ) {
                   tifferr( "treadIFD is unable to allocate memory(3).");
                   return( -3 ); }

                /*  check if IFDvalues is a value or an offset */
                if( (ifd1[i].length * sizes[ ifd1[i].type ]) > 4 ) {
                   if( ifd1[i].type == 3 ) {  /* type = short */
                      tread( stempc, 2, nstrips, ifd1[i].value );
                      for( j=nstrips-1; j>=0; j--)
                           StripByteCounts[j] = (uslong) stempc[j]; }
                   else if( ifd1[i].type == 4 )   /* type = long */
                      tread( StripByteCounts, 4, nstrips, ifd1[i].value );
                   else { ierr = -279; break; } }
                else {
                   if( ifd1[i].type == 3 ) {  /* type = short */
                        ltemp = ifd1[i].value;
                        StripByteCounts[0] = stemp[0]; }
                   else if( ifd1[i].type == 4 )   /* type = long */
                        StripByteCounts[0] = ifd1[i].value;
                   else { ierr = -279; break; } }
                break;
           case 280:                       /* MinSampleValue tag */
                if( (ifd1[i].type != 3) || (ifd1[i].length != 1 ) )
                        { ierr = -280; break; }
                ltemp = ifd1[i].value;
                MinSampleValue =  stemp[0];
                break;
           case 281:                       /* MaxSampleValue tag */
                if( (ifd1[i].type != 3) || (ifd1[i].length != 1 ) )
                        { ierr = -281; break; }
                ltemp = ifd1[i].value;
                MaxSampleValue =  stemp[0];
                break;
           case 282:                       /* XResolution tag */
                if( (ifd1[i].length != 1 ) ||
                    (ifd1[i].type != 5 )   ) {
                        ierr = -282; break; }
                tread( &XResolution[0], 4, 1, ifd1[i].value );
                tread( &XResolution[1], 4, 1, -1L );
                break;
           case 283:                       /* YResolution tag */
                if( (ifd1[i].length != 1 ) ||
                    (ifd1[i].type != 5 )   ) {
                        ierr = -283; break; }
                tread( &YResolution[0], 4, 1, ifd1[i].value );
                tread( &YResolution[1], 4, 1, -1L );
                break;
           case 284:                       /* PlanarConfiguration tag */
                if( (ifd1[i].type != 3) || 
                    (ifd1[i].length != 1 ) ) {
                        ierr = -284; 
                        ss << "PlanarConfiguration type= " << ifd1[i].type << ", length= " 
                                << ifd1[i].length << endl;
                        s = ss.str();
                        break; }
                ltemp = ifd1[i].value;
                PlanarConfiguration = stemp[0];
                break;
           case 286:                       /* XPosition tag */
                if( (ifd1[i].length != 1 ) ||
                    (ifd1[i].type != 5 )   ) {
                        ierr = -286; break; }
                tread( &XPosition[0], 4, 1, ifd1[i].value );
                tread( &XPosition[1], 4, 1, -1L );
                break;
           case 287:                       /* YPosition tag */
                if( (ifd1[i].length != 1 ) ||
                    (ifd1[i].type != 5 )   ) {
                        ierr = -287; break; }
                tread( &YPosition[0], 4, 1, ifd1[i].value );
                tread( &YPosition[1], 4, 1, -1L );
                break;
           case 296:             /* ResolutionUnit tag */
                if( (ifd1[i].type != 3) || 
                    (ifd1[i].length != 1 ) ) {
                        ierr = -296; break; }
                ltemp = ifd1[i].value;
                ResolutionUnit = stemp[0];
                break;
           case 306:             /* (ASCII) DateTime tag */
                if( ifd1[i].type != 2 ) {
                        ierr = -296; break; }
                char ctemp[100];
                //tread( DateTime, 1, 20, ifd1[i].value);   //  ????
                tread( ctemp, 1, 20, ifd1[i].value);
                ctemp[20] = 0; //  null terminator
                DateTime = string( ctemp );
                break;
           case 317:                       /* Predictor tag */
                if( (ifd1[i].type != 3) || 
                    (ifd1[i].length != 1 ) ) {
                        ierr = -317; break; }
                ltemp = ifd1[i].value;
                Predictor = stemp[0];
                break;
           case 339:             /* SampleFormat tag */
                if( (ifd1[i].type != 3) || 
                    (ifd1[i].length != 1 ) ) {
                        ierr = -339; break; }
                ltemp = ifd1[i].value;
                SampleFormat = stemp[0];
                break;
           default:
                ierr = -1000;
                ss << "warning, treadIFD unknown TIFF tag # " << ifd1[i].tag << endl;
                s = ss.str();
                break;
           }
          if ( ierr != 0 ) {
                tifferr( s );
                if( ierr > -999 ) return( ierr ); }

        }  /* end for(i=0; i<nIFD; i++) */

        /* if RowsPerStrip is missing make a guess */
        if( RowsPerStrip == 0 ) RowsPerStrip = ImageLength;

        return( 1 );

}  /* end treadIFD */



/*--------------------- tsetByteOrder() -------------------*/
/*
   determine the byte order (big-endian or little-endian)
   of the computer this is running on

   also double check the data type lengths to be sure
   this is a standard 32 bit machine
*/
void floatTIFF::tsetByteOrder()
{
   usshort us;
   unsigned char *uc;
   stringstream ss;

   int sizec, sizes, sizel, sizef, sized;

   /*  find byte ordering and save results */
   uc = (unsigned char*) &us;
   uc[0] = 1;
   uc[1] = 0;
   if( us > 1 ) CompByteOrder = tBigEndian;
          else  CompByteOrder = tLittleEndian;

   /* test data type lengths */
   sizec = sizeof(char);
   sizes = sizeof(short);
   sizel = sizeof(long32);
   sizef = sizeof(float);
   sized = sizeof(double);
   if( (sizec!=1) || (sizes!=2) || (sizel!=4)
    || (sizef!=4) || (sized!=8) ) {
    ss << "Sorry, can't continue...\n"
       << "This TIFF library requires that the C data types\n"
       << "char, short, long, float and double \n"
       << "have lengths (in bytes) of 1, 2, 4, 4, 8\n"
       <<"but this computer has lengths of "
        << sizec << ", " << sizes << ", "  << sizel << ", " 
        << sizef<< ", "  << sized << endl;
    tifferr( ss.str() );
    exit( 0 );
   }


}  /* end tsetByteOrder() */

/*--------------------- write() ----------------*/
/*
    create a 32 bit floating point image file with an
    8 bit greyscale image in the front in TIFF format
    (the first image is a standard 8 bit greyscale image
    for a standard TIFF reader and the 2nd image is a
    32 bit floating point image for number crunching
     in extended TIFF format - both represent the same pix )

    WARNING: This conforms to the TIFF standard but
    is an uncommon usage.  Also it assumes that the
    computer it is running on uses IEE floating point.

  the data, param's, nx, ny, npix must be set before calling
  this function

  arguments used to scale the 8 bit image:
        rmax  = param[1] = max real part
        imax  = param[2] = max imag part
        rmin  = param[3] = min real part
        imin  = param[4] = min imag part
    dxi   = param[14]
        dyi   = param[15] = determine the aspect ratio

   started 18-mar-1996 ejk
   made 8-bit image pixels be square via interpolation/expansion
        14-apr-1996 ejk
   remove extra free()'s of non-existant memory 19-feb-1997 ejk
   convert to class version 18-nov-2003 ejk
   updated 18-mar-2012 ejk
   convert buffers to new/delete and convert to fstream IO 2-feb-2014 ejk
     (nwrite usage still as for fwrite() not streams so it looks funny)
*/

int floatTIFF::write( const char *file, float rmin, float rmax, float imin, float imax,
    float dxi, float dyi )
{
#define  ntags 14

   unsigned char uctemp, *StripBuf;
   usshort *stemp, st;
   short version;
   int i, strip;
   uslong offset, offsetIFD1, offsetIFD2, offsetIFD1to2,
      offsetIFD3, offsetIFD2to3;
   long32 ltemp, StripsPerImage, xy[2], nwrite, nx2, ny2;
   long32 *Offsets, *ByteCounts, ix, iy, nrows, count;
   IFD ifd2[ntags];
   ofstream tfp2;
   float scaler, scalei, dx, dy, *FStripBuf;
   double x, y, z;

   int nx, ny;

   time_t caltime;
   struct tm *mytime;

    tsetByteOrder();

   /* test for valid arguments */

    scaler = rmax - rmin;   //???  old param[1] - param[3];
    if( scaler < 0.0F ) scaler = -scaler;
    if( scaler  < 1.0E-20 ) return( -1 );

    if( (npix<1) || (npix>2) ) return( -2 );

   /* first write the 8 bit image */
   /* expand to get square pixels so simple viewers will work
    nx2, ny2 = new output size (in pixels)  */

    nx = int(nxl);  // to use the old code from tiffsubs.c
    ny = int(nyl);

    dx = dxi;  //???  old  param[14];
    if( dx <= 0.0F ) dx = 1.0F;
    dy = dyi;  //???  old param[15];
    if( dy <= 0.0F ) dy = dx;
    if( dx > dy ) {
        nx2 = (long32) ( nx * dx/dy + 0.5F);
        ny2 = ny;
    } else if( dx < dy ) {
        nx2 = nx;
        ny2 = (long32) ( ny * dy/dx + 0.5F);;
    } else {
        nx2 = nx;
        ny2 = ny;
    }
    
    /* rows per strip, make strips about 8k */  
    nrows = 8192 / (nx2) ;
    if( nrows <= 0 ) nrows = 1;
    if( nrows > ny2 ) nrows = ny2;  /* if whole image < 1 strip */

    for( i=0; i<ntags; i++) ifd2[i].length = 1;  /* init tags */

    /* set up easy tags */

    ifd2[0].tag = 254;  /* NewSubfileType */
    ifd2[0].type = 4;
    ifd2[0].value = 0;

    ifd2[1].tag = 256;  /* ImageWidth */
    ifd2[1].type = 4;
    ifd2[1].value = (long32) nx2;

    ifd2[2].tag = 257;  /* ImageLength */
    ifd2[2].type = 4;
    ifd2[2].value = (long32) ny2;

    ifd2[3].tag = 258;  /* BitsPerSample */
    ifd2[3].type = 3;
    ifd2[3].value = 8;

    ifd2[4].tag = 259;  /* Compression */
    ifd2[4].type = 3;
    ifd2[4].value = 1;  /* no compression */

    ifd2[5].tag = 262;  /* PhotmetricInterpretation */
    ifd2[5].type = 3;
    ifd2[5].value = 1;  /* black is zero */

/*  ifd2[6].tag = 273;   - StripOffsets - see below */

    ifd2[7].tag = 277;  /* SamplesPerPixel */
    ifd2[7].type = 3;
    ifd2[7].value = 1;

    ifd2[8].tag = 278;  /* RowsPerStrip */
    ifd2[8].type = 4;
    ifd2[8].value = nrows;

/*  ifd2[9].tag = 279;  - StripByteCounts - see below */

/*  ifd2[10].tag = 282; - XResolution - see below */

/*  ifd2[11].tag = 283; - YResolution - see below */

    ifd2[12].tag = 296; /* ResolutionUnit */
    ifd2[12].type = 3;
    ifd2[12].value = 2; /* inches */

    ifd2[13].tag = 339; /* SampleFormat */
    ifd2[13].type = 3;
    ifd2[13].value = 1; /* unsigned integer data */

    { /* word swap short values in IFD*/
    if( CompByteOrder == tBigEndian )
        for( i=0; i<ntags; i++) 
           if( (ifd2[i].type == 3) && (ifd2[i].length <= 1 ) ) {
            stemp = (usshort*) &ifd2[i].value;
            st = stemp[0];  stemp[0] = stemp[1];
            stemp[1] = st;
        }
    }

    StripsPerImage = ( ny2 + nrows - 1 ) / nrows;

    Offsets = new long32[ StripsPerImage ];
    if( NULL == Offsets ) {
        tifferr( "floatTIFF::write cannot allocate memory.");
        return( -3 ); }
    ByteCounts = new long32[ StripsPerImage ];
    if( NULL == ByteCounts ) {
        delete [] Offsets;
        tifferr( "floatTIFF::write cannot allocate memory.");
        return( -4 ); }

    count = nrows * nx2;    /* pixel count for StripBuf */
    StripBuf = new unsigned char[ count ];
    if( NULL == StripBuf ) {
        delete [] Offsets;  delete [] ByteCounts;
        tifferr( "floatTIFF::write cannot allocate memory.");
        return( -5 ); }

    tfp2.open( file, ios_base::binary );
    if( tfp2.fail()) {
        delete [] Offsets;  delete [] ByteCounts; delete [] StripBuf;
        tifferr( "floatTIFF::write cannot open file.");
        return( -6 ); }

    offset = 8; /* start writing after the header */
            /* i.e. fill in header later */
    tfp2.seekp( offset, ios_base::beg);
    if( tfp2.fail() ) {
        delete [] Offsets;  delete [] ByteCounts; delete [] StripBuf;
        tifferr( "Bad file seek in floatTIFF::write().");
        return( -6 ); }

    /*  write image raster data --
        note this will interpolate inbetween real and imaginary
            images in a strange way */
    strip = 0;
    count = 0;
    //???  old rmin = param[3];
    //???  old imin = param[4];
    scaler = 255.0F/(rmax - rmin);  //???  old (param[1] - rmin);
    if( npix > 1 ) scalei = 255.0F/(imax - imin);  //??? (param[2] - imin);
    for( iy=0; iy<ny2; iy++){
       y = (ny-1) * ((double)iy)/((double)(ny2-1)) ;
       for( ix=0; ix<nx2; ix++) {
        x = (nx-1) * ((double)ix)/((double)(nx2-1)) ;
        z = tinter( nx, ny, x, y);
        if( ix < nx2/npix ) ltemp = (long32) ( scaler*(z-rmin) + 0.5 );
                   else ltemp = (long32) ( scalei*(z-imin) + 0.5 );
        if( ltemp < 0 ) ltemp = 0;
        if( ltemp > 255 ) ltemp = 255;
        StripBuf[count+ix] = (unsigned char) (ltemp & 0x00FF); 
       }
       count += nx2;

       if( ( nrows*((iy+1)/nrows) == (iy+1) ) || ( (iy+1)==ny2 ) ) { 
          tfp2.write( (char*) StripBuf, count);
         if( tfp2.fail() ) {
           delete [] Offsets;  delete [] ByteCounts; delete [] StripBuf;
           tifferr("floatTIFF::write cannot write image data");
           return( -8 ); }
        Offsets[ strip ] = offset;
        offset += count;
        ByteCounts[strip] = count;
        count = 0;  strip += 1;
       } /* end if nrows */
    }  /* end for iy loop */

    nwrite = 0;
    ifd2[6].tag = 273;  /* StripOffsets */
    ifd2[6].type = 4;
    ifd2[6].length = strip;
    if( strip > 1 ) {
        ifd2[6].value = offset;
        tfp2.write( (char*) Offsets, strip*sizeof(long32) );
        nwrite += strip;
        offset += strip*sizeof(long32);
    } else {
        ifd2[6].value = Offsets[0];
        nwrite += 1;
    }

    ifd2[9].tag = 279;  /* StripByteCounts */
    ifd2[9].type = 4;
    ifd2[9].length = strip;
    if( strip > 1 ) {
        ifd2[9].value = offset;
        tfp2.write( (char*) ByteCounts, strip*sizeof(long32) );
        nwrite += strip;
        offset += strip*sizeof(long32);
    } else {
        ifd2[9].value = ByteCounts[0];
        nwrite += 1;
    }

    xy[0] = (long32) ( ((double) nx2) / 0.04 );  /* 4 inch wide image */
    xy[1] = 100;
    tfp2.write( (char*) xy, 2*sizeof(long32) );
    nwrite += 2;
    ifd2[10].tag = 282; /* XResolution */
    ifd2[10].type = 5;
    ifd2[10].value = offset;
    offset += 2*sizeof(long32);

    /* remember that the aspect ratio is 1 for first 8 bit pix */
    tfp2.write( (char*) xy, 2*sizeof(long32) );
    nwrite += 2;
    ifd2[11].tag = 283; /* YResolution */
    ifd2[11].type = 5;
    ifd2[11].value = offset;
    offset += 2*sizeof(long32);

    /* remember to store IFD on a word boundary */
    if( (offset/2)*2 != offset ) {
        uctemp = 0;
        tfp2.write( (char*) (&uctemp), 1 );
        offset += 1; }

    offsetIFD1 = offset;    /* pointer to first IFD */

    version = ntags;    /*  write the first IFD */
    tfp2.write( (char*) (&version), sizeof(short) );
    nwrite += 1;
    tfp2.write( (char*) ifd2, ntags * sizeof(IFD) );
    nwrite += ntags;

    /* save location of pointer to 2nd IFD, and write 0 for now */
    offset += sizeof(short) + ntags*sizeof(IFD);
    offsetIFD1to2 = offset;
    ltemp = 0;
    tfp2.write( (char*) (&ltemp), sizeof(long32) );
    nwrite += 1;
    if( tfp2.fail() || (nwrite != ( 2*strip + ntags + 6 ) ) )
       {    tifferr("floatTIFF::write cannot write IFD.");
        delete [] Offsets;  delete [] ByteCounts; delete [] StripBuf;
        return( -9 ); }
    offset += sizeof(long32);

    delete [] StripBuf;
 
   /* next write the 32 bit floating point image */

    /* rows per strip, make strips about 8k */  
    nrows = 8192 / ( nx * sizeof(float) );
    if( nrows <= 0 ) nrows = 1;
    if( nrows > ny ) nrows = ny;  /* if whole image < 1 strip */

    /* change appropriate tags - most are the same */

    ifd2[1].tag = 256;  /* ImageWidth */
    ifd2[1].type = 4;
    ifd2[1].value = (long32) nx;

    ifd2[2].tag = 257;  /* ImageLength */
    ifd2[2].type = 4;
    ifd2[2].value = (long32) ny;

    ifd2[3].tag = 258;  /* BitsPerSample */
    ifd2[3].type = 3;
    ifd2[3].value = 8 * sizeof(float);

    ifd2[4].tag = 259;  /* Compression */
    ifd2[4].type = 3;
    ifd2[4].value = 1;  /* no compression */

    ifd2[5].tag = 262;  /* PhotmetricInterpretation */
    ifd2[5].type = 3;
    ifd2[5].value = 1;  /* black is zero */

/*  ifd2[6].tag = 273;   - StripOffsets - see below */

    ifd2[7].tag = 277;  /* SamplesPerPixel */
    ifd2[7].type = 3;
    ifd2[7].value = 1;

    ifd2[8].tag = 278;  /* RowsPerStrip */
    ifd2[8].type = 4;
    ifd2[8].value = nrows;

/*  ifd2[9].tag = 279;  - StripByteCounts - see below */

/*  ifd2[10].tag = 282; - XResolution - see below */

/*  ifd2[11].tag = 283; - YResolution - see below */

    ifd2[12].tag = 296; /* ResolutionUnit */
    ifd2[12].type = 3;
    ifd2[12].value = 2; /* inches */

    ifd2[13].tag = 339; /* SampleFormat */
    ifd2[13].type = 3;
    ifd2[13].value = 3; /* IEEE floating point */

    if( CompByteOrder == tBigEndian ) { /* word swap short values in IFD*/
       for( i=3; i<ntags; i++) {
          if( (ifd2[i].type == 3) && (ifd2[i].length <= 1 ) ) {
        stemp = (usshort*) &ifd2[i].value;
        st = stemp[0];  stemp[0] = stemp[1];
        stemp[1] = st;
          }
       }
    }

    StripsPerImage = ( ny + nrows - 1 ) / nrows;

    delete [] Offsets;
    Offsets = new long32[ StripsPerImage ];
    if( NULL == Offsets ) {
        tifferr( "floatTIFF::write cannot allocate memory.");
        return( -10 ); }
    delete [] ByteCounts;
    ByteCounts = new long32[ StripsPerImage ];
    if( NULL == ByteCounts ) {
        delete [] Offsets;
        tifferr( "floatTIFF::write cannot allocate memory.");
        return( -11 ); }

    count = nrows * nx ;    /* pixel count for StripBuf */
    FStripBuf = new float[count];
    if( NULL == FStripBuf ) {
        delete [] Offsets;  delete [] ByteCounts;
        tifferr( "floatTIFF::write cannot allocate memory.");
        return( -12 ); }

    /* continue writing where we left off - after the first IFD */
    /* it should be positioned OK from 8 bit write above */

    strip = 0;      /*  write image raster data */
    count = 0;
    for( iy=0; iy<ny; iy++){
       for( ix=0; ix<nx; ix++) {
        FStripBuf[count+ix] = operator()(ix,iy);
       }
       count += nx;

       if( ( nrows*((iy+1)/nrows) == (iy+1) ) || ( (iy+1)==ny ) ) { 
         tfp2.write( (char*) FStripBuf, count * sizeof(float) );
         if( tfp2.fail() ) {
           delete [] Offsets;  delete [] ByteCounts; delete [] FStripBuf;
           tifferr("floatTIFF::write() cannot write image data");
           return( -13 ); }
        Offsets[ strip ] = offset;
        offset += count*sizeof(float);
        ByteCounts[strip] = count*sizeof(float);
        count = 0;  strip += 1;
       } /* end if nrows */
    }  /* end for iy loop */

    delete [] FStripBuf;
    nwrite = 0;
    ifd2[6].tag = 273;  /* StripOffsets */
    ifd2[6].type = 4;
    ifd2[6].length = strip;
    if( strip > 1 ) {
        ifd2[6].value = offset;
        tfp2.write( (char*) Offsets, strip*sizeof(long32) );
        nwrite += strip;
        offset += strip*sizeof(long32);
    } else {
        ifd2[6].value = Offsets[0];
        nwrite += 1;
    }

    ifd2[9].tag = 279;  /* StripByteCounts */
    ifd2[9].type = 4;
    ifd2[9].length = strip;
    if( strip > 1 ) {
        ifd2[9].value = offset;
        tfp2.write( (char*) ByteCounts, strip*sizeof(long32) );
        nwrite += strip;
        offset += strip*sizeof(long32);
    } else {
        ifd2[9].value = ByteCounts[0];
        nwrite += 1;
    }

    xy[0] = (long32) ( (double) nx / 0.04 );  /* 4 inch wide image */
    xy[1] = 100;
    tfp2.write( (char*) xy, 2*sizeof(long32) );
    nwrite += 2;
    ifd2[10].tag = 282; /* XResolution */
    ifd2[10].type = 5;
    ifd2[10].value = offset;
    offset += 2*sizeof(long32);

    //dx = param[14];  /* determine the aspect ratio */  ????
    //dy = param[15];
    if( (dx > 0.0F) && (dy > 0.0F) )
        xy[0] = (long32) ( (dy/dx)* ((double)nx) / 0.04 );
    xy[1] = 100;
    tfp2.write( (char*) xy, 2*sizeof(long32) );
    nwrite += 2;
    ifd2[11].tag = 283; /* YResolution */
    ifd2[11].type = 5;
    ifd2[11].value = offset;
    offset += 2*sizeof(long32);

    /* remember to store IFD on a word boundary */
    if( (offset/2)*2 != offset ) {
        uctemp = 0;
        tfp2.write( (char*) (&uctemp), 1 );
        offset += 1; }

    offsetIFD2 = offset;
    version = ntags;    /*  write the second IFD */
    tfp2.write( (char*) (&version), sizeof(short) );
    nwrite += 1;
    tfp2.write( (char*) ifd2, ntags * sizeof(IFD) );
    nwrite += ntags;
    offset += sizeof(short) + ntags*sizeof(IFD);
    offsetIFD2to3 = offset;
    ltemp = 0;
    tfp2.write( (char*) (&ltemp), sizeof(long32) );
    nwrite += 1;
    if( tfp2.fail() || (nwrite != ( 2*strip + ntags + 6 ) ) )
       {    tifferr("floatTIFF::write cannot write IFD2.");
        delete [] Offsets;  delete [] ByteCounts;
        return( -14 ); }
    offset += sizeof(long32);

  /* next write the 32 bit floating point parameters as a 
    one line image */

    /* one strip and one row per strip */   
    nrows = 1;

    /* set up easy tags */

    ifd2[1].tag = 256;  /* ImageWidth */
    ifd2[1].type = 4;
    ifd2[1].value = (long) PMAX;

    ifd2[2].tag = 257;  /* ImageLength */
    ifd2[2].type = 4;
    ifd2[2].value = 1L;

    ifd2[3].tag = 258;  /* BitsPerSample */
    ifd2[3].type = 3;
    ifd2[3].value = 8 * sizeof(float);

    ifd2[8].tag = 278;  /* RowsPerStrip */
    ifd2[8].type = 4;
    ifd2[8].value = nrows;

    if( CompByteOrder == tBigEndian ) { /* word swap short values in IFD*/
       i = 3;
       if( (ifd2[i].type == 3) && (ifd2[i].length <= 1 ) ) {
        stemp = (usshort*) &ifd2[i].value;
        st = stemp[0];  stemp[0] = stemp[1];
        stemp[1] = st;
       }
    }

    StripsPerImage = 1;

    /* continue writing where we left off - after the second IFD */
    /* it should be positioned OK from 8 bit write above */

    count = PMAX;
    strip = 0;
    param[0] = (float) npix;    /* make sure this is right */
    tfp2.write( (char*) param, count*sizeof(float) );
    if( tfp2.fail() ) {
       delete [] Offsets;  delete [] ByteCounts;
       tifferr("floatTIFF::write cannot write image data");
       return( -15 ); }
    Offsets[ 0 ] = offset;
    offset += count*sizeof(float);
    ByteCounts[ 0 ] = count*sizeof(float);

    ifd2[6].tag = 273;  /* StripOffsets */
    ifd2[6].type = 4;
    ifd2[6].length = 1;
    ifd2[6].value = Offsets[0];

    ifd2[9].tag = 279;  /* StripByteCounts */
    ifd2[9].type = 4;
    ifd2[9].length = 1;
    ifd2[9].value = ByteCounts[0];

    nwrite = 0;
    xy[0] = (long32) ( (double) nx / 0.04 );  /* 4 inch wide image */
    xy[1] = 100;
    tfp2.write( (char*) xy, 2*sizeof(long32) );
    nwrite += 2;
    ifd2[10].tag = 282; /* XResolution */
    ifd2[10].type = 5;
    ifd2[10].value = offset;
    offset += 2*sizeof(long32);

    //??? dx = param[14];  /* determine the aspect ratio */
    //??? dy = param[15];
    if( (dx > 0.0F) && (dy > 0.0F) )
        xy[0] = (long32) ( (dy/dx)* ((double)nx) / 0.04 );
    xy[1] = 100;
    tfp2.write( (char*) xy, 2*sizeof(long32) );
    nwrite += 2;
    ifd2[11].tag = 283; /* YResolution */
    ifd2[11].type = 5;
    ifd2[11].value = offset;
    offset += 2*sizeof(long32);

    /* use the ResolutionUnit tag for the date/time */
    caltime = time( NULL );
    mytime = localtime( &caltime );
    char ctemp[100];
    strftime( ctemp, 20, "%Y:%m:%d %H:%M:%S", mytime );
    DateTime = string(ctemp);
    ifd2[12].tag = 306; /* DateTime */
    ifd2[12].type = 2;  /* ASCII */
    ifd2[12].value = offset;
    tfp2.write( ctemp, 20 );
    nwrite += 20;
    offset += 20*sizeof(char);

    /* remember to store IFD on a word boundary */
    if( (offset/2)*2 != offset ) {
        uctemp = 0;
        tfp2.write( (char*) (&uctemp), 1 );
        offset += 1; }

    offsetIFD3 = offset;
    version = ntags;    /*  write the second IFD */
    tfp2.write( (char*) (&version), sizeof(short) );
    nwrite += 1;
    tfp2.write( (char*) ifd2, ntags * sizeof(IFD) );
    nwrite += ntags;
    ltemp = 0;
    tfp2.write( (char*) (&ltemp), sizeof(long32) );
    nwrite += 1;
    if( tfp2.fail() || ( nwrite != ( ntags + 26 ) ) )
       {    tifferr("floatTIFF::write cannot write IFD3.");
        delete[] Offsets;  delete [] ByteCounts;
        return( -16 ); }

  /*  all done write the header */

    tfp2.seekp( 0, ios_base::beg);
    nwrite = 0;
    tfp2.write( (char*) (&CompByteOrder), sizeof(short) );
    nwrite += 1;
    version = 42;
    tfp2.write( (char*) (&version), sizeof(short) );
    nwrite += 1;
    tfp2.write( (char*) (&offsetIFD1), sizeof(long32) );
    nwrite += 1;
    if( tfp2.fail() || ( nwrite !=  3) )
       {    tifferr("floatTIFF::write cannot write header.");
        delete[] Offsets;  delete [] ByteCounts;
        return( -17 ); }

   /*  fix up the IFD pointers */

    tfp2.seekp( offsetIFD1to2, ios_base::beg);
    tfp2.write( (char*) (&offsetIFD2), sizeof(long32) );
    if( tfp2.fail() )
       {    tifferr("floatTIFF::write cannot write header.");
        delete [] Offsets;  delete [] ByteCounts;
        return( -18 ); }

    tfp2.seekp( offsetIFD2to3, ios_base::beg);
    tfp2.write( (char*) (&offsetIFD3), sizeof(long32) );
    if( tfp2.fail() )
       {    tifferr("floatTIFF::write cannot write header.");
        delete [] Offsets;  delete [] ByteCounts;
        return( -19 ); }

   /* close and exit */
    tfp2.close();
    if(  tfp2.fail() ) {
        tifferr( "floatTIFF::write cannot close file.");
        delete [] Offsets;  delete [] ByteCounts;
        return( -20 ); }

    delete [] Offsets;  delete [] ByteCounts;
    return( 1 );

#undef ntags

}  /* end of floatTIFF::write() */

//--------------------- zeroParam() ----------------------------------

void floatTIFF::zeroParam( )
{
    for( int i=0; i<PMAX; i++) param[i] = 0.0F;
    return;

}  // end floatTIFF::zeroParam()
