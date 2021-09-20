/*              *** floatTIFF.hpp ***

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

   header file for floatTIFF.cpp

   C++ class to read/write floating point images as TIFF format files (in ANSI-C++)
   and manage the memory storage.  There are three (3) images in each file.
   The first image is a simple 8 bit image with square pixels for display only.
   The second image is a 32 bit floating point image with possibly rectangular
   pixels for calculations.  The third image is a 2048 by 1 "image" of
   32 bit floating point parameters related t the calculation.

   This is mainly an easy to use high level front end derived from the old
   tiffsubs.c non-class library.

   The computer is assumed to have either big-endian or little-endian
   byte ordering and to be a 32 bit machine (i.e. a long must be 32 bits).

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

maxParam()      : return maximum nu,ber of parameters

nx()            : return horz. size of current image (in pixels)
ny()            : return vert. size of current image (in pixels)

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
   work on write(), min()/max() subroutines 11-mar-2012 ejk
   add maxParam() member 21-mar-2012
   fix typedef IFD to be long32 for 64 bit gcc 31-may-2012
   convert sprintf()+char[] to strings 1-feb-2014 ejk
   convert to streams IO 1-feb-2014 ejk
   convert remaining buffers to use new/delete 2-feb-2014 ejk
   update wxGUI error message in tifferr() 10-feb-2014 ejk
   last modified 10-feb-2014 ejk

*/

#ifndef FLOATTIFF_HPP   // only include this file if its not already

#define FLOATTIFF_HPP   // remember that this has been included

#include <string>   // STD string class
#include <sstream>  // string streams
#include <fstream>  // STD file IO streams
#include <iostream>

using namespace std;


/*  choose one of each so that a long32 is a 32 bit integer */
/* #define long32 long;     for 16/32 bit machines - rarely used now */
#define long32 int  /*  for 32/64 bit machines */

// define the following symbol to enable bound checking
//     define for debugging, and undefine for final run
//  because these subroutines are not called often it probably doesn't
//  hurt to leave this on all the time
#define floatTIFF_BOUNDS_CHECK


// ---- short hand for some data types
#ifndef TIFF_TYPDEF   // allow TIFFimage and floatTIFF to coexist
#define TIFF_TYPDEF
typedef unsigned short usshort;
typedef unsigned long32 uslong;
typedef struct { unsigned short tag, type;
                unsigned long32 length, value; } IFD;
#endif
        
//------------------------------------------------------------------
class floatTIFF{

public:
    
    // constructor functions
    floatTIFF();            // blank image
    
    ~floatTIFF();    //  destructor function

    int tFloatTest( const char *filename );

    inline int getnpix() const { return( npix ); }
    void setnpix( int i );

    float getParam( int i );
    void setParam( const int i, const float val );
    void zeroParam();
    inline string getDateTime( ) { return DateTime; }

    float max( int ipix );
    float min( int ipix );

    inline int nx() const { return( int(nxl) ); }
    inline int ny() const { return( int(nyl) ); }
    inline int maxParam() const { return( PMAX ); }

    //  remember: operator[] only allows one argument so can't be used for > 1D
    inline float& operator()( const int i1, const int i2 )
    {
        #ifdef floatTIFF_BOUNDS_CHECK
            if( (i1<0) || (i1>=int(nxl)) ||
                (i2<0) || (i2>=int(nyl)) ){
                printf( "out of bounds index in floatTIFF\n"
                    "  size = %d x %d\n  access = (%d , %d)\n", 
                    nxl, nyl, i1, i2 );
                exit( EXIT_FAILURE );
            }
        #endif
        // both should work but one may be faster 
        //   for different operations
        //  return *(data + i2 + i1*nyl);
            return *(data + i1 + i2*nxl);
    }

    int read( const char *filname );

    int resize( const int nx, const int ny );

    int write( const char *file, float rmin, float rmax, float imin, float imax,
    float dxi, float dyi );

private:

    short CompByteOrder;    //  will get the byte order type for this computer
    ifstream tfp;           // pointer to TIFF file
    short FileByteOrder;    // byte order of the TIFF file
    short Version;          // TIFF 'version', must be 42
    uslong NextIFD;         // pointer to next Image File Directory
    uslong CurrentIFD;      // pointer to current Image File Directory
    short nIFD;             // number of entries in IFD
    int nstrips;            // number of strips in file
    int nstripcounts;       // number of strip byte counts in file
    IFD* ifd1;              // IFD data 

    uslong nxl, nyl;        // current stored image size
    float *data;            // current image data buffer
    int npix;               // number of images (1 for real, and 2 for complex)
    float *param;           // buffer for image parameters
    int PMAX;               // max. number of parameters

    usshort* stempo;    // temporary short arrays for IFDread 
    usshort* stempc;

    //------------- TIFF parameters ------------------------------------

    usshort BitsPerSample[3], Compression, FillOrder,
        NewSubFileType,
        MaxSampleValue, MinSampleValue,
        Orientation, PhotometricInterpretation,
        PlanarConfiguration, Predictor, ResolutionUnit,
        SamplesPerPixel, StripByteCount,
        SubfileType, SampleFormat;

    uslong ImageLength, ImageWidth, NewSubfileType,
        RowsPerStrip, XResolution[2], YResolution[2],
        *StripByteCounts, *StripOffsets,
        XPosition[2], YPosition[2];

    string DateTime;

    //----------------- private functions ---------------------
    int tclose( );
    
    void tifferr( const string &error_text );

    int tifftest( const char *filename );

    float tinter( long32 nx, long32 ny, double x, double y );

    int topenFloat( const char *filename );

    int tread( void *bufptr, int size, int n, long32 offset );

    int treadFloatPix( );

    int treadIFD( uslong IFDoffset );
        
    void tsetByteOrder();

};

#endif  // FLOATTIFF_HPP
