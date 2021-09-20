/*
        *** newD.hpp ***

------------------------------------------------------------------------
Copyright 2017 Earl J. Kirkland

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

    template library for allocating/deallocating simple multidimensional arrays
    (should use stl vector<> for 1D arrays)
    
    delete2D()  : free 2D array allocated with malloc2D()
    delete3D()  : free 3D array allocated with malloc3D()
    new2D()     : allocate memory for 2D arrays
    new3D()     : allocate memory for 3D arrays

    this file is formatted for a tab size of 4 characters

    started 27-aug-2017 E. Kirkland
    change messageME() to inline otherwise it may be multiply defined
       if this is included more than once  24-sep-2017 ejk
    change memory.hpp to newD.hpp because C++11 already has a memory
        header  12-nov-2017 ejk
*/

#ifndef MEMORY_HPP   // only include this file if its not already

#define MEMORY_HPP   // remember that this has been included


//
#define USE_SLICELIB   //  do not define to run without slicelib.cpp

#ifdef USE_SLICELIB
#include "slicelib.hpp"    // misc. routines for multislice 
#else
#include <string>
#include <iostream>  //  C++ stream IO
#include <iomanip>   //  to format the output

using namespace std;
#endif

#ifdef  ALTERNATE_MACRO   //  save just in case inline doesn't work somewhere
//Remember:  if this is a subroutine it will produce multiply
//      defined symbols messageME() if included in >1 place!
#ifdef USE_SLICELIB
#define messageME(smsg,level)  messageSL( smsg.c_str(), level )  //  just call slicelib version for now
#else
#define messageME(smsg,level)   cout << smsg << ", level= " << level << endl
#endif
#endif


/*--------------- delete2D() ----------------------------*/
/*   
    free a 2D array f[][] dimensioned as f[ix][iy]
    
    0 < ix < (nx-1)
    0 < iy < (ny-1)   <not used>
    
*/
template< class T > 
void delete2D( T ** f, int nx )
{
    int ix;
    
    for( ix=0; ix<nx; ix++) delete [] f[ix];
    delete [] f ;
    
    return;
    
}  // end delete2D() 

/*--------------- delete3D() ----------------------------*/
/*   
    free a 2D array f[][] dimensioned as f[ix][iy]
    
    0 < ix < (nx-1)
    0 < iy < (ny-1)
    0 < iz < (nz-1) <not used>
    
*/
template< class T > 
void delete3D( T ***f, int nx, int ny )
{
    int ix, iy;
    
    for( ix=0; ix<nx; ix++) {
        for( iy=0; iy<ny; iy++) delete [] f[ix][iy];
        delete [] f[ix];
    }
    delete [] f;
    
    return;
    
}  // end delete3D() 

/* -------------------  messageME() -------------------

Remember:  if this is not "inline" it will produce multiply
        defined symbols messageME() if included in >1 place!

   message output
   direct all output message here to redirect to the command line
   or a GUI status line or message box when appropriate

   msg[] = character string with message to disply
   level = level of seriousness
        0 = simple status message
        1 = significant warning
        2 = possibly fatal error
*/
inline void messageME( std::string &smsg,  int level )
{
#ifdef USE_SLICELIB
        messageSL( smsg.c_str(), level );  //  just call slicelib version for now

#else
        cout << smsg << ", level= " << level << endl;
#endif

}  // end autoslic::messageME()

/*---------------------------- new2D() -------------------------------*/
/*
    2D array allocator for type T (float, double...)
    make space for m[0...(nx-1)][0..(ny-1)]

    nx,ny = number of elements
    T = data type = double, float, int etc.

*/
template< class T > 
T **new2D( int nx, int ny, const char *message )
{   T **m;
    int i;
    std::string stemp;

    m = new T* [nx];
    if( NULL == m ) {
        stemp= "out of memory in new2D(), size=" + toString(nx) +
               ", "+toString(ny)+": "+message;
        messageME( stemp, 2 );
        exit(0);
    }

    for (i=0; i<nx; i++){
        m[i] = new T[ny];
        if( NULL == m[i] ){
            stemp= "out of memory in new2D(), size=" + toString(nx) +
                   ", "+toString(ny)+": "+message;
            messageME( stemp, 2 );
            exit(0);
        }
    }

    return m;

}  // end new2D() 

/*---------------------------- new3D() -------------------------------*/
/*
    3D array allocator for numeric data
    make space for m[0...(nx-1)][0..(ny-1)][0..(nz-1)]

    nx,ny,nz = number of elements
    T = data type = double, float, int etc.

*/
template< class T > 
T ***new3D( int nx, int ny, int nz, const char *message )
{   T ***m;
    int i, j;
    std::string stemp;

    m = new T**[ nx ];
    if( NULL == m ) {
        stemp= "out of memory in new3D(), size=" + toString(nx) + ", " +
                toString(ny) + ", " + toString(nz) + ": " + message;
        messageME( stemp, 2 );
        exit(0);
    }

    for (i=0; i<nx; i++){
        m[i] = new T*[ny];
        if( NULL == m[i] ){
            stemp= "out of memory in new3D(), size=" + toString(nx) + ", " +
                    toString(ny) + ", " + toString(nz) + ": " + message;
            messageME( stemp, 2 );
            exit(0);
        }
        for (j=0; j<ny; j++){
            m[i][j] = new T[nz];
            if( NULL == m[i][j] ){
                stemp= "out of memory in new3D(), size=" + toString(nx) + ", " +
                       toString(ny) + ", " + toString(nz) + ": " + message;
                messageME( stemp, 2 );
                exit(0);
            }
        }
    }

    return m;

}  // end new3D()

#endif


