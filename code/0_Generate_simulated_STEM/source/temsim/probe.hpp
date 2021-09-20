/*      probe.hpp

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

  header file for probe.cpp
  class with member subroutines to
  calculate STEM coherent probe wavefunction

  put the partial cross section (integrated over the ADF detector
  angles) at a single pixels (position of the corresponding atom)
  and convolve with the point spread function (focused probe intensity)

  reference:

  [1] E. Kirkland, "Advanced Computing in Electron Microscopy",
        Plenum 1998, 2nd edit. Springer 2010
  
  this file is formatted for a tab size of 4 char

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

#ifndef PROBE_HPP   // only include this file if its not already

#define PROBE_HPP   // remember that this has been included


#include <cstdio>  // standard ANSI libraries
#include <cstdlib>
#include <cmath>

#include <string>
#include <vector>

using namespace std;

#include "cfpix.hpp"       // complex image handler with FFT
#include "slicelib.hpp"    // misc. routines for multislice

//------------------------------------------------------------------
class probe{

public:
    
    probe( );         // constructor functions
    
    ~probe();        //  destructor function

    void abbPhase2D( cfpix &ab2D,vectorf &param, int multiMode );

    int makeProbe(  cfpix &cpix, int nx, int ny, double xp, double yp,
        vectorf &p, double wavlen, double k2max, double pixel, int multiMode,
        int ismoth, vectorf &kx, vectorf &kx2, vectorf &ky, vectorf &ky2 );

    int makeProbeIntensity( cfpix &pix, vectorf &param, int multiMode,
        int initFFT );

    double prbSize( cfpix &pixsq, int nx, int ny,
        double xp, double yp, double ax, double by );

private:

    void messagePR(  std::string &smsg,  int level = 0 );  // common error message handler

}; // end incostem::

#endif
