/*      incostem.hpp

------------------------------------------------------------------------
Copyright 1998-2018 Earl J. Kirkland

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

  header file for incostem.cpp
  class with member subroutines to
  calculate images in the incoherent STEM approximation

  put the partial cross section (integrated over the ADF detector
  angles) at a single pixels (position of the corresponding atom)
  and convolve with the point spread function (focused probe intensity)

  reference:

  [1] E. Kirkland, "Advanced Computing in Electron Microscopy",
        Plenum 1998, 2nd edit. Springer 2010
  
  this file is formatted for a tab size of 4 char

   started separate incostem class .hpp file  c20-apr-2013 E. Kirkland
   consolodate makeProbe() and prbSize() in probe.cpp and call it from
      here to avoid duplicating code 05-jul-2013 ejk
  convert message() to use string data 5-sep-2013 ejk
  add Pnoise() 29-oct-2014 ejk
  switch to my ranPoisson() for compilers without c++11 on 27-sep-2015 ejk
  convert to use probe::makeProbeIntensity()  which is more generally 
      useful in other places 26-jan-2016 ejk
  convert malloc1D() to vector<> 5-jul-2016 ejk
  add parameter verbose to turn off unneccesary
        terminal echo 26-nov-2016 ejk
  add probe position=(0,0) in param[] for new probe::makeProbeINtensity()
       31-jan-2018 ejk
  change how probe::makeProbeIntensity() calculates the probe size
        3-feb-2018 ejk
*/

#ifndef INCOSTEM_HPP   // only include this file if its not already

#define INCOSTEM_HPP   // remember that this has been included

#include <cstdio>  /* standard ANSI libraries */
#include <cstdlib>
#include <cmath>
#include <ctime>   // to init RNG seed
#include <vector>

//#include <random> // for Poission RNG: requires -std=c++11 option in gcc/g++

using namespace std;

#include "cfpix.hpp"       // complex image handler with FFT
#include "slicelib.hpp"    // misc. routines for multislice

//------------------------------------------------------------------
class incostem{

public:
    
    incostem( );        // constructor functions
    ~incostem();        //  destructor function

    //void calculate( cfpix &pix, float param[], int multiMode, int natom,
    //    int Znum[], float x[], float y[], float z[], float occ[], float wobble[] );

    void calculate2D( cfpix &pix, vectorf &param, int multiMode, int natom,
        vectori &Znum, vectorf &x, vectorf &y, vectorf &occ );
    
    int addNoise( cfpix &pix, int nx, int ny, double probeI, double dwellTime );

    void quiet() { echo = 0 ; };          // turn off unneeded terminal output
    void verbose() { echo = 1; };         // turn on unneeded terminal output

private:

    int NZMAX, FCNatomf, FCNfemr, FCNfemi, echo;

    double twopi, sigmae, wavl;

    double atomf( double t, double p[] );

    double atomsignal( int zatom, double keV, double thetamin, double thetamax );

    double BJ0( double x );

    void feMoliere( double k, int zatom, double *rfe, double *ife );

    double femi( double r, double p[] );

    double femr( double r, double p[] );

    double fint( int FCN, double r, double p[] );

    //double integrate45( double (*fint)(double x, double p[]), double p[],
    double integrate45( int fcn, double p[],
        double xmin, double xmax, double maxerror, int maxsteps );

    void invert2D( float** pix, long nx, long ny );

    void messageIN( std::string &smsg, int level = 0 );  // common error message handler

    std::string sbuff;
    
    int Pnoise( cfpix &pix, cfpix &pixout, int nx, int ny, int imean );

    unsigned long iseedp;

}; // end incostem::

#endif
