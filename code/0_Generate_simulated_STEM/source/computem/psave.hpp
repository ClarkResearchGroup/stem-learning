/*
      *** psave.cpp ***

------------------------------------------------------------------------
Copyright 2014-2019 Earl J. Kirkland

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

   started 12-jan-2014 ejk
   both save and load working 8-feb-2014 ejk
   add dose parameters pPROBEI, pPROBEDT, pDOSE 10-oct-2015 ejk
   convert param[] to vector<> 23-jul-2016 ejk
   add pPMINDET,pPMAXDET for seg. STEM detector 12-may-2018
   last modified 12-may-2018 ejk

*/

#ifndef PSAVE_HPP   // only include this file if its not already
#define PSAVE_HPP   // remember that this has been included


#include "slicelib.hpp"    // misc. routines for multislice

#include <cstdio>  /* ANSI C libraries */
#include <cstdlib>
#include <cstring>
#include <cmath>

#include <string>   // STD string class
#include <vector>   // STD vector class
#include <map>      // STD map class
#include <fstream>  // STD file IO streams
#include <iostream>

using namespace std;


//------------------------------------------------------------------
class psave{

public:
        
    psave();    // constructor function
    ~psave();   // destructor function

    int save( const char fileName[], vectorf &p, int n );  //  save params to file
    int load( const char fileName[], vectorf &p, int n );  //  load params from file

private:
    int pnumMin, pnumMax;
    vector<string>  pnames;
    vector<int> pnumbers;

    map<string, int> name2pvalue;

    void setup( const char name[], int pValue );

    void messagePS( const string &smsg, int level );

}; // end psave::

#endif