/*
      *** psave.cpp ***

------------------------------------------------------------------------
Copyright 2014-2018 Earl J. Kirkland

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


#include "psave.hpp"    // class definition + inline functions here


//========================  public functions ===========================

//------------------ constructor --------------------------------
psave::psave()
{
    pnumMin = pnumMax = pCS;   // start with a known value in range

    //  list of parameter names and their indicies in param[]
    //  only include param. not generated from the spec. coord. 
    //  and not calculated from other param.

    setup( "xBtilt", pXBTILT );
    setup( "yBtilt", pYBTILT );

    setup( "xCtilt", pXCTILT );
    setup( "yCtilt", pYCTILT );

    setup( "defocus", pDEFOCUS );
    setup( "astig", pASTIG );
    setup( "theta", pTHETA );

    setup( "kev", pENERGY );

    setup( "oapert", pOAPERT );

    setup( "Cs3", pCS );

    setup( "capert", pCAPERT );

    setup( "ddf", pDDF );

    setup( "deltaz", pDELTAZ );

    setup( "minDet", pMINDET );
    setup( "maxDet", pMAXDET );

    setup( "phiminDet", pPMINDET );
    setup( "phimaxDet", pPMAXDET );

    setup( "Cs5", pCS5 );

    setup( "C12a", pC12a );
    setup( "C12b", pC12b );

    setup( "C21a", pC21a );
    setup( "C21b", pC21b );
    setup( "C23a", pC23a );
    setup( "C23b", pC23b );

    setup( "C32a", pC32a );
    setup( "C32b", pC32b );
    setup( "C34a", pC34a );
    setup( "C34b", pC34b );

    setup( "C41a", pC41a );
    setup( "C41b", pC41b );
    setup( "C43a", pC43a );
    setup( "C43b", pC43b );
    setup( "C45a", pC45a );
    setup( "C45b", pC45b );

    setup( "C52a", pC52a );
    setup( "C52b", pC52b );
    setup( "C54a", pC54a );
    setup( "C54b", pC54b );
    setup( "C56a", pC56a );
    setup( "C56b", pC56b );

    setup( "Nx", pNX );
    setup( "Ny", pNY );

    setup( "Nxprb", pNXPRB );
    setup( "Nyprb", pNYPRB );

    setup( "Nxout", pNXOUT );
    setup( "Nyout", pNYOUT );

    setup( "Nwobble", pNWOBBLE );

    setup( "temp", pTEMPER );

    setup( "capertMin", pCAPERTMIN );

    setup( "source", pSOURCE );

    setup( "prbPosX", pPPOSX );
    setup( "prbPosY", pPPOSY );

    setup( "oapertMin", pOAPMIN );

    setup( "probeI", pPROBEI );
    setup( "probeDt", pPROBEDT );
    setup( "dose", pDOSE );

    setup( "cdf", pCDF );

    setup( "c2astig", pCDFA2 );
    setup( "c2phi", pCDFA2PHI );

    setup( "c3astig", pCDFA3 );
    setup( "c3phi", pCDFA3PHI );

    setup( "cCs3", pCCS3 );
    setup( "cCs5", pCCS5 );

    setup( "ccapertMin", pCCAPMIN );
    setup( "ccapertMax", pCCAPMAX );

    return;

}  // end psave()

//------------------ destructor ---------------------------------
psave::~psave()
{  
    pnames.clear(); //  should be automatic
    pnumbers.clear();
    name2pvalue.clear();

}  //  end ~psave()

//------------------ load() ---------------------------------
//   load the parameters from a file
int psave::load( const char fileName[], vectorf &p, int np )
{ 
    int i;
    float x;
    string s, serr;
    ifstream fp;

    fp.open( fileName );
    if( fp.fail() )  return( -1 );

    do{ 
        fp >> s >> x ;
        i = name2pvalue[ s ];
        if( (i>0) && (i<np) && (i>= pnumMin) && (i<= pnumMax) )
            p[i] = x;
        else {
            serr= "unrecognized psave keyword= "+ s;
            messagePS( serr, 0 );
        }
    } while ( !fp.eof() );   //  go till the "end of file"

    fp.close();

    return( +1 );

}  //  end load()

//------------------ messagePS() ---------------------------------
//   common error message handler
void psave::messagePS( const string &smsg, int level )
{ 
    messageSL( smsg.c_str(), level );  //  just call slicelib version for now
    return;

}  //  end messagePS()


//------------------ save() ---------------------------------
//   save the parameters to a file
int psave::save( const char fileName[], vectorf &p, int np )
{ 
    int i, ip, n;
    ofstream fp;

    fp.open( fileName );
    if( fp.fail() )  return( -1 );

    n = (int) pnames.size();
    for( i=0; i<n; i++) {
        ip = pnumbers[i];
        if( (ip > 0 ) && (ip < np ) )
            fp << pnames[i] << "    " << p[ip] << endl;
    }

    fp.close( );
    if( fp.fail() )  return( -2 );

    return( +1 );

}  //  end save()

//------------------ setup() ---------------------------------
//   setup param labels and table
//   call once per entry
void psave::setup( const char name[], int pValue )
{   
    pnames.push_back( string(name) );
    pnumbers.push_back( pValue );

    name2pvalue[ name ] = pValue;

    if( pValue < pnumMin ) pnumMin = pValue;
    if( pValue > pnumMax ) pnumMax = pValue;

    return;

}  //  end setup();