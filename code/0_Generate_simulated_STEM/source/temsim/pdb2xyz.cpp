/*  pdb2xyz.cpp
    
    short program to convert PDB (protein data base) files
    to temsim/computem xyz format file - NOT guaranteed to work on all PDB files

------------------------------------------------------------------------
Copyright 2009-2019 Earl J. Kirkland

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
        
    started 2-nov-2009 E. Kirkland
    add carbon support with sub. from biocells.cpp 8-nov-2009
    convert to strings and streams 18-nov-2014 ejk
    convert to vector<> 24-nov-2014 ejk
    add rotation 6-may-2015 ejk
    add option to make just a-carbon film 4-sep-2016 ejk
    convert malloc2D() to vector<>  9-sep-2016 ejk
    change RMIN from 1.0 to 1.5 Ang (C-C bond) 10-sep-2016 ejk
    minor changes to make gcc happy 11-may-2018 ejk
    fix small error in scaling after rotation 7-oct-2018 ejk
    add trap for bad number C atoms 8-oct-2018 ejk
    start option for ICE=H2O  1-may-2019, working 13-may-2019  ejk
	fix minor formatting issue with progress echo 15-may-2019 ejk
	change to ICE_FILL_FRACTION=0.8 on 2-jun-2019 ejk
	put back ICE_FILL_FRACTION=1.0 - was not the problem 6-jun-2019 ejk
    last updated 6-jun-2019 ejk
*/

#include <cstdlib>
#include <cmath>
#include <ctype.h>   // for toupper() 
#include <string.h>  // for strings
#include <time.h>

#include <string>
#include <iostream>  //  C++ stream IO
#include <sstream>   // string streams
#include <fstream>
#include <iomanip>   //  to format the output
#include <vector>

using namespace std;

#include "slicelib.hpp"

const int TRUE=1;
const int FALSE=0;

//  CFILL_FRACTION=1 makes flat support with almost no structure
//      make <1 to get a little surface roughness
const double CWEIGHT=12.01115;       // molec. weight of carbon
const double CDENSITY=2.0;           // density in gm/cm^3 approx. for amorphous C
const double CFILL_FRACTION=0.9;     // filling fraction for this density

const double H2O_WEIGHT=18.0152;     // molecular weight for water (in gm/mole) 
const double H2O_DENSITY=1.0;        // density in gm/cm^3 

//  bio molecules are not visible with ICE_FILL_FRACTION=1 and Scherzer focus
//    must make extremely large defocus (FF=0.8 doesn't help)
const double ICE_WEIGHT=18.0152;     // molecular weight for water (in gm/mole) 
const double ICE_DENSITY= 0.9167;    //  density in gm/cm^3
const double ICE_FILL_FRACTION=1.0;  // filling fraction for this density

const double NAV=6.0225e23;     // Avagadro's number (#/mole)
const double RMIN=1.4;          // min separation distance of C (in Angstroms)

// index for atomic coordinates (x,y,z)= coord, occ=occupancy, znum= atomic number
enum{ OX=0, OY=1, OZ=2, Oocc=3, Oznum=4 };
const int NVAL=5;

unsigned long iseed;        // random number generator seed

/*  subroutines define at the end of the file */
void newcoord( vector<double> &c, 
    const double xmax, const double ymax, const double zmax );
void insertcoord( vector<double> &ctest, vector< vector<double> > &coord, 
        const int ncoord, const int pos );
int testcoord( vector<double> &ctest, vector< vector<double> > &coord, 
    const int ncoord, int *pos, const double rmin );
int fillSolid( vector<vector<double> > &coord, const int ncoord, 
    double ax, double by, double cz, double rmin );

//  more subroutines for ice
int fillIce( vector< vector<double> > &coord, const int nMolIce, const int np,
    double ax, double by, double cz, 
    vector<double> &pdbmin, vector<double> &pdbmax, double rmin );
void sortByZd( vectord &x, vectord &y, vectord &z, vectord &occ,
    vectori &Znum, int natom );
int testInside( vector<double> &ctest,
    vector<double> &pdbmin, vector<double> &pdbmax );
void hydrocoord( vector<double> &ctest, vector<double> &coord );

// modes of operation
enum{ C_PDB_MODE=0, ICE_PDB_MODE=1, C_ONLY_MODE=2, PDB_ONLY_MODE=3 };
int mode=0;


int main()
{
    string cline, filin, filout, symb;
    int i, ic, iz,  np, nlen, ntotal, nAtomsIce, nMolIce, natom, ctop, lreadPDB;
    long ltime;
    double xpos, ypos, zpos, occ, mytime, x0, y0, z0, cr,sr,ct,st;
    double xmin, xmax, ymin, ymax, zmin, zmax, rotat, tilt, pi;
    double ax, by, cz, cthick, xoff, yoff, zoff, wobble, density, vol;

    double iceThick;

    vector<int> znum;
    vector<double> xp, yp, zp, oc;      //  will get atom positions

    vector<double> pdbmin(NVAL), pdbmax(NVAL);

    ifstream fpin;
    ofstream fpout;
    
    /*  the following are the chemical symbols for the periodic table */
    char symbol[] = {
                " HHeLiBe B C N O FNeNaMgAlSi P SCl"
                "Ar KCaScTi VCrMnFeCoNiCuZnGaGeAsSeBr"
                "KrRbSr YZrNbMoTcRuRhPdAgCdInSnSbTe"
                " IXeCsBaLaCePrNdPmSmEuGdTbDyHoErTm"
                "YbLuHfTa WReOsIrPtAuHgTlPbBiPoAtRn"
                "FrRaAcThPa UNpPuAmCmBkCfEsFmMdNoLr"
        };

    //------- what this is ---------- 
    
    cout << "pdb2xyz version dated 6-jun-2019" << endl;
    cout << "  convert PDB data file to xyz format, with a-carbon or a-ice support\n" << endl;
    cout << "  (may not work on all PDB files)\n" << endl;

    //------- convert symbols to upper case for PDB comparison ---------- 
    
    nlen = (int) strlen( symbol );
    for( i=0; i<nlen; i++) symbol[i] = toupper( symbol[i] );
    
    // set defaults
    lreadPDB = true;
    mode = ICE_PDB_MODE;

    //----------- get file names and open ---- 

    cout << "Type name of input file with PDB data ('none' for just a-carbon):" << endl;
    cin >> filin;
    if( filin.compare( "none" ) == 0  ) {
        lreadPDB = FALSE;
        mode = C_ONLY_MODE;
    }

    if( TRUE == lreadPDB ) {
        fpin.open( filin.c_str() );
        if( fpin.fail() ) { cout << "Can't open file "<< filin << endl; exit(0); }
    }

    cout << "Type name of output file to get xyz data:" << endl;
    cin >>  filout;
    fpout.open( filout.c_str() );
    if( fpout.bad() ) { cout << "Can't open file" << filout<<endl; exit(0); }

    if( TRUE == lreadPDB ) {
        cout << "Type azimuthal and polar rotation angles (in degrees):" << endl;
        cin >> rotat >> tilt;
        pi = 4.0 * atan( 1.0 );
        rotat = rotat * pi/180.0;
        tilt  = tilt * pi/180.0;
    }

    //----------- get carbon support or ice support info ---- 
    //    a little awkward - should do better UI sometime
    if( TRUE == lreadPDB ) {
        cout << "Type thickness of carbon support (<0 to disable):" << endl;
        cin >> cthick;
        if( cthick > 0.0 ){
            cthick = fabs( cthick);
            if( (cthick > 0.0) && (TRUE == lreadPDB) ) 
                ctop = askYN( "Do you want the carbon support on the entrance instead of exit");
            iceThick = 0.0;
            mode = C_PDB_MODE;
        } else{
            cout << "Type thickness of ice support (<0 to disable):" << endl;
            cin >> iceThick;
            if( iceThick > 0.0 ) mode = ICE_PDB_MODE;
            else mode = PDB_ONLY_MODE;
        }   
    }

    if( C_ONLY_MODE == mode ) {
        cout << "Type size of amorphous carbon x,y,z (in Ang.):" << endl;
        cin >> ax >> by >> cthick;

        if((ax < 0) || (by < 0) || ( cthick < 0) ) {
            cout << "need > 0 if not reading PDB file, cannot continue...." << endl;
            return( EXIT_SUCCESS );
        }
    }

    cout << "mode = " << mode << endl;  // diagnostic ?????

    mytime = cputim();   //  get CPU time for fun

    // ---- initialize random number generator seed ---- 
    if( (cthick > 0.0) || (iceThick > 0.0) ) {
        ltime = (long) time( NULL );
        iseed = (unsigned) ltime;
        if( ltime == -1 ) {
            cout << "Type initial seed for random number generator:" << endl;
            cin >> iseed;
        } else {
            cout <<  "Random number seed initialized to " << iseed << endl;
        }
    }

    //--------- read data from pdb file in complicated format ------------------------------------
    //  separate read operation so PDB may be used in different modes

    if( TRUE == lreadPDB ) {
        //   read in atoms and get total range and number of atoms 
        cout << "read input PDB file " << filin << endl;
        np = 0; 
        do {
            getline( fpin, cline ); // read a whole line

            //  select lines begin. with ATOM or HETAM  with atom coord.
            if( ( cline.find( "ATOM") == 0 ) ||
                ( cline.find( "HETATM") == 0 )  )  {
            
                //----  read x,y,z coordinates
                istringstream isbuf( cline.substr(30) );
                isbuf >> xpos >> ypos >> zpos >> occ;

                if( 0 == np ) {
                    xmin = xmax = xpos;
                    ymin = ymax = ypos;
                    zmin = zmax = zpos;
                } else {
                    if( xpos < xmin ) xmin = xpos;  // coord. range
                    if( xpos > xmax ) xmax = xpos;
                    if( ypos < ymin ) ymin = ypos;
                    if( ypos > ymax ) ymax = ypos;
                    if( zpos < zmin ) zmin = zpos;
                    if( zpos > zmax ) zmax = zpos;
                }
                xp.push_back( xpos );   //  save coord.
                yp.push_back( ypos );
                zp.push_back( zpos );
                oc.push_back( occ );

                //  find atomic number
                symb= cline.substr(76,77);  // get chemical symbol
                for( i=0; i<nlen; i+=2) {
                    iz = 1 + i/2;
                    if( strncmp( &symbol[i], symb.c_str(), 2 ) == 0 ) break;
                }
                znum.push_back( iz );

                np++;
            }

            // be careful; names like HENDERSON match "END" (!?)
            //   look for "END" in col 1
        }  while( cline.find( "END") != 0 );  //  ???? == string::npos );
        fpin.close( );

        if( np < 1 ) {
            cout << "no atoms read, cannot continue...." << endl;
            return( EXIT_SUCCESS);
        }
    
        cout <<  "Total number of atoms = " << np << endl;
        cout <<  "  with x range " << xmin << " to " << xmax << endl;
        cout <<  "   and y range " << ymin << " to " << ymax << endl;
        cout <<  "   and z range " << zmin << " to " << zmax << endl;

        //----------  rotate to get better view perhaps (from slicview.cpp)  -------------------
        if( ( fabs( rotat) > 1.0e-6 ) || ( fabs( tilt ) > 1.0e-6 ) ) {
            //  Move to center of molecule and rotate
            xoff = 0.5F*( xmax + xmin );
            yoff = 0.5F*( ymax + ymin );
            zoff = 0.5F*( zmax + zmin );

            cout << "rotate about x,y,z= " << xpos << ", " << ypos << ", " << zpos << endl;

            /*  Calculate misc constants  */
            cr = cos( rotat );
            sr = sin( rotat );
            ct = cos( tilt );
            st = sin( tilt );

            for( i=0; i<np; i++) {
                xp[i] = xp[i] - xoff;   // translate to center
                yp[i] = yp[i] - yoff;
                zp[i] = zp[i] - zoff;

                x0 = xp[i];      //  Rotation
                y0 = yp[i];
                xp[i] =  xpos = cr*x0 - sr*y0;
                yp[i] =         sr*x0 + cr*y0;

                y0 = yp[i];      //  Tilt
                z0 = zp[i];
                yp[i] = ypos =  ct*y0 + st*z0;
                zp[i] = zpos = -st*y0 + ct*z0;

                //  find new range
                if( 0 == i ) {
                    xmin = xmax = xpos;
                    ymin = ymax = ypos;
                    zmin = zmax = zpos;
                } else {
                    if( xpos < xmin ) xmin = xpos;
                    if( xpos > xmax ) xmax = xpos;
                    if( ypos < ymin ) ymin = ypos;
                    if( ypos > ymax ) ymax = ypos;
                    if( zpos < zmin ) zmin = zpos;
                    if( zpos > zmax ) zmax = zpos;
                }
            }

            cout <<  "New range of atom coord. after rotation and tilt" << endl;
            cout <<  "   x range " << xmin << " to " << xmax << endl;
            cout <<  "   y range " << ymin << " to " << ymax << endl;
            cout <<  "   z range " << zmin << " to " << zmax << endl;

        }  // end rotate section

    }  //  end if( TRUE == lreadPDB) ...

    //-----------------------------  carbon support modes ----------------
    //  leave old modes alone (they were working) - but separate this section
    //    maybe combine intelligently later

    if( (C_PDB_MODE == mode) || (C_ONLY_MODE == mode) || (PDB_ONLY_MODE == mode) ) {

        if( TRUE == lreadPDB ) {

            cout << "write PDB data to xyz file " << filout << endl;

            //--------- write data to xyz file with offset -------
            ax = 1.5*( xmax - xmin );
            by = 1.5*( ymax - ymin );
            //ax = 1.4*( xmax - xmin );
            //by = 1.4*( ymax - ymin );
            if( by > ax ) ax = by;      // make it square to look right 
            else by = ax;
            if( cthick > 0.0 ) cz = (zmax - zmin) + cthick;
            else cz = zmax - zmin;
            fpout << "pdb2xyz translation of " << filin << endl;
            fpout << setw(16) << ax << setw(16) << by << setw(16) << cz << endl; // cell size 

            xoff = 0.5*ax - 0.5*(xmax+xmin);   // move molecule to the center 
            yoff = 0.5*by - 0.5*(ymax+ymin);

            if( (cthick>0.0) && (ctop==1) ) zoff = -zmin + cthick;
            else zoff = -zmin;
    
            wobble = 0.0;   //  Debye-Waller factor 
            for( i=0; i<np; i++ ) {
                fpout << setw(5) << znum[i] << setw(14) << xp[i]+xoff << setw(14)
                    << yp[i]+yoff << setw(14) << zp[i]+zoff << setw(14) << oc[i] 
                    //<< zp[i]+yoff << setw(14) << yp[i]+zoff << setw(14) << oc[i] 
                    << setw(14) << wobble << endl;
            }

        } else {  //  end if( TRUE == lreadPDB) ...

            //  write header for a-carbon only
            fpout << "pdb2xyz a-carbon" << endl;
            fpout << setw(16) << ax << setw(16) << by << setw(16) << cthick << endl; // cell size 
        }

        //----------- generate carbon support if requested -------------- 

        if( cthick > 0.0 ) {

            cout << "generate random coord. for carbon support"  << endl;
            density = CFILL_FRACTION*NAV*CDENSITY*(1.0e-24)/CWEIGHT; //  # atoms/Angs^3 
            cout << "average density = "<< CFILL_FRACTION*CDENSITY << " gm/cm^3 = "
                << density << " atoms/Angstrom^3" << endl;
            cout << "minimum allowed separation = " << RMIN << " Angstroms" << endl;

            cout << "calculate coord. in a vol. of "<< ax << " x "<< by << 
                " x "<< cthick << " Angstroms" << endl;
            ntotal = (int) ( ax*by*cthick  * density );
            cout << "Total number of carbon atoms = " << ntotal << endl;

            if( ntotal > 0 ) {
                vector<double> tmp( NVAL );
                vector< vector<double> > coord( ntotal, tmp );

                //  fill in random coord.
                ic = fillSolid( coord, ntotal, ax, by, cthick, RMIN );

                iz = 6;   //  atomic number of carbon 
                occ = 1.0;
                wobble = 0.0;
                if( ctop == 1 ) zoff = 0.0; //  on top  
                else zoff = zmax - zmin;    // on bottom 
                for( i=0; i<ic; i++) 
                    fpout << setw(5) << iz << setw(14) << coord[i][0] << setw(14)
                        << coord[i][1] << setw(14) << zoff + coord[i][2] << setw(14) << occ 
                        << setw(14) << wobble << endl;
            }  /// end if( ntotal > 0 )
        } // end if( cthick > 0 )

        //-----------  write end of file mark -------------------------------- 
        iz = -1;
        fpout << setw(5) << iz << endl;  //  end of data 
        fpout.close( );

    // end if ( C_ONLY_MODE of C_PDB_MODE

    } else if( ICE_PDB_MODE == mode ) {

        cout << "start ice mode " << endl;

        //--------- write data to xyz file with offset -------
        ax = 1.5*( xmax - xmin );
        by = 1.5*( ymax - ymin );
        //ax = 1.4*( xmax - xmin );
        //by = 1.4*( ymax - ymin );
        if( by > ax ) ax = by;      // make it square to look right 
        else by = ax;
        cz = iceThick;
        if( iceThick < (zmax-zmin) ) {
            cout << "ice thickness = " << iceThick << " is < PDB structure of " << (zmax-zmin) << endl;
            cout << "cannot continue, exiting...." << endl;
            exit( 0 );  //  rude but what else can we do...
        }
        fpout << "pdb2xyz translation of " << filin << endl;
        fpout << setw(16) << ax << setw(16) << by << setw(16) << cz << endl; // cell size 

        xoff = 0.5*ax - 0.5*(xmax+xmin);   // move molecule to the center of the ice
        yoff = 0.5*by - 0.5*(ymax+ymin);
        zoff = 0.5*cz - 0.5*(zmax+zmin);

        for( i=0; i<np; i++ ) {
            xp[i] += xoff;
            yp[i] += yoff;
            zp[i] += zoff;
            if( i == 0 ){
                xmin = xmax = xp[0];
                ymin = ymax = yp[0];
                zmin = zmax = zp[0];
            }
            if( xp[i] < xmin ) xmin = xp[i];   //  find new range
            if( xp[i] > xmax ) xmax = xp[i];
            if( yp[i] < ymin ) ymin = yp[i];
            if( yp[i] > ymax ) ymax = yp[i];
            if( zp[i] < zmin ) zmin = zp[i];
            if( zp[i] > zmax ) zmax = zp[i];
        }

        cout <<  "New range of atom coord. after move to ice center:" << endl;
        cout <<  "   x range " << xmin << " to " << xmax << endl;
        cout <<  "   y range " << ymin << " to " << ymax << endl;
        cout <<  "   z range " << zmin << " to " << zmax << endl;

        cout << "generate random coord. for ice support"  << endl;
        density = ICE_FILL_FRACTION*NAV*ICE_DENSITY*(1.0e-24)/ICE_WEIGHT; //  # atoms/Angs^3 
        cout << "average density = "<< ICE_FILL_FRACTION*ICE_DENSITY << " gm/cm^3 = "
                << density << " molecules/Angstrom^3" << endl;
        cout << "minimum allowed separation = " << RMIN << " Angstroms" << endl;
        cout << "fill fraction = " << ICE_FILL_FRACTION << endl;

        cout << "calculate coord. in a vol. of "<< ax << " x "<< by << 
                " x "<< cz << " Angstroms" << endl;

        cout << "sort PDB by z..." << endl;
        sortByZd( xp, yp, zp, oc, znum, np );

        //  allocate memory for much more than needed - enough for all ice plus PDB
        //    not efficient but easier to program (sorry?)
        nMolIce = (int) ( ax*by*cz  * density );  //  total number of ice molecules without PDB
        nAtomsIce = 3 * nMolIce;  //  total number of ice atoms without PDB
        cout << "add of order " << nMolIce << " molecules of surrounding ice" << endl;
        ntotal = nAtomsIce + np;  //  should be much more than needed  
        vector<double> tmp( NVAL );
        vector< vector<double> > coord( ntotal, tmp ); // bigger than needed 

        //  copy PDB coord. into common format to use non-overlap subroutines 
        for( i=0; i<np; i++ ) {
            coord[i][OX] = xp[i];
            coord[i][OY] = yp[i];
            coord[i][OZ] = zp[i];
            coord[i][Oocc] = oc[i];
            coord[i][Oznum] = znum[i];
        }

        //  fill the whole vol. with ice excluding the PDB aotm positions
        //    by adding H2O to existing PDB coord.
        //  ice may go into the holes in the PDB so only the density
        //  outside the PDB vol can be calculated.

        vol = ax*by*cz - (zmax-zmin)*(ymax-ymin)*(xmax-xmin);  // vol excluding PDB bounding box
        nMolIce = (int) ( vol  * density );  //  total number of ice mol outside PDB bound. box

        pdbmin[OX] = xmin;
        pdbmin[OY] = ymin;
        pdbmin[OZ] = zmin;
        pdbmax[OX] = xmax;
        pdbmax[OY] = ymax;
        pdbmax[OZ] = zmax;

        //  now do all the work 
        natom = fillIce( coord, nMolIce, np, ax, by, cz, pdbmin, pdbmax, RMIN );

        cout << natom << " total atoms in ice plus PDB structure" << endl;

        // -----   output to file  
        wobble = 0.0;   //  Debye-Waller factor 
        for( i=0; i<natom; i++ ) {
            fpout << setw(5) << coord[i][Oznum] << setw(14) << coord[i][OX] << setw(14)
                << coord[i][OY] << setw(14) << coord[i][OZ] << setw(14) << coord[i][Oocc] 
                << setw(14) << wobble << endl;
        }

        //-----------  write end of file mark -------------------------------- 
        iz = -1;
        fpout << setw(5) << iz << endl;  //  end of data 
        fpout.close( );

    }
    // ------- echo CPU time just for fun ------------ 

    mytime = cputim() - mytime;
    cout << "\ntotal CPU time = " << mytime << " sec." << endl;

    return( EXIT_SUCCESS );

}  // end main 


/*------------------------ fillIce() ------------------------*/
/*
    fill a vol. with amorphous ice as random coord.
    in the range (0,0,0) to (ax,by,cz)

    add ice mol. everywhere avoiding occupied positions
    until density/number outside PDB bounding box is correct
    (i.e. can't calc. ice density around rough edges of PDB structure
    but assume its the same as in vac. if every position is equally 
    added to)

    start with existing PDB coord. and work around it

  coord[][] = start with existing list of PDB atoms sorted by z = coord[][OZ]
                will merge in a set of ICE coord to fill rest of vol.
                dimensions ncoord x NVAL
  nMolIce = number of ice mol. coordinates outside of PDB vol to generate
  np = number of start coord from PDB atoms
  rmin = minimum separation distance

  ax*by*cz = total vol
  (xmin,ymin,zmin) -> (xmax,ymax,zmax) = pdbmin -> pdbmax
        = vol of PDB to exclude for ice density calculation 

  its hard to calculate vol of PDB with random edges so exclude
  a bounding box enclosing the PDB

  assumed globals
     NVAL

*/
int fillIce( vector< vector<double> > &coord, const int nMolIce, const int np,
    double ax, double by, double cz, 
    vector<double> &pdbmin, vector<double> &pdbmax, double rmin )
{
    int ic, i, ncoord;
    vector<double> ctest(NVAL), hcoord(NVAL);

    // initial number of coord. = number in PDB structure
    ncoord = np;    //  total number of atoms PDB+ice
    ic = 1;         // number of ice mol. added
    do {
        // get new coordinate of ice mol (actually only the O atom)
        newcoord( ctest, ax, by, cz );
        ctest[Oocc]  = 1.0; //  occupancy = 1.0
        ctest[Oznum] = 8;   // Z=8 for oxygen
        //cout << "new coord= " << ctest[0] << ", " << 
        //           ctest[1] << ", " <<ctest[2] << endl;   // ??? diagnostic

        //  add it to the list if its OK - this also sorts
        if( testcoord( ctest, coord, ncoord, &i, RMIN ) == TRUE ) {  
            insertcoord( ctest, coord, ncoord, i);
            ncoord++;

            // count mol. outside PDB bounding box
            if( 0 == testInside( ctest, pdbmin, pdbmax) ) { 
                ic++;
                if( ic%1000 == 0 ) cout << "\r nmol= " << ic << flush;
            }

            //  add 2 x H atoms at random positions near O atom
            //  should constrain bond angle but ignore for simplicity
            //   insert is slow - might be faster to insert 3 at one time (more code needed)
            hcoord[Oocc]  = 1.0;    //  occupancy = 1.0
            hcoord[Oznum] = 1;      //  Z=1 for hydrogen
            hydrocoord( hcoord, ctest);
            insertcoord( hcoord, coord, ncoord, i);
            ncoord++;
            hydrocoord( hcoord, ctest);
            insertcoord( hcoord, coord, ncoord, i);
            ncoord++;
        }
    }  while( ic < nMolIce );

    cout << endl;
    return ncoord;

}  // end fillIce()


/*------------------------ fillSolid() ------------------------*/
/*
    fill a vol. with an amorphous solid with random coord.
    in the range (0,0,0) to (ax,by,cz)

  coord[][] = to get list of existing sorted by z = coord[][OZ]
                dimensions ncoord x NVAL
  ncoord = number of coordinates to generate
  rmin = minimum separation distance

  assumed globals
     NVAL

*/
int fillSolid( vector< vector<double> > &coord, const int ncoord, 
    double ax, double by, double cz, double rmin )
{
    int ic, i;
    vector<double> ctest(NVAL);

    newcoord( coord[0], ax, by, cz );   // start coord. list
    ic = 1;
    do {
        // get new coordinate 
        newcoord( ctest, ax, by, cz );
        //cout << "new coord= ", << ctest[0] << ", " << 
        //           ctest[1] << ", " <<ctest[2] << endl;

        //  add it to the list if its OK - this also sorts
        if( testcoord( ctest, coord, ic, &i, RMIN ) == TRUE ) {
            insertcoord( ctest, coord, ic, i);
            ic++;
            if( ic%1000 == 0 ) cout << "\r natom= " << ic << flush;
        }
    }  while( ic < ncoord );

    return ic;

}  // end fillSolid()

/*------------------------ hydrocoord() ------------------------*/
/*
    generate coordinates for Hydrogen around Oxygen atom

  NOTE: this ignores possible overlap of H atoms 
        and the bonding angle of H atoms

  H bond length is < min. sep. distance so overlap unlikely

  ctest[] = to get H coord.
  coord[] = position of O atom

  assumed global:
     iseed = RNG seed
*/
void hydrocoord( vector<double> &ctest, vector<double> &coord )
{
    long i;
    double dx, d=0.0;
    
    /* Hydrogen-Oxygen bond length in Angstroms */
    static const double bondlength=0.958;

    //  get random x,y,z offset at spec. bond length
    for( i=0; i<3; i++) {
        ctest[i] = dx = ranflat( &iseed ) - 0.5 ;
        d += dx*dx;
    }
    d = bondlength / sqrt( d );
    
    for( i=0; i<3; i++) {
        ctest[i] = coord[i] + ctest[i]*d;
    }

    //  probably should generate a 2nd H coord at fixed bond angle here (???)

    return;
} /* end hydrocoord() */

/*------------------------ insertcoord() ------------------------*/
/*
    insert a new coordinate in the list

  ctest[] = new coord to insert
  coord[][] = list of existing coordinates sorted by z = coord[][OZ]
  ncoord = number of coordinates
  pos = position to insert coord. at

  insert new coord at index pos and move all the rest down one

  assumed globals
     NVAL

     might use vector<> insert here sometime - needs an iterator though
*/
void insertcoord( vector<double> &ctest, vector<vector<double> > &coord, 
        const int ncoord, const int pos )
{
    int i, j, n;

    //  check that there is enough memory left
    n = (int) coord.size();
    if( ncoord >= n ) {
        cout << "ncoord = " << ncoord << " too big in insertcoord(), exit...." << endl;
        exit( 0 );
    }
    
    for( i=ncoord; i>pos; i--) {
        for( j=0; j<NVAL; j++)
            coord[i][j] = coord[i-1][j];
    }

    for( j=0; j<NVAL; j++)
            coord[pos][j] = ctest[j];

} // end insertcoord() 


/*------------------------ newcoord() ------------------------*/
/*
    generate a new random coordinate inside the required volume

    xmax, ymax, zmax = volume size

    assumed global
        iseed
*/
void newcoord( vector<double> &c, 
    const double xmax, const double ymax, const double zmax )
{
    c[OX] = xmax * ranflat( &iseed );
    c[OY] = ymax * ranflat( &iseed );
    c[OZ] = zmax * ranflat( &iseed );
    return;

} // end newcoord() 

/*----------------- sortByZd() ------------------------------

from slicelib.cpp but converted from float to double

    improved Shell sort modeled after prog. 6.5 (pg. 274) of
    R. Sedgewick, "Algorithms in C", 3rd edit. Addison-Wesley 1998
    
    x[], y[], z[]   = atom coordinates 
    occ[]           = occupancy of each atom
    Znum[]          = atomic number of each atom
    natom           = number of atoms
*/
void sortByZd( vectord &x, vectord &y, vectord &z, vectord &occ,
    vectori &Znum, int natom )
{
    int i, j, h, Znum2;
    //float x2, y2, z2, occ2;
    double x2, y2, z2, occ2;

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


/*------------------------ testcoord() ------------------------*/
/*
    test the new coordinate to see if its too close to an existing
    coordinate (closer than RMIN)

    A straight search of the whole list is very slow (proportional to N^2)
    and was not pratical for a large set of coordinates (i.e. took
    an absurd amount of CPU time).

    This routines assumes that the list of existing coordinates is
    sorted wrt one coord (z in this case).  Once the position of the new
    coordinate is located in the list (using a binary search), then
    only the small range of coord. near this point need to be
    tested.  This makes the test dramatically faster.

  ctest[] = new coord to test
  coord[][] = list of existing coordinates sorted by z = coord[][OZ]
  ncoord = number of coordinates
  pos = returned position to insert coord.

  assumed globals
     OX, OY, OZ, RMIN

*/
int testcoord( vector<double> &ctest, vector< vector<double> > &coord, 
    const int ncoord, int *pos, const double rmin )
{
    /* this is the more -sophisticated version sorted by Z */
    long i, j, k;
    double d, dx, dy, dz, dz2, rmin2, z, range;

    /* ---- find postion of ctest[] in coord[][] using a binary search  --- 
         i should get the position to insert (between i and i+1)
         this assumes coord[][] is sorted by coord OZ
    */
//  printf("testcoord() top, ncoord= %d\n", ncoord );

    z = ctest[OZ];
    if( z <= coord[0][OZ] ) i = j = 0;
    else if( z >= coord[ncoord-1][OZ] ) i = j = ncoord;
    else { 
        i = 0;
        j = ncoord-1;
        do{ k = ( i + j ) / 2 ;
            if( z < coord[k][OZ] )  j = k;
            else if( z >=  coord[k][OZ] ) i = k;
        } while ( (j-i) > 1 );
    }

    /* now that we have the position of the new point
       we only have to explore within RMIN of this point */

    rmin2 = rmin*rmin;
    range = 2.0*rmin2;   /* add a little safety margin */
    k = j;
    while( k >= 0 ) {
            dx = ctest[OX] - coord[k][OX];
            dy = ctest[OY] - coord[k][OY];
            dz = ctest[OZ] - coord[k][OZ];
            dz2 = dz*dz;
            d = dx*dx + dy*dy + dz2;
            if( d <= rmin2 ) return( FALSE );
            if( dz2 > range ) break;
            k--;
    }
    k = j;
    while( k < ncoord ) {
            dx = ctest[OX] - coord[k][OX];
            dy = ctest[OY] - coord[k][OY];
            dz = ctest[OZ] - coord[k][OZ];
            dz2 = dz*dz;
            d = dx*dx + dy*dy + dz2;
            if( d <= rmin2 ) return( FALSE );
            if( dz2 > range ) break;
            k++;
    }

    /* if it gets to here then this is a good point */

    *pos = j;

    return( TRUE );

} // end testcoord() 

/*------------------------ testInside() ------------------------*/
/*
  test if position ctest=(x,y,z) is inside
 (xmin,ymin,zmin) -> (xmax,ymax,zmax) = pdbmin -> pdbmax
        = vol of PDB to exclude for ice density calculation 

*/
int testInside( vector<double> &ctest,
    vector<double> &pdbmin, vector<double> &pdbmax )
{
    if( (    ctest[OX] > pdbmin[OX] ) && ( ctest[OX] < pdbmax[OX] ) 
        && ( ctest[OY] > pdbmin[OY] ) && ( ctest[OY] < pdbmax[OY] )  
        && ( ctest[OZ] > pdbmin[OZ] ) && ( ctest[OZ] < pdbmax[OZ] )  ) return (1);
    else return( 0 );
}
