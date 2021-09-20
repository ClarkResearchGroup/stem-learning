/*
      *** someAtoms.hpp ***

------------------------------------------------------------------------
Copyright 2016 Earl J. Kirkland

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

function to initialize atomic coordiantes etc. for a few sample specimens
(as would be read in from an xyz file)

inputs:

	nselect = int to select specific specimen
	ncellx,ncelly, ncellz = int number of unit cells in supercell
	ax,by,cz = float/double unit cell size (in Ang.)
	Znum[] = int atomic number
	x[],y[],z[] = float/double atomic coordinats (in Ang.)
	occ[] = float/double occupancies
	wobble[] = float/double thermal vibration amplitude (in Ang.)
	xyzTitle = string to get spec. description


	all but nselect are overwritten

	formatted for a TAB size of 4 char
------------------------------------------------------------------------

   started 15-may-2016 ejk
   working 28-may-2016 ejk
   last modified 26-may-2016 ejk

*/

// parameterize the data type (T= float or double) just in case.....
template< class T >
int someAtoms(  int nselect, int &ncellx, int &ncelly, int &ncellz,
	T &ax, T &by, T &cz, vector<int> &Znum, vector<T> &x, vector<T> &y, vector<T> &z, 
	vector<T> &occ, vector<T> &wobble, string &xyzTitle )
{
	int natoms;

	//  do this here so only one .resize() needed instead of every specimen
	if( 0 == nselect ) natoms = 4;	        // Si110
	else if ( 1 == nselect ) natoms = 5;	// SrTiO
	else if ( 2 == nselect ) natoms = 4;	// graphene
	else if ( 3 == nselect ) natoms = 13;	// ybco

	Znum.resize( natoms );
	x.resize( natoms );
	y.resize( natoms );
	z.resize( natoms );
	occ.resize( natoms );
	wobble.resize( natoms );

	//  really need to find a more compact way to enter xyz data files
	//          for large unit cells
    if( 0 == nselect ) {	// Si 110

            xyzTitle = "Si 110 built-in";
			//natoms = 4;

            ax = (T) 3.83959;  //  unit cell size
            by = (T) 5.4300;
            cz = (T) 3.83959;

            ncellx = 3;
            ncelly = 2;
            ncellz = 3;

            Znum[0] = 14;	//  atom 1
            x[0] = (T) 0.0;
            y[0] = (T) 0.0;
            z[0] = (T) 0.0;
            occ[0] = (T) 1.0;
            wobble[0] = (T) 0.08;

            Znum[1] = 14;	//  atom 2
            x[1] = (T) 0.0;
            y[1] = (T) 1.3575;
            z[1] = (T) 1.91979;
            occ[1] = (T) 1.0;
            wobble[1] = (T) 0.08;

            Znum[2] = 14;	//  atom 3
            x[2] = (T) 1.91979;
            y[2] = (T) 2.7150;
            z[2] = (T) 1.91979;
            occ[2] = (T) 1.0;
            wobble[2] = (T) 0.08;

            Znum[3] = 14;	//  atom 4
            x[3] = (T) 1.91979;
            y[3] = (T) 4.0725;
            z[3] = (T) 3.83959;
            occ[3] = (T) 1.0;
            wobble[3] = (T) 0.08;

    }else if ( 1 == nselect ) {	// SrTiO - wobble prob. not right - used Si

            xyzTitle = "SrTiO built-in";
			//natoms = 5;

            ax = (T) 3.9051;  //  unit cell sizel
            by = (T) 3.9051;
            cz = (T) 3.9051;

            ncellx = 3;
            ncelly = 3;
            ncellz = 3;

            Znum[0] = 38;	//  atom 1
            x[0] = (T) 0.0;
            y[0] = (T) 0.0;
            z[0] = (T) 0.0;
            occ[0] = (T) 1.0;
            wobble[0] = (T) 0.08;

            Znum[1] = 8;	//  atom 2
            x[1] = (T) 1.95255;
            y[1] = (T) 1.95255;
            z[1] = (T) 0.0;
            occ[1] = (T) 1.0;
            wobble[1] = (T) 0.08;

            Znum[2] = 22;	//  atom 3
            x[2] = (T) 1.95255;
            y[2] = (T) 1.95255;
            z[2] = (T) 1.95255;
            occ[2] = (T) 1.0;
            wobble[2] = (T) 0.08;

            Znum[3] = 8;	//  atom 4
            x[3] = (T) 1.95255;
            y[3] = (T) 0.0;
            z[3] = (T) 1.95255;
            occ[3] = (T) 1.0;
            wobble[3] = (T) 0.08;

            Znum[4] = 8;	//  atom 5
            x[4] = (T) 0.0;
            y[4] = (T) 1.95255;
            z[4] = (T) 1.95255;
            occ[4] = (T) 1.0;
            wobble[4] = (T) 0.08;

    }else if ( 2 == nselect ) {	// graphene - wobble prob. not right - used Si

            xyzTitle = "graphene built-in";
			//natoms = 4;

            ax = (T) 2.456;  //  unit cell sizel
            by = (T) 4.2539168;
            cz = (T) 3.3480;

            ncellx = 3;
            ncelly = 2;
            ncellz = 1;

            Znum[0] = 6;	//  atom 1
            x[0] = (T) 0.0;
            y[0] = (T) 0.0;
            z[0] = (T) 1.6740;
            occ[0] = (T) 1.0;
            wobble[0] = (T) 0.08;

            Znum[1] = 6;	//  atom 2
            x[1] = (T) 1.2280;
            y[1] = (T) 0.70898613;
            z[1] = (T) 1.6740;
            occ[1] = (T) 1.0;
            wobble[1] = (T) 0.08;

            Znum[2] = 6;	//  atom 3
            x[2] = (T) 1.2280;
            y[2] = (T) 2.12695839;
            z[2] = (T) 1.6740;
            occ[2] = (T) 1.0;
            wobble[2] = (T) 0.08;

            Znum[3] = 6;	//  atom 4
            x[3] = (T) 2.4560;
            y[3] = (T) 2.83594452;
            z[3] = (T) 1.6740;
            occ[3] = (T) 1.0;
            wobble[3] = (T) 0.08;

    }else if ( 3 == nselect ) {	// YBCO

            xyzTitle = "YBa2Cu3O7 down the b axis, built-in, Beno et al APL 51(1987)57";
			//natoms = 13;

            ax = (T) 3.8231;  //  unit cell sizel
            by = (T) 11.6807;
            cz = (T) 3.8864;

            ncellx = 3;
            ncelly = 1;
            ncellz = 3;

			//39    1.91155    5.84035   1.9432       1.0     0.108  
            Znum[0] = 39;	//  atom 1
            x[0] = (T) 1.91155;
            y[0] = (T) 5.84035;
            z[0] = (T) 1.9432;
            occ[0] = (T) 1.0;
            wobble[0] = (T) 0.108;

			//56    1.91155    2.15275   1.9432       1.0     0.117
            Znum[1] = 56;	//  atom 2
            x[1] = (T) 1.91155;
            y[1] = (T) 2.15275;
            z[1] = (T) 1.9432;
            occ[1] = (T) 1.0;
            wobble[1] = (T) 0.117;

			//29    0.0        0.0       0.0          1.0     0.113
            Znum[2] = 29;	//  atom 3
            x[2] = (T) 0.0;
            y[2] = (T) 0.0;
            z[2] = (T) 0.0;
            occ[2] = (T) 1.0;
            wobble[2] = (T) 0.113;

			//29    0.0        4.15365   0.0          1.0     0.086
            Znum[3] = 29;	//  atom 4
            x[3] = (T) 0.0;
            y[3] = (T) 4.15365;
            z[3] = (T) 0.0;
            occ[3] = (T) 1.0;
            wobble[3] =  (T) 0.086;

			//8     0.0        0.0       1.9432       1.0     0.185
            Znum[4] = 8;	//  atom 5
            x[4] = (T) 0.0;
            y[4] = (T) 0.0;
            z[4] = (T) 1.9432;
            occ[4] = (T) 1.0;
            wobble[4] =  (T) 0.185;

			//8     1.91155    4.40713   0.0          1.0     0.119
            Znum[5] = 8;	//  atom 6
            x[5] = (T) 1.91155;
            y[5] = (T) 4.40713;
            z[5] = (T) 0.0;
            occ[5] = (T) 1.0;
            wobble[5] =  (T) 0.119;

			//8     0.0        4.42582   1.9432       1.0     0.097
            Znum[6] = 8;	//  atom 7
            x[6] = (T) 0.0;
            y[6] = (T) 4.42582;
            z[6] = (T) 1.9432;
            occ[6] = (T) 1.0;
            wobble[6] =  (T) 0.097;

			//8     0.0        1.8502    0.0          1.0     0.130 
            Znum[7] = 8;	//  atom 8
            x[7] = (T) 0.0;
            y[7] = (T) 1.8502;
            z[7] = (T) 0.0;
            occ[7] = (T) 1.0;
            wobble[7] =  (T) 0.130;

			//56    1.91155    9.52795   1.9432       1.0     0.117
            Znum[8] = 56;	//  atom 9
            x[8] = (T) 1.91155;
            y[8] = (T) 9.52795;
            z[8] = (T) 1.9432;
            occ[8] = (T) 1.0;
            wobble[8] =  (T) 0.117;

			//29    0.0        7.52705   0.0          1.0     0.086
            Znum[9] = 29;	//  atom 10
            x[9] = (T) 0.0;
            y[9] = (T) 7.52705;
            z[9] = (T) 0.0;
            occ[9] = (T) 1.0;
            wobble[9] =  (T) 0.086;

			//8     1.91155    7.27357   0.0          1.0     0.119
            Znum[10] = 8;	//  atom 11
            x[10] = (T) 1.91155;
            y[10] = (T) 7.27357;
            z[10] = (T) 0.0;
            occ[10] = (T) 1.0;
            wobble[10] =  (T) 0.119;

			//8     0.0        7.25488   1.9432       1.0     0.097
            Znum[11] = 8;	//  atom 12
            x[11] = (T) 0.0;
            y[11] = (T) 7.25488;
            z[11] = (T) 1.9432;
            occ[11] = (T) 1.0;
            wobble[11] =  (T) 0.097;

			//8     0.0        9.8287    0.0          1.0     0.130
            Znum[12] = 8;	//  atom 13
            x[12] = (T) 0.0;
            y[12] = (T) 9.8287;
            z[12] = (T) 0.0;
            occ[12] = (T) 1.0;
            wobble[12] =  (T) 0.130;
		}

		return( natoms );

}  //  end someAtoms()