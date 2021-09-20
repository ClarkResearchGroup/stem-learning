/*      
    *** cudaSlice.cu ***   

------------------------------------------------------------------------
Copyright 2018-2019 Earl J. Kirkland


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

-----------------------------------------------------------------------------

  cuda subroutines:

	cmplPixMul()     : complex pix mul with shift
	cmplVecMul()     : complex vector mul
	cuAtompot()      : calculate atomic potential of one slice
	cuBWlimit()      : bandwidth limit
	cuFreq()         : calculate FFT frequencies
	cuPhasegrating() : phase grating
	integCBED()      : integrate ADF detector
	magSqPix()       : form sq. magnitude of pix
	probeShift()     : shift probe in FFT space
	zeroDbleArray()  : set double array to zero

  started from autostem_cuda.cu 18-aug-2018 ejk
  add probeShift(), zeroDbleArray(), integCBED() 18-aug-2018 ejk
  add cuAtompot(), cuBWlimit(), cuPhasegrating()  27-sep-2018 ejk 
  add full cache (for all NZMAX) to cuAtompot() 15-sep-2019 ejk

    this file is formatted for a TAB size of 4 characters 
*/


#ifndef CUDASUBS_INCLUDED
#define CUDASUBS_INCLUDED

//==================  extra CUDA stuff ===================


//---------------  cmplPixMul() --------------
//
// CUDA kernel definition for 2D pix mul with trans larger than probe
//     probe is mul. by a subset of trans
//  perform operation probe = probe * trans with offset
//  probe = nxprobe x nyprobe (no bigger than trans)
//  trans = nx x ny  (may be bigger than probe)
//  ixoff, iyoff = offset of probe inside trans ; edges will wrap around
//
__global__ void cmplPixMul( cufftComplex *trans, cufftComplex *probe, int nx, int ny,
     int nxprobe, int nyprobe, int ixoff, int iyoff) 
{
    // 2D index into probe array
    int ix = blockDim.x*blockIdx.x + threadIdx.x;  //  range 0 to (nxprobe-1)
    int iy = blockDim.y*blockIdx.y + threadIdx.y;  //  range 0 to (nyprobe-1)

    //  use scratch var so probe[] can be overwriten
    if( (ix < nxprobe) && (iy < nyprobe) ) {

        //  calculate 2D index into trans array
        int ixt = ix + ixoff;
        if( ixt >= nx ) ixt = ixt - nx;
        else if( ixt < 0 ) ixt = ixt + nx;

        int iyt = iy + iyoff;
        if( iyt >= ny ) iyt = iyt - ny;
        else if( iyt < 0 ) iyt = iyt + ny;

        int it = iyt + ixt*ny;
        int ip = iy + ix*nyprobe;

        float ar, ai, br, bi;
 
        ar = trans[it].x;  // real
        ai = trans[it].y;  // imag
        br = probe[ip].x;
        bi = probe[ip].y;
        probe[ip].x = ar*br - ai*bi;  // real
        probe[ip].y = ar*bi + ai*br;  // imag
    }
}   //  end complPixMul()

//---------------  cmplVecMul() --------------
//
// cuda kernel definition for complex vector mul
//   c = a * b (element by element)
//
__global__ void cmplVecMul( cufftComplex *a, cufftComplex *b, cufftComplex *c, int nmax) 
{
    int i = blockDim.x*blockIdx.x + threadIdx.x;

    //  use scratch var so c[] can overwrite one of inputs
    if( i < nmax ) {
        float ar, ai, br, bi;
        ar = a[i].x;  // real
        ai = a[i].y;  // imag
        br = b[i].x;
        bi = b[i].y;
        c[i].x = ar*br - ai*bi;
        c[i].y = ar*bi + ai*br;
    }
}  //  end cmplPixMul()

/*---------------  cuAtompot() --------------

  CUDA kernel definition to calculate single layer projected atomic potential

  this is actually no faster than doing the potential in real space on the host
	but save anyway just in case!

  calculate the summation over atoms at one point (kx,ky) in reciprocal space
  
  It is better to sum in reciprocal space to use fine-grain parallelism
  as on a GPU. Every point can then run in parallel without trying to
  access the same point. This is very different than the openMP version
  of trlayer() in autostem.cpp which summs the atomic potential in real space.

  potn[] = nx x ny  output array = half of complex plane for C2R FFT

  x[],y[] = real array of atomic coordinates
  occ[]   = real array of occupancies
  Znum[]  = array of atomic numbers

  spec[ k + 4*iatom] = packed array of x,y,occ,Znum (min. GPU transfers)
						(k=0,1,2,3 for x,y,occ,Znum)
  istart  = starting index of atom coord.
  natom   = number of atoms
  ax, by  = size of transmission function in Angstroms
  kev     = beam energy in keV
  trans   = 2D array to get complex specimen
        transmission function
  nx, ny  = dimensions of transmission functions
  *phirms = average phase shift of projected atomic potential
  *nbeams = will get number of Fourier coefficients
  k2max   = square of max k = bandwidth limit
  fparams[] = scattering factor parameters

  scale = mm0 * wavelength (put here for comparison to original trlayer()

	repeat scaling from mulslice.cpp
		mm0 = 1.0F + v0/511.0F;
		wavlen = (float) wavelength( v0 );
		scale = wavlen * mm0;

*/

__global__ void cuAtompot( cufftComplex *potn, 
	 float spec[], int natom,  int istart,
    const float ax, const float by, const float kev,
    const int nx, const int ny,
	float kx[], float ky[], float kx2[], float ky2[],
    const float k2max, double fparams[], const float scale ) 
{

	// for the atomic scattering factor tables
	const int NPMAX=   12;  // number of parameters for each Z
	const int NZMIN=   1;   // min Z 
	const int NZMAX=   103; // max Z 
	const int nl=3, ng=3;   // number of Lorenzians and Gaussians

	int i, j, iatom, Z, nymid;
	float k2, x, y, occ;
	double sumf, sumr, sumi, w;
	double fe[NZMAX+1];  // scatt. factor
	const double twopi = 6.283185307;

    // 2D index into trans array
    int ix = blockDim.x*blockIdx.x + threadIdx.x;  //  range 0 to (nxprobe-1)
    int iy = blockDim.y*blockIdx.y + threadIdx.y;  //  range 0 to (nyprobe-1)

	nymid = 1 + ny/2;

    //  valid trans[] index
    if( (ix < nx) && (iy < nymid) ) {

		// init to something in case it fails later
		potn[ iy + ix*nymid].x = 0.0F;  // real part
		potn[ iy + ix*nymid].y = 0.0F;  // imag part

		k2 = kx2[ix]+ky2[iy];

		if( k2 < k2max ) {
            
            // init scatt. factor to <0 to indicate not yet valid
            //   - should always be > 0.0 for real values
            for( j=0; j<(NZMAX+1); j++) fe[j] = -100.0;

			sumr = sumi = 0.0;
			for( iatom=istart; iatom<(istart+natom); iatom++) {

				Z = (int) ( spec[3 + 4*iatom] + 0.2); // round up to avoid truncation
				if( (Z<NZMIN) || (Z>NZMAX) ) return;

				// save old values in a look up table for repeated Z 
				//  - don't repeat sum - speeds thing up a lot  
				if( fe[Z] <  0.0 ) {
					sumf = 0.0;

					// Lorenzians - from slicelib.cpp  
					for( i=0; i<2*nl; i+=2 )
						sumf += fparams[i + Z*NPMAX]/( k2 + fparams[i+1 + Z*NPMAX] );
							 // mimic fparams[Z][i]/( k2 + fparams[Z][i+1] );

					// Gaussians - from slicelib.cpp 
					for( i=2*nl; i<2*(nl+ng); i+=2 )
						sumf += fparams[i + Z*NPMAX]*exp( - k2 * fparams[i+1 + Z*NPMAX] );
							// mimic fparams[Z][i]*exp( - k2 * fparams[Z][i+1] );

					fe[Z] = sumf;  // save scattering factor for this k
				}
				x = spec[ 0 + 4*iatom];
				y = spec[ 1 + 4*iatom];
				occ = spec[ 2 + 4*iatom];
				w = twopi * ( kx[ix]*x + ky[iy]*y );

				sumr += fe[Z] * cos( -w ) * occ;
				sumi += fe[Z] * sin( -w ) * occ;

			}  //end for( iatom...)

			potn[ iy + ix*nymid].x = scale * sumr;  // real part
			potn[ iy + ix*nymid].y = scale * sumi;  // imag part

		}  // end if( k2 < ....

	}  //  end if((ix<nx)....

	return;

}   //  end cuAtompot()

/*---------------  cuBWlimit() --------------

  bandwidth limit tran[] - assumed to be in reciprocal space
  and add FFT scale
  
  kx2[],ky2[] = spatial freq. sq.
  k2max = max spatial freq.

*/

__global__ void cuBWlimit( cufftComplex *trans, 
	float *kx2, float *ky2, float k2max, const int nx, const int ny ) 
{

    // 2D index into trans array
    int ix = blockDim.x*blockIdx.x + threadIdx.x;  //  range 0 to (nxprobe-1)
    int iy = blockDim.y*blockIdx.y + threadIdx.y;  //  range 0 to (nyprobe-1)

    if( (ix < nx) && (iy < ny) ) {
		float scale = 1.0F/( (float) (nx*ny) );

        float k2 = kx2[ix] + ky2[iy];

        if( k2 > k2max ) {
        	trans[iy + ix*ny].x = 0.0F;
        	trans[iy + ix*ny].y = 0.0F;
        } else {
        	trans[iy + ix*ny].x *= scale;  //  fix FFT scale
        	trans[iy + ix*ny].y *= scale;
		}
        
    }  // end if( ix< nx...
    
}  //  end cuBWlimit()

/*---------------  cuFreq() --------------
//
// cuda kernel definition to calculate spatial freq.
//   
    ko[n]  = real array to get spatial frequencies
    ko2[n] = real array to get k[i]*k[i] 
    nk     = integer number of pixels
    ak     = real full scale size of image in pixels
*/
__global__ void cuFreq( float ko[], float ko2[], int nk, float ak ) 
{
    int i = blockDim.x*blockIdx.x + threadIdx.x;

    //  check for valid thread index
    if( i < nk) {

		int imid = (int) ( nk/2.0 + 0.5);   /* when nk may not be 2^m */

        if ( i > imid ) {
            ko[i]  = ((float)(i-nk)) / ((float)ak);
        } else {
            ko[i]  = ((float)i) / ((float)ak);
        }
        ko2[i] = ko[i] * ko[i];
    }
}  //  end cuFreq()

/*---------------  cuPhasegrating() --------------

  Start with the atomic potential from cuAtompot() after inv. FFT
  and convert to the transmission function as in a phase grating calculation
  - assume its scaled to a phase already
*/

__global__ void cuPhasegrating( float * potnR, cufftComplex *trans, 
    const int nx, const int ny ) 
{

    // 2D index into trans array
    int ix = blockDim.x*blockIdx.x + threadIdx.x;  //  range 0 to (nxprobe-1)
    int iy = blockDim.y*blockIdx.y + threadIdx.y;  //  range 0 to (nyprobe-1)

    if( (ix < nx) && (iy < ny) ) {
         
		//  atomic potential
		double vz= potnR[iy + ix*ny];
		trans[iy + ix*ny].x = (float) cos(vz);  // real
        trans[iy + ix*ny].y = (float) sin(vz);  // imag

    }  // end if( (ix<nx)....
    
}   //  end cuPhasegrating()


//---------------  integCBED() --------------
//
// CUDA kernel definition to integrate STEM detector active regions
//
//  remember:
//      [1] many threads cannot access the same sumation variable at
//           one time so sum along only one direction at a time (into a 1D array)
//				- complete the last sum in 1D on the host
//      [2] many points will not be on the active portion of the detector
//             so there is less competition among threads than it might seem
//
//  cbed = input nx x ny float CBED pix = |cpix|^2
//  sums = oout double[nx]  to get sum |cpix|^2 along iy
//  nx, ny = size of cbed
//  collectorMode = detector type
//  kxp[],kyp[] = spatial freq.
//  kxp2[],kyp2[] = spatial freq. sq.
//  k2min, k2max = detector range in polar direction
//  phimin, phimax = detector range in azimuthal direction
//
__global__ void integCBED( double *sums, float *cbed, int nx, int ny,
    int collectorMode, float *kxp, float *kyp, float *kxp2, float *kyp2,
    float k2min, float k2max, float phiMin, float phiMax ) 
{
    // modes of collector
    enum{ ADF=0, CONFOCAL=1, ADF_SEG=2, TOTAL=3};  //  no confocal here yet
    
    // 2D index into fpix array
    int ix = blockDim.x*blockIdx.x + threadIdx.x;  //  range 0 to (nx-1)
    int iy = blockDim.y*blockIdx.y + threadIdx.y;  //  range 0 to (ny-1)

    if( (ix < nx) && (iy < ny) ) {

        //  calculate 2D index into cbed array
        int it = iy + ix*ny;
        float k2, phi;
        
        k2 = kxp2[ix] + kyp2[iy];

        //  use atomicSum() so only one thread can access sums[] at a time
		//    - need nvcc compiler option "-arch=sm_61" for atomicAdd()
        if( ADF == collectorMode ) {
            if( (k2 >= k2min ) && (k2 < k2max ) ) atomicAdd( &(sums[iy]), (double)cbed[it] );
        } else if( ADF_SEG == collectorMode ) {
            phi = atan2( kyp[iy], kxp[ix] ); 
            if( (k2 >= k2min ) && (k2 < k2max)  
				&& (phi >= phiMin) && (phi< phiMax) )
						atomicAdd( &(sums[iy]), (double)cbed[it] );
        } else if( TOTAL == collectorMode ) {
				atomicAdd( &(sums[iy]), (double)cbed[it] );
        }

    }  // end if( (ix < nx)....
    
}   //  end integCBED()


//---------------  magSqPix() --------------
//
// CUDA kernel definition for 2D pix complex to magnitude
//     take square magnitude on GPU to
// 
//  cpix = nx x ny complex
//  fpix = nx x ny  float = |cpix|^2
//  nx, ny = size of both pix
//
__global__ void magSqPix( float *fpix, cufftComplex *cpix, int nx, int ny ) 
{
    // 2D index into both arrays
    int ix = blockDim.x*blockIdx.x + threadIdx.x;  //  range 0 to (nx-1)
    int iy = blockDim.y*blockIdx.y + threadIdx.y;  //  range 0 to (ny-1)

    if( (ix < nx) && (iy < ny) ) {

        int it = iy + ix*ny;  //  calculate 2D index into trans array
        float ar, ai;
 
        ar = cpix[it].x;  // real
        ai = cpix[it].y;  // imag
        fpix[it] = ar*ar + ai*ai;
    }
}   //  end magSqPix()

//---------------  probeShift() --------------
//
//  CUDA kernel definition for 2D probe shift in FT space
//  perform operation probe *= exp( 2*pi*i * x * k) with offset
//
//  prb0 = input nx x ny complex
//  prbs = output gets prb0 shifted by (xs,ys)
//  xs,ys = amount to shift
//  nx, ny = size of both pix
//  kx[], ky[] = arrays of spatial frequencies
//
__global__ void probeShift( cufftComplex *prbs, cufftComplex *prb0, int nx, int ny,
    float xs, float ys, float *kx, float *ky ) 
{
    // 2D index into both arrays
    int ix = blockDim.x*blockIdx.x + threadIdx.x;  //  range 0 to (nx-1)
    int iy = blockDim.y*blockIdx.y + threadIdx.y;  //  range 0 to (ny-1)

    if( (ix < nx) && (iy < ny) ) {

        int it = iy + ix*ny;  //  calculate 2D index into arrays
        
        double ar, ai, br, bi;
        double w = 6.283185307*(xs*kx[ix] + ys*ky[iy]);
 
        ar = prb0[it].x;  // real
        ai = prb0[it].y;  // imag
        br = cos( w );
        bi = sin( w );
        
        prbs[it].x = (float) (ar*br - ai*bi);  // real
        prbs[it].y = (float) (ar*bi + ai*br);  // imag
    }
}   //  end probeShift()


//---------------  zeroDbleArray() --------------
//
// CUDA kernel definition to zero a double array
// 
//  a[nmax] = double array
//  nmax = size of array
//
__global__ void zeroDbleArray( double *a, int nmax ) 
{
    // index into array
    int i = blockDim.x*blockIdx.x + threadIdx.x;  //  range 0 to (nmax-1)

    if( i < nmax ) a[i] = 0.0;

}   //  end zeroDbleArray()

#endif
