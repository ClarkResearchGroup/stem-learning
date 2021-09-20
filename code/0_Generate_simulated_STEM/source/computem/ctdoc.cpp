/*   ctdoc.cpp

------------------------------------------------------------------------
Copyright 2013-2019 Earl J. Kirkland


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

computem document class

started 31-dec-2012 ejk
change a few things around to work with wx 2.9.4 on 18-apr-2013 ejk
move aber and samp dialog into here 4-may-2013 ejk
add expandCell() 12-may-2013 ejk
start adding autostem 1-sep-2013, working 23-sep-2013ejk
fix output of autostem pix (different nx,ny size) 24-sep-2013 ejk
add code to get ctiltx,y from general dialog 26-sep-2013 ejk
add code to get cpertmin,max 29-sep-2013 ejk
add param[ pMODE ] values to mark which program made the data 
    19-oct-2013 ejk
update parameter checking in autostem mode 26-oct-2013 ejk
add OnNewDocument() with built in specimens (Si110, SrTiO, graphene)
    9-nov-2013 ejk
add save pix as text 16-nov-2013 ejk
add model3d mode 28-nov-2013 ejk
fix scaling in model3d/sphere with nx != ny 7-dec-2013 ejk
add load/save param menu+subroutine 8-feb-2014 ejk
add separate sparam[] with copy of last used param[] for output
    in case param[] gets editted after calcluation 11-feb-2014 ejk
update ReadXYZcoord() argument list for new version of slicelib.cpp
    15-may-2014 ejk
add OnEditCoord() 18 to 27-may-2014  ejk
add DoSaveCoord() 29-may-2014 ejk
add 2D abb. phase error calc. 31-aug-2014 ejk
add OnViewxx() to scale image 11 to 18-sep-2014 ejk
add option to add poisson noise 8-oct-2015 ejk
add YBCO built-in specimen option for new  14-may-2016 ejk
switch to someAtoms.hpp  22-may-2016 ejk
convert malloc1D() etc to vector<> 21-jul-2016 ejk
convert to autoslic.cpp with multithread over phonon TDS config
    28-aug-2016 ejk
add view/fft/real space  9-apr-2017 ejk
fix bug in exit wave calc. when called first 16-apr-2017 ejk
change default STEM image size to 64x64 6-may-2017 ejk
add diffraction 9,15-jul-2017 ejk
remove malloc2D(), malloc3D(), free2D(), free3D()
   and convert to new2D(), new3D(), delete2D(), delete3D() 25-sep-2017 ejk
last modified 25-sep-2017 ejk

*/
#include <sstream>   // STD string streams
#include <fstream>   // STD file IO streams
#include <iomanip>   // STD to format the output

#include "computem.hpp"  //  header for this file
#include "ctdoc.hpp"     //  my document class
#include "ctview.hpp"    //  my view class

#include "aberDialog.hpp"     // custom dialog for aberrations
#include "generalDialog.hpp"  // custom dialog for general parameters
#include "rangeDialog.hpp"    // custom dialog for range of greyscale
#include "sampleDialog.hpp"   // custom dialog for sampling
#include "textDialog.hpp"     // custom dialog for text edit

#include "incostem.hpp"		// incoherent STEM calculation
#include "autoslic.hpp"		// CTEM calculation
#include "autostem.hpp"		// ADF-STEM calculation
#include "psave.hpp"		// to load/save parameters
#include "probe.hpp"        //  probe calculation
#include "someAtoms.hpp"	//  init sample specimens

#include <wx/time.h>		// time to init RNG seed
#include <wx/progdlg.h>		// progress dialog
#include <wx/choicdlg.h>	// choices dialog
#include <wx/textctrl.h>	// test edit dialog

//#include <wx/tokenzr.h>		//  for wxStringTokenizer
//#include <wx/txtstrm.h>		//  for wxTextInputStream


//--------------------------------------------------------------------------------

IMPLEMENT_DYNAMIC_CLASS( ctDoc, wxDocument)

BEGIN_EVENT_TABLE( ctDoc, wxDocument )
        EVT_MENU( ID_SETUP_GENERAL, ctDoc::OnSetGeneral )
        EVT_MENU( ID_SETUP_CNM, ctDoc::OnSetAber )
        EVT_MENU( ID_SETUP_SAMPLE, ctDoc::OnSetSample )
        EVT_MENU( ID_EDIT_COORD, ctDoc::OnEditCoord )
        EVT_MENU( ID_GREY_FULL, ctDoc::OnGreyFull )
        EVT_MENU( ID_GREY_ONE, ctDoc::OnGreyOne )
        EVT_MENU( ID_GREY_SET, ctDoc::OnGreySetRange )
        EVT_MENU( ID_GREY_LOG, ctDoc::OnGreyLog )
        EVT_MENU( ID_VIEW_FFT, ctDoc::OnViewFFT )
        EVT_MENU( ID_VIEW_REALSPACE, ctDoc::OnViewRealSpace )
END_EVENT_TABLE()


//--------- ctDoc() constructor/destructor -----------------
ctDoc::ctDoc()
{
        long  ltime;

        Nparams = myFile.maxParam();   // make place to hold image parameters
        param.resize( Nparams );
        sparam.resize( Nparams );
        sparamSet = 0;

        for(int i=0; i<Nparams; i++) param[i] = 0.0F;

        natom = natoms = 0;
        nx = ny = npix = 0;
        nxb = nyb = 0;
        nxProbe = nyProbe = 0;
        ncellx = ncelly = ncellz = 0;

        xiSTEM = xfSTEM = yiSTEM = yfSTEM = 0.0;
        nThick = ndetect = 0;
        n1pixr = n2pixr = n3pixr = 0;

        npix = 0;
        docType = noDoc;

        sphere3d = 0.12F;

		GreyLogConstant = 100.0;  // to scale log values

        //  start with a different RNG seed each time this is run
        ltime = wxGetLocalTime();	// number of seconds since Jan 1, 1970
        iseed = (unsigned) ltime;	// init global seed

		iViewFFT = 0;  // 0 for real space, 1 for FFT
}

ctDoc::~ctDoc()
{
	//  remember:  vector<> will free its memory on exit so we don't have to

	if( n1pixr > 0 ) delete3D<float>( pixr, n1pixr, n2pixr );

 };  //  end ~ctDoc()

//--------- ctDoc() open file -----------------
bool ctDoc::DoOpenDocument(const wxString &  filenam )
{
		//  this needs to match MODE's enum in slicelib.hpp
        char mode[12][32]={ "unknown", "atompot", "mulslice", "image",
                "probe", "stemslice", "autoslic", "autostem", "incostem", "model3d",
				"cbed", "abb-phase" };
        int ix, iy;

        wxString ws;   //  to format messages
        int itest;

        //--- for testing only
        //wxMessageBox( wxT("in ctDoc::DoOpenDocument"), wxT("computem"),
        //wxMessageBox( GetFilename(), wxT("computem"),
        //		wxOK | wxICON_INFORMATION, wxGetApp().GetTopFrame() );

        if( filenam.EndsWith( wxT(".tif" ) ) ||  filenam.EndsWith( wxT(".TIF" ) )  )
        {
                itest = myFile.tFloatTest( filenam.char_str() );
                if( itest > 0 ) {
                        myFile.read( filenam );
                        nx = myFile.nx();
                        ny = myFile.ny();
                        ix = myFile.getParam( pMODE);
                        if( (ix<0) || (ix>11) ) ix = 0;
                        docMessage.Printf(wxT("open file size %d by %d , mode= %s"), nx, ny, mode[ix] );
                        wxLogStatus( docMessage );
                        wxString dt = "written on " + myFile.getDateTime();
                        wxLogStatus( dt );

                        // remember; nx is 2X if complex in this file
                        //  and npix = 1 or 2 for real or complex
                        npix = myFile.getnpix();
                        nx = nx/npix;
                        for(int i=0; i<Nparams; i++) param[i] = myFile.getParam(i);

                        pix.resize( nx, ny );
                        for(iy=0; iy<ny; iy++) for(ix=0; ix<nx; ix++) {
                                pix.re(ix,iy) = myFile( ix, iy );
                                if( npix > 1 ) pix.im(ix,iy) = myFile(ix+nx, iy );
                                else pix.im(ix,iy) = 0.0F;
                        }
                        myFile.resize( 1, 1 );  // free up memory buffer

                        makeBitmap( );
                        docType = tiffDoc;

                } else  wxLogStatus( wxT("Bad file" ) );

        } else if( filenam.EndsWith( wxT(".xyz" ) ) ||  filenam.EndsWith( wxT(".XYZ" ) )  )
        {
                //wxMessageBox( wxT("read xyz file") );  // for testing

                ncellx = ncelly = ncellz = 1;

                natom = ReadXYZcoord( filenam.char_str(), ncellx, ncelly, ncellz,
                        &ax, &by, &cz,
                        Znum, x, y, z, occ, wobble, xyzTitle );

                for(int i=0; i<Nparams; i++) param[i] = 0.0F;

                if( natom > 0 ) {
                        docType = xyzDoc;
                        nx = ny = 256;		//  set default values
                        nxProbe = nyProbe = 256;
                        ncellx = ncelly = ncellz = 3;
                        param[pENERGY] = 100.0F;
                        param[pOAPERT] = 0.020F;
                        param[pDELTAZ] = 2.0F;
                        param[ pMINDET ] = 0.060F;
                        param[ pMAXDET ] = 0.200F;
                        param[ pSOURCE ] = 0.5F;

                        param[ pNWOBBLE ] = 0.0F;
                        param[ pTEMPER ] = 300.0F;

						param[ pPPOSX ] = 1.0F;
						param[ pPPOSY ] = 1.0F;

                        xiSTEM = yiSTEM = 0.0;
                        xfSTEM = yfSTEM = 10.0;

                        nxiSTEM = nyiSTEM = 64;

                        wxLogStatus( ws.Format(wxT("%d atomic coordinates"), natom) );
                }
        } else wxMessageBox( wxT("not a computem file") );

        return true;

}  //  end ctDoc::DoOpenDocument()

// ---------------------------------------------------------------------------
// output file with result
bool ctDoc::DoSaveDocument(const wxString &filenam )
{
    int i, ix, iy, nxout, nyout;
    float rmin, rmax, aimin, aimax, fdx, fdy;

    if( filenam.EndsWith( wxT(".tif" ) ) ||  filenam.EndsWith( wxT(".TIF" ) )  )
    {
        nxout = pix.nx();
        nyout = pix.ny();
        if( (nxout<1) || (nyout<1) ){
                wxMessageBox( wxT("no available data to output") );
                return true;
        }

        //wxMessageBox( ws.Format( wxT("save file nx,ny= %d, %d, npix= %d"), nx, ny, npix ));

    //  remember to use sparam[] (save at calculation) and NOT
    //     param[] which may have been changed after calculation
        for( i=0; i<Nparams; i++ ) myFile.setParam( i, sparam[i]);
        myFile.resize( nxout*npix, nyout );
        myFile.setnpix( npix );
        rmax  = sparam[pRMAX];
        rmin  = sparam[pRMIN];
        aimax = sparam[pIMAX];
        aimin = sparam[pIMIN];
        fdx = sparam[pDX];
        fdy = sparam[pDY];

        //wxMessageBox( ws.Format( wxT("save file range= %f, %f; %f, %f"), rmin, rmax, aimin, aimax ));
        //wxMessageBox( ws.Format(wxT("output file %s" ), filenam ) );

        for( ix=0; ix<nxout; ix++) for( iy=0; iy<nyout; iy++) {
                myFile(ix,iy) = pix.re(ix,iy); // save real part for output
                if( npix > 1 )
                        myFile(ix+nxout,iy) = pix.im(ix,iy); // save imag part for output
        }
        i = myFile.write( filenam, rmin, rmax, aimin, aimax, fdx, fdy );

        if( i != 1 ) wxMessageBox(
            ws.Format(wxT("computem cannot write TIF file %s" ), filenam ) );

    } else wxMessageBox( 
            ws.Format(wxT("can only output tif image files; %s"), filenam ) );

    return true;

};    // end DoSaveDocument()

// ---------------------------------------------------------------------------
//    print current image as raw text
//    mimic display.cpp
//    imag = 0 for real part and imag=1 for imaginary part
//
void ctDoc::DoSavePixText( const wxString &file, int imag )
{
        int ix, iy, nxp, nyp;
        FILE *fp;

        //wxMessageBox( ws.Format(wxT("output file name; %s"), file ) );  // testing

        nxp = pix.nx();		//  need actual calculated size, NOT nx,ny
        nyp = pix.ny();

        if( (nxp<1) && (nyp<1) ) return;  // nothing to save....

        fp = fopen( file.c_str(), "w+" );
        if( NULL == fp ) {
                wxMessageBox( ws.Format(wxT("cannot open file %s for writing"), file ) );
                return;
        }

        for( iy=(nyp-1); iy>=0; iy--) {
                for( ix=0; ix<nxp; ix++)
                    if( imag == 0 ) fprintf( fp, "%12.4g ", pix.re(ix,iy));
                    else fprintf( fp, "%12.4g ", pix.im(ix,iy));
                fprintf( fp, "\n" );
        }

        fclose( fp );

        return;

} //  end ctDoc::DoSavePixText()

// ---------------------------------------------------------------------------
//    save xyz coord as raw text file
//
void ctDoc::DoSaveCoord( const wxString &file )
{
    int i;

    std::ofstream fp;
    std::string sfile( file.char_str() );  // need to convert to regular char[] in .open()

    fp.open( sfile.c_str() );
    if( fp.fail() ) {
        wxMessageBox( ws.Format(wxT("cannot open file %s"), file ) );
    } else {

        fp << xyzTitle << endl;
        //fp << setw(14) << ax << " " << setw(14) << by << " " << setw(14) << cz << endl;
        // below looks a little better
        wxString ws1 = ws.Format(wxT("%12.6f %12.6f %12.6f \n"), ax, by, cz);
        std::string stt1(ws1.char_str() );
        fp << stt1;

        for( i=0; i<natom; i++) {
            //fp << setw(4) << Znum[i] << " " << setw(14) << x[i] 
            //    << " " << setw(14) << y[i] << " " << setw(14)  << z[i]
            //    << " " << setw(14)  << occ[i] 
            //    << " " << setw(14) << wobble[i] << endl;
            //  below looks better
            wxString ws2 = ws.Format(wxT("%4d %12.6f %12.6f %12.6f %12.6f %12.6f"),
                Znum[i], x[i], y[i], z[i], occ[i], wobble[i]);
            std::string stt2(ws2.char_str() );
            fp << stt2 << endl;
        }
        i = -1;
        fp << setw(4) << i << endl;
    }

    return;
}  // end DoSaveCoord()

// ---------------------------------------------------------------------------
//    save param as raw text file
//
void ctDoc::DoSaveParam( const wxString &file )
{
        psave ps;
        ps.save( file.c_str(), param, Nparams );
        return;
}
// ---------------------------------------------------------------------------
//    load param from raw text file
//
void ctDoc::DoLoadParam( const wxString &file )
{
        psave ps;
        ps.load( file.c_str(), param, Nparams );
        return;
}

//--------- ctDoc() new file -----------------
//   include coord. of some specimens for lazy users
bool ctDoc::OnNewDocument()
{
        int ia, nselect, ntotal;

		// temp. expandable containers to get atomi coord. data
		vector<float> xv, yv, zv, occv, wobblev;
		vector<int> Znumv;

        wxArrayString choices;
        choices.Add( wxT("Si110") );
        choices.Add( wxT("SrTiO") );
        choices.Add( wxT("graphene") );
        choices.Add( wxT("YBCO") );

        wxSingleChoiceDialog dlg( wxGetApp().GetTopFrame(), wxT("available specimen coordinates"),
                wxT("select built in specimen descriptions"), choices );

        dlg.SetSelection( 0 );
        if( dlg.ShowModal() == wxID_OK) {
                //wxMessageBox( dlg.GetStringSelection() );   // testing ????

                nselect = dlg.GetSelection();
                SetDocumentName( choices.Item(nselect) );  // none of these change the window
                SetFilename( choices.Item(nselect) );      // name (why ???)
                //SetTitle( choices.Item(nselect) );  // does not work here ????
				//  change title in "new" window
				((wxTopLevelWindow*)GetFirstView()->GetFrame())->SetTitle( choices.Item(nselect) );

                wxLogStatus( ws.Format(wxT("new specimen using: %s"), choices.Item(nselect)) );

                docType = xyzDoc;
                nx = ny = 256;		//  set default values
                nxProbe = nyProbe = 256;
                ncellx = ncelly = ncellz = 3;
                param[pENERGY] = 100.0F;
                param[pOAPERT] = 0.020F;
                param[pDELTAZ] = 2.0F;
                param[ pMINDET ] = 0.060F;
                param[ pMAXDET ] = 0.200F;
                param[ pSOURCE ] = 0.5F;

                param[ pNWOBBLE ] = 0.0F;
                param[ pTEMPER ] = 300.0F;

				param[ pPPOSX ] = 1.0F;
				param[ pPPOSY ] = 1.0F;

                xiSTEM = yiSTEM = 0.0;
                xfSTEM = yfSTEM = 10.0;

                nxiSTEM = nyiSTEM = 64;

				//  put this long dataset in a separate file so this is easier to follow
				natom = someAtoms<float>( nselect, ncellx, ncelly, ncellz,
					ax, by, cz, Znumv, xv, yv, zv, 
					occv, wobblev, xyzTitle );

                wxLogStatus( ws.Format(wxT("%d atomic coordinates"), natom) );

                ntotal = natom;
                x.resize( ntotal );
                y.resize( ntotal );
                z.resize( ntotal );
                occ.resize( ntotal );
                wobble.resize( ntotal );
                Znum.resize( ntotal );

				for( ia=0; ia< natom; ia++) {  //  copy new atom coord. to good place
					x[ia] = xv[ia];
					y[ia] = yv[ia];
					z[ia] = zv[ia];
					occ[ia] = occv[ia];
					wobble[ia] = wobblev[ia];
					Znum[ia] = Znumv[ia];
				}   // end for( ia....)

        }  // end if( dlg.ShowModal() ... )

    return true;
}  //  end OnNewDocument()

//--------- expandCell -----------------
//  expand unit cell into a bigger supercell
void ctDoc::expandCell( )
{
        int i, j, k, na;

        if( (ncellx < 1) || (ncelly < 1) || (ncellz < 1 ) ) return;

        na = natom * ncellx * ncelly * ncellz;

        //  delete the old arrays and make new ones (if not same size)
        if( natoms != na ) {
            xs.resize( na );
            ys.resize( na );
            zs.resize( na );
            occs.resize( na );
            wobbles.resize( na );
            Znums.resize( na );
            natoms = na;
        }

        //  copy unit cell
        for( i=0; i<natom; i++) {
                xs[i] = x[i];
                ys[i] = y[i];
                zs[i] = z[i];
                occs[i] = occ[i];
                wobbles[i] = wobble[i];
                Znums[i] = Znum[i];
        }
        axs = ax * ncellx;
        bys = by * ncelly;
        czs = cz * ncellz;
        
        //   expand one direction at a time
        na = natom;
        if( ncellx > 1 ) {
            for( i=1; i<ncellx; i++)
            for( j=0; j<na; j++) {
                k = j + na*i;
                xs[ k ]  = xs[j] + i*ax;
                ys[ k ]  = ys[j];
                zs[ k ]  = zs[j];
                occs[ k ]     = occs[j];
                wobbles[ k ]  = wobbles[j];
                Znums[ k ]    = Znums[j];
            }
            na = na * ncellx;
        }

        if( ncelly > 1 ) {
            for( i=1; i<ncelly; i++)
            for( j=0; j<na; j++) {
                k = j + na*i;
                xs[ k ]  = xs[j];
                ys[ k ]  = ys[j] + i*by;
                zs[ k ]  = zs[j];
                occs[ k ]     = occs[j];
                wobbles[ k ]  = wobbles[j];
                Znums[ k ]    = Znums[j];
            }
            na = na * ncelly;
        }

        if( ncellz > 1 ) {
            for( i=1; i<ncellz; i++)
            for( j=0; j<na; j++) {
                k = j + na*i;
                xs[ k ]  = xs[j];
                ys[ k ]  = ys[j];
                zs[ k ]  = zs[j] + i*cz;
                occs[ k ]     = occs[j];
                wobbles[ k ]  = wobbles[j];
                Znums[ k ]    = Znums[j];
            }
            na = na * ncellz;
        }

        return;
};  //  end expandCell( )


// ---------------------------------------------------------------------------
//  make a bitmap that can be sent to the screen
//  from internal calculated image
//  also sets globals nxb, nyb = size of bitmap after interpolation
//
//  here nxp is for the complex portion (not doubled for real + imag)
//
//  remember nx,ny are the transmission function size but
//  STEM image will be different
//
//  LogGreyScale sets log (1) or linear (0) greyscale
//
//  assumed global pix, pixFFT, iViewFFT
//
void ctDoc::makeBitmap( int LogGreyScale )
{
    int ix, iy, ltemp, j, nxp, nyp;
    float scaler, scalei, scaleLog, dx, dy, rmin, rmax, imin, imax;
    double x, y, z, c;

    //  interpolate into a bitmap with square pixels for display 
    //  modify code from floatTIFF to 

    if( (npix<1) || (npix>2)) return;   // just in case

    nxp = pix.nx();
    nyp = pix.ny();

	cfpix px;					/// ?????  - copy orig argument 15-apr-2017 ejk
	px.resize(nxp, nyp);
	if( 0 == iViewFFT ) px = pix;
	else px = pixFFT;

    dx = param[ pDX ];  // pixel size
    dy = param[ pDY ];
	if( iViewFFT !=  0 ){
		dx = 1.0F/( nxp * param[ pDX ] );
		dy = 1.0F/( nyp * param[ pDY ] );
	}

    if( dx <= 0.0F ) dx = 1.0F;
    if( dy <= 0.0F ) dy = dx;
    if( dx > dy ) {
        nxb = (int) ( nxp * dx/dy + 0.5F);
        nyb = nyp;
    } else if( dx < dy ) {
        nxb = nxp;
        nyb = (int) ( nyp * dy/dx + 0.5F);;
    } else {
        nxb = nxp;
        nyb = nyp;
    }

    nxb = nxb * npix;    //  npix=1 for real and 2 for complex, put side by side

    //wxMessageBox( ws.Format(wxT("nxb, nyb = %d, %d"), nxb, nyb ) );  // testing
    //wxMessageBox( ws.Format(wxT("dx, dy = %f, %f"), dx, dy ) );

    // first make a device independent image and then copy to a bitmap
    // remember: can't access data in a raw bitmap, only a wxImage
    ctIMG.Create( nxb, nyb );
    pixData = ctIMG.GetData();

    rmin = param[ pRMIN ];  // image range
    rmax = param[ pRMAX ]; 
    imin = param[ pIMIN ]; 
    imax = param[ pIMAX ];

    if( (rmax-rmin) > 1.0e-25 ) scaler = 255.0F/(rmax-rmin);
    else scaler= 1.0f;
    if( (npix > 1) && ((imax-imin)>1.0e-20) ) scalei = 255.0F/(imax-imin);
    else scalei = 1.0F;

	c = GreyLogConstant;  //  log scale constant
	scaleLog = (float) ( 255.0/log(1.0 + c*255.0) );

    for( iy=0; iy<nyb; iy++){
       y = (nyp-1) * ((double)iy)/((double)(nyb-1)) ;
       for( ix=0; ix<nxb/npix; ix++) {
          x = (nxp-1) * ((double)ix)/((double)(nxb/npix-1)) ;
          z = tinter( px, x, y, 0);
		  if( 0 == LogGreyScale ) 
			ltemp = (int) ( scaler*(z-rmin) + 0.5 );
		  else if ( 1 == LogGreyScale ) 
			ltemp = (int)( scaleLog*(log(1.0+ c*fabs(scaler*(z-rmin)) ) ) );
          if( ltemp < 0 ) ltemp = 0;
          if( ltemp > 255 ) ltemp = 255;
          j = 3*( ix + iy*nxb);  // index of red pixel for real part
          pixData[j] = (unsigned char) (ltemp & 0x00FF);	//red
          pixData[j+1] = pixData[j];				// green
          pixData[j+2] = pixData[j];				// blue
          if( npix >1 ) {
                  z = tinter( px, x, y, 1);
                  ltemp = (int) ( scalei*(z-imin) + 0.5 );
                  if( ltemp < 0 ) ltemp = 0;
                  if( ltemp > 255 ) ltemp = 255;
                  j = 3*( ix +nxb/npix + iy*nxb);  // index of red pixel for imag part
                  pixData[j] = (unsigned char) (ltemp & 0x00FF);	//red
                  pixData[j+1] = pixData[j];				// green
                  pixData[j+2] = pixData[j];				// blue
          }
       }  // end for ix loop

    }  // end for iy loop

    ctBMP = new wxBitmap( ctIMG );   //  make bitmap from image and save it for redrawing
    ctIMG.Destroy();  //  free up memory in image

}  // end makeBitmap()

//--------- OnEditCoord -----------------
void ctDoc::OnEditCoord( wxCommandEvent& event )
{
    int nl, nav, znv;
    float axv, byv, czv, xx, yy, zz, oo, ww;

    std::vector<int> Znumv;
    std::vector<float> xv, yv, zv, occv, wobblev;

    textDialog td( NULL );  //  text dialog

    td.tx.push_back( xyzTitle );    //  fill in the text to edit

    std::stringstream strs;
    //strs <<  setw(14) << ax << " " << setw(14) << by << " " << setw(14) << cz;
    //td.tx.push_back( strs.str() );
    //  below format looks a little better
    wxString ws1 = ws.Format(wxT("%12.6f %12.6f %12.6f"), ax, by, cz);
    std::string stt1(ws1.char_str() );
    td.tx.push_back( stt1 );

    for( int i=0; i<natom; i++) {
        //strs.str("");
        //strs <<  setw(4) << Znum[i] << " " << setw(14) << x[i] 
        //        << " " << setw(14) << y[i] << " " << setw(14)  << z[i]
        //        << " " << setw(14)  << occ[i] 
        //        << " " << setw(14) << wobble[i];
        //td.tx.push_back( strs.str() );
        //  above gives inconsistent spacing, below is better
        wxString ws2 = ws.Format(wxT("%4d %12.6f %12.6f %12.6f %12.6f %12.6f"),
            Znum[i], x[i], y[i], z[i], occ[i], wobble[i]);
        std::string stt2(ws2.char_str() );
        td.tx.push_back( stt2 );
    }

    if( wxID_OK == td.ShowModal() ) { //  show the dialog box

        nl = td.nlines;
        if( (nl > 0) && (td.modified > 0) ) {
            xyzTitle = td.tx[0];
            //std::istringstream isbuf( td.tx[1] );
            //isbuf >> axv >> byv >> czv;
            strs.str( td.tx[1] );
            strs >> axv >> byv >> czv;

            //  parse rest of text lines and count number of coord.
            //     remember; MSVC 2010 does funny things with istringstream
            //     so make it this way to get it too work
            nav = 0;
            for( int i=2; i<nl; i++) {
                if( td.tx[i].length() < 3 ) break;  // skip blank line- usually at end
                std::istringstream isbuf( td.tx[i] );
                //isbuf.str( td.tx[i] );  //  doesn't work here ??? 
                znv = 0;  // just in case the next line fails
                xx = yy = zz = oo = ww = 0;
                isbuf >> znv >> xx >> yy >> zz >> oo >> ww;
                //wxLogStatus( ws.Format(wxT("read =%d %g %g %g"), znv, xx, yy, zz) ); // testing
                if( (znv >= 1) && (znv <= 103) ) {
                    nav += 1;
                    Znumv.push_back( znv );
                    xv.push_back( xx );
                    yv.push_back( yy );
                    zv.push_back( zz );
                    occv.push_back( oo );
                    wobblev.push_back( ww );
               }  else break;  //  read until negative Znum
            }  //  end for( int i=2.... 
            wxLogStatus( ws.Format(wxT("%d atomic coordinates edited"), nav) );

            //  copy back to x,y,z coordinates 
			//   ????  don't need two copies here anymore - should fix this sometime
            if( nav > 0 ){
                x.resize( nav );
                y.resize( nav );
                z.resize( nav );
                occ.resize( nav );
                wobble.resize( nav );
                Znum.resize( nav );
                for( int j=0; j< nav; j++ ) {
                    x[j] = xv[j];
                    y[j] = yv[j];
                    z[j] = zv[j];
                    occ[j] = occv[j];
                    wobble[j] = wobblev[j];
                    Znum[j] = Znumv[j];
                }
                natom = nav;
                ax = axv;
                by = byv;
                cz = czv;
            }  // end if( nav ...

        } // end if( nl > 0 ....

    } // end if(... ShowModal() )

    return;

}  // end OnEditCoord()

//--------- OnSetAber -----------------
void ctDoc::OnSetAber( wxCommandEvent& event )
{
        double keV, wavl;

        aberDialog dlg( NULL );

        //  transfer current parameters to the dialog for display/edit
        dlg.df = param[ pDEFOCUS ];
        dlg.ddf = param[ pDDF ];
        dlg.C12a = param[ pC12a ];
        dlg.C12b = param[ pC12b ];

        dlg.C21a = param[ pC21a ];
        dlg.C21b = param[ pC21b ];
        dlg.C23a = param[ pC23a ];
        dlg.C23b = param[ pC23b ];

        dlg.C30  = param[ pCS ];
        dlg.C32a = param[ pC32a ];
        dlg.C32b = param[ pC32b ];
        dlg.C34a = param[ pC34a ];
        dlg.C34b = param[ pC34b ];

        dlg.C41a = param[ pC41a ];
        dlg.C41b = param[ pC41b ];
        dlg.C43a = param[ pC43a ];
        dlg.C43b = param[ pC43b ];
        dlg.C45a = param[ pC45a ];
        dlg.C45b = param[ pC45b ];

        dlg.C50  = param[ pCS5 ];
        dlg.C52a = param[ pC52a ];
        dlg.C52b = param[ pC52b ];
        dlg.C54a = param[ pC54a ];
        dlg.C54b = param[ pC54b ];
        dlg.C56a = param[ pC56a ];
        dlg.C56b = param[ pC56b ];

		//   below need to calc. random tuning errors
		dlg.iseed = &iseed;   // awkward but seems the only way (?)
		dlg.objApert = param[ pOAPERT ];
		keV = param[pENERGY];
		wavl = param[ pWAVEL ] = wavelength( keV );
		dlg.wavl = param[ pWAVEL ];

        //wxLogStatus( ws.Format(wxT("orig df= %f"), param[ pDEFOCUS ])  );  // ???? testing

        if( wxID_OK == dlg.ShowModal() )  //  show the dialog box
        {
                if( dlg.nerror > 0 )   //  if some values were not numbers
                        wxMessageBox( ws.Format(wxT("%d error(s) in Cnm values"), dlg.nerror),
                                wxT("setup aberrations") );
                if( docType == xyzDoc ) {
                        //  get new values back from the dialog
                        param[ pDEFOCUS ] = dlg.df;
                        param[ pDDF ] = dlg.ddf;
                        param[ pC12a ] = dlg.C12a;
                        param[ pC12b ] = dlg.C12b;

                        param[ pC21a ] = dlg.C21a;
                        param[ pC21b ] = dlg.C21b;
                        param[ pC23a ] = dlg.C23a;
                        param[ pC23b ] = dlg.C23b;

                        param[ pCS ]  = dlg.C30;
                        param[ pC32a ] = dlg.C32a;
                        param[ pC32b ] = dlg.C32b;
                        param[ pC34a ] = dlg.C34a;
                        param[ pC34b ] = dlg.C34b;

                        param[ pC41a ] = dlg.C41a;
                        param[ pC41b ] = dlg.C41b;
                        param[ pC43a ] = dlg.C43a;
                        param[ pC43b ] = dlg.C43b;
                        param[ pC45a ] = dlg.C45a;
                        param[ pC45b ] = dlg.C45b;

                        param[ pCS5 ]  = dlg.C50;
                        param[ pC52a ] = dlg.C52a;
                        param[ pC52b ] = dlg.C52b;
                        param[ pC54a ] = dlg.C54a;
                        param[ pC54b ] = dlg.C54b;
                        param[ pC56a ] = dlg.C56a;
                        param[ pC56b ] = dlg.C56b;
                }
        } // end if(... ShowModal() )

        return;
};  // end OnSetAber()

//--------- OnSetSample -----------------
void ctDoc::OnSetGeneral( wxCommandEvent& event )
{
        generalDialog dlg( NULL );

        dlg.kev = param[ pENERGY ];
        dlg.objAp = param[ pOAPERT ];
        dlg.capertmin = param[ pCAPERTMIN ];
        dlg.capertmax = param[ pCAPERT ];
        dlg.detmin = param[ pMINDET ];
        dlg.detmax = param[ pMAXDET ];
        dlg.dsource = param[ pSOURCE ];
        dlg.probeI = param[ pPROBEI ];
        dlg.probeDt = param[ pPROBEDT ];
        dlg.ctiltx = param[ pXCTILT ];
        dlg.ctilty = param[ pYCTILT ];
        dlg.probePosx = param[ pPPOSX ];
        dlg.probePosy = param[ pPPOSY ];
        dlg.sphere3d = sphere3d;
		dlg.GreyLogConstant = GreyLogConstant;

        if( wxID_OK == dlg.ShowModal() )  //  show the dialog box
        {		
            if( dlg.nerror > 0 )   //  if some values were not numbers
                wxMessageBox( ws.Format(wxT("%d error(s) in general values"), dlg.nerror),
                                wxT("general sampling") );
				GreyLogConstant = dlg.GreyLogConstant;
                //  don't allow changes to old data files- just new xyz files to calculate
                if( docType == xyzDoc ) {
                        param[ pENERGY ] = dlg.kev;
                        param[ pOAPERT ] = dlg.objAp;
                        param[ pCAPERTMIN ] = dlg.capertmin;
                        param[ pCAPERT ] = dlg.capertmax;
                        param[ pMINDET ] = dlg.detmin;
                        param[ pMAXDET ] = dlg.detmax;
                        param[ pSOURCE ] = dlg.dsource;
						param[ pPROBEI ] = dlg.probeI;
						param[ pPROBEDT ] = dlg.probeDt;
                        param[ pXCTILT ] = dlg.ctiltx;
                        param[ pYCTILT ] = dlg.ctilty;
						param[ pPPOSX ] = dlg.probePosx;
						param[ pPPOSY ] = dlg.probePosy;
                        sphere3d = dlg.sphere3d;
                }
        }

        return;
}  // end OnSetGeneral()

//--------- OnSetSample -----------------
void ctDoc::OnSetSample( wxCommandEvent& event )
{
        sampleDialog dlg( NULL );

        dlg.nx = nx;
        dlg.ny = ny;

        dlg.nxProbe = nxProbe;
        dlg.nyProbe = nyProbe;
        dlg.xi = xiSTEM;
        dlg.xf = xfSTEM;
        dlg.yi = yiSTEM;
        dlg.yf = yfSTEM;
        dlg.nxi = nxiSTEM;
        dlg.nyi = nyiSTEM;

        dlg.ncellx = ncellx;
        dlg.ncelly = ncelly;
        dlg.ncellz = ncellz;

        dlg.deltaz = param[pDELTAZ];

        dlg.nwobble = ToInt( param[pNWOBBLE] );
        dlg.temp = param[pTEMPER];

        //wxLogStatus( ws.Format(wxT("orig df= %f"), param[ pDEFOCUS ])  );  // ???? testing

        if( wxID_OK == dlg.ShowModal() )  //  show the dialog box
        {
                if( dlg.nerror > 0 )   //  if some values were not numbers
                        wxMessageBox( ws.Format(wxT("%d error(s) in sampling values"), dlg.nerror),
                                wxT("setup sampling") );

                //  don't allow changes to old data files- just new xyz files to calculate
                if( docType == xyzDoc ) {
                        nx = dlg.nx;
                        ny = dlg.ny;
                        nxProbe = dlg.nxProbe;
                        nyProbe = dlg.nyProbe;
                        ncellx = dlg.ncellx;
                        ncelly = dlg.ncelly;
                        ncellz = dlg.ncellz;
                        param[pDELTAZ] = dlg.deltaz;
                        xiSTEM = dlg.xi;
                        xfSTEM = dlg.xf;
                        yiSTEM = dlg.yi;
                        yfSTEM = dlg.yf;
                        nxiSTEM = dlg.nxi;
                        nyiSTEM = dlg.nyi;

                        param[pNWOBBLE] = dlg.nwobble;
                        param[pTEMPER] = dlg.temp;
                }
        }

        return;
};  //  end OnSetSample()

//--------- OnGreyFull -----------------
//   set greyscale to full range
//  assumed global iViewFFT
void ctDoc::OnGreyFull( wxCommandEvent& event )
{
    float rmin, rmax, aimin, aimax;

    if( (npix<1) || (npix>2)) return;   // just in case

	if( 0 == iViewFFT ) pix.findRange( rmin, rmax, aimin, aimax );
	else pixFFT.findRange( rmin, rmax, aimin, aimax );
	param[ pRMIN ] = rmin;
	param[ pRMAX ] = rmax;

    wxLogStatus( ws.Format(wxT("rescale real pix to %g to %g"), rmin, rmax)  );

	if( npix > 1 ) {   //  if there is an imaginary part
		param[ pIMIN ] = aimin;
		param[ pIMAX ] = aimax;
		wxLogStatus( ws.Format(wxT("rescale imag. pix to %g to %g"), aimin, aimax)  );
	}

	makeBitmap( );
	UpdateAllViews();

    return;
};  //  end OnGreyFull()


//--------- OnGreyOne -----------------
//   expand greyscale to mean +/- one stnd. dev.
//  assumed global iViewFFT
void ctDoc::OnGreyOne( wxCommandEvent& event )
{
    int ix, iy, nxp, nyp;
    float rmin, rmax, imin, imax, fr, fi;
	double r, r2, ai, ai2, scale;

    if( (npix<1) || (npix>2)) return;   // just in case

    nxp = pix.nx();
    nyp = pix.ny();

	//  find the mean and stnd. dev.
	r = r2 = 0.0;
	for( iy=0; iy<nyp; iy++)
	for( ix=0; ix<nxp; ix++) {	// start with real part
		if( 0 == iViewFFT ) fr = pix.re(ix,iy);
		else fr = pixFFT.re(ix,iy);
		r += fr;    r2 += fr*fr;
	}
	scale = 1.0/(nxp*nyp);
	r *= scale;					//  mean of real part
	r2 = sqrt(r2*scale - r*r);	//  stnd dev. of real part
	param[pRMAX] = rmax = r + r2;
	param[pRMIN] = rmin = r - r2;

    wxLogStatus( ws.Format(wxT("rescale real pix to %f to %f"), rmin, rmax)  );

	if( npix > 1 ) {   //  if there is an imaginary part
		ai = ai2 = 0.0;
		for( iy=0; iy<nyp; iy++)
		for( ix=0; ix<nxp; ix++) {
			if( 0 == iViewFFT ) fi = pix.im(ix,iy);
			else fi = pixFFT.im(ix,iy);
			ai += fi;    ai2 += fi*fi;
		}
		ai *= scale;					//  mean of imag part
		ai2 = sqrt(ai2*scale -ai*ai);   //  stnd dev of imag. part
		param[pIMAX] = imin = ai + ai2;
		param[pIMIN] = imax = ai - ai2;

		wxLogStatus( ws.Format(wxT("rescale imag. pix to %f to %f"), imin, imax)  );
	}

	makeBitmap( );
	UpdateAllViews();

    return;
};  //  end OnGreyOne()

//--------- OnGreySetRange -----------------
//   set greyscale to manually requested range
void ctDoc::OnGreySetRange( wxCommandEvent& event )
{
    float rmin, rmax, aimin, aimax;

    rangeDialog dlg( NULL );

    if( (npix<1) || (npix>2)) return;   // just in case

	dlg.npix = npix;  //  real or complex

	dlg.rmin = param[pRMIN];
	dlg.rmax = param[ pRMAX ];
	dlg.aimin = param[ pIMIN ];
	dlg.aimax = param[ pIMAX ];

    if( wxID_OK == dlg.ShowModal() )  //  show the dialog box
    {
            if( dlg.nerror > 0 )   //  if some values were not numbers
                wxMessageBox( ws.Format(wxT("%d error(s) in range values"), dlg.nerror),
                        wxT("setup sampling") );
			if( dlg.rmax > dlg.rmin ) {
				rmin = dlg.rmin;
				rmax = dlg.rmax;
				wxLogStatus( ws.Format(wxT("rescale real pix to %g to %g"), rmin, rmax)  );
				param[ pRMIN ] = rmin;
				param[ pRMAX ] = rmax;
			}
			if( (npix>1) && ( dlg.aimax > dlg.aimin) ) {
				aimin = dlg.aimin;
				aimax = dlg.aimax;
				wxLogStatus( ws.Format(wxT("rescale imag. pix to %g to %g"), aimin, aimax)  );
				param[ pIMIN ] = aimin;
				param[ pIMAX ] = aimax;
			}
	}

	makeBitmap( );
	UpdateAllViews();

    return;
};  //  end OnGreySetRange()


//--------- OnGreyLog -----------------
//   convert greyscale logarithmic
void ctDoc::OnGreyLog( wxCommandEvent& event )
{
    int nxp, nyp;
    float rmin, rmax, imin, imax;

    if( (npix<1) || (npix>2)) return;   // just in case

    nxp = pix.nx();
    nyp = pix.ny();

	rmax = param[pRMAX];
	rmin = param[pRMIN];
	imin = param[pIMAX];
	imax = param[pIMIN];

	if( (npix != 1 ) || (rmin < 0.0F) ){
		wxMessageBox( wxT("cannot display complex or negative images on log scale"),
				wxT("computem") );
		return;
	}

    wxLogStatus( ws.Format(wxT("display real pix on log scale of %f to %f"), rmin, rmax)  );

	makeBitmap( 1 );   // use log scale
	UpdateAllViews();

    return;
};  //  end OnGreyLog()

//--------- OnViewFFT -----------------
//   set to FFT
void ctDoc::OnViewFFT( wxCommandEvent& event )
{
	int nx, ny, ix, iy, mode;
    float rmin, rmax, aimin, aimax, xr, yi;

    if( (npix<1) || (npix>2)) return;   // just in case

	mode = ToInt( param[pMODE] );
	if( !((mode==mMULSLICE)||(mode==mIMAGE)
		||(mode==mPROBE)||(mode==mSTEMSLICE)||(mode==mAUTOSLICE)
		||(mode==mAUTOSTEM)||(mode==mINCOSTEM) )  ){
		wxLogStatus( wxT("sorry can't fft this image type") );
		return;
	}

	nx = pix.nx();
	ny = pix.ny();
	pixFFT.resize( nx, ny );   // need a permenant copy of this to change greyscale

	pixFFT.init( 1 );  // if this is an old file read back there won't be an init
	pixFFT = pix;

	pixFFT.fft();    //  calcualte mag of FFT
	pixFFT.invert2D();
	if( npix == 1 ) {   //  if there is no imaginary part
		for( ix=0; ix< nx; ix++) for( iy=0; iy<ny; iy++) {
			xr = pixFFT.re(ix,iy);  yi= pixFFT.im(ix,iy);
			pixFFT.re(ix,iy) = log( 0.1 + sqrt( xr*xr + yi*yi ) );  // arbitray scale (????)
			pixFFT.im(ix,iy) = 0;
		}
	}

	pixFFT.findRange( rmin, rmax, aimin, aimax );
	param[ pRMIN ] = rmin;
	param[ pRMAX ] = rmax;
	if( npix > 1 ) {   //  if there is an imaginary part
		param[ pIMIN ] = aimin;
		param[ pIMAX ] = aimax;
	}

	wxLogStatus(wxT("display log(0.1+|fft|);") );
    wxLogStatus( ws.Format(wxT("fft pix scale %g to %g"), rmin, rmax)  );

	iViewFFT = 1;    // probably should be an arg. in makeBitmap() 
	makeBitmap( );
	UpdateAllViews();

    return;
};  //  end OnViewFFT()

//--------- OnViewRealSpace -----------------
//   set view to real space 
void ctDoc::OnViewRealSpace( wxCommandEvent& event )
{
    float rmin, rmax, aimin, aimax;

	iViewFFT = 0;  //  reset and repeat fullGreyScale()

    if( (npix<1) || (npix>2)) return;   // just in case

	pix.findRange( rmin, rmax, aimin, aimax );
	param[ pRMIN ] = rmin;
	param[ pRMAX ] = rmax;

    wxLogStatus( ws.Format(wxT("rescale real pix to %g to %g"), rmin, rmax)  );

	if( npix > 1 ) {   //  if there is an imaginary part
		param[ pIMIN ] = aimin;
		param[ pIMAX ] = aimax;
		wxLogStatus( ws.Format(wxT("rescale imag. pix to %g to %g"), aimin, aimax)  );
	}

	makeBitmap( );
	UpdateAllViews();

    return;
};  //  end OnViewRealSpace()


//--------- RUNabbPhase -----------------
void ctDoc::RUNabbPhase( )
{
        int lstart=0, lbeams=0, lcross=0;
        int nbout = 0;  // number of beams to save
        int nwobble=0;

        float rmin, rmax, aimin, aimax;
        float t1, t2, BW=(2.0F/3.0F);

        double keV, wavl;
        int multiMode = 1;

        autoslic aslice;

        if( docType != xyzDoc ) wxMessageBox( wxT("can only calculate from xyz files") );
        else if( (nx<64) && (ny<64) ){
                wxMessageBox( wxT("nx,ny must be larger for abbPhase") );
                return;
		} else {

            wxLogStatus( wxT("calculate the abber. phase in 2D wrapped into -pi to +pi")  );

            expandCell();  //  make supercell from one unit cell

            keV = param[pENERGY];
            wavl = wavelength( keV );

            //   set calculation parameters (some already set in dialogs)
            param[pAX] = axs;		// supercell size
            param[pBY] = bys;
            param[pDX] = (float) (1.0/axs);  // ( axs/((float)nx) );
            param[pDY] = (float) (1.0/bys);  // ( bys/((float)ny) );
            param[ pNX ] = (float) nx;
            param[ pNY ] = (float) ny;
            param[ pMODE ] = mABBPHASE;  // mode = abbPhase

            param[ pWAVEL ] = wavl;

            // repeat calculation max sampling angles for echo 
            t1 = nx/(2.0F*axs);
            t2 = ny/(2.0F*bys);
            if( t2 < t1 ) t1 = t2;
            t1 = wavl * BW * t1;  //  

            wxLogStatus( ws.Format(wxT("max angle =  %f mrad"), t1*1000.0F)  );

            aslice.abbPhase2D( pix, param, multiMode );

            pix.findRange( rmin, rmax, aimin, aimax );
            param[ pRMIN ] = rmin;
            param[ pRMAX ] = rmax;
            param[ pIMIN ] = aimin;
            param[ pIMAX ] = aimax;

            npix = 1;	// real for partial coherence mode
            docMessage.Printf(wxT("phase abb. func.** range %f to %f"), rmin, rmax);
            wxLogStatus( docMessage );

            makeBitmap( );
            UpdateAllViews();

            //   save the parameters actually used
            for( int i=0; i<Nparams; i++) sparam[i] = param[i];
            sparamSet = 1;

            wxLogStatus( wxT("***phase abb. func. finished")  );
        }
        return;

}  //  end RUNabbPhase(  )

//--------- RUNautoslic -----------------
//  CTEM calculation
//  can calculate CTEM image, exit wave or diffraction
//  flags are;
//  lpartl = partial coherent image
//  lwobble = phonons
//  
void ctDoc::RUNautoslic( int lpartl, int lwobble )
{
        int lstart=0, lbeams=0, lcross=0;
        int nbout = 0;  // number of beams to save
        int nwobble=0;
		int iverbose;

        vectori hbeam(10), kbeam(10);  //  not really used but need variables 

        float rmin, rmax, aimin, aimax;

        double keV, wavl;
        int multiMode = 1;

        autoslic aslice;

        cfpix beams, wave0, depthpix;   // not really used here but required

        if( docType != xyzDoc ) wxMessageBox( wxT("can only calculate from xyz files") );
        else if( (nx<64) && (ny<64) ){
                wxMessageBox( wxT("nx,ny must be larger for autoslic") );
                return;

        } else if( param[ pCAPERT ] < param[ pCAPERTMIN ] ){
                wxMessageBox( wxT("bad condencer aperture size in autoslic") );
                return;

        } else {
                //wxMessageBox( wxT("calculate autoslic") );  // testing

            expandCell();  //  make supercell from one unit cell

            wxLogStatus( ws.Format(wxT("autoslic with supercell of %d atoms"), natoms)  );

            keV = param[pENERGY];
            wavl = wavelength( keV );

            nwobble = ToInt( param[pNWOBBLE] );

            /* ---------  setup calculation ----- */
            //   set calculation flags
            aslice.lbeams = lbeams;
            aslice.lcross = lcross;
            aslice.lpartl = lpartl;
            aslice.lstart = lstart;
            if( nwobble > 0 ) aslice.lwobble = 1;
            else  aslice.lwobble = 0;

            // need to get these from somewhere ????
            double ycross = 0.0;	// not used ye
            double dfdelt = 2.0;   // defocus integration step size

            //   set calculation parameters (some already set in dialogs)
            param[pAX] = axs;		// supercell size
            param[pBY] = bys;
            param[pDX] = (float) ( axs/((float)nx) );
            param[pDY] = (float) ( bys/((float)ny) );
            param[ pNX ] = (float) nx;
            param[ pNY ] = (float) ny;
            param[ pMODE ] = mAUTOSLICE;  // mode = autoslic

            param[ pWAVEL ] = wavl;		//  probably recal. in autoslice::calculate()

            //if ( lpartl == 1 ) {
                    //rx = 1.0F/ax;
                    //ry = 1.0F/by;
                    //printf("Illumination angle sampling (in mrad) = %f, %f\n\n",
                    //	1000.*rx*wavlen, 1000.*ry*wavlen);
            //}

            // ------- iterate the multislice algorithm proper -----------
			//
            if ( 0 == lpartl ) {
				iverbose = 1;  // echo steps
				aslice.initAS( param, Znum, natom );
				pix.resize(nx,ny);
				pix.init();  // required if this is the first operation
				aslice.calculate( pix, wave0, depthpix, param, multiMode, natoms,
                        Znums,xs,ys,zs,occs, beams, hbeam, kbeam, nbout,
						ycross, iverbose );
			} else if ( 1 == lpartl ) {
				nbout = 0;
				aslice.calculatePartial( pix, param, multiMode, natoms, &iseed,
						Znums, xs,ys,zs,occs,wobbles, dfdelt );
			}
 
            pix.findRange( rmin, rmax, aimin, aimax );
            param[ pRMIN ] = rmin;
            param[ pRMAX ] = rmax;
            param[ pIMIN ] = aimin;
            param[ pIMAX ] = aimax;

            if( 0 == lpartl ) {
				npix = 2;       // complex for exit wave
				docMessage.Printf(wxT("**exit wave** range %f to %f and %f to %f"), 
							rmin, rmax, aimin, aimax ) ;
				wxLogStatus( docMessage );
            } else {
                npix = 1;	// real for partial coherence mode
                docMessage.Printf(wxT("**ctem** range %f to %f"), rmin, rmax);
                wxLogStatus( docMessage );
            }
            makeBitmap( );
            UpdateAllViews();

            //   save the parameters actually used
            for( int i=0; i<Nparams; i++) sparam[i] = param[i];
            sparamSet = 1;
        }

        return;
};  // end RUNautoslic( )


//--------- RUNautostem -----------------
//  full STEM calculation
void ctDoc::RUNautostem( int lwobble, ctCanvas* canvas )
{
        int multiMode = 1;
        int idetect, it, ix, iy, nxout, nbeamp, nbeampo, nwobble;
        vectori collectorMode;  // for confocal mode

        bool bupdate;

        float  **rmin, **rmax, ***pixrline, rmin0, rmax0, rn;
        float res, thetamax;

        double keV, wavl, dx, dy, x0;
        vectord almin, almax, ThickSave, phiMin, phiMax;

        autostem ast;
        
        nwobble = ToInt( param[pNWOBBLE] );  // number phonon config.

        //   set autostem modes
        ast.l1d = 1;		// l1d;
        ast.lpacbed = 0;	// position averaged CBED
        ast.lxzimage = 0;	//  lxzimage;
        ast.lverbose = 0;	//  for debugging

        //  remember: nwobble must be 1 for non TDS and NOT 0 !!
        if( nwobble > 1 ) ast.lwobble = 1;
        else  ast.lwobble = 0;

        if( docType != xyzDoc ) {
                wxMessageBox( wxT("can only calculate from xyz files") );
                return;

        } else if( (nx<64) && (ny<64) ){
                wxMessageBox( wxT("nx,ny must be larger for autostem") );
                return;

        } else if( (nxiSTEM<1) && (nyiSTEM<1) ){
                wxMessageBox( wxT("nxi,nyi must be larger for autostem") );
                return;

        } else {
            expandCell();  //  make supercell from one unit cell

            wxLogStatus( ws.Format(wxT("autostem with supercell of %d atoms"), natoms)  );

            keV = param[pENERGY];
            wavl = wavelength( keV );

            //   set calculation parameters (some already set in dialogs)
            param[pAX] = axs;			// supercell size
            param[pBY] = bys;
            param[ pNX ] = (float) nx;
            param[ pNY ] = (float) ny;
            param[ pMODE ] = mAUTOSTEM;       // mode = autostem

            param[ pWAVEL ] = wavl;		//  probably recal. autoslice::calculate()

            param[ pNXPRB ] = (float) nxProbe;	// probe size in pixels
            param[ pNYPRB ] = (float) nyProbe;

            param[ pNXOUT ] = nxiSTEM;  //  save output sizel
            param[ pNYOUT ] = nyiSTEM;

            param[ pIMIN ] = 0.0F;	//  no imaginary image here
            param[ pIMAX ] = 0.0F;
            npix = 1;

            ast.CountBeams( param, nbeamp, nbeampo, res, thetamax );

            //  transmission function sampling
            wxLogStatus( ws.Format(
                wxT("Bandwidth limited to a real space resolution of %f Angstroms"), res) );
            wxLogStatus( ws.Format(
                wxT("   (= %.2f mrad)  for symmetrical anti-aliasing."), thetamax*1000.0F) );
        
            //  probe sampling
            wxLogStatus( ws.Format(
                    wxT("Number of symmetrical anti-aliasing beams in probe = %d"), nbeamp) );
            wxLogStatus( ws.Format(
                    wxT("Number of beams in probe aperture = %d"), nbeampo) );

            //--  test for adequate sampling -----
            if( nbeamp < 200 ) { wxMessageBox( 
                    wxT("WARNING: the probe is under sampled, this is a bad idea...") );
                return;
            }
            if( nbeampo < 100 ) { wxMessageBox(
                    wxT("WARNING: the probe aperture is under sampled, this is a bad idea...") );
                return;
            }

            ndetect = 1;   //  ??? need to get this from somewhere
            nThick = 1;

            collectorMode.resize( ndetect );
            for(idetect=0; idetect<ndetect; idetect++) collectorMode[idetect] = ADF;

            phiMin.resize( ndetect );
            phiMax.resize( ndetect );
            almin.resize( ndetect );
            almax.resize( ndetect );
            almin[0] = param[ pMINDET ];
            almax[0] = param[ pMAXDET ];   // detector angles in radians

            /* init the min/max record of total integrated intensity */
            rmin = new2D<float>( nThick, ndetect, "rmin" );
            rmax = new2D<float>( nThick, ndetect, "rmax" );

            //   for now just save the full thickness (should add more sometime ????)
            ThickSave.resize( nThick);
            ThickSave[0] = czs;

            pix.resize( nxiSTEM, nyiSTEM );  //  put image here to make bitmap

            //   memory for full image
            if( (n1pixr !=  ndetect*nThick) || (n2pixr != nxiSTEM) || (n3pixr != nyiSTEM )) {
                if( n1pixr > 0 ) delete3D<float>( pixr, n1pixr, n2pixr  );
                // double up first index to mimic a 4D array
                pixr = new3D<float>( ndetect*nThick, nxiSTEM, nyiSTEM, "pixr"  );
                n1pixr = ndetect*nThick;
                n2pixr = nxiSTEM;
                n3pixr = nyiSTEM;
            }

            //  memory for one line of image
            pixrline = new3D<float>(ndetect*nThick, 1, nyiSTEM, "pixrline"  );

            if( nxiSTEM > 1 ) dx = (xfSTEM - xiSTEM)/ (nxiSTEM-1);
            else dx = 1.0;
            param[pDX] = (float) dx;
            if( nyiSTEM > 1 ) dy = (yfSTEM - yiSTEM)/ (nyiSTEM-1);
            else dy = 1.0;
            param[pDY] = (float) dy;
            
            wxProgressDialog dlg("autostem progress",
                    wxT("calculating..."),
                    nxiSTEM, canvas,
                    //wxPD_CAN_ABORT | wxPD_APP_MODAL | wxPD_ELAPSED_TIME |
                    wxPD_CAN_ABORT | wxPD_ELAPSED_TIME | wxPD_AUTO_HIDE |
                    wxPD_ESTIMATED_TIME | wxPD_REMAINING_TIME );

            //  do one line at a time and paint each line on screen in real time
            //    it seems that the progrees bar blocks screen painting so doesn't quit work yet
            for( ix=0; ix<nxiSTEM; ix++) {

                x0 = xiSTEM + ix*dx;
                nxout = 1;

                //  update progress and check for cancel
                bupdate = dlg.Update( ix, 
                        ws.Format(wxT("calculating line %d"), ix) );
                if( !bupdate ) break;

                //  do one line of the autostem calculation
                ast.calculate( param, multiMode, natoms, &iseed, 
                    Znums, xs,ys,zs, occs, wobbles,
                    x0, x0, yiSTEM, yfSTEM, nxout, nyiSTEM,
                    ThickSave, nThick,
                    almin, almax, collectorMode, ndetect,
                    phiMin, phiMax,
                    pixrline, rmin, rmax, pacbedPix );

                //  save this line in the main buffer
                for( it=0; it<nThick; it++){
                    for( idetect=0; idetect<ndetect; idetect++) 
                        for( iy=0;iy<nyiSTEM; iy++) pixr[idetect+it*ndetect][ix][iy] 
                                = pixrline[idetect+it*ndetect][0][iy];
                }

                if( 0 == ix ) {
                    rmin0 = pixrline[0][0][0];  //pixr[0][ix][0];
                    rmax0 = rmin0;
                    pix = rmin0;
                }
                for( iy=0; iy<nyiSTEM; iy++) {
                    rn = pixrline[0][0][iy];  //pixr[0][ix][iy];
                    if( rn < rmin0 ) rmin0 = rn;
                    if( rn > rmax0 ) rmax0 = rn;
                    pix.re(ix,iy) = pixr[0][ix][iy];
                }
                param[ pRMIN ] = rmin0;
                param[ pRMAX ] = rmax0;

                docMessage.Printf( ws.Format(wxT("**autostem** range %f to %f"), rmin0, rmax0)  );
                wxLogStatus( docMessage );

            }  //  end for( ix; .... )

            //  these don't work inside the main loop for some reason
            //     would be nice to have a live view (?)
            makeBitmap( );
            UpdateAllViews();
			//UpdateDisplay();  //  not defined here ????
			//Update();   //  not defined here????

            //   save the parameters actually used
            for( int i=0; i<Nparams; i++) sparam[i] = param[i];
            sparamSet = 1;

			//  remember; vector<> will free its memory so we don't have to
 
            delete3D<float>( pixrline, ndetect*nThick, 1 );
            delete2D<float>( rmin, nThick );
            delete2D<float>( rmax, nThick );
        }

        return;
};  // end RUNautostem( )

//--------- RUNcbed -----------------
//  lcbed = 1 for CBED
//  lcbed = 1 for parallel illum. (normal diffraction)
void ctDoc::RUNcbed( ctCanvas* canvas, int lcbed )
{
	//--- autoslic parameters 
	int nbout = 0;  // number of beams to save
    int nwobble=0;

	//--- main variables
    int  ismoth;
    vectorf kx, ky, xpos, ypos, kx2, ky2;
    float rmin, rmax, aimin, aimax;

    double k2max, keV, wavl, rx, ry,
        rx2, ry2, xp, yp, pixel, apert;

    int multiMode = 1;   // enable multipole aberrations

    probe prb;
    autoslic aslice;

    if( docType != xyzDoc ) {
        wxMessageBox( wxT("can only calculate from xyz files") );
        return;

    } else if( (nx<64) && (ny<64) ){
        wxMessageBox( wxT("nx,ny must be larger for cbed") );
        return;
	} else {

        expandCell();  //  make supercell from one unit cell

		if( 1 == lcbed ) {  // CBED
			wxLogStatus( ws.Format(wxT("cbed with supercell of %d atoms"), natoms)  );

		} else if( 0 == lcbed ) {  // normal parallel diffraction
			wxLogStatus( ws.Format(wxT("diffract with supercell of %d atoms"), natoms)  );
		}

        keV = param[pENERGY];     // electron energy in keV
        wavl = wavelength( keV ); // electron wavelength in Ang.
		apert = param[pOAPERT];   // obj. aperture in radians

        nwobble = ToInt( param[pNWOBBLE] );

        /* ---------  setup calculation ----- */
        //   set autoslic calculation flags
        aslice.lbeams = 0;
        aslice.lcross = 0;
        aslice.lpartl = 0;
        aslice.lstart = 1;
        if( nwobble > 0 ) aslice.lwobble = 1;
        else  aslice.lwobble = 0;

		if( 0 == lcbed ) {  // normal parallel diffraction
			aslice.lcbed = 0;
		} else  {  // CBED
			aslice.lcbed = 1;
		}

        //   set calculation parameters (some already set in dialogs)
        param[pAX] = axs;		// supercell size
        param[pBY] = bys;
        param[pDX] = (float) ( axs/((float)nx) );
        param[pDY] = (float) ( bys/((float)ny) );
        param[ pNX ] = (float) nx;
        param[ pNY ] = (float) ny;
        param[ pMODE ] = mCBED;  // mode = cbed 

        param[ pWAVEL ] = wavl;		//  probably recal. in autoslice::calculate()

		k2max = apert/wavl;
		k2max = k2max * k2max;

		xp = param[ pPPOSX ];	//   probe position in Angst.
		yp = param[ pPPOSY ];

		if( 0 != lcbed ) {
			docMessage.Printf(wxT("cbed probe at %f, %f Ang."), xp, yp);
			wxLogStatus( docMessage );
		}

		//   probe needs these to be calculated separately
		kx.resize( nx );
		kx2.resize( nx );
		xpos.resize( nx );
		freqn( kx, kx2, xpos, nx, axs );

		ky.resize( ny );
		ky2.resize( ny );
		ypos.resize( ny );
		freqn( ky, ky2, ypos, ny, bys );
   
		rx  = 1.0/axs;
		rx2 = rx * rx;
		ry  = 1.0/bys;
		ry2 = ry * ry;
		pixel = ( rx2 + ry2 );
		ismoth = 0;

        // not used here but need to get these from somewhere ????
        double ycross = 0.0;
        double dfdelt = 2.0;   // defocus integration step size

		pix.resize( nx, ny );
		pix.init( 0 );

		aslice.calculateCBED_TDS( pix, param, multiMode, natoms, &iseed,
                Znums, xs,ys,zs,occs,wobbles );

        pix.findRange( rmin, rmax, aimin, aimax );
        param[ pRMIN ] = rmin;
        param[ pRMAX ] = rmax;
        param[ pIMIN ] = 0;  // remember to output only the real part
        param[ pIMAX ] = 0;

		if( (axs>0.0) && (bys>0.0) ){  // avoid div by 0 just in case
			param[pDX] = (float) ( 1.0/axs );  //  for FFT space !!
			param[pDY] = (float) ( 1.0/bys );
		}
        
		if( 0 != lcbed ) {
			docMessage.Printf(wxT("**cbed** range %f to %f"), rmin, rmax);
		} else {
			docMessage.Printf(wxT("**diffraction** range %f to %f "), 
						rmin, rmax ) ;
		}
		wxLogStatus( docMessage );

        npix = 1;

        makeBitmap( );
        UpdateAllViews();

        //   save the parameters actually used
        for( int i=0; i<Nparams; i++) sparam[i] = param[i];
        sparamSet = 1;

		//---  exit
	}

    return;

}  //  end RUNcbed(  )

//--------- RUNincostem -----------------
void ctDoc::RUNincostem( )
{
		int np;
        float rmin, rmax, aimin, aimax;

        double keV, wavl, probeI, dwellTime, scalep;
        int multiMode;

        incostem inc;

        if( docType != xyzDoc ) wxMessageBox( wxT("can only calculate from xyz files") );
        else if( (nx<64) && (ny<64) ){
                wxMessageBox( wxT("nx,ny must be larger for incostem") );
                return;
        } else if( (param[pMINDET] > param[pMAXDET]) || (param[pMINDET]<param[pOAPERT]) ) {
                wxMessageBox( wxT("bad detector angles for incostem, cannot calculate") );
                return;
        } else {

            expandCell();  //  make supercell from one unit cell
            wxLogStatus( ws.Format(wxT("incostem with supercell of %d atoms"), natoms)  );

            // ---  main calculation ---------------

            keV = param[pENERGY];
            wavl = wavelength( keV );

            multiMode = 1;

            param[pWAVEL] = (float) wavl;
            param[pAX] = axs;
            param[pBY] = bys;
            param[pDX] = (float) ( axs/((float)nx) );
            param[pDY] = (float) ( bys/((float)ny) );
            param[pNX] = (float) nx;
            param[pNY] = (float) ny;
            param[ pMODE ] = mINCOSTEM;  // mode = incostem

			probeI = param[ pPROBEI ] / 1.0e-12;   // probe current in pA
			dwellTime = param[ pPROBEDT ] / 1.0e-6; // probe dwell time in microSec

            // wxMessageBox( ws.Format(wxT("kev = %f"), param[pENERGY] ) );  // testing

            //  perform incostem on supercell
            inc.calculate2D( pix,  param, multiMode, natoms, Znums, xs, ys, occs ); 

            // wxMessageBox( ws.Format(wxT("incostem nx,ny = %d, %d"), nx, ny ) );  // testing

			// remember; 1 pAmp = 6.24146 MHz of electrons
			// incident current in electrons
			scalep = fabs(probeI * dwellTime * 6.24146);  // number of electrons/position
			if( (scalep>0.5) && (dwellTime>0.0) ) {  // scale to number of electrons
				wxLogStatus( "add electron counting noise" );
				np = inc.addNoise( pix, nx, ny, probeI, dwellTime );
			}

            pix.findRange( rmin, rmax, aimin, aimax );
            param[ pRMIN ] = rmin;
            param[ pRMAX ] = rmax;
            param[ pIMIN ] = 0;  // remember to output only the real part
            param[ pIMAX ] = 0;
            
            docMessage.Printf(wxT("**incostem** range %f to %f"), rmin, rmax);
            wxLogStatus( docMessage );

            npix = 1;
            makeBitmap( );
            UpdateAllViews();

            //   save the parameters actually used
            for( int i=0; i<Nparams; i++) sparam[i] = param[i];
            sparamSet = 1;
        }

        return;
}  // end RUNincostem( )

//--------- RUNmodel3d -----------------
void ctDoc::RUNmodel3d( )
{
        int i;
		float d, size0, scalepix;
        float rmin, rmax, aimin, aimax;
        double rotate, tilt;

		vectori icolor;
        vectorf xt, yt, zt, s;

        if( docType != xyzDoc ) {
                wxLogStatus( wxT("can only do model 3D of xyz coord") );
                return;
        }

        expandCell();  //  make supercell from one unit cell
        wxLogStatus( ws.Format(wxT("model 3D of %d atoms"), natoms)  );

        param[pDX] = 1.0F;   // drawing size NOT supercell size;
        param[pDY] = 1.0F;
        param[pNX] = (float) nx;
        param[pNY] = (float) ny;
        param[ pMODE ] = mMODEL3D;  // mode = model3D

        pix.resize( nx, ny );

        // make scratch arrays to work with and fill in extra required data
        //   will be modified by model3d()
        xt.resize( natoms );
        yt.resize( natoms );
        zt.resize( natoms );
        s.resize( natoms);
        icolor.resize( natoms );
        for( i=0; i<natoms; i++) {
                xt[i] = xs[i];
                yt[i] = ys[i];
                zt[i] = zs[i];
                s[i] = 1.0F;
                icolor[i] = 1;
        }

        rotate = 0.0;		// view head on
        tilt = 0.0;
        d = 5.0* fabs( czs );	// viewing distance (guess)
        if( d < 0.1 ) d = 1.0;		// for monolayer specimens like graphene

        // sphere size range 0.0->1.0 (guess)
        if( (sphere3d < 0.0F) || (sphere3d > 1.0F) ) size0 =  0.12F;
        else size0 = sphere3d;

        //  generate a 3D perspective view in pix
        model3( xt,yt,zt, axs, bys, czs, s, icolor, natoms, d, size0,
                rotate, tilt, pix, &scalepix );

		//  remember; vector<> will free its memory

        pix.findRange( rmin, rmax, aimin, aimax );
        param[ pRMIN ] = rmin;
        param[ pRMAX ] = rmax;
        param[ pIMIN ] = 0;  // remember to output only the real part
        param[ pIMAX ] = 0;

        npix = 1;	//  real part only
        makeBitmap( );
        UpdateAllViews();

        //   save the parameters actually used
        for( int i=0; i<Nparams; i++) sparam[i] = param[i];
        sparamSet = 1;

        wxLogStatus( wxT("***model 3D finished")  );

        return;

}  // end RUNmodel3d( )

// ---------------------------------------------------------------------------
/*--------------------- tinter() ----------------------------*/
/*
  Bilinear interpolation from data array

  nx,ny    = integer dimension of the image
  x,y      = real coord of image point to interpolate at
        0->(nx-1) and 0->(ny-1)

  ir = 0 for real part and 1 for imag.
*/
float ctDoc::tinter( cfpix &px, double x, double y, int ir )
{
    int ix, iy, ix2, iy2;
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

    if( 0 == ir ) {
        d = px.re(ix,iy)  - px.re(ix,iy2) - px.re(ix2,iy)
                 + px.re(ix2,iy2);
        c = px.re(ix,iy2) - px.re(ix,iy)  - d*x1;
        b = px.re(ix2,iy) - px.re(ix,iy)  - d*y1;
        a = px.re(ix,iy)  - b*x1 - c*y1 - d*x1*y1;
    } else if( 1 == ir ) {
        d = px.im(ix,iy)  - px.im(ix,iy2) - px.im(ix2,iy)
                + px.im(ix2,iy2);
        c = px.im(ix,iy2) - px.im(ix,iy)  - d*x1;
        b = px.im(ix2,iy) - px.im(ix,iy)  - d*y1;
        a = px.im(ix,iy)  - b*x1 - c*y1 - d*x1*y1;
    }

    ans = (float) ( a + b*x + c*y + d*x*y );
    
    return( ans );

}  /* end tinter() */


//-----------------------  below from sliceview ---------------------------------
//-----------------------  reorganized for computem gui -------------------------

/*-------------------- model3() -----------------*/
/*
  Subroutine to draw 3D perspective of molecules
  into a 2D image with hidden surfaces.
  Each atom will be represented by a shaded solid sphere.

  The CRT is in xy plane with rotation and tilt
  instead of azimuthal and polar angles.

  NOTE: this assumes that nx = ny which may not be the case!! ???

  NOTE: Color is not currently supported. However,
    leave the color code in for compatibility with
    old data files (from when there was color) and
    in case I ever figure out how to do it here.

   rewritten in C 22-feb-1996 E. Kirkland
   added EPS mode (kind of a kluge) 1-july-1997 ejk
   remove EPS mode and write into a pix array instead of TIFF
       file 28-nov-2013 ejk
   change coord. range to be the super cell cube not
      range of coord.  30-nov-2013 ejk
   switch to faster sort (better stride in shell sort) 1-dec-2013 ejk
   switch to sortByZ() from slicelib 2-dec-2013 ejk
   fix scaling with nx != ny 7-dec-2013 ejk

  This subroutine calls sphere()

* x[npts], y[npts], y[npts] = real valued x,y,z coordinates
                                of each atom
  ax, by, cz   = super cell size (range of atomic coord.)
* s[npts]      = real size of  each atom
* icolor[npts] = integer color code for each atom
                        (not currently used)
  npts  = integer number of atoms (or points)
  d     = real viewing distance
  size0 = real initial size of each sphere
          (relative to full size of the 2D image)
  rotat = real rotation angle (in radians)
  tilt  = real tilt angle (in radians)
* pix   = complex float array to get 2d perspective image
                (imag. part not used)
* scalepix       = real value to get the lateral scale of pix[][]

  * = these arguments are changed

  NOTE 2: on exit X,Y,Z,S,ICOLOR will have the scaled and sorted
      coordinates (the original values are overwritten)
*/
void ctDoc::model3( vectorf &x, vectorf &y, vectorf &z, 
        float ax, float by, float cz, vectorf &s, vectori &icolor, int npts,
        float d, float size0, double rotat, double tilt, 
        cfpix &pix, float *scalepix )
{
        int  ix, iy, i, nxout, nyout;
        float xmin, xmax, ymin, ymax,zmin,zmax, smax,
                cr, sr, st, ct, x0,y0,z0, scale, fs, fsx, fsy;
        double  size, pi, size2;

        std::string sbuffer;

        //  initialize graphics image buffer
        nxout = pix.nx();
        nyout = pix.ny();
        for( ix=0; ix<nxout; ix++)
            for( iy=0; iy<nyout; iy++){
                pix.re(ix,iy) = pix.im(ix,iy) = 0;
        }

        //  Calculate misc constants 
        cr = (float) cos( rotat );
        sr = (float) sin( rotat );
        ct = (float) cos( tilt );
        st = (float) sin( tilt );
        pi = 4.0 * atan( 1.0 );

        //  find range of atomic coordinates  */
        smax= -1.0e30F;
        for( i=0; i<npts; i++) {
           if( s[i] > smax ) smax = s[i];
        }
        xmin = 0;
        xmax = ax;
        ymin = 0;
        ymax = by;
        zmin = 0;
        zmax = cz;

        //  Move to center of molecule and rotate  */
        xmin = 0.5F*( xmax + xmin );
        ymin = 0.5F*( ymax + ymin );
        zmin = 0.5F*( zmax + zmin );

        /*  Translate  
                - don't translate in xy version
                rotate about the axis
        */
        for( i=0; i<npts; i++) {
                x[i] = x[i] - xmin;
                y[i] = y[i] - ymin;
                z[i] = z[i] - zmin;

                x0 = x[i];      /*  Rotation  */
                y0 = y[i];
                x[i] =  cr*x0 - sr*y0;
                y[i] =  sr*x0 + cr*y0;

                y0 = y[i];      /*  Tilt  */
                z0 = z[i];
                y[i] =  ct*y0 + st*z0;
                z[i] = -st*y0 + ct*z0;
        }

        //  Sort by depth  ( Shell sort )
        //  need to sort similar arrays as autoslic so use existing sort subroutine
        //    s[], icolor[], npts ->  occ[], Znum[], natoms
        sortByZ( x,y,z, s, icolor, npts );

        /*  Test sort routine -- DELETE this after awhile  */
        for( i=1; i<npts; i++) 
            if( z[i-1] > z[i] ) {
                sbuffer = "Bad sort in model3d !";
                messageSL( sbuffer.c_str() );	// use slicelib message handler for now
            }

        //  Form 2D perspective projection and find scale of x,y 

        xmin = ymin = 1.0e30F;
        xmax = ymax = -xmin;

        for( i=0; i<npts; i++) {
                x[i] = x[i] /( d-z[i] );
                y[i] = y[i] /( d-z[i] );
                if( x[i] < xmin ) xmin = x[i];
                if( x[i] > xmax ) xmax = x[i];
                if( y[i] < ymin ) ymin = y[i];
                if( y[i] > ymax ) ymax = y[i];
        }
        scale = xmax - xmin;
        if( (ymax-ymin) > scale ) scale = ymax-ymin;
        scale = 1.10F * scale;
        scale = ( 1.0F - size0 ) / scale;

        /* Scale coord. to ( 0.0 - 1.0 ) and display molecule
           SIZE0 must be in the 0.0-1.0 range */

        size = size0*( d - z[npts-1] );

        if( nxout < nyout ) fs = (float) nxout; 
                        else fs = (float) nyout;
        fsx = ((float)nxout) / fs;
        fsy = ((float)nyout) / fs;
        x0 = 0.5F* ( fsx - (xmax-xmin)* scale );   // the center (0.0 to 1.0)
        y0 = 0.5F* ( fsy - (ymax-ymin)* scale );
        for( i=0; i<npts; i++) {
                x[i] = ( x[i] - xmin ) * scale  + x0;
                y[i] = ( y[i] - ymin ) * scale  + y0;
                size2 = (s[i]/smax) * size / ( d - z[i] );
                sphere( x[i], y[i], size2, icolor[i], 
                        pix, nxout, nyout );
         }

        //  Save the scale 
        if( (xmax-xmin) > (ymax-ymin) ) scale = (xmax-xmin);
                else scale = (ymax-ymin);
        *scalepix = scale;

        return;

}  // end model3d()


/*------------------------ sphere() ------------------------*/
/*
   Generate shaded sphere at IX0,IY0 with
   radius IRAD  (in units of pixels)
   modified for cfpix images 30-nov-2013 ejk
   fix scaling with nxout != nyout 7-dec-2013 ejk

  x0, y0        = coordinates of sphere range= 0.0 to 1.0
  size          = diameter of sphere range 0.0 to 1.0
                  (this is a percentage of the full scale)
  icolor        =  0 red,  1 green,   2 blue
                        (not currently implemented!)
  pix           = float array to get image
                  (remember that it ix,iy indexes are reversed )
  nxout, nyout  = size of output image in pixels
*/
void ctDoc::sphere( double x0, double y0, double size, int icolor, 
                cfpix &pix, int nxout, int nyout )
{
        int ix, iy, ixd, iyd, ir2, ir, iy1, iy2,
                iyd2, ix1, ix2, ival1, ix0, iy0, irad, nm1;
        float val2, anout;

        //  scale coord. to max screen size of 0-NOUT 
        nm1 = nxout;
        if( nyout < nm1 ) nm1 = nyout;
        nm1 -= 1;
        anout = (double) nm1 ;
        ix0 = (int) (x0 * anout + 0.5);  // center in pixels
        iy0 = (int) (y0 * anout + 0.5);
        irad = (int) (0.5 * anout * size + 0.5);

        //  set up intensity scale 0-> 255
        ival1 = 255;
        val2  = 200.0F;

        //  Calculate sphere size
        iy1 = iy0 - irad;
        iy2 = iy0 + irad;

        if( (nyout-1) < iy1 ) iy1 = (nyout-1);
        if( 0 > iy1 ) iy1 = 0;

        if( (nyout-1) < iy2 ) iy2 = (nyout-1);
        if( 0 > iy2 ) iy2 = 0;

        ir2 = irad*irad;

        //  Generate a shaded circle that looks like a 3D sphere

        for( iy=iy1; iy<=iy2; iy++) {
           iyd  = iy - iy0;
           iyd2 = iyd*iyd;

           if( iyd2 < ir2 ) {
                ixd  = (int) (sqrt( ( (double)(ir2-iyd2) ) ) + 0.5);
                ix1  = ix0 - ixd;
                ix2  = ix0 + ixd;

                if( (nxout-1) < ix1 ) ix1 = (nxout-1);
                if( 0 > ix1 ) ix1 = 0;
                if( (nxout-1) < ix2 ) ix2 = (nxout-1);
                if( 0 > ix2 ) ix2 = 0;

                for( ix=ix1; ix<=ix2; ix++) {
                 ir = ix - ix0;
                 ir = iyd2 + ir*ir;
                  if( ir < ir2 )
                        pix.re(ix,iy) = ival1 - (val2 * ir) / ir2 ;
                }
           }
        }

        return;

}  // end sphere()



