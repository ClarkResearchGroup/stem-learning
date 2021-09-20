/*   ctview.cpp

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

computem view class

started 01-jan-2013 ejk
can show floating point images 24-feb-2013 ejk
add setup menu 3-mar-2013 ejk
start adding aber dialog 10-mar-2013 ejk
add onCalculate() 6-may-2013 ejk
add zoom function 28-may-2013 ejk
add exitWave mode 29-jun-2013 ejk
add OnActivateView() to update status bar 30-jun-2013 ejk
update status messages (move to doc) 20-oct-2013 ejk
changed new window name to new from Child 
   and add help/gettingStarted 10-nov-2013 ejk
add edit/copy = OnCopyPix() 14-nov-2013 ejk
add save pix as text 16-nov-2013 ejk
add model 3d mode 28-nov-2013 ejk
add load/save param menu+subroutine 8-feb-2014 ejk
add save xyz coord menu+subroutine 29-may-2014 ejk
add 2D abb. phase error calc. 31-aug-2014 ejk
add diffraction 9-jul-2017 ejk
change view/RealSpace/FFT to radioItem's 1-may-2018 ejk
add funny help message 29-sep-2019 - change 4-oct-2019 ejk
last modified 4-oct-2019 ejk
*/

//
#include "wx/wx.h"  // regular headers
//#include "wx/wxprec.h"  // for precompiled headers

#include "wx/mdi.h"     // MDI headers
#include "wx/docmdi.h"  // MDI headers for doc/view mode
#include "wx/docview.h"
#include "wx/image.h"
#include "wx/clipbrd.h"
#include <wx/filedlg.h>		// file dialog
#include "wx/string.h"

#include "computem.hpp"  //  header for this file
#include "ctdoc.hpp"     //  my document class
#include "ctview.hpp"    //  my view class

#include "slicelib.hpp"  // for RNG in onOpenDoor()

IMPLEMENT_DYNAMIC_CLASS( ctView, wxView)

BEGIN_EVENT_TABLE( ctView, wxView )
        EVT_MENU( ID_CALCULATE, ctView::OnCalculate )
        EVT_MENU( wxID_COPY, ctView::OnCopyPix )
        EVT_MENU( ID_SAVE_PIXASTEXT, ctView::OnSavePixText )
        EVT_MENU( ID_SAVE_PARAM, ctView::OnSaveParam )
        EVT_MENU( ID_LOAD_PARAM, ctView::OnLoadParam )
        EVT_MENU( ID_SAVE_COORD, ctView::OnSaveCoord )
        EVT_MENU( ID_ZOOM_NORMAL, ctView::OnZoomNormal )
        EVT_MENU( ID_ZOOM_OUT, ctView::OnZoomOut )
        EVT_MENU( ID_ZOOM_IN, ctView::OnZoomIn )
        EVT_MENU( ID_JUST_GUESS, ctView::OnGuess)
END_EVENT_TABLE()


// ---------------------------------------------------------------------------
//   creator
 ctView::ctView()
{
        canvas = (ctCanvas *) NULL;
        subFrame = (wxDocMDIChildFrame *) NULL;

        zoom = 1.0;  // normal zoom to start

        //  start with a different RNG seed each time this is run
		//   -only used for onOpenDoor
        long ltime = wxGetLocalTime();	// number of seconds since Jan 1, 1970
        iseed = (unsigned) ltime;	// init global seed
}

// ---------------------------------------------------------------------------
//
wxDocMDIChildFrame* ctView::CreateChildFrame(wxDocument *doc, wxView *view, bool isCanvas)
{
  //// Make a child frame
  wxDocMDIChildFrame *subframe =
      new wxDocMDIChildFrame(doc, view, wxGetApp().GetTopFrame(), wxID_ANY, wxT("new"),
                             wxPoint(10, 10), wxSize(300, 300),
                             wxDEFAULT_FRAME_STYLE |
                             wxNO_FULL_REPAINT_ON_RESIZE);

#ifdef __WXMSW__
  subframe->SetIcon(wxString(isCanvas ? wxT("chart") : wxT("notepad")));
#endif
#ifdef __X__
  subframe->SetIcon(wxIcon(_T("doc.xbm")));
#endif

  // ------- Make a file menu ------------
  wxMenu *file_menu = new wxMenu;

  file_menu->Append(wxID_NEW, wxT("&New..."));
  file_menu->Append(wxID_OPEN, wxT("&Open..."));
  file_menu->Append(wxID_CLOSE, wxT("&Close"));
  file_menu->Append(wxID_SAVE, wxT("&Save"));
  file_menu->Append(wxID_SAVEAS, wxT("Save &As..."));

  file_menu->AppendSeparator();
  file_menu->Append(ID_SAVE_PIXASTEXT, wxT("&Save pix as text"));
  file_menu->Append(ID_SAVE_PARAM, wxT("&Save param as..."));
  file_menu->Append(ID_LOAD_PARAM, wxT("&Load param from..."));
  file_menu->Append(ID_SAVE_COORD, wxT("&Save xyz coord. as..."));
 
  if (isCanvas)
  {
    file_menu->AppendSeparator();
    file_menu->Append(wxID_PRINT, wxT("&Print..."));
    file_menu->Append(wxID_PRINT_SETUP, wxT("Print &Setup..."));
    file_menu->Append(wxID_PREVIEW, wxT("Print Pre&view"));
  }

  file_menu->AppendSeparator();
  file_menu->Append(wxID_EXIT, _T("E&xit"));

  // ------- Make an edit menu ------------
  wxMenu *edit_menu = (wxMenu *) NULL;
  if (isCanvas)  {
    edit_menu = new wxMenu;
    edit_menu->Append( ID_EDIT_COORD, wxT("&atomic coord."));
    edit_menu->Append(wxID_COPY, wxT("&copy pix\tAlt-C"));
  }

  // ------- Make a help menu ------------
  wxMenu *help_menu = new wxMenu;
  help_menu->Append( ID_GET_START, wxT("&getting started\tF1"));
  help_menu->Append(wxID_ABOUT, wxT("&about\tF2"));
  help_menu->Append( ID_JUST_GUESS, wxT("&guess a number between 0 and 10"));

  // ------- Make a mode menu ------------
  calc_menu = new wxMenu;
  calc_menu->Append(ID_CALCULATE, wxT("&run\tAlt-R"));
  calc_menu->AppendSeparator();
  calc_menu->AppendRadioItem(ID_SETUP_AUTOCTEM, wxT("&ctem"));
  calc_menu->AppendRadioItem(ID_SETUP_CBED, wxT("&cbed"));
  calc_menu->AppendRadioItem(ID_SETUP_EXITWAVE, wxT("&exit wave"));
  calc_menu->AppendRadioItem(ID_SETUP_DIFFRACT, wxT("&diffraction"));
  calc_menu->AppendRadioItem(ID_SETUP_INCOSTEM, wxT("&incostem"));
  calc_menu->AppendRadioItem(ID_SETUP_AUTOSTEM, wxT("&stem"));
  calc_menu->AppendRadioItem(ID_SETUP_MODEL3D, wxT("&model 3d"));
  calc_menu->AppendRadioItem(ID_SETUP_ABBPHASE, wxT("&phase abb. func."));
 
  // ------- Make a setup menu ------------
  wxMenu *setup_menu = new wxMenu;
  setup_menu->Append(ID_SETUP_GENERAL, wxT("&general"));
  setup_menu->Append(ID_SETUP_CNM, wxT("&aberrations"));
  setup_menu->Append(ID_SETUP_SAMPLE, wxT("&Nx,Ny sampling"));
  
  // ------- Make a zoom menu ------------
  wxMenu *zoom_menu = new wxMenu;
  zoom_menu->Append(ID_ZOOM_NORMAL, wxT("&normal"));
  zoom_menu->Append(ID_ZOOM_OUT, wxT("&out\tAlt-+"));
  zoom_menu->Append(ID_ZOOM_IN, wxT("&in\tAlt--"));
  
  // ------- Make a greyscale menu ------------
  wxMenu *grey_menu = new wxMenu;
  grey_menu->Append(ID_GREY_FULL, wxT("&full range"));
  grey_menu->Append(ID_GREY_ONE, wxT("&one stnd dev"));
  grey_menu->Append(ID_GREY_SET, wxT("&set manually"));
  grey_menu->Append(ID_GREY_LOG, wxT("&log"));

  // ------- Make a view menu ------------
  wxMenu *view_menu = new wxMenu;
  view_menu->AppendRadioItem(ID_VIEW_REALSPACE, wxT("&real space"));
  view_menu->AppendRadioItem(ID_VIEW_FFT, wxT("&fft"));

  // ------- Make a menu bar ------------
  wxMenuBar *menu_bar = new wxMenuBar;

  menu_bar->Append(file_menu, wxT("&File"));
  if (isCanvas)
    menu_bar->Append(edit_menu, wxT("&Edit"));
  menu_bar->Append(calc_menu, wxT("&Calculate"));
  menu_bar->Append(setup_menu, wxT("&Setup"));
  menu_bar->Append(grey_menu, wxT("&Greyscale"));
  menu_bar->Append(zoom_menu, wxT("&Zoom"));
  menu_bar->Append(view_menu, wxT("&View"));
  menu_bar->Append(help_menu, wxT("&Help"));

  //// Associate the menu bar with the frame
  subframe->SetMenuBar(menu_bar);

  return subframe;

}  //  end ctView::CreateChildFrame()


// ---------------------------------------------------------------------------
//   update the status bar to explain active window
void ctView::OnActivateView( bool activate, wxView *activeView, wxView *deactivateView )
{
#ifdef NO_LOG_WINDOW	// should only do this if a log winodw has not been initiated

        //  get the image data from the doc
        ctDoc *ctd = (ctDoc*) GetDocument();

        wxLogStatus( ctd->docMessage );
#endif
};


// ---------------------------------------------------------------------------
bool ctView::OnCreate( wxDocument *doc, long WXUNUSED(flags) )
{
        int width=100, height=100;

        // Make a child frame
        //subFrame = wxGetApp().CreateChildFrame( doc, this, true );  // not necessary- can be here
        subFrame = CreateChildFrame( doc, this, true );

        //  make a scrolled window inside the frame to draw in
        //  remember can't draw in a frame
        canvas = new ctCanvas(this, subFrame, wxPoint(0, 0), wxSize(width, height), 0);
        canvas->SetCursor(wxCursor(wxCURSOR_PENCIL));
        canvas->SetScrollbars(20, 20, 50, 50);// Give it scrollbars

        //  connect to appropriate frame for this view
        SetFrame(subFrame);
        //???subFrame->GetCanvas()->SetView(this);

        subFrame->Show(true);

        // make sure the document knows this is the current view
        Activate(true);

        // init edit menu 
        //???doc->GetCommandProcessor()->SetEditMenu(subFrame->GetEditMenu());
        //???doc->GetCommandProcessor()->Initialize();

        nx = ny = nxb = nyb = npix = natom = 0;  //  remember there is no data yet

        count = 0;

        return true;

}   // end ctView::OnCreate()

// ---------------------------------------------------------------------------
// Clean up windows
bool ctView::OnClose(bool deleteWindow)
{
        if (!GetDocument()->Close())
                        return false;

        // Clear the canvas in case we're in single-window mode,
        // and the canvas stays.
        canvas->ClearBackground();
        canvas->view = (wxView *) NULL;
        canvas = (ctCanvas *) NULL;

        wxString s(wxTheApp->GetAppName());
        if (subFrame) subFrame->SetTitle(s);

        SetFrame((wxFrame*)NULL);

        Activate(false);

        if (deleteWindow)
        {
                delete subFrame;
                return true;
        }
        return true;

}  // end ctView::OnClose()

// ---------------------------------------------------------------------------
void ctView::OnUpdate(wxView* sender, wxObject* hint)
{
        wxView::OnUpdate(sender, hint);	//  call the default original first

        //  Refresh() is in wxScrolledWindow inherited from wxWindow()
        //  ask the canvas to update which will then call ctView::OnDraw with a DC
        canvas->Refresh();
        return;

} //  end ctView::OnUpdate()

// ---------------------------------------------------------------------------
//    copy the current bitmap image onto the clipboard for easy transfer to
//     other applications
void ctView::OnCopyPix( wxCommandEvent& event )
{
        ctDoc *ctd = (ctDoc*) GetDocument();
        nxb = ctd->nxb;
        nyb = ctd->nyb;
        if( !( (nxb>0) && (nyb>0) ) ) return;   // if image does not exist

	if( wxTheClipboard->Open() ) {
		wxTheClipboard->SetData( new wxBitmapDataObject( *ctd->ctBMP ) );
		wxTheClipboard->Close();
	}
        return;

} //  end ctView::OnCopyPix()


// ---------------------------------------------------------------------------
void ctView::OnDraw(wxDC *dc)
{
        int xsrc, ysrc, xdest, ydest, nxbz, nybz;

        //  get the image data from the doc
        ctDoc *ctd = (ctDoc*) GetDocument();
        nxb = ctd->nxb;
        nyb = ctd->nyb;
        nx = ctd->nx;
        ny = ctd->ny;
        npix = ctd->npix;

        natom = ctd->natom;
        natoms = ctd->natoms;

        if( (nxb>0) && (nyb>0) ) {  // if image exists

                //dc->DrawBitmap( *ctd->ctBMP, 0, 0, false ); // old no-zoom version
                //  wxDC::StretchBlit() to do zoom
                wxMemoryDC memDC(*ctd->ctBMP);
                xsrc = ysrc = 0;	// source coord.
                xdest = ydest = 0;	// destination coord.
                if( fabs(zoom-1.0) > 0.2 ) {
                        nxbz = int( nxb*zoom + 0.5 );
                        nybz = int( nyb*zoom + 0.5 );
                } else {
                        nxbz = nxb;
                        nybz = nyb;
                }
                canvas->SetVirtualSize( nxbz, nybz); //  size after interpolation+zoom
                dc->StretchBlit( xdest,ydest, nxbz,nybz, &memDC, xsrc,ysrc, nxb,nyb );
                //if( fabs(zoom-1.0) > 0.2 ) {  //  does NOT work ????
				//		// zoom about the center not a corner ???
                //        canvas->Scroll( nxbz/8, nybz/8); 
				//}

                //wxLogStatus( ctd->docMessage );  //  only need this if there is no LogWindow
        }

};   // end ctView::OnDraw()


//--------- OnCalculate -----------------
//  remember that this needs to be here to find setup_menu
//     and query its checked values
void ctView::OnCalculate( wxCommandEvent& event )
{
        //  get the doc
        ctDoc *ctd = (ctDoc*) GetDocument();

        if( calc_menu->IsChecked( ID_SETUP_INCOSTEM ) ) {
                ctd->RUNincostem();

        } else if( calc_menu->IsChecked( ID_SETUP_AUTOSTEM ) ) {
                ctd->RUNautostem( 0, canvas );

        } else if( calc_menu->IsChecked( ID_SETUP_AUTOCTEM ) ) {
                ctd->RUNautoslic( 1, 0 );

        } else if( calc_menu->IsChecked( ID_SETUP_EXITWAVE ) ) {
                ctd->RUNautoslic( 0, 0 );

        } else if( calc_menu->IsChecked( ID_SETUP_DIFFRACT ) ) {
                ctd->RUNcbed( canvas, 0 );

        } else if( calc_menu->IsChecked( ID_SETUP_MODEL3D ) ) {
                ctd->RUNmodel3d();

        } else if( calc_menu->IsChecked( ID_SETUP_ABBPHASE ) ) {
                ctd->RUNabbPhase();

        } else if( calc_menu->IsChecked( ID_SETUP_CBED ) ) {
                ctd->RUNcbed( canvas, 1 );
        }

        return;
};//  end ctView::OnCalculate()

//--------- OnGuess (just for fun) -----------------
void ctView::OnGuess( wxCommandEvent& event )
{
	int nmax=10;
	int i = (int) ( ((double)(nmax+1)) * ranflat( &iseed) ); 
	if( i < 0 ) i = 0;
	if( i > nmax ) i = nmax;
	wxString mes = "If you guessed "+toString(i)+", this might be your lucky day.";
    wxMessageBox( mes );

};  // end ctView::OnGuess

// ---------------------------------------------------------------------------
//    get filename to save text as
//     - must be here (not doc) because FileDialog needs window pointer
//        but xyz coord. are in the doc
void ctView::OnSaveCoord( wxCommandEvent& event )
{
	wxString caption, wildcard, defaultDir, defaultFilename;

    //  get the doc
    ctDoc *ctd = (ctDoc*) GetDocument();

	//   find where to save the xyz coord
	caption= wxT("select output text file to get xyz coord.");
	wildcard = wxT( "xyz files (*.xyz)|*.xyz");
	defaultDir = wxEmptyString;
	defaultFilename = wxEmptyString;

	wxFileDialog dlg( subFrame, caption, defaultDir, defaultFilename,
		wildcard, wxFD_SAVE|wxFD_OVERWRITE_PROMPT );
	if( dlg.ShowModal() == wxID_OK )
		ctd->DoSaveCoord( dlg.GetPath() );

        return;

} //  end ctView::OnSaveCoord()

// ---------------------------------------------------------------------------
//    get filename to save text as
//     - must be here (not doc) because FileDialog needs window pointer
void ctView::OnSavePixText( wxCommandEvent& event )
{
	int nxx, nyy, npixx;
	wxString caption, wildcard, defaultDir, defaultFilename;

        //  get the doc
        ctDoc *ctd = (ctDoc*) GetDocument();

        nxx = ctd->nxb;		//  bitmap size is only non-zero if something has been calculated
        nyy = ctd->nyb;
        npixx = ctd->npix;

	if( (nxx<1) || (nyy<1) || (npixx<1) || (npixx>2) ) {
		wxMessageBox( wxT("no valid image to save...") );
		return;
	}

	//   output real part (if real or complex)
	if( npixx == 1 ) caption= wxT("select output text file");
	else if( npixx == 2 ) caption= wxT("select real output text file");
	wildcard = wxT( "text files (*.txt)|*.txt");
	defaultDir = wxEmptyString;
	defaultFilename = wxEmptyString;

	wxFileDialog dlg( subFrame, caption, defaultDir, defaultFilename,
		wildcard, wxFD_SAVE|wxFD_OVERWRITE_PROMPT );
	if( dlg.ShowModal() == wxID_OK )
		ctd->DoSavePixText( dlg.GetPath(), 0 );

	//   output imaginary part if complex
	if( npixx == 2 ) {
		caption= wxT("select imaginary output text file");
		dlg.SetMessage( caption );
		if( dlg.ShowModal() == wxID_OK )
			ctd->DoSavePixText( dlg.GetPath(), 1 );
	}

        return;

} //  end ctView::OnSavePixText()

// ---------------------------------------------------------------------------
//    get filename to save text as
//     - must be here (not doc) because FileDialog needs window pointer
//        but param[] is in the doc
void ctView::OnSaveParam( wxCommandEvent& event )
{
	wxString caption, wildcard, defaultDir, defaultFilename;

        //  get the doc
        ctDoc *ctd = (ctDoc*) GetDocument();

	//   find where to save the parameters
	caption= wxT("select output text file to get param");
	wildcard = wxT( "text files (*.txt)|*.txt");
	defaultDir = wxEmptyString;
	defaultFilename = wxEmptyString;

	wxFileDialog dlg( subFrame, caption, defaultDir, defaultFilename,
		wildcard, wxFD_SAVE|wxFD_OVERWRITE_PROMPT );
	if( dlg.ShowModal() == wxID_OK )
		ctd->DoSaveParam( dlg.GetPath() );

        return;

} //  end ctView::OnSaveParam()


// ---------------------------------------------------------------------------
//    get filename to save text as
//     - must be here (not doc) because FileDialog needs window pointer
//        but param[] is in the doc
void ctView::OnLoadParam( wxCommandEvent& event )
{
	wxString caption, wildcard, defaultDir, defaultFilename;

        //  get the doc
        ctDoc *ctd = (ctDoc*) GetDocument();

 	//   find where to load the parameters
	caption= wxT("select input text file with param");
	wildcard = wxT( "text files (*.txt)|*.txt");
	defaultDir = wxEmptyString;
	defaultFilename = wxEmptyString;

	wxFileDialog dlg( subFrame, caption, defaultDir, defaultFilename,
		wildcard, wxFD_OPEN|wxFD_FILE_MUST_EXIST );
	if( dlg.ShowModal() == wxID_OK )
		ctd->DoLoadParam( dlg.GetPath() );

        return;

} //  end ctView::OnLoadParam()

//--------- OnZoomNormal -----------------
void ctView::OnZoomNormal( wxCommandEvent& event )
{
        zoom = 1.0;
        canvas->Refresh();
        return;
};//  end ctView::OnZoomNormal()


//--------- OnZoomOut -----------------
void ctView::OnZoomOut( wxCommandEvent& event )
{
        zoom *= 0.5;
        canvas->Refresh();
        return;
};//  end ctView::OnZoomOut()


//--------- OnZoomIn -----------------
void ctView::OnZoomIn( wxCommandEvent& event )
{
        zoom *= 2.0;
        canvas->Refresh();
        return;
};//  end ctView::OnZoomIn()


// ---------------------------------------------------------------------------
// Canvas
// ---------------------------------------------------------------------------
//   remember that you can't draw in a frame - need container to draw in
// Define a constructor for my canvas
ctCanvas::ctCanvas(wxView *v, wxDocMDIChildFrame *parent, const wxPoint& pos, const wxSize& size, long style):
 wxScrolledWindow(parent, wxID_ANY, pos, size, style)
{
  view = v;
}

// Define the repainting behaviour
void ctCanvas::OnDraw(wxDC& dc)
{
  if (view)
    view->OnDraw(& dc);
}
