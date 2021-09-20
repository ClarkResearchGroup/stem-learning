/*   computem.cpp

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

  temsim GUI using wxWidgets

  In wxWidgets first make an application class
  which makes a parent wxFrame window.
  The frame then makes a child frame.  You can't draw
  directly in the child frame but need a scrolled window in t
  child frame to draw in (much harder than it should be!)

-------------------------------------------------
  references:
  
  1. www.wxwidgets.org
  
  2. Julian Smart and Kevin Hock (with Stefan Csomor), "Cross-Platform GUI
  Programming with wxWidgets", Prentice Hall, 2006
     see chp. 19 for document/view widows style used here
     this is NOT complete - still need online help doc!!

  3. Syd Logan, "Cross-Platform Development in C++", Addison-Wesley 2008
 -------------------------------------------------
 Mac OSX 10.8
 
 in /build-cocoa-static run both of these
 
 ./wx-config --cxxflags
 ./wx-config --cxxflags
 
 then copy the settings produced into Xcode compiler and build options
 and add
 -I/opt/local/include  ; to compiler for fftw3.h
 -L/opt/local/lib  -lfftw3f -lfftw3f_threads  ; to build for fftw3 itsself
 
 also set compilre to ;  LLVC GCC 4.2
 
-------------------------------------------------
  MSVS 2010 settings (WXWIN=C:\wxWidgets-2.9.4 or C:\wxWidgets-3.0.0):

  1. precompiler:
     WIN32;__WXMSW__;NDEBUG;_WINDOWS;NOPCH

 2. additional include directories:
    "$(WXWIN)\lib\vc_lib\msw";"$(WXWIN)\include"  - version 2.8.12
    "$(WXWIN)\lib\vc_lib\mswu";"$(WXWIN)\include" - version 2.9.4
    "$(WXWIN)\lib\vc_x64_lib\mswu";"$(WXWIN)\include" - version 3.0.0, 3.1.0

 3. additional linker directories:
    $(WXWIN)\lib\vc_lib      - version 2.9.4
    $(WXWIN)\lib\vc_x64_lib  - version 3.0.0, 3.1.0

 4. link with:

   - version 2.9.4
   wxmsw29u_core.lib wxbase29u.lib wxtiff.lib wxjpeg.lib wxpng.lib
   wxzlib.lib wxregexu.lib wxexpat.lib winmm.lib comctl32.lib rpcrt4.lib 
   wsock32.lib odbc32.lib

   - version 3.0.0
   wxmsw30u_core.lib wxbase30u.lib wxtiff.lib wxjpeg.lib wxpng.lib
   wxzlib.lib wxregexu.lib wxexpat.lib winmm.lib comctl32.lib rpcrt4.lib 
   wsock32.lib odbc32.lib
   
   - version 3.1.0
   wxmsw31u_core.lib wxbase31u.lib wxtiff.lib wxjpeg.lib wxpng.lib
   wxzlib.lib wxregexu.lib wxexpat.lib winmm.lib comctl32.lib rpcrt4.lib 
   wsock32.lib odbc32.lib

 5. char set = "not set" (don't use unicode) or use the "u" versions of
 libraries (for example wxbase30u.lib); must build extra set
-------------------------------------------------

  started from wxpix.cpp 2-jul-2012 ejk
  add FileHistoryUseMenu() 18-apr-2013 ejk
  add New menu 20-apr-2013 ejk
  add wxLogWindow  1-sep-2013 ejk
  make smaller parent frame on mac osx because it doesn't work
      19-oct-2013 ejk
  add help/gettingStarted 9-nov-2013 ejk
  add "do not type here" top message in log window 11-nov-2013 ejk
  update "get started" message 1-dec-2013, 7-dec-2013 ejk
  update about/date 9-feb-2014 ejk
  update to wxWidgets 3.0.0 (almost same) 22-apr-2014 ejk
  update version 28-sep-2015 ejk
  update to wx 3.1.0  may 2016 ejk
  update about info 9-jul-2017 ejk
  update version 1-may-2018 ejk
  update version 17-aug-2019 ejk
  add funny help message 29-sep-2019 ejk
  last updated 29-sep-2019  ejk

*/

//
#include "wx/wx.h"  // regular headers
//#include "wx/wxprec.h"  // for precompiled headers

#include "wx/mdi.h"     // MDI headers
#include "wx/docmdi.h"  // MDI headers for doc/view mode
#include "wx/config.h"  // to get file history
#include "wx/image.h"   // init wxImage for clipboard routines

#include "mondrian.xpm"  // frame icon (from wx demo)

#include "computem.hpp"  //  header for this file
#include "ctdoc.hpp"     //  my document class
#include "ctview.hpp"    //  my view class

//  define to display file use history
#define USE_FILE_HISTORY

// ---------------------------------------------------------------------------
// MyApp
// ---------------------------------------------------------------------------

IMPLEMENT_APP( MyApp ) // Give wxWidgets the means to create a MyApp object

MyApp::MyApp()
{
        m_docManager = NULL;
}

int MyApp::OnExit()
{
#ifdef USE_FILE_HISTORY
        //  save file open history till next run of this application
        //  still not right - this never forgets them!!!!
        m_docManager->FileHistorySave(*wxConfig::Get());
#endif
        delete m_docManager;
        return 0;
}

// --------- Initialize the application -----------------
bool MyApp::OnInit()
{
        int nxSize=800, nySize=600;   // default size

        // create a document manager
        m_docManager = new wxDocManager;

        // create a template relating drawing documents to their views
        (void) new wxDocTemplate( m_docManager, wxT("computem"), wxT("*.tif;*.xyz;*.TIF;*.XYZ"),
                wxT(""), wxT("tif;xyz;TIF;XYZ"), wxT("computem doc"), wxT("computem view"),
                CLASSINFO(ctDoc), CLASSINFO(ctView) );

        // register the drawing document type on Mac
#if defined( __WXMAC__ )  && wxOSX_USE_CARBON
        wxFileName::MacRegisterDefaultTypeAndCreator( wxT("tif"),
                'WXMB', 'WXMA' );
#endif

#if defined( __WXMAC__ ) 
        nxSize= 300;        // parent frames doesn't work on Mac so make it small
        nySize= 100;
#endif

        // max number of doc open at same time - will start to save and del after this
        //      set to a large number !!!
        m_docManager->SetMaxDocsOpen(1000);

        // Create the main application window
        topFrame = new MyFrame( m_docManager, NULL, wxID_ANY, wxT("computem main"),
                wxPoint(0,0), wxSize(nxSize,nySize),
                wxDEFAULT_FRAME_STYLE | wxHSCROLL | wxVSCROLL);
        topFrame->Centre(wxBOTH);
        topFrame->Show(true); 	// Show it
        SetTopWindow( topFrame );

        // for ubuntu to copy bitmap to clipboard
        //wxImage::AddHandler( new wxPNGHandler);  // for ubuntu to copy bitmap to clipboard
        //wxImage::InitStandardHandlers();  // does not work
        wxInitAllImageHandlers();    // runs but no image on ubuntu

        return true;		// Start the event loop
}


// ---------------------------------------------------------------------------
// MyFrame
// ---------------------------------------------------------------------------

IMPLEMENT_CLASS(MyFrame, wxDocMDIParentFrame)

//--------- Event table for MyFrame -----------------
BEGIN_EVENT_TABLE( MyFrame, wxDocMDIParentFrame)
//	    EVT_MENU( wxID_OPEN, MyFrame::OnOpen)
        EVT_MENU( wxID_ABOUT, MyFrame::OnAbout)
        EVT_MENU( ID_GET_START, MyFrame::OnGetStart)
        EVT_MENU( wxID_EXIT, MyFrame::OnQuit)
        EVT_SIZE(MyFrame::OnSize)
END_EVENT_TABLE()

//--------- MyFrame() constructor -----------------
MyFrame::MyFrame( wxDocManager *manager, wxFrame *parent, wxWindowID id,
                const wxString& title, const wxPoint& pos, const wxSize& size, long type)
        :wxDocMDIParentFrame(manager, parent, id, title, pos, size, type ) 
{
        SetIcon( wxIcon(mondrian_xpm)); 	// Set the frame icon
        wxMenu *fileMenu = new wxMenu;		// Create a menu bar

        // The "About" item should be in the help menu
        wxMenu *helpMenu = new wxMenu;
        helpMenu->Append( ID_GET_START, wxT("&getting started...\tF1"),
                wxT("Show getting started dialog"));
        helpMenu->Append( wxID_ABOUT, wxT("&about...\tF2"),
                wxT("Show about dialog"));
        
        fileMenu->Append(wxID_NEW, wxT("&New\tAlt-N"),
                wxT("open new window") );
        fileMenu->Append(wxID_OPEN, wxT("O&pen file\tAlt-O"),
                wxT("opens a file"));
        fileMenu->Append(wxID_EXIT, wxT("E&xit\tAlt-X"),
                wxT("Quit this program"));

        // Now append the freshly created menu to the menu bar
        wxMenuBar *menuBar = new wxMenuBar();
        menuBar->Append(fileMenu, wxT("&File"));
        menuBar->Append(helpMenu, wxT("&Help"));

        // ... and attach this menu bar to the frame
        SetMenuBar( menuBar );

        // remember recent files opened this session
        manager->FileHistoryUseMenu(fileMenu);

#ifdef USE_FILE_HISTORY
        // remember files open in next run of this application
        //   must save history in OnExit()
        manager->FileHistoryLoad(*wxConfig::Get());
#endif

        //  Create a status bar just for fun - does not work on Mac OSX
        //  - wxLogWindow works better
        //CreateStatusBar();
        //SetStatusText(wxT("Welcome to computem"));

        //  make a separate log window for messages
        //   remember: MacOSX does NOT show a status bar in main window
        myLogWind = new wxLogWindow( this, wxT("computem log output") );
        wxLogStatus( wxT("computem log output only: do not type here") );

	//---- wx should used 4999 to 5999 but can check this here if needed
	//wxString ws;   //  to format messages
        //wxLogStatus( ws.Format(wxT("wxID_LOWEST= %d, wxID_HIGHEST= %d"), 
        //        wxID_LOWEST, wxID_HIGHEST)  );  // for testing

        SetSizeHints( 100, 100 );

}
//--------- MyFrame() about -----------------
void MyFrame::OnAbout( wxCommandEvent& event )
{
        wxMessageBox( wxT("computem version 4-oct-2019\n")
                      wxT("AS-IS with NO WARRANTY OR GUARANTEE\n")
                      wxT("under the GNU public license\n")
                      wxT("Copyright (C) 2013-2019 Earl J. Kirkland\n"), 
                wxT("About computem"),
                wxOK | wxICON_INFORMATION, this );
}


//--------- MyFrame() getting started -----------------
void MyFrame::OnGetStart( wxCommandEvent& event )
{
        wxMessageBox( wxT("computem is based on the theory in 'Adv. \n")
                      wxT(" Computing in Electron Microscopy', Springer 2010.\n")
                      wxT("1. Open an .xyz file of atomic coord. or an exiting\n")
                      wxT("   floating point .tif file of previous results\n")
                      wxT("   (view only). Some built in specimens coord. are\n")
                      wxT("   under file/new.\n")
                      wxT("2. Set parameters under Setup.\n")
                      wxT("3. Under Calculate select a mode and run.\n")
                      wxT("4. New .xyz files may be prepared with a plain text\n")
                      wxT("   editor (8 bit char. not unicode)."),
                wxT("getting started..."),
                wxOK | wxICON_INFORMATION, this );
}

//--------- MyFrame() quit -----------------
void MyFrame::OnQuit( wxCommandEvent& event)
{
        Close();	// Destroy the frame
}
