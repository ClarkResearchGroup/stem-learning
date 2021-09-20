/*   ctview.hpp

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

started 1-jan-2013 ejk
can show floating point images 24-feb-2013 ejk
start adding aber dialog 10-mar-2013 ejk
move aber and sampling dialog into ctdoc 4-may-2013 ejk
add onCalculate() 6-may-2013 ejk
add OnActivateView() to update status bar 30-jun-2013 ejk
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
add funny help message 29-sep-2019 - change 3-oct-2019 ejk
last modified 3-oct-2019 ejk
*/

#ifndef CTVIEW_HPP_INCLUDED	// only include this file if its not already
#define CTVIEW_HPP_INCLUDED	// remember that this has been included

#include "computem.hpp"  //  my application
#include "ctdoc.hpp"     //  my document class

#include "slicelib.hpp"   // my multislice subroutine library

enum{
	ID_ZOOM_NORMAL = 2201,
	ID_ZOOM_OUT = 2202,
	ID_ZOOM_IN = 2203,
	ID_SAVE_ROWASTEXT = 2204,
	ID_SAVE_COLASTEXT = 2205,
	ID_SAVE_PIXASTEXT = 2206,
	ID_SAVE_PARAM = 2207,
	ID_LOAD_PARAM = 2208,
    ID_SAVE_COORD = 2209,
	ID_VIEW_FFT = 2210,
	ID_VIEW_REALSPACE = 2211,
    ID_JUST_GUESS = 2212
};


// ---------------------------------------------------------------------------
//--------- canvas class -----------------

//BEGIN_EVENT_TABLE(ctCanvas, wxScrolledWindow)
//    EVT_MOUSE_EVENTS(ctCanvas::OnMouseEvent)
//END_EVENT_TABLE()

class ctCanvas : public wxScrolledWindow
{
public:
    wxView *view;

    ctCanvas(wxView *v, wxDocMDIChildFrame *parent, const wxPoint& pos, const wxSize& size,
		long style);

    virtual void OnDraw(wxDC& dc);
    //void OnEvent(wxMouseEvent& event);
    //??? void OnPaint(  wxPaintEvent &event );


private:

    //DECLARE_EVENT_TABLE()
};

// ---------------------------------------------------------------------------
//--------- ctView class -----------------

class ctView : public wxView
{
public:
        ctView();  // creator

        virtual bool OnCreate( wxDocument *doc, long flags );
        virtual void OnDraw( wxDC *dc);
        virtual void OnUpdate( wxView *sender, wxObject *hint=NULL );
        virtual bool OnClose( bool deleteWindow=true);
        void OnCopyPix( wxCommandEvent& event );
        virtual void OnActivateView( bool activate, wxView *activeView, wxView *deactivateView );

        wxDocMDIChildFrame* CreateChildFrame(wxDocument *doc, wxView *view, bool isCanvas);

private:

        wxDocMDIChildFrame *subFrame;

        ctCanvas* canvas;

        wxBitmap *ctBMP;  // a device dependent bitmap that can be sent to screen
        wxImage ctIMG;   // wx device independent image (cannot be drawn to screen!!)
        double zoom;

        int nx, ny;	// floating point image size (in pixels)
        int nxb, nyb;	// bitmap image size
        int npix;	// number of pix; 1 for real and 2 for complex
        int natom, natoms;	// number of atoms from xyz file

        int count;
        wxString ws;   //  to format messages

        unsigned long iseed;   // for random number generator - only used in onGuess()

        void OnCalculate( wxCommandEvent& WXUNUSED(event) );
        void OnGuess( wxCommandEvent& WXUNUSED(event) );
        void OnZoomNormal( wxCommandEvent& WXUNUSED(event) );
        void OnZoomOut( wxCommandEvent& WXUNUSED(event) );
        void OnZoomIn( wxCommandEvent& WXUNUSED(event) );

        void OnSaveCoord( wxCommandEvent& WXUNUSED(event) );

        void OnSavePixText( wxCommandEvent& WXUNUSED(event) );
        void OnSaveParam( wxCommandEvent& WXUNUSED(event) );
        void OnLoadParam( wxCommandEvent& WXUNUSED(event) );

        wxMenu *calc_menu;	//  save to read checked status

        DECLARE_DYNAMIC_CLASS( ctView )
        DECLARE_EVENT_TABLE()
};

#endif
