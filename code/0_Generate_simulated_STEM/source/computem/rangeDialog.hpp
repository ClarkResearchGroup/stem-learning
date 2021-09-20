/*        rangeDialog.hpp 

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

range dialog to go with computem (manually scale the greyscale)

-------------------------------------------------
  references:

  1. Julian Smart and Kevin Hock (with Stefan Csomor), "Cross-Platform GUI
  Programming with wxWidgets", Prentice Hall, 2006
     see chp. 9  AND appendix J = whole code example

-------------------------------------------------

started from sampleDialog 11-sep-2014 ejk
fix close_box (add to style) 5-jun-2016 ejk
last updated 5-jun-2016 ejk

*/

#ifndef RANGEDIALOG_H
#define RANGEDIALOG_H

//
#include "wx/wx.h"  // regular headers
//#include "wx/wxprec.h"  // for precompiled headers

//#include "wx/spinctrl.h"
//#include "wx/gbsizer.h"


//  remember that 4999 thru 5999 are reserved
enum{
        ID_RANGE_DIALOG=	3400,
        ID_RANGE_RMIN=		3401,
        ID_RANGE_RMAX=		3402,
        ID_RANGE_AIMIN=		3403,
        ID_RANGE_AIMAX=	    3404
};

class rangeDialog: public wxDialog
{
        DECLARE_CLASS( rangeDialog )
public:
    rangeDialog( );
    rangeDialog(wxWindow *parent,
        wxWindowID id = ID_RANGE_DIALOG,
        const wxString& caption = wxT("pix range"),
        const wxPoint& pos = wxDefaultPosition,
        const wxSize& size = wxDefaultSize,
        long style = wxDEFAULT_DIALOG_STYLE); //  all but resize

    //  Initialize our variable
    void Init();

    // Creation
    bool Create( wxWindow* parent,
        wxWindowID id = ID_RANGE_DIALOG,
        const wxString& caption = wxT("pix range"),
        const wxPoint& pos = wxDefaultPosition,
        const wxSize& size = wxDefaultSize,
        long style = wxDEFAULT_DIALOG_STYLE); //  all but resize

    //  create controls and sizers
    void CreateControls();

    bool TransferDataToWindow( );
    bool TransferDataFromWindow( );

    //void OnButton(wxCommandEvent& event);
    //void OnClose(wxCloseEvent& event);

    //----  variables to read in ----
	//   these must be initialized before shpwing the dialog box
    double rmin, rmax, aimin, aimax;		//  range of greyscale
	int npix;   // 1 for real image and 2 for complex image

    int nerror;  // error count

private:
        //  keep these here to refer to in several places
        wxTextCtrl *rminVal, *rmaxVal, *aiminVal, *aimaxVal;

//	DECLARE_EVENT_TABLE()
};


#endif