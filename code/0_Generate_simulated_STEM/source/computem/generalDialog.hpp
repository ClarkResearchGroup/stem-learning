/*        generalDialog.hpp 

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

general dialog to go with computem (display/edit Nx,Ny etc. values)

-------------------------------------------------
  references:

  1. Julian Smart and Kevin Hock (with Stefan Csomor), "Cross-Platform GUI
  Programming with wxWidgets", Prentice Hall, 2006
     see chp. 9  AND appendix J = whole code example

-------------------------------------------------

started from aberrDialog 7-mar-2013 ejk
add ncellx,y 12-may-2013 ejk
add ctiltx,y 25-sep-2013 ejk
add capertmin,max 29-sep-2013 ejk
add model 3D sphere size 3-dec-2013 ejk
add CBED probe position 14-sep-2014 ejk
add STEM probe current and dwell time 29-sep-2015 ejk
fix close_box (fix create style) 5-jun-2016 ejk
last updated 5-jun-2016 ejk
change pA -> pAmp and microS -> microSec 18-aug-2019 ejk

formatted for a TAB size of 4 char

*/

#ifndef GENERALDIALOG_H
#define GENERALDIALOG_H

//
#include "wx/wx.h"  // regular headers
//#include "wx/wxprec.h"  // for precompiled headers

//#include "wx/spinctrl.h"
//#include "wx/gbsizer.h"


//  remember that 4999 thru 5999 are reserved
enum{
        ID_GENERAL_DIALOG=  3500,
        ID_GENERAL_KEV=	    3501,
        ID_GENERAL_OBJ=	    3502,
        ID_GENERAL_CAPERT1=	3503,
        ID_GENERAL_CAPERT2=	3504,
        ID_GENERAL_SRC=		3505,
        ID_GENERAL_DMIN=	3506,
        ID_GENERAL_DMAX=	3507,
        ID_GENERAL_CTILTX=	3508,
        ID_GENERAL_CTILTY=	3509,
        ID_GENERAL_SPHER3D=	3510
};

class generalDialog: public wxDialog
{
        DECLARE_CLASS( generalDialog )
public:
    generalDialog( );
    generalDialog(wxWindow *parent,
        wxWindowID id = ID_GENERAL_DIALOG,
        const wxString& caption = wxT("general"),
        const wxPoint& pos = wxDefaultPosition,
        const wxSize& size = wxDefaultSize,
        long style = wxDEFAULT_DIALOG_STYLE); //  all but resize

    //  Initialize our variable
    void Init();

    // Creation
    bool Create( wxWindow* parent,
        wxWindowID id = ID_GENERAL_DIALOG,
        const wxString& caption = wxT("general"),
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
    double kev, objAp;			// beam energy and objective aperture
    double ctiltx, ctilty;		// crystal tilt in radians
    double capertmin, capertmax;	// condencer angles in radians

    double detmin, detmax;			//  ADF detector size
	double probePosx, probePosy;	//  CBED probe position
    double dsource;					//  ADF source size
    double sphere3d;				//  model 3D rel. sphere size
	double GreyLogConstant;			//  log greyscale scaling constant

	double probeI, probeDt;			//  STEM probe current and dwell time

    int nerror;  // error count

private:
        //  keep these here to refer to in several places
        wxTextCtrl *kevVal, *objVal, *capertminVal, *capertmaxVal, *sphere3dVal;

        wxTextCtrl *dminVal, *dmaxVal, *dsrcVal, *ctiltxVal, *ctiltyVal;

		wxTextCtrl *probePosxVal, *probePosyVal, *GreyLogVal;

		wxTextCtrl *prbIval, *prbDTval;

//	DECLARE_EVENT_TABLE()
};


#endif