/*        aberDialog.hpp 

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

aberration dialog to go with computem (display/edit Cnm values)

-------------------------------------------------
  references:

  1. Julian Smart and Kevin Hock (with Stefan Csomor), "Cross-Platform GUI
  Programming with wxWidgets", Prentice Hall, 2006
     see chp. 9  AND appendix J = whole code example

-------------------------------------------------

started 7-mar-2013 ejk
add ddf = defocus spread 20-may-2013 ejk
fix close_box (fix create style) 5-jun-2016 ejk
add option to add pi/4 random tuning errors 11-jun-2016 ejk
last updated 11-jun-2016 ejk

formatted for a TAB size of 4 char*/

#ifndef ABERRDIALOG_H
#define ABERRDIALOG_H

//
#include "wx/wx.h"  // regular headers
//#include "wx/wxprec.h"  // for precompiled headers

//#include "wx/spinctrl.h"
//#include "wx/gbsizer.h"


//  remember that 4999 thru 5999 are reserved
enum{
        ID_ABERR_DIALOG=3000,
        ID_ABERR_DF=	3001,
        ID_ABERR_DDF=	3002,
        ID_ABERR_C12A=	3003,
        ID_ABERR_C12B=	3004,
        ID_ABERR_C21A=	3005,
        ID_ABERR_C21B=	3006,
        ID_ABERR_C23A=	3007,
        ID_ABERR_C23B=	3008,
        ID_ABERR_C30=   3009,
        ID_ABERR_C32A=	3010,
        ID_ABERR_C32B=	3011,
        ID_ABERR_C34A=	3012,
        ID_ABERR_C34B=	3013,
        ID_ABERR_C41A=	3014,
        ID_ABERR_C41B=	3015,
        ID_ABERR_C43A=	3016,
        ID_ABERR_C43B=	3017,
        ID_ABERR_C45A=	3018,
        ID_ABERR_C45B=	3019,
        ID_ABERR_C50=   3020,
        ID_ABERR_C52A=	3021,
        ID_ABERR_C52B=	3022,
        ID_ABERR_C54A=	3023,
        ID_ABERR_C54B=	3024,
        ID_ABERR_C56A=	3025,
        ID_ABERR_C56B=	3026,
        ID_ABERR_ZERO=  3027,
        ID_ADD_ERRORS=  3028
};

class aberDialog: public wxDialog
{
        DECLARE_CLASS( aberDialog )
public:
    aberDialog( );
    aberDialog(wxWindow *parent,
        wxWindowID id = ID_ABERR_DIALOG,
        const wxString& caption = wxT("aberrations"),
        const wxPoint& pos = wxDefaultPosition,
        const wxSize& size = wxDefaultSize,
        long style = wxDEFAULT_DIALOG_STYLE );

    //  Initialize our variable
    void Init();

    // Creation
    bool Create( wxWindow* parent,
        wxWindowID id = ID_ABERR_DIALOG,
        const wxString& caption = wxT("aberrations"),
        const wxPoint& pos = wxDefaultPosition,
        const wxSize& size = wxDefaultSize,
        long style = wxDEFAULT_DIALOG_STYLE );

    //  create controls and sizers
    void CreateControls();

    bool TransferDataToWindow( );
    bool TransferDataFromWindow( );
    void OnZero( wxCommandEvent &event );
	void OnAddErrors( wxCommandEvent &event );

    //void OnButton(wxCommandEvent& event);
    //void OnClose(wxCloseEvent& event);

    //----  variables to read in ----
    //   all in Angstroms (do scaling here)
    double df, ddf, C12a, C12b, C21a, C21b, C23a, C23b;
    double C30, C32a, C32b, C34a, C34b;
    double C41a, C41b, C43a, C43b, C45a, C45b;
    double C50, C52a, C52b, C54a, C54b, C56a, C56b;

    int nerror;  // error count

	unsigned long *iseed;  //  to get reference to main iseed for RNG
	double objApert;	//  to get objective aperture in radians
	double wavl;		//  wavelength in Angstroms

private:
        //  keep these here to refer to in several places
        wxTextCtrl *dfValue, *ddfValue, *C12aVal, *C12bVal, *C21aVal, *C21bVal,
                *C23aVal, *C23bVal, *C30Val, *C32aVal, *C32bVal, *C34aVal, *C34bVal,
                *C41aVal, *C41bVal, *C43aVal, *C43bVal, *C45aVal, *C45bVal,
                *C50Val, *C52aVal, *C52bVal, *C54aVal, *C54bVal, *C56aVal, *C56bVal;

        //  parameter scaling 
        //  divide by these values to convert from Ang. to displayed units
        double scale1, scale2, scale3, scale4, scale5;

        DECLARE_EVENT_TABLE()
};


#endif