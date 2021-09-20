/*        sampleDialog.hpp 

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

sampling dialog to go with computem (display/edit Nx,Ny etc. values)

-------------------------------------------------
  references:

  1. Julian Smart and Kevin Hock (with Stefan Csomor), "Cross-Platform GUI
  Programming with wxWidgets", Prentice Hall, 2006
     see chp. 9  AND appendix J = whole code example

-------------------------------------------------

started from aberrDialog 7-mar-2013 ejk
add ncellx,y 12-may-2013 ejk
add nwobble and temp. 28-sep-2013 ejk
fix close_box (add to style) 5-jun-2016 ejk
last updated 5-jun-2016 ejk

*/

#ifndef SAMPLEDIALOG_H
#define SAMPLEDIALOG_H

//
#include "wx/wx.h"  // regular headers
//#include "wx/wxprec.h"  // for precompiled headers

//#include "wx/spinctrl.h"
//#include "wx/gbsizer.h"


//  remember that 4999 thru 5999 are reserved
enum{
        ID_SAMPLE_DIALOG=	3200,
        ID_SAMPLE_NX=		3201,
        ID_SAMPLE_NY=		3202,
        ID_SAMPLE_NXP=		3203,
        ID_SAMPLE_NYP=	        3204,
        ID_SAMPLE_NCELLX=	3205,
        ID_SAMPLE_NCELLY=	3206,
        ID_SAMPLE_NCELLZ=	3207,
        ID_SAMPLE_DELTAZ=	3208,
        ID_SAMPLE_NWOBBLE=	3209,
        ID_SAMPLE_TEMPERAT=	3210,
        ID_SAMPLE_XINIT=	3211,
        ID_SAMPLE_XFINAL=	3212,
        ID_SAMPLE_YINIT=	3213,
        ID_SAMPLE_YFINAL=	3214,
        ID_SAMPLE_NXI=		3215,
        ID_SAMPLE_NYI=		3216
};

class sampleDialog: public wxDialog
{
        DECLARE_CLASS( sampleDialog )
public:
    sampleDialog( );
    sampleDialog(wxWindow *parent,
        wxWindowID id = ID_SAMPLE_DIALOG,
        const wxString& caption = wxT("sampling"),
        const wxPoint& pos = wxDefaultPosition,
        const wxSize& size = wxDefaultSize,
        long style = wxDEFAULT_DIALOG_STYLE);

    //  Initialize our variable
    void Init();

    // Creation
    bool Create( wxWindow* parent,
        wxWindowID id = ID_SAMPLE_DIALOG,
        const wxString& caption = wxT("aberrations"),
        const wxPoint& pos = wxDefaultPosition,
        const wxSize& size = wxDefaultSize,
        long style = wxDEFAULT_DIALOG_STYLE);

    //  create controls and sizers
    void CreateControls();

    bool TransferDataToWindow( );
    bool TransferDataFromWindow( );

    //void OnButton(wxCommandEvent& event);
    //void OnClose(wxCloseEvent& event);

    //----  variables to read in ----
    int nx, ny, nxProbe, nyProbe;	//  transmission function and prbe size
    int ncellx, ncelly, ncellz;		//  super cell size (expand unit cell by this)
    int nxi, nyi;			//  image sixe in pixels
    int nwobble;			//  number of phonon config.
    double xi,xf, yi,yf;		//  initial and final x,y positions
    double temp;			//  temperature for phonon configurations

    double deltaz;			// slice thickness (in Ang.)

    int nerror;  // error count

private:
        //  keep these here to refer to in several places
        wxTextCtrl *nxVal, *nyVal, *nxpVal, *nypVal, *ncellxVal, *ncellyVal, *ncellzVal,
                *deltazVal, *nwobbleVal, *tempVal,
                *xiVal, *xfVal, *yiVal, *yfVal, *nxiVal, *nyiVal;


//	DECLARE_EVENT_TABLE()
};


#endif