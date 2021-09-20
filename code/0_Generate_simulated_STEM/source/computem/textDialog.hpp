/*        textDialog.hpp 

------------------------------------------------------------------------
Copyright 2014-2019 Earl J. Kirkland


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

general text dialog to go with computem (to edit xyz values)

-------------------------------------------------
  references:

  1. Julian Smart and Kevin Hock (with Stefan Csomor), "Cross-Platform GUI
  Programming with wxWidgets", Prentice Hall, 2006
     see chp. 9  AND appendix J = whole code example

started from generalDialog 19-may-2014 ejk
fix close_box (add to style) 4-jun-2016 ejk
last updated 4-jun-2016 ejk

*/

#ifndef TEXTDIALOG_H
#define TEXTDIALOG_H

//
#include "wx/wx.h"  // regular headers
//#include "wx/wxprec.h"  // for precompiled headers

//#include "wx/spinctrl.h"
//#include "wx/gbsizer.h"

#include <string>   // STD string class
#include <vector>   // STD vector class

//  remember that 4999 thru 5999 are reserved
enum{
        ID_TEXT_DIALOG=    3300,
        ID_TEXT_TEXT=      3301
};

class textDialog: public wxDialog
{
        DECLARE_CLASS( textDialog )
public:
    textDialog( );
    textDialog(wxWindow *parent,
        wxWindowID id = ID_TEXT_DIALOG,
        const wxString& caption = wxT("text"),
        const wxPoint& pos = wxDefaultPosition,
        const wxSize& size = wxDefaultSize,
        long style = wxDEFAULT_DIALOG_STYLE); //  all but resize
        //long style = wxCLOSE_BOX|wxCAPTION|wxRESIZE_BORDER|wxSYSTEM_MENU );

    //  Initialize our variable
    void Init();

    // Creation
    bool Create( wxWindow* parent,
        wxWindowID id = ID_TEXT_DIALOG,
        const wxString& caption = wxT("text"),
        const wxPoint& pos = wxDefaultPosition,
        const wxSize& size = wxDefaultSize,
        long style = wxDEFAULT_DIALOG_STYLE); //  short version
        //long style = wxCLOSE_BOX|wxCAPTION|wxRESIZE_BORDER|wxSYSTEM_MENU );

    //  create controls and sizers
    void CreateControls();

    bool TransferDataToWindow( );
    bool TransferDataFromWindow( );

    //void OnButton(wxCommandEvent& event);
    //void OnClose( wxCloseEvent& event );

    //----  transfer test back and forth one line per string ----
    std::vector<std::string> tx;

    int nlines;
    int modified;

private:
    //  keep these here to refer to in several places
    wxTextCtrl *TextArea;

	//DECLARE_EVENT_TABLE()

};


#endif