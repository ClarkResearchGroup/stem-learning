/*        textDialog.cpp 

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

#include "textDialog.hpp"

IMPLEMENT_CLASS( textDialog, wxDialog )

//-------- constructors -----------------------
textDialog::textDialog() { Init(); }

textDialog::textDialog(wxWindow *parent,
        wxWindowID id,
        const wxString& caption,
        const wxPoint& pos,
        const wxSize& size,
        long style )
{
        Init();
        Create( parent, id, caption, pos, size, style);
};

 //-------- Create() -----------------------
 bool textDialog::Create( wxWindow* parent,
        wxWindowID id,
        const wxString& caption,
        const wxPoint& pos,
        const wxSize& size,
        long style )
 {
        //  first the built-in wx Create() before ours
        if( !wxDialog::Create( parent, id, caption, pos, size, style))
               return false;
 
        CreateControls();
        //SetDialogHelp();
        //SetDialogValidators();

        GetSizer()->Fit(this);
        GetSizer()->SetSizeHints(this);

        Centre();

        return true;
 };


//-------- Init() -----------------------
 void textDialog::Init( )
 {
        return;
 };

//-------- CreateControls() -----------------------
 void textDialog::CreateControls( )
 {
        wxBoxSizer *topSizer, *boxSizer;
        wxString ws;

        //  ---- top-level sizer -----
        topSizer = new wxBoxSizer( wxVERTICAL );
        this->SetSizer( topSizer );

        // ----- second box sizer to give more space around controls
        boxSizer = new wxBoxSizer( wxVERTICAL );
        topSizer->Add( boxSizer, 0, wxALIGN_CENTER_HORIZONTAL|wxALL, 5);

        //========================================================

        // ------ big text edit area
        // ---------------------------------------------

        TextArea = new wxTextCtrl( this, ID_TEXT_TEXT, wxT(""),
                wxDefaultPosition, wxSize(700, 500), wxTE_MULTILINE);

        boxSizer->Add( TextArea, 0,  wxALIGN_LEFT|wxALL, 5);

        //========================================================
        //---------------- bottom row of action buttons ----------------------
        //---- horz sizer box for reset,OK,Cancel,Help
        wxBoxSizer* okCancelBox = new wxBoxSizer( wxHORIZONTAL );
        boxSizer->Add(okCancelBox, 0, wxALIGN_CENTER_HORIZONTAL|wxALL, 5);

        //---- OK button
        wxButton* ok = new wxButton( this, wxID_OK, wxT("&OK"),
                wxDefaultPosition, wxDefaultSize, 0);
        okCancelBox->Add(ok, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

        //---- Cancel button
        wxButton* cancel = new wxButton( this, wxID_CANCEL, wxT("&Cancel"),
                wxDefaultPosition, wxDefaultSize, 0);
        okCancelBox->Add(cancel, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);
 };

 //========================================================
//-------- TransferDataToWindow() -----------------------
 //   update values on screen when first started
 //    CreateControls() doesn't have right values for some reason
 //
 bool textDialog::TransferDataToWindow( )
 {
        //wxString ws;
        //kevVal->SetValue(ws.Format(wxT("%10.3f"), kev) );

    int i, n;
    n = tx.size();

    if( tx.size() > 0 ) {
        for( i=0; i<n; i++)
            *TextArea << tx[i] << "\n";
    }

    modified = 0;

    return true;
 };

 //-------- TransferDataFromWindow() -----------------------
 bool textDialog::TransferDataFromWindow( )
 {
    int i;
    wxString wxs;

    //if( TextArea->IsModified() ){
	//  IsModified() doesn't seem to work in OSX so always do this
	//   - have to synch on OK vs. Cancel only (ejk 31-may-2014)

    nlines = TextArea->GetNumberOfLines();
    modified = 1;
    tx.clear(); 
    for( i=0; i<nlines; i++) {
        wxs = TextArea->GetLineText(i);
        tx.push_back( wxs.ToStdString() );
    }

    //} else modified = 0;
 
    return true;
 }

 