/*        rangeDialog.cpp 

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

range dialog to go with computem (manually set the greyscale)

-------------------------------------------------
  references:

  1. Julian Smart and Kevin Hock (with Stefan Csomor), "Cross-Platform GUI
  Programming with wxWidgets", Prentice Hall, 2006
     see chp. 9  AND appendix J = whole code example

-------------------------------------------------

started from sampleDialog 11-sep-2014 ejk
fix some horz. spacing for wxWidgets 3.1.0/macosx/ubuntu  12-may-2016 ejk
fix close_box (add to style) 5-jun-2016 ejk
last updated 5-jun-2016 ejk

  formatted for a TAB size of 4 char

*/

#include "rangeDialog.hpp"

IMPLEMENT_CLASS( rangeDialog, wxDialog )

//BEGIN_EVENT_TABLE( aberDialog, wxDialog )
//	EVT_BUTTON( ID_ABERR_ZERO, aberDialog::OnZero )
//END_EVENT_TABLE()

//-------- constructors -----------------------
rangeDialog::rangeDialog() { Init(); }

rangeDialog::rangeDialog(wxWindow *parent,
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
 bool rangeDialog::Create( wxWindow* parent,
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
 void rangeDialog::Init( )
 {
        return;
 };

 //-------- CreateControls() -----------------------
 void rangeDialog::CreateControls( )
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

        // ------ start boxes
        // ---------------------------------------------

        // ----- range of real part ----------------------------------------
        //
        wxStaticBox* staticBox1 = new wxStaticBox( this, wxID_ANY, wxT("real pix range") );
        wxStaticBoxSizer* staticSizer1 = new wxStaticBoxSizer( staticBox1, wxHORIZONTAL );

        wxStaticText *rminName= new wxStaticText( this, wxID_STATIC, wxT("&min ="),
                wxDefaultPosition, wxDefaultSize, 0);
        rminVal = new wxTextCtrl( this, ID_RANGE_RMIN, wxT("0"),
                wxDefaultPosition, wxDefaultSize, 0);
        staticSizer1->Add( rminName, 0, wxALL, 5 );
        staticSizer1->Add( rminVal, 0, wxALL, 5 );

        wxStaticText* rmaxName= new wxStaticText( this, wxID_STATIC, wxT("&max ="),
                wxDefaultPosition, wxDefaultSize, 0);
        rmaxVal = new wxTextCtrl( this, ID_RANGE_RMAX, wxT("0"),
                wxDefaultPosition, wxDefaultSize, 0);
        staticSizer1->Add( rmaxName, 0, wxALL, 5 );
        staticSizer1->Add( rmaxVal, 0, wxALL, 5 );

        boxSizer->Add( staticSizer1, 0,  wxALIGN_LEFT|wxALL, 5);

		//  remember that it seems this is called at creation before parameters
		//  have been intialized so can't just check npix
		//    have to disable(grey out) in TransferDataToWindow() below
		//if( npix > 1 ) {   ///   does NOT work ?????

			// ----- range of imag part ----------------------------------------
			//
			wxStaticBox* staticBox2 = new wxStaticBox( this, wxID_ANY, wxT("imag. pix range") );
			wxStaticBoxSizer* staticSizer2 = new wxStaticBoxSizer( staticBox2, wxHORIZONTAL );

			wxStaticText *aiminName= new wxStaticText( this, wxID_STATIC, wxT("&min ="),
					wxDefaultPosition, wxDefaultSize, 0);
			aiminVal = new wxTextCtrl( this, ID_RANGE_AIMIN, wxT("0"),
					wxDefaultPosition, wxDefaultSize, 0);
			staticSizer2->Add( aiminName, 0, wxALL, 5 );
			staticSizer2->Add( aiminVal, 0, wxALL, 5 );

			wxStaticText* aimaxName= new wxStaticText( this, wxID_STATIC, wxT("&max ="),
					wxDefaultPosition, wxDefaultSize, 0);
			aimaxVal = new wxTextCtrl( this, ID_RANGE_AIMAX, wxT("0"),
					wxDefaultPosition, wxDefaultSize, 0);
			staticSizer2->Add( aimaxName, 0, wxALL, 5 );
			staticSizer2->Add( aimaxVal, 0, wxALL, 5 );

			boxSizer->Add( staticSizer2, 0,  wxALIGN_LEFT|wxALL, 5);

		//}  ///  ????  does not work here

        //========================================================
        //---------------- bottom row of action buttons ----------------------
        //---- horz sizer box for reset,OK,Cancel,Help
        wxBoxSizer* okCancelBox = new wxBoxSizer( wxHORIZONTAL );
        boxSizer->Add(okCancelBox, 0, wxALIGN_CENTER_HORIZONTAL|wxALL, 5);

        //---- OK button
        wxButton* ok = new wxButton( this, wxID_OK, wxT("&OK"),
                wxDefaultPosition, wxDefaultSize, 0);
		ok->SetDefault();
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
 bool rangeDialog::TransferDataToWindow( )
 {
        wxString ws;

       rminVal->SetValue(ws.Format(wxT("%g"), rmin) );
       rmaxVal->SetValue(ws.Format(wxT("%g"), rmax) );

	   if( npix > 1 ) {
           aiminVal->SetValue(ws.Format(wxT("%g"), aimin) );
           aimaxVal->SetValue(ws.Format(wxT("%g"), aimax) );

	   } else {
		   //  grey out imag part if its real
		   //  remember; .Enable() is inherited from wxWindow = base class for wxTextCtrl
		   aiminVal->Enable( false );
		   aimaxVal->Enable( false );
	   }

        return true;
 };

 //-------- TransferDataFromWindow() -----------------------
 bool rangeDialog::TransferDataFromWindow( )
 {
        double xx;

        nerror = 0;  // track number of read errors

        //  remember that GetValue() is from wxTextEntry which wxTextCrtl is derived from
        //

        if( rminVal->GetValue().ToDouble(&xx) ) rmin = xx;
                else nerror++;
        if( rmaxVal->GetValue().ToDouble(&xx) ) rmax = xx;
                else nerror++;

		// only get imaginary scle if its complex
	    if( npix > 1 ) {
			if( aiminVal->GetValue().ToDouble(&xx) ) aimin = xx;
					else nerror++;
			if( aimaxVal->GetValue().ToDouble(&xx) ) aimax = xx;
					else nerror++;
	   }

        return true;
 };

 
