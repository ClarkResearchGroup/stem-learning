/*	  aberDialog.cpp 

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
add scaling constants for Haider aber. symbols 9-may-2016 ejk
fix some horz. spacing for wxWidgets 3.1.0/macosx/ubuntu  11-may-2016 ejk
fix close_box (fix create style) 5-jun-2016 ejk
add option to add pi/4 random tuning errors 11-jun-2016 ejk
last updated 11-jun-2016 ejk

formatted for a TAB size of 4 char

*/

#include "aberDialog.hpp"

#include "slicelib.hpp"    // for rangauss()

IMPLEMENT_CLASS( aberDialog, wxDialog )

BEGIN_EVENT_TABLE( aberDialog, wxDialog )
	EVT_BUTTON( ID_ABERR_ZERO, aberDialog::OnZero )
	EVT_BUTTON( ID_ADD_ERRORS, aberDialog::OnAddErrors )
END_EVENT_TABLE()

//-------- constructors -----------------------
aberDialog::aberDialog() { Init(); }

aberDialog::aberDialog(wxWindow *parent,
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
 bool aberDialog::Create( wxWindow* parent,
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
 void aberDialog::Init( )
 {
	scale1 = 10.0;	 // order 1 in nm
	scale2 = 10.0;	 // order 1 in nm
	scale3 = 1.0e4;   // order 1 in um
	scale4 = 1.0e4;   // order 1 in um
	scale5 = 1.0e7;   // order 1 in mm
 };

 //-------- CreateControls() -----------------------
 void aberDialog::CreateControls( )
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
	//  put Cnm in lots of little boxes

	wxStaticBox* staticBox1 = new wxStaticBox( this, wxID_ANY, wxT("first order (nm)") );
	wxStaticBoxSizer* staticSizer1 = new wxStaticBoxSizer( staticBox1, wxHORIZONTAL );
	boxSizer->Add( staticSizer1, 0, wxALIGN_LEFT|wxALL, 5);

	// ------ start boxes for aberration values
	// ---------------------------------------------
	// ----- C1x ----------------------------------------
	//  remember that wxSpinCtrl() only does int values NOT floating point!!!
	//
	wxStaticText *dfName= new wxStaticText( this, wxID_STATIC, wxT("&(-C10)df="),
		wxDefaultPosition, wxDefaultSize, 0);
	dfValue = new wxTextCtrl( this, ID_ABERR_DF, wxT("0.00"),
		wxDefaultPosition, wxDefaultSize, 0);
	staticSizer1->Add( dfName, 0, wxALL, 5);
	staticSizer1->Add( dfValue, 0, wxALL, 5);

	wxStaticText *ddfName= new wxStaticText( this, wxID_STATIC, wxT("&ddf="),
		wxDefaultPosition, wxDefaultSize, 0);
	ddfValue = new wxTextCtrl( this, ID_ABERR_DDF, wxT("0.00"),
		wxDefaultPosition, wxDefaultSize, 0);
	staticSizer1->Add( ddfName, 0, wxALL, 5);
	staticSizer1->Add( ddfValue, 0, wxALL, 5);

	wxStaticText* C12Name= new wxStaticText( this, wxID_STATIC, wxT("&A1=C12ab="),
		wxDefaultPosition, wxDefaultSize, 0);
	C12aVal = new wxTextCtrl( this, ID_ABERR_C12A, wxT("0.00"),
		wxDefaultPosition, wxDefaultSize, 0);
	C12bVal = new wxTextCtrl( this, ID_ABERR_C12B, wxT("0.00"),
		wxDefaultPosition, wxDefaultSize, 0);
	staticSizer1->Add( C12Name, 0, wxALL, 5);
	staticSizer1->Add( C12aVal, 0, wxALL, 5);
	staticSizer1->Add( C12bVal, 0, wxALL, 5);
	//staticSizer1->AddSpacer(5);

	// ---------------------------------------------
	// ----- C2x ----------------------------------------
	wxStaticBox* staticBox2 = new wxStaticBox( this, wxID_ANY, wxT("second order (nm)") );
	wxStaticBoxSizer* staticSizer2 = new wxStaticBoxSizer( staticBox2, wxHORIZONTAL );
	boxSizer->Add( staticSizer2, 0,  wxALIGN_LEFT|wxALL, 5);

	wxStaticText* C21Name= new wxStaticText( this, wxID_STATIC, wxT("&3*B2=C21ab="),
		wxDefaultPosition, wxDefaultSize, 0);
	C21aVal = new wxTextCtrl( this, ID_ABERR_C21A, wxT("0.00"),
		wxDefaultPosition, wxDefaultSize, 0);
	C21bVal = new wxTextCtrl( this, ID_ABERR_C21B, wxT("0.00"),
		wxDefaultPosition, wxDefaultSize, 0);
	staticSizer2->Add( C21Name, 0, wxALL, 5);
	staticSizer2->Add( C21aVal, 0, wxALL, 5);
	staticSizer2->Add( C21bVal, 0, wxALL, 5);

	wxStaticText* C23Name= new wxStaticText( this, wxID_STATIC, wxT("&A2=C23ab="),
		wxDefaultPosition, wxDefaultSize, 0);
	C23aVal = new wxTextCtrl( this, ID_ABERR_C23A, wxT("0.00"),
		wxDefaultPosition, wxDefaultSize, 0);
	C23bVal = new wxTextCtrl( this, ID_ABERR_C23B, wxT("0.00"),
		wxDefaultPosition, wxDefaultSize, 0); 
		//wxDefaultPosition, wxDefaultSize, wxTE_RIGHT);   //  right justified ????
	staticSizer2->Add( C23Name, 0, wxALL, 5);
	staticSizer2->Add( C23aVal, 0, wxALL, 5);
	staticSizer2->Add( C23bVal, 0, wxALL, 5);

	// ---------------------------------------------
	// ----- C3x ----------------------------------------
	wxStaticBox* staticBox3 = new wxStaticBox( this, wxID_ANY, wxT("third order (um)") );
	wxStaticBoxSizer* staticSizer3 = new wxStaticBoxSizer( staticBox3, wxHORIZONTAL );
	boxSizer->Add( staticSizer3, 0,  wxALIGN_LEFT|wxALL, 5);

	wxStaticText* C30Name= new wxStaticText( this, wxID_STATIC, wxT("&C3=C30="),
		wxDefaultPosition, wxDefaultSize, 0);
	C30Val = new wxTextCtrl( this, ID_ABERR_C30, wxT("0.00"),
		wxDefaultPosition, wxDefaultSize, 0);
	staticSizer3->Add( C30Name, 0, wxALL, 5);
	staticSizer3->Add( C30Val,  0, wxALL, 5);

	wxStaticText* C32Name= new wxStaticText( this, wxID_STATIC, wxT("&4*S3=C32ab="),
		wxDefaultPosition, wxDefaultSize, 0);
	C32aVal = new wxTextCtrl( this, ID_ABERR_C32A, wxT("0.00"),
		wxDefaultPosition, wxDefaultSize, 0);
	C32bVal = new wxTextCtrl( this, ID_ABERR_C32B, wxT("0.00"),
		wxDefaultPosition, wxDefaultSize, 0);
	staticSizer3->Add( C32Name, 0, wxALL, 5);
	staticSizer3->Add( C32aVal, 0, wxALL, 5);
	staticSizer3->Add( C32bVal,  0, wxALL, 5);

	wxStaticText* C34Name= new wxStaticText( this, wxID_STATIC, wxT("&A3=C34ab="),
		wxDefaultPosition, wxDefaultSize, 0);
	C34aVal = new wxTextCtrl( this, ID_ABERR_C34A, wxT("0.00"),
		wxDefaultPosition, wxDefaultSize, 0);
	C34bVal = new wxTextCtrl( this, ID_ABERR_C34B, wxT("0.00"),
		wxDefaultPosition, wxDefaultSize, 0);
	staticSizer3->Add( C34Name, 0, wxALL, 5);
	staticSizer3->Add( C34aVal, 0, wxALL, 5);
	staticSizer3->Add( C34bVal,  0, wxALL, 5);

	// ---------------------------------------------
	// ----- C4x ----------------------------------------
	wxStaticBox* staticBox4 = new wxStaticBox( this, wxID_ANY, wxT("forth order (um)") );
	wxStaticBoxSizer* staticSizer4 = new wxStaticBoxSizer( staticBox4, wxHORIZONTAL );
	boxSizer->Add( staticSizer4, 0,  wxALIGN_LEFT|wxALL, 5);

	wxStaticText* C41Name= new wxStaticText( this, wxID_STATIC, wxT("&4*B4=C41ab ="),
		wxDefaultPosition, wxDefaultSize, 0);
	C41aVal = new wxTextCtrl( this, ID_ABERR_C41A, wxT("0.00"),
		wxDefaultPosition, wxDefaultSize, 0);
	C41bVal = new wxTextCtrl( this, ID_ABERR_C41B, wxT("0.00"),
		wxDefaultPosition, wxDefaultSize, 0);
	staticSizer4->Add( C41Name, 0, wxALL, 5);
	staticSizer4->Add( C41aVal, 0, wxALL, 5);
	staticSizer4->Add( C41bVal,  0, wxALL, 5);

	wxStaticText* C43Name= new wxStaticText( this, wxID_STATIC, wxT("&4*D4=C43ab="),
		wxDefaultPosition, wxDefaultSize, 0);
	C43aVal = new wxTextCtrl( this, ID_ABERR_C43A, wxT("0.00"),
		wxDefaultPosition, wxDefaultSize, 0);
	C43bVal = new wxTextCtrl( this, ID_ABERR_C43B, wxT("0.00"),
		wxDefaultPosition, wxDefaultSize, 0);
	staticSizer4->Add( C43Name, 0, wxALL, 5);
	staticSizer4->Add( C43aVal, 0, wxALL, 5);
	staticSizer4->Add( C43bVal,  0, wxALL, 5);

	wxStaticText* C45Name= new wxStaticText( this, wxID_STATIC, wxT("&A4=C45ab="),
		wxDefaultPosition, wxDefaultSize, 0);
	C45aVal = new wxTextCtrl( this, ID_ABERR_C45A, wxT("0.00"),
		wxDefaultPosition, wxDefaultSize, 0);
	C45bVal = new wxTextCtrl( this, ID_ABERR_C45B, wxT("0.00"),
		wxDefaultPosition, wxDefaultSize, 0);
	staticSizer4->Add( C45Name, 0, wxALL, 5);
	staticSizer4->Add( C45aVal, 0, wxALL, 5);
	staticSizer4->Add( C45bVal,  0, wxALL, 5);

	// ---------------------------------------------
	// ----- C5x ----------------------------------------

	wxStaticBox* staticBox5 = new wxStaticBox( this, wxID_ANY, wxT("fifth order (mm)") );
	wxStaticBoxSizer* staticSizer5 = new wxStaticBoxSizer( staticBox5, wxHORIZONTAL );
	boxSizer->Add( staticSizer5, 0,  wxALIGN_LEFT|wxALL, 5);

	//   too many boxes so make a 2D grid to go inside the box
	wxGridSizer *gridSizer5 = new wxGridSizer( 2, 8, 1, 1 );

	wxStaticText* C50Name= new wxStaticText( this, wxID_STATIC, wxT("&C5=C50="),
		wxDefaultPosition, wxDefaultSize, 0);
	C50Val = new wxTextCtrl( this, ID_ABERR_C50, wxT("0.00"),
		wxDefaultPosition, wxDefaultSize, 0);
	gridSizer5->Add( C50Name, 0, wxALL, 5);
	gridSizer5->Add( C50Val,  0, wxALL, 5);

	wxStaticText* C52Name= new wxStaticText( this, wxID_STATIC, wxT("&6*S5=C52ab="),
		wxDefaultPosition, wxDefaultSize, 0);
	C52aVal = new wxTextCtrl( this, ID_ABERR_C52A, wxT("0.00"),
		wxDefaultPosition, wxDefaultSize, 0);
	C52bVal = new wxTextCtrl( this, ID_ABERR_C52B, wxT("0.00"),
		wxDefaultPosition, wxDefaultSize, 0);
	gridSizer5->Add( C52Name, 0, wxALL, 5);
	gridSizer5->Add( C52aVal, 0, wxALL, 5);
	gridSizer5->Add( C52bVal, 0, wxALL, 5);

	wxStaticText* C54Name= new wxStaticText( this, wxID_STATIC, wxT("&6*R5=C54ab="),
		wxDefaultPosition, wxDefaultSize, 0);
	C54aVal = new wxTextCtrl( this, ID_ABERR_C54A, wxT("0.00"),
		wxDefaultPosition, wxDefaultSize, 0);
	C54bVal = new wxTextCtrl( this, ID_ABERR_C54B, wxT("0.00"),
		wxDefaultPosition, wxDefaultSize, 0);
	gridSizer5->Add( C54Name, 0, wxALL, 5);
	gridSizer5->Add( C54aVal, 0, wxALL, 5);
	gridSizer5->Add( C54bVal, 0, wxALL, 5);

	//  remember; it makes all col. the same width so this doesn't space correctly
	//     but it will have to do
	gridSizer5->AddSpacer( 1 );   // start second row
	gridSizer5->AddSpacer( 1 );   //   and space over C50

	wxStaticText* C56Name= new wxStaticText( this, wxID_STATIC, wxT("&A5=C56ab="),
		wxDefaultPosition, wxDefaultSize, 0);
	C56aVal = new wxTextCtrl( this, ID_ABERR_C56A, wxT("0.00"),
		wxDefaultPosition, wxDefaultSize, 0);
	C56bVal = new wxTextCtrl( this, ID_ABERR_C56B, wxT("0.00"),
		wxDefaultPosition, wxDefaultSize, 0);
	gridSizer5->Add( C56Name, 0, wxALL, 5);
	gridSizer5->Add( C56aVal, 0, wxALL, 5);
	gridSizer5->Add( C56bVal, 0, wxALL, 5);

	//  and finally put the grid inside its box
	staticSizer5->Add( gridSizer5, 0, wxALL, 5);  //  , wxALIGN_CENTER_HORIZONTAL|wxALL, 5);

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

	//---- reset button
	wxButton* reset = new wxButton( this, ID_ABERR_ZERO, wxT("&Zero All"),
		wxDefaultPosition, wxDefaultSize, 0);
	okCancelBox->Add(reset, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);
 
	//---- random errors button
	wxButton* addErrors = new wxButton( this, ID_ADD_ERRORS,
		wxT("&add random pi/4 tuning errors"),
		wxDefaultPosition, wxDefaultSize, 0);
	okCancelBox->Add(addErrors, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);
};

 //========================================================
//-------- TransferDataToWindow() -----------------------
 //   update values on screen when first started
 //    CreateControls() doesn't have right values for some reason
 //
 bool aberDialog::TransferDataToWindow( )
 {
	wxString ws;

	dfValue->SetValue(ws.Format(wxT("%f"), df/scale1) );
	ddfValue->SetValue(ws.Format(wxT("%f"), ddf/scale1) );

	C12aVal->SetValue(ws.Format(wxT("%f"), C12a/scale1) );
	C12bVal->SetValue(ws.Format(wxT("%f"), C12b/scale1) );

	C21aVal->SetValue(ws.Format(wxT("%f"), C21a/scale2) );
	C21bVal->SetValue(ws.Format(wxT("%f"), C21b/scale2) );

	C23aVal->SetValue(ws.Format(wxT("%f"), C23a/scale2) );
	C23bVal->SetValue(ws.Format(wxT("%f"), C23b/scale2) );

	C30Val->SetValue(ws.Format(wxT("%f"), C30/scale3) );
	C32aVal->SetValue(ws.Format(wxT("%f"), C32a/scale3) );
	C32bVal->SetValue(ws.Format(wxT("%f"), C32b/scale3) );
	C34aVal->SetValue(ws.Format(wxT("%f"), C34a/scale3) );
	C34bVal->SetValue(ws.Format(wxT("%f"), C34b/scale3) );

	C41aVal->SetValue(ws.Format(wxT("%f"), C41a/scale4) );
	C41bVal->SetValue(ws.Format(wxT("%f"), C41b/scale4) );
	C43aVal->SetValue(ws.Format(wxT("%f"), C43a/scale4) );
	C43bVal->SetValue(ws.Format(wxT("%f"), C43b/scale4) );
	C45aVal->SetValue(ws.Format(wxT("%f"), C45a/scale4) );
	C45bVal->SetValue(ws.Format(wxT("%f"), C45b/scale4) );

	C50Val->SetValue(ws.Format(wxT("%f"), C50/scale5) );
	C52aVal->SetValue(ws.Format(wxT("%f"), C52a/scale5) );
	C52bVal->SetValue(ws.Format(wxT("%f"), C52b/scale5) );
	C54aVal->SetValue(ws.Format(wxT("%f"), C54a/scale5) );
	C54bVal->SetValue(ws.Format(wxT("%f"), C54b/scale5) );
	C56aVal->SetValue(ws.Format(wxT("%f"), C56a/scale5) );
	C56bVal->SetValue(ws.Format(wxT("%f"), C56b/scale5) );

	return true;
 };

 //-------- TransferDataFromWindow() -----------------------
 bool aberDialog::TransferDataFromWindow( )
 {
	double x;

	nerror = 0;  // track number of read errors

	//--------------  C1x -------------------------------
	if( dfValue->GetValue().ToDouble(&x) ) df = x * scale1; //  defocus
		else nerror++;
	if( ddfValue->GetValue().ToDouble(&x) ) ddf = x * scale1; // defocus spread	
		else nerror++;

	if( C12aVal->GetValue().ToDouble(&x) ) C12a = x * scale1;
		else nerror++;
	if( C12bVal->GetValue().ToDouble(&x) ) C12b = x * scale1;
		else nerror++;

	//--------------  C2x -------------------------------
	if( C21aVal->GetValue().ToDouble(&x) ) C21a = x * scale2;
		else nerror++;
	if( C21bVal->GetValue().ToDouble(&x) ) C21b = x * scale2;
		else nerror++;
	if( C23aVal->GetValue().ToDouble(&x) ) C23a = x * scale2;
		else nerror++;
	if( C23bVal->GetValue().ToDouble(&x) ) C23b = x * scale2;
		else nerror++;

	//--------------  C3x -------------------------------
	if( C30Val->GetValue().ToDouble(&x) ) C30 = x * scale3;
		else nerror++;
	if( C32aVal->GetValue().ToDouble(&x) ) C32a = x * scale3;
		else nerror++;
	if( C32bVal->GetValue().ToDouble(&x) ) C32b = x * scale3;
		else nerror++;
	if( C34aVal->GetValue().ToDouble(&x) ) C34a = x * scale3;
		else nerror++;
	if( C34bVal->GetValue().ToDouble(&x) ) C34b = x * scale3;
		else nerror++;

	//--------------  C4x -------------------------------
	if( C41aVal->GetValue().ToDouble(&x) ) C41a = x * scale4;
		else nerror++;
	if( C41bVal->GetValue().ToDouble(&x) ) C41b = x * scale4;
		else nerror++;
	if( C43aVal->GetValue().ToDouble(&x) ) C43a = x * scale4;
		else nerror++;
	if( C43bVal->GetValue().ToDouble(&x) ) C43b = x * scale4;
		else nerror++;
	if( C45aVal->GetValue().ToDouble(&x) ) C45a = x * scale4;
		else nerror++;
	if( C45bVal->GetValue().ToDouble(&x) ) C45b = x * scale4;
		else nerror++;

	//--------------  C5x -------------------------------
	if( C50Val->GetValue().ToDouble(&x) ) C50 = x * scale5;
		else nerror++;
	if( C52aVal->GetValue().ToDouble(&x) ) C52a = x * scale5;
		else nerror++;
	if( C52bVal->GetValue().ToDouble(&x) ) C52b = x * scale5;
		else nerror++;
	if( C54aVal->GetValue().ToDouble(&x) ) C54a = x * scale5;
		else nerror++;
	if( C54bVal->GetValue().ToDouble(&x) ) C54b = x * scale5;
		else nerror++;
	if( C56aVal->GetValue().ToDouble(&x) ) C56a = x * scale5;
		else nerror++;
	if( C56bVal->GetValue().ToDouble(&x) ) C56b = x * scale5;
		else nerror++;

	return true;
 };

//-------- ZeroAll() -----------------------
//   reset all aberrations to zero
 void aberDialog::OnZero( wxCommandEvent &event )
 {
	df = 0.0;
	ddf = 0.0;

	C12a = 0;
	C12b = 0;

	C21a = 0;
	C21b = 0;
	C23a = 0;
	C23b = 0;

	C30 = 0;
	C32a = 0;
	C32b = 0;
	C34a = 0;
	C34b = 0;

	C41a = 0;
	C41b = 0;
	C43a = 0;
	C43b = 0;
	C45a = 0;
	C45b = 0;

	C50 = 0;
	C52a = 0;
	C52b = 0;
	C54a = 0;
	C54b = 0;
	C56a = 0;
	C56b = 0;

	 //Init();  // prob. not needed here
	 TransferDataToWindow();
 };

//-------- OnAddErrors() -----------------------
//   add random pi/4 random tuning errors
//   to all but defocus and astigmatism
//
 void aberDialog::OnAddErrors( wxCommandEvent &event )
 {
	double oa2, dchiMax2, dchiMax3, dchiMax4, dchiMax5;

	wxString ws;   //  to format messages

	df = 0.0;
	ddf = 0.0;

	C12a = 0;
	C12b = 0;

	//------ calculate max aberr. error for pi/4 phase at max obj. angle
	oa2 = objApert * objApert * objApert;
	dchiMax2 = (3.0*wavl)/( 8.0 * oa2 );
	oa2 *= objApert;
	dchiMax3 = (4.0*wavl)/( 8.0 * oa2);
	oa2 *= objApert;
	dchiMax4 = (5.0*wavl)/( 8.0 * oa2);
	oa2 *= objApert;
	dchiMax5 = (6.0*wavl)/( 8.0 * oa2);

    wxLogStatus( ws.Format(
		wxT("add random tuning errors for %g mrad and wavelength %g Ang."),
			objApert*1000.0, wavl)  );

	C21a += dchiMax2*(2*ranflat(iseed) - 1);  // up to +/- err 
	C21b += dchiMax2*(2*ranflat(iseed) - 1);
	C23a += dchiMax2*(2*ranflat(iseed) - 1);
	C23b += dchiMax2*(2*ranflat(iseed) - 1);

	C30 += dchiMax3*(2*ranflat(iseed) - 1);
	C32a += dchiMax3*(2*ranflat(iseed) - 1);
	C32b += dchiMax3*(2*ranflat(iseed) - 1);
	C34a += dchiMax3*(2*ranflat(iseed) - 1);
	C34b += dchiMax3*(2*ranflat(iseed) - 1);

	C41a += dchiMax4*(2*ranflat(iseed) - 1);
	C41b += dchiMax4*(2*ranflat(iseed) - 1);
	C43a += dchiMax4*(2*ranflat(iseed) - 1);
	C43b += dchiMax4*(2*ranflat(iseed) - 1);
	C45a += dchiMax4*(2*ranflat(iseed) - 1);
	C45b += dchiMax4*(2*ranflat(iseed) - 1);

	C50 += dchiMax5*(2*ranflat(iseed) - 1);
	C52a += dchiMax5*(2*ranflat(iseed) - 1);
	C52b += dchiMax5*(2*ranflat(iseed) - 1);
	C54a += dchiMax5*(2*ranflat(iseed) - 1);
	C54b += dchiMax5*(2*ranflat(iseed) - 1);
	C56a += dchiMax5*(2*ranflat(iseed) - 1);
	C56b += dchiMax5*(2*ranflat(iseed) - 1);

	 TransferDataToWindow();
 };