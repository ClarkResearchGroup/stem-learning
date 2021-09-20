/*        generalDialog.cpp 

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

general setup dialog to go with computem (display/edit keV etc. values)

-------------------------------------------------
  references:

  1. Julian Smart and Kevin Hock (with Stefan Csomor), "Cross-Platform GUI
  Programming with wxWidgets", Prentice Hall, 2006
     see chp. 9  AND appendix J = whole code example

-------------------------------------------------

started from sampleDialog 21-may-2013 ejk
add ncellx,y 12-may-2013 ejk
add ctiltx,y 25-sep-2013 ejk
add capertmin,max 29-sep-2013 ejk
add model 3D sphere size 3-dec-2013 ejk
add CBED probe position 14-sep-2014 ejk
add STEM probe current and dwell time 29-sep-2015 ejk
fix some horz. spacing for wxWidgets 3.1.0/macosx/ubuntu  11-may-2016 ejk
fix close_box (fix create style) 5-jun-2016 ejk
last updated 5-jun-2016 ejk
change pA -> pAmp and microS -> microSec 18-aug-2019 ejk

formatted for a TAB size of 4 char

*/

#include "generalDialog.hpp"

IMPLEMENT_CLASS( generalDialog, wxDialog )


//-------- constructors -----------------------
generalDialog::generalDialog() { Init(); }

generalDialog::generalDialog(wxWindow *parent,
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
 bool generalDialog::Create( wxWindow* parent,
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
 void generalDialog::Init( )
 {
        return;
 };

 //-------- CreateControls() -----------------------
 void generalDialog::CreateControls( )
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
        // ----- energy ----------------------------------------
        //
        wxGridSizer *gridSizer5 = new wxGridSizer( 1, 2, 1, 1 );

        wxStaticBox* staticBox1 = new wxStaticBox( this, wxID_ANY, wxT("beam energy") );
        wxStaticBoxSizer* staticSizer1 = new wxStaticBoxSizer( staticBox1, wxHORIZONTAL );
     
        wxStaticText *kevName= new wxStaticText( this, wxID_STATIC, wxT("&kev"),
                wxDefaultPosition, wxDefaultSize, 0);
        kevVal = new wxTextCtrl( this, ID_GENERAL_KEV, wxT("0"),
                wxDefaultPosition, wxDefaultSize, 0);
        staticSizer1->Add( kevVal, 0, wxALL, 5 );
        staticSizer1->Add( kevName, 0, wxALL, 5 );
     
        // ---------------------------------------------
        // ----- objective aperture ----------------------------------------
        wxStaticBox* staticBox2 = new wxStaticBox( this, wxID_ANY, wxT("obj. apert.") );
        wxStaticBoxSizer* staticSizer2 = new wxStaticBoxSizer( staticBox2, wxHORIZONTAL );
     
        wxStaticText* objName= new wxStaticText( this, wxID_STATIC, wxT("&mrad."),
                wxDefaultPosition, wxDefaultSize, 0);
        objVal = new wxTextCtrl( this, ID_GENERAL_OBJ, wxT("00"),
                wxDefaultPosition, wxDefaultSize, 0);
        staticSizer2->Add( objVal, 0, wxALL, 5 );
        staticSizer2->Add( objName, 0, wxALL, 5 );

        //----  put energy and obj. apert. side by side
        gridSizer5->Add( staticSizer1, 0, wxALIGN_LEFT|wxALL, 5);
        gridSizer5->Add( staticSizer2, 0, wxALIGN_LEFT|wxALL, 5);
        boxSizer->Add( gridSizer5, 0, wxALIGN_CENTER_HORIZONTAL|wxALL, 5);

        // ---------------------------------------------
        // ----- crystal tilt ----------------------------------------
        wxStaticBox* staticBox3 = new wxStaticBox( this, wxID_ANY, wxT("crystal tilt (mrad)") );
        wxStaticBoxSizer* staticSizer3 = new wxStaticBoxSizer( staticBox3, wxHORIZONTAL );
        boxSizer->Add( staticSizer3, 0,  wxALIGN_CENTER_HORIZONTAL|wxALL, 5);
     
        wxStaticText* ctiltxName= new wxStaticText( this, wxID_STATIC, wxT("&ctiltx ="),
                wxDefaultPosition, wxDefaultSize, 0);
        ctiltxVal = new wxTextCtrl( this, ID_GENERAL_CTILTX, wxT("00"),
                wxDefaultPosition, wxDefaultSize, 0);
        staticSizer3->Add( ctiltxName, 0, wxALL, 5 );
        staticSizer3->Add( ctiltxVal, 0, wxALL, 5 );

        wxStaticText* ctiltyName= new wxStaticText( this, wxID_STATIC, wxT("&ctilty ="),
                wxDefaultPosition, wxDefaultSize, 0);
        ctiltyVal = new wxTextCtrl( this, ID_GENERAL_CTILTY, wxT("00"),
                wxDefaultPosition, wxDefaultSize, 0);
        staticSizer3->Add( ctiltyName, 0, wxALL, 5 );
        staticSizer3->Add( ctiltyVal, 0, wxALL, 5 );

        // ---------------------------------------------
        // ----- CTEM condencer angle ----------------------------------------
        wxStaticBox* staticBox4 = new wxStaticBox( this, wxID_ANY, wxT("CTEM condencer angle (mrad)") );
        wxStaticBoxSizer* staticSizer4 = new wxStaticBoxSizer( staticBox4, wxHORIZONTAL );
        boxSizer->Add( staticSizer4, 0,  wxALIGN_CENTER_HORIZONTAL|wxALL, 5);

        wxStaticText* capertminName= new wxStaticText( this, wxID_STATIC, wxT("&min. ="),
                wxDefaultPosition, wxDefaultSize, 0);
        capertminVal = new wxTextCtrl( this, ID_GENERAL_CAPERT1, wxT("00"),
                wxDefaultPosition, wxDefaultSize, 0);
        staticSizer4->Add( capertminName, 0, wxALL, 5 );
        staticSizer4->Add( capertminVal, 0, wxALL, 5 );

        wxStaticText* capertmaxName= new wxStaticText( this, wxID_STATIC, wxT("&max. ="),
                wxDefaultPosition, wxDefaultSize, 0);
        capertmaxVal = new wxTextCtrl( this, ID_GENERAL_CAPERT2, wxT("00"),
                wxDefaultPosition, wxDefaultSize, 0);
        staticSizer4->Add( capertmaxName, 0, wxALL, 5 );
        staticSizer4->Add( capertmaxVal, 0, wxALL, 5 );
     
        // ---------------------------------------------
        // ----- STEM stuff ----------------------------------------
        wxFlexGridSizer *gridSizer6 = new wxFlexGridSizer( 2, 1, 1 );

        // ----- source size ----------------------------------------
        wxStaticBox* staticBox5 = new wxStaticBox( this, wxID_ANY, wxT("STEM source size (incostem)") );
        wxStaticBoxSizer* staticSizer5 = new wxStaticBoxSizer( staticBox5, wxHORIZONTAL );
        //boxSizer->Add( staticSizer5, 0,  wxALIGN_LEFT|wxALL, 5);

        wxStaticText* dsrcName= new wxStaticText( this, wxID_STATIC, wxT("&nm."),
                wxDefaultPosition, wxDefaultSize, 0);
        dsrcVal = new wxTextCtrl( this, ID_GENERAL_SRC, wxT("00"),
                wxDefaultPosition, wxDefaultSize, 0);
        staticSizer5->Add( dsrcVal, 0, wxALL, 5 );
        staticSizer5->Add( dsrcName, 0, wxALL, 5 );


        // ----- probe current, dwell time ----------------------------------------
        wxStaticBox* staticBox6 = new wxStaticBox( this, wxID_ANY, wxT("STEM probe current (incostem)") );
        wxStaticBoxSizer* staticSizer6 = new wxStaticBoxSizer( staticBox6, wxHORIZONTAL );
        //boxSizer->Add( staticSizer5b, 0,  wxALIGN_LEFT|wxALL, 5);

        wxStaticText* prbIName= new wxStaticText( this, wxID_STATIC, wxT("&pAmp"),
                wxDefaultPosition, wxDefaultSize, 0);
        prbIval = new wxTextCtrl( this, ID_GENERAL_SRC, wxT("00"),
                wxDefaultPosition, wxDefaultSize, 0);
        staticSizer6->Add( prbIval, 0, wxALL, 5 );
        staticSizer6->Add( prbIName, 0, wxALL, 5 );

        wxStaticText* prbTName= new wxStaticText( this, wxID_STATIC, wxT("&microSec"),
                wxDefaultPosition, wxDefaultSize, 0);
        prbDTval = new wxTextCtrl( this, ID_GENERAL_SRC, wxT("00"),
                wxDefaultPosition, wxDefaultSize, 0);
        staticSizer6->Add( prbDTval, 0, wxALL, 5 );
        staticSizer6->Add( prbTName, 0, wxALL, 5 );

        //----  put probe stuff side by side
        gridSizer6->Add( staticSizer5, 0, wxALIGN_LEFT|wxALL, 5);
        gridSizer6->Add( staticSizer6, 0, wxALIGN_LEFT|wxALL, 5);
        boxSizer->Add( gridSizer6, 0, wxALIGN_LEFT|wxALL, 5);

        // ---------------------------------------------
        // ----- ADF detector ----------------------------------------
        wxStaticBox* staticBox7 = new wxStaticBox( this, wxID_ANY, wxT("STEM ADF detector (mrad)") );
        wxStaticBoxSizer* staticSizer7 = new wxStaticBoxSizer( staticBox7, wxHORIZONTAL );
        boxSizer->Add( staticSizer7, 0,  wxALIGN_CENTER_HORIZONTAL|wxALL, 5);

        wxStaticText* dminName= new wxStaticText( this, wxID_STATIC, wxT("&min. ="),
                wxDefaultPosition, wxDefaultSize, 0);
        dminVal = new wxTextCtrl( this, ID_GENERAL_DMIN, wxT("00"),
                wxDefaultPosition, wxDefaultSize, 0);
        staticSizer7->Add( dminName, 0, wxALL, 5 );
        staticSizer7->Add( dminVal, 0, wxALL, 5 );

        wxStaticText* dmaxName= new wxStaticText( this, wxID_STATIC, wxT("&max. ="),
                wxDefaultPosition, wxDefaultSize, 0);
        dmaxVal = new wxTextCtrl( this, ID_GENERAL_DMAX, wxT("00"),
                wxDefaultPosition, wxDefaultSize, 0);
        staticSizer7->Add( dmaxName, 0, wxALL, 5 );
        staticSizer7->Add( dmaxVal, 0, wxALL, 5 );

        // ---------------------------------------------
        // ----- CBED probe position ----------------------------------------
        wxStaticBox* staticBox8 = new wxStaticBox( this, wxID_ANY, wxT("CBED probe position (nm)") );
        wxStaticBoxSizer* staticSizer8 = new wxStaticBoxSizer( staticBox8, wxHORIZONTAL );
        boxSizer->Add( staticSizer8, 0,  wxALIGN_CENTER_HORIZONTAL|wxALL, 5);

        wxStaticText* probexName= new wxStaticText( this, wxID_STATIC, wxT("&x ="),
                wxDefaultPosition, wxDefaultSize, 0);
        probePosxVal = new wxTextCtrl( this, ID_GENERAL_DMIN, wxT("00"),
                wxDefaultPosition, wxDefaultSize, 0);
        staticSizer8->Add( probexName, 0, wxALL, 5 );
        staticSizer8->Add( probePosxVal, 0, wxALL, 5 );

        wxStaticText* probeyName= new wxStaticText( this, wxID_STATIC, wxT("&y ="),
                wxDefaultPosition, wxDefaultSize, 0);
        probePosyVal = new wxTextCtrl( this, ID_GENERAL_DMAX, wxT("00"),
                wxDefaultPosition, wxDefaultSize, 0);
        staticSizer8->Add( probeyName, 0, wxALL, 5 );
        staticSizer8->Add( probePosyVal, 0, wxALL, 5 );

        // ---------------------------------------------
        // ----- model 3D ----------------------------------------
        wxFlexGridSizer *gridSizer7 = new wxFlexGridSizer( 2, 1, 1 );
        //wxFlexGridSizer *gridSizer7 = new wxFlexGridSizer( 2 );

		wxStaticBox* staticBox9 = new wxStaticBox( this, wxID_ANY, wxT("3D model") );
        wxStaticBoxSizer* staticSizer9 = new wxStaticBoxSizer( staticBox9, wxHORIZONTAL );
        //boxSizer->Add( staticSizer8, 0,  wxALIGN_LEFT|wxALL, 5);

        wxStaticText* sphere3dName= new wxStaticText( this, wxID_STATIC, wxT("&rel. sphere size="),
                wxDefaultPosition, wxDefaultSize, 0);
        sphere3dVal = new wxTextCtrl( this, ID_GENERAL_SPHER3D, wxT("00"),
                wxDefaultPosition, wxDefaultSize, 0);
        staticSizer9->Add( sphere3dName, 0, wxALL, 5 );
        staticSizer9->Add( sphere3dVal, 0, wxALL, 5 );

        // ---------------------------------------------
        // ----- log greyscale constant ----------------------------------------
        wxStaticBox* staticBox10 = new wxStaticBox( this, wxID_ANY, wxT("log grey constant") );
        wxStaticBoxSizer* staticSizer10 = new wxStaticBoxSizer( staticBox10, wxHORIZONTAL );
     
        GreyLogVal = new wxTextCtrl( this, ID_GENERAL_OBJ, wxT("00"),
                wxDefaultPosition, wxDefaultSize, 0);
        staticSizer10->Add( GreyLogVal, 0, wxALL, 5 );

        //----  put energy and obj. apert. side by side
        gridSizer7->Add( staticSizer9, 1, wxALIGN_CENTER_HORIZONTAL, 5);
        gridSizer7->Add( staticSizer10, 1, wxALIGN_CENTER_HORIZONTAL, 5);
        boxSizer->Add( gridSizer7, 1, wxALIGN_CENTER_HORIZONTAL, 5);

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
 bool generalDialog::TransferDataToWindow( )
 {
        wxString ws;

        kevVal->SetValue(ws.Format(wxT("%10.3f"), kev) );

        objVal->SetValue(ws.Format(wxT("%10.3f"), objAp*1000.0) );	// in mrad.

        capertminVal->SetValue(ws.Format(wxT("%10.3f"), capertmin*1000.0) );	// in mrad.
        capertmaxVal->SetValue(ws.Format(wxT("%10.3f"), capertmax*1000.0) );

        dsrcVal->SetValue(ws.Format(wxT("%10.3f"), dsource/10.0) );	// in nm.
        prbIval->SetValue(ws.Format(wxT("%10.3f"), probeI/1.0e-12) );	// in pA
        prbDTval->SetValue(ws.Format(wxT("%10.3f"), probeDt/1.0e-6) );	// in microSec

        dminVal->SetValue(ws.Format(wxT("%10.3f"), detmin*1000.0) );	// in mrad.
        dmaxVal->SetValue(ws.Format(wxT("%10.3f"), detmax*1000.0) );

        ctiltxVal->SetValue(ws.Format(wxT("%10.3f"), ctiltx*1000.0) );	// in mrad.
        ctiltyVal->SetValue(ws.Format(wxT("%10.3f"), ctilty*1000.0) );

        probePosxVal->SetValue(ws.Format(wxT("%10.3f"), probePosx/10.0) );	// in nm.
        probePosyVal->SetValue(ws.Format(wxT("%10.3f"), probePosy/10.0) );

        sphere3dVal->SetValue(ws.Format(wxT("%10.3f"), sphere3d) );

		GreyLogVal->SetValue(ws.Format(wxT("%10.3f"), GreyLogConstant) );
		 
        return true;
 };

 //-------- TransferDataFromWindow() -----------------------
 bool generalDialog::TransferDataFromWindow( )
 {
        double xx;

        nerror = 0;  // track number of read errors

        //  remember that GetValue() is from wxTextEntry which wxTextCrtl is derived from
        //
        //-----   all of these must be positive so force absolute value 
        //
        if( kevVal->GetValue().ToDouble(&xx) ) kev = fabs(xx);
                else nerror++;

        if( objVal->GetValue().ToDouble(&xx) ) objAp = fabs(xx/1000.0);  // in rad.
                else nerror++;

        if( capertminVal->GetValue().ToDouble(&xx) ) capertmin = fabs(xx/1000.0);  // in rad.
                else nerror++;
        if( capertmaxVal->GetValue().ToDouble(&xx) ) capertmax = fabs(xx/1000.0);  // in rad.
                else nerror++;

        if( dsrcVal->GetValue().ToDouble(&xx) ) dsource = fabs(xx*10.0);  // in Angstroms
                else nerror++;
        if( prbIval->GetValue().ToDouble(&xx) ) probeI = fabs(xx*1.0e-12);  // in Amp
                else nerror++;
        if( prbDTval->GetValue().ToDouble(&xx) ) probeDt = xx*1.0e-6;  // in microSec
                else nerror++;

        if( dminVal->GetValue().ToDouble(&xx) ) detmin = fabs(xx/1000.0);  // in rad.
                else nerror++;
        if( dmaxVal->GetValue().ToDouble(&xx) ) detmax = fabs(xx/1000.0);  // in rad.
                else nerror++;

        //  these can be pos. or neg.
        if( ctiltxVal->GetValue().ToDouble(&xx) ) ctiltx = xx/1000.0;  // in rad.
                else nerror++;
        if( ctiltyVal->GetValue().ToDouble(&xx) ) ctilty = xx/1000.0;  // in rad.
                else nerror++;

        //  these can be pos. or neg.
        if( probePosxVal->GetValue().ToDouble(&xx) ) probePosx = xx*10.0;  // in Ang.
                else nerror++;
        if( probePosyVal->GetValue().ToDouble(&xx) ) probePosy = xx*10.0;  // in Ang.
                else nerror++;

		//  must be positive
		if( sphere3dVal->GetValue().ToDouble(&xx) ) sphere3d = fabs(xx);  // in rad.
                else nerror++;

		//  must be positive
		if( GreyLogVal->GetValue().ToDouble(&xx) ) if(xx>0.0) GreyLogConstant = fabs(xx);  
                else nerror++;


        return true;
 };

 
