/*        sampleDialog.cpp 

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
update parm readback 7-may-2013 ejk
add ncellx,y 12-may-2013 ejk
add nwobble and temp. 28-sep-2013 ejk
fix some horz. spacing for wxWidgets 3.1.0/macosx/ubuntu  11-may-2016 ejk
fix close_box (add to style) 5-jun-2016 ejk
last updated 5-jun-2016 ejk

*/

#include "sampleDialog.hpp"

IMPLEMENT_CLASS( sampleDialog, wxDialog )

//BEGIN_EVENT_TABLE( aberDialog, wxDialog )
//	EVT_BUTTON( ID_ABERR_ZERO, aberDialog::OnZero )
//END_EVENT_TABLE()

//-------- constructors -----------------------
sampleDialog::sampleDialog() { Init(); }

sampleDialog::sampleDialog(wxWindow *parent,
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
 bool sampleDialog::Create( wxWindow* parent,
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
 void sampleDialog::Init( )
 {
        return;
 };

 //-------- CreateControls() -----------------------
 void sampleDialog::CreateControls( )
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

        //  put the first two together
        //  remember; all col will be same width so make both same width
        wxGridSizer *gridSizer5 = new wxGridSizer( 1, 2, 1, 1 );

        // ----- wave function size ----------------------------------------
        //
        wxStaticBox* staticBox1 = new wxStaticBox( this, wxID_ANY, 
		wxT("transmission function size (pixels)") );
        wxStaticBoxSizer* staticSizer1 = new wxStaticBoxSizer( staticBox1, wxHORIZONTAL );

        wxStaticText *nxName= new wxStaticText( this, wxID_STATIC, wxT("&nx ="),
                wxDefaultPosition, wxDefaultSize, 0);
        nxVal = new wxTextCtrl( this, ID_SAMPLE_NX, wxT("0"),
                wxDefaultPosition, wxDefaultSize, 0);
        staticSizer1->Add( nxName, 0, wxALL, 5 );
        staticSizer1->Add( nxVal, 0, wxALL, 5 );

        wxStaticText* nyName= new wxStaticText( this, wxID_STATIC, wxT("&ny ="),
                wxDefaultPosition, wxDefaultSize, 0);
        nyVal = new wxTextCtrl( this, ID_SAMPLE_NY, wxT("0"),
                wxDefaultPosition, wxDefaultSize, 0);
        staticSizer1->Add( nyName, 0, wxALL, 5 );
        staticSizer1->Add( nyVal, 0, wxALL, 5 );


        // ---------------------------------------------
        // ----- slice thickness ----------------------------------------
        wxStaticBox* staticBox2 = new wxStaticBox( this, wxID_ANY, wxT("slice thickness (in nm)") );
        wxStaticBoxSizer* staticSizer2 = new wxStaticBoxSizer( staticBox2, wxHORIZONTAL );
        //boxSizer->Add( staticSizer2, 0,  wxALIGN_LEFT|wxALL, 5);

        wxStaticText* deltazName= new wxStaticText( this, wxID_STATIC, wxT("&deltaz ="),
                wxDefaultPosition, wxDefaultSize, 0);
        deltazVal = new wxTextCtrl( this, ID_SAMPLE_DELTAZ, wxT("00"),
                wxDefaultPosition, wxDefaultSize, 0);
        staticSizer2->Add( deltazName, 0, wxALL, 5 );
        staticSizer2->Add( deltazVal, 0, wxALL, 5 );

        //----  put energy and obj. apert. side by side
        gridSizer5->Add( staticSizer1, 0, wxALL, 5 );
        gridSizer5->Add( staticSizer2, 0, wxALL, 5 );
        boxSizer->Add( gridSizer5, 0, wxALL, 5 );

        // ---------------------------------------------
        // ----- unit cell size ----------------------------------------
        wxStaticBox* staticBox3 = new wxStaticBox( this, wxID_ANY, 
		wxT("super cell size (in unit cells)") );
        wxStaticBoxSizer* staticSizer3 = new wxStaticBoxSizer( staticBox3, wxHORIZONTAL );
        boxSizer->Add( staticSizer3, 0,  wxALIGN_LEFT|wxALL, 5);

        wxStaticText* ncellxName= new wxStaticText( this, wxID_STATIC, wxT("&ncellx ="),
                wxDefaultPosition, wxDefaultSize, 0);
        ncellxVal = new wxTextCtrl( this, ID_SAMPLE_NCELLX, wxT("00"),
                wxDefaultPosition, wxDefaultSize, 0);
        staticSizer3->Add( ncellxName, 0, wxALL, 5 );
        staticSizer3->Add( ncellxVal, 0, wxALL, 5 );

        wxStaticText* ncellyName= new wxStaticText( this, wxID_STATIC, wxT("&ncelly ="),
                wxDefaultPosition, wxDefaultSize, 0);
        ncellyVal = new wxTextCtrl( this, ID_SAMPLE_NCELLY, wxT("00"),
                wxDefaultPosition, wxDefaultSize, 0);
        staticSizer3->Add( ncellyName, 0, wxALL, 5 );
        staticSizer3->Add( ncellyVal, 0, wxALL, 5 );

        wxStaticText* ncellzName= new wxStaticText( this, wxID_STATIC, wxT("&ncellz ="),
                wxDefaultPosition, wxDefaultSize, 0);
        ncellzVal = new wxTextCtrl( this, ID_SAMPLE_NCELLZ, wxT("00"),
                wxDefaultPosition, wxDefaultSize, 0);
        staticSizer3->Add( ncellzName, 0, wxALL, 5 );
        staticSizer3->Add( ncellzVal, 0, wxALL, 5 );

        // ----- phonon configurations ----------------------------------------
        //
        wxStaticBox* staticBox11 = new wxStaticBox( this, wxID_ANY, wxT("frozen phonon configurations") );
        wxStaticBoxSizer* staticSizer11 = new wxStaticBoxSizer( staticBox11, wxHORIZONTAL );
        boxSizer->Add( staticSizer11, 0,  wxALIGN_LEFT|wxALL, 5);

        wxStaticText *nwobbleName= new wxStaticText( this, wxID_STATIC, wxT("&nwobble ="),
                wxDefaultPosition, wxDefaultSize, 0);
        nwobbleVal = new wxTextCtrl( this, ID_SAMPLE_NWOBBLE, wxT("0"),
                wxDefaultPosition, wxDefaultSize, 0);
        staticSizer11->Add( nwobbleName, 0, wxALL, 5 );
        staticSizer11->Add( nwobbleVal, 0, wxALL, 5 );

        wxStaticText *tempName= new wxStaticText( this, wxID_STATIC, wxT("&temperature (deg. K) ="),
                wxDefaultPosition, wxDefaultSize, 0);
        tempVal = new wxTextCtrl( this, ID_SAMPLE_TEMPERAT, wxT("0"),
                wxDefaultPosition, wxDefaultSize, 0);
        staticSizer11->Add( tempName, 0, wxALL, 5 );
        staticSizer11->Add( tempVal, 0, wxALL, 5 );

        // ---------------------------------------------
        // ----- probe size ----------------------------------------
        wxStaticBox* staticBox5 = new wxStaticBox( this, wxID_ANY, wxT("STEM probe size (pixels)") );
        wxStaticBoxSizer* staticSizer5 = new wxStaticBoxSizer( staticBox5, wxHORIZONTAL );
        boxSizer->Add( staticSizer5, 0,  wxALIGN_LEFT|wxALL, 5);

        wxStaticText* nxpName= new wxStaticText( this, wxID_STATIC, wxT("&nxp ="),
                wxDefaultPosition, wxDefaultSize, 0);
        nxpVal = new wxTextCtrl( this, ID_SAMPLE_NXP, wxT("00"),
                wxDefaultPosition, wxDefaultSize, 0);
        staticSizer5->Add( nxpName, 0, wxALL, 5 );
        staticSizer5->Add( nxpVal, 0, wxALL, 5 );

        wxStaticText* nypName= new wxStaticText( this, wxID_STATIC, wxT("&nyp ="),
                wxDefaultPosition, wxDefaultSize, 0);
        nypVal = new wxTextCtrl( this, ID_SAMPLE_NYP, wxT("00"),
                wxDefaultPosition, wxDefaultSize, 0);
        staticSizer5->Add( nypName, 0, wxALL, 5 );
        staticSizer5->Add( nypVal, 0, wxALL, 5 );

        // ---------------------------------------------
        // ----- STEM image size (nm)----------------------------------------
        wxStaticBox* staticBox6 = new wxStaticBox( this, wxID_ANY, wxT("STEM image size (in nm)") );
        wxStaticBoxSizer* staticSizer6 = new wxStaticBoxSizer( staticBox6, wxHORIZONTAL );
        boxSizer->Add( staticSizer6, 0,  wxALIGN_LEFT|wxALL, 5);

        wxStaticText* xiName= new wxStaticText( this, wxID_STATIC, wxT("&xinit ="),
                wxDefaultPosition, wxDefaultSize, 0);
        xiVal = new wxTextCtrl( this, ID_SAMPLE_XINIT, wxT("00"),
                wxDefaultPosition, wxDefaultSize, 0);
        staticSizer6->Add( xiName, 0, wxALL, 5 );
        staticSizer6->Add( xiVal, 0, wxALL, 5 );

        wxStaticText* xfName= new wxStaticText( this, wxID_STATIC, wxT("&xfinal ="),
                wxDefaultPosition, wxDefaultSize, 0);
        xfVal = new wxTextCtrl( this, ID_SAMPLE_XFINAL, wxT("00"),
                wxDefaultPosition, wxDefaultSize, 0);
        staticSizer6->Add( xfName, 0, wxALL, 5 );
        staticSizer6->Add( xfVal, 0, wxALL, 5 );

        wxStaticText* yiName= new wxStaticText( this, wxID_STATIC, wxT("&yinit ="),
                wxDefaultPosition, wxDefaultSize, 0);
        yiVal = new wxTextCtrl( this, ID_SAMPLE_YINIT, wxT("00"),
                wxDefaultPosition, wxDefaultSize, 0);
        staticSizer6->Add( yiName, 0, wxALL, 5 );
        staticSizer6->Add( yiVal, 0, wxALL, 5 );

        wxStaticText* yfName= new wxStaticText( this, wxID_STATIC, wxT("&yfinal ="),
                wxDefaultPosition, wxDefaultSize, 0);
        yfVal = new wxTextCtrl( this, ID_SAMPLE_YFINAL, wxT("00"),
                wxDefaultPosition, wxDefaultSize, 0);
        staticSizer6->Add( yfName, 0, wxALL, 5 );
        staticSizer6->Add( yfVal, 0, wxALL, 5 );

        // ---------------------------------------------
        // ----- STEM image size (pixels)----------------------------------------
        wxStaticBox* staticBox7 = new wxStaticBox( this, wxID_ANY, wxT("STEM image size (pixels)") );
        wxStaticBoxSizer* staticSizer7 = new wxStaticBoxSizer( staticBox7, wxHORIZONTAL );
        boxSizer->Add( staticSizer7, 0,  wxALIGN_LEFT|wxALL, 5);

        wxStaticText* nxiName= new wxStaticText( this, wxID_STATIC, wxT("&nxi ="),
                wxDefaultPosition, wxDefaultSize, 0);
        nxiVal = new wxTextCtrl( this, ID_SAMPLE_NXI, wxT("00"),
                wxDefaultPosition, wxDefaultSize, 0);
        staticSizer7->Add( nxiName, 0, wxALL, 5 );
        staticSizer7->Add( nxiVal, 0, wxALL, 5 );

        wxStaticText* nyiName= new wxStaticText( this, wxID_STATIC, wxT("&nyi ="),
                wxDefaultPosition, wxDefaultSize, 0);
        nyiVal = new wxTextCtrl( this, ID_SAMPLE_NYI, wxT("00"),
                wxDefaultPosition, wxDefaultSize, 0);
        staticSizer7->Add( nyiName, 0, wxALL, 5 );
        staticSizer7->Add( nyiVal, 0, wxALL, 5 );


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
 bool sampleDialog::TransferDataToWindow( )
 {
        wxString ws;

        nxVal->SetValue(ws.Format(wxT("%d"), nx) );
        nyVal->SetValue(ws.Format(wxT("%d"), ny) );

        nxpVal->SetValue(ws.Format(wxT("%d"), nxProbe) );
        nypVal->SetValue(ws.Format(wxT("%d"), nyProbe) );

        ncellxVal->SetValue(ws.Format(wxT("%d"), ncellx) );
        ncellyVal->SetValue(ws.Format(wxT("%d"), ncelly) );
        ncellzVal->SetValue(ws.Format(wxT("%d"), ncellz) );

        deltazVal->SetValue(ws.Format(wxT("%10.3f"), deltaz/10.) );

        nwobbleVal->SetValue(ws.Format(wxT("%d"), nwobble) );
        tempVal->SetValue(ws.Format(wxT("%10.3f"), temp) );

        xiVal->SetValue(ws.Format(wxT("%10.3f"), xi/10.) );
        xfVal->SetValue(ws.Format(wxT("%10.3f"), xf/10.) );
        yiVal->SetValue(ws.Format(wxT("%10.3f"), yi/10.) );
        yfVal->SetValue(ws.Format(wxT("%10.3f"), yf/10.) );

        nxiVal->SetValue(ws.Format(wxT("%d"), nxi) );
        nyiVal->SetValue(ws.Format(wxT("%d"), nyi) );

        return true;
 };

 //-------- TransferDataFromWindow() -----------------------
 bool sampleDialog::TransferDataFromWindow( )
 {
        long ii;
        double xx;

        nerror = 0;  // track number of read errors

        //  remember that GetValue() is from wxTextEntry which wxTextCrtl is derived from
        //
        if( nxVal->GetValue().ToLong(&ii) ) nx = (int) ii;
                else nerror++;
        if( nyVal->GetValue().ToLong(&ii) ) ny = (int) ii;
                else nerror++;

        if( nxpVal->GetValue().ToLong(&ii) ) nxProbe = (int) ii;
                else nerror++;
        if( nypVal->GetValue().ToLong(&ii) ) nyProbe = (int) ii;
                else nerror++;

        if( ncellxVal->GetValue().ToLong(&ii) ) ncellx = (int) ii;
                else nerror++;
        if( ncellyVal->GetValue().ToLong(&ii) ) ncelly = (int) ii;
                else nerror++;
        if( ncellzVal->GetValue().ToLong(&ii) ) ncellz = (int) ii;
                else nerror++;

        if( deltazVal->GetValue().ToDouble(&xx) ) deltaz = xx*10.0;
                else nerror++;

        if( nwobbleVal->GetValue().ToLong(&ii) ) nwobble = (int) ii;
                else nerror++;
        if( tempVal->GetValue().ToDouble(&xx) ) temp = xx;
                else nerror++;

        if( xiVal->GetValue().ToDouble(&xx) ) xi = xx*10.0;
                else nerror++;
        if( xfVal->GetValue().ToDouble(&xx) ) xf = xx*10.0;
                else nerror++;
        if( yiVal->GetValue().ToDouble(&xx) ) yi = xx*10.0;
                else nerror++;
        if( yfVal->GetValue().ToDouble(&xx) ) yf = xx*10.0;
                else nerror++;

        if( nxiVal->GetValue().ToLong(&ii) ) nxi = (int) ii;
                else nerror++;
        if( nyiVal->GetValue().ToLong(&ii) ) nyi = (int) ii;
                else nerror++;

        return true;
 };


 
