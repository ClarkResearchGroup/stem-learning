/*   computem.hpp

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

  header file for computem.cpp
  temsim GUI using wxWidgets

-------------------------------------------------

  started separate header file 05-jan-2013 ejk ejk
  add help/gettingStarted 9-nov-2013 ejk
  add "do not type here" top message in olg window 11-nov-2013 ejk
  update "get started" message 1-dec-2013, 7-dec-2013 ejk
  update about/date 9-feb-2014 ejk
  update about info 9-jul-2017 ejk
  last updated 9-jul-2017  ejk
*/

//
#include "wx/wx.h"  // regular headers
//#include "wx/wxprec.h"  // for precompiled headers

#include "wx/mdi.h"     // MDI headers
#include "wx/docmdi.h"  // MDI headers for doc/view mode

#ifndef COMPUTEM_HPP_INCLUDED	// only include this file if its not already
#define COMPUTEM_HPP_INCLUDED	// remember that this has been included

#include "ctdoc.hpp"  //  header for this file
#include "ctview.hpp"    //  my view class

enum{
        ID_GET_START = 2501
};

// ---------------------------------------------------------------------------
//--------- Parent Frame Class -----------------
//
class MyFrame : public wxDocMDIParentFrame
{
        DECLARE_CLASS(MyFrame)
        DECLARE_EVENT_TABLE() // this class handles events
public:
        // constructor
        MyFrame( wxDocManager *manager, wxFrame *frame, wxWindowID id,
                const wxString& title, const wxPoint& pos, const wxSize& size, long type);

        // Event handlers
        void OnQuit( wxCommandEvent& event );
        void OnAbout( wxCommandEvent& event );
        void OnGetStart( wxCommandEvent& event );
        //void OnOpen( wxCommandEvent& event);
        wxTextCtrl *textWindow;

private:

        wxLogWindow *myLogWind;

};  // end MyFrame()


// ---------------------------------------------------------------------------
// MyApp
// ---------------------------------------------------------------------------
// --------- Declare the application -----------------
class MyApp : public wxApp
{
public:
        MyApp();

        virtual bool OnInit();
        virtual int OnExit();

        MyFrame* GetTopFrame()  { return topFrame; }

private:
        wxDocManager *m_docManager;
        MyFrame* topFrame;

};  // end MyApp()

DECLARE_APP( MyApp )   // Implements MyApp& GetApp()


#ifdef FOOF
// ---------------------------------------------------------------------------
//--------- Declare our child class -----------------
class MyChild: public wxDocMDIChildFrame
{
public:

    MyChild(wxMDIParentFrame *parent, const wxString& title);
    ~MyChild();

    void OnActivate(wxActivateEvent& event);
    void OnQuit(wxCommandEvent& event);
    void OnClose(wxCloseEvent& event);
    MyCanvas *canvas;

private:
        // this class handles events
        DECLARE_EVENT_TABLE()
};
#endif

#endif

