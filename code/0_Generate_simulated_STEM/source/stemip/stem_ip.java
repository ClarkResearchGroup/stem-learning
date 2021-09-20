/******************************************************************************
              stem_ip.java
               
  ImageJ plugin to perform image deconvolution for ADF STEM images
  
  This MUST be in ImageJ\Plugins to run (will compile in other dir. but NOT run)
 
------------------------------------------------------------------------
Copyright 2015-2017 Earl J. Kirkland

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

------------------------------------------------------------------------

references:

---Wiener Filter, and Adaptive Median Filter ----
1. R. C. Gonzalez and R. E. Woods, "Digital Image Processing",
	Prentice Hall, 2008 (3rd edit.)

---the RL method ---
2. W. H. Richardson, J. Opt. Soc. Amer. 62 (1972) p.55-59.
3. L. B. Lucy, The Astronomical Journal, 79 (1974) p.745-754.
4. L. A. Shepp and Y. Vardi, IEEE Trans. on Medical Imaging,
    1, (1982), p.113-122.
5. 	G. M. van Kempen and L. J. van Vliet and P. J. Verveer and
	H. T. M. van der Voort, Journal of Microscopy, 185, (1997), p.354-365.
 
   To Do:
     [1] h(-x) in RL
     [2] copy properties from input
     
  this source code is formatted for a TAB size of 4 char.

  started from jcstem.java with IP_Demo.java example 3-jul-2015 Earl J. Kirkland
  Wiener Filter working 17-jul-2015 ejk
  most of Richardson-Lucy working 18-jul-2015 ejk
  both working 21-jul-2015 ejk
  start adding adaptive median filter 3-dec-2015 ejk
  fix bug so AdaptiveMedian filter no longer requires power of 2 image size
  	20-apr-2017 ejk
  last modified 20-apr-2017 ejk ejk

******************************************************************************/

import java.awt.*;
import java.awt.event.*;
import java.io.*;
import java.util.StringTokenizer;  // for ReadParams() parse params

import ij.plugin.frame.*;
import ij.*;
import ij.process.*;
import ij.gui.*;

/////////////////////////////////////////////////////////////////////////

public class stem_ip extends PlugInFrame implements ActionListener
{
	String VersionDate = "20-apr-2017";
	
	private static final int sizeMax=7;  //  max size of adaptive median filter
	
	// --- parameters: df, Cs3, keV, ddf, beta, alpha, Cs5
	//    set initial default values here
	private static final int NPARAM = 14;
	public double param[] = { 100, 100.0, 100.0, 0.007, 10, 0.0, 50, 30.0, 0.5,
				0, 0, 0, 20.0, 15.0 };
	private double ADFmtfParam[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0  };
	private double wav, keV, fovx, fovy, Cs, Cs5, df, ddf, objApert, sigmae, dsource;
	private int iKEV, iFOVx, iFOVy, iCs3, iCs5, iDF, iDDF, iOBJAP, iADFsource,
			iASTIGa, iASTIGb, iBACK, iSNR, iNITER;
	private int iWIENER, iRL, iADAPTMEDIAN, NXmtf, NYmtf, Niter;
	
	private static final double pi=Math.PI, twopi=2.0*Math.PI;

	/*  absiccas and weights for Gauss-Hermite Quadrature 
	   with exp(-x*x) weighting of integrand 
	   from Abramowitz and Stegun, and Numerical Recipes */
	private static final int NGH=9;   /* number of Gauss-Hermete coeff. to use  */
	private static final double xGH[]={ 3.190993201781528, 2.266580584531843,
		1.468553289216668,
	    0.723551018752838, 0.000000000000000, -0.723551018752838,
	    -1.468553289216668,-2.266580584531843,-3.190993201781528};
	private static final double wGH[]={3.960697726326e-005, 4.943624275537e-003 ,
		8.847452739438e-002,
	    4.326515590026e-001, 7.202352156061e-001, 4.326515590026e-001,
	    8.847452739438e-002, 4.943624275537e-003, 3.960697726326e-005};

	
	//  input boxes for parameters
	private Choice METHODchoice;

	//  parameter descriptions etc.
	private TextField paramText, fovText, dfText, Cs3Text, kevText, ddfText,
		oaText, Cs5Text, ADFsourceText, astigText, backText, SNRText, NiterText;
	private Label titleLabel, fovLabel, dfLabel, Cs3Label, kevLabel, ddfLabel,
		oaLabel, Cs5Label, TEMLabel, ADFsourceLabel, astigLabel, backLabel,
		SNRLabel, NiterLabel;

	// action buttons
	private Button calcButton, aboutButton;
	String aboutString;
	
	// organize input boxes
	private Panel parPanel, topPanel;
	private GridBagConstraints grid;

	//-----------------   java stuff-----------------------
	private Panel panel;
	private int previousID;
	private static Frame instance;

	//-------  constructor ----------------------------------------
	public stem_ip() {
		super("stem_ip");
		if (instance!=null) {
			instance.toFront();
			return;
		}
		instance = this;
		addKeyListener(IJ.getInstance());

		setLayout( new GridBagLayout() );		// to manually position controls
		grid = new GridBagConstraints();
		grid.fill = GridBagConstraints.NONE;	// don't grow components
		grid.anchor = GridBagConstraints.WEST;	// put on left
		grid.insets = new Insets(5,5,5,5);		// 5 pixel margin on all sides
		
		//---- put title+ in a panel to keep them together
		topPanel = new Panel();
		topPanel.setLayout( new GridBagLayout() );
		
		titleLabel = new Label( "ADF STEM image restoration");
		calcButton = new Button( "run" );
		calcButton.addActionListener( this );
		aboutButton = new Button( "about" );
		aboutButton.addActionListener( this );
		
		grid.gridheight = 1;
		grid.gridy = 0;
		grid.weightx = 1;
		grid.weighty = 1;
		grid.gridwidth  = 4;
		grid.gridx = 0;
		topPanel.add( titleLabel, grid );
		grid.gridwidth  = 1;
		grid.gridx = 4;
		topPanel.add( calcButton, grid );
		grid.gridwidth  = 1;
		grid.gridx = 5;
		topPanel.add( aboutButton, grid );
		
		grid.gridheight = 1;
		grid.gridwidth  = 4;
		grid.gridx = 0;
		grid.gridy = 0;
		grid.weightx = 0;
		grid.weighty = 0;
		add( topPanel, grid );
		
		//--------- about box ---------------------------------------
		aboutString = new String(
               " stem_ip, "+VersionDate+", ejk\n"
			+ " \n"
            + "THIS PROGRAM IS SUPPLIED AS-IS\n"
			+ "WITH NO WARRANTY OR GUARANTEE\n"
			+ " \n"
			+ "Restore (deconvolve) a high \n"
			+ "resolution scanning transmission \n"
			+ "electron microscope (ADF-STEM) \n"
			+ "image modeled as a very simple \n"
			+ "incoherent image approx. valid for \n"
			+ "thin specimens. The image parameters \n"
			+ "need to be accurate for good  \n"
			+ "results. The adaptive median filter\n"
			+ "just removes impulsive noise. \n"
			+ "See Kirkland, Acta Cryst, A72, 2016.");
		
		//--------- optical parameters ---------------------------------------
		METHODchoice = new Choice();
		METHODchoice.addItem( "Wiener filter" );
		iWIENER = 0;
		METHODchoice.addItem( "Richardson-Lucy" );
		iRL = 1;
		METHODchoice.addItem( "Adapt. Median" );
		iADAPTMEDIAN = 2;
		METHODchoice.addItemListener( new ItemListener() {
			public void itemStateChanged( ItemEvent e) {}
		});
		TEMLabel = new Label( "parameters" );
		
		//  parameter index values
		iKEV = 0;
		iFOVx = 1;
		iFOVy = 2;
		iCs3 = 3;
		iCs5 = 4;
		iDF = 5;
		iDDF = 6;
		iOBJAP = 7;
		iADFsource = 8;
		iASTIGa = 9;
		iASTIGb = 10;
		iBACK = 11;
		iSNR = 12;
		iNITER = 13;
		
		NXmtf = 0;	//  to remember what has been calculated
		NYmtf = 0;
		
		kevLabel = new Label( "keV" );
		kevText = new TextField(6);
		kevText.setText( Double.toString( param[iKEV] ) );
		kevText.addActionListener( this );
				
		fovLabel = new Label( "FOV x,y (in Ang)" );
		fovText = new TextField(6);
		fovText.setText( Double.toString( param[iFOVx] ) + " " +
			Double.toString( param[iFOVy] ));
		fovText.addActionListener( this );

		Cs3Label = new Label( "Cs (in mm)" );
		Cs3Text = new TextField(6);
		Cs3Text.setText( Double.toString( param[iCs3] ) );
		Cs3Text.addActionListener( this );

		Cs5Label = new Label( "Cs5 (in mm)" );
		Cs5Text = new TextField(6);
		Cs5Text.setText( Double.toString( param[iCs5] ) );
		Cs5Text.addActionListener( this );

		dfLabel = new Label( "defocus (in Ang)" );
		dfText = new TextField(6);
		dfText.setText( Double.toString( param[iDF] ) );
		dfText.addActionListener( this );

		ddfLabel = new Label( "defocus spread (A)" );
		ddfText = new TextField(6);
		ddfText.setText( Double.toString( param[iDDF] ) );
		ddfText.addActionListener( this );

		oaLabel = new Label( "obj. angle (mrad)" );
		oaText = new TextField(6);
		oaText.setText( Double.toString( param[iOBJAP] ) );
		oaText.addActionListener( this );

		ADFsourceLabel = new Label( "ADF source (in A)" );
		ADFsourceText = new TextField(6);
		ADFsourceText.setText( Double.toString( param[iADFsource] ) );
		ADFsourceText.addActionListener( this );
				
		astigLabel = new Label( "astig a,b (in A)" );
		astigText = new TextField(6);
		astigText.setText( Double.toString( param[iASTIGa] ) + " " +
			Double.toString( param[iASTIGb] ));
		astigText.addActionListener( this );
				
		backLabel = new Label( "background" );
		backText = new TextField(6);
		backText.setText( Double.toString( param[iBACK] ) );;
		backText.addActionListener( this );

		SNRLabel = new Label( "SNR (W)" );
		SNRText = new TextField(6);
		SNRText.setText( Double.toString( param[iSNR] ) );
		SNRText.addActionListener( this );

		NiterLabel = new Label( "N iter (RL)" );
		NiterText = new TextField(6);
		NiterText.setText( Double.toString( param[iNITER] ) );
		NiterText.addActionListener( this );

		//---- put all parameters in a panel to keep them together
		parPanel = new Panel();
		parPanel.setLayout( new GridLayout( 0, 2, 5, 5) );
	
		parPanel.add(METHODchoice);	// make a menu-like thing
		parPanel.add(TEMLabel);
		parPanel.add(kevLabel);		//  parameters
		parPanel.add(kevText);
		parPanel.add(fovLabel);
		parPanel.add(fovText);
		parPanel.add(Cs3Label);
		parPanel.add(Cs3Text);
		parPanel.add(Cs5Label);
		parPanel.add(Cs5Text);
		parPanel.add(dfLabel);
		parPanel.add(dfText);
		parPanel.add(ddfLabel);
		parPanel.add(ddfText);
		parPanel.add(oaLabel);
		parPanel.add(oaText);
		parPanel.add(ADFsourceLabel);
		parPanel.add(ADFsourceText);
		parPanel.add(astigLabel);
		parPanel.add(astigText);
		parPanel.add(backLabel);
		parPanel.add(backText);
		
		parPanel.add(SNRLabel);
		parPanel.add(SNRText);
		parPanel.add(NiterLabel);
		parPanel.add(NiterText);
		
		grid.gridheight = 14;
		grid.gridwidth  = 2;
		grid.gridx = 0;
		grid.gridy = 1;
		grid.weightx = 0;
		grid.weighty = 0;
		add( parPanel, grid );
		
		pack();
		GUI.center(this);
		setVisible(true);
		
	}  //  end stem_ip() constructor

	//-------------------------------------------------------------
	//       action
	//-------------------------------------------------------------
	
	public void actionPerformed(ActionEvent e) {
		
		//------------------ about button ---------------------------
		if( e.getSource() == aboutButton ) {
			IJ.showMessage("about stem_ip ",aboutString);
			return;
		}

		//--------- verify that there is an image to work on --------------
		ImagePlus imp = WindowManager.getCurrentImage();
		if (imp==null) {
			IJ.beep();
			IJ.showStatus("No image");
			previousID = 0;
			return;
		}
		
		int nx = imp.getWidth();
		int ny = imp.getHeight();
		int nx2 = Fft.powerof2( nx );
		int ny2 = Fft.powerof2( ny );
		if(  (METHODchoice.getSelectedIndex() != iADAPTMEDIAN  )
			&&	( ( (nx2!=nx) || (ny2!=ny) ) ) ) {
			IJ.beep();
			IJ.showStatus("Nx, Ny must be a power of 2");
			previousID = 0;
			imp.unlock();
			return;
		}
		int nstack = imp.getStackSize();
		if( nstack > 1  ) {
			IJ.beep();
			IJ.showStatus("cannot handle image stacks");
			previousID = 0;
			imp.unlock();
			return;
		}
		
		if (!imp.lock())
			{previousID = 0; return;}
		
		//----  make sure this is a greyscale pix ------------------------
		int it = imp.getType();
		if( !( (ImagePlus.GRAY8==it) || (ImagePlus.GRAY16==it) || (ImagePlus.GRAY32==it) ) ){
			IJ.beep();
			IJ.showStatus("stem_ip cannot handle RGB pix");
			previousID = 0;
			imp.unlock();
			return;
		}
			
		//-------  get parameter values from UI and do simple verification ------------
		ReadParams( );
		if( (param[iKEV]<=0) ||(param[iFOVx]<=0) || (param[iFOVy]<=0) || (param[iOBJAP]<=0)
				|| (param[iADFsource]<0) || (param[iSNR]<0) || (param[iNITER]<0) ){
			IJ.beep();
			IJ.showStatus("parameters seem wrong");
			previousID = 0;
			imp.unlock();
			return;
		}
		
		// model after example IP_DEMO.java			
		int id = imp.getID();
		if (id!=previousID)
			imp.getProcessor().snapshot();
		previousID = id;
		String label = e.getActionCommand();
		if (label==null)
			return;
		new Runner(label, imp);
	}

	public void processWindowEvent(WindowEvent e) {
		super.processWindowEvent(e);
		if (e.getID()==WindowEvent.WINDOW_CLOSING) {
			instance = null;	
		}
	}

	//-------------------------------------------------------------
	//       run action
	//-------------------------------------------------------------
	class Runner extends Thread { // inner class
		private String command;
		private ImagePlus imp;
	
		Runner(String command, ImagePlus imp) {
			super(command);
			this.command = command;
			this.imp = imp;
			setPriority(Math.max(getPriority()-2, MIN_PRIORITY));
			start();
		}
	
		public void run() {
			try {
				runCommand(command, imp);
			} catch(OutOfMemoryError e) {
				IJ.outOfMemory(command);
				if (imp!=null) imp.unlock();
			} catch(Exception e) {
				CharArrayWriter caw = new CharArrayWriter();
				PrintWriter pw = new PrintWriter(caw);
				e.printStackTrace(pw);
				IJ.log(caw.toString());
				IJ.showStatus("");
				if (imp!=null) imp.unlock();
			}
		}
	
		void runCommand(String command, ImagePlus imp) {
			ImageProcessor ip = imp.getProcessor();
			IJ.showStatus(command + "...");
			long startTime = System.currentTimeMillis();
			String title;

			int nx = imp.getWidth();
			int ny = imp.getHeight();
			//IJ.showStatus("size= "+(nx)+" "+ny);   //  ???? debugging
			
			//  calculate the MTF
			float mtfr[][]  = new float[nx][ny]; //  real, imag part of MTF
			float mtfi[][]  = new float[nx][ny]; //  doesn't really need to be complex(?)
			if(  METHODchoice.getSelectedIndex() != iADAPTMEDIAN  ) 
				makeADFmtf( mtfr, mtfi, nx, ny );
			
			float pix[][] = ip.getFloatArray();   //  get input image as a float (new copy)
			imp.unlock();  // don't need to access original any more
			
			//---- do the real deconvolution here ---------
			if( METHODchoice.getSelectedIndex() == iWIENER ) {

				WienerFilter( pix, mtfr, mtfi, nx, ny );
				title = "stem ip Wiener filter";
				
			} else if(  METHODchoice.getSelectedIndex() == iRL  ) {
				
				RichLucy( pix, mtfr, mtfi, nx, ny );
				title = "stem ip Richardson-Lucy";
				
			} else if(  METHODchoice.getSelectedIndex() == iADAPTMEDIAN  ) {
				
				//  use mtfr[][] as scratch array
				AdaptMedian( pix, mtfr, nx, ny, sizeMax );
				title = "stem ip Adaptive Median";
				for( int ix=0; ix<nx; ix++)
				for( int iy=0; iy<ny; iy++) pix[ix][iy] = mtfr[ix][iy];
				
			} else title = "bad mode choice";			
			
			double pmin, pmax, pt;
			pmin = pmax = pix[0][0];
			for( int ix=0; ix<nx; ix++)
			for( int iy=0; iy<ny; iy++) {
				pt = pix[ix][iy];
				if( pt > pmax ) pmax = pt;
				if( pt < pmin ) pmin = pt;
			}
			
			// make new 32 bit float image (leave original unchanged for comparison)			
			ImagePlus imp2 = NewImage.createFloatImage(title,
						nx, ny, 1, NewImage.FILL_RAMP );
			ImageProcessor ip2 = imp2.getProcessor();
			
			//  store the new image and set greyscale range
			ip2.setFloatArray( pix );
			imp2.setDisplayRange( pmin, pmax );
			imp2.updateAndRepaintWindow();	//  put results in new window
			imp2.show();
			IJ.showStatus((System.currentTimeMillis()-startTime)+" milliseconds");
		}
	
	} // Runner inner class

	
	//-------------------------------------------------------------
	//       read all parameters manually (not from an action)
	//-------------------------------------------------------------
	public int ReadParams( )
	{
		Double d;
		Integer i;
		int n;

		try{
			d = Double.valueOf( kevText.getText() );
			param[iKEV] =  d.doubleValue();
	
			StringTokenizer stf = new StringTokenizer( fovText.getText(), " " );
        	n = stf.countTokens();
        	if( n != 2 ) {
	        	param[iFOVx] = param[iFOVy] = 0.0;
				astigText.setText( Double.toString( param[iFOVx] ) + " " +
					Double.toString( param[iFOVy] ));
			} else {
				d = Double.valueOf( stf.nextToken() );
				param[iFOVx] =  d.doubleValue();
				d = Double.valueOf( stf.nextToken() );
				param[iFOVy] =  d.doubleValue();
			}
	
			d = Double.valueOf( Cs3Text.getText() );
			param[iCs3] =  d.doubleValue();
	
			d = Double.valueOf( Cs5Text.getText() );
			param[iCs5] =  d.doubleValue();
	
			d = Double.valueOf( dfText.getText() );		// optical parameters
			param[iDF] =  d.doubleValue();
	
			d = Double.valueOf( ddfText.getText() );
			param[iDDF] =  d.doubleValue();
	
			d = Double.valueOf( oaText.getText() );
			param[iOBJAP] =  d.doubleValue();
	
			StringTokenizer sta = new StringTokenizer( astigText.getText(), " " );
        	n = sta.countTokens();
        	if( n != 2 ) {
	        	param[iASTIGa] = param[iASTIGb] = 0.0;
				astigText.setText( Double.toString( param[iASTIGa] ) + " " +
					Double.toString( param[iASTIGb] ));
			} else {
				d = Double.valueOf( sta.nextToken() );
				param[iASTIGa] =  d.doubleValue();
				d = Double.valueOf( sta.nextToken() );
				param[iASTIGb] =  d.doubleValue();
			}
	
			d = Double.valueOf( ADFsourceText.getText() );
			param[iADFsource] =  d.doubleValue();	        	
	
			d = Double.valueOf( backText.getText() );
			param[iBACK] =  d.doubleValue();	        	
			
			d = Double.valueOf( SNRText.getText() );
			param[iSNR] =  d.doubleValue();
	
			d = Double.valueOf( NiterText.getText() );
			param[iNITER] =  d.doubleValue();
		}
		
		catch( NumberFormatException e ) { return( -1 ); }
		
		return( +1 );
			
	}  // end ReadParams()
	
	
	//--------------------------------------------------------
	//  ADF STEM subroutine group (below)
	//--------------------------------------------------------
	/*------------- makeADFmtf() -------------------------------
	
	calculate new ADF-STEM probe intensity and
	approx. the MTF in a simple incoherent model
	
	leave as full 2D calculation so the multipole aberrations
	can be added at some time if needed
	  
	kx,ky = real scattering vectors
	
	global variables used:
	  	ADFmtfParam[] = parameter array from previous calculation
	  	param[] = parameter array
	  	pi = constant
	
	*/
	private void makeADFmtf( float ADFmtfr[][], float ADFmtfi[][], int nx, int ny )
	{
		int i, ix, iy, ixc, iyc, idf, ndf;
		boolean changed;
		double k2max, mm0;
		double k, dk, k2, mtf, w1, w2, w5, wr, dx, dy, sum, ds;
		double ddf, ddf2, df0, astiga, astigb, astig;
		float ax, by;
		float  tr, ti, hr, hi, ds2, weight;
		float x[]   = new float[nx];
		float kx[]  = new float[nx];
		float kx2[] = new float[nx];
		float y[]   = new float[ny];
		float ky[]  = new float[ny];
		float ky2[] = new float[ny];
		float ctf[][]  = new float[nx][ny]; 
		
		//----  if parameters have not changed then 
		changed = false;
		if( (nx != NXmtf) && (ny != NYmtf) ) changed = true;
		NXmtf = nx;
		NYmtf = ny;
		for( i=0; i<NPARAM; i++) {  // scan all but last 2 = SNR and RL-iter
			if( (Math.abs( param[i] - ADFmtfParam[i] ) > 
				1.0e-4 * Math.abs( param[i] + ADFmtfParam[i] )) && i<(NPARAM-2) )
					changed = true;
			ADFmtfParam[i] = param[i];
		}

		//??? force calc. every time -> if( !changed ) return;
		// recalculate image for new parameters
		IJ.showStatus( "Calculate ADF mtf size= "+(nx)+" "+ny );
		
		// spatial frequencies  (FOVx,y = field of view = ax, by
		ax = (float) param[iFOVx];
		by = (float) param[iFOVy];
		freqn( kx, kx2, x, nx, (float) param[iFOVx] );
		freqn( ky, ky2, y, ny, (float) param[iFOVy] );
		
		keV = ADFmtfParam[iKEV];
		//wav = 12.3986/Math.sqrt( (2*511.0+keV)*keV ); // wavelength
		wav = wavelength( keV ); // wavelength
		mm0 = 1.0F + keV/511.0F;		 // relativistic mass ratio
		sigmae = sigma( keV )/ 1000.0;   // global
		
		k2max = 0.001 * ADFmtfParam[iOBJAP]/wav;		//  obj. apert. 
		k2max = k2max * k2max;
				
		//--- misc constants for transfer function		
		Cs = ADFmtfParam[iCs3] * 1.0e7;	// Cs3, convert mm to Ang.
		w5 = wav * wav * wav;
		w1 = 0.5 * pi * Cs * w5;
		df0 = ADFmtfParam[iDF];			// defocus in Ang.
		astiga = pi * wav * ADFmtfParam[iASTIGa];
		astigb = 2.0 * pi * wav * ADFmtfParam[iASTIGb];
		Cs5 = ADFmtfParam[iCs5] * 1.0e7;	// convert mm to Ang.
		w5 = pi * Cs5 * wav * wav * w5 /3.0;

		ddf = ADFmtfParam[iDDF];
	    if( ddf > 1.0 ) {       //  df integration parameters
	        ndf = NGH;
	        /* ddf2 = Math.sqrt(Math.log(2.0)/(ddf*ddf));   convert from HWHM */
	        ddf2 = Math.sqrt(Math.log(2.0)/(ddf*ddf/4.0));  /* convert from FWHM */
	    } else {
	        ndf = 1;
	        ddf2 = 0.0;
	    }

	    for( ix=0; ix<nx; ix++)
	    for( iy=0; iy<ny; iy++) ctf[ix][iy] = 0;
    
	    //---- integrate over defocus spread if ddf is large enough ----
	    // use Gauss-Hermite quadrature and convert exp(-a^2df^2) to exp(-u^2)  
	    for( idf=0; idf<ndf; idf++) {
	
			if( ndf > 1 ){
				df = df0 + xGH[idf]/ddf2;
				weight = (float) wGH[idf];
			}  else {
				df = df0;
				weight = 1;
			}
			
			w2 = pi * wav * df;
			
			//--- calculate the aperture function ---
			for( ix=0; ix<nx; ix++ )
			for( iy=0; iy<ny; iy++ ) {
				k2 = kx2[ix] + ky2[iy];
				if( k2 <= k2max ) {
					astig = astiga*( kx2[ix] - ky2[iy] ) + astigb*kx[ix]*ky[iy];
					wr = ( (w5*k2 + w1)*k2 - w2 )*k2 + astig;
					ADFmtfr[ix][iy] = (float) ( Math.cos( wr ) );	// transfer function
					ADFmtfi[ix][iy] = (float) (-Math.sin( wr ) );
				} else {
					ADFmtfr[ix][iy] = 0;
					ADFmtfi[ix][iy] = 0;
				}
			} // end for(iy...)
			
			//--- inv. FT and square magnitude to get psf ----		
			Fft.fft2d( ADFmtfr, ADFmtfi, nx, ny, -1 );
			
			for( ix=0; ix<nx; ix++ )
			for( iy=0; iy<ny; iy++ ) {
	            tr = ADFmtfr[ix][iy];
	            ti = ADFmtfi[ix][iy];
	            ctf[ix][iy] += weight* (tr*tr + ti*ti);
	        }
			
		}
		sum = 0;
		for( ix=0; ix<nx; ix++ )
		for( iy=0; iy<ny; iy++ ) sum += ctf[ix][iy];

		tr = (float) (1.0 / (sum) );
		for( ix=0; ix<nx; ix++ )		// normalize to unity integrated itensity
		for( iy=0; iy<ny; iy++ ) {
			ADFmtfr[ix][iy] = tr * ctf[ix][iy];  // probe intensity
			ADFmtfi[ix][iy] = 0;
		}
		// -- this is real so could be done faster - fix sometime ????
		Fft.fft2d( ADFmtfr, ADFmtfi, nx, ny, +1 ); // fwd trans. to get MTF
		
		// ---  add the source size
		dsource = 0.5 * ADFmtfParam[iADFsource];  // convert diameter to radius
		ds = pi*pi * dsource*dsource/Math.log(2);	// source size factor- convert to FWHM
		for( ix=0; ix<nx; ix++ )
		for( iy=0; iy<ny; iy++ ) {
			k2 = kx2[ix] + ky2[iy];
			if( k2 <= 4.0*k2max ) {
				ds2 = (float) Math.exp( -ds*k2 ) ;
				ADFmtfr[ix][iy] *= ds2;
				ADFmtfr[ix][iy] *= ds2;
			} else {
				ADFmtfr[ix][iy] = 0;	// outside the aperture
				ADFmtfr[ix][iy] = 0;
			}
		}
		
		//  leave as FFT for transfer function

	} // end makeADFmtf()
	
	/*------------- WienerFilter() -------------------------------
	
	  perform Wiener Filter
	  
	  arguments:
	  	float pix[][] = real valued input pix in real space
	  					will be overwritten on exit
	  	snr = SNR to use
	  	float ADFmtrr,i [][] = real,imag MTF in FT space
	  	
	  global variables used:
	  		ADFmtfParam[]
	
	*/
	private void WienerFilter( float pix[][],
		float ADFmtfr[][], float ADFmtfi[][],  int nx, int ny )
	{
		int ix, iy;
		float hr, hi, gr, gi, h2, eta, snr;
				
		snr = (float) ADFmtfParam[iSNR];
		
		float pixr[][]  = new float[nx][ny]; //  real, imag part of MTF
		float pixi[][]  = new float[nx][ny]; //  doesn't really need to be complex(?)
		
		for( ix=0; ix<nx; ix++)
		for( iy=0; iy<ny; iy++){
			pixr[ix][iy] = pix[ix][iy];
			pixi[ix][iy] = 0;
		}
		
		Fft.fft2d( pixr, pixi, nx, ny, +1 );
		
		eta = 1.0F/snr;
		for( ix=0; ix<nx; ix++)
		for( iy=0; iy<ny; iy++){
			hr = ADFmtfr[ix][iy];	// transfer function
			hi = ADFmtfi[ix][iy];
			h2 = eta + hr*hr + hi*hi;
			gr = pixr[ix][iy];
			gi = pixi[ix][iy];
			pixr[ix][iy] = ( gr*hr + gi*hi)/h2;	// H*G the long way
			pixi[ix][iy] = (-gr*hi + gi*hr)/h2;
		}
		
		Fft.fft2d( pixr, pixi, nx, ny, -1 );
		
		for( ix=0; ix<nx; ix++)
		for( iy=0; iy<ny; iy++){ pix[ix][iy] = pixr[ix][iy]; }
		
		return;
		
	}  //  end WienerFilter()
	
	/*------------- RichLucy() -------------------------------
	
	  perform Richardson-Lucy
	  
	  I should figure out how to do this with a real to complex FFT sometime
	  and see if it runs faster
	  
	  arguments:
	  	float pix[][] = real valued input pix in real space
	  					will be overwritten on exit

	  	float ADFmtrr,i [][] = real,imag MTF in FT space
	  	nx,ny = size of image
	  	
	  global variables used:
	  	ADFmtfParam[]
	
	*/
	private void RichLucy( float pix[][],
		float ADFmtfr[][], float ADFmtfi[][],  int nx, int ny )
	{
		int ix, iy, it, niter;
		float hr, hi, gr, gi, fr, fi, h2, back;
		double kev, wav, k2max, mtf, ko, k2, pref;
		
		float x[]   = new float[nx];		// for prefilter
		float kx[]  = new float[nx];
		float kx2[] = new float[nx];
		float y[]   = new float[ny];
		float ky[]  = new float[ny];
		float ky2[] = new float[ny];
		
		float gpr[][]  = new float[nx][ny]; //  real, imag part of MTF
		float gpi[][]  = new float[nx][ny]; //  doesn't really need to be complex(?)
		float fpr[][]  = new float[nx][ny]; //  real, imag part of pix
		float fpi[][]  = new float[nx][ny];
		float tpr[][]  = new float[nx][ny]; //  real, imag part of pix
		float tpi[][]  = new float[nx][ny];
		
		niter = (int) (ADFmtfParam[iNITER] + 0.5);  // number of iterations
		back = (float) ADFmtfParam[iBACK];			// background to subtract
		
		freqn( kx, kx2, x, nx, (float) param[iFOVx] );
		freqn( ky, ky2, y, ny, (float) param[iFOVy] );
		
		for( ix=0; ix<nx; ix++)		// copy original to scratch complex pix
		for( iy=0; iy<ny; iy++){
			gpr[ix][iy] = pix[ix][iy];
			gpi[ix][iy] = 0;
		}
		
		//-----  apply a prefilter to get rid of high freq. noise that cannot be real
		//   and include in full transfer function to be correct
		keV = ADFmtfParam[iKEV];
		wav = wavelength( keV ); 					// wavelength in Ang.
		k2max = 0.001 * ADFmtfParam[iOBJAP]/wav;	//  obj. apert. size in 1/Ang 
		k2max = 2.00 * k2max;  // go a little further out for max feq.
		k2max = k2max * k2max;
		pref = 0.50;  // prefilter range - arbitrary value adjust to give good results

		Fft.fft2d( gpr, gpi, nx, ny, +1 );
		for( ix=0; ix<nx; ix++) {
			ko = kx2[ix];
			for( iy=0; iy<ny; iy++) {
				k2 = ko + ky2[iy];				
				mtf = k2/k2max;  // Gaussian
				h2 = (float) Math.exp( -pref*mtf );
				gpr[ix][iy] *= h2;
				gpi[ix][iy] *= h2;
			}
		}
		Fft.fft2d( gpr, gpi, nx, ny, -1 );
		//-- end prefilter
		
		//-- copy to start value and subtract background
		for( ix=0; ix<nx; ix++)
		for( iy=0; iy<ny; iy++){
			hr = gpr[ix][iy] - back;
			if( hr < 0 ) hr = 0;      //  enforce positivity
			gpr[ix][iy] = hr;
			gpi[ix][iy] = 0;
			fpr[ix][iy] = hr;
			fpi[ix][iy] = 0;
		}

		//----- main RL loop --------------------
		for( it=0; it<niter; it++) {
			
			IJ.showStatus("RL iteration "+(it+1)); 
						
			// 1. convolve current pix with h(x) into temp pix t
			for( ix=0; ix<nx; ix++)
			for( iy=0; iy<ny; iy++){
				tpr[ix][iy] = fpr[ix][iy];
				tpi[ix][iy] = fpi[ix][iy];
			}
			Fft.fft2d( tpr, tpi, nx, ny, +1 );
			for( ix=0; ix<nx; ix++) {
				ko = kx2[ix];
				for( iy=0; iy<ny; iy++){
					k2 = ko + ky2[iy];
					mtf = k2/k2max;
					h2 = (float) Math.exp( -pref*mtf );
					// include prefilter
					hr = ADFmtfr[ix][iy] * h2;
					hi = ADFmtfi[ix][iy] * h2;
					fr = tpr[ix][iy];
					fi = tpi[ix][iy];
					tpr[ix][iy] = ( hr*fr - hi*fi);
					tpi[ix][iy] = ( hi*fr + hr*fi);
				}
			}
			Fft.fft2d( tpr, tpi, nx, ny, -1 );
			
			// 2. divide original by current 
			for( ix=0; ix<nx; ix++)
			for( iy=0; iy<ny; iy++){
				tpr[ix][iy] = gpr[ix][iy]/tpr[ix][iy];
				tpi[ix][iy] = 0;
			}
			
			// 3. convolve with h(-x) (mult. by c.c. H)
			Fft.fft2d( tpr, tpi, nx, ny, +1 );
			for( ix=0; ix<nx; ix++) {
				ko = kx2[ix];
				for( iy=0; iy<ny; iy++){
					k2 = ko + ky2[iy];
					mtf = k2/k2max;
					h2 = (float) Math.exp( -pref*mtf );
					hr = ADFmtfr[ix][iy] * h2;
					hi = ADFmtfi[ix][iy] * h2;
					fr = tpr[ix][iy];
					fi = tpi[ix][iy];
					tpr[ix][iy] = ( hr*fr + hi*fi);
					tpi[ix][iy] = (-hi*fr + hr*fi);
				}
			}
			Fft.fft2d( tpr, tpi, nx, ny, -1 );
			
			// 4. update current pix
			for( ix=0; ix<nx; ix++)
			for( iy=0; iy<ny; iy++){
				fpr[ix][iy] *= tpr[ix][iy];
				tpi[ix][iy] = 0;
			}
			
		}  // end for(it...)
		
		//------  store results ---------
		for( ix=0; ix<nx; ix++)
		for( iy=0; iy<ny; iy++){
			pix[ix][iy] = fpr[ix][iy];
		}
		return;
		
	}  //  end RichLucy()
	
	
	//--------------------------------------------------------
	//  from slicelib.java below  (just a few needed so don't need whole class)
	//--------------------------------------------------------
		
	/*------------------------ freqn() ------------------------*/
	/*
		Calculate spatial frequencies for use with fft's
		NOTE: zero freg is in the bottom left corner and
			expands into all other corners - not in the center
			this is required for fft - don't waste time rearanging
	
		This routine must be called once for each direction
	
		ko[n]  = real array to get spatial frequencies
		ko2[n] = real array to get k[i]*k[i]
		xo[n]  = real array to get positions 
		nk     = integer number of pixels
		ak     = real full scale size of image in pixels
	*/
	public static void freqn( float ko[], float ko2[], float xo[], int nk, float ak )
	{
		int i, imid;
	
		imid = nk/2;
	
		for( i=0; i<nk; i++) {
			xo[i] = ((float) (i * ak) ) / ((float)(nk-1));
			if ( i > imid ) {
				ko[i]  = ((float)(i-nk)) / ((float)ak);
			} else {
				ko[i]  = ((float)i) / ((float)ak);
			}
			ko2[i] = ko[i] * ko[i];
		}
	
	}  /*  end freqn() */
	
	
	/*--------------------- sigma() -----------------------------------*/
	/*
		return the interaction parameter sigma in radians/(kv-Angstroms)
		keep this is one place so I don't have to keep typing in these
		constants (that I can never remember anyhow)
	
		ref: Physics Vade Mecum, 2nd edit, edit. H. L. Anderson
			(The American Institute of Physics, New York) 1989
			page 4.
	
		kev = electron energy in keV
	*/
	public static double sigma( double kev )
	{
		double s, pi, wavl, x;
		final double emass=510.99906; /* electron rest mass in keV */
	
		x = ( emass + kev ) / ( 2.0*emass + kev);
		wavl = wavelength( kev );
		pi = 4.0 * Math.atan( 1.0 );
		
		s = 2.0 * pi * x / (wavl*kev);
		
		return( s );
	
	}  /* end sigma() */
	
	/*--------------------- wavelength() -----------------------------------*/
	/*
		return the electron wavelength (in Angstroms)
		keep this is one place so I don't have to keep typing in these
		constants (that I can never remember anyhow)
	
		ref: Physics Vade Mecum, 2nd edit, edit. H. L. Anderson
			(The American Institute of Physics, New York) 1989
			page 4.
	
		kev = electron energy in keV
	*/
	public static double wavelength( double kev )
	{
		double w;
		final double emass=510.99906; /* electron rest mass in keV */
		final double hc=12.3984244; /* Planck's const x speed of light*/
	
		/* electron wavelength in Angstroms */
		w = hc/Math.sqrt( kev * ( 2*emass + kev ) );
	
		return( w );
	
	}  // end wavelength()
	
	//--------------------------------------------------------
	//  some general image processing
	//--------------------------------------------------------
		
	/*------------- AdaptMedian() -------------------------------
	
	  perform Adaptive Median Filter
	  
	  arguments:
	  	float pix[][] = real valued input pix in real space
	  					will be overwritten on exit
	  	snr = SNR to use
	  	float ADFmtrr,i [][] = real,imag MTF in FT space
	  	
	  global variables used:
	  		ADFmtfParam[]
	
	*/
	private void AdaptMedian( float inpix[][], float rpix[][], int nx, int ny, int sizeMax )
	{
		int ix, iy, nx1,ny1, nx2,ny2, ix2, iy2;
		int wsize, wso2, n, doStageB;
	
		float zmin, zmax, zmed, zxy, znew;
	
		float p[] = new float[ sizeMax*sizeMax ];
	
		// test for impossible situation
		if( (nx < sizeMax) || ( ny < sizeMax) || ( sizeMax < 3 ) ) return;
	
		//  loop over all pixels in the image
		for( ix=0; ix<nx; ix++) for( iy=0; iy<ny; iy++) {
	
			wsize = 3;  // initial window size
	
			//--- stage A ------
			do{ 
				//p.clear();
				doStageB = 0;
				znew = 0;
	
				// make a window of size wsize * wsize 
				//   offset to inside pix on edges
				wso2 = wsize/2;
				nx1 = ix - wso2;	//  small corner of window
				ny1 = iy - wso2;
				if( nx1 < 0 ) nx1 = 0;
				if( ny1 < 0 ) ny1 = 0;
	
				nx2 = nx1 + wsize;	//  large corner of window
				ny2 = ny1 + wsize;
				if( nx2 > nx ) { nx2 = nx;     nx1 = nx2 - wsize; }
				if( ny2 > ny ) { ny2 = ny;     ny1 = ny2 - wsize; }
	
				n = 0;
				for( ix2=nx1; ix2<nx2; ix2++) for( iy2=ny1; iy2<ny2; iy2++) {
					p[n] = inpix[ix2][iy2];
					n++;
				}
				shellSort( p, n );
	
				zmax = p[0];	//  max value
				zmin = p[n-1];	//  min value
				zmed = p[n/2];  //  median value
				if( (zmed > zmin) && (zmed < zmax) ) {
					doStageB = 1;
					break;
				}
	
				wsize++;
				znew = zmed;
	
			} while ( wsize < sizeMax );
	
			//---- stage B ------
			if( doStageB != 0 ) {
				zxy = rpix[ix][iy];
				if( (zxy > zmin) && (zxy < zmax ) ) znew = zxy;
				else znew = zmed;
			}
	
			rpix[ix][iy] = znew;
	
		} //  end for(ix... )
	
		return;
			
	}  //  end AdaptMedian()

	/*----------------- shellSort() ------------------------------
	
	    Shell sort modeled after prog. 6.5 (pg. 274) of
	    R. Sedgewick, "Algorithms in C", 3rd edit. Addison-Wesley 1998
	    
	    z[]		= numbers to sort
	    n       = number of atoms
	*/
	void shellSort( float z[], int n )
	{
	    int i, j, h;
	    float z2;
	
	    for( h=1; h<=(n-1)/9; h=3*h+1);
	
	    for( ; h>0; h/=3 )
	        for( i=h; i<n; i++) {
	            j = i;
	            z2 = z[i];
	            while( (j >= h) && ( z2 < z[j-h]) ) {
	                z[j] = z[j-h];
	                j -= h;
	            }
	            z[j] = z2;
	        }
	 
	}  // end shellSort()


} // end stem_ip class
