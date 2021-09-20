#
# makefile to generate computem multislice package in Ubuntu 13.10 Linux
# most likely you will need to edit this on every different linux distro
#
# Put this file in the same directory as the computem
# C/C++ source files and type "make computem" from a command line,
# to compile
#
# start ubuntu version 14-oct-2013 ejk
# switch to ubuntu 13.10 with wxwidgets 2.9.5 
#    (2.6 and 2.8 won't work) on 18-oct-2013 ejk
# add psave.cpp 10-feb-2014 ejk
# switch to ubuntu 14.04 with wxwidgets 3.0.0
#    6-may-2014 ejk
# add rangeDialog.cpp 29-sep-2014 ejk
# switch to ubuntu 16.04 with wxwidgets 3.1.0
#    9-may-2016 ejk
#

# define libraries
MYLIBS =  wxslicelib.o wxfloatTIFF.o cfpix.o psave.o
CTLIBS = autoslic.o autostem.o incostem.o probe.o

CTSUBS =  ctdoc.cpp ctview.cpp 
CTDIAL = aberDialog.cpp generalDialog.cpp sampleDialog.cpp textDialog.cpp rangeDialog.cpp

#  will probably be different on every computer
#  from 
#  wx-config --cxxflags
#  wx-config --libs

CT_CXXFLAGS = -I/usr/local/lib/wx/include/gtk2-unicode-static-3.1 -I/usr/local/include/wx-3.1 -D_FILE_OFFSET_BITS=64 -D__WXGTK__ -pthread

CT_LIB = -L/usr/local/lib -pthread   /usr/local/lib/libwx_gtk2u_xrc-3.1.a /usr/local/lib/libwx_gtk2u_qa-3.1.a /usr/local/lib/libwx_baseu_net-3.1.a /usr/local/lib/libwx_gtk2u_html-3.1.a /usr/local/lib/libwx_gtk2u_adv-3.1.a /usr/local/lib/libwx_gtk2u_core-3.1.a /usr/local/lib/libwx_baseu_xml-3.1.a /usr/local/lib/libwx_baseu-3.1.a -lgthread-2.0 -pthread -lX11 -lSM -lgtk-x11-2.0 -lgdk-x11-2.0 -lpangocairo-1.0 -latk-1.0 -lcairo -lgdk_pixbuf-2.0 -lgio-2.0 -lpangoft2-1.0 -lpango-1.0 -lgobject-2.0 -lglib-2.0 -lfontconfig -lfreetype -lpng -lexpat -lwxregexu-3.1 -lwxtiff-3.1 -lwxjpeg-3.1 -lz -ldl -lm 

#####
CC = gcc -O
CPP = g++ -O
DEL = rm

#
#  entry point to remove compiled files
#
remove:
	${DEL}$ autoslic.o
	${DEL}$ autostem.o
	${DEL}$ cfpix.o
	${DEL}$ incostem.o
	${DEL}$ probe.o
	${DEL}$ psave.o
	${DEL}$ wxfloatTIFF.o
	${DEL}$ wxslicelib.o
#
#  main programs
#

computem: computem.cpp  ${CTSUBS} ${CTDIAL} ${MYLIBS} ${CTLIBS}
	${CPP} ${CT_CXXFLAGS} -o computem computem.cpp \
	${CTSUBS} ${CTDIAL} ${MYLIBS} ${CTLIBS} \
    	 ${CT_LIB} -lfftw3f -lfftw3f_threads


#
# define low level subroutine libraries
#
cfpix.o: cfpix.cpp
	${CPP} ${CT_CXXFLAGS}  -c cfpix.cpp

psave.o: psave.cpp
	${CPP} ${CT_CXXFLAGS}  -c psave.cpp

wxslicelib.o: wxslicelib.cpp
	${CPP} ${CT_CXXFLAGS} -c wxslicelib.cpp

wxfloatTIFF.o: wxfloatTIFF.cpp
	${CPP} ${CT_CXXFLAGS} -c wxfloatTIFF.cpp

#  for testing
sampleDialog.o: sampleDialog.cpp
	${CPP} ${CT_CXXFLAGS} -c sampleDialog.cpp

#
#  define high level subroutine
#
autoslic.o: autoslic.cpp
	${CPP} ${CT_CXXFLAGS} -c autoslic.cpp

autostem.o: autostem.cpp
	${CPP} ${CT_CXXFLAGS} -c autostem.cpp

incostem.o: incostem.cpp
	${CPP} ${CT_CXXFLAGS} -c incostem.cpp

probe.o: probe.cpp
	${CPP} ${CT_CXXFLAGS} -c probe.cpp

