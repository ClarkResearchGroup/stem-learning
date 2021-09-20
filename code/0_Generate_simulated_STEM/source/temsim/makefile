#
# makefile to generate TEMSIM multislice package under ubuntu.
# assumes that FFTW already installed (see www.fftw.org)
#
# Put this file in the same directory as the TEMSIM
# C source files and type "make all" from a command line,
# to compile all of the programs.
#
#  type "make all" to compile everthing
#  type "make remove" to remove all of the compiled files.
#
# Each file name is assumed to be all lower case.
#
# convert to C++ 29-may-2012 ejk
# update for **cmd.cpp 7-may-2014 ejk
# change options for incostem - -std=c++11 
#       for random 10-feb-2015 ejk
# update autoslic to include probe.cpp 24-apr-2017 ejk
# add pdb2xyz, ftiff2raw and remove sumpix,slicview 10-may-2018 ejk
# add autostem_cuda jun 2018 ejk
# add autoslic_cuda oct 2018 ejk
# last modified 24-oct-2018 ejk
#

# define compiler with optimize flag
#CC = gcc -O
CC = g++ -O3
#DEL = del  # windows/mingw - doesn't work without .exe in file name
DEL = rm  # unix

# define libraries
MYLIBS = slicelib.o floatTIFF.o cfpix.o 
LIBS = ${MYLIBS}$
WLIBS = slicelib.o floatTIFF.o cfpix.o -lfftw3f_threads -lfftw3f

#
#  entry point to build everything
#
.PHONY : all
all:
	$(MAKE) atompot
	$(MAKE) autoslic
	$(MAKE) autostem
	$(MAKE) ftiff2raw
	$(MAKE) image
	$(MAKE) incostem
	$(MAKE) mulslice
	$(MAKE) pdb2xyz
	$(MAKE) probe
	$(MAKE) stemslic
	$(MAKE) sumpix

#
#  entry point to remove compiled files
#
.PHONY : remove
remove:
	${DEL}$ atompot
	${DEL}$ autoslic
	${DEL}$ autostem
	${DEL}$ ftiff2raw
	${DEL}$ image
	${DEL}$ incostem
	${DEL}$ mulslice
	${DEL}$ pdb2xyz
	${DEL}$ probe
	${DEL}$ stemslic
	${DEL}$ sumpix
	${DEL}$ cfpix.o
	${DEL}$ slicelib.o
	${DEL}$ floatTIFF.o

#
#  main programs
#

atompot: atompot.cpp  ${MYLIBS}
	${CC} -o atompot atompot.cpp ${WLIBS}

autoslic: autoslic.cpp autosliccmd.cpp ${MYLIBS}
	${CC} -fopenmp -o autoslic autosliccmd.cpp autoslic.cpp probe.cpp ${WLIBS}

autoslic_cuda: autoslic_cuda.cu autosliccmd_cuda.cpp cudaSlice.cu ${MYLIBS}
	nvcc -O3 -o autoslic_cuda autosliccmd_cuda.cpp autoslic_cuda.cu probe.cpp ${WLIBS} -arch=sm_61 -lcuda -lcufft

#  remember that cuda is not thread safe so cannot use use both cuda and openmp 
#     and nvcc has no -fopenmp option(!)
autostem_cuda: autostem_cuda.cu autostemcmd_cuda.cpp cudaSlice.cu ${MYLIBS}
#	nvcc -g -G -o autostem_cuda autostemcmd_cuda.cpp autostem_cuda.cu ${WLIBS} -arch=sm_61 -lcuda -lcufft 
	nvcc -O3 -o autostem_cuda autostemcmd_cuda.cpp autostem_cuda.cu ${WLIBS} -arch=sm_61 -lcuda -lcufft 

autostem: autostem.cpp autostemcmd.cpp ${MYLIBS}
	${CC} -fopenmp -o autostem autostemcmd.cpp autostem.cpp ${WLIBS} 

display: display.cpp  ${MYLIBS}
	${CC} -o display display.cpp ${WLIBS}

ftiff2raw: ftiff2raw.cpp  psave.cpp ${MYLIBS}
	${CC} -o ftiff2raw ftiff2raw.cpp psave.cpp ${WLIBS}

image: image.cpp  ${MYLIBS}
	${CC} -o image image.cpp ${WLIBS}

incostem: incostem.cpp incostemcmd.cpp ${MYLIBS}
	${CC} -std=c++11 -o incostem incostemcmd.cpp incostem.cpp probe.cpp ${WLIBS}

mulslice: mulslice.cpp ${MYLIBS}
	${CC} -o mulslice mulslice.cpp ${WLIBS}

pdb2xyz: pdb2xyz.cpp ${MYLIBS}
	${CC} -o pdb2xyz pdb2xyz.cpp ${WLIBS}

probe: probe.cpp probecmd.cpp ${MYLIBS}
	${CC} -o probe probecmd.cpp probe.cpp ${WLIBS}

stemslic: stemslic.cpp ${MYLIBS}
	${CC} -o stemslic stemslic.cpp ${WLIBS}$

sumpix: sumpix.cpp ${MYLIBS}
	${CC} -o sumpix sumpix.cpp ${WLIBS}

#
# define subroutine library
#

cfpix.o: cfpix.cpp
	${CC} -c cfpix.cpp

slicelib.o: slicelib.cpp
	${CC} -c slicelib.cpp

floatTIFF.o: floatTIFF.cpp
	${CC} -c floatTIFF.cpp

