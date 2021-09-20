
autostem_cuda (in directory source) uses an Nvidia graphics card (GPU) installed in
a computer. It has only been tested on a core i7 computer with Ubuntu 18.04 (Linux),
cuda toolkit 10.1 and a RTX 2080 card. Other operating systems and hardware may or
may not work.  

Most Linux distributions are different enough so that programs cannot be easily
distributed in executable form and must be compiled on the machine in use. To use the
cuda based program the Nvidia cuda driver and the cuda toolkit must also be installed
(obtain separately from Nvidia) in addition to fftw etc. Then (after renaming the
file makefile.ubuntu to makefile), type "make autostem_cuda" from the command line
to compile the program before running it.

6-oct-2019