##  A MAKEFILE FOR tRNA GENE FAMILY SIMULATOR. TO USE, RUN THE COMMAND "make" VIA COMMAND LINE ##
export PATH := bin:$(PATH)
TCFLAGS = -ltcmalloc 
GSL_LIBS = $(shell /usr/local/bin/gsl-config --libs)
CFLAGS = $(shell /usr/local/bin/gsl-config --cflags)
GFLAG = -g 

all:
## IF TCMALLOC IS NOT INSTALLED, COMMENT OUT THE FOLLOWING LINE 
## BY PLACING A # IN FRONT OF THE FOLLOWING LINE AND REMOVE '#' ON THE NEXT LINE
#	$(LINK.cc) -std=c++11  tRNA.cpp $(GFLAG) $(TCFLAGS) $(GSL_LIBS) $(CFLAGS) -o tRNA
	#$(LINK.cc) -O0 -std=c++11 Gene.cpp Individual.cpp Population.cpp CommandLine.cpp Simulation.cpp Driver.cpp $(GFLAG) $(TCFLAGS) $(GSL_LIBS) $(CFLAGS) -o JCTEST
	$(LINK.cc) -O3 -std=c++11 Gene.cpp Individual.cpp Population.cpp CommandLine.cpp Simulation.cpp Driver.cpp $(GSL_LIBS) $(CFLAGS) -o JCTEST


## In many cases, tcmalloc will substantially decrease runtime of SELAM
## we therefore strongly recommend using this library 
## on mac osx it can be installed using homebrew: brew install google-perftools
## on ubuntu distributions you may try: sudo apt-get install google-perftools libgoogle-perftools-dev
