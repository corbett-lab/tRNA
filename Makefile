##  A MAKEFILE FOR tRNA GENE FAMILY SIMULATOR. TO USE, RUN THE COMMAND "make" VIA COMMAND LINE ##

TCFLAGS = -ltcmalloc 
GSL_LIBS = $(shell gsl-config --libs)
CFLAGS = $(shell gsl-config --cflags) 

all:
## IF TCMALLOC IS NOT INSTALLED, COMMENT OUT THE FOLLOWING LINE 
## BY PLACING A # IN FRONT OF THE FOLLOWING LINE AND REMOVE '#' ON THE NEXT LINE
	$(LINK.cc) -std=c++11 -O3 tRNA.cpp $(TCFLAGS) $(GSL_LIBS) $(CFLAGS) -o tRNA
#	$(LINK.cc) -std=c++11 -O3 tRNA.cpp $(GSL_LIBS) $(CFLAGS) -o tRNA


## In many cases, tcmalloc will substantially decrease runtime
## we therefore strongly recommend using this library 
## on mac osx it can be installed using homebrew: brew install google-perftools
## on ubuntu distributions you may try: sudo apt-get install google-perftools libgoogle-perftools-dev
