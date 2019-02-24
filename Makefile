##  A MAKEFILE FOR tRNA GENE FAMILY SIMULATOR. TO USE, RUN THE COMMAND "make" VIA COMMAND LINE ##

TCFLAGS = -ltcmalloc 
GSL_LIBS = $(shell gsl-config --libs)
CFLAGS = $(shell gsl-config --cflags) 
GFLAGS = -g
OPTIMIZER_FLAGS=-O0

clean:
	rm -f tRNA

all:
## IF TCMALLOC IS NOT INSTALLED, COMMENT OUT THE FOLLOWING LINE 
## BY PLACING A # IN FRONT OF THE FOLLOWING LINE AND REMOVE '#' ON THE NEXT LINE
#	$(LINK.cc) -std=c++11 -O3 tRNA.cpp $(TCFLAGS) $(GSL_LIBS) $(CFLAGS) -o tRNA
#	$(LINK.cc) -std=c++11 -O3 population.cpp $(GFLAGS) $(GSL_LIBS) $(CFLAGS) -o population
	#$(LINK.cc) -std=c++11 -O3 gene.cpp tRNA.cpp $(GFLAGS) $(GSL_LIBS) $(CFLAGS) $(OPTIMIZER_FLAGS) -o tRNA
#	$(LINK.cc) -std=c++11 -O3 gene.cpp Individual.cpp tRNA.cpp $(GFLAGS) $(GSL_LIBS) $(CFLAGS) $(OPTIMIZER_FLAGS) -o tRNA
	$(LINK.cc) -std=c++11 -O3 gene.cpp Individual.cpp Population.cpp tRNA.cpp $(GFLAGS) $(GSL_LIBS) $(CFLAGS) $(OPTIMIZER_FLAGS) -o tRNA
#	$(LINK.cc) -std=c++11 -O3 Population.cpp tRNA.cpp $(GFLAGS) $(GSL_LIBS) $(CFLAGS) $(OPTIMIZER_FLAGS) -o tRNA
#	$(LINK.cc) -std=c++11 -O3 tRNA.cpp $(GFLAGS) $(GSL_LIBS) $(CFLAGS) $(OPTIMIZER_FLAGS) -o tRNA


## In many cases, tcmalloc will substantially decrease runtime of SELAM
## we therefore strongly recommend using this library 
## on mac osx it can be installed using homebrew: brew install google-perftools
## on ubuntu distributions you may try: sudo apt-get install google-perftools libgoogle-perftools-dev