/*
 
Need to rework a little for memmory efficiency
 >can stop storing knockout, just draw on the fly during fitness calcs
 >store all of the tRNAs w/ pointers and/or ints pointing towards vector positions
    > create tRNAs with an additional non-functional version one up on the list/vector
 
*/

/// standard headers
#include <string>
#include <string.h>
#include <time.h>
#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <math.h>
#include <map>
#include <set>
#include <utility>
#include <list>
#include <random>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
using namespace std ;

/// declare our gsl
const gsl_rng *rng ;

/// our headers
#include "cmd_line.h"
#include "Gene.h"
#include "Individual.h"
#include "initialize_population.h"
#include "assign_function.h"
#include "neighborhood.h"
#include "mutate.h"
#include "function_to_fitness.h"
#include "fitness.h"
#include "reproduce.h" 
#include "stats.h"
#include "Population.h"


/// main 
int main ( int argc, char **argv ) {

    //// for calculating runtime of code
    clock_t tStart = clock() ;

    /// trna counter to give new name to every trna that arises
    int trna_counter = 0 ;

    //// read command line options
    cmd_line options ;
    options.read_cmd_line( argc, argv ) ;

    // initialize rng for gsl lookup table
    rng = gsl_rng_alloc( gsl_rng_taus2 ) ;
    gsl_rng_set( rng, (long) options.seed ) ;
    // std::uniform_real_distribution<double> distribution(-7,-4);

    /// for now, setting these to follow distributions of min/max:
    // germline: baseline is 1.45e-8 (narismahan et al), active tRNAs ~10x higher, current range 1e-8 - 1e-5
    // somatic: should be ~10x higher than germline (milholland et al). range for somatic will be between 1 - 100x higher than germline.
    // deletion: 
    // duplication: 
    //
    // bergman sees ~1 / mil years in drosophila. we have ~90 species-specific tRNAs in humans over 7 million years (maybe overestimate)
    // this is a proxy for deletion and duplication rates though
    // how many duplications did we lose quickly? 
    // look up some more studies on estimates and see what we can get


    // TODO make this interface very easy for the end-user -- hide all these implementation details but expose the stuff
    // that the developer should need to specify beyond the config options

    /// trna bank
    list<Gene*> trna_bank ;

    // TODO lifespan should be inside the gene object?  or outside if gene gets reaped?
    /// list of tRNA lifespans
    list<float> trna_lifespans ;
    
    // TODO create population class to contain individuals
    // create population of size n with two tRNAs of equivalent function
    initialize_population ( options, trna_bank, trna_counter ) ; 

    /// now copy to population of size n
    Population population = Population(options.n);
    //vector<Individual> population(options.n);

    for ( int i = 0 ; i < population.getSize(); i ++ ) {
    //for ( int i = 0 ; i < population.size() ; i ++ ) {
    	for ( auto t : trna_bank ) { 
    		population.getIndividuals()[i].getMaternal_trnas().push_back( t ) ;
    		population.getIndividuals()[i].getPaternal_trnas().push_back( t ) ;
       		//population[i].getMaternal_trnas().push_back( t ) ;
       		//population[i].getPaternal_trnas().push_back( t ) ;

    	}
    }

    // TODO fitness should be member of individual?
    // fitness vector
    double fitness [options.n]  ;

    // TODO total function should be member of individual?
    double total_function_vector [options.n] ;
    
    // evolve the population forward in time

    // TODO: print header
    //cout << "\tt.first.name_t.first.locus_t.first.function_t.first.neighborhood_t.first.birth_.first.progenitor_t.second" ;

    for ( int g = 1 ; g < options.generations ; g ++ ) {
        
    	// TODO refactor to be a new population object
        /// vector to swap with
        Population new_population = Population(options.n);
    	//vector<Individual> new_population(options.n);


        // TODO: containment refactor -- population.evolve ==> individual.mutate() ==> gene.mutate() ==> individual.fitness() ==> individual.reproduce()
        /// somatic and germline mutations
        vector<Individual> individuals = population.getIndividuals();
    	vector<Individual> new_individuals = new_population.getIndividuals();

        mutate( individuals, options, trna_bank, g, trna_counter ) ;
    	//mutate( population, options, trna_bank, g, trna_counter ) ;

        /// compute fitness
        compute_fitness( total_function_vector, fitness, individuals, options ) ;
        //compute_fitness( total_function_vector, fitness, population, options ) ;

        /// reproduce w/ fitness + recombination
        reproduce( fitness, individuals, new_individuals, options ) ;
        //reproduce( fitness, population, new_population, options ) ;

        /// swap populations
        swap( individuals, new_individuals ) ;
        //swap( population, new_population ) ;

		/// print stats
        if ( g == 1 ){
            cout << "somatic = " << options.somatic_rate ;
            cout << ", germline = " << options.germline_rate ;
            cout << ", dup = " << options.duplication_rate ;
            cout << ", del = " << options.deletion_rate << endl ;
        }

        print_stats( fitness, individuals, g, trna_bank, trna_lifespans, options ) ;
        //print_stats( fitness, population, g, trna_bank, trna_lifespans, options ) ;

    }

    if ( options.output_frequencies == true ){
        std::string frequency_log = "frequencyLog" ;
        frequency_log += std::to_string(options.run_num) ;
        frequency_log += ".txt" ;
        ofstream myfile;
        myfile.open( frequency_log, fstream::app ) ;
        for ( auto t : trna_bank ) {
            myfile << t->getName() << "\t" << t->getBirth() << "\t" << t->getProgenitor() << "\t" ;
            for ( int i = 1 ; i < t->getFrequency().size() ; i ++ ) {
                myfile << t->getFrequency()[i] << "," ;
            }
            myfile << endl ;
        }
        myfile.close();
    }

    printf("Total time: %.2f seconds. Total generations: ", (double)(clock() - tStart)/CLOCKS_PER_SEC) ;
    cout << options.generations << "." << endl ;
    return(0) ;
}

