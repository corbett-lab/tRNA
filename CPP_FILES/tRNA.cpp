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
#include <iterator>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
using namespace std ;

/// declare our gsl
const gsl_rng *rng ;

/// our headers
#include "cmd_line.h"
#include "trna.h"
#include "individual.h"
#include "initialize_population.h"
#include "assign_sequence.h"
#include "nonlocal.h"
#include "local.h"
#include "segdup.h"
#include "mutate.h"
#include "gaussian_fitness.h"
#include "reproduce.h" 
#include "stats.h"


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

    /////////// for now, setting these to follow distributions of min/max:
    // germline: baseline is 1.45e-8 (narismahan et al), active tRNAs ~10x higher, current range 1e-8 - 1e-5
    // somatic: should be ~10x higher than germline (milholland et al). range for somatic will be between 1 - 100x higher than germline.
    // deletion: no idea but maybe give a 4 order of magnitude range in narrowing down
    // duplication: no idea but maybe give a 4 order of magnitude range in narrowing down
    /////////////
    // bergman sees ~1 / mil years in drosophila. we have ~90 species-specific tRNAs in humans over 7 million years (maybe overestimate)
    // this is a proxy for deletion and duplication rates though
    // how many duplications did we lose quickly? 
    // look up some more studies on estimates and see what we can get

    /// OLD :load in all bit score penalties from separate file:
    // std::ifstream is(options.path + "allPenaltiesPct.txt") ;
    // std::istream_iterator<double> start(is), end ;
    // std::vector<double> mutation_penalties(start, end) ;

    // load in different distributions of functions of tRNAs with given number of mutations:
    std::map<int, vector<double>> mutations_to_function ;
    for ( int i = 1 ; i < 11 ; i ++ ){
        std::string myNum = std::to_string(i) ;
        std::ifstream is(options.path + "functionDists/functionDists"+myNum+".txt") ;
        std::istream_iterator<double> start(is), end ;
        std::vector<double> mutation_penalties(start, end) ;
        mutations_to_function[i] = mutation_penalties ;
    }

    /// trna bank
    list<gene*> trna_bank ; 

    /// list of tRNA lifespans
    list<float> trna_lifespans ;
    
    // create population of size n with two tRNAs of equivalent function
    initialize_population ( options, trna_bank, trna_counter ) ; 

    /// now copy to population of size n
    vector<individual> population ( options.n ) ;
    for ( int i = 0 ; i < population.size() ; i ++ ) { 
    	for ( auto t : trna_bank ) { 
    		population[i].maternal_trnas.push_back( t ) ;
    		population[i].paternal_trnas.push_back( t ) ;
    	}
    }

    // fitness vector
    double fitness [options.n]  ;

    // optimal fitness based on inputs
    double opt_fit = (1 / sqrt( 2 * 3.14159265358979323846 * pow(options.fitness_sd, 2) )) ;
    
    // evolve the population forward in time
    for ( int g = 1 ; g <= options.generations ; g ++ ) {
        
        /// vector to swap with
        vector<individual> new_population ( options.n ) ;
                
        /// somatic and germline mutations
        mutate( population, options, trna_bank, g, trna_counter, mutations_to_function ) ;

        /// compute fitness
        compute_fitness( fitness, population, mutations_to_function, opt_fit, options ) ;
        
        /// reproduce w/ fitness + recombination
        reproduce( fitness, population, new_population, options ) ;

        /// swap populations
        swap( population, new_population ) ;

		/////// print stats
        // for first generation give user a quick reminder of what they did:
        if ( g == 1 ){
            cout << "somatic = " << options.somatic_rate ;
            cout << ", germline = " << options.germline_rate ;
            cout << ", dup = " << options.duplication_rate ;
            cout << ", del = " << options.deletion_rate << endl ;
        }

        print_stats( fitness, population, g, trna_bank, trna_lifespans, options ) ;
    }

    if ( options.output_frequencies == true ){
        std::string frequency_log = "frequencyLog" ;
        frequency_log += std::to_string(options.run_num) ;
        frequency_log += ".txt" ;
        ofstream myfile;
        myfile.open( frequency_log, fstream::app ) ;
        for ( auto t : trna_bank ) {
            myfile << t->name << "\t" << t->birth << "\t" << t->progenitor << "\t" ;
            for ( int i = 1 ; i < t->frequency.size() ; i ++ ) {
                myfile << t->frequency[i] << "," ;
            }
            myfile << endl ;
        }
        myfile.close();
    }

    printf("Total time: %.2f seconds. Total generations: ", (double)(clock() - tStart)/CLOCKS_PER_SEC) ;
    cout << options.generations << "." << endl ;
    return(0) ;
}

