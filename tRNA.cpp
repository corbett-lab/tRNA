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
#include "assign_genotype.h"
#include "nonlocal.h"
#include "local.h"
#include "segdup.h"
#include "mutate.h"
#include "gaussian_fitness.h"
#include "exp_fitness.h"
#include "model_fitness.h"
#include "redund_fitness.h"
#include "gene_conversion.h"
#include "reproduce.h" 
#include "stats.h"
#include "sample.h"
#include "sample_all.h"

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



    // for random fitness drawing:
    std::map<int, vector<double>> mutations_to_function ;
    // for pathway-specific fitness drawing:
    std::map<string, double> genotype_to_fitness ;
    std::map<string, vector<string>> genotype_to_genotypes ;
    std::map<string, vector<double>> genotype_to_fitnesses ;

    // load in different distributions of functions of tRNAs with given number of mutations:
    if ( options.mutation_pathways == false ){
        for ( int i = 1 ; i < 11 ; ++i ){
            std::string myNum = std::to_string(i) ;
            std::ifstream is(options.path + "functionDists/functionDists"+myNum+".txt") ;
            std::istream_iterator<double> start(is), end ;
            std::vector<double> mutation_penalties(start, end) ;
            mutations_to_function[i] = mutation_penalties ;
            // if ( i == 1 ){
            //     for ( int m = 1 ; m < mutation_penalties.size() ; ++m ){
            //         cout << i << "\t" << (mutations_to_function[i])[m] << endl ;
            //     }
            // }
        }
    }

    // load in mutation pathways from yeast data:
    else {
        ifstream input(options.path + "functionDists/yeastGenotypePathways.txt") ;
        char const row_delim = '\n' ;
        char const field_delim = '\t' ;
        char const entry_delim = ',' ;
        for (string row; getline(input, row, row_delim); ) {
            istringstream iss(row) ;
            vector<string> tokens{istream_iterator<string>{iss}, istream_iterator<string>{}} ;

            // get genotype
            string genotype ;
            std::istringstream line_stream(row) ;
            line_stream >> genotype ;

            // if mutation chain ends here, add fitness of genotype + 0 as only option after
            if ( tokens.size() == 2 ) {
                double fitness ;
                while (line_stream >> fitness) {
                     genotype_to_fitness[genotype] = fitness ;
                     vector<string> genotypes ;
                     genotypes.push_back( "x" );
                     vector<double> fitnesses ;
                     fitnesses.push_back( 0.0 ) ;
                     genotype_to_genotypes[genotype] = genotypes ;
                     genotype_to_fitnesses[genotype] = fitnesses ;
                 }
            }

            // if mutation chain keeps going, add index 1 as fitness
            // then create vectors of posssible options for future genotypes and fitnesses
            else if ( tokens.size() > 2 ){
                string temp_fitness = tokens.at( 1 ) ;
                double fitness = std::stod( temp_fitness ) ;
                genotype_to_fitness[genotype] = fitness ;

                vector<string> genotypes ;
                vector<double> fitnesses ;
                string temp_genotypes = tokens.at( 2 ) ;
                string temp_fitnesses = tokens.at( 3 ) ;
                int g_count = std::count(temp_genotypes.begin(), temp_genotypes.end(), ',') ;
                int f_count = std::count(temp_fitnesses.begin(), temp_fitnesses.end(), ',') ;
                // if only one option, create empty vector and add them in
                if ( g_count == 0 ){
                    genotypes.push_back( temp_genotypes ) ;
                    double temp_f = std::stod( temp_fitnesses ) ;
                    fitnesses.push_back( temp_f ) ;
                    genotype_to_genotypes[genotype] = genotypes ;
                    genotype_to_fitnesses[genotype] = fitnesses ;
                }
                // if multiple options, split strings by comma and add all elements to vectors
                else {
                    istringstream gss(temp_genotypes) ;
                    for (string temp_g; getline(gss, temp_g, entry_delim); ) {
                        genotypes.push_back( temp_g ) ;
                    }
                    istringstream fss(temp_fitnesses) ;
                    for (string temp_f; getline(fss, temp_f, entry_delim); ) {
                        double temp_f_each = std::stod( temp_f ) ;
                        fitnesses.push_back( temp_f_each ) ;
                    }
                    genotype_to_genotypes[genotype] = genotypes ;
                    genotype_to_fitnesses[genotype] = fitnesses ;
                }
            }
        }
    }


    if ( options.sample == true ){
        std::string sampling_out = std::to_string(options.run_num) + "_sample.txt" ;
        ofstream stream( sampling_out ) ;
        stream << "" ;
    }

    /// trna bank
    list<gene*> trna_bank ; 
    cout << "TRNA_BANK MEMORY ADDRESS: " << &trna_bank << endl ;

    /// reusable pointers
    list<gene*> reusable_pointers ; 

    /// map of tRNA lifespans to count number of tRNAs that lived that long
    std::map<double, int> loci_to_lifespans ;
    std::map<int, int> lifespan_to_count ;
    
    // create population of size n with two tRNAs of equivalent function
    initialize_population ( options, trna_bank, trna_counter ) ; 

    /// now copy to population of size n
    vector<individual> population ( options.n ) ;
    for ( int i = 0 ; i < population.size() ; ++i ) { 
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

        // for first generation give user a quick reminder of what they did:
        if ( g == 1 ){
            cout << "somatic = " << options.somatic_rate ;
            cout << ", germline = " << options.germline_rate ;
            cout << ", dup = " << options.duplication_rate ;
            cout << ", del = " << options.deletion_rate ;
            cout << ", fitness function = " << options.fitness_func ;
            if ( options.fitness_func == "exp" ){
                cout << ", lambda = " << options.fitness_lambda ;
            }
            else if ( options.fitness_func == "gaussian" ) {
                cout << ", mean = " << options.fitness_mean << ", sd = " << options.fitness_sd ;
            }
            else if ( options.fitness_func == "redundant" ) {
                cout << ", minimum = " << options.min_fitness ;
            }
            cout << ", gene conversion rate = " << options.gene_conversion_rate << endl ;
        }
        
        /// vector to swap with
        vector<individual> new_population ( options.n ) ;

        /// somatic and germline mutations
        mutate( population, options, trna_bank, g, trna_counter, mutations_to_function, genotype_to_fitness, genotype_to_genotypes, genotype_to_fitnesses ) ;

        /// compute fitness
        if ( ( options.fitness_func == "exp" ) or ( options.fitness_func == "exponential" ) ) {
            compute_exp_fitness( fitness, population, mutations_to_function, genotype_to_fitness, genotype_to_genotypes, genotype_to_fitnesses, options ) ;
        }
        else if ( options.fitness_func == "model" ) {
            compute_model_fitness( fitness, population, mutations_to_function, genotype_to_fitness, genotype_to_genotypes, genotype_to_fitnesses, options ) ;
        }
        else if ( options.fitness_func == "redundant" ) {
            compute_redund_fitness( fitness, population, mutations_to_function, genotype_to_fitness, genotype_to_genotypes, genotype_to_fitnesses, options ) ;
        }
        else {
            compute_gaussian_fitness( fitness, population, mutations_to_function, genotype_to_fitness, genotype_to_genotypes, genotype_to_fitnesses, opt_fit, options ) ;
        }

        /// non-allelic gene conversion
        gene_conversion ( population, options, trna_bank, g, trna_counter ) ;
        
        /// reproduce w/ fitness + recombination
        reproduce( fitness, population, new_population, options ) ;

        /// swap populations
        swap( population, new_population ) ;

		/////// print stats
        print_stats( fitness, population, g, trna_bank, loci_to_lifespans, lifespan_to_count, options ) ;

        if (( options.sample_all == true ) and ( g >= options.burn_in ) and ( g % options.sampling_frequency == 0 )){
            sample_all( g, population, options ) ;
        }

        else if (( options.sample == true ) and ( g >= options.burn_in ) and ( g % options.sampling_frequency == 0 )){
            sample_individuals( g, population, options ) ;
        }

    }

    printf("Total time: %.2f seconds. Total generations: ", (double)(clock() - tStart)/CLOCKS_PER_SEC) ;
    cout << options.generations << "." << endl ;
    return(0) ;
}

