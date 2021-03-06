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
#include "reduce_ne.h"
#include "assign_genotype.h"
#include "duplication.h"
#include "mutate.h"
#include "fitness.h"
#include "gene_conversion.h"
#include "reproduce.h" 
#include "stats.h"
#include "sample.h"
#include "final_stats.h"
#include "transition.h"
#include "final_vectors.h"

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


    ///////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////// MUTATION PATHWAYS //////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////

    /**
    In this section:
    - if dual_rates == FALSE and mutation_pathways == FALSE:
    -- each SNP incurred by a tRNA gene will result in a change to its fitness
    -- the resulting fitness is drawn from the functionDists folder, where each numbered
    file represents the fitness of a yeast tRNA that has incurred that number of SNPs
    
    - if dual_rates == FALSE and mutation_pathways == TRUE:
    -- each SNP incurred by a tRNA gene will result in a change to its fitness, stored
    as its genotype in a string
    -- the fitness associated with that string is drawn from yeastGenotypePathways.txt,
    which contains the actual SNPs incurred by the tRNA genes in the same experiment.
    -- going forward, a tRNA can only get mutations for which there is a fitness associated
    with the resulting genotype.
    -- for example, if the first mutation is A1, from there it can only become A1-G4, 
    A1-A10, etc..

    - if dual_rates == true, we are using either the gamma model or the model4 model
    */

    // for random fitness drawing:
    std::map<int, vector<double>> mutations_to_function ;
    // for pathway-specific fitness drawing:
    std::map<string, double> genotype_to_fitness ;
    std::map<string, vector<string>> genotype_to_genotypes ;
    std::map<string, vector<double>> genotype_to_fitnesses ;

    // load in different distributions of functions of tRNAs with given number of mutations:
    if ( options.dual_rates == false ){
        if ( options.mutation_pathways == false ){
            for ( int i = 1 ; i < 11 ; ++i ){
                std::string myNum = std::to_string(i) ;
                std::ifstream is(options.path + "functionDists/functionDists"+myNum+".txt") ;
                std::istream_iterator<double> start(is), end ;
                std::vector<double> mutation_penalties(start, end) ;
                mutations_to_function[i] = mutation_penalties ;
            }
        }

        // assign each SNP to associated sequence score
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
    }

    ///////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////// DEMOGRAPHY FILE ////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////

    /**
    In this section, we read in a demography file.
    -- These files are formatted as such:
    node1 node2 branch_length ne1 ne2
    -- all branches making up a tree must be present
    -- there must also be a root as node1 for one branch
    -- before the root is the burn-in phase
    -- the branch lengths and effective population sizes may be scaled if the scaling_factor 
    flag is used.
    **/

    std::map<int, int> branch_to_length ;
    std::map<int, string> branch_to_node1 ;
    std::map<int, string> branch_to_node2 ;
    std::map<string, int> node_to_Ne ;
    std::map<string,vector<double>> node_to_final_active_loci ;
    std::map<string,vector<double>> node_to_final_inactive_loci ;
    std::map<string,map<string,int>> node_to_final_genotypes ;
    std::map<string, vector<individual>> node_to_population ;
    std::map<string, vector<gene*>> node_to_trna_bank ;
    int last_branch = 0 ;
    // add burn-in step, with Ne val from the root
    branch_to_length[0] = options.burn_in ;
    branch_to_node1[0] = "burn_in" ;
    branch_to_node2[0] = "root" ;
    if ( options.demography != "" ) {
        ifstream file( options.demography , ios::in );
        string node1, node2 ;
        int branchLength, node1Ne, node2Ne ;
        while ( file >> node1 >> node2 >> branchLength >> node1Ne >> node2Ne ) {
            last_branch ++ ;
            branch_to_length[last_branch] = branchLength / options.scaling_factor ;
            branch_to_node1[last_branch] = node1 ;
            branch_to_node2[last_branch] = node2 ;
            node_to_Ne[node1] = node1Ne / options.scaling_factor ;
            node_to_Ne[node2] = node2Ne / options.scaling_factor ;
            // cout << last_branch << "\t" << branchLength << "\t" << node1 << "\t" << node2 << endl ;
            // cout << last_branch << "\t" << branch_to_length[last_branch] << "\t" << node_to_Ne[node1] << "\t" << node_to_Ne[node2] << endl ;
        }
        file.close();
    }

    ///////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////// SIMULATION //////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////

    /**
    If the sampling flag is used, after X generations, an XXX_sample.txt file will be written
    with a sampling from the population at that specific time.
    **/

    // enable sampling 
    if ( options.sample == true ){
        std::string sampling_out = std::to_string(options.run_num) + "_sample.txt" ;
        ofstream stream( sampling_out ) ;
        stream << "" ;
    }

    // optimal fitness based on inputs (for use only if fitness is gaussian)
    double opt_fit = (1 / sqrt( 2 * 3.14159265358979323846 * pow(options.fitness_sd, 2) )) ;

    /////////////////////////////////////////
    /////////////////////////////////////////
    /////// NOT USING DEMOGRAPHY FILE ///////
    /////////////////////////////////////////
    /////////////////////////////////////////

    /**
    If a demography file is not used, will simulate one continuous population, with
    no changes to the population size, for the input number of generations. This is
    the most straightforward use of the simulation framework.

    The simulation main loop is as follows:
    -- mutate (go through all genes in all individuals and add mutations)
    -- gene_conversion (at user input rate, allow gene conversion within individuals)
    -- compute_fitness (induce somatic mutations, and, based on user input function, 
    compute fitness according to one of several functions)
    -- reproduce (using fitness including somatic mutations, as well as Haldane's map 
    function for recombination events, create next generation)
    -- swap (built-in function to save new population)
    -- print_stats (give the user )
    **/

    if ( options.demography == "" ) {
        // evolve the population forward in time
        // create population of size n with two tRNAs of equivalent function
        vector<gene*> trna_bank ; 
        initialize_population ( options, trna_bank, trna_counter ) ; 
        /// map of tRNA lifespans to count number of tRNAs that lived that long
        std::map<double, int> loci_to_lifespans ;
        std::map<int, int> lifespan_to_count ;

        // fitness vector
        double fitness [options.n]  ;

        /// now copy to population of size n
        vector<individual> population ( options.n ) ;
        for ( int i = 0 ; i < population.size() ; ++i ) { 
            for ( auto t : trna_bank ) { 
                population[i].maternal_trnas.push_back( t ) ;
                population[i].paternal_trnas.push_back( t ) ;
            }
        }

        /// run simulation for g generationså
        for ( int g = 1 ; g <= options.generations ; g ++ ) {

            // for first generation give user a quick reminder of what they did:
            if ( g == 1 ){
                cout << "somatic = " << options.somatic_rate << ", germline = " << options.germline_rate << ", sdup = " << options.duplication_rate << ", del = "  ;
                cout << options.deletion_rate << ", fitness function = " << options.fitness_func << ", gene conversion rate = " << options.gene_conversion_rate << endl ;
            }
            
            /// vector to swap with
            vector<individual> new_population ( options.n ) ;

            mutate( population, options, trna_bank, g, trna_counter, mutations_to_function, genotype_to_fitness, genotype_to_genotypes, genotype_to_fitnesses ) ;
            gene_conversion ( population, options, trna_bank, g, trna_counter ) ;
            compute_fitness( fitness, population, mutations_to_function, genotype_to_fitness, genotype_to_genotypes, genotype_to_fitnesses, opt_fit, options ) ;
            reproduce( fitness, population, new_population, options ) ;
            swap( population, new_population ) ;
            print_stats( fitness, population, g, options.generations, "default", trna_bank, loci_to_lifespans, lifespan_to_count, options ) ;
            // cout << "TEST: " << options.sample_all << "\t" << g << "\t" << options.burn_in << "\t" << g % options.sampling_frequency << endl ;

            if (( options.sample_all == true ) and ( g >= options.burn_in ) and ( g % options.sampling_frequency == 0 )){
                std::string sampling_out = std::to_string(options.run_num) + "_sample_all.txt" ;
                sample_all( g, population, sampling_out, options ) ;
            }

            else if (( options.sample == true ) and ( g >= options.burn_in ) and ( g % options.sampling_frequency == 0 )){
                sample_individuals( g, population, options ) ;
            }
        }
    }

    /////////////////////////////////////
    /////////////////////////////////////
    /////// USING DEMOGRAPHY FILE ///////
    /////////////////////////////////////
    /////////////////////////////////////

    /**
    If a demography file IS used, we need to simulate each branch, and then store 
    the population at the end of each branch (so that if there are multiple branches
    coming from a single node, we can simulate both from the same starting point).
    For the burn-in, we start from default settings.
    **/

    else {
        for ( int branch = 0 ; branch <= last_branch ; branch ++ ) {

            vector<individual> population ( node_to_Ne[branch_to_node2[branch]] ) ;
            double fitness [node_to_Ne[branch_to_node2[branch]]] ;
            vector<gene*> trna_bank ; 

            /// map of tRNA lifespans to count number of tRNAs that lived that long
            std::map<double, int> loci_to_lifespans ;
            std::map<int, int> lifespan_to_count ;

            // cout << "BRANCH: " << branch << ", length: " << branch_to_length[branch] << ", Ne: " << node_to_Ne[branch_to_node2[branch]] ;
            // cout << ", node1: " << branch_to_node1[branch] << ", node2: " << branch_to_node2[branch] << endl ;
            
            // if no saved population, start from default settings
            if ( !node_to_population.count( branch_to_node1[branch] ) ) {
                initialize_population ( options, trna_bank, trna_counter ) ; 
                for ( int i = 0 ; i < population.size() ; ++i ) { 
                    for ( auto t : trna_bank ) { 
                        population[i].maternal_trnas.push_back( t ) ;
                        population[i].paternal_trnas.push_back( t ) ;
                    }
                }
            }

            // if population we already built, take it and then subsample to get one fitting the new ne
            else {
                // cout << "\n\n" << population.size() << "\n\n" ;
                std::map<gene*, gene*> old_trna_to_new_trna ;
                reduce_ne( options, node_to_population[branch_to_node1[branch]], population, node_to_trna_bank[branch_to_node1[branch]], trna_bank, old_trna_to_new_trna ) ;
            }
            cout << "length: " << branch_to_length[branch] << endl ;

            // evolve the population forward in time
            for ( int g = 1 ; g <= branch_to_length[branch] ; g ++ ) {
                
                // for first generation give user a quick reminder of what they did:
                if ( g == 1 ){
                    cout << "somatic = " << options.somatic_rate << ", germline = " << options.germline_rate << ", sdup = " << options.duplication_rate << ", del = "  ;
                    cout << options.deletion_rate << ", fitness function = " << options.fitness_func << ", gene conversion rate = " << options.gene_conversion_rate << endl ;
                }
                
                /// vector to swap with
                vector<individual> new_population ( population.size() ) ;

                mutate( population, options, trna_bank, g, trna_counter, mutations_to_function, genotype_to_fitness, genotype_to_genotypes, genotype_to_fitnesses ) ;
                gene_conversion ( population, options, trna_bank, g, trna_counter ) ;
                compute_fitness( fitness, population, mutations_to_function, genotype_to_fitness, genotype_to_genotypes, genotype_to_fitnesses, opt_fit, options ) ;
                reproduce( fitness, population, new_population, options ) ;
                swap( population, new_population ) ;
                print_stats( fitness, population, g, branch_to_length[branch], branch_to_node2[branch], trna_bank, loci_to_lifespans, lifespan_to_count, options ) ;

                if (( options.sample_all == true ) and ( g % options.sampling_frequency == 0 )){
                    std::string sampling_out = std::to_string(options.run_num) + "_" + branch_to_node1[branch] + "_to_" + branch_to_node2[branch] + "_sample_all.txt" ;
                    sample_all( g, population, sampling_out, options ) ;
                }
                else if (( options.sample == true ) and ( g % options.sampling_frequency == 0 )){
                    sample_individuals( g, population, options ) ;
                }
            }
            node_to_population[branch_to_node2[branch]] = population ;
            cout << branch_to_node2[branch] << ", " << population.size() << endl ;
            node_to_trna_bank[branch_to_node2[branch]] = trna_bank ;
            // std::string sampling_out = std::to_string(options.run_num) + "_" + branch_to_node1[branch] + "_to_" + branch_to_node2[branch] + "_final_population.txt" ;
            // get_final_stats( population, sampling_out, options ) ;
            update_found( population, branch_to_node2[branch], node_to_final_active_loci, node_to_final_inactive_loci, node_to_final_genotypes, options ) ;
        }
        cout << "running final vectors" << endl ;
        std::string vector_out = std::to_string(options.run_num) + "_final_vector.txt" ;
        final_vectors( node_to_final_active_loci, node_to_final_inactive_loci, node_to_final_genotypes, node_to_Ne, vector_out, options ) ;   
    }

    printf("Total time: %.2f seconds. ", (double)(clock() - tStart)/CLOCKS_PER_SEC) ;
    if ( options.demography == "" ){
        cout << "Total generations: " << options.generations << "." ;
    }
    cout << "\n" ;
    return(0) ;
}

