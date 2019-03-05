/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   CommandLine.h
 * Author: jcasaletto
 *
 * Created on March 3, 2019, 2:52 PM
 */

#ifndef COMMANDLINE_H
#define COMMANDLINE_H

#include <string>
#include <string.h>

using namespace std;

class CommandLine {
public:
    CommandLine();
    CommandLine(const CommandLine& orig);
    virtual ~CommandLine();
        
     // mean proportion of cells in which a portion of the genome is active
    double mean_neighborhood ;
    
    /// population size, constant for now
    int n ;
    
    /// number of generations
    int generations ;
    
    /// number of genes at start
    int start_count ;
    
    /// mutation rates
    double germline_rate ;
    double somatic_rate ;
    
    /// duplicate and deletion rates
    double deletion_rate ;
    double duplication_rate ;

    // start with pseudogene or no
    bool pseudogene ;
    
    /// probability of cluster
    double prob_cluster ; 

    /// to be used in mapping sequence to fitness
    double lambda_seq ;

    // output frequencies at the end or no
    bool output_frequencies ;

    // output lifespans for all tRNAs
    bool output_lifespans ;

    /// rng seed
    int seed ;

    /// map length 
    double map_length ;

    // print every __ generations
    int print_count ;

    // for running many simulations and keeping track of output
    int run_num ;

    // number of generations to be used for burn_in
    // won't start printing results until after this generation
    int burn_in ;

    // path to allPenaltiesPct.txt vector file
    string path ;

    // flags to replicate models from paper
    bool model1 ;
    bool model2 ;
    bool model4 ;
    int model4_count ;
    float model4_deverr ;
    
    /// read relevant information
    void read_cmd_line ( int argc, char *argv[] ) ;
    
private:

};

#endif /* COMMANDLINE_H */

