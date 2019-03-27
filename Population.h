/*
 * Population.h
 *
 *  Created on: Mar 1, 2019
 *      Author: jcasaletto
 */

#ifndef POPULATION_H_
#define POPULATION_H_
#include "Individual.h"
#include "CommandLine.h"
#include <list>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <iostream>
#include <math.h>
#include <algorithm>

class Population {
private:
    int size;
    string name;
    vector<Individual> individuals;

public:
	Population();
        Population(string n);
        Population(string n, int s, vector<Individual> i);
        Population(const Population& orig);

	virtual ~Population();
        
        string getName();
        
        void setName(string n);
        
        void pushback(Individual i);
        
        vector<Individual>& getIndividuals();
        
        void mutate(CommandLine &options, list<Gene*> &trna_bank, int current_gen, int &trna_counter, const gsl_rng rng) ;
        
        void update_gene(Gene &new_trna, vector<Gene*> &trnas, int &trna_counter, list<Gene*> &trna_bank, CommandLine &options, int g, int current_gen, const gsl_rng rng);
       
        void update_gene_duplicate(Gene &new_trna, vector<Gene*> &trnas, int &trna_counter, list<Gene*> &trna_bank, CommandLine &options, int g, int current_gen, const gsl_rng rng);
        
        void assign_function(Gene* old_trna, Gene new_trna, CommandLine &options, const gsl_rng rng);
        
        void neighborhood( Gene* old_trna, Gene new_trna, CommandLine &options, const gsl_rng rng );
        
        void compute_fitness( double fitness[], CommandLine &options, const gsl_rng rng );

        static bool sortByLocus(Gene* a, Gene* b);

        
 

};

#endif /* POPULATION_H_ */
