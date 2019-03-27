/*
 * Population.cpp
 *
 *  Created on: Mar 1, 2019
 *      Author: jcasaletto
 */

#include "Population.h"
using namespace std;

Population::Population() {
	// TODO Auto-generated constructor stub
}

Population::Population(string n) {
    this->name = n;
}



Population::Population(string n, int s, vector<Individual> i) : name(n), individuals(move(i)){
    this->name = n;
    for(int j=0; j<s; j++) 
        this->individuals.push_back(Individual());
}

Population::Population(const Population& orig) {
    
}


Population::~Population() {
	// TODO Auto-generated destructor stub
}

void Population::setName(string n) {
    this->name = n;
}

string Population::getName() {
    return this->name;
}


void Population::pushback(Individual i) {
    this->individuals.push_back(i);
}

vector<Individual>& Population::getIndividuals() {
    return this->individuals;
}

bool Population::sortByLocus(Gene* a, Gene* b) {
	return (a->getLocus() < b->getLocus());
}


void Population::mutate(CommandLine &options, list<Gene*> &trna_bank, int current_gen, int &trna_counter, const gsl_rng rng) {
       // add germline mutations

	/// obviously we will want these mutations to have a range of possible effects on function
	/// bryan expect to modify this. 
	/// the current version simply makes the function 0 if there's a mutation 
    //
    // BRYAN NOTES:
    // - added assign_function.h, which assigns new function to tRNAs with mutations
    // - added neighborhood.h to give new attributes to tRNAs from jumping duplications
    // TODO: add in distribution of mutations
    // - changed locus difference on duplications to +1e-3 for visibility in printing
    // - changed sort function to hopefully resolve bug in recombination
    // - mutations that make things NOT a tRNA despite bit score change? advantageous mutations?
    //
    // TODO: update fitness function such
    //
    //
    // ISSUES:
    // - fitness must have something against total duplicate genes, otherwise will just keep growing
    // - new class of pseudogenes? tRNAs that are there but do not add to fitness. important for comparing to real data!
    // ^ could just set function to zero when it's below a certain point and make fitness function take that into account
    // - *** find consensus fitness function or come up with good one
    // - locus should never go negative or above maximum! must account for this
    // - each duplicate gene should add some fitness, to an extent
    // - how to model without biasing some ideal number of tRNAs to have
    // 

    for ( int i = 0 ; i < this->getIndividuals().size() ; i ++ ) {
        for ( int g = 0 ; g < this->getIndividuals()[i].maternal_trnas.size() ; g ++ ) {
        	if ( gsl_ran_bernoulli(&rng, this->getIndividuals()[i].maternal_trnas[g]->getGermline())) {
            	Gene* new_trna = new Gene;
            	update_gene(*new_trna, this->getIndividuals()[i].maternal_trnas, trna_counter, trna_bank, options, g, current_gen, rng);

            }
        }

        for ( int g = 0 ; g < this->getIndividuals()[i].paternal_trnas.size() ; g ++ ) {
            if ( gsl_ran_bernoulli( &rng, this->getIndividuals()[i].paternal_trnas[g]->getGermline() ) ) {
            	Gene* new_trna = new Gene;
            	update_gene(*new_trna, this->getIndividuals()[i].paternal_trnas, trna_counter, trna_bank, options, g, current_gen, rng);
            }
        }
    }

    // duplicate tRN
    for ( int i = 0 ; i < this->getIndividuals().size() ; i ++ ) {
        for ( int g = 0 ; g < this->getIndividuals()[i].maternal_trnas.size() ; g ++ ) {
            if ( gsl_ran_bernoulli( &rng, options.duplication_rate ) ) {
                /// copy current gene copy
                Gene* new_trna = new Gene;
            	update_gene_duplicate(*new_trna, this->getIndividuals()[i].maternal_trnas, trna_counter, trna_bank, options, g, current_gen, rng);
            }
        }

        for ( int g = 0 ; g < this->getIndividuals()[i].paternal_trnas.size() ; g ++ ) {
            if ( gsl_ran_bernoulli( &rng, options.duplication_rate ) ) {
                /// copy current gene copy
                Gene* new_trna = new Gene;
            	update_gene_duplicate(*new_trna, this->getIndividuals()[i].paternal_trnas, trna_counter, trna_bank, options, g, current_gen, rng);


            }
        }
    }

    // delete tRNAs
    for ( int i = 0 ; i < this->getIndividuals().size() ; i ++ ) {
        for ( int g = this->getIndividuals()[i].maternal_trnas.size() -1 ; g > -1 ; g -- ) {
            if ( gsl_ran_bernoulli( &rng, options.deletion_rate ) ) {
                this->getIndividuals()[i].maternal_trnas.erase( this->getIndividuals()[i].maternal_trnas.begin() + g ) ;
            }
        }
        for ( int g = this->getIndividuals()[i].paternal_trnas.size() -1 ; g > -1 ; g -- ) {
            if ( gsl_ran_bernoulli( &rng, options.deletion_rate ) ) {
                this->getIndividuals()[i].paternal_trnas.erase( this->getIndividuals()[i].paternal_trnas.begin() + g ) ;

            }
        }
    }



    for ( int i = 0 ; i < this->getIndividuals().size() ; i ++ ) {
        // NOW sort, and by locus, so this might fix the recombination issue:
        // NO, next is fitness so we are going to re-sort by function anyway
        /*sort(this->getIndividuals()[i].maternal_trnas.begin() , this->getIndividuals()[i].maternal_trnas.end());
        sort(this->getIndividuals()[i].paternal_trnas.begin() , this->getIndividuals()[i].paternal_trnas.end()); */
        std::sort(this->getIndividuals()[i].maternal_trnas.begin() , this->getIndividuals()[i].maternal_trnas.end(), Population::sortByLocus);
        std::sort(this->getIndividuals()[i].paternal_trnas.begin() , this->getIndividuals()[i].paternal_trnas.end(), Population::sortByLocus);

    }

}

void Population::update_gene(Gene &new_trna, vector<Gene*> &trnas, int &trna_counter, list<Gene*> &trna_bank, CommandLine &options, int g, int current_gen, const gsl_rng rng) {
    new_trna.setLocus(trnas[g]->getLocus()) ;
    new_trna.setSomatic(trnas[g]->getSomatic()) ;
    new_trna.setGermline(trnas[g]->getGermline()) ;
    new_trna.setExpression(trnas[g]->getExpression()) ;
    new_trna.setBirth(current_gen);
    new_trna.getFrequency().push_back( 0 ) ;
    new_trna.setProgenitor(trnas[g]->getName()) ;
    assign_function(trnas[g], new_trna, options, rng) ;
    trna_counter ++ ;
    new_trna.setName(trna_counter);
    trna_bank.push_back( &new_trna ) ;
    trnas[g] = &new_trna ;
    //trnas.push_back( &new_trna);
}



void Population::update_gene_duplicate(Gene &new_trna, vector<Gene*> &trnas, int &trna_counter, list<Gene*> &trna_bank, CommandLine &options, int g, int current_gen, const gsl_rng rng) {
	if ( gsl_ran_bernoulli( &rng, 0.5 ) or (trnas[g]->getLocus() + 1 > options.map_length) ) {
		 new_trna.setLocus(trnas[g]->getLocus() - gsl_rng_uniform( &rng )) ;
	 }
	 else{
		 new_trna.setLocus(trnas[g]->getLocus() + gsl_rng_uniform( &rng )) ;
	 }
	 new_trna.setSomatic(trnas[g]->getSomatic());
	 new_trna.setGermline(trnas[g]->getGermline());
	 new_trna.setSequence(trnas[g]->getSequence());
	 new_trna.setExpression(trnas[g]->getExpression());
	 new_trna.setBirth(current_gen);
	 new_trna.getFrequency().push_back( 0 ) ;
	 new_trna.setProgenitor(trnas[g]->getName());

	 /// if the trna jumps, need to reassigned neighborhood etc
	 if ( gsl_ran_bernoulli( &rng, 1 - options.prob_cluster ) ) {
		 neighborhood(trnas[g], new_trna, options, rng ) ;
	 }
	 trna_counter ++ ;
	 new_trna.setName(trna_counter);
	 trna_bank.push_back( &new_trna ) ;
	 trnas.push_back( &new_trna ) ;

}


void Population::assign_function( Gene* old_trna, Gene new_trna, CommandLine &options, const gsl_rng rng) {

	// this tRNA has a mutation
	// need to draw mutation's fitness effect,
	// calculate fitness accounting for mutation and tRNA's function prior to mutation

	// bit score distribution looks like two gaussians. 2/3 of the time, a mutation
	// will have a moderate effect and 1/3 of the time it has a larger effect.
	// draw from uniform to determine distribution, then draw from that distribution.

	// for the distribution of more 

	//// in the models, all mutations are assumed to be completely inactivating!

	if ( ( options.model1 == true ) or ( options.model2 == true ) or ( options.model4 == true ) ) {
		new_trna.setSequence(0 );
	}



	// mutations should get rid of some function but not affect expression most likely
	// from the classifier, expression levels are pretty static!
	// TODO: replace the arbitrary 0.05 penalty with actual bit score penalties

        // TODO does this reflect the paragraph above?

	else {
		if ( old_trna->getSequence() < 0.05 ){
			new_trna.setSequence(0) ;
		}
		else {
			if ( gsl_rng_uniform( &rng ) < 0.99 ) {
				new_trna.setSequence(( old_trna->getSequence() - ( 0.05 ) ) * new_trna.getExpression()) ;
			}
			else {
				new_trna.setSequence(( old_trna->getSequence() - ( 0.001 ) ) * new_trna.getExpression());
			} 
		}
	}
}

void Population::neighborhood( Gene* old_trna, Gene new_trna, CommandLine &options, const gsl_rng rng ) {

	// for tRNAs that jump to a new portion of the genome, give new tRNA attributes
	// called in mutate.h in non-tandem tRNA duplications

	// make it so that it isn't at either end of the chromosome!


    new_trna.setLocus(( options.map_length * 0.2 ) + ( gsl_rng_uniform( &rng ) * ( options.map_length * 0.6 ) )) ;
	new_trna.setExpression(gsl_ran_exponential( &rng, options.mean_neighborhood )) ;
	new_trna.setSomatic(options.somatic_rate * new_trna.getExpression());
	new_trna.setGermline(options.germline_rate * new_trna.getExpression() );
	new_trna.setSequence(old_trna->getSequence() * new_trna.getExpression() );

}

void Population::compute_fitness(double fitness[],  CommandLine &options, const gsl_rng rng ) {

    // alternative fitness function described in model4 paper:
    // if at least 1 functional tRNA, fitness = 1 - somatic rate ^ number of tRNAs
    if ( options.model4 == true ) {
        for ( int i = 0 ; i < this->getIndividuals().size() ; i ++ ) {
            double mom_function = 0 ;
            double dad_function = 0 ;
            for ( int g = 0 ; g < this->getIndividuals()[i].getMaternal_trnas().size() ; g ++ ) {
                //mom_function += (*population[i]->getMaternal_trnas()[g]).getFunction() ;
                mom_function += this->getIndividuals()[i].getMaternal_trnas()[g]->getSequence() ;

                }
            for ( int g = 0 ; g < this->getIndividuals()[i].getPaternal_trnas().size() ; g ++ ) {
                //dad_function += (*population[i]->getPaternal_trnas()[g]).getFunction() ;
                dad_function += this->getIndividuals()[i].getPaternal_trnas()[g]->getSequence() ;

            }
            if ( (mom_function == 0) and (dad_function == 0) ){
                fitness[i] = 0.0 ;
            }
            else if ( mom_function >= dad_function ){
                fitness[i] = 1.0 - ( pow(options.model4_deverr, mom_function) ) ;
                //cout << options.model4_deverr << "\t" << mom_function << "\t" << pow(options.model4_deverr, mom_function) << endl ;
            }
            else {
                fitness[i] = 1.0 - ( pow(options.model4_deverr, dad_function) ) ;
                //cout << options.model4_deverr << "\t" << dad_function << "\t" << pow(options.model4_deverr, dad_function) << endl ;
            }
            //cout << fitness[i] << endl ;
        }
    }

    // standard fitness function used in all other models:
    // your fitness is equivalent to the highest function of your tRNA genes
    // no inherent penalty in having extra copies
    // (though need to read on if selection acting directly against redundancy exists!)

    else {
        double agg_function = 0 ;
        double mean_function = 0 ;
        double std_dev = 0 ;
        for ( int i = 0 ; i < this->getIndividuals().size() ; i ++ ) {

            double max_function = 0 ;

            //// update with specific functional fitness model for mutations
            for ( int g = 0 ; g < this->getIndividuals()[i].getMaternal_trnas().size() ; g ++ ) {
                if ( gsl_ran_bernoulli( &rng, this->getIndividuals()[i].getMaternal_trnas()[g]->getSomatic() ) ) {
                    if (this->getIndividuals()[i].getMaternal_trnas()[g]->getSequence() - 0.05 > max_function ){
                        max_function = this->getIndividuals()[i].getMaternal_trnas()[g]->getSequence() - 0.05 ;
                    } 
                }
                else {
                    if (this->getIndividuals()[i].getMaternal_trnas()[g]->getSequence() > max_function){
                        max_function = this->getIndividuals()[i].getMaternal_trnas()[g]->getSequence() ;
                    }
                }

            }
            for ( int g = 0 ; g < this->getIndividuals()[i].getPaternal_trnas().size() ; g ++ ) {
                if ( gsl_ran_bernoulli( &rng, this->getIndividuals()[i].getPaternal_trnas()[g]->getSomatic() ) ) {
                    if (this->getIndividuals()[i].getPaternal_trnas()[g]->getSequence() - 0.05 > max_function ){
                        max_function = this->getIndividuals()[i].getPaternal_trnas()[g]->getSequence() - 0.05 ;
                    }
                }
                else {
                    if (this->getIndividuals()[i].getPaternal_trnas()[g]->getSequence() > max_function){
                        max_function = this->getIndividuals()[i].getPaternal_trnas()[g]->getSequence() ;
                    }
                }
            }

            fitness[i] = max_function ;
        }
    }
 
 
}
