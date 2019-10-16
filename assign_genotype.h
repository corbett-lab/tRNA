#ifndef __ASSIGN_GENOTYPE_H
#define __ASSIGN_GENOTYPE_H

std::string random_string( int length ) {
    auto randchar = []() -> char{
        const char charset[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz1234567890" ;
        return charset[ rand() % 62 ] ;
    } ;
    std::string str(length,0);
    std::generate_n( str.begin(), length, randchar );
    return str;
}


void assign_genotype_gamma( gene* old_trna, gene* new_trna, cmd_line &options ) {

	// this is called in the germline block of mutate.h
	// we operate here under the idea that some proportion of mutations destroy a tRNA,
	// and the rest have effects drawn from a gamma DFE.

	new_trna->muts = old_trna->muts + 1 ;
	if ( new_trna->muts < options.max_mutations ){
		if ( gsl_ran_bernoulli( rng, options.prop_destroy ) ) {
			new_trna->genotype = random_string(2) ;
			new_trna->sequence = 0.0 ;
		}
		else{
			new_trna->sequence = old_trna->sequence - (1.0 - gsl_ran_gamma( rng, options.gamma_shape, options.gamma_scale )) ;
			new_trna->genotype = random_string(2) ;
			if ( new_trna->sequence < 0.0 ){
				new_trna->sequence = 0.0 ;
			}
			else if ( new_trna->sequence > 1.0 ){
				new_trna->sequence = 1.0 ;
			}
		}
	}
	else{
		new_trna->sequence = 0.0 ;
	}
}

void assign_genotype_pathways( gene* old_trna, gene* new_trna, std::map<string,double> &genotype_to_fitness, std::map<string,vector<string>> &genotype_to_genotypes, std::map<string,vector<double>> &genotype_to_fitnesses, cmd_line &options ) {

	// this is called in the germline block of mutate.h
	// we read in mutational_penalties in tRNA.cpp
	// here we apply them to the sequence scores of the mutated genes
	////////
	// mutations should get rid of some sequence but not affect expression most likely
	// from the classifier, expression levels are pretty static!
	// TODO: eventually add in mutations affecting expression, but these are so rare

	new_trna->muts = old_trna->muts + 1 ;
	if ( new_trna->muts < options.max_mutations ) {
		int random_index = rand() % (genotype_to_genotypes[old_trna->genotype]).size() ;
		new_trna->genotype = (genotype_to_genotypes[old_trna->genotype])[random_index] ;
		new_trna->sequence = (genotype_to_fitnesses[old_trna->genotype])[random_index] ;
	    if ( new_trna->sequence > 1.0 ) {
	    	new_trna->sequence = 1.0 ;
	    }
	    if ( new_trna->sequence < 0.0 ) {
	    	new_trna->sequence = 0.0 ;
	    }
	}
	else {
		new_trna->sequence = 0.0 ;
	}
}

void assign_genotype_model( gene* old_trna, gene* new_trna, std::map<int, vector<double>> &mutations_to_function, cmd_line &options ) {

	// this is called in the germline block of mutate.h
	// we read in mutational_penalties in tRNA.cpp
	// here we apply them to the sequence scores of the mutated genes
	////////
	// mutations should get rid of some sequence but not affect expression most likely
	// from the classifier, expression levels are pretty static!
	// TODO: eventually add in mutations affecting expression, but these are so rare

	// in the models, all mutations are assumed to be completely inactivating!
	if ( ( options.model1 == true ) or ( options.model2 == true ) or ( options.model4 == true ) ) {
		new_trna->sequence = 0.0 ;
		new_trna->muts = old_trna->muts + 1 ;
	}

	// otherwise, we look in the distribution for that tRNA's number of mutations and grab a new sequence function assignment
	else {
		new_trna->muts = old_trna->muts + 1 ;
		if ( new_trna->muts < options.max_mutations ) {
			int random_index = rand() % (mutations_to_function[new_trna->muts]).size() ;
			if ( ((mutations_to_function[new_trna->muts])[random_index]) <= old_trna->sequence ){
				new_trna->sequence = ((mutations_to_function[new_trna->muts])[random_index]) ;
			}
			else {
				new_trna->sequence = ( old_trna->sequence - ( gsl_rng_uniform( rng ) / 100.0) );
			}
		    if ( new_trna->sequence > 1.0 ) {
		    	new_trna->sequence = 1.0 ;
		    }
		    if ( new_trna->sequence < 0.0 ) {
		    	new_trna->sequence = 0.0 ;
		    }
		}
		else {
			new_trna->sequence = 0.0 ;
		}
	}
}


#endif