#ifndef __ASSIGN_SEQUENCE_H
#define __ASSIGN_SEQUENCE_H

void assign_sequence( gene* old_trna, gene* new_trna, std::map<int, vector<double>> &mutations_to_function, cmd_line &options ) {

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