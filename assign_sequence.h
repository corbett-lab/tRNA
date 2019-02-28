#ifndef __ASSIGN_SEQUENCE_H
#define __ASSIGN_SEQUENCE_H

void assign_sequence( gene* old_trna, gene* new_trna, vector<double> &mutation_penalties, cmd_line &options ) {

	// this is called in the germline block of mutate.h
	// we read in mutational_penalties in tRNA.cpp
	// here we apply them to the sequence scores of the mutated genes
	////////
	// mutations should get rid of some sequence but not affect expression most likely
	// from the classifier, expression levels are pretty static!
	// TODO: eventually add in mutations affecting expression, but these are so rare

	// in the models, all mutations are assumed to be completely inactivating!
	if ( ( options.model1 == true ) or ( options.model2 == true ) or ( options.model4 == true ) ) {
		new_trna->sequence = 0 ;
	}

	// otherwise, we draw a random mutation penalty from the read-in distribution
	else {
		int random_index = rand() % mutation_penalties.size() ;
		new_trna->sequence = (old_trna->sequence + mutation_penalties[random_index]) ;
		if ( new_trna->sequence > 1.0 ) {
			new_trna->sequence = 1.0 ;
		}
		else if ( new_trna->sequence < 0.0 ) {
			new_trna->sequence = 0.0 ;
		}
	}
}

#endif