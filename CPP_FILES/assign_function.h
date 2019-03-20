#ifndef __ASSIGN_SEQUENCE_H
#define __ASSIGN_SEQUENCE_H

void assign_sequence( gene* old_trna, gene* new_trna, vector<double> &mutation_penalties, cmd_line &options ) {

	// this tRNA has a mutation
	// need to draw mutation's fitness effect,
	// calculate fitness accounting for mutation and tRNA's sequence prior to mutation

	// bit score distribution looks like two gaussians. 2/3 of the time, a mutation
	// will have a moderate effect and 1/3 of the time it has a larger effect.
	// draw from uniform to determine distribution, then draw from that distribution.

	// for the distribution of more 

	//// in the models, all mutations are assumed to be completely inactivating!

	if ( ( options.model1 == true ) or ( options.model2 == true ) or ( options.model4 == true ) ) {
		new_trna->sequence = 0 ;
	}


	// mutations should get rid of some sequence but not affect expression most likely
	// from the classifier, expression levels are pretty static!
	// TODO: replace the arbitrary 0.05 penalty with actual bit score penalties

	else {
		int random_index = rand() % mutation_penalties.size() ;
		new_trna->sequence = (old_trna->sequence + mutation_penalties[random_index]) ;
	}
}

#endif