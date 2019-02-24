#ifndef __ASSIGN_FUNCTION_H
#define __ASSIGN_FUNCTION_H

void assign_function( Gene* old_trna, Gene* new_trna, cmd_line &options ) {

	// this tRNA has a mutation
	// need to draw mutation's fitness effect,
	// calculate fitness accounting for mutation and tRNA's function prior to mutation

	// bit score distribution looks like two gaussians. 2/3 of the time, a mutation
	// will have a moderate effect and 1/3 of the time it has a larger effect.
	// draw from uniform to determine distribution, then draw from that distribution.

	// for the distribution of more 

	//// in the models, all mutations are assumed to be completely inactivating!

	/*if ( ( options.model1 == true ) or ( options.model2 == true ) or ( options.model4 == true ) ) {
		new_trna->function = 0 ;
	}*/

	if ( ( options.model == 1 ) or ( options.model == 2 ) or ( options.model == 4 ) ) {
			new_trna->setFunction(0 );
		}

	// mutations should get rid of some function but not affect expression most likely
	// from the classifier, expression levels are pretty static!
	// TODO: replace the arbitrary 0.05 penalty with actual bit score penalties

	else {
		if ( old_trna->getFunction() < 0.05 ){
			new_trna->setFunction(0) ;
		}
		else {
			if ( gsl_rng_uniform( rng ) < 0.99 ) {
				new_trna->setFunction(( old_trna->getFunction() - ( 0.05 ) ) * new_trna->getNeighborhood()) ;
			}
			else {
				new_trna->setFunction(( old_trna->getFunction() - ( 0.001 ) ) * new_trna->getNeighborhood());
			} 
		}
	}
}

#endif
