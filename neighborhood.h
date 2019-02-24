#ifndef __NEIGHBORHOOD_H
#define __NEIGHBORHOOD_H

// TODO where does "neighborhood" term come from?

void neighborhood( Gene* old_trna, Gene* new_trna, cmd_line &options ) {

	// for tRNAs that jump to a new portion of the genome, give new tRNA attributes
	// called in mutate.h in non-tandem tRNA duplications

	// make it so that it isn't at either end of the chromosome!


    new_trna->setLocus(( options.map_length * 0.2 ) + ( gsl_rng_uniform( rng ) * ( options.map_length * 0.6 ) )) ;
	new_trna->setNeighborhood(gsl_ran_exponential( rng, options.mean_neighborhood )) ;
	new_trna->setSomatic(options.somatic_rate * new_trna->getNeighborhood());
	new_trna->setGermline(options.germline_rate * new_trna->getNeighborhood() );
	new_trna->setFunction(old_trna->getFunction() * new_trna->getNeighborhood() );

}

#endif


