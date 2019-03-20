#ifndef __NEIGHBORHOOD_H
#define __NEIGHBORHOOD_H

void neighborhood( gene* old_trna, gene* new_trna, cmd_line &options ) {

	// for tRNAs that jump to a new portion of the genome, give new tRNA attributes
	// called in mutate.h in non-tandem tRNA duplications

	// make it so that it isn't at either end of the chromosome!


    new_trna->locus = ( options.map_length * 0.2 ) + ( gsl_rng_uniform( rng ) * ( options.map_length * 0.6 ) ) ;
	new_trna->neighborhood = gsl_ran_beta( rng, 0.220, 0.287 ) ;
	//cout << new_trna->neighborhood << endl ;
	new_trna->somatic = options.somatic_rate * ((13.0 * (pow(new_trna->neighborhood, 0.42))) + 1.0) ;
	new_trna->germline = options.germline_rate * ((13.0 * (pow(new_trna->neighborhood, 0.42))) + 1.0) ;
	// cout << new_trna->neighborhood << "\t" << new_trna->somatic << "\t" << new_trna->germline << endl ;
	new_trna->function = old_trna->function ;

}

#endif


