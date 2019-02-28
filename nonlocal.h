#ifndef __NONLOCAL_H
#define __NONLOCAL_H

void nonlocal( gene* old_trna, gene* new_trna, cmd_line &options ) {

	// for tRNAs that jump to a new portion of the genome, give new tRNA attributes
	// called in mutate.h in non-tandem tRNA duplications

	// new location is just random place in the genome
    new_trna->locus = ( options.map_length * 0.2 ) + ( gsl_rng_uniform( rng ) * ( options.map_length * 0.6 ) ) ;
	new_trna->expression = gsl_ran_beta( rng, 0.18411666, 0.34907231 ) ;
	new_trna->somatic = options.somatic_rate * ((13.0 * (pow(new_trna->expression, 0.42))) + 1.0) ;
	new_trna->germline = options.germline_rate * ((13.0 * (pow(new_trna->expression, 0.42))) + 1.0) ;
	new_trna->sequence = old_trna->sequence ;

	// duplication rate is set to 0 for nowak models so no need to account for them here
}

#endif


