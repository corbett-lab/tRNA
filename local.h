#ifndef __LOCAL_H
#define __LOCAL_H

void local( gene* old_trna, gene* new_trna, cmd_line &options ) {

	// for tRNAs that are close by, share some similarities to the progenitor
	// called in mutate.h for local tRNA duplications

	// local duplications are pretty close in location to their progenitors
	if ( gsl_ran_bernoulli( rng, 0.5 ) or (old_trna->locus + 1 > options.map_length) ) {
        new_trna->locus = old_trna->locus - gsl_rng_uniform( rng ) ;
    }
    else{
        new_trna->locus = old_trna->locus + gsl_rng_uniform( rng ) ;
    }

    // draw from normal distribution with sd ~.15
	new_trna->expression = old_trna->expression + gsl_ran_gaussian( rng, 0.15856 ) ;
	// if over 1.0, re-set to 1.0
	if (new_trna->expression > 1.0) {
		new_trna->expression = 1.0 ;
	}
	new_trna->somatic = options.somatic_rate * ((13.0 * (pow(new_trna->expression, 0.42))) + 1.0) ;
	new_trna->germline = options.germline_rate * ((13.0 * (pow(new_trna->expression, 0.42))) + 1.0) ;
	new_trna->sequence = old_trna->sequence ;

	// duplication rate is set to 0 for nowak models so no need to account for them here
}

#endif