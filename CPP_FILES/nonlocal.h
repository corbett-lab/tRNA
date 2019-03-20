#ifndef __NONLOCAL_H
#define __NONLOCAL_H

void nonlocal( gene* old_trna, gene* new_trna, cmd_line &options ) {

	// for tRNAs that jump to a new portion of the genome, give new tRNA attributes
	// called in mutate.h in non-tandem tRNA duplications

	// new location is just random place in the genome
    new_trna->locus = ( options.map_length * 0.2 ) + ( gsl_rng_uniform( rng ) * ( options.map_length * 0.6 ) ) ;

    /// note: we are now using beta parameters for expression distribution of tRNAs in STEM CELLS ONLY:

	new_trna->expression = gsl_ran_beta( rng, 0.16780505, 0.17865490 ) ;

	// we are also using the 5' mutation rate mapping equation from polyFit.R and plotSlidingWindowSCOnly.py (STEM CELLS ONLY):

	new_trna->somatic = options.somatic_rate * ((11.8898 * (pow(new_trna->expression, 0.7415))) + 1.3932) ;
	new_trna->germline = options.germline_rate * ((11.8898 * (pow(new_trna->expression, 0.7415))) + 1.3932) ;
	new_trna->sequence = old_trna->sequence ;
	new_trna->muts = old_trna->muts ;

	// duplication rate is set to 0 for nowak models so no need to account for them here
}

#endif


