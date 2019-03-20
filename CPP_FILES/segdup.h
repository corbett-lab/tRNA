#ifndef __SEGDUP_H
#define __SEGDUP_H

void segdup( gene* old_trna, gene* new_trna, cmd_line &options ) {

	// for tRNAs that jump to a new portion of the genome, give new tRNA attributes
	// called in mutate.h in non-tandem tRNA duplications

	// segmental duplications are very close in location to their progenitors
	if ( ( gsl_ran_bernoulli( rng, 0.5 )) or ( old_trna->locus + 1 > options.map_length ) ) {
        new_trna->locus = old_trna->locus - gsl_rng_uniform( rng ) ;
    }
    else {
        new_trna->locus = old_trna->locus + gsl_rng_uniform( rng ) ;
    }

	new_trna->expression = old_trna->expression ;
	new_trna->somatic = old_trna->somatic ;
	new_trna->germline = old_trna->germline ;
	new_trna->sequence = old_trna->sequence ;
	new_trna->muts = old_trna->muts ;

	// duplication rate is set to 0 for nowak models so no need to account for them here
}

#endif
