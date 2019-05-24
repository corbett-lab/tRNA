#ifndef __LOCAL_H
#define __LOCAL_H

void local( gene* old_trna, gene* new_trna, std::map<float, int> &temp_loci, cmd_line &options ) {

	// for tRNAs that are close by, share some similarities to the progenitor
	// called in mutate.h for local tRNA duplications

	// local duplications are pretty close in location to their progenitors
	if ( ( gsl_ran_bernoulli( rng, 0.5 )) or ( old_trna->locus + 2 > options.map_length ) ) {
        new_trna->locus = old_trna->locus - 0.5 - gsl_rng_uniform( rng ) ;
        while ( temp_loci.count( new_trna->locus ) ) {
            new_trna->locus	= old_trna->locus - 0.5 - gsl_rng_uniform( rng ) ;
        }
    }
    else{
        new_trna->locus = old_trna->locus + 0.5 + gsl_rng_uniform( rng ) ;
        while ( temp_loci.count( new_trna->locus ) ) {
            new_trna->locus	= old_trna->locus + 0.5 + gsl_rng_uniform( rng ) ;
        }
    }

    // draw from normal distribution with sd ~.15
	new_trna->expression = old_trna->expression + gsl_ran_gaussian( rng, 0.15856 ) ;

	// if over 1.0, reset to 1.0 / if less than 0, reset to 0
	if (new_trna->expression > 1.0) {
		new_trna->expression = 1.0 ;
	}
	else if (new_trna->expression < 0.0) {
		new_trna->expression = 0.0 ;
	}
	new_trna->somatic = options.somatic_rate * ((options.somatic_coefficient * (pow(new_trna->expression, 0.7415))) + 1.3932) ;
	new_trna->germline = options.germline_rate * ((11.8898 * (pow(new_trna->expression, 0.7415))) + 1.3932) ;
	new_trna->sequence = old_trna->sequence ;
	new_trna->muts = old_trna->muts ;

	// duplication rate is set to 0 for nowak models so no need to account for them here
}

#endif