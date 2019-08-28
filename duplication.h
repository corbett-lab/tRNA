#ifndef __DUPLICATION_H
#define __DUPLICATION_H

void local( gene* old_trna, gene* new_trna, std::map<double, int> &temp_loci, cmd_line &options ) {

	// for tRNAs that are close by, share some similarities to the progenitor
	// called in mutate.h for local tRNA duplications

	// local duplications are pretty close in location to their progenitors
	if ( ( gsl_ran_bernoulli( rng, 0.5 )) or ( old_trna->locus + (options.map_length / 100.0) > options.map_length ) ) {
        new_trna->locus = old_trna->locus - (options.map_length / 10000.0) - ( gsl_rng_uniform( rng ) * (options.map_length / 10000.0) ) ;
        while ( temp_loci.count( new_trna->locus ) ) {
            new_trna->locus	= old_trna->locus - (options.map_length / 10000.0) - ( gsl_rng_uniform( rng ) * (options.map_length / 10000.0) ) ;
        }
    }
    else{
        new_trna->locus = old_trna->locus + (options.map_length / 10000.0) + ( gsl_rng_uniform( rng ) * (options.map_length / 10000.0) ) ;
        while ( temp_loci.count( new_trna->locus ) ) {
        	cout << new_trna->locus ;
            new_trna->locus	= old_trna->locus + (options.map_length / 10000.0) + ( gsl_rng_uniform( rng ) * (options.map_length / 10000.0) ) ;
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
	new_trna->genotype = old_trna->genotype ;
	new_trna->muts = old_trna->muts ;

	// duplication rate is set to 0 for nowak models so no need to account for them here
}

void nonlocal( gene* old_trna, gene* new_trna, std::map<double, int> &temp_loci, cmd_line &options ) {

	// for tRNAs that jump to a new portion of the genome, give new tRNA attributes
	// called in mutate.h in non-tandem tRNA duplications

	// new location is just random place in the genome
    new_trna->locus = ( options.map_length * 0.1 ) + ( gsl_rng_uniform( rng ) * ( options.map_length * 0.85 ) ) ;
    while ( temp_loci.count( new_trna->locus ) ){
    	cout << new_trna->locus ;
        new_trna->locus = ( options.map_length * 0.1 ) + ( gsl_rng_uniform( rng ) * ( options.map_length * 0.85 ) ) ;
    }

    /// note: we are now using beta parameters for expression distribution of tRNAs in STEM CELLS ONLY:
    // 1) beta params for stem cells: 0.16780505, 0.17865490
    // 2) beta params for ALL GENES: 0.1116601893, 0.5867051118
    // 3) beta params for WHOLE GENOME (ignoring gene annotations completely): 0.1313023314, 2.1285553248
    // ^ for now i like option 2 best -- unlikely that all sites in the genome are equally likely for insertion.
    // genes are a good proxy for insertable sites -- if a gene can be inserted, it probably has at some point,
    // and if not, then the distribution of possible sites for insertion is probably not that different.

	new_trna->expression = gsl_ran_beta( rng, 0.1116601893, 0.5867051118 ) ;

	// we are also using the 5' mutation rate mapping equation from polyFit.R and plotSlidingWindowSCOnly.py (STEM CELLS ONLY):

	new_trna->somatic = options.somatic_rate * ((options.somatic_coefficient * (pow(new_trna->expression, 0.7415))) + 1.3932) ;
	new_trna->germline = options.germline_rate * ((11.8898 * (pow(new_trna->expression, 0.7415))) + 1.3932) ;
	new_trna->sequence = old_trna->sequence ;
	new_trna->genotype = old_trna->genotype ;
	new_trna->muts = old_trna->muts ;

	// duplication rate is set to 0 for nowak models so no need to account for them here
}

void segdup( gene* old_trna, gene* new_trna, std::map<double, int> &temp_loci, cmd_line &options ) {

	// for tRNAs that jump to a new portion of the genome, give new tRNA attributes
	// called in mutate.h in non-tandem tRNA duplications

	// segmental duplications are very close in location to their progenitors
	if ( ( gsl_ran_bernoulli( rng, 0.5 )) or ( old_trna->locus + (options.map_length / 100.0) > options.map_length ) ) {
        new_trna->locus = old_trna->locus - (options.map_length / 100000.0) - ( gsl_rng_uniform( rng ) * (options.map_length / 100000.0) ) ;
        while ( temp_loci.count( new_trna->locus ) ) {
            cout << new_trna->locus ;
            new_trna->locus	= old_trna->locus - (options.map_length / 100000.0) - ( gsl_rng_uniform( rng ) * (options.map_length / 100000.0) ) ;
        }
    }
    else {
        new_trna->locus = old_trna->locus + (options.map_length / 100000.0) + ( gsl_rng_uniform( rng ) * (options.map_length / 100000.0) ) ;
        while ( temp_loci.count( new_trna->locus ) ) {
            cout << new_trna->locus ;
            new_trna->locus = old_trna->locus + (options.map_length / 100000.0) + ( gsl_rng_uniform( rng ) * (options.map_length / 100000.0) ) ;
        }
    }

	new_trna->expression = old_trna->expression ;
	new_trna->somatic = old_trna->somatic ;
	new_trna->germline = old_trna->germline ;
	new_trna->sequence = old_trna->sequence ;
    new_trna->genotype = old_trna->genotype ;
	new_trna->muts = old_trna->muts ;

	// duplication rate is set to 0 for nowak models so no need to account for them here
}

#endif