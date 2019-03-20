#ifndef __GENE_CONVERSION_H
#define __GENE_CONVERSION_H

void gene_conversion( gene* first_trna, gene* second_trna, cmd_line &options ) {

	// half the time, first gene takes sequence of second
	// other half, second gene takes sequence of first



	if ( gsl_ran_bernoulli( rng, 0.5 ) ) {
        first_trna->sequence = second_trna->sequence ;
    }
    else{
        second_trna->sequence = first_trna->sequence ;
    }

}

#endif
