#ifndef __GENE_CONVERSION_H
#define __GENE_CONVERSION_H


// function that takes in genes and actually does the conversion event

void non_allelic( gene* new_trna, gene* locus_trna, gene* sequence_trna ) {

	// give new gene chromatin environment characteristics of locus_trna
	new_trna->locus = locus_trna->locus ;
    new_trna->birth = locus_trna->birth ;
	new_trna->expression = locus_trna->expression ;
	new_trna->somatic = locus_trna->somatic ;
	new_trna->germline = locus_trna->germline ;

	// give new gene sequence characteristics of sequence_trna
	new_trna->muts = sequence_trna->muts ;
	new_trna->sequence = sequence_trna->sequence ;
    new_trna->genotype = sequence_trna->genotype ;
	new_trna->progenitor = sequence_trna->name ;
    new_trna->birth_mode = 'c' ;

}

// parse through all individuals and do gene conversion

void gene_conversion( vector<individual> &population, cmd_line &options, list<gene*> &trna_bank, int current_gen, int &trna_counter ) {

	// go through all individuals
	for ( int i = 0 ; i < population.size() ; ++i ) {

		// sort that individual's tRNAs on each chromosome by locus
        std::sort(population[i].maternal_trnas.begin() , population[i].maternal_trnas.end(), sortByLocus);
        std::sort(population[i].paternal_trnas.begin() , population[i].paternal_trnas.end(), sortByLocus);

        // for now, only allow one gene conversion event per individual per generation (should be realistic).
        // i have read a bunch of papers and there seems to be very conflicting evidence on relationship between
        // paralogous genomic distance and gene conversion frequency. so, let's pick two random tRNAs.

        // maternal block
        if ( population[i].maternal_trnas.size() > 1 ) {
        	if ( gsl_ran_bernoulli( rng, options.gene_conversion_rate ) ) {

        		// if this is the case, we have to make a new gene, which will
        		// take the locus of one extant gene and the sequence characteristics of another.
        		// we will also erase the locus gene from that chromosome.

        		gene* new_trna = ::new gene ; 
                new_trna->birth = current_gen ;
                trna_counter ++ ;
                new_trna->name = trna_counter ;

        		if ( population[i].maternal_trnas.size() == 2 ){

        			if ( gsl_ran_bernoulli( rng, 0.5 ) ) {
        			    non_allelic( new_trna, population[i].maternal_trnas[0], population[i].maternal_trnas[1] ) ;
        			    population[i].maternal_trnas.erase( population[i].maternal_trnas.begin() ) ;
        			}
        			else {
        			    non_allelic( new_trna, population[i].maternal_trnas[1], population[i].maternal_trnas[0] ) ;
        			    population[i].maternal_trnas.erase( population[i].maternal_trnas.begin() + 1 ) ;
        			}       			    
        		}
        		else {
        			int random_index = rand() % population[i].maternal_trnas.size() ;
        			if ( random_index != 0 ){
        				non_allelic( new_trna, population[i].maternal_trnas[random_index - 1], population[i].maternal_trnas[random_index] ) ;
        				population[i].maternal_trnas.erase( population[i].maternal_trnas.begin() + (random_index - 1) ) ;
        			}
        			else {
        				non_allelic( new_trna, population[i].maternal_trnas[0], population[i].maternal_trnas[1] ) ;
        				population[i].maternal_trnas.erase( population[i].maternal_trnas.begin() ) ;
        			}
        		}

        		// have to do this at the end or else the new gene gets included in the loop!
        		// new_trna->frequency.push_back( 0 ) ;
            	trna_bank.push_back( new_trna ) ; 
                (population[i].maternal_trnas).push_back( new_trna ) ; 
        	}
        }

        // paternal block
        if ( population[i].paternal_trnas.size() > 1 ) {
        	if ( gsl_ran_bernoulli( rng, options.gene_conversion_rate ) ) {

        		// if this is the case, we have to make a new gene, which will
        		// take the locus of one extant gene and the sequence characteristics of another.
        		// we will also erase the locus gene from that chromosome.

        		gene* new_trna = ::new gene ; 
                trna_counter ++ ;
                new_trna->name = trna_counter ;

        		if ( population[i].paternal_trnas.size() == 2 ){

        			if ( gsl_ran_bernoulli( rng, 0.5 ) ) {
        			    non_allelic( new_trna, population[i].paternal_trnas[0], population[i].paternal_trnas[1] ) ;
        			    population[i].paternal_trnas.erase( population[i].paternal_trnas.begin() ) ;
        			}
        			else {
        			    non_allelic( new_trna, population[i].paternal_trnas[1], population[i].paternal_trnas[0] ) ;
        			    population[i].paternal_trnas.erase( population[i].paternal_trnas.begin() + 1 ) ;
        			}       			    
        		}
        		else {
        			int random_index = rand() % population[i].paternal_trnas.size() ;
        			if ( random_index != 0 ){
        				non_allelic( new_trna, population[i].paternal_trnas[random_index - 1], population[i].paternal_trnas[random_index] ) ;
        				population[i].paternal_trnas.erase( population[i].paternal_trnas.begin() + (random_index - 1) ) ;
        			}
        			else {
        				non_allelic( new_trna, population[i].paternal_trnas[0], population[i].paternal_trnas[1] ) ;
        				population[i].paternal_trnas.erase( population[i].paternal_trnas.begin() ) ;
        			}
        		}

        		// have to do this at the end or else the new gene gets included in the loop!
        		// new_trna->frequency.push_back( 0 ) ;
            	trna_bank.push_back( new_trna ) ; 
                (population[i].paternal_trnas).push_back( new_trna ) ; 
        	}
        }

        // sort that individual's tRNAs on each chromosome by locus
        std::sort(population[i].maternal_trnas.begin() , population[i].maternal_trnas.end(), sortByLocus);
        std::sort(population[i].paternal_trnas.begin() , population[i].paternal_trnas.end(), sortByLocus);
    }
}

#endif
