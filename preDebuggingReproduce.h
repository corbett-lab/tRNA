#ifndef __REPRODUCE_H
#define __REPRODUCE_H


// sort container by locus
bool sortByLocus(gene* a, gene* b) { return (a->locus < b->locus); }


/// create chromosomes for a parent using haldane's map function 


void transmit_chromosome ( individual &parent, vector<gene*> &new_chromosome, vector<float> &my_indiv_loci ) {

	/// go through by position 
	int it_mom = 0 ; 
	int it_dad = 0 ; 

	// if mom == 1, we are currently on mom's chromosome
	float position = 0 ; 
	bool mom = gsl_ran_bernoulli( rng, 0.5 ) ;

	/// while the positiion is less than both 
	while ( it_mom < parent.maternal_trnas.size() && it_dad < parent.paternal_trnas.size() ) {

		/// both have it at the same site 
		if ( (*parent.maternal_trnas[it_mom]).locus == (*parent.paternal_trnas[it_dad]).locus ) {
			if ( gsl_ran_bernoulli( rng, ( 1 - exp( -2*(*parent.maternal_trnas[it_mom]).locus - position ) ) / 2 ) ) mom = !mom ; 
			if ( mom ) { 
				new_chromosome.push_back( parent.maternal_trnas[it_mom] ) ; 
				//cout << ( parent.maternal_trnas[it_mom]->locus ) << "\t" << 31 << endl ;
				if (std::find(my_indiv_loci.begin(), my_indiv_loci.end(), parent.maternal_trnas[it_mom]->locus ) != my_indiv_loci.end()){
					cout << "ERROR LINE 31" << endl ;
				}
				my_indiv_loci.push_back( parent.maternal_trnas[it_mom]->locus ) ;
			}
			else {
				new_chromosome.push_back( parent.paternal_trnas[it_dad] ) ; 
				//cout << ( parent.paternal_trnas[it_dad]->locus ) << "\t" << 38 << endl ;
				if (std::find(my_indiv_loci.begin(), my_indiv_loci.end(), parent.paternal_trnas[it_dad]->locus ) != my_indiv_loci.end()){
					cout << "ERROR LINE 38" << endl ;
				}
				my_indiv_loci.push_back( parent.paternal_trnas[it_dad]->locus ) ;
			}
			position = (*parent.maternal_trnas[it_mom]).locus ; 
			it_mom ++ ;
			it_dad ++ ; 
		} 

		/// if it's dad that's next
		else if ( (*parent.maternal_trnas[it_mom]).locus > (*parent.paternal_trnas[it_dad]).locus ) { 
			if ( gsl_ran_bernoulli( rng, ( 1 - exp( -2*(*parent.paternal_trnas[it_dad]).locus - position ) ) / 2 ) ) mom = !mom ; 
			if ( !mom ) { 
				new_chromosome.push_back( parent.paternal_trnas[it_dad] ) ; 
				cout << ( parent.paternal_trnas[it_dad]->locus ) << "\t" << 53 << endl ;
				if (std::find(my_indiv_loci.begin(), my_indiv_loci.end(), parent.paternal_trnas[it_dad]->locus ) != my_indiv_loci.end()){
					cout << "ERROR LINE 53" << endl ;
				}
				my_indiv_loci.push_back( parent.paternal_trnas[it_dad]->locus ) ;
			}
			position = (*parent.paternal_trnas[it_dad]).locus ;
			it_dad ++ ; 
		}

		/// otherwise, mom is next 
		else  { 
			if ( gsl_ran_bernoulli( rng, ( 1 - exp( -2*(*parent.maternal_trnas[it_mom]).locus - position ) ) / 2 ) ) mom = !mom ; 
			if ( mom ) { 
				new_chromosome.push_back( parent.maternal_trnas[it_mom] ) ; 
				cout << ( parent.maternal_trnas[it_mom]->locus ) << "\t" << 67 << endl ;
				if (std::find(my_indiv_loci.begin(), my_indiv_loci.end(), parent.maternal_trnas[it_mom]->locus ) != my_indiv_loci.end()){
					cout << "ERROR LINE 67" << endl ;
				}
				my_indiv_loci.push_back( parent.maternal_trnas[it_mom]->locus ) ;
			}
			position = (*parent.maternal_trnas[it_mom]).locus ;
			it_mom ++ ; 
		}
	}

	/// if it's dad that's last
	while ( it_dad < parent.paternal_trnas.size() ) { 
		if ( gsl_ran_bernoulli( rng, ( 1 - exp( -2*(*parent.paternal_trnas[it_dad]).locus - position ) ) / 2 ) ) mom = !mom ; 
		if ( !mom ) { 
			new_chromosome.push_back( parent.paternal_trnas[it_dad] ) ; 
			//cout << (parent.paternal_trnas[it_dad]->locus) << "\t" << 83 << endl ;
			if (std::find(my_indiv_loci.begin(), my_indiv_loci.end(), parent.paternal_trnas[it_dad]->locus ) != my_indiv_loci.end()){
					cout << (parent.paternal_trnas[it_dad]->locus) << "\t" << "ERROR LINE 83" << endl ;
					cout << it_mom << "\t" << parent.maternal_trnas.size() << endl ;
					for ( int k = 0 ; k < parent.maternal_trnas.size() ; k ++ ){
						cout << k << "\t" << parent.maternal_trnas[k]-> locus << "\tMOM" << endl;
					}
					cout << it_dad << "\t" << parent.paternal_trnas.size() << endl ;
					for ( int k = 0 ; k < parent.paternal_trnas.size() ; k ++ ){
						cout << k << "\t" << parent.paternal_trnas[k]-> locus << "\tDAD" << endl;
					}
				}
			my_indiv_loci.push_back( parent.paternal_trnas[it_dad]->locus ) ;
		}
		position = (*parent.paternal_trnas[it_dad]).locus ;
		it_dad ++ ; 
	}

	/// otherwise, mom is next 
	while ( it_mom < parent.maternal_trnas.size() ) { 
		if ( gsl_ran_bernoulli( rng, ( 1 - exp( -2*(*parent.maternal_trnas[it_mom]).locus - position ) ) / 2 ) ) mom = !mom ; 
		if ( mom ) { 
			new_chromosome.push_back( parent.maternal_trnas[it_mom] ) ;
			//cout << (parent.maternal_trnas[it_mom]->locus) << "\t" << 95 << endl ;
			if (std::find(my_indiv_loci.begin(), my_indiv_loci.end(), parent.maternal_trnas[it_mom]->locus ) != my_indiv_loci.end()){
					cout << (parent.maternal_trnas[it_mom]->locus) << "\t" << "ERROR LINE 95" << endl ;
					cout << it_mom << "\t" << parent.maternal_trnas.size() << endl ;
					for ( int k = 0 ; k < parent.maternal_trnas.size() ; k ++ ){
						cout << k << "\t" << parent.maternal_trnas[k]-> locus << "\tMOM" << endl;
					}
					cout << it_dad << "\t" << parent.paternal_trnas.size() << endl ;
					for ( int k = 0 ; k < parent.paternal_trnas.size() ; k ++ ){
						cout << k << "\t" << parent.paternal_trnas[k]-> locus << "\tDAD" << endl;
					}
				}
			my_indiv_loci.push_back( parent.maternal_trnas[it_mom]->locus ) ;
		}
		position = (*parent.maternal_trnas[it_mom]).locus ;
		it_mom ++ ; 
	}

}

/// compute fitness
void reproduce( const double *fitness, vector<individual> &population, vector<individual> &new_population, cmd_line &options ) {
    
    // populate parent multinomial samplings
    const gsl_ran_discrete_t *g = gsl_ran_discrete_preproc( population.size(), fitness ) ;
    
    /// iterate through all individuals and draw parents + recomb
    for ( int i = 0 ; i < population.size() ; i ++ ) {
    
        /// new individual
        individual new_ind ;
        
        /// grab mom and dad
        int mom = gsl_ran_discrete( rng, g ) ;
        int dad = gsl_ran_discrete( rng, g ) ;

        // NOW sort, and by locus, so this might fix the recombination issue:
        std::sort(population[i].maternal_trnas.begin() , population[i].maternal_trnas.end(), sortByLocus);
        std::sort(population[i].paternal_trnas.begin() , population[i].paternal_trnas.end(), sortByLocus);

        vector<float> my_indiv_loci ;
               
        // get their chromosomes 
        transmit_chromosome( population[mom], new_ind.maternal_trnas, my_indiv_loci ) ; 

        vector<float> my_indiv_loci2 ;

        transmit_chromosome( population[dad], new_ind.paternal_trnas, my_indiv_loci2 ) ; 

        // swap individual in
        swap( new_population[i], new_ind ) ;
    }
}

#endif
