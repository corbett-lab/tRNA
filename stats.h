#ifndef __STATS_H
#define __STATS_H

void print_stats ( double fitness[], const vector<individual> &population, int g, list<gene*> &trna_bank, std::map<int, int> &lifespan_to_count, cmd_line &options ) {
    
    float trna_count = 0 ;
    float trna_pseudogenes = 0 ;
    float mean_fitness = 0.0 ;  
    map<gene*,int> found ; 

	/// go through and find all extant tRNAs
	for ( int p = 0 ; p < population.size() ; p ++ ) { 

		/// basic counts + fitness stats
		mean_fitness += fitness[p] ;
        trna_count += population[p].maternal_trnas.size() + population[p].paternal_trnas.size() ; 

        /// trna specific stats
		for ( int t = 0 ; t < population[p].maternal_trnas.size() ; t ++ ) { 
			found[population[p].maternal_trnas[t]] ++ ; 
		}
		for ( int t = 0 ; t < population[p].paternal_trnas.size() ; t ++ ) { 
			found[population[p].paternal_trnas[t]] ++ ; 
		}
	}

	/// go through and delete all non-existing trnas
    /// once a tRNA dies, save its lifespan

    /// here iterating through tRNA bank. allows access to frequency (attr of tRNA in bank), but not t.second
	for ( auto t = trna_bank.begin(); t != trna_bank.end() ; ) { 
		if ( !found[*t] ) { 
            /// if not this lifespan already in the dictionary, add it in
            if ( !lifespan_to_count.count( g - (*t)->birth ) ){
                lifespan_to_count.insert( std::pair<int,int> ( g-(*t)->birth, 0 ) ) ;
            }
            lifespan_to_count[ g - (*t)->birth ] += 1 ;
            // (*t)->frequency.push_back( 0 ) ;
			delete *t ;
            // reusable_pointers.push_back( *t ) ;
			t = trna_bank.erase(t) ; 			
		}
		else {
            // (*t)->frequency.push_back(found[*t]) ;
			t ++ ; 
		}
	}

    if ( g == 1 or ( g > options.burn_in and g % options.print_count == 0 ) or ( g == options.generations ) ) {

        cout << g << "\t" ;
        cout << mean_fitness/population.size() << "\t" ; 
        cout << trna_count/population.size() << "\t" ;
        cout << trna_bank.size() ;

        if ( options.quiet == false ){
            for ( auto t : found ) { 
                if ( t.second > 0 ) {
                    cout.precision(10) ;
                    cout << "\t" << (*t.first).name << "_" << (*t.first).locus << "_" << (*t.first).sequence << "_" << (*t.first).expression << "_" << (*t.first).muts ;
                    cout << "_" << (*t.first).birth << "_" << (*t.first).progenitor << "_" << (*t.first).birth_mode << "_" << t.second  ;    
                    if (t.second > options.n*2) {
                        cout << "\t" << "A tRNA HAS BEEN FOUND ON MORE CHROMOSOMES THAN THERE ARE CHROMOSOMES. EXITING PROGRAM." << endl ;
                        // a tRNA can NOT be found more times than there are chromosomes. Exit if this happens.
                        exit(0);
                    }
                }   
            }
        }
        cout << endl ; 
    }

    /// OUTPUT LIFESPAN DISTRIBUTION
    if ( g == options.generations ){

        // first get extant tRNAs at the end
        for ( auto t : found ) { 
            if ( !lifespan_to_count.count( g - (*t.first).birth ) ){
                lifespan_to_count.insert( std::pair<int,int> ( g - (*t.first).birth, 0 ) ) ;
            }
            lifespan_to_count[ g - (*t.first).birth ] += 1 ;
        }

        // then output your lifespan map to a file
        if ( options.output_lifespans == true ){
            std::string lifespan_log = "lifespanLog" + std::to_string(options.run_num) + ".txt" ;
            ofstream stream( lifespan_log ) ;
            for ( auto t : lifespan_to_count ) {
                stream << t.first << "\t" << t.second << "\n" ;
            }
            stream.close() ;
        }
    }

}

#endif
