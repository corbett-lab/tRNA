#ifndef __STATS_H
#define __STATS_H

void print_stats ( double fitness[], const vector<individual> &population, int g, list<gene*> &trna_bank, std::map<double, int> &loci_to_lifespans, std::map<int, int> &lifespan_to_count, cmd_line &options ) {
    
    float trna_count = 0 ;
    float trna_pseudogenes = 0 ;
    float mean_fitness = 0.0 ;  
    map<gene*,int> found ; 
    map<double,int> found_loci ;

	/// go through and find all extant tRNAs
	for ( int p = 0 ; p < population.size() ; p ++ ) { 

		/// basic counts + fitness stats
		mean_fitness += fitness[p] ;
        trna_count += population[p].maternal_trnas.size() + population[p].paternal_trnas.size() ; 

        /// get extant trnas and update lifespan map
		for ( int t = 0 ; t < population[p].maternal_trnas.size() ; t ++ ) { 
			found[population[p].maternal_trnas[t]] ++ ; 
            found_loci[population[p].maternal_trnas[t]->locus] ++ ;
		}
		for ( int t = 0 ; t < population[p].paternal_trnas.size() ; t ++ ) { 
			found[population[p].paternal_trnas[t]] ++ ; 
            found_loci[population[p].paternal_trnas[t]->locus] ++ ;
		}
	}

    if ( options.output_lifespans ){
        // have collected all extant loci in found_loci in map so no double-counting
        // if locus in found_loci but not in loci_to_lifespans, add to loci_to_lifespans + 1
        for ( auto f : found_loci ){
            loci_to_lifespans[f.first] ++ ;
        }

        // if locus in loci_to_lifespans but not in found_loci, that means this locus's lifespan is over
        // - add lifespan length to lifespan_to_count as a key ++
        // - erase from loci_to_lifespans
        map<double, int>::iterator f ; 
        for ( f = loci_to_lifespans.begin() ; f != loci_to_lifespans.end() ; f ++ ){
            if ( !found_loci.count(f->first) ){
                lifespan_to_count[f->second] ++ ;
                loci_to_lifespans.erase( f->first ) ;
            }
        }
    }

	/// go through and delete all non-existing trnas
	for ( auto t = trna_bank.begin(); t != trna_bank.end() ; ) { 
		if ( !found[*t] ) { 
			delete *t ;
			t = trna_bank.erase(t) ; 
		}
		else {
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
                    cout << "_" << (*t.first).genotype << "_" << (*t.first).birth << "_" << (*t.first).progenitor << "_" << (*t.first).birth_mode << "_" << t.second  ;    
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

    // OUTPUT LIFESPAN DISTRIBUTION
    if ( g == options.generations ){

        map<double, int>::iterator f ; 
        for ( f = loci_to_lifespans.begin() ; f != loci_to_lifespans.end() ; f ++ ){
            lifespan_to_count[f->second] ++ ;
        }

        // output lifespan map to a file
        if ( options.output_lifespans == true ){
            std::string lifespan_log = "lifespanLog" + std::to_string(options.run_num) + ".txt" ;
            ofstream stream( lifespan_log ) ;
            for ( auto l : lifespan_to_count ) {
                stream << l.first << "\t" << l.second << "\n" ;
            }
            stream.close() ;
        }
    }

}

#endif
