#ifndef __STATS_H
#define __STATS_H

void print_stats ( double fitness[], const vector<individual> &population, int g, list<gene*> &trna_bank, list<float> &trna_lifespans, cmd_line &options ) {
    
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
            trna_lifespans.push_back( g - (*t)->birth ) ;
            (*t)->frequency.push_back( 0 ) ;

            if ( (options.output_lifespans == true) and (g > options.burn_in) ){
                std::string lifespan_log = "lifespanLog";
                lifespan_log += std::to_string(options.run_num);
                lifespan_log += ".txt" ;
                ofstream myfile;
                myfile.open( lifespan_log, fstream::app ) ;
                myfile << (*t)->name << "\t" << (*t)->birth << "\t" << g << "\t" << g - (*t)->birth << "\t" << (*t)->progenitor << "\t" ;
                for ( int i = 1 ; i < (*t)->frequency.size() ; i ++ ) {
                        myfile << (*t)->frequency[i] << "," ;
                    }
                    myfile << endl ;
                myfile.close();
            }
			delete *t ;
			t = trna_bank.erase(t) ; 			
		}
		else {
            (*t)->frequency.push_back(found[*t]) ;
			t ++ ; 
		}
	}

    if ( g == 1 or ( g > options.burn_in and g % options.print_count == 0 ) ) {
        cout << g << "\t" ;
        cout << mean_fitness/population.size() << "\t" ; 
        cout << trna_count/population.size() << "\t" ;
        cout << trna_bank.size() ;

        for ( auto t : found ) { 
            if ( t.second > 0 ) {
                cout << "\t" << (*t.first).name << "_" << (*t.first).locus << "_" << (*t.first).function << "_" ;
                cout << (*t.first).neighborhood << "_" << (*t.first).birth << "_" << (*t.first).progenitor << "_" << t.second  ;    
                if (t.second > options.n*2) {
                    cout << "\t" << "ERROR\nERROR\nERROR" ;
                    // a tRNA can NOT be found more times than there are chromosomes. Exit if this happens.
                    exit(0);
                }
            }   
        }
        cout << endl ; 

    }
    
    

    // // DEBUGGING:


}

#endif
