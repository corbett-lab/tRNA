#ifndef __STATS_H
#define __STATS_H

//void print_stats ( double fitness[], const vector<individual*> &population, int g, list<gene*> &trna_bank, list<float> &trna_lifespans, cmd_line &options ) {
void print_stats ( double fitness[],   vector<Individual> &population, int g, list<Gene*> &trna_bank, list<float> &trna_lifespans, cmd_line &options ) {

    float trna_count = 0 ;
    float trna_pseudogenes = 0 ;
    float mean_fitness = 0.0 ;  
    map<Gene*,int> found ; 

	/// go through and find all extant tRNAs
	for ( int p = 0 ; p < population.size() ; p ++ ) { 

		/// basic counts + fitness stats
		mean_fitness += fitness[p] ;
        trna_count += population[p].getMaternal_trnas().size() + population[p].getPaternal_trnas().size() ;


        /// trna specific stats

		for ( int t = 0 ; t < population[p].getMaternal_trnas().size() ; t ++ ) {
			found[population[p].getMaternal_trnas()[t]] ++ ;
		}
		for ( int t = 0 ; t < population[p].getPaternal_trnas().size() ; t ++ ) {
			found[population[p].getPaternal_trnas()[t]] ++ ;
		}

	}

	/// go through and delete all non-existing trnas
    /// once a tRNA dies, save its lifespan

    /// here iterating through tRNA bank. allows access to frequency (attr of tRNA in bank), but not t.second
	for ( auto t = trna_bank.begin(); t != trna_bank.end() ; ) { 
		if ( !found[*t] ) { 
            trna_lifespans.push_back( g - (*t)->getBirth() ) ;
            (*t)->getFrequency().push_back(0);

            if ( (options.output_lifespans == true) and (g > options.burn_in) ){
                std::string lifespan_log = "lifespanLog";
                lifespan_log += std::to_string(options.run_num);
                lifespan_log += ".txt" ;
                ofstream myfile;
                myfile.open( lifespan_log, fstream::app ) ;
                myfile << (*t)->getName() << "\t" << (*t)->getBirth() << "\t" << g << "\t" << g - (*t)->getBirth() << "\t" << (*t)->getProgenitor() << "\t" ;
                for ( int i = 1 ; i < (*t)->getFrequency().size() ; i ++ ) {
                        myfile << (*t)->getFrequency()[i] << "," ;
                    }
                    myfile << endl ;
                myfile.close();
            }
			delete *t ;
			t = trna_bank.erase(t) ; 			
		}
		else {
			/*vector<int> f = (*t)->getFrequency();
			f.push_back(found[*t]);
            (*t)->setFrequency(f) ;*/
			(*t)->getFrequency().push_back(found[*t]);
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
                cout << "\t" << (*t.first).getName() << "_" << (*t.first).getLocus() << "_" << (*t.first).getFunction() << "_" ;
                cout << (*t.first).getNeighborhood() << "_" << (*t.first).getBirth() << "_" << (*t.first).getProgenitor() << "_" << t.second  ;
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
