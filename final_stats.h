#ifndef __FINAL_STATS_H
#define __FINAL_STATS_H

void get_final_stats( vector<individual> &population, std::string sampling_out, cmd_line &options ) {

	int actives_ind = 0 ;
	int actives_total = 0 ;
	vector<int> actives_vector ( population.size() ) ;
	for ( int i = 0; i < population.size(); i++ ){
		for ( int m = 0; m < population[i].maternal_trnas.size(); m++ ){
			if ( population[i].maternal_trnas[m]->sequence > 0.05 and population[i].maternal_trnas[m]->expression > 0.05 ){
				actives_ind ++ ;
			}
		}
		for ( int p = 0; p < population[i].paternal_trnas.size(); p++ ){
			if ( population[i].paternal_trnas[p]->sequence > 0.05 and population[i].paternal_trnas[p]->expression > 0.05 ){
				actives_ind ++ ;
			}
		}
		actives_total += actives_ind ;
		actives_vector[i] = actives_ind ;
		actives_ind = 0 ;
	}

	fstream stream( sampling_out, std::fstream::out ) ;
	for ( int i = 0 ; i < actives_vector.size() ; ++i ){
		stream << actives_vector[i] << "\t" ;
	}
	stream << "\n" ;
	stream.close() ;
}

#endif