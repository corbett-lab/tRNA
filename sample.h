#ifndef __SAMPLE_H
#define __SAMPLE_H

void sample_individuals( int g, const vector<individual> &population, cmd_line &options ) {
	std::string sampling_out = std::to_string(options.run_num) + "_sample.txt" ;
    fstream stream( sampling_out, std::fstream::in | std::fstream::out | std::fstream::app ) ;
	for ( int i = 0 ; i < options.sampling_count ; ++i ){
		int random_index = rand() % population.size() ;
		stream << g ;
		for ( int m = 0 ; m < population[random_index].maternal_trnas.size() ; ++m ) {
			stream.precision(15) ;
			stream << "\t" << (*population[random_index].maternal_trnas[m]).name << "_" << (*population[random_index].maternal_trnas[m]).locus << "_" ;
			stream << (*population[random_index].maternal_trnas[m]).sequence << "_" << (*population[random_index].maternal_trnas[m]).expression ;
            stream << "_" << (*population[random_index].maternal_trnas[m]).muts << "_" << (*population[random_index].maternal_trnas[m]).genotype << "_" ;
            stream << (*population[random_index].maternal_trnas[m]).birth << "_" ;
            stream << (*population[random_index].maternal_trnas[m]).progenitor << "_" << (*population[random_index].maternal_trnas[m]).birth_mode  ;
		}
		for ( int p = 0 ; p < population[random_index].paternal_trnas.size() ; ++p ) {
			stream.precision(15) ;
			stream << "\t" << (*population[random_index].paternal_trnas[p]).name << "_" << (*population[random_index].paternal_trnas[p]).locus << "_" ;
			stream << (*population[random_index].paternal_trnas[p]).sequence << "_" << (*population[random_index].paternal_trnas[p]).expression << "_" ;
            stream << (*population[random_index].paternal_trnas[p]).muts << "_" << (*population[random_index].paternal_trnas[p]).genotype << "_" ;
            stream << (*population[random_index].paternal_trnas[p]).birth << "_" ;
            stream << (*population[random_index].paternal_trnas[p]).progenitor << "_" << (*population[random_index].paternal_trnas[p]).birth_mode  ;
		}
		stream << "\n" ;
	}
	stream.close() ;
}

void sample_all( int g, const vector<individual> &population, std::string sampling_out, cmd_line &options ) {
    fstream stream( sampling_out, std::fstream::in | std::fstream::out | std::fstream::app ) ;
	for ( int i = 0 ; i < population.size() ; ++i ){
		stream << g ;
		for ( int m = 0 ; m < population[i].maternal_trnas.size() ; ++m ) {
			stream.precision(15) ;
			stream << "\t" << (*population[i].maternal_trnas[m]).name << "_" << (*population[i].maternal_trnas[m]).locus << "_" ;
			stream << (*population[i].maternal_trnas[m]).sequence << "_" << (*population[i].maternal_trnas[m]).expression ;
            stream << "_" << (*population[i].maternal_trnas[m]).muts << "_" << (*population[i].maternal_trnas[m]).genotype << "_" ;
            stream << (*population[i].maternal_trnas[m]).birth << "_" ;
            stream << (*population[i].maternal_trnas[m]).progenitor << "_" << (*population[i].maternal_trnas[m]).birth_mode  ;
		}
		for ( int p = 0 ; p < population[i].paternal_trnas.size() ; ++p ) {
			stream.precision(15) ;
			stream << "\t" << (*population[i].paternal_trnas[p]).name << "_" << (*population[i].paternal_trnas[p]).locus << "_" ;
			stream << (*population[i].paternal_trnas[p]).sequence << "_" << (*population[i].paternal_trnas[p]).expression ;
            stream << "_" << (*population[i].paternal_trnas[p]).muts << "_" << (*population[i].paternal_trnas[p]).genotype << "_" ;
            stream << (*population[i].paternal_trnas[p]).birth << "_" ;
            stream << (*population[i].paternal_trnas[p]).progenitor << "_" << (*population[i].paternal_trnas[p]).birth_mode  ;
		}
		stream << "\n" ;
	}
	stream.close() ;
}

#endif