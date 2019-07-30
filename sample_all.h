#ifndef __SAMPLE_ALL_H
#define __SAMPLE_ALL_H

void sample_all( int g, const vector<individual> &population, std::string sampling_out, cmd_line &options ) {
    fstream stream( sampling_out, std::fstream::in | std::fstream::out | std::fstream::app ) ;
	for ( int i = 0 ; i < population.size() ; ++i ){
		stream << g ;
		for ( int m = 0 ; m < population[i].maternal_trnas.size() ; ++m ) {
			stream.precision(8) ;
			stream << "\t" << (*population[i].maternal_trnas[m]).name << "_" << (*population[i].maternal_trnas[m]).locus << "_" ;
			stream << (*population[i].maternal_trnas[m]).sequence << "_" << (*population[i].maternal_trnas[m]).expression ;
            stream << "_" << (*population[i].maternal_trnas[m]).muts << "_" << (*population[i].maternal_trnas[m]).genotype << "_" ;
            stream << (*population[i].maternal_trnas[m]).birth << "_" ;
            stream << (*population[i].maternal_trnas[m]).progenitor << "_" << (*population[i].maternal_trnas[m]).birth_mode  ;
		}
		for ( int p = 0 ; p < population[i].paternal_trnas.size() ; ++p ) {
			stream.precision(8) ;
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