#ifndef __UPDATE_GENOTYPES_H
#define __UPDATE_GENOTYPES_H

void update_genotypes( vector<individual> &population, string node, std::map<string,map<string,int>> &node_to_final_genotypes, cmd_line &options ) {

	// node_to_final_genotypes has key node that ties it to a map of loci, whose values are the number of 
	// chromosomes in the final population for that species that have a tRNA with that genotype

	/// go through and find all extant tRNAs
	for ( int p = 0 ; p < population.size() ; p ++ ) { 
        /// get extant trnas
		for ( int t = 0 ; t < population[p].maternal_trnas.size() ; t ++ ) { 
            (node_to_final_genotypes[node])[population[p].maternal_trnas[t]->genotype] ++ ;
		}
		for ( int t = 0 ; t < population[p].paternal_trnas.size() ; t ++ ) { 
            (node_to_final_genotypes[node])[population[p].paternal_trnas[t]->genotype] ++ ;
		}
	}
}

#endif