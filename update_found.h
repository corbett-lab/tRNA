#ifndef __UPDATE_FOUND_H
#define __UPDATE_FOUND_H

void update_found( vector<individual> &population, string node, std::map<string,map<double,int>> &node_to_final_found, cmd_line &options ) {

	// node_to_final_found has key node that ties it to a map of loci, whose values are the number of 
	// chromosomes in the final population for that species that have a tRNA at that locus
	// can we do the same thing with genotypes? why not?

	/// go through and find all extant tRNAs
	for ( int p = 0 ; p < population.size() ; p ++ ) { 
        /// get extant trnas
		for ( int t = 0 ; t < population[p].maternal_trnas.size() ; t ++ ) { 
            (node_to_final_found[node])[population[p].maternal_trnas[t]->locus] ++ ;
		}
		for ( int t = 0 ; t < population[p].paternal_trnas.size() ; t ++ ) { 
            (node_to_final_found[node])[population[p].paternal_trnas[t]->locus] ++ ;
		}
	}
}

#endif