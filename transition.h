#ifndef __TRANSITION_H
#define __TRANSITION_H

void update_found( vector<individual> &population, string node, std::map<string,vector<double>> &node_to_final_active_loci, std::map<string,vector<double>> &node_to_final_inactive_loci, std::map<string,map<string,int>> &node_to_final_genotypes, cmd_line &options ) {

	// node_to_final_found has key node that ties it to a map of loci, whose values are the number of 
	// chromosomes in the final population for that species that have a tRNA at that locus
	////
	/// OLD: if a locus has a tRNA at it in majority of genomes, that species has that tRNA
	/// (every tRNA that can be considered orthologous or """the same""" has the same locus)
	/// need to separate into active and inactive separately
	/// keep current part true, but if most genes at that locus are active in that species,
	/// add it to the active count, else add it to the inactive count
	//
	/// *: separate locus_to_nodes into two dictionaries
	// node to final found is a double map: node to locus to count
	// let's change it to: node to locus to [active,inactive]

	// sample a random individual from the population and 


	/// go through and find all extant tRNAs in random individual
	int p = rand() % population.size() ;
	if ( gsl_ran_bernoulli( rng, 0.5 ) ) {
		for ( int t = 0 ; t < population[p].maternal_trnas.size() ; t ++ ) { 
			if ( population[p].maternal_trnas[t]->expression >= 0.05 and population[p].maternal_trnas[t]->sequence >= 0.05 ){
			    node_to_final_active_loci[node].push_back(population[p].maternal_trnas[t]->locus) ;
			}
			else {
				node_to_final_inactive_loci[node].push_back(population[p].maternal_trnas[t]->locus) ;
			}
			(node_to_final_genotypes[node])[population[p].maternal_trnas[t]->genotype] ++ ;
		}
	}
	else {
        for ( int t = 0 ; t < population[p].paternal_trnas.size() ; t ++ ) { 
			if ( population[p].paternal_trnas[t]->expression >= 0.05 and population[p].paternal_trnas[t]->sequence >= 0.05 ){
			    node_to_final_active_loci[node].push_back(population[p].paternal_trnas[t]->locus) ;
			}
			else {
				node_to_final_inactive_loci[node].push_back(population[p].paternal_trnas[t]->locus) ;
			}
			(node_to_final_genotypes[node])[population[p].paternal_trnas[t]->genotype] ++ ;
		}
	}
}

#endif