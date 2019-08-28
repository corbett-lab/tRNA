#ifndef __TRANSITION_H
#define __TRANSITION_H

void update_found( vector<individual> &population, string node, std::map<string,map<double,int>> &node_to_final_found_active, std::map<string,map<double,int>> &node_to_final_found_inactive, cmd_line &options ) {

	// node_to_final_found has key node that ties it to a map of loci, whose values are the number of 
	// chromosomes in the final population for that species that have a tRNA at that locus


	/// currently: if a locus has a tRNA at it in majority of genomes, that species has that tRNA
	/// need to separate into active and inactive separately
	/// keep current part true, but if most genes at that locus are active in that species,
	/// add it to the active count, else add it to the inactive count
	//
	/// *: separate locus_to_nodes into two dictionaries
	// node to final found is a double map: node to locus to count
	// let's change it to: node to locus to [active,inactive]

	map<double,int> locus_to_total ;
	map<double,int> locus_to_active ;
	map<double,int> locus_to_inactive ;

	/// go through and find all extant tRNAs
	for ( int p = 0 ; p < population.size() ; p ++ ) { 
        /// get extant trnas
		for ( int t = 0 ; t < population[p].maternal_trnas.size() ; t ++ ) { 
			if ( population[p].maternal_trnas[t]->expression >= 0.05 and population[p].maternal_trnas[t]->sequence >= 0.05 ){
				(locus_to_active[population[p].maternal_trnas[t]->locus]) ++ ;
			}
			else{
				(locus_to_inactive[population[p].maternal_trnas[t]->locus]) ++ ;
			}
			(locus_to_total[population[p].maternal_trnas[t]->locus]) ++ ;
		}
		for ( int t = 0 ; t < population[p].paternal_trnas.size() ; t ++ ) { 
            if ( population[p].paternal_trnas[t]->expression >= 0.05 and population[p].paternal_trnas[t]->sequence >= 0.05 ){
				(locus_to_active[population[p].paternal_trnas[t]->locus]) ++ ;
			}
			else{
				(locus_to_inactive[population[p].paternal_trnas[t]->locus]) ++ ;
			}
			(locus_to_total[population[p].paternal_trnas[t]->locus]) ++ ;
		}
	}

	for ( auto a : locus_to_total ){
		if ( a.second >= population.size() ){
			if ( !locus_to_inactive.count(a.first) or ( locus_to_active[a.first] >= locus_to_inactive[a.first] ) ){
				(node_to_final_found_active[node])[a.first] = a.second ;
			}
			else {
				(node_to_final_found_inactive[node])[a.first] = a.second ;
			}
			cout << a.first << "\t" << a.second << endl ;
		}
	}
}

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