#ifndef __FINAL_VECTORS_H
#define __FINAL_VECTORS_H

void final_vectors( std::map<string,vector<double>> &node_to_final_active_loci, std::map<string,vector<double>> &node_to_final_inactive_loci, std::map<string,map<string,int>> &node_to_final_genotypes, std::map<string, int> &node_to_Ne, std::string vector_out, cmd_line &options ) {

	map<string, int>::iterator p ;
	map<double, vector<string>> locus_to_nodes_active ;
	map<double, vector<string>> locus_to_nodes_inactive ;
	map<string, vector<string>> genotype_to_nodes ;
	vector<string> my_nodes ;

	// get all species that are at tips of the phylogeny
	// create dictionary for each locus where the key is what species it is in
	// if a locus is found in at least half of chromosomes, it is in that species
	for ( p = node_to_Ne.begin() ; p != node_to_Ne.end() ; p ++ ) {
		if ( ( (p->first).find("anc") == std::string::npos ) and ( (p->first).find("root") == std::string::npos ) and ( (p->first).find("int") == std::string::npos ) ) {
			my_nodes.push_back( p->first ) ;
			for ( int t = 0 ; t < node_to_final_active_loci[p->first].size() ; t ++ ) { 
				locus_to_nodes_active[node_to_final_active_loci[p->first][t]].push_back( p->first ) ;
			}
			for ( int t = 0 ; t < node_to_final_inactive_loci[p->first].size() ; t ++ ) { 
				locus_to_nodes_active[node_to_final_inactive_loci[p->first][t]].push_back( p->first ) ;
			}
			for ( auto l : node_to_final_genotypes[p->first] ){
				genotype_to_nodes[l.first].push_back( p->first ) ;
			}
		}
	}

	// now make dictionaries connecting all combinations of species to a unique index
	map<int, vector<string>> index_to_nodes ;
	map<vector<string>,int> nodes_to_index ;
	int x = 0 ;
	for ( int r = 0 ; r <= my_nodes.size() ; r ++ ){
		std::vector<bool> v(my_nodes.size());
		std::fill(v.begin(), v.begin() + r, true);
		vector<string> temp_vector ;
		do {
	        for (int i = 0; i < my_nodes.size(); ++i) {
	            if (v[i]) {
	                index_to_nodes[x].push_back(my_nodes[i]) ;
	            }
	        }
	        nodes_to_index[index_to_nodes[x]] = x ;
	        x ++ ;
	    } while (std::prev_permutation(v.begin(), v.end()));
	}

	// now go through and count the number of loci corresponding to each combination of species
	// output them as a vector that can easily be built into exchangeable matrix with other simulations
	vector<int> final_vector_active ( x ) ;
	vector<int> final_vector_inactive ( x ) ;
	map<double, vector<string>>::iterator l ;
	for ( l = locus_to_nodes_active.begin() ; l != locus_to_nodes_active.end() ; l ++ ) {
		cout << (l->first) << "\t" ;
		for ( int i = 0 ; i < (l->second).size() ; i ++ ){
			cout << l->second[i] << "\t" ; 
		}
		cout << endl ;
		if ( nodes_to_index.count( l->second ) ) {
			final_vector_active[nodes_to_index[l->second]] ++ ;
		}
	}

	for ( l = locus_to_nodes_inactive.begin() ; l != locus_to_nodes_inactive.end() ; l ++ ) {
		cout << (l->first) << "\t" ;
		for ( int i = 0 ; i < (l->second).size() ; i ++ ){
			cout << l->second[i] << "\t" ; 
		}
		cout << endl ;
		if ( nodes_to_index.count( l->second ) ) {
			final_vector_inactive[nodes_to_index[l->second]] ++ ;
		}
	}

	vector<int> final_g_vector ( x ) ;
	map<string, vector<string>>::iterator g ;
	for ( g = genotype_to_nodes.begin() ; g != genotype_to_nodes.end() ; g ++ ) {
		cout << (g->first) << "\t" ;
		for ( int i = 0 ; i < (g->second).size() ; i ++ ){
			cout << g->second[i] << "\t" ;
		}
		cout << endl ;
		if ( nodes_to_index.count( g->second ) ) {
			final_g_vector[nodes_to_index[g->second]] ++ ;
		}
	}

	fstream stream( vector_out, std::fstream::out ) ;
	stream << "ACTIVE_VECTOR:\n" ;
	for ( int i = 1 ; i < final_vector_active.size() ; i ++ ){
		stream << i-1 << "\t" << final_vector_active[i] << "\n" ;
	}
	stream << "INACTIVE_VECTOR:\n" ;
	for ( int i = 1 ; i < final_vector_inactive.size() ; i ++ ){
		stream << i-1 << "\t" << final_vector_inactive[i] << "\n" ;
	}
	stream << "GENOTYPE_VECTOR:\n" ;
	for ( int i = 1 ; i < final_g_vector.size() ; i ++ ){
		stream << i-1 << "\t" << final_g_vector[i] << "\n" ;
	}

}


#endif