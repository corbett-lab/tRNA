#ifndef __ASSIGN_GENOTYPE_H
#define __ASSIGN_GENOTYPE_H

void assign_genotype( gene* old_trna, gene* new_trna, std::map<string,double> &genotype_to_fitness, std::map<string,vector<string>> &genotype_to_genotypes, std::map<string,vector<double>> &genotype_to_fitnesses, cmd_line &options ) {

	// this is called in the germline block of mutate.h
	// we read in mutational_penalties in tRNA.cpp
	// here we apply them to the sequence scores of the mutated genes
	////////
	// mutations should get rid of some sequence but not affect expression most likely
	// from the classifier, expression levels are pretty static!
	// TODO: eventually add in mutations affecting expression, but these are so rare

	new_trna->muts = old_trna->muts + 1 ;
	if ( new_trna->muts < options.max_mutations ) {
		int random_index = rand() % (genotype_to_genotypes[old_trna->genotype]).size() ;
		new_trna->genotype = (genotype_to_genotypes[old_trna->genotype])[random_index] ;
		new_trna->sequence = (genotype_to_fitnesses[old_trna->genotype])[random_index] ;
	    if ( new_trna->sequence > 1.0 ) {
	    	new_trna->sequence = 1.0 ;
	    }
	    if ( new_trna->sequence < 0.0 ) {
	    	new_trna->sequence = 0.0 ;
	    }
	}
	else {
		new_trna->sequence = 0.0 ;
	}
}

#endif