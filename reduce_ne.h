#ifndef __REDUCE_NE_H
#define __REDUCE_NE_H

void reduce_ne( cmd_line &options, vector<individual> &old_population, vector<individual> &new_population, vector<gene*> &old_trna_bank, vector<gene*> &new_trna_bank, std::map<gene*,gene*> &old_trna_to_new_trna ) {

	for ( int i = 0 ; i < old_trna_bank.size() ; ++i ) {
		gene* new_trna = ::new gene ;
		new_trna->name = (old_trna_bank[i])->name ;
		new_trna->locus = (old_trna_bank[i])->locus ;
		new_trna->sequence = (old_trna_bank[i])->sequence ;
		new_trna->expression = (old_trna_bank[i])->expression ;
		new_trna->muts = (old_trna_bank[i])->muts ;
		new_trna->genotype = (old_trna_bank[i])->genotype ;
		new_trna->birth = (old_trna_bank[i])->birth ;
		new_trna->progenitor = (old_trna_bank[i])->progenitor ;
		new_trna->birth_mode = (old_trna_bank[i])->birth_mode ;
		old_trna_to_new_trna[ old_trna_bank[i] ] = new_trna ;
		new_trna_bank.push_back( new_trna ) ;
	}

	for ( int i = 0 ; i < new_population.size() ; ++i ){
		int random_index = rand() % old_population.size() ;
		for ( int m = 0 ; m < old_population[random_index].maternal_trnas.size() ; ++m ) {
			new_population[i].maternal_trnas.push_back( old_trna_to_new_trna[old_population[random_index].maternal_trnas[m]] ) ;
		}
		for ( int p = 0 ; p < old_population[random_index].paternal_trnas.size() ; ++p ) {
			new_population[i].paternal_trnas.push_back( old_trna_to_new_trna[old_population[random_index].paternal_trnas[p]] ) ;
		}
	}
}

#endif
    