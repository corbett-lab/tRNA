#ifndef __INDIVIDUAL_CPP
#define __INDIVIDUAL_CPP

#include "Gene.h"
#include "Individual.h"
#include <vector>
using namespace std ;


// TODO refactor, but how?  here he is only calling out the "collection" of maternal and paternal trna genes, but
// in general the user may wish to model more than just one gene type --> how about chromosome?
Individual::Individual() {
}

Individual::~Individual() {
}

vector<Gene*> Individual::getMaternal_trnas() {
	return this->maternal_trnas;
}

vector<Gene*> Individual::getPaternal_trnas() {
	return this->paternal_trnas;
}

#endif
