/*
 * population.h
 *
 *  Created on: Feb 20, 2019
 *      Author: jcasaletto
 */

#ifndef POPULATION_CPP_
#define POPULATION_CPP_

#include "Individual.h"
#include "Population.h"
#include <vector>
#include <string>
using namespace std ;

Population::Population() {
	size = 0;
	label = "";
}
Population::Population(int n) {
	size = n;
	label = "";
	individuals = vector<Individual>(n);
}
Population::Population(int n, string l)  {
	size = n;
	label = l;
	individuals.resize(n);
	for(int i=0; i< n; i++) {
		individuals[i] = Individual();
	}

}
Population::~Population() {
}

vector<Individual> Population::getIndividuals() {
	return this->individuals;
}

void Population::setIndividuals(vector<Individual> individuals) {
	this->individuals = individuals;
	}

int Population::getSize() {
	return size;
}

void Population::setSize(int n) {
	this->individuals.resize(n);
}


#endif /* POPULATION_CPP_ */
