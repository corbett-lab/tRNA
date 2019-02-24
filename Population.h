/*
 * population.h
 *
 */

#ifndef POPULATION_H_
#define POPULATION_H_

#include "Individual.h"
#include <string>
#include <vector>
using namespace std ;


class Population {
private:
	string label;
	int size;
	vector<Individual> individuals;

public:
	Population(); 

	Population(int n);

	Population(int n, string l);

	virtual ~Population();

	vector<Individual> getIndividuals();

	void setIndividuals(vector<Individual> individuals);

	int getSize();

	void setSize(int n);
};

#endif /* POPULATION_H_ */
