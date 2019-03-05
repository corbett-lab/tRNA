/*
 * Individual.h
 *
 *  Created on: Feb 21, 2019
 *      Author: jcasaletto
 */

#ifndef INDIVIDUAL_H_
#define INDIVIDUAL_H_
#include "Gene.h"
#include <vector>

class Individual {
public:
	vector<Gene*> genes;
	int size;
	string name;
        vector<Gene*> maternal_trnas ; 
        vector<Gene*> paternal_trnas ;

	Individual();     
        Individual(string s, int n);
        Individual(const Individual& orig);
	virtual ~Individual();

	int getSize();
	void setSize(int s);
        
	string getName();
	void setName(string n);
        
	vector<Gene*>& getGenes();
        void pushback(Gene* g);
        
        vector<Gene*>& getMaternal_trnas();
        void maternal_pushback(Gene* g);
        
        vector<Gene*>& getPaternal_trnas();
        void paternal_pushback(Gene* g);
        
};

#endif /* INDIVIDUAL_H_ */
