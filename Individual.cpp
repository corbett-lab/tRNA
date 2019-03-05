/*
 * Individual.cpp
 *
 *  Created on: Feb 21, 2019
 *      Author: jcasaletto
 */

#include "Individual.h"

Individual::Individual() {
    //genes = vector<Gene>();
}

Individual::Individual(string n, int s) {
    this->name = n;
    this->size = s;
}

Individual::Individual(const Individual& orig) {
    
}


Individual::~Individual() {
	// TODO Auto-generated destructor stub


}

int Individual::getSize() {
	return this->size;
}

void Individual::setSize(int s) {
	this->size = s;
}

string Individual::getName() {
	return this->name;
}

void Individual::setName(string n) {
	this->name = n;
}

vector<Gene*>& Individual::getGenes() {
	return this->genes;
}

void Individual::pushback(Gene * g) {
    this->genes.push_back(g);
}

void Individual::maternal_pushback(Gene * g) {
    this->maternal_trnas.push_back(g);
}

void Individual::paternal_pushback(Gene * g) {
    this->paternal_trnas.push_back(g);
} 

vector<Gene*>& Individual::getMaternal_trnas() {
    return this->maternal_trnas;
}

vector<Gene*>& Individual::getPaternal_trnas() {
    return this->paternal_trnas;
}



