#ifndef __GENE_CPP
#define __GENE_CPP


// need a way to record mutations and their impacts
// would also be good to record parent of each gene
#include "Gene.h"
#include <vector>

Gene::Gene() {

}

Gene::~Gene() {

}

float Gene::getLocus() {
	return this->locus;
}

void Gene::setLocus(float l) {
	this->locus = l;
}

int Gene::getName() {
	return this->name;
}

void Gene::setName(int n) {
	this->name = n;
}

float Gene::getFunction() {
	return this->function;
}

void Gene::setFunction(float f) {
	this->function = f;
}

float Gene::getNeighborhood() {
	return this->neighborhood;
}

void Gene::setNeighborhood(float n) {
	this->neighborhood = n;
}

vector<int> Gene::getFrequency() {
	return this->frequency;
}

/*void Gene::setFrequency(vector<int> f) {
	this->frequency = f;
}*/

float Gene::getBirth() {
	return this->birth;
}

void Gene::setBirth(float b) {
	this->birth = b;
}

float Gene::getSomatic() {
	return this->somatic;
}

void Gene::setSomatic(float s) {
	this->somatic = s;
}

float Gene::getGermline() {
	return this->germline;
}

void Gene::setGermline(float g) {
	this->germline = g;
}

int Gene::getProgenitor() {
	return this->progenitor;
}

void Gene::setProgenitor(int p) {
	this->progenitor = p;
}

// sort function
bool Gene::operator <(const Gene &g1 ) {
	return g1.locus > locus ;
}

#endif
