/*
 * Gene.cpp
 *
 *  Created on: Feb 21, 2019
 *      Author: jcasaletto
 */

#include "Gene.h"

Gene::Gene() {

}

Gene::Gene(const Gene& orig) {
    
}

Gene::Gene(int n) {
	// TODO Auto-generated constructor stub
	this->name = n;

}

Gene::~Gene() {
	// TODO Auto-generated destructor stub
}

int Gene::getName() {
	return this->name;
}

void Gene::setName(int n) {
	this->name = n;
}

float Gene::getBirth() {
    return this->birth;
}

void Gene::setBirth(float b) {
	this->birth = b;
}

float Gene::getLocus() {
    return this->locus;
}

void Gene::setLocus(float l) {
	this->locus = l;
}

float Gene::getExpression() {
    return this->expression;
}

void Gene::setExpression(float e) {
	this->expression = e;
}

    
vector<int>& Gene::getFrequency() {
    return this->frequency;
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
    

bool Gene::operator <(const Gene g1 ) const {
        return g1.locus > this->locus ;
}

void Gene::setSequence(float s) {
	this->sequence = s;
}
    
float Gene::getSequence() {
    return this->sequence;
}

