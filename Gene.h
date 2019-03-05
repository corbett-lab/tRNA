/*
 * Gene.h
 *
 *  Created on: Feb 21, 2019
 *      Author: jcasaletto
 */

#ifndef GENE_H_
#define GENE_H_

#include <vector>
#include <string>
using namespace std ;

class Gene {
private:

       /// position of the tRNA
    float locus ;

    /// the name of the tRNA - just a number unique for that run
    int name ;
    
    /// sequence
    float sequence ;

    // expression
    float expression ;

    // frequency over time ;
    vector<int> frequency ;

    // generation when gene was first seen
    float birth ;  

    /// u rates per locus
    float somatic ;
    float germline ;

    // progenitor 
    int progenitor ;
    
    // sort function
    
    
    

public:
    Gene();
    Gene(const Gene& orig);
    Gene(int n);
    virtual ~Gene();
    
    void setLocus(float l);
    float getLocus();
    
    void setName(int n);
    int getName();
    
    void setExpression(float l);
    float getExpression();
    
    vector<int>& getFrequency();
    
    void setBirth(float b);
    float getBirth();
    
    void setSomatic(float s);
    float getSomatic();
    
    void setGermline(float g);
    float getGermline();
    
    void setProgenitor(int p);
    int getProgenitor();
    
    void setSequence(float s);
    float getSequence();
    
    
    bool operator <(const Gene g1 ) const;
};

#endif /* GENE_H_ */
