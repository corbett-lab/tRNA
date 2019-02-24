#ifndef __GENE_H
#define __GENE_H


// need a way to record mutations and their impacts
// would also be good to record parent of each gene
#include <vector>
using namespace std;

class Gene {
public:
    Gene();
    ~Gene();

    float getLocus();
    void setLocus(float l);
    int getName();
    void setName(int n);
    float getFunction();
    void setFunction(float f);
    float getNeighborhood();
    void setNeighborhood(float n);
    vector<int> getFrequency();
    //void setFrequency(vector<int> f);
    float getBirth();
    void setBirth(float b);
    float getSomatic();
    void setSomatic(float s);
    float getGermline();
    void setGermline(float g);
    int getProgenitor();
    void setProgenitor(int p);

    // sort function
    bool operator <(const Gene &g1 );

private:

    /// position of the tRNA
    float locus ;

    /// the name of the tRNA - just a number unique for that run
    int name ;
    
    /// function
    float function ;

    // neighborhood impact
    float neighborhood ;

    // frequency over time ;
    vector<int> frequency ;

    // generation when gene was first seen
    float birth ;  

    /// u rates per locus
    float somatic ;
    float germline ;

    // progenitor 
    int progenitor ;
    


} ;

#endif
