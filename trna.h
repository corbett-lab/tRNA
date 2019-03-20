#ifndef __TRNA_H
#define __TRNA_H

// need a way to record mutations and their impacts
// would also be good to record parent of each gene

class gene {
public:

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
    bool operator <(const gene &g1 ) {
        return g1.locus > locus ;
    }

} ;

#endif