#ifndef __INDIVIDUAL_H
#define __INDIVIDUAL_H

#include "Gene.h"
#include <vector>
using namespace std ;


// TODO refactor, but how?  here he is only calling out the "collection" of maternal and paternal trna genes, but
// in general the user may wish to model more than just one gene type --> how about chromosome?
class Individual {
private:
    vector<Gene*> maternal_trnas ;
    vector<Gene*> paternal_trnas ;

public:
    Individual();
    virtual ~Individual();
    vector<Gene*> getMaternal_trnas();
    vector<Gene*> getPaternal_trnas();


} ;

#endif
