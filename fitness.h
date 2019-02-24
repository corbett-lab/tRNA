#ifndef __FITNESS_H
#define __FITNESS_H

/// compute fitness
/// obviously this needs a complete rework

// TODO what is the 0.05 fudge?  amount of degredation in function based on somatic mutation rate?  shouldn't it be in mutate.h then?


//void compute_fitness( double total_function_vector[], double fitness[], vector<individual*> &population, cmd_line &options ) {
void compute_fitness( double total_function_vector[], double fitness[], vector<Individual> &population, cmd_line &options ) {

    // alternative fitness function described in model4 paper:
    // if at least 1 functional tRNA, fitness = 1 - somatic rate ^ number of tRNAs
    if ( options.model4 == true ) {
        for ( int i = 0 ; i < population.size() ; i ++ ) {
            double mom_function = 0 ;
            double dad_function = 0 ;
            for ( int g = 0 ; g < population[i].getMaternal_trnas().size() ; g ++ ) {
                mom_function += (*population[i].getMaternal_trnas()[g]).getFunction() ;
                }
            for ( int g = 0 ; g < population[i].getPaternal_trnas().size() ; g ++ ) {
                dad_function += (*population[i].getPaternal_trnas()[g]).getFunction() ;
                }
            if ( (mom_function == 0) and (dad_function == 0) ){
                fitness[i] = 0.0 ;
            }
            else if ( mom_function >= dad_function ){
                fitness[i] = 1.0 - ( pow(options.model4_deverr, mom_function) ) ;
                //cout << options.model4_deverr << "\t" << mom_function << "\t" << pow(options.model4_deverr, mom_function) << endl ;
            }
            else {
                fitness[i] = 1.0 - ( pow(options.model4_deverr, dad_function) ) ;
                //cout << options.model4_deverr << "\t" << dad_function << "\t" << pow(options.model4_deverr, dad_function) << endl ;
            }
            //cout << fitness[i] << endl ;
        }
    }

    // standard fitness function used in all other models:
    // your fitness is equivalent to the highest function of your tRNA genes
    // no inherent penalty in having extra copies
    // (though need to read on if selection acting directly against redundancy exists!)

    else {
        double agg_function = 0 ;
        double mean_function = 0 ;
        double std_dev = 0 ;
        for ( int i = 0 ; i < population.size() ; i ++ ) {

            double max_function = 0 ;

            //// update with specific functional fitness model for mutations
            for ( int g = 0 ; g < population[i].getMaternal_trnas().size() ; g ++ ) {
                if ( gsl_ran_bernoulli( rng, (*population[i].getMaternal_trnas()[g]).getSomatic() ) ) {
                    if ((*population[i].getMaternal_trnas()[g]).getFunction() - 0.05 > max_function ){
                        max_function = (*population[i].getMaternal_trnas()[g]).getFunction() - 0.05 ;
                    } 
                }
                else {
                    if ((*population[i].getMaternal_trnas()[g]).getFunction() > max_function){
                        max_function = (*population[i].getMaternal_trnas()[g]).getFunction() ;
                    }
                }

            }
            for ( int g = 0 ; g < population[i].getPaternal_trnas().size() ; g ++ ) {
                if ( gsl_ran_bernoulli( rng, (*population[i].getPaternal_trnas()[g]).getSomatic() ) ) {
                    if ((*population[i].getPaternal_trnas()[g]).getFunction() - 0.05 > max_function ){
                        max_function = (*population[i].getPaternal_trnas()[g]).getFunction() - 0.05 ;
                    }
                }
                else {
                    if ((*population[i].getPaternal_trnas()[g]).getFunction() > max_function){
                        max_function = (*population[i].getPaternal_trnas()[g]).getFunction() ;
                    }
                }
            }

            fitness[i] = max_function ;
        }
    }
 
    // mean_function = agg_function/(population.size()) ;
    // cout << "mean_function:\t" << mean_function << endl ;

    // for ( int i = 0 ; i < population.size() ; i ++ ) {
    //     std_dev += pow( total_function_vector[i] - mean_function, 2 ) ;
    // }
    // // std_dev = sqrt(std_dev / population.size() ) ;
    //cout << std_dev << endl ;

    //function_to_fitness( total_function_vector, fitness, mean_function, population, options ) ;
}

#endif
