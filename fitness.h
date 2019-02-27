#ifndef __FITNESS_H
#define __FITNESS_H

/// compute fitness
void compute_fitness( double fitness[], vector<individual> &population, vector<double> &mutation_penalties, cmd_line &options ) {

    /////////////////////
    ////// MODEL 4 //////
    /////////////////////
    // alternative fitness function described in model4 paper:
    // if at least 1 functional tRNA, fitness = 1 - somatic rate ^ number of tRNAs
    // entirely dependent on sequence -- expression assumed the same
    // linear mapping of function to fitness

    if ( options.model4 == true ) {
        for ( int i = 0 ; i < population.size() ; i ++ ) {

            // get total function of all of an individual's tRNAs
            double mom_function = 0 ;
            double dad_function = 0 ;

            // maternal fitness block
            for ( int g = 0 ; g < population[i].maternal_trnas.size() ; g ++ ) {
                mom_function += (*population[i].maternal_trnas[g]).sequence ;
                }

            // paternal fitness block
            for ( int g = 0 ; g < population[i].paternal_trnas.size() ; g ++ ) {
                dad_function += (*population[i].paternal_trnas[g]).sequence ;
                }

            // if they have no tRNAs, set their fitness to 0
            if ( (mom_function == 0) and (dad_function == 0) ){
                fitness[i] = 0.0 ;
            }

            // take whichever chromosome has more total function and take error rate raised to that power
            // (total function == number of tRNAs here because they're all 0.0 and 1.0 in this model)
            // (also originally a haploid model so this our adaptation to diploid)
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

    //////////////////////////
    ////// MODELS 1 & 2 //////
    //////////////////////////
    // in models 1 and 2, just take the best functioning tRNA and that's your fitness
    else if (( options.model1 == true ) or ( options.model2 == true )) {

        for ( int i = 0 ; i < population.size() ; i ++ ) {
            double max_function = 0 ;

            // maternal fitness block
            for ( int g = 0 ; g < population[i].maternal_trnas.size() ; g ++ ) {
                if ((*population[i].maternal_trnas[g]).sequence > max_function){
                    max_function = (*population[i].maternal_trnas[g]).sequence ;
                }
            }

            // paternal fitness block
            for ( int g = 0 ; g < population[i].paternal_trnas.size() ; g ++ ) {
                if ((*population[i].paternal_trnas[g]).sequence > max_function){
                    max_function = (*population[i].paternal_trnas[g]).sequence ;
                }
            }
            fitness[i] = max_function ;
        }
    }


    //////////////////////////////////
    //////////////////////////////////
    ///////// REGULAR METHOD ///////// 
    //////////////////////////////////
    //////////////////////////////////
    //
    // tRNA MODEL BASED ON OBSERVED DATA:
    // fitness = fitness_seq * fitness_exp
    // fitness_seq = 1 - e^((options.lambda_seq) * x)
    // fitness_exp = 1 - e^(-10.5 * x)

    else {
        for ( int i = 0 ; i < population.size() ; i ++ ) {
            double max_function = 0.0 ;

            // maternal fitness block
            for ( int g = 0 ; g < population[i].maternal_trnas.size() ; g ++ ) {
                if ( gsl_ran_bernoulli( rng, (*population[i].maternal_trnas[g]).somatic ) ) {
                    int rand_index = rand() % mutation_penalties.size() ;
                    if ((1.0 - exp((options.lambda_seq) * ((*population[i].maternal_trnas[g]).sequence + mutation_penalties[rand_index]))) * (1.0 - exp(-10.5 * (*population[i].maternal_trnas[g]).expression)) > max_function ){
                        max_function = (1.0 - exp((options.lambda_seq) * ((*population[i].maternal_trnas[g]).sequence + mutation_penalties[rand_index]))) * (1.0 - exp(-10.5 * (*population[i].maternal_trnas[g]).expression)) ;
                    } 
                }
                else {
                    if ((1.0 - exp((options.lambda_seq) * ((*population[i].maternal_trnas[g]).sequence))) * (1.0 - exp(-10.5 * (*population[i].maternal_trnas[g]).expression)) > max_function){
                        max_function = (1.0 - exp((options.lambda_seq) * ((*population[i].maternal_trnas[g]).sequence))) * (1.0 - exp(-10.5 * (*population[i].maternal_trnas[g]).expression)) ;
                    }
                }

            }

            // paternal fitness block
            for ( int g = 0 ; g < population[i].paternal_trnas.size() ; g ++ ) {
                if ( gsl_ran_bernoulli( rng, (*population[i].paternal_trnas[g]).somatic ) ) {
                    int rand_index = rand() % mutation_penalties.size() ;
                    if ((1.0 - exp((options.lambda_seq) * ((*population[i].paternal_trnas[g]).sequence + mutation_penalties[rand_index]))) * (1.0 - exp(-10.5 * (*population[i].paternal_trnas[g]).expression)) > max_function ){
                        max_function = (1.0 - exp((options.lambda_seq) * ((*population[i].paternal_trnas[g]).sequence + mutation_penalties[rand_index]))) * (1.0 - exp(-10.5 * (*population[i].paternal_trnas[g]).expression)) ;
                    }
                }
                else {
                    if ((1.0 - exp((options.lambda_seq) * ((*population[i].paternal_trnas[g]).sequence))) * (1.0 - exp(-10.5 * (*population[i].paternal_trnas[g]).expression)) > max_function){
                        max_function = (1.0 - exp((options.lambda_seq) * ((*population[i].paternal_trnas[g]).sequence))) * (1.0 - exp(-10.5 * (*population[i].paternal_trnas[g]).expression)) ;
                    }
                }
            }

            fitness[i] = max_function ;
            //cout << max_function << endl ;
        }
    }
}

#endif