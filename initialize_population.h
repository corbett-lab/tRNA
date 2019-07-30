#ifndef __INITIALIZE_POPULATION_H
#define __INITIALIZE_POPULATION_H


void initialize_population( cmd_line &options, vector<gene*> &trna_bank, int &trna_counter ) {


    ///////////////
    /// MODEL 1 ///
    ///////////////
    // model1: two genes, equal function, equal mut rates, no dup/del, no somatic

    if ( options.model1 == true ){

        // BOTH GENES THE SAME:
        for ( int t = 0 ; t < 2 ; ++t ) { 
            gene* new_trna = ::new gene ; 
            /// map position of the initial tRNA 
            // make it so that it isn't at either end of the chromosome!
            new_trna->locus = ( options.map_length * 0.1 ) + ( gsl_rng_uniform( rng ) * ( options.map_length * 0.85 ) ) ;
            new_trna->sequence = 1.0 ; 
            new_trna->expression = 1.0 ;
            new_trna->progenitor = 0 ;
            new_trna->birth = 0 ;
            new_trna->muts = 0 ;
            new_trna->genotype = "A0" ;
            // new_trna->frequency.push_back( 0 ) ;
            /// mutation rates -- no somatic under model 1
            new_trna->somatic = 0 ;
            new_trna->germline = options.germline_rate ;  
            new_trna->birth_mode = 'f' ;
            // store full info in our vector of trnas
            trna_counter ++ ;
            new_trna->name = trna_counter ;
            trna_bank.push_back( new_trna ) ; 
        }
    }

    ///////////////
    /// MODEL 2 ///
    ///////////////
    // model2: two genes, one with lower function but also lower mut rates, no dup/del, no somatic

    else if ( options.model2 == true ){

        // FIRST GENE:
        gene* trna1 = ::new gene ; 
        /// map position of the initial tRNA 
        // make it so that it isn't at either end of the chromosome!
        trna1->locus = ( options.map_length * 0.1 ) + ( gsl_rng_uniform( rng ) * ( options.map_length * 0.85 ) ) ;
        /// first gene has higher function but also higher mutation rate (will reach an equilibrium)
        trna1->sequence = 1.0 ; 
        trna1->expression = 1.0 ;
        trna1->progenitor = 0 ;
        trna1->birth = 0 ;
        trna1->muts = 0 ;
        trna1->genotype = "A0" ;
        // trna1->frequency.push_back( 0 ) ;
        /// mutation rates -- no somatic under model 2
        trna1->somatic = 0 ;
        trna1->germline = options.germline_rate ;  
        trna1->birth_mode = 'f' ;
        // store full info in our vector of trnas
        trna_counter ++ ;
        trna1->name = trna_counter ;
        trna_bank.push_back( trna1 ) ; 

        // SECOND GENE:
        gene* trna2 = ::new gene ; 
        /// map position of the initial tRNA 
        // make it so that it isn't at either end of the chromosome!
        trna2->locus = ( options.map_length * 0.1 ) + ( gsl_rng_uniform( rng ) * ( options.map_length * 0.85 ) ) ;
        /// second gene has lower function but also lower mutation rate (will reach an equilibrium)
        trna2->sequence = 0.9 ; 
        trna2->expression = 1.0 ;
        trna2->progenitor = 0 ;
        trna2->birth = 0 ;
        trna2->muts = 0 ;
        trna2->genotype = "A0" ;
        // trna2->frequency.push_back( 0 ) ;
        /// mutation rates -- no somatic under model 2
        trna2->somatic = 0 ;
        trna2->germline = options.germline_rate / 10 ;  
        trna2->birth_mode = 'f' ;
        // store full info in our vector of trnas
        trna_counter ++ ;
        trna2->name = trna_counter ;
        trna_bank.push_back( trna2 ) ; 
    }

    ///////////////
    /// MODEL 4 ///
    ///////////////
    // model4: arbitrary number of genes with some germline rate lower than deverr rate
    // deverr is separate from somatic because of the way it is modeled in the paper as direct penalty off fitness

    else if ( options.model4 == true ){

        // ALL GENES THE SAME:
        for ( int t = 0 ; t < options.model4_count ; ++t ) { 
            gene* new_trna = ::new gene ; 
            /// map position of the initial tRNA 
            // make it so that it isn't at either end of the chromosome!
            new_trna->locus = ( options.map_length * 0.1 ) + ( gsl_rng_uniform( rng ) * ( options.map_length * 0.85 ) ) ;
            new_trna->sequence = 1.0 ; 
            new_trna->expression = 1.0 ;
            new_trna->progenitor = 0 ;
            new_trna->birth = 0 ;
            new_trna->muts = 0 ;
            new_trna->genotype = "A0" ;
            // new_trna->frequency.push_back( 0 ) ;
            /// mutation rates -- somatic modeled differently under model 4 so set to 0 here
            new_trna->somatic = 0 ;
            new_trna->germline = options.germline_rate ; 
            new_trna->birth_mode = 'f' ; 
            // store full info in our vector of trnas
            trna_counter ++ ;
            new_trna->name = trna_counter ;
            trna_bank.push_back( new_trna ) ; 
        }
    }

    //////////////////////////////////
    //////////////////////////////////
    ///////// REGULAR METHOD ///////// 
    //////////////////////////////////
    //////////////////////////////////

    else {
        for ( int t = 0 ; t < options.start_count ; ++t ) { 

            gene* new_trna = ::new gene ; 

            /// map position of the initial tRNA 
            // make it so that it isn't at either end of the chromosome!
            new_trna->locus = ( options.map_length * 0.1 ) + ( gsl_rng_uniform( rng ) * ( options.map_length * 0.85 ) ) ;

            /// starting conditions:
            new_trna->sequence = 1.0 ; 
            new_trna->expression = 1.0 ;
            new_trna->progenitor = 0 ;
            new_trna->birth = 0 ;
            new_trna->muts = 0 ;
            new_trna->genotype = "A0" ;
            // new_trna->frequency.push_back( 0 ) ;

            /// mutation rates :
            new_trna->somatic = options.somatic_rate * ((options.somatic_coefficient * (pow(new_trna->expression, 0.7415))) + 1.3932) ;
            new_trna->germline = options.germline_rate * ((11.8898 * (pow(new_trna->expression, 0.7415))) + 1.3932) ;  
            new_trna->birth_mode = 'f' ;

            // store full info in our vector of trnas
            trna_counter ++ ;
            new_trna->name = trna_counter ;
            trna_bank.push_back( new_trna ) ; 
        }
    }



    //////////////////////
    /// ADD PSEUDOGENE ///
    //////////////////////

    if ( options.pseudogene == true ){
        for ( int t = 0 ; t < options.start_count ; ++t ) { 

            gene* new_trna = ::new gene ; 

            /// map position of the initial tRNA 
            // make it so that it isn't at either end of the chromosome!
            new_trna->locus = ( options.map_length * 0.1 ) + ( gsl_rng_uniform( rng ) * ( options.map_length * 0.85 ) ) ;

            /// this should probably be defined by some starting conditions
            new_trna->sequence = 0.0 ; 
            new_trna->expression = 0.0 ;
            new_trna->progenitor = 0 ;
            new_trna->birth = 0 ;
            new_trna->muts = 0 ;
            new_trna->genotype = "x" ;
            // new_trna->frequency.push_back( 0 ) ;

            /// mutation rates 
            new_trna->somatic = options.somatic_rate * ((options.somatic_coefficient * (pow(new_trna->expression, 0.7415))) + 1.3932) ;
            new_trna->germline = options.germline_rate * ((11.8898 * (pow(new_trna->expression, 0.7415))) + 1.3932) ;
            new_trna->birth_mode = 'f' ;  

            // store full info in our vector of trnas
            trna_counter ++ ;
            new_trna->name = trna_counter ;
            trna_bank.push_back( new_trna ) ; 
        }
    }
}

#endif
