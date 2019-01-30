#ifndef __INITIALIZE_POPULATION_H
#define __INITIALIZE_POPULATION_H


void initialize_population( cmd_line &options, list<gene*> &trna_bank, int &trna_counter ) {

    //////////////////////
    /// REGULAR METHOD ///
    //////////////////////

    for ( int t = 0 ; t < options.start_count ; t ++ ) { 

        gene* new_trna = new gene ; 

        /// map position of the initial tRNA 
        
        // make it so that it isn't at either end of the chromosome!
        new_trna->locus = ( options.map_length * 0.2 ) + ( gsl_rng_uniform( rng ) * ( options.map_length * 0.6 ) ) ;

        /// starting conditions:
        new_trna->function = 1 ; 
        new_trna->neighborhood = 1 ;
        new_trna->progenitor = 0 ;
        new_trna->birth = 0 ;
        new_trna->frequency.push_back( 0 ) ;

        /// mutation rates :
        new_trna->somatic = options.somatic_rate ;
        new_trna->germline = options.germline_rate ;  

        // store full info in our vector of trnas
        trna_counter ++ ;
        new_trna->name = trna_counter ;
        trna_bank.push_back( new_trna ) ; 
    }

    ///////////////
    /// MODEL 1 ///
    ///////////////
    // model1: two genes, equal function, equal mut rates, no dup/del, no somatic

    if ( options.model1 == true ){
        for ( int t = 0 ; t < options.start_count ; t ++ ) { 

            gene* new_trna = new gene ; 

            /// map position of the initial tRNA 
            // make it so that it isn't at either end of the chromosome!
            new_trna->locus = ( options.map_length * 0.2 ) + ( gsl_rng_uniform( rng ) * ( options.map_length * 0.6 ) ) ;

            /// second gene is identical to the first except higher mutation rate (should go extinct)
            new_trna->function = 1 ; 
            new_trna->neighborhood = 1 ;
            new_trna->progenitor = 0 ;
            new_trna->birth = 0 ;
            new_trna->frequency.push_back( 0 ) ;

            /// mutation rates 
            new_trna->somatic = 0 ;
            new_trna->germline = options.germline_rate ;  

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

        for ( int t = 0 ; t < options.start_count ; t ++ ) { 

            gene* new_trna = new gene ; 

            /// map position of the initial tRNA 
            // make it so that it isn't at either end of the chromosome!
            new_trna->locus = ( options.map_length * 0.2 ) + ( gsl_rng_uniform( rng ) * ( options.map_length * 0.6 ) ) ;

            /// second gene has lower function but also lower mutation rate (will reach an equilibrium)
            new_trna->function = 0.8 ; 
            new_trna->neighborhood = 1 ;
            new_trna->progenitor = 0 ;
            new_trna->birth = 0 ;
            new_trna->frequency.push_back( 0 ) ;

            /// mutation rates 
            new_trna->somatic = 0 ;
            new_trna->germline = options.germline_rate / 100 ;  

            // store full info in our vector of trnas
            trna_counter ++ ;
            new_trna->name = trna_counter ;
            trna_bank.push_back( new_trna ) ; 
        }
    }

    ///////////////
    /// MODEL 4 ///
    ///////////////
    // model4: arbitrary number of genes with some germline rate lower than deverr rate
    // deverr is separate from somatic because of the way it is modeled in the paper as direct penalty off fitness

    else if ( options.model4 == true ){
        for ( int t = 0 ; t < options.model4_count ; t ++ ) { 

            gene* new_trna = new gene ; 

            /// map position of the initial tRNA 
            
            // make it so that it isn't at either end of the chromosome!
            new_trna->locus = ( options.map_length * 0.2 ) + ( gsl_rng_uniform( rng ) * ( options.map_length * 0.6 ) ) ;

            /// second gene has lower function but also lower mutation rate (will reach an equilibrium)
            new_trna->function = 1 ; 
            new_trna->neighborhood = 1 ;
            new_trna->progenitor = 0 ;
            new_trna->birth = 0 ;
            new_trna->frequency.push_back( 0 ) ;

            /// mutation rates 
            new_trna->somatic = 0 ;
            new_trna->germline = options.germline_rate ;  

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
        for ( int t = 0 ; t < options.start_count ; t ++ ) { 

            gene* new_trna = new gene ; 

            /// map position of the initial tRNA 
            
            // make it so that it isn't at either end of the chromosome!
            new_trna->locus = ( options.map_length * 0.2 ) + ( gsl_rng_uniform( rng ) * ( options.map_length * 0.6 ) ) ;

            /// this should probably be defined by some starting conditions
            new_trna->function = 0 ; 
            new_trna->neighborhood = 0 ;
            new_trna->progenitor = 0 ;
            new_trna->birth = 0 ;
            new_trna->frequency.push_back( 0 ) ;

            /// mutation rates 
            new_trna->somatic = options.somatic_rate / 10 ;
            new_trna->germline = options.germline_rate / 10 ;  

            // store full info in our vector of trnas
            trna_counter ++ ;
            new_trna->name = trna_counter ;
            trna_bank.push_back( new_trna ) ; 
        }
    }
}

#endif
