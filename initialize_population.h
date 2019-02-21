#ifndef __INITIALIZE_POPULATION_H
#define __INITIALIZE_POPULATION_H

void store_info(gene* g, int &counter, list<gene*> &bank) {
    counter ++ ;
    g->name = counter ;
    bank.push_back( g ) ;

}

void initialize_gene(cmd_line &options, gene* g, double function, int neighborhood, int progenitor, int birth, double somatic, double germline) {
	g->locus = ( options.map_length * 0.2 ) + ( gsl_rng_uniform( rng ) * ( options.map_length * 0.6 ) ) ;
	g->function = function;
	g->neighborhood = neighborhood;
	g->progenitor = progenitor;
	g->birth = birth;
	g->somatic = somatic;
	g->germline = germline;
	g->frequency.push_back(0);
}

void initialize_population( cmd_line &options, list<gene*> &trna_bank, int &trna_counter ) {


    for ( int t = 0 ; t < options.start_count ; t ++ ) { 
        //////////////////////
        /// REGULAR METHOD ///
        //////////////////////
        gene* regular_new_trna = new gene ;
        initialize_gene(options, regular_new_trna, 1, 1, 0, 0, options.somatic_rate, options.germline_rate);
        store_info(regular_new_trna, trna_counter, trna_bank);


        //////////////
        /// MODELS ///
        //////////////
        gene* model_new_trna = new gene ;

        switch(options.model) {
        case 1:
        	///////////////
        	/// MODEL 1 ///
        	///////////////
        	// model1: two genes, equal function, equal mut rates, no dup/del, no somatic
        	/// second gene is identical to the first except higher mutation rate (should go extinct)
        	initialize_gene(options, model_new_trna, 1, 1, 0, 0, 0.0, options.germline_rate);
        	break;

        case 2:
        	///////////////
        	/// MODEL 2 ///
        	///////////////
        	// model2: two genes, one with lower function but also lower mut rates, no dup/del, no somatic
        	// map position of the initial tRNA
        	// make it so that it isn't at either end of the chromosome!
        	initialize_gene(options, model_new_trna, 0.8, 1, 0, 0, 0.0, options.germline_rate/100);
        	break;

        case 4:
        	///////////////
        	/// MODEL 4 ///
        	///////////////
        	// model4: arbitrary number of genes with some germline rate lower than deverr rate
        	// deverr is separate from somatic because of the way it is modeled in the paper as direct penalty off fitness
        	// make it so that it isn't at either end of the chromosome!
        	initialize_gene(options, model_new_trna, 1, 1, 0, 0, 0.0, options.germline_rate);
        	break;

        }

        // store full info in our vector of trnas
        store_info(model_new_trna, trna_counter, trna_bank);


        //////////////////////
        /// ADD PSEUDOGENE ///
        //////////////////////

        if ( options.pseudogene == true ){
        	gene* pseudogene_new_trna = new gene ;
        	initialize_gene(options, pseudogene_new_trna, 0, 0, 0, 0, options.somatic_rate / 10, options.germline_rate / 10);
        	// store full info in our vector of trnas
            store_info(pseudogene_new_trna, trna_counter, trna_bank);

        }
    }

}

#endif
