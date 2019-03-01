#ifndef __CMD_LINE_H
#define __CMD_LINE_H

/// command line information and global parameters
class cmd_line {
public:
    
    /// population size, constant for now
    int n ;
    
    /// number of generations
    int generations ;
    
    /// number of genes at start
    int start_count ;
    
    /// mutation rates
    double germline_rate ;
    double somatic_rate ;
    
    /// duplicate and deletion rates
    double deletion_rate ;
    double duplication_rate ;

    // start with pseudogene or no
    bool pseudogene ;
    
    /// probability of cluster
    double prob_cluster ; 

    /// to be used in mapping sequence to fitness
    double lambda_seq ;

    // output frequencies at the end or no
    bool output_frequencies ;

    // output lifespans for all tRNAs
    bool output_lifespans ;

    /// rng seed
    int seed ;

    /// map length 
    double map_length ;

    // print every __ generations
    int print_count ;

    // for running many simulations and keeping track of output
    int run_num ;

    // number of generations to be used for burn_in
    // won't start printing results until after this generation
    int burn_in ;

    // path to allPenaltiesPct.txt vector file
    string path ;

    // flags to replicate models from paper
    bool model1 ;
    bool model2 ;
    bool model4 ;
    int model4_count ;
    float model4_deverr ;
    
    /// read relevant information
    void read_cmd_line ( int argc, char *argv[] ) ;

} ;

void cmd_line::read_cmd_line ( int argc, char *argv[] ) {
    
    /// set parameter defaults
    n = 1000 ;
    generations = 1000000 ;

    /// for now, setting these to follow distributions of min/max:
    // germline: baseline is 1.45e-8 (narismahan et al), active tRNAs ~10x higher, current range 1e-8 - 1e-5
    // somatic: should be ~10x higher than germline (milholland et al). range for somatic will be between 1 - 100x higher than germline.
    // deletion: 
    // duplication: 
    //
    // bergman sees ~1 / mil years in drosophila. we have ~90 species-specific tRNAs in humans over 7 million years (maybe overestimate)
    // this is a proxy for deletion and duplication rates though
    // how many duplications did we lose quickly? 
    // look up some more studies on estimates and see what we can get

    // for now, let's do 1e-6 to 1e-4 for both

    germline_rate = 1e-6 ;
    somatic_rate = 1e-5 ;
    deletion_rate = 1e-5 ; // + (gsl_rng_uniform( rng ) * 100) ;
    duplication_rate = 1e-5 ; // + (gsl_rng_uniform( rng ) * 100) ;
    prob_cluster = 0.5 ; 
    lambda_seq = -5.0 ;

    seed = time(NULL) ;
    start_count = 1 ; 
    pseudogene = false ;
    map_length = 30 ; 
    print_count = 1000 ;
    burn_in = 50000 ;
    path = "" ;
    output_frequencies = false ;
    output_lifespans = false ;
    run_num = 0 ;

    /// flags to replicate models in http://ped.fas.harvard.edu/files/ped/files/nature97_0.pdf
    //
    // model1: two genes, equal function, equal mut rates, no dup/del, no somatic
    // model2: two genes, one with lower function but also lower mut rates, no dup/del, no somatic
    // model3 involves two genes performing different functions -- not redundancy, not testing this (at least for now)
    // model4: two genes with two different germline and somatic mutation rates,
    //   such that g1 < s2 and g2 < s1
    //   call somatic here separate because it's modeled as part of fitness, not part of function

    model1 = false ;
    model2 = false ;
    model4 = false ;
    model4_count = 1 ;
    model4_deverr = 1e-4 ;
    
    /// accept command line parameters
    for (int i=1; i<argc; i++) {

        if ( strcmp(argv[i],"-n") == 0 ) {
            n = atoi(argv[++i]) ;
        }
        if ( strcmp(argv[i],"-g") == 0 ) {
            generations = atoi(argv[++i]) ;
        }
        if ( strcmp(argv[i],"--ug") == 0 ) {
            germline_rate = atof(argv[++i]) ;
        }
        if ( strcmp(argv[i],"--us") == 0 ) {
            somatic_rate = atof(argv[++i]) ;
        }
        if ( strcmp(argv[i],"--del") == 0 ) { 
            deletion_rate = atof(argv[++i]) ;
        }
        if ( strcmp(argv[i],"--dup") == 0 ) { 
            duplication_rate = atof(argv[++i]) ;
        }
        if ( strcmp(argv[i],"-s") == 0 ) {
            seed = atoi(argv[++i]) ;
        }
        if ( strcmp(argv[i],"--start") == 0 ) {
            start_count = atoi(argv[++i]) ;
        }
        if ( strcmp(argv[i],"--pseudo") == 0 ) {
            pseudogene = true ;
        }
        if ( strcmp(argv[i],"-c") == 0 ) {
            prob_cluster = atof(argv[++i]) ;
        }
        if ( strcmp(argv[i],"--lambda") == 0 ) {
            lambda_seq = atof(argv[++i]) ;
        }
        if ( strcmp(argv[i],"-m") == 0 ) {
            map_length = atof(argv[++i]) ;
        }
        if ( strcmp(argv[i],"--print") == 0 ) {
            print_count = atoi(argv[++i]) ;
        }
        if (strcmp(argv[i],"-b") == 0 ) {
            burn_in = atoi(argv[++i]) ;
        }
        if (strcmp(argv[i],"--run") == 0 ) {
            run_num = atof(argv[++i]) ;
        }
        if (strcmp(argv[i],"--path") == 0 ) {
            path = argv[++i] ;
        }
        if ( strcmp(argv[i],"--output-frequencies") == 0 ) {
            output_frequencies = true ;
            std::string frequency_log = "frequencyLog" ;
            frequency_log += std::to_string(run_num) ;
            frequency_log += ".txt" ;
            ofstream myfile;
            myfile.open( frequency_log ) ;
            myfile << "tRNA\tbirth\tprogenitor\tfrequencies\n" ;
            myfile.close();
        }
        if ( strcmp(argv[i],"--output-lifespans") == 0 ) {
            output_lifespans = true ;
            std::string lifespan_log = "lifespanLog" ;
            lifespan_log += std::to_string(run_num) ;
            lifespan_log += ".txt" ;
            ofstream myfile;
            myfile.open( lifespan_log ) ;
            myfile << "tRNA\tbirth\tdeath\tlifespan\tprogenitor\tfrequencies\n" ;
            myfile.close();
        }

        /// replicating models requires changes to other parameters
        // this is last, so model flags overrule everything else!

        // model1: two genes, equal function, equal mut rates, no dup/del, no somatic
        if ( strcmp(argv[i],"--model1") == 0 ) {
            model1 = true ;
            somatic_rate = 0 ;
            duplication_rate = 0 ;
            deletion_rate = 0 ;
            start_count = 2 ;
        }

        // model2: two genes, one with lower function but also lower mut rates, no dup/del, no somatic
        if ( strcmp(argv[i],"--model2") == 0 ) {
            model2 = true ;
            somatic_rate = 0 ;
            duplication_rate = 0 ;
            deletion_rate = 0 ;
            start_count = 2 ;
        }

        // model4: 
        if ( strcmp(argv[i],"--model4") == 0 ) {
            model4 = true ;
            somatic_rate = 0 ;
            duplication_rate = 0 ;
            deletion_rate = 0 ;
            start_count = 2 ;
        }
        if ( strcmp(argv[i],"--model4-count") == 0 ) {
            model4_count = atoi(argv[++i]) ;
        }
        if ( strcmp(argv[i],"--model4-deverr") == 0 ) {
            model4_deverr = atof(argv[++i]) ;
        }

    }
    
    return ;
}


#endif
