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

    // flag to use gaussian fitness function (default) or exponential
    string fitness_func ;

    // fitness redundant fitting parameters
    float min_fitness ;

    // fitness gaussian fitting parameters
    float fitness_mean ;
    float fitness_sd ;

    // fitness exponential fitting parameters
    float fitness_lambda ;
    
    /// mutation rates
    double germline_rate ;
    double somatic_rate ;

    // gamma / desroy model
    bool dual_rates ;
    double gamma_shape ;
    double gamma_scale ;
    double prop_destroy ;
    
    /// duplicate and deletion rates
    double deletion_rate ;
    double duplication_rate ;

    /// somatic cnv
    double somatic_del ;
    double somatic_dup ;

    /// multiple segmental duplications
    double segdup_exp ;

    /// somatic dup as a fraction of somatic del
    double somatic_dup_mult ;

    // start with pseudogene or no
    bool pseudogene ;
    
    /// proportion of duplications that are local
    double prob_cluster ; 

    /// rate at which gene conversion events occur
    double gene_conversion_rate ; 

    // output lifespans for all tRNAs
    bool output_lifespans ;

    // output sample data
    bool sample ;

    // sample whole population
    bool sample_all ;

    // sample how many generations
    int sampling_frequency ;

    // sample how many individuals each time
    int sampling_count ;

    /// rng seed
    int seed ;

    /// map length 
    double map_length ;

    // print every __ generations
    int print_count ;

    // quiet mode (don't print all individual tRNAs every time)
    bool quiet ;

    // coefficient for relating tRNA gene expression to somatic mutation rate
    double somatic_coefficient ;

    // max number of mutations a tRNA is allowed to have before it is no longer a tRNA
    int max_mutations ;

    // for running many simulations and keeping track of output
    int run_num ;

    // number of generations to be used for burn_in
    // won't start printing results until after this generation
    int burn_in ;

    // amount 
    int scaling_factor ;

    // path to fitness values for each number of mutations files
    string path ;

    // path to demography file (user must specify whole thing otherwise will assume none)
    string demography ;

    // mutation pathways bool
    bool mutation_pathways ;

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
    n = 10000 ;
    generations = 10000000 ;
    scaling_factor = 1 ;

    ////////////////////
    ///// PRIMATES /////
    ////////////////////
    //
    // germline: baseline is 1.45e-8 (narasimhan et al) +- 0.05*1e-8 [1.4e-8..1.5e-8] per base pair per generation
    //              but we can scale this for each tRNA by multiplying by 70: [9.8e-7..1.05e-6]
    //              use 1e-6
    //
    // somatic:  baseline is 2.8e-7 (milholland et al) 
    //              scaled the same way: 1.96e-5
    //
    // deletion: calculated using DGV data, 95% confidence interval = [0.00709368, 0.01064052] for theta
    //              divide by 4*Ne (40,000) to get [1.77342e-7, 2.66013e-7]
    //              use 2e-7
    // duplication: calculated using DGV data, 95% confidence interval = [0.00771959, 0.01074484] for theta
    //              divide by 4*Ne (40,000) to get [1.9298975e-7, 2.68621e-7]
    //              use 2e-7
    //
    //
    /////////////////////
    ///// MOUSE/RAT /////
    /////////////////////
    //
    // germline: baseline is 5.3e-9 (milholland et al)
    //              * 70 = 3.7e-7
    //
    // somatic:  baseline is 4.4e-7 (milholland et al) 
    //              scaled the same way: 3.08e-5
    //
    //
    //

    // TO FIT TO:
    //      - extremely low dup/del rates per million years (0.008, 0.009, etc.)
    //      - a proper fit would mean very little change overall over the course of the simulation from the start
    //
    // TO DO:
    //      - identify parameter sets that result in steady states with little to no variation
    //      - should end up being slightly different depending on how much expression you really want
    //      -- vary start counts and try both fitness functions? parameters unlikely to vary much by gene family 
    //      -- so maybe just min_fitness or even just variation of the model function??

    //
    // bergman sees ~1 / mil years in drosophila. we have ~90 species-specific tRNAs in humans over 7 million years (maybe overestimate)
    // this is a proxy for deletion and duplication rates though
    // how many duplications did we lose quickly? 
    // look up some more studies on estimates and see what we can get

    // for now, let's do 1e-6 to 1e-4 for both

    germline_rate = 1e-6 ;
    somatic_rate = 1.96e-5 ;
    deletion_rate = 3.7e-6 ; 
    duplication_rate = 3.7e-6 ; 
    somatic_del = 3.7e-6 ;

    // gamma model:
    dual_rates = false ;
    gamma_shape = 0.2 ;
    gamma_scale = 0.035 ;
    prop_destroy = 0.5 ;

    /// somatic deletion rate is almost certainly higher than somatic del
    /// so let's map somatic dup as some fraction of somatic del (from 0.0 to 1.0x sdel)
    somatic_dup_mult = 1.0 ;
    segdup_exp = 1/5763.17444  ;

    gene_conversion_rate = 2.5e-7 ;

    prob_cluster = 0.6 ; 

    sampling_frequency = 10000 ;
    sampling_count = 10 ;

    sample = false ;
    sample_all = false ;

    // mutation pathways flag
    mutation_pathways = false ;

    // quiet flag
    quiet = false ;

    // coefficient for relating tRNA gene expression to somatic mutation rate
    // default is the same as germline
    somatic_coefficient = 11.8898 ;

    // maximum number of SNPs a tRNA is allowed to have before it is no longer a tRNA gene
    max_mutations = 7 ;

    // fitness fitting parameters -- redundant
    min_fitness = 4.0 ;

    // fitness fitting parameters -- gaussian
    fitness_mean = 10.52110756 ;
    fitness_sd = fitness_mean/4.0 ;

    // fitness fitting parameters -- exponential
    fitness_lambda = -15.0 ;

    fitness_func = "redundant" ;

    seed = time(NULL) ;
    start_count = 1 ; 
    pseudogene = false ;
    map_length = 30 ; 
    print_count = 1000 ;
    burn_in = 50000 ;
    path = "/Users/Bryan/Desktop/simulator_github/" ;
    demography = "" ;
    // output_frequencies = false ;
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
        if ( strcmp(argv[i],"--dual-rates") == 0 ) {
            dual_rates = true ;
        }
        if ( strcmp(argv[i],"--gamma-shape") == 0 ) {
            gamma_shape = atof(argv[++i]) ;
        }
        if ( strcmp(argv[i],"--gamma-scale") == 0 ) {
            gamma_scale = atof(argv[++i]) ;
        }
        if ( strcmp(argv[i],"--pd") == 0 ) {
            prop_destroy = atof(argv[++i]) ;
        }
        if ( strcmp(argv[i],"--del") == 0 ) { 
            deletion_rate = atof(argv[++i]) ;
        }
        if ( strcmp(argv[i],"--dup") == 0 ) { 
            duplication_rate = atof(argv[++i]) ;
        }
        if ( strcmp(argv[i],"--sdel") == 0 ) { 
            somatic_del = atof(argv[++i]) ;
        }
        if ( strcmp(argv[i],"--sdup") == 0 ) { 
            somatic_dup_mult = atof(argv[++i]) ;
        }
        if ( strcmp(argv[i],"--seed") == 0 ) {
            seed = atoi(argv[++i]) ;
        }
        if ( strcmp(argv[i],"--scale") == 0 ) {
            scaling_factor = atoi(argv[++i]) ;
        }
        if ( strcmp(argv[i],"--start") == 0 ) {
            start_count = atoi(argv[++i]) ;
        }
        if ( strcmp(argv[i],"--pseudo") == 0 ) {
            pseudogene = true ;
        }
        if ( strcmp(argv[i],"--local") == 0 ) {
            prob_cluster = atof(argv[++i]) ;
        }
        if ( strcmp(argv[i],"--geneconv") == 0 ) {
            gene_conversion_rate = atof(argv[++i]) ;
        }
        if ( strcmp(argv[i],"--sample") == 0 ) {
            sample = true ;
        }
        if ( strcmp(argv[i],"--sample-all") == 0 ) {
            sample_all = true ;
        }
        if ( strcmp(argv[i],"--sample-freq") == 0 ) {
            sampling_frequency = atof(argv[++i]) ;
        }
        if ( strcmp(argv[i],"--sample-count") == 0 ) {
            sampling_count = atof(argv[++i]) ;
        }
        if ( strcmp(argv[i],"--min-fitness") == 0 ) {
            min_fitness = atof(argv[++i]) ;
        }
        if ( strcmp(argv[i],"--fitmean") == 0 ) {
            fitness_mean = atof(argv[++i]) ;
        }
        if ( strcmp(argv[i],"--fitsd") == 0 ) {
            fitness_sd = atof(argv[++i]) ;
        }
        if ( strcmp(argv[i],"--fitlambda") == 0 ) {
            fitness_lambda = atof(argv[++i]) ;
        }
        if ( strcmp(argv[i],"-m") == 0 ) {
            map_length = atof(argv[++i]) ;
        }
        if ( strcmp(argv[i],"--print") == 0 ) {
            print_count = atoi(argv[++i]) ;
        }
        if ( strcmp(argv[i],"--quiet") == 0 ) {
            quiet = true ;
        }
        if ( strcmp(argv[i],"--somatic-coefficient") == 0 ) {
            somatic_coefficient = atof(argv[++i]) ;
        }
        if ( strcmp(argv[i],"--max-mutations") == 0 ) {
            max_mutations = atof(argv[++i]) ;
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
        if (strcmp(argv[i],"--demography") == 0 ) {
            demography = argv[++i] ;
        }
        if (strcmp(argv[i],"--function") == 0 ) {
            fitness_func = argv[++i] ;
        }
        if ( strcmp(argv[i],"--output-lifespans") == 0 ) {
            output_lifespans = true ;
        }
        if ( strcmp(argv[i],"--mutation-pathways") == 0 ) {
            mutation_pathways = true ;
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
            gene_conversion_rate = 0 ;
            fitness_func = "model" ;
        }

        // model2: two genes, one with lower function but also lower mut rates, no dup/del, no somatic
        if ( strcmp(argv[i],"--model2") == 0 ) {
            model2 = true ;
            somatic_rate = 0 ;
            duplication_rate = 0 ;
            deletion_rate = 0 ;
            start_count = 2 ;
            gene_conversion_rate = 0 ;
            fitness_func = "model" ;
        }

        // model4: 
        if ( strcmp(argv[i],"--model4") == 0 ) {
            model4 = true ;
            somatic_rate = 0 ;
            duplication_rate = 0 ;
            deletion_rate = 0 ;
            start_count = 2 ;
            gene_conversion_rate = 0 ;
            fitness_func = "model" ;
        }
        if ( strcmp(argv[i],"--model4-count") == 0 ) {
            model4_count = atoi(argv[++i]) ;
        }
        if ( strcmp(argv[i],"--model4-deverr") == 0 ) {
            model4_deverr = atof(argv[++i]) ;
        }
    }

    somatic_dup = somatic_dup_mult * somatic_del ;

    if ( scaling_factor != 1 ){
        germline_rate = germline_rate * scaling_factor ;
        somatic_rate = somatic_rate * scaling_factor ;
        deletion_rate = deletion_rate * scaling_factor ;
        duplication_rate = duplication_rate * scaling_factor ;
        somatic_del = somatic_del * scaling_factor ;
        somatic_dup = somatic_dup * scaling_factor ;
        gene_conversion_rate = gene_conversion_rate * scaling_factor ;
        map_length = map_length / scaling_factor ;
        n = n / scaling_factor ;
        generations = generations / scaling_factor ;
        sampling_frequency = sampling_frequency / scaling_factor ;
        sampling_count = sampling_count / scaling_factor ;
        burn_in = burn_in / scaling_factor ;
    }
    
    return ;
}


#endif
