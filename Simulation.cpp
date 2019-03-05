/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Simulation.cpp
 * Author: jcasaletto
 * 
 * Created on March 2, 2019, 2:31 PM
 */



/// standard headers
#include <string>
#include <string.h>
#include <time.h>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <math.h>
#include <set>
#include <utility>
#include <list>
#include <random>
#include <iterator>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
using namespace std ;

/// declare our gsl
const gsl_rng *rng ;


#include "CommandLine.h"
#include "Simulation.h"
#include "Individual.h"


Simulation::Simulation(map<string, Population>& m_of_p) : populations(move(m_of_p)) {
}

Simulation::Simulation(map<string, Population>& m_of_p, string n) : populations(move(m_of_p)), name(n) {
    this->name = n;
}

Simulation::Simulation(map<string, Population>& m_of_p, string n, list<Gene*>& b) : populations(move(m_of_p)), name(n), trna_bank(b) {
}

Simulation::Simulation(map<string, Population>& m_of_p, string n, list<Gene*>& b, CommandLine& o) : populations(move(m_of_p)), name(n), trna_bank(move(b)), options(move(o)) {
    this->options = o;
}

/*Simulation::Simulation(const Simulation& orig) {
} */

Simulation::~Simulation() {
}

/*map<string, Population>& Simulation::getPopulations() {
    return this->populations;
}*/

Population& Simulation::getPopulation(string s) {
    return this->populations.at(s);
}

string& Simulation::getName() {
    return this->name;
}

void Simulation::addPopulation(string s , Population p) {
    //this->populations.insert(make_pair(s, p));
    this->populations[s] = p;
    cout << "just added population " << this->populations.at("s").getName() << endl;
}


void Simulation::storeInfo(Gene* g, int &counter) {
    counter ++ ;
    g->setName(counter) ;
    //g->name = counter;

    this->trna_bank.push_back( g ) ;
}

void Simulation::initializeGene(CommandLine &options, Gene* g, double function, int neighborhood, int progenitor, int birth, double somatic, double germline, const gsl_rng rng) {
        g->setLocus(( options.map_length * 0.2 ) + ( gsl_rng_uniform( &rng ) * ( options.map_length * 0.6 ) )) ;
        g->setSequence(function);
        g->setExpression(neighborhood);
        g->setProgenitor(progenitor);
        g->setBirth(birth);
        g->setSomatic(somatic);
        g->setGermline(germline);
        g->getFrequency().push_back(0);

}

void Simulation::initialize( CommandLine &options, int &trna_counter, const gsl_rng rng ) {

    //////////////////////
    /// REGULAR METHOD ///
    //////////////////////

    for ( int t = 0 ; t < options.start_count ; t ++ ) { 

        Gene* new_trna = new Gene ;
        initializeGene(options, new_trna, 1, 1, 0, 0, options.somatic_rate, options.germline_rate, rng);
        // store full info in our vector of trnas
        storeInfo(new_trna, trna_counter);

    }

    ///////////////
    /// MODEL 1 ///
    ///////////////
    // model1: two genes, equal function, equal mut rates, no dup/del, no somatic

    if ( options.model1 == true ){
        for ( int t = 0 ; t < options.start_count ; t ++ ) { 

            Gene* new_trna =  new Gene ;
            /// second gene is identical to the first except higher mutation rate (should go extinct)
            initializeGene(options, new_trna, 1, 1, 0, 0, 0.0, options.germline_rate, rng);
            // store full info in our vector of trnas
            storeInfo(new_trna, trna_counter);

        }
    }

    ///////////////
    /// MODEL 2 ///
    ///////////////
    // model2: two genes, one with lower function but also lower mut rates, no dup/del, no somatic

    else if ( options.model2 == true ){

        for ( int t = 0 ; t < options.start_count ; t ++ ) { 

            Gene* new_trna = new Gene;
            /// second gene has lower function but also lower mutation rate (will reach an equilibrium)
            initializeGene(options, new_trna, 0.8, 1, 0, 0, 0.0, options.germline_rate/100, rng);
            // store full info in our vector of trnas
            storeInfo(new_trna, trna_counter);


        }
    }

    ///////////////
    /// MODEL 4 ///
    ///////////////
    // model4: arbitrary number of genes with some germline rate lower than deverr rate
    // deverr is separate from somatic because of the way it is modeled in the paper as direct penalty off fitness

    else if ( options.model4 == true ){
        for ( int t = 0 ; t < options.model4_count ; t ++ ) { 

            Gene* new_trna =  new Gene ;
            /// second gene has lower function but also lower mutation rate (will reach an equilibrium)
            initializeGene(options, new_trna, 1, 1, 0, 0, 0.0, options.germline_rate, rng);
            // store full info in our vector of trnas
            storeInfo(new_trna, trna_counter);

        }
    }



    //////////////////////
    /// ADD PSEUDOGENE ///
    //////////////////////

    if ( options.pseudogene == true ){
        for ( int t = 0 ; t < options.start_count ; t ++ ) { 

            Gene* new_trna =  new Gene ;
            /// this should probably be defined by some starting conditions
            initializeGene(options, new_trna, 0, 0, 0, 0, options.somatic_rate / 10, options.germline_rate / 10, rng);
            // store full info in our vector of trnas
            storeInfo(new_trna, trna_counter);
        }
    }

}

void Simulation::transmit_chromosome ( Individual &parent, vector<Gene*> &new_chromosome, const gsl_rng rng ) {

    /// go through by position 
    int it_mom = 0 ; 
    int it_dad = 0 ; 

    // if mom == 1, we are currently on mom's chromosome
    float position = 0 ; 
    bool mom = gsl_ran_bernoulli( &rng, 0.5 ) ; 

    /// while the positiion is lt both 
    while ( it_mom < parent.getMaternal_trnas().size() && it_dad < parent.getPaternal_trnas().size() ) {

        /// both have it at the same site 
        if ( parent.getMaternal_trnas()[it_mom]->getLocus() == parent.getPaternal_trnas()[it_dad]->getLocus() ) {

            if ( gsl_ran_bernoulli( &rng, ( 1 - exp( -2*parent.getMaternal_trnas()[it_mom]->getLocus() - position ) ) / 2 ) ) mom = !mom ;
            if ( mom ) { 
                new_chromosome.push_back( parent.getMaternal_trnas()[it_mom] ) ;
            }
            else {
                new_chromosome.push_back( parent.getPaternal_trnas()[it_dad] ) ;
            }
            position = parent.getMaternal_trnas()[it_mom]->getLocus() ;
            it_mom ++ ;
            it_dad ++ ; 
        } 
        /// if it's dad that's next
        else if ( parent.getMaternal_trnas()[it_mom]->getLocus() > parent.getPaternal_trnas()[it_dad]->getLocus() ) {
            if ( gsl_ran_bernoulli( &rng, ( 1 - exp( -2*parent.getPaternal_trnas()[it_dad]->getLocus() - position ) ) / 2 ) ) mom = !mom ;
            if ( !mom ) { 
                new_chromosome.push_back( parent.getPaternal_trnas()[it_dad] ) ;
            }
            position = parent.getPaternal_trnas()[it_dad]->getLocus() ;
            it_dad ++ ; 
        }
        /// otherwise, mom is next 
        else  { 
            if ( gsl_ran_bernoulli( &rng, ( 1 - exp( -2*parent.getMaternal_trnas()[it_mom]->getLocus() - position ) ) / 2 ) ) mom = !mom ;
            if ( mom ) { 
                new_chromosome.push_back( parent.getMaternal_trnas()[it_mom] ) ;
            }
            position = parent.getMaternal_trnas()[it_mom]->getLocus() ;

            it_mom ++ ; 
        }
    }

    /// if it's dad that's last
    while ( it_dad < parent.getPaternal_trnas().size() ) {
        if ( gsl_ran_bernoulli( &rng, ( 1 - exp( -2*parent.getPaternal_trnas()[it_dad]->getLocus() - position ) ) / 2 ) ) mom = !mom ;
        if ( !mom ) { 
            new_chromosome.push_back( parent.getPaternal_trnas()[it_dad] ) ;
        }
        position = parent.getPaternal_trnas()[it_dad]->getLocus() ;
        it_dad ++ ; 
    }
    /// otherwise, mom is next 
    while ( it_mom < parent.getMaternal_trnas().size() ) {
        if ( gsl_ran_bernoulli( &rng, ( 1 - exp( -2*parent.getMaternal_trnas()[it_mom]->getLocus() - position ) ) / 2 ) ) mom = !mom ;
        if ( mom ) { 
            new_chromosome.push_back( parent.getMaternal_trnas()[it_mom] ) ;
        }
        position = parent.getMaternal_trnas()[it_mom]->getLocus() ;
        it_mom ++ ; 
    }

}


void Simulation::reproduce(const double* fitness, const gsl_rng rng) {
    
    // populate parent multinomial samplings
    gsl_ran_discrete_t *g = gsl_ran_discrete_preproc( this->getPopulation("current").getIndividuals().size(), fitness ) ;

    /// iterate through all individuals and draw parents + recomb
    for ( int i = 0 ; i < this->getPopulation("current").getIndividuals().size() ; i ++ ) {
    
        /// new individual
        Individual new_ind = Individual();
        
        /// grab mom and dad
        int mom = gsl_ran_discrete( &rng, g ) ;
        int dad = gsl_ran_discrete( &rng, g ) ;

        // NOW sort, and by locus, so this might fix the recombination issue:
        vector<Gene*> momTrnas = this->getPopulation("current").getIndividuals()[i].getMaternal_trnas();
        vector<Gene*> dadTrnas = this->getPopulation("current").getIndividuals()[i].getPaternal_trnas();
        std::sort(momTrnas.begin() , momTrnas.end());
        std::sort(dadTrnas.begin() , dadTrnas.end());
        
               
        // get their chromosomes 
       /* transmit_chromosome( *population[mom], new_ind.getMaternal_trnas() ) ;
        *
        transmit_chromosome( *population[dad], new_ind.getPaternal_trnas() ) ; */

        transmit_chromosome( this->getPopulation("current").getIndividuals()[mom], new_ind.getMaternal_trnas() , rng) ;
        transmit_chromosome( this->getPopulation("current").getIndividuals()[dad], new_ind.getPaternal_trnas() , rng) ;

        // swap individual in
        //swap( *new_population[i], new_ind ) ;
        swap( this->getPopulation("new").getIndividuals()[i], new_ind ) ;
        

    }

    // de-allocate space for multinomial samplings
    if (g) {
        gsl_ran_discrete_free( g ) ;
    }
}


void Simulation::myswap(vector<Individual>& p1, vector<Individual>& p2) {
	vector<Individual> temp = p1;
	p1 = p2;
	p2 = temp;
}
 

void Simulation::print_stats ( double fitness[],  int g, list<Gene*> &trna_bank, list<float> &trna_lifespans, CommandLine &options ) {

    float trna_count = 0 ;
    float trna_pseudogenes = 0 ;
    float mean_fitness = 0.0 ;  
    map<Gene*,int> found ;
    vector<Individual> population = this->getPopulation("current").getIndividuals();

	/// go through and find all extant tRNAs
	for ( int p = 0 ; p < population.size() ; p ++ ) { 

		/// basic counts + fitness stats
		mean_fitness += fitness[p] ;
        trna_count += population[p].getMaternal_trnas().size() + population[p].getPaternal_trnas().size() ;


        /// trna specific stats
/*
		for ( int t = 0 ; t < population[p]->getMaternal_trnas().size() ; t ++ ) {
			found[population[p]->getMaternal_trnas()[t]] ++ ;
		}
		for ( int t = 0 ; t < population[p]->getPaternal_trnas().size() ; t ++ ) {
			found[population[p]->getPaternal_trnas()[t]] ++ ;
		}

*/
	}

	/// go through and delete all non-existing trnas
    /// once a tRNA dies, save its lifespan

    /// here iterating through tRNA bank. allows access to frequency (attr of tRNA in bank), but not t.second   
    for ( auto t1 = trna_bank.begin(); t1 != trna_bank.end() ; ) {
        if ( !found[(*t1)] ) {
            trna_lifespans.push_back( g - (*t1)->getBirth() ) ;
            (*t1)->getFrequency().push_back(0);

            if ( (options.output_lifespans == true) and (g > options.burn_in) ){
                std::string lifespan_log = "lifespanLog";
                lifespan_log += std::to_string(options.run_num);
                lifespan_log += ".txt" ;
                ofstream myfile;
                myfile.open( lifespan_log, fstream::app ) ;
                myfile << (*t1)->getName() << "\t" << (*t1)->getBirth() << "\t" << g << "\t" << g - (*t1)->getBirth() << "\t" << (*t1)->getProgenitor() << "\t" ;
                for ( int i = 1 ; i < (*t1)->getFrequency().size() ; i ++ ) {
                        myfile << (*t1)->getFrequency()[i] << "," ;
                }
                myfile << endl ;
                myfile.close();
            }

            delete *t1;
            t1 = trna_bank.erase(t1) ;
        }
        else {
            (*t1)->getFrequency().push_back(found[*t1]);
            t1 ++ ;
        }   
    }

    if ( g == 1 or ( g > options.burn_in and g % options.print_count == 0 ) ) {
        cout << g << "\t" ;
        cout << mean_fitness/population.size() << "\t" ; 
        cout << trna_count/population.size() << "\t" ;
        cout << trna_bank.size() ;

        map<Gene*, int>::iterator t;

       // for ( auto const& t : found ) {
        for (t = found.begin(); t != found.end(); t++) {
            if ( t->second > 0 ) {
            	Gene* bob = t->first;
                cout << "\t" << bob->getName() << "_" << bob->getLocus() << "_" <<  bob->getSequence() << "_" ;
                cout << bob->getExpression() << "_" << bob->getBirth() << "_" << bob->getProgenitor() << "_" << t->second  ;
                if (t->second > options.n*2) {
                    cout << "\t" << "ERROR\nERROR\nERROR" ;
                    // a tRNA can NOT be found more times than there are chromosomes. Exit if this happens.
                    exit(0);
                }
            }   
        }
        cout << endl ; 

    }
    
    

    // // DEBUGGING:


}



int Simulation::run() {
    
    /*for(int i=0; i< 50; i++) {
	// instantiate individual on the stack
        //simulation.getPopulation("newPop").pushback(Individual("bob", i));
        this->getPopulation("current").pushback(Individual());
        this->getPopulation("current").getIndividuals()[i].setName("bob");
        this->getPopulation("current").getIndividuals()[i].setSize(i);
        for(int j=0; j< this->getPopulation("current").getIndividuals()[i].getSize(); j++){
            this->getPopulation("current").getIndividuals()[i].pushback(new Gene(j));
            cout << this->getPopulation("current").getIndividuals()[i].getGenes().back()->getName() << endl;
        }
        cout << this->getPopulation("current").getIndividuals()[i].getName() << "\n";
        cout << this->getPopulation("current").getIndividuals()[i].getSize() << "\n";     
    }*/
    
    
        //// for calculating runtime of code
    clock_t tStart = clock() ;

    /// trna counter to give new name to every trna that arises
    int trna_counter = 0 ;

    //// read command line options


    // initialize rng for gsl lookup table
    const gsl_rng *rng = gsl_rng_alloc( gsl_rng_taus2 ) ;
    gsl_rng_set( rng, (long) options.seed ) ;
    cout << options.path + "allPenaltiesPct.txt" << endl;
    std::ifstream is(options.path + "allPenaltiesPct.txt") ;
    std::istream_iterator<double> start(is), end ;
    std::vector<double> mutation_penalties(start, end) ;

    /// list of tRNA lifespans
    list<float> trna_lifespans ;
    
    // create population of size n
    //Population population ("current", options.n ) ;
   

    // initialize trna_bank
    this->initialize(options, trna_counter, *rng);
    /// now copy to population of size n

    for ( int i = 0 ; i < this->getPopulation("current").getIndividuals().size() ; i ++ ) { 
    	for ( auto t : this->trna_bank ) { 
            this->getPopulation("current").getIndividuals()[i].getMaternal_trnas().push_back(t);
            this->getPopulation("current").getIndividuals()[i].getPaternal_trnas().push_back(t);
    	}
    } 

    
    // fitness vector
    double fitness [options.n]  ;
    
    // evolve the population forward in time
    for ( int g = 1 ; g < options.generations ; g ++ ) {
        
        /// somatic and germline mutations
        this->getPopulation("current").mutate(options, trna_bank, g, trna_counter, *rng ) ;
        
        
        /// compute fitness
        this->getPopulation("current").compute_fitness( fitness, options, *rng ) ;
        
        /// reproduce w/ fitness + recombination
        this->reproduce(fitness, *rng) ;

      
        /// swap populations
        swap( this->getPopulation("current").getIndividuals(), this->getPopulation("new").getIndividuals() ) ;

		/////// print stats
        // for first generation give user a quick reminder of what they did:
        if ( g == 1 ){
            cout << "somatic = " << options.somatic_rate ;
            cout << ", germline = " << options.germline_rate ;
            cout << ", dup = " << options.duplication_rate ;
            cout << ", del = " << options.deletion_rate << endl ;
        }

        print_stats( fitness, g, trna_bank, trna_lifespans, options ) ;
    }
/*
    if ( options.output_frequencies == true ){
        std::string frequency_log = "frequencyLog" ;
        frequency_log += std::to_string(options.run_num) ;
        frequency_log += ".txt" ;
        ofstream myfile;
        myfile.open( frequency_log, fstream::app ) ;
        for ( auto t : trna_bank ) {
            myfile << t->getName() << "\t" << t->birth << "\t" << t->progenitor << "\t" ;
            for ( int i = 1 ; i < t->getFrequency().size() ; i ++ ) {
                myfile << t->getFrequency[i] << "," ;
            }
            myfile << endl ;
        }
        myfile.close(); 
    }*/

    printf("Total time: %.2f seconds. Total generations: ", (double)(clock() - tStart)/CLOCKS_PER_SEC) ;
    cout << options.generations << "." << endl ; 
    cout << "trna bank size is " << trna_bank.size() << endl;
    return(0) ;

}


