#include <iostream>
#include <string>
#include <vector>
#include "Gene.h"
#include "Individual.h"
#include "Simulation.h"
using namespace std;

int main ( int argc, char **argv ) {
    
    CommandLine options = CommandLine();
    options.read_cmd_line( argc, argv ) ;
    map<string, Population>  populations ;
    vector<Individual> currentIndividuals, newIndividuals;
    populations["current"] = Population("current", options.n, currentIndividuals);
    populations["new"] = Population("new", options.n, newIndividuals);
    list<Gene*> trna_bank;
    Simulation simulation = Simulation(populations, "test", trna_bank, options);

    simulation.run();

    return(0);

}

