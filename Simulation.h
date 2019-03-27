/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Simulation.h
 * Author: jcasaletto
 *
 * Created on March 2, 2019, 2:31 PM
 */

#ifndef SIMULATION_H
#define SIMULATION_H

#include <vector>
#include <string>
#include <map>
#include <iostream>
#include "Population.h"
using namespace std;

class Simulation {
public:
    Simulation(map<string, Population>& m_of_p);
    Simulation(map<string, Population>& m_of_p, string n);
    Simulation(map<string, Population>& m_of_p, string n, list<Gene*>& t);
    Simulation(map<string, Population>& m_of_p, string n, list<Gene*>& t, CommandLine& o);

    //Simulation(const Simulation& orig);
    virtual ~Simulation();
    Population& getPopulation(string s);
    string& getName();
    int run();
    void addPopulation(string s, Population p);
    void initialize( CommandLine &options, int &trna_counter, const gsl_rng rng );
    void initializeGene(CommandLine &options, Gene* g, double function, int neighborhood, int progenitor, int birth, double somatic, double germline, const gsl_rng rng);
    void storeInfo(Gene* g, int &counter);
    void reproduce( const double* fitness, Population &newPopulation, const gsl_rng rng) ;
    void transmit_chromosome ( Individual &parent, vector<Gene*> &new_chromosome, const gsl_rng rng);
    void myswap(vector<Individual>& p1, vector<Individual>& p2);
    void print_stats ( double fitness[], int g, list<Gene*> &trna_bank, list<float> &trna_lifespans, CommandLine &options );


    
private:
    string name;
    map<string, Population> populations;
    list<Gene*> trna_bank;
    CommandLine options;



};

#endif /* SIMULATION_H */

