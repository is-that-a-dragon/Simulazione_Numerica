#pragma once
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cassert>
#include <iomanip>
#include <string>

#include "random.h"

using namespace std;

void wake_up_geffrey(Random &);

class Travel_Salesman{
    public:
        Travel_Salesman() {};
        Travel_Salesman(unsigned int N, vector<int> starter, vector<vector<double>> pos, Random& geffrey) {
            N_cities = N;
            m_start = starter;
            m_city_pos = pos;
            m_geffrey = geffrey; 
        };
        ~Travel_Salesman() {};

        void check();
        void add_path(vector<int> percorso) {m_population.push_back(percorso);};
        void set_pop(vector<vector<int>> pop) {m_population = pop;};
        vector<vector<int>> copy_pop() {return m_population;};
        int get_pop() {return m_population.size();};
        int get_city(int chromo_index, int city_index) { return m_population[chromo_index][city_index]; };

        void set_crossover_prob(double P) { crossover_prob = P;};
        void set_mutation_prob(double P) { mutation_prob = P;};
        void set_elite(bool E) { m_elite = E; };
        
        double tomtom(int, int);
        vector<double> fit();
        double fit(vector<int>);

        void set_ave_fit(int);
        vector<double> get_ave_fit(){ return average_fit; };
        vector<double> get_best_fit() { return best_fit; };

        void print_best(string);
        vector<int> get_best() { return m_population[0]; };

        void order();
        vector<vector<int>> select(vector<vector<int>>);
        vector<vector<int>> crossover(vector<vector<int>>);
        vector<int> mutate(vector <int>);
        void evolve();
    
    protected:
        unsigned int N_cities;
        vector<int> m_start;
        vector< vector <double> > m_city_pos;
        vector< vector <int> > m_population;
        Random m_geffrey;
        
        double crossover_prob;
        double mutation_prob;
        bool m_elite;
        vector<double> average_fit;
        vector<double> best_fit;
};

Travel_Salesman create_population();

Travel_Salesman US_population();

void read_map(vector<double>& , vector<double>& , vector<int>&);