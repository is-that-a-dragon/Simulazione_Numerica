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

void wake_up_geffrey(Random &, string);

class city{

    public:
        city() {};
        city(unsigned int indice, double x, double y){
            m_index = indice;
            coordinates.push_back(x);
            coordinates.push_back(y);
        };
        city(unsigned int indice, vector<double> coord){
            m_index = indice;
            coordinates = coord;
        }
        ~city() {};

        bool operator== (const city& rhs){
            return (this->m_index == rhs.m_index) and (this->coordinates == rhs.coordinates);
        };

        int get_index() const {return m_index;};
        int get_dimension() const {return coordinates.size();};
        double get_coord(int i) const {return coordinates[i];};
        void print_city(){
            cout << setw(10) << get_index() << setw(10) << get_coord(0) << setw(10) << get_coord(1) << endl;
        }
    
    protected:
        unsigned int m_index;
        vector<double> coordinates; 
};

double distance(city& A, city& B);

class path{
    public:
        path (){};
        path (vector<city> percorso) {
            //m_cities = percorso;
            for(unsigned int i{}; i<percorso.size(); i++)
                m_cities.push_back( percorso[i] );
            
            for(unsigned int i{}; i<m_cities.size()-1; i++)
                m_fit += distance( m_cities[i], m_cities[i+1] );
            m_fit += distance( m_cities[ m_cities.size()-1 ] , m_cities[0] ); 
        };
        ~path() {};

        int get_size() const {return m_cities.size();};
        int get_city_dim() const {return m_cities[0].get_dimension();};
        double get_city_index(int posizione) const { return m_cities[posizione].get_index(); };
        double get_city_coord(int posizione, int i) const { return m_cities[posizione].get_coord(i); };
        double get_fit() const {return m_fit;};
        void print_path();

        city get_city (int i) const { return m_cities[i]; };
        vector<city> get_cities() const { return m_cities; };
        void set_fit(){
            m_fit = 0. ;
            for(unsigned int i{}; i<m_cities.size()-1; i++)
                m_fit += distance( m_cities[i], m_cities[i+1] );
            m_fit += distance( m_cities[ m_cities.size()-1 ] , m_cities[0] );
        }
        void set_cities(vector<city> A) {
            m_cities = A;
            set_fit();
            };

    protected:
        vector<city> m_cities;
        double m_fit;
};

class Travel_Salesman{
    public:
        Travel_Salesman() {};
        Travel_Salesman(unsigned int N, vector<city> starter, Random& geffrey) {
            N_cities = N;
            m_start = starter;
            m_geffrey = geffrey; 
        };
        ~Travel_Salesman() {};

        void check();
        void print_best(string);

        void set_mutation_prob(double P) {mutation_prob = P;};
        void set_crossover_prob(double P) {crossover_prob = P;};
        void set_elite(bool E) {m_elite = E;};

        int get_pop() const {return m_population.size();};
        double get_fit(int posizione) { return m_population[posizione].get_fit(); };
        vector<city> access_start() {return m_start;};

        void add_path(path A){ m_population.push_back(A); };
        void set_path(path A, int posizione){ m_population[posizione] = A;};
               
        void path_fit(){
            for(unsigned int j{}; j < m_start.size(); j++)
                m_population[j].set_fit();
        };

        void set_ave_fit(int);
        void set_best_fit() { best_fit.push_back( get_fit(0) ); };
        vector<double> get_ave_fit() const {return average_fit;};
        vector<double> get_best_fit() const {return best_fit;};

        void order();
        vector<path> select(vector<path>);
        vector<path> crossover(vector<path>);
        path mutate(path gene);
        void evolve();
    
    protected:
        unsigned int N_cities;
        vector<city> m_start;
        vector<path> m_population;
        Random m_geffrey;
        
        double crossover_prob;
        double mutation_prob;
        bool m_elite;
        vector<double> average_fit;
        vector<double> best_fit;
};

Travel_Salesman create_population();