#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include <random>

#include "city.h"
#include "mpi.h"

using namespace std;

/*
    READ ME:
        as it is, this main computes the 50 US capitals
        to change that you have to:
        > change the n_cities value (line 41)
        > change the "read_map" function to read from circle_town or square_town
        > change the names of the output files so they do not overlap
*/


int main(int argc , char* argv[]){

    int size , rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size( MPI_COMM_WORLD, &size );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );

    int n_migr{20};        // migration frequency
    int n_send{10};        // migrating chromosomes
    int n_print{1};        // print best frequency
    int pop_size{1000};
    int evolve_time{200};

    double c_prob{0.8};
    double m_prob{0.3};
    bool eli = false; 

    //set variables to initialize population in other ranks
    unsigned int n_cities{50};
    vector<int> starter (n_cities); //vector of n_cities dimension filled with 0
    vector<double> x_pos (n_cities);
    vector<double> y_pos (n_cities);
    Random geffrey;
    wake_up_geffrey(geffrey);

    //rank 0 reads positions and broadcasts them to all other nodes
    //vectors are guaranteed to have contiguously allocated memory:
    //  I can treat them like C arrays by getting a pointer to the first element. 
    if(rank==0)
        read_map(x_pos, y_pos, starter);
   
    if(size>1){
        MPI_Bcast( x_pos.data(), x_pos.size(), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD);
        MPI_Bcast( y_pos.data(), y_pos.size(), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD);
        MPI_Bcast( starter.data(), starter.size(), MPI_INTEGER, 0, MPI_COMM_WORLD);
        //vector::data() returns a direct pointer to the memory array
        //          used internally by the vector to store its owned elements
    }
    //create position matrix
    vector<vector<double>> pos;
    pos.push_back(x_pos);
    pos.push_back(y_pos);

    //create first generation, each node uses a different seed
    vector<vector<int>> first_gen;
    int seed = pow(rank+1, 2);
    vector<int> copy = starter;
    for(int i{}; i < pop_size; i++){
        shuffle(copy.begin()+1, copy.end(), std::default_random_engine(seed));
        first_gen.push_back(copy);
    }

    //create population for every node, check and order
    Travel_Salesman popolazione(n_cities, starter, pos, geffrey);
    popolazione.set_pop( first_gen );
    popolazione.set_crossover_prob( c_prob );
    popolazione.set_mutation_prob( m_prob );
    popolazione.set_elite( eli );
    popolazione.check();
    popolazione.order();

    //print cost and best path to file
    ofstream out;
    vector<int> best_path = popolazione.get_best();
    out.open("silent_rank_costs/cost"+to_string(rank)+".dat");
    out << popolazione.fit( best_path ) << "\n";
    out.close();

    //actually useful appo variables
    vector<int> best_chromos( n_send*n_cities );
    vector<int> migrated_chromos( n_send*n_cities*size );
    vector<vector<int>> appo_pop;
    vector<int> chromosome(n_cities);

    //yay! evolution time!
    for(int i{}; i<evolve_time; i++){
        //print to show evolution best path fit evolution
        if(i%n_print == 0 && i!=0){
            best_path = popolazione.get_best();
            out.open("silent_rank_costs/cost"+to_string(rank)+".dat", ios::app);
            out << popolazione.fit( best_path ) << "\n";
            out.close();
        }

        popolazione.evolve();
        if(rank==0){
            if(i*20%evolve_time == 0)
                cout << "\033[1m\033[37m" << "evolution at:\t" << (double) i*100./evolve_time << "%\n" << "\r\033[F";
            if(i==evolve_time-1)
                cout << "\033[1m\033[37m" << "evolution at:\t" << "100%\n" << "\033[0m";
        }

        //send best chromosomes
        
        if(i%n_migr == 0 && i!=0 && size>1){
            popolazione.order();
            appo_pop = popolazione.copy_pop();
            //all n_send chromosomes are put in one vector (best_chromo)
            for (int j{}; j < n_send; j++){ //chromo index
                for (unsigned int h{}; h < n_cities; h++) //city index in chromo
                    best_chromos[h+j*n_cities] = popolazione.get_city(j,h);
            }
            MPI_Allgather(best_chromos.data(), n_send*n_cities, MPI_INTEGER, migrated_chromos.data(), n_send*n_cities, MPI_INTEGER, MPI_COMM_WORLD);
                        //starting address of sender //#n elem to send //type //address of receivee //#n elem received //type //communicator
            //update population: each nodes receives n_send*size chromosomes from the other nodes
            //now it substitutes its worst elements
            for (int j{}; j < n_send*size; j++) { //number of sent chromosomes
                for (unsigned int k{}; k < n_cities; k++) //city in a chromosome
                    chromosome[k] = migrated_chromos[k+j*n_cities];
                appo_pop[ pop_size-1 -j ] = chromosome;
                fill(chromosome.begin(), chromosome.end(), 0);
            }
            popolazione.set_pop( appo_pop );
        }
    }
    //print
    popolazione.order();
    out.open("rank_costs/cost"+to_string(rank)+".dat", ios::app);
    best_path = popolazione.get_best();
    out << popolazione.fit( best_path ) << "\n";
    out.close();

    string title = "rank_paths/best_chromo"+to_string(rank)+".dat";
    popolazione.print_best( title );

    MPI_Finalize();
    
    return 0;
}