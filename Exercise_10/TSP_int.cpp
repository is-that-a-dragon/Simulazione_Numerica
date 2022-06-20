#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>

#include "city.h"

using namespace std;

int main() {
    int evolve_time = 1000;

    Travel_Salesman popolazione = US_population();
    for(int i{}; i<evolve_time; i++){
        popolazione.evolve();
        if(i*10%evolve_time == 0)
            cout << "\033[1m\033[37m" << "evolution at:\t" << (double) i*100./evolve_time << "%\n" << "\r\033[F";
        if(i==evolve_time-1)
            cout << "\033[1m\033[37m" << "evolution at:\t" << "100%\n" << "\033[0m";
    }
    
    popolazione.print_best("best_path_small.dat");

    vector<double> costi = popolazione.get_ave_fit();
    ofstream A;
    A.open("avefit._smalldat");
    for(unsigned int i{}; i<costi.size(); i++)
        A << i << "," << costi[i] << endl;
    A.close();

    vector<double> best_cost = popolazione.get_best_fit();
    ofstream B;
    B.open("best_fit_US_small.dat");
    for(unsigned int i{}; i<best_cost.size(); i++)
        B << i << "," << best_cost[i] << endl;
    
    cout << "\n\nprinted\n\n";

    return 0;
}