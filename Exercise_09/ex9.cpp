#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>

#include "city.h"

using namespace std;

int main() {
    int evolve_time = pow(10,3);

    Travel_Salesman popolazione = create_population();

    popolazione.order();

    for(int i{}; i<evolve_time; i++){
        popolazione.evolve();
        if(i%100 ==0) cout << i/10 <<"%" <<endl;
    }
    
    popolazione.print_best("prova.dat");

    vector<double> costi = popolazione.get_ave_fit();
    ofstream A;
    A.open("average_fit_data/average_fit_square.dat");
    for(unsigned int i{}; i<costi.size(); i++)
        A << i << "," << costi[i] << endl;

    vector<double> best_cost = popolazione.get_best_fit();
    ofstream B;
    B.open("best_fit_data/b_fit_square.dat");
    for(unsigned int i{}; i<best_cost.size(); i++)
        B << i << "," << best_cost[i] << endl;

    cout << "\n\nprinted\n\n";
    A.close();
    return 0;
}