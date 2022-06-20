#include "libre.h"
#include "funzioni.h"
#include "IntegralMC.h"
#include <ctime>

using namespace std;

int main(){

    Random geffrey;
    wake_up_geffrey(geffrey);
    IntegralMC calculator(geffrey);                                    // integral class

    double integral_test{};                                            // sum of elements in a block 

    vector<double> prog_sum_unif(N, 0.);
    vector<double> prog_sum_unif2(N, 0.);
    vector<double> prog_error_unif(N, 0.);
    Integranda coseno;

    vector<double> prog_sum_simp(N, 0.);
    vector<double> prog_sum_simp2(N, 0.);
    vector<double> prog_error_simp(N, 0.);
    Integranda_simp sviluppo;

    vector<double> prog_sum_badsimp(N, 0.);
    vector<double> prog_sum_badsimp2(N, 0.);
    vector<double> prog_error_badsimp(N, 0.);
    Integranda_IS_bad cosx2;

    ofstream out;
    out.open("Data1.csv");
    ofstream costiout;
    costiout.open("Costi.csv");

        
    double delta_t_unif{}, delta_t_simp{};
    double unif_cost{}, simp_cost{};
    clock_t t_start_unif, t_final_unif, t_start_simp, t_final_simp;


    //uniform distribution sampling
    t_start_unif = clock();
    for(int i=0; i<N; i++){                                                       //cycle on N blocks                         
        integral_test = calculator.IntegralAVE(0., 1., coseno, L);                           //block element = integral evaluation with L points
        
        if(i==0){           //mean and square mean of block 0
            prog_sum_unif[i] = integral_test;
            prog_sum_unif2[i] = pow(integral_test,2);
        }
        else{               //add new mean and mean squared to form next block
            prog_sum_unif[i] = prog_sum_unif[i-1] + integral_test;
            prog_sum_unif2[i] = prog_sum_unif2[i-1] + pow(integral_test,2);       
        }
        integral_test = 0;                                                        //reset test to 0 for new block element
    }

    for(int i=0; i<N; i++){
        prog_sum_unif[i] /= i+1;                                             //progressive mean  (divide per number of block added)
        prog_sum_unif2[i] /= i+1;                                            //progressive square mean (divide per number of block added)
        prog_error_unif[i] = error(prog_sum_unif, prog_sum_unif2, i);        //error per new added block -- it's the std
    }

    t_final_unif = clock();
    delta_t_unif = double(t_final_unif - t_start_unif)/(CLOCKS_PER_SEC);
    cout << "tempo metodo media:\t" << delta_t_unif << "\n";

    //importance sampling
    t_start_simp = clock();
    for(int i=0; i<N; i++){                                                                 //cycle on N blocks                         
        integral_test = calculator.IntegralSIMP(0., 1., 1.5, sviluppo, L);              //block element = integral evaluation with L points
        
        if(i==0){           //mean and square mean of block 0
            prog_sum_simp[i] = integral_test;
            prog_sum_simp2[i] = pow(integral_test,2);
        }
        else{               //add new mean and mean squared to form next block
            prog_sum_simp[i] = prog_sum_simp[i-1] + integral_test;
            prog_sum_simp2[i] = prog_sum_simp2[i-1] + pow(integral_test,2);       
        }
        integral_test = 0;                                                        //reset test to 0 for new block element
    }

    for(int i=0; i<N; i++){
        prog_sum_simp[i] /= i+1;                                             //progressive mean  (divide per number of block added)
        prog_sum_simp2[i] /= i+1;                                            //progressive square mean (divide per number of block added)
        prog_error_simp[i] = error(prog_sum_simp, prog_sum_simp2, i);        //error per new added block -- it's the std
    }
    t_final_simp = clock();
    delta_t_simp = double(t_final_simp - t_start_simp)/(CLOCKS_PER_SEC);
    cout << "tempo importance sampling:\t" << delta_t_simp << "\n";
    unif_cost = delta_t_unif*prog_error_unif[prog_error_unif.size()-1];
    simp_cost = delta_t_simp*prog_error_simp[prog_error_simp.size()-1];
    cout << scientific << "\nmean method cost:\t" << unif_cost << "\n";
    cout << "importance sampling cost:\t" << simp_cost << "\n";
    costiout << scientific << unif_cost << "," << simp_cost << "\n";


    //importance sampling with different (and not so good) p(x)
    for(int i=0; i<N; i++){                                                                 //cycle on N blocks                         
        integral_test = calculator.IntegraleIS_bad(0., 1., cosx2, L);                                  //block element = integral evaluation with L points
        
        if(i==0){           //mean and square mean of block 0
            prog_sum_badsimp[i] = integral_test;
            prog_sum_badsimp2[i] = pow(integral_test,2);
        }
        else{               //add new mean and mean squared to form next block
            prog_sum_badsimp[i] = prog_sum_badsimp[i-1] + integral_test;
            prog_sum_badsimp2[i] = prog_sum_badsimp2[i-1] + pow(integral_test,2);       
        }
        integral_test = 0;                                                        //reset test to 0 for new block element
    }

    for(int i=0; i<N; i++){
        prog_sum_badsimp[i] /= i+1;                                             //progressive mean  (divide per number of block added)
        prog_sum_badsimp2[i] /= i+1;                                            //progressive square mean (divide per number of block added)
        prog_error_badsimp[i] = error(prog_sum_badsimp, prog_sum_badsimp2, i);  //error per new added block -- it's the std
    }

    for(int i=0; i<N; i++){
        out << prog_sum_unif[i] << "," << prog_error_unif[i] << "," << prog_sum_badsimp[i] << "," << prog_error_badsimp[i] << "," << prog_sum_simp[i] << "," << prog_error_simp[i] << "\n";
    }
    cout << "Data created. All praise geffrey." << endl;

    out.close();
    costiout.close();

    return 0;
}

    