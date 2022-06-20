#include "libre.h"

using namespace std;

int main (){
    double T{1}, K{100}, r{0.1}, sigma{0.25};
    double S0{100};
    double S_T{};

    Random geffrey;
    wake_up_geffrey(geffrey);
    double call_sum{};
    double put_sum{};                                    
    vector<double> prog_call_sum(N, 0.);
    vector<double> prog_call_sum2(N, 0.);
    vector<double> prog_call_error(N, 0.);
    vector<double> prog_put_sum(N, 0.);
    vector<double> prog_put_sum2(N, 0.);
    vector<double> prog_put_error(N, 0.);

    for(int i=0; i<N; i++){                                             //cycle on N blocks                         
        for(int j=0; j<L; j++){                                         //cycle on L elements of i-th block
            S_T = generate_S_T(S0, r, sigma, T, geffrey);               //generate price at time T in one step
            call_sum += claim_call(r, T, S_T, K);                       //calculate call price as of now
            put_sum += claim_put(r, T, S_T, K);                         //calculate put price as of now
        }       //it's ok to use the same S(T) for call and put since both take into accout just the price at time T
        if(i==0){
            prog_call_sum[i] = call_sum/double(L);                                //mean of block 0
            prog_call_sum2[i] = pow(call_sum/double(L),2);                        //square mean of block 0
            prog_put_sum[i] = put_sum/double(L);
            prog_put_sum2[i] = pow(put_sum/double(L),2);
        }
        else{
            prog_call_sum[i] = prog_call_sum[i-1] + call_sum/double(L);
            prog_call_sum2[i] = prog_call_sum2[i-1] + pow(call_sum/double(L),2);
            prog_put_sum[i] = prog_put_sum[i-1] + put_sum/double(L);
            prog_put_sum2[i] = prog_put_sum2[i-1] + pow(put_sum/double(L),2);       
        }
        call_sum = 0.;
        put_sum = 0.;                                                        //reset sum to 0 for new block
    }
    for(int i=0; i<N; i++){
        prog_call_sum[i] /= i+1;                                             //progressive mean  (divide per number of block added)
        prog_call_sum2[i] /= i+1;                                            //progressive square mean (divide per number of block added)
        prog_call_error[i] = error(prog_call_sum, prog_call_sum2, i);        //error per new added block
        prog_put_sum[i] /= i+1;
        prog_put_sum2[i] /= i+1;
        prog_put_error[i] = error(prog_put_sum, prog_put_sum2, i);
    }

    Print_2(prog_call_sum, prog_call_error, "Direct_call.csv");
    Print_2(prog_put_sum, prog_put_error, "Direct_put.csv");



    //part 2: now it's progressive baby
    fill(prog_call_sum.begin(), prog_call_sum.end(), 0.);               //reset all vectors to zero
    fill(prog_call_sum2.begin(), prog_call_sum2.end(), 0.);
    fill(prog_call_error.begin(), prog_call_error.end(), 0.);
    fill(prog_put_sum.begin(), prog_put_sum.end(), 0.);
    fill(prog_put_sum2.begin(), prog_put_sum2.end(), 0.);
    fill(prog_put_error.begin(), prog_put_error.end(), 0.);

    for(int i=0; i<N; i++){                                              //cycle on N blocks                         
        for(int j=0; j<L; j++){                                          //cycle on the L elements of i-th block
            S_T = generate_S_T(S0, r, sigma, T, geffrey, 100);           //generate price at time T in 100 steps
            call_sum += claim_call(r, T, S_T, K);
            put_sum += claim_put(r, T, S_T, K);
        }
        if(i==0){
            prog_call_sum[i] = call_sum/double(L);                                //mean of block 0
            prog_call_sum2[i] = pow(call_sum/double(L),2);                        //square mean of block 0
            prog_put_sum[i] = put_sum/double(L);
            prog_put_sum2[i] = pow(put_sum/double(L),2);
        }
        else{
            prog_call_sum[i] = prog_call_sum[i-1] + call_sum/double(L);
            prog_call_sum2[i] = prog_call_sum2[i-1] + pow(call_sum/double(L),2);
            prog_put_sum[i] = prog_put_sum[i-1] + put_sum/double(L);
            prog_put_sum2[i] = prog_put_sum2[i-1] + pow(put_sum/double(L),2);       
        }
        call_sum = 0.;
        put_sum = 0.;                                                        //reset sum to 0 for new block
    }
    for(int i=0; i<N; i++){
        prog_call_sum[i] /= i+1;                                             //progressive mean  (divide per number of block added)
        prog_call_sum2[i] /= i+1;                                            //progressive square mean (divide per number of block added)
        prog_call_error[i] = error(prog_call_sum, prog_call_sum2, i);                  //error per new added block
        prog_put_sum[i] /= i+1;
        prog_put_sum2[i] /= i+1;
        prog_put_error[i] = error(prog_put_sum, prog_put_sum2, i);
    }

    Print_2(prog_call_sum, prog_call_error, "Progressive_call.csv");
    Print_2(prog_put_sum, prog_put_error, "Progressive_put.csv");

    return 0;
}