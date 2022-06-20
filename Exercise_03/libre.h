#pragma once

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>

#include "random.h"

using namespace std;
#define _USE_MATH_DEFINES

int M = pow(10,6);               //throws
int N = 100;                     //block number
int L = int(M/N);                //elements per block



double error(const vector <double> ave, const vector <double> ave2, const int n){
    if(n==0)
        return 0;
    else
        return sqrt((ave2[n] - pow(ave[n],2))/n);
}


void wake_up_geffrey(Random &geffrey, string sing=""){
    ifstream inprimes("Primes");
    ifstream inseed("seed.in");
    int p1{}, p2{};
    int seed[4];
    string check;

    if(!inprimes)
        cerr << "Primes file not opening!\n";
    else{
        inprimes >> p1 >> p2;
    }

    if(!inseed)
        cerr << "Seed file not opening!\n";
    else{
        while(!inseed.eof()){
            inseed >> check;
            if(check == "RANDOMSEED"){
                inseed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            }
        }
        geffrey.SetRandom(seed, p1, p2);
    }

    inprimes.close();
    inseed.close();
    cout << "geffrey has awoken\n";

    if(sing=="sing!"){
        cout << "\033[1m\033[37m" <<"♫     ♪\n";
        cout << "\033[1m\033[34m" << "CEO, entrepreneur\nBorn in 1964\nJeffrey\nJeffrey Bezos\n";
        cout << "\033[1m\033[36m" << "CEO, entrepreneur\nBorn in 1964\nJeffrey\nJeffrey Bezos\n";
        cout << "\033[1m\033[32m" << "Come on, Jeffrey, you can do it\nPave the way, put your back into it\n";
        cout << "\033[1m\033[33m" << "Tell us why\nShow us how\nLook at where you came from\nLook at you now\n";
        cout << "\033[1m\033[31m" << "Zuckerberg and Gates and Buffet\nAmateurs can fucking suck it\n";
        cout << "\033[1m\033[35m" << "Fuck their wives, drink their blood\nCome on, Jeff, get 'em!\n";
        cout << "\033[1m\033[37m" <<"♫     ♪\n" << "\033[0m";
    }
}

template <typename T> void Print(const vector<T> &v, const char *filename) {
    ofstream outf(filename); 
    for (T elem : v){
        outf << elem << "\n";}
    outf.close();
    cout << filename << " created successfully, all praise geffrey.\n";
}

template <typename T> void Print_2(const vector<T> &v1, const vector<T> &v2, const char *filename) {
    ofstream outf(filename);
    if(v1.size() != v2.size())
        cerr << "sizes of given vectors are not equal";
    for (int i=0; i<int(v1.size()); i++){
        outf << v1[i] << "," << v2[i] << "\n" ;    
    }
    outf.close();
    cout << filename << " created successfully, all praise geffrey.\n";
}

double generate_S_T(double S_0, double r, double sigma, double T_final, Random& geffrey, double N=1){
    vector<double> S(N,0.);
    double Z{};
    for(int i{0}; i<N; i++){        //if no N is specified, N is put to 1 and i can only be 0
        Z = geffrey.Gauss(0.,1.);
        if(i==0)                    //if N == 1, T_final is not divided and S is generateed in one step
            S[i] = (S_0*exp( ( r-pow(sigma,2)/2 )*(T_final/N) + sigma*Z*sqrt(T_final/N) ));
        else
            S[i] = (S[i-1]*exp( ( r-pow(sigma,2)/2 )*(T_final/N) + sigma*Z*sqrt(T_final/N) ));  
    }
    return S[S.size()-1];           //return last element of vector which is either S[0] if N==1 or S[99]
}

double claim_call(double r, double T_final, double S, double K){
        return exp(-r*T_final)*max(0., S-K);
}

double claim_put(double r, double T_final, double S, double K){
    return exp(-r*T_final)*max(0., K-S);
}