#pragma once

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>

#include "random.h"

using namespace std;

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

template <typename T> int signum(const T val){
    return ((T(0)<val) - (val<T(0))); 
}

bool stick_intersection(double &l, Random& geffrey){
    //importat: x coordinates are not useful since the problem is x invariant
    //if a stick doesn't cross an horizontal line, translating it horizontally won't do the trick either
    double y0 = geffrey.Rannyu();                   //generate stick starting point in [0;1)
    double v = geffrey.Rannyu();
    double u = geffrey.Rannyu(-1, 1);               //generate 2 points in rectangle 2x1

    while(pow(u,2) + pow(v,2) > 1){                 //accept - reject: are the points in the radius?               
        v = geffrey.Rannyu();
        u = geffrey.Rannyu(-1, 1);
    }

    double sin_theta{(2*u*v) / (pow(u,2)+pow(v,2))};

    double y = y0 + l*sin_theta ;                //find y coordinate of stick end point

    if(y0 == 0 || y > 1 || y < 0)               // since the problem is y invariant i only need 2 lines to simulate
        return true;                            // an entire plane
    else
        return false;
    //see this two articles about how to extract sin from two random numbers
    // 1. https://math.stackexchange.com/questions/3183253/sine-cosine-of-random-angle-from-0-to-2-pi
    // 2. https://pdg.lbl.gov/2012/reviews/rpp2012-rev-monte-carlo-techniques.pdf (section 37.4.3)
}
