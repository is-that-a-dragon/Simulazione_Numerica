#pragma once

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>

#include "random.h"

using namespace std;
#define _USE_MATH_DEFINES

int M = pow(10,8);               //throws
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

    if(sing=="sing"){
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

vector<int> make_trajectory(Random &geffrey, int L){
    vector<int> traj (3,0);            //vector with x,y,z

    if(L==0)
        return traj;

    for(int i = 0 ; i < L ; i++){                               //L: RW lenght
        if(geffrey.Rannyu() < 0.5)                              //right or left step
            traj[floor(geffrey.Rannyu(0., 3.))] += 1;           //choose which coordinate to increase/decrease
        else
            traj[floor(geffrey.Rannyu(0., 3.))] -= 1;
    }

return traj;
}

vector<double> make_trajectory_continous( Random &geffrey , int L ){ //L is number of steps to take
    vector<double> traj (3,0);       //vector with x,y,z
    if(L==0)
        return traj;

    double phi{}, cos_theta{}, sin_theta{}; //phi angolo polare, theta azimutale

    for( int i{}; i < L; i++){
        cos_theta = geffrey.Rannyu(-1.,1.);
        sin_theta = sqrt(1-pow(cos_theta,2));      // sin must be in [0,1] bc theta is in [0, pi] but that's fine since it's a sqrt
        phi = geffrey.Rannyu(0, 2*M_PI);

        traj[0] += cos(phi)*sin_theta;             //make step
        traj[1] += sin(phi)*sin_theta;             //convert from spherical to cartesian, R = 1
        traj[2] += cos_theta;
    }

    return traj;
}

template <typename T> double dist_from_origin_square(vector<T> A){
    double dist2{};
    for(unsigned int i{}; i < A.size(); i++) //A.size is 3, one for each component in the 3D space
        dist2 += pow(A[i],2.); //distance squared from the origin = (0-X_final)^2 + (0-Y_final)^2 + (0-Z_final)^2
    // remember that you have to calculate sqrt( pow( distance,2 ) ), the sqrt will be taken later
    return dist2;
}