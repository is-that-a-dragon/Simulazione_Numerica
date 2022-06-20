#include "integral.h"
#include "random.h"
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <vector>

using namespace std;

const int n_punti=10000, n_eq=10;                   //n points in integral, equilibration points
const double x_0=0., dx=1.;                         //starting point, metropolis step,

const double sigma0=1., mu0=1.;                     //sigma and mu starting values
double dsigma=0.3, dmu=0.3;                         //metropolis step for sigma and mu
const double nstep=200;                             //metropolis test per beta value
const double beta=1., dbeta=0.5, bmax=200+dbeta;    //beta start, metro step beta, beta max  

int M{10000};            // Total number of throws
int N{50};               // Number of blocks
int L{int(M/N)};         // Number of throws in each block

double psi_squared(vector <double> v){
    double x=v[0];
    double sigma=v[1];
    double mu=v[2];
    double s1=pow((x+mu)/sigma,2);
    double s2=pow((x-mu)/sigma,2);
    double e1=exp(-s1/2);
    double e2=exp(-s2/2);
    
    return pow(e1+e2,2);
}

double V(double x){
    return pow(x,4)-5./2.*pow(x,2);
    //return pow(x,2)/2; //armonic potential test
}

double E(vector <double> v){
    double x=v[0];
    double sigma=v[1];
    double mu=v[2];
    double s1=pow((x+mu)/sigma,2);
    double s2=pow((x-mu)/sigma,2);
    double e1=exp(-s1/2);
    double e2=exp(-s2/2);
    
    return -(s2*e2+s1*e1-e1-e2)/(e1+e2)/(2*sigma*sigma)+V(x);
}

int main(){
    Random geffrey;
    wake_up_geffrey(geffrey);
    
    vector <double> variables(3);
    Integral cicciopizzo(E);
    
    variables[0]=x_0;
    variables[1]=sigma0;
    variables[2]=mu0;
    
    double old_integral = cicciopizzo.Metropolis(n_punti, psi_squared, variables, dx, n_eq);
    double sigma_old    = variables[1];
    double mu_old       = variables[2];
    
    cout<<endl<<"beta =\t\t"<<beta<<"\nmu =\t\t"<<variables[2]<<"\nsigma =\t\t"<<variables[1]<<"\nI =\t\t"<<old_integral<<endl;
    
    ofstream output1;
    output1.open("integral_beta.dat");
    
    for(double b=beta; b<bmax; b+=dbeta){
        double sum{};
        double sum2{};
        for(int i=0; i<nstep; i++){
            //choose new mu and ne sigma
            variables[1] = fabs(geffrey.Rannyu(variables[1]-dsigma, variables[1]+dsigma));
            variables[2] = fabs(geffrey.Rannyu(variables[2]-dmu, variables[2]+dmu));
            //evaluate integral with new sigma and mu
            double r = cicciopizzo.Metropolis(n_punti, psi_squared, variables, dx, n_eq);

            if(r < old_integral){         //is it better?
                old_integral = r;
                sigma_old    = variables[1];
                mu_old       = variables[2];
            }
            else{                       //random change   
                if( geffrey.Rannyu() < exp(b*(old_integral-r)) ){
                    old_integral = r;
                    sigma_old    = variables[1];
                    mu_old       = variables[2];
                }   
                else {
                    variables[1] = sigma_old;
                    variables[2] = mu_old;
                }
            } 
            sum += old_integral;
            sum2 += old_integral*old_integral;
        }
        output1<<b<<","<<variables[1]<<","<<variables[2]<<","<<sum/nstep<<","<<sqrt(sum2/nstep-pow(sum/nstep,2))<<endl;
        cout<<endl<<"beta =\t\t"<<b<<"\nmu =\t\t"<<variables[2]<<"\nsigma =\t\t"<<variables[1]<<"\nI =\t\t"<<old_integral<<endl;

        //dsigma = dsigma - dsigma/20.; //change sigma and mu stepsize
        //dmu = dmu - dmu/20.;
    }

    cout<<"\n\n\033[1m\033[33mFinal values:\033[0m"<<endl;
    cout<<"mu = "<<variables[2]<<endl;
    cout<<"sigma = "<<variables[1]<<endl<<endl;

    //once mu and sigma are found
    ofstream output;
    output.open("block_ave.dat");
    double ave{};
    double av2{};
    
    for(int i=0; i<N; i++){                      //N blocks
        double sum{};  
        for (int j=0; j<L; j++)                  //L integral evaluation per block
            sum += cicciopizzo.Metropolis(n_punti,psi_squared,variables,dx,n_eq);
        output<<i+1;                             // block number
        ave += (double)sum/L;                    // single block average
        output<<","<<(double)ave/(i+1);          // printing block average
        av2 += pow(sum/L,2);                     // (single block average)^2
        if(i==0)
            output<<","<<0;
        else
            output<<","<< sqrt(((double)av2/(i+1) - pow((double)ave/(i+1),2))/(i)) ; //error
        output<<endl;
    }
    output.close();
    cout << "\n\033[1m\033[33mMean value of H for best parameters sigma and mu:\n\033[0m";
    cout << "<H> =\t\t" << (double)ave/N << endl << endl; 
    return 0;
}