#include "libre.h"

using namespace std;

int main (){

    double sum{};                                                       //sum of elements in a block                                     
    vector<double> prog_sum(N, 0.);
    vector<double> prog_sum2(N, 0.);
    vector<double> prog_error(N, 0.);

    int Mchi = 100;                                                     //chi squared, ex 1.3
    int nchi = pow(10,4);
    int Lchi = int(nchi/Mchi);
    int rep = 100;
    vector<int> conteggi(Mchi, 0);
    vector<double> chi_squared(rep, 0.);

    Random geffrey;                                                     //random generator
    wake_up_geffrey(geffrey);                                           //initialize

    for(int i=0; i<N; i++){                                             //cycle on N blocks                         
        for(int j=0; j<L; j++)                                          //cycle on the L elements of i-th block
            sum += geffrey.Rannyu();                                    //sum of L numbers
        if(i==0){
            prog_sum[i] = sum/double(L);                                //mean of block 0
            prog_sum2[i] = pow(sum/double(L),2);                        //square mean of block 0
        }
        else{
            prog_sum[i] = prog_sum[i-1] + sum/double(L);                //add new mean and mean squared to form next block
            prog_sum2[i] = prog_sum2[i-1] + pow(sum/double(L),2);       
        }
        sum = 0;                                                        //reset sum to 0 for new block
    } 
    for(int i=0; i<N; i++){
        prog_sum[i] /= i+1;                                             //progressive mean (divide per number of block added)
        prog_sum2[i] /= i+1;                                            //progressive square mean (divide per number of block added)
        prog_error[i] = error(prog_sum, prog_sum2, i);                  //error per new added block
    }

    Print_2(prog_sum, prog_error, "Data1_1");

    //ex 1.2
    for(int i=0; i<N; i++){                                             //cycle on N blocks                         
        for(int j=0; j<L; j++)                                          //cycle on the L elements of i-th block
            sum += pow(geffrey.Rannyu()-0.5, 2);                        //sum of (random number-0.5)^2 
        if(i==0){
            prog_sum[i] = sum/double(L);                                //mean of block 0
            prog_sum2[i] = pow(sum/double(L),2);                        //square mean of block 0
        }
        else{
            prog_sum[i] = prog_sum[i-1] + sum/double(L);                //add new mean and mean squared to form next block
            prog_sum2[i] = prog_sum2[i-1] + pow(sum/double(L),2);       
        }
        sum = 0;                                                        //reset sum to 0 for new block
    }
    for(int i=0; i<N; i++){          
        prog_sum[i] /= i+1;                                             //progressive mean  (divide per number of block added)
        prog_sum2[i] /= i+1;                                            //progressive square mean (divide per number of block added)
        prog_error[i] = error(prog_sum, prog_sum2, i);                  //error per new added block
    }

    Print_2(prog_sum, prog_error, "Data1_2");

    //ex 1.3
    for(int j=0; j<rep; j++){                                           // j-th chi squared
        for(int i=0; i<nchi; i++)
            conteggi[floor(geffrey.Rannyu() * 100)] ++;                 // generate nchi random number, increment relative position in vector
        for(int k=0; k<Mchi; k++)
            chi_squared[j] += pow((conteggi[k]-Lchi),2)/Lchi;           // j-th chi squared using relative conteggi vector
        fill(conteggi.begin(), conteggi.end(), 0);                      // empty conteggi
    }

    Print(chi_squared, "Data1_3");

    return 0;
}