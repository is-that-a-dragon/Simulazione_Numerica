#include "libre.h"

using namespace std;

int main(){

    Random geffrey;                                                     //random generator
    wake_up_geffrey(geffrey, "sing!");                                           //initialize

    double L_stick{0.8};                                                //stick lenght
    double d{1.};                                                       //line distance
    M *= 10;                                                            // 10^7 total throws
    //N *= 10;                                                          // 1000 number of blocks
    int Nhit{};                                                         //number of intersections

    vector<double> pie(N, 0.);
    vector<double> pie2(N, 0.);
    vector<double> prog_error(N, 0.);

    for(int i=0; i<N; i++){                                             //cycle on N blocks                         
        for(int j=0; j<L; j++)                                          //cycle on the L elements of i-th block
            if(stick_intersection(L_stick, geffrey))
                Nhit ++;                                                //number of intersection per block
        if(i==0){
            pie[i] = 2*L_stick*L/(Nhit*d);                              //pi of block 0
            pie2[i] = pow(2*L_stick*L/(Nhit*d),2);                      //square pi of block 0
        }
        else{
            pie[i] = pie[i-1] + 2*L_stick*L/(Nhit*d);                   //add new pi value and pi squared value
            pie2[i] = pie2[i-1] + pow(2*L_stick*L/(Nhit*d),2);          //to form next block
        }
        Nhit = 0;                                                       //reset Nhit to 0 for new block
    }

    for(int i=0; i<N; i++){ 
        pie[i] /= i+1;                                                  //progressive value of pi  (divide per number of block added)
        pie2[i] /= i+1;                                                 //progressive square value of pi (divide per number of block added)
        prog_error[i] = error(pie, pie2, i);                            //error per new added block
    }

    Print_2(pie, prog_error, "Data3");

    return 0;
}