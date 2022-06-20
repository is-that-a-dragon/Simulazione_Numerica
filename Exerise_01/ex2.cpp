#include "libre.h"

using namespace std;

int main(){

    vector<int> N {1,2,10,100};                   //how many elements in Sn per different cycle

    int Sn_tries = pow(10,4);
    double sunif{};
    double sexp{};
    double slor{};

    Random geffrey;
    wake_up_geffrey(geffrey);

    ofstream out;
    out.open("Data2");
    if(!out)
        cerr << "Error: can't open file!\n";

    for(unsigned int i=0; i<N.size(); i++){                                   //cicle on the 4 possibilities of N
        for(int j=0; j<Sn_tries; j++){                                        //number of Sn
            for(int k=0; k<N[i] ; k++){                                       //add N[i] elements in Sn
                sunif += geffrey.Rannyu();
                sexp += geffrey.Exp(1.);
                slor += geffrey.Lorentz(0.,1.);
            }
            out << sunif/double(N[i]) << "," << sexp/double(N[i]) << "," << slor/double(N[i]) << "\n"; //actual Sn
            sunif = 0.;                                                       //reset values
            sexp = 0.;
            slor = 0. ;
        }
    }

    cout << "Data2 created successfully, all praise geffrey.\n";
    out.close();
    return 0;
}