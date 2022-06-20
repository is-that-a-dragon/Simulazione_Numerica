#include "random.h"
#include "libre.h"

using namespace std;

int main(){

    Random geffrey;
    wake_up_geffrey(geffrey);

    M = pow(10,5);                 //number of simulations
    L = 100;                       //RW  lenght and number of blocks

    ofstream out("Data2.csv");

    if(!out){
        cout << "Error in file opening!\n";
        return -1;
    }

    vector<int> trajectory;
    vector<double> trajectory_cont;

    vector<double> sum_prog (L,0);
    vector<double> sum_prog_2 (L,0);
    vector<double> errors (L,0);
    vector<double> sum_continous (L,0);
    vector<double> sum_continous_2 (L,0);
    vector<double> errors_continous (L,0);

    sum_prog[0]=0;                //first step is always at 0
    sum_prog_2[0]=0;
    errors[0]=0;
    sum_continous[0]=0;
    sum_continous_2[0]=0;
    errors_continous[0]=0 ;

    out << 0 << "," << 0 << "," << 0 << "," << 0 << endl;

    double lattice_dist{};
    double continous_dist{};

    for(int i = 1 ; i <= L; i++){                                   // each iteration increase the number of steps in the trajectory
        for(int j{}; j < L ; j++){                                  // cycle on blocks
            for(int k{}; k<int(M/L) ; k++){                         // make M/L trajectories (aka random walks), each made of i steps
            trajectory = make_trajectory(geffrey, i);                       
            trajectory_cont = make_trajectory_continous(geffrey, i);

            lattice_dist += dist_from_origin_square(trajectory);    // sum the squared distance from origin of each traj, will be averaged later
            continous_dist += dist_from_origin_square(trajectory_cont);
        }
        // each block contains the average squared distances from the origin
        // one block element is the average of M/L squared distances
        // so for each step I average the distance on L * M/L = M random walks
        // repeat for L = 100 steps
        sum_prog[j] = sum_prog[j-1] + lattice_dist/(M/L);
        sum_prog_2[j] = sum_prog_2[j-1] + pow(lattice_dist/(M/L),2);  
        sum_continous[j] = sum_continous[j-1] + continous_dist/(M/L);
        sum_continous_2[j] = sum_continous_2[j-1] + pow(continous_dist/(M/L),2);
        lattice_dist = 0.;                      //reset step variable to zero
        continous_dist = 0.;
        }

        for(int s{}; s < L; s++){                //mean average of squared distance (progressive)
            sum_prog[s] /= (s+1);
            sum_prog_2[s] /= (s+1);
            errors[s] = error(sum_prog,sum_prog_2,s) / ( 2*sqrt(sum_prog[s]) ); //error propagation
            sum_continous[s] /= (s+1);
            sum_continous_2[s] /= (s+1);
            errors_continous[s] = error(sum_continous,sum_continous_2,s) / ( 2*sqrt(sum_continous[s]) ); //error propagation
        }

        //error propagation: sum_prog has d^2 in it, I want the error on d = sqrt(d^2)
        //error = (error on d^2) * (derivative on sqrt(sum_prog) ) = (error on d^2) / (2*sqrt(sum_prog))

        //print the last value, which is the most precise bc it averages on M random walks
        //i goes from 0 to L and increase each cycle
        //so going down the column of printed values the number of step per random walk increases from 0 to 100
        out << sqrt(sum_prog[L-1]) << "," << errors[L-1] << "," << sqrt(sum_continous[L-1]) << "," << errors_continous[L-1] << endl;

        fill(sum_prog.begin(), sum_prog.end(), 0.);           //reset vectors to zero
        fill(sum_prog_2.begin(), sum_prog_2.end(), 0.);
        fill(errors.begin(), errors.end(), 0.);
        fill(sum_continous.begin(), sum_continous.end(), 0.);
        fill(sum_continous_2.begin(), sum_continous_2.end(), 0.);
        fill(errors_continous.begin(), errors_continous.end(), 0.);
    }

    out.close();
    cout << "Data2.csv created. All praise geffrey.\n"; 

    return 0;
}