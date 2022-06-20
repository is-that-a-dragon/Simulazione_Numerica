/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include <limits>
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

int main()
{ 
  Input(); //Inizialization
  for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
  {
    Reset(iblk);   //Reset block averages
    for(int istep=1; istep <= nstep; ++istep)
    {
      Move(metro);
      Measure();
      Accumulate(); //Update block averages
    }
    Averages(iblk, nblk);   //Print results for current block, if iblk == nblk print to final
  }

  ConfFinal(); //Write final configuration

  return 0;
}




void Input(void)
{
  ifstream ReadInput;

  cout << "Classic 1D Ising model             " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Nearest neighbour interaction      " << endl << endl;
  cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
  cout << "The program uses k_B=1 and mu_B=1 units " << endl;

  //Read seed for random numbers
  int p1, p2;
  ifstream Primes("Primes");
  for(int i=0; i < 16; ++i){
        Primes.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
    }
  Primes >> p1 >> p2;
  Primes.close();

  string start;
  cout << "\033[1m\033[37m" << "\n\nFirst time, man?\t(y/n)\n" << "\033[0m";
  cin >> start;

  if(start == "y"){           //se non è la prima volta uso i semi di uscita per riproducibilità
    ifstream inseed("seed.in");
    inseed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
    rnd.SetRandom(seed,p1,p2);
    inseed.close();
    equilibration = true;
  }
  else{
    ifstream inseed("seed.out");
    if(!inseed){
      cout << "Error opening file 'seed.out'\t I'm a failure, end the program\n";
    }
    inseed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
    rnd.SetRandom(seed,p1,p2);
    inseed.close();
    equilibration = false;
  }
  
 //Read input informations
  ReadInput.open("input.dat");

  ReadInput >> temp;
  beta = 1.0/temp;
  cout << "Temperature = " << temp << endl;

  ReadInput >> nspin;
  cout << "Number of spins = " << nspin << endl;

  ReadInput >> J;
  cout << "Exchange interaction = " << J << endl;

  ReadInput >> h;
  cout << "External field = " << h << endl << endl;
    
  ReadInput >> metro; // if=1 Metropolis else Gibbs

  ReadInput >> nblk;

  ReadInput >> nstep;

  if(metro==1) cout << "The program perform Metropolis moves" << endl;
  else cout << "The program perform Gibbs moves" << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();


 //Prepare arrays for measurements
  iu = 0; //Energy
  ic = 1; //Heat capacity
  im = 2; //Magnetization
  ix = 3; //Magnetic susceptibility
 
  n_props = 4; //Number of observables

 //initial configuration
  if(start == "y"){           //set random spin if it's the first time 
    for (int i=0; i<nspin; ++i){
      if(rnd.Rannyu() >= 0.5) s[i] = 1;
      else s[i] = -1;
    }
  }
  else{                       //read file config.out if it's not the first time
    ifstream Readconf("config.final");
    for(int i=0; i<nspin; i++)
      Readconf >> s[i];
  }  
 //Evaluate energy etc. of the initial configuration
  Measure();

 //Print initial values for the potential energy and virial
  cout << "Initial energy = " << walker[iu]/(double)nspin << endl;
}


void Move(int metro)
{
  int o;
  double p, energy_old, energy_new, sm;
  double deltaE;
  double energy_up, energy_down;

  for(int i=0; i<nspin; ++i)
  {
  //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
    o = (int)(rnd.Rannyu()*nspin);

    if(metro==1){ //Metropolis
      energy_old = Boltzmann(s[o], o);
      s[o] *= -1;
      energy_new = Boltzmann(s[o], o);
      deltaE = energy_new - energy_old;

      accepted ++; //se ho rifiutato la mossa accepted rimane allo stesso valore di prima
      attempted ++;

      if( deltaE > 0){
        if( rnd.Rannyu() > exp( -beta*deltaE ) ){
          s[o] *= -1; //ovvero non varia nulla, rifiuta la mossa
          accepted --;
        }
      }
    }

    else{ //Gibbs sampling   
        energy_up = Boltzmann(1, o);
        energy_down = Boltzmann(-1,o);
        deltaE = energy_down - energy_up;
        p = 1./(1. + exp(-beta*deltaE));

        if( rnd.Rannyu() < p)
          s[o] = 1;
        else
          s[o] = -1;
        
        accepted ++;
        attempted ++;
    }
  }
}

double Boltzmann(int sm, int ip) //valore spin e posizione dello spin
{
  double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm; //energia sui primi vicini
  return ene;
}

void Measure()
{
  int bin;
  double u = 0.0, m = 0.0;

 //cycle over spins
  for (int i=0; i<nspin; ++i) {
    u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
    m += s[i];
  }
  walker[iu] = u;
  walker[ic] = u*u; // qui mi basta mettere u^2, poi nelle medie lo userò di conseguenza
  walker[im] = m;
  walker[ix] = m*m;
}


void Reset(int iblk) //Reset block averages
{
   
   if(iblk == 1)
   {
       for(int i=0; i<n_props; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}


void Accumulate(void) //Update block averages
{

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] += walker[i];
   }
   blk_norm = blk_norm + 1.0;
}


void Averages(int iblk, int nblk) //Print results for current block
{
   ofstream Ene, Heat, Mag, Chi, Equil;
   const int wd=12;
    
    cout << "Block number " << iblk << endl;
    cout << "Acceptance rate " << accepted/attempted << endl << endl;

    //compute averages
    stima_u = blk_av[iu]/blk_norm/(double)nspin; //Energy
    glob_av[iu]  += stima_u;
    glob_av2[iu] += stima_u*stima_u;
    err_u = Error(glob_av[iu],glob_av2[iu],iblk);

    stima_c = ( beta*beta*( blk_av[ic]/blk_norm - pow(blk_av[iu]/blk_norm,2) ) ) /(double)nspin; //Heat capacity
    glob_av[ic]  += stima_c;
    glob_av2[ic] += stima_c*stima_c;
    err_c = Error(glob_av[ic],glob_av2[ic],iblk);

    stima_m = blk_av[im]/blk_norm/(double)nspin; //Magnetisation
    glob_av[im] += stima_m;
    glob_av2[im] += stima_m*stima_m;
    err_m = Error(glob_av[im], glob_av2[im], iblk);

    stima_x = ( beta*blk_av[ix] ) /blk_norm/(double)nspin; //Magnetic susceptibility
    glob_av[ix]  += stima_x;
    glob_av2[ix] += stima_x*stima_x;
    err_x = Error(glob_av[ix],glob_av2[ix],iblk);


    cout << "----------------------------" << endl << endl;
    
    //print file for equilibration
    if(equilibration){
      Equil.open("output.equil.dat",ios::app);
      // ith-step == i-th block --- energy 
      Equil << iblk << setw(wd) << stima_u << endl;
    }

    else {
      Ene.open("output.ene.0",ios::app);
      Heat.open("output.heat.0",ios::app);
      Mag.open("output.mag.0",ios::app);
      Chi.open("output.suscet.0",ios::app);

      // block number --- mean on ith block --- progressive mean --- progressive error
      Ene << setw(wd) << iblk <<  setw(wd) << stima_u << setw(wd) << glob_av[iu]/(double)iblk << setw(wd) << err_u << endl;
      Heat << setw(wd) << iblk <<  setw(wd) << stima_c << setw(wd) << glob_av[ic]/(double)iblk << setw(wd) << err_c << endl;
      Mag << setw(wd) << iblk <<  setw(wd) << stima_m << setw(wd) << glob_av[im]/(double)iblk << setw(wd) << err_m << endl;
      Chi << setw(wd) << iblk <<  setw(wd) << stima_x << setw(wd) << glob_av[ix]/(double)iblk << setw(wd) << err_x << endl;

      //print final values
      if(iblk == nblk){
        ofstream Final;
        Final.open("gibbs_wh.dat",ios::app);
        // temperature --- energy --- energy_error --- heat capacity --- hc_error --- magnetisation --- mag_error --- mag_susceptibility --- mag_susceptibility_error
        Final << temp << "," << glob_av[iu]/(double)iblk << "," << err_u  << "," << glob_av[ic]/(double)iblk << "," << err_c  << "," << glob_av[im]/(double)iblk << "," << err_m << "," << glob_av[ix]/(double)iblk << "," << err_x << endl;
        Final.close();
        cout << "\nFinal values printend in gibbs_wh.dat\n\n";
      }
    }
    Ene.close();
    Heat.close();
    Mag.close();
    Chi.close();
    Equil.close();
}


void ConfFinal(void)
{
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");
  for (int i=0; i<nspin; ++i)
  {
    WriteConf << s[i] << endl;
  }
  WriteConf.close();

  rnd.SaveSeed();
  if(equilibration)
    cout << "\033[1m\033[37m" << "\nEquilibration finished, have a nice day and restart the program!\n\n" << "\033[0m";
}

int Pbc(int i)  //Algorithm for periodic boundary conditions
{
    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
}

double Error(double sum, double sum2, int iblk)
{
    if(iblk==1) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
