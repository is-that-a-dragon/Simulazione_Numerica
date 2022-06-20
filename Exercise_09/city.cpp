#include "city.h"

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include <random>

using namespace std;

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
        cout << "\033[1m\033[37m" <<"♫\t\t♪\n";
        cout << "\033[1m\033[34m" << "CEO, entrepreneur\nBorn in 1964\nJeffrey\nJeffrey Bezos\n";
        cout << "\033[1m\033[36m" << "CEO, entrepreneur\nBorn in 1964\nJeffrey\nJeffrey Bezos\n";
        cout << "\033[1m\033[32m" << "Come on, Jeffrey, you can do it\nPave the way, put your back into it\n";
        cout << "\033[1m\033[33m" << "Tell us why\nShow us how\nLook at where you came from\nLook at you now\n";
        cout << "\033[1m\033[31m" << "Zuckerberg and Gates and Buffet\nAmateurs can fucking suck it\n";
        cout << "\033[1m\033[35m" << "Fuck their wives, drink their blood\nCome on, Jeff, get 'em!\n";
        cout << "\033[1m\033[37m" <<"♫\t\t♪\n" << "\033[0m";
    }
}

//funzione costo come il modulo quadro
double distance(city& A, city& B){
    double dist{};
    assert( A.get_dimension() == B.get_dimension() );
    for(int i{}; i < A.get_dimension(); i++ )
        dist += pow(A.get_coord(i) - B.get_coord(i), 2); 
    return dist;
}

void path::print_path(){
    cout << "printing cities in order:\n";
    for(unsigned int i{}; i < m_cities.size(); i++){
        cout << setw(10) << get_city_index(i) << setw(10) << get_city_coord(i,0) << setw(10) << get_city_coord(i,1);
        cout << endl;
    }
    cout << "\npath fit:\t\t" << m_fit << endl;
}

void Travel_Salesman::check(){
    for(unsigned int i{}; i < m_population.size(); i++){
        if( m_population[i].get_size() != static_cast<int>(N_cities) ){
            cout << "one path in population has more than " << N_cities << " cities\n";
            exit(-1);
        }
        if( m_population[i].get_city_index(0) != m_start[0].get_index() ){
            cout << "one path has a different starting city\n";
            exit(-1);
        }
        for(unsigned int j{}; j < m_start.size(); j++ ){
            for(unsigned int k=j+1; k < m_start.size(); k++ ){
                if( m_population[i].get_city_index(j) == m_population[i].get_city_index(k) ){
                    cout << "one path has the same city two times\n";
                    exit(-1);
                }
            }
        }
    }
    //cout << "\033[1m\033[32m" << "\nchecked\n" << "\033[0m";
}

void Travel_Salesman::print_best(string filename){
    ofstream best;
    best.open(filename, ios::app);
    order();
    for(int j{}; j < m_population[0].get_size(); j++){
        best << m_population[0].get_city_index(j) << ",";
        for(int i{}; i < m_population[0].get_city_dim(); i++)
            best << m_population[0].get_city_coord(j, i) << ",";
        best << endl;
    }
    cout << "path fit:\t\t" << get_fit(0);
    best.close();
}

void Travel_Salesman::set_ave_fit(int index){
    double ave_fit{};
    order();
    for(int i{}; i < index; i++){
        ave_fit += get_fit(i);
    }
    average_fit.push_back( (double) ave_fit / index );
}

void Travel_Salesman::order(){
    vector<double> messy_cost( get_pop(), 0. );         //unsorted costs
    for(int i{}; i<get_pop(); i++)
        messy_cost[i] = m_population[i].get_fit();

    vector<double> sorted_cost = messy_cost;
    sort(sorted_cost.begin() , sorted_cost.end());      //now it's sorted, ascending order

    vector<int> sorted_indexes;                            
    //sorted_indexes will contain the indexes of the elements in messy_cost as if it was sorted
    for(unsigned int j{}; j<messy_cost.size(); j++){
        auto e = find( messy_cost.begin(), messy_cost.end(), sorted_cost.at(j) );
        //e is an iterator to found element, if not found e is equal to messy_cost.end()
        if( e != messy_cost.end() )     //aka if the element is found
            sorted_indexes.push_back( distance( messy_cost.begin(), e ) );
    }

    vector<path> copy_pop = m_population;
    for(int i{}; i<get_pop(); i++){
        m_population[i] = copy_pop[ sorted_indexes[i] ];
    }

    check(); 
}

vector<path> Travel_Salesman::select(vector<path> parents){
    double r1 = m_geffrey.Rannyu();
    double r2 = m_geffrey.Rannyu();

    int index_father = static_cast<int> ( get_pop() * pow(r1,4) );
    int index_mother = static_cast<int> ( get_pop() * pow(r2,4) );
    // r is a number between 0 and <1, r^4 is more likely to be closer to 0 than 1
    // since the population will be sorted from more fit to less fit, the index is more likely to choose
    // a path with a good fit

    if(index_father == get_pop())           //just to be sure
        index_father --;
    if(index_mother == get_pop())
        index_mother --;
    
    parents.push_back( m_population[index_father] );
    parents.push_back( m_population[index_mother] );
    return parents;
}

vector<path> Travel_Salesman::crossover(vector<path> parents){
    unsigned int index_cut = static_cast<int>(m_geffrey.Rannyu() * N_cities);
    while (index_cut == 0 or index_cut == N_cities)
        index_cut = static_cast<int>(m_geffrey.Rannyu() * N_cities);

    vector<city> father_head;
    vector<city> mother_head;
    vector<city> father_tail;
    vector<city> mother_tail;
    for(unsigned int i{}; i<N_cities-index_cut; i++){
        father_head.push_back( parents[0].get_city(i) );
        mother_head.push_back( parents[1].get_city(i) );
    }

    for(unsigned int i{}; i < N_cities; i++){
        auto a = find( father_head.begin(), father_head.end(), parents[1].get_city(i) );
        if(a == father_head.end() ) //aka if the city in mother is not present in father_head
            father_tail.push_back( parents[1].get_city(i) );

        auto b = find( mother_head.begin(), mother_head.end(), parents[0].get_city(i) );
        if(b == mother_head.end() ) //aka if the city in father is not present in mother_head
            mother_tail.push_back( parents[0].get_city(i) );
    }

    //merge head and tail
    for(unsigned int i{}; i < father_tail.size(); i++){
        father_head.push_back( father_tail[i] );
        mother_head.push_back( mother_tail[i] );
    }
    parents[0].set_cities( father_head );
    parents[1].set_cities( mother_head );

    return parents;
}

path Travel_Salesman::mutate(path gene){
    vector<city> mutated = gene.get_cities();
    //mutation 1: swap adjacent cities
    if( m_geffrey.Rannyu() <= mutation_prob ){
        int index = floor( m_geffrey.Rannyu(1, N_cities-2) );
        swap( mutated.at(index), mutated.at(index+1) );
    }
    //mutation 2: swap between random cities
    if( m_geffrey.Rannyu() <= mutation_prob ){
        int index1 = floor( m_geffrey.Rannyu(1, N_cities-1) );
        int index2 = floor( m_geffrey.Rannyu(1, N_cities-1) );

        while(index1 == index2)
            index2 = floor( m_geffrey.Rannyu(1, N_cities-1) );
        
        swap( mutated.at(index1), mutated.at(index2) );
    }
    //mutation 3: shift m cities (starting at n) of k
    if( m_geffrey.Rannyu() <= mutation_prob ){
        int n = floor( m_geffrey.Rannyu(1, N_cities-1) );       //starting point
        //N_cities-1 is excluded by Rannyu so, thanks to floor, n is between 1 and N-2
        int m = floor( m_geffrey.Rannyu(0, N_cities-1-n) );     //end point
        //same as above: m is between 0 and N_cities-2-n (last element)
        int k = floor( m_geffrey.Rannyu(1, N_cities-m-n) );     //shift of k
        //same as above: k is a number between 1 and m* - n
        //  (last - start position, a shift of exactly m*-n +1 return every value to its startin position) 
        //where m* is the position indicated by m but starting from the beginning of the array
        //  (m start from the end since I use the iterator .end)
        //now: m* = N_cities -m -1 and m*-n = N_cities -m -1 -n
        //the use of Rannyu and floor lets me skip th -1 in the formula

        rotate( mutated.begin()+n , mutated.begin()+n+k , mutated.end()-m );
        // beginning of original range, which element should go at the beginning after rotation, end of original rage
        // it's a shift to the left
    }
    //mutation 4: swap group of contiguos cities with other group
    if( m_geffrey.Rannyu() <= mutation_prob ){
        int n = floor( m_geffrey.Rannyu( 1, N_cities/2) );          // n < N/2, group 1 starting index
        int k = floor( m_geffrey.Rannyu( 0, (N_cities/2 - n) ) );   // k+1 elements in a group
        int m = floor( m_geffrey.Rannyu( N_cities/2 , N_cities-k ) );
        // group 2 starting index, last value (max N_cities) is excluded by Rannyu and floor -- at max m < N_cities

        for(int i{}; i <= k; i++)
            swap( mutated.at(n+i), mutated.at(m+i) );
    }
    //mutation 5: invert order
    if( m_geffrey.Rannyu() <= mutation_prob ){
        int n = floor( m_geffrey.Rannyu(1, N_cities-1) );       //starting point
        //N_cities-1 is excluded by Rannyu so, thanks to floor, n is between 1 and N-2
        int m = floor( m_geffrey.Rannyu(0, N_cities-1-n) );     //end point
        //same as above: m is between 0 and N_cities-2-n (last element)

        reverse( mutated.begin()+n , mutated.end()-m );
    }
    
    gene.set_cities(mutated);
    gene.set_fit();

    return gene;
}

void Travel_Salesman::evolve(){
    if(m_elite){
        vector<path> new_gen;
        order();
        
        int keep = static_cast<int>(get_pop()/2.);
        for (int i{};  i< get_pop()-keep; i+=2){
            vector<path> parents;
            vector<path> offspring;

            parents = select(parents);
            //faccio comunque il select su tutta la popolazione cosi da non cadere in minimi

            if(m_geffrey.Rannyu() < crossover_prob)
                offspring = crossover(parents);
            else
                for(unsigned int j{}; j < parents.size(); j++)
                    offspring.push_back( parents[j] );

            offspring[0] = mutate( offspring[0] );
            offspring[1] = mutate( offspring[1] );

            new_gen.push_back( offspring[0] );
            new_gen.push_back( offspring[1] );
        }

        for(int i{keep}; i < get_pop(); i++)
            set_path( new_gen[i-keep], i );
        order();
        set_ave_fit( round( get_pop()/2 ) );
        best_fit.push_back( get_fit(0) );
    }

    else{
        vector<path> new_gen;
        order();
        for (int i{};  i< get_pop(); i+=2){
            vector<path> parents;
            vector<path> offspring;

            parents = select(parents);

            if(m_geffrey.Rannyu() < crossover_prob)
                offspring = crossover(parents);
            else
                for(unsigned int j{}; j < parents.size(); j++)
                    offspring.push_back( parents[j] );
        
            offspring[0] = mutate( offspring[0] );
            offspring[1] = mutate( offspring[1] );

            new_gen.push_back( offspring[0] );
            new_gen.push_back( offspring[1] );
        }
        //cout << "new generation created!\n\n";

        for(int i{}; i < get_pop(); i++)
            set_path( new_gen[i], i );

        order();
        set_ave_fit( round( get_pop()/2 ) );
        best_fit.push_back( get_fit(0) );      
    }
}



Travel_Salesman create_population(){
    cout << "\033[1m\033[36m" << "\n... Booting up ...\n" << "\033[0m";

    Random citypos;
    Random geffrey;
    wake_up_geffrey(citypos);
    wake_up_geffrey(geffrey);

    bool circle, elite;
    unsigned int N, city_dim;
    double prob1, prob2;

    ifstream Readinput;
    Readinput.open("input.in");
    if(!Readinput){
        cout << "Error in file opening!\n";
        exit (-1);
    }

    Readinput >> circle;
    if(circle)
        cout << "Creating cities fixed on a circle\n";
    else
        cout << "Creating cities inside a square\n";
    
    Readinput >> N >> city_dim;
    cout << N << " cities in " << city_dim << "D" << endl;
    Readinput >> prob1 >> prob2;
    cout << "Crossover probability:\t\t" << prob1 << endl;
    cout << "Mutation probability:\t\t" << prob2 << endl;
    Readinput >> elite;
    if(elite)
        cout << "\nElite Genetic Algorithm\n";

    vector<city> starter;
    
    if(circle){             //circle of radius 1
        //ofstream circle_towns;
        //circle_towns.open("circle_towns.dat", ios::app);
        for(unsigned int i{}; i < N; i++){
            double theta = citypos.Rannyu( 0, 2*M_PI );
            double x = cos(theta);
            double y = sin(theta);
            city town(i+1, x, y);
            starter.push_back(town);
            //circle_towns << i+1 << "," << x << "," << y << endl;
        }
    }
    else{                   //inside a square of side 1
        //ofstream square_towns;
        //square_towns.open("square_towns.dat", ios::app);
        for(unsigned int i{}; i < N; i++){
            double x = citypos.Rannyu();
            double y = citypos.Rannyu();
            city town(i+1, x, y);
            starter.push_back(town);  
            //square_towns << i+1 << "," << x << "," << y << endl;  
        }
    }

    Travel_Salesman popolazione(N, starter, geffrey);
    popolazione.set_crossover_prob(prob1);
    popolazione.set_mutation_prob(prob2);
    popolazione.set_elite(elite);

    vector<city> copy = starter;

    for(unsigned int i{}; i < N*N; i++){    //34*34 = 1156 << 33! quindi improbabile avere stesso path due volte
        shuffle(copy.begin()+1, copy.end(), std::default_random_engine(1));
        path gene(copy); //automatically compute the fit
        popolazione.add_path(gene);
    }

    cout << "\nPopulation is made up of:\t\t" << popolazione.get_pop() << " paths\n";
    popolazione.check();
    cout << "\nCheck passed successfully\n";

    Readinput.close();
    cout << "\033[1m\033[36m" << "\n... Bootalicious ...\n\n" << "\033[0m";

    return popolazione;
}