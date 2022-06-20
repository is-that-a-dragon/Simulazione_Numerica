#include "city.h"

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include <random>

using namespace std;

void wake_up_geffrey(Random &geffrey){
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
}

void Travel_Salesman::check(){
    for(unsigned int i{}; i < m_population.size(); i++){
        if( m_population[i].size() != N_cities ){
            cout << "\033[1m\033[31m" << "one path in population has more than " << N_cities << " cities\n"<<"\033[0m";
            exit(-1);
        }
        if( m_population[i][0] != m_start[0] ){
            cout << "\033[1m\033[31m"<< "one path has a different starting city\n"<< "\033[0m";
            exit(-1);
        }
        for(unsigned int j{}; j < m_start.size(); j++ ){
            for(unsigned int k=j+1; k < m_start.size(); k++ ){
                if( m_population[i][j] == m_population[i][k] ){
                    cout << "\033[1m\033[31m"<< "one path has the same city two times\n"<< "\033[0m";
                    exit(-1);
                }
            }
        }
    }
    //cout << "\033[1m\033[32m" << "\nchecked\n" << "\033[0m";
}

double Travel_Salesman::tomtom(int A, int B){
    double dist{};
    for(int i{}; i<2; i++)
        dist += pow( m_city_pos[i][A-1] - m_city_pos[i][B-1] , 2 );
    return dist; 
}

vector<double> Travel_Salesman::fit(){
    vector<double> costi;
    double path_fit;
    for(int i{}; i<get_pop(); i++){
        for(unsigned int j{}; j<N_cities-1; j++)
            path_fit += tomtom( m_population[i][j], m_population[i][j+1] );       //distance between cities
        path_fit += tomtom( m_population[i][N_cities-1], m_population[i][0] );    //add first and last

        costi.push_back( path_fit );
        path_fit = 0. ;     //reset to zero for new path
    }
    return costi;
}

double Travel_Salesman::fit(vector<int> path){
    //cannot say vector<int> path = m_population[0] bc than both fit function are w/out arguments
    double costo{}; 
    for(unsigned int i{}; i<N_cities-1; i++)
        costo += tomtom( path[i], path[i+1]);
    costo += tomtom( path[N_cities-1], path[0] );
    return costo;
}

void Travel_Salesman::set_ave_fit(int index){
    double ave_fit{};
    order();
    vector<double> costi = fit();
    ave_fit = accumulate( costi.begin(), costi.begin()+index, 0 );
    average_fit.push_back( (double) ave_fit / index );
}

void Travel_Salesman::order(){
    vector<double> messy_cost = fit();         //unsorted costs
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

    vector< vector <int> > copy_pop = m_population;
    for(int i{}; i<get_pop(); i++)
        m_population[i] = copy_pop[ sorted_indexes[i] ];

    check(); 
}

void Travel_Salesman::print_best(string filename){
    ofstream best;
    best.open(filename, ios::app);
    order();
    int index;
    for(unsigned int i{}; i < N_cities; i++){
        index = m_population[0][i];
        best << m_population[0][i] << "," << m_city_pos[0][index-1] << "," << m_city_pos[1][index-1] << endl;
    }
    vector<double> costo = fit();
    cout << "\033[1m\033[37m" << "best path fit:\t\t" << costo[0] << "\033[0m" << endl;
    best.close();
}

vector<vector<int>> Travel_Salesman::select(vector<vector<int>> parents){
    double r1 = m_geffrey.Rannyu();
    double r2 = m_geffrey.Rannyu();

    order();

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

vector<vector<int>> Travel_Salesman::crossover(vector<vector<int>> parents){
    unsigned int index_cut = static_cast<int>(m_geffrey.Rannyu() * N_cities);
    while (index_cut == 0 or index_cut == N_cities)
        index_cut = static_cast<int>(m_geffrey.Rannyu() * N_cities);

    vector<int> father_head;
    vector<int> mother_head;
    vector<int> father_tail;
    vector<int> mother_tail;
    for(unsigned int i{}; i<N_cities-index_cut; i++){
        father_head.push_back( parents[0][i] );
        mother_head.push_back( parents[1][i] );
    }

    for(unsigned int i{}; i < N_cities; i++){
        auto a = find( father_head.begin(), father_head.end(), parents[1][i] );
        if(a == father_head.end() ) //aka if the city in mother is not present in father_head
            father_tail.push_back( parents[1][i] );

        auto b = find( mother_head.begin(), mother_head.end(), parents[0][i] );
        if(b == mother_head.end() ) //aka if the city in father is not present in mother_head
            mother_tail.push_back( parents[0][i] );
    }

    //merge head and tail
    for(unsigned int i{}; i < father_tail.size(); i++){
        father_head.push_back( father_tail[i] );
        mother_head.push_back( mother_tail[i] );
    }
    parents[0] = father_head;
    parents[1] = mother_head;

    return parents;
}

vector<int> Travel_Salesman::mutate(vector <int> gene){
    vector<int> mutated = gene;
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
    
    return mutated;
}

void Travel_Salesman::evolve(){
    if(m_elite){
        vector<vector<int>> new_gen;
        order();
        int keep = static_cast<int>(get_pop()/2.);
        for (int i{};  i< get_pop()-keep; i+=2){
            vector<vector<int>> parents;
            vector<vector<int>> offspring;

            parents = select(parents);
            //faccio comunque il select su tutta la popolazione cosi da non cadere in minimi

            if(m_geffrey.Rannyu() < crossover_prob)
                offspring = crossover(parents);
            else
                offspring = parents;

            offspring[0] = mutate( offspring[0] );
            offspring[1] = mutate( offspring[1] );

            new_gen.push_back( offspring[0] );
            new_gen.push_back( offspring[1] );
        }

        for(int i{keep}; i < get_pop(); i++)
            m_population[i] = new_gen[i-keep];
        order();
        set_ave_fit( round( get_pop()/2. ) );
        best_fit.push_back( fit( m_population[0] ) ); 
    }

    else{
        vector<vector<int>> new_gen;
        order();
        for (int i{};  i< get_pop(); i+=2){
            vector<vector<int>> parents;
            vector<vector<int>> offspring;

            parents = select(parents);

            if(m_geffrey.Rannyu() < crossover_prob)
                offspring = crossover(parents);
            else
                offspring = parents;
        
            offspring[0] = mutate( offspring[0] );
            offspring[1] = mutate( offspring[1] );

            new_gen.push_back( offspring[0] );
            new_gen.push_back( offspring[1] );
        }
        m_population = new_gen;
        order();
        set_ave_fit( round( get_pop()/2. ) );
        best_fit.push_back( fit( m_population[0] ) );      
    }
}



Travel_Salesman create_population(){
    cout << "\033[1m\033[36m" << "\n... Booting up ...\n" << "\033[0m";

    Random citypos;
    Random geffrey;
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

    vector<int> starter(N, 0);
    vector< vector <double> > posizioni( N, vector<double>(2,0.) );
    int s;
    double p1, p2;
    
    if(circle){             //circle of radius 1
        ifstream in;
        in.open("circle_towns.dat");
        int i{};
        while(in >> s >> p1 >> p2){
            starter[i] = s;
            posizioni[i][0] = p1;
            posizioni[i][1] = p2;
            i++;
        }
        in.close();
    }
    else{                   //inside a square of side 1
        ifstream in;
        in.open("square_towns.dat");
        int i{};
        while( in >> s >> p1 >> p2 ){
            starter[i] = s;
            posizioni[i][0] = p1;
            posizioni[i][1] = p2;
            i++;
        }
        in.close();
    }
    Travel_Salesman popolazione(N, starter, posizioni, geffrey);
    popolazione.set_crossover_prob(prob1);
    popolazione.set_mutation_prob(prob2);
    popolazione.set_elite(elite);

    vector<int> copy = starter;
    for(unsigned int i{}; i < N*N; i++){    //34*34 = 1156 << 33! quindi improbabile avere stesso path due volte
        shuffle(copy.begin()+1, copy.end(), std::default_random_engine(1));
        popolazione.add_path(copy);
    }

    cout << "\nPopulation is made up of:\t\t" << popolazione.get_pop() << " paths\n";
    popolazione.order();
    popolazione.check();
    cout << "\nCheck passed successfully\n";
    Readinput.close();
    cout << "\033[1m\033[36m" << "\n... Bootalicious ...\n\n" << "\033[0m";

    return popolazione;
}

Travel_Salesman US_population(){
    cout <<"\033[1m\033[37m"<<"\n... Bo"<<"\033[1m\033[34m"<<"oting"<<"\033[1m\033[31m"<<"up ...\n"<< "\033[0m";

    Random citypos;
    Random geffrey;
    wake_up_geffrey(geffrey);

    bool elite;
    unsigned int N{50};
    unsigned int city_dim{2};
    double prob1, prob2;
    unsigned int pop_size{300};

    ifstream Readinput;
    Readinput.open("US_input.in");
    if(!Readinput){
        cout << "Error in opening input!\n";
        exit (-1);
    }

    cout << N << " cities in " << city_dim << "D" << endl;
    Readinput >> prob1 >> prob2;
    cout << "Crossover probability:\t\t" << prob1 << endl;
    cout << "Mutation probability:\t\t" << prob2 << endl;
    Readinput >> elite;
    if(elite)
        cout << "\nElite Genetic Algorithm\n";
    Readinput.close();

    Readinput.open("US_capitals.dat");
    if(!Readinput){
        cout << "Error in opening positions!\n";
        exit (-1);
    }

    vector<int> starter(N, 0);
    vector< vector <double> > posizioni( 2, vector<double>(N,0.) );
    int s{1}, i{};
    double p1, p2;
    string state;
    string capital;

    while(Readinput >> state >> capital >> p1 >> p2){
        starter[i] = s;
        posizioni[0][i] = p1;
        posizioni[1][i] = p2;
        i ++;
        s ++;
    }
    Readinput.close();

    Travel_Salesman popolazione(N, starter, posizioni, geffrey);
    popolazione.set_crossover_prob(prob1);
    popolazione.set_mutation_prob(prob2);
    popolazione.set_elite(elite);

    vector<int> copy = starter;
    for(unsigned int i{}; i < pop_size; i++){    //50*50 = 2500 << 50! quindi improbabile avere stesso path due volte
        shuffle(copy.begin()+1, copy.end(), std::default_random_engine(1));
        popolazione.add_path(copy);
    }

    cout << "\nPopulation is made up of:\t\t" << popolazione.get_pop() << " paths\n";
    popolazione.order();
    popolazione.check();
    cout << "\nCheck passed successfully\n";
    cout <<"\033[1m\033[37m"<<"\n... Boo"<<"\033[1m\033[34m"<<"talic"<<"\033[1m\033[31m"<<"ious ...\n\n"<< "\033[0m";

    return popolazione;
}

void read_map( vector<double>& x, vector<double>& y, vector<int>& starter){
    ifstream map;
    map.open("US_capitals.dat");
    if(!map){
        cout << "Error in opening positions!\n";
        exit (-1);
    }
    
    string state, capital;
    double p1, p2;
    int s{1};
    int i{};
    while(map >> state >> capital >> p1 >> p2){
        starter[i] = s;
        x[i] = p1;
        y[i] = p2;
        s++;
        i++;
    }
    map.close();
}