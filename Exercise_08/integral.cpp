#include "integral.h"

using namespace std;

Integral::Integral(const function <double(vector<double>)> f){
    m_f = f;
    m_prec = 1.E-6;
    wake_up_geffrey(m_geffrey);
}

Integral::Integral(double a, double b, const function <double(vector<double>)> f){
    m_f = f;
    m_a = min(a,b);
    m_b = max(a,b);
    if ( a > b) m_sign = -1;
    else m_sign = 1;
    m_prec = 1.E-6;
    wake_up_geffrey(m_geffrey);
}

Integral::Integral(double a, double b, double prec, const function <double(vector<double>)> f){
    m_f = f;
    m_a = min(a,b);
    m_b = max(a,b);
    if ( a > b) m_sign = -1;
    else m_sign = 1;
    m_prec = prec;
    wake_up_geffrey(m_geffrey);
}


double Integral::Metropolis(int punti,const function <double(vector<double>)> f, vector <double> var, double dx, int neq){
    m_sum = 0.;
    attempted = 0;
    accepted = 0;
    //equilibration before integral
    for (int i=0; i<neq; i++){
        double x1 = m_geffrey.Rannyu(var[0]-dx, var[0]+dx);     //new point to test
        vector<double> var1(var);                               //copy of variables
        var1[0] = x1;                                           //set new starting point
        if( f(var) < f(var1) )                                  //evaluate if better
            var[0] = x1;
        else                                                    //random change
            if ( m_geffrey.Rannyu() < ( (double)f(var1)/f(var) ) )
                var[0]=x1;
    }
    //integral
    for(int i=0; i<punti; i++){                                 
        m_sum += m_f(var);                                      //evaluate integral
        double x1 = m_geffrey.Rannyu(var[0]-dx, var[0]+dx);     //metropolis change as before
        vector <double> var1(var);
        var1[0] = x1;
        attempted += 1;
        if( f(var) < f(var1) ) {
            var[0] = x1;
            accepted += 1;
        }
        else{
            if ( m_geffrey.Rannyu() < ( (double)f(var1)/f(var) ) ) {
                var[0] = x1;
                accepted += 1;
            }
        }

    } 
    m_integral = (double)m_sum/punti;
    //cout<<"Acceptance rate: "<<(double)accepted/(double)attempted<<endl;
    //cout << "-------------------------------------------------\n\n";

    return m_integral;
}