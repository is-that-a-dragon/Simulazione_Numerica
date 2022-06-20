#pragma once

#include "random.h"
#include "float.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <functional>
#include <vector>

using namespace std;

class Integral{

    public:
        Integral(double a, double b, const function <double(vector <double>)> f);
        Integral(const function <double(vector<double>)> f);
        Integral(double a, double b, double prec, const function <double(vector<double>)> f);
        
        void SetPrec (double const prec) {
                if(prec >= DBL_EPSILON)
                    m_prec = prec;
                else
                    m_prec = DBL_EPSILON;
        };
        double GetPrec() const {return m_prec;};

        void SetH(double const h) {m_h = h;};
        double GetH() const {return m_h;};
        double Metropolis (int,const function <double(vector <double>)>, vector <double>, double, int);

    private:
        double m_a, m_b;                                        //Estremi di integrazione
        double m_sum;                                           //Integrale senza segno
        double m_h;                                             //Larghezza intervalli
        int m_sign;                                             //Se estremi non sono in ordine
        double m_integral;                                      //Valore dell'integrale
        function <double(vector<double>)> m_f;                  //Funzione integranda
        double m_prec;                                          //Precisione
        Random m_geffrey;
        int accepted, attempted;
};