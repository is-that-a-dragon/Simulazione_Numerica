#pragma once
#include <cmath>
#include <vector>

using namespace std;


class FunzioneBase {

    public:
        virtual double Eval(double x) const = 0;
};


class Integranda : public FunzioneBase {

    public:
        Integranda () {};
        ~Integranda() {};

        double Eval(double x)  const override {return M_PI/2. * cos(M_PI*x/2.);};
};

class Integranda_simp : public FunzioneBase {

    public:
        Integranda_simp () {};
        ~Integranda_simp () {};
        // g(x) = f(x) / p(x)
        double Eval(double x) const override {return (M_PI/3.) * cos(M_PI*x/2.)/(1-pow(x,2.)); } ;
};


class Integranda_IS_bad : public FunzioneBase {
    public:
        Integranda_IS_bad () {};
        ~Integranda_IS_bad () {};

        double Eval( double x ) const override { return M_PI/(4. * cos(M_PI*x/2.)); } ;
};