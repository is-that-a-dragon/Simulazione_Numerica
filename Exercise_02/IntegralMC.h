#pragma once

#include "funzioni.h"
#include "random.h"

class IntegralMC {

public:

	IntegralMC(Random &numbergen) : geffrey{numbergen} {};
	~IntegralMC() {};

	double IntegralHoM(double, double, double, const FunzioneBase*, int);
    	double IntegralHoM(double, double, double, const FunzioneBase&, int);
	double IntegralAVE(double, double, const FunzioneBase*, int);
	double IntegralAVE(double, double, const FunzioneBase&, int);
	double IntegralSIMP(double, double, double, const FunzioneBase*, int);
	double IntegralSIMP(double, double, double, const FunzioneBase&, int);
	double IntegraleIS_bad(double, double, FunzioneBase*, int);
	double IntegraleIS_bad(double, double, FunzioneBase&, int);

private:

	Random geffrey;

};
