#include "IntegralMC.h"

//METODO HIT OR MISS: //in ingresso estremi di integrazione, max funzione nell'intervallo, funzione da integrare, numero punti da usare
double IntegralMC::IntegralHoM(double xmin, double xmax, double fmax, const FunzioneBase *f, int npunti) {

	int Ntot{0};
	int Nhit{0};
	double x{};
	double y{};

	for ( int i{0}; i < npunti; ++i ) {
		x = geffrey.Rannyu(xmin, xmax);
		y = geffrey.Rannyu(0, fmax);

		if ( y <= f->Eval(x) ) Nhit++;
		Ntot++;
	}
	return (xmax - xmin) * fmax * (double) Nhit / Ntot;
}
//by reference
double IntegralMC::IntegralHoM(double xmin, double xmax, double fmax, const FunzioneBase &f, int npunti) {

	int Ntot{0};
	int Nhit{0};
	double x{};
	double y{};

	for ( int i{0}; i < npunti; ++i ) {
		x = geffrey.Rannyu(xmin, xmax);
		y = geffrey.Rannyu(0, fmax);

		if ( y <= f.Eval(x) ) Nhit++;
		Ntot++;
	}
	return (xmax - xmin) * fmax * (double) Nhit / Ntot;
}

//METODO DELLA MEDIA: in ingresso estremi di integrazione, funzione da integrare, numero di punti da usare
double IntegralMC::IntegralAVE(double xmin, double xmax, const FunzioneBase *f, int npunti) {

	double sum{0};

	for ( int i{0}; i < npunti; i++ )
		sum += f->Eval(geffrey.Rannyu(xmin, xmax));

	return (xmax - xmin) * sum / (double) npunti;
}
//by reference
double IntegralMC::IntegralAVE(double xmin, double xmax, const FunzioneBase &f, int npunti) {

	double sum{0};

	for ( int i{0}; i < npunti; i++ )
		sum += f.Eval(geffrey.Rannyu(xmin, xmax));
    
	return (xmax - xmin) * sum / (double) npunti;
}

//IMPORTANCE SAMPLING
double IntegralMC::IntegralSIMP(double xmin, double xmax, double fmax, const FunzioneBase *f, int npunti){
    
	double sum{0};

	for ( int i{0}; i < npunti; i++ )
		sum += f->Eval(geffrey.simp(xmin, xmax, fmax));

	return sum/(double) npunti;
}
//by reference
double IntegralMC::IntegralSIMP(double xmin, double xmax, double fmax, const FunzioneBase &f, int npunti){
    
	double sum{0};

	for ( int i{0}; i < npunti; i++ )
		sum += f.Eval(geffrey.simp(xmin, xmax, fmax));

	return sum/(double) npunti;
}

double IntegralMC::IntegraleIS_bad(double xmin, double xmax, FunzioneBase *f, int npunti) {

	double sum = 0;

	for ( int i{0}; i < npunti ; i++ )
		sum += f->Eval(geffrey.simp_bad(xmin, xmax, 1.)); //1 it's the distribution max

	return sum/(double)npunti;
}

double IntegralMC::IntegraleIS_bad(double xmin, double xmax, FunzioneBase &f, int npunti) {

	double sum = 0;

	for ( int i{0}; i < npunti ; i++ )
		sum += f.Eval(geffrey.simp_bad(xmin, xmax, 1.)); //1 it's the distribution max

	return sum/(double)npunti;
}