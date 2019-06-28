#ifndef EVOLUTIVOS_H
#define EVOLUTIVOS_H

#include "Funciones.h"

const int TAM_PBL = 30;   // Tamaño población cromosomas
const int TAM_MEM = 10;		// Tamaño población meméticos
const double ALPHA_AGG = 0.3;
const double PB_CRUCE = 0.7;
const double PB_MUT = 0.001;

struct Cromosoma{
	vector<double> w;
	double pts;
};

vector<int> getListBestCromosoma(const vector<Cromosoma>& poblacion);
pair<Cromosoma, Cromosoma> cruceBLX(const Cromosoma& p1, const Cromosoma& p2, default_random_engine& gen);
Cromosoma cruceArit(const Cromosoma &p1, const Cromosoma &p2);
int binaryTournament(const vector<Cromosoma>& poblacion, default_random_engine& gen);
int getBestCromosoma(const vector<Cromosoma>& poblacion);
int getWorstCromosoma(const vector<Cromosoma>& poblacion);

#endif
