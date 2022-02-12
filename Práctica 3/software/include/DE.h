#ifndef DE_H
#define DE_H

#include "KNN.h"

const int TAM_PBL = 50;   // Tamaño población cromosomas
const double CR = 0.5;
const double F = 0.5;

struct Cromosoma{
	vector<double> w;
	double pts;
};

vector<int> select_3_parents(default_random_engine& gen);
vector<int> select_2_parents(default_random_engine& gen);
void DE_rand(const vector<FicheroCSV*>& train, vector<double>& w);
void DE_best(const vector<FicheroCSV*>& train, vector<double>& w);
Resultados ejecutarDE_rand(vector<vector<FicheroCSV*>>& particion, int num_part, const int& SEED);
Resultados ejecutarDE_best(vector<vector<FicheroCSV*>>& particion, int num_part, const int& SEED);


#endif
