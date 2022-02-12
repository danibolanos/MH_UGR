#ifndef MEMETICOS_H
#define MEMETICOS_H

#include "KNN.h"
#include "Evolutivos.h"

const double PB_PLS = 0.1;

void BL_MEM(const vector<FicheroCSV*>& train, Cromosoma& h, int& num_eval);
void AM_BLX_10_1_0(const vector<FicheroCSV*>& train, vector<double>& w);
void AM_BLX_10_0_1(const vector<FicheroCSV*>& train, vector<double>& w);
void AM_BLX_10_0_1_mej(const vector<FicheroCSV*>& train, vector<double>& w);
Resultados ejecutarAM_BLX_10_1_0(vector<vector<FicheroCSV*>>& particion, int num_part, const int& SEED);
Resultados ejecutarAM_BLX_10_0_1(vector<vector<FicheroCSV*>>& particion, int num_part, const int& SEED);
Resultados ejecutarAM_BLX_10_0_1_mej(vector<vector<FicheroCSV*>>& particion, int num_part, const int& SEED);

#endif
