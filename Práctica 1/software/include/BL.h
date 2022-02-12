#ifndef BL_H
#define BL_H

#include "KNN.h"

const double SIGMA = 0.3;
const int MAX_ITER = 15000;
const int MAX_NEIGHBOUR = 20;

void BL(const vector<FicheroCSV*>& train, vector<double>& w, const int& SEED);
Resultados ejecutarBL(vector<vector<FicheroCSV*>>& particion, int num_part, const int& SEED);


#endif
