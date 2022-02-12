#ifndef BL_H
#define BL_H

#include "KNN.h"

const int MAX_NEIGHBOUR = 20;
const int MAX_ITER_ILS = 1000;

void BL(const vector<FicheroCSV*>& train, vector<double>& w);
Resultados ejecutarBL(vector<vector<FicheroCSV*>>& particion, int num_part, const int& SEED);


#endif
