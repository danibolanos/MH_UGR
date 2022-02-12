#ifndef BL_H
#define BL_H

#include "KNN.h"

const int MAX_NEIGHBOUR = 20;

void BL(const vector<FicheroCSV*>& train, vector<double>& w);
Resultados ejecutarBL(vector<vector<FicheroCSV*>>& particion, int num_part, const int& SEED);


#endif
