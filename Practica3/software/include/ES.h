#ifndef ES_H
#define ES_H

#include "KNN.h"

const double MU = 0.3;
const double PHI = 0.3;
const double VAL_TF = 0.001;

void ES(const vector<FicheroCSV*>& train, vector<double>& w);
Resultados ejecutarES(vector<vector<FicheroCSV*>>& particion, int num_part, const int& SEED);


#endif
