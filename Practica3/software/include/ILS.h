#ifndef ILS_H
#define ILS_H

#include "BL.h"

const double SIGMA_ILS = 0.4;

void ILS(const vector<FicheroCSV*>& train, vector<double>& w);
Resultados ejecutarILS(vector<vector<FicheroCSV*>>& particion, int num_part, const int& SEED);


#endif
