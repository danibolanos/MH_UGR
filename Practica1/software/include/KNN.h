#ifndef KNN_H
#define KNN_H

#include "Funciones.h"

const double ALPHA = 0.5;

int nearestNeighbour(const vector<FicheroCSV*>& train, const FicheroCSV& actual, const vector<double>& w);
Resultados KNN(const vector<FicheroCSV*>& train, const vector<FicheroCSV*>& test, const vector<double>& w);
Resultados ejecutarKNN(vector<vector<FicheroCSV*>>& particion, int num_part);
double agregado(const double& t_clas, const double& t_red);

#endif
