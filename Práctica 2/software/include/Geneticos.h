#ifndef GENETICOS_H
#define GENETICOS_H

#include "KNN.h"
#include "Evolutivos.h"

void AGG_BLX(const vector<FicheroCSV*>& train, vector<double>& w);
void AGE_BLX(const vector<FicheroCSV*>& train, vector<double>& w);
void AGG_CA(const vector<FicheroCSV*>& train, vector<double>& w);
void AGE_CA(const vector<FicheroCSV*>& train, vector<double>& w);
Resultados ejecutarAGG_BLX(vector<vector<FicheroCSV*>>& particion, int num_part, const int& SEED);
Resultados ejecutarAGE_BLX(vector<vector<FicheroCSV*>>& particion, int num_part, const int& SEED);
Resultados ejecutarAGG_CA(vector<vector<FicheroCSV*>>& particion, int num_part, const int& SEED);
Resultados ejecutarAGE_CA(vector<vector<FicheroCSV*>>& particion, int num_part, const int& SEED);

#endif
