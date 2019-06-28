#ifndef RELIEF_H
#define RELIEF_H

#include "KNN.h"

void nearestFriendEnemy(const vector<FicheroCSV*>& train, const FicheroCSV& actual, int& posFriend, int& posEnemy);
void relief(const vector<FicheroCSV*>& train, vector<double>& w);
Resultados ejecutarRelief(vector<vector<FicheroCSV*>>& particion, int num_part);

#endif
