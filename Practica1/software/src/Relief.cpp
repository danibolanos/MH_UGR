#include "Relief.h"

void nearestFriendEnemy(const vector<FicheroCSV*>& train, const FicheroCSV& actual, int& posFriend, int& posEnemy){
  double distancia_actual;
  double mejor_distancia_a, mejor_distancia_e;
  mejor_distancia_a = mejor_distancia_e = train[0]->traits.size() + 1;

  for (int i=0; i < train.size(); i++) {
    if (&actual != train[i]){
      distancia_actual = euclideanDistance(train[i]->traits, actual.traits);
      if(train[i]->category == actual.category){
        if (distancia_actual < mejor_distancia_a){
          mejor_distancia_a = distancia_actual;
          posFriend = i;
        }
      }
      else
        if (distancia_actual < mejor_distancia_e){
          mejor_distancia_e = distancia_actual;
          posEnemy = i;
        }
    }
  }
}

void relief(const vector<FicheroCSV*>& train, vector<double>& w){
  int pos_nearest_friend, pos_nearest_enemy;
  double w_max = -1.0;

  for (int i=0; i < train.size(); i++){
    nearestFriendEnemy(train, *train[i], pos_nearest_friend, pos_nearest_enemy);
    for (int j=0; j < w.size(); j++){
      w[j] += fabs(train[i]->traits[j]-train[pos_nearest_enemy]->traits[j])
      - fabs(train[i]->traits[j]-train[pos_nearest_friend]->traits[j]);
      if(w_max < w[j])
        w_max = w[j];
    }
  }

  for (int i=0; i < w.size(); i++){
    if(w[i] < 0.0)
      w[i] = 0.0;
    else if (w_max != 0)
      w[i] /= w_max;
  }
}

Resultados ejecutarRelief(vector<vector<FicheroCSV*>>& particion, int num_part) {
  vector<FicheroCSV*> train, test;
  vector<string> clasificacion;
  vector<double> w;
  Resultados results;

  test = particion[num_part];

  for(int i=0; i < particion.size(); i++){
		if(i!=num_part){
			for(int j=0; j < particion[i].size(); j++)
				train.push_back(particion[i][j]);
		}
	}

  w.resize(train[0]->n);
  for (int i=0; i < train[0]->n; i++){
      w[i] = 0.0;
  }

  start_timers();
  relief(train, w);
  results = KNN(train, test, w);
  results.tiempo = elapsed_time(REAL);

  cout << "Partición " << num_part+1 << " %_clas: " << results.tasa_clas << " %_red: " << results.tasa_red
  << " Agr: " << agregado(results.tasa_clas, results.tasa_red) << " T: " << results.tiempo << endl;

  return results;
}
