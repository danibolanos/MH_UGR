#include "KNN.h"

double agregado(const double& t_clas, const double& t_red){
  return ALPHA*t_clas+(1.0-ALPHA)*t_red;
}

int nearestNeighbour_LOO(const vector<FicheroCSV*>& train, const FicheroCSV& actual, const vector<double>& w){
  int pos_nearest_neighbour = -1;
  double distancia_actual;
  double mejor_distancia = train[0]->traits.size() + 1;

  for (int i=0; i < train.size(); i++){
    if (&actual != train[i]){
      distancia_actual = euclideanDistance(train[i]->traits, actual.traits, w);
          if (distancia_actual < mejor_distancia){
              mejor_distancia = distancia_actual;
              pos_nearest_neighbour = i;
          }
    }
  }

  return pos_nearest_neighbour;
}

int nearestNeighbour(const vector<FicheroCSV*>& train, const FicheroCSV& actual, const vector<double>& w){
  int pos_nearest_neighbour = -1;
  double distancia_actual;
  double mejor_distancia = train[0]->traits.size() + 1;

  for (int i=0; i < train.size(); i++){
      distancia_actual = euclideanDistance(train[i]->traits, actual.traits, w);
      if (distancia_actual < mejor_distancia){
          mejor_distancia = distancia_actual;
          pos_nearest_neighbour = i;
      }
  }

  return pos_nearest_neighbour;
}

Resultados KNN_LOO(const vector<FicheroCSV*>& train, const vector<FicheroCSV*>& test, const vector<double>& w){
  int pos_nearest_neighbour = -1;
  Resultados results;
  int num_w_menor = 0;
  results.aciertos = 0;
  results.tasa_red = 0;

  for (int i=0; i < test.size(); i++){
      pos_nearest_neighbour = nearestNeighbour_LOO(train, *test[i], w);
      if (train[pos_nearest_neighbour]->category == test[i]->category)
          results.aciertos++;
  }

  for(int i=0; i < w.size(); i++){
    if(w[i] < 0.2)
      num_w_menor++;
  }

  results.tasa_clas = 100.0*((1.0*results.aciertos) / (1.0*test.size()));
  results.tasa_red = 100.0*((1.0*num_w_menor) / (1.0*w.size()));

  return results;
}

Resultados KNN(const vector<FicheroCSV*>& train, const vector<FicheroCSV*>& test, const vector<double>& w){
  int pos_nearest_neighbour = -1;
  Resultados results;
  int num_w_menor = 0;
  results.aciertos = 0;
  results.tasa_red = 0;

  for (int i=0; i < test.size(); i++){
      pos_nearest_neighbour = nearestNeighbour(train, *test[i], w);
      if (train[pos_nearest_neighbour]->category == test[i]->category)
          results.aciertos++;
  }

  for(int i=0; i < w.size(); i++){
    if(w[i] < 0.2)
      num_w_menor++;
  }

  results.tasa_clas = 100.0*((1.0*results.aciertos) / (1.0*test.size()));
  results.tasa_red = 100.0*((1.0*num_w_menor) / (1.0*w.size()));

  return results;
}

Resultados ejecutarKNN(vector<vector<FicheroCSV*>>& particion, int num_part) {
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
  for (int i=0; i < train[0]->n; i++) {
      w[i] = 1.0;
  }

  start_timers();
  results = KNN(train, test, w);
  results.tiempo = elapsed_time(REAL);

  cout << "Partición " << num_part+1 << " %_clas: " << results.tasa_clas << " %_red: " << results.tasa_red
  << " Agr: " << agregado(results.tasa_clas, results.tasa_red) << " T: " << results.tiempo << endl;

  return results;
}
