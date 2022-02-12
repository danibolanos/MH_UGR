#include "BL.h"

void BL(const vector<FicheroCSV*>& train, vector<double>& w, const int& SEED){
  normal_distribution<double> normal(0.0, SIGMA);
  uniform_real_distribution<double> uniform(0.0,1.0);
  default_random_engine generator(SEED);
  vector<int> index;
  vector<double> w_mut;
  double agr_ant, agr_new;
  Resultados antiguo, nuevo;
  int iter = 0, neighbour = 0, aux = 0;

  index.resize(w.size());

  for(int i=0; i < w.size(); i++){
    index[i] = i;
    w[i] = uniform(generator);
  }

  shuffle(index.begin(), index.end(), generator);

  antiguo = KNN(train, train, w);
  agr_ant = agregado(antiguo.tasa_clas, antiguo.tasa_red);

  while (iter < MAX_ITER && neighbour < w.size()*MAX_NEIGHBOUR){
      aux = index[iter % w.size()];
      w_mut = w;
      w_mut[aux] += normal(generator);
      if (w_mut[aux] > 1)
        w_mut[aux] = 1;
      else if (w_mut[aux] < 0)
        w_mut[aux] = 0;
      nuevo = KNN(train, train, w_mut);
      agr_new = agregado(nuevo.tasa_clas, nuevo.tasa_red);
      iter++;
      if (agr_new > agr_ant) {
        w = w_mut;
        agr_ant = agr_new;
        neighbour = 0;
      }
      else
        neighbour++;
      if(iter % w.size()==0)
        shuffle(index.begin(), index.end(), generator);
  }
}

Resultados ejecutarBL(vector<vector<FicheroCSV*>>& particion, int num_part, const int& SEED){
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
  BL(train, w, SEED);
  results = KNN(train, test, w);
  results.tiempo = elapsed_time(REAL);

  cout << "Partición " << num_part+1 << " %_clas: " << results.tasa_clas << " %_red: " << results.tasa_red
  << " Agr: " << agregado(results.tasa_clas, results.tasa_red) << " T: " << results.tiempo << endl;

  return results;
}
