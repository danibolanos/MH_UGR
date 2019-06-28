#include "ILS.h"

normal_distribution<double> normal_ILS(0.0, SIGMA_ILS);
uniform_real_distribution<double> uniform_ILS(0.0,1.0);
default_random_engine generator_ILS;

void ILS(const vector<FicheroCSV*>& train, vector<double>& w){
  uniform_int_distribution<int> random_int(0, w.size()-1);
  vector<double> w_mejor;
  Resultados antiguo, nuevo;
  double agr_new, agr_mejor;
  int num_mutaciones = ceil(0.1*w.size());

  // Evaluamos la solución inicial
  for(int i=0; i < w.size(); i++)
    w[i] = uniform_ILS(generator_ILS);
  // Aplicamos la BL y tomamos el agr actual como mejor
  BL(train, w);
  antiguo = KNN_LOO(train, train, w);
  agr_mejor = agregado(antiguo.tasa_clas, antiguo.tasa_red);
  w_mejor = w;

  // Hacemos 14 iteraciones
  for(int i=1; i < 15; i++){
    // Mutamos la mejor solución
    for(int j=0; j < num_mutaciones; j++){
      int aux = random_int(generator_ILS);
      w[aux] += normal_ILS(generator_ILS);
      if (w[aux] > 1.0)
        w[aux] = 1.0;
      else if (w[aux] < 0.0)
        w[aux] = 0.0;
    }
    BL(train, w);
    nuevo = KNN_LOO(train, train, w);
    agr_new = agregado(nuevo.tasa_clas, nuevo.tasa_red);
    // Si la nueva solución es mejor se toma esa
    if(agr_new > agr_mejor){
      w_mejor = w;
      agr_mejor = agr_new;
    }
    w = w_mejor;
  }
}

Resultados ejecutarILS(vector<vector<FicheroCSV*>>& particion, int num_part, const int& SEED){
	generator_ILS = default_random_engine(SEED);

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
  ILS(train, w);
  results = KNN(train, test, w);
  results.tiempo = elapsed_time(REAL);

  cout << "Partición " << num_part+1 << " %_clas: " << results.tasa_clas << " %_red: " << results.tasa_red
  << " Agr: " << agregado(results.tasa_clas, results.tasa_red) << " T: " << results.tiempo << endl;

  return results;
}
