#include "ES.h"

normal_distribution<double> normal_ES(0.0, SIGMA);
uniform_real_distribution<double> uniform_ES(0.0,1.0);
default_random_engine generator_ES;

void ES(const vector<FicheroCSV*>& train, vector<double>& w){
  uniform_int_distribution<int> random_int(0, w.size()-1);
  int max_exitos, max_vecinos, M;
  int num_eval=0, neighbour=0, num_exitos;
  vector<double> w_mut, w_mejor;
  double agr_ant, agr_new, agr_mejor;
  Resultados antiguo, nuevo;
  double To, Tf, beta, Tk;

  max_vecinos = 10*w.size();
  max_exitos = 0.1*max_vecinos;
  M = MAX_ITER / max_vecinos;
  num_exitos = max_exitos;

  // Evaluamos la solución inicial
  for(int i=0; i < w.size(); i++)
    w[i] = uniform_ES(generator_ES);

  antiguo = KNN_LOO(train, train, w);
  agr_ant = agregado(antiguo.tasa_clas, antiguo.tasa_red);
  num_eval++;

  w_mejor = w;
  agr_mejor = agr_ant;

  // MU = PHI = 0.3 en las ejecuciones
  To = (MU*agr_mejor)/(-1.0*log(PHI));
  Tf = VAL_TF;
  beta = (To-Tf)/(M*To*Tf);
  Tk = To;

  //Comprobamos que Tf es menor que To, si no se cumpliese multiplicamos por 0.1 la final
  while (Tf > To)
    Tf = Tf*0.1;


  while(num_eval < MAX_ITER && num_exitos != 0){
    num_exitos = 0;
    neighbour = 0;
    while(num_eval < MAX_ITER && num_exitos < max_exitos && neighbour < max_vecinos){
      // Mutar y truncar el vector de pesos
      int aux = random_int(generator_ES);
      w_mut = w;
      w_mut[aux] += normal_ES(generator_ES);
      if (w_mut[aux] > 1.0)
        w_mut[aux] = 1.0;
      else if (w_mut[aux] < 0.0)
        w_mut[aux] = 0.0;
      // Calculo la tasa de agregado del vector mutado
      nuevo = KNN_LOO(train, train, w_mut);
      agr_new = agregado(nuevo.tasa_clas, nuevo.tasa_red);
      num_eval++;
      neighbour++;
      // Cálculo del incremento de la diff
      double diff = agr_ant-agr_new;
      if(diff < 0 || uniform_ES(generator_ES) <= exp(-1.0*diff/Tk)){
        w = w_mut;
        agr_ant = agr_new;
        num_exitos++;
        if(agr_new > agr_mejor){
          w_mejor = w_mut;
          agr_mejor = agr_new;
        }
      }
    }
    // Enfriamiento (Esquema de Cauchy)
    Tk = Tk/(1.0+beta*Tk);
  }

  w = w_mejor;
}

Resultados ejecutarES(vector<vector<FicheroCSV*>>& particion, int num_part, const int& SEED){
	generator_ES = default_random_engine(SEED);

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
  ES(train, w);
  results = KNN(train, test, w);
  results.tiempo = elapsed_time(REAL);

  cout << "Partición " << num_part+1 << " %_clas: " << results.tasa_clas << " %_red: " << results.tasa_red
  << " Agr: " << agregado(results.tasa_clas, results.tasa_red) << " T: " << results.tiempo << endl;

  return results;
}
