#include "DE.h"

normal_distribution<double> normal_DE(0.0, SIGMA);
uniform_real_distribution<double> uniform_DE(0.0,1.0);
uniform_int_distribution<int> random_int(0, TAM_PBL-1);
default_random_engine generator_DE;

vector<int> select_3_parents(default_random_engine& gen){
  vector<int> index(3);

  index[0] = random_int(generator_DE);
  do{
    index[1] = random_int(generator_DE);
  }while(index[0]==index[1]);
  do{
    index[2] = random_int(generator_DE);
  }while(index[2]==index[0] || index[2]==index[1]);

  return index;
}

vector<int> select_2_parents(default_random_engine& gen){
  vector<int> index(2);

  index[0] = random_int(generator_DE);
  do{
    index[1] = random_int(generator_DE);
  }while(index[0]==index[1]);

  return index;
}

void DE_rand(const vector<FicheroCSV*>& train, vector<double>& w){
  uniform_int_distribution<int> random_int_mut(0, w.size()-1);
  vector<Cromosoma> offspring(TAM_PBL), poblacion(TAM_PBL), padres(3);
  vector<double> x_pesos;
  vector<int> index(3);
  int num_eval=0, mejor=0;
  Resultados results;

  x_pesos.resize(w.size());

  // Inicializamos la población inicial y guardamos la pos del mejor
  for(int i=0; i < TAM_PBL; i++){
    poblacion[i].w.resize(w.size());
    for(int j=0; j < w.size(); j++){
      poblacion[i].w[j] = uniform_DE(generator_DE);
    }
    results = KNN_LOO(train, train, poblacion[i].w);
    poblacion[i].pts = agregado(results.tasa_clas, results.tasa_red);
    if(poblacion[i].pts > poblacion[mejor].pts)
      mejor = i;
    num_eval++;
  }

  // Condición de parada
  while(num_eval < MAX_ITER){
    for(int i=0; i < TAM_PBL; i++){
      // Generamos los 3 padres aleatorios excluyentes
      index = select_3_parents(generator_DE);
      padres[0] = poblacion[index[0]];
      padres[1] = poblacion[index[1]];
      padres[2] = poblacion[index[2]];

      // Seleccionamos una caracteristica random
      int random = random_int_mut(generator_DE);

      for(int j=0; j < w.size(); j++){
        // Recombinación binomial
        if(uniform_DE(generator_DE) < CR || random == j){
          double val = padres[0].w[j]+F*(padres[1].w[j]-padres[2].w[j]);
          if (val > 1.0)
            val = 1.0;
          else if (val < 0.0)
            val = 0.0;
          x_pesos[j] = val;
        }
        else
          x_pesos[j] = poblacion[i].w[j];
      }
      offspring[i].w = x_pesos;
    }

    // Evaluamos los nuevos hijos
    for (int i=0; i < TAM_PBL; i++){
      results = KNN_LOO(train, train, offspring[i].w);
      offspring[i].pts = agregado(results.tasa_clas, results.tasa_red);
      num_eval++;
    }

    // Sustituyo los hijos mejores en la población
    for(int i=0; i < TAM_PBL; i++){
      if(offspring[i].pts > poblacion[i].pts)
        poblacion[i] = offspring[i];
    }

    //Obtengo el mejor de la nueva poblacion
    for(int i=0; i < TAM_PBL; i++){
      if(poblacion[i].pts > poblacion[mejor].pts){
        mejor = i;
      }
    }
  }

  w = poblacion[mejor].w;
}

void DE_best(const vector<FicheroCSV*>& train, vector<double>& w){
  uniform_int_distribution<int> random_int_mut(0, w.size()-1);
  vector<Cromosoma> offspring(TAM_PBL), poblacion(TAM_PBL), padres(2);
  vector<double> x_pesos;
  vector<int> index(2);
  int num_eval=0, mejor=0;
  Resultados results;

  x_pesos.resize(w.size());

  // Inicializamos la población inicial y guardamos la pos del mejor
  for(int i=0; i < TAM_PBL; i++){
    poblacion[i].w.resize(w.size());
    for(int j=0; j < w.size(); j++){
      poblacion[i].w[j] = uniform_DE(generator_DE);
    }
    results = KNN_LOO(train, train, poblacion[i].w);
    poblacion[i].pts = agregado(results.tasa_clas, results.tasa_red);
    if(poblacion[i].pts > poblacion[mejor].pts)
      mejor = i;
    num_eval++;
  }

  // Condición de parada
  while(num_eval < MAX_ITER){
    for(int i=0; i < TAM_PBL; i++){
      // Generamos los 3 padres aleatorios excluyentes
      index = select_2_parents(generator_DE);
      padres[0] = poblacion[index[0]];
      padres[1] = poblacion[index[1]];

      // Seleccionamos una caracteristica random
      int random = random_int_mut(generator_DE);

      for(int j=0; j < w.size(); j++){
        // Recombinación binomial
        if(uniform_DE(generator_DE) < CR || random == j){
          double val = poblacion[i].w[j]+F*(poblacion[mejor].w[j]-poblacion[i].w[j])+F*(padres[0].w[j]-padres[1].w[j]);
          if (val > 1.0)
            val = 1.0;
          else if (val < 0.0)
            val = 0.0;
          x_pesos[j] = val;
        }
        else
          x_pesos[j] = poblacion[i].w[j];
      }
      offspring[i].w = x_pesos;
    }

    // Evaluamos los nuevos hijos
    for (int i=0; i < TAM_PBL; i++){
      results = KNN_LOO(train, train, offspring[i].w);
      offspring[i].pts = agregado(results.tasa_clas, results.tasa_red);
      num_eval++;
    }

    // Sustituyo los hijos mejores en la población
    for(int i=0; i < TAM_PBL; i++){
      if(offspring[i].pts > poblacion[i].pts)
        poblacion[i] = offspring[i];
    }

    //Obtengo el mejor de la nueva poblacion
    for(int i=0; i < TAM_PBL; i++){
      if(poblacion[i].pts > poblacion[mejor].pts){
        mejor = i;
      }
    }
  }

  w = poblacion[mejor].w;
}

Resultados ejecutarDE_rand(vector<vector<FicheroCSV*>>& particion, int num_part, const int& SEED){
	generator_DE = default_random_engine(SEED);

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
  DE_rand(train, w);
  results = KNN(train, test, w);
  results.tiempo = elapsed_time(REAL);

  cout << "Partición " << num_part+1 << " %_clas: " << results.tasa_clas << " %_red: " << results.tasa_red
  << " Agr: " << agregado(results.tasa_clas, results.tasa_red) << " T: " << results.tiempo << endl;

  return results;
}

Resultados ejecutarDE_best(vector<vector<FicheroCSV*>>& particion, int num_part, const int& SEED){
	generator_DE = default_random_engine(SEED);

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
  DE_best(train, w);
  results = KNN(train, test, w);
  results.tiempo = elapsed_time(REAL);

  cout << "Partición " << num_part+1 << " %_clas: " << results.tasa_clas << " %_red: " << results.tasa_red
  << " Agr: " << agregado(results.tasa_clas, results.tasa_red) << " T: " << results.tiempo << endl;

  return results;
}
