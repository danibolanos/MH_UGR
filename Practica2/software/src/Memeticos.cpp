#include "Memeticos.h"

normal_distribution<double> normal_MEM(0.0, SIGMA);
uniform_real_distribution<double> uniform_MEM(0.0,1.0);
default_random_engine generator_MEM;

void BL_MEM(const vector<FicheroCSV*>& train, Cromosoma& h, int& num_eval){
  vector<int> index;
  vector<double> w_mut;
  double agr_ant, agr_new;
  Resultados antiguo, nuevo;
  int iter = 0, neighbour = 0, aux = 0;

  index.resize(h.w.size());

  for(int i=0; i < h.w.size(); i++)
    index[i] = i;

  shuffle(index.begin(), index.end(), generator_MEM);

  agr_ant = h.pts;

  // Condición de parada
  while (neighbour < h.w.size()*2){
      aux = index[iter % h.w.size()];
      w_mut = h.w;
      w_mut[aux] += normal_MEM(generator_MEM);
      if (w_mut[aux] > 1.0)
        w_mut[aux] = 1.0;
      else if (w_mut[aux] < 0.0)
        w_mut[aux] = 0.0;
      nuevo = KNN_LOO(train, train, w_mut);
      agr_new = agregado(nuevo.tasa_clas, nuevo.tasa_red);
      iter++;
      num_eval++;
      if (agr_new > agr_ant) {
        h.w = w_mut;
        agr_ant = agr_new;
      }
      neighbour++;
      if(iter % h.w.size()==0)
        shuffle(index.begin(), index.end(), generator_MEM);
  }
  h.pts = agr_ant;
}

void AM_BLX_10_1_0(const vector<FicheroCSV*>& train, vector<double>& w) {
  vector<Cromosoma> poblacion(TAM_MEM), padres(TAM_MEM), poblacion_intermedia(TAM_MEM);
  pair<Cromosoma, Cromosoma> hijosBLX;
  int mejor = 0, peorNueva = 0, mejorNueva = 0, select = 0, num_cruces, num_mutaciones;
  int genes_por_generacion, generacion = 0, t = 0;
  int cromosoma_mutar, gen_mutar, valor_mut;
  double aux_mut;
  Resultados results = {0.0};

  num_cruces = floor(PB_CRUCE * (TAM_MEM/2))*2;
  genes_por_generacion = num_cruces * w.size();
  num_mutaciones = ceil(PB_MUT * genes_por_generacion);

  uniform_int_distribution<int> random_mutaciones(0, genes_por_generacion-1);

  for(int i=0; i < TAM_MEM; i++){
    poblacion[i].w.resize(w.size());
    for(int j=0; j < w.size(); j++){
      poblacion[i].w[j] = uniform_MEM(generator_MEM);
    }
    results = KNN_LOO(train, train, poblacion[i].w);
    poblacion[i].pts = agregado(results.tasa_clas, results.tasa_red);
    if(poblacion[i].pts > poblacion[mejor].pts)
      mejor = i;
    t++;
  }

  while(t < MAX_ITER){

    generacion++;

    // Seleccionamos los padres con el torneo binario
    for (int i=0; i < TAM_MEM; i++){
      select = binaryTournament(poblacion, generator_MEM);
      padres[i] = poblacion[select];
    }

    // Hacemos el cruce BLX para obtener los hijos
    for (int i=0; i < num_cruces; i+=2) {
      hijosBLX = cruceBLX(padres[i], padres[i+1], generator_MEM);
      poblacion_intermedia[i] = hijosBLX.first;
      poblacion_intermedia[i+1] = hijosBLX.second;
    }

    // Hacemos el número de mutaciones esperadas
    for(int i=0; i < num_mutaciones; i++){
      valor_mut = random_mutaciones(generator_MEM);
      cromosoma_mutar = valor_mut / w.size();
      gen_mutar = valor_mut % w.size();
      aux_mut = poblacion_intermedia[cromosoma_mutar].w[gen_mutar] + normal_MEM(generator_MEM);
      if (aux_mut > 1)
        aux_mut = 1.0;
      else if (aux_mut < 0)
        aux_mut = 0.0;
      poblacion_intermedia[cromosoma_mutar].w[gen_mutar] = aux_mut;
    }

    // Introducimos los últimos padres en la población
    for (int i=num_cruces; i < TAM_MEM; i++)
      poblacion_intermedia[i] = padres[i];

    // Aplicamos Local Search a toda la población
    if (generacion % 10 == 0){
      for (int i=0; i < TAM_MEM; i++){
        BL_MEM(train, poblacion_intermedia[i], t);
      }
    }

    // Calculamos el mejor y el peor Cromosoma de la nueva población
    mejorNueva = getBestCromosoma(poblacion_intermedia);
    peorNueva = getWorstCromosoma(poblacion_intermedia);

    // Elitismo
    if (poblacion[mejor].pts > poblacion_intermedia[mejorNueva].pts){
      poblacion_intermedia[peorNueva] = poblacion[mejor];
      mejor = peorNueva;
    }
    else
      mejor = mejorNueva;

    // Sustituimos la población
    poblacion = poblacion_intermedia;
  }

  w = poblacion[mejor].w;
}

void AM_BLX_10_0_1(const vector<FicheroCSV*>& train, vector<double>& w) {
  vector<Cromosoma> poblacion(TAM_MEM), padres(TAM_MEM), poblacion_intermedia(TAM_MEM);
  pair<Cromosoma, Cromosoma> hijosBLX;
  int mejor = 0, peorNueva = 0, mejorNueva = 0, select = 0, num_cruces, num_mutaciones, num_bl_cromosoma;
  int genes_por_generacion, generacion = 0, t = 0;
  int cromosoma_mutar, gen_mutar, valor_mut;
  double aux_mut;
  Resultados results = {0.0};

  num_cruces = floor(PB_CRUCE * (TAM_MEM/2))*2;
  genes_por_generacion = num_cruces * w.size();
  num_mutaciones = ceil(PB_MUT * genes_por_generacion);
  num_bl_cromosoma = ceil(PB_PLS * TAM_MEM);

  uniform_int_distribution<int> random_mutaciones(0, genes_por_generacion-1);
  uniform_int_distribution<int> random_BL(0, TAM_MEM-1);

  for(int i=0; i < TAM_MEM; i++){
    poblacion[i].w.resize(w.size());
    for(int j=0; j < w.size(); j++){
      poblacion[i].w[j] = uniform_MEM(generator_MEM);
    }
    results = KNN_LOO(train, train, poblacion[i].w);
    poblacion[i].pts = agregado(results.tasa_clas, results.tasa_red);
    if(poblacion[i].pts > poblacion[mejor].pts)
      mejor = i;
    t++;
  }

  while(t < MAX_ITER){

    generacion++;

    // Seleccionamos los padres con el torneo binario
    for (int i=0; i < TAM_MEM; i++){
      select = binaryTournament(poblacion, generator_MEM);
      padres[i] = poblacion[select];
    }

    // Hacemos el cruce BLX para obtener los hijos
    for (int i=0; i < num_cruces; i+=2) {
      hijosBLX = cruceBLX(padres[i], padres[i+1], generator_MEM);
      poblacion_intermedia[i] = hijosBLX.first;
      poblacion_intermedia[i+1] = hijosBLX.second;
    }

    // Hacemos el número de mutaciones esperadas
    for(int i=0; i < num_mutaciones; i++){
      valor_mut = random_mutaciones(generator_MEM);
      cromosoma_mutar = valor_mut / w.size();
      gen_mutar = valor_mut % w.size();
      aux_mut = poblacion_intermedia[cromosoma_mutar].w[gen_mutar] + normal_MEM(generator_MEM);
      if (aux_mut > 1)
        aux_mut = 1.0;
      else if (aux_mut < 0)
        aux_mut = 0.0;
      poblacion_intermedia[cromosoma_mutar].w[gen_mutar] = aux_mut;
    }

    // Introducimos los últimos padres en la población
    for (int i=num_cruces; i < TAM_MEM; i++)
      poblacion_intermedia[i] = padres[i];

    // Evaluamos de nuevo la poblacion
    for (int i=0; i < num_cruces; i++){
        results = KNN_LOO(train, train, poblacion_intermedia[i].w);
        poblacion_intermedia[i].pts = agregado(results.tasa_clas, results.tasa_red);
        t++;
    }

    // Aplicamos Local Search a un cromosoma aleatorio de la poblacion
    if (generacion % 10 == 0){
      for (int i=0; i < num_bl_cromosoma; i++){
        select = random_BL(generator_MEM);
        BL_MEM(train, poblacion_intermedia[select], t);
      }
    }

    // Calculamos el mejor y el peor Cromosoma de la nueva población
    mejorNueva = getBestCromosoma(poblacion_intermedia);
    peorNueva = getWorstCromosoma(poblacion_intermedia);

    // Elitismo
    if (poblacion[mejor].pts > poblacion_intermedia[mejorNueva].pts){
      poblacion_intermedia[peorNueva] = poblacion[mejor];
      mejor = peorNueva;
    }
    else
      mejor = mejorNueva;

    // Sustituimos la población
    poblacion = poblacion_intermedia;
  }

  w = poblacion[mejor].w;
}

void AM_BLX_10_0_1_mej(const vector<FicheroCSV*>& train, vector<double>& w) {
  vector<Cromosoma> poblacion(TAM_MEM), padres(TAM_MEM), poblacion_intermedia(TAM_MEM);
  pair<Cromosoma, Cromosoma> hijosBLX;
  vector<int> lista;
  int mejor = 0, peorNueva = 0, mejorNueva = 0, select = 0, num_cruces, num_mutaciones, num_bl_cromosoma;
  int genes_por_generacion, generacion = 0, t = 0;
  int cromosoma_mutar, gen_mutar, valor_mut;
  double aux_mut;
  Resultados results = {0.0};

  num_cruces = floor(PB_CRUCE * (TAM_MEM/2))*2;
  genes_por_generacion = num_cruces * w.size();
  num_mutaciones = ceil(PB_MUT * genes_por_generacion);
  num_bl_cromosoma = ceil(PB_PLS * TAM_MEM);

  uniform_int_distribution<int> random_mutaciones(0, genes_por_generacion-1);
  uniform_int_distribution<int> random_BL(0, TAM_MEM-1);

  for(int i=0; i < TAM_MEM; i++){
    poblacion[i].w.resize(w.size());
    for(int j=0; j < w.size(); j++){
      poblacion[i].w[j] = uniform_MEM(generator_MEM);
    }
    results = KNN_LOO(train, train, poblacion[i].w);
    poblacion[i].pts = agregado(results.tasa_clas, results.tasa_red);
    if(poblacion[i].pts > poblacion[mejor].pts)
      mejor = i;
    t++;
  }

  while(t < MAX_ITER){

    generacion++;

    // Seleccionamos los padres con el torneo binario
    for (int i=0; i < TAM_MEM; i++){
      select = binaryTournament(poblacion, generator_MEM);
      padres[i] = poblacion[select];
    }

    // Hacemos el cruce BLX para obtener los hijos
    for (int i=0; i < num_cruces; i+=2) {
      hijosBLX = cruceBLX(padres[i], padres[i+1], generator_MEM);
      poblacion_intermedia[i] = hijosBLX.first;
      poblacion_intermedia[i+1] = hijosBLX.second;
    }

    // Hacemos el número de mutaciones esperadas
    for(int i=0; i < num_mutaciones; i++){
      valor_mut = random_mutaciones(generator_MEM);
      cromosoma_mutar = valor_mut / w.size();
      gen_mutar = valor_mut % w.size();
      aux_mut = poblacion_intermedia[cromosoma_mutar].w[gen_mutar] + normal_MEM(generator_MEM);
      if (aux_mut > 1)
        aux_mut = 1.0;
      else if (aux_mut < 0)
        aux_mut = 0.0;
      poblacion_intermedia[cromosoma_mutar].w[gen_mutar] = aux_mut;
    }

    // Introducimos los últimos padres en la población
    for (int i=num_cruces; i < TAM_MEM; i++)
      poblacion_intermedia[i] = padres[i];

    // Evaluamos la poblacion
    for (int i=0; i < num_cruces; i++){
        results = KNN_LOO(train, train, poblacion_intermedia[i].w);
        poblacion_intermedia[i].pts = agregado(results.tasa_clas, results.tasa_red);
        t++;
    }

    lista = getListBestCromosoma(poblacion_intermedia);

    // Aplicamos Local Search sobre el mejor cromosoma de la nueva población
    if (generacion % 10 == 0){
      for (int i=0; i < num_bl_cromosoma; i++){
        select = lista[i];
        BL_MEM(train, poblacion_intermedia[select], t);
      }
    }

    // Calculamos el peor y mejor Cromosomas de la nueva población
    mejorNueva = getBestCromosoma(poblacion_intermedia);
    peorNueva = getWorstCromosoma(poblacion_intermedia);

    // Elitismo
    if (poblacion[mejor].pts > poblacion_intermedia[mejorNueva].pts){
      poblacion_intermedia[peorNueva] = poblacion[mejor];
      mejor = peorNueva;
    }
    else
      mejor = mejorNueva;

    // Sustituimos la población
    poblacion = poblacion_intermedia;
  }

  w = poblacion[mejor].w;
}

Resultados ejecutarAM_BLX_10_1_0(vector<vector<FicheroCSV*>>& particion, int num_part, const int& SEED){
	generator_MEM = default_random_engine(SEED);

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
  AM_BLX_10_1_0(train, w);
  results = KNN(train, test, w);
  results.tiempo = elapsed_time(REAL);

  cout << "Partición " << num_part+1 << " %_clas: " << results.tasa_clas << " %_red: " << results.tasa_red
  << " Agr: " << agregado(results.tasa_clas, results.tasa_red) << " T: " << results.tiempo << endl;

  return results;
}

Resultados ejecutarAM_BLX_10_0_1(vector<vector<FicheroCSV*>>& particion, int num_part, const int& SEED){
	generator_MEM = default_random_engine(SEED);

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
  AM_BLX_10_0_1(train, w);
  results = KNN(train, test, w);
  results.tiempo = elapsed_time(REAL);

  cout << "Partición " << num_part+1 << " %_clas: " << results.tasa_clas << " %_red: " << results.tasa_red
  << " Agr: " << agregado(results.tasa_clas, results.tasa_red) << " T: " << results.tiempo << endl;

  return results;
}

Resultados ejecutarAM_BLX_10_0_1_mej(vector<vector<FicheroCSV*>>& particion, int num_part, const int& SEED){
	generator_MEM = default_random_engine(SEED);

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
  AM_BLX_10_0_1_mej(train, w);
  results = KNN(train, test, w);
  results.tiempo = elapsed_time(REAL);

  cout << "Partición " << num_part+1 << " %_clas: " << results.tasa_clas << " %_red: " << results.tasa_red
  << " Agr: " << agregado(results.tasa_clas, results.tasa_red) << " T: " << results.tiempo << endl;

  return results;
}
