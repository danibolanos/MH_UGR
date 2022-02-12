#include "Geneticos.h"

normal_distribution<double> normal_GEN(0.0, SIGMA);
uniform_real_distribution<double> uniform_GEN(0.0,1.0);
default_random_engine generator_GEN;

void AGG_BLX(const vector<FicheroCSV*>& train, vector<double>& w) {
  vector<Cromosoma> poblacion(TAM_PBL), padres(TAM_PBL), poblacion_intermedia(TAM_PBL);
  pair<Cromosoma, Cromosoma> hijosBLX;
  int mejor = 0, peorNueva = 0, mejorNueva = 0, select = 0, num_cruces, num_mutaciones;
  int genes_por_generacion, t = 0;
  int cromosoma_mutar, gen_mutar, valor_mut;
  double aux_mut;
  Resultados results = {0.0};

  num_cruces = floor(PB_CRUCE * (TAM_PBL/2))*2;
  genes_por_generacion = num_cruces * w.size();
  num_mutaciones = ceil(PB_MUT * genes_por_generacion);

  uniform_int_distribution<int> random_mutaciones(0, genes_por_generacion-1);

  for(int i=0; i < TAM_PBL; i++){
    poblacion[i].w.resize(w.size());
    for(int j=0; j < w.size(); j++){
      poblacion[i].w[j] = uniform_GEN(generator_GEN);
    }
    results = KNN_LOO(train, train, poblacion[i].w);
    poblacion[i].pts = agregado(results.tasa_clas, results.tasa_red);
    if(poblacion[i].pts > poblacion[mejor].pts)
      mejor = i;
    t++;
  }

  while(t < MAX_ITER){

    // Seleccionamos los padres con el torneo binario
    for (int i=0; i < TAM_PBL; i++){
      select = binaryTournament(poblacion, generator_GEN);
      padres[i] = poblacion[select];
    }

    // Hacemos el cruce BLX para obtener los hijos
    for (int i=0; i < num_cruces; i+=2) {
      hijosBLX = cruceBLX(padres[i], padres[i+1], generator_GEN);
      poblacion_intermedia[i] = hijosBLX.first;
      poblacion_intermedia[i+1] = hijosBLX.second;
    }

    // Hacemos el número de mutaciones esperadas
    for(int i=0; i < num_mutaciones; i++){
      valor_mut = random_mutaciones(generator_GEN);
      cromosoma_mutar = valor_mut / w.size();
      gen_mutar = valor_mut % w.size();
      aux_mut = poblacion_intermedia[cromosoma_mutar].w[gen_mutar] + normal_GEN(generator_GEN);
      if (aux_mut > 1)
        aux_mut = 1.0;
      else if (aux_mut < 0)
        aux_mut = 0.0;
      poblacion_intermedia[cromosoma_mutar].w[gen_mutar] = aux_mut;
    }

    // Introducimos los últimos padres en la población
    for (int i=num_cruces; i < TAM_PBL; i++)
      poblacion_intermedia[i] = padres[i];

    // Evaluamos los nuevos hijos
    for (int i=0; i < num_cruces; i++){
      results = KNN_LOO(train, train, poblacion_intermedia[i].w);
      poblacion_intermedia[i].pts = agregado(results.tasa_clas, results.tasa_red);
      t++;
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

void AGG_CA(const vector<FicheroCSV*>& train, vector<double>& w) {
  vector<Cromosoma> poblacion(TAM_PBL), padres, poblacion_intermedia(TAM_PBL);
  Cromosoma hijoCA;
  int mejor = 0, peorNueva = 0, mejorNueva = 0, select = 0, num_cruces, num_mutaciones;
  int genes_por_generacion, t = 0;
  int cromosoma_mutar, gen_mutar, valor_mut, cruces_dobles;
  double aux_mut;
  Resultados results = {0.0};

  num_cruces = floor(PB_CRUCE * (TAM_PBL/2))*2;
  genes_por_generacion = num_cruces * w.size();
  num_mutaciones = ceil(PB_MUT * genes_por_generacion);
  cruces_dobles = num_cruces+TAM_PBL;

  padres.resize(cruces_dobles);

  uniform_int_distribution<int> random_mutaciones(0, genes_por_generacion-1);

  for(int i=0; i < TAM_PBL; i++){
    poblacion[i].w.resize(w.size());
    for(int j=0; j < w.size(); j++){
      poblacion[i].w[j] = uniform_GEN(generator_GEN);
    }
    results = KNN_LOO(train, train, poblacion[i].w);
    poblacion[i].pts = agregado(results.tasa_clas, results.tasa_red);
    if(poblacion[i].pts > poblacion[mejor].pts)
      mejor = i;
    t++;
  }

  while(t < MAX_ITER){

    // Seleccionamos los padres con el torneo binario
    for (int i=0; i < cruces_dobles; i++){
      select = binaryTournament(poblacion, generator_GEN);
      padres[i] = poblacion[select];
    }

    // Hacemos el cruce Aritmetico para obtener los hijos
    for (int i=0; i < num_cruces*2; i+=2) {
      hijoCA = cruceArit(padres[i], padres[i+1]);
      poblacion_intermedia[i/2] = hijoCA;
    }

    // Hacemos el número de mutaciones esperadas
    for(int i=0; i < num_mutaciones; i++){
      valor_mut = random_mutaciones(generator_GEN);
      cromosoma_mutar = valor_mut / w.size();
      gen_mutar = valor_mut % w.size();
      aux_mut = poblacion_intermedia[cromosoma_mutar].w[gen_mutar] + normal_GEN(generator_GEN);
      if (aux_mut > 1)
        aux_mut = 1.0;
      else if (aux_mut < 0)
        aux_mut = 0.0;
      poblacion_intermedia[cromosoma_mutar].w[gen_mutar] = aux_mut;
    }

    // Introducimos los últimos padres en la población
    for (int i=num_cruces*2; i < cruces_dobles; i++)
      poblacion_intermedia[i-num_cruces] = padres[i];

    // Evaluamos los nuevos hijos
    for (int i=0; i < num_cruces; i++) {
      results = KNN_LOO(train, train, poblacion_intermedia[i].w);
      poblacion_intermedia[i].pts = agregado(results.tasa_clas, results.tasa_red);
      t++;
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

void AGE_BLX(const vector<FicheroCSV*>& train, vector<double>& w) {
  vector<Cromosoma> poblacion(TAM_PBL), poblacion_intermedia(2);
  pair<Cromosoma, Cromosoma> hijosBLX;
  Cromosoma p1, p2;
  int peor = 0, mejorNueva = 0, peorNueva = 0, select = 0, num_cruces, num_mutaciones;
  int genes_por_generacion, t = 0;
  int cromosoma_mutar, gen_mutar, valor_mut;
  double aux_mut;
  Resultados results = {0.0};

  num_cruces = floor(PB_CRUCE * (TAM_PBL/2))*2;
  genes_por_generacion = 2*w.size();
  num_mutaciones = ceil(PB_MUT * genes_por_generacion);

  uniform_int_distribution<int> random_mutaciones(0, genes_por_generacion-1);

  for(int i=0; i < TAM_PBL; i++){
    poblacion[i].w.resize(w.size());
    for(int j=0; j < w.size(); j++){
      poblacion[i].w[j] = uniform_GEN(generator_GEN);
    }
    results = KNN_LOO(train, train, poblacion[i].w);
    poblacion[i].pts = agregado(results.tasa_clas, results.tasa_red);
    if(poblacion[i].pts < poblacion[peor].pts)
      peor = i;
    t++;
  }

  while(t < MAX_ITER){

    // Seleccionamos los padres con el torneo binario
    select = binaryTournament(poblacion, generator_GEN);
    p1 = poblacion[select];
    select = binaryTournament(poblacion, generator_GEN);
    p2 = poblacion[select];

    // Hacemos el cruce BLX para obtener los hijos
    hijosBLX = cruceBLX(p1, p2, generator_GEN);
    poblacion_intermedia[0] = hijosBLX.first;
    poblacion_intermedia[1] = hijosBLX.second;

    // Hacemos el número de mutaciones esperadas
    for(int i=0; i < num_mutaciones; i++){
      valor_mut = random_mutaciones(generator_GEN);
      cromosoma_mutar = valor_mut / w.size();
      gen_mutar = valor_mut % w.size();
      aux_mut = poblacion_intermedia[cromosoma_mutar].w[gen_mutar] + normal_GEN(generator_GEN);
      if (aux_mut > 1)
        aux_mut = 1.0;
      else if (aux_mut < 0)
        aux_mut = 0.0;
      poblacion_intermedia[cromosoma_mutar].w[gen_mutar] = aux_mut;
    }

    // Evaluamos los nuevos hijos
    for(int i=0; i < 2; i++){
      results = KNN_LOO(train, train, poblacion_intermedia[i].w);
      poblacion_intermedia[i].pts = agregado(results.tasa_clas, results.tasa_red);
      t++;
    }

    if(poblacion_intermedia[0].pts > poblacion_intermedia[1].pts){
      mejorNueva = 0;
      peorNueva = 1;
    }
    else{
      mejorNueva = 1;
      peorNueva = 0;
    }

    if(poblacion_intermedia[mejorNueva].pts > poblacion[peor].pts){
      poblacion[peor] = poblacion_intermedia[mejorNueva];
      peor = getWorstCromosoma(poblacion);
      if(poblacion_intermedia[peorNueva].pts > poblacion[peor].pts){
        poblacion[peor] = poblacion_intermedia[peorNueva];
        peor = getWorstCromosoma(poblacion);
      }
    }

  }

  w = poblacion[getBestCromosoma(poblacion)].w;
}

void AGE_CA(const vector<FicheroCSV*>& train, vector<double>& w) {
  vector<Cromosoma> poblacion(TAM_PBL), poblacion_intermedia(2);
  Cromosoma p1, p2;
  int peor = 0, mejorNueva = 0, peorNueva = 0, select = 0, num_cruces, num_mutaciones;
  int genes_por_generacion, t = 0;
  int cromosoma_mutar, gen_mutar, valor_mut;
  double aux_mut;
  Resultados results = {0.0};

  num_cruces = floor(PB_CRUCE * (TAM_PBL/2))*2;
  genes_por_generacion = 2*w.size();
  num_mutaciones = ceil(PB_MUT * genes_por_generacion);

  uniform_int_distribution<int> random_mutaciones(0, genes_por_generacion-1);

  for(int i=0; i < TAM_PBL; i++){
    poblacion[i].w.resize(w.size());
    for(int j=0; j < w.size(); j++){
      poblacion[i].w[j] = uniform_GEN(generator_GEN);
    }
    results = KNN_LOO(train, train, poblacion[i].w);
    poblacion[i].pts = agregado(results.tasa_clas, results.tasa_red);
    if(poblacion[i].pts < poblacion[peor].pts)
      peor = i;
    t++;
  }

  while(t < MAX_ITER){

    // Seleccionamos los padres con el torneo binario
    select = binaryTournament(poblacion, generator_GEN);
    p1 = poblacion[select];
    select = binaryTournament(poblacion, generator_GEN);
    p2 = poblacion[select];

    // Hacemos el cruce Aritmetico para obtener los hijos
    poblacion_intermedia[0] = cruceArit(p1, p2);

    select = binaryTournament(poblacion, generator_GEN);
    p1 = poblacion[select];
    select = binaryTournament(poblacion, generator_GEN);
    p2 = poblacion[select];

    poblacion_intermedia[1] = cruceArit(p1, p2);

    // Hacemos el número de mutaciones esperadas
    for(int i=0; i < num_mutaciones; i++){
      valor_mut = random_mutaciones(generator_GEN);
      cromosoma_mutar = valor_mut / w.size();
      gen_mutar = valor_mut % w.size();
      aux_mut = poblacion_intermedia[cromosoma_mutar].w[gen_mutar] + normal_GEN(generator_GEN);
      if (aux_mut > 1)
        aux_mut = 1.0;
      else if (aux_mut < 0)
        aux_mut = 0.0;
      poblacion_intermedia[cromosoma_mutar].w[gen_mutar] = aux_mut;
    }

    // Evaluamos los nuevos hijos
    for(int i=0; i < 2; i++){
      results = KNN_LOO(train, train, poblacion_intermedia[i].w);
      poblacion_intermedia[i].pts = agregado(results.tasa_clas, results.tasa_red);
      t++;
    }

    if(poblacion_intermedia[0].pts > poblacion_intermedia[1].pts){
      mejorNueva = 0;
      peorNueva = 1;
    }
    else{
      mejorNueva = 1;
      peorNueva = 0;
    }

    if(poblacion_intermedia[mejorNueva].pts > poblacion[peor].pts){
      poblacion[peor] = poblacion_intermedia[mejorNueva];
      peor = getWorstCromosoma(poblacion);
      if(poblacion_intermedia[peorNueva].pts > poblacion[peor].pts){
        poblacion[peor] = poblacion_intermedia[peorNueva];
        peor = getWorstCromosoma(poblacion);
      }
    }

  }

  w = poblacion[getBestCromosoma(poblacion)].w;
}

Resultados ejecutarAGG_BLX(vector<vector<FicheroCSV*>>& particion, int num_part, const int& SEED){
	generator_GEN = default_random_engine(SEED);

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
  AGG_BLX(train, w);
  results = KNN(train, test, w);
  results.tiempo = elapsed_time(REAL);

  cout << "Partición " << num_part+1 << " %_clas: " << results.tasa_clas << " %_red: " << results.tasa_red
  << " Agr: " << agregado(results.tasa_clas, results.tasa_red) << " T: " << results.tiempo << endl;

  return results;
}

Resultados ejecutarAGE_BLX(vector<vector<FicheroCSV*>>& particion, int num_part, const int& SEED){
	generator_GEN = default_random_engine(SEED);

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
  AGE_BLX(train, w);
  results = KNN(train, test, w);
  results.tiempo = elapsed_time(REAL);

  cout << "Partición " << num_part+1 << " %_clas: " << results.tasa_clas << " %_red: " << results.tasa_red
  << " Agr: " << agregado(results.tasa_clas, results.tasa_red) << " T: " << results.tiempo << endl;

  return results;
}

Resultados ejecutarAGG_CA(vector<vector<FicheroCSV*>>& particion, int num_part, const int& SEED){
	generator_GEN = default_random_engine(SEED);

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
  AGG_CA(train, w);
  results = KNN(train, test, w);
  results.tiempo = elapsed_time(REAL);

  cout << "Partición " << num_part+1 << " %_clas: " << results.tasa_clas << " %_red: " << results.tasa_red
  << " Agr: " << agregado(results.tasa_clas, results.tasa_red) << " T: " << results.tiempo << endl;

  return results;
}

Resultados ejecutarAGE_CA(vector<vector<FicheroCSV*>>& particion, int num_part, const int& SEED){
	generator_GEN = default_random_engine(SEED);

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
  AGE_CA(train, w);
  results = KNN(train, test, w);
  results.tiempo = elapsed_time(REAL);

  cout << "Partición " << num_part+1 << " %_clas: " << results.tasa_clas << " %_red: " << results.tasa_red
  << " Agr: " << agregado(results.tasa_clas, results.tasa_red) << " T: " << results.tiempo << endl;

  return results;
}
