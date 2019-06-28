#include "Evolutivos.h"

pair<Cromosoma, Cromosoma> cruceBLX(const Cromosoma& p1, const Cromosoma& p2, default_random_engine& gen){
  Cromosoma h1, h2;
  h1.w.resize(p1.w.size());
  h2.w.resize(p1.w.size());
  double cmax, cmin, diff, random;

  for (int i=0; i < p1.w.size(); i++){

    cmax = max(p1.w[i], p2.w[i]);
    cmin = min(p1.w[i], p2.w[i]);

    diff = cmax - cmin;

    uniform_real_distribution<double> distribucion(cmin-diff*ALPHA_AGG, cmax+diff*ALPHA_AGG);

    random = distribucion(gen);
    if (random > 1)
      random = 1;
    else if (random < 0)
      random = 0;
    h1.w[i] = random;

    random = distribucion(gen);
    if (random > 1)
      random = 1.0;
    else if (random < 0)
      random = 0.0;
    h2.w[i] = random;
  }

  return make_pair(h1, h2);
}

Cromosoma cruceArit(const Cromosoma &p1, const Cromosoma &p2){
  Cromosoma h;
  h.w.resize(p1.w.size());

  for (int i=0; i < p1.w.size(); i++)
    h.w[i] = (p1.w[i]+p2.w[i]) / 2.0;

  return h;
}

int binaryTournament(const vector<Cromosoma>& poblacion, default_random_engine& gen){
  int num1, num2;
  int max = 0;
  uniform_int_distribution<int> random_int(0, poblacion.size()-1);

  num1 = random_int(gen);
  do {
    num2 = random_int(gen);
  } while (num1 == num2);

  if (poblacion[num1].pts > poblacion[num2].pts)
    max = num1;
  else
    max = num2;

  return max;
}

vector<int> getListBestCromosoma(const vector<Cromosoma>& poblacion){
  vector<int> index(poblacion.size());
  for(int i=0; i < poblacion.size(); i++)
    index[i] = i;

  sort(index.begin(), index.end(), [&](int i1, int i2) {
    return poblacion[i1].pts > poblacion[i2].pts; });

  return index;
}

int getBestCromosoma(const vector<Cromosoma>& poblacion){
  int max = 0;
  double pts_max = 0.0;
  for(int i=0; i < poblacion.size(); i++){
    if(poblacion[i].pts > pts_max){
      max = i;
      pts_max = poblacion[i].pts;
    }
  }
  return max;
}

int getWorstCromosoma(const vector<Cromosoma>& poblacion){
  int min = 0;
  double pts_min = 100.0;
  for(int i=0; i < poblacion.size(); i++){
    if(poblacion[i].pts < pts_min){
      min = i;
      pts_min = poblacion[i].pts;
    }
  }
  return min;
}
