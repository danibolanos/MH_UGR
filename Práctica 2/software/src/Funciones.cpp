#include "Funciones.h"

void read_csv(string filename, vector<FicheroCSV>& result){
	ifstream file(filename);
	string line;
	vector<double> max, min;
	max.clear();
	min.clear();

	getline(file, line);

	while (getline(file, line)) {
		FicheroCSV e;
	  istringstream s(line);
	  string field;

	  while (getline(s, field, ',')) {
	   	if (s.peek() != '\n' && s.peek() != EOF)
	   		e.traits.push_back(stod(field));
	    }

	    e.category = field;
	    e.n = e.traits.size();
	    result.push_back(e);

		if(max.empty()){
			max = e.traits;
			min = e.traits;
		}
		else{
			for (int i=0; i < e.n; i++){
				if( max[i] < e.traits[i])
					max[i] = e.traits[i];
				if( min[i] > e.traits[i])
					min[i] = e.traits[i];
			}
		}
	}

	int atributos = result[0].n;
	for (int i=0; i < result.size(); i++)
		for (int j=0; j < atributos; j++)
			if (max[j] != min [j])
				result[i].traits[j] = (result[i].traits[j]-min[j])/(max[j]-min[j]);
}

int makePartitions(vector<FicheroCSV>& datos, vector<vector<FicheroCSV*>>& particion){
	 set<string> categorias;
   string clase;
	 int contador = 0;
   bool stop = false;

	 particion.resize(5);

	 for (int i=0; i < datos.size() && !stop; i++){
		 clase = "\0";
   	 for (int j=i; j < datos.size(); j++){
			 if (categorias.count(datos[j].category)==0 &&  clase == "\0"){
				 categorias.insert(datos[j].category);
   			 clase = datos[j].category;
   			 particion[contador].push_back(&datos[j]);
         contador = (++contador)%5;
   		 }
			 else if (clase == datos[j].category){
   		 	 particion[contador].push_back(&datos[j]);
         contador = (++contador)%5;
       }
   	}
    if(clase == "\0")
    	stop = true;
    }
	return categorias.size();
}

double euclideanDistance(const vector<double>& v1, const vector<double>& v2, const vector<double>& w){
	double dist=0.0;
	for (int i=0; i < v1.size(); i++)
		if (w[i] >= 0.2)
			dist += w[i]*(v2[i]-v1[i])*(v2[i]-v1[i]);

	return dist;
}

double euclideanDistance(const vector<double>& v1, const vector<double>& v2){
	double dist=0.0;
	for (int i=0; i < v1.size(); i++)
			dist += (v2[i]-v1[i])*(v2[i]-v1[i]);
	return dist;
}
