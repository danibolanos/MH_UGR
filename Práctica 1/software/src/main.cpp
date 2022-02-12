#include "Relief.h"
#include "BL.h"

int main(int argc, char *argv[]) {
	int seed;
	double tam = 5.0;

	if(argc==3){
		vector<FicheroCSV> v;
		vector<vector<FicheroCSV*>> parts;
		int casos = 0;
		int categorias = 0;
		Resultados resultados;
		Resultados media = {0.0};
		seed = stoi(argv[2]);

		read_csv(argv[1], v);
		categorias = makePartitions(v, parts);
		cout << "/-*****************************************************************-/" << endl;
		cout << "\t\t\t         1-NN         " << endl;
		cout << "/-*****************************************************************-/" << endl;
		for(int i=0; i < parts.size(); i++){
			resultados = ejecutarKNN(parts, i);
			media.tasa_clas += resultados.tasa_clas;
			media.tasa_red += resultados.tasa_red;
			media.tiempo += resultados.tiempo;
		}
		cout << "/-------------------------------------------------------------------/" << endl;
		cout << "Media:" << " %_clas: " << media.tasa_clas/tam << " %_red: " << media.tasa_red/tam
		<< " Agr: " << agregado(media.tasa_clas/tam, media.tasa_red/tam)
		<< " T: " << media.tiempo/tam << endl << endl;

		media = {0.0};

		cout << "/-*****************************************************************-/" << endl;
		cout << "\t\t\t    Greedy-Relief    " << endl;
		cout << "/-*****************************************************************-/" << endl;
		for(int i=0; i < parts.size(); i++){
			resultados = ejecutarRelief(parts, i);
			media.tasa_clas += resultados.tasa_clas;
			media.tasa_red += resultados.tasa_red;
			media.tiempo += resultados.tiempo;
		}
		cout << "/-------------------------------------------------------------------/" << endl;
		cout << "Media:" << " %_clas: " << media.tasa_clas/tam << " %_red: " << media.tasa_red/tam
		<< " Agr: " << agregado(media.tasa_clas/tam, media.tasa_red/tam)
		<< " T: " << media.tiempo/tam << endl << endl;

		media = {0.0};

		cout << "/-*****************************************************************-/" << endl;
		cout << "\t\t\t   Búsqueda Local   " << endl;
		cout << "/-*****************************************************************-/" << endl;
		for(int i=0; i < parts.size(); i++){
			resultados = ejecutarBL(parts, i, seed);
			media.tasa_clas += resultados.tasa_clas;
			media.tasa_red += resultados.tasa_red;
			media.tiempo += resultados.tiempo;
		}
		cout << "/-------------------------------------------------------------------/" << endl;
		cout << "Media:" << " %_clas: " << media.tasa_clas/tam << " %_red: " << media.tasa_red/tam
		<< " Agr: " << agregado(media.tasa_clas/tam, media.tasa_red/tam)
		<< " T: " << media.tiempo/tam << endl << endl;

		for(int i=0; i < parts.size(); ++i)
			casos += parts[i].size();

		cout << "Había " << casos << " casos divididos en " << categorias << " categorias." << endl;
	}

	else if(argc==2){
		vector<FicheroCSV> colposcopy, ionosphere, texture;
		vector<vector<FicheroCSV*>> partsC, partsI, partsT;
		double tiempo = 0.0;
		Resultados resultados;
		seed = stoi(argv[1]);
		string tablas = "./data/tablas_" + to_string(seed) + ".csv";
		Resultados mediaC = {0.0}, mediaI = {0.0}, mediaT = {0.0};

		cout << "Creando tablas...(El proceso puede tardar varios minutos)" << endl;

		ofstream myfile;
		myfile.open(tablas);
		myfile << " ,COLPOSCOPY,,,,IONOSPHERE,,,,TEXTURE,,,,\n";
		myfile << ",1-NN,,,,1-NN,,,,1-NN\n";
		myfile << " , %tasa_clas,%tasa_red,Agr.,T(seg),%tasa_clas,%tasa_red,Agr.,T(seg),%tasa_clas,%tasa_red,Agr.,T(seg)\n";

		read_csv("./data/colposcopy.csv", colposcopy);
		makePartitions(colposcopy, partsC);
		read_csv("./data/ionosphere.csv", ionosphere);
		makePartitions(ionosphere, partsI);
		read_csv("./data/texture.csv", texture);
		makePartitions(texture, partsT);

		for(int i=0; i < partsC.size(); i++){
				resultados = ejecutarKNN(partsC, i);
				myfile << "Partición " << i+1 << "," << resultados.tasa_clas << "," << resultados.tasa_red << "," <<
				agregado(resultados.tasa_clas, resultados.tasa_red) << "," << resultados.tiempo;
				mediaC.tasa_clas += resultados.tasa_clas;
				mediaC.tasa_red += resultados.tasa_red;
				mediaC.tiempo += resultados.tiempo;
				resultados = ejecutarKNN(partsI, i);
				myfile << "," << resultados.tasa_clas << "," << resultados.tasa_red << "," <<
				agregado(resultados.tasa_clas, resultados.tasa_red) << "," << resultados.tiempo;
				mediaI.tasa_clas += resultados.tasa_clas;
				mediaI.tasa_red += resultados.tasa_red;
				mediaI.tiempo += resultados.tiempo;
				resultados = ejecutarKNN(partsT, i);
				myfile << "," << resultados.tasa_clas << "," << resultados.tasa_red << "," <<
				agregado(resultados.tasa_clas, resultados.tasa_red) << "," << resultados.tiempo << endl;
				mediaT.tasa_clas += resultados.tasa_clas;
				mediaT.tasa_red += resultados.tasa_red;
				mediaT.tiempo += resultados.tiempo;
		}

		myfile << "Media:" << "," << mediaC.tasa_clas/tam << "," << mediaC.tasa_red/tam
		<< "," << agregado(mediaC.tasa_clas/tam, mediaC.tasa_red/tam)
		<< "," << mediaC.tiempo/tam;
		myfile << "," << mediaI.tasa_clas/tam << "," << mediaI.tasa_red/tam
		<< "," << agregado(mediaI.tasa_clas/tam, mediaI.tasa_red/tam)
		<< "," << mediaI.tiempo/tam;
		myfile << "," << mediaT.tasa_clas/tam << "," << mediaT.tasa_red/tam
		<< "," << agregado(mediaT.tasa_clas/tam, mediaT.tasa_red/tam)
		<< "," << mediaT.tiempo/tam << endl;

		mediaC = mediaI = mediaT = {0.0};

		myfile << ",Greedy-Relief,,,,Greedy-Relief,,,,Greedy-Relief\n";

		for(int i=0; i < partsC.size(); i++){
				resultados = ejecutarRelief(partsC, i);
				myfile << "Partición " << i+1 << "," << resultados.tasa_clas << "," << resultados.tasa_red << "," <<
				agregado(resultados.tasa_clas, resultados.tasa_red) << "," << resultados.tiempo;
				mediaC.tasa_clas += resultados.tasa_clas;
				mediaC.tasa_red += resultados.tasa_red;
				mediaC.tiempo += resultados.tiempo;
				resultados = ejecutarRelief(partsI, i);
				myfile << "," << resultados.tasa_clas << "," << resultados.tasa_red << "," <<
				agregado(resultados.tasa_clas, resultados.tasa_red) << "," << resultados.tiempo;
				mediaI.tasa_clas += resultados.tasa_clas;
				mediaI.tasa_red += resultados.tasa_red;
				mediaI.tiempo += resultados.tiempo;
				resultados = ejecutarRelief(partsT, i);
				myfile << "," << resultados.tasa_clas << "," << resultados.tasa_red << "," <<
				agregado(resultados.tasa_clas, resultados.tasa_red) << "," << resultados.tiempo << endl;
				mediaT.tasa_clas += resultados.tasa_clas;
				mediaT.tasa_red += resultados.tasa_red;
				mediaT.tiempo += resultados.tiempo;
		}

		myfile << "Media:" << "," << mediaC.tasa_clas/tam << "," << mediaC.tasa_red/tam
		<< "," << agregado(mediaC.tasa_clas/tam, mediaC.tasa_red/tam)
		<< "," << mediaC.tiempo/tam;
		myfile << "," << mediaI.tasa_clas/tam << "," << mediaI.tasa_red/tam
		<< "," << agregado(mediaI.tasa_clas/tam, mediaI.tasa_red/tam)
		<< "," << mediaI.tiempo/tam;
		myfile << "," << mediaT.tasa_clas/tam << "," << mediaT.tasa_red/tam
		<< "," << agregado(mediaT.tasa_clas/tam, mediaT.tasa_red/tam)
		<< "," << mediaT.tiempo/tam << endl;

		mediaC = mediaI = mediaT = {0.0};

		myfile << ",BL,,,,BL,,,,BL\n";

		for(int i=0; i < partsC.size(); i++){
				resultados = ejecutarBL(partsC, i, seed);
				myfile << "Partición " << i+1 << "," << resultados.tasa_clas << "," << resultados.tasa_red << "," <<
				agregado(resultados.tasa_clas, resultados.tasa_red) << "," << resultados.tiempo;
				mediaC.tasa_clas += resultados.tasa_clas;
				mediaC.tasa_red += resultados.tasa_red;
				mediaC.tiempo += resultados.tiempo;
				resultados = ejecutarBL(partsI, i, seed);
				myfile << "," << resultados.tasa_clas << "," << resultados.tasa_red << "," <<
				agregado(resultados.tasa_clas, resultados.tasa_red) << "," << resultados.tiempo;
				mediaI.tasa_clas += resultados.tasa_clas;
				mediaI.tasa_red += resultados.tasa_red;
				mediaI.tiempo += resultados.tiempo;
				resultados = ejecutarBL(partsT, i, seed);
				myfile << "," << resultados.tasa_clas << "," << resultados.tasa_red << "," <<
				agregado(resultados.tasa_clas, resultados.tasa_red) << "," << resultados.tiempo << endl;
				mediaT.tasa_clas += resultados.tasa_clas;
				mediaT.tasa_red += resultados.tasa_red;
				mediaT.tiempo += resultados.tiempo;
		}

		myfile << "Media:" << "," << mediaC.tasa_clas/tam << "," << mediaC.tasa_red/tam
		<< "," << agregado(mediaC.tasa_clas/tam, mediaC.tasa_red/tam)
		<< "," << mediaC.tiempo/tam;
		myfile << "," << mediaI.tasa_clas/tam << "," << mediaI.tasa_red/tam
		<< "," << agregado(mediaI.tasa_clas/tam, mediaI.tasa_red/tam)
		<< "," << mediaI.tiempo/tam;
		myfile << "," << mediaT.tasa_clas/tam << "," << mediaT.tasa_red/tam
		<< "," << agregado(mediaT.tasa_clas/tam, mediaT.tasa_red/tam)
		<< "," << mediaT.tiempo/tam << endl;

		myfile.close();
		cout << "Proceso terminado." << endl;
	}

  return 0;
}
