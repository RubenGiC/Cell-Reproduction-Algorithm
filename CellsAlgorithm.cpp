//============================================================================
// Name        : Celulas.cpp
// Author      : Ruben Girela Castellón
// Version     : 0.0
// Copyright   : Ruben Girela Castellón
// Description : Cells Algotihm in C++, Ansi-style
// Created     : 13 jun. 2021
//============================================================================

#include <iostream>
#include <vector>
#include <random> // generación de valores aleatorios
#include <cmath> // calculos matematicos
#include <stdlib.h> //necesario para los includes de abajo
#include <algorithm>// para barajar el vector
#include <ctime>// libreria de time para el shuffle aleatorio de vectores

extern "C"{
	#include "cec17.h"
}

using namespace std;

//algoritmo Solis-Wet
void clip(vector<double> &sol, int lower, int upper) {
  for (auto &val : sol) {
    if (val < lower) {
      val = lower;
    }
    else if (val > upper) {
      val = upper;
    }
  }
}

void increm_bias(vector<double> &bias, vector<double> dif) {
  for (unsigned i = 0; i < bias.size(); i++) {
    bias[i] = 0.2*bias[i]+0.4*(dif[i]+bias[i]);
  }
}

void decrement_bias(vector<double> &bias, vector<double> dif) {
  for (unsigned i = 0; i < bias.size(); i++) {
    bias[i] = bias[i]-0.4*(dif[i]+bias[i]);
  }
}

/**
 * Aplica el Solis Wets
 *
 * @param  sol solucion a mejorar.
 * @param fitness fitness de la solución.
 */
template <class Random>
void soliswets(vector<double> &sol, double &fitness, double delta, int maxevals, int lower, int upper, Random &random, int &eval_actual, int max_global) {
  const size_t dim = sol.size();
  vector<double> bias (dim), dif (dim), newsol (dim);
  double newfit;
  size_t i;

  int evals = 0;
  int num_success = 0;
  int num_failed = 0;

  while (evals < maxevals && eval_actual < max_global) {
    std::uniform_real_distribution<double> distribution(0.0, delta);

    for (i = 0; i < dim; i++) {
      dif[i] = distribution(random);
      newsol[i] = sol[i] + dif[i] + bias[i];
    }

    clip(newsol, lower, upper);
    newfit = cec17_fitness(&newsol[0]);
    evals += 1;
    ++eval_actual;

    if (newfit < fitness) {
      sol = newsol;
      fitness = newfit;
      increm_bias(bias, dif);
      num_success += 1;
      num_failed = 0;
    }
    else if (evals < maxevals && eval_actual < max_global) {

      for (i = 0; i < dim; i++) {
        newsol[i] = sol[i] - dif[i] - bias[i];
      }

      clip(newsol, lower, upper);
      newfit = cec17_fitness(&newsol[0]);
      evals += 1;
      ++eval_actual;

      if (newfit < fitness) {
        sol = newsol;
        fitness = newfit;
        decrement_bias(bias, dif);
        num_success += 1;
        num_failed = 0;
      }
      else {
        for (i = 0; i < dim; i++) {
          bias[i] /= 2;
        }

        num_success = 0;
        num_failed += 1;
      }
    }

    if (num_success >= 5) {
      num_success = 0;
      delta *= 2;
    }
    else if (num_failed >= 3) {
      num_failed = 0;
      delta /= 2;
    }
  }

}
//-----------------------------------------------------------------------

//hace una explotación suave
vector<pair<double, vector<double>>> explotacion(vector<pair<double, vector<double>>> cells, int dim, int &eval, int max_eval){
	vector<pair<double,vector<double>>> new_cells=cells;
	vector<pair<double,vector<double>>> actual_cells=cells;
	int max_nodo=2, n_nodo=0, cell=0, nodo=0;
	double val, fitness;

	//mientras no supere el numero máximo de evaluaciones
	for(int i=0; i<(1000*dim) && eval < max_eval; ++i){

		
		//calculo el nuevo valor del nodo de esa celula
		val = new_cells[cell].second[nodo] * (rand() % 10 + 2);

		//lo modularizo entre -100 y 100
		if(val < 0){
			val = (-1)*fmod(abs(val),100.000000001);
		}else{
			val = fmod(val,100.000000001);
		}

		//y si el valor aleatorio es > a 5 cambia de signo el valor
		if((rand() % 10 + 1) > 5)
			val = val *(-1);

		//reemplazamos el valor
		new_cells[cell].second[nodo] = val;	

		//calculamos su fitness
		fitness = cec17_fitness(&new_cells[cell].second[0]);
		++eval;

		//cout << fitness << endl;

		//si es mejor que la que tiene actualmente actualiza
		if(fitness < actual_cells[cell].first){

			actual_cells[cell].first = fitness;
			actual_cells[cell].second[nodo] = new_cells[cell].second[nodo];
			new_cells[cell].first = fitness;

			//pasamos al siguiente nodo o nodo de la siguiente celula si hemos
			//recorrido todos los nodos de esa celula
			++nodo;
			n_nodo = 0;

			//si hemos recorrido todos los nodos de esa celula pasamos a la siguiente celula
			if(nodo >= dim){
				nodo = 0;
				++cell;
				//si hemos recorrido todas las celulas hace otra pasada 
				if(cell >= cells.size())
					cell = 0;
			}
		}else{
			//restauramos el valor
			new_cells[cell].second[nodo] = actual_cells[cell].second[nodo];

			++n_nodo;
			//si supera el numero maximo de cambios de cada nodo
			if(n_nodo>=max_nodo){
				//pasamos al siguiente nodo o nodo de la siguiente celula si hemos
				//recorrido todos los nodos de esa celula
				++nodo;
				n_nodo = 0;

				//si hemos recorrido todos los nodos de esa celula pasamos a la siguiente celula
				if(nodo >= dim){
					nodo = 0;
					++cell;
					//si hemos recorrido todas las celulas hace otra pasada 
					if(cell >= cells.size())
						cell = 0;
				}
			}
		}
	}

	return actual_cells;
}

//en la mejora solo he cambiado el numero de evaluaciones a evaluar el resto es igual
vector<pair<double, vector<double>>> explotacionMejorado(vector<pair<double, vector<double>>> cells, int dim, int &eval, int max_eval, int seed){
	vector<pair<double,vector<double>>> new_cells=cells;
	vector<pair<double,vector<double>>> actual_cells=cells;
	int max_nodo=2, n_nodo=0, cell=0, nodo=0;
	double val, fitness;

	uniform_real_distribution<> dis(-100.0, 100.0);

	mt19937 gen(seed);

	//mientras no supere el numero máximo de evaluaciones
	for(int i=0; i<(1000*dim) && eval < max_eval; ++i){

		//reemplazamos el valor por un valor aleatorio entre [-100, 100]
		new_cells[cell].second[nodo] = dis(gen);	

		//calculamos su fitness
		fitness = cec17_fitness(&new_cells[cell].second[0]);
		++eval;

		//cout << fitness << endl;

		//si es mejor que la que tiene actualmente actualiza
		if(fitness < actual_cells[cell].first){

			actual_cells[cell].first = fitness;
			actual_cells[cell].second[nodo] = new_cells[cell].second[nodo];
			new_cells[cell].first = fitness;

			//pasamos al siguiente nodo o nodo de la siguiente celula si hemos
			//recorrido todos los nodos de esa celula
			++nodo;
			n_nodo = 0;

			//si hemos recorrido todos los nodos de esa celula pasamos a la siguiente celula
			if(nodo >= dim){
				nodo = 0;
				++cell;
				//si hemos recorrido todas las celulas hace otra pasada 
				if(cell >= cells.size())
					cell = 0;
			}
		}else{
			//restauramos el valor
			new_cells[cell].second[nodo] = actual_cells[cell].second[nodo];

			++n_nodo;
			//si supera el numero maximo de cambios de cada nodo
			if(n_nodo>=max_nodo){
				//pasamos al siguiente nodo o nodo de la siguiente celula si hemos
				//recorrido todos los nodos de esa celula
				++nodo;
				n_nodo = 0;

				//si hemos recorrido todos los nodos de esa celula pasamos a la siguiente celula
				if(nodo >= dim){
					nodo = 0;
					++cell;
					//si hemos recorrido todas las celulas hace otra pasada 
					if(cell >= cells.size())
						cell = 0;
				}
			}
		}
	}

	return actual_cells;
}

vector<pair<double, vector<double>>> Meiosis(vector<pair<double, vector<double>>> cells, int &eval, int max_eval, vector<int> indice, int dim, int decrece, int n_cells, int seed){

	//para los valores aleatorios
	uniform_real_distribution<> dis(-100.0, 100.0);
	mt19937 gen(seed);

	//guarda el nuevo conjunto de celulas
	vector<pair<double, vector<double>>> new_cells;

	//vector de indices de las celulas para elegir las parejas aleatoriamente y sin repetir
	vector<int> select_cell;

	for(unsigned int i=0; i< cells.size(); ++i){
		select_cell.push_back(i);
	}

	//cada pareja se fragmenta en 4 partes
	vector<vector<double>> fragmentos;
	vector<double> frag;
	vector<double> frag2;
	vector<double> cell;

	//contador de celulas creadas por cada pareja
	int count=0;
	//calculamos el numero maximo de nodos por fragmento
	int max = ceil((float)dim/4);

	//creamos el vector de indices para seleccionar los fragmentos aleatorios
	vector<int> ind_frag(8);

	for(int i=0; i<8; ++i){
		ind_frag[i]=i;
	}

	random_shuffle(ind_frag.begin(), ind_frag.end());//barajo el vector

	random_shuffle(select_cell.begin(), select_cell.end());//barajo el vector

	int l = 0;

	//si no decrementa el conjunto de celulas, sino que aumenta
	if(decrece == 0){

		//recorre cada 2 celulas
		for(unsigned int i=0; i < select_cell.size() && eval < max_eval; ++i){

			++l;
			if(l>=select_cell.size()) l = 0;

			//y de cada pareja fragmentamos cada celula en 4 fragmentos, en total 8
			for(int e=0; e<dim; ++e){
				frag.push_back(cells[select_cell[i]].second[e]);
				frag2.push_back(cells[select_cell[l]].second[e]);
				++count;

				//si ya tiene un fragmento lo añade
				if(count == max || e+1==dim){
					//añado los fragmentos y limpio esos fragmentos
					fragmentos.push_back(frag);
					frag.clear();
					fragmentos.push_back(frag2);
					frag2.clear();
					//reseteo el contador
					count = 0;
				}
			}

			//una vez fragmentados elegimos 3 fragmentos aleatorios y que no sean del mismo padre los 3
			for(int e = 0; e < n_cells && eval < max_eval; ++e){

				//añadimos los 3 fragmentos aleatoriamente a la celula
				cell = fragmentos[ind_frag[1]];
				cell.insert(cell.end(), fragmentos[ind_frag[2]].begin(),fragmentos[ind_frag[2]].end());
				cell.insert(cell.end(), fragmentos[ind_frag[3]].begin(),fragmentos[ind_frag[3]].end());

				unsigned int size = cell.size();

				//el resto de nodos faltantes se generan aleatoriamente
				do{
					cell.push_back(dis(gen));
					++size;
				}while(size < dim);

				//añadimos la nueva celula con su fitness
				new_cells.push_back(pair<double,vector<double>>(cec17_fitness(&cell[0]),cell));
				++eval;
				cell.clear();

				//rebarajamos la lista de indices
				random_shuffle(ind_frag.begin(), ind_frag.end());//barajo el vector
			}
			//y limpiamos la lista de fragmentos
			for(unsigned int e = 0; e < fragmentos.size(); ++e)
				fragmentos[e].clear();
			fragmentos.clear();
		}
	}else{
		//recorre cada 2 celulas
		for(unsigned int i=0; i < (select_cell.size()/n_cells) && eval < max_eval; ++i){

			++l;

			if(l>=select_cell.size()) l = 0;

			//y de cada pareja fragmentamos cada celula en 4 fragmentos, en total 8
			for(int e=0; e<dim; ++e){
				frag.push_back(cells[select_cell[i]].second[e]);
				frag2.push_back(cells[select_cell[l]].second[e]);
				++count;

				if(count == max){
					//añado los fragmentos y limpio esos fragmentos
					fragmentos.push_back(frag);
					frag.clear();
					fragmentos.push_back(frag2);
					frag2.clear();
					//reseteo el contador
					count = 0;
				}
			}

			//una vez fragmentados elegimos 3 fragmentos aleatorios
			cell = fragmentos[ind_frag[1]];
			cell.insert(cell.end(), fragmentos[ind_frag[2]].begin(),fragmentos[ind_frag[2]].end());
			cell.insert(cell.end(), fragmentos[ind_frag[3]].begin(),fragmentos[ind_frag[3]].end());

			unsigned int size = cell.size();

			//el resto de nodos faltantes se generan aleatoriamente
			do{
				cell.push_back(dis(gen));
				++size;
			}while(size < dim);

			new_cells.push_back(pair<double,vector<double>>(cec17_fitness(&cell[0]),cell));
			++eval;
			cell.clear();

			random_shuffle(ind_frag.begin(), ind_frag.end());//barajo el vector

			//y limpiamos la lista de fragmentos
			for(unsigned int e = 0; e < fragmentos.size(); ++e)
				fragmentos[e].clear();
			fragmentos.clear();
		}
	}
	return new_cells;
}


vector<pair<double, vector<double>>> Mitosis(vector<pair<double, vector<double>>> cells, int &eval, int max_eval, vector<int> indice, int dim, int decrece, int seed){
	
	//para los valores aleatorios
	uniform_real_distribution<> dis(-100.0, 100.0);
	mt19937 gen(seed);

	vector<pair<double, vector<double>>> new_cells;
	vector<double> cell1(dim);
	vector<double> cell2(dim);

	//crece el numero de celulas
	if(decrece == 0){
		//recorro cada celula
		for(unsigned int i=0; i < cells.size() && eval < max_eval; ++i){
			//para cada celula se recorre en orden distinto
			random_shuffle(indice.begin(), indice.end());//barajo el vector
			//recorro cada nodo de la celula
			for(unsigned int e = 0; e < indice.size(); ++e){
				//copio la mitad de los datos de la celula en la nueva celula
				if(e < indice.size()/2){
					cell1[indice[e]] = cells[i].second[indice[e]];
					cell2[indice[e]] = dis(gen);
				//y el resto se generan aleatoriamente
				}else{
					cell1[indice[e]] = dis(gen);
					cell2[indice[e]] = cells[i].second[indice[e]];
				}
			}

			//evaluo el fitness y las guardo
			new_cells.push_back(pair<double,vector<double>>(cec17_fitness(&cell1[0]),cell1));
			new_cells.push_back(pair<double,vector<double>>(cec17_fitness(&cell2[0]),cell2));
			eval += 2;
		}

	//decrece el numero de celulas
	}else{
		//celulas que se seleccionaran aleatoriamente
		vector<int> select_cell;

		for(unsigned int i=0; i< cells.size(); ++i){
			select_cell.push_back(i);
		}

		random_shuffle(select_cell.begin(), select_cell.end());//barajo el vector

		//recorro cada celula
		for(unsigned int i=0; i < select_cell.size()/2 && eval < max_eval; ++i){
			//para cada celula se recorre en orden distinto
			random_shuffle(indice.begin(), indice.end());//barajo el vector
			//recorro cada nodo de la celula
			for(unsigned int e = 0; e < indice.size(); ++e){
				//copio la mitad de los datos de la celula en la nueva celula
				if(e < indice.size()/2){
					cell1[indice[e]] = cells[select_cell[i]].second[indice[e]];
				//y el resto se generan aleatoriamente
				}else{
					cell1[indice[e]] = dis(gen);
				}
			}
			//evaluo el fitness y las guardo
			new_cells.push_back(pair<double,vector<double>>(cec17_fitness(&cell1[0]),cell1));
			++eval;
		}
	}

	return new_cells;
}

double Cells(vector<vector<double>> celulas, int dim, int max_eval, int seed){
	
	double best_f=-1, f=-1;
	vector<double> queen_cell;
	vector<pair<double,vector<double>>> cells1;
	vector<pair<double,vector<double>>> cells2;

	int eval=0, grupo = -1;

	int max_size = 0;

	int decrece = 0, decrece2 = 0;
	int constant_meiosis = 4;

	vector<int> indice(dim);

	for(int i=0; i < dim; ++i)
		indice[i]=i;
	
	//creamos los 2 grupos de celulas
	for(unsigned int i=0; i < celulas.size(); i+=2){

		cells1.push_back(pair<double,vector<double>>(cec17_fitness(&celulas[i][0]),celulas[i]));
		cells2.push_back(pair<double,vector<double>>(cec17_fitness(&celulas[i][0]),celulas[i+1]));
		eval +=2;
	}

	while(eval < max_eval){

		//aplicamos la explotación
		cells1 = explotacion(cells1, dim, eval, max_eval);
		cells2 = explotacion(cells2, dim, eval, max_eval);

		//calculamos el tamaño maximo de los 2 grupos para ahorrar iteraciones inecesarias
		max_size=cells1.size();
		if(cells2.size() > max_size)
			max_size=cells2.size();

		//guardamos el que tenga mejor fitness
		for(unsigned int i=0; i<max_size; ++i){

			//para el grupo 1 de celulas
			if(i < cells1.size()){
				if(cells1[i].first < best_f || best_f==-1){
					best_f = cells1[i].first;
					queen_cell = cells1[i].second;
					grupo = 1;
				}
			}
			//para el grupo 2 de celulas
			if(i < cells2.size()){
				if(cells2[i].first < best_f || best_f==-1){
					best_f = cells2[i].first;
					queen_cell = cells2[i].second;
					grupo = 2;
				}
			}
		}
		//Propagación o reproducción
		//cells1 = Meiosis(cells1, eval, max_eval, indice, dim, decrece2, constant_meiosis, seed);
		cells1 = Mitosis(cells1, eval, max_eval, indice, dim, decrece, seed);
		cells2 = Meiosis(cells2, eval, max_eval, indice, dim, decrece2, constant_meiosis, seed);

		//controlo que la constante este entre 4 y 2
		if(decrece2==0){
			--constant_meiosis;
		}else{
			++constant_meiosis;
		}

		//controlo que los rangos en la Meiosis este entre [5 y 120]
		if(constant_meiosis == 1){
			decrece2 = 1;
			++constant_meiosis;
		}
		else if(constant_meiosis == 5){
			decrece2 = 0;
			--constant_meiosis;
		}

		//controlo que los rangos en la Mitosis este entre [10 y 160]
		if(cells1.size() >= 160)
			decrece = 1;
		else if(cells1.size() <= 10)
			decrece = 0;

		//guardamos el que tenga mejor fitness
		for(unsigned int i=0; i<max_size; ++i){

			if(i < cells1.size()){
				if(cells1[i].first < best_f || best_f==-1){
					best_f = cells1[i].first;
					queen_cell = cells1[i].second;
					grupo = 1;
				}
			}

			if(i < cells2.size()){
				if(cells2[i].first < best_f || best_f==-1){
					best_f = cells2[i].first;
					queen_cell = cells2[i].second;
					grupo = 2;
				}
			}
		}

	}

	return best_f;
}

//mejora de celulas
//en la mejora solo he cambiado el numero de evaluaciones a evaluar en la explotación el resto es igual
double CellsMejorado(vector<vector<double>> celulas, int dim, int max_eval, int seed){
	
	double best_f=-1, f=-1;
	vector<double> queen_cell;
	vector<pair<double,vector<double>>> cells1;

	int eval=0, decrece = 0, constant_meiosis = 4;

	vector<int> indice(dim);

	for(int i=0; i < dim; ++i)
		indice[i]=i;
	
	//creamos los 2 grupos de celulas
	for(unsigned int i=0; i < celulas.size(); ++i){

		cells1.push_back(pair<double,vector<double>>(cec17_fitness(&celulas[i][0]),celulas[i]));
		++eval;
	}

	while(eval < max_eval){

		//aplicamos la explotación
		cells1 = explotacionMejorado(cells1, dim, eval, max_eval, seed);

		//guardamos el que tenga mejor fitness
		for(unsigned int i=0; i<cells1.size(); ++i){

			//para el grupo 1 de celulas
			if(i < cells1.size()){
				if(cells1[i].first < best_f || best_f==-1){
					best_f = cells1[i].first;
					queen_cell = cells1[i].second;
				}
			}
		}
		//Propagación o reproducción
		cells1 = Meiosis(cells1, eval, max_eval, indice, dim, decrece, constant_meiosis, seed);

		//controlo que la constante este entre 4 y 2
		if(decrece==0){
			--constant_meiosis;
		}else{
			++constant_meiosis;
		}

		//controlo que los rangos en la Meiosis este entre [5 y 120]
		if(constant_meiosis == 1){
			decrece = 1;
			++constant_meiosis;
		}
		else if(constant_meiosis == 5){
			decrece = 0;
			--constant_meiosis;
		}

		//guardamos el que tenga mejor fitness
		for(unsigned int i=0; i<max_size; ++i){

			if(i < cells1.size()){
				if(cells1[i].first < best_f || best_f==-1){
					best_f = cells1[i].first;
					queen_cell = cells1[i].second;
				}
			}
		}

	}

	return best_f;
}

//hibridación entre el algoritmo de las celulas y Solis-Wet
double CellsWithSW(vector<vector<double>> celulas, int dim, int max_eval, int seed){
	
	double best_f=-1, f=-1;
	vector<double> queen_cell;
	vector<pair<double,vector<double>>> cells1;
	vector<pair<double,vector<double>>> cells2;

	int eval=0, grupo = -1;

	int max_size = 0;

	int decrece = 0, decrece2 = 0;
	int constant_meiosis = 4, maxtimes=5;

	vector<int> indice(dim);

	mt19937 gen(seed);

	for(int i=0; i < dim; ++i)
		indice[i]=i;
	
	//creamos los 2 grupos de celulas
	for(unsigned int i=0; i < celulas.size(); i+=2){

		cells1.push_back(pair<double,vector<double>>(cec17_fitness(&celulas[i][0]),celulas[i]));
		cells2.push_back(pair<double,vector<double>>(cec17_fitness(&celulas[i][0]),celulas[i+1]));
		eval +=2;
	}

	while(eval < max_eval){

		//calculamos el tamaño maximo de los 2 grupos para ahorrar iteraciones inecesarias
		max_size=cells1.size();
		if(cells2.size() > max_size)
			max_size=cells2.size();

		//guardamos el que tenga mejor fitness
		for(unsigned int i=0; i<max_size; ++i){

			//para el grupo 1 de celulas
			if(i < cells1.size()){
				//aplicamos la explotación
				soliswets(cells1[i].second, cells1[i].first, 0.2, 100000/maxtimes-1, -100, 100, gen, eval, max_eval);

				if(cells1[i].first < best_f || best_f==-1){
					best_f = cells1[i].first;
					queen_cell = cells1[i].second;
					grupo = 1;
				}
			}
			//para el grupo 2 de celulas
			if(i < cells2.size()){
				//aplicamos la explotacíón
				soliswets(cells2[i].second, cells2[i].first, 0.2, 100000/maxtimes-1, -100, 100, gen, eval, max_eval);
				if(cells2[i].first < best_f || best_f==-1){
					best_f = cells2[i].first;
					queen_cell = cells2[i].second;
					grupo = 2;
				}
			}
		}
		//Propagación o reproducción
		//cells1 = Meiosis(cells1, eval, max_eval, indice, dim, decrece2, constant_meiosis);
		cells1 = Mitosis(cells1, eval, max_eval, indice, dim, decrece, seed);
		cells2 = Meiosis(cells2, eval, max_eval, indice, dim, decrece2, constant_meiosis, seed);

		//controlo que la constante este entre 4 y 2
		if(decrece2==0){
			--constant_meiosis;
		}else{
			++constant_meiosis;
		}

		//controlo que los rangos en la Meiosis este entre [5 y 120]
		if(constant_meiosis == 1){
			decrece2 = 1;
			++constant_meiosis;
		}
		else if(constant_meiosis == 5){
			decrece2 = 0;
			--constant_meiosis;
		}

		//controlo que los rangos en la Mitosis este entre [10 y 160]
		if(cells1.size() >= 160)
			decrece = 1;
		else if(cells1.size() <= 10)
			decrece = 0;

		//guardamos el que tenga mejor fitness
		for(unsigned int i=0; i<max_size; ++i){

			if(i < cells1.size()){
				if(cells1[i].first < best_f || best_f==-1){
					best_f = cells1[i].first;
					queen_cell = cells1[i].second;
					grupo = 1;
				}
			}

			if(i < cells2.size()){
				if(cells2[i].first < best_f || best_f==-1){
					best_f = cells2[i].first;
					queen_cell = cells2[i].second;
					grupo = 2;
				}
			}
		}

	}

	return best_f;
}


int main(int argc, char** argv) {

	if(argc>1){

		string alg = argv[1];

		if(alg == "cells" or alg == "cellsSW" or alg == "cellsMejorado"){

			//int dim = 10;
			vector<int> seed(10);
			seed.push_back(42);
			seed.push_back(5);
			seed.push_back(100);
			seed.push_back(-1);
			seed.push_back(50);
			seed.push_back(1995);
			seed.push_back(19);
			seed.push_back(21);
			seed.push_back(80);
			seed.push_back(7);
			
			double fitness;
			double fitness_best=-1;

			vector<vector<double>> solutions;

			double f=0;

			clock_t end_part;
			float elapsed_part;
			clock_t start_part;

			if(alg == "cells")
				cout << "USANDO EL ALGORITMO CELLS" << endl;
			else if (alg == "cellsMejorado")
				cout << "USANDO EL ALGORITMO CELLS MEJORADO" << endl;
			else
				cout << "USANDO EL ALGORITMO CELLS CON SOLIS-WET" << endl;

			clock_t start = clock();
			//recorre las dimensiones 10 30 50
			for(int dim=10; dim <=50; dim+=20){

				vector<double> sol(dim);
				vector<double> best(dim);

				start_part = clock();
				//recorre cada problema
				for(int j=0; j < 30; ++j){

					//repite el problema 10 veces
					for(int k=0; k<10; ++k){

						uniform_real_distribution<> dis(-100.0, 100.0);

						mt19937 gen(seed[k]);

						cec17_init(alg.c_str(), j+1, dim);

						for(int e = 0; e < 10; ++e){
							for(int i=0; i < dim; ++i){
								sol[i] = dis(gen);

							}
							solutions.push_back(sol);

							
						}
						if(alg == "cells")
							f += Cells(solutions, dim, 10000*dim, seed[k]);
						else if(alg == "cellsMejorado")
							f += CellsMejorado(solutions, dim, 10000*dim, seed[k]);
						else
							f += CellsWithSW(solutions, dim, 10000*dim, seed[k]);

						for(unsigned int i =0; i < solutions.size(); ++i)
							solutions[i].clear();

						solutions.clear();
					}
					f = f/10;
					cout << "Media Best cell[" << j+1 << "] con dimension (" << dim << "): " << scientific << cec17_error(f) << endl;
					f =0;
				}
				end_part = clock();
				elapsed_part = float(end_part - start_part)/CLOCKS_PER_SEC;

				cout << "Elapse con dimensión (" << dim << "): " << (elapsed_part/60.0) << " (minutos)" << endl;

			}
			clock_t end = clock();
			float elapsed = float(end - start)/CLOCKS_PER_SEC;

			cout << "En total tarda: " << (elapsed/60.0) << " (minutos)" << endl;
		}else{
			cout << "Error solo usa los algoritmos (cells, cellsSW o cellsMejorado), escriba bien los algoritomos" << endl;
		}
	}else{
		cout << "Error necesita especificar el nombre del algotitmo ./cells <cells, cellsSW o cellsMejorado>" << endl;
	}

	return 0;
}


