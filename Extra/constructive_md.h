/****************************************************************
 *
 * Candidate list
 *
 ****************************************************************/

#ifndef CONSTRUCTIVENEW_H
#define CONSTRUCTIVENEW_H

#include <stdio.h>
#include "solution.h"
// Gabriel - alteração 30/06/2008	- início 
#include "elite_simples.h"
#include <vector>
#include <iostream>
#include <set>
// Gabriel - alteração 30/06/2008	- fim

class PMConstructiveMD {
	private:
		enum {RPG, RGREEDY, PGREEDY, PWORST, PREVERSE, MST, SAMPLE, MIX, NMETHODS};
		char method_name [256];  //name of the current method
		double param [NMETHODS]; //default parameters for all methdos
		int method;              //current method

		static void fatal (const char *func, const char *msg);
		bool tryMethod (char *m, const char *name, int code, int np);
		// Gabriel - alteração 30/06/2008	- início 
		int seed;
		int construcaoInicial;
		int tam_elite_simples;
		int npadroes;
		int supmin;
		std::vector< std::set<int> > ti_padroes;
		std::vector< std::set<int> >::iterator ti_padroes_iter;
		// Gabriel - alteração 30/06/2008	- fim
	public:
		PMConstructiveMD () {setMethod ("sample:2");construcaoInicial=0;} //constructor; default is random solution

		//--------------------------------------------------
		// routines dealing the 'current' or default method
		//--------------------------------------------------
		void setMethod (char *m);      //define current method
		void printMethod (FILE *file); //print name of current method
		// Gabriel - alteração 30/06/2008	- início 
		// Adicionado o parametro "es" no método run
		void run(PMSolution *s);       //run current method
// Gabriel - alteração 06/09/2008 - início 
		// void run (PMSolution *s, PMEliteSimple *es);       //run current method
		void setSeed (int s) {seed = s;}
		void setTamEliteSimples (int es) {tam_elite_simples = es;}
		void setNumeroPadroes (int np) {npadroes = np;}
		void setSuporteMinimo (int sup) {supmin = sup;}
// Gabriel - alteração 06/09/2008 - fim 
		// Gabriel - alteração 30/06/2008	- fim
		//-------------------------
		// constructive algorithms
		//-------------------------
		static void addRandom (PMSolution *s, int k);
		static void pworst    (PMSolution *s, int *candlist=NULL, double exponent=1, bool stats=false);
		static void greedy    (PMSolution *s, bool proportional=false, double alpha=0, int *candlist=NULL);
		static void rpg       (PMSolution *s, double rfrac, bool stats=false);
		static void preverse  (PMSolution *s, bool stats=false);
		static void mst       (PMSolution *s, bool stats=false);
		void sample    (PMSolution *s, double frac, bool stats=false);
};

#endif
