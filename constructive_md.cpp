#include "constructive_md.h"
#include "bossa_stack.h"
#include "bossa_timer.h"
#include "bossa_random.h"
#include "rfw_stats.h"
#include "bossa_heap.h"
#include <stdlib.h>
#include <string.h>
// Gabriel - alteração 13/09/2008 - início 
#include <sstream>
using namespace std;
// Gabriel - alteração 13/09/2008 - fim
void PMConstructiveMD::fatal (const char *func, const char *msg) {
	fprintf (stderr, "PMConstructiveMD::%s: %s.\n", func, msg);
	exit(-1);
}



/***************************************************
 *
 * print the name of the current method into 'file'
 *
 ***************************************************/

void PMConstructiveMD::printMethod (FILE *file) {
	fprintf (file, method_name);
}


/***************************************************
 *
 * parse string that specifies a method; returns
 * true if the string is correctly parsed (in that
 * case, all relevant local variables are changed)
 *
 ***************************************************/

bool PMConstructiveMD::tryMethod (char *m, const char *name, int code, int np) {
	if (np == 0) {
		if (strcmp (m, name) == 0) {
			sprintf (method_name, m);
			method = code;
			return true;		
		}
	} else if (np == 1) {
		char format[256], buffer[256];
		sprintf (format, "%s:%%s", name);
		if (sscanf (m, format, buffer) == 1) {
			sprintf (method_name, m);
			method = code;
			param[method] = atof (buffer);
			return true;
		}
	}
	return false;
}


/************************************************
 *
 * define which method will be executed by 'run'
 * - m is a string with the name of the methods
 *   followed by its parameters (if any) --- the
 *   separator is a colon (:).
 *
 ************************************************/

void PMConstructiveMD::setMethod (char *m) {
	if (tryMethod (m, "rgreedy", RGREEDY, 1)) return;
	if (tryMethod (m, "rpg", RPG, 1)) return;
	if (tryMethod (m, "pgreedy", PGREEDY, 1)) return;
	if (tryMethod (m, "pworst",  PWORST, 1)) return;
	if (tryMethod (m, "preverse", PREVERSE, 0)) return;
	if (tryMethod (m, "mix", MIX, 0)) return;
	if (tryMethod (m, "sample", SAMPLE, 1)) return;
	if (tryMethod (m, "mst", MST, 0)) return;

	fprintf (stderr, "WARNING: '%s' is not a valid constructive heuristic.\n", m);
}



/***********************************
 * 
 * run method defined by setMethod
 *
 ***********************************/

void PMConstructiveMD::run (PMSolution *s) {
	int m;
	if (method != MIX) {m = method;}
	else {m = BossaRandom::getInteger(0,MIX-1);}

	// Gabriel - alteração 30/06/2008	- início 
	// AQUI DEVE SER REALIZADA A MINERAÇÃO
	if (construcaoInicial == 0)
	{
		ostringstream buffer("");		
		ostringstream arq_padroes("");
		//buffer << MD_PATH << "./fpmax_hnmp 1 1 " << MD_PATH << "elite.txt 16 " << "1" << " 99 " << MD_PATH << "padroes.txt";

		int i = supmin;
		fprintf (stdout, "supmin %d\n", i);
		int id_arq_tmp = i;
		BossaTimer mine(false); //Isabel - 10/09/27
		arq_padroes << "padroes-" << s->getM() << "-" << time (NULL) << ".txt";
		buffer << "./fpmax_hnmp "  // minerador de padroes frequentes
			   << seed << " "                                          // seed
			   << id_arq_tmp                                           // id de arq temporario
			   << " elite.txt "  									   // base a ser minerada
			   << tam_elite_simples << " "                             // tamanho da base
			   << i << " "                                             // suporte minimo
			   << npadroes << " "                                      // quantidade de padroes a minerar
			   << arq_padroes.str().c_str(); // arquivo onde serão armazenados os padroes
		mine.resume();
		system(buffer.str().c_str());
		mine.pause(); 					
		buffer.str("");
		fprintf (stdout, "mine: %.5f\n", mine.getTime());		

		FILE *padroes = fopen(arq_padroes.str().c_str(), "r");
		if (padroes != NULL)
		{
			char linha[90000];
			while (fgets(linha,90000,padroes))
			{
				istringstream facilities;
				int facility;
				set<int> wa_padroes;

				facilities.str(linha);
				while (facilities >> facility)
				{
					wa_padroes.insert(facility);
				}
				ti_padroes.push_back(wa_padroes);
				wa_padroes.clear();
			}
		}
		// Passar o resultado da mineracao para uma tabela interna
		ti_padroes_iter = ti_padroes.begin();
		construcaoInicial = 1;
		remove(arq_padroes.str().c_str());	
		remove("elite.txt");
	}
	else
	{
	// Caso hajam menos de 16 padroes, faça um uso circular dos mesmos
		ti_padroes_iter++;
		if (ti_padroes_iter == ti_padroes.end())
			ti_padroes_iter = ti_padroes.begin();
	}

	// Gabriel - alteração 30/06/2008	- fim
	double v = param[m];
	switch (m) {
		case MST     : mst      (s, false); break;
		case RPG     : rpg      (s, v, false); break;
		case RGREEDY : greedy   (s, false, v, NULL); break;
		case PGREEDY : greedy   (s, true, v, NULL); break;
		case PWORST  : pworst   (s, NULL, v, false); break;
		case PREVERSE: preverse (s, false); break;
		case SAMPLE  : sample   (s, v, false); break;
		default: fatal ("run", "invalid method");
	}
}




/***********************************
 *
 * add k new elements to solution s
 *
 ***********************************/

void PMConstructiveMD::addRandom (PMSolution *s, int k) {
	int i;
	int p = s->getP();    //number of facilities in the solution
	int m = s->getM();    //potential facilities
	int count = m-p;      //number of candidates
	if (k>count) fatal ("addRandom", "too many facilities");

	int *v = new int[m+1];                                  
	for (i=1; i<=count; i++) v[i]=s->getFacility(i+p); //list of candidates
	partialShuffle (v, 1, count, k);                   //pick first k at random
	for (i=1; i<=k; i++) s->add(v[i],true);            //add them to the solution
	delete [] v;
}



/*******************************************************
 *
 * Random plus greedy; add a fraction of the facilities
 * at random, complete with the greedy algorithm
 *
 *******************************************************/

void PMConstructiveMD::rpg (PMSolution *s, double rfrac, bool stats) {
	if (rfrac<0 || rfrac>1) rfrac = BossaRandom::getDouble();

	BossaTimer t(true);

	int r = (int) ceil (rfrac * (double)s->getInstance()->getP()); //number of random facilities

	s->reset();             //start with an empty solution
	addRandom (s, r);       //add r elements at random
	greedy(s,false,0,NULL); //complete greedily

	//print statistics if necessary
	if (stats) {
		fprintf (stdout, "rpg %.2f\n", s->calcCost());
		fprintf (stdout, "rpgtime %.2f\n", t.getTime());
	}
}



/************************************************************
 *
 * greedy: adds facilities to s one at a time
 *
 *   running time: O(nk), where k is the number of times
 *   cities 'change' hands during the execution of the 
 *   algorithm. One can find trivial upper and lower bounds
 *   for k: n <= k <= np (n is a lower bound if we start
 *   from an empty solution; if we start with something else
 *   we may have k < n).
 *
 ************************************************************/

void PMConstructiveMD::greedy (PMSolution *s, bool proportional, double alpha, int *initlist) {
	int i, f;
	int verbose = 0;
	int *candlist;
	double *savelist = NULL;

	bool random_alpha = (alpha<0 || alpha>1); // fatal ("greedy", "invalid value of alpha");

	//-----------------------
	// get basic information
	//-----------------------
	PMInstance *inst = s->getInstance();	
	int n = inst->getN(); //number of cities
	int m = inst->getM(); //number of potential facilities
	int p = inst->getP(); //number of facilities that must be in the solution
	if (s->getP() == p) return; //may we're done already

	//------------------------------
	// check the list of candidates
	//------------------------------
	bool list_provided = (initlist!=NULL);
	if (!list_provided) {
		candlist = new int [m+1];
		candlist[0] = 0;             //position 0 contains the size of the array
		for (i=s->getP()+1; i<=m; i++) {
			candlist[0]++;
			candlist[candlist[0]] = s->getFacility(i);
		}
	} else candlist = initlist;

	//---------------------------------
	// determine closeness information
    //---------------------------------
	double *save = new double [m+1]; //save[i] = how much we save if we add i to the solution
	savelist = new double [m+1];

	save[0] = -1; //we'll use this in the future
	for (i=1; i<=m; i++) save[i]=0;
	s->ensureConsistency();

	//-------------------------------------
	// calculate 'save' for all facilities
    // this takes O(nm) time
	//-------------------------------------
	for (i=1; i<=n; i++) {             //for every city i
		double d = s->getDistClosest(i); //distance to the closest facility
		for (int j=s->getP()+1; j<=m; j++) {
			f = s->getFacility(j);
			double profit = d - inst->getDist(i,f);
			if (profit>0) save[f] += profit;
		}
	}

	if (verbose>1) fprintf (stderr, "Initial amounts computed.\n");

	//----------------------------------
	// we now add facilities one by one
	//----------------------------------
	while (s->getP() < p) {
		if (verbose>1) fprintf (stderr, "Adding facility %d.\n", s->getP()+1);

		//--------------------------------------
		// determine best facility to add:
		//--------------------------------------
		int ncand = candlist[0];
		int bestj, bestf;
	
		//-------------------------------------------
		// proportional method: facilities that save
		// the most have higher probabilitis
		//-------------------------------------------
		for (int j=1; j<=ncand; j++) {savelist[j]=save[candlist[j]];} //list of candidates
		
		if (proportional) {
			double factor = alpha;
			if (factor == 0) {
				factor = 2 * BossaRandom::getDouble();
			}
			double min = RFWStats::min (savelist, 1, ncand);
			double max = RFWStats::max (savelist, 1, ncand);
			if (min == max) {
				bestj = BossaRandom::getInteger (1, ncand);
			} else {
				double x = BossaRandom::getDouble();
				for (int j=1; j<=ncand; j++) {savelist[j]-=min;}
				bestj = RFW::getWeightedIndex (savelist, 1, ncand, x, factor);
			}
		} 
		
		//--------------------------------------------
		// normal method: select at random among the
		// best candidates (fraction alpha)
		//-------------------------------------------
		else {
			if (random_alpha) alpha = BossaRandom::getDouble();
			int c = (int) ceil (alpha*candlist[0]);                  //number of candidates to consider
			RFW::bound (c, 1, candlist[0]);                          //avoiding rounding errors
			bestj = RFW::randomLargest(savelist, 1, candlist[0], c); //get one of those c elements
		}
		bestf = candlist[bestj];

		//remove best from the list
		candlist[bestj] = candlist[candlist[0]];
		candlist[0] --;

		//finally, add facility
		s->add (bestf, true); //closeness info is automatically updated

		//---------------------------------------------------------
		// make all vertices 'change hands' for the next iteration
		// (if there is a next iteration, of course) 		
		//---------------------------------------------------------
		if (s->getP() < p) {
			for (int i=1; i<=n; i++) {
				if (s->getClosest(i)==bestf) { 
					double d1 = s->getDistClosest(i);  //current closest
					double d2 = s->getDistClosest2(i); //previous closest (has to be)
					for (int j=s->getP()+1; j<=m; j++) { //update 'save' for facilities not in solution
						int f = s->getFacility(j);
						double dif = inst->getDist(i,f);
						if (dif<d2) {
							save[f] -= (d2-dif);             //we won't save w.r.t. previous closest
							if (dif<d1) save[f] += (d1-dif); //but we may still save w.r.t. current closest
						}
					}
				}
			}
		}
	}

	if (verbose>1) fprintf (stderr, "All iterations executed.\n");

	//----------
	// clean up
	//----------
	if (!list_provided) {
		delete [] candlist;
		candlist = NULL;
	}
	delete [] save;
	if (savelist!=NULL) delete [] savelist;

	if (verbose>1) fprintf (stderr, "Greedy algorithm executed successfully.\n\n");
}


/*****************************************************************
 *
 * Proportional worst: 
 * - starts with a random facility
 * - in each further iteration, chooses costumer with probability
 *   proportional to how much thed'd save if chosen; adds closest
 *   facility to customer
 *
 * s (output): solution (input/output);
 * candlist (input): list of allowed facilities (if NULL, all)
 *                   (candlist[0] must be the length of the list)
 * stats (input): print statistics?
 *
 *****************************************************************/

void PMConstructiveMD::pworst (PMSolution *s, int *candlist, double exponent, bool stats) {
	BossaTimer st;
	st.start();

	int bestf, f, i, j;
	PMInstance *inst = s->getInstance();
	int p = inst->getP();
	int n = inst->getN();
	int m = inst->getM();

	if (p==0) fatal ("pworst", "number of facilities cannot be zero");

	double *w  = new double [n+1]; //weights used to compute probabilities
	double *dc = new double [n+1]; //distance to closest candidate

	bool list_given = (candlist!=NULL); //list of candidates may be given as input
	if (!list_given) {                  //if it isn't, assume every facility is a candidate
		candlist = new int [m+1];  
		candlist[0] = m;
		for (i=1; i<=m; i++) candlist[i] = i;
	}

	if (candlist[0]<p) fatal ("pworst", "not enough candidates");
	
	//---------------------------------------------------------------
	// find distance from each user to the closest facility
	// THIS IS NOT ALWAYS RELEVANT --- IN MOST CASES, IT'S JUST ZERO
	//---------------------------------------------------------------
	if (n==m) {
		for (i=1; i<=n; i++) dc[i] = inst->getDist(i,i); //this deals with the trivial case
	} else {
		for (i=1; i<=n; i++) dc[i] = INFINITY; //this is the general case
	}

	for (i=1; i<=n; i++) {
		if (dc[i]==0) continue;
		for (j=1; j<=candlist[0]; j++) {
			f = candlist[j];
			double d = inst->getDist(i,f);
			if (d<dc[i]) dc[i] = d;
		}
	}

	//-----------------------------------
	// start by adding a random facility
	//-----------------------------------
	s->reset();
	bestf = candlist[BossaRandom::getInteger(1,candlist[0])];
	s->add (bestf, true);
	
	//--------------
	// initialize w
	//--------------
	for (i=1; i<=n; i++) {
		w[i] = inst->getDist(i,bestf) - dc[i];
	}

	if (exponent == 0) exponent = 2 * BossaRandom::getDouble();
	fprintf (stderr, "<%.02f> ", exponent);

	//----------------------------------
	// main loop: add facilities
	//----------------------------------
	while (s->getP() < p) {
		//select a user at random with relative probabilities given by w
		int u = RFW::getWeightedIndex (w, 1, n, BossaRandom::getDouble(), exponent);

		//find closest facility not in the solution
		bestf = 0;
		for (j=1; j<=candlist[0]; j++) {
			f = candlist[j];
			if (s->contains(f)) continue;
			if (inst->getDist(u,f) < inst->getDist(u,bestf)) bestf = f;
		}
		
		//add facility to the solution
		s->add(bestf, true);

		//update weights
		for (u=1; u<=n; u++) {
			w[u] = inst->getDist(u,s->getClosest(u)) - dc[u];
		}
	}

	delete [] w;
	delete [] dc;
	if (!list_given) delete [] candlist;

	if (stats) {
		double result = s->calcCost();
		fprintf (stdout, "srandom %.2f\n", result);
		fprintf (stdout, "srandomtime %.3f\n", st.getTime());
	}
}








/**********************************
 *
 * This thing is not general
 *
 **********************************/

void PMConstructiveMD::mst (PMSolution *s, bool stats) {
	BossaVertexHeap <double> *heap;
	PMInstance *inst = s->getInstance();
	int m = inst->getM();
	//int n = inst->getN();
	int p = inst->getP();
	int i;
	BossaTimer t;
	t.start();

	s->reset();
	int *pred = new int [m+1];
	double *predval = new double [m+1];
	bool *intree = new bool [m+1];

	for (i=1; i<=m; i++) {
		pred[i] = 0;
		predval[i] = 0.0;
		intree[i] = false;
	}

	heap = new BossaVertexHeap <double> (m);

	int r = BossaRandom::getInteger (1, m);
	heap->insert (r, 0);
	while (!heap->isEmpty()) {
		double v;
		heap->removeFirst (i, v);
		intree[i] = true;

		for (int j=1; j<=m; j++) {
			if (intree[j]) continue;
			v = inst->getFacDist(i,j) * (1 + 0.1 * BossaRandom::getDouble()); //BossaRandom::getInteger (1,2);
			if (heap->insert (j, v)) {
				//fprintf (stderr, "Inserted %d with value %f.\n", j, v);
				pred[j] = i;
				predval[j] = v;
			}
		}
	}

	//calc cost
	double cost = 0;
	for (i=1; i<=m; i++) {
		if (i==r) continue;
		cost += inst->getFacDist(i,pred[i]);
	}
	//fprintf (stderr, "Total cost is %f.\n", cost);

	//-----------------------
	// remove heaviest edges
	//-----------------------
	IntDouble *list = new IntDouble [m+1];
	for (i=1; i<=m; i++) {
		double v;
		if (i==r) v = -INFINITY;
		else v = predval[i]; //inst->getFacDist(i,pred[i]);
		list[i-1].id = i;
		list[i-1].value = v;
	}

	//fprintf (stderr, "Sorting...\n");
	qsort (list, m, sizeof (class IntDouble), IntDouble::compare);
	for (i=m-1; i>m-p; i--) {intree[list[i].id] = false;} //remove p-1 edges
	intree[r] = false; //root is also in a component by itself
	delete [] list;

	//------------------------------------
	// reset colors
	//------------------------------------
	int *color = new int[m+1];
	//int *root  = new int[p+1];
	int c = 0;
	for (i=1; i<=m; i++) {
		if (!intree[i]) {c++; color[i] = c;}
		else color[i] = 0;
	}

	//fprintf (stderr, "Number of colors: %d\n", c);

	//-----------------
	// find components
	//-----------------
	BossaStack <int> *stack = new BossaStack <int> (m);
	for (i=1; i<=m; i++) {
		if (color[i]>0) {continue;}
		int j = i;
		while (color[j] == 0) {stack->push(j); j = pred[j];} //find root
		while (!stack->isEmpty()) {color[stack->pop()] = color[j];} //assign its color to everyoen
	}
	delete stack;
	
	//-----------------------------------------
	// choose a random element from each color
	//-----------------------------------------

	int *count = new int [p+1]; //how many elements in each color

	for (i=1; i<=p; i++) {count[i] = 0;}
	//fprintf (stderr, "reset...");
	for (i=1; i<=m; i++) {count[color[i]]++;}
	//fprintf (stderr, "counted...\n");
	//for (i=1; i<=p; i++) {
	//	fprintf (stderr, "%d ", count[i]);
	//}
	//fprintf (stderr, "\n");
	
	for (i=1; i<=p; i++) {count[i] = BossaRandom::getInteger(1,count[i]);}
	for (i=1; i<=m; i++) {
		c = color[i];
		if (count[c]>0) {
			count[c]--;
			if (count[c]==0) s->add(i, true);
		}
	}
	delete [] count;
	
	//fprintf (stderr, "Solution cost is %f.\n", s->calcCost());

	//delete [] root;
	delete [] color;
	delete [] intree;
	delete [] pred;
	delete [] predval;
	delete heap;

	if (stats) {
		fprintf (stdout, "mst %.2f\n", s->calcCost());
		fprintf (stdout, "msttime %.3f\n", t.getTime());
	}
}



void PMConstructiveMD::preverse (PMSolution *s, bool stats) {
	BossaTimer t;
	t.start();

	PMInstance *inst = s->getInstance();
	int p = inst->getP(); //number of facilities that must remain open
	int n = inst->getN(); //total number of users
	int m = inst->getM();
	int i, j, f;
	double *w    = new double [m+1];
	double *loss = new double [m+1];

	//--------------------------------
	// add all facilities to solution
	//--------------------------------
	for (f=1; f<=m; f++) {s->add(f,false);}
	s->ensureConsistency();

	//----------------------------------
	// now remove facilities one by one
	//----------------------------------
	
	while (s->getP() > p ) {
		int cp = s->getP();

		for (f=1; f<=m; f++) loss[f] = 0;
		for (i=1; i<=n; i++) {
			int c1 = s->getClosest(i);
			int c2 = s->getClosest2(i);
			double diff = inst->getDist(i,c2) - inst->getDist(i,c1);
			loss[s->getPosition(c1)] += diff;
		}

		double max = RFWStats::max (loss, 1, cp); //largest loss
		double min = RFWStats::min (loss, 1, cp); //smallest lost

		int best;
		if (max == min) {
			best = BossaRandom::getInteger (1, cp);
		} else {
		    for (j=1; j<=cp; j++) {w[j] = max - loss[j];}
			best = RFW::getWeightedIndex (w, 1, cp, BossaRandom::getDouble());
		}
		int bestf = s->getFacility(best);

		s->remove(bestf, true);
	}


	s->ensureConsistency();
	
	if (stats) {
		fprintf (stdout, "preverse %.2f\n", s->calcCost());
		fprintf (stdout, "preversetime %.3f\n", t.getTime());
	}
	
	delete [] loss;
	delete [] w;
}











/************************************************************
 *
 * sample: greedy, but sampling only 
 *
 ************************************************************/

void PMConstructiveMD::sample (PMSolution *s, double fraction, bool stats) {
	BossaTimer t;
	t.start();
	int m, n, p, verbose, i;

	PMInstance *inst = s->getInstance();	
	verbose = 0;
	n = inst->getN(); //number of cities
	m = inst->getM(); //number of potential facilities
	p = inst->getP(); //number of facilities that must be in the solution

	int sample_size;
	if (fraction<0) fraction = BossaRandom::getDouble();
	if (fraction>1) {
		sample_size = (int) ceil ((log ((double)m/(double)p))/LN2);
		
		//double a = 1.0/(double)p;
		//double b = (log ((double)m) / LN2) / (double)m; 
		//fraction = (a<b) ? a : b;
	} else {
		sample_size = (int)ceil (fraction * (double)m);
	}
	if (sample_size <= 0) sample_size = 1;
	if (verbose > 1) {fprintf (stderr, "Fraction:%f, absolute %d.\n", fraction, sample_size);}

	//--------------------------
	// build list of candidates
	//--------------------------
	int *candlist = new int [m+1];
	candlist[0] = 0;             
	for (i=s->getP()+1; i<=m; i++) {
		candlist[0]++;
		candlist[candlist[0]] = s->getFacility(i);
	}

	//RFW::shuffle (candlist, candlist[0]);

	//---------------------------------
	// Running time: O(p m fraction n)
	//---------------------------------
	s->reset(); //start with an empty solution
	// Gabriel - alteração 22/07/2008	- início 
	// AQUI DEVE SER !APLICADA! A MINERAÇÃO!!!!
	std::set<int> wa_padroes = *ti_padroes_iter;
	set<int>::iterator wa_padroes_iter; 
	for (wa_padroes_iter = wa_padroes.begin(); wa_padroes_iter != wa_padroes.end(); wa_padroes_iter++)
	{
		s->add(*wa_padroes_iter,true);
	}
	//s->add(fac, true);
	// Gabriel - alteração 22/07/2008	- fim
	int count = 0;
	while (s->getP() < p) {
		count++;
		if (sample_size > candlist[0]) sample_size = candlist[0];

		//RFW::shuffle(candlist, candlist[0]);
		partialShuffle (candlist, 1, candlist[0], sample_size);

		//find best option in the sample
		int bestc = -1; /* WARNING, THIS CHANGED */
		if (sample_size==1) bestc = 1;
		else {
			double bestsave = -1;
			for (int c=1; c<=sample_size; c++) {
				int f = candlist[c];

				//how much do we save by adding f?
				double save = 0; //s->calcSave (f); //0;
				for (int i=1; i<=n; i++) {
					double diff = s->getDistClosest(i) - inst->getDist(i,f); //(inst->getDist(i,s->getClosest(i)) - inst->getDist(i,f)); 
					if (diff > 0) save += diff;	
				}	
				if (save > bestsave) {bestsave = save; bestc = c;}
			}
		}
		int bestf = candlist[bestc];
		s->add (bestf, true);

		//remove facility from list of candidates
		candlist[bestc] = candlist[candlist[0]];
		candlist[0]--;
	}

	s->ensureConsistency();
		
	delete [] candlist;

	if (stats) {
		fprintf (stdout, "sample %.2f\n", s->calcCost());
		fprintf (stdout, "samplef %.3f\n", fraction);
		fprintf (stdout, "sampletime %.3f\n", t.getTime());
	}
// Gabriel - alteração 07/09/2008 - início - contando tempo de construcao
//	fprintf (stdout, "sampletime %.3f\n", t.getTime());
// Gabriel - alteração 07/09/2008 - fim
	if (verbose>0) fprintf (stderr, "(%.2f)", t.getTime());
}
