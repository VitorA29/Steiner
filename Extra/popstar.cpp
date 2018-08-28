/****************************************************************
 *
 * program: pmedian
 *          provides heuristic solutions to the pmedian problem
 * author: Renato Werneck (rwerneck@princeton.edu)
 * log: May 29, 2002 - file created.
 *
 ****************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "basics.h"
#include "bossa_timer.h"
#include "config.h"
#include "elite.h"
#include "instance.h"
#include "search.h"
#include "solution.h"
#include "bossa_random.h"
#include "bossa_heap.h"
#include "rfw_stats.h"
#include "search_tables.h"
#include "constructive.h"
#include "distance.h"
#include "matrix_instance.h"
#include "euclidean_instance.h"
#include "graph_instance.h"
#include "geo_instance.h"
#include "hybrid_instance.h"
// Gabriel - alteração 30/06/2008	- início 
#ifdef MINERACAO
#include "elite_simples.h"
#include "constructive_md.h"

static int g_seed;
static int g_tam_elite_simples;
static int g_npadroes;
static int g_supmin;

#endif
// Gabriel - alteração 30/06/2008	- fim

class Partition {
	private:
		bool *used;   //used[i] = is name i used?
		int *prev;    //prev[i] = element that precedes i in its group
		int *next;    //next[i] = element that succeeds i in its group
		int *group;   //group[i]: group i belongs to
		int g; //number of groups
		int n; //number of elements

	public:
		



};





/******************************************************
 *
 * auxiliary class; stores the best value obtained
 * by an iteractive algorithm --- with solution value,
 * iteration number and time
 *
 ******************************************************/

class PMRecorder {
	public: 
		int best_iteration;
		double best_time;
		double best_value;

		PMRecorder () {
			best_time = 0;
			best_iteration = 0;
			best_value = INFINITY;
		}

		bool report (int iteration, double time, double value) {
			if (value < best_value) {
				best_value = value;
				best_time = time;
				best_iteration = iteration;
				return true;
			}
			return false;
		}
};



int main (int argc, char *argv[]);

double grasp (
  PMSolution *s,   //output: will hold the best solution found
  PMConstructive *constructive, //input: used to run the appropriate constructive heuristic
  int it,          //input: number of iterations to be performed
  PMElite *elite,  //input/output: will hold best solutions found and be used for path-relinking (can be NULL)
  char *lsmethod,   //local search method to be used
  char *combmethod,
  bool stats,      //input: print statistics to stdout?
  bool partials    //input: print partial information to stdout?
);

// Isabel - 30/08/2010
double grasp_target (
  PMSolution *s,   //output: will hold the best solution found
  PMConstructive *constructive, //input: used to run the appropriate constructive heuristic    
  int it, 	   //input: number of iterations to be performed
  double target, 
  PMElite *elite,  //input/output: will hold best solutions found and be used for path-relinking (can be NULL)
  char *lsmethod,  //local search method to be used
  char *combmethod,
  int *flag,
  bool stats,      //input: print statistics to stdout?
  bool partials    //input: print partial information to stdout?
);


//double solve2 (PMSolution *s, int *best1, int *best2, int *assign, double *acost);
void vns (PMSolution *opt_s, int kmax, int rmax, bool rvns, PMElite *elite, bool stats = false);
void runLocalSearch (PMSolution *s, char *lsmethod);

int upcount = 0; //number of call to "update"

/*************************************************
 *
 * fatal: prints error message and exists program
 *
 *************************************************/

void fatal (const char *func, const char *msg) {
	fprintf (stderr, "%s: %s.\n", func, msg);
	exit (-1);
}


/***********************************************************
 *
 * given a candidate for insertion (in), finds a the best
 * candidate for removal (out), together with the variation
 * in solution cost (w) associated with the resulting
 * swap (w)
 *
 ***********************************************************/

void move (PMSolution *s, int in, int *out, double *w) {
	int i;
	int p = s->getP();
	int n = s->getN();
	PMInstance *inst = s->getInstance();
	if (s->contains(in)) fatal ("move", "invalid insertion facility");
	
	double *v = new double [n+1];

	//initialization
	*w = 0;  //change in the objective function value obtained 
		 //by the best interchange; always negative
	for (i=1; i<=p; i++) v[s->getFacility(i)] = 0;
	for (i=1; i<=n; i++) v[i] = 0;
	
	//best deletion
	for (i=1; i<=n; i++) {
		double din = inst->getDist (i, in);
		//double d1  = inst->getDist (i, s->getClosest(i));
		//double d2  = inst->getDist (i, s->getClosest2(i));
		double d1 = s->getDistClosest(i);
		double d2 = s->getDistClosest2(i);

		if (din < d1) {   //WARNING: THIS LINE WAS WRONG IN THE PAPER (originally, din < d2)
			*w += (din - d1); //in will be the new closest, we'll save something
		} else {
			//if we delete closest[i], we'll save d1, but we'll 
			//have to spend either din or d2
			double min = (din<d2) ? din : d2;
			v[s->getClosest(i)] += (min-d1); //we'll have to pay for the difference 
		}
	}	

	//find best facility to remove
	int best = s->getFacility(1);
	for (i=1; i<=p; i++) {
		if (v[s->getFacility(i)]<v[best]) best = s->getFacility(i);
	}

	//fprintf (stderr, " (o:%d i:%d w:%8.2f v:%8.2f)", *out, in, *w, v[best]);
	*out = best;
	*w += v[best];

	delete [] v;
}


/************************************************************
 *
 * Given a filename, creates an instances of the appropriate
 * type (based on the extension), reads the instance, and
 * returns a pointer to the allocated instance; the instance
 * must be deallocated elsewhere.
 *
 ************************************************************/

PMInstance *readInstance (const char *filename, int p) {
	PMInstance *instance = NULL;
	bool dimacs, tsp, pmi, pmm, geo;
	tsp = pmi = pmm = geo = false;
	dimacs = true;
	if (strlen(filename)>=4) {
		if (strcmp (&filename[strlen(filename)-4], ".tsp") == 0) tsp = true;
		else if (strcmp (&filename[strlen(filename)-4], ".pmm") == 0) pmm = true;
		else if (strcmp (&filename[strlen(filename)-4], ".pmi") == 0) pmi = true;
		//else if (strcmp (&filename[strlen(filename)-4], ".dat") == 0) dat = true;
		else if (strcmp (&filename[strlen(filename)-4], ".geo") == 0) geo = true;
	}
	
	FILE *input = fopen (filename, "r");
	if (input == NULL) fatal ("readInstance", "could not open input file");
	if (pmm) {
		instance = new PMMatrixInstance ();
		((PMMatrixInstance*)instance)->readPMM (input, p);
	} else if (tsp) {
		instance = new PMEuclideanInstance ();
		((PMEuclideanInstance*)instance)->readTSP (input, p);
	} else if (pmi) {
		instance = new PMEuclideanInstance ();
		((PMEuclideanInstance*)instance)->readPMI (input, p);
	} else if (geo) {
		instance = new PMGeoInstance ();
		((PMGeoInstance*)instance)->readGEO(input, p);
	} else if (dimacs) {
		instance = new PMGraphInstance ();
		((PMGraphInstance*)instance)->readDimacs(input, p);
	} else {
		fatal ("readInstance", "format not supported");
	}

	fclose(input);
	return instance;
}



/****************************************************************************
 *
 * updateClosest: finds the closest and the second closest facilities
 *                with respect to u (considering only the facilities in 
 *   solution s). Both arrays must be previously initialized with something 
 *   that makes sense (for instance, with zeros).
 * 
 *   running time: O(p)
 *
 ***************************************************************************/

void updateClosest (PMSolution *s, int u, int *closest, int *closest2) {
	int f, i, p;
	PMInstance *inst = s->getInstance();
	
	p = s->getP(); //number of facilities

	if (closest2 != NULL) {
		closest[u] = closest2[u] = 0;
		for (i=1; i<=p; i++) {
			f = s->getFacility(i); //f is the current facility
			double duf = inst->getDist(u,f);
			if (duf < inst->getDist(u,closest2[u])) {
				if (duf < inst->getDist(u,closest[u])) {
					closest2[u] = closest[u];
					closest[u] = f;
				} else closest2[u] = f;
			}
		} 
		if (closest[u] == closest2[u]) {
			fprintf (stderr, "[%d,%d,%8.2f]!\n", u, closest[u], inst->getDist(u,closest[u]));
			exit (-1);
		}
	} else {
		closest[u]=0;
		for (i=1; i<=p; i++) {
			f = s->getFacility(i);
			if (inst->getDist(u,f) < inst->getDist(u,closest[u])) closest[u] = f;
		}
	}	
}



void initializeClosest (PMSolution *s, int *closest, int *closest2) {
	int i;
	int n = s->getN();
	
	for (i=0; i<=n; i++) closest[i] = 0;
	if (closest2 != NULL) for (i=0; i<=n; i++) closest2[i] = 0;	

	for (i=1; i<=n; i++) updateClosest (s, i, closest, closest2);
}	






/*********************************************
 *
 * remove c random elements from the solution
 *
 *********************************************/

void removeRandom (PMSolution *s, int c) {
	int i, p, *v;
	p = s->getP();
	if (c>=p) s->reset(); 
	else {
		v = new int[p+1];
		for (i=1; i<=p; i++) v[i] = s->getFacility(i);
		RFW::shuffle (v, p);
		for (i=1; i<=c; i++) s->remove(v[i], false);	
		delete [] v;
	}
}




/***************************************************************
 *
 * updateArrays: updates the contents of save, lose, and extra
 * - checks only the neighborhood of vertex j
 * - uses closest and closet2, which must be up-to-date
 * - this function is a subroutine of localSearch
 *
 ***************************************************************/

void updateArrays (
  PMSolution *s,     //input: current solution
  bool *affected,    //input: list of vertices whose neighborhood we need to check
  SearchTables *tab, //input/output: these are the arrays that will be updated
  bool add           //input: if true, add contributions; if false, subtract them
) {
	double *save = tab->save;
	double *lose = tab->lose;
	PMInstance *inst = s->getInstance();
	int n = s->getN();
	//int m = s->getM();			 //number of potential facilities
	
	for (int j=1; j<=n; j++) {
		if (!affected[j]) continue;
		int r = s->getClosest(j);           //element that is closest to j
		double d1 = s->getDistClosest (j);  //distance to the closest
		double d2 = s->getDistClosest2 (j); //distance to the second closest
		tab->setCurrent(r);
		int *closelist = inst->getCloser(j, d2); //case 1

		if (add) {
			lose[r] += (d2-d1);
			for (int x=1; x<=closelist[0]; x++) { 
				int i = closelist[x];
				if (s->contains(i)) continue; //not in the solution	
				double dji = inst->getDist(j,i); 
				double profit = d1-dji;   
				double correction = d2;
				if (profit>0) {correction-=d1; save[i]+=profit;} //case 3
				         else {correction-=dji;}                 //case 2
				tab->addExtra(i,correction);
			}
		} else {
			lose[r] -= (d2-d1);
			for (int x=1; x<=closelist[0]; x++) { 
				int i = closelist[x];
				if (s->contains(i)) continue; 
				double dji = inst->getDist(j,i);			
				double profit = d1-dji;
				double correction=d2;
				if (profit>0) {correction-=d1; save[i]-=profit;} //case 3
				         else {correction-=dji;}                 //case 2
				tab->addExtra(i,-correction);
			}
		}
	}
}
















/*************************************************************
 *
 * bestNeighbor: use save, lose, and extra to find the best
 *               neighbor in O(p(m-p)) time. Returns the profit
 * that would be obtained if bestr was replaced with besti.
 * This function is used by localSearch.
 * 
 *************************************************************/

double bestNeighbor (
	PMSolution *s,     //input
	SearchTables *tab, //input: save, lose, extra...
	int *rcand,        //input (list of candidates for removal)
	int *icand,        //input (list of candidates for insertion)
	int *bestr,        //output: vertex to be removed
	int *besti,        //output: vertex to be inserted
	bool first         //input: use first improvement?
) {
	double *save = tab->save;
	double *lose = tab->lose;
	double **extra = tab->extra;

	double bestprofit = -INFINITY;

	if (first) RFW::shuffle (rcand, rcand[0]);

	*bestr = 0;
	*besti = 0;
	
	if (tab->usingList()) {
		//fprintf (stderr, "Finding best neighbor (%d,%d)... ", rcand[0], icand[0]);

		int j;
		//------------------------------------------------
		// build incidence vector of insertion candidates
		//------------------------------------------------
		int m = s->getInstance()->getM();
		bool *iscand= new bool [m+1];
		for (j=0; j<=m; j++) iscand[j] = false;
		for (j=icand[0]; j>0; j--) iscand[icand[j]] = true;
		
		//----------------------------------------------------
		// calculate maxsave and minlose; this will give the
		// best profit assuming all extras are zero; positive
		// extras can only improve the result
		//----------------------------------------------------
		double minlose = INFINITY;
		for (j=rcand[0]; j>0; j--) {
			int r = rcand[j];
			if (lose[r] < minlose) {
				*bestr = r;
				minlose = lose[r];
			}
		}
		double maxsave = -1;
		for (j=icand[0]; j>0; j--) {
			int i = icand[j];
			if (save[i] > maxsave) {
				*besti = i;
				maxsave = save[i];
			}
		}
		bestprofit = save[*besti] - lose[*bestr];

		//fprintf (stderr, "P(%d,%d) ", *bestr, *besti);

		//-------------------------------------------------
		// find the best moves considering "extra" factors
		//-------------------------------------------------
		for (j=rcand[0]; j>0; j--) {
			int r = rcand[j];
			ExtraNode *current = tab->list[r]->next;
			while (current->node <= m) {
				int i = current->node;
				if (iscand[i]) {
					double profit = save[i] - lose[r] + current->value;//extra[r][i];
					if (profit > bestprofit) {
						bestprofit = profit;
						*bestr = r; *besti = i;
					}
				}
				current = current->next;
			}
		}
		delete [] iscand;
		//fprintf (stderr, "done.");
	} else {
		for (int j=rcand[0]; j>0; j--) {
			int r = rcand[j];                //facility considered for removal
			for (int k=icand[0]; k>0; k--) {
				int i = icand[k];        //candidate for insertion
				double profit = save[i] - lose[r] + extra[r][i];
				if (profit > bestprofit) {
					bestprofit = profit;
					*bestr = r; *besti = i;
				}			
			}
			if (first && bestprofit > EPSILON) break;
		}
	} 
	return bestprofit;
}


/*************************************************************
 *
 * search: in each iteration, removes one facility from
 *         the solution and inserts another one to keep
 * its place. The pair of facilities involved is the one that
 * minimizes the total profit. Each iteration runs in O(n^2)
 * time. The stop criterion depends on sm.
 *
 *************************************************************/

double search (PMSolution *s, bool best_improvement, PMSearch *sm, PMElite *elite, bool use_list) {
	bool debug = false;
	
	enum {DOWN, SAME, UP};
	bool in_valley = false;

	int dir;

	//-----------------------
	// variable declarations
    //-----------------------
	BossaTimer tc (false);
	BossaTimer tb (false);
	BossaTimer tu (false);
	//BossaTimer t (true);
	BossaTimer tr (true);

	SearchTables *tab;
	bool *affected;         //vertices affected (incidence list)
	int *icand;             //candidates for insertion
	int *rcand;             //candidates for removal
	int i, n, m, p; 
	double start_cost = -1; //cost of the original solution (only if needed)
	double tot_profit = 0;  //highest profit in history
	double cur_profit = 0;  //current profit (w.r.t. original solution)

 	PMInstance *inst = s->getInstance(); //the instance we're dealing with
	n = s->getN(); //total number of cities
	m = s->getM(); //number of potential facilities
	p = s->getP(); //number of facilities to select

	//------------
	// allocation 
	//------------
	affected = new bool    [n+1];
	tab      = new SearchTables (m, p, use_list);
	icand    = new int     [m+1];
	rcand    = new int     [m+1];
	if (!tab->usingList()) {
		for (i=1; i<=m; i++) //allocate the rows we need in extra (only p rows are needed)
			tab->extra[i] = (s->contains(i) ? new double[m+1] : NULL);
	}

	dir = UP;

	//----------------
	// initialization
	//----------------
	for (i=1; i<=n; i++) affected[i] = true; //everything must be initialized

	for (i=1; i<=m; i++) {
		tab->save[i] = tab->lose[i] = 0; 
		if (!tab->usingList()) {
			if (tab->extra[i] != NULL) 
				for (int j=1; j<=m; j++) tab->extra[i][j] = 0;
		}
	}

	s->ensureConsistency();

	start_cost = s->calcCost();


	bool success = true;
	while (success) {

		//------------------------------
		// update save, lose, and extra	
		//------------------------------		
		tr.pause(); tu.resume();
		updateArrays (s, affected, tab, true);
		tu.pause();	tr.resume();

		//------------------------------------------------------------
		// build list of candidates for removal and insertion
		// the set depends on what search method we're dealing with
		//------------------------------------------------------------
		sm->getCandidates (rcand, icand);

		//--------------------------
		// find the best candidates
		//--------------------------
		int bestr, besti;
		double profit;
		tr.pause(); tb.resume();
		profit = bestNeighbor (s, tab, rcand, icand, &bestr, &besti, sm->firstImprovement());
		tb.pause(); tr.resume();
		
		//------------------------------------
		// check if that's indeed a good move
		//------------------------------------
		success = sm->reportMove (bestr, besti, profit);

		/*
		if (profit > 0) { //improving
			dir = DOWN;
		} else { //not improving
			if (dir==DOWN && elite!=NULL) { //if we were improving before, we've found a local optimum
				elite->add(s, s->calcCost());
			}
			if (profit < 0) {dir = UP;}
			else {dir = SAME;}
		}*/

		if (elite!=NULL) {
			if (profit>0) { //just went down
				in_valley = true;
			} else if (profit < 0) { //went up --- no longer in value
				if (in_valley) elite->add(s, s->calcCost()); //we have a local optimum in our hands
				in_valley = false;
			}
		}

		//----------------
		// compute profit
		//----------------
		cur_profit += profit;
		if (cur_profit>tot_profit) tot_profit = cur_profit;


		//--------------------------------------------------------
		// update things for the next iteration (if there is one)
		//--------------------------------------------------------
		if (success) {
			//------------------------------------------------------------------
			// we'll have to determine all vertices that were somehow affected
			//------------------------------------------------------------------
			for (i=1; i<=n; i++) {
				affected[i] = 
					(s->getClosest(i)==bestr)  || //closest vertex removed
					(s->getClosest2(i)==bestr) || //second closest vertex removed
					(s->getDistClosest2(i) > inst->getDist(i,besti)); // new vertex is very close
			}

			//------------------------------------------------------
			// for every vertex that was affected, we must subtract
			// its contributions to save, lose, and extra. New 
			// contributions will be computed by updateArray later
			//------------------------------------------------------
			tr.pause(); tu.resume();
			updateArrays (s, affected, tab, false);
			tu.pause(); tr.resume();
			
			//------------------------------------
			// now we finally change the solution
			//------------------------------------
			tr.pause(); tc.resume();
			s->swap (bestr, besti, true);
			tc.pause(); tr.resume();

			//------------------------------
			// add the solution to the pool
			//------------------------------ 
			//if (elite!=NULL) elite->add (s, start_cost - cur_profit);

			//----------------------------------------------------
			// lose[besti] and save[bestr] were undefined before,
			// we have to initialize them now. The actual values
			// will be calculated later in by updateArray
			//----------------------------------------------------
			tab->lose[besti] = tab->save[bestr] = 0;
			
			//-----------------------------------------------------------
			// the search table must know it's dealing with a different
			// set of vertices now
			//-----------------------------------------------------------
			tab->reportMove(besti, bestr);
		}
	}

	if (debug) fprintf (stderr, "[u:%.2f b:%.2f r:%.2f c:%.2f] ", tu.getTime(), tb.getTime(), tr.getTime(), tc.getTime());

	//-------------------
	// deallocate arrays
	//-------------------
	delete [] icand;
	delete [] rcand;
	delete [] affected;
	if (!tab->usingList()) {
		for (i=1; i<=m; i++) 
			if (tab->extra[i] != NULL) delete [] tab->extra[i];
	}
	delete tab;

	//---------------------------
	// return how much was saved
	//---------------------------
	return tot_profit;
}



/*******************************************
 *
 * consider each star to be a cluster and
 * calculate its 1-median
 *
 *******************************************/

 bool calc1medians (PMSolution *s) {
	PMInstance *inst = s->getInstance();
	int n = inst->getN();
	int m = inst->getM();
	int p = inst->getP();
	int i;
	int *color = new int [n+1];
	int *map   = new int [m+1]; //map[x] = rank of x among all facilities

	double cost = s->calcCost();

	fprintf (stderr, "%f -> ", cost);

	//initially, every customer painted with color of closest element
	for (i=1; i<=n; i++) {color[i] = s->getClosest(i);}

	//closest element will be mapped to its position in list
	for (i=1; i<=m; i++) {map[s->getFacility(i)] = i;}

	//now every customer will have a color in interval [1,p]
	for (i=1; i<=m; i++) color[i] = map[color[i]];

	//now calculate the 1-median of each cluster
	double *value = new double [m+1];
	value[0] = INFINITY;
	for (i=1; i<=m; i++) value[i] = 0;

	//calculate how much each center would cost
	for (i=1; i<=m; i++) {
		int c = color[i];
		for (int j=1; j<=m; j++) { //this is really ugly
			if (color[j]!=c) continue;	
			value[i] += inst->getDist(j,i);
		}
	}

	//determine best center for each facility
	int *best = new int [p+1];
	for (i=1; i<=p; i++) best[i] = 0;
	for (i=1; i<=m; i++) {
		int c = color[i];
		if (value[i] < value[best[c]]) best[c] = i;
	}

	s->reset();

	for (i=1; i<=p; i++) {
		s->add(best[i], true);
	}

	if (s->getP() < p) PMConstructive::addRandom (s, p - s->getP());

	fprintf (stderr, "%f\n", s->calcCost());

	//remember there may be fewer than n guys
	delete [] best;
	delete [] value;
	delete [] color;
	delete [] map;

	return (cost > s->calcCost() + EPSILON);

}

bool clusterSearch (PMSolution *s) {
	do {} while (calc1medians(s)); 
	return true;

	PMInstance *inst = s->getInstance();
	int n = inst->getN();
	int m = inst->getM();
	int p = inst->getP();
	int i;
	int *color = new int [n+1];
	int *map   = new int [m+1]; //map[x] = rank of x among all facilities

	double cost = s->calcCost();

	fprintf (stderr, "%f -> ", cost);

	//initially, every customer painted with color of closest element
	for (i=1; i<=n; i++) {color[i] = s->getClosest(i);}

	//closest element will be mapped to its position in list
	for (i=1; i<=m; i++) {map[s->getFacility(i)] = i;}

	//now every customer will have a color in interval [1,p]
	for (i=1; i<=m; i++) color[i] = map[color[i]];

	//now calculate the 1-median of each cluster
	double *value = new double [m+1];
	value[0] = INFINITY;
	for (i=1; i<=m; i++) value[i] = 0;

	//calculate how much each center would cost
	for (i=1; i<=m; i++) {
		int c = color[i];
		for (int j=1; j<=m; j++) { //this is really ugly
			if (color[j]!=c) continue;	
			value[i] += inst->getDist(j,i);
		}
	}

	//determine best center for each facility
	int *best = new int [p+1];
	for (i=1; i<=p; i++) best[i] = 0;
	for (i=1; i<=m; i++) {
		int c = color[i];
		if (value[i] < value[best[c]]) best[c] = i;
	}

	s->reset();

	for (i=1; i<=p; i++) {
		s->add(best[i], true);
	}

	if (s->getP() < p) PMConstructive::addRandom (s, p - s->getP());

	fprintf (stderr, "%f\n", s->calcCost());

	//remember there may be fewer than n guys
	delete [] best;
	delete [] value;
	delete [] color;
	delete [] map;

	return (cost > s->calcCost() + EPSILON);
 }








void localSearch (PMSolution *s, bool best_improvement, bool uselist, int *itcount=NULL, bool stats=false) {
	int it;
	BossaTimer t;
	bool verbose = false;

	t.start();
	if (verbose) fprintf (stderr, "Allocating local search object...");
	PMSearchLocal *sm = new PMSearchLocal (s);
	if (verbose) fprintf (stderr, "Running local search...");
	search (s, best_improvement, sm, NULL, uselist);
	if (verbose) fprintf (stderr, "end of local search.\n");
	it = sm->getIterations();
	if (itcount != NULL) *itcount = it;
	if (verbose) fprintf (stderr, "[%d]  ", it);

	delete sm;
	if (stats) {
		fprintf (stdout, "ls %.2f\n", s->calcCost());
		fprintf (stdout, "lstime %.2f\n", t.getTime()); 
		fprintf (stdout, "lsit %d\n", it);
	}
	if (verbose) fprintf (stderr, "[ls:%.2f] ", t.getTime());
}




/************************************************************
 *
 * performs path relinking from s1 to s2; result goes to t
 *
 ************************************************************/

double combinePR (PMSolution *t, PMSolution *s1, PMSolution *s2) {
	//BossaTimer timer(true);
	PMInstance *inst = t->getInstance();
	PMSearchPathRelinking *sm = new PMSearchPathRelinking (inst);
	PMElite *single_elite = new PMElite (inst, 1);

	sm->setSolutions (s1, s2);  //set both solutions
	single_elite->reset();      //reset set of elite solutions
	t->copyFrom (s1);
	search (t, true, sm, single_elite, true);       //perform path-relinking

	double cost;
	if (single_elite->countConsistent() == 0) { //no local minimum was found
		if (BossaRandom::getInteger(1,2) == 1) {
			t->copyFrom(s1);
			cost = s1->calcCost();
		} else {
			t->copyFrom(s2);
			cost = s2->calcCost();
		}
	} else {
		cost = single_elite->getSolution (1, t); 
	}
		
	//deallocate structures
	delete single_elite;
	delete sm;

	//fprintf (stderr, "<%.2f>", timer.getTime());
	return cost;

}


/********************************************************
 *
 * Combines s1 and s2, creating a new solution t, by
 * solving a subproblem containing only facilities that
 * are different in both cases
 *
 ********************************************************/

double combineSub (PMSolution *t, PMSolution *s1, PMSolution *s2) {	
	double t_init, t_solve, t_recover;

	BossaTimer timer(true);
	PMHybridInstance *subinst = new PMHybridInstance (s1, s2);
	PMSolution       *subsol  = new PMSolution (subinst);
	PMConstructive constructive;
	t_init = timer.getTime();

	timer.start();
	constructive.setMethod ("rpg:1");
	//PMElite *subelite = new PMElite (subinst, 4);
	grasp (subsol, &constructive, 4, NULL, "s", "r:u", false, false);
	t_solve = timer.getTime();
	
	timer.start();
	subinst->restoreSolution (subsol, t);

	//delete subelite;
	delete subsol;
	delete subinst;

	double cost = t->calcCost();
	t_recover = timer.getTime();

	fprintf (stderr, "[%.2f:%.2f:%.2f] ", t_init, t_solve, t_recover);
	return cost;
}


double combineSubPR (PMSolution *t, PMSolution *s1, PMSolution *s2) {
	double t_init, t_initsub, t_build, t_solve, t_recover;

	BossaTimer timer(true);
	bool verbose = false;

	
	//PMInstance *inst = t->getInstance();
	PMHybridInstance *subinst = new PMHybridInstance (s1, s2);
	t_init = timer.start();
	
	PMSolution *sub_s1 = new PMSolution (subinst);
	PMSolution *sub_s2 = new PMSolution (subinst);
	PMSolution *sub_t  = new PMSolution (subinst);
	t_initsub = timer.start();

	subinst->original2sub (s1, sub_s1);
	subinst->original2sub (s2, sub_s2);
	t_build = timer.start();

	combinePR (sub_t, sub_s1, sub_s2);
	t_solve = timer.start();
	subinst->restoreSolution (sub_t, t);
	t_recover = timer.getTime();
	delete sub_t;
	delete sub_s2;
	delete sub_s1;
	delete subinst;

	double cost = t->calcCost();
	if (verbose) fprintf (stderr, "Solution value is %.2f.\n", cost);
	fprintf (stderr, "[%.2f:%.2f:%.2f:%.2f:%.2f] ", t_init, t_initsub, t_build, t_solve, t_recover);
	return cost;
}






double combinePartial (PMSolution *t, PMSolution *s1, PMSolution *s2) {
	PMInstance *inst  = t->getInstance();
	PMSolution *start = new PMSolution (inst);
	PMSolution *best  = new PMSolution (inst);
	double best_value;
	int p = inst->getP();
	//int n = inst->getN();


	int i;
	//----------------------------------------
	// build list of candidates for insertion
	//----------------------------------------
	int *cand = new int [p+1];
	cand[0] = 0;
	for (i=1; i<=s2->getP(); i++) {
		int f = s2->getFacility(i);
		if (s1->contains(f)) continue;
		cand[0]++;
		cand[cand[0]] = f;
	}

	best->copyFrom (s1);
	best_value = INFINITY; //	best->calcCost();
	
	int maxdiff = 3;
	if (maxdiff > cand[0]/2) maxdiff = cand[0]/2;
	int min = 0;

	for (int it=1; it<=5; it++) {
		int x = BossaRandom::getInteger (min, maxdiff);
		if (x==0) min = 1;
		start->copyFrom (s1);
		
		partialShuffle (cand, 1, cand[0], x);
		for (int i=1; i<=x; i++) {
			int in = cand[i];
			int out;
			double profit;
			move (start, in, &out, &profit);
			start->swap(out, in, true);
			//fprintf (stderr, "(%d->%d) ", out, in);
		}
		double v = combinePR (t, start, s2);
		fprintf (stderr, "[%.2f]", v);
		if (v < best_value) {
			best_value = v;
			best->copyFrom (t);
		}
	}

	t->copyFrom(best);
	delete [] cand;
	delete start;
	delete best;
	return t->calcCost();

}


/******************************************
 *
 * Combine two solutions to find a new one
 *
 ******************************************/

double combine_r (PMSolution *t, PMSolution *s1, PMSolution *s2, const char* method) {
	if (method[0] != 'r') fatal ("combine_r", "invalid method");
	if (method[1] != ':') fatal ("combine_r", "wrong format");

	if (method[2] == 'b') {
		double v1 = combinePR (t, s1, s2);
		PMSolution *tcopy = new PMSolution (t->getInstance());
		tcopy->copyFrom(t);
		double v2 = combinePR (t, s2, s1);
		if (v2 > v1) t->copyFrom (tcopy);
		delete tcopy;
		return t->calcCost();
	} else {
		bool direct = true; //one to two?
		switch (method[2]) {
			case 'u': direct = (s1->calcCost()<=s2->calcCost()); break;
			case 'd': direct = (s1->calcCost()>s2->calcCost()); break;
			case 'r': direct = (BossaRandom::getInteger(0,1) == 1); break;
			case 's': direct = true; break;
			case 'i': direct = false; break;
			default: fatal ("combine_r", "invalid parameter for method 'r'");
		}
		if (direct) return combinePR (t, s1, s2);
		else return combinePR (t, s2, s1);
	}
}



double combine (PMSolution *t, PMSolution *s1, PMSolution *s2, const char* method) {
	if (s1->computeDifference(s2) < 2) {
		t->copyFrom(s1);
		return t->calcCost();
	}
	if (method==NULL) fatal ("combine", "invalid method");

	switch (method[0]) {
		case 's': return combineSub (t, s1, s2); break;
		case 'p': return combinePartial (t, s1, s2); break;
		case 'r': return combine_r (t, s1, s2, method); break;
		case 'X': return combineSubPR (t, s1, s2); break;
		case 'n': t->copyFrom(s1); return t->calcCost(); break;

		default: 
			fatal ("combine", "unrecognized method"); 
	}
	return 0; //never reached; here just to avoid complaints from the compiler
}
	


/**************************************************
 *
 * combine: merges two solutions (s1, s2) creating
 *          a third one.
 *
 **************************************************/

double combine2 (PMSolution *t, PMSolution *s1, PMSolution *s2) {
	PMInstance *inst = t->getInstance();
	int ins;
	int i;
	int p = inst->getP();
	int *candidates = new int [p+1];

	fprintf (stderr, "Combining solutions...\n");

	//--------------------------------------------
	// add p/2 facilities from the first solution
	//--------------------------------------------
	for (i=1; i<=p; i++) candidates[i] = s1->getFacility(i);
	RFW::shuffle (candidates, p); //we'll select the facilities at random

	t->reset();
	ins = 0;
	while (ins < p/2) {
		ins++;
		t->add(candidates[ins], true);
	}

	//---------------------------------------------
	// add p/2 facilities from the second solution
	//---------------------------------------------
	for (i=1; i<=p; i++) candidates[i] = s2->getFacility(i);
	RFW::shuffle (candidates, p);

	i = 0;
	while (ins < p) {
		i++;
		int c = candidates[i];
		if (t->contains(c)) continue;
		ins++;
		t->add(c, true);
	}

	delete [] candidates;
 
	fatal ("combine2", "deprecated function");
	//localSearch (t, true, NULL);

	fprintf (stderr, "done combining.\n");

	return t->calcCost();
}



/***************************************************************
 *
 * pathRelinking: execute multiple-generation path-relinking, 
 *     considering all pairs of solutions in 'elite'. The best 
 *     solution found in the process is returned in best_s. The 
 *     contents of 'elite' may change.
 *
 ***************************************************************/

void pathRelinking (PMSolution *best_s, PMElite *elite, const char *method, bool stats = false) {
	if (elite == NULL) fatal ("pathRelinking", "list of elite solutions cannot be NULL");

	BossaTimer timer;
	timer.start();
	
	bool advance = true;
	bool debug = false;
	
	PMElite *local_elite;
	PMSolution *s, *source, *target;
	PMInstance *inst = elite->getInstance();

	elite->getSolution (1, best_s);
	int capacity = elite->getCapacity();

	local_elite  = new PMElite (inst, capacity); //will hold the next generation
	s      = new PMSolution (inst); //termporary solution
	source = new PMSolution (inst); //first solution in the path
	target = new PMSolution (inst); //last solution in the path

	double vsource, vtarget; //, profit;
	double best = elite->getSolution (1, source);
	double oldbest = best + 1;

	int generation = 0;

	//save original best solution
	elite->getSolution(1, best_s);

	while (advance) {
		advance = false;
		generation ++;
		int count = elite->countConsistent();
		fprintf (stderr, "GENERATION %d\n", generation);
		oldbest = best;
		local_elite->reset();
		//for (int i=1; i<count; i++) {	
		
		for (int i=count; i>1; i--) {
			//get source
			vsource = elite->getSolution (i, source);
			if (debug) elite->checkConsistency (i);

			for (int j=i-1; j>=1; j--) {
			//for (int j=i+1; j<=count; j++) {
				//if (i==j) continue;

				if (debug) {
					elite->checkConsistency(j);
					fprintf (stderr, "Combining solutions (%d:%d); there are %d consistent solutions..\n", i, j, count);
				}
				//------------
				// get target
				//------------
				vtarget = elite->getSolution (j, target);
				PMSolution *r;
				double rvalue1, rvalue2; //before and after local search

				int diff = source->computeDifference(target);


				//--------------------------------------------------------------------------
				// perform path-relinking: best solution will be in r, its value in rvalue2
				//--------------------------------------------------------------------------
				rvalue1 = rvalue2 = combine (s, source, target, method);
				r = s;
				if (r->computeDifference(source)>1 && r->computeDifference(target)>1) {
					localSearch (r, true, true);  //try to improve solution with local search
					rvalue2 = r->calcCost();      //get the value
				}

				//-------------------
				// update statistics
				//-------------------
				//current = vsource - profit;
				if (rvalue2 < best - EPSILON) {
					best = rvalue2;
					fprintf (stderr, "> ");
				} else {
					fprintf (stderr, "  ");
				}
				fprintf (stderr, "(%02d:%.2f) x (%02d:%.2f) [%d]: %.2f", i, vsource, j, vtarget, diff, rvalue1);
				if (rvalue2 < rvalue1 - EPSILON) fprintf (stderr, " -> %.2f", rvalue2);
				fprintf (stderr, "\n");

				//------------------------------------
				// store solution for next generation
				//------------------------------------
				local_elite->add (r, rvalue2);
			}	
		}
		
		//---------------------------------------------------
		// if there is improvement, prepare a new generation
		//---------------------------------------------------
		if (best < oldbest-EPSILON) {
			advance = true;	                         //there is a new generation
			local_elite->getSolution (1, best_s);    //best_s holds the best solution now
			elite->reset();                          //remove old elite solutions
			count = local_elite->countConsistent();  
			for (int i=1; i<=count; i++) {           //add new set of elite solutions
				vsource = local_elite->getSolution (i, source);
				elite->add (source, vsource);
			}
		}
	}
			

	//------------------
	// print statistics
	//------------------
	if (stats) {
		fprintf (stdout, "generations %d\n", generation);
		fprintf (stdout, "prmethod %s\n", method); 
		fprintf (stdout, "prelite %d\n", elite->getCapacity());
		fprintf (stdout, "prtime %.3f\n", timer.getTime());
		fprintf (stdout, "pr %.2f\n", best);
	}

	//----------
	// clean up 
	//----------
	delete local_elite;
	delete s;
	delete source;
	delete target;
}



/***************************************************************
 *
 * pathRelinking_target: execute multiple-generation path-relinking, 
 *     considering all pairs of solutions in 'elite'. The best 
 *     solution found in the process is returned in best_s. The 
 *     contents of 'elite' may change.
 *
 ***************************************************************/
void pathRelinking_target (PMSolution *best_s, PMElite *elite, const char *method, double alvo, bool stats = false) {
	if (elite == NULL) fatal ("pathRelinking", "list of elite solutions cannot be NULL");

	BossaTimer timer;
	timer.start();
	
	bool advance = true;
	bool debug = false;
	
	PMElite *local_elite;
	PMSolution *s, *source, *target;
	PMInstance *inst = elite->getInstance();

	elite->getSolution (1, best_s);
	int capacity = elite->getCapacity();

	local_elite  = new PMElite (inst, capacity); //will hold the next generation
	s      = new PMSolution (inst); //termporary solution
	source = new PMSolution (inst); //first solution in the path
	target = new PMSolution (inst); //last solution in the path

	double vsource, vtarget; //, profit;
	double best = elite->getSolution (1, source);
	double oldbest = best + 1;

	int generation = 0;

	//save original best solution
	elite->getSolution(1, best_s);

	while ((best > alvo) && (elite->countConsistent() != 0)) {
		//advance = false;
		generation ++;
		int count = elite->countConsistent();
		fprintf (stderr, "GENERATION %d\n", generation);
		oldbest = best;
		local_elite->reset();
		//for (int i=1; i<count; i++) {	
		
		for (int i=count; i>1; i--) {
			//get source
			vsource = elite->getSolution (i, source);
			if (debug) elite->checkConsistency (i);

			for (int j=i-1; j>=1; j--) {
			//for (int j=i+1; j<=count; j++) {
				//if (i==j) continue;

				if (debug) {
					elite->checkConsistency(j);
					fprintf (stderr, "Combining solutions (%d:%d); there are %d consistent solutions..\n", i, j, count);
				}
				//------------
				// get target
				//------------
				vtarget = elite->getSolution (j, target);
				PMSolution *r;
				double rvalue1, rvalue2; //before and after local search

				int diff = source->computeDifference(target);


				//--------------------------------------------------------------------------
				// perform path-relinking: best solution will be in r, its value in rvalue2
				//--------------------------------------------------------------------------
				rvalue1 = rvalue2 = combine (s, source, target, method);
				r = s;
				if (r->computeDifference(source)>1 && r->computeDifference(target)>1) {
					localSearch (r, true, true);  //try to improve solution with local search
					rvalue2 = r->calcCost();      //get the value
				}

				//-------------------
				// update statistics
				//-------------------
				//current = vsource - profit;
				if (rvalue2 < best - EPSILON) {
					best = rvalue2;
					fprintf (stderr, "> ");
				} else {
					fprintf (stderr, "  ");
				}
				fprintf (stderr, "(%02d:%.2f) x (%02d:%.2f) [%d]: %.2f", i, vsource, j, vtarget, diff, rvalue1);
				if (rvalue2 < rvalue1 - EPSILON) fprintf (stderr, " -> %.2f", rvalue2);
				fprintf (stderr, "\n");

				//------------------------------------
				// store solution for next generation
				//------------------------------------
				local_elite->add (r, rvalue2);

				//Isabel - 17/03/11
				if(best <= alvo){ 
					if (best < oldbest-EPSILON) local_elite->getSolution (1, best_s);
					fprintf(stdout, "break 1\n"); 
					break; 
				}
			}
			//Isabel - 17/03/11
			if(best <= alvo){ 
				if (best < oldbest-EPSILON) local_elite->getSolution (1, best_s);
				fprintf(stdout, "break 2\n"); 
				break; 
			}
		}
		//Isabel - 17/03/11
		if(best <= alvo){ 
			if (best < oldbest-EPSILON) local_elite->getSolution (1, best_s);
			fprintf(stdout, "break 3\n"); 
			break; 
		}

		
		//---------------------------------------------------
		// if there is improvement, prepare a new generation
		//---------------------------------------------------
		//if (best < oldbest-EPSILON) {			 //criação da nova geração sempre 
			//advance = true; Isabel 17/03/11	 //there is a new generation
			local_elite->getSolution (1, best_s);    //best_s holds the best solution now
			elite->reset();                          //remove old elite solutions
			count = local_elite->countConsistent();  
			for (int i=1; i<=count; i++) {           //add new set of elite solutions
				vsource = local_elite->getSolution (i, source);
				elite->add (source, vsource);
			}
		//}
	}
			

	//------------------
	// print statistics
	//------------------
	if (stats) {
		fprintf (stdout, "generations %d\n", generation);
		fprintf (stdout, "prmethod %s\n", method); 
		fprintf (stdout, "prelite %d\n", elite->getCapacity());
		fprintf (stdout, "prtime %.3f\n", timer.getTime());
		fprintf (stdout, "pr %.2f\n", best);
	}

	//----------
	// clean up 
	//----------
	delete local_elite;
	delete s;
	delete source;
	delete target;
}



/*************************************************************
 * 
 * VNDS, as proposed by Hansen and Mladenovic
 *
 *************************************************************/

double vnds (PMSolution *s, PMElite *elite, int kmax, int bmax, bool verb) {
	fprintf (stderr, "Running vnds now, from a solution whose value is %.2f...\n", s->calcCost());

	//-----------------------
	// variable declarations
	//-----------------------
	bool stop;
	bool *insubfac;   
	BossaTimer t, uptimer, subtimer, consttimer;
	BossaVertexHeap <double> *heap;
	double best_cost, best_time; 
	int best_it, i, it, j, k, m, n, p;
	int *subfac, *subusers; //, *closest;   
	int **faclist; //faclist[i][j]: j-th closest facility to i
	PMInstance *inst; //, *subproblem;
	//PMSolution *subsol;
	PMConstructive constructive;
	//constructive.useProportionalWorst();
	constructive.setMethod ("pworst");

	uptimer.start();
	uptimer.pause();
	subtimer.start();
	subtimer.pause();
	consttimer.start();
	consttimer.pause();
	t.start();

	//--------------------------------------
	// define how much information to print
	//--------------------------------------
	bool stats = verb;
	int verbose;
	if (verb) verbose=1; 
	else verbose = 0;

	//------------------------------
	// initialize basic information
	//------------------------------
	s->ensureConsistency();
	inst = s->getInstance(); //current instance
	m = s->getM();           //number of potential facilities
	n = s->getN();           //number of users
	p = s->getP();           //number of facilities that must be selected 

	
	
	//int rmax  = 1000;  //number of RVNS iterations
	//int kmax  = (int)(ceil((double)p/log((double)p)));     //number of VNDS iterations
	
	if (kmax == -1) kmax = p;

	//int kmaxp = 2;     //kmax prime (for RVNS)
	//int kmaxpp = 5;    //kmax double prime (for VNS)
	//int b = 400;       //max size of the subproblem
	//int maxtries = (int)ceil ((double)log((double)p));
	
	//int maxtries = (int) ceil ((double)p / (double)kmax);
	
	subfac    = new int  [m+1]; //list of facilities in the subproblem (subfac[0] is the length)
	subusers  = new int  [n+1]; //list of users in the subproblem (subusers[0] is the length)
	insubfac  = new bool [m+1]; //insubfac[i] -> is facility 'i' in the subproblem?
	
	//initializeClosest (s, closest, NULL);



	//---------------------------------------------
	// allocate and initialize list of facilities
	// we're using heapsort here (just because we 
	// have a heap already)
	//---------------------------------------------
	if (verbose>1) fprintf (stderr, "Creating list of facilities...\n");
	heap = new BossaVertexHeap <double> (m); //will 
	faclist = new int *[m+1];
	for (i=1; i<=m; i++) {
		faclist[i] = new int [m+1];
		heap->reset();
		for (j=1; j<=m; j++) {
			if (i==j) continue;
			heap->insert (j, inst->getFacDist(i,j));
		}
		for (j=1; j<m; j++) {
			int v;
			double d;
			heap->removeFirst (v, d);
			faclist[i][j] = v;
		}
	}
	if (verbose>1) fprintf (stderr, "done.\n");
	delete heap;

	//------------------
	// initialize stats
	//------------------
	best_time = 0; //time in which the best solution was found
	best_it   = 0; //iteration in which the best solution was found
	best_cost = 0; //cost of the best solution found
	best_cost = s->calcCost();
	
	//------------
	// main loop
	//------------
	it = 0;
	stop = false;
	while (!stop) {
		k = 1; //int tries = 1;

		fatal ("vnds", "feature no longer implemented");
		/*
		while (k<=kmax) {
			it ++;
			if (verbose>0) {
				if (k==1) fprintf (stderr, "\n%.2f (best:%.2fs sub:%.2fs const:%.2fs up:%.2fs)", 
					best_cost, best_time, subtimer.getTime(), consttimer.getTime(),uptimer.getTime());
				fprintf (stderr, "%c", (k%10==0) ? ':' : '.');
			}

			//-------------------------------
			// choose a random facility (rf)
			//-------------------------------
			int rf = s->getFacility (BossaRandom::getInteger(1,p)); 

			//--------------------------------------------------------------------
			// create a list containing this facility and the k-1 closest medians
			//--------------------------------------------------------------------
			subfac[0] = 0;                             //nobody in this list...
			for (i=1; i<=m; i++) insubfac[i] = false;  //...in the beginning

			//-----------------------
			// add facilty rf itself
			//-----------------------
			subfac[0] = 1;
			subfac[1] = rf;
			insubfac[rf] = true;

			//-------------------------------------------------------------------
			// add other facilities (closest k-1 facilities in current solution)
			//-------------------------------------------------------------------
			i = 0;
			while (subfac[0]<k) {
				i++;
				int g = faclist[rf][i];
				if (!s->contains(g)) continue;
				subfac[0]++;
				subfac[subfac[0]] = g;
				insubfac[g] = true;
			}

			//------------------------------------------------
			// find the set of users whose closest facilities
			// belong to subfac
			//------------------------------------------------
			double original_cost = 0;
			subusers[0] = 0; //no user initially 
			for (i=1; i<=n; i++) { //run through the list of users, find those that are close enough
				if (insubfac[s->getClosest(i)]) {
					original_cost += inst->getDist(i,s->getClosest(i));
					subusers[0]++;
					subusers[subusers[0]] = i;
				}
			}

			//-----------------------------
			// create and solve subproblem
			//-----------------------------
			bool improved = false;
			if (subusers[0]>k) {

				//------------------------------------
				// create new problem an new solution
				//------------------------------------
				consttimer.resume();
				subproblem = new PMInstance (inst, subusers, subusers);
				subproblem->setP(k);
				subsol = new PMSolution (subproblem);
				consttimer.pause();

				//-------------------
				// solve the problem
				//-------------------
				subtimer.resume();
				if (k>1) {
					subsol->reset();
					for (i=1; i<=k; i++) {subsol->add(i, true);}

					if (subusers[0] > bmax) {
						//smartRandom (subsol, NULL, false);
						//vns (subsol, 2, subusers[0], true, NULL, false);
						grasp (subsol, &constructive, 1, NULL, false, false, false);
						  //fprintf (stderr, "(n=%d p1=%d p2=%d m=%d)", subproblem->getN(), subproblem->getP(), subsol->getP(), subproblem->getM());
						//fprintf (stderr, "Original cost: %f\n", subsol->calcCost());
						//vnds (subsol, NULL, kmax, bmax / 2, false);
					} else {
						//PMElite *elite = new PMElite (subproblem, 2);
						grasp (subsol, &constructive, 1, NULL, false, false, false);
						//greedy (subsol, NULL);
						delete elite;
					}
				} else {
					PMConstructive::greedy (subsol, false, 0.0, NULL);
				}
				subtimer.pause();

				double new_cost = subsol->calcCost();
				

				if (verbose>1) fprintf (stderr, "New solution cost: %.2f.\n", new_cost);

				//----------------------------------------------------
				// if there is an improvement, go to the new solution
				//----------------------------------------------------
				if (new_cost < original_cost - EPSILON) {
					uptimer.resume();
					improved = true;

					//-------------------------------------------------
					// remove all facilities that were in the solution
					//-------------------------------------------------
					for (i=1; i<=subfac[0]; i++) s->remove(subfac[i], true);

					//----------------------------------
					// insert the new set of facilities
					//----------------------------------
					for (i=1; i<=subsol->getP(); i++) {
						int fsub = subsol->getFacility(i); //facility number as seen by the subproblem
						int f = subusers[fsub];            //actual facility number
						if (s->contains(f)) fatal ("vnds", "invalid insertion");
						s->add(f, true);
					}

					//---------------------------------------------------
					// update closeness information for the new solution
					// (THIS CAN BE MADE MORE EFFICIENT)
					//---------------------------------------------------
					//for (i=1; i<=n; i++) closest[i] = 0;
					//for (i=1; i<=n; i++) updateClosest (s, i, closest, NULL);

					//-------------------
					// update statistics
					//-------------------
					//best_cost = 0;
					//for (i=1; i<=n; i++) best_cost+=inst->getDist(i,closest[i]);
					best_cost = s->calcCost();
					best_time = t.getTime();
					best_it = it;

					//----------------------------------------
					// include in the list of elite solutions
					//----------------------------------------
					if (elite!=NULL) elite->add(s,best_cost);
					uptimer.pause();
				}
			
				delete subsol;
				delete subproblem;
			} else {
				if (verbose>1) fprintf (stderr, "Trivial subproblem ignored.\n");
			}
			if (improved) {
				k = 1;
				tries = 1;
				improved = false;
			} else {
				tries ++;
				if (tries > maxtries) {
					tries = 1;
					k++;
				}
			}

		}
		*/
		stop = true;
	}


	if (stats) {
		fprintf (stdout, "vnds %.2f\n", best_cost);
		fprintf (stdout, "vndstime %.3f\n", t.getTime());
		fprintf (stdout, "vndsbesttime %.3f\n", best_time);
		fprintf (stdout, "vndsk %d\n", kmax);
		fprintf (stdout, "vndsb %d\n", bmax);
		fprintf (stdout, "vndsbestit %d\n", best_it);
	}

	//--------------------------------
	// deallocate lists of facilities
	//--------------------------------
	for (i=1; i<=m; i++) delete [] faclist[i];
	delete [] faclist;
	delete [] subfac;
	delete [] insubfac;
	delete [] subusers;

	return best_cost;

}




/****************************************************************
 *
 * Tabu Search: performs tabu search starting with a solution
 *              obtained by a constructive heuristic. Input
 * parameters are tenure for insertions and tenure for removals.
 *
 ****************************************************************/

void tabu (PMSolution *s, int it, int itenure, int rtenure, int nelite) {
	BossaTimer t;
	t.start();
	PMInstance *inst = s->getInstance();
	
	PMElite *elite;
	elite = (nelite>0) ? new PMElite (inst, nelite) : NULL;

	s->reset();
	PMConstructive::greedy (s);
	double cost = s->calcCost();
	PMSearchTabu *sm = new PMSearchTabu (s, it);
	sm->setInsertTenure (itenure);
	sm->setRemoveTenure (rtenure);
	double profit = search (s, true, sm, elite, true);

	delete sm;

	pathRelinking (s, elite, "r:u");
	delete elite;
	fprintf (stdout, "tabu %.2f\n", cost - profit);
	fprintf (stdout, "tabutime %.3f\n", t.getTime());
}



/*********************************************
 *
 * update routine for the FI heuristic, as
 * described in HM97
 *
 *********************************************/

void update (PMSolution *s, int in, int out) { //, int *closest, int *closest2) {
	PMInstance *inst = s->getInstance();
	int n = s->getN();
	int i;

	upcount++;
	
	if (!s->contains(in)) fatal ("update", "invalid insertion");
	if (s->contains(out)) fatal ("update", "invalid removal"); 
	
	bool *affected = new bool [n+1];
	for (i=1; i<=n; i++) affected[i] = false;

	for (i=1; i<=n; i++) {
		double din = inst->getDist(i,in);
		double d1 = inst->getDist (i, s->getClosest(i));
		double d2 = inst->getDist (i, s->getClosest2(i));
		
		//case a: we're removing the closest element
		if (s->getClosest(i) == out) {
			if (din <= d2) { //first case: in is the new closest
				s->setClosest(i, in); //closest2[i] remains unaltered, it isn't out
			} else { //second case: in is not the closest
				s->setClosest(i, s->getClosest2(i));
				s->setClosest2(i, 0);
				affected[i] = true; 
			}
		}
		//we're removing an arbitrary element
		else { 
			//in is the new closest
			if (din < d1) {
				s->setClosest2(i, s->getClosest(i));
				s->setClosest(i, in);
			} else {
				//in is closer than the second closest
				if (din < d2) {
					s->setClosest2(i,in); //in becomes new second closest
				} else { //last case: second closest was affected
					if (s->getClosest2(i)==out) {affected[i] = true;}
				}
			}
		}
	}	
	
	// now we just have to find the second 
	// closest for affected elements
	int count = 0;
	for (i = 1; i<=n; i++) {
		if (s->getClosest(i) == s->getClosest2(i)) fatal ("update", "something is wrong");
		if (affected[i]) {
			count++;
			s->updateClosest (i);
		}
	}

	delete [] affected;
}



/****************************************************
 *
 * fastInterchange: executes the fast interchange
 * local search procedure, implemented as suggested
 * in HM97. 'it' (an output parameter) is the number
 * of iterations actually performed.
 *
 ****************************************************/

void fastInterchange (PMSolution *s, int *itcount = NULL) {
	int i, it;
	int p = s->getP();
	int n = s->getN();
	
	it = 0;
	//PMInstance *inst = s->getInstance();
	int verbose = 0;

	double cur_cost = s->calcCost();
	bool success = true;

	BossaTimer t_timer;
	BossaTimer m_timer;

	t_timer.start();
	m_timer.start();
	m_timer.pause();
	
	while (success) {
		double best_delta = INFINITY;
		int best_in = 0;
		int best_out = 0;
		
		m_timer.resume();
		for (i = p+1; i<=n; i++) {
			int in = s->getFacility(i);
			int out;
			double delta;
			move (s, in,&out, &delta);
			if (delta < best_delta) {
				best_delta = delta;
				best_in = in;
				best_out = out;
			}
		}
		m_timer.pause();

		it++;
		
		success = (best_delta < -EPSILON);
	    
		if (success) {
			cur_cost += best_delta;
			s->remove(best_out, false);
			s->add(best_in, false);
			update (s, best_in, best_out);
		}
	}	
	if (verbose) fprintf (stderr, "[%d] ", it);

	if (itcount!=NULL) *itcount = it;
}


/*********************************************
 *
 * VNS, as suggested by Hansen and Mladenovic
 *
 *********************************************/
				
void vns (PMSolution *opt_s, int kmax, int rmax, bool rvns, PMElite *elite, bool stats) {
	int p = opt_s->getP();
	//int n = opt_s->getN();
	int m = opt_s->getM();
	int i, k;
	int best_it = 0;
	int bestup = 0;
	int *candidates = new int [m+1];
	int *quitters   = new int [m+1];
	bool *removed   = new bool [m+1];
	
	bool ls;
	bool stop = false;

	int verbose;
	if (stats) verbose = 2;
	else verbose = 0;
	ls = !rvns;

	int rcount = 0;

	if (kmax<=0) fatal ("vns", "kmax out of range");
	if (verbose>1) fprintf (stderr, "Running VNS...\n");

	BossaTimer t;
	t.start();
	
	PMInstance *inst = opt_s->getInstance();

	PMSolution *cur_s = new PMSolution (inst);
	opt_s->ensureConsistency(); //find closest and second closest for everybody
	cur_s->copyFrom (opt_s);
	
	
	//------------------
	// initialize costs
	//------------------
	double cur_cost = cur_s->calcCost();
	double opt_cost = opt_s->calcCost();
	if (cur_cost != opt_cost) fatal ("vns", "solutions should be identical");
	
	if (verbose>1) fprintf (stderr, "Initial solution value is %.2f.\n", opt_cost);
	
	//-----------
	// main loop
	//-----------
	int it = 0;
	int itmax = 5000;
	k = 1;
	

	while (!stop) {		
		if (verbose>0) fprintf (stderr, "it:%d  best:%.2f\n", it, opt_cost);
		k = 1;
		it++;
		while (k<=kmax && !stop) {
			//---------	
			// shaking
			//---------			
			for (i=1; i<=m-p; i++) candidates[i] = cur_s->getFacility(i+p);
			RFW::shuffle (candidates, m-p);
			
			if (verbose>2) {
				for (int i=1; i<=k; i++) fprintf (stderr, "%d ", candidates[i]);
				fprintf (stderr, "\n");
			}
			
			if (m-p < k) fatal ("vns", "not enough facilities to insert");

			//------------
			// move stuff
			//------------
			for (int j=1; j<=k; j++) {
				int goin, goout;
			    double delta;

				goin = candidates[j];
				move (cur_s, goin, &goout, &delta);
				quitters[j] = goout;
				cur_cost += delta;
				cur_s->add   (goin, true);
				cur_s->remove(goout, true); //closest and closest2 automatically updated
			}
			
			//--------------
			// local search
			//--------------
			if (ls) {
				int fit;
				localSearch (cur_s, true, true, &fit);
				it += fit;
				cur_cost = cur_s->calcCost();
			}
	
			if (verbose > 1) {
				fprintf (stderr, "it=%d(%d) up=%d(%d) k=%d(%d) cost=%.2f(%.2f) tpi=%.4f r=%d/%d\n", 
					it, best_it, upcount, bestup, k, kmax, cur_cost, opt_cost, t.getTime()/(double)it, rcount, rmax);
			}
			//-------------
			// move or not
			//-------------
			if (cur_cost < opt_cost - EPSILON) {  //success
				opt_s->copyFrom (cur_s); //move to current solution	
				opt_cost = cur_cost;     //save current cost
				best_it = it;            //this is the best iteration
				bestup = upcount;	
				k = 1;                   //k=1 again
			} else {  //case 2: failure
				k++;                     //try a farther neighborhoood
				cur_s->copyFrom (opt_s); //go back to the best known solution
				cur_cost = opt_cost;     //get old cost
				rcount ++;               //one more failure
				if (rcount>=rmax && rmax>0) stop = true; //maybe we've had enough
			}
		}
		if (it>=itmax) stop = true;
	}
	
	
	//------------------
	// print statistics
	//------------------
	if (stats) {
		if (rvns) {
			fprintf (stdout, "rvnsit %d\n", it);
			fprintf (stdout, "rvnsbestit %d\n", best_it);
			fprintf (stdout, "rvns %.2f\n", opt_cost);
			fprintf (stdout, "rvnstime %.3f\n", t.getTime());
			fprintf (stdout, "rvnstimeit %.3f\n", t.getTime() / it);
		} else {
			fprintf (stdout, "vnsit %d\n", it);
			fprintf (stdout, "vnsbestit %d\n", best_it);
			fprintf (stdout, "vns %.2f\n", opt_cost);
			fprintf (stdout, "vnstime %.3f\n", t.getTime());
			fprintf (stdout, "vnstimeit %.3f\n", t.getTime() / it);
		}
	}

	//----------
	// clean up
	//----------
	delete [] candidates;
	delete [] quitters;
	delete [] removed;
}	



/*****************************************************
 *
 * GRASP: executes a GRASP with n interations and
 *        possibly with path-relinking (if nelite 
 * is greater than 0). Uses fastInterchange if fi 
 * is true, otherwise uses localSearch. The best 
 * solution found (after path-relinking) is returned 
 * in s.
 *  
 *****************************************************/

double grasp (
  PMSolution *s,   //output: will hold the best solution found
  PMConstructive *constructive, //input: used to run the appropriate constructive heuristic
  int it,          //input: number of iterations to be performed
  PMElite *elite,  //input/output: will hold best solutions found and be used for path-relinking (can be NULL)
  char *lsmethod,   //local search method to be used
  char *combmethod,
  bool stats,      //input: print statistics to stdout?
  bool partials    //input: print partial information to stdout?
) {
	int verbose = stats ? 1 : 0;

	PMRecorder recorder;
	BossaTimer t(true);          //overall
	BossaTimer ls_timer(false);  //local search
	BossaTimer con_timer(false); //constructive
	BossaTimer comb_timer(false); //combination time

	PMInstance *inst = s->getInstance();
	PMElite *elite_best = new PMElite (inst, 1);
	// Gabriel - alteração 30/06/2008	- início 	
#ifdef MINERACAO
	PMEliteSimple *elite_simples;
	PMConstructiveMD constructive_md;
	elite_simples = new PMEliteSimple(inst, g_tam_elite_simples); // Elite de tamanho 4
	//constructive_md = new PMConstructiveMD;
	constructive_md.setSeed(g_seed);
	constructive_md.setMethod("sample:2");
	constructive_md.setTamEliteSimples(g_tam_elite_simples);
	constructive_md.setNumeroPadroes(g_npadroes);
	constructive_md.setSuporteMinimo(g_supmin);
#endif	
	// Gabriel - alteração 30/06/2008	- fim
	PMSolution *r  = new PMSolution (inst);
	PMSolution *rs = new PMSolution (inst); //relink solution
	
	int p = inst->getP();
	int n = inst->getN();
	int m = inst->getM();
	
	if (verbose > 1) fprintf (stderr, "GRASP: m=%d n=%d p=%d\n", m, n, p);
	if (verbose) fprintf (stderr, "Running grasp...\n");

		for (int i=1; i<=it; i++) {
		//-------------------------------
		// build a solution from scratch
		//-------------------------------
		con_timer.resume();
		s->reset();
		// Gabriel - alteração 30/06/2008 - início 
#ifdef MINERACAO
		if (i <= (it/2)) 
			constructive->run(s);
		else
		{
			if (i == ((it/2)+1))
			{
//				FILE *elite_file = fopen ("c:\\ProjetoFinal\\popstar2002\\Debug\\elite.txt", "w");
				FILE *elite_file = fopen ("elite.txt", "w");
				if (elite_file == NULL) 
					fatal ("elite_simples", "could not open input file");
				elite_simples->output(elite_file,true);
				fclose(elite_file); 					
			}
				constructive_md.run(s);
		}
#endif
#ifndef MINERACAO
		constructive->run(s);
#endif
		// Gabriel - alteração 30/06/2008 - fim  
		double original = s->calcCost();
		con_timer.pause();
		
		//----------------------
		// perform local search
		//----------------------
		ls_timer.resume();
		runLocalSearch (s, lsmethod);
		double result = s->calcCost();
		ls_timer.pause();

		if (verbose) fprintf (stderr, "%5d: [%.2f] %.2f -> %.2f", i, elite_best->getSolutionCost(1), original, result);

		unsigned int seed = BossaRandom::getNonZero();

		if (elite!=NULL) {
			comb_timer.resume();
			int d = elite->countConsistent();
			bool insert = true;

			//if there are elite solutions, relink
			if (d>=1) { 
				int e = elite->getRandomDifferent (s, 1, d); //a probably different solution
				
				if (e!=0) { //solution is not identical to some solution in the pool
					double elite_cost = elite->getSolution (e, r); //r is the solution
					double oldcost = (elite_cost < result) ? elite_cost : result; //get best of two original

					if (!r->equals(s)) { //if they are not equal, try to relink them
						//--------------------------------
						// combine both solutions into rs
						//--------------------------------
						double newcost = combine (rs, s, r, combmethod); //s is the one produced by local search, r is the elite one
						if (verbose) fprintf (stderr, " x %.2f -> %.2f", elite_cost, newcost);
						if (rs->computeDifference(r)>1 && rs->computeDifference(s)>1) {
							localSearch (rs, true, true); //try to improve by local search if necssary
							newcost = rs->calcCost();
							if (verbose) fprintf (stderr, " -> %.2f", newcost);
						}

						//-----------------------------------------------------------
						// if there is an improvement, try to save 'hybrid' solution
						//-----------------------------------------------------------
						if (newcost < oldcost - EPSILON) {
							elite->add (rs, newcost);      //add to the set of elite solutions
							elite_best->add (rs, newcost); //save it in the local list of elite solutions
							if (recorder.report (i, t.getTime(), newcost)) {
								if (partials) fprintf (stdout, "partial %.3f %.3f\n", recorder.best_time, recorder.best_value);
							}
						}
					}
				} else insert = false;
			}
			if (insert) elite->add (s, result);
			comb_timer.pause();
		}

		// Gabriel - alteração 30/06/2008 - início 
#ifdef MINERACAO
		if (i <= 16)
		{
			elite_simples->add(s, result);
			//elite_simples->outputScreen();
		}
#endif
		// Gabriel - alteração 30/06/2008 - fim

		//------------------------
		// save original solution
		//------------------------
		elite_best->add(s, result); 
		if (recorder.report (i, t.getTime(), result)) {
			if (partials) {fprintf (stdout, "partial %.3f %.3f\n", recorder.best_time, recorder.best_value);}
		}
		
		if (verbose) {
			if (recorder.best_iteration==i) fprintf (stderr, " <");
			fprintf (stderr, "\n");
		}
		BossaRandom::randomize(seed);		
	}
	
	//----------------------------------------------------
	// copy best solution to s (and its value to solcost)
	//----------------------------------------------------
	double solcost = elite_best->getSolution (1, s);
	
	//-----------------------
	// deallocate structures
	//-----------------------
	delete elite_best;
	delete r;
	delete rs;

	//--------------
	// print report
	//--------------
	if (stats) {
		fprintf (stdout, "grasp %.2f\n", solcost);
		fprintf (stdout, "graspconst ");
		constructive->printMethod(stdout);
		fprintf (stdout, "\n");
		fprintf (stdout, "lsmethod %s\n", lsmethod);
		fprintf (stdout, "graspcomb %s\n", combmethod);
		fprintf (stdout, "grasplstimeavg %.2f\n", ls_timer.getTime()/(double)it);
		fprintf (stdout, "graspconsttimeavg %.2f\n", con_timer.getTime()/(double)it);
		fprintf (stdout, "grasplstime %.2f\n", ls_timer.getTime());
		fprintf (stdout, "graspconsttime %.2f\n", con_timer.getTime());
		fprintf (stdout, "grasptime %.2f\n", t.getTime());
		fprintf (stdout, "graspit %d\n", it);
		fprintf (stdout, "graspbestit %d\n", recorder.best_iteration);
		fprintf (stdout, "graspbesttime %.2f\n", recorder.best_time); 
	}

	return solcost;
}


// Isabel - 30/08/2010
/*****************************************************
 *
 * GRASP_target: executes a GRASP with n interations and
 *        possibly with path-relinking (if nelite 
 * is greater than 0). Uses fastInterchange if fi 
 * is true, otherwise uses localSearch. The best 
 * solution found (after path-relinking) is returned 
 * in s.
 *  
 *****************************************************/

double grasp_target (
  PMSolution *s,   //output: will hold the best solution found
  PMConstructive *constructive, //input: used to run the appropriate constructive heuristic
  int it,          //input: number of iterations to be performed
  double target, 
  PMElite *elite,  //input/output: will hold best solutions found and be used for path-relinking (can be NULL)
  char *lsmethod,   //local search method to be used
  char *combmethod,
  int *flag,
  bool stats,      //input: print statistics to stdout?
  bool partials    //input: print partial information to stdout?
) {
	int verbose = stats ? 1 : 0;

	PMRecorder recorder;
	BossaTimer t(true);          //overall
	BossaTimer ls_timer(false);  //local search
	BossaTimer con_timer(false); //constructive
	BossaTimer comb_timer(false); //combination time

	PMInstance *inst = s->getInstance();
	PMElite *elite_best = new PMElite (inst, 1);
	// Gabriel - alteração 30/06/2008	- início 	
#ifdef MINERACAO
	PMEliteSimple *elite_simples;
	PMConstructiveMD constructive_md;
	elite_simples = new PMEliteSimple(inst, g_tam_elite_simples); // Elite de tamanho 4
	//constructive_md = new PMConstructiveMD;
	constructive_md.setSeed(g_seed);
	constructive_md.setMethod("sample:2");
	constructive_md.setTamEliteSimples(g_tam_elite_simples);
	constructive_md.setNumeroPadroes(g_npadroes);
	constructive_md.setSuporteMinimo(g_supmin);
#endif	
	// Gabriel - alteração 30/06/2008	- fim
	PMSolution *r  = new PMSolution (inst);
	PMSolution *rs = new PMSolution (inst); //relink solution
	
	int p = inst->getP();
	int n = inst->getN();
	int m = inst->getM();
	
	if (verbose > 1) fprintf (stderr, "GRASP: m=%d n=%d p=%d\n", m, n, p);
	if (verbose) fprintf (stderr, "Running grasp...\n");

        int i;
        double solcost;
	for(i = 1; i <= it; i++){
		//i++;
		//-------------------------------
		// build a solution from scratch
		//-------------------------------
		con_timer.resume();
		s->reset();
		// Gabriel - alteração 30/06/2008 - início 
#ifdef MINERACAO
		if (i <= (it/2)) 
			constructive->run(s);
		else
		{
			if (i == ((it/2) + 1))
			{
//				FILE *elite_file = fopen ("c:\\ProjetoFinal\\popstar2002\\Debug\\elite.txt", "w");
				FILE *elite_file = fopen ("elite.txt", "w");
				if (elite_file == NULL) 
					fatal ("elite_simples", "could not open input file");
				elite_simples->output(elite_file,true);
				fclose(elite_file); 					
			}
				constructive_md.run(s);
		}
#endif
#ifndef MINERACAO
		constructive->run(s);
#endif
		// Gabriel - alteração 30/06/2008 - fim  
		double original = s->calcCost();
		con_timer.pause();
                if(original <= target){
			solcost = original;
			fprintf (stdout, "solution %.2f\n", solcost);
                        (*flag) = 1;
               		break;
		}
		
		//----------------------
		// perform local search
		//----------------------
		ls_timer.resume();
		runLocalSearch (s, lsmethod);
		double result = s->calcCost();
		ls_timer.pause();
                if(result <= target){
			solcost = result;
			fprintf (stdout, "solution %.2f\n", solcost);
                        (*flag) = 1;
               		break;
		}

		if (verbose) fprintf (stderr, "%5d: [%.2f] %.2f -> %.2f", i, elite_best->getSolutionCost(1), original, result);

		unsigned int seed = BossaRandom::getNonZero();

		if (elite!=NULL) {
			comb_timer.resume();
			int d = elite->countConsistent();
			bool insert = true;

			//if there are elite solutions, relink
			if (d>=1) { 
				int e = elite->getRandomDifferent (s, 1, d); //a probably different solution
				
				if (e!=0) { //solution is not identical to some solution in the pool
					double elite_cost = elite->getSolution (e, r); //r is the solution
					double oldcost = (elite_cost < result) ? elite_cost : result; //get best of two original

					if (!r->equals(s)) { //if they are not equal, try to relink them
						//--------------------------------
						// combine both solutions into rs
						//--------------------------------
						double newcost = combine (rs, s, r, combmethod); //s is the one produced by local search, r is the elite one
						if (verbose) fprintf (stderr, " x %.2f -> %.2f", elite_cost, newcost);
						if (rs->computeDifference(r)>1 && rs->computeDifference(s)>1) {
							localSearch (rs, true, true); //try to improve by local search if necssary
							newcost = rs->calcCost();
							if (verbose) fprintf (stderr, " -> %.2f", newcost);
							if(newcost <= target){
								solcost = newcost;
								fprintf (stdout, "solution %.2f\n", solcost);
					                        (*flag) = 1;
	               						break;
							} //if(newcost <= target)
						} //if (rs->computeDifference(r)>1 && rs->computeDifference(s)>1)

						//-----------------------------------------------------------
						// if there is an improvement, try to save 'hybrid' solution
						//-----------------------------------------------------------
						if (newcost < oldcost - EPSILON) {
							elite->add (rs, newcost);      //add to the set of elite solutions
							elite_best->add (rs, newcost); //save it in the local list of elite solutions
							if (recorder.report (i, t.getTime(), newcost)) {
								if (partials) fprintf (stdout, "partial %.3f %.3f\n", recorder.best_time, recorder.best_value);
							} //if (recorder.report (i, t.getTime(), newcost))
						} //if (newcost < oldcost - EPSILON)


					} //if (!r->equals(s))
				} else insert = false;
			} //if (d>=1)
			if (insert) elite->add (s, result);
			comb_timer.pause();
		} //if (elite!=NULL)

		// Gabriel - alteração 30/06/2008 - início 
#ifdef MINERACAO
		if (i <= (it/2))
		{
			elite_simples->add(s, result);
			//elite_simples->outputScreen();
		} //if (i <= 16)
#endif
		// Gabriel - alteração 30/06/2008 - fim

		//------------------------
		// save original solution
		//------------------------
		elite_best->add(s, result); 
		if (recorder.report (i, t.getTime(), result)) {
			if (partials) {fprintf (stdout, "partial %.3f %.3f\n", recorder.best_time, recorder.best_value);}
		} //if (recorder.report (i, t.getTime(), result))
		
		if (verbose) {
			if (recorder.best_iteration==i) fprintf (stderr, " <");
			fprintf (stderr, "\n");
		} //if(verbose)
		BossaRandom::randomize(seed);		
	}
	
	//----------------------------------------------------
	// copy best solution to s (and its value to solcost)
	//----------------------------------------------------
	if((*flag) == 0) solcost = elite_best->getSolution (1, s);
	
	//-----------------------
	// deallocate structures
	//-----------------------
	delete elite_best;
	delete r;
	delete rs;

	//--------------
	// print report
	//--------------
	fprintf (stdout, "i %d\n", i);
	fprintf (stdout, "grasp target %.2f\n", target);
	fprintf (stdout, "grasp best solution %.2f\n", solcost); 
	fprintf (stdout, "grasptime %.2f\n", t.getTime());

	//if (stats) {
		//fprintf (stdout, "i %d\n", i);
		//fprintf (stdout, "grasp target %d\n", target); 
		//fprintf (stdout, "grasp %.2f\n", solcost);
		//fprintf (stdout, "grasptime %.2f\n", t.getTime());
		//fprintf (stdout, "graspconst ");
		//constructive->printMethod(stdout);
		//fprintf (stdout, "\n");
		//fprintf (stdout, "lsmethod %s\n", lsmethod);
		//fprintf (stdout, "graspcomb %s\n", combmethod);
		//fprintf (stdout, "grasplstimeavg %.2f\n", ls_timer.getTime()/(double)it);
		//fprintf (stdout, "graspconsttimeavg %.2f\n", con_timer.getTime()/(double)it);
		//fprintf (stdout, "grasplstime %.2f\n", ls_timer.getTime());
		//fprintf (stdout, "graspconsttime %.2f\n", con_timer.getTime());
		//fprintf (stdout, "graspbestit %d\n", recorder.best_iteration);
		//fprintf (stdout, "graspbesttime %.2f\n", recorder.best_time); 
	//}

	return solcost;
}





/************************************************
 *
 * Show usage: prints message explaining how the 
 * program should be used and exits
 *
 ************************************************/

void showUsage (char *argv[]) {
	fprintf (stderr, "\nUsage: %s <file> [options]\n", argv[0]);
	//fprintf (stderr, "    Valid options:\n");
	//fprintf (stderr, "    -grasp       : execute grasp\n");
	fprintf (stderr, "    -graspit <n> : number of grasp iterations [default=32]\n");
	fprintf (stderr, "    -elite <x>   : size of the pool of elite solutions [default=10]\n");
	//fprintf (stderr, "  -greedy      : execute greedy algorithm\n"); 
	fprintf (stderr, "    -p <x>       : sets number of facilities (p) to x; overrides value given in the input file, if any\n");
	fprintf (stderr, "    -output <x>  : output solution to file 'x'\n");
	//fprintf (stderr, "  -ls <met>    : use <met> as the local search method (default is m:-1)\n");
	//fprintf (stderr, "  -ch <met>    : use <met> as the constructive heuristic (default is rpg:1)");
	//fprintf (stderr, "  -vnds      : run vnds\n");
	fprintf (stderr, "\n");
	exit (-1);
}



/*******************************
 *
 * test constructive heuristics
 *
 *******************************/

void testConstructive (PMInstance *inst, PMConstructive *tc, PMConfig *config) {
	
	int i,j,k,p;
	k = config->tc; //number or iterations
	double *v, *lit;

	p = inst->getP();

	PMSolution **s;
	
	s = new PMSolution* [k+1];
	v = new double [k+1];
	lit = new double [k+1];
	for (i=1; i<=k; i++) s[i] = new PMSolution (inst);

	fprintf (stderr, "tcmethod ");
	tc->printMethod(stderr);
	fprintf (stderr, "\n");


	bool run_ls = true;

	BossaTimer t;
	t.start();
	for (i=1; i<=k; i++) {
		BossaRandom::randomize(config->seed - 1 + i);
		tc->run(s[i]);

		lit[i] = 0.0;
		if (run_ls) {
			int it;
			localSearch (s[i], true, true, &it);
			lit[i] = (double)it;
		}
		v[i] = s[i]->calcCost();
	}
	fprintf (stdout, "tctime %.2f\n", t.getTime() / (double)k);
	fprintf (stdout, "tccount %.2f\n", (double)k);
	fprintf (stdout, "tcmin %.2f\n", RFWStats::min(v,1,k));
	fprintf (stdout, "tcavg %.2f\n", RFWStats::average(v,1,k));
	fprintf (stdout, "tcstddev %.2f\n", RFWStats::stddev(v,1,k));
	fprintf (stdout, "tcmed %.2f\n", RFWStats::median(v,1,k));
	fprintf (stdout, "tcmax %.2f\n", RFWStats::max(v,1,k));

	fprintf (stdout, "tcminit %.2f\n", RFWStats::min(lit,1,k));
	fprintf (stdout, "tcavgit %.2f\n", RFWStats::average(lit,1,k));
	fprintf (stdout, "tcstddevit %.2f\n", RFWStats::stddev(lit,1,k));
	fprintf (stdout, "tcmedit %.2f\n", RFWStats::median(lit,1,k));
	fprintf (stdout, "tcmaxit %.2f\n", RFWStats::max(lit,1,k));



	int npairs = k * (k-1) / 2;
	double *diff = new double [npairs+1];

	int c = 0;
	for (i=1; i<k; i++) {
		for (j=i+1; j<=k; j++) {
			c++;
			diff[c] = (double) (s[i]->computeDifference(s[j]));
		}
	}

	fprintf (stdout, "tcmindiff %.3f\n", RFWStats::min(diff,1,npairs) / (double)p);
	fprintf (stdout, "tcavgdiff %.3f\n", RFWStats::average(diff,1,npairs)/ (double)p);
	fprintf (stdout, "tcmeddiff %.3f\n", RFWStats::median(diff,1,npairs)/ (double)p);
	fprintf (stdout, "tcmaxdiff %.3f\n", RFWStats::max(diff,1,npairs)/ (double)p);
	
	for (i=1; i<=k; i++) delete s[i];
	delete [] s;
	delete [] v;
}



void testPathRelinking (PMInstance *inst, PMConstructive *c, PMElite *elite, const char *method) {
	if (elite==NULL) fatal ("testPathRelinking", "number of elite solutions cannot be zero");
	fprintf (stderr, "Resetting list of elite solutions...\n" );
	PMElite *single_elite = new PMElite (inst, 1);
	
	elite->reset();
	fprintf (stderr, "Creating new solution...\n");
	PMSolution *s1 = new PMSolution (inst);
	PMSolution *s2 = new PMSolution (inst);
	PMSolution *s  = new PMSolution (inst);
	PMSearchPathRelinking *sm = new PMSearchPathRelinking (inst);


	int i, j, t;
	int e = elite->getCapacity();
	int maxtries = 2 * e;

	//----------------------------------
	// fill the pool of elite solutions
	//----------------------------------
	fprintf (stderr, "Entering main loop...\n");
	for (t=1; t<=maxtries; t++) {
		fprintf (stderr, "%d ", t);
		s1->reset();
		c->run(s1);                         //constructive algorithm
		
		//localSearch(s1, true, NULL, false); //local search
		localSearch(s1, true, true, false, false); 
		elite->add(s1, s1->calcCost());
		fprintf (stderr, "%d/%d\n", t, elite->countConsistent());
		if (elite->countConsistent()>=e) break;
	}

	enum {DOWN, HOR, UP};

	int m = inst->getM(); //number of candidate facilities
	
	
	double *imp  [3];
	int    *count[3]; 
	double *sum  [3];

	for (i=0; i<3; i++) {
		count[i] = new int [m+1];
		imp  [i] = new double [m+1];
		sum  [i] = new double [m+1];
		for (j=0; j<=m; j++) {
			count[i][j] = 0;
			imp  [i][j] = 0;
			sum  [i][j] = 0;
		}
	}

	int ns = elite->countConsistent();

	int eliteplus = 10;

	int bestuniform = 0; //number of times uniform won
	int bestbiased  = 0; //number of times biased won
	int ties        = 0; //number of times there was a tie
	double sumuniform = 0;
	double sumbiased  = 0;

	for (int i2=eliteplus+1; i2<=ns; i2++) {
		int totalcount = 0;
		int rightcount = 0;
		int totaldiff = 0;
		int rightdiff = 0;
		for (int i1=1; i1<=eliteplus; i1++) {

			if (i1==i2) continue;
			double v1 = elite->getSolution (i1, s1);
			double v2 = elite->getSolution (i2, s2);
			int diff = s1->computeDifference(s2);

			int dir = HOR;
			if (v1<v2) dir = UP;
			else if (v1 > v2) dir = DOWN;

			combine (s, s1, s2, method);

			/*
			sm->setSolutions (s1, s2);           //define the extremes of the path
			single_elite->reset();               //the best solution in the path will be here
			s->copyFrom (s1);              
			search (s, true, sm, single_elite, true);  //excecute path-relinking
			*/
			//double vsol = single_elite->getSolution (1, s);   //copy solution to s and get value
			double vsol = s->calcCost();
			if (s->computeDifference(s1)>1 && s->computeDifference(s2)>1) {
				localSearch (s, true, true);  //try to improve solution with local search
				vsol = s->calcCost();   //get the value
			}

			//calculate profit
			double profit;
			if (v1<v2) profit = v1 - vsol;
			else profit = v2 - vsol;
			if (profit<0) profit = 0;

			//bool success = (profit > 0);

			fprintf (stderr, "[%d:%d] %f x %f (%d,%d) -> %f (%f)\n", i1, i2, v1, v2, dir, diff, vsol, profit);

			count[dir][diff] ++;
			imp  [dir][diff] += profit;
			sum  [dir][diff] += v1+v2;

			totalcount ++;
			totaldiff  += diff;

			if (profit) {
				rightcount ++;
				rightdiff += diff;
			}

		}
		double probcount = (double)rightcount / (double)totalcount;
		double probdiff  = (double)rightdiff  / (double)totaldiff;

		double ratio;
		if (probcount>0) ratio = probdiff / probcount;
		else ratio = 1;
		
		fprintf (stdout, "tpd %d %f %f %f\n", i2, probdiff, probcount, ratio);
		fprintf (stderr, "tpd %d %f %f %f\n", i2, probdiff, probcount, ratio);

		double mindiff = 0.001;
		if (ratio > 1 + mindiff) {
			bestbiased ++;	
		} else if (ratio < 1 - mindiff) {
			bestuniform ++;
		} else {
			ties ++;
		}

		sumbiased += probdiff;
		sumuniform += probcount;
	}


	{
		int x = bestbiased + bestuniform + ties;
		fprintf (stdout, "trt %f\n", (double)ties/(double)x);
		fprintf (stdout, "tru %f\n", (double)bestuniform/(double)x);
		fprintf (stdout, "trb %f\n", (double)bestbiased/(double)x);
		fprintf (stdout, "trpu %f\n", (double)sumuniform/(double)x);
		fprintf (stdout, "trpb %f\n", (double)sumbiased/(double)x);
		fprintf (stdout, "trpr %f\n", (double)sumbiased/(double)sumuniform);


	}

	for (i=0; i<3; i++) {
		for (j=0; j<=m; j++) {
			if (count[i][j]==0) continue;
			fprintf (stdout, "tpr %d %d %d %f\n", i, j, count[i][j], imp[i][j]); 
		}
	}

	for (j=0; j<=m; j++) {
		int c = count[0][j] + count[1][j] + count[2][j];
		double v = imp[0][j] + imp[1][j] + imp[2][j];
		double a = (sum[0][j] + sum[1][j] + sum[2][j]) / (2.0 * (double)c);
		double im = v / (double)c;
		if (c>0) {
			fprintf (stdout, "tpa %d %d %f %f %f \n", j, c, a, im, a-im);
		}
	}

	delete single_elite;
	delete sm;
	delete s1;
	delete s2;
	delete s;
	for (i=0; i<3; i++) {
		delete [] count[i];
		delete [] imp[i];
		delete [] sum[i];
	}
}




/***********************************************
 *
 * run the local search method specified
 * by lsmethod; changes the original solution s
 *
 ***********************************************/

void runLocalSearch (PMSolution *s, char *lsmethod) {
	int i=0;
	while (true) {
		switch (lsmethod[i]) {
			case '\0': return;
			case 'n': break;                                               //none
			case 'c': clusterSearch (s); break;
			case 'f': fastInterchange (s); break;                          //fast interchange
			case 'm': localSearch (s, true, false, NULL, false); break;    //matrix
			case 's': localSearch (s, true, true, NULL, false); break;     //sparse matrix
			default: fatal ("runLocalSearch", "invalid local search method");
		}
		i++;
	}
}




/****************************************************************
 *
 * procedure to test (time) local search methods. First, 
 * builds a solution using the method specified by constr. Then,
 * executes local search on the resulting solution; if the 
 * minimum time (given mt), restores the original solution and
 * does it again --- until mt is reached. Prints the average 
 * time to run everything.
 *
 ****************************************************************/

void testLocalSearch (PMInstance *inst, PMConstructive *constr, char *method, double mt) {
	BossaTimer timer(1);
	PMSolution *greedy = new PMSolution (inst);
	constr->run (greedy);
	fprintf (stdout, "lsconsttime %.2f\n", timer.getTime());
	fprintf (stdout, "lsconstmethod ");
	constr->printMethod(stdout);
	fprintf (stdout, "\n");
	fprintf (stdout, "lstestmaxtime %f\n", mt);
	fprintf (stdout, "lstestmethod %s\n", method);
	fprintf (stdout, "lstestcf %f\n", PMDistanceOracle::cache_factor);
	fprintf (stdout, "lstestmethodcf %s:%.1f\n", method, PMDistanceOracle::cache_factor);

	PMSolution *temp = new PMSolution (inst);
	
	double t;
	timer.start();
	int runs = 0;
		
	while (1) {
		temp->copyFrom(greedy);
		runLocalSearch (temp, method);
		t = timer.getTime();
		runs++;
		if (t > mt) break;
	}

	fprintf (stdout, "lstestconstr %.2f\n", greedy->calcCost());
	fprintf (stdout, "lstestvalue %.2f\n", temp->calcCost());
	fprintf (stdout, "lstestruns %d\n", runs);
	fprintf (stdout, "lstesttime %.6f\n", t / (double)runs);

	delete greedy;
	delete temp;
}





/****************
 *
 * Main program
 *
 ****************/

int main (int argc, char *argv[]) {
	//testPartition ();

	PMInstance *inst;
	BossaTimer t;
	PMConfig c;

	//-------------------------
	// read input parameters
	//-------------------------	
	if (argc<2) showUsage (argv);
	bool success = c.readParameters (argc, argv);
	if (!success) showUsage (argv);	

	PMDistanceOracle::cache_factor = c.cf;
// Gabriel - alteração 07/09/2008 - início - parametro inútil
	//fprintf (stdout, "cf %.3f\n", PMDistanceOracle::cache_factor);
// Gabriel - alteração 07/09/2008 - fim

	//-------------------------
	// read instance
	//-------------------------
	t.start();
	inst = readInstance (argv[1], c.pvalue);
	if (inst->getP() <= 0) fatal ("main", "invalid number of centers");

	fprintf (stdout, "inittime %.2f\n", t.getTime());
	fprintf (stdout, "oracletime %.3f\n", inst->getOracleTime());


	
	//---------------------------
	// print instance attributes
	//---------------------------
// Gabriel - alteração 07/09/2008 - início
#ifdef MINERACAO
	g_seed = c.seed;
	g_tam_elite_simples = c.elite_simples;
	g_npadroes = c.npadroes;
	g_supmin = c.supmin;
#endif
// Gabriel - alteração 07/09/2008 - fim
	BossaRandom::randomize(c.seed);

	char *filename = argv[1];
	char *instance = new char [strlen(filename)+1];
	strcpy (instance, filename);
	RFW::stripPath (instance);
	RFW::stripExtension (instance);

	fprintf (stdout, "file %s\n", filename);
	fprintf (stdout, "instance %s\n", instance);
	fprintf (stdout, "instancep %s-%d\n", instance, inst->getP());
	fprintf (stdout, "vertices %d\n", inst->getN());
	fprintf (stdout, "centers %d\n", inst->getP());
	fprintf (stdout, "seed %d\n", c.seed);
	fprintf (stdout, "target %f\n", c.target);

	//----------------------------------------
	// create a variable to hold the solution
	//----------------------------------------
	PMSolution *s = new PMSolution(inst);      //this variable will hold the solution
	PMConstructive::addRandom (s, inst->getP());	//starts with a random solutoin

	//---------------------------------------
	// set constructive algorithm parameters
	//---------------------------------------
	PMConstructive constructive;
	constructive.setMethod (c.ch);


	//--------------------------------
	// create pool of elite solutions
	//--------------------------------
	PMElite *elite = (c.elite>0) ? new PMElite (inst, c.elite) : NULL; 
	
	//--------------
	// start timing
	//--------------
	t.start();
	
	//-----------------------------------------------------
	// run grasp (updates s and populates elite)
	//-----------------------------------------------------
	// Isabel - 30/08/2010
        int flag = 0;
	if (c.run_grasp) grasp_target (s, &constructive, c.graspit, c.target, elite, c.ls, c.cmg, &flag, true, c.partials);
	
        /*
	//------------------------------------------
	// run rvns (updates s and populates elite)
	//------------------------------------------
	if (c.run_rvns) {vns (s, c.rvns_k, c.rvns_r, true, elite, true);}
	
	if (c.run_vns) {
		if (c.vns_k <= 0) vns (s, inst->getP(), -1, false, elite, true);	
		else vns (s, c.vns_k, -1, false, elite, true);
	}


	if (c.run_vnds) vnds (s, elite, c.vnds_k, c.vnds_b, true);

	//-------------------------------------
	// run tabu search with path-relinking
	//-------------------------------------
	if (c.run_tabu) tabu (s, c.tabuit, 
				(int)floor(0.5+inst->getP()*c.itenure), 
				(int)floor(0.5+inst->getP()*c.rtenure), 
				c.elite
			);

	*/
	//--------------------
	// run path-relinking
	//--------------------
        fprintf (stdout, "flag = %d\n", flag);
	if (elite!=NULL && !c.test_pr && c.postopt &&!flag) {
		pathRelinking_target (s, elite, c.cmp, c.target, true);
		fprintf (stdout, "postopt 1\n");
		delete elite;
	} else fprintf (stdout, "postopt 0\n");

	fprintf (stdout, "usepr %d\n", (elite!=NULL) ? 1 : 0);

	//--------------------------
	// test individual elements
	//--------------------------
	if (c.tc>0) {testConstructive (inst, &constructive, &c);}
	if (c.tls)  {testLocalSearch (inst, &constructive, c.ls, c.tlstime);}
	if (c.test_pr) {testPathRelinking (inst, &constructive, elite, c.cmg);}


	//------------------------------
	// output best solution to file 
	//------------------------------
	if (c.output) {
		t.pause();
		FILE *outfile = fopen (c.outname, "w");
		if (outfile == NULL) fprintf (stderr, "Could not open %s for writing.\n", c.outname);
		else {
			s->printSolution(outfile, instance);
			fclose (outfile);
		}
		t.resume();
	}
	
	/*
	fprintf (stdout, "totaltime %.2f\n", totaltime);
	totaltime += inst->getOracleTime();
	fprintf (stdout, "totaltimeo %.2f\n", totaltime);
	totaltime += inst->getFloydTime();
	fprintf (stdout, "totaltimeof %.2f\n", totaltime);	*/

	double totaltime = t.getTime();	fprintf (stdout, "cputime %.2f\n", totaltime);
	fprintf (stdout, "bestsol %.2f\n", s->calcCost());

	//----------
	// clean up
	//----------
	delete [] instance;
	delete s;	
	delete inst;

	return 0;
}

/*

void move (PMSolution *s, int in, int *closest, int *closest2, int *out, double *w) {
	int i;
	int p = s->getP();
	int n = s->getN();
	PMInstance *inst = s->getInstance();
	if (s->contains(in)) fatal ("move", "invalid insertion facility");
	
	double *v = new double [n+1];

	//initialization
	*w = 0;  //change in the objective function value obtained 
		 //by the best interchange; always negative
	for (i=1; i<=p; i++) v[s->getFacility(i)] = 0;
	for (i=1; i<=n; i++) v[i] = 0;
	
	//best deletion
	for (i=1; i<=n; i++) {
		double din = inst->getDist (i, in);
		double d1  = inst->getDist (i, closest[i]);
		double d2  = inst->getDist (i, closest2[i]);
		
		if (din < d1) {   //WARNING: THIS LINE WAS WRONG IN THE PAPER (originally, din < d2)
			*w += (din - d1); //in will be the new closest, we'll save something
		} else {
			//if we delete closest[i], we'll save d1, but we'll 
			//have to spend either din or d2
			double min = (din<d2) ? din : d2;
			v[closest[i]] += (min-d1); //we'll have to pay for the difference 
		}
	}	

	//find best facility to remove
	int best = s->getFacility(1);
	for (i=1; i<=p; i++) {
		if (v[s->getFacility(i)]<v[best]) best = s->getFacility(i);
	}

	//fprintf (stderr, " (o:%d i:%d w:%8.2f v:%8.2f)", *out, in, *w, v[best]);
	*out = best;
	*w += v[best];

	delete [] v;
}
*/

/*
double solve2(PMSolution *s, int *fac1, int *fac2, int *assign, double *acost) {
	PMInstance *inst = s->getInstance();
	int n = inst->getN();
	int m = inst->getM();

	int f1, f2, c;
	double bestprofit = -INFINITY;
	double baseprofit = 0;
	
	int best1, best2;
	best1 = best2 = 0;

	for (f1=1; f1<m; f1++) {
		for (f2 = f1+1; f2<=m; f2++) {
			double profit = 0;
			for (c=1; c<=n; c++) {
				double d1 = inst->getDist(c,f1);
				double d2 = inst->getDist(c,f2);
				double d = (d1<d2) ? d1 : d2;
				if (d<acost[c]) profit += (acost[c]-d);
			}	
			if (profit > bestprofit + EPSILON) {
				if (best1 == 0) baseprofit = profit; 
				//fprintf (stderr, "{%d,%d:%8.2f}\n", f1, f2, profit - baseprofit);
				best1 = f1;
				best2 = f2;
				bestprofit = profit;
			}
		}
	}
	*fac1 = best1;
	*fac2 = best2;
	return (bestprofit - baseprofit);	
}*/



