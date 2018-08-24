#include "elite_simples.h"

void PMEliteSimple::fatal (const char *func, const char *msg) {
	fprintf (stderr, "PMEliteSimple::%s: %s.\n", func, msg);
	exit (-1);
}

		
/***********************************************************
 *
 * promote an instance from the i-th position to its proper
 * position (tries only positions below i)
 *
 ***********************************************************/

void PMEliteSimple::promote(int i) {
	while (i>1) {
		if (v[i-1] <= v[i]) break;
		PMSolution *temps;
		double tempv;
		temps = s[i]; tempv = v[i];
		s[i] = s[i-1]; v[i] = v[i-1];
		s[i-1] = temps; v[i-1] = tempv;
		i--;
	}			
}

/*******************************************
 *
 * Sort list of elite solutions by solution
 * value. The algorithm used is insertion
 * sort.
 *
 *******************************************/

void PMEliteSimple::sort() {
	for (int i=2; i<=consistent; i++) promote(i);
}


/************************************
 *
 * Copy into 'sol' the i-th solution
 *
 ************************************/

double PMEliteSimple::getSolution (int i, PMSolution *sol) {
	sol->copyFrom (s[i]); 
	return v[i];
}


/**********************************************
 *
 * returns the number of the solution that 
 * is most similar to t; returns 0 if there is
 * a solution that is too similar
 *
 **********************************************/

int PMEliteSimple::getMostSimilar (PMSolution *t, int p1, int p2) {
	int diff, mindiff, sub;
	mindiff = t->getP() + 1;
	sub = 0;

	if (p1>p2) fatal ("getMostSimilar", "invalid range");

	for (int i=p1; i<=p2; i++) {
		diff = t->computeDifference(s[i]);
		//if (diff<=2) return 0; //there's not point in adding solution is too 
		if (diff==0) fatal ("getMostSimilar", "solution is not unique");
		if (diff<mindiff) {
			mindiff = diff;
			sub = i;
		}
	}
	return sub;
}


/************************************************
 *
 * checks if t is unique: returns true if no
 * other solution is equal to t, false otherwise
 *
 ************************************************/
	
bool PMEliteSimple::isUnique(PMSolution *t) {
	int c = countConsistent();
	bool unique = true;
	for (int i=1; i<=c; i++) {
		if (t->equals(s[i])) {unique=false; break;}
	}
	return unique;		
}


/************************************************
 *
 * select an element in the pool at random to
 * be combined with s; the probability of a
 * given element to be selected is proportional
 * to how different it is from s; returns 0 if
 * a near-identical solution is found
 *
 ************************************************/

int PMEliteSimple::getRandomDifferent(PMSolution *t, int p1, int p2) {
	int diff, sumdiff, sub;
	sub = diff = sumdiff = 0;
			
	int p = t->getP();

	//check if the solution is identical to some other solution
	for (int i=p1; i<=p2; i++) {
		if (s[i]->getP() != p) fatal ("getRandomDifferent", "invalid solution in pool");

		//if very similar to an existing solution, there's a problem
		diff = t->computeDifference(s[i]);
		if (diff == 0) return 0;

		//among those that cost more, select the most similar
		sumdiff += diff;
		if (BossaRandom::getInteger(1,sumdiff) <= diff) sub = i;
	}
	return sub;
}


/*************************************
 *
 * output elite solutions to files
 *
 *************************************/

void PMEliteSimple::outputToFiles (const char *prefix) {
	FILE *file;
	char name [256];
	PMSolution *x;

	for (int i=1; i<=cap; i++) {
		x = s[i];
		sprintf (name, "%s%02d.out", prefix, i);
		fprintf (stderr, "Writing file %s... ", name);
		file = fopen (name, "w");
		x->printSolution(file, NULL);
		fclose (file);
		fprintf (stderr, "done.");
	}
}


/*****************************************************************
 *
 * for every pair of instances, output their symmetric difference
 *
 *****************************************************************/

void PMEliteSimple::outputDifferences (FILE *file) {
	for (int i=1; i<=cap; i++) {
		for (int j=1; j<=cap; j++) {
			fprintf (file, "(%d:%f x %d:%f): %5d\n", i, v[i], j, v[j], s[i]->computeDifference(s[j]));
		}
	}
}



void PMEliteSimple::output (FILE *file, bool complete) {
//	fprintf (file, "(%d/%d) ", consistent, cap);
//	for (int i=1; i<=cap; i++) {
//		if (complete) {
//			fprintf (file, "%d: %f ->", i, v[i]);
//			PMSolution *t = s[i];
//			for (int f=1; f<=t->getP(); f++) {
//				fprintf (file, " %d", t->getFacility(f));
//			}
//			fprintf (file, "\n");
//		} else {
//			fprintf (file, "%f ", v[i]);
//		}
//	}
//	if (!complete) fprintf (file, "\n");

	for (int i=1; i<=cap; i++) {
			PMSolution *t = s[i];
			for (int f=1; f<=t->getP(); f++) {
				fprintf (file, "%d ", t->getFacility(f));
			}
			fprintf (file, "\n");
	}
}		

void PMEliteSimple::outputScreen () {
//	fprintf (file, "(%d/%d) ", consistent, cap);
//	for (int i=1; i<=cap; i++) {
//		if (complete) {
//			fprintf (file, "%d: %f ->", i, v[i]);
//			PMSolution *t = s[i];
//			for (int f=1; f<=t->getP(); f++) {
//				fprintf (file, " %d", t->getFacility(f));
//			}
//			fprintf (file, "\n");
//		} else {
//			fprintf (file, "%f ", v[i]);
//		}
//	}
//	if (!complete) fprintf (file, "\n");

	for (int i=1; i<=cap; i++) {
			printf ("\n");
			PMSolution *t = s[i];
			for (int f=1; f<=t->getP(); f++) {
				printf ("%d ", t->getFacility(f));
			}
			printf ("\n");
	}
}		

/**************************************************
 *
 * compute difference from every valid solution in
 * the pool to t; returns the results in diff
 *
 **************************************************/

void PMEliteSimple::computeDifferences (int *diff, PMSolution *t) {
	int i;
	for (i=1; i<=consistent; i++) {diff[i] = t->computeDifference(s[i]);}
	for (i=consistent+1; i<=cap; i++) {diff[i] = t->getM()+1;
	}
}



/*******************************************************
 *
 * Tries to add a solution to the pool; the following 
 * criteria must be observed fro the insertion to be
 * made:
 * - solution must have distance at least 3 to every
 *   solution already in the pool
 * - 
 *
 *
 *******************************************************/

bool PMEliteSimple::add (PMSolution *t, double c) { 
	bool debug = true;
	//bool verbose = false;
	// Gabriel - alteração 30/06/2008	- início
	// Removendo restrições de adição do conjunto elite	
	// const int THRESHOLD = 4;
	const int THRESHOLD = -1;
	// Gabriel - alteração 30/06/2008	- Fim 
	if (c>v[cap]) return false; //not cheap enough

	int bestpos  = 0;
	int bestdiff = 2 * t->getM() + THRESHOLD + 1; //a big number

	for (int i=1; i<=cap; i++) {
		if (i>consistent) { //not empty yet
			if (bestdiff>=THRESHOLD) bestpos=i;  //add in last position if there's nobody very similar
			break;
		} 
		int diff = t->computeDifference(s[i]);
		if (v[i]<c) {
			if (diff < THRESHOLD) {bestpos=0; break;} //too similar --- don't add
		} else {
			if (diff<bestdiff) {
				bestpos=i;
				bestdiff=diff;
			}
		}
	}

	if (bestpos>0) {
		if (v[bestpos]<c) fatal ("add", "invalid position for insertion");
		s[bestpos]->copyFrom(t);
		if (bestpos>consistent) consistent = bestpos;
		v[bestpos] = c;
		if (debug) checkConsistency (bestpos);
		sort(); //this is insertion sort, so it takes O(cap) time to update
	} 

	return (bestpos>0);
}


/********************************************
 *
 * get solution that differs the most from t
 *
 ********************************************/

int PMEliteSimple::getMostDifferent (PMSolution *t) {
	int c = countConsistent();
	if (c < 1) fatal ("getMostDifferent", "pool is empty");
			
	int besti = 1;
	int bestdiff = t->computeDifference(s[1]);
	for (int i=2; i<=c; i++) {
		int diff = t->computeDifference(s[i]);
		if (diff > bestdiff) {
			besti = i;
			bestdiff = diff;
		}
	}
	return besti;
}


/**************
 *
 * Constructor
 *
 **************/

PMEliteSimple::PMEliteSimple (PMInstance *_inst, int _cap) {
	cap = _cap;
	inst = _inst;

	s = new PMSolution* [cap+1];
	v = new double [cap+1];

	for (int i=1; i<=cap; i++) {
		s[i] = new PMSolution (inst);
	}
	reset();
}


/**************************
 *
 * reset list of solutions
 *
 **************************/

void PMEliteSimple::reset() {
	consistent = 0;
	for (int i=1; i<=cap; i++) {
		s[i]->reset();
		v[i] = s[i]->calcCost(); //this should be a very big number
	}
}


/***********************
 *
 * destructor
 *
 ***********************/

PMEliteSimple::~PMEliteSimple() {
	for (int i=1; i<=cap; i++) delete s[i];
	delete [] s;
	delete [] v;
}


/************************************
 *
 * check the consistency of the node
 *
 ************************************/

void PMEliteSimple::checkConsistency (int i) {
	double diff = v[i] - s[i]->calcCost();
	if (diff>EPSILON || diff<-EPSILON) {
		fprintf (stderr, "(saved:%f calculated:%f)\n", v[i], s[i]->calcCost());
		fatal ("checkConsistency", "inconsistent solution");
	}
}

