#ifndef ELITESIMPLE_H
#define ELITESIMPLE_H
#include "solution.h"
#include "basics.h"
#include "bossa_random.h"

class PMEliteSimple {
	private:
		int cap;          //total number of solutions allowed
		PMSolution **s;   //list of solutions
		double *v;        //list of solution values
		PMInstance *inst; //instance to which the solutions... are solutions
		int consistent;   //number of consistent solutions

		void fatal (const char *func, const char *msg);
		void promote(int i);
		void sort();
		int getMostSimilar (PMSolution *t, int p1, int p2);	
		bool isUnique(PMSolution *t);

	public:
		PMEliteSimple (PMInstance *_inst, int _cap);
		int getRandomDifferent(PMSolution *t, int p1, int p2);
		void outputToFiles (const char *prefix);
		void outputDifferences (FILE *file);
		void output (FILE *file, bool complete = false);
		double getSolution (int i, PMSolution *sol);
		bool add (PMSolution *t, double c); //, int pos=0);
		void computeDifferences (int *diff, PMSolution *t);
		int getMostDifferent (PMSolution *t);
		void reset();
		void checkConsistency (int i);
		~PMEliteSimple();
// Gabriel - alteração 07/09/2008 - início
		void outputScreen();
// Gabriel - alteração 07/09/2008 - fim

		inline int count () {return cap;}
		inline int countConsistent () {return consistent;}
		inline PMInstance* getInstance() {return inst;}
		inline int getCapacity() {return cap;}
		inline double getSolutionCost (int i) {return v[i];}

};

#endif