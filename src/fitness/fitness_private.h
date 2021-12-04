#ifndef _FITNESS_PRIVATE_H
#define _FITNESS_PRIVATE_H

/** \file fitness_private.h Private header file solely for usage by fitness files. */

/* Header to be included just by files in the fitness/ directory
 */
#include <chaininghp.h>
#include <numtrd.h>

/**********************************
 *    FitnessCalc Procedures      *
 **********************************/

#define MAX_MEMORY ((long int) 4*1E9) // Max total size of memory allocated

/** Structure that holds resources to be reused throughout calls to functions. */
typedef struct FitnessCalc_ {
	const HPElem * chaininghp;
	int hpSize;
	void *space3d;
	int axisSize;
	double maxGyration;
} FitnessCalc;

/** Holds a triple of double values. */
typedef struct {
	double x;
	double y;
	double z;
} DPoint;

/** Holds a pair of double values. */
typedef struct {
	double first;
	double second;
} DPair;

/** Holds a pair of integer values. */
typedef struct {
	int first;
	int second;
} IPair;

/** Holds a triple of integer values. */
typedef struct {
	int first;
	int second;
	int third;
} ITriple;

/** Used to return many measures to the caller. */
typedef struct {
	int hh, pp, hp, hb, pb, bb;
	int collisions;
} BeadMeasures;

FitnessCalc FitnessCalc_get(); // Returns the FIT_BUNDLE of the protein being assessed.
BeadMeasures proteinMeasures(const numtrd *BBbeads, const numtrd *SCbeads, const HPElem *chaininghp, int hpSize);

#endif
