#ifndef ABC_ALG_H
#define ABC_ALG_H

/** \file abc_alg.h Public header file with routine for performing prediction of a protein structure. */

#include <migrch.h>
#include <chaininghp.h>
#include <config.h>

#include <solution/solution.h>

/** Structure for returning prediction results to the user. */
typedef struct PredResults_ {
	double fitness;    /**< Fitness of the predicted protein */
	int contactsH;     /**< Number of H contacts */
	int collisions;    /**< Number of collisions among beads */
	double bbGyration; /**< Gyration radius for the backbone beads */
} PredResults;

/** Given a protein in the HPElem * format, searches the 3D conformation with minimal energy.
 * 'nCycles' is the number of cycles desired for the algorithm to run.
 *
 * Returns a shiftmel *, which is a sequence of movements of the backbone of the protein,
 *   and also the movements of the side-chain beads relative to the backbone.
 *
 * If 'results' is not NULL, it receives additional values about the predicted
 *   protein (see the structure itself for what values are given).
 */
Solution ABC_predict_structure(const HPElem * chaininghp, int hpSize, int nCycles, PredResults *results);

#endif // ABC_ALG_H
