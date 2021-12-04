#ifndef FITNESS_H
#define FITNESS_H

/** \file fitness.h Public header file with routines for calculating the fitness of a protein represented as a chain of relative movements + hp chain. */

#include <numtrd.h>
#include <chaininghp.h>
#include <migrch.h>
#include <config.h>

/* Initialize the resources needed for calling the functions in this library.
 *   of this library can be called in parallel among threads.
 */
void FitnessCalc_initialize(const HPElem * chaininghp, int hpSize);

/* Cleans up resources allocated.
 */
void FitnessCalc_cleanup();


/* Returns the fitness for a protein already registered with FitnessCalc_initialize,
 *   considering that the protein has its 3d coordinates in coordsBB and coordsSC.
 */
double FitnessCalc_run(const numtrd *coordsBB, const numtrd *coordsSC);

/* Returns the fitness for a protein already registered with FitnessCalc_initialize,
 *   considering that the protein has movement chain 'chain'.
 */
double FitnessCalc_run2(const shiftmel * chain);

/* Returns measures for a given movement chain.
 * chain    - the movement chain from which to extract measures
 *
 * For all the other arguments, if they are NULL they are ignored, if they aren't null, the memory pointed to receives the desired value.
 * Hcontacts_p  - the number of hidrophobic contacts
 * collisions_p - the number of collisions among beads
 * bbGyration_p - the gyration radius for the backbone beads
 */
void FitnessCalc_measures(const shiftmel *chain, int *Hcontacts_p, int *collisions_p, double *bbGyration_p);

#endif // FITNESS_H

