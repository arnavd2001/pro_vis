#ifndef _CUDA_HEADER_
#define _CUDA_HEADER_

/** \file CUDA_header.h CUDA routines for calculating collisions and contacts among a beads. */

/** A "key" that can be used to fetch results of non-blocking computations sent to the GPU. */
struct CollisionCountPromise {  // Vectors of # of collisions
	int *d_toReduce;
	int *d_reduced;
};

/** Holds a triple of floats. */
typedef struct {
	float x;
	float y;
	float z;
} ElfFloat3d;

/** Holds a triple of integer numbers. */
typedef struct {
	int x;
	int y;
	int z;
} Elfnumtrd;

/*
Launches the GPU procedure for counting collisions in 'vector' which has size 'size'.

This function does not wait until the GPU procedure is finished.
It returns a "Promise" structure which represents a promise for a future return value.
The return value can be fetched with the _fetch corresponding function.
*/
struct CollisionCountPromise
	count_collisions_launch(ElfFloat3d *vector, int size);

/*
Returns the number of collisions associated with the given "Promise" structure.

The "Promise" structure is a promise for a future return value, which is returned
  by the non-blocking _launch function.
*/
int count_collisions_fetch(struct CollisionCountPromise promise);


// The functions below have similar behavior as the above ones.
struct CollisionCountPromise
	count_contacts_launch(ElfFloat3d *vector, int size);

int count_contacts_fetch(struct CollisionCountPromise promise);

#endif
