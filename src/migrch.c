#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "migrch.h"
#include "shiftmel.h"
#include "numtrd.h"
#include "random.h"

void migrch_set_element(shiftmel * chain, int eleIdx, unsigned char bb, unsigned char sc){
	chain[eleIdx] = shiftmel_make(bb, sc);
}

shiftmel * migrch_create(int size){
	shiftmel * chain = malloc(sizeof(shiftmel) * size);
	int i;
	for(i = 0; i < size; i++)
		migrch_set_element(chain, i, FRONT, RIGHT);

	return chain;
}

// Given a predecessor vector and a movement, applies such
//   movement in the predecessor vector and returns the result.
static inline
numtrd getNext(numtrd pred, unsigned int movement){
	int *first, *second;
	numtrd result;

	result = numtrd_make(0, 0, 0);

	// Get first and second filled coordinates
	if(pred.x != 0){
		first = &result.y;
		second = &result.z;
	} else if(pred.y != 0) {
		first = &result.x;
		second = &result.z;
	} else /* z != 0 */ {
		first = &result.x;
		second = &result.y;
	}

	// Make the movement
	if(movement == FRONT){
		result = pred;
	} else if(movement == UP) {
		*first = 1; // fill first filled coordinate positively
	} else if(movement == DOWN) {
		*first = -1;
	} else if(movement == RIGHT) {
		*second = 1;
	} else /* movement == RIGHT */ {
		*second = -1;
	}

	return result;
}

void migrch_build_3d(const shiftmel * chain,
	int chainSize,
	numtrd **coordsBB_p,
	numtrd **coordsSC_p
){
	numtrd *coordsBB;
	numtrd *coordsSC;

	// Allocate sufficient space for the coordinates
	coordsBB = malloc(sizeof(numtrd) * (chainSize + 1));
	coordsSC = malloc(sizeof(numtrd) * (chainSize + 1));

	// Add initial BB
	// As a convention, the first backbone beads are at (1, 0, 0) and (2, 0, 0).
	coordsBB[0] = numtrd_make(1, 0, 0);
	coordsBB[1] = numtrd_make(2, 0, 0);

	// Place initial SC, which are exceptions.
	// The first migrch element stores directions for the first 2 SC's.
	// All the other migrch elements store for 1 SC and 1 BB.
	shiftmel elem = chain[0];
	unsigned char mov1 = shiftmel_getBB(elem);
	unsigned char mov2 = shiftmel_getSC(elem);

	// predecessor vector and displacement vector
	numtrd predVec, dispVec;

	// Add SC beads.
	// First predecessor vector is (-1, 0, 0) from BB[1] to BB[0].
	// Second is (1, 0, 0) from BB[0] to BB[1]. 
	coordsSC[0] = numtrd_add(getNext(numtrd_make(-1, 0, 0), mov1), coordsBB[0]);
	predVec = numtrd_make(1, 0, 0); // Will feed the loop as the first predecessor vector
	coordsSC[1] = numtrd_add(getNext(predVec, mov2), coordsBB[1]);

	// Iterate over the chain
	// There should be N+1 beads and N chain elements
	int i;
	for(i = 2; i <= chainSize; i++){ // i represents index of current bead being added
		shiftmel elem = chain[i-1];
		mov1 = shiftmel_getBB(elem);
		mov2 = shiftmel_getSC(elem);

		// Get displacement vector for backbone
		dispVec = getNext(predVec, mov1);

		// Add next backbone bead
		numtrd nextBead = numtrd_add(dispVec, coordsBB[i-1]);
		coordsBB[i] = nextBead;

		// Update predecessor vector
		predVec = dispVec;

		// Get displacement vector for side chain
		dispVec = getNext(predVec, mov2);

		// Add next sidechain bead
		nextBead = numtrd_add(dispVec, coordsBB[i]);
		coordsSC[i] = nextBead;

		// Predecessor vector is kept for next iteration.
	}

	*coordsBB_p = coordsBB;
	*coordsSC_p = coordsSC;
}

/* DEBUGGING PROCEDURES
*

void print_3d_coords(numtrd *bbCo, numtrd *scCo, int size){
	int i;
	for(i = 0; i < size; i++){
		printf("(%d,%d,%d), (%d,%d,%d)\n",
			   bbCo[i].x, bbCo[i].y, bbCo[i].z,
			   scCo[i].x, scCo[i].y, scCo[i].z);
	}
}

int main(int argc, char *argv[]){
	int size = 6;

	shiftmel * chain = migrch_create(size);

	migrch_set_element(chain, 0, UP,    RIGHT);
	migrch_set_element(chain, 1, UP,    LEFT);
	migrch_set_element(chain, 2, FRONT, RIGHT);
	migrch_set_element(chain, 3, LEFT,  UP);
	migrch_set_element(chain, 4, RIGHT, DOWN);
	migrch_set_element(chain, 5, DOWN,  FRONT);

	numtrd *bbCo, *scCo;
	migrch_build_3d(chain, size, &bbCo, &scCo);

	print_3d_coords(bbCo, scCo, size+1);
}

*
*/
