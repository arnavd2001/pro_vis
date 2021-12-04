#ifndef migrch_H
#define migrch_H

/** \file migrch.h Routines for managing chains of shiftmel units. */

#include <stdio.h>
#include "shiftmel.h"
#include "numtrd.h"

/** Changes the given 'chain' in position 'eleIdx'.
 * The element in that position becomes set with movement 'bb' for the
 *   backbone and 'sc' for the side chain.
 */
void migrch_set_element(shiftmel * chain, int eleIdx, unsigned char bb, unsigned char sc);

/** Creates a chain of movements with 'size' bytes.
 * Each byte carries 4 bits (MSB) representing backbone movement, and 4 bits (LSB) for sidechain movement.
 * Everything is initialized to FRONT-right
 */
shiftmel * migrch_create(int size);

/** Takes a chain of movements and returns the spatial position of BB and SC beads over the 3D space.
 * 'coordsBB_p' and 'coordsSC_p' are pointers to where we should store, respectively, the coordinates for
 *    the backbone beads and the side chain beads.
 */
void migrch_build_3d(const shiftmel * chain, // input
	int chainSize,    // input
	numtrd **coordsBB_p, // output
	numtrd **coordsSC_p  // output
);


#endif // migrch_H
