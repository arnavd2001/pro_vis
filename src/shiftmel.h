#ifndef _shiftmel_H_
#define _shiftmel_H_

/** \file shiftmel.h Routines for manipulating shiftmel units, which represent relative movements within a protein.
 *
 * We decided to encode 2 relative movements in a each byte (char), each occupying 4 bits.
 * The most significant bits hold a relative movement for the backbone chain.
 */

#include <stdio.h>
#include <stdlib.h>
#include "random.h"

/** Type that holds 2 movements, one for the backbone and one for the side chain.
 * The backbone movement is on the most significant 4 bits of this type.
 * The side chain movement is on the least significant 4 bits.
 * The type's size is supposed to be 1 byte.
 *
 * migrch is how we refer to an array of shiftmel.
 */
typedef unsigned char shiftmel;

/** Constants for making shiftmel elements */
enum _MovUnits {
	FRONT = 0,
	LEFT  = 1,
	RIGHT = 2,
	UP    = 3,
	DOWN  = 4,
};

#ifndef shiftmel_SOURCE_FILE
	#define shiftmel_INLINE inline
#else
	#define shiftmel_INLINE extern inline
#endif

/** Creates a shiftmel.
 *
 * \param bb The backbone element (e.g. FRONT, LEFT)
 * \param sc The side-chain element (e.g. FRONT, LEFT)
 */
shiftmel_INLINE
shiftmel shiftmel_make(unsigned char bb, unsigned char sc){
	return (bb << 4) | sc;
}

/** Returns a uniformly random shiftmel. */
shiftmel_INLINE
shiftmel shiftmel_random(){
	return shiftmel_make(urandom_max(DOWN+1), urandom_max(DOWN+1));
}

/** Returns the movement for the backbone (BB) stored in a shiftmel. */
shiftmel_INLINE
unsigned char shiftmel_getBB(shiftmel elem){
	return elem >> 4;
}

/** Returns the movement for the side chain (SC) stored in a shiftmel. */
shiftmel_INLINE
unsigned char shiftmel_getSC(shiftmel elem){
	return elem & 0x0F;
}

/** Convert a shiftmel to a number.
 * Each possible shiftmel map into a single number in interval [0, 24].
 */
shiftmel_INLINE
unsigned char shiftmel_to_number(shiftmel elem){
	unsigned char bb = shiftmel_getBB(elem);
	unsigned char sc = shiftmel_getSC(elem);
	return bb * 5 + sc;
}

/** Convert a number to shiftmel.
 * Each possible shiftmel map into a single number in interval [0, 24].
 */
shiftmel_INLINE
shiftmel shiftmel_from_number(unsigned char num){
	return shiftmel_make(num / 5, num % 5);
}

/** Prints a movement in the format "%c,%c", without leading/trailing spaces.
 * Prints to file 'fp'.
 */
shiftmel_INLINE
void shiftmel_print(shiftmel elem, FILE *fp){
	static const char chars[] = {'F','L','R','U','D'};
	unsigned char bb = elem >> 4;
	unsigned char sc = elem & 0x0F;
	fprintf(fp, "%c,%c", chars[bb], chars[sc]);
}

#endif
