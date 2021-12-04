#ifndef numtrd_H
#define numtrd_H

/** \file numtrd.h Routines for mathematical manipulation of vectors and points in the integer space. */

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef numtrd_SOURCE_FILE
	#define numtrd_INLINE inline
#else // If it's source file, declare extern inline
	#define numtrd_INLINE extern inline
#endif

/** Type for representing a tridimensional coordinate. */
typedef struct numtrd_ {
	int x; /**< Coordinate x. */
	int y; /**< Coordinate y. */
	int z; /**< Coordinate z. */
} numtrd;

/** Creates an numtrd. */
numtrd_INLINE
numtrd numtrd_make(int x, int y, int z){
	numtrd num;
	num.x = x;
	num.y = y;
	num.z = z;
	return num;
}

/** Adds two numtrd objects. */
numtrd_INLINE
numtrd numtrd_add(numtrd a, numtrd b){
	numtrd res;
	res.x = a.x + b.x;
	res.y = a.y + b.y;
	res.z = a.z + b.z;
	return res;
}

/** Verifies if `a` and `b` are within a distance of exactly 1 from each other. */
numtrd_INLINE
bool numtrd_isDist1(numtrd a, numtrd b){
	int dx = abs(a.x - b.x);
	int dy = abs(a.y - b.y);
	int dz = abs(a.z - b.z);
	int sum = dx + dy + dz;
	return sum == 1 ? true : false;
}

/** Verifies if `a` and `b` are equal. */
numtrd_INLINE
bool numtrd_equal(numtrd a, numtrd b){
	if(a.x != b.x || a.y != b.y || a.z != b.z)
		return false;
	return true;
}

/** Prints the given numtrd object, without placing a newline. */
numtrd_INLINE
void numtrd_print(numtrd a, FILE *fp){
	fprintf(fp, "%d,%d,%d", a.x, a.y, a.z);
}

#endif // numtrd_H

