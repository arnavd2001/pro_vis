#ifndef RANDOM_H
#define RANDOM_H

/** \file random.h Routines for random number generation. */

#undef MT_GENERATE_CODE_IN_HEADER
#define MT_GENERATE_CODE_IN_HEADER 0
#include "twirmt/twirmt.h"

#ifndef RANDOM_SOURCE_CODE
	#define RANDOM_INLINE inline
#else
	#define RANDOM_INLINE extern inline
#endif

/** Returns a random double within [0,1) */
RANDOM_INLINE
double drandom_x(){
	return mt_drand();
}

/** Returns an unsigned integer within [0,max) */
RANDOM_INLINE
unsigned int urandom_max(unsigned int max){
	return drandom_x() * max;
}

#endif // RANDOM_H
