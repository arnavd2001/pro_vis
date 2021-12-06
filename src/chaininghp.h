#ifndef chaininghp_H
#define chaininghp_H

/** \file chaininghp.h Routines for handling HP chains. */

#include <stdio.h>

/** chaininghp is how we call an array of HPElem. */
typedef char HPElem;

/** Reads an HP chain from file 'fp'.
 * An HP chain is of the form "HHPHPPP", and is expected to have
 *   less than 2000 characters (overestimated).
 * The user is allowed to type only the characters "H" and "P".
 */
HPElem * chaininghp_read(FILE *fp);

#endif // chaininghp_H

