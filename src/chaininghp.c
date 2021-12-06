#include <stdio.h>
#include <string.h>

#include "chaininghp.h"

HPElem *chaininghp_read(FILE *fp){
	char buffer[2001];
	char *chain;
	int retval;
	
	retval = fscanf(fp, " %2000[HP]", buffer);

	if(retval == 1){
		chain = strdup(buffer);
	} else {
		chain = NULL;
	}

	return chain;
}
