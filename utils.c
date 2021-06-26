// Eduardo Vudala Senoski GRR20195689

#include "utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void remove_eol(char *line){
	int i = strlen(line) - 1;
	while (line[i] == '\n' || line[i] == '\r') {
		line[i] = 0;
		i--;
	}
}

void must_alloc(void *ptr, const char *desc){
	if (!ptr){
		fprintf(stderr, "Malloc failure: %s", desc);
		exit(FAILURE);
	}
}
