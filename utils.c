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


void **alloc_matrix(unsigned int m, unsigned int n, unsigned int element_size)
{
	void **matrix = malloc(sizeof(void*) * m);
	must_alloc(matrix, "matrix");

	for (int i = 0; i < m; i++)
	{
		matrix[i] = malloc(element_size * n);
		must_alloc(matrix[i], "matrix[i]");
	}

	return matrix;
}


void free_matrix(void **matrix, unsigned int m)
{
	for (int i = 0; i < m; i++)
		free(matrix[i]);
	free(matrix);
}
