#ifndef __UTILS_H__
#define __UTILS_H__

#define FAILURE 1

// Certifica que um ponteiro recebeu uma alocação adequada
void must_alloc(void *ptr, const char *desc);

// Remove os fins de linha de uma string
void remove_eol(char *line);

// Aloca uma matriz não contígua de elementos de tamanho qualquer
void **alloc_matrix(unsigned int m, unsigned int n, unsigned int element_size);

// Libera o espaço ocupado por uma matriz não contígua qualquer
void free_matrix(void **matrix, unsigned int m);

#endif