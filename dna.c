#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>

#include "utils.h"

#define NUM_QUERIES 10000 // Quantas queries serão lidas
#define MAX_WORDSIZE 100 // Tamanho máximo de uma string auxiliar
#define MAX_SUBSTRING 1001 // Tamanho máximo de uma substring
#define MAX_DATABASE 25000 // Tamanho da database
#define DNA_SECTIONS 10 // Quantas seções de dna serão analisadas[
#define STD_TAG 0

// MAX char table (ASCII)
#define MAX 256

// Boyers-Moore-Hospool-Sunday algorithm for string matching
int bmhs(char *string, char *substr) {
	int d[MAX];
	int n = strlen(string);
	int m = strlen(substr);

	int i, j, k;
	for (j = 0; j < MAX; j++)
		d[j] = m + 1;

	for (j = 0; j < m; j++)
		d[(int) substr[j]] = m - j;

	for (i = m - 1; i < n;){
		k = i;
		j = m - 1;
		while ((j >= 0) && (string[k] == substr[j])) {
			j--;
			k--;
		}
		if (j < 0)
			return k + 1;

		i += d[(int) string[i + 1]];
	}

	return -1;
}


// Lê as queries que serão buscadas
void read_queries(FILE *file, char **queries, char **queries_descs){
    char str[MAX_SUBSTRING];
    int query_id = 0;

    while (fgets(str, MAX_SUBSTRING, file) && query_id < NUM_QUERIES) {
        remove_eol(str);
        if (str[0] == '>'){
            queries_descs[query_id] = (char*) malloc(sizeof(char) * MAX_WORDSIZE);
            strcpy(queries_descs[query_id], str);
            
            !fgets(str, MAX_SUBSTRING, file);
            remove_eol(str);
            queries[query_id] = (char*) malloc(sizeof(char) * MAX_SUBSTRING);
            strcpy(queries[query_id], str);
            
            query_id++;
        }
    }
}


// Lê a database e a separa em databases menores
void read_database(FILE *file, char **bases, char **descs){
    char *line = (char*) malloc(sizeof(char) * MAX_SUBSTRING);
	must_alloc(line, "line");
	int base_id = 0;

    while (fgets(line, MAX_WORDSIZE, file)) {
        remove_eol(line);
        if (line[0] == '>'){
            bases[base_id][0] = 0;
            strcpy(descs[base_id], line);
            base_id++;
        }
        else 
            strcat(bases[base_id-1], line);
    }

	free(line);
}


void resolve_chunk(char **bases, char **descs, char **queries, unsigned int n_queries, char **query_descs, char **query_results){
	int result, found;
	char aux[MAX_WORDSIZE];

	for (int i = 0; i < n_queries; i++)
	{
		found = 0;
		sprintf(query_results[i], "%s\n", query_descs[i]);

		for (int j = 0; j < DNA_SECTIONS; j++) {
			result = bmhs(bases[j], queries[i]);
			if (result > 0) {
				sprintf(aux, "%s\n%i\n", descs[j], result);
				strcat(query_results[i], aux);
				found = 1;
			}
		}
		if (!found) {
			sprintf(aux, "NOT FOUND\n");
			strcat(query_results[i], aux);
		}
	}
}

void **alloc_matrix(unsigned int m, unsigned int n, unsigned int element_size){
	void **matrix = malloc(sizeof(void*) * m);
	must_alloc(matrix, "matrix");

	for (int i = 0; i < m; i++){
		matrix[i] = malloc(element_size * n);
		must_alloc(matrix[i], "matrix[i]");
	}

	return matrix;
}


int main(int argc, char **argv) {

	int rank, n_procs;
	MPI_Status status;

	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &n_procs);

	char **bases = (char**) alloc_matrix(DNA_SECTIONS, MAX_DATABASE, sizeof(char));
	char **descs = (char**) alloc_matrix(DNA_SECTIONS, MAX_WORDSIZE, sizeof(char));

	int slice = NUM_QUERIES / n_procs;
	int limit = rank == n_procs - 1 ? NUM_QUERIES - (n_procs - 1) * slice : (rank == 0 ? NUM_QUERIES : slice);

	char **queries = (char**) alloc_matrix(limit, MAX_SUBSTRING, sizeof(char));
	char **query_descs = (char**) alloc_matrix(limit, MAX_WORDSIZE, sizeof(char));
	char **query_results = (char**) alloc_matrix(limit, MAX_SUBSTRING, sizeof(char));

	if (rank == 0) {
		FILE *fdatabase = fopen("dna.in", "r");
		must_alloc(fdatabase, "fdatabase");

		FILE *fquery = fopen("query.in", "r");
		must_alloc(fquery, "fquery");

		FILE *fout = fopen("dna.out", "w");
		must_alloc(fout, "fout");

		read_database(fdatabase, bases, descs);
		read_queries(fquery, queries, query_descs);

		for (int i = 1; i < n_procs; i++){
			for (int j = 0; j < DNA_SECTIONS; j++){
				MPI_Send(bases[j], MAX_DATABASE, MPI_CHAR, i, STD_TAG, MPI_COMM_WORLD);
				MPI_Send(descs[j], MAX_WORDSIZE, MPI_CHAR, i, STD_TAG, MPI_COMM_WORLD);
			}
			for (int j = i * slice; j < (i + 1) * slice; j++){
				MPI_Send(queries[j], MAX_SUBSTRING, MPI_CHAR, i, STD_TAG, MPI_COMM_WORLD);
				MPI_Send(query_descs[j], MAX_WORDSIZE, MPI_CHAR, i, STD_TAG, MPI_COMM_WORLD);
			}
		}
		for (int i = slice * (n_procs - 1); i < NUM_QUERIES; i++){
			MPI_Send(queries[i], MAX_SUBSTRING, MPI_CHAR, n_procs - 1, STD_TAG, MPI_COMM_WORLD);
			MPI_Send(query_descs[i], MAX_WORDSIZE, MPI_CHAR, n_procs - 1, STD_TAG, MPI_COMM_WORLD);
		}

		fclose(fdatabase);
		fclose(fquery);

		resolve_chunk(bases, descs, queries, slice, query_descs, query_results);

		for (int i = 1; i < n_procs; i++)
			for (int j = i * slice; j < (i + 1) * slice; j++)
				MPI_Recv(query_results[j], MAX_SUBSTRING, MPI_CHAR, i, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

		for (int i = 0; i < NUM_QUERIES; i++)
			fprintf(fout, "%s", query_results[i]);

		fclose(fout);
	}
	else {
		for (int j = 0; j < DNA_SECTIONS; j++){
			MPI_Recv(bases[j], MAX_DATABASE, MPI_CHAR, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			MPI_Recv(descs[j], MAX_WORDSIZE, MPI_CHAR, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		}
		for (int j = 0; j < limit; j++){
			MPI_Recv(queries[j], MAX_SUBSTRING, MPI_CHAR, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			MPI_Recv(query_descs[j], MAX_WORDSIZE, MPI_CHAR, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		}

		int result, found;
		char aux[MAX_WORDSIZE];

		for (int i = 0; i < limit; i++)
		{
			found = 0;
			sprintf(query_results[i], "%s\n", query_descs[i]);

			for (int j = 0; j < DNA_SECTIONS; j++) {
				result = bmhs(bases[j], queries[i]);
				if (result > 0) {
					sprintf(aux, "%s\n%i\n", descs[j], result);
					strcat(query_results[i], aux);
					found = 1;
				}
			}
			if (!found) {
				sprintf(aux, "NOT FOUND\n");
				strcat(query_results[i], aux);
			}
		}

		for (int i = 0; i < limit; i++)
			MPI_Send(query_results[i], MAX_SUBSTRING, MPI_CHAR, 0, STD_TAG, MPI_COMM_WORLD);
	}

	MPI_Finalize();

	// char *query_results[NUM_QUERIES];

	// int result, found;
	// char aux[MAX_WORDSIZE];
	// for (int i = rank; i < NUM_QUERIES; i += n_procs)
	// {
	// 	found = 0;
	// 	query_results[i] = (char*) malloc(sizeof(char) * 1000);
	// 	sprintf(query_results[i], "%s\n", queries_descs[i]);

	// 	for (int j = 0; j < DNA_SECTIONS; j++) {
	// 		result = bmhs(bases[j], queries[i]);
	// 		if (result > 0) {
	// 			sprintf(aux, "%s\n%i\n", descs[j], result);
	// 			strcat(query_results[i], aux);
	// 			found = 1;
	// 		}
	// 	}
	// 	if (!found) {
	// 		sprintf(aux, "NOT FOUND\n");
	// 		strcat(query_results[i], aux);
	// 	}
	// }
	

	// for (int i = 0; i < NUM_QUERIES; i++)
	// 	fprintf(fout, "%s", query_results[i]);
	
	
	
	return EXIT_SUCCESS;
}

// int main(int argc, char **argv) {

// 	MPI_Init(&argc, &argv);

// 	int rank;
// 	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

// 	FILE *fdatabase = fopen("dna.in", "r");
// 	must_alloc(fdatabase, "fdatabase");

// 	FILE *fquery = fopen("query.in", "r");
// 	must_alloc(fquery, "fquery");

// 	FILE *fout = fopen("dna.out", "w");
// 	must_alloc(fout, "fout");

// 	char *bases[DNA_SECTIONS];
// 	char *descs[DNA_SECTIONS];

// 	char *queries[NUM_QUERIES];
// 	char *queries_descs[NUM_QUERIES];

// 	read_database(fdatabase, bases, descs);
// 	read_queries(fquery, queries, queries_descs);

// 	char *query_results[NUM_QUERIES];


// 	int result, found;
// 	char aux[MAX_WORDSIZE];
// 	for (int i = 0; i < NUM_QUERIES; i++)
// 	{
// 		found = 0;
// 		query_results[i] = (char*) malloc(sizeof(char) * 1000);
// 		sprintf(query_results[i], "%s\n", queries_descs[i]);

// 		for (int j = 0; j < DNA_SECTIONS; j++) {
// 			result = bmhs(bases[j], queries[i]);
// 			if (result > 0) {
// 				sprintf(aux, "%s\n%i\n", descs[j], result);
// 				strcat(query_results[i], aux);
// 				found = 1;
// 			}
// 		}
// 		if (!found) {
// 			sprintf(aux, "NOT FOUND\n");
// 			strcat(query_results[i], aux);
// 		}
// 	}
	

// 	for (int i = 0; i < NUM_QUERIES; i++)
// 		fprintf(fout, "%s", query_results[i]);
	
// 	for (int i = 0; i < NUM_QUERIES; i++){
// 		free(queries[i]);
// 		free(queries_descs[i]);
// 		free(query_results[i]);
// 	}

// 	for (int i = 0; i < DNA_SECTIONS; i++){
// 		free(bases[i]);
// 		free(descs[i]);
// 	}

// 	fclose(fdatabase);
// 	fclose(fquery);
// 	fclose(fout);

// 	MPI_Finalize();
	
// 	return EXIT_SUCCESS;
// }
