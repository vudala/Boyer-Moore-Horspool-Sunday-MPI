// Eduardo Vudala Senoski GRR20195689

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>

#include "utils.h"

#define NUM_QUERIES 80000 // Quantas queries serão lidas
#define MAX_WORDSIZE 100 // Tamanho máximo de uma string auxiliar
#define MAX_SUBSTRING 1001 // Tamanho máximo de uma substring
#define MAX_DATABASE 25000 // Tamanho da database
#define DNA_SECTIONS 10 // Quantas seções de dna serão analisadas
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
    char *str = (char*) malloc(sizeof(char) * MAX_SUBSTRING);
	must_alloc(str, "str");
    int query_id = 0;

    while (fgets(str, MAX_SUBSTRING, file) && query_id < NUM_QUERIES) {
        remove_eol(str);
        if (str[0] == '>'){
            strcpy(queries_descs[query_id], str);
            !fgets(str, MAX_SUBSTRING, file);
            remove_eol(str);
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


void solve_chunk(char **bases, char **descs, char **queries, unsigned int n_queries, char **query_descs, char **query_results){
	int result, found;
	char *aux = (char*) malloc(sizeof(char) * MAX_WORDSIZE);
	must_alloc(aux, "aux");

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

// determina o tamanho da fatia de queries baseado no rank
int rank_limit(unsigned int rank, unsigned int n_procs, unsigned int slice){
	return rank == n_procs - 1 ? NUM_QUERIES - (n_procs - 1) * slice : (rank == 0 ? NUM_QUERIES : slice);
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
	int limit = rank_limit(rank, n_procs, slice);

	char **queries = (char**) alloc_matrix(limit, MAX_SUBSTRING, sizeof(char));
	char **query_descs = (char**) alloc_matrix(limit, MAX_WORDSIZE, sizeof(char));
	
	char **query_results = (char**) alloc_matrix(limit, MAX_SUBSTRING, sizeof(char));

	int buff_size = DNA_SECTIONS * (MAX_WORDSIZE + MAX_DATABASE) + (rank_limit(n_procs - 1, n_procs, slice)) * (MAX_SUBSTRING + MAX_WORDSIZE);
	char *buffer = (char*) malloc(sizeof(char) * buff_size);
	must_alloc(buffer, "buffer");

	if (rank == 0) {
		FILE *fdatabase = fopen("dna.in", "r");
		must_alloc(fdatabase, "fdatabase");

		FILE *fquery = fopen("query.in", "r");
		must_alloc(fquery, "fquery");

		read_database(fdatabase, bases, descs);
		read_queries(fquery, queries, query_descs);

		fclose(fdatabase);
		fclose(fquery);

		MPI_Request req;
		int pos;
		for (int i = 1; i < n_procs; i++){
			int r_limit = rank_limit(i, n_procs, slice);
			pos = 0;

			for (int j = 0; j < DNA_SECTIONS; j++){
				MPI_Pack(bases[j], MAX_DATABASE, MPI_CHAR, buffer, buff_size, &pos, MPI_COMM_WORLD);
				MPI_Pack(descs[j], MAX_WORDSIZE, MPI_CHAR, buffer, buff_size, &pos, MPI_COMM_WORLD);
			}
			for (int j = i * slice; j < i * slice + r_limit; j++){
				MPI_Pack(queries[j], MAX_SUBSTRING, MPI_CHAR, buffer, buff_size, &pos, MPI_COMM_WORLD);
				MPI_Pack(query_descs[j], MAX_WORDSIZE, MPI_CHAR, buffer, buff_size, &pos, MPI_COMM_WORLD);
			}
		
			MPI_Isend(buffer, buff_size, MPI_PACKED, i, STD_TAG, MPI_COMM_WORLD, &req);
		}

		solve_chunk(bases, descs, queries, slice, query_descs, query_results);
	
		for (int i = 1; i < n_procs; i++){
			pos = 0;
			int r_limit = rank_limit(i, n_procs, slice);
			MPI_Recv(buffer, buff_size, MPI_PACKED, i, STD_TAG, MPI_COMM_WORLD, &status);
			for (int j = i * slice; j < i * slice + r_limit; j++)
				MPI_Unpack(buffer, buff_size, &pos, query_results[j], MAX_SUBSTRING, MPI_CHAR, MPI_COMM_WORLD);
		}

		FILE *fout = fopen("dna.out", "w");
		must_alloc(fout, "fout");

		for (int i = 0; i < NUM_QUERIES; i++)
			fprintf(fout, "%s", query_results[i]);
		
		fclose(fout);
	}
	else {
		MPI_Recv(buffer, buff_size, MPI_PACKED, 0, STD_TAG, MPI_COMM_WORLD, &status);

		int pos = 0;
		for (int j = 0; j < DNA_SECTIONS; j++){
			MPI_Unpack(buffer, buff_size, &pos, bases[j], MAX_DATABASE, MPI_CHAR, MPI_COMM_WORLD);
			MPI_Unpack(buffer, buff_size, &pos, descs[j], MAX_WORDSIZE, MPI_CHAR, MPI_COMM_WORLD);
		}
		for (int j = 0; j < limit; j++){
			MPI_Unpack(buffer, buff_size, &pos, queries[j], MAX_SUBSTRING, MPI_CHAR, MPI_COMM_WORLD);
			MPI_Unpack(buffer, buff_size, &pos, query_descs[j], MAX_WORDSIZE, MPI_CHAR, MPI_COMM_WORLD);
		}

		solve_chunk(bases, descs, queries, limit, query_descs, query_results);
		
		int res_buff_size = limit * MAX_SUBSTRING;
		char *res_buffer = (char*) malloc(sizeof(char) * res_buff_size);
		must_alloc(res_buffer, "res_buffer");
		pos = 0;
		for (int i = 0; i < limit; i++)
			MPI_Pack(query_results[i], MAX_SUBSTRING, MPI_CHAR, res_buffer, res_buff_size, &pos, MPI_COMM_WORLD);

		MPI_Request req;
		MPI_Isend(res_buffer, res_buff_size, MPI_PACKED, 0, STD_TAG, MPI_COMM_WORLD, &req);
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