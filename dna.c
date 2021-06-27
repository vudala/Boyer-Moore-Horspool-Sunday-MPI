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
#define ROOT 0

// MAX char table (ASCII)
#define MAX 256

// Boyers-Moore-Hospool-Sunday algorithm for string matching
int bmhs(char *string, char *substr)
{
	int d[MAX];
	int n = strlen(string);
	int m = strlen(substr);

	int i, j, k;
	for (j = 0; j < MAX; j++)
		d[j] = m + 1;

	for (j = 0; j < m; j++)
		d[(int) substr[j]] = m - j;

	for (i = m - 1; i < n;)
	{
		k = i;
		j = m - 1;
		while ((j >= 0) && (string[k] == substr[j]))
		{
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
void read_queries(FILE *file, char **queries, char **queries_descs)
{
    char *str = (char*) malloc(sizeof(char) * MAX_SUBSTRING);
	must_alloc(str, "str");
    int query_id = 0;

    while (fgets(str, MAX_SUBSTRING, file) && query_id < NUM_QUERIES)
	{
        remove_eol(str);
        if (str[0] == '>')
		{
            strcpy(queries_descs[query_id], str);
            !fgets(str, MAX_SUBSTRING, file);
            remove_eol(str);
            strcpy(queries[query_id], str);
            
            query_id++;
        }
    }

	free(str);
	str = NULL;
}


// Lê a database e a separa em databases menores
void read_database(FILE *file, char **bases, char **descs)
{
    char *line = (char*) malloc(sizeof(char) * MAX_SUBSTRING);
	must_alloc(line, "line");
	int base_id = 0;

    while (fgets(line, MAX_WORDSIZE, file))
	{
        remove_eol(line);
        if (line[0] == '>')
		{
            bases[base_id][0] = 0;
            strcpy(descs[base_id], line);
            base_id++;
        }
        else 
            strcat(bases[base_id-1], line);
    }

	free(line);
	line = NULL;
}


// Aplica o BMHS sobre um subvetor de queries de tamanho n
void solve_chunk(char **bases, char **descs, char **queries, unsigned int n_queries, char **query_descs, char **query_results)
{
	int result, found;
	char *aux = (char*) malloc(sizeof(char) * MAX_WORDSIZE);
	must_alloc(aux, "aux");

	for (int i = 0; i < n_queries; i++)
	{
		found = 0;
		sprintf(query_results[i], "%s\n", query_descs[i]);
		for (int j = 0; j < DNA_SECTIONS; j++)
		{
			result = bmhs(bases[j], queries[i]);
			if (result > 0)
			{
				sprintf(aux, "%s\n%i\n", descs[j], result);
				strcat(query_results[i], aux);
				found = 1;
			}
		}
		if (!found)
		{
			sprintf(aux, "NOT FOUND\n");
			strcat(query_results[i], aux);
		}
	}

	free(aux);
	aux = NULL;
}


// Determina o tamanho da fatia de queries baseado no rank
int rank_limit(unsigned int rank, unsigned int n_procs, unsigned int slice)
{
	return rank == n_procs - 1 ? NUM_QUERIES - (n_procs - 1) * slice : (rank == 0 ? NUM_QUERIES : slice);
}


int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);

	int rank, n_procs;
	MPI_Status status;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &n_procs);

	char **bases = (char**) alloc_matrix(DNA_SECTIONS, MAX_DATABASE, sizeof(char));
	char **descs = (char**) alloc_matrix(DNA_SECTIONS, MAX_WORDSIZE, sizeof(char));
	int bases_buff_size = DNA_SECTIONS * (MAX_WORDSIZE + MAX_DATABASE);
	char *bases_buff = (char*) malloc(sizeof(char) * bases_buff_size);
	must_alloc(bases_buff, "bases_buff");
	
	if (rank == ROOT)
	{
		FILE *fdatabase = fopen("dna.in", "r");
		must_alloc(fdatabase, "fdatabase");
		read_database(fdatabase, bases, descs);
		fclose(fdatabase);

		int pos = 0;
		for (int j = 0; j < DNA_SECTIONS; j++)
		{
			MPI_Pack(bases[j], MAX_DATABASE, MPI_CHAR, bases_buff, bases_buff_size, &pos, MPI_COMM_WORLD);
			MPI_Pack(descs[j], MAX_WORDSIZE, MPI_CHAR, bases_buff, bases_buff_size, &pos, MPI_COMM_WORLD);
		}
	}
	MPI_Bcast(bases_buff, bases_buff_size, MPI_PACKED, ROOT, MPI_COMM_WORLD);

	int slice = NUM_QUERIES / n_procs; // Tamanho da fatia
	int limit = rank_limit(rank, n_procs, slice); // Tamanho da fatia baseado no rank

	char **queries = (char**) alloc_matrix(limit, MAX_SUBSTRING, sizeof(char));
	char **query_descs = (char**) alloc_matrix(limit, MAX_WORDSIZE, sizeof(char));
	char **query_results = (char**) alloc_matrix(limit, MAX_SUBSTRING, sizeof(char));

	int query_buff_size = (rank_limit(n_procs - 1, n_procs, slice)) * (MAX_SUBSTRING + MAX_WORDSIZE);
	char *query_buff = (char*) malloc(sizeof(char) * query_buff_size);
	must_alloc(query_buff, "query_buff");

	if (rank == ROOT)
	{
		FILE *fquery = fopen("query.in", "r");
		must_alloc(fquery, "fquery");

		// Lê as queries
		read_queries(fquery, queries, query_descs);

		fclose(fquery);

		MPI_Request req;
		int pos;

		// Fatia e empacota as queries e as bases para serem enviadas a todos os outros processos
		for (int i = 1; i < n_procs; i++)
		{
			int r_limit = rank_limit(i, n_procs, slice);
			pos = 0;

			for (int j = i * slice; j < i * slice + r_limit; j++)
			{
				MPI_Pack(queries[j], MAX_SUBSTRING, MPI_CHAR, query_buff, query_buff_size, &pos, MPI_COMM_WORLD);
				MPI_Pack(query_descs[j], MAX_WORDSIZE, MPI_CHAR, query_buff, query_buff_size, &pos, MPI_COMM_WORLD);
			}
		
			MPI_Send(query_buff, query_buff_size, MPI_PACKED, i, STD_TAG, MPI_COMM_WORLD);
		}

		// Aplica o BMHS sobre a primeira fatia de queries
		solve_chunk(bases, descs, queries, slice, query_descs, query_results);
	
		// Desempacota os dados processados pelos outros processos
		for (int i = 1; i < n_procs; i++)
		{
			pos = 0;
			int r_limit = rank_limit(i, n_procs, slice);
			MPI_Recv(query_buff, query_buff_size, MPI_PACKED, i, STD_TAG, MPI_COMM_WORLD, &status);
			for (int j = i * slice; j < i * slice + r_limit; j++)
				MPI_Unpack(query_buff, query_buff_size, &pos, query_results[j], MAX_SUBSTRING, MPI_CHAR, MPI_COMM_WORLD);
		}
		
		// Escreve os resultados em fout
		FILE *fout = fopen("dna.out", "w");
		must_alloc(fout, "fout");

		for (int i = 0; i < NUM_QUERIES; i++)
			fprintf(fout, "%s", query_results[i]);
		
		fclose(fout);
	}
	else
	{
		int pos = 0;
		for (int j = 0; j < DNA_SECTIONS; j++)
		{
			MPI_Unpack(bases_buff, bases_buff_size, &pos, bases[j], MAX_DATABASE, MPI_CHAR, MPI_COMM_WORLD);
			MPI_Unpack(bases_buff, bases_buff_size, &pos, descs[j], MAX_WORDSIZE, MPI_CHAR, MPI_COMM_WORLD);
		}

		// Desempacota as informações enviadas pelo processo root
		MPI_Recv(query_buff, query_buff_size, MPI_PACKED, 0, STD_TAG, MPI_COMM_WORLD, &status);
		pos = 0;
		for (int j = 0; j < limit; j++)
		{
			MPI_Unpack(query_buff, query_buff_size, &pos, queries[j], MAX_SUBSTRING, MPI_CHAR, MPI_COMM_WORLD);
			MPI_Unpack(query_buff, query_buff_size, &pos, query_descs[j], MAX_WORDSIZE, MPI_CHAR, MPI_COMM_WORLD);
		}

		// Aplica o BMHS utilizando queries recebidas do root
		solve_chunk(bases, descs, queries, limit, query_descs, query_results);

		// Envia o resultado do processo de volta ao root
		int res_buff_size = limit * MAX_SUBSTRING;
		char *res_buffer = (char*) malloc(sizeof(char) * res_buff_size);
		must_alloc(res_buffer, "res_buffer");
		pos = 0;

		// Empacota em res_buffer as informações pra enviar de volta ao root
		for (int i = 0; i < limit; i++)
			MPI_Pack(query_results[i], MAX_SUBSTRING, MPI_CHAR, res_buffer, res_buff_size, &pos, MPI_COMM_WORLD);

		MPI_Send(res_buffer, res_buff_size, MPI_PACKED, 0, STD_TAG, MPI_COMM_WORLD);

		free(res_buffer);
		res_buffer = NULL;
	}

	// libera e desponta os espaços de memória da aplicação
	free_matrix((void**) bases, DNA_SECTIONS);
	free_matrix((void**) descs, DNA_SECTIONS);
	free_matrix((void**) queries, limit);
	free_matrix((void**) query_descs, limit);
	free_matrix((void**) query_results, limit);
	free(query_buff);
	free(bases_buff);

	bases = NULL;
	descs = NULL;
	bases_buff = NULL;
	queries = NULL;
	query_descs = NULL;
	query_results = NULL;
	query_buff = NULL;

	MPI_Finalize();
	
	return EXIT_SUCCESS;
}