
#pragma warning(disable:4996)
#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
// include Intellisense
#include "cuda.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
// includes, project
#include <cutil_inline.h>
#include <shrQATest.h>

extern "C"{
#include "mirandaVar.h"
#include "miranda.h"
}
#include "mirandaCPP.h"

/* invalid nucleotide matrix */
/* '!' indicates an illegal character */

char nt_invalid[256];

int initialize_nt_invalid() {
	int i;
	for (i = 0; i < 256; i++) nt_invalid[i] = '\1';
	nt_invalid['a'] = '\0';
	nt_invalid['A'] = '\0';
	nt_invalid['c'] = '\0';
	nt_invalid['C'] = '\0';
	nt_invalid['g'] = '\0';
	nt_invalid['G'] = '\0';
	nt_invalid['u'] = '\0';
	nt_invalid['U'] = '\0';
	nt_invalid['t'] = '\0';
	nt_invalid['T'] = '\0';
	nt_invalid['n'] = '\0';
	nt_invalid['N'] = '\0';
	nt_invalid['r'] = '\0';
	nt_invalid['R'] = '\0';
	nt_invalid['y'] = '\0';
	nt_invalid['Y'] = '\0';
	nt_invalid['k'] = '\0';
	nt_invalid['K'] = '\0';
	nt_invalid['m'] = '\0';
	nt_invalid['M'] = '\0';
	nt_invalid['s'] = '\0';
	nt_invalid['S'] = '\0';
	nt_invalid['w'] = '\0';
	nt_invalid['W'] = '\0';
	nt_invalid['b'] = '\0';
	nt_invalid['B'] = '\0';
	nt_invalid['d'] = '\0';
	nt_invalid['D'] = '\0';
	nt_invalid['h'] = '\0';
	nt_invalid['H'] = '\0';
	nt_invalid['v'] = '\0';
	nt_invalid['V'] = '\0';
	nt_invalid['x'] = '\0';
	nt_invalid['X'] = '\0';
	return 0;
}

ExpString *sequence_buffer = 0;
ExpString *description_buffer = 0;
ExpString *name_buffer = 0;

void initialize_seqio_buffers() {
	initialize_nt_invalid();
	create_ExpString(&sequence_buffer);
	create_ExpString(&description_buffer);
	create_ExpString(&name_buffer);
}

void destroy_seqio_buffers() {
	destroy_ExpString(&sequence_buffer);
	destroy_ExpString(&description_buffer);
	destroy_ExpString(&name_buffer);
}

/* returns 0 if EOF is reached, otherwise 1 */
int skip_leading_whitespace(FILE* fp) {
	for (;;) {
		int cbuf = fgetc(fp);
		switch (cbuf) {
			case '\n': case '\r': case '\t': case ' ': break;
			case EOF: return 0;
			default: ungetc((char)cbuf, fp); return 1;
		}
	}
	return 0; /* never reached */
}

int skip_leading_space(FILE* fp) {
	for (;;) {
		int cbuf = fgetc(fp);
		switch (cbuf) {
			case '\r': case '\t': case ' ': break;
			case EOF: return 0;
			default: ungetc((char)cbuf, fp); return 1;
		}
	}
	return 0; /* never reached */
}

/* requires file position to be at '>' */
void read_fasta_name_into_buffer(FILE* fp, ExpString *name_buf) {
	int cbuf = fgetc(fp);
	if (cbuf != '>') {
		fprintf(stderr, "Error: illegal format in sequence file (record does not start with '>')\n");
		exit(1);
	}
	clear_ExpString(name_buf);
	for (;;) {
		cbuf = fgetc(fp);
		if (cbuf == EOF) {
			fprintf(stderr, "Error: illegal format in sequence file (EOF while reading name)\n");
			exit(1);
		}
		if (cbuf == '\n' || cbuf == '\r' || cbuf == '\t' || cbuf == ' ') {
			ungetc((char)cbuf, fp);
			break;
		}
		append_char_ExpString(name_buf, (char)cbuf);
	}
}

void read_fasta_description_into_buffer(FILE* fp, ExpString *description_buf) {
	if (skip_leading_space(fp) == 0) {
		fprintf(stderr, "Error: illegal format in sequence file (EOF while reading description)\n");
		exit(1);
	}
	clear_ExpString(description_buf);
	for (;;) {
		int cbuf = fgetc(fp);
		if (cbuf == EOF) {
			fprintf(stderr, "Error: illegal format in sequence file (EOF while reading description)\n");
			exit(1);
		}
		if (cbuf == '\n') return;
		append_char_ExpString(description_buf, (char)cbuf);
	}
}

void read_fasta_sequence_into_buffer(FILE* fp, ExpString *sequence_buf, ExpString *name_buf) {
	clear_ExpString(sequence_buf);
	for (;;) {
		int cbuf = fgetc(fp);
		if (cbuf == EOF) return;
		if (cbuf == '\n' || cbuf == '\r' || cbuf == '\t' || cbuf == ' ') continue;
		if (cbuf == '>') {
			ungetc((char)cbuf,fp);
			return;
		}
		if (nt_invalid[cbuf]) {
			fprintf(stderr, "Error: illegal character '%c' in sequence file in sequence: %s\n", (char)cbuf, access_ExpString(name_buf));
			exit(1);
		}
		append_char_ExpString(sequence_buf, (char)cbuf);
	}
}

int readinseq(FILE* fp, char** sequence, char** description, char** name, int* seqlen) {
	if (skip_leading_whitespace(fp) == 0) return 0; /* end of file reached */
	read_fasta_name_into_buffer(fp, name_buffer);
	read_fasta_description_into_buffer(fp, description_buffer);
	read_fasta_sequence_into_buffer(fp, sequence_buffer, name_buffer);
	if (length_ExpString(name_buffer) == 0) {
		fprintf(stderr, "Error: sequence file contains a record with a length 0 name\n");
		exit(1);
	}
	if (length_ExpString(sequence_buffer) == 0) {
		fprintf(stderr, "Error: sequence file contains a record with a length 0 sequence\n");
		exit(1);
	}
	*name = (char *)realloc(*name, (length_ExpString(name_buffer) + 1) * sizeof(char));
	*description = (char *)realloc(*description, (length_ExpString(description_buffer) + 1) * sizeof(char));
	*sequence = (char *)realloc(*sequence, (length_ExpString(sequence_buffer) + 1) * sizeof(char));
	if (*name == 0 || *description == 0 || *sequence == 0) {
		fprintf(stderr,"Error: memory allocation failed\n");
		exit(1);
	}
	strcpy(*name, access_ExpString(name_buffer));
	strcpy(*description, access_ExpString(description_buffer));
	strcpy(*sequence, access_ExpString(sequence_buffer));
	*seqlen = length_ExpString(sequence_buffer);
	return 1;
}

/* Stefan & Shuang*/

input_type* readinallseq(FILE* fp) {
	int i = 0;
	input_type* array = (input_type *) malloc(sizeof(input_type));
	array->sequence = (char **)malloc(sizeof(char *));
	array->sequence1D = NULL;
	array->name = (char **)malloc(sizeof(char *));
	array->description = (char **)malloc(sizeof(char *));	
	array->seqlen = (int *)malloc(sizeof(int));
	array->size = 0;
	array->seqlen_max = 0;
	array->seqlen_min = 100000000;
	array->seqlen_mean = 0;
	while (skip_leading_whitespace(fp) > 0) {
		read_fasta_name_into_buffer(fp, name_buffer);
		read_fasta_description_into_buffer(fp, description_buffer);
		read_fasta_sequence_into_buffer(fp, sequence_buffer, name_buffer);
		if (length_ExpString(name_buffer) == 0) {
			fprintf(stderr, "Error: sequence file contains a record with a length 0 name\n");
			exit(1);
		}
		if (length_ExpString(sequence_buffer) == 0) {
			fprintf(stderr, "Error: sequence file contains a record with a length 0 sequence\n");
			exit(1);
		}

		array->name[array->size] = (char *)malloc((length_ExpString(name_buffer) + 1) * sizeof(char));
		array->description[array->size] = (char *)malloc((length_ExpString(description_buffer) + 1) * sizeof(char));
		array->sequence[array->size] = (char *)malloc((length_ExpString(sequence_buffer) + 1) * sizeof(char));

		if (array->name[array->size] == 0 || array->description[array->size] == 0 || array->sequence[array->size] == 0) {
			fprintf(stderr,"Error: memory allocation failed\n");
			exit(1);
		}

		strcpy(array->name[array->size], access_ExpString(name_buffer));
		strcpy(array->description[array->size], access_ExpString(description_buffer));
		strcpy(array->sequence[array->size], access_ExpString(sequence_buffer));
		array->seqlen[array->size] = length_ExpString(sequence_buffer);		
		array->seqlen_mean += array->seqlen[array->size];
		if(array->seqlen[array->size]> array->seqlen_max) array->seqlen_max = array->seqlen[array->size];
		if(array->seqlen[array->size] < array->seqlen_min) array->seqlen_min = array->seqlen[array->size];		
		array->size = array->size + 1;
		array->name = (char **)realloc(array->name, (array->size + 1) * sizeof(char*));
		array->description = (char **)realloc(array->description, (array->size + 1) * sizeof(char*));
		array->sequence = (char **)realloc(array->sequence, (array->size + 1) * sizeof(char*));
		array->seqlen = (int *)realloc(array->seqlen, (array->size + 1) * sizeof(int));
	}
	array->seqlen_mean /= array->size;
	array->seqlen_sd = 0;
	for(i = 0; i < array->size; i++){		
		array->seqlen_sd += (FloatType) pow(array->seqlen[i] - array->seqlen_mean, 2);
	}
	array->seqlen_sd = (FloatType) sqrt(array->seqlen_sd/array->size);

	return array;
}

/*Shuang */
input_type* loadQuery(FILE* query_fp){	
	int i = 0, j = 0;
	int q;
	intXY * sortedID;
	// Read the query sequence(s) (microRNA(s)) from a FASTA file
	input_type* query_array, *query_array_tmp;
	query_array_tmp = readinallseq(query_fp);
	if(((query_array_tmp->size - 1)/WARP_SIZE + 1)*WARP_SIZE > MAX_NUM_QUERY_FOR_CONSTANT_MEMORY){
		fprintf(stderr, "Error: %d queries given, at most %d are supported per run; split the query file\n", query_array_tmp->size, MAX_NUM_QUERY_FOR_CONSTANT_MEMORY);
		exit(1);
	}
	query_array_tmp->sequence1D_reshape = NULL;
	query_array_tmp->sequence1D_reshape_raw = NULL;
	sortedID = (intXY *)malloc(query_array_tmp->size*sizeof(intXY));
	for(i = 0; i < query_array_tmp->size; i++){
		sortedID[i].val = query_array_tmp->seqlen[i];
		sortedID[i].id = i;
	}
	qsort(sortedID, query_array_tmp->size, sizeof(intXY), cmpsIntXY);
	query_array = (input_type*) calloc(1, sizeof(input_type));
	query_array->seqlen_max  = query_array_tmp->seqlen_max;	
	query_array->seqlen_min  = query_array_tmp->seqlen_min;
	query_array->seqlen_mean = query_array_tmp->seqlen_mean;
	query_array->seqlen_sd   = query_array_tmp->seqlen_sd;
	query_array->size        = query_array_tmp->size;
	query_array->sequence = (char **) calloc(query_array->size, sizeof(char*));
	query_array->description = (char **) calloc(query_array->size, sizeof(char*));
	query_array->name = (char **) calloc(query_array->size, sizeof(char*));
	query_array->seqlen = (int *) calloc(MAX_NUM_QUERY_FOR_CONSTANT_MEMORY, sizeof(int));	
	query_array->sequence1D = (char *) calloc(query_array->size * MAX_QUERY_LENGHT, sizeof(char));
	query_array->size_patched = ((query_array_tmp->size - 1)/WARP_SIZE + 1)*WARP_SIZE;
	query_array->sequence1D_reshape = (char4 *) calloc(query_array->size_patched * MAX_QUERY_LENGHT/4, sizeof(char4));
	query_array->sequence1D_reshape_raw = (char2 *) calloc(query_array->size_patched * MAX_QUERY_LENGHT, sizeof(char2));
	query_array->total_seqlen = 0;
	for(i = 0; i < query_array->size; i++){
		q = sortedID[i].id;
		query_array->seqlen[i] = query_array_tmp->seqlen[q];
		query_array->sequence[i] = query_array->sequence1D + i * MAX_QUERY_LENGHT;
		memcpyCharWithFixedMaxLength(query_array->sequence[i], query_array_tmp->sequence[q], query_array_tmp->seqlen[q]+1, MAX_QUERY_LENGHT);
        query_array->name[i]  = memcpyChar(query_array_tmp->name[q], (int)strlen(query_array_tmp->name[q])+1);
		query_array->description[i] = memcpyChar(query_array_tmp->description[q], (int)strlen(query_array_tmp->description[q])+1);
	}
	// Reverse the query (microRNA) sequence for aligning
	length_3p_for_weighting_array = (int*)calloc(MAX_NUM_QUERY_FOR_CONSTANT_MEMORY, sizeof(int));
	for(i = 0; i < query_array->size; i++){
		query_array->total_seqlen += query_array->seqlen[i];
		revstring(query_array->sequence[i]);
		length_3p_for_weighting_array[i] = query_array->seqlen[i] - length_5p_for_weighting;
	}
	for(i = 0; i < MAX_QUERY_LENGHT; i++){
		for(j = 0; j < query_array->size; j ++){
			query_array->sequence1D_reshape_raw[i*query_array->size_patched + j].x = query_array->sequence[j][i];
			query_array->sequence1D_reshape_raw[i*query_array->size_patched + j].y = bases[query_array->sequence[j][i]];
		}
	}
	for(i = 0; i < MAX_QUERY_LENGHT/4; i++){
		for(j = 0; j < query_array->size; j ++){
			query_array->sequence1D_reshape[i*query_array->size_patched + j].x = bases[query_array->sequence[j][i * 4]];
			query_array->sequence1D_reshape[i*query_array->size_patched + j].y = bases[query_array->sequence[j][i * 4 + 1]];
			query_array->sequence1D_reshape[i*query_array->size_patched + j].z = bases[query_array->sequence[j][i * 4 + 2]];
			query_array->sequence1D_reshape[i*query_array->size_patched + j].w = bases[query_array->sequence[j][i * 4 + 3]];
		}
	}
	if(query_array->seqlen_max > MAX_QUERY_LENGHT){
		printf("Error lenght of a signle query is larger than %d\n", MAX_QUERY_LENGHT);
		exit(0);
	}
	freeInputType(query_array_tmp);
	free(sortedID);
	return query_array;
}

/*Shuang */
input_type* loadRef(char* filename){	
	FILE* reference_fp = 0;
	input_type* reference_array;
	int r, i;
	// Read the reference sequence(s) (microRNA(s)) from a FASTA file
	if ((reference_fp = fopen(filename, "r")) == NULL) {
		fprintf(stderr, "Error: Cannot open file %s\n", filename);
		exit(1);
	}
	reference_array = readinallseq(reference_fp);
	if (truncated && reference_array->seqlen_max > truncated){
		reference_array->seqlen_max = truncated;
	}
	if (truncated && reference_array->seqlen_min > truncated){
		reference_array->seqlen_min = truncated;
	}
	reference_array->seqlen_mean = 0;
	for(r = 0; r < reference_array->size; r++) {
		if (truncated && reference_array->seqlen[r] > truncated) {
			reference_array->sequence[r][truncated] = '\0';
			reference_array->seqlen[r] = truncated;
			
		}		
		reference_array->seqlen_mean += reference_array->seqlen[r];
	}
	reference_array->seqlen_mean /= reference_array->size;
	reference_array->seqlen_sd = 0;
	for(i = 0; i < reference_array->size; i++){		
		reference_array->seqlen_sd += (FloatType) pow(reference_array->seqlen[i] - reference_array->seqlen_mean, 2);
	}
	reference_array->seqlen_sd = (FloatType) sqrt(reference_array->seqlen_sd/reference_array->size);

	fclose(reference_fp);
	return reference_array;
	//*/
}
