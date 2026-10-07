/* pairs.c */
/* -------------------------------------------------------------------
 * CUDA miRanda- An miRNA target scanner with CUDA speedup
  *
 * Copyright (C) (2012) UCSD DBMI Shuang Wang
 *
 * Distributed under the GNU Public License (GPL)
 * See the files 'COPYING' and 'LICENSE' for details
 *
 * Authors: Shuang Wang, */
#pragma warning(disable:4244)
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

/* Load in a flat-file containing miRNA / UTR pairs to scan*/
int readinpairs(FILE* fp , pair_struct** pairs) {
	char BUFFER[400];
	int n_pairs = 0;
	int N_pairs = 16;
	if (!(*pairs = (pair_struct *) malloc(N_pairs * sizeof(pair_struct)))) {
		fprintf(stderr, "Cannot obtain %ld pair structs\n", (long)N_pairs);
		exit(1);
	}
	while (fgets(BUFFER, 400, fp)) {
		if (n_pairs >= N_pairs) {
			N_pairs *= 1.44;
			if (!(*pairs = (pair_struct *) realloc(*pairs, N_pairs * sizeof(pair_struct)))) {
				fprintf(stderr, "Cannot realloc to %ld pair structs\n", (long)N_pairs);
				exit(1);
			}
		}
		sscanf(BUFFER, "%s\t%s", (*pairs)[n_pairs].identifier1, (*pairs)[n_pairs].identifier2);
		n_pairs++;
	}
	return n_pairs;
}


int sort_pairs(const void* s1, const void* s2) {
	pair_struct* a = (pair_struct*)s1;
	pair_struct* b = (pair_struct*)s2;
	if (strcmp(a->identifier1, b->identifier1) == 0) {
		return strcmp(a->identifier2, b->identifier2);
	} else {
		return strcmp(a->identifier1, b->identifier1);
	}
}

int load_pairs(FILE* fp_pairs, pair_struct** pairs) {
	int total_pairs = readinpairs(fp_pairs, pairs);
	qsort(*pairs, total_pairs, sizeof(pair_struct), sort_pairs);
	return total_pairs;
}

int find_pair(char* ident1, char* ident2, int total_pairs, pair_struct* pairs) {
	pair_struct search_key;
	strcpy(search_key.identifier1, ident1);
	strcpy(search_key.identifier2, ident2);
	return bsearch(&search_key, pairs, total_pairs, sizeof(pair_struct), sort_pairs) != 0;
}
