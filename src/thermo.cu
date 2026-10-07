/* thermo.c */
/* -------------------------------------------------------------------
 * CUDA miRanda- An miRNA target scanner with CUDA speedup
  *
 * Copyright (C) (2012) UCSD DBMI Shuang Wang
 *
 * Distributed under the GNU Public License (GPL)
 * See the files 'COPYING' and 'LICENSE' for details
 *
 * Authors: Shuang Wang, */

/*
 * This file includes headers and links to the RNAlib library of
 * Ivo Hofackers Vienna Package
 */
//#pragma warning(disable:4018)
//#pragma warning(disable:4267)
#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
// include Intellisense
#include "cuda.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
// includes, project
#include <cutil_inline.h>
#include <shrQATest.h>

extern "C"{
#include "./H/fold_vars.h"
#include "./H/RNAlib_utils.h"
#include "./H/fold.h"
#include "mirandaVar.h"
#include "miranda.h"
}
#include "mirandaCPP.h"

FloatType vfold(char* sequence) {
	void* struct1;
	FloatType e1;
	struct1 = (char*)space(sizeof(char) * ((unsigned int)strlen(sequence) + 1));
	temperature = 30;
	initialize_fold((int)strlen(sequence));
	e1 = fold(sequence, (char*)struct1);
	free_arrays();
	free(struct1);
	return e1;
}

FloatType get_energy(hit_struct* hit) {
	FloatType energy = 0;
	int i = 0;
	int j = 0;
	char foldsequence[5000];
	revstring(hit->alignment[0]);
	revstring(hit->rest[1]);
	revstring(hit->rest[4]);
	foldsequence[0] = '\0';
	for (i = 0; i < strlen(hit->rest[0]); i++) {
		foldsequence[j] = hit->rest[0][i];
		j++;
	}
	for (i = 0; i < strlen(hit->alignment[0]); i++) {
		if (hit->alignment[0][i] != '-') {
			foldsequence[j] = hit->alignment[0][i];
			j++;
		}
	}
	for (i = 0; i < strlen(hit->rest[3]); i++) {
		foldsequence[j] = hit->rest[3][i];
		j++;
	}
	for (i = 0; i < 7; i++) {
		foldsequence[j] = 'X';
		j++;
	}
	for (i = 0; i < strlen(hit->rest[4]); i++) {
		foldsequence[j] = hit->rest[4][i];
		j++;
	}
	for (i = 0; i < strlen(hit->alignment[2]); i++) {
		if (hit->alignment[2][i] != '-') {
			foldsequence[j] = hit->alignment[2][i];
			j++;
		}
	}
	for (i = 0; i < strlen(hit->rest[1]); i++) {
		foldsequence[j] = hit->rest[1][i];
		j++;
	}
	foldsequence[j] = '\0';
	/* printf("FOLD: %s\n", foldsequence);*/
	energy = vfold(foldsequence);
	revstring(hit->alignment[0]);
	revstring(hit->rest[1]);
	revstring(hit->rest[4]);
	return energy;
}
