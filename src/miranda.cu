/* -------------------------------------------------------------------
 * CUDA miRanda- An miRNA target scanner with CUDA speedup
  *
 * Copyright (C) (2012) UCSD DBMI Shuang Wang
 *
 * Distributed under the GNU Public License (GPL)
 * See the files 'COPYING' and 'LICENSE' for details
 *
 * Authors: Shuang Wang, */


// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

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

FloatType scale;
int strict;
int debug;
FloatType gap_open;
FloatType gap_extend;
FloatType score_threshold;
FloatType energy_threshold;
int length_5p_for_weighting;
int* length_3p_for_weighting_array;
int key_value_pairs;
int no_energy;
int verbosity;
int outfile;
int truncated;
int restricted;
int ref_batch_size;
int nGPUs;
FILE* fpout = stdout;
char filename2[200];
FILE* query_fp = 0;
pair_struct* pairs = 0;
int total_pairs = 0;

void initialize_globals() {
	initialize_bases(); /* Prepare the generic base lookup array*/
	initialize_scores();
	initialize_seqio_buffers();
}

void destroy_globals() {
	destroy_seqio_buffers();
}

int loadParametersFromArguments(int argc, char** argv){
	char filename1[200];	
	char pairs_file[200];
	char fileout[200];		
	FILE* reference_fp = 0;
	FILE* fp_pairs = 0;	
	
	/* Set Default Parameter Values*/
	length_5p_for_weighting = 8;	/* The 5' sequence length to be weighed  except for the last residue*/
	scale = 4.0;			/* The 5' miRNA scaling parameter*/
	strict = 0;			/* Strict seed model on/off*/
	debug = 0;			/* Debugging mode on/off*/
	key_value_pairs = 0;
	gap_open = -9.0;		/* Gap-open Penalty*/
	gap_extend = -4.0;		/* Gap-extend Penalty*/
	score_threshold = 140.0;	/* SW Score Threshold for reporting hits*/
	energy_threshold = 1.0;		/* Energy Threshold (DG) for reporting hits*/
	verbosity = 1;			/* Verbose mode on/off*/
	outfile = 0;			/* Dump to file on/off*/
	truncated = 0;			/* Truncate sequences on/off*/
	no_energy = 0;			/* Turn off Vienna Energy Calcs - FASTER*/
	restricted = 0;			/* Perform restricted search space*/
	ref_batch_size = 1;     /* the size of a group of ref sequence which would be processed in parallel*/
	nGPUs = 1;
	parse_command_line(argc, argv, filename1, filename2, fileout, pairs_file);
	if (gap_open > 0.0 || gap_extend > 0.0) { 
		fprintf(stderr, "Error: gap penalties may not be greater than 0\n");
		return 1;
	}
	if (truncated < 0) {
		fprintf(stderr, "Error: negative value give for UTR truncation\n");
		return 1;
	}
	if ((query_fp = fopen(filename1, "r")) == NULL) {
		fprintf(stderr, "Error: Cannot open file %s\n", filename1);
		return 1;
	} 
	if ((reference_fp = fopen(filename2, "r")) == NULL) {
		fprintf(stderr, "Error: Cannot open file %s\n", filename2);
		return 1;
	}
	fclose(reference_fp);
	if ((outfile) && ((fpout = fopen(fileout, "w")) == NULL)) {
		fprintf(stderr, "Error: Cannot create output file %s\n", fileout);
		return 1;
	}
	if (restricted) {
		if ((fp_pairs = fopen(pairs_file, "r")) == NULL) {
			fprintf(stderr, "Error: Cannot open restrict pairs file %s\n", pairs_file);
			return 1;
		}
		/* Initialize the pairs list for restriced searches*/
		total_pairs = load_pairs(fp_pairs, &pairs);
		fclose(fp_pairs);
	}
	initialize_globals();
	print_parameters(filename1, filename2, fpout);
	if (restricted && verbosity) {
		printf("Performing Restricted Scan on:%d pairs\n", total_pairs);
	}	
	return 1;   
}
    
void RUN_GPU(int argc, char **argv){
	int res = loadParametersFromArguments(argc, argv);	
	find_targets(query_fp, fpout, pairs, total_pairs, filename2);
	destroy_globals();
	if (outfile) fclose(fpout);
	fclose(query_fp);	 
}  
////////////////////////////////////////////////////////////////////////////////
// Program main
////////////////////////////////////////////////////////////////////////////////
int main( int argc, char** argv) 
{	
	
	RUN_GPU(argc, argv);	
	cutilDeviceReset(); 
}

