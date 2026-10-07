/* utils.c */
/* -------------------------------------------------------------------
 * CUDA miRanda- An miRNA target scanner with CUDA speedup
  *
 * Copyright (C) (2012) UCSD DBMI Shuang Wang
 *
 * Distributed under the GNU Public License (GPL)
 * See the files 'COPYING' and 'LICENSE' for details
 *
 * Authors: Shuang Wang, */
//#pragma warning(disable:4267)
#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
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

int cmpsIntXY(const void* p1, const void* p2){
	intXY * s1 = (intXY*) p1;
	intXY * s2 = (intXY*) p2;
	return (s2->val - s1->val);
}
int cmpscores(const void* p1, const void* p2) {
	ushort4* s1 = (ushort4*)p1;
	ushort4* s2 = (ushort4*)p2;
	/* primary key: score */
	if (s1->x < s2->x) return 1;
	if (s1->x > s2->x) return -1;
	/* secondary key: earlier end in reference */
	if (s1->w > s2->w) return 1;
	if (s1->w < s2->w) return -1;
	return 0;
}

void clear_hit(hit_struct* hit, int query_length, int reference_length) {
	hit->score = 0;
	hit->query_start = 0;
	hit->query_end = 0;
	hit->ref_start = 0;
	hit->ref_end = 0;
	memset(hit->alignment[0], '\0', MAX_QUERY_LENGHT);
	memset(hit->alignment[1], '\0', MAX_QUERY_LENGHT);
	memset(hit->alignment[2], '\0', MAX_QUERY_LENGHT);
	memset(hit->rest[0], '\0', query_length + 10);
	memset(hit->rest[1], '\0', query_length + 10);
	memset(hit->rest[2], '\0', query_length + 10);
	memset(hit->rest[3], '\0', query_length + 10);
	memset(hit->rest[4], '\0', query_length + 10);
	memset(hit->rest[5], '\0', query_length + 10);
}

void revstring(char s[]) {
	int c, i, j;
	for (i = 0, j = (int) (strlen(s) - 1); i < j; i++, j--) {
		c = s[i];
		s[i] = s[j];
		s[j] = c;
	}
}

void string_toupper(char *s) {
	while (*s) {
		*s = toupper(*s);
		s++;
	}
}

void freeInputType(input_type * in){
	int i = 0;
	
	if(in->sequence1D_reshape_raw != NULL){
		free(in->sequence1D_reshape_raw);
	}
	
	if(in->sequence1D_reshape != NULL){
		free(in->sequence1D_reshape);
	}
	
	if(in->sequence1D != NULL){
	
		free(in->sequence1D);
	}else{
	
		for (i = 0; i < in->size; ++i) {
	
			free(in->sequence[i]);
		}
	
	}
	
	for (i = 0; i < in->size; ++i) {		
		free(in->description[i]);
		free(in->name[i]);
	}
	
	free(in->sequence);	
	free(in->description);
	free(in->name);
	free(in->seqlen);
	free(in);
}

void free2dFloatType(FloatType ** in, int len){
	int i = 0;
	for(i = 0; i< len; i++){
		free(in[i]);
	}
	free(in);
}

FloatType ** new2dFloatType(int size1, int size2){
	FloatType ** out;
	int i;
	out = (FloatType**) calloc(size1, sizeof(FloatType*));
	for(i = 0; i < size1; i++){
		out[i] = (FloatType *) calloc(size2, sizeof(FloatType));
	}
	return out;
}

void free2dHitSummary(HitSummary ** in, int size1, int size2){
	int i = 0;
	int j = 0;
	for(i = 0; i< size1; i++){
		for(j = 0; j < size2; j++){
			destroy_ExpString(&(in[i][j].position_list));
		}
		free(in[i]);
	}
	free(in);
}

HitSummary ** new2dHitSummary(int size1, int size2){
	HitSummary ** out;
	int i;
	int j;
	out = (HitSummary**) calloc(size1, sizeof(HitSummary*));
	for(i = 0; i < size1; i++){
		out[i] = (HitSummary *) calloc(size2, sizeof(HitSummary));
		for(j = 0; j < size2; j++){
			out[i][j].position_list = 0;
			create_ExpString(&(out[i][j].position_list));
			out[i][j].no_hits = 0;
			out[i][j].max_hit = 0;
			out[i][j].max_score = 0;
			out[i][j].scan_score = 0;
			out[i][j].total_score = 0;
			clear_ExpString(out[i][j].position_list);
		}
	}
	return out;
}

void free2dInt(int ** in, int len){
	int i = 0;
	for(i = 0; i< len; i++){
		free(in[i]);
	}
	free(in);
}

int ** new2dInt(int size1, int size2){
	int ** out;
	int i;
	out = (int**) calloc(size1, sizeof(int*));
	for(i = 0; i < size1; i++){
		out[i] = (int *) calloc(size2, sizeof(int));
	}
	return out;
}

void clear2dInt(int **in, int size1, int size2){	
	int i;	
	for(i = 0; i < size1; i++){		
			memset(in[i], 0, size2*sizeof(int));
	}
	
}

void free3dInt(int *** in, int size1, int size2){
	int i, j;	
	for(i = 0; i< size1; i++){
		for(j = 0; j< size2; j++){
			free(in[i][j]);
		}
		free(in[i]);
	}
	free(in);
}

int *** new3dInt(int size1, int size2, int size3){
	int *** out;
	int i, j;
	out = (int***) calloc(size1, sizeof(int**));
	for(i = 0; i < size1; i++){
		out[i] = (int **) calloc(size2, sizeof(int*));
		for(j = 0; j < size2; j++){
			out[i][j] = (int *) calloc(size3, sizeof(int));
		}
	}
	return out;
}

hit_struct * newHit_struct(int query_length, int reference_length){
	int i;
	/*Struct to store hit information*/	
	hit_struct* hit = (hit_struct*) calloc(1,sizeof(hit_struct));
	/* Initialize the hit / alignment constructs for this sequence*/
	for(i = 0; i < 3; i++){
		//hit->alignment[i] = (char *)calloc(query_length + reference_length, sizeof(char));
		hit->alignment[i] = (char *)calloc(MAX_QUERY_LENGHT, sizeof(char));
	}
	for(i = 0; i < 6; i++){
		hit->rest[i] = (char *)calloc(query_length + 10, sizeof(char));		
	}
	return hit;
}

void freeHit_struct(hit_struct * hit){
	int i;
	for(i = 0; i < 3; i++){
		free(hit->alignment[i]);
	}
	for(i = 0; i < 6; i++){
		free(hit->rest[i]);
	}
}

void free4dInt(int **** in, int size1, int size2, int size3){
	int i, j, k;	
	for(i = 0; i< size1; i++){
		for(j = 0; j< size2; j++){
			for(k = 0; k < size3; k++){
				free(in[i][j][k]);
			}
			free(in[i][j]);
		}
		free(in[i]);
	}
	free(in);
}

int **** new4dInt(int size1, int size2, int size3, int size4){
	int **** out;
	int i, j, k;
	out = (int****) calloc(size1, sizeof(int***));
	for(i = 0; i < size1; i++){
		out[i] = (int ***) calloc(size2, sizeof(int**));
		for(j = 0; j < size2; j++){
			out[i][j] = (int **) calloc(size3, sizeof(int*));
			for(k = 0; k < size3; k++){
				out[i][j][k] = (int *) calloc(size4, sizeof(int));
			}
		}
	}
	return out;
}

void free2dScore_struct(int4 ** in, int size1){
	int i;
	for(i = 0; i < size1; i++){
		free(in[i]);
	}
	free(in);
}
int4 ** new2dScore_struct(int size1, int* size2, int factor){
	int i;
	int4** scores = (int4**)calloc(size1, sizeof(int4*));
	for(i = 0; i < size1; i++){
		scores[i] = (int4*)calloc(size2[i]*factor, sizeof(int4));
	}
	return scores;
}

void clear2dScore_struct(int4** scores, int size1, int* size2, int factor){
	int i;	 
	for(i = 0; i < size1; i++){
		memset(scores[i], 0, size2[i]*factor*sizeof(int4));
	}	
}

void free3dScore_struct(int4 *** in, int size1, int size2){
	int i, j;
	for(i = 0; i < size1; i++){
		for(j = 0; j < size2; j++){
			free(in[i][j]);
		}
		free(in[i]);
	}
	free(in);
}
int4 *** new3dScore_struct(int size1, int size2, int* size3, int factor){
	int i, j;
	int4*** scores = (int4***)calloc(size1, sizeof(int4**));
	for(i = 0; i < size1; i++){
		scores[i] = (int4**)calloc(size2, sizeof(int4*));
		for(j = 0; j < size2; j++){
			scores[i][j] = (int4*)calloc(size3[j]*factor, sizeof(int4));
		}
	}
	return scores;
}

void clear3dScore_struct(int4*** scores, int size1, int size2, int* size3, int factor){
	int i, j;	 
	for(i = 0; i < size1; i++){
		for(j = 0; j < size2; j++){
			memset(scores[i][j], 0, size3[j]*factor*sizeof(int4));
		}
	}	
}



void free4dTrack(int **** in, int size1, int size2, int* size3){
	int i, j, k;	
	for(i = 0; i< size1; i++){
		for(j = 0; j< size2; j++){
			for(k = 0; k < size3[i] + 1; k++){
				free(in[i][j][k]);
			}
			free(in[i][j]);
		}
		free(in[i]);
	}
	free(in);
}

int **** new4dTrack(int size1, int size2, int* size3, int size4){
	int **** out;
	int i, j, k;
	out = (int****) calloc(size1, sizeof(int***));
	for(i = 0; i < size1; i++){
		out[i] = (int ***) calloc(size2, sizeof(int**));
		for(j = 0; j < size2; j++){
			out[i][j] = (int **) calloc(size3[i] + 1, sizeof(int*));
			for(k = 0; k < size3[i] + 1; k++){
				out[i][j][k] = (int *) calloc(size4 + 1, sizeof(int));
			}
		}
	}
	return out;
}

void clear4dTrack(int **** out, int size1, int size2, int* size3, int size4){	
	int i, j, k;	
	for(i = 0; i < size1; i++){		
		for(j = 0; j < size2; j++){			
			for(k = 0; k < size3[i] + 1; k++){
				memset(out[i][j][k], 0, (size4 + 1)*sizeof(int));
			}
		}
	}	
}

void free5dTrack(char4 **** in, int size1, int size2, int* size4){
	int i, j, k;	
	for(i = 0; i< size1; i++){
		for(j = 0; j< size2; j++){			
			for(k = 0; k < size4[j] + 1; k++){
				free(in[i][j][k]);
			}
			free(in[i][j]);			
		}
		free(in[i]);
	}
	free(in);
}

char4 **** new5dTrack(int size1, int size2, int* size4, int size5){
	char4 **** out;
	int i, j, k;
	out = (char4****) calloc(size1, sizeof(char4***));
	for(i = 0; i < size1; i++){
		out[i] = (char4 ***) calloc(size2, sizeof(char4**));
		for(j = 0; j < size2; j++){			
			out[i][j] = (char4 **) calloc(size4[j] + 1, sizeof(char4*));			
			for(k = 0; k < size4[j] + 1; k++){				
				out[i][j][k] = (char4 *) calloc(size5 + 1, sizeof(char4));
			}			
		}
	}
	return out;
}

void clear5dTrack(char4 **** out, int size1, int size2, int* size4, int size5){	
	int i, j, k;	
	for(i = 0; i < size1; i++){		
		for(j = 0; j < size2; j++){				
			for(k = 0; k < size4[j] + 1; k++){
				memset(out[i][j][k], 0, (size5 + 1)*sizeof(char4));
			}			
		}
	}	
}


void free3dTable_s(int *** in, int size1, int size2){
	int i, j;	
	for(i = 0; i< size1; i++){
		for(j = 0; j< size2 + 1; j++){
			free(in[i][j]);
		}
		free(in[i]);
	}
	free(in);
}

int *** new3dTable_s(int size1, int size2, int size3){
	int *** out;
	int i, j;
	out = (int***) calloc(size1, sizeof(int**));
	for(i = 0; i < size1; i++){
		out[i] = (int **) calloc(size2 + 1, sizeof(int*));
		for(j = 0; j < size2 + 1; j++){
			out[i][j] = (int *) calloc(size3 + 1, sizeof(int));
		}
	}
	return out;
}

void clear3dTable_s(int *** out, int size1, int size2, int size3){	
	int i, j;	
	for(i = 0; i < size1; i++){		
		for(j = 0; j < size2 + 1; j++){
			memset(out[i][j], 0, (size3 + 1)*sizeof(int));
		}
	}

}


void free3dTable(int *** in, int size1, int* size2){
	int i, j;	
	for(i = 0; i< size1; i++){
		for(j = 0; j< size2[i] + 1; j++){
			free(in[i][j]);
		}
		free(in[i]);
	}
	free(in);
}

int *** new3dTable(int size1, int* size2, int size3){
	int *** out;
	int i, j;
	out = (int***) calloc(size1, sizeof(int**));
	for(i = 0; i < size1; i++){
		out[i] = (int **) calloc(size2[i] + 1, sizeof(int*));
		for(j = 0; j < size2[i] + 1; j++){
			out[i][j] = (int *) calloc(size3 + 1, sizeof(int));
		}
	}
	return out;
}

void clear3dTable(int *** out, int size1, int* size2, int size3){	
	int i, j;	
	for(i = 0; i < size1; i++){		
		for(j = 0; j < size2[i] + 1; j++){
			memset(out[i][j], 0, (size3 + 1)*sizeof(int));
		}
	}

}

void free4dTable(int **** in, int size1, int size2, int* size3){
	int i, j, k;	
	for(i = 0; i< size1; i++){
		for(j = 0; j< size2; j++){
			for(k = 0; k< size3[j] + 1; k++){
				free(in[i][j][k]);
			}
		}
		free(in[i]);
	}
	free(in);
}

int **** new4dTable(int size1, int size2, int* size3, int size4){
	int **** out;
	int i, j, k;
	out = (int****) calloc(size1, sizeof(int***));
	for(i = 0; i < size1; i++){
		out[i] = (int ***) calloc(size2, sizeof(int**));
		for(j = 0; j < size2; j++){
			out[i][j] = (int **) calloc(size3[j] + 1, sizeof(int*));
			for(k = 0; k < size3[j] + 1; k++){
				out[i][j][k] = (int *) calloc(size4 + 1, sizeof(int));
			}
		}
	}
	return out;
}

void clear4dTable(int **** out, int size1, int size2, int* size3, int size4){	
	int i, j, k;	
	for(i = 0; i < size1; i++){		
		for(j = 0; j < size2; j++){
			for(k = 0; k < size3[j] + 1; k++){
				memset(out[i][j][k], 0, (size4 + 1)*sizeof(int));
			}
		}
	}

}

int * new1dInt(int size1){	
	return (int*) calloc(size1, sizeof(int));
}

void clear1dInt(int* in, int size1){	
	memset(in, 0, size1*sizeof(int));
}

int ** initializeRestricted_Array(Reference_batch *reference_array, input_type *  query_array, pair_struct* pairs, int total_pairs){
	int ** restricted_Array; 
	int r, q;
	restricted_Array = new2dInt(reference_array->size, query_array->size);
	for(r = 0; r < reference_array->size; r++){
		for(q = 0; q < query_array->size; q++){					
			if (restricted && !find_pair(query_array->name[q], reference_array->name[r], total_pairs, pairs)) {
				restricted_Array[r][q] = 1;
			}
		}
	}
	return restricted_Array;
}

void memcpyCharWithFixedMaxLength(char* dist, char * source, int size, int patch_size){	
	if(size > patch_size){
		fprintf(stderr, "The current CUDA version does not support a query length %d > %d.\n", size - 1, patch_size - 1);
        exit(EXIT_FAILURE); /* indicate failure.*/
	}	
	memcpy(dist, source, size*sizeof(char));	
}

char * memcpyChar(char * source, int size){
	char * dist;
	dist = (char *) calloc(size, sizeof(char));
	memcpy(dist, source, size*sizeof(char));
	return dist;
}

int * memcpyInt(int * source, int size){
	int * dist;
	dist = (int *) calloc(size, sizeof(int));
	memcpy(dist, source, size*sizeof(int));
	return dist;
}

Ref_batch* getRef_batch(char * filename, input_type * query_array, pair_struct* pairs, int total_pairs){
	int i, j, r, k;
	intXY * sortedID;
	int count = 0;
	Ref_batch* ref_batch;
	input_type* reference_array;	
	reference_array = loadRef(filename); // array for storing all references;		
	reference_array->sequence1D_reshape = NULL;
	reference_array->sequence1D_reshape_raw = NULL;
	reference_array->sequence1D = NULL;
	ref_batch = (Ref_batch *) calloc(1, sizeof(Ref_batch)); 
	ref_batch->ref_batch_size = ref_batch_size;
	//ref_batch->ref_batch_size_patched = ((ref_batch_size - 1)/WARP_SIZE + 1)*WARP_SIZE;
	ref_batch->batch_length = (reference_array->size - 1)/(ref_batch_size*nGPUs) + 1;	
	ref_batch->overall_reference_size = reference_array->size;
	ref_batch->overall_seqlen_max = reference_array->seqlen_max;
	ref_batch->overall_seqlen_max_patched = ((ref_batch->overall_seqlen_max + 1 - 1)/WARP_SIZE + 1)*WARP_SIZE;
	ref_batch->overall_seqlen_min = reference_array->seqlen_min;
	ref_batch->overall_mean = reference_array->seqlen_mean;
	ref_batch->overall_sd = reference_array->seqlen_sd; // standard deviation; 
	int length_in_last_batch = (reference_array->size - (ref_batch->batch_length - 1)*ref_batch_size*nGPUs);
	int num_GPUS_in_last_batch = (length_in_last_batch - 1)/ref_batch_size + 1;
	int length_in_last_batch_and_last_GPU = (length_in_last_batch - (num_GPUS_in_last_batch - 1)*ref_batch_size);
	sortedID = (intXY *)malloc(reference_array->size*sizeof(intXY));
	
	if(ref_batch->ref_batch_size > MAX_RefBatchSize){
		printf("ref_batch_size_patched should be less than 64\n");
		exit(0);
	}

	for(i = 0; i < reference_array->size; i++){
		sortedID[i].val = reference_array->seqlen[i];
		sortedID[i].id = i;
	}
	qsort(sortedID, reference_array->size, sizeof(intXY), cmpsIntXY);
	
	ref_batch->reference_batch = (Reference_batch **) calloc(ref_batch->batch_length, sizeof(Reference_batch*));	
	for(i = 0; i < ref_batch->batch_length; i++){
		ref_batch->reference_batch[i] = (Reference_batch *) calloc(nGPUs, sizeof(Reference_batch));
		for(int g = 0; g < nGPUs; g++){
			if(i ==  (ref_batch->batch_length - 1)){
				if(g < num_GPUS_in_last_batch - 1){
					ref_batch->reference_batch[i][g].size = ref_batch->ref_batch_size;		
				}else if(g == num_GPUS_in_last_batch - 1){
					ref_batch->reference_batch[i][g].size = length_in_last_batch_and_last_GPU;
				}else{
					break;
				}
			}else{
				ref_batch->reference_batch[i][g].size = ref_batch->ref_batch_size;		
			}			
			ref_batch->reference_batch[i][g].sequence = (char **) calloc(ref_batch->reference_batch[i][g].size, sizeof(char*));
			ref_batch->reference_batch[i][g].description = (char **) calloc(ref_batch->reference_batch[i][g].size, sizeof(char*));
			ref_batch->reference_batch[i][g].name = (char **) calloc(ref_batch->reference_batch[i][g].size, sizeof(char*));
			ref_batch->reference_batch[i][g].seqlen = (int *) calloc(MAX_RefBatchSize, sizeof(int));		
			ref_batch->reference_batch[i][g].end_score = new2dFloatType(ref_batch->reference_batch[i][g].size, query_array->size);	
			ref_batch->reference_batch[i][g].hit_summary = new2dHitSummary(ref_batch->reference_batch[i][g].size, query_array->size);	 /*Summary of best hits for reporting*/	
			ref_batch->reference_batch[i][g].seqlen_max = 0;			
			for(j = 0; j < ref_batch->reference_batch[i][g].size; j++){			
				r = sortedID[count].id;	
				ref_batch->reference_batch[i][g].seqlen[j]      = reference_array->seqlen[r];
				if(ref_batch->reference_batch[i][g].seqlen_max < ref_batch->reference_batch[i][g].seqlen[j])
					ref_batch->reference_batch[i][g].seqlen_max = ref_batch->reference_batch[i][g].seqlen[j];				
				ref_batch->reference_batch[i][g].description[j] = memcpyChar(reference_array->description[r], (int)strlen(reference_array->description[r])+1);
				ref_batch->reference_batch[i][g].name[j]        = memcpyChar(reference_array->name[r], (int)strlen(reference_array->name[r])+1);
				count++;
			}
			ref_batch->reference_batch[i][g].seqlen_max_patched = ((ref_batch->reference_batch[i][g].seqlen_max + 1 - 1)/WARP_SIZE + 1)*WARP_SIZE; // + 1 for "\0";			
			ref_batch->reference_batch[i][g].seqlen_max_patched_char4 = ref_batch->reference_batch[i][g].seqlen_max_patched/4;			
			ref_batch->reference_batch[i][g].restricted_Array = initializeRestricted_Array(&(ref_batch->reference_batch[i][g]), query_array, pairs, total_pairs);
			//*/
		}
	}
	
	assert(count == reference_array->size);	
	count = 0;	
	for(i = 0; i < ref_batch->batch_length; i++){
		for(int g = 0; g < nGPUs; g++){
			if(ref_batch->reference_batch[i][g].size == 0){
				break;
			}
			ref_batch->reference_batch[i][g].sequence1D = (char *) calloc(ref_batch->reference_batch[i][g].size * ref_batch->reference_batch[i][g].seqlen_max_patched, sizeof(char));				
			ref_batch->reference_batch[i][g].sequence1D_revalued = (char4 *) calloc(ref_batch->reference_batch[i][g].size * ref_batch->reference_batch[i][g].seqlen_max_patched_char4, sizeof(char4));				
			ref_batch->reference_batch[i][g].sequence1D_revalued_raw = (char2 *) calloc(ref_batch->reference_batch[i][g].size * ref_batch->reference_batch[i][g].seqlen_max_patched, sizeof(char2));				
			for(j = 0; j < ref_batch->reference_batch[i][g].size; j++){
				r = sortedID[count].id;	
				ref_batch->reference_batch[i][g].sequence[j] = ref_batch->reference_batch[i][g].sequence1D + j * ref_batch->reference_batch[i][g].seqlen_max_patched;			
				memcpyCharWithFixedMaxLength(ref_batch->reference_batch[i][g].sequence[j], reference_array->sequence[r], reference_array->seqlen[r]+1, ref_batch->reference_batch[i][g].seqlen_max_patched);
				
				count++;
			}
			for(k = 0; k < ref_batch->reference_batch[i][g].size * ref_batch->reference_batch[i][g].seqlen_max_patched; k++){			
				ref_batch->reference_batch[i][g].sequence1D_revalued_raw[k].x = ref_batch->reference_batch[i][g].sequence1D[k];
				ref_batch->reference_batch[i][g].sequence1D_revalued_raw[k].y = bases[ref_batch->reference_batch[i][g].sequence1D[k]];
			}
			for(k = 0; k < ref_batch->reference_batch[i][g].size * ref_batch->reference_batch[i][g].seqlen_max_patched/4; k++){							
					ref_batch->reference_batch[i][g].sequence1D_revalued[k].x = bases[ref_batch->reference_batch[i][g].sequence1D[k*4]];								
					ref_batch->reference_batch[i][g].sequence1D_revalued[k].y = bases[ref_batch->reference_batch[i][g].sequence1D[k*4 + 1]];
					ref_batch->reference_batch[i][g].sequence1D_revalued[k].z = bases[ref_batch->reference_batch[i][g].sequence1D[k*4 + 2]];
					ref_batch->reference_batch[i][g].sequence1D_revalued[k].w = bases[ref_batch->reference_batch[i][g].sequence1D[k*4 + 3]];				
			}
		}
	}
	
	assert(count == reference_array->size);	
	free(sortedID);	
	freeInputType(reference_array);	
	return ref_batch;
}

void freeReference_batch_partal(Reference_batch * in, int query_array_size){
	int i = 0;
	free(in->sequence1D);
	free(in->sequence1D_revalued);
	free(in->sequence1D_revalued_raw);
	for (i = 0; i < in->size; ++i) {
		//free(in->sequence[i]);
		free(in->description[i]);
		free(in->name[i]);
	}
	free(in->sequence);
	free(in->description);
	free(in->name);
	free(in->seqlen);
	free2dHitSummary(in->hit_summary,in->size, query_array_size);	
	free2dInt(in->restricted_Array, in->size);
	free2dFloatType(in->end_score, in->size);
	//free(in);
}
void freeRef_batch(Ref_batch* ref_batch, int query_array_size){
	int i;
	for(i = 0; i < ref_batch->batch_length; i++){
		for(int g = 0; g < nGPUs; g++){
			freeReference_batch_partal(&ref_batch->reference_batch[i][g], query_array_size);
		}
	}
	free(ref_batch->reference_batch);
	free(ref_batch);
}


void recalloc3dScore_struct(int4*** scores, int size1, int size2, int* size3, int factor){
	int i, j;	 
	for(i = 0; i < size1; i++){
		for(j = 0; j < size2; j++){			
			free(scores[i][j]);	scores[i][j] = (int4 *)calloc(size3[j]*factor, sizeof(int4));	
			
		}
	}	
}

void recalloc5dTrack(char4 ***** out, int size1, int size2, int size3, int* size4, int size5){	
	int i, j, k, l;	
	for(i = 0; i < size1; i++){		
		for(j = 0; j < size2; j++){	
			for(k = 0; k < size3; k++){
				for(l = 0; l < size4[j] + 1; l++){					
					free(out[i][j][k][l]);	out[i][j][k][l] = (char4 *)calloc(size5 + 1, sizeof(char4));					
				}
			}
		}
	}	
}

void recalloc4dTable(int **** out, int size1, int size2, int* size3, int size4){	
	int i, j, k;	
	for(i = 0; i < size1; i++){		
		for(j = 0; j < size2; j++){
			for(k = 0; k < size3[j] + 1; k++){				
				free(out[i][j][k]);	out[i][j][k] = (int *)calloc(size4 + 1, sizeof(int));	
			}
		}
	}

}