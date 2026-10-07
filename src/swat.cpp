/* swat.c */
/* -------------------------------------------------------------------
 * CUDA miRanda- An miRNA target scanner with CUDA speedup
  *
 * Copyright (C) (2012) UCSD DBMI Shuang Wang
 *
 * Distributed under the GNU Public License (GPL)
 * See the files 'COPYING' and 'LICENSE' for details
 *
 * Authors: Shuang Wang, */
//#pragma warning(disable:4244)
#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <limits.h>
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

/* Master list of Nucleotide Bases*/
int bases[256];

void initialize_bases() {
	int i;
	for (i = 0; i < 256; i++) bases[i] = INT_MAX;
	bases['c'] = 0;
	bases['C'] = 0;
	bases['g'] = 1;
	bases['G'] = 1;
	bases['a'] = 2;
	bases['A'] = 2;
	bases['u'] = 3;
	bases['U'] = 3;
	bases['t'] = 3;
	bases['T'] = 3;
	bases['n'] = 4;
	bases['N'] = 4;
	bases['r'] = 4;
	bases['R'] = 4;
	bases['y'] = 4;
	bases['Y'] = 4;
	bases['k'] = 4;
	bases['K'] = 4;
	bases['m'] = 4;
	bases['M'] = 4;
	bases['s'] = 4;
	bases['S'] = 4;
	bases['w'] = 4;
	bases['W'] = 4;
	bases['b'] = 4;
	bases['B'] = 4;
	bases['d'] = 4;
	bases['D'] = 4;
	bases['h'] = 4;
	bases['H'] = 4;
	bases['v'] = 4;
	bases['V'] = 4;
	bases['x'] = 4;
	bases['X'] = 4;
}

const int ALIGNMENT_MATCH_SCORE = 5;
const int ALIGNMENT_MISMATCH_SCORE = -3;
const int ALIGNMENT_WOBBLE_SCORE = 1;
const int ALIGNMENT_MASKED_SCORE = -1;
const char ALIGNMENT_MATCH_REPRESENTATION = '|';
const char ALIGNMENT_MISMATCH_REPRESENTATION = ' ';
const char ALIGNMENT_WOBBLE_REPRESENTATION = ':';
const char ALIGNMENT_MASKED_REPRESENTATION = ' ';

/* Scoring Matrix - Only Takes Normal Nucleotides*/

int ali_score[32]; //[5][5];
int ali_score5p[32]; //[5][5];
char ali_rep[32]; //[5][5];

int bases_are_complementary(char nt1, char nt2) {
	if (nt1 == 'C') return nt2 == 'G';
	if (nt1 == 'G') return nt2 == 'C';
	if (nt1 == 'A') return nt2 == 'U';
	if (nt1 == 'U') return nt2 == 'A';
	return 0;
}

int bases_are_wobble_pair(char nt1, char nt2) {
	if (nt1 == 'G') return nt2 == 'U';
	if (nt1 == 'U') return nt2 == 'G';
	return 0;
}

void initialize_scores() {
	char scored_bases[5] = {'C','G','A','U','X'};
	int i;
	int j;
	for (i = 0; i < sizeof(scored_bases) / sizeof(char); i++) {
		for (j = 0; j < sizeof(scored_bases) / sizeof(char); j++) {
			if (scored_bases[i] == 'X' || scored_bases[j] == 'X') {
				ali_rep[i * ali_patch + j] = ALIGNMENT_MASKED_REPRESENTATION;
				ali_score[i * ali_patch + j] = ALIGNMENT_MASKED_SCORE;
				ali_score5p[i * ali_patch + j] = (int)(scale * ALIGNMENT_MASKED_SCORE);
				continue;
			}
			if (bases_are_complementary(scored_bases[i],scored_bases[j])) {
				ali_rep[i * ali_patch + j] = ALIGNMENT_MATCH_REPRESENTATION;
				ali_score[i * ali_patch + j] = ALIGNMENT_MATCH_SCORE;
				ali_score5p[i * ali_patch + j] = (int)(scale * ALIGNMENT_MATCH_SCORE);
				continue;
			}
			if (bases_are_wobble_pair(scored_bases[i],scored_bases[j])) {
				ali_rep[i * ali_patch + j] = ALIGNMENT_WOBBLE_REPRESENTATION;
				ali_score[i * ali_patch + j] = ALIGNMENT_WOBBLE_SCORE;
				ali_score5p[i * ali_patch + j] = (int)(scale * ALIGNMENT_WOBBLE_SCORE);
				continue;
			}
			ali_rep[i * ali_patch + j] = ALIGNMENT_MISMATCH_REPRESENTATION;
			ali_score[i * ali_patch + j] = ALIGNMENT_MISMATCH_SCORE;
			ali_score5p[i * ali_patch + j] = (int)(scale * ALIGNMENT_MISMATCH_SCORE);
		}
	}
}

inline int score(char nt1, char nt2) {
	//return ali_score[bases[(int)nt1] * ali_patch + bases[(int)nt2]];
	return ali_score[nt1 * ali_patch + nt2];
}

inline int score5p(char nt1, char nt2) {
	//return ali_score5p[bases[(int)nt1] * ali_patch + bases[(int)nt2]];
	return ali_score5p[nt1 * ali_patch + nt2];
}


void nullify_overlaps(int total_elements, ushort4* scores) {
	int index;
	int outer_ref_end;
	int inner_index;
	for (index = 0; index < total_elements; index++) {
		if (scores[index].x == 0) {continue;}		
		/* use query end to figure where ref offset relative to fictional query offset 0 would be */
		outer_ref_end = scores[index].w - scores[index].z;
		for (inner_index = total_elements; inner_index > index; inner_index--) {
			int inner_ref_end = scores[inner_index].w -
					scores[inner_index].z;
			int diff;
			diff = abs(inner_ref_end - outer_ref_end);
			if (diff <= OVERLAP_LIMIT) {
				if(scores[index].x < scores[inner_index].x){  // remover small score
					scores[index].x = 0;					
					break;
				}else if(scores[index].x > scores[inner_index].x){ // keep large score
					scores[inner_index].x = 0;
				}else{
					if(scores[inner_index].w < scores[index].w){ // for the same score, keep the earlier one;
						scores[index].x = 0;						
						break;
					}else{
						scores[inner_index].x = 0;
					}
				}
			}
		}
	}
}

void old_nullify_overlaps(int total_elements, ushort4* scores, int *scores_length) {
	int index;
	int outer_ref_end;
	int inner_index;
	for (index = 0; index <= total_elements; index++) {
		if (scores[index].x == 0) {continue;}
		(*scores_length)++;
		/* use query end to figure where ref offset relative to fictional query offset 0 would be */
		outer_ref_end = scores[index].w - scores[index].z;
		for (inner_index = total_elements; inner_index > index; inner_index--) {
			int inner_ref_end = scores[inner_index].w -
					scores[inner_index].z;
			int diff;
			diff = abs(inner_ref_end - outer_ref_end);
			if (diff <= OVERLAP_LIMIT) {
				scores[inner_index].x = 0;
			}
		}
	}
}

inline void max_finder_and_track_twostates_with_0_bit(int zero, int a, int* max , signed char * track){	
	*max = a - (((a - zero) >> (INTBITS32 - 1))&(a - zero));
	*track = (signed char) (3 - *max + (*max&~3));
	*max >>= 2;
}

inline void max_finder_and_track_fourstates_bit(int a, int b, int c, int* max, signed char * track, int tmp_score){
	*max = a - ((( a - b) >> (INTBITS32 - 1))&(a - b));
	*max = *max - (((*max - c) >> (INTBITS32 - 1))&(*max - c));
	*track =(signed char)( 3 - *max + (*max&~3));
	*max >>= 2;
	max_finder_and_track_twostates_with_0_bit(3, (((*max + tmp_score)<<2)+ 3 - *track), max, track);
}

inline void max_finder_and_track_fourstates_norm(int a, int b, int c, int * max, signed char *track, int tmp_score) {	
	if (a >= b && a >= c) {
		*max = a;
		*track = 1;
	} else if (b >= c) {
		*max = b;
		*track = 2;
	} else {
		*max = c;
		*track = 3;
	}	
	*max += tmp_score;
	if(*max <= 0){
		*max = 0;
		*track = 0;
	}
}

inline void max_finder_and_track_fourstates(int a, int b, int c, int * max, signed char *track, int tmp_score) {	
	max_finder_and_track_fourstates_bit((a<<2)+2, (b<<2)+1, (c<<2), max, track, tmp_score);
	//max_finder_and_track_fourstates_norm(a, b, c, max, track, tmp_score);
}

inline void max_finder_and_track_twostates_bit(int a, int b, int* max, signed char * track){
	*max = a - ((( a - b) >> (INTBITS32 - 1))&(a - b));	
	*track =(signed char)( 3 - *max + (*max&~3));
	*max >>= 2;	
}


inline void max_finder_and_track_twostates_norm(int tmp_a, int tmp_b, int * max, signed char *track){
	if (tmp_b > tmp_a) {
		*max = tmp_b;
		*track = 2;
	} else {
		*max = tmp_a;
		*track = 1;
	}
}

inline void max_finder_and_track_twostates_ab(int a, int b, int * max, signed char *track) {	
	max_finder_and_track_twostates_bit((a<<2)+2, (b<<2)+1, max, track);
	//max_finder_and_track_twostates_norm(a, b, max, track);
}

inline void max_finder_and_track_twostates_ac(int a, int c, int * max, signed char *track) {	
	max_finder_and_track_twostates_bit((a<<2)+2, (c<<2), max, track);
	//max_finder_and_track_twostates_norm(a, c, max, track);
}

inline void max_finder_and_track_threestates_bit(int zero, int a, int b, int c, int* max, signed char * track){
	*max = a - ((( a - b) >> (INTBITS32 - 1))&(a - b));
	*max = *max - (((*max - c) >> (INTBITS32 - 1))&(*max - c));
	*max = *max - (((*max - zero) >> (INTBITS32 - 1))&(*max - zero));
	*track = (signed char) (3 - *max + (*max&~3));
	*max >>= 2;
}

inline void max_finder_and_track_threestates_norm(int a, int b, int c, int * max, signed char *track) {	
	if (a >= b && a >= c) {
		*max = a;
		*track = 1;
	} else if (b >= c) {
		*max = b;
		*track = 2;
	} else {
		*max = c;
		*track = 3;
	}
	if (*max <= 0) {
		*max = 0;
		*track = 0;
	}	
}

inline void max_finder_and_track_threestates(int a, int b, int c, int* max, signed char *track) {
	max_finder_and_track_threestates_bit(3, (a<<2)+2, (b<<2)+1, (c<<2), max, track);
	//max_finder_and_track_threestates_norm(a, b, c, max, track);
}

inline void memcpy_d(int query_array_size_patched, int **a_nt_nt_dist, int**a_nt_nt_src){
	for(int i = 0; i < MAX_QUERY_LENGHT; i++){
		memcpy(a_nt_nt_dist[i],a_nt_nt_src[i], sizeof(int) * query_array_size_patched);
	}
}


void build_matrix(int shift_offset_ref, input_type * query_array, char* ref_seq, int* best, char4* track,
		          int reference_length, int* length_3p_for_weighting_array, char3 *open_ext_in_seed, ushort4* scores, int *good_count) {
	int** a_nt_nt0  = new2dInt(MAX_QUERY_LENGHT, query_array->size_patched);
	int** b_gap_nt0 = new2dInt(MAX_QUERY_LENGHT, query_array->size_patched);
	int** c_nt_gap0 = new2dInt(MAX_QUERY_LENGHT, query_array->size_patched);
	int** a_nt_nt1  = new2dInt(MAX_QUERY_LENGHT, query_array->size_patched);
	int** b_gap_nt1 = new2dInt(MAX_QUERY_LENGHT, query_array->size_patched);
	int** c_nt_gap1 = new2dInt(MAX_QUERY_LENGHT, query_array->size_patched);
	int** a_nt_nt2  = new2dInt(MAX_QUERY_LENGHT, query_array->size_patched);
	int** b_gap_nt2 = new2dInt(MAX_QUERY_LENGHT, query_array->size_patched);
	int** c_nt_gap2 = new2dInt(MAX_QUERY_LENGHT, query_array->size_patched);		
	for (int j = 2; j <= (reference_length + MAX_QUERY_LENGHT); j++) {	
		for (int i = 1; i <= MAX_QUERY_LENGHT; i++) {	
				for(int q = 0; q < query_array->size; q++){
					//if( i > query_array->seqlen[q]) continue;
					//if( j > (reference_length + query_array->seqlen[q])) continue;
					if( i <= query_array->seqlen[q] &&
						j <= (reference_length + query_array->seqlen[q]) &&
						i < j &&
						!((j > reference_length) && (i < (j - reference_length)))){
							int shift_ref_id = j * shift_offset_ref;

							//if(i >= j) break;
							//if((j > reference_length) && (i < (j - reference_length))) continue;	
							int shift_query_id = shift_ref_id + i*query_array->size_patched;
							int tmp_score = 0;							
							if(i > 2){	
								char query_seq_val = 0;
								//int query_seq_val = query_array->sequence1D_reshape[(i-1)*query_array->size_patched + q];
								int pp = (i - 1)/4;
								switch((i-1)%4){
								case 0:
										query_seq_val = query_array->sequence1D_reshape[pp*query_array->size_patched + q].x;
										break;
								case 1:
										query_seq_val = query_array->sequence1D_reshape[pp*query_array->size_patched + q].y;
										break;
								case 2:
										query_seq_val = query_array->sequence1D_reshape[pp*query_array->size_patched + q].z;
										break;
								case 3:
										query_seq_val = query_array->sequence1D_reshape[pp*query_array->size_patched + q].w;
										break;
								}
								if(i <= length_3p_for_weighting_array[q]){
									tmp_score = score(query_seq_val, ref_seq[j - i - 1]);
								}else if(i < query_array->seqlen[q]){														
									tmp_score = score5p(query_seq_val, ref_seq[j - i - 1]);
								}
							}			
							int shift_id = i * query_array->size_patched;
							//max_finder_and_track_fourstates(a_nt_nt[i-1][j-1], b_gap_nt[i-1][j-1], c_nt_gap[i-1][j-1], &(a_nt_nt[i][j]), &(track[1][i][j]), tmp_score);										
							//max_finder_and_track_fourstates(a_nt_nt0[i-1][q], b_gap_nt0[i-1][q], c_nt_gap0[i-1][q], &(a_nt_nt2[i][q]), &(track[q][i][j - i].y), tmp_score);						
							max_finder_and_track_fourstates(a_nt_nt0[i-1][q], b_gap_nt0[i-1][q], c_nt_gap0[i-1][q], &(a_nt_nt2[i][q]), &(track[shift_query_id + q].y), tmp_score);						
							//max_finder_and_track_twostates_ab((a_nt_nt[i][j-1] + open), (b_gap_nt[i][j-1] + ext), &b_gap_nt[i][j], &track[2][i][j]);				
							//max_finder_and_track_twostates_ab((a_nt_nt1[i][q] + open[q + shift_id]), (b_gap_nt1[i][q] + ext[q + shift_id]), &b_gap_nt2[i][q], &track[q][i][j - i].z);
							max_finder_and_track_twostates_ab((a_nt_nt1[i][q] + open_ext_in_seed[q + shift_id].x), (b_gap_nt1[i][q] + open_ext_in_seed[q + shift_id].y), &b_gap_nt2[i][q], &track[shift_query_id + q].z);
							if (!open_ext_in_seed[q + shift_id].z) {
								//max_finder_and_track_twostates_ac((a_nt_nt[i-1][j] + open), (c_nt_gap[i-1][j] + ext), &c_nt_gap[i][j], &track[3][i][j]);				
								//max_finder_and_track_twostates_ac((a_nt_nt1[i-1][q] + open[q + shift_id]), (c_nt_gap1[i-1][q] + ext[q + shift_id]), &c_nt_gap2[i][q], &track[q][i][j - i].w);				
								max_finder_and_track_twostates_ac((a_nt_nt1[i-1][q] + open_ext_in_seed[q + shift_id].x), (c_nt_gap1[i-1][q] + open_ext_in_seed[q + shift_id].y), &c_nt_gap2[i][q], &track[shift_query_id + q].w);				
							} else {
								//c_nt_gap[i][j] = -1;
								//track[3][i][j] = 1;
								c_nt_gap2[i][q] = -1;
								//track[q][i][j - i].w = 1;
								track[shift_query_id + q].w = 1;
							}
							//max_finder_and_track_threestates(a_nt_nt[i][j], b_gap_nt[i][j], c_nt_gap[i][j], &(best[i][j]), &(track[0][i][j]));				
							//max_finder_and_track_threestates(a_nt_nt2[i][q], b_gap_nt2[i][q], c_nt_gap2[i][q], &(best[q][i][j - i]), &(track[q][i][j - i].x));
							max_finder_and_track_threestates(a_nt_nt2[i][q], b_gap_nt2[i][q], c_nt_gap2[i][q], &(best[shift_query_id + q]), &(track[shift_query_id + q].x));

							if (best[shift_query_id + q] >= score_threshold && track[shift_query_id + q].x != 2 && good_count[q] < MAX_REPORTS_OF_FIRST_HIT) {				
								// reject alignments ending with an unaligned reference base				
								int shift_id2 = good_count[q] * query_array->size_patched;// q * MAX_REPORTS_OF_FIRST_HIT;
								scores[q + shift_id2].x = best[shift_query_id + q];
								scores[q + shift_id2].y = good_count[q];
								scores[q + shift_id2].z = i;
								scores[q + shift_id2].w = j - i;
								good_count[q]++;				
							}
					}
				}// for(int q = 0; q < query_array->size; q++){				
		} // for i
		memcpy_d(query_array->size_patched, a_nt_nt0, a_nt_nt1);memcpy_d(query_array->size_patched, a_nt_nt1, a_nt_nt2);
		memcpy_d(query_array->size_patched, b_gap_nt0, b_gap_nt1);memcpy_d(query_array->size_patched, b_gap_nt1, b_gap_nt2);
		memcpy_d(query_array->size_patched, c_nt_gap0, c_nt_gap1);memcpy_d(query_array->size_patched, c_nt_gap1, c_nt_gap2);	
	} // for (int j = 1; j <= reference_length; j++) {	

	free2dInt(a_nt_nt0,  MAX_QUERY_LENGHT);
	free2dInt(b_gap_nt0, MAX_QUERY_LENGHT); 
	free2dInt(c_nt_gap0, MAX_QUERY_LENGHT);
	free2dInt(a_nt_nt1,  MAX_QUERY_LENGHT);
	free2dInt(b_gap_nt1, MAX_QUERY_LENGHT); 
	free2dInt(c_nt_gap1, MAX_QUERY_LENGHT);
	free2dInt(a_nt_nt2,  MAX_QUERY_LENGHT);
	free2dInt(b_gap_nt2, MAX_QUERY_LENGHT); 
	free2dInt(c_nt_gap2, MAX_QUERY_LENGHT);
}
/* good
void build_matrix(int shift_offset_ref, int shift_offset_score, input_type * query_array, char* ref_seq, int* best, char4* track,
		          int reference_length, int* length_3p_for_weighting_array, int * open, int* ext, int* in_seed, int4* scores, int *good_count) {
	int** a_nt_nt0  = new2dInt(MAX_QUERY_LENGHT, query_array->size);
	int** b_gap_nt0 = new2dInt(MAX_QUERY_LENGHT, query_array->size);
	int** c_nt_gap0 = new2dInt(MAX_QUERY_LENGHT, query_array->size);
	int** a_nt_nt1  = new2dInt(MAX_QUERY_LENGHT, query_array->size);
	int** b_gap_nt1 = new2dInt(MAX_QUERY_LENGHT, query_array->size);
	int** c_nt_gap1 = new2dInt(MAX_QUERY_LENGHT, query_array->size);
	int** a_nt_nt2  = new2dInt(MAX_QUERY_LENGHT, query_array->size);
	int** b_gap_nt2 = new2dInt(MAX_QUERY_LENGHT, query_array->size);
	int** c_nt_gap2 = new2dInt(MAX_QUERY_LENGHT, query_array->size);

	for(int q = 0; q < query_array->size; q++){	
		int tmp_score;			
		for (int j = 2; j <= (reference_length + MAX_QUERY_LENGHT); j++) {	
			if( j > (reference_length + query_array->seqlen[q])) continue;
			int shift_ref_id = j * shift_offset_ref;
			for (int i = 1; i <= query_array->seqlen[q]; i++) {	
				if(i >= j) break;
				if((j > reference_length) && (i < (j - reference_length))) continue;	
				int shift_query_id = shift_ref_id + i*query_array->size_patched;
				if(i <= 2 || i == query_array->seqlen[q]){
					tmp_score = 0;
				}else{
					if(i <= length_3p_for_weighting_array[q]){
						tmp_score = score(query_array->sequence1D_reshape[(i-1)*query_array->size_patched + q], ref_seq[j - i - 1]);
					}else{
						tmp_score = score5p(query_array->sequence1D_reshape[(i-1)*query_array->size_patched + q], ref_seq[j - i - 1]);
					}
				}				
				int shift_id = i * query_array->size_patched;
				//max_finder_and_track_fourstates(a_nt_nt[i-1][j-1], b_gap_nt[i-1][j-1], c_nt_gap[i-1][j-1], &(a_nt_nt[i][j]), &(track[1][i][j]), tmp_score);										
				//max_finder_and_track_fourstates(a_nt_nt0[i-1][q], b_gap_nt0[i-1][q], c_nt_gap0[i-1][q], &(a_nt_nt2[i][q]), &(track[q][i][j - i].y), tmp_score);						
				max_finder_and_track_fourstates(a_nt_nt0[i-1][q], b_gap_nt0[i-1][q], c_nt_gap0[i-1][q], &(a_nt_nt2[i][q]), &(track[shift_query_id + q].y), tmp_score);						
				//max_finder_and_track_twostates_ab((a_nt_nt[i][j-1] + open), (b_gap_nt[i][j-1] + ext), &b_gap_nt[i][j], &track[2][i][j]);				
				//max_finder_and_track_twostates_ab((a_nt_nt1[i][q] + open[q + shift_id]), (b_gap_nt1[i][q] + ext[q + shift_id]), &b_gap_nt2[i][q], &track[q][i][j - i].z);
				max_finder_and_track_twostates_ab((a_nt_nt1[i][q] + open[q + shift_id]), (b_gap_nt1[i][q] + ext[q + shift_id]), &b_gap_nt2[i][q], &track[shift_query_id + q].z);
				if (!in_seed[q + shift_id]) {
					//max_finder_and_track_twostates_ac((a_nt_nt[i-1][j] + open), (c_nt_gap[i-1][j] + ext), &c_nt_gap[i][j], &track[3][i][j]);				
					//max_finder_and_track_twostates_ac((a_nt_nt1[i-1][q] + open[q + shift_id]), (c_nt_gap1[i-1][q] + ext[q + shift_id]), &c_nt_gap2[i][q], &track[q][i][j - i].w);				
					max_finder_and_track_twostates_ac((a_nt_nt1[i-1][q] + open[q + shift_id]), (c_nt_gap1[i-1][q] + ext[q + shift_id]), &c_nt_gap2[i][q], &track[shift_query_id + q].w);				
				} else {
					//c_nt_gap[i][j] = -1;
					//track[3][i][j] = 1;
					c_nt_gap2[i][q] = -1;
					//track[q][i][j - i].w = 1;
					track[shift_query_id + q].w = 1;
				}
				//max_finder_and_track_threestates(a_nt_nt[i][j], b_gap_nt[i][j], c_nt_gap[i][j], &(best[i][j]), &(track[0][i][j]));				
				//max_finder_and_track_threestates(a_nt_nt2[i][q], b_gap_nt2[i][q], c_nt_gap2[i][q], &(best[q][i][j - i]), &(track[q][i][j - i].x));
				max_finder_and_track_threestates(a_nt_nt2[i][q], b_gap_nt2[i][q], c_nt_gap2[i][q], &(best[shift_query_id + q]), &(track[shift_query_id + q].x));

				if (best[shift_query_id + q] >= score_threshold && track[shift_query_id + q].x != 2) {				
				// reject alignments ending with an unaligned reference base				
					int shift_id2 = q * shift_offset_score;
					scores[good_count[q] + shift_id2].x = best[shift_query_id + q];
					scores[good_count[q] + shift_id2].y = good_count[q];
					scores[good_count[q] + shift_id2].z = i;
					scores[good_count[q] + shift_id2].w = j - i;
					good_count[q]++;				
				}
			}
			memcpy_d(q, a_nt_nt0, a_nt_nt1);memcpy_d(q, a_nt_nt1, a_nt_nt2);
			memcpy_d(q, b_gap_nt0, b_gap_nt1);memcpy_d(q, b_gap_nt1, b_gap_nt2);
			memcpy_d(q, c_nt_gap0, c_nt_gap1);memcpy_d(q, c_nt_gap1, c_nt_gap2);			
		}	// for (int j = 1; j <= reference_length; j++) {			
	} // for(int q = 0; q < query_array->size; q++){	

	free2dInt(a_nt_nt0,  MAX_QUERY_LENGHT);
	free2dInt(b_gap_nt0, MAX_QUERY_LENGHT); 
	free2dInt(c_nt_gap0, MAX_QUERY_LENGHT);
	free2dInt(a_nt_nt1,  MAX_QUERY_LENGHT);
	free2dInt(b_gap_nt1, MAX_QUERY_LENGHT); 
	free2dInt(c_nt_gap1, MAX_QUERY_LENGHT);
	free2dInt(a_nt_nt2,  MAX_QUERY_LENGHT);
	free2dInt(b_gap_nt2, MAX_QUERY_LENGHT); 
	free2dInt(c_nt_gap2, MAX_QUERY_LENGHT);
}*/
/*
void build_matrix(input_type * query_array, char* ref_seq, int*** best, char4*** track,
		          int reference_length, int* length_3p_for_weighting_array, int ** open, int** ext, int** in_seed) {
	for(int q = 0; q < query_array->size; q++){	
		int good_count = 0;		
		int** a_nt_nt  = new2dInt(2, reference_length + 1);
		int** b_gap_nt = new2dInt(2, reference_length + 1);
		int** c_nt_gap = new2dInt(2, reference_length + 1);
		int tmp_score;		
		for (int i = 1; i <= query_array->seqlen[q]; i++) {			
			for (int j = 1; j <= reference_length; j++) {	
				if(i <= 2 || i == query_array->seqlen[q]){
					tmp_score = 0;
				}else{
					if(i <= length_3p_for_weighting_array[q]){
						tmp_score = score(query_array->sequence[q][i-1], ref_seq[j - 1]);
					}else{
						tmp_score = score5p(query_array->sequence[q][i-1], ref_seq[j - 1]);
					}
				}			
				//max_finder_and_track_fourstates(a_nt_nt[i-1][j-1], b_gap_nt[i-1][j-1], c_nt_gap[i-1][j-1], &(a_nt_nt[i][j]), &(track[1][i][j]), tmp_score);						
				max_finder_and_track_fourstates(a_nt_nt[0][j-1], b_gap_nt[0][j-1], c_nt_gap[0][j-1], &(a_nt_nt[1][j]), &(track[q][i][j].y), tmp_score);						
				//max_finder_and_track_twostates_ab((a_nt_nt[i][j-1] + open), (b_gap_nt[i][j-1] + ext), &b_gap_nt[i][j], &track[2][i][j]);
				max_finder_and_track_twostates_ab((a_nt_nt[1][j-1] + open[q][i]), (b_gap_nt[1][j-1] + ext[q][i]), &b_gap_nt[1][j], &track[q][i][j].z);
				if (!in_seed[q][i]) {
					//max_finder_and_track_twostates_ac((a_nt_nt[i-1][j] + open), (c_nt_gap[i-1][j] + ext), &c_nt_gap[i][j], &track[3][i][j]);				
					max_finder_and_track_twostates_ac((a_nt_nt[0][j] + open[q][i]), (c_nt_gap[0][j] + ext[q][i]), &c_nt_gap[1][j], &track[q][i][j].w);				
				} else {
					//c_nt_gap[i][j] = -1;
					//track[3][i][j] = 1;
					c_nt_gap[1][j] = -1;
					track[q][i][j].w = 1;
				}
				//max_finder_and_track_threestates(a_nt_nt[i][j], b_gap_nt[i][j], c_nt_gap[i][j], &(best[i][j]), &(track[0][i][j]));
				max_finder_and_track_threestates(a_nt_nt[1][j], b_gap_nt[1][j], c_nt_gap[1][j], &(best[q][i][j]), &(track[q][i][j].x));
			}
			memcpy(a_nt_nt[0], a_nt_nt[1], sizeof(int) * (reference_length + 1));
			memcpy(b_gap_nt[0], b_gap_nt[1], sizeof(int) * (reference_length + 1));
			memcpy(c_nt_gap[0], c_nt_gap[1], sizeof(int) * (reference_length + 1));
		}	
		free2dInt(a_nt_nt, 2);
		free2dInt(b_gap_nt, 2); 
		free2dInt(c_nt_gap, 2);
	}	
} //*/


void updateScores(ushort4* scores, int good_count) {		
	//qsort(scores, good_count, sizeof(int4), cmpscores);
	nullify_overlaps(good_count, scores);
	//qsort(scores, good_count, sizeof(int4), cmpscores);//*/

}

void traceback_cpy(hit_struct* hit_ptr, S_HIT_struct s_hit_struct, char4 * s_hit_struct_alignment) {

	hit_ptr->query_end = s_hit_struct.query_end;
	hit_ptr->ref_end = s_hit_struct.ref_end;
	hit_ptr->score = s_hit_struct.score;
	hit_ptr->query_start = s_hit_struct.query_start;
	hit_ptr->ref_start = s_hit_struct.ref_start;

	for(int i = 0; i < MAX_QUERY_LENGHT+MAX_QUERY_APPENDIX; i++){		
		hit_ptr->alignment[0][i] = s_hit_struct_alignment[i].x;
		hit_ptr->alignment[1][i] = s_hit_struct_alignment[i].y;
		hit_ptr->alignment[2][i] = s_hit_struct_alignment[i].z;
		if(s_hit_struct_alignment[i].x == '\0') break;
	}
	
}

void traceback(int query_array_size_patched, int q, int* best, char4* track, char* sequence1, char* sequence2, int i, int j,
		unsigned short hit_score, int path, int* shift_offset_level2, int *shift_offset_level3, S_HIT_struct *s_hit_struct, char4 * s_hit_struct_alignment) {
	int length, track_array;
	length = 0;
	s_hit_struct->path = path;
	s_hit_struct->query_end = i;
	s_hit_struct->ref_end = j;
	s_hit_struct->score = hit_score;
	int shift_id = q + shift_offset_level3[i] +  shift_offset_level2[j + i];
	track_array = track[shift_id].x;
	while (best[shift_id] > 0) {
		if (track_array == 1) {
			length++;
			s_hit_struct_alignment[length - 1].x = sequence1[i - 1];			
			s_hit_struct_alignment[length - 1].y = ali_rep[(int)bases[(int)sequence1[i - 1]] * ali_patch + bases[(int)sequence2[j - 1]]];
			s_hit_struct_alignment[length - 1].z = sequence2[j - 1];
			track_array = track[shift_id].y;
			i--;
			j--;
			shift_id = q + shift_offset_level3[i] +  shift_offset_level2[j + i];
		} else if (track_array == 2) {
			length++;
			s_hit_struct_alignment[length - 1].x = '-';			
			s_hit_struct_alignment[length - 1].y = ' ';
			s_hit_struct_alignment[length - 1].z = sequence2[j - 1];
			track_array = track[shift_id].z;
			j--;
			shift_id = q + shift_offset_level3[i] +  shift_offset_level2[j + i];
		} else if (track_array == 3) {
			length++;
			s_hit_struct_alignment[length - 1].x = sequence1[i - 1];			
			s_hit_struct_alignment[length - 1].y = ' ';
			s_hit_struct_alignment[length - 1].z = '-';
			track_array = track[shift_id].w;
			i--;
			shift_id = q + shift_offset_level3[i] +  shift_offset_level2[j + i];
		} else {
			s_hit_struct_alignment[length - 1].x = '\0';			
			s_hit_struct_alignment[length - 1].y = '\0';
			s_hit_struct_alignment[length - 1].z = '\0';
			break;
		}
	}
	s_hit_struct->query_start = i;
	s_hit_struct->ref_start = j;
}

void old_traceback(int query_array_size_patched, int q, int* best, char4* track, char* sequence1, char* sequence2, int i, int j,
		hit_struct* hit_ptr, FloatType hit_score, int* shift_offset_level2, int *shift_offset_level3) {
	int length, track_array;
	length = 0;
	hit_ptr->query_end = i;
	hit_ptr->ref_end = j;
	hit_ptr->score = hit_score;
	int shift_id = q + shift_offset_level3[i] +  shift_offset_level2[j + i];
	track_array = track[shift_id].x;
	while (best[shift_id] > 0) {
		if (track_array == 1) {
			length++;
			hit_ptr->alignment[0][length - 1] = sequence1[i - 1];
			hit_ptr->alignment[2][length - 1] = sequence2[j - 1];
			hit_ptr->alignment[1][length - 1] = ali_rep[(int)bases[(int)sequence1[i - 1]] * ali_patch + bases[(int)sequence2[j - 1]]];
			track_array = track[shift_id].y;
			i--;
			j--;
			shift_id = q + shift_offset_level3[i] +  shift_offset_level2[j + i];
		} else if (track_array == 2) {
			length++;
			hit_ptr->alignment[0][length - 1] = '-';
			hit_ptr->alignment[2][length - 1] = sequence2[j - 1];
			hit_ptr->alignment[1][length - 1] = ' ';
			track_array = track[shift_id].z;
			j--;
			shift_id = q + shift_offset_level3[i] +  shift_offset_level2[j + i];
		} else if (track_array == 3) {
			length++;
			hit_ptr->alignment[0][length - 1] = sequence1[i - 1];
			hit_ptr->alignment[2][length - 1] = '-';
			hit_ptr->alignment[1][length - 1] = ' ';
			track_array = track[shift_id].w;
			i--;
			shift_id = q + shift_offset_level3[i] +  shift_offset_level2[j + i];
		} else {
			hit_ptr->alignment[0][length] = '\0';
			hit_ptr->alignment[1][length] = '\0';
			hit_ptr->alignment[2][length] = '\0';
			break;
		}
	}
	hit_ptr->query_start = i;
	hit_ptr->ref_start = j;
}
/*
void old_traceback(int query_array_size_patched, int q, int* best, char4* track, char* sequence1, char* sequence2, int i, int j,
		hit_struct* hit_ptr, FloatType hit_score, int* shift_offset_level2, int *shift_offset_level3) {
	int length, track_array;
	length = 0;
	hit_ptr->query_end = i;
	hit_ptr->ref_end = j;
	hit_ptr->score = hit_score;
	int shift_id = q + shift_offset_level3[i] +  shift_offset_level2[j + i];
	track_array = track[shift_id].x;
	while (best[shift_id] > 0) {
		if (track_array == 1) {
			length++;
			hit_ptr->alignment[0][length - 1] = sequence1[i - 1];
			hit_ptr->alignment[2][length - 1] = sequence2[j - 1];
			hit_ptr->alignment[1][length - 1] = ali_rep[(int)bases[(int)sequence1[i - 1]] * ali_patch + bases[(int)sequence2[j - 1]]];
			track_array = track[shift_id].y;
			i--;
			j--;
			shift_id = q + shift_offset_level3[i] +  shift_offset_level2[j + i];
		} else if (track_array == 2) {
			length++;
			hit_ptr->alignment[0][length - 1] = '-';
			hit_ptr->alignment[2][length - 1] = sequence2[j - 1];
			hit_ptr->alignment[1][length - 1] = ' ';
			track_array = track[shift_id].z;
			j--;
			shift_id = q + shift_offset_level3[i] +  shift_offset_level2[j + i];
		} else if (track_array == 3) {
			length++;
			hit_ptr->alignment[0][length - 1] = sequence1[i - 1];
			hit_ptr->alignment[2][length - 1] = '-';
			hit_ptr->alignment[1][length - 1] = ' ';
			track_array = track[shift_id].w;
			i--;
			shift_id = q + shift_offset_level3[i] +  shift_offset_level2[j + i];
		} else {
			hit_ptr->alignment[0][length] = '\0';
			hit_ptr->alignment[1][length] = '\0';
			hit_ptr->alignment[2][length] = '\0';
			break;
		}
	}
	hit_ptr->query_start = i;
	hit_ptr->ref_start = j;
}//*/

int overlap(int start_a, int end_a, int start_b, int end_b) {
	int max_start = start_a > start_b ? start_a : start_b;
	int min_end = end_a < end_b ? end_a : end_b;
	return min_end - max_start + 1;
}

int testfor_overlap(int* good_ones_starts_j, int* good_ones_ends_j, int* good_ones_count, int test_start, int test_end) {
	int index, good_call;
	good_call = 1;
	if (*good_ones_count < 0) {
		good_call = 1;
		*good_ones_count = *good_ones_count + 1;
		return good_call;
	}
	for (index = 0; index <= *good_ones_count; index++) {
		if (overlap(good_ones_starts_j[index],good_ones_ends_j[index],test_start,test_end) > 6) {
			good_call = 0;
			break;
		}
	}
	if (good_call == 1) {*good_ones_count = *good_ones_count + 1;}
	/*
	 printf("ref_s ref_e test_s test_e min_end max_start good_call good_ones_count %d %d %d %d %d %d %d %d\n",
			 ref_start, ref_end, test_start, test_end, min_end, max_start, good_call, *good_ones_count);
	 */
	return good_call;
}
