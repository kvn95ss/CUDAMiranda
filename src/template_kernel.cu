
__constant__ int d_c_ali_score[32];
__constant__ char d_c_ali_rep[32];
__constant__ int d_c_ali_score5p[32];
__constant__ int d_c_length_3p_for_weighting_array[MAX_NUM_QUERY_FOR_CONSTANT_MEMORY];
__constant__ int d_c_query_seqlen[MAX_NUM_QUERY_FOR_CONSTANT_MEMORY];
__constant__ char2 d_c_open_ext[MAX_NUM_QUERY_FOR_CONSTANT_MEMORY*MAX_QUERY_LENGHT];
__constant__ unsigned int d_c_in_seed[MAX_NUM_QUERY_FOR_CONSTANT_MEMORY];
__constant__ int d_c_reference_batch_seqlen[MAX_RefBatchSize];

texture<char2, cudaTextureType2D, cudaReadModeElementType> d_tex_reference_batch_sequence_raw; // Texture;
texture<char4, cudaTextureType2D, cudaReadModeElementType> d_tex_reference_batch_sequence; // Texture;
texture<char2, cudaTextureType2D, cudaReadModeElementType> d_tex_query_sequence_raw;
texture<char4, cudaTextureType2D, cudaReadModeElementType> d_tex_query_sequence;

#define Block_Size_X 32
#define ModBase32 31

#define RefBufferFactor 4
#define RefBufferLength 128
#define ModBaseRefBuffer 127

#define ModBase4 3

#define score(nt1, nt2)      d_c_ali_score[nt1 * ali_patch + nt2]
#define score5p(nt1, nt2)      d_c_ali_score5p[nt1 * ali_patch + nt2]

#define getX(val) (val&3)
#define getW(val) ((val>>2) & 1)
#define getZ(val) ((val>>3) & 1)
#define getI(val, pos) ((val>>pos) & 1)

__forceinline__ __device__ unsigned char getNextTrack(unsigned short val, int pos){
	return ((val >> ((3 - pos)*4)) & 15);
}//*/

__global__ void d_build_matrix(int d_g_shift_offset_batch_ref, int d_g_shift_offset_ref, int d_g_query_array_size_patched, 
	                           int d_g_shift_offset_socre_batch_ref, int d_g_score_threshold, 
							   int query_array_size,
	                           unsigned short* track, ushort4 * scores, int *good_count){

	__shared__ short  s_a_nt_nt1[MAX_QUERY_LENGHT][WARP_SIZE];  // MAX_QUERY_LENGHT by query_array->size_patched
	__shared__ short s_b_gap_nt1[MAX_QUERY_LENGHT][WARP_SIZE];
	__shared__ short s_c_nt_gap1[MAX_QUERY_LENGHT][WARP_SIZE];
	__shared__ short  s_a_nt_nt2[MAX_QUERY_LENGHT][WARP_SIZE];  // MAX_QUERY_LENGHT by query_array->size_patched
	__shared__ short s_b_gap_nt2[MAX_QUERY_LENGHT][WARP_SIZE];
	__shared__ short s_c_nt_gap2[MAX_QUERY_LENGHT][WARP_SIZE];  // 6 * 2k = 12288 used so far	

	                      	
	__shared__ int s_good_count[WARP_SIZE];                   
	__shared__ char s_ref_seq[WARP_SIZE];                      
	__shared__ char s_ref_seq_buffer[RefBufferLength];  
	__shared__ unsigned int s_in_seed[MAX_QUERY_LENGHT];         // 8*32 + 1 * 128 = 288  
	__shared__ char s_query_seq[MAX_QUERY_LENGHT][WARP_SIZE];  
	__shared__ char2 s_open_ext[MAX_QUERY_LENGHT][WARP_SIZE];  // 32*32*3 = 3072
	
	
	__shared__ unsigned short tmp_best2_track[MAX_QUERY_LENGTH_AFTER_SPLIT][WARP_SIZE];  // y z w x are represented by 2 bit in a char type with 8 bits.
	__shared__ unsigned char tmp_best2[MAX_QUERY_LENGHT][WARP_SIZE];       // 2048
	
	__shared__ unsigned char tmp_best1[MAX_QUERY_LENGHT][WARP_SIZE];       
	__shared__ unsigned char tmp_best0[MAX_QUERY_LENGHT][WARP_SIZE];       // 2048 		
	
	int query_array_size_patched = d_g_query_array_size_patched;
	int shift_offset_socre_batch_ref = d_g_shift_offset_socre_batch_ref;
	int d_score_threshold = d_g_score_threshold;

	int q = blockIdx.x * blockDim.x + threadIdx.x;	
	int blockIdx_z_shift_offset_batch_ref = blockIdx.z * d_g_shift_offset_batch_ref;
	// load 4 char in one time
	if(threadIdx.y < 8){
		char4 seq_tmp = tex2D(d_tex_query_sequence, q, threadIdx.y);
		int threadIdx_y_4 = threadIdx.y*4;
		s_query_seq[threadIdx_y_4][threadIdx.x] = seq_tmp.x;
		s_query_seq[threadIdx_y_4 + 1][threadIdx.x] = seq_tmp.y;
		s_query_seq[threadIdx_y_4 + 2][threadIdx.x] = seq_tmp.z;
		s_query_seq[threadIdx_y_4 + 3][threadIdx.x] = seq_tmp.w;
	}
	__syncthreads();
#pragma unroll	 
	for(int ii = 0; ii < KERNEL_SPLIT_FACTOR; ii++){
		int i = threadIdx.y* KERNEL_SPLIT_FACTOR + ii;		
		s_a_nt_nt1[i][threadIdx.x] = 0;
		s_a_nt_nt2[i][threadIdx.x] = 0;		
		s_b_gap_nt1[i][threadIdx.x] = 0;
		s_b_gap_nt2[i][threadIdx.x] = 0;		
		s_c_nt_gap1[i][threadIdx.x] = 0;
		s_c_nt_gap2[i][threadIdx.x] = 0;
		// load query onto shared memory;
		s_open_ext[i][threadIdx.x] = d_c_open_ext[q + i * query_array_size_patched];		
		tmp_best0[i][threadIdx.x] = 0;
		tmp_best1[i][threadIdx.x] = 0;
	}
	s_in_seed[threadIdx.x] = d_c_in_seed[q];
	 
	if(threadIdx.y == 0){		
		// load first 32 ref;	
		if(threadIdx.x < 8){
			char4 ref_tmp = tex2D(d_tex_reference_batch_sequence, threadIdx.x, blockIdx.z);
			s_ref_seq[threadIdx.x*4] = ref_tmp.x;
			s_ref_seq[threadIdx.x*4 + 1] = ref_tmp.y;
			s_ref_seq[threadIdx.x*4 + 2] = ref_tmp.z;
			s_ref_seq[threadIdx.x*4 + 3] = ref_tmp.w;
		}
		s_good_count[threadIdx.x] = 0;
	}	
	__syncthreads();	

	for (int j = 2; j <= (d_c_reference_batch_seqlen[blockIdx.z] + MAX_QUERY_LENGHT); j++) {	
			tmp_best2_track[threadIdx.y][threadIdx.x] = 0;
			for(int ii = 0; ii < KERNEL_SPLIT_FACTOR; ii++){
				int i = threadIdx.y * KERNEL_SPLIT_FACTOR + ii;
				tmp_best2[i][threadIdx.x] = 0;
				tmp_best2_track[threadIdx.y][threadIdx.x] <<= 4;
				unsigned short tmp_best2_track_var = 0;
				if( (i > 0) &&
					(i <= d_c_query_seqlen[q]) && 
					(j <= (d_c_reference_batch_seqlen[blockIdx.z] + d_c_query_seqlen[q])) && 
					(i < j) && 
					!((j > d_c_reference_batch_seqlen[blockIdx.z]) && (i < (j - d_c_reference_batch_seqlen[blockIdx.z])))
					){
						int tmp_score = 0;					
						if(i > 2){						
							if(i <= d_c_length_3p_for_weighting_array[q]){							
								tmp_score = score(s_query_seq[i - 1][threadIdx.x],s_ref_seq[(j - i - 1) & ModBase32]);								
							}else if(i < d_c_query_seqlen[q]){							
								tmp_score = score5p(s_query_seq[i - 1][threadIdx.x],s_ref_seq[(j - i - 1) & ModBase32]);								
							}							
						}						
						s_a_nt_nt2[i][threadIdx.x] = max(0, tmp_best0[i-1][threadIdx.x] + tmp_score);
						short a = s_a_nt_nt1[i][threadIdx.x] + s_open_ext[i][threadIdx.x].x;
						s_b_gap_nt2[i][threadIdx.x] = max(a, (s_b_gap_nt1[i][threadIdx.x] + s_open_ext[i][threadIdx.x].y));
						tmp_best2_track_var |= (s_b_gap_nt2[i][threadIdx.x] == a); // this is for z;						
						tmp_best2_track_var <<= 1;
						if (getI(s_in_seed[threadIdx.x], i)) {						
							s_c_nt_gap2[i][threadIdx.x] = -1;						
							tmp_best2_track_var |= 1; // for  w;
						}else{
							a = s_a_nt_nt1[i - 1][threadIdx.x] + s_open_ext[i][threadIdx.x].x;
							s_c_nt_gap2[i][threadIdx.x] = max(a, (s_c_nt_gap1[i - 1][threadIdx.x] + s_open_ext[i][threadIdx.x].y));
							tmp_best2_track_var |= (s_c_nt_gap2[i][threadIdx.x] == a);	// for w															
						}
						tmp_best2_track_var <<= 2;
						tmp_best2[i][threadIdx.x] = max(max(max(s_a_nt_nt2[i][threadIdx.x], s_b_gap_nt2[i][threadIdx.x]), s_c_nt_gap2[i][threadIdx.x]), 0);	
						tmp_best2_track_var |= (((3^(tmp_best2[i][threadIdx.x] == s_b_gap_nt2[i][threadIdx.x]))>>(tmp_best2[i][threadIdx.x] == s_a_nt_nt2[i][threadIdx.x])) >> (tmp_best2[i][threadIdx.x] == 0));
				} // for if( (i > 0) &&		
				tmp_best2_track[threadIdx.y][threadIdx.x] |= tmp_best2_track_var;
			} // for (ii)
		//} // if(q < query_array_size)
		// update array
		__syncthreads();
		if(threadIdx.y == 0){
#pragma unroll 1 // do not unroll this;
			for(int i = 0; i < MAX_QUERY_LENGHT; i++){
				if (tmp_best2[i][threadIdx.x] >= d_score_threshold && getX(getNextTrack(tmp_best2_track[i/4][threadIdx.x], (i&3))) != 2 && s_good_count[threadIdx.x] < MAX_REPORTS_OF_FIRST_HIT-1) {				
					// reject alignments ending with an unaligned reference base										
					ushort4 tmp_scores;						
					tmp_scores.x = tmp_best2[i][threadIdx.x];
					tmp_scores.y = s_good_count[threadIdx.x];
					tmp_scores.z = i;
					tmp_scores.w = j - i;
					scores[blockIdx.z * shift_offset_socre_batch_ref + s_good_count[threadIdx.x]*query_array_size_patched + q] = tmp_scores;
					//printf("%d, %d,%d,%d,%d,%d\n", q, tmp_best[i][threadIdx.x], s_good_count[threadIdx.x],threadIdx.x, i, j-i);	
					atomicAdd(&s_good_count[threadIdx.x], 1);	//*/
				}
			} //*/
			int jj = j - 2;// 2 -> 0; 33 -> 31; 34 - > 32;								
			if((jj & ModBaseRefBuffer) == 0){				
				int threadIdx_x4 = threadIdx.x*4;					
				char4 ref_tmp = tex2D(d_tex_reference_batch_sequence, (jj/4) + threadIdx.x, blockIdx.z);
				s_ref_seq_buffer[threadIdx_x4] = ref_tmp.x; 
				s_ref_seq_buffer[threadIdx_x4 + 1] = ref_tmp.y;
				s_ref_seq_buffer[threadIdx_x4 + 2] = ref_tmp.z;
				s_ref_seq_buffer[threadIdx_x4 + 3] = ref_tmp.w;
			}
			if(threadIdx.x == 0){					
				s_ref_seq[jj & ModBase32] = s_ref_seq_buffer[jj & ModBaseRefBuffer];					
			}
		}
		__syncthreads();	
		int threadIdx_y4 = threadIdx.y * KERNEL_SPLIT_FACTOR;
		int shift_id3 = blockIdx_z_shift_offset_batch_ref + j * d_g_shift_offset_ref + threadIdx.y*query_array_size_patched + q;
		track[shift_id3] = tmp_best2_track[threadIdx.y][threadIdx.x];
#pragma unroll	
		for(int ii = 0; ii < KERNEL_SPLIT_FACTOR; ii++){
			int i = threadIdx_y4 + ii;
			tmp_best0[i][threadIdx.x] =  tmp_best1[i][threadIdx.x];
			tmp_best1[i][threadIdx.x] = tmp_best2[i][threadIdx.x];

			s_a_nt_nt1[i][threadIdx.x] =  s_a_nt_nt2[i][threadIdx.x];			
			s_b_gap_nt1[i][threadIdx.x] = s_b_gap_nt2[i][threadIdx.x];			
			s_c_nt_gap1[i][threadIdx.x] = s_c_nt_gap2[i][threadIdx.x];	
		}		
		__syncthreads();
	} // for (int j = 1; j <= reference_length; j++) {	
	if(threadIdx.y == 0){
		good_count[blockIdx.z * query_array_size_patched + q] = s_good_count[threadIdx.x];	
	}
}

__global__ void d_updateScores(int shift_offset_score_level1, int query_array_size_patched, ushort4* d_scores, int* d_good_count){		
		__shared__ int s_good_count[WARP_SIZE];    
		// the max reports is bounded by the limits of shared memeory
		// (48 * 1024 - 4 * 32) /(2*4)/32 = 191
		__shared__ ushort4 s_scores[MAX_REPORTS_OF_FIRST_HIT][WARP_SIZE];		
		
		int j_shift_offset_score_level1 = blockIdx.z * shift_offset_score_level1;		
		int q = blockIdx.x * blockDim.x + threadIdx.x;	
		s_good_count[threadIdx.x] = d_good_count[blockIdx.z  * query_array_size_patched + q];
		// load data onto shared memory;
		for (int index = 0; index < s_good_count[threadIdx.x]; index++) {
			s_scores[index][threadIdx.x] = d_scores[j_shift_offset_score_level1 + index * query_array_size_patched + q];			
		}
		__syncthreads();
		// bubble sorting		
		for(int i = 1; i < s_good_count[threadIdx.x]; i++)
		{
			for (int j = 0; j < i; j++)
			{
				if (s_scores[i][threadIdx.x].x > s_scores[j][threadIdx.x].x || (s_scores[i][threadIdx.x].x == s_scores[j][threadIdx.x].x && s_scores[i][threadIdx.x].w <= s_scores[j][threadIdx.x].w))      // ascending order simply changes to <
				{ 
					ushort4 tmp = s_scores[j][threadIdx.x];	
					s_scores[j][threadIdx.x] = s_scores[i][threadIdx.x];
					s_scores[i][threadIdx.x] = tmp;
				}
			}
		}
		__syncthreads();		
		for (int index = 0; index < s_good_count[threadIdx.x]; index++) {				
			if (s_scores[index][threadIdx.x].x == 0) {continue;}		
			/* use query end to figure where ref offset relative to fictional query offset 0 would be */
			int outer_ref_end = s_scores[index][threadIdx.x].w - s_scores[index][threadIdx.x].z;
			for (int inner_index = (s_good_count[threadIdx.x]-1); inner_index > index; inner_index--) {					
				int inner_ref_end = s_scores[inner_index][threadIdx.x].w - s_scores[inner_index][threadIdx.x].z;										
				if (abs(inner_ref_end - outer_ref_end) <= OVERLAP_LIMIT) {
					s_scores[inner_index][threadIdx.x].x = 0;					
				}
			}
		} // for (index = 0; index < total_elements; index++) 			
		__syncthreads();
		for (int index = 0; index < s_good_count[threadIdx.x]; index++) {
			d_scores[j_shift_offset_score_level1 + index * query_array_size_patched + q] = s_scores[index][threadIdx.x];
		}
}

__global__ void d_traceback(int score_threshold, int d_shift_offset_level1,
							int d_shift_offset_level1_hit_alignment, int d_shift_offset_level2_hit_alignment, 
							int d_shift_offset_level2,
	                        int d_shift_offset_level1_score, int query_array_size_patched, 
							int * good_count, unsigned short* track, ushort4* scores,
							S_HIT_struct *s_hit_struct, char4 * s_hit_struct_alignment){
		__shared__ int s_good_count[WARP_SIZE];		
		__shared__ char s_ali_rep[32];  
		s_ali_rep[threadIdx.x] = d_c_ali_rep[threadIdx.x];
		__syncthreads();	
		int jj = blockIdx.z;
		int shift_id_count = jj * query_array_size_patched;		
		//for(int q = 0; q < query_array->size; q++){			    			
		int q = blockIdx.x * blockDim.x + threadIdx.x;		
		int d_shift_offset_level1_score_jj = d_shift_offset_level1_score * jj;
		int d_shift_offset_level1_hit_alignment_jj = d_shift_offset_level1_hit_alignment *jj;
		int d_shift_offset_level3_hit_alignment_q = (MAX_QUERY_LENGHT +MAX_QUERY_APPENDIX)*q;
		int d_shift_offset_level1_jj = d_shift_offset_level1 * jj;
		s_good_count[threadIdx.x] = good_count[shift_id_count + q];
		for (int ii = 0; ii < s_good_count[threadIdx.x]; ii++) {						
			int tmp_shift_id = d_shift_offset_level1_score_jj + query_array_size_patched *ii + q;		
			int tmp_shift_id_hit_aligment = d_shift_offset_level1_hit_alignment_jj + d_shift_offset_level2_hit_alignment *ii + d_shift_offset_level3_hit_alignment_q ;
			if (scores[tmp_shift_id].x >= score_threshold) {					
				int length = 0;	
				int i = scores[tmp_shift_id].z;
				int j = scores[tmp_shift_id].w;
				s_hit_struct[tmp_shift_id].path = scores[tmp_shift_id].y;
				s_hit_struct[tmp_shift_id].query_end = i;
				s_hit_struct[tmp_shift_id].ref_end = j;
				s_hit_struct[tmp_shift_id].score = scores[tmp_shift_id].x;
				int shift_id = d_shift_offset_level1_jj + q + query_array_size_patched * (i/4) +  d_shift_offset_level2 * (j + i);				
				unsigned char tmp_track = getNextTrack(track[shift_id], (i&3));
				unsigned char track_array = getX(tmp_track);
				while (getX(tmp_track) > 0){ //(best[shift_id] > 0) {
					if(length >= (MAX_QUERY_LENGHT+MAX_QUERY_APPENDIX)){
						break;
					}
					if (track_array == 1) {
						length++;
						char2 sequence1 =  tex2D(d_tex_query_sequence_raw, q, (i - 1)); // query_array->sequence1D_reshape_raw[(i - 1)*query_array->size_patched + q];
						char2 sequence2 = tex2D(d_tex_reference_batch_sequence_raw, j - 1, jj); //reference_array->sequence1D_revalued_raw[(jj*reference_array->seqlen_max_patched) + j - 1];							
						s_hit_struct_alignment[tmp_shift_id_hit_aligment + length - 1].x = sequence1.x;			
						s_hit_struct_alignment[tmp_shift_id_hit_aligment + length - 1].y = s_ali_rep[sequence1.y * ali_patch + sequence2.y]; // sequence2_base
						s_hit_struct_alignment[tmp_shift_id_hit_aligment + length - 1].z = sequence2.x;
						
						i--;
						j--;
						shift_id = d_shift_offset_level1_jj + q + query_array_size_patched * (i/4) +  d_shift_offset_level2 * (j + i);
						tmp_track = getNextTrack(track[shift_id], (i&3));
						track_array = getX(tmp_track);						
					} else if (track_array == 2) {
						length++;
						char2 sequence2 = tex2D(d_tex_reference_batch_sequence_raw, j - 1, jj); //reference_array->sequence1D_revalued_raw[(jj*reference_array->seqlen_max_patched) + j - 1];
						s_hit_struct_alignment[tmp_shift_id_hit_aligment + length - 1].x = '-';			
						s_hit_struct_alignment[tmp_shift_id_hit_aligment + length - 1].y = ' ';
						s_hit_struct_alignment[tmp_shift_id_hit_aligment + length - 1].z = sequence2.x;
						track_array = (2 >> (getZ(tmp_track)));
						j--;
						shift_id = d_shift_offset_level1_jj + q + query_array_size_patched * (i/4) +  d_shift_offset_level2 * (j + i);
						tmp_track = getNextTrack(track[shift_id], (i&3));						
					} else if (track_array == 3) {
						length++;
						char2 sequence1 = tex2D(d_tex_query_sequence_raw, q, (i - 1)); //query_array->sequence1D_reshape_raw[(i - 1)*query_array->size_patched + q];
						s_hit_struct_alignment[tmp_shift_id_hit_aligment + length - 1].x = sequence1.x;			
						s_hit_struct_alignment[tmp_shift_id_hit_aligment + length - 1].y = ' ';
						s_hit_struct_alignment[tmp_shift_id_hit_aligment + length - 1].z = '-';
						track_array = (3 >> (getW(tmp_track)));
						i--;
						shift_id = d_shift_offset_level1_jj + q + query_array_size_patched * (i/4) +  d_shift_offset_level2 * (j + i);
						tmp_track = getNextTrack(track[shift_id], (i&3));						
					} else {
						s_hit_struct_alignment[tmp_shift_id_hit_aligment + length].x = '\0';	// -1		
						s_hit_struct_alignment[tmp_shift_id_hit_aligment + length].y = '\0';    // -1
						s_hit_struct_alignment[tmp_shift_id_hit_aligment + length].z = '\0';    // -1
						break;
					}
										
				}
				s_hit_struct[tmp_shift_id].query_start = i;
				s_hit_struct[tmp_shift_id].ref_start = j;
			}
		}
}