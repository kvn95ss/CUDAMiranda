struct ExpStringT;
typedef struct ExpStringT ExpString;
/* summary of best hits for reporting - for convenient argument passing */
typedef struct HitSummaryT
{
	FloatType scan_score;
	int no_hits;
	FloatType max_hit;
	FloatType max_score;
	FloatType total_score;
	ExpString *position_list;
} HitSummary;

typedef struct intXY{
	int id;
	int val;
}intXY;

typedef struct input_type
{	
	char** sequence;  // the pointer here indicate the 2d location in sequence1D;
	char* sequence1D; // all the sequences will be flatten into 1D for CUDA;
	char4* sequence1D_reshape; // reshape the 1D sequence so that it can fit the CUDA memory alignment. replaced with base[] array
	char2* sequence1D_reshape_raw; // reshape the 1D sequence so that it can fit the CUDA memory alignment. without replaced base[] array  .x is raw, .y is replaced data
	char** description;
	char** name;
	int* seqlen;	
	int    size;
	int    size_patched;
	int	   seqlen_max;
	int	   seqlen_min;
	FloatType seqlen_mean;
	FloatType seqlen_sd; // standard deviation;
	int total_seqlen; // sum over the length of all query;
}input_type;

typedef struct Reference_batch{
	char** sequence;
	char* sequence1D; // all the sequences will be flatten into 1D for CUDA;
	char4* sequence1D_revalued; // all the sequences will be flatten into 1D for CUDA; replaced with base[] array
	char2* sequence1D_revalued_raw; // all the sequences will be flatten into 1D for CUDA; without replaced base[] array
	char** description;
	char** name;
	int*   seqlen;
	int    size;
	int	   seqlen_max;
	int    seqlen_max_patched; // make sure seqlen_max is a multiply of 32;
	int    seqlen_max_patched_char4; // make sure seqlen_max is a multiply of 32;	
	HitSummary** hit_summary;		/*Summary of best hits for reporting*/
	FloatType** end_score;
	int ** restricted_Array;
}Reference_batch;

typedef struct Ref_batch{
	Reference_batch **reference_batch;
	int    batch_length;	
	int    ref_batch_size;
	//int    ref_batch_size_patched;

	int    overall_reference_size;
	int	   overall_seqlen_max;
	int    overall_seqlen_max_patched;
	int	   overall_seqlen_min;
	FloatType overall_mean;
	FloatType overall_sd; // standard deviation; 
}Ref_batch;

/* Structure to store individual hit information*/
typedef struct hit_struct
{
	FloatType score;
	int query_start;
	int query_end;
	int ref_start;
	int ref_end;
	char* alignment[3];
	char* rest[6];
}hit_struct;

/* Structure to store individual hit information*/
typedef struct S_HIT_struct
{
	unsigned short score;	
	unsigned short ref_start;
	unsigned short ref_end;
	unsigned short path; // y;
	unsigned char query_start;
	unsigned char query_end;
}s_HIT_struct;

typedef struct TempData{
	short d_a_nt_nt1;	
	short d_c_nt_gap1;	
	unsigned char d_tmp_best0;
	unsigned char d_tmp_best0_1_track_x;
}TempData;

typedef struct MultiCPUPlan{
	S_HIT_struct *s_hit_struct;
	char4 * s_hit_struct_alignment;
}MultiCPUPlan;

typedef struct MultiGPUPlan{
	// initialize CUDA Array;
	cudaArray * d_reference_batch_sequence; // Texture;
	cudaArray * d_reference_batch_sequence_raw; // Texture;
	cudaArray * d_query_sequence; // Texture;	
	cudaArray * d_query_sequence_raw; // Texture;	 

	unsigned short * d_track; // global
	//short * d_best;    // global
	ushort4* d_scores;	 // global	
	int* d_good_count;// global
	//int* d_reference_batch_seqlen; // global
	S_HIT_struct *d_s_hit_struct;
	char4 * d_s_hit_struct_alignment;
	TempData *d_tempData;
	//Stream for asynchronous command execution
    cudaStream_t stream;
}MultiGPUPlan;

/* Score structure for non-optimal path detection*/
//typedef struct score_struct
//{
//	int score; //x
//	int path;  //y
//	int query_trace_end; //z
//	int reference_trace_end; //w
//} score_struct; //

/* Structure for pair-wise restriction*/
typedef struct pair_struct
{
	char identifier1[200];
	char identifier2[200];
} pair_struct;



