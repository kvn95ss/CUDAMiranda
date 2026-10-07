/* Adjustable Algorithm Parameters*/
#define INTBITS32 32
// (3+1)*32N + 8N + 3*32*4 < 64000  N = 448 = 14 * 32;
#define MAX_NUM_QUERY_FOR_CONSTANT_MEMORY 480
// to achieve the maximum performance we set Max_Query_length = 32;
#define MAX_QUERY_LENGHT 32 
#define MAX_QUERY_LENGHT_HALF 16
#define KERNEL_SPLIT_FACTOR 4
#define MAX_QUERY_LENGTH_AFTER_SPLIT 8
#define MAX_RefBatchSize 2048
// make sure there are enough room for alignment score;  (max_query_length + max_query_appendix)*char4 , it should be char3 however, char4 save more space.
#define MAX_QUERY_APPENDIX 8
// warp size;
#define WARP_SIZE 32
// the max reports is bounded by the limits of shared memeory
// (48 * 1024 - 4 * 32) /(2*4+1)/32 = 191
#define MAX_REPORTS_OF_FIRST_HIT 191 
#define OVERLAP_LIMIT 6

typedef float FloatType;
typedef char4 char4;
extern FloatType scale;
extern int strict;
extern int debug;
extern FloatType gap_open;
extern FloatType gap_extend;
extern FloatType score_threshold;
extern FloatType energy_threshold;
extern int length_5p_for_weighting;
extern int* length_3p_for_weighting_array;
extern int key_value_pairs;
extern int no_energy;
extern int verbosity;
extern int outfile;
extern int truncated;
extern int restricted;
extern int ref_batch_size;
extern int nGPUs;
//extern const int MAX_QUERY_LENGHT;
//extern const int WARP_SIZE;
extern int bases[256];
#define ali_patch 5
extern int ali_score[32]; //[5][5];
extern int ali_score5p[32]; //[5][5];
extern char ali_rep[32]; //[5][5];