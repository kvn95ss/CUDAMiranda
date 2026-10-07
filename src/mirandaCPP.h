/* Expandable string type (allocates memory as needed) */
void create_ExpString(ExpString **es);
void destroy_ExpString(ExpString **es);
void clear_ExpString(ExpString *es);
int length_ExpString(ExpString *es);
void append_char_ExpString(ExpString *es, char c);
void append_string_ExpString(ExpString *es, char *c);
void append_int_ExpString(ExpString *es, int i);
char *access_ExpString(ExpString *es); 

/* Functions Declarations */
/*   in scan.c */
int find_targets(FILE*, FILE*, pair_struct*, int, char*);
/*   in swat.c */
void traceback(int query_array_size_patched, int q, int*, char4*, char*, char*, int, int, unsigned short, int path, int* shift_offset_level2, int *shift_offset_level3,S_HIT_struct *s_hit_struct, char4 * s_hit_struct_alignment);
void traceback_cpy(hit_struct* hit_ptr, S_HIT_struct s_hit_struct, char4 * s_hit_struct_alignment);
void old_traceback(int query_array_size_patched, int q, int*, char4*, char*, char*, int, int, hit_struct*, FloatType, int* shift_offset_level2, int *shift_offset_level3);
int testfor_overlap(int*, int*, int*, int, int);
void build_matrix(int shift_offset_ref, input_type * query_array, char* ref_seq, int* best, char4* track,
		          int reference_length, int* length_3p_for_weighting_array, char3 *open_ext_in_seed, ushort4* scores, int *good_count);
void updateScores(ushort4* scores, int good_count);
void initialize_bases();
void initialize_scores();
/*   in output.c */
int parse_command_line(int, char**, char*, char*, char*, char*);
void printhit(char*, int, char*, hit_struct*, double, int, FILE*);
void print_parameters(char*, char*, FILE*);


/*   in seqio.c */
input_type* loadRef(char* filename);
input_type* loadQuery(FILE* query_fp);
input_type *readinallseq(FILE* fp);
int readinseq(FILE*, char**, char**, char**, int *);
void initialize_seqio_buffers();
void destroy_seqio_buffers(); //*/
/*   in utils.c */
void revstring(char s[]);
void clear_hit(hit_struct*, int, int);
int cmpscores(const void*, const void*);
void string_toupper(char*);
void freeInputType(input_type * in);
int ** new2dInt(int size1, int size2);
void free2dInt(int ** in, int len);
void free2dFloatType(FloatType ** in, int len);
FloatType ** new2dFloatType(int size1, int size2);
HitSummary ** new2dHitSummary(int size1, int size2);
void free2dHitSummary(HitSummary ** in, int size1, int size2);
int *** new3dInt(int size1, int size2, int size3);
void free3dInt(int *** in, int size1, int size2);
hit_struct * newHit_struct(int query_length, int reference_length);
void freeHit_struct(hit_struct * hit);
int **** new4dInt(int size1, int size2, int size3, int size4);
void free4dInt(int **** in, int size1, int size2, int size3);
int4 ** new2dScore_struct(int size1, int* size2, int factor);
void free2dScore_struct(int4 ** in, int size1);
int **** new4dTrack(int size1, int size2, int* size3, int size4);
void free4dTrack(int **** in, int size1, int size2, int* size3);
int *** new3dTable(int size1, int* size2, int size3);
int *** new3dTable_s(int size1, int size2, int size3);
void free3dTable(int *** in, int size1, int* size2);
void free3dTable_s(int *** in, int size1, int size2);
int * new1dInt(int size1);
void clear1dInt(int* in, int size1);
void clear4dTrack(int **** out, int size1, int size2, int* size3, int size4);
void clear3dTable(int *** out, int size1, int* size2, int size3);
void clear3dTable_s(int *** out, int size1, int size2, int size3);
void clear2dScore_struct(int4** scores, int size1, int* size2, int factor);
void free4dTable(int **** in, int size1, int size2, int * size3);
int **** new4dTable(int size1, int size2, int* size3, int size4);
void clear4dTable(int **** out, int size1, int size2, int* size3, int size4);
void free5dTrack(char4 **** in, int size1, int size2, int* size4);
char4 **** new5dTrack(int size1, int size2, int* size4, int size5);
void clear5dTrack(char4 **** out, int size1, int size2, int* size4, int size5);
int4 *** new3dScore_struct(int size1, int size2, int* size3, int factor);
void clear3dScore_struct(int4*** scores, int size1, int size2, int* size3, int factor);
void free3dScore_struct(int4 *** in, int size1, int size2);
void clear2dInt(int **in, int size1, int size2);
Ref_batch* getRef_batch(char * filename, input_type * query_array, pair_struct* pairs, int total_pairs);
void freeReference_batch(Reference_batch * in, int query_array_size);
void freeRef_batch(Ref_batch* ref_batch, int query_array_size);
void recalloc3dScore_struct(int4*** scores, int size1, int size2, int* size3, int factor);
void recalloc5dTrack(char4 ***** out, int size1, int size2, int size3, int* size4, int size5);
void recalloc4dTable(int **** out, int size1, int size2, int* size3, int size4);
/*   in pairs.c */
int load_pairs(FILE*, pair_struct**);
int find_pair(char*, char*, int, pair_struct*);
/*   in thermo.c */
FloatType get_energy(hit_struct*);
int cmpsIntXY(const void* p1, const void* p2);
int * memcpyInt(int * source, int size);
char * memcpyChar(char * source, int size);
void memcpyCharWithFixedMaxLength(char* dist, char * source, int size, int patch_size);