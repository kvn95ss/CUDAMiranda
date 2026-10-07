/* Minimal replacement for the long-gone CUDA SDK cutil helpers used by cudaMiranda. */
#ifndef CUTIL_INLINE_COMPAT_H
#define CUTIL_INLINE_COMPAT_H
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <cuda_runtime.h>

#define cutilSafeCall(call) do { cudaError_t e_ = (call); if (e_ != cudaSuccess) { \
    fprintf(stderr, "CUDA error %s:%d: %s\n", __FILE__, __LINE__, cudaGetErrorString(e_)); exit(EXIT_FAILURE); } } while (0)
#define cutilCheckMsg(msg) do { cudaError_t e_ = cudaGetLastError(); if (e_ != cudaSuccess) { \
    fprintf(stderr, "%s:%d: %s: %s\n", __FILE__, __LINE__, (msg), cudaGetErrorString(e_)); exit(EXIT_FAILURE); } } while (0)
#define cutilCheckError(x) (void)(x)
#define cutilDeviceReset() cudaDeviceReset()

/* Simple wall-clock stopwatches addressed by integer handle (values in ms). */
#define CUT_MAX_TIMERS 64
static double cut_start_[CUT_MAX_TIMERS], cut_total_[CUT_MAX_TIMERS];
static int cut_used_[CUT_MAX_TIMERS];
static inline double cut_now_ms_(void) { struct timeval tv; gettimeofday(&tv, 0); return tv.tv_sec * 1e3 + tv.tv_usec / 1e3; }
static inline int cutCreateTimer(unsigned int *h) {
    for (unsigned int i = 1; i < CUT_MAX_TIMERS; i++) if (!cut_used_[i]) { cut_used_[i] = 1; cut_total_[i] = 0; *h = i; return 1; }
    return 0;
}
static inline int cutDeleteTimer(unsigned int h) { cut_used_[h] = 0; return 1; }
static inline int cutStartTimer(unsigned int h) { cut_start_[h] = cut_now_ms_(); return 1; }
static inline int cutStopTimer(unsigned int h) { cut_total_[h] += cut_now_ms_() - cut_start_[h]; return 1; }
static inline float cutGetTimerValue(unsigned int h) { return (float)cut_total_[h]; }
#endif
