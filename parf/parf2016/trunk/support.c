#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <signal.h>
#include <time.h>
#include <sys/time.h>

#define BENCHMARK_BYTES (8*1024*1024)
#define BENCHMARK_DURATION 2
#define BENCHMARK_TYPE double
#define BENCHMARK_SIZE (BENCHMARK_BYTES / sizeof(BENCHMARK_TYPE))

#define INT8 long long
#define INT4 long
#define REAL4 float

INT8 microtime_() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return ((INT8) tv.tv_sec) * 1000000 + tv.tv_usec;
}

INT4 timeseed_() {
  return (INT4) microtime_();
}

INT4 ismissing_(float *number) {
  return isnan(*number);
}

int benchmark_waiting;
void benchmark_done(int dummy) {
  benchmark_waiting = 0;
}

unsigned long int benchmark_loc() {
  static unsigned long int loc;
  loc += random();
  return (loc % BENCHMARK_SIZE);
}

REAL4 benchmark_() {
  void (*oldalarm)(int);
  unsigned long int iter;
  BENCHMARK_TYPE *matrix;
  BENCHMARK_TYPE temp, sum = 0;

  benchmark_waiting = -1;
  oldalarm = signal(SIGALRM, benchmark_done);
  alarm(BENCHMARK_DURATION);
  
  matrix = (BENCHMARK_TYPE *) malloc(BENCHMARK_SIZE * sizeof(BENCHMARK_TYPE));
  for (iter = 0; benchmark_waiting; iter++) {
    temp = matrix[benchmark_loc()];
    if (!temp) temp = 1;
    sum += matrix[benchmark_loc()] / temp;
  }
  free(matrix);
  signal(SIGALRM, oldalarm);
  return ((REAL4) 1.0)/iter;
}

INT8 (*microtime__)() = microtime_;
INT4 (*timeseed__)() = timeseed_;
INT4 (*ismissing__)(REAL4 *) = ismissing_;
REAL4 (*benchmark__)() = benchmark_;
