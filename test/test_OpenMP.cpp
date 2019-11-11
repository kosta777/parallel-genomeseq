#ifdef USEOMP
#include <gtest/gtest.h>
#include <cstdio>
#include <omp.h>

TEST(OpenMP, Usability) {
  #pragma omp parallel for default(none) num_threads(3)
  for(int n=0; n<6; ++n) printf("Thread: %d, Num: %d\n", omp_get_thread_num(), n);
}
#endif
