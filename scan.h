#include <math.h>

#include "parlaylib/include/parlay/parallel.h"
using namespace parlay;

/**
 * Inplace block scan DIRECTLY without filter (from coding assignment 2)
 */
template <typename T>
T scan(T *A, size_t n) {
  /**
   * REPLACE THE BODY OF THE FOLLOWING
   * FUNCTION WITH YOUR PARALLEL IMPLEMENTATION
  */
  const size_t k  = std::sqrt(n);  // chunks to break into
  const size_t sn = (n+k-1)/k;     // chunk size
  T* S = (T*)malloc( k * sizeof(T) );

  // compute S[0:k]
  parallel_for(0, k, [&](size_t i) {
    S[i] = 0;
    const size_t j_start = i * sn;
    const size_t j_end = (i+1)*sn < n ? (i+1)*sn : n;
    for (size_t j = j_start; j < j_end; ++j) {
      S[i] += A[j];
    }
  });

  // calculate prefix sum of S on S
  T s_total = 0;
  for (size_t i = 0; i < k; ++i) {
    T tmp = S[i];
    S[i] = s_total;
    s_total += tmp;
  }

  // calculate B using A and S
  T a_total;
  parallel_for(0, k, [&](size_t i) {
    T total = S[i];
    const size_t j_start = i * sn;
    const size_t j_end   = (i+1)*sn < n ? (i+1)*sn : n;
    for (size_t j = j_start; j < j_end; ++j) {
      T tmp = A[j];
      A[j] = total;
      total += tmp;
    }
    if (i == k-1) a_total = total;
  });
  return a_total;
}
