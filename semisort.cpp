#include "semisort.h"

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <map>

#include "parlaylib/include/parlay/parallel.h"
#include "parlaylib/include/parlay/internal/get_time.h"

using Type = long long;

inline uint64_t hash64(uint64_t u) {
  uint64_t v = u * 3935559000370003845ul + 2691343689449507681ul;
  v ^= v >> 21;
  v ^= v << 37;
  v ^= v >> 4;
  v *= 4768777513237032717ul;
  v ^= v << 20;
  v ^= v >> 41;
  v ^= v << 5;
  return v;
}

int main(int argc, char* argv[]) {
  size_t n = 1e9;
  int num_rounds = 3;
  int random_seed = 1;
  if (argc >= 2) n = atoll(argv[1]);
  if (argc >= 3) num_rounds = atoi(argv[2]);
  if (argc >= 4) random_seed = atoi(argv[3]);  

  Type* A = (Type*)malloc(n * sizeof(Type));
  parlay::parallel_for(0, n, [&](size_t j) { A[j] = hash64(j * random_seed) % 1000; });

  // [Checker] get counts
  std::map<int, int> counts;
  for (size_t i = 0; i < n; ++i) {
    ++counts[A[i]];
  }

  double total_time = 0;
  for (int i = 0; i <= num_rounds; i++) {
    // Generate random arrays

    parlay::internal::timer t;
    semisort(A, n);
    t.stop();

    // [Checker] check counts and semisorted (when it changes, there should be no more left)
    int curr = -1;
    auto remaining = counts;
    for (size_t i = 0; i < n; ++i) {
      const auto x = A[i];
      if (x != curr) {
        if (remaining[curr]) {
          std::cout << "Output is not semisorted\n";
          std::cout << "Segment with key=" << curr << " ended with " << remaining[x] << " missing elements\n";
          exit(0);
        }
        curr = x;
      }
      --remaining[x];
    }


    if (i == 0) {
      std::cout << "Warmup round running time: " << t.total_time() << std::endl;
    } else {
      std::cout << "Round " << i << " running time: " << t.total_time()
                << std::endl;
      total_time += t.total_time();
    }
  }
  std::cout << "Average running time: " << total_time / num_rounds << std::endl;
  return 0;
}