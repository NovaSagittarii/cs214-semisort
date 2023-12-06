#include "semisort.h"

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <map>
#include <utility>

#include "parlaylib/include/parlay/parallel.h"
#include "parlaylib/include/parlay/internal/get_time.h"

#include "hash64.h"

using Type = std::pair<long long, long long>;

void run_test(size_t n, size_t keys, int num_rounds, int random_seed) {
  Type* A = (Type*)malloc(n * sizeof(Type));
  parlay::parallel_for(0, n, [&](size_t j) { A[j] = {hash64(j * random_seed) % keys, hash64(j * random_seed)}; });
  
  // Note: checker doesn't run well with few duplicates.
  const bool run_checker = false && (n <= 1e7 || keys <= 1e4);
  // [Checker] get counts
  std::map<long long, int> counts;
  if (run_checker) {
    std::cout << "checker is enabled" << std::endl;
    for (size_t i = 0; i < n; ++i) {
      ++counts[A[i].first];
    }
  }

  double total_time = 0;
  for (int i = 0; i <= num_rounds; i++) {
    // Generate random arrays
    parlay::parallel_for(0, n, [&](size_t j) { A[j] = {hash64(j * random_seed) % keys, hash64(j * random_seed)}; });

    parlay::internal::timer t;
    semisort(A, n);
    t.stop();

    // [Checker] check counts and semisorted (when it changes, there should be no more left)
    if (run_checker) {
      long long curr = -1;
      int remaining = 0;
      for (size_t i = 0; i < n; ++i) {
        const auto x = A[i].first;
        if (x != curr) {
          if (remaining) {
            std::cout << "Output is not semisorted\n";
            std::cout << "Segment with key=" << curr << " ended with " << remaining << " missing elements\n";
            std::cout << "Found " << (counts[x]-remaining) << " elements but expected " << counts[x] << " elements\n";
            exit(0);
          }
          curr = x;
          remaining = counts[x];
        }
        --remaining;
      }
    }

    if (i == 0) {
      std::cout << "Warmup round running time: " << t.total_time() << std::endl;
    } else {
      std::cout << "Round " << i << " running time: " << t.total_time()
                << std::endl;
      total_time += t.total_time();
    }
  }
  free(A);
  std::cout << "Average running time: " << total_time / num_rounds << std::endl;
}

int main(int argc, char* argv[]) {
  size_t n = 1e9;
  int num_rounds = 3;
  int random_seed = 1;
  if (argc >= 2) n = atoll(argv[1]);
  if (argc >= 3) num_rounds = atoi(argv[2]);
  if (argc >= 4) random_seed = atoi(argv[3]);  

  for (auto element_types : {1e1, 1e3, 1e5, 1e7, 1e9}) {
    std::cout << "Keys=" << element_types << std::endl;
    run_test(n, element_types, num_rounds, random_seed);
    std::cout << std::endl;
  }
  
  return 0;
}
