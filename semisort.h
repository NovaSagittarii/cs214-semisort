#include <algorithm>
#include <map>
#include <math.h>
#include <utility>

#include "parlaylib/include/parlay/parallel.h"

#include "hash64.h"
#include "scan.h"

// parameters (move to function later?)
// #define nL    1024      // number of light buckets
#define alpha 32*1024   // base case size

size_t GetKey(const std::pair<long long, long long> &el) {
  return el.first;
}
// hash function
template <class T>
size_t h(T x) {
  // return x;
  return hash64(x);
  // return hash64(x) & (((size_t)1 << 63) - 1);
  // https://graphics.stanford.edu/~seander/bithacks.html#IntegerAbs
  // auto v = hash64(x);
  // size_t const mask = (v >> sizeof(size_t)) - 1;
  // return (v + mask) ^ mask;
}

#define print(x) std::cout<< x
#define println(x) print(x << std::endl)
#define print_array(A, n) for (size_t i = 0; i < n; ++i) {print(A[i] << " ");} println("");

size_t GetBucketId(const size_t k, const std::map<size_t, size_t>& H, const size_t nL) {
  if (H.count(k)) return H.at(k);
  return k % nL;
}

template <class T>
void ssort(T *A, T *B, size_t *C, size_t *SC, const size_t An, const size_t n, const bool in_place) {
  if (n <= alpha) {
    std::sort(A, A + n);
    if (in_place) {
      parlay::parallel_for(0, n, [&](size_t i) { // copy step; TODO: remove later
        B[i] = A[i];
      });
    }
    return;
  }
  if (n <= 1) return;

  // Parameters
  const int target_nL = std::pow(n, 0.5)/2;
  const size_t nL = std::min(1024, target_nL);
  const size_t l = n / target_nL; // subarray size
  const size_t m = (n + l - 1) / l; // subarrays

  // println("n'=" << n << " l=" << l << " subarrays=" << m);
  // print_array(A, n);

  // Sampling and Bucketing
  const size_t Sn = nL * std::log(n);
  std::map<size_t, size_t> S;
  for (size_t i = 0; i < Sn; ++i) {
    ++S[h(GetKey(A[hash64(i + Sn) % Sn]))];
  }

  std::map<size_t, size_t> H;
  {
    size_t id = nL;
    for (const auto [k, v] : S) {
      if (v >= std::log(n)) {
        // Assign bucket id to heavy key k
        H[k] = id;
        ++id;
      }
    }
  }
  const size_t nH = H.size();
  const size_t nB = nL + nH; // bucket count (needed for matrix dimensions)
  // print(nH << " heavy bucket(s): ");
  // for (auto [k,v]: H) print(k<<" ("<<GetBucketId(k, H)<<") ");
  // println("");
  const bool needs_extra = nB * m > n;
  if (needs_extra) {
    std::cout << "parameters messed up, using more memory than allocated" << std::endl;
    std::cout << "requires " << nB*m << "B but has " << n << "B" << std::endl;
  }
  //  nB * m = n
  // 2nL * m = n
  //       m = n / (2nL)

  // Blocked Distributing
    // (nL+nH) * (n'/l) <= 2nL * (n/l)
    // since l >> 2nL, O(n) extra space is sufficient
    // 2nL * n/l <= n
    // 2nL <= l
  parlay::parallel_for(0, n, [&](size_t i) { // initialize matrix C
    C[i] = 0;
  });
  parlay::parallel_for(0, m, [&](size_t i) { // to determine bucket size
    for (size_t j = i*l; j < std::min((i+1)*l, n); ++j) {
      const size_t id = GetBucketId(h(GetKey(A[j])), H, nL);
      // C[i][id] = C[id * m + i] = records falling into bucket id in subarray i
      ++C[id * m + i];
    }
  });
    // println("C = ");
    // for (size_t i = 0; i < nB; ++i) {
    //   auto D = C+m*i;
    //   print_array(D, m);
    // }
    // println("===");
  C[m * nB] = scan(C, SC, m * nB); // computing the prefixsum to get exact bucket sizes
  // println("C[-1]="<< C[m * nB]);
    // println("X = ");
    // for (size_t i = 0; i < nB; ++i) {
    //   auto D = C+m*i;
    //   print_array(D, m);
    // }
    // println("===");

  // C[0][id] = C[id * m] is the offset of each bucket
  // size_t* offsets = C + (l * (nL + nH)) + 1; // subproblem partition messes up this data
  size_t* offsets = (size_t*)malloc((nL + 1) * sizeof(size_t));
  parlay::parallel_for(0, nL + 1, [&](size_t i) {
    offsets[i] = C[i * m];
  });
  // print_array(offsets, nL+nH+1);
  // println("aux= " << nL+nH+1 + (l*(nL+nH)) + 1 << " n=" << n);

  // initialize B to make mapping fails more obvious
  // parlay::parallel_for(0, n, [&](size_t i) {
  //   B[i] = -1;
  // });

  parlay::parallel_for(0, m, [&](size_t i) { // distributing items into buckets
    for (size_t j = i*l; j < std::min((i+1)*l, n); ++j) {
      const size_t id = GetBucketId(h(GetKey(A[j])), H, nL);
      B[C[id * m + i]] = A[j];
      ++C[id * m + i];
    }
  });

  // print_array(A, n);
  // print_array(B, n);
  // print_array(offsets, nL+1);
  if (!in_place) {
    parlay::parallel_for(offsets[nL], n, [&](size_t i) { // write heavy buckets back
      A[i] = B[i];
    });
  }

  // Local Refining
  // const size_t memory_block = mem / nL;
  parlay::parallel_for(0, nL, [&](size_t i) {
    const size_t o = offsets[i];
    ssort(B + o, A + o, C + o, SC + o, An, offsets[i+1] - o, !in_place);
  });
  free(offsets);
}

template <class T>
void semisort(T *A, size_t n) {
  T* B = (T*)malloc(n * sizeof(T)); // output buffer
  size_t* C = (size_t*)malloc(n * sizeof(size_t)); // subarray counts
  size_t* SC = (size_t*)malloc(n * sizeof(size_t)); // scan buffer
  // print_array(A, n);
  ssort(A, B, C, SC, n, n, false);
  // print_array(A, n);
  free(C);
  free(B);
}