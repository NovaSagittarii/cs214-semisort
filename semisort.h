#include <algorithm>
#include <math.h>

#include "parlaylib/include/parlay/parallel.h"

// hash function
template <class T>
T h(T x) {
    return x;
}

template <class T>
void semisort(T *A, size_t n) {
    std::sort(A, A+n);
}