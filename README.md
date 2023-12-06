# CS214 Project 2
> Parallel Semisort

## Results
This was run on a 48-core (Intel(R) Xeon(R) Silver 4214R CPU @ 2.40GHz) server with 755GB RAM. A test is run 4 times and the performance is the average of the last 3 runs (first run is a warm-up round). There are `1e9` elements in every test.

### = (equality with hashing)
| Types of Keys | paper | mine |
| --- | --- | --- |
| 1e9 | 2.47 | 4.55 |
| 1e7 | 2.20 | 4.93 |
| 1e5 | 1.59 | 5.68 |
| 1e3 | 1.56 | 6.53 |
| 1e1 | 1.35 | 2.56 |

### =-i (equality without hashing, using integer keys)
| Types of Keys | paper | mine |
| --- | --- | --- |
| 1e9 | 1.85 | 4.21 |
| 1e7 | 1.80 | 4.51 |
| 1e5 | 1.30 | 5.38 |
| 1e3 | 1.30 | 6.38 |
| 1e1 | 1.06 | 2.39 |

## Running
```sh
make && ./semisort
```

## References
Xiaojun Dong, Yunshu Wu, Zhongqi Wang, Laxman Dhulipala, Yan Gu, Yihan Sun. High-Performance and Flexible Parallel Algorithms for Semisort and Related Problems. *arXiv preprint: 2304.10078*, 2023.