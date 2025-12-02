/*

Demonstrate SIMD vectorization, once by using AVX2 intrinsics frob_vec() and once using std::simd in frob_vec2().

We need to enable FMA support, native optimizations, and use a new g++:
g++ -O3 -march=native -mfma simd.cpp

Sadly, The performance using std::simd is significantly slower (N=2048):

frob_old:  0.7e6 MFlops
frob:      5.0e6 MFlops  (~ 7x, uses auto vec with SSE2, 2 doubles)
frob_vec:  7.7e6 MFlops  (~ 10x, vec with AVX2 and FMA, 4 doubles)
frob_vec2: 5.1e6 MFlops


Look at assembly:
g++ -O1 -march=native -mfma ../simd.cpp -S simd.s

*/

#include <cmath>
#include <iostream>
#include <vector>
#include <experimental/simd>

#include <immintrin.h> // for AVX2 intrinsics


#include "timer.h"

// A very simple matrix class
class matrix
{
public:
  matrix(int N)
    : N(N)
  {
    data.resize(N * N);
  }

  double &
  operator()(int i, int j)
  {
    return data[j * N + i];
  }

  int                 N;
  std::vector<double> data;
};

void
fill(matrix &mat)
{
  for (int i = 0; i < mat.N; ++i)
    for (int j = 0; j < mat.N; ++j)
      {
        mat(i, j) = 1.0 * (i + j);
      }
}

// compute the Frobenius norm
double
frob_old(matrix &mat)
{
  // Manual loop (i,j)

  double result = 0.0;
  for (int i = 0; i < mat.N; ++i)
    for (int j = 0; j < mat.N; ++j)
      {
        // try what happens when you do mat(j,i) instead
        // (i,j) is slow, (j,i) is fast
        //double t = std::abs(mat(j,i));
        double t = mat(i, j);
        result += t * t;
      }
  return std::sqrt(result);
}

double
frob(matrix &mat)
{
  // Unrolled to a single loop over i

  double result = 0.0;
  int limit = mat.N*mat.N;

  for (int i = 0; i < limit; ++i)
      {
        // try what happens when you do mat(j,i) instead
        // (i,j) is slow, (j,i) is fast
        //double t = std::abs(mat(j,i));
        double t = mat.data[i];
        result += t * t;
      }
  return std::sqrt(result);
}

double
frob_vec(matrix &mat)
{
  // Manually vectorize with AVX2 and use fused multiply add

  double result = 0.0;
  const int N = mat.N;
  const int block = 4; // AVX2 processes 4 doubles at a time (256 bits)
  int N2 = N * N;
  int limit = (N2 / block) * block;

  __m256d sum = _mm256_setzero_pd();

  double* data = mat.data.data();

  int i = 0;
  for (; i < limit; i += block)
    {
      __m256d v = _mm256_loadu_pd(&data[i]);
      sum = _mm256_fmadd_pd(v, v, sum); // sum += v*v
    }

  // Horizontal sum of the vector register
  double partial[4];
  _mm256_storeu_pd(partial, sum);
  result = partial[0] + partial[1] + partial[2] + partial[3];

  // tail elements (if N*N is not divisible by 4)
  for (; i < N2; ++i)
    {
      double t = data[i];
      result += t * t;
    }

  return std::sqrt(result);
}


double
frob_vec2(matrix &mat)
{
  // Now using std::simd, sadly, it does not give the performance as expected

  using namespace std::experimental;
  constexpr std::size_t simd_size = simd<double>::size();

  double result = 0.0;
  double* data = mat.data.data();
  const int N2 = mat.N * mat.N;
  const int limit = (N2 / simd_size) * simd_size;

  simd<double> v; 
  simd<double> sum = 0.0;

  int i = 0;
  for (; i < limit; i += simd_size)
    {
      v.copy_from(&data[i], element_aligned);
      //or: sum += v * v;
      sum = fma(v,v, sum);
    }

  result = reduce(sum);

  // tail elements
  for (; i < N2; ++i)
    {
      double t = data[i];
      result += t * t;
    }

  return std::sqrt(result);
}

void
test(int size)
{
  matrix mat(size);
  fill(mat);
  Timer<> clock; // Timer is defined in timer.h in the same folder

  // Let's time it:
  clock.tick();
  int    runs = 1000000 / size;
  double r    = 1.0;
  for (int i = 0; i < runs; ++i)
    r += frob(mat);
   //r += frob_vec(mat);
    // r += frob_vec2(mat);
  clock.tock();

  const double secs  = clock.duration().count() / 1000.0 / runs;
  const double bytes = 1.0 * sizeof(double) * size * size;
  const double flops = 3.0 * size * size;

  std::cout << "N=" << size
            << " mem: " << std::round(bytes / 1024 / 1024 * 100) / 100 << "mb"
            << " Mflops/s: " << std::round(flops / secs / 1000000)
            << " mem: " << std::round(bytes / secs / 1024 / 1024) << "mb/s"
            << " took " << secs << "s"
            << " runs: " << runs << " r=" << r << std::endl;
}

int
main(int argc, char *argv[])
{
  if (argc != 2)
    {
      std::cout << "Usage: \n\t" << argv[0] << " [size]\n\t" << argv[0]
                << " all\n";
      return 0;
    }

  if (std::string(argv[1]) == "all")
    {
      for (int i = 1; i < 30; ++i)
        {
          int size = 200 * i;
          test(size);
        }
    }
  else
    {
      int size = atoi(argv[1]);

      test(size);
    }
  return 1;
}
