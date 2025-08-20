/**
 * A C++ example to compute the Frobenius norm of an NxN matrix with
 * double values.
 *
 * Compile with:   g++ matrix-demo.cpp
 * Run with:       ./a.out
 */

#include <cmath>
#include <iostream>
#include <vector>

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
frob(matrix &mat)
{
  double result = 0.0;
  for (int i = 0; i < mat.N; ++i)
    for (int j = 0; j < mat.N; ++j)
      {
        // try what happens when you do mat(j,i) instead
        double t = std::abs(mat(i, j));
        result += t * t;
      }
  return std::sqrt(result);
}

int
main()
{
  matrix mat(8000);
  fill(mat);

  Timer<> clock; // Timer is defined in timer.h in the same folder

  // Let's time it:
  clock.tick();
  double result = frob(mat);
  clock.tock();

  std::cout << result << std::endl;
  std::cout << "took " << clock.duration().count() << " ms." << std::endl;
}
