/* Consider an MPI parallel program running on K ranks with an array of N strings split evenly between
the K ranks (each rank only knows its own part). You may assume there are no duplicates, even
between different ranks.
(a) All ranks call a function “find owner of(X)“ that returns (on all ranks!) the integer (between 0
and K-1) of the unique rank who has this string X in their array. Assume it exists on exactly
one rank. Detail the necessary steps for an efficient implementation of this function including
communication (be precise about who sends/receives what!).
*/

// compile using "mpicxx main.cc" and run with
// "mpirun -n 3 ./main"

#include <mpi.h>

#include <algorithm>
#include <cstddef>
#include <iostream>
#include <numeric>
#include <vector>
#include <random>

bool find(const std::vector<int> &data, int key)
{
  return std::find(data.begin(), data.end(), key) != data.end();
}

void shuffle(std::vector<int> &data)
{
    std::random_device rd;
    std::mt19937 g(rd());
 
    std::shuffle(data.begin(), data.end(), g);
}

int find_owner(const std::vector<int> &data, int key)
{
  bool found = find(data, key);

  return -1;
}


int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int  name_len;
  MPI_Get_processor_name(processor_name, &name_len);

  std::cout << "Hi, I am process " << rank << " of " << size
            << " and I am running on " << processor_name << std::endl;


  const std::size_t N = 1000;
  std::vector<int> data(N/size);
  std::iota(data.begin(), data.end(), 1);
  
  shuffle(data);

  int key = 42;
  int owner = find_owner(data, key);

  if (rank == 0)
    std::cout << "Owner of " << key << " is " << owner << std::endl;


  // if (rank == 0)
  //   {
  //     for (int i = 1; i < size; ++i)
  //       {
  //         double     value = 0.0;
  //         MPI_Status status;

  //         // Use MPI_ANY_SOURCE or i as the source:
  //         MPI_Recv(&value,         // void *buf
  //                  1,              // int count
  //                  MPI_DOUBLE,     // MPI_Datatype datatype,
  //                  MPI_ANY_SOURCE, // int source
  //                  0,              // int tag
  //                  MPI_COMM_WORLD,
  //                  &status); // MPI_Status *status

  //         std::cout << "I got value = " << value << " from "
  //                   << status.MPI_SOURCE << std::endl;
  //       }
  //   }
  // else
  //   {
  //     // std::this_thread::sleep_for(std::chrono::milliseconds(3000));
  //     srand(rank);
  //     double my_value = (rand() % 1000) / 1000.0;

  //     MPI_Send(&my_value,  // void *buf
  //              1,          // int count
  //              MPI_DOUBLE, // MPI_Datatype datatype
  //              0,          // int dest
  //              0,          // int tag
  //              MPI_COMM_WORLD);
  //   }

  MPI_Finalize();
}
