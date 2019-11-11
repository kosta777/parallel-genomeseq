#ifdef USEMPI
#include <gtest/gtest.h>
#include <mpi.h>

TEST(MPI, Usability) {

  // Initialisation
  int *argc={};
  char ***argv={};
  MPI_Init(argc, argv);

  // Reading size and rank
  int size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // Printing
  std::cout << "Hello world, from process #" << rank << std::endl;

  // Finalisation
  MPI_Finalize();

}
#endif

