# Parallel Genome Sequencing
In bioinformatics, algorithms that compare strings for similarity play a decisive role in genome sequencing. Although today’s research is dominated by programs like BLAST that produce approximate solutions on large datasets, they fundamentally rely on dynamic algorithms such as:

1. Smith-Waterman algorithm
2. Needleman-Wunsch algorithm
3. Hirschberg's algorithm

These three closely related algorithms all fall under the category of dynamic programming algorithms for finding the edit distance of two strings (i.e. a measure for their dissimilarity). Our project consists of exploring different methods for parallelization of these fundamental algorithms while using different techniques (MPI, threading, combination etc.).


## Installation

The project is built with CMake. Available options:
1. `-DDEBUG=ON` to compile in debug mode
2. `-DVERBOSE=ON` to turn on verbose output
3. `-DUSEMPI=ON` to use MPI (requires openMPI library installed)
4. `-DUSEOMP=ON` to use OpenMP (requires openMP library installed)

### Building a release version:

```bash
cd build
cmake ..
make
```

### Building a debug version:

```bash
cd build
cmake .. -DDEBUG=ON #-DVERBOSE=ON #if need verbose output
```

### Building MPI version:

```bash
cd build
cmake .. -DUSEMPI=ON
make
cd ../bin
mpiexec parseq
```

## Development
Please wrap anything using MPI by `#ifdef USEMPI ... #endif` and anything using OpenMP by `#ifdef USEOMP ... #endif`.
The rule of thumb is do not break serial code when MPI and openMP is not available.
```C++
//Usage
#include <memory>
{
  auto la = std::make_unique<SWAligner>(str1,str2);
  la->calculateScore();
  \\... = la->get...();
}
//Constructor
class SWAligner : public LocalAligner {
  public:
    SWAligner(std::string_view, std::string_view);
    SWAligner(std::string_view, std::string_view, std::function<double(const char &, const char &)> &&);
}
//API
class LocalAligner {
  public:
    virtual double calculateScore() = 0;
    virtual double getScore() const = 0;
    virtual unsigned int getPos() const = 0;
    virtual std::string_view getConsensus_x() const = 0;
    virtual std::string_view getConsensus_y() const = 0;
    virtual const Similarity_Matrix& getSimilarity_matrix() const =0;
};
```

## Usage

Run the binary `parsequal` in the `bin` directory.

-Generate a modified input file with fields: index, QNAME, SEQ, POS (POS="true" alignment positions extracted from a SAM file), solve the small dataset using Smith-Waterman, and evaluate the alignment result wrt "true" positions using a py script:

```bash
#sw_solve_small
cd build
cmake ..
make sw_solve_small
cd ../py
python reader.py gen_input
cd ..
./bin/sw_solve_small
cd py
python eval.py sw_solve_small
#test
cd build
cmake ..
make test
../bin/tests
```
