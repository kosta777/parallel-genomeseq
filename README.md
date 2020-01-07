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

### Building OpenMP version:

```bash
cd build
cmake .. -DUSEOMP=ON
make
cd ..
./bin/sw_solve_small #sw_solve_small is equipped with openmp-parallellocalaligner
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
  auto la = std::make_unique<SWAligner<Similarity_Matrix>>(str1,str2); //or SWAligner<Similarity_Matrix_Skewed>
  // or use parallel LocalAligner breaking ref into 4 substrings and setting overlaplength = samplelength * 2.0
  auto la = std::make_unique<OMPParallelLocalAligner<Similarity_Matrix_Skewed, SWAligner<Similarity_Matrix_Skewed>>>(row[2],fa_string,4,2.0);
  la->calculateScore();
  //... = la->get...();
}
//Constructor
template <class Similarity_Matrix_Type>
class SWAligner : public LocalAligner <Similarity_Matrix_Type> {
  public:
    SWAligner(std::string_view, std::string_view);
    SWAligner(std::string_view, std::string_view, std::function<float(const char &, const char &)> &&);
}
//API
template <class Similarity_Matrix_Type>
class LocalAligner {
  public:
    virtual float calculateScore() = 0;
    virtual float getScore() const = 0;
    virtual unsigned int getPos() const = 0;
    virtual std::string_view getConsensus_x() const = 0;
    virtual std::string_view getConsensus_y() const = 0;
    virtual const Similarity_Matrix_Type& getSimilarity_matrix() const =0;
};
template <class Similarity_Matrix_Type, class LocalAligner_Type>
class OMPParallelLocalAligner : public ParallelLocalAligner<Similarity_Matrix_Type, LocalAligner_Type> {
  public:
    OMPParallelLocalAligner(std::string_view, std::string_view, int, float);
    OMPParallelLocalAligner(std::string_view, std::string_view, int, float, float);
    OMPParallelLocalAligner(std::string_view, std::string_view, int, float, std::function<float(const char &, const char &)> &&);
    OMPParallelLocalAligner(std::string_view, std::string_view, int, float, std::function<float(const char &, const char &)> &&, float);
    float calculateScore();
    float getScore() const;
    unsigned int getPos() const;
    std::string_view getConsensus_x() const;
    std::string_view getConsensus_y() const;
};
//API for similarity matrix
class Abstract_Similarity_Matrix{
  public:
    virtual void iterate(const std::function<float(const char &, const char &)> &scoring_function,
                         float gap_penalty) = 0;
    virtual std::tuple<Eigen::Index, Eigen::Index, float> find_index_of_maximum() const = 0;
    virtual void print_matrix() const = 0;
    virtual const Eigen::MatrixXf &get_matrix() const = 0;
    virtual float operator()(Eigen::Index row, Eigen::Index col) const = 0;
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
make
../bin/tests

#Coarse-grained MPI
cd build
cmake -DUSEMPI=ON ..
make
#do the next part once to prepare input for MPI
cd ../py
python reader.py mpi_prepare_input
cd ..
mpiexec -np {node_num, ie. 6} ./bin/mpi_sw_solve_small

#MPI Benchmark on UNIPROT
cd build
cmake -DUSEMPI=ON ..
make
#do the next part once to prepare input
cd ../py
python reader.py uniprot_prepare
cd ..
mpiexec -np {node_num, ie. 6} ./bin/mpi_sw_solve_uniprot

#Coarse-grained OMP
cd build
cmake -DUSEOMP=ON ..
make
cd ..
./bin/sw_solve_small 
#create fixed number of subtasks by breaking the reference string into several pieces,
#number of subtasks and length of overlaps are controlled by pLocalAligner constructor:
#auto la = std::make_unique<OMPParallelLocalAligner<Similarity_Matrix_Skewed, SWAligner<Similarity_Matrix_Skewed>>>(row[2],fa_string,4,2.0);
#openmp automatically distribute subtasks to available processors
#number of threads can be controlled by setting env like (execute in bash): `export OMP_NUM_THREADS=16`

#Fine-grained OMP
cd build
cmake -DUSEOMP=ON ..
make
cd ..
bin/omp_sw_solve_small solve_small n_reads n_threads
(e.g.) bin/omp_sw_solve_small solve_small 3 2
```

## Benchmarking Specific Usage
### Fine-grained OMP
```bash
sh benchmark/leonhardsetup.sh

#download chromosome22 of hg19 reference from USCS database
sh benchmark/ompfg/data_get_chr22.sh

module load gcc/4.8.2 python/3.6.1
#Generate a custom reference from a section of the downloaded hg19 chr22 and a set of custom reads
#Default: |ref|=30k, n_reads=100, |read|=100
python py/ompfg_data_prep.py  --option gen_ref_custom
#Reference default output: data/custom_ref_1.fa
python py/ompfg_data_prep.py  --option gen_reads_custom
#Reference default input: data/custom_ref_1.fa
#Reads default output: data/custom_reads_1.csv

#examples:
python py/ompfg_data_prep.py  --option gen_ref_custom --start_pos 0 --ref_len 50000000 --remove_N=true
python py/ompfg_data_prep.py  --option gen_reads_custom --read_len=? --n_reads=?

rm -rf build
mkdir build
cd build
module load gcc
module load cmake

#Original ompfg benchmark
cmake .. -DUSEOMP=ON
#ompfg with multithreaded SIMD benchmark
cmake .. -DUSEOMP=ON -DMTSIMD=ON

make
cd ..

#Run OMP finegrain Smith-Waterman with different n_threads settings, append times to a csv
#arguments to the following command are specified in the sh file. Check for correctness if necessary
#cmd=$project_dir"/bin/omp_sw_solve_small "$code_section" "$n_reads" "$n_threads" "$finegrain_type" "$timing_file_path" "$ref_file_path" "$reads_file_path" "$mt_simd

#Original ompfg benchmark
sh benchmark/ompfg/ompfg_bench.sh
#ompfg with multithreaded SIMD benchmark
sh benchmark/ompfg/ompfg_mtsimd_bench.sh

#Visualization of numerical results in csv
module load gcc/4.8.2 python/3.6.1
python py/eval.py --option ompfg --yaxis gcups --plot_type box_plot --fit false
```
