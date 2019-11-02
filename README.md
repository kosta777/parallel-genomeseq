# Parallel Genome Sequencing
In bioinformatics, algorithms that compare strings for similarity play a decisive role in genome sequencing. Although today’s research is dominated by programs like BLAST that produce approximate solutions on large datasets, they fundamentally rely on dynamic algorithms such as:

1. Smith-Waterman algorithm
2. Needleman-Wunsch algorithm
3. Hirschberg's algorithm

These three closely related algorithms all fall under the category of dynamic programming algorithms for finding the edit distance of two strings (i.e. a measure for their dissimilarity). Our project consists of exploring different methods for parallelization of these fundamental algorithms while using different techniques (MPI, threading, combination etc.).


## Installation

The project is built with CMake.

### Building a release version:

```bash
cd build
cmake ..
make
```

### Building a debug version:

```bash
cd build
cmake .. -DDEBUG=ON
```

## Usage

Run the binary `parsequal` in the `bin` directory.

-Solving the small dataset and evaluate output using py script:

```bash
parseqal sw_solve_small
python eval.py sw_solve_small
```
