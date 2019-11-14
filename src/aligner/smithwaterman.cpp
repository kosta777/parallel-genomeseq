#include <iostream>
#include <Eigen/Dense>
#include <utility>
#include "smithwaterman.h"

template <class SMT>
SWAligner<SMT>::SWAligner(std::string_view first_sequence, std::string_view second_sequence) :
    SWAligner(first_sequence, second_sequence, [](const char &a, const char &b) { return a == b ? 3.0 : -3.0; }, 2.0) {}

template <class SMT>
SWAligner<SMT>::SWAligner(std::string_view first_sequence,
                          std::string_view second_sequence,
                          double gap_penalty) :
    SWAligner(first_sequence, second_sequence, [](const char &a, const char &b) { return a == b ? 3.0 : -3.0; }, gap_penalty) {}

template <class SMT>
SWAligner<SMT>::SWAligner(std::string_view first_sequence,
                          std::string_view second_sequence,
                          std::function<double(const char &, const char &)> &&scoring_function) :
    SWAligner(first_sequence, second_sequence, std::move(scoring_function), 2.0) {}

template <class SMT>
SWAligner<SMT>::SWAligner(std::string_view first_sequence,
                     std::string_view second_sequence,
                     std::function<double(const char &, const char &)> &&scoring_function, double gap_penalty) :
    pos(0),
    max_score(-1),
    gap_penalty(gap_penalty),
    sequence_x(first_sequence),
    sequence_y(second_sequence),
    similarity_matrix(first_sequence, second_sequence),
    scoring_function(std::move(scoring_function)) {}

template <class SMT>
void SWAligner<SMT>::traceback(index_tuple idx, unsigned int &preliminary_pos) {
  // stopping criterion
  auto [index_x, index_y] = idx;
  if (similarity_matrix(index_x - 1, index_y - 1) == 0) {
    sequence_x.insert(index_x - 1, "(");
    sequence_y.insert(index_y - 1, "(");

    pos = index_y;
    return;
  }

  auto n1 = similarity_matrix(index_x - 1, index_y - 1);
  auto n2 = similarity_matrix(index_x, index_y - 1);
  auto n3 = similarity_matrix(index_x - 1, index_y);

  // move north-west
  if ((n1 >= n2) && (n1 >= n3)) {
    index_tuple new_idx(index_x - 1, index_y - 1);
    traceback(new_idx, preliminary_pos);
  }
    // move west
  else if ((n2 >= n1) && (n2 >= n3)) {
    sequence_x.insert(index_x, "-");
    index_tuple new_idx(index_x, index_y - 1);
    traceback(new_idx, preliminary_pos);
  }
    // move north
  else {
    sequence_y.insert(index_y, "-");
    index_tuple new_idx(index_x - 1, index_y);
    traceback(new_idx, preliminary_pos);
  }
}

template <class SMT>
double SWAligner<SMT>::calculateScore() {
  similarity_matrix.iterate(scoring_function, gap_penalty);

  auto [index_x, index_y, max] = similarity_matrix.find_index_of_maximum();
  max_score = max;

  sequence_x.insert(index_x, ")");
  sequence_y.insert(index_y, ")");

  traceback(index_tuple(index_x, index_y), pos);
#ifdef VERBOSE
  raw_matrix.print_matrix();
  std::cout << "POS: " << pos << std::endl;

  std::cout << sequence_x << std::endl;
  std::cout << sequence_y << std::endl;
#endif
  return max;
}

//DO NOT FORGET TO INSTANTIATE the template class
template class SWAligner<Similarity_Matrix>;
template class SWAligner<Similarity_Matrix_Skewed>;