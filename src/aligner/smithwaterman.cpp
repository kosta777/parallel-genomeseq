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
                          float gap_penalty) :
    SWAligner(first_sequence, second_sequence, [](const char &a, const char &b) { return a == b ? 3.0 : -3.0; }, gap_penalty) {}

template <class SMT>
SWAligner<SMT>::SWAligner(std::string_view first_sequence,
                          std::string_view second_sequence,
                          std::function<float(const char &, const char &)> &&scoring_function) :
    SWAligner(first_sequence, second_sequence, std::move(scoring_function), 2.0) {}

template <class SMT>
SWAligner<SMT>::SWAligner(std::string_view first_sequence,
                     std::string_view second_sequence,
                     std::function<float(const char &, const char &)> &&scoring_function, float gap_penalty) :
    sm_timings(2),
    pos(0),
    max_score(-1),
    gap_penalty(gap_penalty),
    sequence_x(first_sequence),
    sequence_y(second_sequence),
    consensus_x(),
    consensus_y(),
    similarity_matrix(first_sequence, second_sequence),
    scoring_function(std::move(scoring_function)) {
  consensus_x.reserve(sequence_x.size());
  consensus_y.reserve(sequence_x.size());
}

template <class SMT>
void SWAligner<SMT>::traceback(index_tuple similarity_matrix_max) {
  // stopping criterion
  auto [index_x, index_y] = similarity_matrix_max;

  while(true) {
    auto n1 = similarity_matrix(index_x - 1, index_y - 1);
    auto n2 = similarity_matrix(index_x, index_y - 1);
    auto n3 = similarity_matrix(index_x - 1, index_y);

    // stopping criterion
    if (n1 == 0 || n2 == 0 || n3 == 0) {
      consensus_x.push_back( sequence_x[index_x - 1] );
      consensus_y.push_back( sequence_y[index_y - 1] );
      pos = index_y;
      break;
    }

    // move north-west
    if ((n1 >= n2) && (n1 >= n3)) {
      consensus_x.push_back( sequence_x[index_x - 1] );
      consensus_y.push_back( sequence_y[index_y - 1] );
      index_x -= 1;
      index_y -= 1;
    }
    // move west
    else if ((n2 >= n1) && (n2 >= n3)) {;
      consensus_x.push_back( '-' );
      consensus_y.push_back( sequence_y[index_y - 1] );
      index_y -= 1;
    }
    // move north
    else {
      consensus_x.push_back( sequence_x[index_x - 1] );
      consensus_y.push_back( '-' );
      index_x -= 1;
    }
  }
}

template <class SMT>
float SWAligner<SMT>::calculateScore() {
#ifdef USEOMP
  similarity_matrix.sm_nthreads = sw_nthreads;
  similarity_matrix.sm_finegrain_type = sw_finegrain_type;
#endif

  similarity_matrix.iterate(scoring_function, gap_penalty);
#ifdef MTSIMD
  sw_mt_simd = similarity_matrix.sm_mt_simd;
#endif
  sm_timings = similarity_matrix.getTimings();

#ifdef USEOMP
//  sw_iter_ad_i_times_sum = (similarity_matrix.sm_iter_ad_i_times).sum();
//  sw_iter_ad_read_time = similarity_matrix.sm_iter_ad_read_time;
#endif
  auto [index_x, index_y, max] = similarity_matrix.find_index_of_maximum();
  max_score = max;
  traceback(index_tuple(index_x, index_y));
#ifdef VERBOSE
  similarity_matrix.print_matrix();
  std::cout << "POS: " << pos << std::endl;

  std::cout << consensus_x << std::endl;
  std::cout << consensus_y << std::endl;
#endif
  return max;
}

//DO NOT FORGET TO INSTANTIATE the template class
template class SWAligner<Similarity_Matrix>;
template class SWAligner<Similarity_Matrix_Skewed>;
