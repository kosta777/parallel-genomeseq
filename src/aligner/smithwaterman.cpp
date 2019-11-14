#include <iostream>
#include <Eigen/Dense>
#include <utility>
#include "smithwaterman.h"

SWAligner::SWAligner(std::string_view first_sequence, std::string_view second_sequence) :
    pos(0),
    max_score(-1),
    sequence_x(first_sequence),
    sequence_y(second_sequence),
    similarity_matrix(first_sequence, second_sequence),
    scoring_function([](const char &a, const char &b) { return a == b ? 3.0 : -3.0; }) {}

SWAligner::SWAligner(std::string_view first_sequence,
                     std::string_view second_sequence,
                     std::function<double(const char &, const char &)> &&scoring_function) :
    pos(0),
    max_score(-1),
    sequence_x(first_sequence),
    sequence_y(second_sequence),
    similarity_matrix(first_sequence, second_sequence),
    scoring_function(std::move(scoring_function)) {}

void SWAligner::traceback(index_tuple similarity_matrix_max) {
  auto matrix = similarity_matrix.get_matrix();

  auto [index_x, index_y] = similarity_matrix_max;

  while(true) {
    auto n1 = matrix(index_x - 1, index_y - 1);
    auto n2 = matrix(index_x, index_y - 1);
    auto n3 = matrix(index_x - 1, index_y);

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

void SWAligner::calculate_similarity_matrix() {
  auto cb = [&sequence_x = sequence_x, &sequence_y = sequence_y, &
      scoring_function = scoring_function](const Eigen::MatrixXd &matrix, index_tuple idx) {
    auto a = sequence_x[idx.first - 1];
    auto b = sequence_y[idx.second - 1];
    auto west = matrix(idx.first, idx.second - 1);
    auto north = matrix(idx.first - 1, idx.second);
    auto north_west = matrix(idx.first - 1, idx.second - 1);
    auto gap_penalty = 2;
    Eigen::Vector4d v{
        north_west + scoring_function(a, b),
        west - gap_penalty,
        north - gap_penalty,
        0
    };
    auto result = v.maxCoeff();
    return result;
  };
  similarity_matrix.iterate_anti_diagonal(cb);
}

double SWAligner::calculateScore() {
  calculate_similarity_matrix();

  auto max_idx = similarity_matrix.find_index_of_maximum();

  traceback(max_idx);
#ifdef VERBOSE
  similarity_matrix.print_matrix();
  std::cout << "POS: " << pos << std::endl;

  std::cout << sequence_x << std::endl;
  std::cout << sequence_y << std::endl;
#endif
  max_score = similarity_matrix.get_matrix()(max_idx.first, max_idx.second);
  return max_score;
}
