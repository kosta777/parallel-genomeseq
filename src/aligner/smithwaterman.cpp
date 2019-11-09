#include <iostream>
#include <Eigen/Dense>
#include <utility>
#include "smithwaterman.h"

SWAligner::SWAligner(std::string_view first_sequence, std::string_view second_sequence) :
    pos(0),
    sequence_x(first_sequence),
    sequence_y(second_sequence),
    similarity_matrix(first_sequence, second_sequence),
    scoring_function([](const char &a, const char &b) { return a == b ? 3.0 : -3.0; }) {}

SWAligner::SWAligner(std::string_view first_sequence,
                     std::string_view second_sequence,
                     std::function<double(const char &, const char &)> &&scoring_function) :
    pos(0),
    sequence_x(first_sequence),
    sequence_y(second_sequence),
    similarity_matrix(first_sequence, second_sequence),
    scoring_function(std::move(scoring_function)) {}

void SWAligner::traceback(index_tuple idx, unsigned int &preliminary_pos) {
  const Eigen::MatrixXd &matrix = similarity_matrix.get_matrix();

  // stopping criterion
  auto [index_x, index_y] = idx;
  if (matrix(index_x - 1, index_y - 1) == 0) {
    sequence_x.insert(index_x - 1, "(");
    sequence_y.insert(index_y - 1, "(");

    pos = index_y;
    return;
  }

  double n1 = matrix(index_x - 1, index_y - 1);
  double n2 = matrix(index_x, index_y - 1);
  double n3 = matrix(index_x - 1, index_y);

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

void SWAligner::calculate_similarity_matrix() {
  auto cb = [&sequence_x = sequence_x, &sequence_y = sequence_y, &
      scoring_function = scoring_function](const Eigen::MatrixXd &matrix, index_tuple idx) {
    const char a = sequence_x[idx.first - 1];
    const char b = sequence_y[idx.second - 1];
    const double west = matrix(idx.first, idx.second - 1);
    const double north = matrix(idx.first - 1, idx.second);
    const double north_west = matrix(idx.first - 1, idx.second - 1);
    const double gap_penalty = 2;
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

  index_tuple max_idx = similarity_matrix.find_index_of_maximum();

  sequence_x.insert(max_idx.first, ")");
  sequence_y.insert(max_idx.second, ")");

  traceback(max_idx, pos);
#ifdef DEBUG
  similarity_matrix.print_matrix();
  std::cout << "POS: " << pos << std::endl;

  std::cout << sequence_x << std::endl;
  std::cout << sequence_y << std::endl;
#endif
  return similarity_matrix.get_matrix()(max_idx.first, max_idx.second);
}