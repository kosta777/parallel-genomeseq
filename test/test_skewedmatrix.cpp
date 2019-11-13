#include <gtest/gtest.h>
#include <iostream>
#include "similaritymatrix.h"

TEST(SimilarityMatrix,SkewedMatrixIndex) {
  std::string sequence_x = "GGTTGACTA";
  std::string sequence_y = "TGTTACG";
  auto len_x = sequence_x.size() + 1;
  auto len_y = sequence_y.size() + 1;
  auto skewed1 = Similarity_Matrix_Skewed(sequence_x, sequence_y);
  auto skewed2 = Similarity_Matrix_Skewed(sequence_y, sequence_x);
  auto size = len_x * len_y;
  Eigen::VectorXd initvec = Eigen::VectorXd::LinSpaced(size, 0., size - 1);
  auto raw_matrix = Eigen::MatrixXd(len_y, len_x);
  auto true_matrix = Eigen::Map<Eigen::MatrixXd>(initvec.data(), len_x, len_y);
  for (size_t j = 0; j < len_y; j++) {
    for (size_t i = 0; i < len_x; i++) {
      auto idx1 = index_tuple(i,j);
      auto idx2 = index_tuple(j,i);
      ASSERT_EQ(idx1,skewed1.rawindex2trueindex(skewed1.trueindex2rawindex(idx1)));
      ASSERT_EQ(idx2,skewed2.rawindex2trueindex(skewed2.trueindex2rawindex(idx2)));
      ASSERT_EQ(idx2,skewed2.trueindex2rawindex(skewed2.rawindex2trueindex(idx2)));
      auto [ri,rj] = skewed1.trueindex2rawindex(idx1);
      raw_matrix(ri,rj) = true_matrix(i,j);
#ifdef VERBOSE
      std::cout << "True index: (" << i <<","<< j << "), Raw index: (" << ri <<","<< rj << ")" << std::endl;
#endif
    }
  }
#ifdef VERBOSE
  std::cout << raw_matrix << std::endl;
#endif
}

TEST(SimilarityMatrix, SkewedMatrixDP) {
  std::string sequence_x = "GGTTGACTA";
  std::string sequence_y = "TGTTACG";
  auto len_x = sequence_x.size() + 1;
  auto len_y = sequence_y.size() + 1;
  auto skewed = Similarity_Matrix_Skewed(sequence_x, sequence_y);
  auto normal = Similarity_Matrix(sequence_x, sequence_y);
  auto skewed2 = Similarity_Matrix_Skewed(sequence_y, sequence_x);
  auto normal2 = Similarity_Matrix(sequence_y, sequence_x);
  auto scoring_function = [] (const char& a, const char& b) {return a == b ? 3.0 : -3.0;};
  skewed.iterate_anti_diagonal(scoring_function, 2.0);
  normal.iterate_anti_diagonal(scoring_function, 2.0);
  skewed2.iterate_anti_diagonal(scoring_function, 2.0);
  normal2.iterate_anti_diagonal(scoring_function, 2.0);
  auto skewed_mat = skewed.get_matrix();
  auto normal_mat = normal.get_matrix();
  auto skewed2_mat = skewed2.get_matrix();
  auto normal2_mat = normal2.get_matrix();
  for (size_t j = 0; j < len_y; j++) {
    for (size_t i = 0; i < len_x; i++) {
      auto [ri, rj] = skewed.trueindex2rawindex(index_tuple(i, j));
      auto [rj2, ri2] = skewed2.trueindex2rawindex(index_tuple(j, i));
      ASSERT_EQ(normal_mat(i, j), skewed_mat(ri, rj));
      ASSERT_EQ(normal2_mat(j, i), skewed2_mat(rj2, ri2));
    }
  }
#ifdef VERBOSE
  std::cout << "Similarity_Matrix result:\n" << normal.get_matrix() << std::endl;
  std::cout << "Similarity_Matrix_Skewed result:\n" ;
  skewed.print_matrix();
#endif
}
