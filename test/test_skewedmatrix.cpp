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
  auto true_matrix = Eigen::MatrixXd(len_x, len_y);
  //auto true_matrix = Eigen::Map<Eigen::MatrixXd>(initvec.data(), len_x, len_y);
  for (size_t j = 0; j < len_y; j++) {
    for (size_t i = 0; i < len_x; i++) {
      auto idx1 = index_tuple(i,j);
      auto idx2 = index_tuple(j,i);
      ASSERT_EQ(idx1,skewed1.rawindex2trueindex(skewed1.trueindex2rawindex(idx1)));
      ASSERT_EQ(idx2,skewed2.rawindex2trueindex(skewed2.trueindex2rawindex(idx2)));
      ASSERT_EQ(idx2,skewed2.trueindex2rawindex(skewed2.rawindex2trueindex(idx2)));
      //raw_matrix(ri,rj) = true_matrix(i,j);
#ifdef VERBOSE
      auto [ri,rj] = skewed1.trueindex2rawindex(idx1);
      true_matrix(i,j) = i + j;
      raw_matrix(ri,rj) = i + j;
      //std::cout << "True index: (" << i <<","<< j << "), Raw index: (" << ri <<","<< rj << ")" << std::endl;
#endif
    }
  }
#ifdef VERBOSE
  std::cout << "Before Transform:\n" << true_matrix << std::endl;
  std::cout << "After Transform:\n" <<raw_matrix << std::endl;
#endif
}
