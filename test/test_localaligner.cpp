#include <gtest/gtest.h>
#include <memory>
#include "localaligner.h"
#include "smithwaterman.h"
#include "similaritymatrix.h"

namespace {
  class SWAligner_Test : public testing::Test {
    protected:
      void SetUp() override {
        // Wikipedia Waterman-Smith example
        // https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm#/media/File:Smith-Waterman-Algorithm-Example-Step2.png
        std::string sequence_x = "GGTTGACTA";
        std::string sequence_y = "TGTTACGG";

        la = std::make_unique<SWAligner<Similarity_Matrix_Skewed>>(sequence_x, sequence_y);
        la->calculateScore();
      }

      std::unique_ptr<SWAligner<Similarity_Matrix_Skewed>> la;
  };


  TEST_F(SWAligner_Test, Example_small_sequence_alignment) {
    ASSERT_EQ(la->getScore(), 13);
    ASSERT_EQ(la->getPos(), 2);
  }

  /*
  TEST_F(SWAligner_Test, Example_output) {
    Eigen::MatrixXf refrence_m(10, 9);
    refrence_m <<
      0,  0,  0,  0,  0,  0,  0,  0,  0,
      0,  0,  3,  1,  0,  0,  0,  3,  3,
      0,  0,  3,  1,  0,  0,  0,  3,  6,
      0,  3,  1,  6,  4,  2,  0,  1,  4,
      0,  3,  1,  4,  9,  7,  5,  3,  2,
      0,  1,  6,  4,  7,  6,  4,  8,  6,
      0,  0,  4,  3,  5, 10,  8,  6,  5,
      0,  0,  2,  1,  3,  8, 13, 11,  9,
      0,  3,  1,  5,  4,  6, 11, 10,  8,
      0,  1,  0,  3,  2,  7,  9,  8,  7;


    ASSERT_PRED2(
        [](const Eigen::MatrixXf &lhs, const Eigen::MatrixXf &rhs) { return lhs.isApprox(rhs, 1e-4);},
        refrence_m,
        la->getSimilarity_matrix().get_matrix()
    );
  }
   */

  TEST_F(SWAligner_Test, Verify_consensus_strings) {
    std::string_view expected_consensus_x = "CAGTTG";
    std::string_view expected_consensus_y = "CA-TTG";

    ASSERT_EQ(la->getConsensus_x(), expected_consensus_x);
    ASSERT_EQ(la->getConsensus_y(), expected_consensus_y);
  }
}
