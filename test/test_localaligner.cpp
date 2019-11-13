#include <gtest/gtest.h>
#include <memory>
#include "localaligner.h"
#include "smithwaterman.h"
#include "similaritymatrix.h"

TEST(LocalAligner,SmallSeqAlignment) {
  // https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm#/media/File:Smith-Waterman-Algorithm-Example-Step2.png
  std::string sequence_x = "GGTTGACTA";
  std::string sequence_y = "TGTTACGG";
  {
    auto la = std::make_unique<SWAligner<Similarity_Matrix_Skewed>>(sequence_x, sequence_y);
    auto maxval = la->calculateScore();
    auto pos = la->getPos();
    ASSERT_EQ(maxval,13);
    ASSERT_EQ(pos,2);
  }
}
