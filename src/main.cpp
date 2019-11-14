#include <iostream>
#include <string>
#include <algorithm>
#include <vector>
#include <memory>

#include "localaligner.h"
#include "smithwaterman.h"
#include "similaritymatrix.h"

int main() {
  // https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm#/media/File:Smith-Waterman-Algorithm-Example-Step2.png
  std::string sequence_x = "GGTTGACTA";
  std::string sequence_y = "TGTTACGG";
  {
    auto la = std::make_unique<SWAligner<Similarity_Matrix_Skewed>>(sequence_x, sequence_y);
    la->calculateScore();
  }
  return 0;
}
