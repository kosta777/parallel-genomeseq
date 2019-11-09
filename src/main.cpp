#include <iostream>
#include <string>
#include <algorithm> 
#include <vector>

#include "aligner/localaligner.h"
#include "aligner/smithwatermanlinear.h"


int main() {
  // https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm#/media/File:Smith-Waterman-Algorithm-Example-Step2.png
  std::string sequence_x = "GGTTGACTA";
  std::string sequence_y = "TGTTACGG";

  LocalAligner *la = new SWAligner(sequence_x, sequence_y);
  la->calculateScore();

  return 0;
}
