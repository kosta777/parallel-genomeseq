#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <cmath>
#include <Eigen/Dense>

#include "swaligner.h"
#include "localaligner.h"

using namespace std;

void SWAligner::setFirstSequence(std::string seq) {
  this->seqA = seq;
}

void SWAligner::setSecondSequence(std::string seq) {
  this->seqB = seq;
}

void SWAligner::calculateScore() {
  // initialize some variables
  int lengthSeqA = seqA.length();
  int lengthSeqB = seqB.length();

  // initialize matrix
  auto matrix = Eigen::MatrixXd(lengthSeqA + 1, lengthSeqB + 1);
  matrix.setZero();

  auto I_i = Eigen::MatrixXi(lengthSeqA + 1, lengthSeqB + 1);
  auto I_j = Eigen::MatrixXi(lengthSeqA + 1, lengthSeqB + 1);
  auto traceback = Eigen::VectorXd(4);

  int ind; //max index of traceback

  //start populating matrix
  //https://github.com/ding-lab/pindel2/blob/master/src/smith_waterman_alignment.cpp
  for (int i = 1; i <= lengthSeqA; i++) {
    for (int j = 1; j <= lengthSeqB; j++) {
#ifdef DEBUG
      cout << i << " " << j << endl;
#endif
      traceback(0) = matrix(i - 1, j - 1) + similarityScore(seqA[i - 1], seqB[j - 1]);
      traceback(1) = matrix(i - 1, j) + penalty;
      traceback(2) = matrix(i, j - 1) + penalty;
      traceback(3) = 0;
      matrix(i, j) = traceback.maxCoeff(&ind);
      switch (ind) {
        case 0:I_i(i, j) = i - 1;
          I_j(i, j) = j - 1;
          break;
        case 1:I_i(i, j) = i - 1;
          I_j(i, j) = j;
          break;
        case 2:I_i(i, j) = i;
          I_j(i, j) = j - 1;
          break;
        case 3:I_i(i, j) = i;
          I_j(i, j) = j;
          break;
        default:throw; //ind<=3
      }
    }
  }

#ifdef DEBUG
  // print the scoring matrix to console
  cout << matrix << endl;
#endif
  // find the max score in the matrix
  int i_max = 0, j_max = 0;
  matrix_max = matrix.maxCoeff(&i_max, &j_max);

#ifdef DEBUG
  cout << "Max score in the matrix is " << matrix_max << endl;
#endif

  // traceback

  int current_i = i_max, current_j = j_max;
  int next_i = I_i(current_i, current_j);
  int next_j = I_j(current_i, current_j);
  int tick = 0;
  char consensus_a[lengthSeqA + lengthSeqB + 2], consensus_b[lengthSeqA + lengthSeqB + 2];

  while (((current_i != next_i) || (current_j != next_j)) && (next_j != 0) && (next_i != 0)) {

    if (next_i == current_i) consensus_a[tick] = '-';  // deletion in A
    else consensus_a[tick] = seqA[current_i - 1];      // match/mismatch in A

    if (next_j == current_j) consensus_b[tick] = '-';  // deletion in B
    else consensus_b[tick] = seqB[current_j - 1];      // match/mismatch in B

    current_i = next_i;
    current_j = next_j;
    next_i = I_i(current_i, current_j);
    next_j = I_j(current_i, current_j);
    tick++;
  }
  pos_max = current_j;

#ifdef DEBUG
  //print the consensus sequences
  cout << endl << " " << endl;
  cout << "Alignment:" << endl << endl;
  for (int i = 0; i < lengthSeqA; i++) { cout << seqA[i]; };
  cout << "  and" << endl;
  for (int i = 0; i < lengthSeqB; i++) { cout << seqB[i]; };
  cout << endl << endl;
  for (int i = tick - 1; i >= 0; i--) cout << consensus_a[i];
  cout << endl;
  for (int j = tick - 1; j >= 0; j--) cout << consensus_b[j];
  cout << endl;
#endif

}

double SWAligner::similarityScore(char a, char b) {
  double result;
  if (a == b) {
    result = 1;
  } else {
    result = penalty;
  }
  return result;
}

double SWAligner::getScore() const {
  return matrix_max;
}

int SWAligner::getPos() const {
  return pos_max;
}
