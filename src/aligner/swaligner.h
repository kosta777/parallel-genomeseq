#ifndef _SW_ALIGNER_H_
#define _SW_ALIGNER_H_

#include <string>
#include "localaligner.h"

class SWAligner : public LocalAligner {
  public:
    void setFirstSequence(std::string seq) override;
    void setSecondSequence(std::string seq) override;
    void calculateScore() override;
    double getScore() const override;
  private:
    int ind;
    double penalty = -4;
    std::string seqA;
    std::string seqB;
    double matrix_max;

    double similarityScore(char a, char b);
    double findMax(const double array[], int length);
};
#endif
