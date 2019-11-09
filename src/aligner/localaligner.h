#ifndef _LOCAL_ALIGNER_H_
#define _LOCAL_ALIGNER_H_

#include<string>
#include"similaritymatrix.h"

class LocalAligner {
  public:
    virtual double calculateScore() = 0;
    virtual double getScore() const = 0;
    virtual unsigned int getPos() const = 0;
    virtual std::string_view getConsensus_x() const = 0;
    virtual std::string_view getConsensus_y() const = 0;
    virtual const Similarity_Matrix& getSimilarity_matrix() const =0;
};

#endif
