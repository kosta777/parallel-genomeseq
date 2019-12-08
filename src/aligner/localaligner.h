#ifndef _LOCAL_ALIGNER_H_
#define _LOCAL_ALIGNER_H_

#include<string>
#include"similaritymatrix.h"

template <class Similarity_Matrix_Type>
class LocalAligner {
  public:
    virtual float calculateScore() = 0;
    virtual float getScore() const = 0;
    virtual unsigned int getPos() const = 0;
    virtual std::string_view getConsensus_x() const = 0;
    virtual std::string_view getConsensus_y() const = 0;
    virtual const Similarity_Matrix_Type& getSimilarity_matrix() const =0;
    virtual Eigen::VectorXf getTimings() const = 0;
};

template <class Similarity_Matrix_Type, class LocalAligner_Type>
class ParallelLocalAligner {
  public:
    virtual float calculateScore() = 0;
    virtual float getScore() const = 0;
    virtual unsigned int getPos() const = 0;
    virtual std::string_view getConsensus_x() const = 0;
    virtual std::string_view getConsensus_y() const = 0;
    virtual Eigen::VectorXf getTimings() const = 0;
};
#endif
