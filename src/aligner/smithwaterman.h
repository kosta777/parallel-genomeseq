#ifndef _SMITHWATERMAN_LINEAR_H_
#define _SMITHWATERMAN_LINEAR_H_

#include <functional>
#include <Eigen/Dense>
#include <string_view>
#include "localaligner.h"

typedef std::pair<Eigen::Index, Eigen::Index> index_tuple;

template <class Similarity_Matrix_Type>
class SWAligner : public LocalAligner<Similarity_Matrix_Type> {
  public:
    SWAligner(std::string_view, std::string_view);
    SWAligner(std::string_view, std::string_view, std::function<double(const char &, const char &)> &&);
    /*
     * Maximum score in similarity matrix.
     */
    double calculateScore() override;

    double getScore() const override {return max_score;}

    /*
     * Position of the first matched position in sequence_x. Corresponds to POS in SAM file.
     *
     * @return 1-based (as opposed to standard CS 0-based) counter.
     */
    unsigned int getPos() const override {return pos;}
    std::string_view getConsensus_x() const override {return sequence_x;}
    std::string_view getConsensus_y() const override {return sequence_y;}
    const Similarity_Matrix_Type& getSimilarity_matrix() const override { return similarity_matrix;}

  private:
    unsigned int pos;
    double max_score;
    std::string sequence_x;
    std::string sequence_y;
    Similarity_Matrix_Type similarity_matrix;
    void traceback(index_tuple, unsigned int &);
    std::function<double(const char &, const char &)> scoring_function;
};

#endif