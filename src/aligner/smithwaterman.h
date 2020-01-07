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
    SWAligner(std::string_view, std::string_view, float);
    SWAligner(std::string_view, std::string_view, std::function<float(const char &, const char &)> &&);
    SWAligner(std::string_view, std::string_view, std::function<float(const char &, const char &)> &&, float);
    /*
     * Calculate similarity matrix entries and traceback.
     *
     * @return Maximum score in the similarity matrix.
     */
    float calculateScore() override;

    float getScore() const override {return max_score;}

    /*
     * Position of the first matched position in sequence_x. Corresponds to POS in SAM file.
     *
     * @return 1-based (as opposed to standard CS 0-based) counter.
     */
    unsigned int getPos() const override { return pos; }
    std::string_view getConsensus_x() const override { return consensus_x; }
    std::string_view getConsensus_y() const override { return consensus_y; }
    const Similarity_Matrix_Type& getSimilarity_matrix() const override { return similarity_matrix; }
    Eigen::VectorXf getTimings() const override { return sm_timings; }
#ifdef USEOMP
    int sw_nthreads;
    int sw_finegrain_type;
    int sw_mt_simd;
#endif
  private:
    Eigen::VectorXf sm_timings;
    void traceback(index_tuple);
    unsigned int pos;
    float max_score;
    float gap_penalty;
    std::string_view sequence_x;
    std::string_view sequence_y;
    /*
     * Portions of sequences x and y found by the algorithm as best match.
     * In reverse order. Gaps are denoted by a dash symbol "-".
     */
    std::string consensus_x;
    std::string consensus_y;
    Similarity_Matrix_Type similarity_matrix;
    std::function<float(const char &, const char &)> scoring_function;
};

#endif
