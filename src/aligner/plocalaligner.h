#ifndef _PLOCALALIGNER_H_
#define _PLOCALALIGNER_H_
#include "localaligner.h"

template <class Similarity_Matrix_Type, class LocalAligner_Type>
class OMPParallelLocalAligner : public ParallelLocalAligner<Similarity_Matrix_Type, LocalAligner_Type> {
  public:
    OMPParallelLocalAligner(std::string_view, std::string_view, int, double);
    OMPParallelLocalAligner(std::string_view, std::string_view, int, double, double);
    OMPParallelLocalAligner(std::string_view, std::string_view, int, double, std::function<double(const char &, const char &)> &&);
    OMPParallelLocalAligner(std::string_view, std::string_view, int, double, std::function<double(const char &, const char &)> &&, double);
    double calculateScore() override;
    double getScore() const override { return max_score; };
    unsigned int getPos() const override { return pos; };
    std::string_view getConsensus_x() const override { return consensus_x; }
    std::string_view getConsensus_y() const override { return consensus_y; }
  private:
    unsigned int pos;
    double max_score;
    double gap_penalty;
    double overlap_ratio; // heuristics of maximum len(consensus_x)/len(sequence_x)
    int npiece; // number of pieces that sequence_y is broken into
    std::string consensus_x;
    std::string consensus_y;
    std::string sequence_x;
    std::string sequence_y;
    std::function<double(const char &, const char &)> scoring_function;
};

#endif //_PLOCALALIGNER_H_
