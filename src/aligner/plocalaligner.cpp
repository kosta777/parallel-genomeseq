#include "plocalaligner.h"
#include "smithwaterman.h"
#include <string_view>
#include <memory>
#include <vector>
#include <assert.h>

template<class SMT, class LAT>
OMPParallelLocalAligner<SMT, LAT>::OMPParallelLocalAligner(std::string_view first_sequence, std::string_view second_sequence, int npiece, double overlap_ratio)
    :
    OMPParallelLocalAligner(first_sequence, second_sequence, npiece, overlap_ratio, [](const char &a, const char &b) { return a == b ? 3.0 : -3.0; }, 2.0) {}

template<class SMT, class LAT>
OMPParallelLocalAligner<SMT, LAT>::OMPParallelLocalAligner(std::string_view first_sequence, std::string_view second_sequence, int npiece, double overlap_ratio, double gap_penalty) :
    OMPParallelLocalAligner(first_sequence, second_sequence, npiece, overlap_ratio, [](const char &a, const char &b) { return a == b ? 3.0 : -3.0; }, gap_penalty) {}

template<class SMT, class LAT>
OMPParallelLocalAligner<SMT, LAT>::OMPParallelLocalAligner(std::string_view first_sequence, std::string_view second_sequence, int npiece, double overlap_ratio,
                                                      std::function<double(const char &, const char &)> &&scoring_function) :
    OMPParallelLocalAligner(first_sequence, second_sequence, npiece, overlap_ratio, std::move(scoring_function), 2.0) {}

template<class SMT, class LAT>
OMPParallelLocalAligner<SMT, LAT>::OMPParallelLocalAligner(std::string_view first_sequence, std::string_view second_sequence, int npiece, double overlap_ratio,
                                                      std::function<double(const char &, const char &)> &&scoring_function,
                                                      double gap_penalty) :
    sm_timings(),
    pos(0),
    max_score(-1),
    gap_penalty(gap_penalty),
    overlap_ratio(overlap_ratio),
    npiece(npiece),
    consensus_x(),
    consensus_y(),
    sequence_x(first_sequence),
    sequence_y(second_sequence),
    scoring_function(std::move(scoring_function)) {
  consensus_x.reserve(sequence_x.size());
  consensus_y.reserve(sequence_x.size());
}

std::vector<std::pair<size_t, size_t>> _make_string_range(int npiece, size_t shortstringlength ,size_t longstringlength, double overlap_ratio) {
  auto overlaplength = (size_t)(shortstringlength * overlap_ratio);
  std::vector<std::pair<size_t, size_t>> string_ranges;
  auto piecelength = (longstringlength + (npiece - 1) * overlaplength) / npiece;
  assert(overlaplength <= piecelength);
  size_t left = 0;
  size_t right = piecelength;
  string_ranges.emplace_back(left, right);
  while ((double)(longstringlength - right + overlaplength) > piecelength ) {
    left = std::max((size_t)0, right - overlaplength);
    right = std::min(left + piecelength, longstringlength);
    if (longstringlength - right < std::min((size_t)20, piecelength - overlaplength)) right = longstringlength;
    string_ranges.emplace_back(left, right);
  }
  if (right < longstringlength) string_ranges.emplace_back(std::max((size_t)0, right - overlaplength), longstringlength);
  assert(string_ranges.size() == (size_t)npiece);
  return string_ranges;
}

template<class SMT, class LAT>
double OMPParallelLocalAligner<SMT, LAT>::calculateScore() {
  auto shortstringlength = sequence_x.size();
  auto longstringlength = sequence_y.size();
  double max_score_l = -1.0;
  int max_score_piece = -1;
  auto string_ranges = _make_string_range(npiece, shortstringlength, longstringlength, overlap_ratio);
#ifdef USEOMP
  #pragma omp parallel for default(none) shared(string_ranges, shortstringlength, longstringlength, max_score_l, max_score_piece, \
                                                overlap_ratio, sequence_x, sequence_y, scoring_function, gap_penalty)
#endif
  for (int i = 0; i < npiece; i++) {
    auto [left, right] = string_ranges[i];
    std::string_view sequence_y_i = sequence_y.substr(left, right - left);
    {
      auto smp = std::make_unique<SMT>(sequence_x, sequence_y_i);
      smp -> iterate(scoring_function, gap_penalty);
      auto [x, y, max] = smp -> find_index_of_maximum();
      if (max > max_score_l) {
        max_score_l = max;
        max_score_piece = i;
      }
    }
  }
  max_score = max_score_l;
  auto [left, right] = string_ranges[max_score_piece];
  {
    std::string_view sequence_y_i = sequence_y.substr(left, right - left);
    auto la = std::make_unique<LAT>(sequence_x, sequence_y_i);
    la -> calculateScore();
    pos = la -> getPos() + left;
    sm_timings = la -> getTimings();
    consensus_x = la -> getConsensus_x();
    consensus_y = la -> getConsensus_y();
  }
  return max_score;
}

template class OMPParallelLocalAligner<Similarity_Matrix, SWAligner<Similarity_Matrix>>;
template class OMPParallelLocalAligner<Similarity_Matrix_Skewed, SWAligner<Similarity_Matrix>>;
template class OMPParallelLocalAligner<Similarity_Matrix, SWAligner<Similarity_Matrix_Skewed>>;
template class OMPParallelLocalAligner<Similarity_Matrix_Skewed, SWAligner<Similarity_Matrix_Skewed>>;