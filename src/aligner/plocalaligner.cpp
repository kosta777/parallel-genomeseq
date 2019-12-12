#include "plocalaligner.h"
#include "smithwaterman.h"
#include <string_view>
#include <memory>
#include <vector>
#include <assert.h>
#include <chrono>

template<class SMT, class LAT>
OMPParallelLocalAligner<SMT, LAT>::OMPParallelLocalAligner(std::string_view first_sequence,
                                                           std::string_view second_sequence,
                                                           int npiece,
                                                           float overlap_ratio)
    :
    OMPParallelLocalAligner(first_sequence,
                            second_sequence,
                            npiece,
                            overlap_ratio,
                            [](const char &a, const char &b) { return a == b ? 3.0 : -3.0; },
                            2.0) {}

template<class SMT, class LAT>
OMPParallelLocalAligner<SMT, LAT>::OMPParallelLocalAligner(std::string_view first_sequence,
                                                           std::string_view second_sequence,
                                                           int npiece,
                                                           float overlap_ratio,
                                                           float gap_penalty) :
    OMPParallelLocalAligner(first_sequence,
                            second_sequence,
                            npiece,
                            overlap_ratio,
                            [](const char &a, const char &b) { return a == b ? 3.0 : -3.0; },
                            gap_penalty) {}

template<class SMT, class LAT>
OMPParallelLocalAligner<SMT, LAT>::OMPParallelLocalAligner(std::string_view first_sequence,
                                                           std::string_view second_sequence,
                                                           int npiece,
                                                           float overlap_ratio,
                                                           std::function<float(const char &,
                                                                               const char &)> &&scoring_function) :
    OMPParallelLocalAligner(first_sequence, second_sequence, npiece, overlap_ratio, std::move(scoring_function), 2.0) {}

std::vector<std::pair<Eigen::Index, Eigen::Index>> _make_string_range(int npiece,
                                                                      Eigen::Index shortstringlength,
                                                                      Eigen::Index longstringlength,
                                                                      float overlap_ratio) {
  auto overlaplength = (Eigen::Index) (shortstringlength * overlap_ratio);
  std::vector<std::pair<Eigen::Index, Eigen::Index>> string_ranges;
  if (npiece == 1) {string_ranges.emplace_back(0, longstringlength); return string_ranges;}
  auto piecelength = (longstringlength + (npiece - 1) * overlaplength) / npiece;
  assert(overlaplength <= piecelength);
  Eigen::Index left = 0;
  Eigen::Index right = piecelength;
  string_ranges.emplace_back(left, right);
  int n_ranges = 1;
  while (n_ranges < npiece - 1) {
    left = std::max((Eigen::Index) 0, right - overlaplength);
    right = std::min(left + piecelength, longstringlength);
    string_ranges.emplace_back(left, right);
    n_ranges += 1;
  }
  assert(right < longstringlength);
  string_ranges.emplace_back(std::max((Eigen::Index) 0, right - overlaplength), longstringlength);
  assert(string_ranges.size() == (size_t) npiece);
  return string_ranges;
}

template<class SMT, class LAT>
OMPParallelLocalAligner<SMT, LAT>::OMPParallelLocalAligner(std::string_view first_sequence,
                                                           std::string_view second_sequence,
                                                           int npiece,
                                                           float overlap_ratio,
                                                           std::function<float(const char &,
                                                                               const char &)> &&scoring_function,
                                                           float gap_penalty) :
    sm_timings(2),
    pos(0),
    max_score(-1),
    gap_penalty(gap_penalty),
    overlap_ratio(overlap_ratio),
    npiece(npiece),
    consensus_x(),
    consensus_y(),
    sequence_x(first_sequence),
    sequence_y(second_sequence),
    scoring_function(std::move(scoring_function)),
    string_ranges(_make_string_range(npiece, sequence_x.size(), sequence_y.size(), overlap_ratio)),
    smptr_vec(){
  consensus_x.reserve(sequence_x.size());
  consensus_y.reserve(sequence_x.size());
#ifdef USEOMP
#pragma omp parallel for default(none) shared(npiece,sequence_x, sequence_y, smptr_vec) num_threads(npiece) schedule(static, 1)
#endif
  for (int i = 0; i < npiece; i++) {
    Eigen::Index left = string_ranges[i].first;
    Eigen::Index right = string_ranges[i].second;
    std::string_view sequence_y_i = sequence_y.substr(left, right - left);
    auto smptr = std::make_unique<SMT>(sequence_x, sequence_y_i);
#pragma omp critical
    smptr_vec.emplace_back(std::move(smptr));
  }
}

template<class SMT, class LAT>
float OMPParallelLocalAligner<SMT, LAT>::calculateScore() {
  float max_score_l = -1.0;
  int max_score_piece = 0;
    auto iter_ad_read_start = std::chrono::high_resolution_clock::now();
#ifdef USEOMP
#pragma omp parallel for default(none) shared(npiece, max_score_l, max_score_piece, scoring_function, gap_penalty, smptr_vec) num_threads(npiece) schedule(static,1)
#endif
    for (int i = 0; i < npiece; i++) {
      smptr_vec[i]->iterate(scoring_function, gap_penalty);
    }
    auto iter_ad_read_end = std::chrono::high_resolution_clock::now();//OVERESTIMATE TIME
    sm_timings[0] = (float) std::chrono::duration_cast<std::chrono::microseconds>(iter_ad_read_end - iter_ad_read_start).count();
    double time_iter = 0.0;
#ifdef USEOMP
#pragma omp parallel for default(none) shared(npiece, max_score_l, max_score_piece, smptr_vec, time_iter) num_threads(npiece) schedule(static, 1)
#endif
    for (int i = 0; i < npiece; i++) {
      auto[x, y, max] = smptr_vec[i]->find_index_of_maximum();
      time_iter += smptr_vec[i]->getTimings()[0];
      if (max > max_score_l) {
        max_score_l = max;
        max_score_piece = i;
      }
    }
  sm_timings[1] = time_iter;

  auto[left, right] = string_ranges[max_score_piece];
  {
    std::string_view sequence_y_i = sequence_y.substr(left, right - left);
    auto la = std::make_unique<LAT>(sequence_x, sequence_y_i);
    la->calculateScore();
    pos = la->getPos() + left;
    max_score = la->getScore();
    consensus_x = la->getConsensus_x();
    consensus_y = la->getConsensus_y();
  }
  return max_score;
}

template
class OMPParallelLocalAligner<Similarity_Matrix, SWAligner<Similarity_Matrix>>;
template
class OMPParallelLocalAligner<Similarity_Matrix_Skewed, SWAligner<Similarity_Matrix>>;
template
class OMPParallelLocalAligner<Similarity_Matrix, SWAligner<Similarity_Matrix_Skewed>>;
template
class OMPParallelLocalAligner<Similarity_Matrix_Skewed, SWAligner<Similarity_Matrix_Skewed>>;
