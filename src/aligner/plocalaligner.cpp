#include "plocalaligner.h"
#include "smithwaterman.h"
#include <string_view>
#include <memory>

template<class SMT, class LAT>
OMPParallelLocalAligner<SMT, LAT>::OMPParallelLocalAligner(std::string_view first_sequence, std::string_view second_sequence)
    :
    OMPParallelLocalAligner(first_sequence, second_sequence, [](const char &a, const char &b) { return a == b ? 3.0 : -3.0; }, 2.0) {}

template<class SMT, class LAT>
OMPParallelLocalAligner<SMT, LAT>::OMPParallelLocalAligner(std::string_view first_sequence, std::string_view second_sequence, double gap_penalty) :
    OMPParallelLocalAligner(first_sequence, second_sequence, [](const char &a, const char &b) { return a == b ? 3.0 : -3.0; }, gap_penalty) {}

template<class SMT, class LAT>
OMPParallelLocalAligner<SMT, LAT>::OMPParallelLocalAligner(std::string_view first_sequence, std::string_view second_sequence,
                                                      std::function<double(const char &, const char &)> &&scoring_function) :
    OMPParallelLocalAligner(first_sequence, second_sequence, std::move(scoring_function), 2.0) {}

template<class SMT, class LAT>
OMPParallelLocalAligner<SMT, LAT>::OMPParallelLocalAligner(std::string_view first_sequence, std::string_view second_sequence,
                                                      std::function<double(const char &, const char &)> &&scoring_function,
                                                      double gap_penalty) :
    pos(0),
    max_score(-1),
    gap_penalty(gap_penalty),
    overlap_ratio(2.0),
    sequence_x(first_sequence),
    sequence_y(second_sequence),
    scoring_function(std::move(scoring_function)) {}

std::pair<size_t, size_t> _make_string_range(int ithread, int nthreads, size_t shortstringlength ,size_t longstringlength, double overlap_ratio) {
  auto overlaplength = (size_t)(shortstringlength * overlap_ratio);
  overlaplength = (overlaplength * 2) > (longstringlength / nthreads) ? longstringlength / (2 * nthreads) : overlaplength;
  auto recurrencelength = (longstringlength - 2 * overlaplength) / nthreads;
  return {ithread * recurrencelength, (ithread + 1) * recurrencelength + overlaplength * 2 };
}

template<class SMT, class LAT>
double OMPParallelLocalAligner<SMT, LAT>::calculateScore() {
#ifdef USEOMP
  auto nthreads = OMP_NUM_THREADS ;
#else
  auto nthreads = 1;
#endif
  auto shortstringlength = sequence_x.size();
  auto longstringlength = sequence_y.size();
  double max_score_l = -1.0;
  int max_score_thread = -1;
#ifdef USEOMP
  #pragma omp parallel for default(none) shared(nthreads, shortstringlength, longstringlength, max_score_l, max_score_thread, \
                                                overlap_ratio, sequence_x, sequence_y, scoring_function, gap_penalty)
#endif
  for (int i = 0; i < nthreads; i++) {
    auto [left, right] = _make_string_range(i, nthreads, shortstringlength, longstringlength, overlap_ratio);
    std::string_view sequence_y_i = sequence_y.substr(left, right - left);
    {
      auto smp = std::make_unique<SMT>(sequence_x, sequence_y_i);
      smp -> iterate(scoring_function, gap_penalty);
      auto [x, y, max] = smp -> find_index_of_maximum();
      if (max > max_score_l) {
        max_score_l = max;
        max_score_thread = i;
      }
    }
  }
  max_score = max_score_l;
  auto [left, right] = _make_string_range(max_score_thread, nthreads, shortstringlength, longstringlength, overlap_ratio);
  {
    std::string_view sequence_y_i = sequence_y.substr(left, right - left);
    auto la = std::make_unique<LAT>(sequence_x, sequence_y_i);
    la -> calculateScore();
    pos = la -> getPos() + left;
  }
  return max_score;
}

template class OMPParallelLocalAligner<Similarity_Matrix, SWAligner<Similarity_Matrix>>;
template class OMPParallelLocalAligner<Similarity_Matrix_Skewed, SWAligner<Similarity_Matrix>>;
template class OMPParallelLocalAligner<Similarity_Matrix, SWAligner<Similarity_Matrix_Skewed>>;
template class OMPParallelLocalAligner<Similarity_Matrix_Skewed, SWAligner<Similarity_Matrix_Skewed>>;