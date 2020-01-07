#define EIGEN_NO_DEBUG //SUPRESS EIGEN ASSERTION: causing LLC missing
#define EIGEN_DONT_PARALLELIZE
#include "similaritymatrix.h"
#include <iostream>
#include <algorithm>
#include <Eigen/Dense>
#include <chrono>
#include <immintrin.h>

#ifdef USEOMP
#include <omp.h>
#endif

#define N_PACK 32

Similarity_Matrix::Similarity_Matrix(std::string_view sequence_x, std::string_view sequence_y) :
    raw_matrix(sequence_x.size() + 1, sequence_y.size() + 1), sequence_x(sequence_x), sequence_y(sequence_y) {
  raw_matrix.setZero();
}

std::tuple<Eigen::Index, Eigen::Index, float> Similarity_Matrix::find_index_of_maximum() const {
  Eigen::Index x = 0, y = 0;
  auto max = raw_matrix.maxCoeff(&x, &y);
#ifdef VERBOSE
  std::cout << "Maximum is " << max << " @ (" << x << ", " << y << ")" << std::endl;
#endif
  return {x, y, max};
}

Eigen::VectorXf Similarity_Matrix::getTimings() const {
  Eigen::VectorXf sm_timings_tmp = Eigen::VectorXf::Zero(2);
  sm_timings_tmp(0) = sm_iter_ad_read_time;
  sm_timings_tmp(1) = sm_iter_ad_i_times.sum();
  return sm_timings_tmp;
}

void Similarity_Matrix::print_matrix() const {
  std::cout << raw_matrix << std::endl;
}

inline float _fast_max(float x, float y) {
  bool flag = x > y;
  return flag * x + (!flag) * y;
}

#ifdef USEOMP
#pragma omp declare simd uniform(gap_penalty)
#endif
float dp_func(float north, float west, float north_west, float score, float gap_penalty) {
  auto x = north_west + score;
  auto y = west - gap_penalty;
  auto z = north - gap_penalty;
  return std::max(std::max(x, y), std::max(z, (float) 0.0));
}

inline Eigen::Array4f dp_func(Eigen::Array4f north,
                              Eigen::Array4f west,
                              Eigen::Array4f north_west,
                              Eigen::Array4f score,
                              float gap_penalty) {
  auto x = north_west + score;
  auto y = west - gap_penalty;
  auto z = north - gap_penalty;
  return x.max(y).max(z.max((float) 0.0));
}

inline __m256 dp_func(__m256 north, __m256 west, __m256 north_west, __m256 score, __m256 gap_penalty_vec) {
  auto zeros = _mm256_set1_ps(0.0);
  auto x = north_west + score;
  auto y = west - gap_penalty_vec;
  auto z = north - gap_penalty_vec;
  return _mm256_max_ps(_mm256_max_ps(x, y), _mm256_max_ps(z, zeros));
}

inline __m256i dp_func(__m256i north, __m256i west, __m256i north_west, __m256i score_p, __m256i score_m, __m256i gap_penalty_vec) {
  auto x = _mm256_adds_epu8(north_west, score_p);
  x = _mm256_subs_epu8(x, score_m);
  auto y = _mm256_subs_epu8(west, gap_penalty_vec);
  auto z = _mm256_subs_epu8(north, gap_penalty_vec);
  return _mm256_max_epu8(_mm256_max_epu8(x, y), z);
}

void Similarity_Matrix::iterate_column(const std::function<float(const char &, const char &)> &scoring_function,
                                       float gap_penalty) {
  auto nrows = raw_matrix.rows();
  auto ncols = raw_matrix.cols();
  for (Eigen::Index j = 1; j < ncols; j++) {
    for (Eigen::Index i = 1; i < nrows; i++) {
      auto west = raw_matrix(i, j - 1);
      auto north = raw_matrix(i - 1, j);
      auto north_west = raw_matrix(i - 1, j - 1);
      auto a = sequence_x[i - 1];
      auto b = sequence_y[j - 1];
      raw_matrix(i, j) = dp_func(north, west, north_west, scoring_function(a, b), gap_penalty);
    }
  }
}

void Similarity_Matrix::iterate(const std::function<float(const char &, const char &)> &scoring_function,
                                float gap_penalty) {
  const unsigned int dim_x = raw_matrix.rows();
  const unsigned int dim_y = raw_matrix.cols();

  auto iter_ad_read_start = std::chrono::high_resolution_clock::now();
  for (Eigen::Index i = 1; i < dim_x + dim_y - 2; ++i) {
    Eigen::Index local_i = i;
    Eigen::Index starting_k = 1;
    Eigen::Index ending_k = i;
    if (local_i > dim_x - 1) {
      local_i = dim_x - 1;
      starting_k = i - local_i + 1;
      ending_k = starting_k + dim_x - 2;
    }
    if (ending_k > dim_y - 1) {
      ending_k = dim_y - 1;
    }

#ifdef USEOMP
    int omp_n_threads2;
    omp_n_threads2 = omp_get_num_threads();
    auto iter_ad_i_start = std::chrono::high_resolution_clock::now();
    if (sm_finegrain_type==0){   //finegrain type 0: first attempt
      auto ad_len = ending_k-starting_k+1; //anti diagonal length
      Eigen::VectorXf k_vec = Eigen::VectorXf::LinSpaced(ad_len,starting_k,ending_k);
      Eigen::VectorXf local_i_vec = Eigen::VectorXf::LinSpaced(ad_len,local_i,local_i-ad_len+1);
      Eigen::Index ad_idx;
#pragma omp parallel for default(none) shared(ad_len, k_vec, local_i_vec, gap_penalty, scoring_function) private(ad_idx)
      for (ad_idx=0; ad_idx<ad_len; ++ad_idx) {
        index_tuple idx(local_i_vec(ad_idx), k_vec(ad_idx));
        auto west = raw_matrix(idx.first, idx.second - 1);
        auto north = raw_matrix(idx.first - 1, idx.second);
        auto north_west = raw_matrix(idx.first - 1, idx.second - 1);
        auto a = sequence_x[idx.first - 1];
        auto b = sequence_y[idx.second - 1];
        raw_matrix( local_i_vec(ad_idx), k_vec(ad_idx) ) = dp_func(north, west, north_west, scoring_function(a, b), gap_penalty);
//        omp_n_threads = omp_get_num_threads();
      }
    }
    else if (sm_finegrain_type==1){   //finegrain type 1: create chunks manually
      auto ad_len = ending_k-starting_k+1; //anti diagonal length
      Eigen::VectorXf k_vec = Eigen::VectorXf::LinSpaced(ad_len,starting_k,ending_k);
      Eigen::VectorXf local_i_vec = Eigen::VectorXf::LinSpaced(ad_len,local_i,local_i-ad_len+1);

      auto chunk_len = ad_len/sm_nthreads; //length per chunk, nchunks=nthreads
      Eigen::VectorXf chunk_start_vec = Eigen::VectorXf::LinSpaced(sm_nthreads, 0, (sm_nthreads-1)*chunk_len);
      Eigen::VectorXf chunk_end_vec = Eigen::VectorXf::LinSpaced(sm_nthreads, chunk_len, sm_nthreads*chunk_len);
      chunk_end_vec(sm_nthreads-1) = ad_len;

#pragma omp parallel for default(none) shared(ad_len, k_vec, local_i_vec, gap_penalty, scoring_function, chunk_start_vec, chunk_end_vec)
      for (Eigen::Index chunk_idx=0; chunk_idx<sm_nthreads; chunk_idx++){
        for (Eigen::Index ad_idx=chunk_start_vec(chunk_idx); ad_idx<chunk_end_vec(chunk_idx); ++ad_idx) {
          index_tuple idx(local_i_vec(ad_idx), k_vec(ad_idx));
          auto west = raw_matrix(idx.first, idx.second - 1);
          auto north = raw_matrix(idx.first - 1, idx.second);
          auto north_west = raw_matrix(idx.first - 1, idx.second - 1);
          auto a = sequence_x[idx.first - 1];
          auto b = sequence_y[idx.second - 1];
          raw_matrix( local_i_vec(ad_idx), k_vec(ad_idx) ) = dp_func(north, west, north_west, scoring_function(a, b), gap_penalty);
//          omp_n_threads = omp_get_num_threads();
        }
      }
    }
    else if (sm_finegrain_type==2){   //finegrain type 2: schedule static
      auto ad_len = ending_k-starting_k+1; //anti diagonal length
      Eigen::VectorXf k_vec = Eigen::VectorXf::LinSpaced(ad_len,starting_k,ending_k);
      Eigen::VectorXf local_i_vec = Eigen::VectorXf::LinSpaced(ad_len,local_i,local_i-ad_len+1);
      Eigen::Index ad_idx;
#pragma omp parallel for default(none) schedule(static) shared(ad_len, k_vec, local_i_vec, gap_penalty, scoring_function) private(ad_idx)
      for (ad_idx=0; ad_idx<ad_len; ++ad_idx) {
        index_tuple idx(local_i_vec(ad_idx), k_vec(ad_idx));
        auto west = raw_matrix(idx.first, idx.second - 1);
        auto north = raw_matrix(idx.first - 1, idx.second);
        auto north_west = raw_matrix(idx.first - 1, idx.second - 1);
        auto a = sequence_x[idx.first - 1];
        auto b = sequence_y[idx.second - 1];
        raw_matrix( local_i_vec(ad_idx), k_vec(ad_idx) ) = dp_func(north, west, north_west, scoring_function(a, b), gap_penalty);
//        omp_n_threads = omp_get_num_threads();
      }
    }
    else if (sm_finegrain_type==3){   //finegrain type 1: create chunks manually, but with scheule static, with array
      int L1_cachesize = 64;
      auto ad_len = ending_k-starting_k+1; //anti diagonal length
      auto ad_len_padded = (ad_len/L1_cachesize+1)*L1_cachesize;
      std::cout<<"ad_len_padded: "<<ad_len_padded<<std::endl;
      Eigen::VectorXf k_vec = Eigen::VectorXf::LinSpaced(ad_len,starting_k,ending_k);
      Eigen::VectorXf local_i_vec = Eigen::VectorXf::LinSpaced(ad_len,local_i,local_i-ad_len+1);

      auto chunk_len = ad_len/sm_nthreads; //length per chunk, nchunks=nthreads
      Eigen::VectorXf chunk_start_vec = Eigen::VectorXf::LinSpaced(sm_nthreads, 0, (sm_nthreads-1)*chunk_len);
      Eigen::VectorXf chunk_end_vec = Eigen::VectorXf::LinSpaced(sm_nthreads, chunk_len, sm_nthreads*chunk_len);
      chunk_end_vec(sm_nthreads-1) = ad_len;

      Eigen::VectorXf ad_vec_tmp = Eigen::VectorXf::Zero(ad_len);

#pragma omp parallel for default(none) schedule(static) shared(ad_vec_tmp, ad_len, k_vec, local_i_vec, gap_penalty, scoring_function, chunk_start_vec, chunk_end_vec)
      for (Eigen::Index chunk_idx=0; chunk_idx<sm_nthreads; chunk_idx++){
        for (Eigen::Index ad_idx=chunk_start_vec(chunk_idx); ad_idx<chunk_end_vec(chunk_idx); ++ad_idx) {
          index_tuple idx(local_i_vec(ad_idx), k_vec(ad_idx));
          auto west = raw_matrix(idx.first, idx.second - 1);
          auto north = raw_matrix(idx.first - 1, idx.second);
          auto north_west = raw_matrix(idx.first - 1, idx.second - 1);
          auto a = sequence_x[idx.first - 1];
          auto b = sequence_y[idx.second - 1];
          ad_vec_tmp(ad_idx) = dp_func(north, west, north_west, scoring_function(a, b), gap_penalty);
        }
      }
      for (Eigen::Index ad_idx=0; ad_idx<ad_len; ad_idx++){
        raw_matrix( local_i_vec(ad_idx), k_vec(ad_idx) ) = ad_vec_tmp(ad_idx);
//      omp_n_threads = omp_get_num_threads();
      }
    }
    else if (sm_finegrain_type==4){   //finegrain type 1: create chunks manually, but with scheule static, with array
      int L1_cachesize = 64;   //getconf LEVEL1_DCACHE_LINESIZE
      //#define CACHE_LINE_SIZE sysconf(_SC_LEVEL1_DCACHE_LINESIZE)
      auto ad_len = ending_k-starting_k+1; //anti diagonal length
      auto ad_len_padded = (ad_len/(L1_cachesize*sm_nthreads)+1)*L1_cachesize*sm_nthreads;
      if (i%1000==0){
        std::cout<<"ad_len_padded: "<<ad_len_padded<<std::endl;
      }
      Eigen::VectorXf k_vec = Eigen::VectorXf::LinSpaced(ad_len,starting_k,ending_k);
      Eigen::VectorXf local_i_vec = Eigen::VectorXf::LinSpaced(ad_len,local_i,local_i-ad_len+1);

      Eigen::VectorXf ad_vec_tmp = Eigen::VectorXf::Zero(ad_len_padded);

#pragma omp parallel for default(none) schedule(static,ad_len_padded/sm_nthreads) shared(ad_vec_tmp, ad_len, ad_len_padded,k_vec, local_i_vec, gap_penalty, scoring_function)
      for (Eigen::Index ad_idx=0; ad_idx<ad_len_padded; ++ad_idx) {
    if (ad_idx<ad_len){
          index_tuple idx(local_i_vec(ad_idx), k_vec(ad_idx));
          auto west = raw_matrix(idx.first, idx.second - 1);
          auto north = raw_matrix(idx.first - 1, idx.second);
          auto north_west = raw_matrix(idx.first - 1, idx.second - 1);
          auto a = sequence_x[idx.first - 1];
          auto b = sequence_y[idx.second - 1];
          ad_vec_tmp(ad_idx) = dp_func(north, west, north_west, scoring_function(a, b), gap_penalty);
    }
      }
      for (Eigen::Index ad_idx=0; ad_idx<ad_len; ad_idx++){
        raw_matrix( local_i_vec(ad_idx), k_vec(ad_idx) ) = ad_vec_tmp(ad_idx);
      }
    }
    auto iter_ad_i_end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(iter_ad_i_end-iter_ad_i_start);
    sm_iter_ad_i_times.resize(i);
    sm_iter_ad_i_times(i-1) = (float) duration.count();

#else
    for (Eigen::Index k = starting_k; k <= ending_k; ++k) {
      index_tuple idx(local_i, k);
      auto west = raw_matrix(idx.first, idx.second - 1);
      auto north = raw_matrix(idx.first - 1, idx.second);
      auto north_west = raw_matrix(idx.first - 1, idx.second - 1);
      auto a = sequence_x[idx.first - 1];
      auto b = sequence_y[idx.second - 1];
      raw_matrix(local_i, k) = dp_func(north, west, north_west, scoring_function(a, b), gap_penalty);
      --local_i;
    }
#endif
  }
  auto iter_ad_read_end = std::chrono::high_resolution_clock::now();
  auto read_duration = std::chrono::duration_cast<std::chrono::microseconds>(iter_ad_read_end - iter_ad_read_start);
  sm_iter_ad_read_time = (float) read_duration.count();
//  std::cout<<"Similarity_Matrix::iterate omp_n_threads: "<<omp_n_threads<<std::endl;
//  std::cout<<"Similarity_Matrix::iterate omp_n_threads2: "<<omp_n_threads2<<std::endl;
}

const Eigen::MatrixXf &Similarity_Matrix::get_matrix() const {
  return raw_matrix;
}

inline Eigen::Index _trueindex2invindex(Eigen::Index ti, Eigen::Index size) {// 1-INDEX based for both input and output
  return size + 1 - ti;
}

Similarity_Matrix_Skewed::Similarity_Matrix_Skewed(std::string_view _sequence_x, std::string_view _sequence_y)
    : sequence_x(_sequence_y.size() + N_PACK),
      inv_sequence_y(_sequence_x.size() + N_PACK),
      sm_iter_ad_read_time(-1.0),
      sm_iter_ad_i_times() {//SWITCHING SEQX AND SEQY
  len_x = _sequence_y.size() + 1;//SWITCHING SEQX AND SEQY
  len_y = _sequence_x.size() + 1;//SWITCHING SEQX AND SEQY
  nrows = std::min(len_x, len_y);
  ncols = std::max(len_x, len_y);
  for (Eigen::Index i = 0; i < len_x - 1; i++) sequence_x(i) = (uint8_t) _sequence_y[i];
  for (Eigen::Index i = 0; i < len_y - 1; i++)
    inv_sequence_y(_trueindex2invindex(i + 1, len_y - 1) - 1) = (uint8_t) _sequence_x[i];
  //std::reverse_copy(sequence_x.begin(), sequence_x.end(), inv_sequence_y.begin());//SWITCHING SEQX AND SEQY
  raw_matrix = MatrixX8u(nrows + N_PACK, ncols);//Avoid blendv/store failure in the end
  raw_matrix.setZero();
}

std::tuple<Eigen::Index, Eigen::Index, float> Similarity_Matrix_Skewed::find_index_of_maximum() const {
  index_tuple maxidx(0, 0);
  auto max = raw_matrix.maxCoeff(&maxidx.first, &maxidx.second);
  auto[x, y] = rawindex2trueindex(maxidx);
#ifdef VERBOSE
  std::cout << "Maximum is " << max << " @( " << maxidx.first << "," << maxidx.second << ")" << std::endl;
#endif
  return {y, x, max}; //SWITCHING SEQX AND SEQY
}

void Similarity_Matrix_Skewed::print_matrix() const {
  Eigen::MatrixXf similarity_matrix(len_y, len_x);
  for (int i = 0; i < len_x; i++) {
    for (int j = 0; j < len_y; j++){
      index_tuple trueidx(i, j);
      auto[ri, rj] = trueindex2rawindex(trueidx);
      similarity_matrix(j, i) = raw_matrix(ri, rj);
    }
  }
  std::cout<< "\n" << similarity_matrix << std::endl;
}

void Similarity_Matrix_Skewed::print_matrix_raw() const {
  Eigen::MatrixXf similarity_matrix(nrows + N_PACK, ncols);
  for (int j = 0; j < ncols; j++) {
    for (int i = 0; i < nrows + N_PACK; i++) {
      similarity_matrix(i, j) = raw_matrix(i, j);
    }
  }
  std::cout << "\n" <<similarity_matrix << std::endl;
}

Eigen::VectorXf Similarity_Matrix_Skewed::getTimings() const {
  Eigen::VectorXf sm_timings_tmp = Eigen::VectorXf::Zero(2);
  sm_timings_tmp(0) = sm_iter_ad_read_time;
  sm_timings_tmp(1) = sm_iter_ad_i_times.sum();
  return sm_timings_tmp;
}

index_tuple _rawindex2trueindex(size_t ri, size_t rj, size_t nrows, size_t ncols, size_t len_x, size_t len_y) {
  if (rj < nrows - 1) {
    if (ri <= rj) {//upper triangular part
      return index_tuple(ri, rj - ri);
    } else {//lower triangular part
      return index_tuple(len_x - nrows + ri,
                         len_y - ri
                             + rj);//index_tuple(len_x - (nrows - 1 - rj) + (ri - rj - 1), len_y - 1 - (ri - rj - 1));
    }
  } else {//equal-length diagonal part
    if (len_x <= len_y) {//diagonal propagates horizontally (+y)
      return index_tuple(ri, rj - ri);
    } else {//diagonal propagates vertically (+x)
      return index_tuple(rj - (nrows - 1) + ri, nrows - 1 - ri);
    }
  }
}

index_tuple Similarity_Matrix_Skewed::rawindex2trueindex(index_tuple raw_index) const {
  auto[ri, rj] = raw_index;
  return _rawindex2trueindex(ri, rj, nrows, ncols, len_x, len_y);//Always have nrows <= ncols
}

index_tuple _trueindex2rawindex(size_t ti, size_t tj, size_t nrows, size_t ncols, size_t len_x, size_t len_y) {
  if (ti + tj < nrows - 1) {//upper triangular part
    return index_tuple(ti, ti + tj);
  } else if (ti + tj > ncols - 1) {//lower triangular part
    return index_tuple(ti - ncols + len_y,
                       ti + tj - (ncols - 1)
                           - 1);//index_tuple(ti + tj - (ncols - 1) + len_y - 1 - tj , ti + tj - (ncols - 1) - 1);
  } else {//equal-length diagonal part
    auto delta_x = len_x <= len_y ? ti : len_y - 1 - tj;
    return index_tuple(delta_x, ti + tj);
  }
}

index_tuple Similarity_Matrix_Skewed::trueindex2rawindex(index_tuple true_index) const {
  auto[ti, tj] = true_index;
  return _trueindex2rawindex(ti, tj, nrows, ncols, len_x, len_y);//Always have nrows <= ncols
}

void _print_string(uint8_t* strp, int size) {
  std::string_view view((char*)strp, size);
  std::cout << view << std::endl;
}

uint8_t _saturate(float a) {
  if (a < 0) {
    return 0;
  } else if (a > 255) {
    return 255;
  } else {
    return (uint8_t) a;
  }
}

void Similarity_Matrix_Skewed::iterate(const std::function<float(const char &, const char &)> &scoring_function,
                                       float gap_penalty) {
  auto flag = len_x < len_y;
  auto match_score = _mm256_set1_epi8(_saturate(scoring_function('A', 'A')));
  auto mismatch_score = _mm256_set1_epi8(_saturate(-scoring_function('A', 'T'))); //ASSUMING negative mismatch score
  auto zero_score = _mm256_set1_epi8(0);
  auto gap_penalty_vec = _mm256_set1_epi8(_saturate(gap_penalty));
  ArrayX8i mask_arr_init{N_PACK}, mask_arr{N_PACK};
  for (int i = 0; i < N_PACK; i++) mask_arr_init[i] = (int8_t) (i - 1);
  auto iter_ad_read_start = std::chrono::high_resolution_clock::now();
#ifdef VERBOSE
#ifdef USEOMP
  int omp_n_threads = omp_get_num_threads();
  std::cout<<"Similarity_Matrix_Skewed::iterate line 205 omp_n_threads: "<<omp_n_threads<<std::endl;
#endif
#endif
  //Phase 1: Upper triangular part
  for (Eigen::Index j = 2; j < nrows; j++) {
    auto i0 = 1;
    auto in = j - 1;
    auto[ti1, tj1] = _rawindex2trueindex(i0, j, nrows, ncols, len_x, len_y);//here tj1 > tj2 while ti1 < ti2
    auto tj1_inv = _trueindex2invindex(tj1, len_y - 1);//now tj1_inv < tj2_inv
    Eigen::Index i;
    for (i = i0; i <= in - (N_PACK - 1); i += N_PACK) {
      __m256i north = _mm256_loadu_si256((__m256i *)&raw_matrix(i - 1, j - 1));
      __m256i west = _mm256_loadu_si256((__m256i *)&raw_matrix(i, j - 1));
      __m256i north_west = _mm256_loadu_si256((__m256i *)&raw_matrix(i - 1, j - 2));
      __m256i seqx = _mm256_loadu_si256((__m256i *) &sequence_x(ti1 - 1 + i - i0));
      __m256i seqy = _mm256_loadu_si256((__m256i *) &inv_sequence_y(tj1_inv - 1 + i - i0));
      __m256i seqmask = _mm256_cmpeq_epi8(seqx, seqy);
      __m256i scores_p = _mm256_blendv_epi8(zero_score, match_score, seqmask);
      __m256i scores_m = _mm256_blendv_epi8(mismatch_score, zero_score, seqmask);
      __m256i dest = dp_func(north, west, north_west, scores_p, scores_m, gap_penalty_vec);
      _mm256_storeu_si256((__m256i *)&raw_matrix(i, j), dest);
    }
    if (i <= in) {
      auto di = (int8_t) (i - in);
      mask_arr = mask_arr_init + di;
      __m256i mask = _mm256_loadu_si256((__m256i *) mask_arr.data());
      __m256i north = _mm256_loadu_si256((__m256i *)&raw_matrix(i - 1, j - 1));
      __m256i west = _mm256_loadu_si256((__m256i *)&raw_matrix(i, j - 1));
      __m256i north_west = _mm256_loadu_si256((__m256i *)&raw_matrix(i - 1, j - 2));
      __m256i seqx = _mm256_loadu_si256((__m256i *) &sequence_x(ti1 - 1 + i - i0));
      __m256i seqy = _mm256_loadu_si256((__m256i *) &inv_sequence_y(tj1_inv - 1 + i - i0));
      __m256i seqmask = _mm256_cmpeq_epi8(seqx, seqy);
      __m256i scores_p = _mm256_blendv_epi8(zero_score, match_score, seqmask);
      __m256i scores_m = _mm256_blendv_epi8(mismatch_score, zero_score, seqmask);
      __m256i dest_ori = _mm256_loadu_si256((__m256i *)&raw_matrix(i, j));
      __m256i dest = _mm256_blendv_epi8(dest_ori, dp_func(north, west, north_west, scores_p, scores_m, gap_penalty_vec), mask);
      _mm256_storeu_si256((__m256i *)&raw_matrix(i, j), dest);
    }
    //_print_string(&sequence_x(ti1 - 1 + i - i0), in - i + 1);
    //_print_string(&inv_sequence_y(tj1_inv - 1 + i - i0), in - i + 1);
  }
  //Phase 2: Equal-length diagonal part
  for (Eigen::Index j = nrows; j < ncols; j++) {
    int i0 = (int) flag;
    auto in = nrows - 1 + i0 - 1;
    auto[ti1, tj1] = _rawindex2trueindex(i0, j, nrows, ncols, len_x, len_y);//here tj1 > tj2 while ti1 < ti2
    auto tj1_inv = _trueindex2invindex(tj1, len_y - 1);//now tj1_inv < tj2_inv
    if (flag) {//Condition 1: diagonal propagate horizontaly (+y)
      Eigen::Index i;
#ifdef MTSIMD
//      std::cout<<"sm_mt_simd activate: h propa "<<std::endl;
      sm_mt_simd = 666;
      #pragma omp parallel for default(none) private(i0,in,j,ti1,tj1_inv,match_score,zero_score,mismatch_score,gap_penalty_vec) schedule(static)
#endif
      for (i = i0; i <= in - (N_PACK - 1); i += N_PACK) {
        __m256i north = _mm256_loadu_si256((__m256i *)&raw_matrix(i - 1, j - 1));
        __m256i west = _mm256_loadu_si256((__m256i *)&raw_matrix(i, j - 1));
        __m256i north_west = _mm256_loadu_si256((__m256i *)&raw_matrix(i - 1, j - 2));
        __m256i seqx = _mm256_loadu_si256((__m256i *) &sequence_x(ti1 - 1 + i - i0));
        __m256i seqy = _mm256_loadu_si256((__m256i *) &inv_sequence_y(tj1_inv - 1 + i - i0));
        __m256i seqmask = _mm256_cmpeq_epi8(seqx, seqy);
        __m256i scores_p = _mm256_blendv_epi8(zero_score, match_score, seqmask);
      __m256i scores_m = _mm256_blendv_epi8(mismatch_score, zero_score, seqmask);
        __m256i dest = dp_func(north, west, north_west, scores_p, scores_m, gap_penalty_vec);
        _mm256_storeu_si256((__m256i *)&raw_matrix(i, j), dest);
      }
      if (i <= in) {
        auto di = (int8_t) (i - in);
        mask_arr = mask_arr_init + di;
        __m256i mask = _mm256_loadu_si256((__m256i *) mask_arr.data());
        __m256i north = _mm256_loadu_si256((__m256i *)&raw_matrix(i - 1, j - 1));//picked whose HIGHEST BIT IS ZERO
        __m256i west = _mm256_loadu_si256((__m256i *)&raw_matrix(i, j - 1));
        __m256i north_west = _mm256_loadu_si256((__m256i *)&raw_matrix(i - 1, j - 2));
        __m256i seqx = _mm256_loadu_si256((__m256i *) &sequence_x(ti1 - 1 + i - i0));
        __m256i seqy = _mm256_loadu_si256((__m256i *) &inv_sequence_y(tj1_inv - 1 + i - i0));
        __m256i seqmask = _mm256_cmpeq_epi8(seqx, seqy);
        __m256i scores_p = _mm256_blendv_epi8(zero_score, match_score, seqmask);
      __m256i scores_m = _mm256_blendv_epi8(mismatch_score, zero_score, seqmask);
        __m256i dest_ori = _mm256_loadu_si256((__m256i *)&raw_matrix(i, j));
        __m256i dest = _mm256_blendv_epi8(dest_ori, dp_func(north, west, north_west, scores_p, scores_m, gap_penalty_vec), mask);
        _mm256_storeu_si256((__m256i *)&raw_matrix(i, j), dest);
      }
    } else {//Condition 2: diagonal propagate vertically (+x)
      auto di_nw = j == nrows ? 0 : 1;
      Eigen::Index i;
#ifdef MTSIMD
//      std::cout<<"sm_mt_simd activate: v propa "<<std::endl;
      sm_mt_simd = 667;
      #pragma omp parallel for default(none) private(i0,in,j,ti1,tj1_inv,match_score,zero_score,mismatch_score,gap_penalty_vec, di_nw) schedule(static)
#endif
      for (i = i0; i <= in - (N_PACK - 1); i += N_PACK) {
        __m256i north = _mm256_loadu_si256((__m256i *)&raw_matrix(i, j - 1));
        __m256i west = _mm256_loadu_si256((__m256i *)&raw_matrix(i + 1, j - 1));
        __m256i north_west = _mm256_loadu_si256((__m256i *)&raw_matrix(i + di_nw, j - 2));
        __m256i seqx = _mm256_loadu_si256((__m256i *) &sequence_x(ti1 - 1 + i - i0));
        __m256i seqy = _mm256_loadu_si256((__m256i *) &inv_sequence_y(tj1_inv - 1 + i - i0));
        __m256i seqmask = _mm256_cmpeq_epi8(seqx, seqy);
        __m256i scores_p = _mm256_blendv_epi8(zero_score, match_score, seqmask);
      __m256i scores_m = _mm256_blendv_epi8(mismatch_score, zero_score, seqmask);
        __m256i dest = dp_func(north, west, north_west, scores_p, scores_m, gap_penalty_vec);
        _mm256_storeu_si256((__m256i *)&raw_matrix(i, j), dest);
      }
      if (i <= in) {
        auto di = (int8_t) (i - in);
        mask_arr = mask_arr_init + di;
        __m256i mask = _mm256_loadu_si256((__m256i *) mask_arr.data());
        __m256i north = _mm256_loadu_si256((__m256i *)&raw_matrix(i, j - 1));//picked whose HIGHEST BIT IS ZERO
        __m256i west = _mm256_loadu_si256((__m256i *)&raw_matrix(i + 1, j - 1));
        __m256i north_west = _mm256_loadu_si256((__m256i *)&raw_matrix(i + di_nw, j - 2));
        __m256i seqx = _mm256_loadu_si256((__m256i *) &sequence_x(ti1 - 1 + i - i0));
        __m256i seqy = _mm256_loadu_si256((__m256i *) &inv_sequence_y(tj1_inv - 1 + i - i0));
        __m256i seqmask = _mm256_cmpeq_epi8(seqx, seqy);
        __m256i scores_p = _mm256_blendv_epi8(zero_score, match_score, seqmask);
      __m256i scores_m = _mm256_blendv_epi8(mismatch_score, zero_score, seqmask);
        __m256i dest_ori = _mm256_loadu_si256((__m256i *)&raw_matrix(i, j));
        __m256i dest = _mm256_blendv_epi8(dest_ori, dp_func(north, west, north_west, scores_p, scores_m, gap_penalty_vec), mask);
        _mm256_storeu_si256((__m256i *)&raw_matrix(i, j), dest);
      }
    }
  }
  //Phase 3: Lower triangular part
  for (Eigen::Index j = 0; j < nrows - 1; j++) {
    auto j_prev = j == 0 ? ncols - 1 : j - 1;
    auto j_prev2 = j <= 1 ? ncols - 2 + j : j - 2;
    auto di_nw = (j == 0) && (!flag) ? 1 : 0;
    auto i0 = j + 1;
    auto in = nrows - 1;
    auto[ti1, tj1] = _rawindex2trueindex(i0, j, nrows, ncols, len_x, len_y);//here tj1 > tj2 while ti1 < ti2
    auto tj1_inv = _trueindex2invindex(tj1, len_y - 1);//now tj1_inv < tj2_inv
    Eigen::Index i;
    for (i = i0; i <= in - (N_PACK - 1); i += N_PACK) {
      __m256i north = _mm256_loadu_si256((__m256i *)&raw_matrix(i - 1, j_prev));
      __m256i west = _mm256_loadu_si256((__m256i *)&raw_matrix(i, j_prev));
      __m256i north_west = _mm256_loadu_si256((__m256i *)&raw_matrix(i - 1 + di_nw, j_prev2));
      __m256i seqx = _mm256_loadu_si256((__m256i *) &sequence_x(ti1 - 1 + i - i0));
      __m256i seqy = _mm256_loadu_si256((__m256i *) &inv_sequence_y(tj1_inv - 1 + i - i0));
      __m256i seqmask = _mm256_cmpeq_epi8(seqx, seqy);
      __m256i scores_p = _mm256_blendv_epi8(zero_score, match_score, seqmask);
      __m256i scores_m = _mm256_blendv_epi8(mismatch_score, zero_score, seqmask);
      __m256i dest = dp_func(north, west, north_west, scores_p, scores_m, gap_penalty_vec);
      _mm256_storeu_si256((__m256i *)&raw_matrix(i, j), dest);
    }
    if (i <= in) {
      auto di = (int8_t) (i - in);
      mask_arr = mask_arr_init + di;
      __m256i mask = _mm256_loadu_si256((__m256i *) mask_arr.data());
      __m256i north = _mm256_loadu_si256((__m256i *)&raw_matrix(i - 1, j_prev));//picked whose HIGHEST BIT IS ZERO
      __m256i west = _mm256_loadu_si256((__m256i *)&raw_matrix(i, j_prev));
      __m256i north_west = _mm256_loadu_si256((__m256i *)&raw_matrix(i - 1 + di_nw, j_prev2));
      __m256i seqx = _mm256_loadu_si256((__m256i *) &sequence_x(ti1 - 1 + i - i0));
      __m256i seqy = _mm256_loadu_si256((__m256i *) &inv_sequence_y(tj1_inv - 1 + i - i0));
      __m256i seqmask = _mm256_cmpeq_epi8(seqx, seqy);
      __m256i scores_p = _mm256_blendv_epi8(zero_score, match_score, seqmask);
      __m256i scores_m = _mm256_blendv_epi8(mismatch_score, zero_score, seqmask);
      __m256i dest_ori = _mm256_loadu_si256((__m256i *)&raw_matrix(i, j));
      __m256i dest = _mm256_blendv_epi8(dest_ori, dp_func(north, west, north_west, scores_p, scores_m, gap_penalty_vec), mask);
      _mm256_storeu_si256((__m256i *)&raw_matrix(i, j), dest);
    }
  }
  auto iter_ad_read_end = std::chrono::high_resolution_clock::now();
  auto read_duration = std::chrono::duration_cast<std::chrono::microseconds>(iter_ad_read_end - iter_ad_read_start);
  sm_iter_ad_read_time = (float) read_duration.count();
}
