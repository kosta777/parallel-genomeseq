#ifndef _SIMILARITY_MATRIX_H_
#define _SIMILARITY_MATRIX_H_

#include <Eigen/Dense>
#include <string_view>
#include <iostream>

typedef std::pair<Eigen::Index, Eigen::Index> index_tuple;
typedef Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> MatrixX8u;
typedef Eigen::Array<uint8_t, Eigen::Dynamic, 1, Eigen::ColMajor> ArrayX8u;
typedef Eigen::Array<int8_t, Eigen::Dynamic, 1, Eigen::ColMajor> ArrayX8i;

class Abstract_Similarity_Matrix{
  public:
    /**
     * Iterates through all matrix entries in a fashion that respects data dependencies of the Smith-Waterman algrithm.
     */
    virtual void iterate(const std::function<float(const char &, const char &)> &scoring_function,
                         float gap_penalty) = 0;
    virtual std::tuple<Eigen::Index, Eigen::Index, float> find_index_of_maximum() const = 0;
    virtual void print_matrix() const = 0;
    virtual float operator()(Eigen::Index row, Eigen::Index col) const = 0;
    virtual Eigen::VectorXf getTimings() const = 0;
};

class Similarity_Matrix : public Abstract_Similarity_Matrix{
  public:
    // unimplemented to prevent implicit construction
    //Similarity_Matrix();

    /**
     * Constructor for generating a zero-initalized similarity matrix.
     *
     * @param sequence_x Sequence along the x-axis of the matrix.
     * @param sequence_y Sequence along the y-axis.
     */
    Similarity_Matrix(std::string_view sequence_x, std::string_view sequence_y);
    void iterate(const std::function<float(const char &, const char &)> &scoring_function,
                 float gap_penalty) override;
    void iterate_column(const std::function<float(const char &, const char &)> &scoring_function,
                        float gap_penalty);
    std::tuple<Eigen::Index, Eigen::Index, float> find_index_of_maximum() const override;
    void print_matrix() const override;
    const Eigen::MatrixXf &get_matrix() const;
    float operator()(Eigen::Index row, Eigen::Index col) const override { return raw_matrix(row, col); };
    Eigen::VectorXf getTimings() const override;
#ifdef USEOMP
    int sm_nthreads;
    int sm_finegrain_type;
#ifdef MTSIMD
    int sm_mt_simd;
#endif
#endif

  private:
    Eigen::VectorXf sm_timings;
    Eigen::MatrixXf raw_matrix; // column major
    std::string_view sequence_x;
    std::string_view sequence_y;
    float sm_iter_ad_read_time;  //anti-diagonal
    Eigen::VectorXf sm_iter_ad_i_times;  //anti-diagonal
};

class Similarity_Matrix_Skewed: public Abstract_Similarity_Matrix {
  public:
    Similarity_Matrix_Skewed(std::string_view sequence_x, std::string_view sequence_y);
    void iterate(const std::function<float(const char &, const char &)> &scoring_function,
                 float gap_penalty) override;
    std::tuple<Eigen::Index, Eigen::Index, float> find_index_of_maximum() const override;
    void print_matrix() const override;
    virtual void print_matrix_raw() const;
    const MatrixX8u &get_matrix() const { return raw_matrix; };
    Eigen::VectorXf getTimings() const override;
    index_tuple rawindex2trueindex(index_tuple raw_index) const;
    index_tuple trueindex2rawindex(index_tuple true_index) const;
    float operator()(Eigen::Index row, Eigen::Index col) const override {
      auto [ri, rj] = trueindex2rawindex(index_tuple(col, row));//SWITCHING SEQX AND SEQY
      return raw_matrix(ri, rj);
    };
#ifdef USEOMP
    int sm_nthreads;
    int sm_finegrain_type;
#ifdef MTSIMD
    int sm_mt_simd;
#endif
#endif

  private:
    MatrixX8u raw_matrix; //column major
    Eigen::Index len_x;
    Eigen::Index len_y;
    Eigen::Index nrows;
    Eigen::Index ncols;
    ArrayX8u sequence_x;
    //Eigen::ArrayXi sequence_y;
    ArrayX8u inv_sequence_y;
    float sm_iter_ad_read_time;  //anti-diagonal
    Eigen::VectorXf sm_iter_ad_i_times;  //anti-diagonal

};
#endif
