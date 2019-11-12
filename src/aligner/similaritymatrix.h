#ifndef _SIMILARITY_MATRIX_H_
#define _SIMILARITY_MATRIX_H_

#include <Eigen/Dense>
#include <string_view>

typedef std::pair<Eigen::Index, Eigen::Index> index_tuple;

class Abstract_Similarity_Matrix{
  public:
    /**
     * Iterates through all matrix entries in a fashion that respects data dependencies of the Smith-Waterman algrithm.
     * Starts in the north-west corner and moves along anti-diagonals to the south-east corner of the matrix.
     *
     * @param callback A callback function that is called for each matrix entry (with a reference to the matrix itself as a first argument).
     */
    virtual void iterate_anti_diagonal(const std::function<double(const char &, const char &)> &scoring_function,
                                       double gap_penalty) = 0;
    virtual index_tuple find_index_of_maximum() const = 0;
    virtual void print_matrix() const = 0;
    virtual const Eigen::MatrixXd &get_matrix() const = 0;
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
    void iterate_anti_diagonal(const std::function<double(const char &, const char &)> &scoring_function,
                               double gap_penalty) override;
    index_tuple find_index_of_maximum() const override;
    void print_matrix() const override;
    const Eigen::MatrixXd &get_matrix() const override;
  private:
    Eigen::MatrixXd similarity_matrix; // column major
    std::string_view sequence_x;
    std::string_view sequence_y;
};

class Similarity_Matrix_Skewed: public Abstract_Similarity_Matrix {
  public:
    Similarity_Matrix_Skewed(std::string_view sequence_x, std::string_view sequence_y);
    void iterate_anti_diagonal(const std::function<double(const char &, const char &)> &scoring_function,
                               double gap_penalty) override {};
    index_tuple find_index_of_maximum() const override;
    void print_matrix() const override;
    const Eigen::MatrixXd &get_matrix() const override { return raw_matrix; };
    index_tuple rawindex2trueindex(index_tuple raw_index) const;
    index_tuple trueindex2rawindex(index_tuple true_index) const;
  private:
    Eigen::MatrixXd raw_matrix; //column major
    Eigen::Index len_x;
    Eigen::Index len_y;
    std::string_view sequence_x;
    std::string_view sequence_y;
};

#endif
