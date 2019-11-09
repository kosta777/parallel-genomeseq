#ifndef _SIMILARITY_MATRIX_H_
#define _SIMILARITY_MATRIX_H_

#include <Eigen/Dense>
#include <string_view>

typedef std::pair<Eigen::Index, Eigen::Index> index_tuple;

class Similarity_Matrix {
    public:
        // unimplemented to prevent implicit construction
        Similarity_Matrix();

        /**
         * Constructor for generating a zero-initalized similarity matrix.
         * 
         * @param sequence_x Sequence along the x-axis of the matrix.
         * @param sequence_y Sequence along the y-axis.
         */
        Similarity_Matrix(std::string_view sequence_x, std::string_view sequence_y);

        /**
         * Iterates through all matrix entries in a fashion that respects data dependencies of the Smith-Waterman algrithm.
         * Starts in the north-west corner and moves along anti-diagonals to the south-east corner of the matrix.
         * 
         * @param callback A callback function that is called for each matrix entry (with a reference to the matrix itself as a first argument).
         */
        void iterate_anti_diagonal(const std::function<double(const Eigen::MatrixXd &, index_tuple)> &callback);
        
        index_tuple find_index_of_maximum() const; 
        void print_matrix() const;
        const Eigen::MatrixXd& get_matrix() const;

    private:
        Eigen::MatrixXd similarity_matrix;
};

#endif