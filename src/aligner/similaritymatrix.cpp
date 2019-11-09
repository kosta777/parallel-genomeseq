#include "similaritymatrix.h"
#include <iostream>
#include <algorithm> 
#include <Eigen/Dense>

Similarity_Matrix::Similarity_Matrix(std::string_view sequence_x, std::string_view sequence_y) :
    similarity_matrix(sequence_x.size() + 1, sequence_y.size() + 1) 
{
    similarity_matrix.setZero();
}
  
index_tuple Similarity_Matrix::find_index_of_maximum() const {
    Eigen::Index x = 0, y = 0;
    similarity_matrix.maxCoeff(&x, &y);

    std::cout << "Maximum is " << similarity_matrix(x, y) << " @ (" << x << ", " << y << ")" << std::endl;

    index_tuple max(x, y);
    return max;
}

void Similarity_Matrix::print_matrix() const {
    std::cout << similarity_matrix << std::endl;
}

void Similarity_Matrix::iterate_anti_diagonal(const std::function<double(const Eigen::MatrixXd &,
                                                                         index_tuple)> &callback) {
    const unsigned int dim_x = similarity_matrix.rows();
    const unsigned int dim_y = similarity_matrix.cols();

    for(Eigen::Index i = 1; i < dim_x + dim_y - 2; ++i) {
        Eigen::Index local_i = i;
        Eigen::Index starting_k = 1;
        Eigen::Index ending_k = i;
        if(local_i > dim_x - 1) {
            local_i = dim_x - 1;
            starting_k = i - local_i + 1;
            ending_k = starting_k + dim_x - 2;
        }
        if(ending_k > dim_y - 1) {
            ending_k = dim_y - 1;
        }
        for(Eigen::Index k = starting_k; k <= ending_k; ++k) {
            index_tuple idx(local_i, k);
            similarity_matrix(local_i, k) = callback(similarity_matrix, idx);

            --local_i;
        }
    }
}

const Eigen::MatrixXd& Similarity_Matrix::get_matrix() const  {
    return similarity_matrix;
}