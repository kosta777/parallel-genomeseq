#include "similaritymatrix.h"

#include <iostream>
#include <string>
#include <algorithm> 

typedef boost::multi_array<double, 2> array_type;
typedef array_type::index index;
typedef boost::array<index, 2> index_tupel;


Similarity_Matrix::Similarity_Matrix(const std::string& sequence_x, const std::string& sequence_y) :
    similarity_matrix(boost::extents[sequence_x.size() + 1][sequence_y.size() + 1]) {}

index_tupel Similarity_Matrix::find_index_of_maximum() const {
    int distance = std::distance(
        similarity_matrix.origin(), 
        std::max_element( 
            similarity_matrix.origin(), 
            similarity_matrix.origin() + similarity_matrix.num_elements()
        )
    );

    const int dim_y = similarity_matrix.shape()[1];
    const unsigned int x = distance/dim_y;
    const unsigned int y = distance - x*dim_y;

    index_tupel idx = {{x, y}};

    std::cout << "Maximum is " << similarity_matrix[x][y] << " @ (" << x << ", " << y << ")" << std::endl;

    return idx;
}

void Similarity_Matrix::print_matrix() const {
    index dim_x = similarity_matrix.shape()[0];
    index dim_y = similarity_matrix.shape()[1];

    for(index x = 0; x < dim_x; x++) {
        for(index y = 0; y < dim_y; y++) {
            std::cout << similarity_matrix[x][y] << "\t";
        }
        std::cout << "\n";
    }
}

void Similarity_Matrix::iterate_anti_diagonal(std::function<double(const array_type& matrix, index local_i, index k)> callback) {
    const unsigned int dim_x = similarity_matrix.shape()[0];
    const unsigned int dim_y = similarity_matrix.shape()[1];

    for(index i = 1; i < dim_x + dim_y - 2; ++i) {
        index local_i = i;
        index starting_k = 1;
        index ending_k = i;
        if(local_i > dim_x - 1) {
            local_i = dim_x - 1;
            starting_k = i - local_i + 1;
            ending_k = starting_k + dim_x - 2;
        }
        if(ending_k > dim_y - 1) {
            ending_k = dim_y - 1;
        }
        for(index k = starting_k; k <= ending_k; ++k) {

            similarity_matrix[local_i][k] = callback(similarity_matrix, local_i, k);

            --local_i;
        }
    }
}

const array_type& Similarity_Matrix::get_matrix() const  {
    return similarity_matrix;
}