#include <iostream>
#include <string>
#include <algorithm> 
#include <vector>
#include "boost/multi_array.hpp"

typedef boost::multi_array<double, 2> array_type;
typedef array_type::index index;

double scoring_function(const char* a, const char* b) {
    if(*a == *b) {
        return 3;
    }
    else {
        return -3;
    }
}

double calculate_scoring_Smith_Waterman_linear(const double* west, const double* north, const double* north_west, const char* a, const char* b) {
    const double gap_penalty = 2;
    std::vector<double> v {
        *north_west + scoring_function(a, b), 
        *west - gap_penalty, 
        *north - gap_penalty, 
        0
    };
    std::vector<double>::iterator result = std::max_element(v.begin(), v.end());
    return *result;
}

void print_matrix(array_type& matrix) {
    for(size_t x = 0; x < matrix.shape()[0]; x++) {
        for(size_t y = 0; y < matrix.shape()[1]; y++) {
            std::cout << matrix[x][y] << "\t";
        }
        std::cout << "\n";
    }
}


int main()
{
    std::string sequence_A = "GGTTGACTA";
    std::string sequence_B = "TGTTACGG";

    int dim_x = sequence_A.size() + 1;
    int dim_y = sequence_B.size() + 1;
    array_type A(boost::extents[dim_x][dim_y]);

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

            char a = sequence_A[local_i-1];
            char b = sequence_B[k-1];

            A[local_i][k] = calculate_scoring_Smith_Waterman_linear(
                &A[local_i][k-1], // west
                &A[local_i-1][k], // north
                &A[local_i-1][k-1], // north west
                &a, 
                &b
            );

            --local_i;
        }
    }

    print_matrix(A);

    std::cout << "Max: " << *std::max_element( A.origin(), A.origin() + A.num_elements());

    return 0;
}
