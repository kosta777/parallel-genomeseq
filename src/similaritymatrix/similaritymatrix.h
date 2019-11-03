#ifndef _SIMILARITY_MATRIX_H_
#define _SIMILARITY_MATRIX_H_

#include "boost/multi_array.hpp"

typedef boost::multi_array<double, 2> array_type;
typedef array_type::index index;
typedef boost::array<index, 2> index_tupel;

class Similarity_Matrix {

    public:
        Similarity_Matrix();
        Similarity_Matrix(const std::string& sequence_x, const std::string& sequence_y);
        index_tupel find_index_of_maximum() const; 
        void print_matrix() const;
        void iterate_anti_diagonal(std::function<double(const array_type& matrix, index local_i, index k)> callback);
        const array_type& getMatrix() const;
    private:
        array_type similarity_matrix;
};

#endif