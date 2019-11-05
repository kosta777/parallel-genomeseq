#ifndef _SIMILARITY_MATRIX_H_
#define _SIMILARITY_MATRIX_H_

#include "boost/multi_array.hpp"

typedef boost::multi_array<double, 2> array_type;
typedef array_type::index index;
typedef boost::array<index, 2> index_tupel;

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
        Similarity_Matrix(const std::string& sequence_x, const std::string& sequence_y);

        /**
         * Iterates through all matrix entries in a fashion that respects data dependencies of the Smith-Waterman algrithm.
         * Starts in the north-west corner and moves along anti-diagonals to the south-east corner of the matrix.
         * 
         * @param callback A callback function that is called for each matrix entry (with a reference to the matrix itself as a first argument).
         */
        void iterate_anti_diagonal(std::function<double(const array_type&, index, index)> callback);
        
        index_tupel find_index_of_maximum() const; 
        void print_matrix() const;
        const array_type& get_matrix() const;

    private:
        array_type similarity_matrix;
};

#endif