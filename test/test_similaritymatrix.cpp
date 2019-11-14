#include <iostream>

#include <Eigen/Dense>

#include "../src/aligner/similaritymatrix.h"
#include "gtest/gtest.h"

namespace {
    double increment_by_one(const Eigen::MatrixXd& matrix, index_tuple idx) {
        return matrix(idx.first-1, idx.second-1) + matrix(idx.first-1, idx.second) + matrix(idx.first, idx.second-1) + 1;
    }

    double do_not_increment(const Eigen::MatrixXd& matrix, index_tuple idx) {
        return matrix(idx.first-1, idx.second-1) + matrix(idx.first-1, idx.second) + matrix(idx.first, idx.second-1);
    }


    TEST(Similarity_Matrix, Iterate_anti_diagonal) {
        Similarity_Matrix sm("FOO", "FOOBAR");
        sm.iterate_anti_diagonal(increment_by_one);

        Eigen::MatrixXd refrence_m(4, 7);
        refrence_m <<   0,  0,  0,  0,  0,  0,  0,
                        0,  1,  2,  3,  4,  5,  6,
                        0,  2,  6,  12, 20, 30, 42,
                        0,  3,  12, 31, 64, 115,188;

        ASSERT_PRED2( 
            [](const Eigen::MatrixXd &lhs, const Eigen::MatrixXd &rhs) { return lhs.isApprox(rhs, 1e-4);},
            refrence_m, 
            sm.get_matrix()
        );

        ASSERT_EQ(sm.find_index_of_maximum().first, 3);
        ASSERT_EQ(sm.find_index_of_maximum().second, 6);
    }

    TEST(Similarity_Matrix, Iterate_anti_diagonal_without_incrementing) {
        Similarity_Matrix sm("FOOBARFOO", "FOOBAR");
        sm.iterate_anti_diagonal(do_not_increment);

        Eigen::MatrixXd refrence_m_zeros(10, 7);
        refrence_m_zeros.setZero();

        ASSERT_PRED2( 
            [](const Eigen::MatrixXd &lhs, const Eigen::MatrixXd &rhs) { return lhs.isApprox(rhs, 1e-4);},
            refrence_m_zeros, 
            sm.get_matrix()
        );
    }


}