#include "smithwatermanlinear.h"

#include <functional>
#include <iostream>

#include "similaritymatrix.h"

SWAligner::SWAligner(const std::string& first_sequence, const std::string& second_sequence) :
    sequence_x(first_sequence), sequence_y(second_sequence), similarity_matrix(first_sequence, second_sequence) {}

void SWAligner::traceback(index_tuple idx, unsigned int& preliminary_pos) {
    const Eigen::MatrixXd& matrix = similarity_matrix.get_matrix();

    // stopping criterion
    Eigen::Index index_x = idx.first;
    Eigen::Index index_y = idx.second;
    if(matrix(index_x-1, index_y-1) == 0) {
        sequence_x.insert(index_x-1, "(");
        sequence_y.insert(index_y-1, "(");

        pos = index_y;
        return;
    }

    double n1 = matrix(index_x-1, index_y-1);
    double n2 = matrix(index_x, index_y-1);
    double n3 = matrix(index_x-1, index_y);

    // move north-west
    if((n1 >= n2) && (n1 >= n3)) {
        index_tuple new_idx(index_x-1, index_y-1); 
        traceback(new_idx, preliminary_pos);
    }
    // move west
    else if ((n2 >= n1) && (n2 >= n3)) {
        sequence_x.insert(index_x, "-");
        index_tuple new_idx(index_x, index_y-1); 
        traceback(new_idx, preliminary_pos);
    }
    // move north
    else {
        sequence_y.insert(index_y, "-");
        index_tuple new_idx(index_x-1, index_y); 
        traceback(new_idx, preliminary_pos);
    }
}

double calculate_scoring_Smith_Waterman_linear(const double& west, const double& north, const double& north_west, const char& a, const char& b) {
    const double gap_penalty = 2;
    auto scoring_function = [](const char& a, const char& b) { return a == b ? 3.0 : -3.0; };
    std::vector<double> v {
      north_west + scoring_function(a, b), 
      west - gap_penalty, 
      north - gap_penalty, 
      0
    };
    std::vector<double>::iterator result = std::max_element(v.begin(), v.end());
    return *result;
}

double SWAligner::on_each_iteration(const Eigen::MatrixXd& matrix, index_tuple idx) {        
    char a = sequence_x[idx.first-1];
    char b = sequence_y[idx.second-1];

    return calculate_scoring_Smith_Waterman_linear(
        matrix(idx.first, idx.second-1), // west
        matrix(idx.first-1, idx.second), // north
        matrix(idx.first-1, idx.second-1), // north west
        a, 
        b
    );
}

void SWAligner::calculate_similarity_matrix() {
    // https://stackoverflow.com/questions/2298242/callback-functions-in-c
    using namespace std::placeholders;

    auto cb = std::bind(&SWAligner::on_each_iteration, this, _1, _2);
    similarity_matrix.iterate_anti_diagonal(cb);
}

double SWAligner::calculateScore() {
    calculate_similarity_matrix();

    similarity_matrix.print_matrix();

    index_tuple max_idx = similarity_matrix.find_index_of_maximum();

    sequence_x.insert(max_idx.first, ")");
    sequence_y.insert(max_idx.second, ")");

    traceback(max_idx, pos);
    std::cout << "POS: " << pos << std::endl;

    std::cout << sequence_x << std::endl;
    std::cout << sequence_y << std::endl;

    return similarity_matrix.get_matrix()(max_idx.first, max_idx.second);
}

unsigned int SWAligner::getPos() const {
    return pos;
}