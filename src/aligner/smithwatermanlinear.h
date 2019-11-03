#ifndef _SMITHWATERMAN_LINEAR_H_
#define _SMITHWATERMAN_LINEAR_H_

#include <string>
#include "boost/multi_array.hpp"
#include "localaligner.h"
#include "../similaritymatrix/similaritymatrix.h"

typedef boost::multi_array<double, 2> array_type;
typedef array_type::index index;
typedef boost::array<index, 2> index_tupel;

class SWAligner : public LocalAligner {
    public:
        SWAligner();
        SWAligner(const std::string& first_sequence, const std::string& second_sequence);

        void calculateScore() override;
        double getScore() const override;
        int getPos() const override;
    private:
        std::string sequence_x;
        std::string sequence_y;
        Similarity_Matrix similarity_matrix;

        void traceback(index_tupel idx);
        void calculate_similarity_matrix();
        double on_each_iteration(const array_type& matrix, index local_i, index k);
};

#endif