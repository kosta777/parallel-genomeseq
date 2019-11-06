#ifndef _SMITHWATERMAN_LINEAR_H_
#define _SMITHWATERMAN_LINEAR_H_

#include <string>

#include "boost/multi_array.hpp"

#include "localaligner.h"
#include "similaritymatrix.h"

typedef boost::multi_array<double, 2> array_type;
typedef array_type::index index;
typedef boost::array<index, 2> index_tupel;

class SWAligner : public LocalAligner {
    public:
        SWAligner();
        SWAligner(const std::string&, const std::string&);

        /*
         * Maximum score in similarity matrix.
         */
        double calculateScore() override;

        /*
         * Position of the first matched position in sequence_x. Corresponds to POS in SAM file.
         * 
         * @return 1-based (as opposed to standard CS 0-based) counter.
         */
        unsigned int getPos() const override;

    private:
        unsigned int pos;

        std::string sequence_x;
        std::string sequence_y;
        
        Similarity_Matrix similarity_matrix;

        void traceback(index_tupel, unsigned int&);
        void calculate_similarity_matrix();
        double on_each_iteration(const array_type&, index, index);
};

#endif