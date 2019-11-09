#ifndef _SMITHWATERMAN_LINEAR_H_
#define _SMITHWATERMAN_LINEAR_H_

#include <string>

#include <Eigen/Dense>

#include "localaligner.h"
#include "similaritymatrix.h"

typedef std::pair<Eigen::Index, Eigen::Index> index_tuple;

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

        void traceback(index_tuple, unsigned int&);
        void calculate_similarity_matrix();
        double on_each_iteration(const Eigen::MatrixXd&, index_tuple);
};

#endif