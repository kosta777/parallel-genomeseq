#include <iostream>
#include <memory>
#include <fstream>
#include <string>
#include <vector>
#include "localaligner.h"
#include "smithwaterman.h"
#include "similaritymatrix.h"
#ifdef USEOMP
#include "plocalaligner.h"
#endif

int main(int argc, char **argv) {
#ifdef USEOMP
  int npiece;
  if (argc < 2) {
    std::cout << "Please specify number of pieces to break e.g: `sw_solve_big 16`" <<std::endl;
    return -1;
  } else {
    npiece = std::stoi(argv[1]);
//#pragma omp parallel default(none) shared(npiece, std::cout)
    std::cout << "[INFO] npiece: " << npiece << std::endl;//<< ", nthreads: " << omp_get_num_threads() << std::endl;
  }
#endif

  const std::string fa_file_path = "data/custom_ref_1.fa"; //fa contains reference
  const std::string input_file_path = "data/custom_reads_1.csv";

  std::ifstream fa;
  fa.open(fa_file_path);
  std::string fa_string;
  std::getline(fa, fa_string);
  fa.close();
  //std::cout <<"Reference string size: "<< fa_string.size() << std::endl;

  std::ifstream align_input;
  align_input.open(input_file_path);
  std::string input_line;
  std::string delimiter = ",";

  std::string input_header_line;
  std::string output_header_line;

  double time_avg = 0.0;
  unsigned long long num_cells = 0;
  int i = 0;
  std::vector<double> GCUPS_vec;
  while (std::getline(align_input, input_line)) {
    std::vector<std::string> row;
    std::string input_line_tmp = input_line;

    input_line_tmp += ",";
    size_t pos = 0;
    std::string token;
    while ((pos = input_line_tmp.find(delimiter)) != std::string::npos) {
      token = input_line_tmp.substr(0, pos);
      input_line_tmp.erase(0, pos + delimiter.length());
      row.push_back(token);
    }

    if (i > 0) {
      {
        auto matsize = row[2].size()*fa_string.size();
        if (i == 1) std::cout << "[INFO] Estimated Memory consumption " << (double)(matsize*sizeof(uint8_t))*1e-9 << "GB" << std::endl;
        //std::cout << "progress: " << i << std::endl;
#ifdef USEOMP
        auto la = std::make_unique<OMPParallelLocalAligner<Similarity_Matrix_Skewed, SWAligner<Similarity_Matrix_Skewed>>>(row[2],fa_string,npiece,2.0);
#else
        auto la = std::make_unique<SWAligner<Similarity_Matrix_Skewed>>(row[2], fa_string);
#endif
        la->calculateScore();
        auto time = la->getTimings()[0];
        time_avg += time;
        GCUPS_vec.emplace_back(matsize/time*1e-3);
        num_cells += matsize;
      }
    }

    i++;
  }
  auto GCUPS = num_cells/time_avg*1e-3;
  time_avg /= (i - 1);
  Eigen::Map<Eigen::ArrayXd> GCUPS_arr(&GCUPS_vec[0],i - 1);
  align_input.close();
  std::cout << "[INFO] Average SW iter_ad_read times: " << time_avg*1e-6 << "s, GCUPS:" << GCUPS << std::endl;
  std::cout << "[INFO] GCUPS avg:" << GCUPS_arr.mean() << "GCUPS std:" << sqrt((GCUPS_arr - GCUPS_arr.mean()).square().mean()) << std::endl;
  std::cout << GCUPS_arr.transpose() << std::endl;
}

