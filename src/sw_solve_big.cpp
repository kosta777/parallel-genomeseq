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
#include <omp.h>
#endif

int main(int argc, char **argv) {
  int nrepeat = 1;
#ifdef USEOMP
  omp_set_dynamic(0);
  int npiece;
  if (argc < 3) {
    std::cout << "Please specify number of pieces to break e.g: `sw_solve_big <npiece> <nrepeat>`" <<std::endl;
    return -1;
  } else {
    npiece = std::stoi(argv[1]);
    nrepeat = std::stoi(argv[2]);
//#pragma omp parallel default(none) shared(npiece, std::cout)
    std::cout << "[INFO] npiece: " << npiece << ", nrepeat:" << nrepeat << std::endl;//<< ", nthreads: " << omp_get_num_threads() << std::endl;
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
  double time_iter_avg = 0.0;
  unsigned long long num_cells = 0;
  int i = 0;
  std::vector<double> GCUPS_vec;
  auto overlaprate = 2.0;
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
        if (i == 1) {
          std::cout << "[INFO] Estimated Memory consumption " << (double)(matsize*sizeof(uint8_t))*1e-9 << "GB" << std::endl;
#ifdef USEOMP
          std::cout << "[INFO] Theoretical GCUPS on Leonhard: " << npiece*4.6/(fa_string.size()+2*(npiece-1)*overlaprate*row[2].size())*fa_string.size() << std::endl;
#endif
        }
        //std::cout << "progress: " << i << std::endl;
#ifdef USEOMP
        auto la = std::make_unique<OMPParallelLocalAligner<Similarity_Matrix_Skewed, SWAligner<Similarity_Matrix_Skewed>>>(row[2],fa_string,npiece*2,overlaprate);
#else
        auto la = std::make_unique<SWAligner<Similarity_Matrix_Skewed>>(row[2], fa_string);
#endif
        auto time_min = 9e20;
	auto time_iter_min = 9e20;
        for (auto j = 0; j < nrepeat; j++) {
          la->calculateScore();
          time_min = std::min(time_min, (double) (la->getTimings()[0]));
          time_iter_min = std::min(time_iter_min, (double) (la->getTimings()[1]));
        }

        time_avg += time_min;
	time_iter_avg += time_iter_min;
        GCUPS_vec.emplace_back(matsize/time_min*1e-3);
        num_cells += matsize;
      }
    }

    i++;
  }
  auto GCUPS = num_cells/time_avg*1e-3;
  auto GCUPS_iter = num_cells/time_iter_avg*1e-3;
  time_avg /= (i - 1);
  Eigen::Map<Eigen::ArrayXd> GCUPS_arr(&GCUPS_vec[0],i - 1);
  align_input.close();
  std::cout << "[INFO] Average SW iter_ad_read times: " << time_avg*1e-6 << "s, GCUPS:" << GCUPS << ", GCPUS per iteration: " << GCUPS_iter << std::endl;
  std::cout << "[INFO] GCUPS avg:" << GCUPS_arr.mean() << ", GCUPS std:" << sqrt((GCUPS_arr - GCUPS_arr.mean()).square().mean()) << std::endl;
  std::cout << GCUPS_arr.transpose() << std::endl;
}

