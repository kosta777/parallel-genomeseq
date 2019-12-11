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

int main() {
  const std::string fa_file_path = "data/data_small/genome.chr22.5K.fa"; //fa contains reference
  const std::string input_file_path = "data/data_small_ground_truth.csv";
  const std::string output_file_path = "data/align_output.csv";

  std::cout << "Hello sw_solve_small" << std::endl;

  std::ifstream fa;
  fa.open(fa_file_path);
  std::string fa_string;
  std::string fa_line;
  int i = 0;
  while (std::getline(fa, fa_line)) {
    if (i > 0) {
      fa_string += fa_line;
    }
    i++;
  }
  fa.close();

#ifdef VERBOSE
  std::cout << "fa stuff: " << std::endl;
    std::cout << fa_string << std::endl;
    std::cout << fa_string.size() << std::endl;
#endif

  std::ifstream align_input;
  align_input.open(input_file_path);
  std::string input_line;
  std::string delimiter = ",";

  std::ofstream align_output;
  align_output.open(output_file_path);
  std::string output_line;

  std::string input_header_line;
  std::string output_header_line;

  float score_tmp;
  size_t pos_pred_tmp;
  double time_avg = 0.0;
  unsigned long long num_cells = 0;
  i = 0;
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

    if (i == 0) {
      input_header_line = input_line;
      output_header_line = input_line;
      output_header_line += ",pos_pred";
      output_header_line += ",score";
      align_output << output_header_line << "\n";
    } else {
#ifdef VERBOSE
      std::cout << "about to call SWAligner" << std::endl;
      std::cout << row[2] << std::endl;
#endif
      {
#ifdef USEOMP
        auto la = std::make_unique<OMPParallelLocalAligner<Similarity_Matrix_Skewed, SWAligner<Similarity_Matrix_Skewed>>>(row[2],fa_string,17,2.0);
#else
        auto la = std::make_unique<SWAligner<Similarity_Matrix_Skewed>>(row[2], fa_string);
#endif
        score_tmp = la->calculateScore();
        pos_pred_tmp = la->getPos();
        time_avg += la->getTimings()[0];
        num_cells += row[2].size()*fa_string.size();
      }
      align_output << input_line << ", "
                   << pos_pred_tmp << ", "
                   << score_tmp << "\n";
    }

    if (i % 50 == 0) {
      std::cout << "progress: " << i << std::endl;
    }

    i++;
  }
  auto GCUPs = num_cells/time_avg*1e-3;
  time_avg /= i;
  align_input.close();
  align_output.close();
  std::cout << "Average SW iter_ad_read times: " << time_avg << "us, GCUP:" << GCUPs << std::endl;
  std::cout << "Done, output file see: " << output_file_path << std::endl;
}

