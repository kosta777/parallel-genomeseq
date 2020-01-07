#include <iostream>
#include <memory>
#include <chrono>
#include <fstream>
#include <string>
#include <vector>
#include "localaligner.h"
#include "smithwaterman.h"
#include <Eigen/Core>

#include <set>
#include <iterator>
#include <algorithm>

#ifdef USEOMP
#include <cstdio>
#include <omp.h>
#endif

class CSVWriter
{
	std::string fileName;
	std::string delimeter;
	int linesCount;
 
public:
	CSVWriter(std::string filename, std::string delm = ",") :
			fileName(filename), delimeter(delm), linesCount(0)
	{}
	/*
	 * Member function to store a range as comma seperated value
	 */
	template<typename T>
	void addDatainRow(T first, T last, int append);
};
 
/*
 * This Function accepts a range and appends all the elements in the range
 * to the last row, seperated by delimeter (Default is comma)
 */
template<typename T>
void CSVWriter::addDatainRow(T first, T last, int append)
{
	std::fstream file;
	// Open the file in truncate mode if first line else in Append Mode
//	file.open(fileName, std::ios::out | (linesCount ? std::ios::app : std::ios::trunc));
	file.open(fileName, std::ios::out | (append ? std::ios::app : std::ios::trunc));


	// Iterate over the range and add each lement to file seperated by delimeter.
	for (; first != last; )
	{
		file << *first;
		if (++first != last)
			file << delimeter;
	}
	file << "\n";
	linesCount++;
 
	// Close the file
	file.close();
}

int main(int argc, char **argv) {
#ifdef USEOMP
  if (argc < 8) exit(-1);
  std::string argv1 = argv[1];
  int arg_nreads = std::stoi(argv[2]);
  int arg_nthreads = std::stoi(argv[3]);
  int arg_finegrain_type = std::stoi(argv[4]);
  std::string arg_timing_file_path = argv[5];
  std::string arg_ref_file_path = argv[6];
  std::string arg_reads_file_path = argv[7];

  if (argv1 == "0") {
    std::cout << "Hello omp" << std::endl;
#pragma omp parallel for default(none) num_threads(3)
    for (int n = 0; n < 6; ++n) printf("Thread: %d, Num: %d\n", omp_get_thread_num(), n);
  } else if (argv1 == "solve_small") {

    const std::string timing_file_path = arg_timing_file_path;

//    const std::string fa_file_path = "data/test_fa_custom1.fa"; //fa contains reference
    const std::string fa_file_path = arg_ref_file_path; //fa contains reference

//    const std::string input_file_path = "data/ground_truth_custom10k.csv";
    const std::string input_file_path = arg_reads_file_path;
    const std::string align_output_file_path = "data/ompfg_align_output.csv";

    int fa_file_has_header = 0;

    std::cout << "Hello sw_solve_small" << std::endl;

    std::ifstream fa;
    fa.open(fa_file_path);
    std::string fa_string;
    std::string fa_line;
    int i = 0;
    while (std::getline(fa, fa_line)) {
      if (i+1 > fa_file_has_header) {
        fa_string += fa_line;
      }
      i++;
    }
    fa.close();
    std::cout<<"omp_sw_solve_small line 99, i:  "<<i<<std::endl;

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
    align_output.open(align_output_file_path);
    std::string output_line;

    std::string input_header_line;
    std::string output_header_line;

    float score_tmp;
    Eigen::Index pos_pred_tmp;
    Eigen::VectorXf SM_timings_tmp;
    i = 0;

    Eigen::VectorXf calculateScore_times = Eigen::VectorXf::Zero(arg_nreads);
    Eigen::VectorXf read_sm_iterate_times = Eigen::VectorXf::Zero(arg_nreads);
    Eigen::VectorXf adsum_sm_iterate_times = Eigen::VectorXf::Zero(arg_nreads);
    while (std::getline(align_input, input_line)) {
      std::cout << "*** i=" << i << " ***" << std::endl;
      if (i > arg_nreads) {
        break;
      }
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
#ifdef MTSIMD
          auto la = std::make_unique<SWAligner<Similarity_Matrix_Skewed>>(row[2], fa_string);
	  la->sw_finegrain_type = -1;
#else
          auto la = std::make_unique<SWAligner<Similarity_Matrix>>(row[2], fa_string);
	  la->sw_finegrain_type = arg_finegrain_type;
#endif

	  la->sw_nthreads = arg_nthreads;
          auto start = std::chrono::high_resolution_clock::now();
          score_tmp = la->calculateScore();
          auto end = std::chrono::high_resolution_clock::now();
          auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
          float duration_val = duration.count();
	  
	  std::cout<<"MT SIMD status code: "<<la->sw_mt_simd<<std::endl;
          std::cout << "i:" << i << ", dur_val:" << duration_val << std::endl;
          calculateScore_times(i - 1) = duration_val;

          pos_pred_tmp = la->getPos();

	  SM_timings_tmp = la->getTimings();
          std::cout<<"sw_iter_ad_read_times: "<<SM_timings_tmp(0)<<"us"<<std::endl;
          read_sm_iterate_times(i-1) = SM_timings_tmp(0);

          std::cout<<"sw_iter_ad_i_times_sum: "<<SM_timings_tmp(1)<<"us"<<std::endl;
          adsum_sm_iterate_times(i-1) = SM_timings_tmp(1);


#ifdef VERBOSE
          std::cout << "OMP test stuff" << std::endl;
          std::cout << la->pub_max_score << std::endl;
#endif
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

    align_input.close();
    align_output.close();
    std::cout << "Done, align output file see: " << align_output_file_path << std::endl;
    
    std::cout << "calculateScore_times (us): " << std::endl;
    std::cout << calculateScore_times << std::endl;

    std::cout << "sm_iterate_times (us): " << std::endl;
    std::cout << read_sm_iterate_times << std::endl;

    std::cout << "adsum sm_iterate_times (us): " << std::endl;
    std::cout << adsum_sm_iterate_times << std::endl;


    std::ifstream timing_file_exist(timing_file_path);
    std::cout<<timing_file_exist.good()<<std::endl;
    CSVWriter timing_writer(timing_file_path);
    std::vector<std::string> header_tmp;
    std::vector<float> row_tmp;
    float arg_nreads_float = (float) arg_nreads;
    float arg_nthreads_float = (float) arg_nthreads;
    float arg_finegrain_type_float = (float) arg_finegrain_type;
    if (timing_file_exist.good() != 1){
//      header_tmp = {"n_reads","n_threads","time1","time2"};
      header_tmp = {"n_reads","n_threads","finegrain_type","avg_t_calcscore","avg_t_adread","avg_t_adisum"};
      timing_writer.addDatainRow(header_tmp.begin(),header_tmp.end(),0);
    }
//    row_tmp = {arg_nreads_float, arg_nthreads_float, 
    row_tmp = {arg_nreads_float, arg_nthreads_float, arg_finegrain_type_float,
	    calculateScore_times.mean(), read_sm_iterate_times.mean(), adsum_sm_iterate_times.mean()};
    timing_writer.addDatainRow(row_tmp.begin(),row_tmp.end(),1);

    /*
    for (i=0; i<arg_nreads; i++){
      row_tmp = {arg_nreads_float, arg_nthreads_float, calculateScore_times(i), read_sm_iterate_times(i) };
      timing_writer.addDatainRow(row_tmp.begin(),row_tmp.end());
    }
    */
    std::cout << "Done, timing file see: " << timing_file_path << std::endl;

    int a = 59/4;
    std::cout<<"testing123: "<<a<<std::endl;
    std::cout<<"END***** omp_sw_solve_small ************************************\n\n\n"<<std::endl;


  } else if (argv1 == "eigen1") {
    std::cout << "Hello omp eigen" << std::endl;

    Eigen::VectorXf test_vec = Eigen::VectorXf::Zero(12);
    std::cout << test_vec << std::endl;
    std::cout << test_vec.size() << std::endl;

    omp_set_num_threads(arg_nthreads);
#pragma omp parallel for default(shared)
    for (int i = 0; i < test_vec.size(); i++) {
      printf("i: %d, 1 Thread no: %d, \n", omp_get_thread_num(), i);
      printf("i: %d, 2 Thread no: %d, \n", omp_get_thread_num(), i);
      test_vec(i) = omp_get_thread_num();
      printf("i:%d, value assigned \n\n", i);
    }
    std::cout << test_vec << std::endl;
  } else if (argv1 == "eigen2") {
    std::cout << "Hello omp eigen" << std::endl;

    Eigen::VectorXf test_vec = Eigen::VectorXf::Zero(12);
    std::cout << test_vec << std::endl;
    std::cout << test_vec.size() << std::endl;

    omp_set_num_threads(arg_nthreads);
#pragma omp parallel for default(shared)
    for (int i = 0; i < test_vec.size(); i++) {
      int cur_thread_num = omp_get_thread_num();
      printf("i: %d, 1 Thread no: %d, \n", cur_thread_num, i);
      printf("i: %d, 2 Thread no: %d, \n", cur_thread_num, i);
      test_vec(i) = cur_thread_num;
      printf("i:%d, value assigned \n\n", i);
    }
    std::cout << test_vec << std::endl;
  } else if (argv1 == "eigen3") {
    std::cout << "Hello omp eigen" << std::endl;

    Eigen::VectorXf test_vec = Eigen::VectorXf::LinSpaced(12, 11, 0);
    std::cout << test_vec << std::endl;
    std::cout << test_vec.size() << std::endl;

    omp_set_num_threads(arg_nthreads);
    int cur_thread_num;
    int i;
    int n = test_vec.size();
#pragma omp parallel for shared(n) private(i, cur_thread_num)
    for (i = 0; i < n; i++) {
      cur_thread_num = omp_get_thread_num();
      printf("i: %d, 1 Thread no: %d, \n", cur_thread_num, i);
      printf("i: %d, 2 Thread no: %d, \n", cur_thread_num, i);
      test_vec(i) = cur_thread_num;
      printf("i:%d, value assigned \n\n", i);
    }
    std::cout << test_vec << std::endl;
  } else if (argv1 == "hello1") {
    omp_set_num_threads(arg_nthreads);
#pragma omp parallel
    {
      int ID = omp_get_thread_num();
      printf("hello(%d)", ID);
      printf("world(%d) \n", ID);
    }
  } else if (argv1 == "hello2") {
    omp_set_num_threads(arg_nthreads);
#pragma omp parallel
    {
      printf("hello(%d)", omp_get_thread_num());
      printf("world(%d) \n", omp_get_thread_num());
    }
  } else if (argv1 == "hello3") {
    omp_set_num_threads(arg_nthreads);
#pragma omp parallel for
    for (int i = 0; i < arg_nthreads; i++) {
      printf("hello(%d)", omp_get_thread_num());
      printf("world(%d) \n", omp_get_thread_num());
    }
  } else if (argv1 == "1") {
    std::cout << "Hello omp" << std::endl;
    omp_set_num_threads(arg_nthreads);
#pragma omp parallel for default(shared)
    for (int n = 0; n < 15; ++n) printf("Thread: %d, Num: %d\n", omp_get_thread_num(), n);
  } else if (argv1 == "2") {
    std::cout << "Hello omp" << std::endl;
#pragma omp parallel for default(none) num_threads(3)
    for (int n = 0; n < 6; ++n) printf("Thread: %d, Num: %d\n", omp_get_thread_num(), n);
  }
#endif
}

