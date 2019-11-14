#include <iostream>
#include <memory>
#include <chrono>
#include <fstream>
#include <string>
#include <vector>
#include "localaligner.h"
#include "smithwaterman.h"
#include <Eigen/Core>

#ifdef USEOMP
#include <cstdio>
#include <omp.h>
#endif

int main(int argc, char** argv) {
#ifdef USEOMP
  std::string argv1 = argv[1];
  int arg_nreads = std::stoi(argv[2]);
  int arg_nthreads = std::stoi(argv[3]);

  if (argv1.compare("0")==0) {
    std::cout<<"Hello omp"<<std::endl;
    #pragma omp parallel for default(none) num_threads(3)
    for(int n=0; n<6; ++n) printf("Thread: %d, Num: %d\n", omp_get_thread_num(), n);
  }

  else if (argv1.compare("solve_small")==0) {
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
  
    double score_tmp;
    int pos_pred_tmp;
    i = 0;

    Eigen::VectorXf calculateScore_times = Eigen::VectorXf::Zero(arg_nreads);
    while (std::getline(align_input, input_line)) {
      std::cout<<"*** i="<<i<<" ***"<<std::endl;
      if (i>arg_nreads){
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
          auto la = std::make_unique<SWAligner>(row[2],fa_string);
	  la->sw_OMP_nthreads = arg_nthreads;
	  std::cout<<"sw_OMP_nthreads: "<<la->sw_OMP_nthreads<<std::endl;

	  auto start = std::chrono::high_resolution_clock::now();
          score_tmp = la->calculateScore();
          auto end = std::chrono::high_resolution_clock::now();
	  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end-start);
	  float duration_val = duration.count();
	  std::cout<<"i:"<<i<<", dur_val:"<<duration_val<<std::endl;
	  calculateScore_times(i-1) = duration_val;

          pos_pred_tmp = la->getPos();
	  std::cout<<"swaligner instance iter method: "<<la->sw_iter_method<<std::endl;
          std::cout<<"sw_iter_ad_read_times: "<<la->sw_iter_ad_read_time<<"us"<<std::endl;

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

    std::cout<<"calculateScore_times (us): "<<std::endl;
    std::cout<<calculateScore_times<<std::endl;

    align_input.close();
    align_output.close();
    std::cout << "Done, output file see: " << output_file_path << std::endl;
  }

  else if (argv1.compare("eigen1")==0) {
    std::cout<<"Hello omp eigen"<<std::endl;

    Eigen::VectorXd test_vec = Eigen::VectorXd::Zero(12);
    std::cout<<test_vec<<std::endl;
    std::cout<<test_vec.size()<<std::endl;

    omp_set_num_threads(arg_nthreads);
    #pragma omp parallel for default(shared)
    for(int i=0; i<test_vec.size(); i++){
      printf("i: %d, 1 Thread no: %d, \n", omp_get_thread_num(), i);
      printf("i: %d, 2 Thread no: %d, \n", omp_get_thread_num(), i);
      test_vec(i) = omp_get_thread_num();
      printf("i:%d, value assigned \n\n",i);
    }
    std::cout<<test_vec<<std::endl;
  }

  else if (argv1.compare("eigen2")==0) {
    std::cout<<"Hello omp eigen"<<std::endl;

    Eigen::VectorXd test_vec = Eigen::VectorXd::Zero(12);
    std::cout<<test_vec<<std::endl;
    std::cout<<test_vec.size()<<std::endl;

    omp_set_num_threads(arg_nthreads);
    #pragma omp parallel for default(shared)
    for(int i=0; i<test_vec.size(); i++){
      int cur_thread_num = omp_get_thread_num();
      printf("i: %d, 1 Thread no: %d, \n", cur_thread_num, i);
      printf("i: %d, 2 Thread no: %d, \n", cur_thread_num, i);
      test_vec(i) = cur_thread_num;
      printf("i:%d, value assigned \n\n",i);
    }
    std::cout<<test_vec<<std::endl;
  }

  else if (argv1.compare("eigen3")==0) {
    std::cout<<"Hello omp eigen"<<std::endl;

    Eigen::VectorXd test_vec = Eigen::VectorXd::LinSpaced(12,11,0);
    std::cout<<test_vec<<std::endl;
    std::cout<<test_vec.size()<<std::endl;

    omp_set_num_threads(arg_nthreads);
    int cur_thread_num;
    int i;
    int n = test_vec.size();
    #pragma omp parallel for shared(n) private(i,cur_thread_num)
    for(i=0; i<n; i++){
      cur_thread_num = omp_get_thread_num();
      printf("i: %d, 1 Thread no: %d, \n", cur_thread_num, i);
      printf("i: %d, 2 Thread no: %d, \n", cur_thread_num, i);
      test_vec(i) = cur_thread_num;
      printf("i:%d, value assigned \n\n",i);
    }
    std::cout<<test_vec<<std::endl;
  }

  else if (argv1.compare("hello1")==0) {
    omp_set_num_threads(arg_nthreads);
    #pragma omp parallel
    {
      int ID = omp_get_thread_num();
      printf("hello(%d)",ID);
      printf("world(%d) \n",ID);
    }
  }

  else if (argv1.compare("hello2")==0) {
    omp_set_num_threads(arg_nthreads);
    #pragma omp parallel
    {
      printf("hello(%d)",omp_get_thread_num());
      printf("world(%d) \n",omp_get_thread_num());
    }
  }

  else if (argv1.compare("hello3")==0) {
    omp_set_num_threads(arg_nthreads);
    #pragma omp parallel for
    for(int i=0;i<arg_nthreads;i++){
      printf("hello(%d)",omp_get_thread_num());
      printf("world(%d) \n",omp_get_thread_num());
    }
  }

  else if (argv1.compare("1")==0) {
    std::cout<<"Hello omp"<<std::endl;
    omp_set_num_threads(arg_nthreads);
    #pragma omp parallel for default(shared)
    for(int n=0; n<15; ++n) printf("Thread: %d, Num: %d\n", omp_get_thread_num(), n);
  }

  else if (argv1.compare("2")==0) {
    std::cout<<"Hello omp"<<std::endl;
    #pragma omp parallel for default(none) num_threads(3)
    for(int n=0; n<6; ++n) printf("Thread: %d, Num: %d\n", omp_get_thread_num(), n);
  }
#endif
}

