#include <iostream>
#include <fstream>
#include "localaligner.h"
#include "swaligner.h"
#include <sstream>
#include <string>
#include <vector>

int main(int argc, char** argv) {
  std::string argv1 = argv[1];
  if (argv1.compare("swa")==0) {
    std::cout << "Hello world" << std::endl;

    std::string s1 = "Montag";
    std::string s2 = "Donnerstag";
  //std::string s2 = "sdkjfdfjksdnsdkfdskndkjsnfdMontaggmdfskmfsdds";

    LocalAligner *la = new SWAligner();
    la->setFirstSequence(s1);
    la->setSecondSequence(s2);
    la->calculateScore();
    std::cout << la->getScore() << std::endl;
    std::cout << la->getPos() << std::endl;
  }

  else if (argv1.compare("sw_solve_small") ==0){
    std::string fa_file_path = "data/data_small/genome.chr22.5K.fa"; //fa contains reference
    std::string input_file_path = "data/data_small_ground_truth.csv";
    std::string output_file_path = "data/align_output.csv";

    std::cout<<"Hello sw_solve_small"<<std::endl;

    std::ifstream fa;
    fa.open(fa_file_path);
    std::string fa_string;
    std::string fa_line;
    int i=0;
    while (std::getline(fa,fa_line)){
      if (i>0){
        fa_string += fa_line;
      }
      i++;
    }
    fa.close();

    #ifdef DEBUG
    std::cout<<"fa stuff: "<<std::endl;
    std::cout<<fa_string<<std::endl;
    std::cout<<fa_string.size()<<std::endl;
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

    int score_tmp;
    int pos_pred_tmp;
    i=0;
    while (std::getline(align_input,input_line)){
      std::vector<std::string> row;
      std::string input_line_tmp=input_line;

      input_line_tmp += ",";
      size_t pos = 0;
      std::string token;
      while ((pos = input_line_tmp.find(delimiter)) != std::string::npos) {
	token = input_line_tmp.substr(0, pos);
	input_line_tmp.erase(0, pos + delimiter.length());
	row.push_back(token);
      }

      if (i==0){
	input_header_line = input_line;
	output_header_line = input_line;
	output_header_line += ",pos_pred";
	output_header_line += ",score";
	align_output << output_header_line << "\n";
      }
      else {
	#ifdef DEBUG
	std::cout<<"about to call SWAligner"<<std::endl;
	std::cout<<row[2]<<std::endl;
        #endif
	LocalAligner *la = new SWAligner();
	la->setFirstSequence(row[2]);
	la->setSecondSequence(fa_string);
	la->calculateScore();

	score_tmp = la->getScore();
	pos_pred_tmp = la->getPos();
	align_output << input_line << ", "
                     << pos_pred_tmp << ", "
                     << score_tmp <<"\n";
      }

      if (i%50==0) {
        std::cout<<"progress: "<<i<<std::endl;
      }

      i++;
    }
    align_input.close();
    align_output.close();
    std::cout<<"Done, output file see: "<<output_file_path<<std::endl;
  }
  return 0;
}
