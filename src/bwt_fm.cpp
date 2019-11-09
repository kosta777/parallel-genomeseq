#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <vector>
#include <map>
#include <unordered_map>

template<class T>
void printVector(std::vector<std::vector<T>> const &mat) {
  for (const std::vector<T> &row: mat) {
    for (T val: row) {
      std::cout << val << " ";
    }
    std::cout << '\n';
  }
}

class BWTFM {
  public:
    std::string reference_genome$;
    std::string Alphabet = "$ACGT";

    //std::string FL_table;
    std::vector<std::pair<std::string, std::string>> FL_table;
    //FM index uses LF mapping
    std::string FM_F; //first: suffix beginning
    std::string FM_L; //last: char preceding suffix
    std::map<char, int> nt_2_index = {
        {'$', 0},
        {'A', 1},
        {'C', 2},
        {'G', 3},
        {'T', 4},
    };
    std::map<char, int> FM_C;
    std::vector<int> FM_C_vec;
    std::map<char, std::vector<int>> FM_Occ;
    std::vector<std::vector<int>> FM_Occ_mat;

    std::string query_pattern;
    int pattern_sp;
    int pattern_ep;

    void create_FL_table() {
      int len_rg$ = reference_genome$.size();
      FL_table = {{reference_genome$.substr(0, len_rg$ - 1), "$"}};

      int i = 0;
      while (i < len_rg$ - 1) {
        //		std::cout<<i<<std::endl;
        //		std::cout<<reference_genome$.substr(i+1,len_rg$-1)<<std::endl;
        //                std::cout<<reference_genome$.at(i)<<std::endl;

        //initialize tmp variable (char->str)
        std::string F_tmp;
        F_tmp = reference_genome$.substr(i + 1, len_rg$ - 1) + reference_genome$.substr(0, i);

        std::string L_tmp;
        L_tmp.push_back(reference_genome$.at(i));

        //append
        FL_table.emplace_back(F_tmp, L_tmp);
        i += 1;
      }

      std::vector<std::pair<std::string, std::string>> to_sort = FL_table;
      std::sort(to_sort.begin(), to_sort.end());

      FL_table = to_sort;

      std::cout << "FL creation finished" << std::endl;
      for (auto p : FL_table) {
        std::cout << p.first << ", " << p.second << std::endl;

        FM_F.push_back(p.first.at(0));
        FM_L.push_back(p.second.at(0));
      }
    }

    void create_C() {
      FM_C['$'] = 0;
      FM_C_vec.push_back(0);
      for (size_t i = 1; i < Alphabet.size(); i++) {
//	    std::cout<<Alphabet[i]<<std::endl;
        auto count =
            FM_C[Alphabet[i - 1]] + std::count(reference_genome$.begin(), reference_genome$.end(), Alphabet[i - 1]);
        FM_C[Alphabet[i]] = count;
        FM_C_vec.push_back(count);
      }
    }

    void create_Occ() {
      std::cout << "Create_Occ CALLED" << std::endl;
      //initialization
      std::cout << FM_L << std::endl;
      std::vector<int> init_vec(FM_L.size(), 0);
      for (char i : Alphabet) {
        FM_Occ[i] = init_vec;
      }
      char tmp_char = FM_L.at(0);
      FM_Occ['$'].at(2) = 499;
      FM_Occ[tmp_char][3] = 500;
      for (char i : Alphabet) {
        FM_Occ[i] = init_vec;
      }
    }

    void create_Occ_mat() {
      std::cout << "Create_Occ_mat CALLED" << std::endl;
      //initialization
      FM_Occ_mat.resize(Alphabet.size(), std::vector<int>(FM_L.size()));

      std::cout << nt_2_index[FM_L.at(0)] << std::endl;
      FM_Occ_mat[nt_2_index[FM_L.at(0)]][0] = 1;

      for (size_t j = 1; j < FM_L.size(); j++) {
        for (size_t i = 0; i < Alphabet.size(); i++) {
          FM_Occ_mat[i][j] = FM_Occ_mat[i][j - 1];
        }
        FM_Occ_mat[nt_2_index[FM_L.at(j)]][j] += 1;
      }

    }

    //Need to keep track of original coordinates for read mapping
    void locate_pattern() {
      std::cout << "locating pattern" << std::endl;
      std::string P = query_pattern;
      auto i = P.size() - 1;
      int sp = FM_C_vec[nt_2_index[P.at(i)]] + 1;
      int ep = FM_C_vec[nt_2_index[P.at(i)] + 1];
      std::cout << "i=" << i << ": " << sp << " , " << ep << std::endl;

      while ((sp <= ep) & (i >= 1)) {
        char c = P.at(i - 1);
        int c_idx = nt_2_index[c];
        sp = FM_C_vec[c_idx] + FM_Occ_mat[c_idx][sp - 2] + 1;
        ep = FM_C_vec[c_idx] + FM_Occ_mat[c_idx][ep - 1];
        i -= 1;
        std::cout << "i=" << i << ": " << sp << " , " << ep << std::endl;
      }
      pattern_sp = sp;
      pattern_ep = ep;
    }

    //use vector of pairs and sort the first element
    //.first=F, .second=L
    //Naive suffix array construction
    //(+) Can improve using Ukkonen's algo
    void create_suffix_array() {
      if (reference_genome$.find('$') != std::string::npos) {
        std::cout << "check reference genome" << std::endl;
      } else {
        reference_genome$ += "$";
        int len_rg$ = reference_genome$.size();
        FL_table = {{reference_genome$.substr(0, len_rg$ - 1), "$"}};

        int i = 0;
        while (i < len_rg$ - 1) {
//		std::cout<<i<<std::endl;
//		std::cout<<reference_genome$.substr(i+1,len_rg$-1)<<std::endl;
//                std::cout<<reference_genome$.at(i)<<std::endl;

          //initialize tmp variable (char->str)
          std::string L_char;
          L_char.push_back(reference_genome$.at(i));

          //append
          FL_table.emplace_back(reference_genome$.substr(i + 1, len_rg$ - 1), L_char);
          i += 1;
        }
        std::cout << "FL creation finished" << std::endl;
        for (const auto &p : FL_table) {
          std::cout << p.first << ", " << p.second << std::endl;
        }
      }
    }
};

int main() {
  std::cout << "Hello world" << std::endl;
  std::string line;
  std::ifstream fa;

  fa.open("../data/data_small/genome.chr22.5K.fa");
  std::string reference_genome;

  int i = 0;
  while (std::getline(fa, line)) {
    std::cout << "line number: " << i << std::endl;
    i++;

    if (line.find_first_not_of("ATCGNX") != std::string::npos) {
      std::cout << "metadata or erroneous data here" << std::endl;
      std::cout << line << std::endl;
    } else {
      reference_genome += line;
    }
  }

  BWTFM bwt1;
  bwt1.reference_genome$ = "TAGAGA$";

  bwt1.create_FL_table();

  std::cout << "*****************" << std::endl;
  auto myvec = bwt1.FL_table;
  for (const auto &p : myvec) {
    std::cout << p.first << ", " << p.second << std::endl;
  }

  std::cout << "*********************" << std::endl;
  auto to_sort = myvec;
  std::sort(to_sort.begin(), to_sort.end());
  for (unsigned int i = 0; i < to_sort.size(); i++) {
    auto p = to_sort[i];
    std::cout << p.first << ", " << p.second << std::endl;
    std::cout << bwt1.FM_F.at(i) << " , " << bwt1.FM_L.at(i) << std::endl;
  }

  bwt1.create_C();
  bwt1.create_Occ_mat();

  for (const auto &t : bwt1.FM_C)
    std::cout << t.first << " "
              << t.second << "\n";

  std::cout << "printing FM_Occ_mat: " << std::endl;
  printVector(bwt1.FM_Occ_mat);

  bwt1.query_pattern = "AGA";
  bwt1.locate_pattern();
  std::cout << bwt1.pattern_sp << " , " << bwt1.pattern_ep << std::endl;

  fa.close();
  return 0;
}
