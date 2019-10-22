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
	for (std::vector<T> row: mat) {
		for (T val: row) {
			std::cout << val << " ";
		}
		std::cout << '\n';
	}
}

class BWTFM
{
    public:
    std::string testmem = "hellogruezi";
    std::string reference_genome$;
    std::string Alphabet = "$ACGT";

//    std::string FL_table;
    std::vector<std::pair<std::string,std::string>> FL_table;
    //FM index uses LF mapping
    std::string FM_F; //first: suffix beginning
    std::string FM_L; //last: char preceding suffix
//    std::map<int,int> FM_C {{0,0}};
//    FM_C.push_back({3,5});
    std::map<char,int> nt_index = {
	    {'$',0},
	    {'A',1},
	    {'C',2},
	    {'G',3},
	    {'T',4},
    };
    std::map<char,int> FM_C;
    std::map<char,std::vector<int>> FM_Occ;
    std::vector<std::vector<int>> FM_Occ_mat;

    void printname()
    {
        std::cout<<"inside BWTFM, testmem is: "<<testmem<<std::endl;
    }

    void easy_day()
    {
        std::cout<<"inside BWTFM, ref is: "<<reference_genome$<<std::endl;
	testmem = "changed by easy_day";
    }

    void create_FL_table() 
    {
        int len_rg$ = reference_genome$.size();
        FL_table = { {reference_genome$.substr(0,len_rg$-1),"$"} };
    
        int i=0;
        while (i<len_rg$-1){
    //		std::cout<<i<<std::endl;
    //		std::cout<<reference_genome$.substr(i+1,len_rg$-1)<<std::endl;
    //                std::cout<<reference_genome$.at(i)<<std::endl;
    
	    //initialize tmp variable (char->str)
	    std::string F_tmp;
	    F_tmp = reference_genome$.substr(i+1,len_rg$-1) + reference_genome$.substr(0,i);
    
	    std::string L_tmp;
	    L_tmp.push_back(reference_genome$.at(i));
    
	    //append
	    FL_table.push_back({F_tmp,L_tmp });
	    i+=1;
        }
    
        std::vector<std::pair<std::string,std::string>> to_sort = FL_table;
        std::sort(to_sort.begin(),to_sort.end());
    
        FL_table = to_sort;
    
        std::cout<<"FL creation finished"<<std::endl;
        for (int i=0;i<FL_table.size();i++){
	    std::pair<std::string,std::string> p = FL_table[i];
	    std::cout<<p.first<<", "<<p.second<<std::endl;
    
	    FM_F.push_back(p.first.at(0));
	    FM_L.push_back(p.second.at(0));
    
	    /*
	    std::string L_tmp;
	    L_tmp.push_back(p.second.at(0));
	    FM_L.push_back(L_tmp);
	    */
    
        }
    }

    void create_C(){
	FM_C['$'] = 0;
	for (int i=1;i<Alphabet.size();i++){
//	    std::cout<<Alphabet[i]<<std::endl;
	    FM_C[Alphabet[i]] = FM_C[Alphabet[i-1]] + std::count(reference_genome$.begin(),reference_genome$.end(),Alphabet[i-1]);
	}
    }

    void create_Occ(){
	std::cout<<"Create_Occ CALLED"<<std::endl;
        //initialization
	std::cout<<FM_L<<std::endl;
	std::vector<int> init_vec(FM_L.size(),0);
        for (int i=0;i<Alphabet.size();i++){
            FM_Occ[Alphabet[i]] = init_vec;
	}
	char tmp_char = FM_L.at(0);
        FM_Occ['$'].at(2) = 499;
        FM_Occ[tmp_char][3] = 500;
        for (int i=0;i<Alphabet.size();i++){
            FM_Occ[Alphabet[i]] = init_vec;
	}
    }

    void create_Occ_mat(){
	std::cout<<"Create_Occ_mat CALLED"<<std::endl;
        //initialization
	FM_Occ_mat.resize(Alphabet.size(), std::vector<int>( FM_L.size() ));

	std::cout<<nt_index[FM_L.at(0)]<<std::endl;
        FM_Occ_mat[ nt_index[FM_L.at(0)] ][0] = 1;

	for (int j=1;j<FM_L.size();j++){
            for (int i=0;i<Alphabet.size();i++){
                FM_Occ_mat[i][j] = FM_Occ_mat[i][j-1];
            }
	    FM_Occ_mat[ nt_index[FM_L.at(j)] ][j] += 1;
        }
//        FM_Occ_mat[1][0] = 1;

/*
//	vector<vector<int>> FM_Occ_mat

	std::cout<<FM_L<<std::endl;
	std::vector<int> init_vec(FM_L.size(),0);
        for (int i=0;i<Alphabet.size();i++){
            FM_Occ[Alphabet[i]] = init_vec;
	}
	char tmp_char = FM_L.at(0);
        FM_Occ['$'].at(2) = 499;
        FM_Occ[tmp_char][3] = 500;
        for (int i=0;i<Alphabet.size();i++){
            FM_Occ[Alphabet[i]] = init_vec;
	}
    */
    }


    
    //use vector of pairs and sort the first element
    //.first=F, .second=L
    //Naive suffix array construction
    //(+) Can improve using Ukkonen's algo
    void create_suffix_array() 
    {
        if (reference_genome$.find("$")!=std::string::npos){
            std::cout<<"check reference genome"<<std::endl;
	}
	else{
	    reference_genome$ += "$";
	    int len_rg$ = reference_genome$.size();
            FL_table = { {reference_genome$.substr(0,len_rg$-1),"$"} };

	    int i=0;
	    while (i<len_rg$-1){
//		std::cout<<i<<std::endl;
//		std::cout<<reference_genome$.substr(i+1,len_rg$-1)<<std::endl;
//                std::cout<<reference_genome$.at(i)<<std::endl;

                //initialize tmp variable (char->str)
		std::string L_char;
		L_char.push_back(reference_genome$.at(i));

		//append
	        FL_table.push_back( {reference_genome$.substr(i+1,len_rg$-1),L_char });
		i+=1;
	    }
	    std::cout<<"FL creation finished"<<std::endl;
            for (int i=0;i<FL_table.size();i++){
	        std::pair<std::string,std::string> p = FL_table[i];
                std::cout<<p.first<<", "<<p.second<<std::endl;
	    }
	}
    }
};

int main()
{
    std::cout<<"Hello world"<<std::endl;
    std::string line;
    std::ifstream fa;

    fa.open("../data/data_small/genome.chr22.5K.fa");
    std::string reference_genome;

    int i = 0;
    while (std::getline(fa,line))
    {
	std::cout<<"line number: "<<i<<std::endl;
	i++;

	if (line.find_first_not_of("ATCGNX") != std::string::npos){
	    std::cout<<"metadata or erroneous data here"<<std::endl;
            std::cout<<line<<std::endl;
	}
	else{
	    reference_genome += line;
	}
    }

    BWTFM bwt1;
    bwt1.reference_genome$ = "TAGAGA$";
    /*
    bwt1.testmem = "yet another test";
    std::cout<<bwt1.testmem<<std::endl;
    std::cout<<"try memfunc: "<<std::endl;

    bwt1.easy_day();

    bwt1.printname();
    std::cout<<bwt1.testmem[0]<<bwt1.testmem[2]<<std::endl;
    bwt1.testmem += "$";
    bwt1.printname();

    std::string to_sort = bwt1.testmem;
    std::sort(to_sort.begin(),to_sort.end());

    std::cout<<to_sort<<std::endl;

    std::cout<<reference_genome$<<std::endl;

    int l = reference_genome.size();
    std::cout<<reference_genome.substr(0,l-1)<<std::endl;
*/

    bwt1.create_FL_table();

    /*
    std::vector<std::pair<std::string,std::string>> myvec = { {"bdsoij","jfuhd"},{"aijofd","lpddkd"} };
    myvec.push_back({"fdjsd","djfdjg"});
    for (int i=0;i<myvec.size();i++){
	std::pair<std::string,std::string> p = myvec[i];
        std::cout<<p.first<<", "<<p.second<<std::endl;
    }
    */
    std::cout<<"line129"<<std::endl;
    std::vector<std::pair<std::string,std::string>> myvec = bwt1.FL_table;
    for (int i=0;i<myvec.size();i++){
	std::pair<std::string,std::string> p = myvec[i];
        std::cout<<p.first<<", "<<p.second<<std::endl;
    }

    std::cout<<"*********************"<<std::endl;
    std::vector<std::pair<std::string,std::string>> to_sort = myvec;
    std::sort(to_sort.begin(),to_sort.end());
    for (int i=0;i<to_sort.size();i++){
	std::pair<std::string,std::string> p = to_sort[i];
        std::cout<<p.first<<", "<<p.second<<std::endl;
	std::cout<<bwt1.FM_F.at(i)<<" , "<<bwt1.FM_L.at(i)<<std::endl;
    }

    
    bwt1.create_C();
    bwt1.create_Occ_mat();


    for (auto& t : bwt1.FM_C)
        std::cout << t.first << " " 
                  << t.second << "\n";

//    for (i=0;i<5;i++)
//    bwt1.FM_Occ['A'][0] = 1;
//   bwt1.FM_Occ['G'][0] = 383;

/*
    for (auto& t : bwt1.FM_Occ)
            std::cout << t.first << " " 
                      << t.second[0] << " " 
                      << t.second[1] << " " 
                      << t.second[2] << " " 
                      << t.second[3] << "\n";
*/
    std::cout<<"printing FM_Occ_mat: "<<std::endl;
    printVector(bwt1.FM_Occ_mat);
//    std::cout<<bwt1.FM_Occ[bwt1.FM_L.at(1)][0]<< std::endl;
/*
    std::vector<int> init_vec(bwt1.FM_L.size(),0);
    for (int i=0;i<bwt1.Alphabet.size();i++){
        bwt1.FM_Occ[bwt1.Alphabet[i]] = init_vec;
    }
    char tmp_char = bwt1.FM_L.at(0);
    bwt1.FM_Occ[tmp_char][0] = 500;
*/
    bwt1.printname();
//    std::cout<<myvec<<std::endl;
//

    int a = 2;
    int b = 4;
    std::vector<std::vector<int>> mat;

    /*
    std::vector<std::vector<int>> mat {
				{ 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0 },
				{ 0, 0, 50, 0, 0 },
				{ 0, 0, 0, 0, 0 }
			};
*/
    mat.resize(a, std::vector<int>(b));
    mat[0][3] = 49;
    printVector(mat);
    fa.close();
    return 0;
}
