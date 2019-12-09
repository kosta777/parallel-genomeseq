#include <iostream>
#include <fstream>
#include <string>
#include <memory>

#include "localaligner.h"
#include "smithwaterman.h"

#ifdef USEMPI
#include <mpi.h>
#endif

const int read_size = 125;
struct read_output{
    char buff[read_size + 1];
    int pos_pred;
    double score;
};

int main(int argc, char* argv[])
{
#ifdef USEMPI
    int size, rank, errclass, resultlen;
    MPI_File fh;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Status status;
    MPI_Offset total_number_of_bytes, number_of_bytes;

    char err_buffer[MPI_MAX_ERROR_STRING];
    char *buff;

    const std::string fa_file_path = "data/query/P02232.fasta";
    const std::string output_file_path = "data/align_output.csv";
    const std::string stats_file_path = "data/uniprot/stats.txt";
    const std::string reffile_file_path = "data/uniprot/";

    //Initialize MPI 
    MPI_Init(&argc, &argv);
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    //Since size is used to determine the amount of reads per reader, assume that we have 1 node less since one node will be in charge of writing data to output
    if(size == 1)
    {
        std::cout<<"Cannot run with only 1 node."<<std::endl;
        MPI_File_close(&fh);
        MPI_Finalize();
        return 0;
    }
   
    //Read reference genome
    std::ifstream fa;
    fa.open(fa_file_path);
    std::string fa_string;
    std::string fa_line;
    int i = 0;
    while (std::getline(fa, fa_line))
    {
        if (i > 0)
            fa_string += fa_line;
        i++;
    }
    fa.close();

    std::ifstream stats;
    stats.open(stats_file_path);
    int filecnt = 0;
    stats >> filecnt;
    stats.close();

    int files_per_rank = filecnt/(size-1);
    std::cout<<files_per_rank<<std::endl;

    struct read_output out;
    MPI_Datatype Outputtype;
    MPI_Datatype type[3] = { MPI_CHAR, MPI_INT, MPI_DOUBLE };
    int blocklen[3] = { read_size+1, 1, 1 };
    MPI_Aint disp[3];
    disp[0] = (char*)&out.buff - (char*)&out;
    disp[1] = (char*)&out.pos_pred - (char*)&out;
    disp[2] = (char*)&out.score - (char*)&out;
    
    MPI_Type_create_struct(3, blocklen, disp, type, &Outputtype);
    MPI_Type_commit(&Outputtype);

    if(rank != size -1)
    {
        //Perform sequencing for every read 
        int ind = 0, off =0;
        i=0;
        std::string delimiter = ",";
        std::ofstream align_output;
        std::string output_line;
        std::string ref_string;
        while(i< files_per_rank)
        {
            std::ifstream reffile;
            reffile.open(reffile_file_path+std::to_string(rank * files_per_rank + i)+".fasta");
            std::cout<<"Rank"<<rank<<" opening" << reffile_file_path+std::to_string(rank * files_per_rank + i)+".fasta"<<std::endl;
            
            std::string input_line = "";
            //extract read from reffile
            std::string pom_line;
            //throw away first line
            std::getline(reffile, pom_line);
            while (std::getline(reffile, pom_line))
            {
               input_line += pom_line; 
            }
            reffile.close();
            std::string input_line_tmp = input_line;
            std::string input_header_line;
            std::string output_header_line;
            double score_tmp;
            int pos_pred_tmp;

            ind++;
            off++;
            {
                auto la = std::make_unique<SWAligner<Similarity_Matrix>>(input_line,fa_string);
                score_tmp = la->calculateScore();
                pos_pred_tmp = la->getPos();
            }


            if (i % 50 == 0) {
              std::cout<<"Rank "<<rank<< " progress: " << i << std::endl;
            }

            i++;
            //Send data to the writer node
            sprintf(out.buff, "%.126s", input_line.c_str());
            out.pos_pred = pos_pred_tmp;
            out.score = score_tmp;

            MPI_Send(&out, 1, Outputtype, size-1, 123, comm);   
            
       }
#ifdef VERBOSE
       std::cout<<"Rank "<<rank<<" done("<<ind<<")"<<std::endl;
#endif
    }
    else
    {      
        int count = 0;
        std::string input_header_line;
        std::string output_header_line;
        std::ofstream align_output;
        align_output.open(output_file_path);
        std::string output_line;
        while(count < files_per_rank+1)
        {
            std::string input_line;
            if(count == 0)
            {
                  output_header_line = "read";
                  output_header_line += ",pos_pred";
                  output_header_line += ",score";
                  align_output << output_header_line << "\n";
            }
            else
            {
                MPI_Status rec_status;
                MPI_Recv(&out, 1, Outputtype, MPI_ANY_SOURCE, 123, comm, &rec_status);
                align_output << out.buff << ", "
                   << out.pos_pred << ", "
                   << out.score << "\n";
            }

            count++;
#ifdef VERBOSE
            std::cout<<"Reader done with "<<count<<std::endl;
#endif
        }

        align_output.close();
#ifdef VERBOSE
        std::cout<<"Reader rank done"<<std::endl;
#endif
        
    }
           
    MPI_Barrier(comm);
    MPI_File_close(&fh);
    free(buff);
  
    MPI_Finalize();

#endif
    return 0;
}




