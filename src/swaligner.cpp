#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <cmath>
#include <sys/time.h>

#include "swaligner.h"
#include "localaligner.h"

using namespace std;


void SWAligner::setFirstSequence(std::string seq)
{
   this->seqA = seq; 
}

void SWAligner::setSecondSequence(std::string seq)
{
    this->seqB = seq;
}

void SWAligner::calculateScore() 
{
    // initialize some variables
    int lengthSeqA = seqA.length();
    int lengthSeqB = seqB.length();
    
    // initialize matrix
    double matrix[lengthSeqA+1][lengthSeqB+1];
    for(int i=0;i<=lengthSeqA;i++)
    {
        for(int j=0;j<=lengthSeqB;j++)
        {
            matrix[i][j]=0;
        }
    }

    double traceback[4];
    int I_i[lengthSeqA+1][lengthSeqB+1];
    int I_j[lengthSeqA+1][lengthSeqB+1];


    //start populating matrix
    for (int i=1;i<=lengthSeqA;i++)
    {
        for(int j=0;j<=lengthSeqB;j++)
                {
                        #ifdef DEBUG
                        cout << i << " " << j << endl;
                        #endif
                        traceback[0] = matrix[i-1][j-1]+similarityScore(seqA[i-1],seqB[j-1]);
                        traceback[1] = matrix[i-1][j]+penalty;
                        traceback[2] = matrix[i][j-1]+penalty;
                        traceback[3] = 0;
                        matrix[i][j] = findMax(traceback,4);
                        switch(ind)
                        {
                            case 0:
                                I_i[i][j] = i-1;
                                I_j[i][j] = j-1;
                                break;
                            case 1:
                                I_i[i][j] = i-1;
                                I_j[i][j] = j;
                                                    break;
                            case 2:
                                I_i[i][j] = i;
                                                    I_j[i][j] = j-1;
                                                    break;
                            case 3:
                                I_i[i][j] = i;
                                                    I_j[i][j] = j;
                                                    break;
                        }
                }
    }

    #ifdef DEBUG
    // print the scoring matrix to console
    for(int i=1;i<lengthSeqA;i++)
    {
        for(int j=1;j<lengthSeqB;j++)
        {
            cout << matrix[i][j] << " ";
        }
        cout << endl;
    }
    #endif
    // find the max score in the matrix
    matrix_max = 0;
    int i_max=0, j_max=0;
    for(int i=1;i<lengthSeqA;i++)
    {
        for(int j=1;j<lengthSeqB;j++)
        {
            if(matrix[i][j]>matrix_max)
            {
                matrix_max = matrix[i][j];
                i_max=i;
                j_max=j;
            }
        }
    }
    
    #ifdef DEBUG
    cout << "Max score in the matrix is " << matrix_max << endl;
    #endif

    // traceback

    int current_i=i_max,current_j=j_max;
    int next_i=I_i[current_i][current_j];
    int next_j=I_j[current_i][current_j];
    int tick=0;
    char consensus_a[lengthSeqA+lengthSeqB+2],consensus_b[lengthSeqA+lengthSeqB+2];

    while(((current_i!=next_i) || (current_j!=next_j)) && (next_j!=0) && (next_i!=0))
    {

        if(next_i==current_i)  consensus_a[tick] = '-';                  // deletion in A
        else                   consensus_a[tick] = seqA[current_i-1];   // match/mismatch in A

        if(next_j==current_j)  consensus_b[tick] = '-';                  // deletion in B
        else                   consensus_b[tick] = seqB[current_j-1];   // match/mismatch in B

        current_i = next_i;
        current_j = next_j;
        next_i = I_i[current_i][current_j];
        next_j = I_j[current_i][current_j];
        tick++;
    }

    #ifdef DEBUG
    //print the consensus sequences
    cout<<endl<<" "<<endl;
    cout<<"Alignment:"<<endl<<endl;
    for(int i=0;i<lengthSeqA;i++){cout<<seqA[i];}; cout<<"  and"<<endl;
    for(int i=0;i<lengthSeqB;i++){cout<<seqB[i];}; cout<<endl<<endl;  
    for(int i=tick-1;i>=0;i--) cout<<consensus_a[i]; 
    cout<<endl;
    for(int j=tick-1;j>=0;j--) cout<<consensus_b[j];
    cout<<endl;
    #endif

}

double SWAligner::similarityScore(char a, char b)
{
    double result;
    if(a==b)
    {
        result=1;
    }
    else
    {
        result=penalty;
    }
    return result;
}

double SWAligner::findMax(double array[], int length)
{
    double max = array[0];
    ind = 0;

    for(int i=1; i<length; i++)
    {
        if(array[i] > max)
        {
            max = array[i];
            ind=i;
        }
    }
    return max;
}

double SWAligner::getScore() const
{
    return matrix_max;
}

