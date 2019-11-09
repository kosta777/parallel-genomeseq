#ifndef _LOCAL_ALIGNER_H_
#define _LOCAL_ALIGNER_H_

#include<string>

class LocalAligner {
  public:
    virtual double calculateScore() = 0;
    virtual unsigned int getPos() const = 0;

};

#endif
