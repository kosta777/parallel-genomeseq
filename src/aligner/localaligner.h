#ifndef _LOCAL_ALIGNER_H_
#define _LOCAL_ALIGNER_H_

#include<string>

class LocalAligner {
  public:
    virtual void calculateScore() = 0;
    virtual double getScore() const = 0;
    virtual int getPos() const = 0;
};

#endif