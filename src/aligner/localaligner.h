#ifndef _LOCAL_ALIGNER_H_
#define _LOCAL_ALIGNER_H_

#include<string>

class LocalAligner {
  public:
    virtual void calculateScore() = 0;
    virtual double getScore() const;
    virtual int getPos() const;
};

#endif