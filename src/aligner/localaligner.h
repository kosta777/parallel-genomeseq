#ifndef _LOCAL_ALIGNER_H_
#define _LOCAL_ALIGNER_H_

#include<string>

class LocalAligner {
  public:
    virtual void setFirstSequence(std::string) = 0;
    virtual void setSecondSequence(std::string) = 0;

    virtual void calculateScore() = 0;
    virtual double getScore() const = 0;
};

#endif
