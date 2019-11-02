#include <iostream>
#include "localaligner.h"
#include "swaligner.h"

#include <string>

int main() {
  std::cout << "Hello world" << std::endl;

  std::string s1 = "Montag";
  std::string s2 = "Donnerstag";

  LocalAligner *la = new SWAligner();
  la->setFirstSequence(s1);
  la->setSecondSequence(s2);
  la->calculateScore();
  std::cout << la->getScore() << std::endl;

  return 0;
}

