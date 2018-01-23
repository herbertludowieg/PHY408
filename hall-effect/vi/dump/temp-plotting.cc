#include <iostream>
#include <fstream>
#include <string>
#include "calc-class.h"

//compile with
//g++ temp-plotting.cc -lboost_iostreams -lboost_system -lboost_filesystem
int main(int argc, char ** argv) {
  if (argc < 2) {
    std::cout<<"ERROR Cannot parse arguments"<<std::endl
             <<"USAGE: "<<argv[0]<<" [file to analyze]"<<std::endl;
  }
  std::ifstream in(argv[1]);
  std::string line;
  Calculate calc;
  while (std::getline(in,line)){
    calc.input(in);
  }
  calc.print();
  return 0;
}
