#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<stdlib.h>
#include"van-der-pauw.h"

int main(int argc, char * argv[]) {
  if (argc < 2) {
    std::cout << "ERROR cannot parse command line argument" << std::endl
              << "USAGE: " << argv[0] << " [Voltage data file]" << std::endl
              << "Voltage data file format:\nMust have a scale labeled scale."
              << "\nMust have a current labeled current." 
              << "Labels format:\nscale\n####\ncurrent\n#####\nV#"
              << "\n########\n# = a number"
              << std::endl;
    return -1;
  }
  std::ifstream in(argv[1]);
  VanDerPauw first;
  std::vector<VanDerPauw> voltages;
  //std::vector<VanDerPauw> * vptr = &voltages;
  first.reset_counter();
  first.PI = 3.141592653589793;
  first.t = 300e-6;
  while (first.input(in)) 
    voltages.push_back(first);
  first.vanderpauw();
}
