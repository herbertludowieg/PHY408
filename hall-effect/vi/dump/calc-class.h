#ifndef calc_h
#define calc_h
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <stdlib.h>

class Calculate {
  public:
  Calculate();
  ~Calculate();
  bool input ( std::istream & in );
  void print ( std::ostream & out = std::cout );
  void least_squares ();
  private:
  std::vector<double> bfield_,bf_voltage_,ae_voltage_;
  std::vector<double> bfparameters_,bfsigma_;
  std::vector<double> aeparameters_,aesigma_;
};

#endif
