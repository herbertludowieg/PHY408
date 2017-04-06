#include "calc-class.h"

Calculate::Calculate () {}
Calculate::~Calculate () {}

bool Calculate::input ( std::istream & in ) {
  std::string line;
  std::getline(in,line,';');
  bfield_.push_back(std::atof(line.c_str()));
  std::getline(in,line,';');
  bf_voltage_.push_back(std::atof(line.c_str()));
  std::getline(in,line,';');
  ae_voltage_.push_back(std::atof(line.c_str()));
  if (line == "")
    return false;
  else
    return true;
}

void Calculate::least_squares () {
  
}
void Calculate::plot () {
  //least_squares ();
  Gnuplot gp;
  gp<<"plot'-'\n";
  gp.send1d(bfield_);
}
