#ifndef VANDERPAUW_H
#define VANDERPAUW_H

#include<string>
#include<iostream>
#include<cstdlib>
#include<fstream>
#include<vector>
#include<cmath>

class VanDerPauw {

  public:
  VanDerPauw();
  ~VanDerPauw();
  double PI;
  double t;
  double scale() const;
  double current() const;
  void reset_counter();
  //double voltages(vector<double>::const_iterator) const;
  void print(std::ostream & out /*= std::cout*/);
  bool input(std::istream & in);
  double find_resistivity(double v1, double v2);
  void vanderpauw(void);
  
  private:
  std::vector<double> voltages_;
  double scale_;
  double current_;
  int counter_;
  std::vector<double> resistivity_;
};
#endif
