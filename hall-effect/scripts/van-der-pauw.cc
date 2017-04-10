#include"van-der-pauw.h"

VanDerPauw::VanDerPauw(){}
VanDerPauw::~VanDerPauw(){}
double VanDerPauw::scale() const {return scale_;}
double VanDerPauw::current() const {return current_;}
void VanDerPauw::reset_counter(){counter_ = 0;}
void VanDerPauw::print(std::ostream & out = std::cout) {
  out << "Scale: " << scale() << " V" << std::endl << "Current: " << current() 
      << " A" << std::endl << "Voltages:" << std::endl;
  for (std::vector<double>::const_iterator i=voltages_.begin();
       i!=voltages_.end();++i) {
    out << *i << std::endl;
  }
}
bool VanDerPauw::input (std::istream & in) {
  std::string line;
  std::getline(in,line);
  if (line[0]=='#') {
    return true;
  } else if (line=="scale" && counter_==0) {
    std::getline(in,line);
    scale_=std::atof(line.c_str());
  } else if (line=="current" && counter_==1) {
    std::getline(in,line);
    current_=std::atof(line.c_str());
  } else if (line[0]=='V' && counter_>1) {
    std::getline(in,line);
    voltages_.push_back(std::atof(line.c_str()));
  } else if (line!=""){
    std::cout << "ERROR missing current/scale/V# tag.\n"
              << "Make sure proper tag format is being used.\n"
              << "To view format type executable name without args."
              << std::endl;
    return false;
  }
  if (line=="") {
    return false;
  } else {
  counter_ += 1;
  return true;
  }
}
double VanDerPauw::find_resistivity(double v1, double v2){
  double resistivity,f,r_val,ln2=std::log(2.0),other_f,res_error;
  r_val = ((v1-v2)/(v1+v2));
  f = 1-std::pow(r_val,2)*(ln2/2.0)-std::pow(r_val,4)*
           ((std::pow(ln2,2)/4.0)-(std::pow(ln2,3)/12.0));
  resistivity = ((PI*t)/(ln2))*(((v1+v2)*scale_)/(2.0*current_))*f;
  other_f = v1/v2;
  res_error = ((PI*t)/(ln2))*((std::sqrt(2.)*(0.005)*scale_)/(2.0*current_))*f;
  std::cout<<"============================"
           <<"\nR_AB,CD = "<<v1
           <<"\nR_BC,AD = "<<v2
           <<"\nVoltage scale = "<<scale_
           <<"\nPI*t/ln2 = "<<PI*t/ln2
           <<"\nResistivity = "<<resistivity
           <<"\nRes Error = "<<res_error
           <<"\nf value = "<<f
           <<"\nR_AB,CD/R_BC,AD = "<<other_f
           <<"\nother stuff: "<<((PI*t)/(ln2))*((v1+v2)/(2.0*current_))
           <<"\nCurrent = "<<current_<<std::endl;
  return resistivity;
}
void VanDerPauw::vanderpauw(void){
  if (voltages_.size() > 4) {
    double resistivity;
    for (unsigned int i=0; i<voltages_.size(); ++i) {
      if (i==0) {
        resistivity = find_resistivity(
                                 voltages_[voltages_.size()-2],voltages_[0]);
      } else if (i==1) {
        resistivity = find_resistivity(
                                 voltages_[1],voltages_[voltages_.size()-1]);
      } else if ((i%2)==0){
        resistivity = find_resistivity(voltages_[i-2],voltages_[i]);
      } else if ((i%2)!=0) {
        resistivity = find_resistivity(
                voltages_[voltages_.size()+2-i],voltages_[voltages_.size()-i]);
      }
      resistivity_.push_back(resistivity);
      //std::cout<<f<<" "<<other_f<< " " <<resistivity<<" "
      //         <<r_val<<" "<<i<<std::endl;
    }
    for (std::vector<double>::const_iterator i=resistivity_.begin();
         i!=resistivity_.end();++i) {
      std::cout << *i << std::endl;
    }
  } else {
    double resistivity;
    for (unsigned int i=0;i<voltages_.size();++i) {
      if (i!=voltages_.size()-1) {
        resistivity = find_resistivity(voltages_[i],voltages_[i+1]);
      } else {
        resistivity = find_resistivity(voltages_[i],voltages_[0]);
      }
      resistivity_.push_back(resistivity);
    }
  }
}
