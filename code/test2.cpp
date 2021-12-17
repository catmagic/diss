 
#include <boost/multiprecision/float128.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <iostream> 
#include"quadmath.h"
using namespace boost::multiprecision;

int main()
{
  
std::cout << std::setprecision(200);
long double ld=3.3333333333333333333333333333333333333333333333333333333333333333Q;
float128 londfloat=3.333333333333333333333333333333333333333333333333333333333333Q;
std::cout<<ld<<'\n';
std::cout<<londfloat<<'\n';
std::cout<<(float128)sqrtq((__float128)2.0Q);
return 0;
}
