 
#include <boost/multiprecision/float128.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <iostream> 
//#include<quadmath>
using namespace boost::multiprecision;

int main()
{
  

float128 **A_ij,**Betta;
A_ij=(float128**)malloc(5*sizeof(float128*));
Betta=(float128**)malloc(5*sizeof(float128*));
for(int i=0;i<5;i++)
{
    A_ij[i]=(float128*)malloc(5*sizeof(float128));
    Betta[i]=(float128*)malloc(5*sizeof(float128));
}
for(int i=0;i<5;i++)
{
    for(int j=0;j<5;j++)
    {
        Betta[i][j]=0.0;
    }
}
for(int i=0;i<5;i++)
{
    for(int j=0;j<5;j++)
    {
        A_ij[i][j]=Betta[i][j];
    }
}
for(int i=0;i<5;i++)
{
    free(A_ij[i]);
    free(Betta[i]);
}
free(A_ij);
free(Betta);
return 0;
}
