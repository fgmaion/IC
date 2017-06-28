#ifndef GAUSS_HPP
#define GAUSS_HPP

#include <iostream>
#include <cmath>
using namespace std;

double pi =3.141592653589793238462643383279;

double Gauss_Kernel(double x, double y,double h)
{
double d  = sqrt(pow(x,2)+pow(y,2));
double A = exp(-d*d/(h*h))/(pi*h*h);
return (A);
}

double D_Gauss_Kernel(double x, double y, double h)
{
  double d  = sqrt(pow(x,2)+pow(y,2));
  double A = -( 2*d/ ( h*h ) ) * exp(-d*d / ( h*h ) )/( pi*h*h );
  return(A);
}

double D2_Gauss_Kernel(double x, double y, double h)
{
  double d = sqrt(pow(x,2)+pow(y,2));
  double A = -( 2/( pi * h*h*h*h ) ) * ( exp( -d*d/( h*h ) ) - (2 * d*d/( h*h ) )* exp( -d*d/ ( h*h ) ) );
  return(A);
}
#endif
