#include <iostream>
#include <cmath>
const double pi = 3.14159265358979323846264338328;

using namespace std;

double GradKernel2DX (double x, double y, double h)
{
	double R = 0;
	double q = sqrt(pow(x,2)+pow(y,2));
	if (q<1 && q>=0)
		{
			R = ((3*x/q)*pow((2-q),2)-(12*x/q)*pow(1-q,2))*(15/(14*pi*pow(h,2)));
		}
	if ( q>=1 && q<=2)
		{
			R = (3*x/q)*pow(2-q,2)*(15/(14*pi*pow(h,2)));
		}
	else
		{
			R = 0;
		}
	return (R);
}

double GradKernel2DY (double x, double y, double h)
{
	double R = 0;
	double q = sqrt(pow(x,2)+pow(y,2));
	if (q<1 && q>=0)
		{
			R = ((3*y/q)*pow((2-q),2)-(12*y/q)*pow(1-q,2))*(15/(14*pi*pow(h,2)));
		}
	if ( q>=1 && q<=2)
		{
			R = (3*y/q)*pow(2-q,2)*(15/(14*pi*pow(h,2)));
		}
	else
		{
			R = 0;
		}
	return (R);
}


double Kernel(double x, double y, double h)
	{	
		double q = sqrt(pow(x,2)+pow(y,2));
		double R = 0;
		if (0<=q && q<1)
			{
				R = (pow((2-q),3)-4*pow((1-q),3))*(15/(14*pi*pow(h,2)));
			}
		if ((1<=q && q<2 ))
			{
				R = (pow((2-q),3))*(15/(14*pi*pow(h,2)));
			}
		if (2<=q)
			{
				R = 0;
			}
		return (R);
	}
