# include <iostream>
#include <cmath>
const double pi = 3.14159265358979323846264338328;
using namespace std;

double GradKernel(double q, double h)
	{
		double R = 0;
		if(0<=q && q<1)
		{
			R = (-3*pow((2-q),2)+12*pow((1-q),2))/(6*pow(h,2));
		}
		if(-1<q && q<=0)
		{
			R = (3*pow((2+q),2)-12*pow((1+q),2))/(6*pow(h,2));
		}
		if (1<=q && q<2)
		{
			R = -3*pow((2-q),2)/(6*pow(h,2));
		}
		if (-2 < q && q<=-1)
		{
			R = 3*pow((2+q),2)/(6*pow(h,2));
		}
		if (q<=-2)
		{
			R=0;
		}
		if (2<=q)
		{
			R=0;
		}
	return (R);
	}
double func(double y)
{
	return sin(y);
}
double GradFunc(double x)
{
	return cos(x);
}
double Aprox(double x, double h)
	{
	double A = 0;
	double dx = h*0.1;
	for ( int i =-2/dx; i < (2/dx); i++)
		{
			cout << A << "  " << i << "  " << func(x+h*i*dx) << "  \n";
			A = A + func(x+h*i*dx)*GradKernel(-i*dx, h)*h*dx;
		}
	return A;
	}
	
int main()
{
cout << "Aproximando a derivada da funcao seno" << endl;
cout << " A aproximacao de cos(pi/4) resultou  " << Aprox(pi/6,0.01) << "\n";
cout << "O valor numerico e de " << GradFunc(pi/6);
}


