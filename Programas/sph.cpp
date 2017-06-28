#include <iostream>
#include <cmath>
#include <fstream>
const double pi = 3.14159265358979323846264338328;
///////////////////////////////
//Codigo para aproximacao de uma funcao pelo metodo SPH

using namespace std;

enum CHOICE {Funcao, Derivada};

//kernel parece com funcao delta de dirac
double Kernel(double q, double h)
	{
		double R = 0;
		if (0<=q && q<1)
			{
				R = (pow((2-q),3)-4*pow((1-q),3))/(6*h);
			}
		if ((1<=q && q<2 ))
			{
				R = (pow((2-q),3))/(6*h);
			}
		if (2<=q)
			{
				R = 0;
			}
			return (R);
	}
	
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
	
//funcao a ser aproximada
double func(double y)
	{
		return (sin(y));
	}

double GradFunc(double x)
{
	return cos(x);
}

//funcao Aprox aproxima a quantidade desejada
double AproxFunc(double x, double h)
	{
	cout << "pegando a derivada";
	double A = 0;
	double dx = h*0.1;
	for ( int i =-2/dx; i < (2/dx); i++)
		{
			A = A + func(x+h*i*dx)*Kernel(abs(i*dx), h)*h*dx;
			cout << i << endl;
		}
	return A;
	}
	
double AproxGradFunc(double x, double h)
{
	double A = 0;
	double dx = h*0.1;
	for ( int i =-2/dx; i < (2/dx); i++)
		{
			A = A + func(x+h*i*dx)*GradKernel(-i*dx, h)*h*dx;
		}
	return(A);
}

int main()
{
cout << "Aproximando a funcao seno, valor em pi/6" << AproxFunc(pi/6,0.01) << " valor numerico " << func(pi/6) << endl;
cout << "Aproximando a derivada  " << AproxGradFunc(pi/3,0.01) << " valor numerico " << GradFunc(pi/3);
}


