#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <time.h>
#include <cmath>
#include </home/francisco/Documents/Fisica/IC-master/Programas/Headers/GaussKernel.hpp>
#include </home/francisco/Documents/Fisica/IC-master/Programas/Headers/Poiseulle_Particle.h>
#include </home/francisco/Documents/Fisica/IC-master/Programas/Headers/Novo.hpp>

using namespace std;

int main()
{
double r_0 = 1;
double rho_0 = 10.0; //colocamos uma densidade alta para observar menos turbulencia
double v_0 = 0.1; //nao sei o que significa isso. Suponho que seja uma velocidade caracteristica
//a qual queremos que seja baixa para minimizar turbulencia.
long double m = 0.00444*rho_0*r_0*r_0;
//long double h = 1.4*sqrt(m/rho_0;
long double h = 0.1;
long double mu = 0.2*r_0*v_0*rho_0;
double c = 2*v_0;
long double lambda = -2*mu*v_0/(rho_0*r_0*r_0);

double dt = 0.001;

double xcell;
double ycell;
int running;

double Reynolds = rho_0*v_0*r_0/(mu);
cout << " h " << h << "\n";
cout << "Reynolds Number " << Reynolds;
clock_t tStart = clock();
////////////////////////////////////////////////////////////////////////////////////////////
ParticleSystem* PS1 = new ParticleSystem;
PS1->Init(h,dt,m,c,mu,lambda,rho_0);
ofstream myfile10;
{
myfile10.open("/home/francisco/Documents/Fisica/Mecflu_IC/Programas/Imagens/Ibagens2/myfile10.txt");
for(int i = 0; i <1000; i++)
	{
	PS1->Run(h,dt,m,c,mu,lambda,rho_0);
	for (int l = 0; l<810; l++)
		{
			myfile10 << PS1->Particles[l].x_t[0] << ";" << PS1->Particles[l].x_t[1] <<  "\n";
		}
		myfile10 << "split\n";
	}
	myfile10.close();
}
ofstream myfile11;
{
myfile11.open("/home/francisco/myfile11.csv");
for(int i = 0; i<810; i++)
{
	myfile11 << PS1->Particles[i].x_t[1] << "," <<  PS1->Particles[i].velo_t[0] << "\n";
}
myfile11.close();
}
//////////////////////////////PARTE DE EXPORTAR OS DADOS EM TXT/////////////////////////////
printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);

/*
//EXPORTA OS DADOS PARA UM ARQUIVO
ofstream myfile2 ("/home/francisco/myfile4.txt");
if (myfile2.is_open())
{
	for (int j = 0; j<1000; j++)
		{
		myfile2 << PS1->Particles[j]->GetX0()<< "," << PS1->Particles[j]->GetX1() << "," << PS1->CalculateDensity(PS1->Particles[j], h) << "\n";
		}
myfile2.close();
}

ofstream myfile3 ("/home/francisco/myfile5.txt");
if (myfile3.is_open())
{
	for (int k = 0; k<1000; k++)
	{
		double r = sqrt(pow(PS1->Particles[k]->GetX0(),2)+pow(PS1->Particles[k]->GetX1(),2));
		myfile3 << r << "," << PS1->CalculateDensity(PS1->Particles[k],h) << "\n";
	}
}
*/

}
