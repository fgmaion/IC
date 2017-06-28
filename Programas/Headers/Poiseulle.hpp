#ifndef POIS_HPP
#define POIS_HPP

#include <cmath>
#include </home/francisco/Documents/Fisica/Mecflu_IC/Programas/Headers/Poiseulle_Particle.h>
#include </home/francisco/Documents/Fisica/Mecflu_IC/Programas/Headers/GaussKernel.hpp>
#include <stdlib.h>

/*
//PARAMETROS DA SIMULACAO
//lambda = 1
//k = 1/4
double c = 1;
double m =1;
double h = 0.093;
int N_part = 810;
double GRID_SIZE = 10;
*/

class ParticleSystem
	{
	public:
		//constructor
		ParticleSystem();
		~ParticleSystem();

		double CalculateDensities(ParticleGrid*, double);//Funcao calcula as densidades

		void Init(double);

		void Run(double, double, double, double);

		Particle* Particles;

		ParticleGrid* PG1;
};

ParticleSystem::ParticleSystem()
{
}
///////////////////////////////////////////////////////////////////////////////////////////
ParticleSystem::~ParticleSystem()
{
}
///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////
void ParticleSystem::Init(double m)
{
srand(time(NULL));

Particles = new Particle[810];
PG1 = new ParticleGrid(360,810);
for (int i = 0; i<810; i++)
{
			Particle* p1 = &Particles[i];
			p1->mass  = (m);
			p1->x_t[0] = ((rand() % 18000)/10000.0);
			p1->x_t[1] = ((rand() % 20000) /10000.0)-1;
			p1->velo_t[0] = 0;
			p1->velo_t[1] = 0;
			p1->a_t[0] = 0;
			p1->a_t[1] = 0;
}
}
///////////////////////////////////////////////////////////////////////////////////////////
//////////////IMPLEMENTACAO DAS FUNCOES DE DENSIDADE //////////////////////////////////////
/*double ParticleSystem::CalculateDensities(ParticleGrid* PG1, double h)
{
for (int i =0; i < 810; i++)
	{
		double A  = 0;
		Particle* p1 = &Particles[i];
		for (int j = 0; j<PG1->size(i); j++)
		{

			A = A + (*PG1)[i][j]->mass
			*Gauss_Kernel(sqrt(pow(p1->x_t[0]-(*PG1)[i][j]->x_t[0],2)
			+pow(p1->x_t[1]-(*PG1)[i][j]->x_t[1],2)),h);

			p1->density = (A);
			}

	}
}
*/
//////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
void ParticleSystem::Run(double h, double dt, double m, double c)
{
PG1->clear();
//boundaries and cell separation
for (int i=0; i<810; i++)
	{

			//Aqui calculamos (x + 1)/2 * tamanho do grid. Somamos 1 para fazer com que x E [0,2], e dividimos por 2 para
			//ter x E [0,1]. Dessa forma temos qual a porcentagem do grid ocupada pela posicao x. Multiplicamos pelo tamanho do grid
			//e tiramos (int) para obter em qual celula ele esta.
			Particle* p1 = &Particles[i];
//			cout << p1->x_t[0] <<"  " << i <<"\n";
			if(p1->x_t[0] > 1.8)
			{
				if(p1->velo_t[0] > 0)
				{
					p1->x_t[0] = 0;
				}
				else{}
			}
			if(p1->x_t[0] < 0)
			{
				if(p1->velo_t[0] < 0)
				{
					p1->x_t[0] = 1.8;
				}
				else{}
			}



			int xcell = ((p1->x_t[0]) / 1.8)*18;
			int ycell = ((p1->x_t[1] + 1 ) / 2)*20;

			if(xcell < 18 && xcell >= 0)
			{
				p1->cell_x = ( (int)(xcell) );
			}
			else
			{
				if ( xcell > 0 )
				{
				p1->cell_x = 17;
				}
				else
				{
				p1->cell_x = 0;
				}
			}

			if(ycell < 20 && ycell >= 0)
			{
			p1->cell_y = ( (int)(ycell) );
			}
			else
			{
				if ( ycell > 0 )
				{
				p1->cell_y = 19;
				}
				else
				{
				p1->cell_y = 0;
				}
			}
			p1->cell = (p1->cell_x + 18*p1->cell_y);
			PG1->push_back(p1->cell,p1);
			//AQUI TODAS AS PARTICULAS TEM A ELAS UMA CELULA ASSOCIADA
		}
//densidade
for (int i = 0; i<810; i++)
	{
	Particle* p1 = &Particles[i];
	double A = 0;
	//cout << p1->cell <<"  "<< i << "\n";
	for(int j  = 0; j < PG1->size(p1->cell); j++)
		{
			A = A + m*Gauss_Kernel(p1->x_t[0]-(*PG1)[p1->cell][j]->x_t[0],p1->x_t[1]-(*PG1)[p1->cell][j]->x_t[1],h);
			cout << (*PG1)[p1->cell][j]->x_t[0] << "\n";
		}
	p1->density = A;
	//cout <<  << "\n";
}
//accelerations
for(int i = 0; i<810; i++)
	{
		Particle* p1 = &Particles[i];
		double A = 0;
		double B = 0;
		for(int j = 0; j < PG1->size(p1->cell); j++)
			{
			double distance = sqrt(pow( (p1->x_t[0]-(*PG1)[p1->cell][j]->x_t[0] ),2)+pow((p1->x_t[1]-(*PG1)[p1->cell][j]->x_t[1]), 2));

			A = A - m*c*c*((p1->density
			+ (*PG1)[p1->cell][j]->density) /(p1->density*(*PG1)[p1->cell][j]->density*distance)
			*D_Gauss_KernelX(p1->x_t[0]-(*PG1)[p1->cell][j]->x_t[0]
			,p1->x_t[1]-(*PG1)[p1->cell][j]->x_t[1],h)
			*(p1->x_t[0]-(*PG1)[p1->cell][j]->x_t[0]));

			B = B -m*c*c*( (p1->density
			+(*PG1)[p1->cell][j]->density) /(p1->density*(*PG1)[p1->cell][j]->density*distance)
			*D_Gauss_KernelY(p1->x_t[0]-(*PG1)[p1->cell][j]->x_t[0]
			,p1->x_t[1]-(*PG1)[p1->cell][j]->x_t[1],h)
			*(p1->x_t[1]-(*PG1)[p1->cell][j]->x_t[1]));

			}
		p1->a_t[0]=(A+0.5);
		p1->a_t[1]=(B);
		p1->x_t[0]=(p1->x_t[0] + p1->velo_t[0]*dt);
		p1->x_t[1]=(p1->x_t[1] + p1->velo_t[1]*dt);
		p1->velo_t[0] = (p1->velo_t[0] + p1->a_t[0]*dt);
		p1->velo_t[1] = (p1->velo_t[1] + p1->a_t[1]*dt);
	}
}
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////

#endif
