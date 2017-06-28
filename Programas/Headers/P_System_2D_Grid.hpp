#ifndef PGRID_HPP
#define PGRID_HPP

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <time.h>
#include <cmath>
#include <vector>
#include </home/francisco/Documents/Fisica/Mecflu_IC/Programas/Headers/GradKernel2D.hpp>
#include </home/francisco/Documents/Fisica/Mecflu_IC/Programas/Headers/Kernel2D.hpp>
#include </home/francisco/Documents/Fisica/Mecflu_IC/Programas/Headers/Particle_Grid.hpp>


#define NPART 1000
#define GRID_SIZE 10 //NUMERO DE CELULAS DO GRID
#define BOX_SIZE 2 // Tamanho da parede da fronteira
#define H 0.1
#define mass pi/2

using namespace std;

class ParticleSystem
{
	public:
		//constructor
		ParticleSystem();
		~ParticleSystem();

		//Inicializa a posicao, velocidade e aceleracao das particulas
		void initialize(double dt, double h);

		//Faz com que o sistema de um passo temporal
		void Step(double dt, double h);
		void Boundaries();
		Particle* Particles[1000];
		///////////////////////////////////
		void fill_grid(Particle* Particles[1000], double h);
		//////////////////////////////////////
		double CalculateDensity(Particle * p1, double h);
		double CalculateDensityGradX(Particle * p1, double h);
		double CalculateDensityGradY(Particle * p1, double h);
};
ParticleSystem::ParticleSystem()
{
}

//FRONTEIRA CHECA O SENTIDO DA VELOCIDADE E TROCA, SE NECESSARIO
void ParticleSystem::Boundaries()
{
	for (int i =0; i < NPART; i++)
	{
				Particle* p1 = Particles[i];
				if (p1->GetX0() > 1)
					{
						if (p1->GetV0() > 0)
						{
							p1->SetV0(-0.6*p1->GetV0());
						}
					}
				if (p1->GetX0()<-1)
					{
						if (p1->GetV0() < 0)
						{
							p1->SetV0(-0.6*p1->GetV0());
						}
					}
				if (p1->GetX1() > 1)
					{
						if (p1->GetV1()>0)
						{
							p1->SetV1(-0.6*p1->GetV1());
						}
					}
				if (p1->GetX1() < -1)
					{
						if (p1->GetV1()<0)
						{
							p1->SetV1(-0.6*p1->GetV1());
						}
					}
	}
}

 ParticleSystem::fill_grid(Particle* Particles[NPART], double h)
{
double xcell;
double ycell;
for (int i=0; i<NPART; i++)
{
Particle* p1 = Particles[i];
//Aqui calculamos (x + 1)/2 * tamanho do grid. Somamos 1 para fazer com que x E [0,2], e dividimos por 2 para
//ter x E [0,1]. Dessa forma temos qual a porcentagem do grid ocupada pela posicao x. Multiplicamos pelo tamanho do grid
//e tiramos (int) para obter em qual celula ele esta.
xcell = ((p1->GetX0() + 1 ) / 2)*GRID_SIZE;
ycell = ((p1->GetX1() + 1 ) / 2 )*GRID_SIZE;
if(xcell < GRID_SIZE && xcell > 0
	&& ycell <GRID_SIZE && ycell >0)
{
p1->cell_x = (int)(xcell);
p1->cell_y = (int)(ycell);
grid[p1->cell_x + 3*p1->cell_y].push_back(p1);
}
else
{
p1->cell_x = GRID_SIZE;
p1->cell_y = GRID_SIZE;
grid[p1->cell_x + 3*p1->cell_y].push_back(p1);//A maneira pela qual as celulas foram numeradas esta descrita no caderno vermelho,
//na entrada de 10/05/2017
}
//AQUI TODAS AS PARTICULAS TEM A ELAS UMA CELULA ASSOCIADA
}
}

double ParticleSystem::CalculateDensity(Particle* p1, double h)
{
double densidade = 0;
int p1_cell = p1->cell_x+3*p1->cell_y;
for (int i = 0; i<grid[p1_cell].size(); i++)
{
	double d = sqrt(pow(grid[p1_cell][i]->GetX0()-p1->GetX0(),2)+pow(grid[p1_cell][i]->GetX1()-p1->GetX1(),2));
		if( d < 2*h)
			{
				densidade = densidade + mass*Kernel2D(d,h);
			}
}
return densidade;
}

double ParticleSystem::CalculateDensityGradX(Particle* p1,double h)
{
	double densidade = 0;
	int p1_cell = p1->cell_x+3*p1->cell_y;
	for (int i = 0; i<grid[p1_cell].size(); i++)
	{
			double d = sqrt(pow(grid[p1_cell][i]->GetX0()-p1->GetX0(),2)+pow(grid[p1_cell][i]->GetX1()-p1->GetX1(),2));
			if( d < 2*h)
				{
					densidade = densidade + mass*GradKernel2DX(grid[p1_cell][i]->GetX0(),grid[p1_cell][i]->GetX1(),h);
				}
	}
	return densidade;
}

double ParticleSystem::CalculateDensityGradY(Particle* p1,double h)
{
	double densidade = 0;
	int p1_cell = p1->cell_x+3*p1->cell_y;
	for (int i = 0; i<grid[p1_cell].size(); i++)
	{
			double d = sqrt(pow(grid[p1_cell][i]->GetX0()-p1->GetX0(),2)+pow(grid[p1_cell][i]->GetX1()-p1->GetX1(),2));
			if( d < 2*h)
				{
					densidade = densidade + mass*GradKernel2DY(grid[p1_cell][i]->GetX0(),grid[p1_cell][i]->GetX1(),h);
				}
	}
	return densidade;
}

void ParticleSystem::initialize(double dt, double h)
{
	srand(time(NULL));
	double m = 1/1000.0;
	for (int i = 0; i<1000; i++)
	{
				Particle* p1 = Particles[i];

				fill_grid(Particle[i],h);

				p1->SetMass(m);
				cout << "rodando 3.3\n";
				p1->SetX0((rand() % 20000)/10000.0 - 1);
				cout << "rodando 3.4\n";
				p1->SetX1((rand() % 20000)/10000.0 - 1);

				p1->SetA0(-(0.5)*CalculateDensityGradX(p1,h) -p1->GetX0());
				cout << "rodando 3.\n";
				p1->SetA1(-(0.5)*CalculateDensityGradY(p1,h) -p1->GetX1());

				p1->SetV0(p1->GetA0()*(dt/2));
				p1->SetV1(p1->GetA1()*(dt/2));
	}
}

void ParticleSystem::Step(double dt, double h)
{
	for (int i = 0; i< 1000; i++)
		{
				Particle* p1 = Particles[i];

				double X0_prev = p1->GetX0();
				double X1_prev = p1->GetX1();

				p1->SetX0(p1->GetX0() + p1->GetV0()*dt);
				p1->SetX1(p1->GetX1() + p1->GetV1()*dt);

				p1->SetA0(-(0.5)*CalculateDensityGradX(p1,h) -X0_prev -0.5*(p1->GetV0()));
				p1->SetA1(-(0.5)*CalculateDensityGradY(p1,h) -X1_prev -0.5*(p1->GetV1()));

				p1->SetV0(p1->GetV0() + p1->GetA0()*dt);
				p1->SetV1(p1->GetV1() + p1->GetA1()*dt);

		}
}
#endif
