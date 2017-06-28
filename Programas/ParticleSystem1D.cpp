#include <iostream>
#include <cmath>
#include <fstream>
#include </home/francisco/Documents/Física/Mecflu_IC/Programas/Headers/Kernel.hpp>
#include </home/francisco/Documents/Física/Mecflu_IC/Programas/Headers/GradKernel.hpp>
using namespace std;

class Particle
	{
	public:
	float x_t;
	float velo_t;
	float a_t;
	float mass;

	Particle();
	};
	
	Particle::Particle()
	{
	}
//////////////////////////////Classe ParticleSystem///////////////////////////////////////////////////////
class ParticleSystem
	{
	public:
		//constructor
		ParticleSystem();
		
		//Inicializa a posicao, velocidade e aceleracao das particulas
		void initialize(float dt);
		
		//Faz com que o sistema de um passo temporal
		void Step(float dt);
		void Boundaries();
		///////////////////////////////////
		double CalculateDensity(Particle * p1, double h);
		double CalculateDensityGrad(Particle * p1, double h);
		Particle Particles[200];
	};
	ParticleSystem::ParticleSystem()
{ 
}

	void ParticleSystem::Boundaries()
	{
		for ( int i =0; i < 200; i++)
		{
			Particle* p1 = &Particles[i];
			if (p1->x_t > 1)
				{
					p1->x_t = p1->x_t -2;
				}
			if ( p1->x_t < -1)
				{
					p1->x_t = p1->x_t +2;
				}
		}
	}
	void ParticleSystem::initialize(float dt)
	{	
		double GridGauge = 2/200.0;
		for (int i =0 ; i < 200 ; i++)
			{	
				Particle* p1 = &Particles[i];
				
				p1->a_t = 0;
				p1->x_t = (i-100)*GridGauge;
				p1->velo_t = 0;
				p1->mass = 4/600.0;
			}
	}
	
double ParticleSystem::CalculateDensityGrad(Particle * p1, double h)
	{
		double A = 0;
		for (unsigned int i =0; i<200; i++)
		{
				if ( abs(p1->x_t - Particles[i].x_t) < 2*h )
				{
					A = A + Particles[i].mass*GradKernel((p1->x_t - Particles[i].x_t)/h, h);
				}			
		}
		return A;

	}

double ParticleSystem::CalculateDensity(Particle * p1, double h)
{
	double A = 0;
	for ( int j = 0; j < 200; j++)
		{	
			if ( abs(p1->x_t-Particles[j].x_t) < 2*h )
			{
				A = A + Particles[j].mass*Kernel(abs(p1->x_t-Particles[j].x_t)/h, h);
			}
		}
	return A;
}

void ParticleSystem::Step( float dt)
{
for (unsigned short int i =0; i< 200; i++)
	{	
		Particle* p1 = &Particles[i];
		
		double A = p1->a_t;
		p1->a_t = -(0.5)*CalculateDensityGrad(&Particles[i],0.12)- p1->x_t-(0.02*p1->velo_t);
		//Calcula v
		p1->velo_t = p1->velo_t + 0.5*(A + p1->a_t)*dt;
		//Calcula a nova posicao
		p1->x_t = p1->x_t + p1->velo_t * dt + 0.5*A*pow(dt,2);
	}
}
int main()
{
ParticleSystem* PS1 = new ParticleSystem;
PS1->initialize(0.1);
for (int j =0; j<50000; j++)
{
PS1->Step(0.01);
PS1->Boundaries();
}


ofstream myfile ("/home/francisco/myfile.txt");	

if (myfile.is_open())
{
	for (int count = 0; count < (200); count++)
	{	
		myfile << PS1->Particles[count].x_t << "," << PS1->CalculateDensity(&PS1->Particles[count],0.12) << "\n";
	}
myfile.close();
}

}
