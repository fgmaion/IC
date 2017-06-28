#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <time.h>
#include <cmath>
#include </home/francisco/Documents/Fisica/Mecflu_IC/Programas/Headers/GradKernel2D.hpp>
#include </home/francisco/Documents/Fisica/Mecflu_IC/Programas/Headers/Kernel2D.hpp>
#include </home/francisco/Documents/Fisica/Mecflu_IC/Programas/Headers/Particle2D.hpp>

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
		///////////////////////////////////
		double CalculateDensity(Particle * p1, double h);
		double CalculateDensityGradX(Particle * p1, double h);
		double CalculateDensityGradY(Particle * p1, double h);
		Particle Particles[1000];
};
ParticleSystem::ParticleSystem()
{
}	

//FRONTEIRA CHECA O SENTIDO DA VELOCIDADE E TROCA, SE NECESSARIO
void ParticleSystem::Boundaries()
{
	for (int i =0; i < 1000; i++)
	{
				Particle* p1 = &Particles[i];
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

//INICIALIZA ALEATORIAMENTE A POSICAO DAS PARTICULAS NUM QUADRADO DE LADO 2		
void ParticleSystem::initialize(double dt, double h)
{
	srand(time(NULL));
	double m = 1/1000.0;
	for (int i = 0; i<1000; i++)
	{
				Particle* p1 = &Particles[i];
				
				p1->SetMass(m);
				
				p1->SetX0((rand() % 20000)/10000.0 - 1);
				p1->SetX1((rand() % 20000)/10000.0 - 1);
				
				p1->SetA0(-(0.5)*CalculateDensityGradX(p1,h) -p1->GetX0());
				p1->SetA1(-(0.5)*CalculateDensityGradY(p1,h) -p1->GetX1());
				
				p1->SetV0(p1->GetA0()*(dt/2));
				p1->SetV1(p1->GetA1()*(dt/2));
	}
}

double ParticleSystem::CalculateDensity(Particle* p1, double h)
{
double A = 0;
for(int i =0; i<1000; i++)
{
	if (Particles[i].GetX0() > (p1->GetX0()-2*h) && Particles[i].GetX0() < (p1->GetX0()+2*h))
			{
				if (Particles[i].GetX1()> p1->GetX1()-2*h && Particles[i].GetX1() < p1->GetX1()+2*h)
				{
						A = A + Particles[i].GetMass()*Kernel2D(sqrt(pow((p1->GetX0()-Particles[i].GetX0()),2)+pow((p1->GetX1()-Particles[i].GetX1()),2))/h,h);
				}
			}
}
return A;
}

double ParticleSystem::CalculateDensityGradX(Particle* p1,double h)
{
double A =0;
for (int i =0; i <1000; i++)
{
	if (Particles[i].GetX0() > (p1->GetX0()-2*h) && Particles[i].GetX0() < (p1->GetX0()+2*h))
			{
				if (Particles[i].GetX1()> p1->GetX1()-2*h && Particles[i].GetX1() < p1->GetX1()+2*h)
				{
					A = A + Particles[i].GetMass()*GradKernel2DX(p1->GetX0()-Particles[i].GetX0(),p1->GetX1()-Particles[i].GetX1(),h);
				}
			}
}
return A;
}

double ParticleSystem::CalculateDensityGradY(Particle* p1,double h)
{
double A =0;
for (int i =0; i <1000; i++)
{
	if (Particles[i].GetX0() > (p1->GetX0()-2*h) && Particles[i].GetX0() < (p1->GetX0()+2*h))
			{
				if (Particles[i].GetX1()> p1->GetX1()-2*h && Particles[i].GetX1() < p1->GetX1()+2*h)
				{
					A = A + Particles[i].GetMass()*GradKernel2DY(p1->GetX0()-Particles[i].GetX0(),p1->GetX1()-Particles[i].GetX1(),h);
					
				}
			}
}
return A;
}

void ParticleSystem::Step(double dt, double h)
{
	for (int i = 0; i< 1000; i++)
		{		
				Particle* p1 = &Particles[i];
				
				double X0_prev = p1->GetX0();
				double X1_prev = p1->GetX1();
				
				p1->SetX0(p1->GetX0() + p1->GetV0()*dt);
				p1->SetX1(p1->GetX1() + p1->GetV1()*dt);
				
				p1->SetA0(-(0.5)*CalculateDensityGradX(p1,h) -X0_prev -0.5*(p1->GetV0()));
				p1->SetA1(-(0.5)*CalculateDensityGradY(p1,h) -X1_prev -0.5*(p1->GetV1()));
				
				p1->SetV0(p1->GetV0() + p1->GetA0()*dt);
				p1->SetV1(p1->GetV1() + p1->GetA1()*dt);
				
				//ITERACAO TEMPORAL NO METODO LEAPFROG FEITO EM UMA UNICA ETAPA
			}
}


int main()
{
double h = 0.1;

clock_t tStart = clock();
ParticleSystem * PS1 = new ParticleSystem;
PS1->initialize(0.004, h);
for (int i = 0; i<5000; i++)
	{
		
		PS1->Step(0.004,h);
		PS1->Boundaries();
	}
printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
	


//EXPORTA OS DADOS PARA UM ARQUIVO
ofstream myfile2 ("/home/francisco/myfile2.txt");
if (myfile2.is_open())
{
	for (int j = 0; j<1000; j++)
		{
		myfile2 << PS1->Particles[j].GetX0()<< "," << PS1->Particles[j].GetX1() <<"\n";
		}
myfile2.close();
}

ofstream myfile3 ("/home/francisco/myfile3.txt");
if (myfile3.is_open())
{
	for (int k = 0; k<1000; k++)
	{
		double r = sqrt(pow(PS1->Particles[k].GetX0(),2)+pow(PS1->Particles[k].GetX1(),2));
		myfile3 << r << "," << PS1->CalculateDensity(&PS1->Particles[k],h) << "\n"; 
	}
}

}
