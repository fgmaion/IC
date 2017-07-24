#ifndef POIS_HPP
#define POIS_HPP

#include <cmath>
#include <math.h>
#include </home/francisco/Documents/Fisica/IC-master/Programas/Headers/Poiseulle_Particle.h>
#include </home/francisco/Documents/Fisica/IC-master/Programas/Headers/GaussKernel.hpp>
#include <stdlib.h>

class ParticleSystem
	{
	public:
		//constructor
		ParticleSystem();
		~ParticleSystem();

		long double CalculateDensities(ParticleGrid*, long double);//Funcao calcula as densidades

		void Init(long double, long double, long double, long double,long double,long double, double, long double, long double);

		void Run(long double, long double, long double, long double,long double,long double, double, long double, long double);

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
void ParticleSystem::Init(long double h, long double dt, long double m, long double c, long double mu, long double lambda, double rho_0, long double v_0, long double r_0)
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

for (int i=0; i<810; i++)
	{
			Particle* p1 = &Particles[i];
			//Cell distribution
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
		}
//densidade
for (int i = 0; i<810; i++)
	{
	Particle* p1 = &Particles[i];
	long double A = 0;
	for(int j  = 0; j < PG1->size(p1->cell); j++)
		{
				A = A + m*Gauss_Kernel(p1->x_t[0]-(*PG1)[p1->cell][j]->x_t[0],p1->x_t[1]-(*PG1)[p1->cell][j]->x_t[1],h);
		}
//Boundaries y over density
	if( p1->cell_y <20 )
		{
			p1->density = A / (1 - erf( sqrt( pow( p1->x_t[1] + 1, 2) ) / h )/2 );
		}
	if( p1->cell_y >=0 )
		{
			p1->density = A / (1 - erf( sqrt( pow( p1->x_t[1] - 1, 2) ) / h )/2 );
		}
	}

//accelerations
for(int i = 0; i<810; i++)
{
		Particle* p1 = &Particles[i];
		long double A = 0;
		long double B = 0;
		for(int k = 0; k<360; k++)
			{
			if(k == p1->cell || k == p1->cell - 18 || k == p1->cell+18)
			{
				for(int l = k-1; l<k+2; l++)
				{
				if(l<360 && l>=0)
				{
				for(int j = 0; j < PG1->size(l); j++)
					{
					long double distance = sqrt(pow( (p1->x_t[0]-(*PG1)[l][j]->x_t[0]),2)+pow((p1->x_t[1]-(*PG1)[l][j]->x_t[1]), 2));
					if(p1->density*(*PG1)[l][j]->density*distance  == 0)
					{}
					else{
							long double Xij = p1->x_t[0] - (*PG1)[l][j]->x_t[0];
							long double Yij = p1->x_t[1] - (*PG1)[l][j]->x_t[1];
							long double VXij = p1->velo_t[0] - (*PG1)[l][j]->velo_t[0];
							long double VYij = p1->velo_t[1] - (*PG1)[l][j]->velo_t[1];
							long double T0i_1 = sqrt( pow( p1->x_t[1] + 1, 2) )/h;
							long double T0i_2 = sqrt( pow( p1->x_t[1] - 1, 2) )/h;

							A =
							//Pressure
							A - m*c*c*( ( ( p1->density
							+ (*PG1)[l][j]->density ) /( p1->density*(*PG1)[l][j]->density*distance ) )
							*D_Gauss_Kernel( Xij , Yij ,h )
							* Xij  )
							//Viscosity
						//	OVER HERE THERS SHIT
							- ( mu*m*m/ ( p1->density*(*PG1)[l][j]->density ) )*( p1->velo_t[0]
							*D_Gauss_Kernel( Xij, Yij ,h )
							*( 7/ (3*distance) ) + ( Xij*(Xij*VXij + Yij*VYij)/3 + VXij*( Xij*Xij + Yij*Yij ) )
							*( -D_Gauss_Kernel(Xij, Yij, h ) / ( distance*distance*distance )
							+ D2_Gauss_Kernel( Xij, Yij, h ) / ( distance*distance ) ) );
							//Boundaries y over pressure

							if( p1->cell_y <20 )
								{
								A = A
								//Pressure boundary over y
								+( lambda/2 )*erf( T0i_1 )
								//Viscosity boundary over y
								-( mu*( 2*T0i_1*T0i_1 + 1 ) / ( p1->density*sqrt(pi)*h*h*T0i_1 ) )
								*( p1->velo_t[0]*exp(-T0i_1*T0i_1 ) );
								}
							if( p1->cell_y >=0 )
								{
								A = A
								//Pressure boundary over y
								+( lambda/2 )*erf( T0i_2 )
								//Viscosity boundary over y
								-( mu*( 2*T0i_2*T0i_2 + 1 ) / ( p1->density*sqrt(pi)*h*h*T0i_2 ) )
								*( p1->velo_t[0]*exp(-T0i_2*T0i_2 ) );
								}

							B =
							//pressure
							B -m*c*c*( ( (p1->density
							+(*PG1)[l][j]->density) /(p1->density*(*PG1)[l][j]->density*distance) )
							*D_Gauss_Kernel( Xij , Yij ,h)
							* Yij )
							//Viscosity
							- ( mu*m*m/ ( p1->density*(*PG1)[l][j]->density ) )*( p1->velo_t[1]
							*D_Gauss_Kernel( Xij, Yij ,h )
							*( 7/ (3*distance) ) + ( Yij*(Xij*VXij + Yij*VYij)/3 + VYij*(Xij*Xij + Yij*Yij ) )
							*( -D_Gauss_Kernel(Xij, Yij, h ) / ( distance*distance*distance )
							+ D2_Gauss_Kernel( Xij, Yij, h ) / ( distance*distance ) ) );


							if( p1->cell_y <20 )
								{
								B = B
								//Pressure boundary over y
								- 2*(p1->density*c*c/rho_0) * (p1->density*c*c/rho_0)*exp( - pow ( T0i_1 , 2 ) )/(sqrt(pi)*h)
								//Viscosity Boundary over y
								-( mu*( 2*T0i_1*T0i_1 + 1 ) / ( p1->density*sqrt(pi)*h*h*T0i_1 ) )
								*( p1->velo_t[1] + p1->velo_t[1]/3 )*exp(-T0i_1*T0i_1 ) ;
								}
							if( p1->cell_y >=0 )
								{
								B = B
								//Pressure boundary over y
								+ 2*(p1->density*c*c/rho_0) * (p1->density*c*c/rho_0)*exp( - pow ( T0i_2 , 2 ) )/(sqrt(pi)*h)
								//viscosity boundary over y
								-( mu*( 2*T0i_2*T0i_2 + 1 ) / ( p1->density*sqrt(pi)*h*h*T0i_2 ) )
								*( p1->velo_t[1] + p1->velo_t[1]/3 )*exp(-T0i_2*T0i_2 );
								}

							}
						}
					}
				}
			}
		}
		p1->a_t[0]=(A+2*mu*v_0/(r_0*r_0));
		p1->a_t[1]=(B);
		p1->velo_t[0] = (p1->velo_t[0] + p1->a_t[0]*dt/2);
		p1->velo_t[1] = (p1->velo_t[1] + p1->a_t[1]*dt/2);
}

}
///////////////////////////////////////////////////////////////////////////////////////////
//////////////IMPLEMENTACAO DAS FUNCOES DE DENSIDADE //////////////////////////////////////
/*long double ParticleSystem::CalculateDensities(ParticleGrid* PG1, long double h)
{
for (int i =0; i < 810; i++)
	{
		long double A  = 0;
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
void ParticleSystem::Run(long double h, long double dt, long double m, long double c, long double mu, long double lambda, double rho_0, long double v_0, long double r_0)
{
PG1->clear();

for (int i=0; i<810; i++)
	{
			Particle* p1 = &Particles[i];
//Boundaries x
			if(p1->x_t[0] > 1.8)
			{
				if(p1->velo_t[0] > 0)
				{
					p1->x_t[0] = 0;
				}
			}
			if(p1->x_t[0] < 0)
			{
				if(p1->velo_t[0] < 0)
				{
					p1->x_t[0] = 1.8;
				}
			}
			//Cell distribution
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
		}
//densidade
for (int i = 0; i<810; i++)
	{
	Particle* p1 = &Particles[i];
	long double A = 0;
	for(int j  = 0; j < PG1->size(p1->cell); j++)
		{
				A = A + m*Gauss_Kernel(p1->x_t[0]-(*PG1)[p1->cell][j]->x_t[0],p1->x_t[1]-(*PG1)[p1->cell][j]->x_t[1],h);
		}
//Boundaries y over density
	if( p1->cell_y <20 )
		{
			p1->density = A / (1 - erf( sqrt( pow( p1->x_t[1] + 1, 2) ) / h )/2 );
		}
	if( p1->cell_y >=0 )
		{
			p1->density = A / ( 1 - erf( sqrt( pow( p1->x_t[1] - 1, 2) ) / h )/2 );
		}
}

//accelerations
for(int i = 0; i<810; i++)
	{
		Particle* p1 = &Particles[i];
		long double A = 0;
		long double B = 0;
		for (int k = 0; k<360; k++)
		{
			if(k == p1->cell || k == p1->cell - 18 || k == p1->cell+18)
			{
				for (int l = k-1; l <k+2; l++)
				{
					if(l<360 && l>=0)
					{
					for(int j = 0; j < PG1->size(l); j++)
					{
					long double distance = sqrt( pow( (p1->x_t[0]-(*PG1)[l][j]->x_t[0]),2)+pow((p1->x_t[1]-(*PG1)[l][j]->x_t[1]), 2));
					if(p1->density*(*PG1)[l][j]->density*distance  == 0)
					{}
					else{

							long double Xij = p1->x_t[0] - (*PG1)[l][j]->x_t[0];
							long double Yij = p1->x_t[1] - (*PG1)[l][j]->x_t[1];
							long double VXij = p1->velo_t[0] - (*PG1)[l][j]->velo_t[0];
							long double VYij = p1->velo_t[1] - (*PG1)[l][j]->velo_t[1];
							long double T0i_1 = sqrt( pow( p1->x_t[1] + 1, 2) )/h;
							long double T0i_2 = sqrt( pow( p1->x_t[1] - 1, 2) )/h;

							A =
							//Pressure
							A - m*c*c*( ( ( p1->density
							+ (*PG1)[l][j]->density ) /( p1->density*(*PG1)[l][j]->density*distance ) )
							*D_Gauss_Kernel( Xij , Yij ,h )
							* Xij  )
							//Viscosity
							- ( mu*m*m/ ( p1->density*(*PG1)[l][j]->density ) )*( p1->velo_t[0]
							*D_Gauss_Kernel( Xij, Yij ,h )
							*( 7/ (3*distance) ) + ( Xij*(Xij*VXij + Yij*VYij)/3 + VXij*(Xij*Xij + Yij*Yij ) )
							*( -D_Gauss_Kernel(Xij, Yij, h ) / ( distance*distance*distance )
							+ D2_Gauss_Kernel( Xij, Yij, h ) / ( distance*distance ) ) );
							//Boundaries y over pressure

							if( p1->cell_y <20 )
								{
								A = A
								//Density boundary over y
								+( lambda/2 )*erf( T0i_1 )
								//Viscosity boundary over y
								-( mu*( 2*T0i_1*T0i_1 + 1 ) / ( p1->density*sqrt(pi)*h*h*T0i_1 ) )
								*( p1->velo_t[0]*exp(-T0i_1*T0i_1 ) );
								}
							if( p1->cell_y >=0 )
								{
								A = A
								//Density boundary over y
								+( lambda/2 )*erf( T0i_2 )
								//Viscosity boundary over y
								-( mu*( 2*T0i_2*T0i_2 + 1 ) / ( p1->density*sqrt(pi)*h*h*T0i_2 ) )
								*( p1->velo_t[0]*exp(-T0i_2*T0i_2 ) );
								}

							B =
							//pressure term
							B -m*c*c*( ( (p1->density
							+(*PG1)[l][j]->density) /(p1->density*(*PG1)[l][j]->density*distance) )
							*D_Gauss_Kernel( Xij , Yij ,h)
							* Yij )
							//Viscosity term
							- ( mu*m*m/ ( p1->density*(*PG1)[l][j]->density ) )*( p1->velo_t[1]
							*D_Gauss_Kernel( Xij, Yij ,h )
							*( 7/ (3*distance) ) + ( Yij*(Xij*VXij + Yij*VYij)/3 + VYij*(Xij*Xij + Yij*Yij ) )
							*( -D_Gauss_Kernel(Xij, Yij, h ) / ( distance*distance*distance )
							+ D2_Gauss_Kernel( Xij, Yij, h ) / ( distance*distance ) ) );


							if( p1->cell_y <20 )
								{
								B = B
								//Pressure boundary over y
								- 2*(p1->density*c*c/rho_0) * (p1->density*c*c/rho_0)*exp( - pow ( T0i_1 , 2 ) )/(sqrt(pi)*h)
								//Viscosity Boundary over y
								-( mu*( 2*T0i_1*T0i_1 + 1 ) / ( p1->density*sqrt(pi)*h*h*T0i_1 ) )
								*( p1->velo_t[1] + p1->velo_t[1]/3 )*exp(-T0i_1*T0i_1 ) ;
								}
							if( p1->cell_y >=0 )
								{
								B = B
								//Pressure boundary over y
								+ 2*(p1->density*c*c/rho_0) * (p1->density*c*c/rho_0)*exp( - pow ( T0i_2 , 2 ) )/(sqrt(pi)*h)
								//viscosity boundary over y
								-( mu*( 2*T0i_2*T0i_2 + 1 ) / ( p1->density*sqrt(pi)*h*h*T0i_2 ) )
								*( p1->velo_t[1] + p1->velo_t[1]/3 )*exp(-T0i_2*T0i_2 );
								}
							}
						}
					}
					}
				}
			}
		p1->a_t[0]=(A+2*mu*v_0/(r_0*r_0));
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
