#ifndef PUSED_HPP
#define PUSED_HPP

#include <iostream>
#include <vector>
class Particle
{
	public:

  const int GetCellx();
  void SetCellx(int x);

  const int GetCelly();
  void SetCelly(int y);

	const double GetMass();
	void SetMass(double m);

	const double GetX0();
	const double GetX1();

	const double GetV0();
	const double GetV1();

	const double GetA0();
	const double GetA1();

	void SetX0(double x0);
	void SetX1(double x1);

	void SetV0(double v0);
	void SetV1(double v1);

	void SetA0(double a0);
	void SetA1(double a1);

	void SetDensity(double d);
	const double GetDensity();

	void SetCell(int);
	const int GetCell();

	Particle();
	~Particle();

	private:
  unsigned int cell;
	unsigned int cellx;
	unsigned int celly;

	double density;
	double mass;

	double x_t[2];
	double velo_t[2];
	double a_t[2];
};


void Particle::SetCellx(int c)
{
cellx = c;
}

void Particle::SetCell(int c)
{
cell = c;
}

const int Particle::GetCell()
{
	return cell;
}

void Particle::SetX0(double x0)
{
		x_t[0] = x0;
}

void Particle::SetX1(double x1)
{
		x_t[1] = x1;
}

void Particle::SetV0(double v0)
{
	velo_t[0] = v0;
}

void Particle::SetDensity(double d)
{
	density = d;
}

const double Particle::GetX0()
{
	return x_t[0];
}

const double Particle::GetX1()
{
	return x_t[1];
}

const double Particle::GetV0()
{
	return velo_t[0];
}

const double Particle::GetV1()
{
	return velo_t[1];
}

const double Particle::GetA0()
{
	return a_t[0];
}

const double Particle::GetA1()
{
	return a_t[1];
}

void Particle::SetMass(double m)
{
	mass = m;
}
const double Particle::GetMass()
{
	return mass;
}

const double Particle::GetDensity()
{
	return density;
}

Particle::~Particle()
{
}

Particle::Particle()
{
}
#endif //TESTE_HPP


class ParticleGrid
{
unsigned int total_size;
unsigned int grid_size;
unsigned int* actual_size;
Particle*** grid;//Eh um array de arrays de particulas
 											//Guarda os reais tamanhos de cada celula do vetor grid
											//isso pois nao queremos tratar dinamicamente o tamanho de
											//cada vetor guardado na celula

public:
ParticleGrid(unsigned int, unsigned int);
~ParticleGrid();
Particle**& operator[] (unsigned int);
int push_back(unsigned int, Particle*);
void clear();
};
/////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
Particle**& ParticleGrid::operator[](unsigned int i)
{
return(grid[i]);
}
//////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
int ParticleGrid::push_back(unsigned int celula, Particle* p1)
{
if(actual_size[celula]< total_size)
	{
		grid[celula][actual_size[celula]] =	p1;
		actual_size[celula]++;
		return (0);
	}
return(1);
}
///////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
ParticleGrid::ParticleGrid(unsigned int size_to_allocate, unsigned int N_particles)
{
total_size = N_particles;
grid_size = size_to_allocate;
grid = new Particle**[grid_size];
for(int i = 0; i<grid_size; i++)
	{
		grid[i] = new Particle*[total_size]; //cria-se um vetor de particulas do tamanho do numero
																							// de particulas no sistema
		for(int j = 0; j<total_size; j++)
			{
				grid[i][j] = NULL;
			}
		actual_size[i] = 0;
	}
}
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
ParticleGrid::~ParticleGrid()
{
for (i = 0; i < grid_size; i++)
	{
	delete [] grid[i];
	}
delete [] grid;
}
/////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
void ParticleGrid::clear()
{
	for(int i =0; i<grid_size; i++)
	{
		actual_size[i]=0;
	}
}
