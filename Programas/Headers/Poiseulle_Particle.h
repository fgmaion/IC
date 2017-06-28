#ifndef PUSED_HPP
#define PUSED_HPP

#include <iostream>
class Particle
{
	public:
	Particle();
	~Particle();

  unsigned int cell;
	unsigned int cell_x;
	unsigned int cell_y;

	long double density;
	long double mass;

	long double x_t[2];
	long double velo_t[2];
	long double a_t[2];
};

Particle::~Particle()
{
}

Particle::Particle()
{
}


class ParticleGrid
{

unsigned int total_size;
unsigned int grid_size;
unsigned int* actual_size;
public:
Particle*** grid;//Eh um array de arrays de particulas
		 											//Guarda os reais tamanhos de cada celula do vetor grid
													//isso pois nao queremos tratar dinamicamente o tamanho de
													//cada vetor guardado na celula
ParticleGrid(unsigned int, unsigned int);
~ParticleGrid();
Particle**& operator[] (unsigned int);
unsigned int push_back(unsigned int, Particle*);
void clear();
unsigned int size(unsigned int);
};

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
unsigned int ParticleGrid::size(unsigned int index)
{
	return(actual_size[index]);
}
/////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
Particle**& ParticleGrid::operator[](unsigned int i)
{
return(grid[i]);
}
//////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
unsigned int ParticleGrid::push_back(unsigned int celula, Particle* p1)
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
actual_size = new unsigned int[grid_size];
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
for (int i = 0; i < grid_size; i++)
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

#endif
