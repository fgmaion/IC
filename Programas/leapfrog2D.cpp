#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

//leap frog method

class Particle
	{
	public:
	float x_t[2];
	float velo_t[2];
	float a_t[2];
	float mass; 
	
	void walk(float dt);
	void initialize(float x[], float v[], float a[], float dt, float m);
	
	Particle();
	~Particle();
	};
	
	Particle::Particle()
	{
	}
	
	void Particle::initialize(float x[], float v[], float a[], float dt, float m)
	{
		cout << "inicializa";
		a_t[0] = a[0];
		a_t[1] = a[1];

		x_t[0] = x[0];
		cout << x_t[0] << "Valor de c iniciado";
		x_t[1] = x[1];
		
		velo_t[0] = v[0]-(a[0]*dt/2);
		velo_t[1] = v[1]-(a[1]*dt/2);
		
		mass = m;
	}
	void Particle::walk( float dt)
	{
	//Raio da trajetoria 
	double R = 2;
	//calcula a
	a_t[0] = -(pow(velo_t[0],2)+pow(velo_t[1],2))*(x_t[0]/sqrt(pow(x_t[0],2)+pow(x_t[1],2)))/R;
	a_t[1] = -(pow(velo_t[0],2)+pow(velo_t[1],2))*(x_t[1]/sqrt(pow(x_t[0],2)+pow(x_t[1],2)))/R;
	
	//Calcula v
	velo_t[0] = (velo_t[0] + a_t[0]*dt);
	velo_t[1] = (velo_t[1] + a_t[1]*dt);
	
	//Calcula a nova posicao
	x_t[0] = (x_t[0] + velo_t[0]*dt);
	x_t[1] = (x_t[1] + velo_t[1]*dt);
	
	}
	
int main()
{

double T;
double dt = 0.01;

double Position1 [5000];
double Position2 [5000];
//double Time [500];

cout << "Voce deseja saber a posicao da particula em que tempo? (No maximo 5. Duas casas decimais.)\n";
cin >> T;

//Cria a "particula"
Particle *P1 = new Particle();
float arr_x[2] = {2,0};
float arr_v[2] = {0,1};
float arr_a[2] = {(-1/2),0};
P1->initialize(arr_x,arr_v,arr_a,dt, 0.5);

// Loop de iteracao

for (int i; i<(T/dt); i++)
	{
		//Faz com que a particula ande
		P1->walk(dt);
		cout  << P1->x_t[1] <<" " << P1->x_t[2] << "\n";
		//Coloca em arrays a posicao e o tempo
		Position1 [i] = P1->x_t[0];
		Position2 [i] = P1->x_t[1];
		//Time [i] = dt*i;
	}
//Imprime os arrays num arquivo de texto
ofstream myfile ("/home/francisco/myfile.txt");	
if (myfile.is_open())
{

	for (int count = 0; count < (T/dt); count++)
	{
		myfile << Position1 [count]<< " " << Position2[count] << " \n";
	}
myfile.close();
}
}
