#include <iostream>
#include <fstream>
using namespace std;

//LEAP FROG METHOD

class Particle
	{
	public:
	double x_t;
	double velo_t;
	double a_t;
	double mass; 
		
	//Funcao que atualiza a posição, velocidade e aceleração da particula no tempo
	void walk(float dt, float k_stiff);
	//dá à partícula os parametros iniciais
	void initialize(float x, float v, float a, float dt, float m);
	
	Particle();
	~Particle();
	};
	
	Particle::Particle()
	{
	}
	
	void Particle::initialize(float x, float v, float a, float dt, float m)
	{
		a_t = a;
		x_t = x;
		//passo atrás de meio dt, faz parte do método Leap Frog.
		velo_t = v-(a*dt/2);
		mass = m;
	}
	
	void Particle::walk( float dt, float k_stiff)
	{
	double x_prev = x_t;
	
	x_t = (x_t + velo_t*dt);
	
	a_t = (-k_stiff*x_prev/mass);
	
	velo_t = velo_t + a_t*dt;
	
	}
	
int main()
{

double T;
double dt = 0.01;

double Position [500];
double Time [500];

cout << "Voce deseja saber a posicao da particula em que tempo? (No maximo 5. Duas casas decimais.)\n";
cin >> T;

//Cria a "particula"
Particle *P1 = new Particle();
P1->initialize(1,0.5,0,dt, 0.5);

// Loop de iteracao

for (int i = 0 ; i<(T/dt); i++)
	{
		//Faz com que a particula ande
		P1->walk(dt, 100);
		cout << P1->x_t << "\n";
		//Coloca em arrays a posicao e o tempo
		Position [i] = P1->x_t;	
		Time [i] = dt*i;
	}
//Imprime os arrays num arquivo de texto
ofstream myfile ("/home/francisco/myfile.txt");	
if (myfile.is_open())
{

	for (int count = 0; count < (T/dt); count++)
	{
		myfile << Time[count]<< " " << Position[count] << " \n";
	}
myfile.close();
}
}
