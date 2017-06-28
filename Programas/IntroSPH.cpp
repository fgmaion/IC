#include <iostream>
#include <vector>
using namespace std;

class Particle{
public:
	Particle();
	float x[2], v[2], a[2];
	float pressure, density;
	float v_prev[2];
	Particle* next;
};

class ParticleSystem{
public:
	ParticleSystem();
	void Step();
	void Draw();
	float Parameters[N_PARAMS];
protected:
	void Boundaries(Particle* a);
	void AddGravity();
	void AddMotionDamping();
	virtual void ComputeAccelerations();
	virtual void InitialPositions();
	std::vector<Particle> Particles;
};

int main()
{
	for (int i; i<Particles.size();i++)
		{
		pcurr = &

