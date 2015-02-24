#ifndef __particle_fluid_h__
#define __particle_fluid_h__

#include <vector>
#include <glm/glm.hpp>

using namespace std;

struct FluidParticle {
	FluidParticle() {
		density = 0;
		nearDensity = 0;
		mass = 0;
		pressure = 0;
		nearPressure = 0;
		memset(pos, 0, sizeof(pos));
		memset(prevPos, 0, sizeof(prevPos));
		memset(relaxedPos, 0, sizeof(relaxedPos));
		memset(v, 0, sizeof(v));
	}
	float pos[2];
	
	float prevPos[2];
	
	float relaxedPos[2];
	float v[2];


	float pressure;
	float nearPressure;

	float density;
	float nearDensity;
	
	float mass;

	FluidParticle* next;
	int gridX;
	int gridY;
};

class ParticleFluid {
private:
	vector<FluidParticle*> _particles;
	const int _maxParticles;
public:
	float particleRadius;
	float stiffness;
	float nearStiffness;
	float linearViscosity;
	float quadraticViscosity;
	float particleMass;
	float damping;
	float densityAtRest;
	float surfaceTension;
	
public:
	ParticleFluid(int maxParticles);
	void addParticles(vector<FluidParticle*> particles);
    void doit(glm::vec3 velo);

	const vector<FluidParticle*>& getParticles() const;

// 	float getStiffness() const;
// 	float getViscosity() const;
// 	float getParticleMass() const;
// 	float getDamping() const;
// 	float getDensityAtRest() const;
// 
// 	void setStiffness(float stiffness);
// 	void setViscosity(float viscosity);
// 	void setParticleMass(float particleMass);
// 	void setDamping(float damping);
// 	void setDensityAtRest(float densityAtRest);
};	

#endif