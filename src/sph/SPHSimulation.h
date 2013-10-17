#ifndef __sphsimulation_h__
#define __sphsimulation_h__

#include "ParticleFluid.h"
#include "Grid.h"

struct Boundary {
	glm::vec2 norm;
	float center;
};

struct Neighbors
{
	FluidParticle* particles[64];
	float r[64];
	size_t count;
};

class SPHSimualation {
private:
	Grid* _grid;
	float _cellSize;
	float _viewWidth;
	float _viewHeight;
	float _kernel;
	float _nearKernel;
	Boundary _walls[4];
	Neighbors neighbors[5000];
public:
	SPHSimualation(float cellSize, int width, int height, glm::vec2 domain);
	void simulate(ParticleFluid& fluid, float dt);
private:
	void advance(ParticleFluid& fluid, float dt);	
	void resolveCollisions(ParticleFluid& fluid, float dt);	
	void calculatePressure(ParticleFluid& fluid);
	void calculateRelaxedPositions(ParticleFluid& fluid, float dt);
	void moveToRelaxedPositions(ParticleFluid& fluid, float dt);

	
private:
	//old stuff

	/*
	void computeDensityAndPressure();
	void computeInternalForces();
	void computeExternalForces();
	void updateParticles(float dt);
	float kernel(glm::vec2 r);
	glm::vec2 kernelGradient(glm::vec2 r);
	float kernelLaplacian(glm::vec2 r);
	glm::vec2 pressureKernelGradient(glm::vec2 r);
	float viscosityKernelLaplacian(glm::vec2 r);
	*/
};

#endif