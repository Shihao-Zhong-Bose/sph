#ifndef __grid_h__
#define __grid_h__

#include "ParticleFluid.h"
#include <vector>
#include <list>
using namespace std;
#define MAX_NEIGHBORS 64
class Grid {
public:
	float cellSize;
	int width;
	int height;
	int cellCount;
	vector<FluidParticle*> grid;
public:	
	Grid(float cellSize, int width, int height);
	~Grid();
	void refresh(const ParticleFluid& fluid);
	std::list<FluidParticle*> getNeighbors(FluidParticle& particle);
};
#endif