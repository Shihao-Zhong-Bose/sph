#include "Grid.h"

Grid::Grid(float cellSize, int width, int height) {
	this->cellSize = cellSize;
	this->width = width;
	this->height = height;
	this->cellCount = width * height;
 	grid.resize(cellCount);
 	std::fill(grid.begin(), grid.end(), nullptr);
}

void Grid::refresh(const ParticleFluid& fluid) {

	// Clear grid
	std::fill(grid.begin(), grid.end(), nullptr);

	vector<FluidParticle*> particles = fluid.getParticles();
	// Add particles to grid
	for (size_t i=0; i<particles.size(); ++i)
	{
		FluidParticle& pi = *particles[i];
		int x = pi.pos[0] / cellSize;
		int y = pi.pos[1] / cellSize;


		if (x < 1)	x = 1;
		else if (x > width - 2) x = width - 2;

		if (y < 1) y = 1;
		else if (y > height - 2) y = height - 2;

		pi.next = grid[x + (y * width)];
		grid[x + (y * width)] = &pi;

		pi.gridX = x;
		pi.gridY = y;
	}

// 	std::fill(_grid.begin(), _grid.end(), nullptr);
// 
// 	vector<FluidParticle*> particles = fluid.getParticles();
// 	for(int i = 0; i < particles.size(); ++i) {
// 		FluidParticle* pi = particles[i];
// 		int x = (int)glm::clamp((pi->pos.x / _cellSize), 1.f, (float)(_width - 2));
// 		int y = (int)glm::clamp((pi->pos.y / _cellSize), 1.f, (float)(_height - 2)); 
// 		
// 		pi->next = _grid[x + y * _width];
// 		_grid[x + y * _width] = pi;
// 		pi->gridX = x;
// 		pi->gridY = y;
// 	}
}

std::list<FluidParticle*> Grid::getNeighbors(FluidParticle& particle) {
	//slow as shit


 	std::list<FluidParticle*> neighbors;
// 	for(int i = -1; i <= 2; ++i) {
// 		if(particle.gridX + i >= _width) continue;
// 		for(int j = -1; j <= 2; ++j) {
// 			if(particle.gridY + j >= _height) continue;
// 			int index = particle.gridX + i + (particle.gridY + j) * _width;
// 			FluidParticle* p = _grid[index];
// 
// 			while(p) {
// 				if(neighbors.size() >= MAX_NEIGHBORS) break;
// 				float r  = glm::length(p->pos - particle.pos);
// 				if(r < 1e-8 || r > _cellSize) { 
// 					p = p->next;
// 					continue;
// 				}
// 				neighbors.push_back(p);
// 				p = p->next;
// 			}
// 			
// 		}
// 	}
 	return neighbors;
}