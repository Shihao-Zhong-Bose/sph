#include "SPHSimulation.h"

SPHSimualation::SPHSimualation(float cellSize, int width, int height, glm::vec2 domain) {
	_grid = new Grid(cellSize, width, height);
	_viewWidth = domain.x;
	_viewHeight = domain.y;
	_cellSize = cellSize;

	//not sure where these kernels came from, but blog post mixed kernels from different papers to get these
	// they work, so whatever
	_kernel = (20.f/(2.f*3.1415926f*_cellSize*_cellSize));
	_nearKernel = (30.f/(2.f*3.1415926f*_cellSize*_cellSize));

	_walls[0].norm = glm::vec2(1, 0);
	_walls[0].center = 0;

	_walls[1].norm = glm::vec2(0, 1);
	_walls[1].center = 0;

	_walls[2].norm = glm::vec2(0, -1);
	_walls[2].center = -_viewHeight;

	_walls[3].norm = glm::vec2(-1, 0);
	_walls[3].center = -_viewWidth;

}

void SPHSimualation::simulate(ParticleFluid& fluid, float dt) {
	for(int i(0); i < fluid.getParticles().size(); ++i) {
		fluid.getParticles()[i]->v[1] += -9.8f * dt;
	}
	advance(fluid, dt);
	_grid->refresh(fluid);
	calculatePressure(fluid);
	calculateRelaxedPositions(fluid, dt);
	moveToRelaxedPositions(fluid, dt);
	resolveCollisions(fluid, dt);
}

void SPHSimualation::advance(ParticleFluid& fluid, float dt) {
	//move particle
	vector<FluidParticle*> particles = fluid.getParticles();
	for(int i = 0; i < particles.size(); ++i) {
		FluidParticle* p = particles[i];
		p->prevPos[0] = p->pos[0];
		p->prevPos[1] = p->pos[1];
		p->pos[0] += dt * p->v[0];		
		p->pos[1] += dt * p->v[1];
	}
}

void SPHSimualation::resolveCollisions(ParticleFluid& fluid, float dt) {
	vector<FluidParticle*> particles = fluid.getParticles();
	float particleRad = fluid.particleRadius;
	for(int i = 0; i < particles.size(); ++i) {
		FluidParticle* p = particles[i];

// 		if(p->pos[1] > _viewHeight) {
// 				p->pos[1] = _viewHeight;
// 				//p->v[1] *= -0.9f;
// 		} else if(p->pos[1] < 0) {
// 			//p->pos[1] = 0;
// 			p->v[1] *= -.9f;
// 		}
// 
// 		if(p->pos[0] > _viewWidth) {
// 			//p->pos[0] = _viewWidth;
// 			p->v[0] *= -0.9f;
// 		} else if(p->pos[0] < 0) {
// 			//p->pos[0] = 0;
// 			p->v[0] *= -0.9f;
// 		}
		
		//this way makes stuff look better

		// if we are hitting a wall, add some velocity in the direction away from the wall
		for (int j = 0; j < 4; ++j) {
			Boundary& wall = _walls[j];
			float dis = wall.norm.x * p->pos[0] + wall.norm.y * p->pos[1] - wall.center;
			if (dis < particleRad) {	
				if (dis < 0) {
					dis = 0;
				}
				p->v[0] += 0.8 * (particleRad - dis) * wall.norm.x / dt;
				p->v[1] += 0.8 * (particleRad - dis) * wall.norm.y / dt;
			}
		}
	}
}

void SPHSimualation::calculatePressure( ParticleFluid& fluid ) {
	vector<FluidParticle*> particles = fluid.getParticles();
	float density, nearDensity;
	for(int i = 0; i < particles.size(); ++i) {
		FluidParticle* pi = particles[i];
		int xIndex = pi->gridX;
		int yIndex = pi->gridY * _grid->width;

		neighbors[i].count = 0;

		float density = 0;
		float nearDensity = 0;

		// so my version of getting neighbors is just way too slow. I found this other....ugly way but is much faster

// 		std::list<FluidParticle*> neighbors = _grid->getNeighbors(*pi);
//  		for(auto it = neighbors.bexIndexn(); it != neighbors.end(); ++it) {
// 		 	FluidParticle& pj = *(*it);
		for (int j = (xIndex - 1); j <= (xIndex + 1); ++j) {
			for (int k = yIndex - _grid->width; k <= yIndex + _grid->width; k += _grid->width) {
				for (FluidParticle* pj = _grid->grid[j+k]; pj != NULL; pj = pj->next) {
					
					//TIL, using glm::vec2 or pow function in here significantly slows down your frame.
					//glm::vec2 dist = pj.pos - pi->pos;
					float dx = pj->pos[0] - pi->pos[0];
					float dy = pj->pos[1] - pi->pos[1];
					float r2 = dx*dx + dy*dy;
					//float r2 = dist.x * dist.x + dist.y * dist.y;
					if (r2 < 1e-8 || r2 > (_cellSize * _cellSize))
						continue;
					
					float r = sqrt(r2);
					float a = 1.f - (r/_cellSize); //kernel from vicoelastic paper
					density += pj->mass * a*a * _kernel;
					nearDensity += pj->mass * a*a*a * _nearKernel;

					//save neighbors for quick access
					if (neighbors[i].count < MAX_NEIGHBORS) {
						neighbors[i].particles[neighbors[i].count] = pj;
						neighbors[i].r[neighbors[i].count] = r;
						neighbors[i].count++;
					}
				}
			}
		}

		pi->density = density;
		pi->nearDensity = nearDensity;
		pi->pressure = fluid.stiffness * (pi->density - fluid.densityAtRest);
		pi->nearPressure = fluid.nearStiffness * pi->nearDensity;
	}
// 	vector<FluidParticle*> particles = fluid.getParticles();
// 	float density, nearDensity;
// 	for(int i = 0; i < particles.size(); ++i) {
// 		FluidParticle* pi = particles[i];
// 		density = 0;
// 		nearDensity = 0;
// 		std::list<FluidParticle*> neighbors = _grid->getNeighbors(*pi);
// 		for(auto it = neighbors.begin(); it != neighbors.end(); ++it) {
// 			FluidParticle* pj = *it;
// 			
// 			glm::vec2 dist = pj->pos - pi->pos;
// 			float len = glm::length(dist);
// 			if(len < 0.f || len > _cellSize) continue;
// 
// 			float a = 1.f - (len / _cellSize);
// 			density += pj->mass * pow(a, 3) * _norm;
// 			nearDensity += pj->mass * pow(a, 3) * _nearNorm;
// 		}
// 
// 		pi->density = density;
// 		pi->nearDensity = nearDensity;
// 		pi->pressure = fluid.stiffness * (pi->density - fluid.densityAtRest);
// 		pi->nearPressure = fluid.nearStiffness * pi->nearDensity;
// 	}

}

void SPHSimualation::calculateRelaxedPositions( ParticleFluid& fluid, float dt) {
	vector<FluidParticle*> particles = fluid.getParticles();
	for(int i = 0; i < particles.size(); ++i) {
		FluidParticle* pi = particles[i];

		float x = pi->pos[0];
		float y = pi->pos[1];

		// so instead of calculating forces, as i failed to do, we calculate relaxed positions based on description in Viscoelastic paper. 
		//That wasnt working, but then I found a blog post that did things a little different and mixed methods from different papers. That worked, so I get to not fully understand where the maths came from.	 
		for (int j = 0; j < neighbors[i].count; ++j) {
			FluidParticle* pj = neighbors[i].particles[j];
			float r = neighbors[i].r[j];
			
			float dx = pj->pos[0] - pi->pos[0];
			float dy = pj->pos[1] - pi->pos[1];

			float a = 1.f - (r/_cellSize);

			//(6)
			float d = 0.5f * dt * dt * 
				((pi->pressure + pj->nearPressure) * a *a * _nearKernel + 
				(pi->pressure + pj->pressure) * a * _kernel);
			
			// get new position based off of pressure
			x -= (d * dx) / (r*pi->mass);
			y -= (d * dy) / (r*pi->mass);

			// update position due to surface tension
			if (pi->mass == pj->mass) {
				float a1 = (fluid.surfaceTension) * a * a * _kernel;
				x += a1 * dx;
				y += a1 * dy;
			}

			// update position due to viscosity impulses
			if(a  < 1) {
				//inward radial velocity
				float du = pi->v[0] - pj->v[0];
				float dv = pi->v[1] - pj->v[1];
				float u = (du * dx) + (dv * dy);
				if (u > 0) {
					//use linear and quadratic impulses
					u /= r;
					float I = 0.5f * dt * a * 
						((fluid.linearViscosity * u) + (fluid.quadraticViscosity * u * u));
					x -= I * dx * dt;
					y -= I * dy * dt;
				}
			}

		}
		
		//set new position
		pi->relaxedPos[0] = x;
		pi->relaxedPos[1] = y;
	}
}

void SPHSimualation::moveToRelaxedPositions( ParticleFluid& fluid, float dt ) {
	vector<FluidParticle*> particles = fluid.getParticles();
	for(int i = 0; i < particles.size(); ++i) {
		FluidParticle* pi = particles[i];
		pi->pos[0] = pi->relaxedPos[0];
		pi->pos[1] = pi->relaxedPos[1];
		pi->v[0] = (pi->pos[0] - pi->prevPos[0]) / dt;
		pi->v[1] = (pi->pos[1] - pi->prevPos[1]) / dt;
	}
}




//// OLD STUFF =========================================

/*
glm::vec2 SPHSimulation::pressureKernelGradient( glm::vec2 r ) {
	float len = glm::length(r);
	float a = -45.f / (Constants::PI * pow(_smoothingKernel, 6));
	float b = pow(_smoothingKernel - len, 2);
	return r * (1.f / len) * a * b;
}

float SPHSimulation::viscosityKernelLaplacian( glm::vec2 r ) {
	float len = glm::length(r);
	float a = 45.f / (Constants::PI * pow(_smoothingKernel, 6));
	float b = (len - len);
	return a*b;
}

glm::vec2 SPHSimulation::kernelGradient( glm::vec2 r ) {
	float a = -945.f / (32.f * Constants::PI * pow(_smoothingKernel, 9));
	float b = pow(pow(_smoothingKernel, 2) - pow(glm::length(r), 2), 2);
	return r * a * b;

}

float SPHSimulation::kernelLaplacian( glm::vec2 r ) {
	float lenSqrd = pow(glm::length(r), 2);
	float a = -945.f / (32.f * Constants::PI * pow(_smoothingKernel, 9));
	float b = pow(_smoothingKernel, 2) - lenSqrd;
	float c = 3.f * pow(_smoothingKernel, 2) - 7.f * lenSqrd;
	return a * b * c;
}


void SPHSimulation::update(float dt) {
	std::cout << "ComputeDensityAndPressure()\n";
	computeDensityAndPressure();
	std::cout << "ComputeInternalForces()\n";
	computeInternalForces();
	std::cout << "ComputeExternalForces()\n";
	computeExternalForces();
	std::cout << "UpdatingParticles()\n";
	updateParticles(dt);
	std::cout << "Done\n";
	//_grid->refresh(_particles);
}

float SPHSimulation::kernel( glm::vec2 r ) {
	float len = glm::length(r);
	if(len > _smoothingKernel) {
		return 0;
	}

	float a = 315.0/(64.0*Constants::PI*pow(_smoothingKernel,9));
	float b = pow(pow(_smoothingKernel, 2) - pow(len,2),3);
	return a*b;   
}

void SPHSimulation::computeDensityAndPressure() {
	for(int i(0); i < _particles.size(); ++i) {
		FluidParticle* pi = _particles[i];

		pi->density = 0.f;
		std::set<FluidParticle*> neighbors = _grid->getNeighborhood(pi->pos);
		for(auto it = neighbors.begin(); it != neighbors.end(); ++it) {
			FluidParticle* pj = *it;
			if(pi == pj) continue;
			pi->density += _particleMass * kernel(pi->pos - pj->pos);
		}

		pi->pressure = _gasStiffness * (pi->density - _densityAtRest);
	}
}

void SPHSimulation::computeInternalForces() {
	for(int i(0); i < _particles.size(); ++i) {
		FluidParticle* pi = _particles[i];
		pi->fpressure = glm::vec2();
		pi->fviscosity = glm::vec2();

		std::set<FluidParticle*> neighbors = _grid->getNeighborhood(pi->pos);
		for(auto it = neighbors.begin(); it != neighbors.end(); ++it) {
			FluidParticle* pj = *it;
			if(pi == pj) continue;

			pi->fpressure += pressureKernelGradient(pi->pos - pj-> pos) 
				* _particleMass 
				* ((pi->pressure / pow(pi->density, 2)) + (pj->pressure / pow(pj->density, 2)));
			pi->fviscosity += ((pj->velocity - pi->velocity) 
				* (_particleMass / pj->density)) 
				* viscosityKernelLaplacian(pi->pos - pj->pos);
		}
		pi->fpressure *= -1.f * pi->density;
		pi->fviscosity *= _viscosity;
		pi->finternal = pi->fpressure + pi->fviscosity;
	}
}

void SPHSimulation::computeExternalForces() {
	for(int i(0); i < _particles.size(); i++) {
		FluidParticle *pi = _particles[i];
		pi->fgravity = pi->density * Constants::GRAVITY;


		pi->surfaceNormal = glm::vec2();
		pi->color = 0.0;
		std::set<FluidParticle*> neighbors = _grid->getNeighborhood(pi->pos);
		for(auto it = neighbors.begin(); it != neighbors.end(); ++it) {
			FluidParticle *pj = *it;
			pi->surfaceNormal += kernelGradient(pi->pos - pj->pos) * (_particleMass / pj->density);
			pi->color += kernel(pi->pos - pj->pos) * (_particleMass / pj->density); 
		}

		pi->fsurface = glm::vec2();
		if(glm::length(pi->surfaceNormal) >= sqrt(_densityAtRest / _avgKernelParticles)) {
			float tmp = 0.f;
			for(auto it = neighbors.begin(); it != neighbors.end(); ++it) {
				FluidParticle *pj = *it;
				if(pi == pj) continue;
				tmp += (_particleMass / pj->density) * kernelLaplacian(pi->pos - pj->pos);
			}
			glm::normalize(pi->surfaceNormal);
			pi->fsurface = pi->surfaceNormal * (-_surfaceTension * tmp);
		}

		pi->fexternal = pi->fgravity + pi->fsurface;
	}
}

void SPHSimulation::updateParticles(float dt) {
	for(int i(0); i < _particles.size(); i++) {
		FluidParticle* pi = _particles[i];
		pi->ftotal = pi->finternal + pi->fexternal;

		pi->accel = pi->ftotal * (1.f / pi->density);
		pi->velocityOld = pi->velocity;
		pi->velocity += dt * pi->accel;
		pi->pos += dt * pi->velocity;


		//Collision Detection
		if(pi->pos.y < -0.9) {
			pi->velocity = pi->velocity - glm::vec2(0.f, 1.f) * (1.f+0.f*(-0.2f-pi->pos.y)/(dt*glm::length(pi->velocity)))*glm::dot(pi->velocity,glm::vec2(0.f, 1.f));
			pi->pos.y = -0.9;     	
		} else if(pi->pos.y > 0.9) {
			pi->pos.y = 0.9;
		}

		pi->velocity = (pi->velocityOld + pi->velocity)*(1.0f/2.0f);
	}   
}
*/