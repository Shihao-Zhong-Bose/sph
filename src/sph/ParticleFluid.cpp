#include "ParticleFluid.h"



ParticleFluid::ParticleFluid(int maxParticles) 
	: _maxParticles(maxParticles) {
	
	stiffness = 0.08;;
	nearStiffness = 0.1;
	linearViscosity = 0.5f;
	quadraticViscosity = 1.f;
	damping = 0.1f;
	densityAtRest = 20.f;
	surfaceTension = 0.0004f;
	particleRadius = 0.05f;

	float dt = 0.1;
}

void ParticleFluid::addParticles(vector<FluidParticle*> particles) {
	int remainingParticles = _maxParticles - _particles.size();
	if(remainingParticles <= 0) return;

	if(remainingParticles >= particles.size()) {
		_particles.insert(_particles.end(), particles.begin(), particles.end()); 
	} else {
		_particles.insert(_particles.end(), particles.begin(), particles.begin() + remainingParticles); 
	}

}

const vector<FluidParticle*>& ParticleFluid::getParticles() const {
	return _particles;
}

// float ParticleFluid::getStiffness() const {
// 	return stiffness;
// }
// 
// float ParticleFluid::getViscosity() const {
// 	return viscosity;
// }
