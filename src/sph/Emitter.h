#ifndef __emitter_h__
#define __emitter_h__

#include <glm/glm.hpp>
#include <vector>
#include "ParticleFluid.h"

using namespace std;

class Emitter {
public:
	glm::vec2 pos;
	glm::vec2 direction;
	float speed;
	float delayMS;
	int particleToEmit;
	int particleEmitted;
	float size;
	int delay;
	float particlesPerSec;
public:
	Emitter(glm::vec2 pos, glm::vec2 direction, float particlesPerSecond, float particlesToEmit, float size);
	vector<FluidParticle*> emitParticles(float dt);
};

#endif

