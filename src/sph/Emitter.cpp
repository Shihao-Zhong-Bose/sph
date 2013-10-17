#include "Emitter.h"

Emitter::Emitter(glm::vec2 pos, glm::vec2 direction, float particlesPerSecond, float particlesToEmit, float size) {
	//this->pos = pos;
	//this->direction = glm::normalize(direction);
	//this->speed = speed;
	//this->delayMS = delayMS;
	this->particleToEmit = particlesToEmit;
	this->particlesPerSec = particlesPerSecond;
	this->size = size;
	this->pos = pos;
	this->direction = direction;
	particleEmitted = 0;
}

vector<FluidParticle*> Emitter::emitParticles(float dt) {
	vector<FluidParticle*> particles;
	
	if(particleEmitted >= particleToEmit) return particles;

	int pc = (int)(particlesPerSec * dt) + 0.5f;

	for(int i = 0; i < pc; ++i) {
		float r = size * ((float)rand()/RAND_MAX - 0.5f);

		FluidParticle* p = new FluidParticle();
		p->pos[0] = pos.x;
		p->pos[1] = pos.y + r;
		p->prevPos[0] = p->pos[0];
		p->prevPos[1] = p->pos[1];
		p->mass = 0.5f;
		p->v[0] = (float)rand()/RAND_MAX * 10 * direction.x;
		p->v[1] = (float)rand()/RAND_MAX * 5 * direction.y;
		particles.push_back(p);
	}

	particleEmitted += pc;

	return particles;
}