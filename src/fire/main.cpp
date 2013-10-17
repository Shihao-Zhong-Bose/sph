#define GLEW_STATIC
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <iostream>
#include <iomanip>
#include <vector>

using namespace std;

struct FluidDesc {
	FluidDesc(string name, float densityAtRest, float dt) : density(densityAtRest), timestep(dt) {};
	float density;
	float timestep;
};

//----- constants
const int screenWidth = 800;
const int screenHeight = 600;
const float viewWidth = 10.f;
const float viewHeight = (screenHeight * viewWidth) / screenWidth;
const float particleRadius = 0.05f;

const float PI = 3.1415926f;
const float epsilon = 1e-8;
const int maxParticles = 5000;

const float ParticlesPerSecond = 2000;
const float tpp = 1.f/(float)ParticlesPerSecond;



//----- variables
GLFWwindow* window;



//----- prototypes
void initWindow();
void update(float dt);
void render();
void createParticles(int numParticles);


//----- functions

struct Particle {
	Particle() {
		memset(pos, 0, sizeof(pos));
		memset(v, 0, sizeof(v));
		lifespan = 0;
		noBurst = true;
	}
	bool noBurst;
	float pos[2];
	float v[2];
	float lifespan;
	glm::vec3 col;
};


Particle* particles[maxParticles];
int particleCount = 0;

void createParticles(int numParticles) {
	for(int i(0); i < numParticles; ++i) {
		Particle* p = new Particle();
		p->lifespan = 2 + ((float)rand()/RAND_MAX * 3);
		p->pos[0] = 5;
		p->pos[1] = 0;
		p->v[0] = 0;
		p->v[1] = 5;

		particles[particleCount++] = p;
	}
}

void init() {
	initWindow();

	//createParticles(10);


}

void initWindow() {
	if (!glfwInit()) {
		std::cout << "Failed to init GLFW\n";
		return;
	}

	window = glfwCreateWindow(1280, 800, "Particles", NULL, NULL);
	if (!window) {
		glfwTerminate();	
		std::cout << "Failed to init GLFW\n";
		return;
	}

	glfwMakeContextCurrent(window);

	glewExperimental = GL_TRUE;
	GLenum err = glewInit();
	if(!err == GL_FALSE) {
		std::cout << "GLEW init failed (Error:" << glewGetErrorString(err) << ")" << std::endl;
		return;
	}

	glClearColor(0.f, 0.f, 0.f, 1.f);
}

void emitParticles(float dt) {
	static float acc = 0;
	acc += dt;
	for(; acc / tpp > 1.f; acc -= tpp) {
		Particle* p = new Particle();
		float r = -2 +  (float)rand()/RAND_MAX * 4;
		p->noBurst = true;
		p->lifespan = 0.5 + ((float)rand()/RAND_MAX * 1.25);
		p->pos[0] = 5 + r;
		p->pos[1] = 0;
		p->v[0] = -3 + 4 * (float)rand()/RAND_MAX;
		p->v[1] = 1 + 3 * (float)rand()/RAND_MAX;

		particles[particleCount++] = p;
	}

	// 	int particlesToEmit = (int)(((float)ParticlesPerSecond * dt) + 0.5f);
	// 
	// 	for(int i = 0; i < particlesToEmit && (i - particleCount) < particlesToEmit && i + particleCount < maxParticles; ++i) {
	// 		Particle* p = new Particle();
	// 		float r = 1 * (float)rand()/RAND_MAX;
	// 		p->noBurst = false;
	// 		p->lifespan = 0.5 + ((float)rand()/RAND_MAX * 1);
	// 		p->pos[0] = 5 + r;
	// 		p->pos[1] = 0;
	// 		p->v[0] = -2 + 4 * (float)rand()/RAND_MAX;
	// 		p->v[1] = 5+ 5 * (float)rand()/RAND_MAX;
	// 
	// 		particles[particleCount++] = p;
	// 	}
}

void burst(Particle* particle) {
	int num = 5 + rand() % 30;
	for(int i(0); i < num; ++i) {
		Particle* p = new Particle();


		p->lifespan = 0.25 +  0.5 * (float)rand()/RAND_MAX;
		p->pos[0] = particle->pos[0];
		p->pos[1] = particle->pos[1];
		p->col.b = 0.5;
		p->col.r = -0.5;
		float angle = ((float)rand())/RAND_MAX*PI*2.f;

		float x = cos(angle) * (float)rand()/RAND_MAX * 10;
		float y = sin(angle) * (float)rand()/RAND_MAX * 10;


		p->v[0] = x;
		p->v[1] = y;

		particles[particleCount++] = p;
	}
}

struct GravityWell {
	GravityWell(float x, float y, float s, float r) {
		pos[0] = x;
		pos[1] = y;
		strength = s;
		rad = r;
	}
	float rad;
	float pos[2];
	float strength;
};

GravityWell well(viewWidth/2.f, viewHeight, 10, 5);

void update(float dt) {
	for(int i(0); i < particleCount; ++i) {
		Particle* p = particles[i];

		p->lifespan -= dt;
		if(p->lifespan <= 0) {
			if(p->noBurst == false) {
				burst(p);
			}
			particles[i] = 0;
			particles[i] = particles[particleCount - 1];
			particles[particleCount - 1] = NULL;
			particleCount--;
			i--;
			delete p;
			continue;
		}

		p->v[0] -= 0;
		p->v[1] -= 0 * dt;
		p->pos[0] += p->v[0] * dt + (-5 + 10* sin((float)rand() / RAND_MAX * PI)) * dt;
		p->pos[1] += p->v[1] * dt;

		float dx = well.pos[0] + (-0.25 + (float)rand()/RAND_MAX * 0.5f) - p->pos[0];
		float dy = well.pos[1] - p->pos[1];
		float r2 = dx * dx + dy * dy;
		//if(r2 < (well.rad * well.rad)) {
			float vdir[2];
			float d = sqrt(r2);
			vdir[0] = dx / d;
			vdir[1] = dy / d;
			p->pos[0] += well.strength * vdir[0] * 1.f/d * dt;
			p->pos[1] += well.strength * vdir[1] * 1.f/d * dt;
		//}

		if(p->pos[1] > viewHeight) {
		//	p->pos[1] = viewHeight;
		//	p->v[1] *= -0.5f;
		} else if(p->pos[1] < 0) {
			p->pos[1] = 0;
			p->v[1] *= -0.5f;
		}

		if(p->pos[0] > viewWidth) {
			p->pos[0] = viewWidth;
			p->v[0] *= -0.5f;
		} else if(p->pos[0] < 0) {
			p->pos[0] = 0;
			p->v[0] *= -0.5f;
		}

	}
}

void render() {
	glClearColor(0.f, 0.f, 0.f, 1);
	glClear(GL_COLOR_BUFFER_BIT);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(0, viewWidth, 0, viewHeight, 0, 1);

	glEnable(GL_POINT_SMOOTH);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);


	vector<glm::vec2> positions;
	vector<glm::vec4> colors;
	for(int i = 0; i < particleCount; ++i) {
		Particle* pi = particles[i];
		positions.push_back(glm::vec2(pi->pos[0], pi->pos[1]));
		glm::vec2 sd = glm::normalize(glm::vec2(pi->v[0], pi->v[1]));



		float l = pi->lifespan > 1 ? 1.f : pi->lifespan;
		float alpha = 1.f;
		
		alpha = pi->lifespan;
		
		colors.push_back(glm::vec4(1.f, pi->lifespan * pi->lifespan , 0.f, alpha));

		
	}

	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_COLOR_ARRAY);

	glPointSize(2.f*particleRadius*screenWidth/viewWidth);

	if(positions.size() > 0) {
		glColorPointer(4, GL_FLOAT, sizeof(glm::vec4), &colors[0]);
		glVertexPointer(2, GL_FLOAT, sizeof(glm::vec2), &positions[0]);
		glDrawArrays(GL_POINTS, 0, positions.size());
	}

	glDisableClientState(GL_COLOR_ARRAY);
	glDisableClientState(GL_VERTEX_ARRAY);


}

#include <sstream>
#include <string>

int main(int argc, char** argv) {
	init();
	float updateTime = 0;
	float renderTime = 0;
	float dt = 0.0;
	int frameCount = 0;
	float startFrame = 0;
	float endFrame = 0;
	float accumulate = 0;
	float frameTime = 0.1;
	while (!glfwWindowShouldClose(window)) {
		if(glfwGetKey(window, GLFW_KEY_ESCAPE)) break;

		startFrame = glfwGetTime();

		emitParticles(frameTime);

		emitParticles(frameTime);
		update(frameTime);


		glClear(GL_COLOR_BUFFER_BIT);


		render();

		glfwSwapBuffers(window);
		glfwPollEvents();
		endFrame = glfwGetTime();
		frameTime = endFrame - startFrame;
		accumulate += frameTime;
		frameCount++;
		if(accumulate >= 1.f) {
			ostringstream os;
			os << "FPS: " << frameCount << " Particles: " << particleCount;
			string s = os.str();
			glfwSetWindowTitle(window, s.c_str());
			accumulate = 0;
			frameCount = 0;
		}
	}

	glfwTerminate();

	return 0;
}