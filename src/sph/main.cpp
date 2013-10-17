#define GLEW_STATIC
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <iostream>
#include <iomanip>
#include <vector>
#include "ParticleFluid.h"
#include "SPHSimulation.h"
#include "Emitter.h"
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

// play with these
const int usePreset = 1; // set to -1 for none;
const bool startWithParticles = true;
const bool useEmitters = true;
const float densityAtRest = 40.f; //82.f
const float stiffness = 0.08f;
const float nearStiffness = 0.01f;
const float surfaceTension = 0.0004;
float timestep = 0.005; 

const float linearViscosity = 0.5f;
const float quadViscosity = 1.f;
const float particleRadius = 0.05f;
//-----

FluidDesc fluidDescs[2] = {
	FluidDesc("good fluid", 40, 0.005),
	FluidDesc("lava lamp", 100, 0.002)
};

const float PI = 3.1415926f;
const float epsilon = 1e-8;
const int maxParticles = 5000;
const float h = 6 * particleRadius;
const float cellSize = h;
const int gridWidth = (int)(viewWidth / cellSize);
const int gridHeight = (int)(viewHeight / cellSize);
const int gridSize = gridWidth * gridHeight;


//----- variables
GLFWwindow* window;
SPHSimualation* sph;
ParticleFluid* fluid;
Emitter* emitter;
Emitter* emitter2;

//----- prototypes
void initWindow();
void update(float dt);
void render();
void createParticles(int numParticles);


//----- functions

void createParticles(int numParticles) {
	std::vector<FluidParticle*> particles;
	for(int i(0); i < numParticles; ++i) {
		FluidParticle* p = new FluidParticle();
		p->pos[0] = (float)rand()/RAND_MAX * 5;
		p->pos[1] = (float)rand()/RAND_MAX * 5;
		p->prevPos[0] = p->pos[0];
		p->prevPos[1] = p->pos[1];
		p->mass = 0.5f;
		p->v[0] = (float)rand()/RAND_MAX * 10;
		p->v[1] = (float)rand()/RAND_MAX * 10;
		particles.push_back(p);
	}
	fluid->addParticles(particles);
}

void init() {
	initWindow();

	//create fluid
	fluid = new ParticleFluid(5000);
	fluid->stiffness = stiffness;
	fluid->surfaceTension = surfaceTension;
	fluid->densityAtRest = densityAtRest;
	fluid->particleRadius = particleRadius;
	fluid->linearViscosity = linearViscosity;
	fluid->quadraticViscosity = quadViscosity;
	fluid->nearStiffness = nearStiffness;

	if(usePreset != -1) {
		fluid->densityAtRest = fluidDescs[usePreset].density;
		timestep = fluidDescs[usePreset].timestep;
	}

	sph = new SPHSimualation(cellSize, gridWidth, gridHeight, glm::vec2(viewWidth, viewHeight));
	emitter = new Emitter(glm::vec2(0.f, 7.f), glm::vec2(1.f, -1.f), 100, 500, .75);
	emitter2 = new Emitter(glm::vec2(10.f, 7.f), glm::vec2(-1.f, -1.f), 25, 500, 0.1);

	if(startWithParticles) {
		createParticles(1500);
	}
	
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
	if(!useEmitters) return;
	// wait 10 seconds before emitting
	static float accumulate = 0;
	accumulate += dt;

	if(accumulate > 10) {
		fluid->addParticles(emitter->emitParticles(dt));
		
	}
	if(accumulate > 12) {
		fluid->addParticles(emitter2->emitParticles(dt));
	}
}

void update(float dt) {
	
	sph->simulate(*fluid, dt);
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

	vector<FluidParticle*> particles = fluid->getParticles();
	vector<glm::vec2> positions;
	vector<glm::vec4> colors;
	for(int i = 0; i < particles.size(); ++i) {
		FluidParticle* pi = particles[i];
		positions.push_back(glm::vec2(pi->pos[0], pi->pos[1]));
		glm::vec2 sd = glm::normalize(glm::vec2(pi->v[0], pi->v[1]));
		colors.push_back(glm::vec4(sd, 1.f, 0.8f));
	}

	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_COLOR_ARRAY);

	glPointSize(2.f*particleRadius*screenWidth/viewWidth);

	glColorPointer(4, GL_FLOAT, sizeof(glm::vec4), &colors[0]);
	glVertexPointer(2, GL_FLOAT, sizeof(glm::vec2), &positions[0]);
	glDrawArrays(GL_POINTS, 0, particles.size());

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
	float frameTime = 0;
	while (!glfwWindowShouldClose(window)) {
		if(glfwGetKey(window, GLFW_KEY_ESCAPE)) break;
		
		startFrame = glfwGetTime();

		emitParticles(frameTime);

		//multiple steps per frame
		dt = glfwGetTime();
		update(timestep);
		update(timestep);
		update(timestep);
		updateTime = glfwGetTime() - dt;
		
		glClear(GL_COLOR_BUFFER_BIT);
		
		dt = glfwGetTime();
		render();
		renderTime = glfwGetTime() - dt;
		cout << "render: " << setw(6)  << renderTime << " | update: " << updateTime << "\n";
		glfwSwapBuffers(window);
		glfwPollEvents();
		endFrame = glfwGetTime();
		frameTime = endFrame - startFrame;
		accumulate += frameTime;
		frameCount++;

		//update caption every second
		if(accumulate >= 1.f) {
			ostringstream os;
			os << "FPS: " << frameCount << " Particles: " << fluid->getParticles().size();
			string s = os.str();
			glfwSetWindowTitle(window, s.c_str());
			accumulate = 0;
			frameCount = 0;
		}
	}

	glfwTerminate();

	return 0;
}