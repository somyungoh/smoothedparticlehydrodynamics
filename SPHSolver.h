#ifndef SM_SPHSOLVER_H_
#define SM_SPHSOLVER_H_

#include "glm\glm.hpp"
#include "Particle.h"
#include <vector>

using namespace glm;

class SPHSolver
{
public:
	SPHSolver();
	SPHSolver(const vec2 &LLC, const vec2 &URC, int initial_num_particles);
	~SPHSolver();

	enum MODE {EULER, LEAPFROG, SIXTH};

	// main simulation method
	void simulate_oneStep(float dt);
	void display();
	void apply_change();

	// getters and setters
	int   getSolverMode()	const;
	int   getParticleNum()  const;
	vec2  getGravity()	const;
	float getD0()		const;
	float getG()		const;
	float getA()		const;
	float getB()		const;
	float getE()		const;
	float getH()		const;
	float getMass()		const;
	float getWS()		const;
	void  setSolverMode(int MODE);
	void  setGravity(const vec2 &g);
	void  setD0(float d0);
	void  setG(float g);
	void  setA(float a);
	void  setB(float B);
	void  setE(float e);
	void  setH(float h);
	void  setMass(float m);
	void  setWS(float ws);
	void  addParticleDump(float cx, float cy, int amount);

private:
	
	// volume
	vec2  LLC;	// lower left corner
	vec2  URC;	// upper right corner
	int   Nx, Ny;	// occupancy grids
	float dx, dy;	// occupancy grid size (approx. to radius)

	// particle set
	std::vector<Particle> particles;
	std::vector<std::vector<size_t>> occupancy_volume;

	// variables
	int TOTAL_PARTICLE;
	int SOLVER_MODE;
	
	// solver control parameters
	vec2  gravity;	 	// gravity (m/s^2)
	float d0;		// initial density  (kg/m^3)
	float g;		// tait exp. factor (gamma)
	float a;		// viscosity constant (alpha)
	float B;		// pressure constant
	float e;		// viscosity epsilon
	float wall_sticky;	// wall stickiness for boundary condition

	// particle control parameters
	float h;		// radius (m)
	float m;		// mass	  (kg)

	// sim methods
	void sim_eulerian(float dt);
	void sim_leapfrog(float dt);
	void sim_sixth(float dt);
	void update_density();
	void update_pressure();
	void update_momentum();
	void update_velocity(float dt);
	void update_position(float dt);
	void boundary_check();
	void update_occupancy();

	// helper methods
	void init_particles();
	void get_nearbyParticles(const vec2 &pos, std::vector<size_t> &list);
	float compute_force(const Particle &p1, const Particle &p2);
	float compute_weight(const vec2 &p1, const vec2 &p2, float radius);
	vec2  compute_gradientWeight(const vec2 &p1, const vec2 &p2, float radius);
};

#endif // !SM_SPHSOLVER_H_
