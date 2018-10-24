#include "SPHSolver.h"
#include <iostream>

/********************************************************
*
*	SPHSolver.cpp
*	
*	Implmenation of SPHSolver Class.
*
*	This class consists all the particle data
*	that is used in simulation as well as different
*	solvers.
*
*	Notes about finding optimal values
*
*	lower d0 & high B -> more dynamics. (easy to get crazy)
*	high a -> help to keep the volume
*	g -> 3 was the best. Others had terrible output.
*	e -> could not find difference in change
*
*********************************************************/



//***************************************//
//				Constructors			 //
//***************************************//

SPHSolver::SPHSolver() {};
SPHSolver::SPHSolver(const vec2 &LLC, const vec2 &URC, int initial_num_particles) : LLC(LLC), URC(URC), TOTAL_PARTICLE(initial_num_particles) {

	SOLVER_MODE = LEAPFROG;	// initial solver mode

	// initialize parameters
	gravity = vec2(0, -9.81); // gravity
	d0 = 200;				// initial density
	g = 3;					// tait exp. factor (gamma)
	a = 6.0;				// viscosity constant
	e = 0.05;				// viscosity epsilon
	B = 0.0009;				// pressure constant
	h = 0.015;				// radius
	m = 0.1;				// mass
	wall_sticky = 0.5;  // wall stickiness for boundary condition


	// these working good!

	//gravity = vec2(0, -9.81); // gravity
	//d0 = 25;				// initial density
	//g = 3;					// tait exp. factor (gamma)
	//a = 6.0;				// viscosity constant
	//e = 0.15;				// viscosity epsilon
	//B = 0.05;				// pressure constant
	//h = 0.15;				// radius
	//m = 0.1;				// mass
	//wall_sticky = 0.5;  // wall stickiness for boundary condition

	Nx = (int)((URC.x - LLC.x) / h + 1);
	Ny = (int)((URC.y - LLC.y) / h + 1);
	dx = (URC.x - LLC.x) / (float)(Nx - 1);
	dy = (URC.y - LLC.y) / (float)(Ny - 1);

	init_particles();
};
SPHSolver::~SPHSolver() {}

void SPHSolver::apply_change() {

#pragma omp parallel for
	for (int i = 0; i < TOTAL_PARTICLE; i++) {
		particles[i].setRadius(h);
	}
	Nx = (int)((URC.x - LLC.x) / h + 1);
	Ny = (int)((URC.y - LLC.y) / h + 1);
	dx = (URC.x - LLC.x) / (float)(Nx - 1);
	dy = (URC.y - LLC.y) / (float)(Ny - 1);

};

// getters and setters
int	  SPHSolver::getSolverMode()	const { return SOLVER_MODE; };
int   SPHSolver::getParticleNum()	const { return TOTAL_PARTICLE; };
vec2  SPHSolver::getGravity()		const { return gravity; };
float SPHSolver::getD0()			const { return d0; };
float SPHSolver::getG()				const { return g; };
float SPHSolver::getA()				const { return a; };
float SPHSolver::getB()				const { return B; };
float SPHSolver::getE()				const { return e; };
float SPHSolver::getH()				const { return h; };
float SPHSolver::getMass()			const { return m; };
float SPHSolver::getWS()			const { return wall_sticky; };
void  SPHSolver::setSolverMode(int MODE) { SOLVER_MODE = MODE; }
void  SPHSolver::setGravity(const vec2 &g) { gravity = g; };
void  SPHSolver::setD0(float _d0)		  { d0 = _d0; };
void  SPHSolver::setG(float _g)			  { g = _g; };
void  SPHSolver::setA(float _a)			  { a = _a; };
void  SPHSolver::setB(float _B)			  { B = _B; };
void  SPHSolver::setE(float _e)			  { e = _e; };
void  SPHSolver::setH(float _h)			  { h = _h; };
void  SPHSolver::setMass(float _m)		  { m = _m; };
void  SPHSolver::setWS(float ws)		  { wall_sticky = ws; };


//***************************************//
//			simulation methods			 //
//***************************************//


void SPHSolver::simulate_oneStep(float dt) {
	
	if (SOLVER_MODE == EULER)		sim_eulerian(dt);
	if (SOLVER_MODE == LEAPFROG)	sim_leapfrog(dt);
	if (SOLVER_MODE == SIXTH)		sim_sixth(dt);

	std::cout << particles[0].getDensity() << std::endl;
};

void SPHSolver::sim_eulerian(float dt) {
	update_density();
	update_pressure();
	update_momentum();
	update_velocity(dt);
	update_position(dt);
	boundary_check();
	update_occupancy();
};

void SPHSolver::sim_leapfrog(float dt) {

	update_position(0.5f * dt);	// half-step
	boundary_check();
	update_occupancy();
	update_density();			
	update_pressure();
	update_momentum();
	update_velocity(dt);
	update_position(0.5f * dt); // another half-step
	boundary_check();
	update_occupancy();
};

void SPHSolver::sim_sixth(float dt) {
	
	float a = 1 / (4 - pow(4, 1 / 3));	// coefficients
	float b = 1 - 4 * a;

	sim_leapfrog(a * dt);
	sim_leapfrog(a * dt);
	sim_leapfrog(b * dt); 
	sim_leapfrog(a * dt);
	sim_leapfrog(a * dt);
};



// update density of all particles
void SPHSolver::update_density() {

	vec2 pos_p1, pos_p2;

	for (int i = 0; i < TOTAL_PARTICLE; i++) {

		float new_density = 0;				// stores new density
		pos_p1 = particles[i].getPosition();
		
		std::vector<size_t> neighbors;		// get neighbor particles (in occupancy volume)
		get_nearbyParticles(pos_p1, neighbors);

		for (int j = 0; j < neighbors.size(); j++) {
			int index = neighbors[j];

			pos_p2 = particles[index].getPosition();
			new_density += particles[index].getMass() * compute_weight(pos_p1, pos_p2, particles[i].getRadius());
		}
		particles[i].setDensity(new_density);
	}
};


// update pressure of all particles
void SPHSolver::update_pressure() {

#pragma omp parallel for
	for (int i = 0; i < TOTAL_PARTICLE; i++) {

		float new_pressure = B * (pow(particles[i].getDensity() / d0, g) - 1);
		particles[i].setPressure(new_pressure);
	}
};


// update momentum(acceleration) of all particles
void SPHSolver::update_momentum() {

#pragma omp parallel for
	for (int i = 0; i < TOTAL_PARTICLE; i++) {
		
		vec2 new_acceleration(0, 0);

		std::vector<size_t> neighbors;		// get neighbor particles (in occupancy volume)
		get_nearbyParticles(particles[i].getPosition(), neighbors);

		for (int j = 0; j < neighbors.size(); j++) {

			int index = neighbors[j];
			if (i == index) continue;		// ignore same index

			// 1. compute force
			float S_force = compute_force(particles[i], particles[index]);
			// 2. compute gradient weight
			vec2 grad = compute_gradientWeight(particles[i].getPosition(), particles[index].getPosition(), particles[i].getRadius());		// gradient

			new_acceleration += particles[index].getMass() * S_force * grad;		// final gathering

		}
		particles[i].setAcceleration(-1.f * new_acceleration + gravity);
	}
};


// update velocity due to acceleration, timestep
void SPHSolver::update_velocity(float dt) {

#pragma omp parallel for
	for (int i = 0; i < TOTAL_PARTICLE; i++) {
		vec2 new_velocity(particles[i].getVelocity() + particles[i].getAcceleration() / particles[i].getMass() * dt);
		particles[i].setVelocity(new_velocity);
	}
};


// update position due to acceleration, timestep
void SPHSolver::update_position(float dt) {

#pragma omp parallel for
	for (int i = 0; i < TOTAL_PARTICLE; i++) {
		vec2 new_position(particles[i].getPosition() + particles[i].getVelocity() * dt);
		particles[i].setPosition(new_position);
	}
};


// checks boundary due to LLC, URC
void SPHSolver::boundary_check() {

	for (int i = 0; i < TOTAL_PARTICLE; i++) {

		vec2 pos = particles[i].getPosition();
		vec2 vel = particles[i].getVelocity();

		if (pos.x < LLC.x) {		// left boundary check
			vec2 newPos(LLC.x + (LLC.x - 1 * pos.x) * wall_sticky, pos.y);
			vec2 newVel(-1.f * vel.x * wall_sticky, vel.y);
			particles[i].setPosition(newPos);
			particles[i].setVelocity(newVel);
		}
		if (pos.x > URC.x) {		// right boundary check
			vec2 newPos(URC.x + (URC.x -1 * pos.x) * wall_sticky, pos.y);
			vec2 newVel(-1.f * vel.x * wall_sticky, vel.y);
			particles[i].setPosition(newPos);
			particles[i].setVelocity(newVel);

		}
		if (pos.y < LLC.y) {		// bottom boundary check
			vec2 newPos(pos.x, LLC.y + (LLC.y - 1 * pos.y) * wall_sticky);
			vec2 newVel(vel.x, -1.f * vel.y * wall_sticky);
			particles[i].setPosition(newPos);
			particles[i].setVelocity(newVel);
		}
		if (pos.y > URC.y) {		// top boundary check
			vec2 newPos(pos.x, URC.y + (URC.y -1 * pos.y) * wall_sticky);
			vec2 newVel(vel.x, -1 * vel.y * wall_sticky);
			particles[i].setPosition(newPos);
			particles[i].setVelocity(newVel);
		}
	}
};


// generate occupancy volume
void SPHSolver::update_occupancy() {
	
	// init. occupancy volume
	occupancy_volume.clear();
	for (int i = 0; i < Nx * Ny; i++) occupancy_volume.push_back(std::vector<size_t>());

	for (int i = 0; i < TOTAL_PARTICLE; i++) {

		// 1. boundary check
		vec2 p = particles[i].getPosition() - LLC;
		if (p.x < 0 || p.y <0 || p.x > URC.x - LLC.x || p.y > URC.y - LLC.y) continue;

		// 2. insert current particle
		//	  to occupancy volume at index
		int ix = (int)(p.x / dx);
		int iy = (int)(p.y / dy);
		int index = ix + iy * Nx;

		if (index < 0 || index >= occupancy_volume.size()) continue;		// index over bound

		occupancy_volume[index].push_back(i);		// add this particle index
	}
};




//***************************************//
//				helper methods			 //
//***************************************//


// display method
void SPHSolver::display() {

	for (int i = 0; i < TOTAL_PARTICLE; i++) { particles[i].draw(); }
}

// called by constructor
// initialize particles due to intial number of particles
void SPHSolver::init_particles() {

	// update particle array
	for (int i = 0; i < TOTAL_PARTICLE; i++) {
		// initial particle attributes

		int Nwidth = sqrt(TOTAL_PARTICLE);

		float px = (i % Nwidth) / (float)Nwidth * 0.45f + 0.1f;
		float py = (int)(i / Nwidth) / (float)Nwidth * 0.65f - 0.3f;

		vec2  position(px, py);
		particles.push_back(Particle(i, position, h, m));
	}
};


void SPHSolver::addParticleDump(float cx, float cy, int amount) {
	
	for (int i = 0; i < amount; i++) {
	
		float x = ((float)rand() / (float)RAND_MAX) * h * 0.5 - 0.25f * h + cx;
		float y = ((float)rand() / (float)RAND_MAX) * h * 0.5 - 0.25f * h + cy;

		Particle new_particle(++TOTAL_PARTICLE, vec2(x,y), h, m);
		new_particle.setColor(vec3(1.0, 0.3647, 0.333));

		particles.push_back(new_particle);
	}
};


// returns list of particle indices that is inside of occupancy
void SPHSolver::get_nearbyParticles(const vec2 &pos, std::vector<size_t> &list) {

	vec2 p = pos - LLC;
	// 1. boundary check
	if (p.x < 0 || p.y <0 || p.x > URC.x - LLC.x || p.y > URC.y - LLC.y) return;
	
	int ix = (int)(p.x / dx);
	int iy = (int)(p.y / dy);

	int index[9];							// get 9 neighboring grids
	index[0] = ix - 1 + (iy - 1) * Nx;
	index[1] = ix - 1 + (iy + 0) * Nx;
	index[2] = ix - 1 + (iy + 1) * Nx;
	index[3] = ix + 0 + (iy - 1) * Nx;
	index[4] = ix + 0 + (iy + 0) * Nx;
	index[5] = ix + 0 + (iy + 1) * Nx;
	index[6] = ix + 1 + (iy - 1) * Nx;
	index[7] = ix + 1 + (iy + 0) * Nx;
	index[8] = ix + 1 + (iy + 1) * Nx;

	// collect particles in occupancy grid
	for (int i = 0; i < 9; i++) {
		int idx = index[i];
		if (idx < 0 || idx >= occupancy_volume.size()) continue;		// index over bound

		for (int j = 0; j < occupancy_volume[idx].size(); j++) {
			list.push_back(occupancy_volume[idx][j]);
		}
	}
};


// called by update momentum()
// calculates force due to pressure + viscosity
float SPHSolver::compute_force(const Particle &p1, const Particle &p2) {

	float S_force, S_pressure, S_viscosity;	// variables that need to be computed.
	double v;
	vec2  Xab, Vab;

	vec2  pos_p1	= p1.getPosition();		// prepare particle data
	vec2  pos_p2	= p2.getPosition();
	vec2  vel_p1	= p1.getVelocity();
	vec2  vel_p2	= p2.getVelocity();
	float press_p1	= p1.getPressure();
	float press_p2	= p2.getPressure();
	float den_p1	= p1.getDensity();
	float den_p2	= p2.getDensity();


	// 1. compute pressure
	S_pressure = (press_p1 / (den_p1 * den_p1)) + (press_p2 / (den_p2 * den_p2));

	// 2. compute viscosity
	Xab = pos_p1 - pos_p2;
	Vab = vel_p1 - vel_p2;
	v = (a * p1.getRadius()) / (den_p1 + den_p2);
	double dotVXab = dot(Vab, Xab);
	double dotXab  = dot(Xab, Xab);

	if (dotVXab < 0)  S_viscosity = -1.0f * v * (dotVXab / (dotXab + e * pow(p1.getRadius(), 2)));
	else			  S_viscosity = 0;

	S_force = S_pressure + S_viscosity;

	// **** DEBUGGING INFINITE FORCES
	if (isinf(S_force)) std::cout << "ERROR::Infinite Pressure Occured!!" << std::endl;
	if (isinf(S_force)) S_force = 1.0;	// kill inifinite forces
	
	return S_force;
};



// returns weight between to given particles
float SPHSolver::compute_weight(const vec2 &p1, const vec2 &p2, float radius) {
	
	float r, rh;
	float w;

	r = length(p1 - p2);
	rh = r / radius;

	if (rh < 1) w = pow(1 - rh, 3) * 10.f / (PI * pow(radius, 3)); 	// distance is less than radius
	else		w = 0;

	return w;
};


// returns gradient weight between to given particles
vec2 SPHSolver::compute_gradientWeight(const vec2 &p1, const vec2 &p2, float radius) {

	float r_len = length(p1 - p2);
	vec2  r_dir  = p1 - p2;
	vec2  grad_w;

	if (abs(r_len) < 1e-10) grad_w = vec2(0, 0);
	if (r_len / radius < 1) grad_w = (float)(- 1.f * pow(1 - r_len / radius, 2) * 30 / (PI * pow(radius, 3))) * (r_dir / r_len);
	else					grad_w = vec2(0, 0);

	return grad_w;
};