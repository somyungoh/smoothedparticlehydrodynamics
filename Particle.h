#ifndef SM_PARTICLE_H_
#define SM_PARTICLE_H_

/****************************************************************
*
*		Particle.h
*
*		Particle Class.
*
****************************************************************/

#include "glm\glm.hpp"
#include "GL\glut.h"
#define PI 3.14159265359

using namespace glm;

class Particle{
public:
	// constructors
	Particle();
	Particle(int id, const vec2 &pos, float radius, float mass);
	~Particle();

	// getter & setters
	int   getID()			const;
	vec2  getPosition()		const;
	vec2  getVelocity()		const;
	vec2  getAcceleration()		const;
	float getRadius()		const;
	float getMass()			const;
	float getDensity()		const;
	float getPressure()		const;
	void  setPosition(const vec2 &p);
	void  setVelocity(const vec2 &v);
	void  setAcceleration(const vec2 &a);
	void  setRadius(float r);
	void  setMass(float m);
	void  setDensity(float d);
	void  setPressure(float p);
	void  setColor(const vec3 &color);

	// display function
	void  draw();

private:
	int   id;
	vec2  position;
	vec2  velocity;
	vec2  acceleration;
	float radius;
	float mass;
	float density;
	float pressure;
	vec3  color;
	bool  using_defaultColor;
};

#endif // !SM_PARTICLE_H_
