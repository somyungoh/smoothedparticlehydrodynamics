#include "Particle.h"

/****************************************************************
*
*		Particle.cpp
*
*		Implementation of Particle Class.
*
****************************************************************/

Particle::Particle() {}
Particle::Particle(int id, const vec2 &pos, float radius, float mass) : id(id), position(pos), radius(radius), mass(mass) {

	// init. attributes
	velocity	 = vec2(0,0);
	acceleration = vec2(0,0);
	density		 = 0;
	pressure	 = 0;
	using_defaultColor = true;

}
Particle::~Particle() {}
 
// getter & setters
int	  Particle::getID()				const { return id; };
vec2  Particle::getPosition()		const { return position; };
vec2  Particle::getVelocity()		const { return velocity; };
vec2  Particle::getAcceleration()	const { return acceleration; };
float Particle::getRadius()			const { return radius; };
float Particle::getMass()			const { return mass; };
float Particle::getDensity()		const { return density; };
float Particle::getPressure()		const { return pressure; };
void  Particle::setPosition(const vec2 &p)		{ position = p; };
void  Particle::setVelocity(const vec2 &v)		{ velocity = v; };
void  Particle::setAcceleration(const vec2 &a)	{ acceleration = a; };
void  Particle::setRadius(float r)				{ radius = r; };
void  Particle::setMass(float m)				{ mass = m; };
void  Particle::setDensity(float d)				{ density = d; };
void  Particle::setPressure(float p)			{ pressure = p; };
void  Particle::setColor(const vec3 &_color)	{ color = _color; using_defaultColor = false; }

// display function
void  Particle::draw() {

	vec3 dark_blue(0.00213, 0.5, 0.8549);
	vec3 light_blue(0.615, 1.0, 0.7254);
	
///	float t = density / 1800.f > 1 ? 1 : density / 1800.f;
	float t = density / 300000.f > 1 ? 1 : density / 300000.f;

	if(using_defaultColor) color = t * dark_blue + (1 -t) * light_blue;

	glColor3f(color.r, color.g, color.b);
	glPointSize(2.0f);
	glBegin(GL_POINTS);
	glVertex2f(position.x, position.y);
	glEnd();
};