/****************************************************************
*
*		Somyung(David) Oh
*
*		Project4 [Smooth Particle Hydrodynamics]
*
*		Advanced Topics in Physically Based Modeling
*		2018 Spring, Texas A&M University
*		Instructor - Jerry Tessendorf
*
*
****************************************************************/

#define _CRT_SECURE_NO_WARNINGS

#include <Windows.h>
#include <omp.h>
#include <iostream>
#include <sstream>
#include <string>
#include "GL\glui.h"
#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image.h"
#include "stb_image_write.h"
#include "SPHSolver.h"

// define constants
#define WIDTH			1280
#define HEIGHT			720
#define INITIALPARTICLE 50000
#define ADD_AMOUNT		6
#define MAP_SIZE		WIDTH * HEIGHT
#define TSTEP_STEP		0.1

using namespace std;

//	MAIN SOLVER VARIABLE	//
SPHSolver solver;
float aspect_ratio = (float)WIDTH / (float)HEIGHT;

//		GLUI variables		//

int		sph_solverType;
bool	isSimulate;
float	time_step;
bool	capture_screen;				// screen capture
int		frame;
int		particles;
string  captured_file_basename;
int		xmouse_prev, ymouse_prev;	// mouse xy

// solver control parameters
float sph_gravity;		// gravity (m/s^2)
float sph_d0;			// initial density  (kg/m^3)
float sph_a;			// viscosity constant
float sph_g;			// tait exp. factor (gamma)
float sph_B;			// pressure constant
float sph_e;			// viscosity epsilon
float sph_ws;			// wall stickiness for boundary condition
float sph_h;			// radius (m)
float sph_m;			// mass	  (kg)

// glui variables
int				main_window;
GLUI			*glui;		
GLUI_Panel		*panel_solverType, *panel_solverControl, *panel_description, *panel_status;
GLUI_RadioGroup	*radio_solverType;
GLUI_EditText	*tb_timestep, *tb_gravity, *tb_d0, *tb_a, *tb_g, *tb_B, *tb_e, *tb_ws, *tb_h, *tb_m;
GLUI_EditText	*tb_frame, *tb_p_num;
GLUI_StaticText *txt_blank, *txt_d0, *txt_a, *txt_g, *txt_B, *txt_e, *txt_ws, *txt_h, *txt_m;
GLUI_Button		*button_apply, *button_run, *button_reset, *button_capture;

#define ID_PANEL_PARAMETER		100
#define ID_RADIO_SOLVER			101
#define ID_BUTTON_APPLY			200
#define ID_BUTTON_RUN			201
#define ID_BUTTON_RESET			202
#define ID_BUTTON_CAPTURE		203




void writeImage(const char* file_name, unsigned char * img) {

	string dir = "Render";		// saving directory
	string format = "jpg";		// save format
	string out_name = dir + "/" + file_name + "." + format;		// full output name

	unsigned char *out_map = new unsigned char[WIDTH * HEIGHT * 3];
																// copy data
	for (int y = 0; y <HEIGHT; y++) {
#pragma omp parallel for
		for (int x = 0; x < WIDTH; x++) {
			int index = (y * WIDTH + x) * 3;
			int out_index = ((HEIGHT - y - 1) * WIDTH + x) * 3;
			out_map[out_index + 0] = img[index + 0];
			out_map[out_index + 1] = img[index + 1];
			out_map[out_index + 2] = img[index + 2];
		}
	}

	stbi_write_jpg(out_name.c_str(), WIDTH, HEIGHT, 3, out_map, 100);

	delete[] out_map;
}


// called when the APPLY button is pressed
void apply_change() {

	solver.setSolverMode(sph_solverType);
	solver.setGravity(vec2(0, sph_gravity));
	solver.setD0(sph_d0);
	solver.setA(sph_a);
	solver.setG(sph_g);
	solver.setB(sph_B);
	solver.setE(sph_e);
	solver.setWS(sph_ws);
	solver.setH(sph_h);
	solver.setMass(sph_m);
	solver.apply_change();
}


//----------------------------------------------------
//
//  GL and GLUT callbacks
//
//----------------------------------------------------

void reshape(int x, int y)
{
	glViewport(0, 0, x, y);

	glutPostRedisplay();
}

void cbDisplay(void)
{
	glClear(GL_COLOR_BUFFER_BIT);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-1.f, 1.f, -(1.f / aspect_ratio), (1.f / aspect_ratio), -1.f, 1.f);
	glMatrixMode(GL_MODELVIEW);

	glClearColor(0.26, 0.26, 0.3, 1);

	solver.display();

	glutSwapBuffers();
	glFlush();

	if (capture_screen)
	{
		std::stringstream os; os << frame;
		string dispframe = os.str();
		if (frame < 1000) { dispframe = "0" + dispframe; }
		if (frame < 100) { dispframe = "0" + dispframe; }
		if (frame < 10) { dispframe = "0" + dispframe; }
		string fname = captured_file_basename + "." + dispframe;

		unsigned char *capture_map = new unsigned char[WIDTH * HEIGHT * 3];
		glReadPixels(0, 0, WIDTH, HEIGHT, GL_RGB, GL_UNSIGNED_BYTE, capture_map);
		writeImage(fname.c_str(), capture_map);
		delete[] capture_map;
		cout << "Frame written to file " << fname << endl;
	}
}

// animate and display new result
void cbIdle()
{
	if (glutGetWindow() != main_window)
		glutSetWindow(main_window);

	// read values
	particles = solver.getParticleNum();

	// run simulation
	if (isSimulate) {
		solver.simulate_oneStep(time_step);
		frame++;
	}
	glutPostRedisplay();
	glui->sync_live();
}

void cbOnKeyboard(unsigned char key, int x, int y)
{
	switch (key)
	{
	case '+': case '=':
		break;

	default:
		break;
	}
}

void cbMouseDown(int button, int state, int x, int y)
{
	if (button != GLUT_LEFT_BUTTON) { return; }
	if (state != GLUT_DOWN) { return; }
	xmouse_prev = x;
	ymouse_prev = y;

	float px = (2 * x - WIDTH) / (float)WIDTH;				// add particles on mouse click!
	float py = (2 * (HEIGHT - y) - HEIGHT) / (float)HEIGHT / aspect_ratio;
	solver.addParticleDump(px, py, ADD_AMOUNT);

	cout << "INPUT::Added Particles\n";
}

void cbMouseMove(int x, int y)
{
	xmouse_prev = x;
	ymouse_prev = y;
}



void PrintUsage()
{
	cout << "==============================================================\n";
	cout << "\n";
	cout << "\t Project4 - Smooth Particle Hydrodynamics\n";
	cout << "\t Made by Somyung (David) Oh.\n";
	cout << "\n";
	cout << "\t Advanced Topics and Physically Based Modeling\n";
	cout << "\t Spring 2018, Texas A&M University.\n";
	cout << "\t Instructor: Dr. Jerry Tessendorf\n";
	cout << "\n";
	cout << "\t Note: Click screen to add particles.\n";
	cout << "\n";
	cout << "==============================================================\n";
}


void init_components() {
	
	vec2 LLC	= vec2(-1, -(1.f / aspect_ratio));
	vec2 URC	= vec2(1, (1.f / aspect_ratio));

	solver		= SPHSolver(LLC, URC, INITIALPARTICLE);
	sph_solverType = solver.getSolverMode();
	sph_gravity = solver.getGravity().y;
	sph_d0		= solver.getD0();
	sph_a		= solver.getA();
	sph_g		= solver.getG();
	sph_B		= solver.getB();
	sph_e		= solver.getE();
	sph_ws		= solver.getWS();
	sph_h		= solver.getH();
	sph_m		= solver.getMass();

	frame = 1;									// screen capture
	capture_screen = false;
	captured_file_basename = "sph";
	isSimulate = false;

	time_step = 0.0005;
}

void buttonCallback(GLUI_Control* control) {

	switch (control->get_id()) {
	case ID_BUTTON_APPLY:
		apply_change();
		cout << "INPUT::Applied changes.\n";
		break;

	case ID_BUTTON_RUN:
		isSimulate = !isSimulate;
		if (isSimulate) button_run->set_name("Pause");
		else button_run->set_name("Resume");
		cout << "INPUT::Simulate " << isSimulate << endl;
		break;

	case ID_BUTTON_RESET:
		init_components();
		cout << "INPUT::Reset Simulation\n";
		break;

	case ID_BUTTON_CAPTURE:
		capture_screen = !capture_screen;
		if (capture_screen) button_capture->set_name("Stop Capture");
		else button_capture->set_name("Capture");
		cout << "INPUT::Capture " << capture_screen << endl;
		break;
	}
}


// initialize GLUI
void initGLUI() {

	// GLUI variables
	glui = GLUI_Master.create_glui("Controller");

	panel_solverControl = new GLUI_Panel(glui, "Solver Controller");
	tb_timestep = new GLUI_EditText(panel_solverControl, "Timestep", GLUI_EDITTEXT_FLOAT, &time_step, -1, buttonCallback);
	tb_gravity = new GLUI_EditText(panel_solverControl, "gravity", GLUI_EDITTEXT_FLOAT, &sph_gravity, -1, buttonCallback);
	tb_d0 = new GLUI_EditText(panel_solverControl, "d0", GLUI_EDITTEXT_FLOAT, &sph_d0, -1, buttonCallback);
	tb_a = new GLUI_EditText(panel_solverControl, "a", GLUI_EDITTEXT_FLOAT, &sph_a, -1, buttonCallback);
	tb_g = new GLUI_EditText(panel_solverControl, "g", GLUI_EDITTEXT_FLOAT, &sph_g, -1, buttonCallback);
	tb_B = new GLUI_EditText(panel_solverControl, "B", GLUI_EDITTEXT_FLOAT, &sph_B, -1, buttonCallback);
	tb_e = new GLUI_EditText(panel_solverControl, "e", GLUI_EDITTEXT_FLOAT, &sph_e, -1, buttonCallback);
	tb_ws = new GLUI_EditText(panel_solverControl, "ws", GLUI_EDITTEXT_FLOAT, &sph_ws, -1, buttonCallback);
	tb_h = new GLUI_EditText(panel_solverControl, "h", GLUI_EDITTEXT_FLOAT, &sph_h, -1, buttonCallback);
	tb_m = new GLUI_EditText(panel_solverControl, "m", GLUI_EDITTEXT_FLOAT, &sph_m, -1, buttonCallback);
	button_apply = new GLUI_Button(panel_solverControl, "Apply", ID_BUTTON_APPLY, buttonCallback);

	glui->add_column(true);		// second column

	panel_description = new GLUI_Panel(glui, "Description");
	txt_blank = new GLUI_StaticText(panel_description, "");
	txt_d0 = new GLUI_StaticText(panel_description, "d0: base density");
	txt_a = new GLUI_StaticText(panel_description, "a:  viscosity_c");
	txt_g = new GLUI_StaticText(panel_description, "g:  tait exp");
	txt_B = new GLUI_StaticText(panel_description, "B:  pressure_c");
	txt_e = new GLUI_StaticText(panel_description, "e:  viscosity_epsilon");
	txt_ws = new GLUI_StaticText(panel_description, "ws: wall sticky");
	txt_h = new GLUI_StaticText(panel_description, "h:  radius");
	txt_m = new GLUI_StaticText(panel_description, "m:  mass");

	glui->add_column(true);

	panel_solverType = new GLUI_Panel(glui, "Solver Type");
	radio_solverType = new GLUI_RadioGroup(panel_solverType, &sph_solverType, ID_RADIO_SOLVER, buttonCallback);
	new GLUI_RadioButton(radio_solverType, "Eulerian");
	new GLUI_RadioButton(radio_solverType, "Leap Frog");
	new GLUI_RadioButton(radio_solverType, "Sixth");

	panel_status = new GLUI_Panel(glui, "Sim Status");
	tb_frame = new GLUI_EditText(panel_status, "frame", GLUI_EDITTEXT_INT, &frame, -1, buttonCallback);
	tb_p_num = new GLUI_EditText(panel_status, "particles", GLUI_EDITTEXT_INT, &particles, -1, buttonCallback);

	button_run = new GLUI_Button(glui, "Run", ID_BUTTON_RUN, buttonCallback);
	button_reset = new GLUI_Button(glui, "Reset", ID_BUTTON_RESET, buttonCallback);
	
	glui->add_separator();
	button_capture = new GLUI_Button(glui, "Capture", ID_BUTTON_CAPTURE, buttonCallback);

	glui->set_main_gfx_window(main_window);
	GLUI_Master.set_glutIdleFunc(cbIdle);

}


// main
int main(int argc, char** argv)
{

	init_components();
	PrintUsage();
	omp_set_num_threads(4);

	// GLUT routines
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
	glutInitWindowSize(WIDTH, HEIGHT);
	main_window = glutCreateWindow("Project4 - Somyung Oh");
	glutDisplayFunc(&cbDisplay);
	glutReshapeFunc(reshape);
	glutKeyboardFunc(&cbOnKeyboard);
	glutMouseFunc(&cbMouseDown);
	glutMotionFunc(&cbMouseMove);

	initGLUI();		// init GLUI

	glutMainLoop();
	return 1;
}
