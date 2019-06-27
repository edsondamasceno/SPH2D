#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include <stdlib.h>
#include <iostream>
#include "sph.hh"

using namespace std;

int winWidth = 512, winHeight = 512;
static int ox, oy;
bool buttondown[3] = {false, false, false};

float mouse_force = 200;
Vector mouse_point(std::numeric_limits<double>::infinity(), 
                   std::numeric_limits<double>::infinity(), 
                   std::numeric_limits<double>::infinity());
std::vector<NeighborData> selected;

struct {
  const char *name;
  Parameters params;
} parameters[] = {
// name          dt    r  gamma K rho0 zeta  nu xi  xsph             g          a 
{"incompressible", {0.001, 2.5, 1, 500, 1, 0.5, 0, 0.01, false, Vector(0, -9.81, 0), 0.3}},
{"compressible", {0.001, 2.5, 1, 20, 1, 0.5, 0, 0.01, false, Vector(0, -9.81, 0), 0.3}},
{"real XSPH", {0.001, 2.5, 1, 500, 1, 0.5, 0, 0.01,  true, Vector(0, -9.81, 0), 0.3}},
{"fake XSPH high", {0.001, 2.5, 1, 500, 1, 0.5, 0, 1, false, Vector(0, -9.81, 0), 0.3}},
{"viscosity too low", {0.001, 2.5, 1, 500, 1, 0.5, 0.001, 0, false, Vector(0, -9.81, 0), 0.3}},
{"viscosity ok (low)", {0.001, 2.5, 1, 500, 1, 0.5, 0.002, 0, false, Vector(0, -9.81, 0), 0.3}},
{"viscosity highest", {0.001, 2.5, 1, 500, 1, 0.5, 0.4, 0, false, Vector(0, -9.81, 0), 0.3}},
{"viscosity too high", {0.001, 2.5, 1, 500, 1, 0.5, 0.5, 0, false, Vector(0, -9.81, 0), 0.3}},
{"low smoothing radius", {0.001, 1.8, 1, 500, 1, 0.5, 0, 0.1, false, Vector(0, -9.81, 0), 0.3}},
{"high smoothing radius", {0.001, 5, 1, 500, 1, 0.5, 0, 0.01, false, Vector(0, -9.81, 0), 0.3}}
};

SPH sph;

bool running = false, drawgrid = false, drawsurface = false, drawhash = false, drawtree = false;

void redraw(void)
{
  glClearColor(0, 0, 0.1, 0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  glBegin(GL_QUADS);
  glColor3f(0, 0, 0.3);
  glVertex2f(0, 0);
  glVertex2f(0, 10);
  glVertex2f(1, 10);
  glVertex2f(1, 0);
  glEnd();  
  
  glColor3f(1, 1, 1);
  glPointSize(5);
  glEnable(GL_POINT_SMOOTH);
  glBegin(GL_POINTS);
  // draw particles
  for (int i = 0; i < sph.n_particles(); ++i) {
    Vector const &p = sph.particle(i).x;
    glVertex2f(p.x, p.y);
  }
  // selected particles in red
  glColor3f(1, 0, 0);
  for (int i = 0; i < selected.size(); ++i) {
    Vector const &p = sph.particle(selected[i].idx).x;
    glVertex2f(p.x, p.y);
  }
  glEnd();
 
  // if requested, draw search data structures
  if (drawhash) {
    sph.hashgrid().draw(mouse_point);
  }
  
  if (drawtree) {
    sph.kdtree().draw(mouse_point);
  }
  
  if (drawgrid) {
    sph.drawSurfaceGrid();
  }
  
  if (drawsurface) {
    sph.drawSurface();
  }
  
  glFinish();
  glutSwapBuffers();
}  

void reshape(int width, int height)
{
  glViewport(0, 0, width, height);
  winWidth = width;
  winHeight = height;
  
  glMatrixMode(GL_PROJECTION);
  
  glLoadIdentity();

  double r = (float)width/height;
  double xw = 1.2 * r;
  double yw = 1.2 / r;
  
  if (width < height)
    glOrtho(-.1, 1.1, -(yw-1)/2, (yw+1)/2, -1, 1);
  else
    glOrtho(-(xw-1)/2, (xw+1)/2, -.1, 1.1, -1, 1);    
    
  glMatrixMode(GL_MODELVIEW);  
  
  glutPostRedisplay();
}

Vector worldSpace(int x, int y) {
  int vp[4];
  double MV[16], P[16];
  double ox, oy, oz;
  glGetIntegerv(GL_VIEWPORT, vp);
  glGetDoublev(GL_MODELVIEW_MATRIX, MV);
  glGetDoublev(GL_PROJECTION_MATRIX, P);
  gluUnProject(x, winHeight-y, 0, MV, P, vp, &ox, &oy, &oz);
  
  oz = 0;
  Vector o(ox, oy, oz);

  return o;
}

void keyboard(unsigned char key, int x, int y)
{
  switch(key) {
    case '\033':
    case 'q':
      exit(0);
    case ' ':
      running = !running;
      break;
    case 't':
      drawtree = !drawtree;
      break;
    case 'h':
      drawhash = !drawhash;
      break;
    case 's':
      drawsurface = !drawsurface;
      break;
    case 'g':
      drawgrid = !drawgrid;
      break;
    case 'r':
      // reset
      sph.init();
      break;
    case '+':
      sph.smoothing_radius(1.1 * sph.smoothing_radius());
      break;
    case '-':
      sph.smoothing_radius(1./1.1 * sph.smoothing_radius());
      break;
    case '0':
    case '1': 
    case '2': 
    case '3': 
    case '4': 
    case '5': 
    case '6': 
    case '7': 
    case '8': 
    case '9': 
    {
      int set = key-'1';
      if (set < 0)
        set = 9;
      // set preset parameter values
      std::cout << "setting parameter set " << set << ": " << parameters[set].name << std::endl;
      sph.parameters(parameters[set].params);
      break;
    }
  }
  glutPostRedisplay();
}

void motion(int x, int y)
{
  Vector point = worldSpace(x,y);

  if(buttondown[GLUT_LEFT_BUTTON])
  { 
    sph.kdtree().neighbors(point, selected);

    // drag fluid
    for (int i = 0; i < selected.size(); ++i) {
      Particle &p = sph.particle(selected[i].idx);
      float k = sph.vkernel.value(selected[i].d, selected[i].d_squared);
      p.a += mouse_force * k * (point - mouse_point);
    }   
  }
  else if (buttondown[GLUT_RIGHT_BUTTON])
  {
  }
  
  mouse_point = point;
  glutPostRedisplay();
}

void button(int b, int state, int x, int y)
{
  ox = x;
  oy = y;
  
  // collect neighhoring particles
  if (state == GLUT_DOWN) {
    mouse_point = worldSpace(x,y);
    sph.kdtree().neighbors(mouse_point, selected);
    buttondown[b] = true;
  } else {
    selected.clear();
    buttondown[b] = false;
  }
  
  glutPostRedisplay();
}

void idle()
{
  if (running) {
    sph.step();
    glutPostRedisplay();
  }
}

int main(int argc, char ** argv) 
{
  glutInitWindowSize(winWidth, winHeight);
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_RGBA | GLUT_DEPTH | GLUT_ALPHA | GLUT_DOUBLE);
  glutCreateWindow("SPH");

  sph.init();

  // GL init
  glEnable(GL_DEPTH_TEST);
  glDepthFunc(GL_LEQUAL);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  
  glutDisplayFunc(redraw);
  glutMotionFunc(motion);
  glutPassiveMotionFunc(motion);
  glutMouseFunc(button);
  glutReshapeFunc(reshape);
  glutKeyboardFunc(keyboard);
  glutIdleFunc(idle);
  
  glutMainLoop();
  
  return 0;
}

