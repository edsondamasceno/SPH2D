#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glew.h>
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glext.h>
#endif

#include <stdlib.h>
#include <iostream>


#include "sph.hh"

#define		REFRESH_DELAY		1		//ms

using namespace std;
using namespace of;


//int winWidth = 800, winHeight = 600;
static int ox, oy;
bool buttondown[3] = {false, false, false};

float mouse_force = 500;
Vector mouse_point(std::numeric_limits<double>::infinity(), 
                   std::numeric_limits<double>::infinity(), 
                   std::numeric_limits<double>::infinity());
std::vector<NeighborData> selected;

struct {
  const char *name;
  Parameters params;
} parameters[] = {
// name          dt    r  gamma K rho0 rho1 zeta  nu xi  xsph             g          a 
{"incompressible", {0.001, 2.5, 1, 500, 1, 0.9, 0.5, 0, 0.01, false, Vector(0, -9.81, 0), 0.3}},
{"compressible", {0.001, 2.5, 1, 20, 1, 0.9, 0.5, 0, 0.01, false, Vector(0, -9.81, 0), 0.3}},
{"real XSPH", {0.001, 2.5, 1, 500, 1, 0.9, 0.5, 0, 0.01,  true, Vector(0, -9.81, 0), 0.3}},
{"fake XSPH high", {0.001, 3.0, 1, 500, 1, 0.9, 0.5, 0, 1, false, Vector(0, -9.81, 0), 0.3}},
{"viscosity too low", {0.001, 3.0, 1, 500, 1, 0.9, 0.5, 0.001, 0, false, Vector(0, -9.81, 0), 0.3}},
{"viscosity ok (low)", {0.001, 3.0, 1, 500, 1, 0.9, 0.5, 0.002, 0, false, Vector(0, -9.81, 0), 0.3}},
{"viscosity highest", {0.001, 3.0, 1, 500, 1, 0.9, 0.5, 0.4, 0, false, Vector(0, -9.81, 0), 0.3}},
{"viscosity too high", {0.001, 3.0, 1, 500, 1, 0.9, 0.5, 0.2, 0, false, Vector(0, -9.81, 0), 0.3}},
{"low smoothing radius", {0.001, 1.8, 1, 500, 1, 0.9, 0.5, 0, 0.1, false, Vector(0, -9.81, 0), 0.3}},
{"high smoothing radius", {0.001, 5, 1, 500, 1, 0.9, 0.5, 0, 0.01, false, Vector(0, -9.81, 0), 0.3}}
};

SPH sph;

bool running = false, drawgrid = false, drawsurface = false, drawhash = false, drawtree = false, obstok=false;
int saveFile =false;
int saveImageFile=false;
void RenderScene(void){ 
  	allCommands->Execute();
   //Print->Faces(white);
   
	if(sph.Tmax==sph.Tmin) 
	{
     	// Muda a cor do fluido
      	//glColor3f(0, 0, 1);
	  	//glColor3d(1,0,0);
      	//glPointSize(5);
      	//glEnable(GL_POINT_SMOOTH);
      	//glBegin(GL_POINTS);
      	// draw particles
   		for (int i = 0; i < sph.n_particles(); ++i) {
        if(i < (sph.maximoParticle/2)) {
          	//glColor3f(0, 0, 1);
          	glColor3f(1.0, 1.0, 0);
        	Vector const &p = sph.particle(i).x;
			//glColor3f(sph.particle(i).color.R,sph.particle(i).color.G,sph.particle(i).color.B);
			//glVertex3f(p.x, p.y,p.z);
			glPushMatrix();
			glTranslated(p.x, p.y, p.z);
			glutSolidSphere(sph.sphereRadius0,6,6);
			glPopMatrix();
        }
        else {
          	//glColor3f(1.0, 1.0, 0);
          	glColor3f(0, 0, 1);
          	Vector const &p = sph.particle(i).x;
			//glColor3f(sph.particle(i).color.R,sph.particle(i).color.G,sph.particle(i).color.B);
			//glVertex3f(p.x, p.y,p.z);
			glPushMatrix();
			glTranslated(p.x, p.y, p.z);
			glutSolidSphere(sph.sphereRadius1,6,6);
			glPopMatrix();
        }
			//Vector const &p = sph.particle(i).x;
			//glColor3f(sph.particle(i).color.R,sph.particle(i).color.G,sph.particle(i).color.B);
			//glVertex3f(p.x, p.y,p.z);
			//glPushMatrix();
			//glTranslated(p.x, p.y, p.z);
			//glutSolidSphere(sph.sphereRadius0,6,6);
			//glPopMatrix(); 
   		}
    //
		glEnd();
    }
    /*
    else
    {
    	//glPointSize(5);
     	//glEnable(GL_POINT_SMOOTH);
     	//glBegin(GL_POINTS);
      	// draw particles
      	for (int i = 0; i < sph.n_particles(); ++i) {
        	Vector const &p = sph.particle(i).x;
        	Particle &pi = sph.particles[i];
        	sph.NormalizedColor(sph.particle(i).T,&(pi.C));
        	glColor4f(sph.particle(i).C.R,sph.particle(i).C.G,sph.particle(i).C.B,sph.particle(i).C.A);
       		// glVertex3f(p.x, p.y,p.z);
	 		glPushMatrix();
	  		glTranslated(p.x, p.y,p.z);
			glutSolidSphere(sph.sphereRadius0,6,6);
			glPopMatrix(); 
      	}
       	//
   
   		//glEnd();
    }
  	// selected particles in red
   	//glPointSize(10);
    //glEnable(GL_POINT_SMOOTH);
    //glBegin(GL_POINTS);
  	glColor3f(1, 0, 0);
  	for (int i = 0; i < selected.size(); ++i) {
    	Vector const &p = sph.particle(selected[i].idx).x;
    	//std::cout << p << std::endl;
    	//glVertex3f(p.x, p.y,p.z);
    	glPushMatrix();
	  	glTranslated(p.x, p.y,p.z);
		glutSolidSphere(sph.sphereRadius0*2.0,6,6);
		glPopMatrix(); 
  	}
  	//glEnd();
  	*/
 
  	// if requested, draw search data structures
  	if (drawhash) {
    	//sph.hashgrid().draw(mouse_point);
  	}
  
  	if (drawtree) {
    	sph.kdtree().draw(mouse_point);
  	}
  
  	if (drawgrid) {
    	sph.surfaceCells(50);
    	sph.drawSurfaceGrid();
  	}
  
  	if (drawsurface) {
    	sph.surfaceCells(50);
    	sph.drawSurface();
  	}
  
  	glFinish();
  	glutSwapBuffers();
}


  

void reshape(int width, int height)
{
  glViewport(0, 0, width, height);
  sph.winWidth = width;
  sph.winHeight = height;
  
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
  gluUnProject(x, sph.winHeight-y, 0, MV, P, vp, &ox, &oy, &oz);
  
  oz = 0;
  Vector o(ox+0.5, oy+0.5, oz);

  return o;
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
      p.T+=0.1;
      float k = sph.vkernel.value(selected[i].d, selected[i].d_squared);
      p.a += mouse_force * k * (point - mouse_point);
      //std::cout << p << std::endl;
      //std::cout <<  " T = " << p.T << std::endl;
    }   
  }
  else if (buttondown[GLUT_RIGHT_BUTTON])
  {
  }
  
  mouse_point = point;
  Interactor->Refresh_List();
  glutPostRedisplay();
}

void button(int b, int state, int x, int y)
{
  ox = x;
  oy = y;
  
  // collect neighhoring particles
  if (state == GLUT_DOWN) {
    mouse_point = worldSpace(x,y);
    //std::cout << mouse_point << std::endl;
    sph.kdtree().neighbors(mouse_point, selected);
    //std::cout << "selected size = " << selected.size() << std::endl;
    buttondown[b] = true;
  } else {
    selected.clear();
    buttondown[b] = false;
  }
  Interactor->Refresh_List();
  glutPostRedisplay();
}


void HandleKeyboard(unsigned char key, int x, int y){	
	
	
	
	double coords[3];
	char *xs[10];
	allCommands->Keyboard(key);
	
	switch (key) {

		case 'e':
			exit(1);
		break;
	
	  case '\033':
    case 'f':
      saveImageFile = !saveImageFile;
      break;
      case 'i':
      sph.WriteScreenImage();
      break;
    case ' ':
      running = !running;
      break;
    case 't':
      drawtree = !drawtree;
      break;
    case 'h':
      drawhash = !drawhash;
      break;
    case 'x':
    {
       std::cout << " Tmin = " << sph.Tmin << " Tmax = " << sph.Tmax << std::endl;
    break;
    }
    case 's':
      drawsurface = !drawsurface;
      break;
    case 'g':
      drawgrid = !drawgrid;
      break;
    case 'r':
      // reset
      sph.init(0.02,0.02);
      break;
    case '.':
      sph.smoothing_radius(1.1 * sph.smoothing_radius());
      break;
    case ',':
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
	Interactor->Refresh_List();
	glutPostRedisplay();
}


void constructBox()
{ 
  double vertex[3];
  int cell[3];
  //face frontal
  
  vertex[0]=0.0;vertex[1]=0.0;vertex[2]=0.5;
  malha ->addVertex(vertex);
  vertex[0]=1.01;vertex[1]=0.0;vertex[2]=0.5;
  malha ->addVertex(vertex);
  vertex[0]=1.01;vertex[1]=1.0;vertex[2]=0.5;
  malha ->addVertex(vertex);
  vertex[0]=0.0;vertex[1]=1.0;vertex[2]=0.5;
  malha ->addVertex(vertex);
  cell[0]=0;cell[1]=1;cell[2]=3;
  //malha->addCell(cell);
  cell[0]=1;cell[1]=3;cell[2]=2	;
  //malha->addCell(cell);
  
  //chao
  
  vertex[0]=0.0;vertex[1]=0.0;vertex[2]=-0.5;
  malha ->addVertex(vertex);
  vertex[0]=1.01;vertex[1]=0.0;vertex[2]=-0.5;
  malha ->addVertex(vertex);
  cell[0]=0;cell[1]=4;cell[2]=1;
  malha->addCell(cell);
  cell[0]=1;cell[1]=4;cell[2]=5	;
  malha->addCell(cell);

//traz
  
  vertex[0]=0.0;vertex[1]=1.0;vertex[2]=-0.5;
  malha ->addVertex(vertex);
  vertex[0]=1.01;vertex[1]=1.0;vertex[2]=-0.5;
  malha ->addVertex(vertex);
  cell[0]=4;cell[1]=5;cell[2]=6;
  malha->addCell(cell);
  cell[0]=5;cell[1]=7;cell[2]=6	;
  malha->addCell(cell);

//lateral esquerda
  
  
  cell[0]=4;cell[1]=0;cell[2]=6;
  malha->addCell(cell);
  cell[0]=6;cell[1]=0;cell[2]=3	;
  malha->addCell(cell);
   
  
  //lateral direita
  
  
  cell[0]=1;cell[1]=5;cell[2]=2;
  malha->addCell(cell);
  cell[0]=2;cell[1]=5;cell[2]=7	;
  malha->addCell(cell);
  
}

void timerEvent(int value) {
	

	if (running) {
	
	  //glutPostRedisplay();
	  if (saveFile)
	     sph.surfaceCells(60);
	  sph.step();
	  if(saveImageFile)
	  {
	    sph.surfaceCells(60);
	    Interactor->WriteScreenImage();
	  }
	  Interactor->Refresh_List();
	  glutPostRedisplay();
	  }
	
glutTimerFunc(REFRESH_DELAY, timerEvent, 0);
}

void idle()
{
  if (running) {
    sph.step();
    if(saveImageFile)
    {
      sph.surfaceCells(50);
      Interactor->WriteScreenImage();
    }
    Interactor->Refresh_List();
    glutPostRedisplay();
  }
}
	
int main(int argc, char ** argv) 
{
  
    Interactor->setDraw(RenderScene);
    glutInit(&argc, argv);
  
  
    sph.init(0.02,0.02);
    
    malha = new TMesh();
   
    meshHandler.Set(malha);
    
    Print = new TPrintOf(meshHandler);
    
   	constructBox();
   
  	allCommands = new TMyCommands(Print, Interactor);
  	double a,x1,x2,y1,y2,z1,z2; 

	of::ofVerticesIterator<TTraits> iv(&meshHandler);

	iv.initialize();
	x1 = x2 = iv->getCoord(0);
	y1 = y2 = iv->getCoord(1);
	z1 = z2 = iv->getCoord(2);

	for(iv.initialize(); iv.notFinish(); ++iv){
		if(iv->getCoord(0) < x1) x1 = a = iv->getCoord(0);
		if(iv->getCoord(0) > x2) x2 = a = iv->getCoord(0);
		if(iv->getCoord(1) < y1) y1 = a = iv->getCoord(1);
		if(iv->getCoord(1) > y2) y2 = a = iv->getCoord(1);
		if(iv->getCoord(2) < z1) z1 = a = iv->getCoord(2);
		if(iv->getCoord(2) > z2) z2 = a = iv->getCoord(2);
	}

	double maxdim;
	maxdim = fabs(x2 - x1);
	if(maxdim < fabs(y2 - y1)) maxdim = fabs(y2 - y1);
	if(maxdim < fabs(z2 - z1)) maxdim = fabs(z2 - z1);

	maxdim *= 0.5;
	
	Point center((x1+x2)/2.0, (y1+y2)/2.0, (z1+z2)/2.0 );
	Interactor->Init(center[0]-maxdim, center[0]+maxdim,
					center[1]-maxdim, center[1]+maxdim,
					center[2]-maxdim, center[2]+maxdim);
	
	

 	AddKeyboard(HandleKeyboard);
 	//allCommands->Help(std::cout);
	std::cout<< std::endl<< "Press \"?\" key for help"<<std::endl<<std::endl;
 
  	glutTimerFunc(REFRESH_DELAY, timerEvent, 0);
  	glutMotionFunc(motion);
  	glutPassiveMotionFunc(motion);
 	glutMouseFunc(button);
  	Init_Interactor();
	// Interactor->Refresh_List();
  
  return 0;
}
