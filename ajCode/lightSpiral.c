#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#if defined(__APPLE__) || defined(MACOSX)
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#include <GL/glu.h>
#include <GL/gl.h>
#endif

/* Function Prototypes */
void init( void );
void display( void );
void keyboard( unsigned char, int, int );
void reshape( int, int );
void drawText( void );
void buildMesh( void );
void drawMesh( void );
void sphereToCart( GLdouble, GLdouble, GLdouble, double *, double *, double * );
void normalize( GLdouble * );
void mouseMotion( int, int );
void mouse( int, int, int , int );
void idle( void );

/* Spine Curve Attributes */
/* Number of points to be calculated */
#define TUBE_N 512
#define TUBE_M 64
#define RADIANS_PER_PIXEL 0.01
#define EPS 0.01
/* indices stuff */
#define INDICES_PER_STRIP (2*TUBE_M+2)
#define NUM_STRIPS TUBE_N
#define NUM_STRIP_INDICES (NUM_STRIPS*INDICES_PER_STRIP)
#define VERTEX_INDEX(i,j) ((i)*TUBE_M + (j))
GLushort tubeStrips[NUM_STRIP_INDICES];
/* remaining */
GLdouble tube[TUBE_N][TUBE_M][3];
GLdouble tubeNormals[TUBE_N][TUBE_M][3];
GLdouble a = 1.0; // Large radii of torus
GLdouble b = 0.45; //small radii of torus
GLdouble R = 0.25; // Tube radii
GLdouble p = 1; //Times the curve winds around the origin
GLdouble q = 3; //number of times the curve coils around the torus
/* Vertex buffer handler */
GLuint vertexBuffer;
GLuint indexBuffer;
GLuint normalBuffer;
/* Eye stuff */
GLdouble eyeRho = 3;
GLdouble eyeTheta = 0;
GLdouble eyePhi = EPS;
int mousex, mousey;

float seconds = 0.0;

void buildTube(){
  static GLboolean firstTime = GL_TRUE;
  if( firstTime ) {
    glGenBuffers(1,&vertexBuffer);
    glGenBuffers(1,&indexBuffer);
    glGenBuffers(1,&normalBuffer);
    /* Fill and buffer index data, since its done only once */
    int n = 0;
    int i, j;
    for ( i = 0; i < TUBE_N; i++ ) {
      for( j = 0; j < TUBE_M; j++ ) {
	tubeStrips[n++] = VERTEX_INDEX((i+1)%TUBE_N,j);
 	tubeStrips[n++] = VERTEX_INDEX(i,j);
      }
      tubeStrips[n++] = VERTEX_INDEX((i+1)%TUBE_N,0);
      tubeStrips[n++] = VERTEX_INDEX(i,0);
    }
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexBuffer);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(tubeStrips), tubeStrips, GL_STATIC_DRAW);

    firstTime = GL_FALSE;
    
  }
  double dt = (2*M_PI)/TUBE_N; 
  double du = (2*M_PI)/TUBE_M;
  /* Calculate Points */
  double t;
  int i;
  for(i = 0, t = 0.0; i < TUBE_N; i++, t += dt ){
    GLdouble C[3], T[3], A[3], B[3], N[3];
    /* x, y, z */
    C[0] = ( a + b*cos(q*t)) * cos(p*t);
    C[1] = ( a + b*cos(q*t)) * sin(p*t);
    C[2] = b * sin(q*t); 
    /* dx, dy, dz */
    T[0] = -p*C[1] - b*q*sin(q*t)*cos(p*t);
    T[1] = p*C[0] - b*q*sin(q*t)*sin(p*t);
    T[2] = b*q*cos(q*t);
    /* ddx, ddy, ddz */
    A[0] = -p*T[1] + b*q*(p*sin(q*t)*sin(p*t) - q*cos(q*t)*cos(p*t));
    A[1] = p*T[0] - b*q*(p*sin(q*t)*cos(p*t) + q*cos(q*t)*sin(p*t));
    A[2] = -(q*q)*b*sin(q*t);
    /* B = T cross A */
    B[0] = T[1]*A[2] - T[2]*A[1];
    B[1] = T[2]*A[0] - T[0]*A[2];
    B[2] = T[0]*A[1] - T[1]*A[0];
    /* Normalize T and B */
    normalize(T); 
    normalize(B);
    /* N = B cross T */
    N[0] = B[1]*T[2] - B[2]*T[1];
    N[1] = B[2]*T[0] - B[0]*T[2];
    N[2] = B[0]*T[1] - B[1]*T[0];
    //normalize(N);
    float u;
    int j;
    for( j = 0,  u = 0.0; j < TUBE_M; j++, u += du ){
      for( int k = 0; k < 3; k++ ){
	tubeNormals[i][j][k] = cos(u)*B[k] + sin(u)*N[k];
	tube[i][j][k] = C[k] + R*tubeNormals[i][j][k];
      }
    }
  }
  /* Load vertex buffer */
  glBindBuffer(GL_ARRAY_BUFFER, normalBuffer);
  glBufferData(GL_ARRAY_BUFFER, sizeof(tubeNormals), tubeNormals, GL_DYNAMIC_DRAW);
  glBindBuffer(GL_ARRAY_BUFFER, vertexBuffer);
  glBufferData(GL_ARRAY_BUFFER, sizeof(tube), tube, GL_DYNAMIC_DRAW);

}

void drawTube() {
  glEnableClientState( GL_VERTEX_ARRAY );
  glEnableClientState( GL_NORMAL_ARRAY );
  glBindBuffer(GL_ARRAY_BUFFER, vertexBuffer);
  glVertexPointer( 3, GL_DOUBLE, 0, (GLvoid*) 0);
  glBindBuffer(GL_ARRAY_BUFFER, normalBuffer);
  glNormalPointer(GL_DOUBLE, 0, (GLvoid*) 0 );
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexBuffer);

  char * offset;
  int i;
  for( i = 0, offset = 0; i < NUM_STRIPS; 
       i++, offset += INDICES_PER_STRIP*sizeof(GLushort))
    glDrawElements(GL_QUAD_STRIP, INDICES_PER_STRIP, GL_UNSIGNED_SHORT,
		   (GLvoid*) offset);
  
}



void display( ) {
  glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  /* Fixed light Pos */
  GLfloat position[] = {15*sin(0.5*seconds), 0 , 15*cos(0.5*seconds)};
  glLightfv(GL_LIGHT0, GL_POSITION, position);
  /* Set up gluLookAt */
  double eyex, eyey, eyez;
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  sphereToCart(eyeTheta, eyePhi, eyeRho, &eyex, &eyey, &eyez);
  gluLookAt(eyex, eyey, eyez, 0,0,0, 0,0,1);

  drawTube();
 
  glutSwapBuffers();
}

void keyboard( unsigned char key, int x, int y ){
  GLint ESC = 27;
  if( key == ESC || key == 'q' || key == 'Q' ) exit(0);
  if( key == 'f' ) { p++; buildTube(); }
  if( key == 'g' ) { p--; buildTube(); }
  if( key == 'v' ) { q++; buildTube(); }
  if( key == 'b' ) { q--; buildTube(); }
  glutPostRedisplay();
}

void mouse( int button, int state, int x, int y){
  mousex = x;
  mousey = y;
}

void mouseMotion( int x, int y ){
  int dx = x - mousex, dy = y - mousey;
  //double eyex, eyey, eyez;
  eyeTheta -= dx*RADIANS_PER_PIXEL;
  eyePhi -= dy*RADIANS_PER_PIXEL;
  if( eyePhi < EPS ) eyePhi = EPS;
  if( eyePhi > (M_PI - EPS) ) eyePhi = M_PI - EPS;
  mousex = x; mousey = y;
  /*  
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  sphereToCart(eyeTheta, eyePhi, eyeRho, &eyex, &eyey, &eyez);
  gluLookAt(eyex, eyey, eyez, 0,0,0, 0,0,1);
  */
  glutPostRedisplay();
}

void reshape( int w, int h ){
  glViewport( 0, 0, w, h );
  glMatrixMode( GL_PROJECTION );
  glLoadIdentity();
  gluPerspective( 80, ((double)(w)/(double)(h)), 0.1, 20 );
}

int main( int argc, char ** argv ){
	
  glutInit( &argc, argv );
  glutInitDisplayMode( GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH );
  glutInitWindowSize( 1000, 600 );
  glutInitWindowPosition( 10, 10 );
  glutCreateWindow("Toroidal Spiral");

  glutMotionFunc( mouseMotion );
  glutMouseFunc( mouse );
  glutDisplayFunc(display);
  glutReshapeFunc(reshape);
  glutKeyboardFunc(keyboard);
  glutIdleFunc( idle );
	
  init();
	
  glutMainLoop();
	
  return 0;
}

void idle () {
  seconds = ((float)glutGet(GLUT_ELAPSED_TIME))/1000.0;
  glutPostRedisplay();
}

GLfloat globalAmbient[] = {0.1, 0.1, 0.1, 1.0 };
GLfloat ambient[] = {0.01, 0.01, 0.01 };
GLfloat diffuse[] = {0.6, 0.6, 0.6};
GLfloat specular[] = {0.5, 0.5, 0.5};
GLfloat shininess[] = {30};

void init(){
  glClearColor(0.0,0.0,0.0,1.0);
  glFrontFace(GL_CW);
  glCullFace(GL_BACK);
  glEnable(GL_CULL_FACE);
  glEnable(GL_DEPTH_TEST);
  glShadeModel(GL_SMOOTH);
  glEnable( GL_LIGHTING);
  glEnable( GL_LIGHT0 );
  
  glLightModelfv(GL_LIGHT_MODEL_AMBIENT, globalAmbient);
  glLightfv(GL_LIGHT0, GL_AMBIENT, ambient);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse);
  glLightfv(GL_LIGHT0, GL_SPECULAR, specular);
  glLightfv(GL_LIGHT0, GL_SHININESS, shininess);
  GLfloat mColor[] = {0.0,1.0,1.0};
  glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, mColor);
  GLfloat spec[] = {0.4,0.4,0.4};
  glMaterialfv(GL_FRONT, GL_SPECULAR, spec);
  glMaterialfv(GL_FRONT, GL_SHININESS, shininess);

  buildTube();
}

void sphereToCart( GLdouble theta, GLdouble phi, GLdouble rho, 
		   double * x, double * y, double * z ){
  *x = rho*sin(phi)*cos(theta);
  *y = rho*sin(phi)*sin(theta);
  *z = rho*cos(phi);
}

void normalize( GLdouble * v ){
  double mag = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
  v[0] /= mag;
  v[1] /= mag;
  v[2] /= mag;
}
