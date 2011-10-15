#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#if defined(__APPLE__) || defined(MACOSX)
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#else
#define GL_GLEXT_PROTOTYPES
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define TUBE_N 64
#define TUBE_M 32

static float a = 100.0, b = 40.0;
static float R = 20.0;
static int p = 1, q = 7;
GLfloat xOf(float t){
    return ((a+b*cos(q*t))*cos(p*t));
}

GLfloat yOf(float t){
    return ((a+b*cos(q*t))*sin(p*t));
}

GLfloat zOf(float t){
    return (b*sin(q*t));
}

GLfloat dxOf(float t){
    return -p * yOf(t) - b*q*sin(q*t)*cos(p*t);
}

GLfloat dyOf(float t){
    return p * xOf(t) - b*q *sin(q*t) *sin(p*t);
}

GLfloat dzOf(float t){
    return b * q * cos(q*t);
}

GLfloat ddxOf(float t){
    return -p * dyOf(t) + b *q * (p * sin(q*t) * sin(p*t) - q*cos(q*t)*cos(p*t));
}

GLfloat ddyOf(float t){
    return p * dxOf(t) - b*q*(p*sin(q*t)*cos(p*t)+q*cos(q*t)*sin(p*t));
}

GLfloat ddzOf(float t){
    return -(q*q)*b*sin(q*t);
}

void spine(){
    static GLboolean first = GL_TRUE;
    static GLuint vertBuffer;
    if (first) {
        GLfloat* verts = malloc(sizeof(GLfloat)*TUBE_N*3);
        for(int i = 0; i < TUBE_N; i++){
            float t = 2.0*M_PI*(float)i / (float)TUBE_N;
            verts[3*i] = xOf(t);
            verts[3*i+1] = yOf(t); 
            verts[3*i+2] = zOf(t); 
            //printf("%f, %f, %f, %f\n", t, verts[3*i], verts[3*i+1], verts[3*i+2]);
        }
        glGenBuffers(1, &vertBuffer);
        glBindBuffer(GL_ARRAY_BUFFER, vertBuffer);
        glBufferData(GL_ARRAY_BUFFER, TUBE_N*sizeof(GLfloat)*3,
                verts, GL_STATIC_DRAW);
        first = GL_FALSE;
    }
    glEnableClientState(GL_VERTEX_ARRAY);
    glBindBuffer(GL_ARRAY_BUFFER, vertBuffer); 
    glVertexPointer(3, GL_FLOAT, 0, (GLvoid*) 0);
    glDrawArrays(GL_LINE_LOOP, 0, TUBE_N);
}

void cross(GLdouble A[3], GLdouble B[3], GLdouble C[3]){
    GLdouble a1 = A[0], a2 = A[1], a3 = A[2];
    GLdouble b1 = B[0], b2 = B[1], b3 = B[2];
    C[0] = a2*b3-a3*b2;
    C[1] = a3*b1-a1*b3;
    C[2] = a1*b2-a2*b1;
}

void normalize(GLdouble A[3]){
    GLdouble mag = sqrt(A[0]*A[0]+A[1]*A[1]+A[2]*A[2]);
    //printf("%f\n", mag);
    A[0] = A[0]/mag;
    A[1] = A[1]/mag;
    A[2] = A[2]/mag;
}

GLdouble tube[TUBE_N][TUBE_M][3];
GLdouble tubeNormals[TUBE_N][TUBE_M][3];
void mesh(void){
    float dt = 2*M_PI/TUBE_N;
    float du = 2*M_PI/TUBE_M;
    float t = 0.0;
    int i;
    for (i = 0, t = 0.0; i < TUBE_N; i++, t += dt) {
        GLdouble C[3], T[3], A[3], B[3], N[3];
        C[0] = xOf(t), C[1] = yOf(t), C[2] = zOf(t);
        T[0] = dxOf(t), T[1] = dyOf(t), T[2] = dzOf(t);
        A[0] = ddxOf(t), A[1] = ddyOf(t), A[2] = ddzOf(t);
        cross(T,A,B);
        normalize(T);
        normalize(B);
        cross(B,T,N);
        int j, k; float u;
        for (j = 0, u = 0.0; j < TUBE_M; j++, u += du){
            for (k = 0; k < 3; k++) {
                tubeNormals[i][j][k] = cos(u)*B[k] + sin(u)*N[k];
                tube[i][j][k] = C[k] + R*tubeNormals[i][j][k];
                //printf("%f ", tubeNormals[i][j][k]);
            }
            //printf("\n");
        }
    }
}

#define INDICES_PER_STRIP (2*TUBE_M+2)
#define NUM_STRIPS TUBE_N
#define NUM_STRIP_INDICES (NUM_STRIPS*INDICES_PER_STRIP)
GLushort tubeStrips[NUM_STRIP_INDICES];
#define VERTEX_INDEX(i,j) ((i)*TUBE_M + (j))
void createMeshStripIndices(void) {
    int n = 0;
    int i,j;
    for (i = 0; i < TUBE_N; i++) {
        for (j = 0; j < TUBE_M; j++) {
            tubeStrips[n++] = VERTEX_INDEX((i+1)%TUBE_N,j);
            tubeStrips[n++] = VERTEX_INDEX(i,j);
        }
        tubeStrips[n++] = VERTEX_INDEX((i+1)%TUBE_N,0);
        tubeStrips[n++] = VERTEX_INDEX(i,0);
    }
}

void drawMesh(){
    static GLuint vertBuffer;
    static GLuint indexBuffer;
    static GLuint normBuffer;

    glGenBuffers(1, &vertBuffer);
    glBindBuffer(GL_ARRAY_BUFFER, vertBuffer);
    glBufferData(GL_ARRAY_BUFFER, sizeof(tube), tube, GL_STATIC_DRAW);

    glGenBuffers(1, &indexBuffer);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexBuffer);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(tubeStrips), tubeStrips, GL_STATIC_DRAW);

    glGenBuffers(1, &normBuffer);
    glBindBuffer(GL_ARRAY_BUFFER, normBuffer);
    glBufferData(GL_ARRAY_BUFFER, sizeof(tubeNormals), tubeNormals, GL_STATIC_DRAW);

    glEnableClientState(GL_NORMAL_ARRAY);
    glBindBuffer(GL_ARRAY_BUFFER, normBuffer); 
    glNormalPointer(GL_FLOAT, 0, (GLvoid*)normBuffer);

    glEnableClientState(GL_VERTEX_ARRAY);
    glBindBuffer(GL_ARRAY_BUFFER, vertBuffer); 
    glVertexPointer(3, GL_DOUBLE, 0, (GLvoid*) 0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexBuffer); 
    int i;
    char *offset;
    for (i = 0, offset = 0; i < NUM_STRIPS; i++, offset += INDICES_PER_STRIP*sizeof(GLushort)){
        glDrawElements(GL_QUAD_STRIP, INDICES_PER_STRIP, GL_UNSIGNED_SHORT, (GLvoid*) offset);
    }
}

void reshape(int w, int h) {
#define HALF_SIZE 18.0 /* window half-size */
    glViewport(0,0, w,h); /* viewport is entire application window */
    //glMatrixMode(GL_PROJECTION); /* set window into world coord system */
    //glLoadIdentity();
    if (w > h)
        gluOrtho2D(-HALF_SIZE, HALF_SIZE, -HALF_SIZE*h/w, HALF_SIZE*h/w);
    else
        gluOrtho2D(-HALF_SIZE*w/h, HALF_SIZE*w/h, -HALF_SIZE, HALF_SIZE);
    //glMatrixMode(GL_MODELVIEW); /* modelview matrix = identity */
    //glLoadIdentity();
}
GLfloat theta = 0.0;
float lmodel_ambient[4]    = {0.1, 0.1, 0.1, 1.0};
GLfloat light0_position[4]   = {0.0, 0.0, 500.0, 0.0};
GLfloat light0_ambient[4]    = {0.3, 0.3, 0.3, 1.0};
GLfloat light0_diffuse[4]    = {1.0, 1.0, 1.0, 1.0};
GLfloat light0_specular[4]   = {1.0, 1.0, 1.0, 1.0};


void display(void) {
    static int is_first = 1;
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
    glColor3f(1, 0, 0);
    glLightfv(GL_LIGHT0, GL_POSITION, light0_position);
    //spine();
    mesh();

    if(is_first){
        createMeshStripIndices();
        is_first = 0;
    }
    drawMesh();
    glutSwapBuffers();
}

GLdouble eye[3] = {1.0, 1.0, 500.0};
GLdouble lookat[3] = {0.0, 0.0, 0.0};
GLdouble up[3] = {0.0, 0.0, 1.0};
GLdouble hither = 0.1, yon = 1000.0;

void sphericalToCartesian(double r, double theta, double phi,
        double *x, double *y, double *z) {
    double sin_phi = sin(phi);
    *x = r*cos(theta)*sin_phi;
    *y = r*sin(theta)*sin_phi;
    *z = r*cos(phi);
}

GLdouble eyeRadius = 500, eyePhi = M_PI/4.0, eyeTheta = M_PI/4.0;
void setView() {
    sphericalToCartesian(eyeRadius, eyeTheta, eyePhi, &eye[0], &eye[1], &eye[2]);
    eye[0] += lookat[0]; eye[1] += lookat[1]; eye[2] += lookat[2];
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glLightfv(GL_LIGHT0, GL_POSITION, light0_position);
    gluLookAt(eye[0], eye[1], eye[2], lookat[0], lookat[1], lookat[2], up[0], up[1], up[2]);
}



GLfloat ambient[4] = {0.24725, 0.1995, 0.0745, 1.0};
GLfloat diffuse[4] = {0.75164, 0.60648, 0.22648, 1.0};
GLfloat specular[4] = {0.628281, 0.555802, 0.366025, 1.0};
GLfloat shininess = 51.2;


void init() {
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lmodel_ambient);
    glLightfv(GL_LIGHT0, GL_AMBIENT, light0_ambient);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, light0_diffuse);
    glLightfv(GL_LIGHT0, GL_SPECULAR, light0_specular);
    glLightfv(GL_LIGHT0, GL_POSITION, light0_position);
    glEnable(GL_LIGHT0);

    glMaterialfv(GL_FRONT, GL_AMBIENT, ambient);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, diffuse);
    glMaterialfv(GL_FRONT, GL_SPECULAR, specular);
    glMaterialfv(GL_FRONT, GL_SHININESS, &shininess);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(40.0, 1, hither, yon);

    setView();

    glEnable(GL_DEPTH_TEST);

    glFrontFace(GL_CW);
    glCullFace(GL_BACK);
    glDisable(GL_CULL_FACE);
    glEnable(GL_CULL_FACE);

    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    //glPolygonMode(GL_FRONT, GL_LINE);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glEnable(GL_LIGHTING);

    glEnable(GL_DEPTH_TEST);
    glShadeModel(GL_SMOOTH);

    glColor3f(1,1,1);

}





bool mouseRotate = false;
int mousex, mousey;
const double epsilon = 0.0000001;

void mouse(int button, int state, int x, int y) {
    if (state == GLUT_DOWN) {
        if (button == GLUT_LEFT_BUTTON) {
            mouseRotate = true;
            mousex = x;
            mousey = y;
        }
    } else if (state == GLUT_UP) {
        mouseRotate = false;
    }
}

void mouseMotion(int x, int y) {
    if (mouseRotate) {
        const double radiansPerPixel = M_PI/(2*90.0);
        int dx = x - mousex, dy = y - mousey;
        eyeTheta -= dx*radiansPerPixel;
        eyePhi -= dy*radiansPerPixel;
        if (eyePhi >= M_PI)
            eyePhi = M_PI - epsilon;
        else if (eyePhi <= 0.0)
            eyePhi = epsilon;
        setView();
        mousex = x;
        mousey = y;
        glutPostRedisplay();
    }
}

int main(int argc, char *argv[]) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH); /* double buffering */
    glEnable(GL_DEPTH_TEST);
    glutInitWindowSize(800, 800);
    glutInitWindowPosition(10, 10);
    glutCreateWindow("Dancing Santa Robot");
    glutMouseFunc(mouse);
    glutMotionFunc(mouseMotion);
    //glutReshapeFunc(reshape);
    glClearColor(1.0f,1.0f,1.0f,1.0f);
    glutDisplayFunc(display);
    init();
    glutMainLoop();
    return 0;
}
