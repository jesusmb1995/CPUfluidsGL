/*
 * Copyright 1993-2015 NVIDIA Corporation.  All rights reserved.
 *
 * Please refer to the NVIDIA end user license agreement (EULA) associated
 * with this source code for terms and conditions that govern your use of
 * this software. Any use, reproduction, disclosure, or distribution of
 * this software and related documentation outside the terms of the EULA
 * is strictly prohibited.
 *
 */

// This is a ported CPU version of the CUDA fluidsGL demo
// originally published by NVidia.
//
// Port made by Jesús Martín Berlanga
//
// Info and documentation at rtfluids.bdevel.org

// OpenGL Graphics includes
#if defined(__APPLE__) || defined(MACOSX)
  #pragma clang diagnostic ignored "-Wdeprecated-declarations"
  #include <GLUT/glut.h>
  #ifndef glutCloseFunc
  #define glutCloseFunc glutWMCloseFunc
  #endif
#else
#include <GL/glew.h>
#include <GL/freeglut.h>
#endif

// Includes
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "defines.h"
#include "fluidsGL_cpu.hpp"

#define MAX_EPSILON_ERROR 1.0f

const char *sSDKname = "fluidsGL";

void cleanup(void);
void reshape(int x, int y);

static cData *vxfield = NULL;
static cData *vyfield = NULL;

cData *hvfield = NULL;
static int wWidth  = MAX(512, DIM);
static int wHeight = MAX(512, DIM);

// FPS counter
static int fpsCount = 0;
static int time = 0;
static int base_time = 0;
static int timeLimit = 1000;
static float fps = 0;

// Variables for the automatic stress test
static bool stress_test = false;
static int stress_test_step_index = 0;
static int stress_test_n_steps = 10000;
static int fps_n_accumulations = 0;
static float fps_accumulation = 0;

static int clicked  = 0;

// Particle data
GLuint vbo = 0;                 // OpenGL vertex buffer object
static cData *particles = NULL; // particle positions in host memory
static int lastx = 0, lasty = 0;

// Texture pitch
size_t tPitch = 0; // Now this is compatible with gcc in 64-bit

char *ref_file         = NULL;
bool g_bQAAddTestForce = true;
int  g_TotalErrors     = 0;

// 10 vueltas de bucle en el test
int  g_dumpFrames = 10;

//  Para el test:
// No recorrer matrizes enteras, 
// nos llevaría demasiado tiempo,
// recorrer 40 veces menos (saltando 
// de 40 en 40).
int  g_compareDivisor  = 1;

bool g_bExitESC = false;

//#define DUMPEABLE
#ifdef DUMPEABLE
FILE* dump_file = NULL;
void dump_title(const char* title) {
    if(dump_file != NULL) {
        fprintf(dump_file, 
"\n\n\n///////////////////////////////////[%s]///////////////////////////////////\n"
        ,title);
    }
}
void dump(const char* tag) {
    if(dump_file != NULL) {
        fprintf(dump_file, 
"\n\n=============[%s]=============\n"
        ,tag);

        // Dump de hvfield
        fprintf(dump_file,"\n----[vfield]----\n");
        for(int x = 0; x < DIM; x = x + g_compareDivisor) {
            for(int y = 0; y < DIM; y = y + g_compareDivisor) {
                int i = y*DIM + x;
                cData iv = hvfield[i];
                fprintf(dump_file,"(%f,%f)",iv.x,iv.y);  
            }
            fprintf(dump_file,"\n");
        }

        // Dump de vxfield
        fprintf(dump_file,"\n----[vxfield memory]----\n");
        for(int x = 0; x < CPADW; x = x + g_compareDivisor) {
            for(int y = 0; y < DIM; y = y + g_compareDivisor) {
                int i = y*CPADW + x;
                cData vx = vxfield[i];
                fprintf(dump_file,"|%f|%f",
                    vx.x,vx.y
                );  
            }
            fprintf(dump_file,"\n");
        }

        // Dump de vyfield
        fprintf(dump_file,"\n----[vyfield memory]----\n");
        for(int x = 0; x < CPADW; x = x + g_compareDivisor) {
            for(int y = 0; y < DIM; y = y + g_compareDivisor) {
                int i = y*CPADW + x;
                cData vy = vyfield[i];
                fprintf(dump_file,"|%f|%f",
                    vy.x,vy.y
                );  
            }
            fprintf(dump_file,"\n");
        }
    }
}
#endif

void gl_print_errors(const char* where) {
	GLenum err;
	while ((err = glGetError()) != GL_NO_ERROR) {
		fprintf(stderr, "OpenGL error at %s: %#04x\n", where, err);
	}
}

cData* glBufferMap;
#ifdef DUMPEABLE
void dumpMapGLBuffer() {
	fprintf(dump_file, "\n''glBufferMap mapGLBuffer()''''''''''''''''''''\n");
        for(int x = 0; x < DIM; x = x + g_compareDivisor) {
            for(int y = 0; y < DIM; y = y + g_compareDivisor) {
                int i = y*DIM + x;
                cData part = glBufferMap[i];
                fprintf(dump_file,"(%f,%f)",part.x,part.y);  
            }
            fprintf(dump_file,"\n");
        }
        fprintf(dump_file, "\n'''''''''''''''''''''''''''''''''''\n");
}
#endif
void mapGlBuffer() {
	glBindBuffer(GL_ARRAY_BUFFER, vbo); // Establish to which buffer we will do the mapping
        // Map buffer with:
        //   void *glMapBufferRange​(GLenum target​, GLintptr offset​,
        //                    GLsizeiptr length​, GLbitfield access​);
        glBufferMap = (cData*) glMapBufferRange(GL_ARRAY_BUFFER, 0,
        		sizeof(cData) * DS, GL_MAP_WRITE_BIT | GL_MAP_READ_BIT);
        if(glBufferMap == NULL) {
	    fprintf(stderr, "Failed to map GL_ARRAY_BUFFER\n");
	    gl_print_errors("glBufferMap");
            exit(EXIT_FAILURE);
        }
	glBindBuffer(GL_ARRAY_BUFFER, 0); // Unbind
   	#ifdef DUMPEABLE 
	dumpMapGLBuffer();
   	#endif
}

void unmapGLBuffer() {
	#ifdef DUMPEABLE 
	dumpMapGLBuffer();
	#endif
	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	// Unmap buffer if we will no longer operate with it
        // in order to evade unexpected errors.
	// For example, if a buffer is used for feedback at the
        // same time it's used for another purpose it will cause
        // an error. e.g. mapped buffer at the same time we call
        // glDrawArrays will cause an error.
	glUnmapBuffer(GL_ARRAY_BUFFER);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
}

extern "C" void addForces(cData *v, int dx, int dy, int spx, int spy, float fx, float fy, int r);
extern "C" void advectVelocity(cData *v, float *vx, float *vy, int dx, int pdx, int dy, float dt);
extern "C" void diffuseProject(cData *vx, cData *vy, int dx, int dy, float dt, float visc);
extern "C" void updateVelocity(cData *v, float *vx, float *vy, int dx, int pdx, int dy);
extern "C" void advectParticles(cData *part, cData *v, int dx, int dy, float dt);

void simulateFluids(void)
{
    // simulate fluid
    #ifdef DUMPEABLE 
    dump_title("simulateFluids");
    #endif

    #ifdef DUMPEABLE 
    dump("advectVelocity (pre)");
    #endif
    advectVelocity(hvfield, (float *)vxfield, (float *)vyfield, DIM, RPADW, DIM, DT);
    #ifdef DUMPEABLE 
    dump("advectVelocity (post)");
    #endif

    #ifdef DUMPEABLE 
    dump("diffuseProject (pre)");
    #endif
    diffuseProject(vxfield, vyfield, CPADW, DIM, DT, VIS);
    #ifdef DUMPEABLE 
    dump("diffuseProject (post)");
    #endif

    #ifdef DUMPEABLE 
    dump("updateVelocity (pre)");
    #endif
    updateVelocity(hvfield, (float *)vxfield, (float *)vyfield, DIM, RPADW, DIM);
    #ifdef DUMPEABLE 
    dump("updateVelocity (post)");
    #endif

    mapGlBuffer();
    #ifdef DUMPEABLE 
    dump("advectParticles (pre)");
    #endif
    advectParticles(glBufferMap, hvfield, DIM, DIM, DT);
    #ifdef DUMPEABLE 
    dump("advectParticles (post)");
    #endif
    unmapGLBuffer();
}

void display(void)
{

    if (!ref_file)
    {
        simulateFluids();
    }

    // render points from vertex buffer
    glClear(GL_COLOR_BUFFER_BIT);
    glColor4f(0,1,0,0.5f);
    glPointSize(1);
    glEnable(GL_POINT_SMOOTH);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnableClientState(GL_VERTEX_ARRAY);
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_CULL_FACE);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glVertexPointer(2, GL_FLOAT, 0, NULL);
    gl_print_errors("glVertexPointer");
    glDrawArrays(GL_POINTS, 0, DS);
    gl_print_errors("glDrawArrays");
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glDisableClientState(GL_VERTEX_ARRAY);
    glDisableClientState(GL_TEXTURE_COORD_ARRAY);
    glDisable(GL_TEXTURE_2D);

    if (ref_file)
    {
        return;
    }

    if(stress_test && stress_test_step_index < stress_test_n_steps) {
        if(stress_test_step_index % 10 == 0) {
            // Add Force (simulate click)
            // each 10 steps
            int x = wWidth/(stress_test_step_index/100+1);
            int y = wHeight/(stress_test_step_index/100+1);
            float fx = (x / (float)wWidth);
            float fy = (y / (float)wHeight);
            int nx = (int)(fx * DIM);
            int ny = (int)(fy * DIM);

            int ddx = 35;
            int ddy = 35;
            fx = ddx / (float)wWidth;
            fy = ddy / (float)wHeight;
            int spy = ny-FR;
            int spx = nx-FR;

            addForces(hvfield, DIM, DIM, spx, spy, FORCE * DT * fx, FORCE * DT * fy, FR);

            lastx = x;
            lasty = y;
        }
        stress_test_step_index++;
    }

    fpsCount++;
    time = glutGet(GLUT_ELAPSED_TIME);
    if ((time - base_time) > timeLimit)
    {
	// multiply by 1000 to convert ms to seconds
        fps=fpsCount*1000.0f/(time - base_time);
        // update last taken fps time (base_time)
	base_time = time;
        // reset fps count
	fpsCount=0;
	char fps_string[256];	
	sprintf(fps_string, "CPU/GL Stable Fluids (%d x %d): %f fps", DIM, DIM, fps);  	
 	glutSetWindowTitle(fps_string);        

    	if(stress_test) {
       		printf("%f fps (%d/%d)\n",fps,stress_test_step_index,stress_test_n_steps);
		fps_accumulation += fps;
		fps_n_accumulations++;
                
                // Exit test just right after
                // the last measure.
                if(stress_test_step_index >= stress_test_n_steps) {
                    printf("fps average: %f\n",fps_accumulation/fps_n_accumulations);
                    g_bExitESC = true;
                    #if defined (__APPLE__) || defined(MACOSX)
                        exit(EXIT_SUCCESS);
                    #else
                        glutDestroyWindow(glutGetWindow());
                        return;
                    #endif
                }
    	}
    }

    // Finish timing before swap buffers to avoid refresh sync
    glutSwapBuffers();

    glutPostRedisplay();
}

/*
 * Escribe en un fichero auxiliar los resultados
 * de cada uno de los pasos del algoritmo.
 */
void dumpTest()
{    
    #ifdef DUMPEABLE
    dump_file = fopen("dump_test/dump.txt", "w");
    #endif
    reshape(wWidth, wHeight);

    for (int count=0; count<g_dumpFrames; count++)
    {
        simulateFluids();

        // add in a little force so the automated testing is interesting.
        if (ref_file)
        {
	    #ifdef DUMPEABLE 
	    dump("addForces (pre)");
	    #endif
            int x = wWidth/(count+1);
            int y = wHeight/(count+1);
            float fx = (x / (float)wWidth);
            float fy = (y / (float)wHeight);
            int nx = (int)(fx * DIM);
            int ny = (int)(fy * DIM);
	    #ifdef DUMPEABLE 
	    fprintf(dump_file,"\n----[dumpTest automatic forces]----\n");
            fprintf(dump_file, "\tx = %d/(%d+1) = %d",wWidth,count,x);
	    fprintf(dump_file, "\ty = %d/(%d+1) = %d",wHeight,count,y);
	    fprintf(dump_file, "\tfx = %d/((float)%d) = %f",x,wWidth,fx);	
	    fprintf(dump_file, "\tfy = %d/((float)%d) = %f",y,wHeight,fy);
            #endif

            int ddx = 35;
            int ddy = 35;
            fx = ddx / (float)wWidth;
            fy = ddy / (float)wHeight;
            int spy = ny-FR;
            int spx = nx-FR;

            #ifdef DUMPEABLE
	    fprintf(dump_file, "\tFORCE * DT * fx = %f * %f * fx = %f",FORCE,DT,FORCE * DT * fx); 
	    fprintf(dump_file, "\tFORCE * DT * fx = %f * %f * fx = %f\n",FORCE,DT,FORCE * DT * fy);
            #endif
            addForces(hvfield, DIM, DIM, spx, spy, FORCE * DT * fx, FORCE * DT * fy, FR);
	    #ifdef DUMPEABLE 
	    dump("addForces (post)");
	    #endif

            lastx = x;
            lasty = y;
        }
    }

    display();
    #ifdef DUMPEABLE
    fclose(dump_file);
    #endif
}

// very simple von neumann middle-square prng.  can't use rand() in -qatest
// mode because its implementation varies across platforms which makes testing
// for consistency in the important parts of this program difficult.
float myrand(void)
{
    static int seed = 72191;
    char sq[22];

    if (ref_file)
    {
        seed *= seed;
        sprintf(sq, "%010d", seed);
        // pull the middle 5 digits out of sq
        sq[8] = 0;
        seed = atoi(&sq[3]);

        return seed/99999.f;
    }
    else
    {
        return rand()/(float)RAND_MAX;
    }
}

void initParticles(cData *p, int dx, int dy)
{
    int i, j;

    for (i = 0; i < dy; i++)
    {
        for (j = 0; j < dx; j++)
        {
            p[i*dx+j].x = (j+0.5f+(myrand() - 0.5f))/dx;
            p[i*dx+j].y = (i+0.5f+(myrand() - 0.5f))/dy;
        }
    }
}

void keyboard(unsigned char key, int x, int y)
{
    // Only enable key functionallity if 
    // the automatic stress test is not 
    // taking place
    if(!stress_test) switch (key)
    {
        case 27:
            g_bExitESC = true;
            #if defined (__APPLE__) || defined(MACOSX)
                exit(EXIT_SUCCESS);
            #else
                glutDestroyWindow(glutGetWindow());
                return;
            #endif
            break;

        case 'r':
            memset(hvfield, 0, sizeof(cData) * DS);
            initParticles(particles, DIM, DIM);

            glBindBuffer(GL_ARRAY_BUFFER, vbo);
            glBufferData(GL_ARRAY_BUFFER, sizeof(cData) * DS,
                            particles, GL_DYNAMIC_DRAW_ARB);
            glBindBuffer(GL_ARRAY_BUFFER, 0);

            break;

        default:
            break;
    }
}

void click(int button, int updown, int x, int y)
{
    // Only enable click functionallity if 
    // the automatic stress test is not 
    // taking place
    if(!stress_test) {
        lastx = x;
        lasty = y;
        clicked = !clicked; 
    }
}

void motion(int x, int y)
{
    // Convert motion coordinates to domain
    float fx = (lastx / (float)wWidth);
    float fy = (lasty / (float)wHeight);
    int nx = (int)(fx * DIM);
    int ny = (int)(fy * DIM);

    if (clicked && nx < DIM-FR && nx > FR-1 && ny < DIM-FR && ny > FR-1)
    {
        int ddx = x - lastx;
        int ddy = y - lasty;
        fx = ddx / (float)wWidth;
        fy = ddy / (float)wHeight;
        int spy = ny-FR;
        int spx = nx-FR;
        addForces(hvfield, DIM, DIM, spx, spy, FORCE * DT * fx, FORCE * DT * fy, FR);
        lastx = x;
        lasty = y;
    }

    glutPostRedisplay();
}

void reshape(int x, int y)
{
    wWidth = x;
    wHeight = y;
    glViewport(0, 0, x, y);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0, 1, 1, 0, 0, 1);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glutPostRedisplay();
}

void cleanup(void)
{
    // Free all host and device resources
    free(hvfield);
    free(particles);

    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glDeleteBuffers(1, &vbo);
}

int initGL(int *argc, char **argv)
{
    glutInit(argc, argv);
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
    glutInitWindowSize(wWidth, wHeight);
    glutCreateWindow("Compute Stable Fluids");
    GLenum err = glewInit();
    if (GLEW_OK != err)
    {
        /* Problem: glewInit failed, something is seriously wrong. */
	fprintf(stderr, "Error: %s\n", glewGetErrorString(err));
	return false;
    }
    fprintf(stdout, "Status: Using GLEW %s\n", glewGetString(GLEW_VERSION));
    glutDisplayFunc(display);
    glutKeyboardFunc(keyboard);
    glutMouseFunc(click);
    glutMotionFunc(motion);
    glutReshapeFunc(reshape);

    return true;
}


int main(int argc, char **argv)
{
    char* lastArg = argv[argc-1];

    if (strcmp(lastArg,"test") == 0)
    {
        // Restablecemos la dimensión del problema por defecto 
        // a algo mas manejable que sirva para comprobar el funcionamiento
        // a simple vista.
        DIM = 8;
	FR = 2; // Dismunimos radio de aplicacion de fuerza 
                // a un tamaño razonable para nuestra pequeña dimensión
        updateVariables();
    }
    else if (strcmp(lastArg,"stress_test") == 0)
    {
        printf("Starting automatic stress test...\n");
        stress_test = true;
    }

#if defined(__linux__)
    setenv ("DISPLAY", ":0", 0);
#endif

    printf("%s Starting...\n\n", sSDKname);

    // First initialize OpenGL context, so we can properly set the GL for CUDA.
    // This is necessary in order to achieve optimal performance with OpenGL/CUDA interop.
    if (false == initGL(&argc, argv))
    {
        exit(EXIT_SUCCESS);
    }

    // Allocate and initialize host data
    GLint bsize;

    hvfield = (cData *)malloc(sizeof(cData) * DS);
    memset(hvfield, 0, sizeof(cData) * DS);

    vxfield = (cData*) malloc(sizeof(cData) * PDS);
    vyfield = (cData*) malloc(sizeof(cData) * PDS);

    // Create particle array
    particles = (cData *)malloc(sizeof(cData) * DS);
    memset(particles, 0, sizeof(cData) * DS);

    initParticles(particles, DIM, DIM);

    glGenBuffers(1, &vbo);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, sizeof(cData) * DS,
                    particles, GL_DYNAMIC_DRAW_ARB);
    // DRAW: Buffer de escritura/lectura

    glGetBufferParameteriv(GL_ARRAY_BUFFER, GL_BUFFER_SIZE, &bsize);

    if (bsize != ((int) (sizeof(cData) * DS)))
        goto EXTERR;

    glBindBuffer(GL_ARRAY_BUFFER, 0);

    if (strcmp(lastArg,"test") == 0)
    {
        ref_file = (char*) 1;
	printf("Starting dump test...\n");
	dumpTest();
        cleanup();
    }
    else
    {
#if defined (__APPLE__) || defined(MACOSX)
        atexit(cleanup);
#else
        glutCloseFunc(cleanup);
#endif
        glutMainLoop();
    }

    if (!ref_file)
    {
        exit(EXIT_SUCCESS);
    }

    return 0;

EXTERR:
    printf("Failed to initialize GL extensions.\n");

    exit(EXIT_FAILURE);
}
