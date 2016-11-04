/** File:		sph_main.cpp
 ** Author:		Dongli Zhang
 ** Contact:	dongli.zhang0129@gmail.com
 **
 ** Copyright (C) Dongli Zhang 2013
 **
 ** This program is free software;  you can redistribute it and/or modify
 ** it under the terms of the GNU General Public License as published by
 ** the Free Software Foundation; either version 2 of the License, or
 ** (at your option) any later version.
 **
 ** This program is distributed in the hope that it will be useful,
 ** but WITHOUT ANY WARRANTY;  without even the implied warranty of
 ** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See
 ** the GNU General Public License for more details.
 **
 ** You should have received a copy of the GNU General Public License
 ** along with this program;  if not, write to the Free Software 
 ** Foundation, 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 */

#include "sph_header.h"
#include "sph_data.h"
#include "sph_system.h"
#include <GL/glew.h>
#include <string.h>

#ifdef __APPLE__
#include <GLUT/glut.h>
#include <OpenGL.h>
#else
#include <GL/glut.h>
#endif

SPHSystem *sph;

GLuint v;
GLuint f;
GLuint p;

void set_shaders()
{
	char *vs=NULL;
	char *fs=NULL;

	vs=(char *)malloc(sizeof(char)*10000);
	fs=(char *)malloc(sizeof(char)*10000);
	memset(vs, 0, sizeof(char)*10000);
	memset(fs, 0, sizeof(char)*10000);

	FILE *fp;
	char c;
	int count;

    fp=fopen("Shader/shader.vs", "r");
	count=0;
	while((c=fgetc(fp)) != EOF)
	{
		vs[count]=c;
		count++;
	}
	fclose(fp);

    fp=fopen("Shader/shader.fs", "r");
	count=0;
	while((c=fgetc(fp)) != EOF)
	{
		fs[count]=c;
		count++;
	}
	fclose(fp);

	v=glCreateShader(GL_VERTEX_SHADER);
	f=glCreateShader(GL_FRAGMENT_SHADER);

	const char *vv;
	const char *ff;
	vv=vs;
	ff=fs;

	glShaderSource(v, 1, &vv, NULL);
	glShaderSource(f, 1, &ff, NULL);

	int success;

	glCompileShader(v);
	glGetShaderiv(v, GL_COMPILE_STATUS, &success);
	if(!success)
	{
		char info_log[5000];
		glGetShaderInfoLog(v, 5000, NULL, info_log);
		printf("Error in vertex shader compilation!\n");
		printf("Info Log: %s\n", info_log);
	}

	glCompileShader(f);
	glGetShaderiv(f, GL_COMPILE_STATUS, &success);
	if(!success)
	{
		char info_log[5000];
		glGetShaderInfoLog(f, 5000, NULL, info_log);
		printf("Error in fragment shader compilation!\n");
		printf("Info Log: %s\n", info_log);
	}

	p=glCreateProgram();
	glAttachShader(p, v);
	glAttachShader(p, f);
	glLinkProgram(p);
	glUseProgram(p);

	free(vs);
	free(fs);
}

void draw_box(float ox, float oy, float oz, float width, float height, float length)
{
    glLineWidth(1.0f);
    glColor3f(0.0f, 0.0f, 0.0f);

    glBegin(GL_LINES);   
		
        glVertex3f(ox, oy, oz);
        glVertex3f(ox+width, oy, oz);

        glVertex3f(ox, oy, oz);
        glVertex3f(ox, oy+height, oz);

        glVertex3f(ox, oy, oz);
        glVertex3f(ox, oy, oz+length);

        glVertex3f(ox+width, oy, oz);
        glVertex3f(ox+width, oy+height, oz);

        glVertex3f(ox+width, oy+height, oz);
        glVertex3f(ox, oy+height, oz);

        glVertex3f(ox, oy+height, oz+length);
        glVertex3f(ox, oy, oz+length);

        glVertex3f(ox, oy+height, oz+length);
        glVertex3f(ox, oy+height, oz);

        glVertex3f(ox+width, oy, oz);
        glVertex3f(ox+width, oy, oz+length);

        glVertex3f(ox, oy, oz+length);
        glVertex3f(ox+width, oy, oz+length);

        glVertex3f(ox+width, oy+height, oz);
        glVertex3f(ox+width, oy+height, oz+length);

        glVertex3f(ox+width, oy+height, oz+length);
        glVertex3f(ox+width, oy, oz+length);

        glVertex3f(ox, oy+height, oz+length);
        glVertex3f(ox+width, oy+height, oz+length);

    glEnd();
}

void init_sph_system()
{
	real_world_origin.x=-10.0f;
	real_world_origin.y=-10.0f;
	real_world_origin.z=-10.0f;

    real_world_side.x=8.0f;
    real_world_side.y=8.0f;
    real_world_side.z=8.0f;

	sph=new SPHSystem();
}

void init()
{
	glewInit();

	glViewport(0, 0, window_width, window_height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

    gluPerspective(45.0, (float)window_width/window_height, 0.1f, 100.0);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

void init_ratio()
{
    sim_ratio.x=real_world_side.x/sph->getWorldSize().x;
    sim_ratio.y=real_world_side.y/sph->getWorldSize().y;
    sim_ratio.z=real_world_side.z/sph->getWorldSize().z;
}

void render_particles()
{
    glPointSize(1.0f);

    for(uint i=0; i<sph->getNumPart(); i++)
    {
        glBegin(GL_POINTS);
        glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glColor3f(sph->mem[i].colour.x, sph->mem[i].colour.y, sph->mem[i].colour.z);
			glVertex3f(sph->mem[i].pos.x*sim_ratio.x+real_world_origin.x, 
						sph->mem[i].pos.y*sim_ratio.y+real_world_origin.y,
                        sph->mem[i].pos.z*sim_ratio.z+real_world_origin.z);
        glEnd();
    }

}

void display_func()
{
    glEnable(GL_POINT_SMOOTH);

    glDepthMask(GL_TRUE);
    glEnable(GL_DEPTH_TEST);

    glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glPushMatrix();

	if(buttonState == 1)
	{
		xRot+=(xRotLength-xRot)*0.1f;
		yRot+=(yRotLength-yRot)*0.1f;
	}

	glTranslatef(xTrans, yTrans, zTrans);
    glRotatef(xRot, 1.0f, 0.0f, 0.0f);
    glRotatef(yRot, 0.0f, 1.0f, 0.0f);

	sph->animation();

    glUseProgram(p);
    render_particles();

	glUseProgram(0);
	draw_box(real_world_origin.x, real_world_origin.y, real_world_origin.z, real_world_side.x, real_world_side.y, real_world_side.z);

	glPopMatrix();

    glutSwapBuffers();
	
    glutSetWindowTitle("SPH Multiple Fluid Prototype");
}

void idle_func()
{
	glutPostRedisplay();
}

void reshape_func(GLint width, GLint height)
{
	window_width=width;
	window_height=height;

	glViewport(0, 0, width, height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

    gluPerspective(45.0, (float)width/height, 0.1, 100.0);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glTranslatef(0.0f, 0.0f, -3.0f);
}

void keyboard_func(unsigned char key, int x, int y)
{
	if(key == ' ')
	{
        sph->setSysRunning(1-sph->getSysRunning());
	}
    //waves
    if(key == '1')
    {
        sph->loadScenario(1);
        //sph->init_system();
    }
    //dam
    if(key == '2' )
    {
        sph->loadScenario(2);
        //sph->damnScenario();
    }
    //drop
    if(key == '3')
    {
        sph->loadScenario(3);
        //sph->dropScenario();
    }
	if(key == 'w')
	{
        zTrans += 1.0f;
	}

	if(key == 's')
	{
        zTrans -= 1.0f;
	}

	if(key == 'a')
	{
        xTrans += 0.5f;
	}

	if(key == 'd')
	{
        xTrans -= 0.5f;
	}

    if(key == 27)
    {
       exit(0);
    }

	glutPostRedisplay();
}

void SpecialInput(int key, int x, int y)
{
    switch(key)
    {
        case GLUT_KEY_UP:
            yTrans -= 0.5f;
        break;
        case GLUT_KEY_DOWN:
            yTrans += 0.5f;
        break;
        case GLUT_KEY_LEFT:
            xTrans += 0.5f;
        break;
        case GLUT_KEY_RIGHT:
            xTrans -= 0.5f;
        break;
    }

    glutPostRedisplay();
}

void mouse_func(int button, int state, int x, int y)
{
	if (state == GLUT_DOWN)
	{
        buttonState = 1;
	}
    else if (state == GLUT_UP)
	{
        buttonState = 0;
	}

    ox = x; oy = y;

    glutPostRedisplay();
}

void motion_func(int x, int y)
{
    float dx, dy;
    dx = (float)(x - ox);
    dy = (float)(y - oy);

	if (buttonState == 1) 
	{
		xRotLength += dy / 5.0f;
		yRotLength += dx / 5.0f;
	}

	ox = x; oy = y;

	glutPostRedisplay();
}

int main(int argc, char **argv)
{
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
    glutInitWindowPosition(0, 0);
    glutInitWindowSize(window_width, window_height);
    glutCreateWindow("SPH Fluid 3D");

    init_sph_system();
	init();
	init_ratio();
	set_shaders();
	glEnable(GL_VERTEX_PROGRAM_POINT_SIZE_NV);
	glEnable(GL_POINT_SPRITE_ARB);
	glTexEnvi(GL_POINT_SPRITE_ARB, GL_COORD_REPLACE_ARB, GL_TRUE);

    glutDisplayFunc(display_func);
	glutReshapeFunc(reshape_func);
	glutKeyboardFunc(keyboard_func);
    glutSpecialFunc(SpecialInput);
	glutMouseFunc(mouse_func);
	glutMotionFunc(motion_func);
	glutIdleFunc(idle_func);

    glutMainLoop();

    return 0;
}
