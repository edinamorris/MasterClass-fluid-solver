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

/// \file SphMain.cpp
/// \brief Drawing code and initial OpenGL and camera setup
/// \author Dongli Zhang modified by Edina Morris

#include <GL/glew.h>
#include <string.h>

#ifdef __APPLE__
#include <GLUT/glut.h>
#include <OpenGL.h>
#else
#include <GL/glut.h>
#endif

#include "SphHeader.h"
#include "SphData.h"
#include "SphSystem.h"

SPHSystem *sph;

GLuint v, f, p;

void setShaders()
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

//drawing boundary using info from sph_data.h
void drawBox(float _ox, float _oy, float _oz, float _width, float _height, float _length)
{
    glLineWidth(1.0f);
    glColor3f(0.0f, 0.0f, 0.0f);

    glBegin(GL_LINES);   
		
        glVertex3f(_ox, _oy, _oz);
        glVertex3f(_ox+_width, _oy, _oz);

        glVertex3f(_ox, _oy, _oz);
        glVertex3f(_ox, _oy+_height, _oz);

        glVertex3f(_ox, _oy, _oz);
        glVertex3f(_ox, _oy, _oz+_length);

        glVertex3f(_ox+_width, _oy, _oz);
        glVertex3f(_ox+_width, _oy+_height, _oz);

        glVertex3f(_ox+_width, _oy+_height, _oz);
        glVertex3f(_ox, _oy+_height, _oz);

        glVertex3f(_ox, _oy+_height, _oz+_length);
        glVertex3f(_ox, _oy, _oz+_length);

        glVertex3f(_ox, _oy+_height, _oz+_length);
        glVertex3f(_ox, _oy+_height, _oz);

        glVertex3f(_ox+_width, _oy, _oz);
        glVertex3f(_ox+_width, _oy, _oz+_length);

        glVertex3f(_ox, _oy, _oz+_length);
        glVertex3f(_ox+_width, _oy, _oz+_length);

        glVertex3f(_ox+_width, _oy+_height, _oz);
        glVertex3f(_ox+_width, _oy+_height, _oz+_length);

        glVertex3f(_ox+_width, _oy+_height, _oz+_length);
        glVertex3f(_ox+_width, _oy, _oz+_length);

        glVertex3f(_ox, _oy+_height, _oz+_length);
        glVertex3f(_ox+_width, _oy+_height, _oz+_length);

    glEnd();
}

void initSphSystem()
{
    realWorldOrigin.m_x=-10.0f;
    realWorldOrigin.m_y=-10.0f;
    realWorldOrigin.m_z=-10.0f;

    realWorldSide.m_x=8.0f;
    realWorldSide.m_y=8.0f;
    realWorldSide.m_z=8.0f;

	sph=new SPHSystem();
    sph->loadScenario(2);
}

//camera setup
void init()
{
	glewInit();

    glViewport(0, 0, windowWidth, windowHeight);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

    gluPerspective(45.0, (float)windowWidth/windowHeight, 0.1f, 100.0);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

void initRatio()
{
    simRatio.m_x=realWorldSide.m_x/sph->getWorldSize().m_x;
    simRatio.m_y=realWorldSide.m_y/sph->getWorldSize().m_y;
    simRatio.m_z=realWorldSide.m_z/sph->getWorldSize().m_z;
}

//using position and colour calculated in sph_system.cpp
void renderParticles()
{
    glPointSize(1.0f);

    for(uint i=0; i<sph->getNumPart(); i++)
    {
        glBegin(GL_POINTS);
        glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glColor3f(sph->mem[i].m_colour.m_x, sph->mem[i].m_colour.m_y, sph->mem[i].m_colour.m_z);
        glVertex3f(sph->mem[i].m_pos.m_x*simRatio.m_x+realWorldOrigin.m_x,
                   sph->mem[i].m_pos.m_y*simRatio.m_y+realWorldOrigin.m_y,
                   sph->mem[i].m_pos.m_z*simRatio.m_z+realWorldOrigin.m_z);
        glEnd();
    }

}

void displayFunc()
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
    renderParticles();

	glUseProgram(0);
    drawBox(realWorldOrigin.m_x, realWorldOrigin.m_y, realWorldOrigin.m_z, realWorldSide.m_x, realWorldSide.m_y, realWorldSide.m_z);

	glPopMatrix();

    glutSwapBuffers();
	
    glutSetWindowTitle("SPH Multiple Fluid Prototype");
}

void idleFunc()
{
	glutPostRedisplay();
}

void reshapeFunc(GLint _width, GLint _height)
{
    windowWidth=_width;
    windowHeight=_height;

    glViewport(0, 0, _width, _height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

    gluPerspective(45.0, (float)_width/_height, 0.1, 100.0);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glTranslatef(0.0f, 0.0f, -3.0f);
}

//scenario keyboard setup, also start/stop button
void keyboardFunc(unsigned char key, int x, int y)
{
	if(key == ' ')
	{
        sph->setSysRunning(1-sph->getSysRunning());
	}
    //waves
    if(key == '1')
    {
        sph->loadScenario(1);
    }
    //dam
    if(key == '2' )
    {
        sph->loadScenario(2);
    }
    //drop
    if(key == '3')
    {
        sph->loadScenario(3);
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

//needed for arrow keys, movement of camera
void SpecialInput(int key, int x, int y)
{
    switch(key)
    {
        case GLUT_KEY_UP: { yTrans -= 0.5f; break; }
        case GLUT_KEY_DOWN: { yTrans += 0.5f; break; }
        case GLUT_KEY_LEFT: { xTrans += 0.5f; break; }
        case GLUT_KEY_RIGHT: { xTrans -= 0.5f; break; }
    }
    glutPostRedisplay();
}

//mouse state readin
void mouseFunc(int button, int state, int x, int y)
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

//mouse dragging for camera rotation
void motionFunc(int _x, int _y)
{
    float dx, dy;
    dx = (float)(_x - ox);
    dy = (float)(_y - oy);

	if (buttonState == 1) 
	{
		xRotLength += dy / 5.0f;
		yRotLength += dx / 5.0f;
	}

    ox = _x; oy = _y;

	glutPostRedisplay();
}

int main(int argc, char **argv)
{
    glutInit(&argc, argv);
                                        //need this extra depth tag on mac osx, otherwise get translucent behaviour
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
    glutInitWindowPosition(0, 0);
    glutInitWindowSize(windowWidth, windowHeight);
    glutCreateWindow("SPH Fluid 3D");

    initSphSystem();
	init();
    initRatio();
    setShaders();
	glEnable(GL_VERTEX_PROGRAM_POINT_SIZE_NV);
	glEnable(GL_POINT_SPRITE_ARB);
	glTexEnvi(GL_POINT_SPRITE_ARB, GL_COORD_REPLACE_ARB, GL_TRUE);

    glutDisplayFunc(displayFunc);
    glutReshapeFunc(reshapeFunc);
    glutKeyboardFunc(keyboardFunc);
    glutSpecialFunc(SpecialInput);
    glutMouseFunc(mouseFunc);
    glutMotionFunc(motionFunc);
    glutIdleFunc(idleFunc);

    glutMainLoop();

    return 0;
}
