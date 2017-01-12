/** File:		SphData.h
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

/// \file SphData.h
/// \brief holds data for camera rotation and world positioning
/// \author Dongli Zhang modified by Edina Morris

#ifndef SPHDATA_H_
#define SPHDATA_H_

#include "SphHeader.h"
#include "SphType.h"
#include "Vec3.h"

float windowWidth=1000;
float windowHeight=1000;

//initial camera values
float xRot = 0.0f;
float yRot = -30.0f;
float xTrans = 0.0;
float yTrans = 3.0;
float zTrans = -25.0;

int ox;
int oy;
int buttonState;
float xRotLength = 0.0f;
float yRotLength = 0.0f;

//used for scaling and positions
vec3 realWorldOrigin;
vec3 realWorldSide;
vec3 simRatio;

float worldWidth;
float worldHeight;
float worldLength;

#endif //SPHDATA_H_
