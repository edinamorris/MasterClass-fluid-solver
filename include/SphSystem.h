/** File:		SphSystem.h
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

/// \file SphSystem.h
/// \brief main system for calculating fluid simulation including spatial partitioning methods, multiple fluid calculations
///        and advecting particles using acceleration and velocity.
/// \author Dongli Zhang modified by Edina Morris

#ifndef SPHSYSTEM_H_
#define SPHSYSTEM_H_

#include <stdio.h>
#include <stdlib.h>
#include "SphType.h"
#include "vector"
#include "Vec3.h"
#include "SphParticle.h"

//system for fluid simulation and multiple fluid calculations
class SPHSystem
{
public:
	SPHSystem();
	~SPHSystem();

    Particle *mem;
    Particle **cell;
    Particle **phases;
    float m_poly6Value;
    float m_lplcPoly6;
    float m_kernel;
    float m_kernel2;

    //accessors
    uint getNumPart() const {return m_numParticle;}
    vec3 getWorldSize() const {return m_worldSize;}
    uint getSysRunning()const {return m_sysRunning;}
    int getNumPhases() const {return m_numberOfPhases;}

    //mutators
    void setSysRunning(const uint _sysRun){m_sysRunning=_sysRun;}

    //operation functions
    void addParticle(const int _phase, const vec3 _pos, const vec3 _vel, const int _misc);
    vec3 spikyGrad(const vec3 _iPos, const vec3 _jPos);

    //different scenarios and animation set up, by defualt fluid simulation is not started, user must start the animation after choosing scenario
	void animation();
    void loadScenario(const int _scenario);
	void init_system();
    void damnScenario();
    void dropScenario();

//variables
private:
    int m_numParticle;
    int m_numberOfPhases;
    vec3 m_worldSize;
    //initially 0, set to 1 when user hits space bar
    uint m_sysRunning;

    float m_spikyValue;
    float m_viscoValue;

    float m_gradPoly6;
    float m_gradSpiky;
    uint m_maxParticle;

    float m_cellSize;
    uint3 m_gridSize;
    uint m_totCell;

    vec3 m_gravity;
    float m_wallDamping;
    float m_gasConstant;
    float m_timeStep;
    float m_surfNorm;
    float m_surfCoe;
    float m_diffuse;

//functions
private:
    //sorting particles into a grid system each time
    void buildTable();
    int3 calcCellPos(const vec3 _p);
    uint calcCellHash(int3 io_cellPos);
    //multiple fluid calculations
    void densityPressure();
    void acceleration();
    void driftVelocity();
    void advection();
    void advectVolumeFractions();
    void correctVolumeFraction();
    void updateColour();
};

#endif //SPHSYSTEM_H_
