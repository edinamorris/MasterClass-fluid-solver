/** File:		SphSystem.cpp
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

/// \file sph_system.cpp
/// \brief main calculations for fluid simulation, multiple fluid drift velocity and volume fraction
/// \author Dongli Zhang modified by Edina Morris

#ifdef __APPLE__
#include <GLUT/glut.h>
#include <OpenGL.h>
#else
#include <GL/glut.h>
#endif

#include "SphSystem.h"
#include "SphHeader.h"
#include "Vec3.h"
#include "SphPhaseData.h"

phaseData updatePhases;

SPHSystem::SPHSystem()
{
    //initial setup of world
    m_maxParticle=30000;
    m_numParticle=0;
    m_kernel=0.04f;
    m_worldSize.m_x=0.64f;
    m_worldSize.m_y=0.64f;
    m_worldSize.m_z=0.64f;
    m_cellSize=m_kernel;
    m_gridSize.x=(uint)ceil(m_worldSize.m_x/m_cellSize);
    m_gridSize.y=(uint)ceil(m_worldSize.m_y/m_cellSize);
    m_gridSize.z=(uint)ceil(m_worldSize.m_z/m_cellSize);
    m_totCell=m_gridSize.x*m_gridSize.y*m_gridSize.z;

    //world constants
    m_gravity=vec3(0.0f,-9.8f,0.0f);
    m_wallDamping=-0.5f;
    m_gasConstant=1.0f;
    m_numberOfPhases=PHASES;
    m_timeStep=0.003f;
    m_surfNorm=6.0f;
    m_surfCoe=0.1f;
    m_poly6Value=315.0f/(64.0f * PI * pow(m_kernel, 9));;
    m_spikyValue=-45.0f/(PI * pow(m_kernel, 6));
    m_viscoValue=45.0f/(PI * pow(m_kernel, 6));
    m_gradPoly6=-945/(32 * PI * pow(m_kernel, 9));
    m_lplcPoly6=-945/(8 * PI * pow(m_kernel, 9));
    m_gradSpiky=-45.0/(PI*pow(m_kernel, 6));
    m_kernel2=m_kernel*m_kernel;

    cell=(Particle **)malloc(sizeof(Particle *)*m_totCell);
    m_sysRunning=0;

    //set to 0 for no diffusion effect
    m_diffuse=1;

	printf("Initialize SPH:\n");
    printf("World Width : %f\n", m_worldSize.m_x);
    printf("World Height: %f\n", m_worldSize.m_y);
    printf("World Length: %f\n", m_worldSize.m_z);
    printf("Cell Size  : %f\n", m_cellSize);
    printf("Grid Width : %u\n", m_gridSize.x);
    printf("Grid Height: %u\n", m_gridSize.y);
    printf("Grid Length: %u\n", m_gridSize.z);
    printf("Total Cell : %u\n", m_totCell);
}

SPHSystem::~SPHSystem()
{
	free(mem);
	free(cell);
}

//if system isn't paused then start calculations
void SPHSystem::animation()
{
    if(m_sysRunning == 0)
	{
        return;
    }

    buildTable();
    densityPressure();
    //Multiple fluid calculations
    //step1
    driftVelocity();
    //step2
    advectVolumeFractions();
    //step3
    correctVolumeFraction();
    //step4
    acceleration();
    //step 5 - advect acceleration and velocity
	advection();
    updateColour();
}

void SPHSystem::loadScenario(const int _scenario)
{
    //remove any particles
    if(m_numParticle>0)
    {
        free(mem);
    }

    //reset variables
    m_numParticle=0;

    //waves
    if(_scenario==1)
    {
        mem=(Particle *)malloc(sizeof(Particle)*m_maxParticle);
        phases=(Particle **)malloc(sizeof(Particle *)*m_numberOfPhases);

        vec3 pos, vel;

        vel=vec3(0,0,0);

        for(pos.m_x=m_worldSize.m_x*0.0f; pos.m_x<m_worldSize.m_x*0.4f; pos.m_x+=(m_kernel*0.5f))
        {
            for(pos.m_y=m_worldSize.m_y*0.0f; pos.m_y<m_worldSize.m_y*0.7f; pos.m_y+=(m_kernel*0.5f))
            {
                for(pos.m_z=m_worldSize.m_z*0.0f; pos.m_z<m_worldSize.m_z*0.4f; pos.m_z+=(m_kernel*0.5f))
                {
                    //phase, position, velocity and miscibility
                    addParticle(0, pos, vel,1);
                }
            }
        }

        for(pos.m_x=m_worldSize.m_x*0.0f; pos.m_x<m_worldSize.m_x*0.4f; pos.m_x+=(m_kernel*0.5f))
        {
            for(pos.m_y=m_worldSize.m_y*0.0f; pos.m_y<m_worldSize.m_y*0.7f; pos.m_y+=(m_kernel*0.5f))
            {
                for(pos.m_z=m_worldSize.m_z*0.0f+0.4; pos.m_z<m_worldSize.m_z*0.4f+0.4; pos.m_z+=(m_kernel*0.5f))
                {
                    //phase, position, velocity and miscibility
                    addParticle(1, pos, vel, 1);
                }
            }
        }

        std::cout<<"number of particles >"<<m_numParticle<<"\n";
    }
    //dam setup
    else if(_scenario==2)
    {
        mem=(Particle *)malloc(sizeof(Particle)*m_maxParticle);
        phases=(Particle **)malloc(sizeof(Particle *)*m_numberOfPhases);

        vec3 pos, vel;

        vel=vec3(0,0,0);

        for(pos.m_x=m_worldSize.m_x*0.0f; pos.m_x<m_worldSize.m_x*0.2f; pos.m_x+=(m_kernel*0.5f))
        {
            for(pos.m_y=m_worldSize.m_y*0.0f; pos.m_y<m_worldSize.m_y*0.7f; pos.m_y+=(m_kernel*0.5f))
            {
                for(pos.m_z=m_worldSize.m_z*0.0f+0.1; pos.m_z<m_worldSize.m_z*0.7f+0.1; pos.m_z+=(m_kernel*0.5f))
                {
                    //phase, position, velocity and miscibility
                    addParticle(0, pos, vel, 1);
                }
            }
        }

        for(pos.m_x=m_worldSize.m_x*0.0f+0.5; pos.m_x<m_worldSize.m_x*0.2f+0.5; pos.m_x+=(m_kernel*0.5f))
        {
            for(pos.m_y=m_worldSize.m_y*0.0f; pos.m_y<m_worldSize.m_y*0.7f; pos.m_y+=(m_kernel*0.5f))
            {
                for(pos.m_z=m_worldSize.m_z*0.0f+0.1; pos.m_z<m_worldSize.m_z*0.7f+0.1; pos.m_z+=(m_kernel*0.5f))
                {
                    //phase, position, velocity and miscibility
                    addParticle(1, pos, vel, 1);
                }
            }
        }

        std::cout<<"number of particles >"<<m_numParticle<<"\n";
    }
    //drop simulation
    else if(_scenario==3)
    {
        mem=(Particle *)malloc(sizeof(Particle)*m_maxParticle);
        phases=(Particle **)malloc(sizeof(Particle *)*m_numberOfPhases);

        vec3 pos, vel;

        vel=vec3(0,0,0);

        for(pos.m_x=m_worldSize.m_x*0.0f+0.2; pos.m_x<m_worldSize.m_x*0.35f+0.25; pos.m_x+=(m_kernel*0.5f))
        {
            for(pos.m_y=m_worldSize.m_y*0.0f+0.2; pos.m_y<m_worldSize.m_y*0.3f+0.2; pos.m_y+=(m_kernel*0.5f))
            {
                for(pos.m_z=m_worldSize.m_z*0.0f+0.2; pos.m_z<m_worldSize.m_z*0.35f+0.25; pos.m_z+=(m_kernel*0.5f))
                {
                    //phase, position, velocity and miscibility
                    addParticle(0, pos, vel, 1);
                }
            }
        }
        for(pos.m_x=m_worldSize.m_x*0.0f; pos.m_x<m_worldSize.m_x*1.0; pos.m_x+=(m_kernel*0.5f))
        {
            for(pos.m_y=m_worldSize.m_y*0.0f; pos.m_y<m_worldSize.m_y*0.2f; pos.m_y+=(m_kernel*0.5f))
            {
                for(pos.m_z=m_worldSize.m_z*0.0f; pos.m_z<m_worldSize.m_z*1.0; pos.m_z+=(m_kernel*0.5f))
                {
                    //phase, position, velocity and miscibility
                    addParticle(1, pos, vel, 1);
                }
            }
        }

        std::cout<<"number of particles >"<<m_numParticle<<"\n";
    }
}

//using pre-set phase(liquid) data to assign for a particle belonging to that phase
void SPHSystem::addParticle(const int _phase, const vec3 _pos, const vec3 _vel, const int _misc)
{
    Particle *p;
    p=&(mem[m_numParticle]);

    p->m_id=m_numParticle;
    p->m_pos=_pos;
    p->m_vel=_vel;
    p->m_misc=_misc;
    //set values based on preset phase data
    p->m_colour=updatePhases.phases[_phase].m_colour;
    p->m_mass=updatePhases.phases[_phase].m_individualMass;
    p->m_visc=updatePhases.phases[_phase].m_individualVisc;
    p->m_selfDens=updatePhases.phases[_phase].m_selfDens;
    p->m_lplcColour=updatePhases.phases[_phase].m_selfLplcColor;
    p->m_restdens=updatePhases.phases[_phase].m_dens;
    p->m_dens=p->m_restdens;
    p->m_phase=_phase;
    //adding to the phase counter
    updatePhases.addParticle(_phase);

    //initial values
    p->m_acc=vec3(0,0,0);
    p->m_ev=vec3(0,0,0);
    p->m_pres=0.0f;
    for(int i=0; i<m_numberOfPhases; i++)
    {
        //drift velocity of particle for first and second liquids
        p->m_driftVelocity[i]=vec3(0,0,0);
        //clear any previous volume fraction values first
        p->m_volumeFraction[i]=0.0f;
    }
    //setting original volume fraction, only sets the volume fraction of the particles initial phase to 1
    p->m_volumeFraction[_phase]=1.0f;

    p->next=NULL;
    //still add one to overall count
    m_numParticle++;
}

//spatial partitioning
void SPHSystem::buildTable()
{
    Particle *p;
    int hash;

    for(int i=0; i<int(m_totCell); i++)
	{
		cell[i]=NULL;
	}

    for(int j=0; j<m_numParticle; j++)
    {
        p=&(mem[j]);
        hash=calcCellHash(calcCellPos(p->m_pos));

        if(cell[hash] == NULL)
        {
            p->next=NULL;
            cell[hash]=p;
        }
        else
        {
            p->next=cell[hash];
            cell[hash]=p;
        }
    }
}

//uses volume fraction of a particle for each liquid to calculate
void SPHSystem::updateColour()
{
    Particle *p;
    for(int i=0; i<m_numParticle; i++)
    {
        p=&(mem[i]);
        vec3 updatedColour=vec3();
        for(int phase=0; phase<m_numberOfPhases; phase++)
        {
            updatedColour+=updatePhases.getColour(phase)*p->m_volumeFraction[phase];

            p->m_colour=updatedColour;
        }
    }
}

//calculating density and pressure and storing a particles updated neighbour list
void SPHSystem::densityPressure()
{
	Particle *p;
	Particle *np;

    int3 cellPos;
    int3 nearPos;
	uint hash;

    vec3 relPos;
	float r2;

    for(int i=0; i<m_numParticle; i++)
	{
		p=&(mem[i]); 
        cellPos=calcCellPos(p->m_pos);

        p->m_dens=0.0f;
        p->m_pres=0.0f;
        //each time clearing vector
        p->neighbour.clear();

		for(int x=-1; x<=1; x++)
		{
			for(int y=-1; y<=1; y++)
			{
				for(int z=-1; z<=1; z++)
				{
                    nearPos.x=cellPos.x+x;
                    nearPos.y=cellPos.y+y;
                    nearPos.z=cellPos.z+z;
                    hash=calcCellHash(nearPos);

					if(hash == 0xffffffff)
					{
						continue;
					}

					np=cell[hash];
					while(np != NULL)
					{
                        relPos=np->m_pos-p->m_pos;
                        r2=relPos.m_x*relPos.m_x+relPos.m_y*relPos.m_y+relPos.m_z*relPos.m_z;

                        if(r2<INF || r2>=m_kernel2)
						{
                            np=np->next;
                            continue;
                        }

                        p->m_dens=p->m_dens + p->m_mass * m_poly6Value * pow(m_kernel2-r2, 3);

                        //check for neighbours
                        if(r2>m_kernel*m_kernel)
                        {
                            np=np->next;
                            continue;
                        }
                        //adding neighbours for each particle
                        p->neighbour.push_back(np);

						np=np->next;
					}
				}
			}
		}
        p->m_dens=p->m_dens+p->m_selfDens;
        p->m_pres=(pow(p->m_dens / p->m_restdens, 7) - 1) *m_gasConstant;
	}
}

//Calculate the acceleration of the mixture particle according to Eq. (8). SPH formulations
//of the related terms are provided in Eqs. (19)–(21).
void SPHSystem::acceleration()
{

    //MULTIPLE FLUID PAPER ACCELERATION - NOT WORKING WITH OTHER PARTS
    /*Particle *p;
    for(int i=0; i<m_numParticle; i++)
    {
        p=&(mem[i]);

        //equation 20 - pressure gradient
        vec3 firstPhase=vec3(0,0,0);
        //equation 21 - viscosity tensor
        vec3 secondPhase=vec3(0,0,0);
        //equation 19 - corrective momentum change
        vec3 thirdPhase=vec3(0,0,0);

        //first phase - equation 20
        vec3 pressureGradient=vec3(0,0,0);
        vec3 smoothingKernel=vec3(0,0,0);
        //neighbourhood particles
        for(int j=0; j<int(p->neighbour.size()); j++)
        {
            Particle *neighbour;
            neighbour=(p->neighbour[j]);

            //for all other calculations involving derivative of smoothing function
            smoothingKernel=spikyGrad(p->pos, neighbour->pos);

            pressureGradient+=neighbour->mass*(p->pres+neighbour->pres/2*neighbour->dens)*smoothingKernel;
        }
        firstPhase=pressureGradient/p->restdens;

        //second phase - equation 21
        vec3 viscTensor=vec3(0,0,0);
        //neighbourhood particles
        for(int j=0; j<int(p->neighbour.size()); j++)
        {
            Particle *neighbour;
            neighbour=(p->neighbour[j]);

            //for all other calculations involving derivative of smoothing function
            smoothingKernel=spikyGrad(p->pos, neighbour->pos);

            vec3 positionDiff;
            positionDiff=neighbour->pos-p->pos;

            viscTensor=neighbour->mass/neighbour->dens*(p->visc+neighbour->visc)*(neighbour->vel-p->vel)*
                    ((positionDiff.dot(smoothingKernel)/(positionDiff.dot(smoothingKernel))));
        }
        secondPhase=viscTensor/p->restdens;

        //third phase - equation 19
        vec3 correctMomentum =vec3(0,0,0);
        //neighbourhood particles
        for(int j=0; j<int(p->neighbour.size()); j++)
        {
            Particle *neighbour;
            neighbour=(p->neighbour[j]);

            //for all other calculations involving derivative of smoothing function
            smoothingKernel=spikyGrad(p->pos, neighbour->pos);

            //resetting for new neighbour particle
            vec3 summationPhase=vec3(0,0,0);

            //go through phases
            for(int phase=0; phase<m_numberOfPhases; phase++)
            {
                summationPhase+=neighbour->restdens*(neighbour->volumeFraction[phase]*neighbour->driftVelocity[phase]*(neighbour->driftVelocity[phase].dot(smoothingKernel))+
                        (p->volumeFraction[phase]*p->driftVelocity[phase]*(p->driftVelocity[phase].dot(smoothingKernel))));
            }

            correctMomentum+=neighbour->mass/neighbour->dens*summationPhase;

        }
        thirdPhase=correctMomentum/p->restdens;

        p->acc=(-1*firstPhase+secondPhase+thirdPhase)+m_gravity;
    }*/

    //This acceleration method works
    Particle *p;
	Particle *np;

    int3 cellPos, nearPos;
	uint hash;

    //relative position and velocity
    vec3 relPos, relVel;
    float r2, r, kernelr, V;
    float presKernel, viscKernel, tempForce;

    vec3 gradColor;
    float lplcColor;

    for(int i=0; i<m_numParticle; i++)
	{
		p=&(mem[i]); 
        cellPos=calcCellPos(p->m_pos);

        //resetting variables
        p->m_acc=vec3(0,0,0);
        gradColor=vec3(0,0,0);
        lplcColor=0.0f;
		
		for(int x=-1; x<=1; x++)
		{
			for(int y=-1; y<=1; y++)
			{
				for(int z=-1; z<=1; z++)
				{
                    nearPos.x=cellPos.x+x;
                    nearPos.y=cellPos.y+y;
                    nearPos.z=cellPos.z+z;
                    hash=calcCellHash(nearPos);

					if(hash == 0xffffffff)
					{
						continue;
					}

					np=cell[hash];
					while(np != NULL)
					{
                        relPos=p->m_pos-np->m_pos;
                        r2=relPos.m_x*relPos.m_x+relPos.m_y*relPos.m_y+relPos.m_z*relPos.m_z;

                        if(r2 < m_kernel2 && r2 > INF)
						{
							r=sqrt(r2);
                            V=p->m_mass/np->m_dens/2;
                            kernelr=m_kernel-r;

                            presKernel=m_spikyValue * kernelr * kernelr;
                            tempForce=V * (p->m_pres+np->m_pres) * presKernel;
                            p->m_acc=p->m_acc-relPos*tempForce/r;

                                         //estimated velocity
                            relVel=np->m_ev-p->m_ev;

                            viscKernel=m_viscoValue*(m_kernel-r);
                            tempForce=V * p->m_visc * viscKernel;
                            p->m_acc=p->m_acc + relVel*tempForce;

                            float temp=(-1) * m_gradPoly6 * V * pow(m_kernel2-r2, 2);
                            gradColor += temp * relPos;
                            lplcColor += m_lplcPoly6 * V * (m_kernel2-r2) * (r2-3/4*(m_kernel2-r2));
						}

						np=np->next;
					}
				}
			}
        }

        lplcColor+=p->m_lplcColour/p->m_dens;
        p->m_surfNorm=sqrt(gradColor.m_x*gradColor.m_x+gradColor.m_y*gradColor.m_y+gradColor.m_z*gradColor.m_z);

        if(p->m_surfNorm > m_surfNorm)
		{
            p->m_acc+=(m_surfCoe * lplcColor * gradColor / p->m_surfNorm);
		}
    }
}

//Compute the drift velocity of each phase according to Eq. (9) using pk calculated from Eq. (11) or Eq. (12).
//The SPH formulation of the gradient terms is given in Eqs. (13) and (14).
void SPHSystem::driftVelocity()
{
    Particle *p;
    //for each phase AND each particle, i.e. each particle will have multiple drift velocities and volume fractions attached to it
    for(int i=0; i<m_numParticle; i++)
    {
        p=&(mem[i]);
        for(int phase=0; phase<m_numberOfPhases; phase++)
        {
            p->m_driftVelocity[phase]=vec3(0,0,0);
            float strengthFac=1e-8;
            float diffuseConst;
            if(m_diffuse==1)
            {
                diffuseConst=1e-5f;
            }
            else
            {
                diffuseConst=0.0f;
            }
            vec3 firstPhase=vec3(0,0,0);
            vec3 secondPhase=vec3(0,0,0);
            vec3 thirdPhase=vec3(0,0,0);

            //FIRST PHASE OF EQUATION - accounts for the inertia effect
            float inertia=0;
            //looping through all other phases j represent k' here
            for(int j=0; j<m_numberOfPhases; j++)
            {
                if(j==phase)
                {
                    continue;
                }
                //sum of all other phases related to current k' phase
                inertia+=((p->m_volumeFraction[j]*updatePhases.getDensity(j))/p->m_restdens)*updatePhases.getDensity(j);
            }

            //rest density of current overall particle phase used here
            firstPhase=strengthFac*(updatePhases.getDensity(phase)-inertia)*(m_gravity-p->m_acc);

            //SECOND PHASE OF EQUATION - accounts for the pressure effect
            vec3 spikyKernelGrad;
            //interpolated pressure for current phase [phase]
            vec3 interpolatedPressureInd=vec3(0,0,0);
            //neighbourhood particles
            for(int j=0; j<int(p->neighbour.size()); j++)
            {
                Particle *neighbour;
                neighbour=(p->neighbour[j]);

                //for all other calculations involving derivative of smoothing function
                spikyKernelGrad=spikyGrad(p->m_pos, neighbour->m_pos);

                //calculating first interpolated pressure (for individual particles) - using the spiky for other derivatives of smoothing function
                if((p->m_misc==1) && (neighbour->m_misc==1))
                {
                    interpolatedPressureInd+= (neighbour->m_mass/neighbour->m_dens)*(neighbour->m_pres*neighbour->m_volumeFraction[phase]-p->m_pres*p->m_volumeFraction[phase])*spikyKernelGrad;
                }
                else
                {
                    interpolatedPressureInd = (neighbour->m_mass/neighbour->m_dens)*(neighbour->m_pres-p->m_pres)*spikyKernelGrad;
                }
            }

            //interpolated pressure accumulated for all phases [j]
            vec3 pressure=vec3(0,0,0);
            for(int j=0; j<m_numberOfPhases; j++)
            {
                if(j==phase)
                {
                    continue;
                }
                for(int neigh=0; neigh<int(p->neighbour.size()); neigh++)
                {
                    Particle *neighbour;
                    neighbour=(p->neighbour[neigh]);

                    //for all other calculations involving derivative of smoothing function
                    spikyKernelGrad=spikyGrad(p->m_pos, neighbour->m_pos);

                    //interpolated pressure of all phases
                    vec3 interpolatedPressureAll;

                    if((p->m_misc==1) && (neighbour->m_misc==1))
                    {
                        interpolatedPressureAll=(neighbour->m_mass/neighbour->m_dens)*(neighbour->m_pres*neighbour->m_volumeFraction[j]-p->m_pres*p->m_volumeFraction[j])*spikyKernelGrad;
                    }
                    else
                    {
                        interpolatedPressureAll=(neighbour->m_mass/neighbour->m_dens)*(neighbour->m_pres-p->m_pres)*spikyKernelGrad;
                    }

                    pressure+=((p->m_volumeFraction[j]*updatePhases.getDensity(j))/p->m_restdens)*interpolatedPressureAll;
                }
            }

                                                       // - sum of other phases values
            secondPhase=strengthFac*(interpolatedPressureInd-(pressure));

            //THIRD PHASE OF EQUATION - accounts for the diffusion effect
            //interpolated volume fraction just for this phase [phase]
            vec3 interpolatedVolFracInd=vec3(0,0,0);
            //neighbourhood particles
            for(int j=0; j<int(p->neighbour.size()); j++)
            {
                Particle *neighbour;
                neighbour=(p->neighbour[j]);

                //for all other calculations involving derivative of smoothing function
                spikyKernelGrad=spikyGrad(p->m_pos, neighbour->m_pos);

                //calculating first interpolated volume fraction (for this phase) - using the spiky for other derivatives of smoothing function
                interpolatedVolFracInd+=(neighbour->m_mass/neighbour->m_dens)*(neighbour->m_volumeFraction[phase]-p->m_volumeFraction[phase])*spikyKernelGrad;
            }

            vec3 diffusion=vec3(0,0,0);
            //j is k' here
            for(int j=0; j<m_numberOfPhases; j++)
            {
                if(j==phase)
                {
                    continue;
                }
                for(int neigh=0; neigh<int(p->neighbour.size()); neigh++)
                {
                    Particle *neighbour;
                    neighbour=(p->neighbour[neigh]);

                    //for all other calculations involving derivative of smoothing function
                    spikyKernelGrad=spikyGrad(p->m_pos, neighbour->m_pos);

                    //interpolated volume fraction of all phases [j]                                volume fraction for k' phase
                    vec3 interpolatedVolumeAll=(neighbour->m_mass/neighbour->m_dens)*(neighbour->m_volumeFraction[j]-p->m_volumeFraction[j])*spikyKernelGrad;

                    vec3 volumeMixVolFrac;
                    //nan if dividing by 0
                    if(p->m_volumeFraction[j]==0)
                    {
                        volumeMixVolFrac=vec3(0,0,0);
                    }
                    else
                    {
                        volumeMixVolFrac=interpolatedVolumeAll/p->m_volumeFraction[j];
                    }

                    //volume fraction here uses k' ie. j
                    diffusion+=((p->m_volumeFraction[j]*updatePhases.getDensity(j)/p->m_restdens)*(volumeMixVolFrac));
                }
            }

            vec3 VolFrac=vec3(0.0f,0.0f,0.0f);

            //makes drift velocity nan if volume fraction is 0, cannot divide by 0
            if(p->m_volumeFraction[phase]!=0)
            {
                VolFrac=interpolatedVolFracInd/p->m_volumeFraction[phase];
            }

                                            // - sum of other phases values
            thirdPhase=diffuseConst*((VolFrac)-diffusion);

            //drift velocity - second phase is where pressure relationship will affect the misc/immiscibility of the mixture
            //each particle will have a drift velocity calculated for however many liquids (phases) there are
            p->m_driftVelocity[phase]=firstPhase-secondPhase-thirdPhase;
        }
    }
}

//Advect volume fraction values according to Eq. (7), where the
//relevant SPH formulations are given in Eqs. (17) and (18).
void SPHSystem::advectVolumeFractions()
{
    Particle *p;
    Particle *neighbour;

    for(int i=0; i<m_numParticle; i++)
    {
        p=&(mem[i]);
        for(int phase=0; phase<m_numberOfPhases; phase++)
        {
            float eqs17=0.0f;
            float eqs18=0.0f;
            for(int j=0; j<int(p->neighbour.size()); j++)
            {
                neighbour=(p->neighbour[j]);
                //dont want to compare same particle
                if(p->m_id==neighbour->m_id)
                {
                    continue;
                }

                vec3 spikyGradient=spikyGrad(p->m_pos, neighbour->m_pos);

                //average volume fraction value between particle and it's neighbours
                float diffvolfrac=neighbour->m_volumeFraction[phase]+p->m_volumeFraction[phase]/2;
                //collation of drift velocities by phase fraction per particle
                vec3 volumeFracDrift=(neighbour->m_volumeFraction[phase]*neighbour->m_driftVelocity[phase]+
                                 p->m_volumeFraction[phase]*p->m_driftVelocity[phase]);

                //Eq. (17) reflects the change of volume fraction due to the aggregate motion of the mixture,
                eqs17+=(neighbour->m_mass/neighbour->m_dens)*((diffvolfrac)*(neighbour->m_vel-p->m_vel).dot(spikyGradient));
                //Eq. (18) reflects the change of volume fraction due to the discrepancy between phase velocities,
                //that is, the difference of drift fluxes αk umk between particles.
                eqs18+=((neighbour->m_mass/neighbour->m_dens)*(volumeFracDrift).dot(spikyGradient));
            }
            p->m_volumeFraction[phase]=p->m_volumeFraction[phase]-eqs17-eqs18;
        }
    }
}

//Check the bound of volume fraction according to Eq.(3) and, if the bound is invalidated, correct the volume
//fraction within the particle and calculate the pressure adjustment accordingly as described in Section 5.1.
//For particles with a corrected volume fraction, update into the pressure term the pressure adjustment as Eq. (31).
void SPHSystem::correctVolumeFraction()
{
    Particle *p;

    for(int i=0; i<m_numParticle; i++)
    {
        float volumeFraction=0.0f;
        p=&(mem[i]);
        for(int phase=0; phase<m_numberOfPhases; phase++)
        {
            //correcting any unforseen errors
            bool notANumberCheck=std::isnan(p->m_volumeFraction[phase]);
            if(notANumberCheck==true)
            {
                //resets the particles volume fraction to original starting values
                if(phase==p->m_phase)
                {
                    p->m_volumeFraction[phase]=1.0f;
                }
                else
                {
                    p->m_volumeFraction[phase]=0.0f;
                }
            }
            if(p->m_volumeFraction[phase]<0)
            {
                p->m_volumeFraction[phase]=0.0f;
            }
            //adding all volume fraction values of a particle for every phase
            volumeFraction+=p->m_volumeFraction[phase];
        }
        if(volumeFraction>=0.0f&&volumeFraction<=0.01)
        {
            p->m_volumeFraction[p->m_phase]=1.0f;
            for(int k=0; k<m_numberOfPhases; k++)
            {
                if(p->m_phase!=k)
                {
                    p->m_volumeFraction[k]=0.0f;
                }
            }
        }
        //checking if values add up
        else if(volumeFraction!=1.0)
        {
            float scale=1.0/volumeFraction;
            for(int phase=0; phase<m_numberOfPhases; phase++)
            {
                p->m_volumeFraction[phase]=p->m_volumeFraction[phase]*scale;
            }
        }
    }
}

//used for calculations involving the smoothing kernel function
vec3 SPHSystem::spikyGrad(const vec3 _iPos, const vec3 _jPos)
{
    vec3 gradient;
    vec3 diffPos=_iPos-_jPos;
    float radius=sqrt(diffPos.dot(diffPos));

    gradient=m_gradSpiky*pow(m_kernel-radius,2)*diffPos/radius;

    return gradient;
}

//Final step, advect particles using their velocities and accelerations, also bounds checking
void SPHSystem::advection()
{
	Particle *p;

    for(int i=0; i<m_numParticle; i++)
	{
		p=&(mem[i]);

        p->m_vel=p->m_vel+p->m_acc*m_timeStep/p->m_dens+m_gravity*m_timeStep;
        p->m_pos=p->m_pos+p->m_vel*m_timeStep;

        if(p->m_pos.m_x >= m_worldSize.m_x-BOUNDARY)
		{
            p->m_vel.m_x=p->m_vel.m_x*m_wallDamping;
            p->m_pos.m_x=m_worldSize.m_x-BOUNDARY;
		}

        if(p->m_pos.m_x < 0.0f)
		{
            p->m_vel.m_x=p->m_vel.m_x*m_wallDamping;
            p->m_pos.m_x=0.0f;
		}

        if(p->m_pos.m_y >= m_worldSize.m_y-BOUNDARY)
		{
            p->m_vel.m_y=p->m_vel.m_y*m_wallDamping;
            p->m_pos.m_y=m_worldSize.m_y-BOUNDARY;
		}

        if(p->m_pos.m_y < 0.0f)
		{
            p->m_vel.m_y=p->m_vel.m_y*m_wallDamping;
            p->m_pos.m_y=0.0f;
		}

        if(p->m_pos.m_z >= m_worldSize.m_z-BOUNDARY)
		{
            p->m_vel.m_z=p->m_vel.m_z*m_wallDamping;
            p->m_pos.m_z=m_worldSize.m_z-BOUNDARY;
		}

        if(p->m_pos.m_z < 0.0f)
		{
            p->m_vel.m_z=p->m_vel.m_z*m_wallDamping;
            p->m_pos.m_z=0.0f;
		}

        p->m_ev=(p->m_ev+p->m_vel)/2;
	}
}

int3 SPHSystem::calcCellPos(const vec3 _p)
{
    int3 cellPos;
    cellPos.x = int(floor((_p.m_x) / m_cellSize));
    cellPos.y = int(floor((_p.m_y) / m_cellSize));
    cellPos.z = int(floor((_p.m_z) / m_cellSize));

    return cellPos;
}

uint SPHSystem::calcCellHash(int3 io_cellPos)
{
    if(io_cellPos.x<0 || io_cellPos.x>=(int)m_gridSize.x || io_cellPos.y<0 || io_cellPos.y>=(int)m_gridSize.y || io_cellPos.z<0 || io_cellPos.z>=(int)m_gridSize.z)
	{
		return (uint)0xffffffff;
	}

    io_cellPos.x = io_cellPos.x & (m_gridSize.x-1);
    io_cellPos.y = io_cellPos.y & (m_gridSize.y-1);
    io_cellPos.z = io_cellPos.z & (m_gridSize.z-1);

    return ((uint)(io_cellPos.z))*m_gridSize.y*m_gridSize.x + ((uint)(io_cellPos.y))*m_gridSize.x + (uint)(io_cellPos.x);
}
