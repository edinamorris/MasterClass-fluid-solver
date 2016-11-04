/** File:		sph_system.cpp
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

#include "sph_system.h"
#include "sph_header.h"
#include "vec3.h"
#include "sph_phasedata.h"

#ifdef __APPLE__
#include <GLUT/glut.h>
#include <OpenGL.h>
#else
#include <GL/glut.h>
#endif

phaseData updatePhases;

SPHSystem::SPHSystem()
{
    max_particle=30000;
	num_particle=0;

    //h
    kernel=0.04f;

    world_size.x=0.64f;
    world_size.y=0.64f;
    world_size.z=0.64f;
	cell_size=kernel;
	grid_size.x=(uint)ceil(world_size.x/cell_size);
	grid_size.y=(uint)ceil(world_size.y/cell_size);
	grid_size.z=(uint)ceil(world_size.z/cell_size);
	tot_cell=grid_size.x*grid_size.y*grid_size.z;

	gravity.x=0.0f; 
    gravity.y=-9.8f;
	gravity.z=0.0f;
	wall_damping=-0.5f;
    gas_constant=1.0f;

    numberOfPhases=2;

	time_step=0.003f;
	surf_norm=6.0f;
	surf_coe=0.1f;

	poly6_value=315.0f/(64.0f * PI * pow(kernel, 9));;
	spiky_value=-45.0f/(PI * pow(kernel, 6));
	visco_value=45.0f/(PI * pow(kernel, 6));

	grad_poly6=-945/(32 * PI * pow(kernel, 9));
	lplc_poly6=-945/(8 * PI * pow(kernel, 9));
    grad_spiky=-45.0/(PI*pow(kernel, 6));

	kernel_2=kernel*kernel;

	cell=(Particle **)malloc(sizeof(Particle *)*tot_cell);

	sys_running=0;

    //miscibility
    misc=1;

	printf("Initialize SPH:\n");
	printf("World Width : %f\n", world_size.x);
	printf("World Height: %f\n", world_size.y);
	printf("World Length: %f\n", world_size.z);
	printf("Cell Size  : %f\n", cell_size);
	printf("Grid Width : %u\n", grid_size.x);
	printf("Grid Height: %u\n", grid_size.y);
	printf("Grid Length: %u\n", grid_size.z);
	printf("Total Cell : %u\n", tot_cell);
	printf("Poly6 Kernel: %f\n", poly6_value);
	printf("Spiky Kernel: %f\n", spiky_value);
	printf("Visco Kernel: %f\n", visco_value);
    //printf("Self Density: %f\n", self_dens);
}

SPHSystem::~SPHSystem()
{
	free(mem);
	free(cell);
}

void SPHSystem::animation()
{
	if(sys_running == 0)
	{
        return;
	}

	build_table();
	comp_dens_pres();
	comp_force_adv();
    //driftVelocity();
	advection();
}

void SPHSystem::loadScenario(int _scenario)
{
    //remove any particles
    if(num_particle>0)
    {
        free(mem);
    }

    //reset variables
    num_particle=0;
    for(int i=0; i<numberOfPhases; i++)
    {
        updatePhases.setNumberOfParticles(0,i);
    }

    //waves
    if(_scenario==1)
    {
        mem=(Particle *)malloc(sizeof(Particle)*max_particle);
        phases=(Particle **)malloc(sizeof(Particle *)*numberOfPhases);

        vec3 pos;
        vec3 vel;

        vel.x=0.0f;
        vel.y=0.0f;
        vel.z=0.0f;

        for(pos.x=world_size.x*0.0f; pos.x<world_size.x*0.4f; pos.x+=(kernel*0.5f))
        {
            for(pos.y=world_size.y*0.0f; pos.y<world_size.y*0.7f; pos.y+=(kernel*0.5f))
            {
                for(pos.z=world_size.z*0.0f; pos.z<world_size.z*0.4f; pos.z+=(kernel*0.5f))
                {
                    //phase, position and velocity
                    add_particle(0, pos, vel);
                }
            }
        }

        for(pos.x=world_size.x*0.0f; pos.x<world_size.x*0.4f; pos.x+=(kernel*0.5f))
        {
            for(pos.y=world_size.y*0.0f; pos.y<world_size.y*0.7f; pos.y+=(kernel*0.5f))
            {
                for(pos.z=world_size.z*0.0f+0.4; pos.z<world_size.z*0.4f+0.4; pos.z+=(kernel*0.5f))
                {
                    //phase, position and velocity
                    add_particle(1, pos, vel);
                }
            }
        }

        for(int i=0; i<numberOfPhases; i++)
        {
            int numParticlePhase = updatePhases.getNumberOfParticles(i);
            int volumeFraction=numParticlePhase/num_particle;
            updatePhases.setVolumeFraction(volumeFraction, i);
        }

        printf("Init Particle: %u\n", num_particle);
    }
    //damn
    else if(_scenario==2)
    {
        mem=(Particle *)malloc(sizeof(Particle)*max_particle);
        phases=(Particle **)malloc(sizeof(Particle *)*numberOfPhases);

        vec3 pos;
        vec3 vel;

        vel.x=0.0f;
        vel.y=0.0f;
        vel.z=0.0f;

        for(pos.x=world_size.x*0.0f; pos.x<world_size.x*0.2f; pos.x+=(kernel*0.5f))
        {
            for(pos.y=world_size.y*0.0f; pos.y<world_size.y*0.7f; pos.y+=(kernel*0.5f))
            {
                for(pos.z=world_size.z*0.0f+0.1; pos.z<world_size.z*0.7f+0.1; pos.z+=(kernel*0.5f))
                {
                    //phase, position and velocity
                    add_particle(0, pos, vel);
                }
            }
        }

        for(pos.x=world_size.x*0.0f+0.5; pos.x<world_size.x*0.2f+0.5; pos.x+=(kernel*0.5f))
        {
            for(pos.y=world_size.y*0.0f; pos.y<world_size.y*0.7f; pos.y+=(kernel*0.5f))
            {
                for(pos.z=world_size.z*0.0f+0.1; pos.z<world_size.z*0.7f+0.1; pos.z+=(kernel*0.5f))
                {
                    //phase, position and velocity
                    add_particle(1, pos, vel);
                }
            }
        }

        for(int i=0; i<numberOfPhases; i++)
        {
            int numParticlePhase = updatePhases.getNumberOfParticles(i);
            int volumeFraction=numParticlePhase/num_particle;
            updatePhases.setVolumeFraction(volumeFraction, i);
        }

        printf("Init Particle: %u\n", num_particle);
    }
    //drop simulation
    else if(_scenario==3)
    {
        mem=(Particle *)malloc(sizeof(Particle)*max_particle);
        phases=(Particle **)malloc(sizeof(Particle *)*numberOfPhases);

        vec3 pos;
        vec3 vel;

        vel.x=0.0f;
        vel.y=0.0f;
        vel.z=0.0f;

        //weird depth test for second set of particles, become completely transparent
        for(pos.x=world_size.x*0.0f+0.2; pos.x<world_size.x*0.4f+0.2; pos.x+=(kernel*0.5f))
        {
            for(pos.y=world_size.y*0.0f+0.3; pos.y<world_size.y*0.5f+0.3; pos.y+=(kernel*0.5f))
            {
                for(pos.z=world_size.z*0.0f+0.2; pos.z<world_size.z*0.4f+0.2; pos.z+=(kernel*0.5f))
                {
                    //phase, position and velocity
                    add_particle(0, pos, vel);
                }
            }
        }
        for(pos.x=world_size.x*0.0f; pos.x<world_size.x*1.0; pos.x+=(kernel*0.5f))
        {
            for(pos.y=world_size.y*0.0f; pos.y<world_size.y*0.25f; pos.y+=(kernel*0.5f))
            {
                for(pos.z=world_size.z*0.0f; pos.z<world_size.z*1.0; pos.z+=(kernel*0.5f))
                {
                    //phase, position and velocity
                    add_particle(1, pos, vel);
                }
            }
        }

        for(int i=0; i<numberOfPhases; i++)
        {
            int numParticlePhase = updatePhases.getNumberOfParticles(i);
            int volumeFraction=numParticlePhase/num_particle;
            updatePhases.setVolumeFraction(volumeFraction, i);
        }

        printf("Init Particle: %u\n", num_particle);
    }
}

//changing sp particles get added into lists
void SPHSystem::add_particle(int _phase, vec3 pos, vec3 vel)
{
    Particle *p;
    p=&(mem[num_particle]);

    //phases[_phase]=NULL;
    p->id=num_particle;

    p->pos=pos;
    p->vel=vel;
    //set values based on preset phase data
    p->colour=updatePhases.phases[_phase].colour;
    p->mass=updatePhases.phases[_phase].individualMass;
    p->visc=updatePhases.phases[_phase].individualVisc;
    p->selfDens=updatePhases.phases[_phase].selfDens;
    p->lplcColour=updatePhases.phases[_phase].self_lplc_color;
    p->restdens=updatePhases.phases[_phase].dens;
    p->dens=p->restdens;
    //adding to the phase counter
    updatePhases.addParticle(_phase);

    p->acc.x=0.0f;
    p->acc.y=0.0f;
    p->acc.z=0.0f;
    p->ev.x=0.0f;
    p->ev.y=0.0f;
    p->ev.z=0.0f;
    p->driftVelocity.x=0.0f;
    p->driftVelocity.y=0.0f;
    p->driftVelocity.z=0.0f;
    p->pres=0.0f;

    p->next=NULL;

    //adding p to the phases list
    //phases[_phase]=p;

    //still add one to overall count
    num_particle++;
}

void SPHSystem::build_table()
{
    Particle *p;
    //Particle **pp;
	uint hash;

	for(uint i=0; i<tot_cell; i++)
	{
		cell[i]=NULL;
	}

        for(uint j=0; j<num_particle; j++)
        {
            p=&(mem[j]);
            //pp[j]->pos;
            hash=calc_cell_hash(calc_cell_pos(p->pos));

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

void SPHSystem::comp_dens_pres()
{
	Particle *p;
	Particle *np;
    Particle **pp;

	int3 cell_pos;
	int3 near_pos;
	uint hash;

    vec3 rel_pos;
	float r2;

	for(uint i=0; i<num_particle; i++)
	{
		p=&(mem[i]); 
		cell_pos=calc_cell_pos(p->pos);

		p->dens=0.0f;
		p->pres=0.0f;
        //each time clearing vector
        p->neighbour.clear();

		for(int x=-1; x<=1; x++)
		{
			for(int y=-1; y<=1; y++)
			{
				for(int z=-1; z<=1; z++)
				{
					near_pos.x=cell_pos.x+x;
					near_pos.y=cell_pos.y+y;
					near_pos.z=cell_pos.z+z;
					hash=calc_cell_hash(near_pos);

					if(hash == 0xffffffff)
					{
						continue;
					}

					np=cell[hash];
					while(np != NULL)
					{
						rel_pos.x=np->pos.x-p->pos.x;
						rel_pos.y=np->pos.y-p->pos.y;
						rel_pos.z=np->pos.z-p->pos.z;
						r2=rel_pos.x*rel_pos.x+rel_pos.y*rel_pos.y+rel_pos.z*rel_pos.z;

						if(r2<INF || r2>=kernel_2)
						{
                            np=np->next;
                            continue;
                        }

                        p->dens=p->dens + p->mass * poly6_value * pow(kernel_2-r2, 3);

                        //check for neighbours
                        if(r2>kernel*kernel)
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
        p->dens=p->dens+p->selfDens;
        p->pres=(pow(p->dens / p->restdens, 7) - 1) *gas_constant;
	}
}

void SPHSystem::comp_force_adv()
{
	Particle *p;
	Particle *np;

	int3 cell_pos;
	int3 near_pos;
	uint hash;

    vec3 rel_pos;
    vec3 rel_vel;

	float r2;
	float r;
	float kernel_r;
	float V;

	float pres_kernel;
	float visc_kernel;
	float temp_force;

    vec3 grad_color;
	float lplc_color;

	for(uint i=0; i<num_particle; i++)
	{
		p=&(mem[i]); 
		cell_pos=calc_cell_pos(p->pos);

		p->acc.x=0.0f;
		p->acc.y=0.0f;
		p->acc.z=0.0f;

		grad_color.x=0.0f;
		grad_color.y=0.0f;
		grad_color.z=0.0f;
		lplc_color=0.0f;
		
		for(int x=-1; x<=1; x++)
		{
			for(int y=-1; y<=1; y++)
			{
				for(int z=-1; z<=1; z++)
				{
					near_pos.x=cell_pos.x+x;
					near_pos.y=cell_pos.y+y;
					near_pos.z=cell_pos.z+z;
					hash=calc_cell_hash(near_pos);

					if(hash == 0xffffffff)
					{
						continue;
					}

					np=cell[hash];
					while(np != NULL)
					{
						rel_pos.x=p->pos.x-np->pos.x;
						rel_pos.y=p->pos.y-np->pos.y;
						rel_pos.z=p->pos.z-np->pos.z;
						r2=rel_pos.x*rel_pos.x+rel_pos.y*rel_pos.y+rel_pos.z*rel_pos.z;

						if(r2 < kernel_2 && r2 > INF)
						{
							r=sqrt(r2);
                            V=p->mass/np->dens/2;
							kernel_r=kernel-r;

							pres_kernel=spiky_value * kernel_r * kernel_r;
                            temp_force=V * (p->pres+np->pres) * pres_kernel;
							p->acc.x=p->acc.x-rel_pos.x*temp_force/r;
							p->acc.y=p->acc.y-rel_pos.y*temp_force/r;
							p->acc.z=p->acc.z-rel_pos.z*temp_force/r;

							rel_vel.x=np->ev.x-p->ev.x;
							rel_vel.y=np->ev.y-p->ev.y;
							rel_vel.z=np->ev.z-p->ev.z;

							visc_kernel=visco_value*(kernel-r);
                            temp_force=V * p->visc * visc_kernel;
							p->acc.x=p->acc.x + rel_vel.x*temp_force; 
							p->acc.y=p->acc.y + rel_vel.y*temp_force; 
							p->acc.z=p->acc.z + rel_vel.z*temp_force; 

							float temp=(-1) * grad_poly6 * V * pow(kernel_2-r2, 2);
							grad_color.x += temp * rel_pos.x;
							grad_color.y += temp * rel_pos.y;
							grad_color.z += temp * rel_pos.z;
							lplc_color += lplc_poly6 * V * (kernel_2-r2) * (r2-3/4*(kernel_2-r2));
						}

						np=np->next;
					}
				}
			}
		}

        lplc_color+=p->lplcColour/p->dens;
		p->surf_norm=sqrt(grad_color.x*grad_color.x+grad_color.y*grad_color.y+grad_color.z*grad_color.z);

		if(p->surf_norm > surf_norm)
		{
			p->acc.x+=surf_coe * lplc_color * grad_color.x / p->surf_norm;
			p->acc.y+=surf_coe * lplc_color * grad_color.y / p->surf_norm;
			p->acc.z+=surf_coe * lplc_color * grad_color.z / p->surf_norm;
		}
	}
}

//needs some work - mainly design for more than two phases
void SPHSystem::driftVelocity()
{
    Particle *p;

    for(uint i=0; i<num_particle; i++)
    {
        p=&(mem[i]);
        //DRIFT VELOCITY - Need to change design so it's not reliant on just two phases
        //user constants
        //float strengthFac=0.000001f;
        float diffuseConst=0.0001f;
        float strengthFac=1.0f;
        vec3 firstPhase;
        vec3 secondPhase;
        vec3 thirdPhase;

         //FIRST PHASE OF EQUATION
        //sigma k, all other phases
        float summation=0;
        float mixDensity=0;
        //looping through other phases (not self)
        for(int i=0; i<numberOfPhases; i++)
        {
            //only when it's not itself
            if(numberOfPhases!=p->phase)
            {
                summation+=(updatePhases.getVolumeFraction(i)*updatePhases.getDensity(i));
                //not sure if this is correct order for summation
                summation=summation*updatePhases.getDensity(i);
            }
            //mixture Density
            mixDensity+=(updatePhases.getVolumeFraction(i)*updatePhases.getDensity(i));

        }

        //changed first part of equation from rest density as that's the phase density, believe it should be individual particle density
        firstPhase.x=strengthFac*(p->dens-(summation/mixDensity))*p->acc.x;
        firstPhase.y=strengthFac*(p->dens-(summation/mixDensity))*p->acc.y;
        firstPhase.z=strengthFac*(p->dens-(summation/mixDensity))*p->acc.y;

        //SECOND PHASE OF EQUATION
        vec3 poly6KernelGrad;
        vec3 spikyKernelGrad;
        vec3 interpolatedPressureInd;
        //reset for each particle
        interpolatedPressureInd.x=0.0f;
        interpolatedPressureInd.y=0.0f;
        interpolatedPressureInd.z=0.0f;
        //neighbourhood particles
        for(int j=0; j<int(p->neighbour.size()); j++)
        {
            Particle *neighbour;
            neighbour=(p->neighbour[j]);

            vec3 interpolatedDensity;
            vec3 calc_1;
            float phasePressureSelf;
            float phasePressureN;
            if(misc==1)
            {
                //mixture pressure * individual particles volume fraction
                phasePressureSelf=p->pres*p->volFrac;
                phasePressureN=neighbour->pres*neighbour->volFrac;
            }
            else
            {
                phasePressureSelf=p->pres;
                phasePressureN=neighbour->pres;
            }

            // /\Wij - for density interpolation
            poly6KernelGrad=poly6Grad(p->pos, neighbour->pos);

            //for all other calculations involving derivative of smoothing function
            spikyKernelGrad=spikyGrad(p->pos, neighbour->pos);

            //following standard method for interpolated density
            interpolatedDensity.x=(neighbour->mass*poly6KernelGrad.x);
            interpolatedDensity.y=(neighbour->mass*poly6KernelGrad.y);
            interpolatedDensity.z=(neighbour->mass*poly6KernelGrad.z);

            //calculating first interpolated pressure (for individual particles) - using the spiky for other derivatives of smoothing function
            calc_1.x = (neighbour->mass/interpolatedDensity.x)*(phasePressureN-phasePressureSelf)*spikyKernelGrad.x;
            calc_1.y = (neighbour->mass/interpolatedDensity.y)*(phasePressureN-phasePressureSelf)*spikyKernelGrad.y;
            calc_1.z = (neighbour->mass/interpolatedDensity.z)*(phasePressureN-phasePressureSelf)*spikyKernelGrad.z;

            //accumulates to form the sum of all this particles neighbours
            interpolatedPressureInd.x+=calc_1.x;
            interpolatedPressureInd.y+=calc_1.y;
            interpolatedPressureInd.z+=calc_1.z;
        }
        //for mixture
        float phasePressureSelf;
        float phasePressurePhase;
        if(misc==1)
        {
            //mixture pressure * phase volume fraction - NEED TO GET PRESSURE OF EACH PHASE, NOT INDIVIDUAL PARTICLE
            phasePressureSelf=p->pres*updatePhases.getVolumeFraction(i);
            phasePressurePhase=p->pres*updatePhases.getVolumeFraction(1);//needs to be automated
        }
        else
        {
            phasePressureSelf=p->pres;
            phasePressurePhase=p->pres;
        }
        //interpolated pressure of other phases (summed of all OTHER phases rather than neighbours) - in this case only 2 phases needs changing
        vec3 interpolatedPressureMix;
        //not sure about this, may not be interpolated as its for the other phases - using poly6 for interpolated density
        vec3 interpolatedDensityMix;
        interpolatedDensityMix.x=updatePhases.getMass(i)*poly6KernelGrad.x;
        interpolatedDensityMix.y=updatePhases.getMass(i)*poly6KernelGrad.y;
        interpolatedDensityMix.z=updatePhases.getMass(i)*poly6KernelGrad.z;

        vec3 massSummation;
        massSummation.x=0;
        massSummation.y=0;
        massSummation.z=0;
        //looping through other phases (not self)
        for(int i=0; i<numberOfPhases; i++)
        {
            //only when it's not itself
            if(numberOfPhases!=p->phase)
            {
                massSummation.x+=(updatePhases.getMass(i)/interpolatedDensityMix.x);
                massSummation.y+=(updatePhases.getMass(i)/interpolatedDensityMix.y);
                massSummation.z+=(updatePhases.getMass(i)/interpolatedDensityMix.z);

            }
        }
                                                                            //pressure of other phase minus pressure of current phase
        interpolatedPressureMix.x=(massSummation.x)*(phasePressurePhase-phasePressureSelf)*spikyKernelGrad.x;
        interpolatedPressureMix.y=(massSummation.y)*(phasePressurePhase-phasePressureSelf)*spikyKernelGrad.y;
        interpolatedPressureMix.z=(massSummation.z)*(phasePressurePhase-phasePressureSelf)*spikyKernelGrad.z;

        float summation2=0;
        float mixDensity2=0;
        //looping through other phases (not self)
        for(int i=0; i<numberOfPhases; i++)
        {
            //only when it's not itself
            if(numberOfPhases!=p->phase)
            {
                summation2+=(updatePhases.getVolumeFraction(i)*updatePhases.getDensity(i));
                //not sure if this is correct order for summation
                summation2=summation2*updatePhases.getDensity(i);
            }
            //mixture Density
            mixDensity2+=(updatePhases.getVolumeFraction(i)*updatePhases.getDensity(i));

        }

                                                              // - sum of other phases values
        secondPhase.x=strengthFac*(interpolatedPressureInd.x-(summation2)/(mixDensity2)*interpolatedPressureMix.x);//may have to multiply these as part of the summation
        secondPhase.y=strengthFac*(interpolatedPressureInd.y-(summation2)/(mixDensity2)*interpolatedPressureMix.y);
        secondPhase.z=strengthFac*(interpolatedPressureInd.z-(summation2)/(mixDensity2)*interpolatedPressureMix.z);

        //THIRD PHASE OF EQUATION
        vec3 interpolatedVolFracInd;
        //reset for each particle
        interpolatedVolFracInd.x=0.0f;
        interpolatedVolFracInd.y=0.0f;
        interpolatedVolFracInd.z=0.0f;
        //neighbourhood particles
        for(int j=0; j<int(p->neighbour.size()); j++)
        {
            Particle *neighbour;
            neighbour=(p->neighbour[j]);

            vec3 interpolatedDensityInd;
            vec3 calc_1;

            // /\Wij - for density interpolation
            poly6KernelGrad=poly6Grad(p->pos, neighbour->pos);

            //for all other calculations involving derivative of smoothing function
            spikyKernelGrad=spikyGrad(p->pos, neighbour->pos);

            //following standard method for interpolated density
            interpolatedDensityInd.x=(neighbour->mass*poly6KernelGrad.x);
            interpolatedDensityInd.y=(neighbour->mass*poly6KernelGrad.y);
            interpolatedDensityInd.z=(neighbour->mass*poly6KernelGrad.z);

            //calculating first interpolated volume fraction (for individual particles) - using the spiky for other derivatives of smoothing function
            calc_1.x = (neighbour->mass/interpolatedDensityInd.x)*(neighbour->volFrac-p->volFrac)*spikyKernelGrad.x;
            calc_1.y = (neighbour->mass/interpolatedDensityInd.y)*(neighbour->volFrac-p->volFrac)*spikyKernelGrad.y;
            calc_1.z = (neighbour->mass/interpolatedDensityInd.z)*(neighbour->volFrac-p->volFrac)*spikyKernelGrad.z;

            //accumulates to form the sum of all this particles neighbours
            interpolatedVolFracInd.x+=calc_1.x;
            interpolatedVolFracInd.y+=calc_1.y;
            interpolatedVolFracInd.z+=calc_1.z;
        }
        //interpolated pressure of other phases (summed of all OTHER phases rather than neighbours) - in this case only 2 phases needs changing
        vec3 interpolatedVolumeMix;
        //not sure about this, may not be interpolated as its for the other phases - using poly6 for interpolated density
        interpolatedDensityMix.x=updatePhases.getMass(i)*poly6KernelGrad.x;
        interpolatedDensityMix.y=updatePhases.getMass(i)*poly6KernelGrad.y;
        interpolatedDensityMix.z=updatePhases.getMass(i)*poly6KernelGrad.z;

        vec3 summationMass;
        summationMass.x-0;
        summationMass.y-0;
        summationMass.z=0;
        float mixVolFrac=0;
        //looping through other phases (not self)
        for(int i=0; i<numberOfPhases; i++)
        {
            //only when it's not itself
            if(numberOfPhases!=p->phase)
            {
                summationMass.x+=(updatePhases.getMass(i)/interpolatedDensityMix.x);
                summationMass.y+=(updatePhases.getMass(i)/interpolatedDensityMix.y);
                summationMass.z+=(updatePhases.getMass(i)/interpolatedDensityMix.z);
                mixVolFrac+=updatePhases.getVolumeFraction(i);
            }

        }
                                                  //other phase volume fraction - self's volume fraction
        interpolatedVolumeMix.x=(summationMass.x)*(mixVolFrac-updatePhases.getVolumeFraction(i))*spikyKernelGrad.x;
        interpolatedVolumeMix.y=(summationMass.y)*(mixVolFrac-updatePhases.getVolumeFraction(i))*spikyKernelGrad.y;
        interpolatedVolumeMix.z=(summationMass.z)*(mixVolFrac-updatePhases.getVolumeFraction(i))*spikyKernelGrad.z;

        float summation3=0;
        float mixDensity3=0;
        float volumeMix=0;
        //looping through other phases (not self)
        for(int i=0; i<numberOfPhases; i++)
        {
            //only when it's not itself
            if(numberOfPhases!=p->phase)
            {
                summation3+=(updatePhases.getVolumeFraction(i)*updatePhases.getDensity(i));
                volumeMix+=(updatePhases.getVolumeFraction(i));
            }
            //mixture Density
            mixDensity3+=(updatePhases.getVolumeFraction(i)*updatePhases.getDensity(i));

        }

                                                                       // - sum of other phases values
        thirdPhase.x=diffuseConst*((interpolatedVolFracInd.x/p->volFrac)-(summation3)/(mixDensity3)
                                  *(interpolatedVolumeMix.x/volumeMix));
        thirdPhase.y=diffuseConst*((interpolatedVolFracInd.y/p->volFrac)-(summation3)/(mixDensity3)
                                  *(interpolatedVolumeMix.y/volumeMix));
        thirdPhase.z=diffuseConst*((interpolatedVolFracInd.z/p->volFrac)-(summation3)/(mixDensity3)
                                  *(interpolatedVolumeMix.z/volumeMix));


        //drift velocity - second phase is where pressure relationship will affect the misc/immiscibility of the mixture
        p->driftVelocity.x=firstPhase.x-secondPhase.x-thirdPhase.x;
        p->driftVelocity.y=firstPhase.y-secondPhase.y-thirdPhase.y;
        p->driftVelocity.z=firstPhase.z-secondPhase.z-thirdPhase.z;
    }
}

//used for interpolated densities
vec3 SPHSystem::poly6Grad(vec3 _iPos, vec3 _jPos)
{
    vec3 gradient;
    vec3 diffPos=_iPos-_jPos;
    float radiusSquared=diffPos.dot(diffPos);

    gradient=grad_poly6*pow((kernel*kernel)-radiusSquared, 2)*diffPos;
    return gradient;
}

//used for all other calculations involving the smoothing kernel function
vec3 SPHSystem::spikyGrad(vec3 _iPos, vec3 _jPos)
{
    vec3 gradient;
    vec3 diffPos=_iPos-_jPos;
    float radius=sqrt(diffPos.dot(diffPos));

    gradient=grad_spiky*pow(kernel-radius,2)*diffPos/radius;
    return gradient;
}

void SPHSystem::advection()
{
	Particle *p;

	for(uint i=0; i<num_particle; i++)
	{
		p=&(mem[i]);

        //std::cout<<"driftVelocity X -> "<<p->driftVelocity.x<<"\n";
        //std::cout<<"driftVelocity Y -> "<<p->driftVelocity.y<<"\n";
        //std::cout<<"driftVelocity Z -> "<<p->driftVelocity.z<<"\n";

        //potentially add drift velocity in with this step or have it replace velocity as a whole
        p->vel.x=p->vel.x+p->driftVelocity.x+p->acc.x*time_step/p->dens+gravity.x*time_step;
        p->vel.y=p->vel.y+p->driftVelocity.y+p->acc.y*time_step/p->dens+gravity.y*time_step;
        p->vel.z=p->vel.z+p->driftVelocity.z+p->acc.z*time_step/p->dens+gravity.z*time_step;

		p->pos.x=p->pos.x+p->vel.x*time_step;
		p->pos.y=p->pos.y+p->vel.y*time_step;
		p->pos.z=p->pos.z+p->vel.z*time_step;

		if(p->pos.x >= world_size.x-BOUNDARY)
		{
			p->vel.x=p->vel.x*wall_damping;
			p->pos.x=world_size.x-BOUNDARY;
		}

		if(p->pos.x < 0.0f)
		{
			p->vel.x=p->vel.x*wall_damping;
			p->pos.x=0.0f;
		}

		if(p->pos.y >= world_size.y-BOUNDARY)
		{
			p->vel.y=p->vel.y*wall_damping;
			p->pos.y=world_size.y-BOUNDARY;
		}

		if(p->pos.y < 0.0f)
		{
			p->vel.y=p->vel.y*wall_damping;
			p->pos.y=0.0f;
		}

		if(p->pos.z >= world_size.z-BOUNDARY)
		{
			p->vel.z=p->vel.z*wall_damping;
			p->pos.z=world_size.z-BOUNDARY;
		}

		if(p->pos.z < 0.0f)
		{
			p->vel.z=p->vel.z*wall_damping;
			p->pos.z=0.0f;
		}

		p->ev.x=(p->ev.x+p->vel.x)/2;
		p->ev.y=(p->ev.y+p->vel.y)/2;
		p->ev.z=(p->ev.z+p->vel.z)/2;
	}
}

int3 SPHSystem::calc_cell_pos(vec3 p)
{
	int3 cell_pos;
	cell_pos.x = int(floor((p.x) / cell_size));
	cell_pos.y = int(floor((p.y) / cell_size));
	cell_pos.z = int(floor((p.z) / cell_size));

    return cell_pos;
}

uint SPHSystem::calc_cell_hash(int3 cell_pos)
{
	if(cell_pos.x<0 || cell_pos.x>=(int)grid_size.x || cell_pos.y<0 || cell_pos.y>=(int)grid_size.y || cell_pos.z<0 || cell_pos.z>=(int)grid_size.z)
	{
		return (uint)0xffffffff;
	}

	cell_pos.x = cell_pos.x & (grid_size.x-1);  
    cell_pos.y = cell_pos.y & (grid_size.y-1);  
	cell_pos.z = cell_pos.z & (grid_size.z-1);  

	return ((uint)(cell_pos.z))*grid_size.y*grid_size.x + ((uint)(cell_pos.y))*grid_size.x + (uint)(cell_pos.x);
}
