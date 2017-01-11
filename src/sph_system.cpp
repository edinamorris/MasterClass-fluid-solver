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

    numberOfPhases=PHASES;

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
    //step1
    driftVelocity();
    //step2
    advectVolumeFractions();
    //step3
    correctVolumeFraction();
    //step4 - acceleration
    comp_force_adv();
    //step 5 - advect acceleration and velocity
	advection();
    updateColour();
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
                    //phase, position, velocity and miscibility
                    add_particle(0, pos, vel,1);
                }
            }
        }

        for(pos.x=world_size.x*0.0f; pos.x<world_size.x*0.4f; pos.x+=(kernel*0.5f))
        {
            for(pos.y=world_size.y*0.0f; pos.y<world_size.y*0.7f; pos.y+=(kernel*0.5f))
            {
                for(pos.z=world_size.z*0.0f+0.4; pos.z<world_size.z*0.4f+0.4; pos.z+=(kernel*0.5f))
                {
                    //phase, position, velocity and miscibility
                    add_particle(1, pos, vel, 1);
                }
            }
        }

        std::cout<<"number of particles >"<<num_particle<<"\n";
    }
    //dam setup
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
                    //phase, position, velocity and miscibility
                    add_particle(0, pos, vel, 1);
                }
            }
        }

        for(pos.x=world_size.x*0.0f+0.5; pos.x<world_size.x*0.2f+0.5; pos.x+=(kernel*0.5f))
        {
            for(pos.y=world_size.y*0.0f; pos.y<world_size.y*0.7f; pos.y+=(kernel*0.5f))
            {
                for(pos.z=world_size.z*0.0f+0.1; pos.z<world_size.z*0.7f+0.1; pos.z+=(kernel*0.5f))
                {
                    //phase, position, velocity and miscibility
                    add_particle(1, pos, vel, 1);
                }
            }
        }

        std::cout<<"number of particles >"<<num_particle<<"\n";
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

        for(pos.x=world_size.x*0.0f+0.2; pos.x<world_size.x*0.4f+0.2; pos.x+=(kernel*0.5f))
        {
            for(pos.y=world_size.y*0.0f+0.3; pos.y<world_size.y*0.5f+0.3; pos.y+=(kernel*0.5f))
            {
                for(pos.z=world_size.z*0.0f+0.2; pos.z<world_size.z*0.4f+0.2; pos.z+=(kernel*0.5f))
                {
                    //phase, position, velocity and miscibility
                    add_particle(0, pos, vel, 1);
                }
            }
        }
        for(pos.x=world_size.x*0.0f; pos.x<world_size.x*1.0; pos.x+=(kernel*0.5f))
        {
            for(pos.y=world_size.y*0.0f; pos.y<world_size.y*0.25f; pos.y+=(kernel*0.5f))
            {
                for(pos.z=world_size.z*0.0f; pos.z<world_size.z*1.0; pos.z+=(kernel*0.5f))
                {
                    //phase, position, velocity and miscibility
                    add_particle(1, pos, vel, 1);
                }
            }
        }

        std::cout<<"number of particles >"<<num_particle<<"\n";
    }
}

//changing sp particles get added into seperate lists
void SPHSystem::add_particle(int _phase, vec3 pos, vec3 vel, int misc)
{
    Particle *p;
    p=&(mem[num_particle]);

    p->id=num_particle;
    p->pos=pos;
    p->vel=vel;
    p->misc=misc;
    //set values based on preset phase data
    p->colour=updatePhases.phases[_phase].colour;
    p->mass=updatePhases.phases[_phase].individualMass;
    p->visc=updatePhases.phases[_phase].individualVisc;
    p->selfDens=updatePhases.phases[_phase].selfDens;
    p->lplcColour=updatePhases.phases[_phase].self_lplc_color;
    p->restdens=updatePhases.phases[_phase].dens;
    p->dens=p->restdens;
    p->phase=_phase;
    //adding to the phase counter
    updatePhases.addParticle(_phase);

    p->acc=vec3(0,0,0);
    p->ev=vec3(0,0,0);
    for(int i=0; i<numberOfPhases; i++)
    {
        //drift velocity of particle for first and second liquids
        p->driftVelocity[i]=vec3(0,0,0);
    }

    //setting original volume fraction and mass
    p->volumeFraction[_phase]=1.0f;

    p->pres=0.0f;

    p->next=NULL;

    //still add one to overall count
    num_particle++;
}

void SPHSystem::build_table()
{
    Particle *p;
    //Particle **pp;
    int hash;

    for(int i=0; i<int(tot_cell); i++)
	{
		cell[i]=NULL;
	}

        for(int j=0; j<num_particle; j++)
        {
            p=&(mem[j]);
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

void SPHSystem::updateColour()
{
    Particle *p;
    for(int i=0; i<num_particle; i++)
    {
        p=&(mem[i]);
        vec3 updatedColour=vec3();
        for(int phase=0; phase<numberOfPhases; phase++)
        {
            updatedColour+=updatePhases.getColour(phase)*p->volumeFraction[phase];

            p->colour=updatedColour;
        }
    }
}

void SPHSystem::comp_dens_pres()
{
	Particle *p;
	Particle *np;

	int3 cell_pos;
	int3 near_pos;
	uint hash;

    vec3 rel_pos;
	float r2;

    for(int i=0; i<num_particle; i++)
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
                        rel_pos=np->pos-p->pos;
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

//comment out acceleration here, needs to be calculated after drift velocity with multiple fluid method
//Calculate the acceleration of the mixture particle according to Eq. (8). SPH formulations
//of the related terms are provided in Eqs. (19)–(21).
void SPHSystem::comp_force_adv()
{

    //MULTIPLE FLUID PAPER EQUATIONS - NOT WORKING WITH OTHER PARTS
    /*Particle *p;
    for(int i=0; i<num_particle; i++)
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
            for(int phase=0; phase<numberOfPhases; phase++)
            {
                summationPhase+=neighbour->restdens*(neighbour->volumeFraction[phase]*neighbour->driftVelocity[phase]*(neighbour->driftVelocity[phase].dot(smoothingKernel))+
                        (p->volumeFraction[phase]*p->driftVelocity[phase]*(p->driftVelocity[phase].dot(smoothingKernel))));
            }

            correctMomentum+=neighbour->mass/neighbour->dens*summationPhase;

        }
        thirdPhase=correctMomentum/p->restdens;

        p->acc=(-1*firstPhase+secondPhase+thirdPhase)+gravity;
        //std::cout<<"acceleration x > "<<p->acc.x<<"\n";
        //std::cout<<"acceleration y > "<<p->acc.y<<"\n";
        //std::cout<<"acceleration z > "<<p->acc.z<<"\n";
    }*/

    //WORKS
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

    for(int i=0; i<num_particle; i++)
	{
		p=&(mem[i]); 
		cell_pos=calc_cell_pos(p->pos);

        p->acc=vec3(0,0,0);

        grad_color=vec3(0,0,0);
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
                        rel_pos=p->pos-np->pos;
						r2=rel_pos.x*rel_pos.x+rel_pos.y*rel_pos.y+rel_pos.z*rel_pos.z;

						if(r2 < kernel_2 && r2 > INF)
						{
							r=sqrt(r2);
                            V=p->mass/np->dens/2;
							kernel_r=kernel-r;

							pres_kernel=spiky_value * kernel_r * kernel_r;
                            temp_force=V * (p->pres+np->pres) * pres_kernel;
                            p->acc=p->acc-rel_pos*temp_force/r;

                                         //estimated velocity
                            rel_vel=np->ev-p->ev;

							visc_kernel=visco_value*(kernel-r);
                            temp_force=V * p->visc * visc_kernel;
                            p->acc=p->acc + rel_vel*temp_force;

							float temp=(-1) * grad_poly6 * V * pow(kernel_2-r2, 2);
                            grad_color += temp * rel_pos;
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
            p->acc+=(surf_coe * lplc_color * grad_color / p->surf_norm);
		}
    }
}

//needs some work - readjusting based on new findings
void SPHSystem::driftVelocity()
{
    Particle *p;
    //for each phase AND each particle, i.e. each particle will have two drift velocities attached to it
    for(int i=0; i<num_particle; i++)
    {
        p=&(mem[i]);
        for(int phase=0; phase<numberOfPhases; phase++)
        {
            p->driftVelocity[phase]=vec3(0,0,0);
            float strengthFac=1e-8;
            //set to 0 for no diffuse effect
            float diffuseConst=1e-5f;
            vec3 firstPhase=vec3(0,0,0);
            vec3 secondPhase=vec3(0,0,0);
            vec3 thirdPhase=vec3(0,0,0);

            //FIRST PHASE OF EQUATION - accounts for the inertia effect
            float inertia=0;
            //looping through all other phases j represent k' here
            for(int j=0; j<numberOfPhases; j++)
            {
                if(j==phase)
                {
                    continue;
                }
                //sum of all other phases (volumefraction of this particle related to current k' phase * phase density)/the mix density and then * k' phase rest density
                inertia+=((p->volumeFraction[j]*updatePhases.getDensity(j))/p->restdens)*updatePhases.getDensity(j);
            }

            //rest density of current looping phase used here
            firstPhase=strengthFac*(updatePhases.getDensity(phase)-inertia)*(gravity-p->acc);

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
                spikyKernelGrad=spikyGrad(p->pos, neighbour->pos);

                //calculating first interpolated pressure (for individual particles) - using the spiky for other derivatives of smoothing function
                if((p->misc==1) && (neighbour->misc==1))
                {
                    interpolatedPressureInd+= (neighbour->mass/neighbour->dens)*(neighbour->pres*neighbour->volumeFraction[phase]-p->pres*p->volumeFraction[phase])*spikyKernelGrad;
                }
                else
                {
                    interpolatedPressureInd = (neighbour->mass/neighbour->dens)*(neighbour->pres-p->pres)*spikyKernelGrad;
                }
            }

            //interpolated pressure accumulated for all phases [j]
            vec3 pressure=vec3(0,0,0);
            for(int j=0; j<numberOfPhases; j++)
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
                    spikyKernelGrad=spikyGrad(p->pos, neighbour->pos);

                    //interpolated pressure of all phases
                    vec3 interpolatedPressureAll;

                    if((p->misc==1) && (neighbour->misc==1))
                    {
                        interpolatedPressureAll=(neighbour->mass/neighbour->dens)*(neighbour->pres*neighbour->volumeFraction[j]-p->pres*p->volumeFraction[j])*spikyKernelGrad;
                    }
                    else
                    {
                        interpolatedPressureAll=(neighbour->mass/neighbour->dens)*(neighbour->pres-p->pres)*spikyKernelGrad;
                    }

                    pressure+=((p->volumeFraction[j]*updatePhases.getDensity(j))/p->restdens)*interpolatedPressureAll;
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
                spikyKernelGrad=spikyGrad(p->pos, neighbour->pos);

                //calculating first interpolated volume fraction (for this phase) - using the spiky for other derivatives of smoothing function
                interpolatedVolFracInd+=(neighbour->mass/neighbour->dens)*(neighbour->volumeFraction[phase]-p->volumeFraction[phase])*spikyKernelGrad;
            }

            vec3 diffusion=vec3(0,0,0);
            //j is k' here
            for(int j=0; j<numberOfPhases; j++)
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
                    spikyKernelGrad=spikyGrad(p->pos, neighbour->pos);

                    //interpolated volume fraction of all phases [j]                                volume fraction for k' phase
                    vec3 interpolatedVolumeAll=(neighbour->mass/neighbour->dens)*(neighbour->volumeFraction[j]-p->volumeFraction[j])*spikyKernelGrad;

                    vec3 volumeMixVolFrac;
                    //nan if dividing by 0
                    if(p->volumeFraction[j]==0)
                    {
                        volumeMixVolFrac=vec3(0,0,0);
                    }
                    else
                    {
                        volumeMixVolFrac=interpolatedVolumeAll/p->volumeFraction[j];
                    }

                    //sum of all other phases (volumefraction of phase * phase density)/the mix density and then multiply by
                    //gradient of volume fraction/volume fraction of k' phase                volume fraction here uses k' ie. j
                    diffusion+=((p->volumeFraction[j]*updatePhases.getDensity(j)/p->restdens)*(volumeMixVolFrac));
                }
            }

            vec3 VolFrac=vec3(0.0f,0.0f,0.0f);

            //makes drift velocity nan if volume fraction is 0, cannot divide by 0
            if(p->volumeFraction[phase]!=0)
            {
                VolFrac=interpolatedVolFracInd/p->volumeFraction[phase];
            }

                                                                           // - sum of other phases values
            thirdPhase=diffuseConst*((VolFrac)-diffusion);

            //drift velocity - second phase is where pressure relationship will affect the misc/immiscibility of the mixture
            //each particle will have a drift velocity calculated for however many liquids (phases) there are
            p->driftVelocity[phase]=firstPhase-secondPhase-thirdPhase;
        }
    }
}

//Advect volume fraction values according to Eq. (7), where the
//relevant SPH formulations are given in Eqs. (17) and (18).
void SPHSystem::advectVolumeFractions()
{
    Particle *p;
    Particle *neighbour;

    for(int i=0; i<num_particle; i++)
    {
        p=&(mem[i]);
        for(int phase=0; phase<numberOfPhases; phase++)
        {
            float eqs17=0.0f;
            float eqs18=0.0f;
            for(int j=0; j<int(p->neighbour.size()); j++)
            {
                neighbour=(p->neighbour[j]);
                //dont want to compare same particle
                if(p->id==neighbour->id)
                {
                    continue;
                }

                vec3 spikyGradient=spikyGrad(p->pos, neighbour->pos);

                //average volume fraction value between particle and it's neighbours
                float diffvolfrac=neighbour->volumeFraction[phase]+p->volumeFraction[phase]/2;
                //collation of drift velocities by phase fraction per particle
                vec3 volumeFracDrift=(neighbour->volumeFraction[phase]*neighbour->driftVelocity[phase]+
                                 p->volumeFraction[phase]*p->driftVelocity[phase]);

                eqs17+=(neighbour->mass/neighbour->dens)*((diffvolfrac)*(neighbour->vel-p->vel).dot(spikyGradient));
                eqs18+=((neighbour->mass/neighbour->dens)*(volumeFracDrift).dot(spikyGradient));
            }
            p->volumeFraction[phase]=p->volumeFraction[phase]-eqs17-eqs18;
        }
    }
}

//Check the bound of volume fraction according to Eq.(3) and, if the bound is invalidated, correct the volume
//fraction within the particle and calculate the pressure adjustment accordingly as described in Section 5.1.
//For particles with a corrected volume fraction, update into the pressure term the pressure adjustment as Eq. (31).
void SPHSystem::correctVolumeFraction()
{
    Particle *p;

    for(int i=0; i<num_particle; i++)
    {
        float volumeFraction=0.0f;
        p=&(mem[i]);
        for(int phase=0; phase<numberOfPhases; phase++)
        {
            //correcting any unforseen errors
            bool notANumberCheck=std::isnan(p->volumeFraction[phase]);
            if(notANumberCheck==true)
            {
                //resets the particles volume fraction to original starting values
                if(phase==p->phase)
                {
                    p->volumeFraction[phase]=1.0f;
                }
                else
                {
                    p->volumeFraction[phase]=0.0f;
                }
            }
            if(p->volumeFraction[phase]<0)
            {
                p->volumeFraction[phase]=0.0f;
            }
            //adding all volume fraction values of a particle for every phase
            volumeFraction+=p->volumeFraction[phase];
        }
        if(volumeFraction>=0.0f&&volumeFraction<=0.01)
        {
            p->volumeFraction[p->phase]=1.0f;
            for(int k=0; k<numberOfPhases; k++)
            {
                if(p->phase!=k)
                {
                    p->volumeFraction[k]=0.0f;
                }
            }
        }
        //checking if values add up
        else if(volumeFraction!=1.0)
        {
            float pressureAdjustment=0.0f;
            float scale=1.0/volumeFraction;
            for(int phase=0; phase<numberOfPhases; phase++)
            {
                p->volumeFraction[phase]=p->volumeFraction[phase]*scale;
                pressureAdjustment+=-gas_constant*updatePhases.getDensity(phase)*(p->volumeFraction[phase]-p->prevVolumeFraction[phase]);
            }
            //CRAY CRAY
            //p->pres=p->pres+pressureAdjustment;
        }
    }
}

//used for calculations involving the smoothing kernel function
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

    for(int i=0; i<num_particle; i++)
	{
		p=&(mem[i]);

        p->vel=p->vel+p->acc*time_step/p->dens+gravity*time_step;

        p->pos=p->pos+p->vel*time_step;

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

        p->ev=(p->ev+p->vel)/2;
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
