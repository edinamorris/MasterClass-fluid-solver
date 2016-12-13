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

    //temporary stroage of phase rest density
    restDensityArray[0]=1300;
    restDensityArray[1]=1000;

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
    updateColour();
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

        std::cout<<"number of particles >"<<num_particle<<"\n";
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

        std::cout<<"number of particles >"<<num_particle<<"\n";
    }
}

//changing sp particles get added into seperate lists
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
    for(int i=0; i<numberOfPhases; i++)
    {
        //drift velocity of particle for first and second liquids
        p->driftVelocity[i].x=0.0f;
        p->driftVelocity[i].y=0.0f;
        p->driftVelocity[i].z=0.0f;
    }

    //TEMP initial volume fraction
    if(_phase==0)
    {
        //if particle starts in liquid 0 then it will be completely in liquid 0, no mixing has occurerd
        p->volumeFraction[0]=0.99f;
        p->volumeFraction[1]=0.01f;
    }
    else if(_phase==1)
    {
        //opposite if starts in liquid 1
        p->volumeFraction[0]=0.01f;
        p->volumeFraction[1]=0.99f;
    }

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
        p->updatedColour=vec3(0,0,0);
    }
    for(int phase=0; phase<numberOfPhases; phase++)
    {
        for(int i=0; i<num_particle; i++)
        {
            p=&(mem[i]);
            p->updatedColour+=updatePhases.getColour(phase)*p->volumeFraction[phase];

            p->colour=p->updatedColour;
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

//comment out acceleration here, needs to be calculated after drift velocity with multiple fluid method
//Calculate the acceleration of the mixture particle according to Eq. (8). SPH formulations
//of the related terms are provided in Eqs. (19)–(21).
/*void SPHSystem::comp_force_adv()
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

    for(int i=0; i<num_particle; i++)
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

                                         //estimated velocity
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
}*/

//Calculate the acceleration of the mixture particle according to Eq. (8). SPH formulations
//of the related terms are provided in Eqs. (19)–(21).
void SPHSystem::comp_force_adv()
{
    Particle *p;
    for(int i=0; i<num_particle; i++)
    {
        p=&(mem[i]);

        //p->acc.x=0.0f;
        //p->acc.y=0.0f;
        //p->acc.z=0.0f;
        vec3 firstPhase;
        vec3 secondPhase;
        vec3 thirdPhase;
        firstPhase.x=0.0f;
        firstPhase.y=0.0f;
        firstPhase.z=0.0f;
        secondPhase.x=0.0f;
        secondPhase.y=0.0f;
        secondPhase.z=0.0f;
        thirdPhase.x=0.0f;
        thirdPhase.y=0.0f;
        thirdPhase.z=0.0f;

        //first phase - equation 20
        vec3 summation;
        summation.x=0.0f;
        summation.y=0.0f;
        summation.z=0.0f;
        float interpolatedDensity=0.0f;
        vec3 smoothingKernel;
        smoothingKernel.x=0.0f;
        smoothingKernel.y=0.0f;
        smoothingKernel.z=0.0f;
        float poly6Kernel;
        //neighbourhood particles
        for(int j=0; j<int(p->neighbour.size()); j++)
        {
            Particle *neighbour;
            neighbour=(p->neighbour[j]);

            // Wij - for density interpolation - changed to scalar after teemu's discovery
            poly6Kernel=poly6(p->pos, neighbour->pos);

            //for all other calculations involving derivative of smoothing function
            smoothingKernel=spikyGrad(p->pos, neighbour->pos);

            //now scalar value
            interpolatedDensity=(neighbour->mass*poly6Kernel);

            summation.x+=neighbour->mass*(p->pres+neighbour->pres/2*interpolatedDensity)*smoothingKernel.x;
            summation.y+=neighbour->mass*(p->pres+neighbour->pres/2*interpolatedDensity)*smoothingKernel.y;
            summation.z+=neighbour->mass*(p->pres+neighbour->pres/2*interpolatedDensity)*smoothingKernel.z;
        }
        firstPhase.x=summation.x/p->restdens;
        firstPhase.y=summation.y/p->restdens;
        firstPhase.z=summation.z/p->restdens;

        //second phase - equation 21
        vec3 summation2;
        summation2.x=0.0f;
        summation2.y=0.0f;
        summation2.z=0.0f;
        //neighbourhood particles
        for(int j=0; j<int(p->neighbour.size()); j++)
        {
            Particle *neighbour;
            neighbour=(p->neighbour[j]);

            // Wij - for density interpolation - changed to scalar after teemu's discovery
            poly6Kernel=poly6(p->pos, neighbour->pos);

            //for all other calculations involving derivative of smoothing function
            smoothingKernel=spikyGrad(p->pos, neighbour->pos);

            //now scalar value
            interpolatedDensity=(neighbour->mass*poly6Kernel);

            //checks
            if(interpolatedDensity==0)
            {
                interpolatedDensity=0.1;
            }
            vec3 positionDiff;
            positionDiff.x=neighbour->pos.x-p->pos.x;
            positionDiff.y=neighbour->pos.y-p->pos.y;
            positionDiff.z=neighbour->pos.z-p->pos.z;
            if(positionDiff.x==0)
            {
                positionDiff.x=0.01;
            }
            else if(positionDiff.y==0)
            {
                positionDiff.y=0.01;
            }
            else if(positionDiff.z==0)
            {
                positionDiff.z=0.01;
            }

            summation2.x=neighbour->mass/interpolatedDensity*(p->visc+neighbour->visc)*(neighbour->vel.x-p->vel.x)*
                    ((neighbour->pos.x-p->pos.x)*(smoothingKernel.x)/(positionDiff.x)*(neighbour->pos.x-p->pos.x));
            summation2.y=neighbour->mass/interpolatedDensity*(p->visc+neighbour->visc)*(neighbour->vel.y-p->vel.y)*
                    ((neighbour->pos.y-p->pos.y)*(smoothingKernel.y)/(positionDiff.y)*(neighbour->pos.y-p->pos.y));
            summation2.z=neighbour->mass/interpolatedDensity*(p->visc+neighbour->visc)*(neighbour->vel.z-p->vel.z)*
                    ((neighbour->pos.z-p->pos.z)*(smoothingKernel.z)/(positionDiff.z)*(neighbour->pos.z-p->pos.z));

            std::cout<<"testing x > "<<positionDiff.x<<"\n";
            std::cout<<"testing y > "<<positionDiff.y<<"\n";
            std::cout<<"testing z > "<<positionDiff.z<<"\n";
        }
        secondPhase.x=summation2.x/p->restdens;
        secondPhase.y=summation2.y/p->restdens;
        secondPhase.z=summation2.z/p->restdens;

        //std::cout<<"testing x > "<<secondPhase.x<<"\n";
        //std::cout<<"testing y > "<<secondPhase.y<<"\n";
        //std::cout<<"testing z > "<<secondPhase.z<<"\n";

        //third phase - equation 19
        vec3 summation3;
        summation3.x=0.0f;
        summation3.y=0.0f;
        summation3.z=0.0f;
        vec3 summationPhase;
        //neighbourhood particles
        for(int j=0; j<int(p->neighbour.size()); j++)
        {
            Particle *neighbour;
            neighbour=(p->neighbour[j]);

            // Wij - for density interpolation - changed to scalar after teemu's discovery
            poly6Kernel=poly6(p->pos, neighbour->pos);

            //for all other calculations involving derivative of smoothing function
            smoothingKernel=spikyGrad(p->pos, neighbour->pos);

            //now scalar value
            interpolatedDensity=(neighbour->mass*poly6Kernel);

            //resetting for new neighbour particle
            summationPhase.x=0.0f;
            summationPhase.y=0.0f;
            summationPhase.z=0.0f;

            //go through phases
            for(int phase=0; phase<numberOfPhases; phase++)
            {
                summationPhase.x+=neighbour->restdens*(neighbour->volumeFraction[phase]*neighbour->driftVelocity[phase].x*(neighbour->driftVelocity[phase].x*smoothingKernel.x))+
                        (p->volumeFraction[phase]*p->driftVelocity[phase].x*(p->driftVelocity[phase].x*smoothingKernel.x));
                summationPhase.y+=neighbour->restdens*(neighbour->volumeFraction[phase]*neighbour->driftVelocity[phase].y*(neighbour->driftVelocity[phase].y*smoothingKernel.z))+
                        (p->volumeFraction[phase]*p->driftVelocity[phase].y*(p->driftVelocity[phase].y*smoothingKernel.y));
                summationPhase.z+=neighbour->restdens*(neighbour->volumeFraction[phase]*neighbour->driftVelocity[phase].z*(neighbour->driftVelocity[phase].y*smoothingKernel.z))+
                        (p->volumeFraction[phase]*p->driftVelocity[phase].z*(p->driftVelocity[phase].z*smoothingKernel.z));
            }

            summation3.x+=neighbour->mass/interpolatedDensity*summationPhase.x;
            summation3.y+=neighbour->mass/interpolatedDensity*summationPhase.y;
            summation3.z+=neighbour->mass/interpolatedDensity*summationPhase.z;

        }
        thirdPhase.x=-summation3.x/p->restdens;
        thirdPhase.y=-summation3.y/p->restdens;
        thirdPhase.z=-summation3.z/p->restdens;

        p->acc.x=gravity.x-1*(-firstPhase.x+gravity.x+secondPhase.x+thirdPhase.x);
        p->acc.y=gravity.y-1*(-firstPhase.y+gravity.y+secondPhase.y+thirdPhase.y);
        p->acc.z=gravity.z-1*(-firstPhase.z+gravity.z+secondPhase.z+thirdPhase.z);
        //std::cout<<"acceleration x > "<<p->acc.x<<"\n";
        //std::cout<<"acceleration x > "<<p->acc.y<<"\n";
        //std::cout<<"acceleration x > "<<p->acc.z<<"\n";
    }
}

vec3 SPHSystem::MixturePressure()
{
    vec3 mixPressure;
    //= sum( volume fraction of phase k * pressure gradient of phase k)
    for(int phase=0; phase<numberOfPhases; phase++)
    {
        mixPressure.x=5.0f;
        mixPressure.y=5.0f;
        mixPressure.z=5.0f;
    }

    return mixPressure;
}

//needs some work - readjusting based on new findings
void SPHSystem::driftVelocity()
{
    vec3 mixturePressure=MixturePressure();
    Particle *p;
    //for each phase AND each particle, i.e. each particle will have two drift velocities attached to it
    for(int i=0; i<num_particle; i++)
    {
        p=&(mem[i]);
        for(int phase=0; phase<numberOfPhases; phase++)
        {
            p->driftVelocity[phase].x=0;
            p->driftVelocity[phase].y=0;
            p->driftVelocity[phase].z=0;
            //set to 0 for immiscibility
            float strengthFac=0.00000001f;
            //set to 0 for no diffuse effect
            //float diffuseConst=0.0f;
            float diffuseConst=0.000001f;
            vec3 firstPhase;
            vec3 secondPhase;
            vec3 thirdPhase;
            firstPhase.x=0;
            firstPhase.y=0;
            firstPhase.z=0;
            secondPhase.x=0;
            secondPhase.y=0;
            secondPhase.z=0;
            thirdPhase.x=0;
            thirdPhase.y=0;
            thirdPhase.z=0;

            //FIRST PHASE OF EQUATION
            //sigma k, all other phases
            float summation=0;
            float mixDensity=0;
            //looping through all phases j represent k' here
            for(int j=0; j<numberOfPhases; j++)
            {
                //mixture Density
                mixDensity+=(p->volumeFraction[j]*restDensityArray[j]);
            }
            for(int j=0; j<numberOfPhases; j++)
            {
                //sum of all other phases (volumefraction of this particle related to current k' phase * phase density)/the mix density and then * k' phase rest density
                summation+=((p->volumeFraction[j]*restDensityArray[j])/mixDensity)*restDensityArray[j];
            }

            //rest density of current looping phase used here
            firstPhase.x=strengthFac*(restDensityArray[phase]-summation)*p->acc.x;
            firstPhase.y=strengthFac*(restDensityArray[phase]-summation)*p->acc.y;
            firstPhase.z=strengthFac*(restDensityArray[phase]-summation)*p->acc.y;

            //SECOND PHASE OF EQUATION
            float poly6Kernel;
            vec3 spikyKernelGrad;
            vec3 interpolatedPressureInd;
            //INTERPOLATED PRESSURE FOR INDIVIDUAL PARTICLE
            interpolatedPressureInd.x=0.0f;
            interpolatedPressureInd.y=0.0f;
            interpolatedPressureInd.z=0.0f;
            //neighbourhood particles
            for(int j=0; j<int(p->neighbour.size()); j++)
            {
                Particle *neighbour;
                neighbour=(p->neighbour[j]);

                float interpolatedDensity;
                vec3 calc_1;
                calc_1.x=0;
                calc_1.y=0;
                calc_1.z=0;
                vec3 phasePressureSelf;
                vec3 phasePressureN;
                //the ability to turn misicibility on or off
                if(misc==1)
                {
                    //mixture pressure * individual particles volume fraction
                    phasePressureSelf.x=mixturePressure.x*p->volumeFraction[phase];
                    phasePressureSelf.y=mixturePressure.y*p->volumeFraction[phase];
                    phasePressureSelf.z=mixturePressure.z*p->volumeFraction[phase];
                    phasePressureN.x=mixturePressure.x*neighbour->volumeFraction[phase];
                    phasePressureN.y=mixturePressure.y*neighbour->volumeFraction[phase];
                    phasePressureN.z=mixturePressure.z*neighbour->volumeFraction[phase];
                }
                //mixture pressure ie. same value so equation turns to 0
                else
                {
                    phasePressureSelf.x=mixturePressure.x;
                    phasePressureSelf.y=mixturePressure.y;
                    phasePressureSelf.z=mixturePressure.z;
                    phasePressureN.x=mixturePressure.x;
                    phasePressureN.y=mixturePressure.y;
                    phasePressureN.z=mixturePressure.z;
                }

                // Wij - for density interpolation - changed to scalar after teemu's discovery
                poly6Kernel=poly6(p->pos, neighbour->pos);

                //for all other calculations involving derivative of smoothing function
                spikyKernelGrad=spikyGrad(p->pos, neighbour->pos);

                //now scalar value
                interpolatedDensity=(neighbour->mass*poly6Kernel);

                //calculating first interpolated pressure (for individual particles) - using the spiky for other derivatives of smoothing function
                if((interpolatedDensity!=0))
                {
                    calc_1.x = (neighbour->mass/interpolatedDensity)*(phasePressureN.x-phasePressureSelf.x)*spikyKernelGrad.x;
                    calc_1.y = (neighbour->mass/interpolatedDensity)*(phasePressureN.y-phasePressureSelf.y)*spikyKernelGrad.y;
                    calc_1.z = (neighbour->mass/interpolatedDensity)*(phasePressureN.z-phasePressureSelf.z)*spikyKernelGrad.z;
                }

                //accumulates to form the sum of all this particles neighbours
                interpolatedPressureInd.x+=calc_1.x;
                interpolatedPressureInd.y+=calc_1.y;
                interpolatedPressureInd.z+=calc_1.z;
            }

            //INTERPOLATED PRESSURE OF ALL PHASES
            //for mixture
            vec3 phasePressureSelf;
            vec3 phasePressureNeigh;

            vec3 summation2;
            summation2.x=0.0f;
            summation2.y=0.0f;
            summation2.z=0.0f;
            float mixDensity2=0.0f;
            //looping through all phases j represents k' here
            for(int j=0; j<numberOfPhases; j++)
            {
                //mixture Density
                mixDensity2+=(p->volumeFraction[j]*restDensityArray[j]);
            }
            for(int j=0; j<numberOfPhases; j++)
            {
                for(int neigh=0; neigh<int(p->neighbour.size()); neigh++)
                {
                    Particle *neighbour_mix;
                    neighbour_mix=(p->neighbour[neigh]);

                    if(misc==1)
                    {
                        //mixture density * volume fraction value of current phase loop
                        phasePressureSelf.x=mixturePressure.x*p->volumeFraction[j];
                        phasePressureSelf.y=mixturePressure.y*p->volumeFraction[j];
                        phasePressureSelf.z=mixturePressure.z*p->volumeFraction[j];//self and neighbour is current phase that's being looped through - k'
                        phasePressureNeigh.x=mixturePressure.x*neighbour_mix->volumeFraction[j];
                        phasePressureNeigh.y=mixturePressure.y*neighbour_mix->volumeFraction[j];
                        phasePressureNeigh.z=mixturePressure.z*neighbour_mix->volumeFraction[j];

                    }
                    else
                    {
                        //mixture pressure ie. same value so equation turns to 0
                        phasePressureSelf.x=mixturePressure.x;
                        phasePressureSelf.y=mixturePressure.y;
                        phasePressureSelf.z=mixturePressure.z;
                        phasePressureNeigh.x=mixturePressure.x;
                        phasePressureNeigh.y=mixturePressure.y;
                        phasePressureNeigh.z=mixturePressure.z;
                    }

                    //Wij - for density interpolation
                    poly6Kernel=poly6(p->pos, neighbour_mix->pos);

                    //for all other calculations involving derivative of smoothing function
                    spikyKernelGrad=spikyGrad(p->pos, neighbour_mix->pos);

                    //interpolated pressure of all phases
                    vec3 interpolatedPressureMix;
                    //changed to be a scalar
                    float interpolatedDensityMix;
                    interpolatedDensityMix=neighbour_mix->mass*poly6Kernel;

                    //pressure of other phase minus pressure of current phase
                    if((interpolatedDensityMix!=0))
                    {
                        interpolatedPressureMix.x=(neighbour_mix->mass/interpolatedDensityMix)*(phasePressureNeigh.x-phasePressureSelf.x)*spikyKernelGrad.x;
                        interpolatedPressureMix.y=(neighbour_mix->mass/interpolatedDensityMix)*(phasePressureNeigh.y-phasePressureSelf.y)*spikyKernelGrad.y;
                        interpolatedPressureMix.z=(neighbour_mix->mass/interpolatedDensityMix)*(phasePressureNeigh.z-phasePressureSelf.z)*spikyKernelGrad.z;
                    }

                    //sum of all phases (volumefraction of k' phase * k' phase rest density)/the mix density and then * phase density
                    summation2.x+=((p->volumeFraction[j]*restDensityArray[j])/mixDensity2)*interpolatedPressureMix.x;
                    summation2.y+=((p->volumeFraction[j]*restDensityArray[j])/mixDensity2)*interpolatedPressureMix.y;
                    summation2.z+=((p->volumeFraction[j]*restDensityArray[j])/mixDensity2)*interpolatedPressureMix.z;
                }
            }

                                                       // - sum of other phases values
            secondPhase.x=strengthFac*(interpolatedPressureInd.x-(summation2.x));
            secondPhase.y=strengthFac*(interpolatedPressureInd.y-(summation2.y));
            secondPhase.z=strengthFac*(interpolatedPressureInd.z-(summation2.z));

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

                float interpolatedDensityInd;
                vec3 calc_1;

                // /\Wij - for density interpolation
                poly6Kernel=poly6(p->pos, neighbour->pos);

                //for all other calculations involving derivative of smoothing function
                spikyKernelGrad=spikyGrad(p->pos, neighbour->pos);

                //following standard method for interpolated density pi=sigma(massj*smoothingkernelfunction)
                interpolatedDensityInd=(neighbour->mass*poly6Kernel);

                //std::cout<<"mass >"<<interpolatedDensityInd.x<<"\n";

                if((interpolatedDensityInd!=0.0))
                {
                    //calculating first interpolated volume fraction (for individual particles) - using the spiky for other derivatives of smoothing function
                                                                          //currently these are the same value
                    calc_1.x = (neighbour->mass/interpolatedDensityInd)*(neighbour->volumeFraction[phase]-p->volumeFraction[phase])*spikyKernelGrad.x;
                    calc_1.y = (neighbour->mass/interpolatedDensityInd)*(neighbour->volumeFraction[phase]-p->volumeFraction[phase])*spikyKernelGrad.y;
                    calc_1.z = (neighbour->mass/interpolatedDensityInd)*(neighbour->volumeFraction[phase]-p->volumeFraction[phase])*spikyKernelGrad.z;
                }

                //accumulates to form the sum of all this particles neighbours
                interpolatedVolFracInd.x+=calc_1.x;
                interpolatedVolFracInd.y+=calc_1.y;
                interpolatedVolFracInd.z+=calc_1.z;
            }

            vec3 summation3;
            summation3.x=0.0f;
            summation3.y=0.0f;
            summation3.z=0.0f;
            float mixDensity3=0;
            //looping through all phases - j is k' here
            for(int j=0; j<numberOfPhases; j++)
            {
                //mixture Density
                mixDensity3+=(p->volumeFraction[j]*restDensityArray[j]);
            }
            //j is k' here
            for(int j=0; j<numberOfPhases; j++)
            {
                for(int neigh=0; neigh<int(p->neighbour.size()); neigh++)
                {
                    Particle *neighbour_mix;
                    neighbour_mix=(p->neighbour[neigh]);

                    //Wij - for density interpolation
                    poly6Kernel=poly6(p->pos, neighbour_mix->pos);

                    //for all other calculations involving derivative of smoothing function
                    spikyKernelGrad=spikyGrad(p->pos, neighbour_mix->pos);

                    //interpolated pressure of other phases (summed of all OTHER phases rather than neighbours) - in this case only 2 phases needs changing
                    vec3 interpolatedVolumeMix;
                    float interpolatedDensityMix;
                    //using poly6 for interpolated density
                    interpolatedDensityMix=neighbour_mix->mass*poly6Kernel;

                    if((interpolatedDensityMix!=0))
                    {
                        //other phase volume fraction - self's volume fraction                              volume fraction for k' phase
                        interpolatedVolumeMix.x=(neighbour_mix->mass/interpolatedDensityMix)*(p->volumeFraction[j]-p->volumeFraction[j])*spikyKernelGrad.x;
                        interpolatedVolumeMix.y=(neighbour_mix->mass/interpolatedDensityMix)*(p->volumeFraction[j]-p->volumeFraction[j])*spikyKernelGrad.y;
                        interpolatedVolumeMix.z=(neighbour_mix->mass/interpolatedDensityMix)*(p->volumeFraction[j]-p->volumeFraction[j])*spikyKernelGrad.z;
                    }

                    vec3 volumeMixVolFrac;
                    //nan if dividing by 0
                    if(p->volumeFraction[j]==0)
                    {
                        volumeMixVolFrac.x=0.0;
                        volumeMixVolFrac.y=0.0;
                        volumeMixVolFrac.z=0.0;
                    }
                    else
                    {
                        volumeMixVolFrac.x=interpolatedVolumeMix.x/p->volumeFraction[j];
                        volumeMixVolFrac.y=interpolatedVolumeMix.y/p->volumeFraction[j];
                        volumeMixVolFrac.z=interpolatedVolumeMix.z/p->volumeFraction[j];
                    }

                    //sum of all other phases (volumefraction of phase * phase density)/the mix density and then multiply by
                    //gradient of volume fraction/volume fraction of k' phase                volume fraction here uses k' ie. j
                    summation3.x+=((p->volumeFraction[j]*restDensityArray[j])/mixDensity3)*(volumeMixVolFrac.x);
                    summation3.y+=((p->volumeFraction[j]*restDensityArray[j])/mixDensity3)*(volumeMixVolFrac.y);
                    summation3.z+=((p->volumeFraction[j]*restDensityArray[j])/mixDensity3)*(volumeMixVolFrac.z);
                }
            }

            vec3 VolFrac;
            //TESTING
            if(p->volumeFraction[phase]==0)
            {
                VolFrac.x=0.0;
                VolFrac.y=0.0;
                VolFrac.z=0.0;
            }
            else
            {
                VolFrac.x=interpolatedVolFracInd.x/p->volumeFraction[phase];
                VolFrac.y=interpolatedVolFracInd.y/p->volumeFraction[phase];
                VolFrac.z=interpolatedVolFracInd.z/p->volumeFraction[phase];
            }

                                                                           // - sum of other phases values
            thirdPhase.x=diffuseConst*((VolFrac.x)-summation3.x);
            thirdPhase.y=diffuseConst*((VolFrac.y)-summation3.y);
            thirdPhase.z=diffuseConst*((VolFrac.z)-summation3.z);

            //drift velocity - second phase is where pressure relationship will affect the misc/immiscibility of the mixture
            //each particle will have a drift velocity calculated for however many liquids (phases) there are
            p->driftVelocity[phase].x=firstPhase.x-secondPhase.x-thirdPhase.x;
            p->driftVelocity[phase].y=firstPhase.y-secondPhase.y-thirdPhase.y;
            p->driftVelocity[phase].z=firstPhase.z-secondPhase.z-thirdPhase.z;

            //std::cout<<"drift velocity x >"<<p->driftVelocity[phase].x<<"\n";
            //std::cout<<"drift velocity y >"<<p->driftVelocity[phase].y<<"\n";
            //std::cout<<"drift velocity z >"<<p->driftVelocity[phase].z<<"\n";
        }
    }
}

//Advect volume fraction values according to Eq. (7), where the
//relevant SPH formulations are given in Eqs. (17) and (18).
void SPHSystem::advectVolumeFractions()
{
    Particle *p;

    for(int i=0; i<num_particle; i++)
    {
        p=&(mem[i]);
        for(int phase=0; phase<numberOfPhases; phase++)
        {
            float eqs17=0.0f;
            float eqs18=0.0f;
            for(int j=0; j<int(p->neighbour.size()); j++)
            {
                Particle *neighbour;
                neighbour=(p->neighbour[j]);

                float interpolateDensity=neighbour->mass*poly6(p->pos, neighbour->pos);
                vec3 spikyGradient=spikyGrad(p->pos, neighbour->pos);
                float diffvolfrac;
                float volumeFracDrift;
                vec3 neighDV=neighbour->driftVelocity[phase];

                diffvolfrac=neighbour->volumeFraction[phase]+p->volumeFraction[phase]/2==0;
                if(diffvolfrac==0)
                {
                    diffvolfrac=0.01;
                }
                volumeFracDrift=(neighbour->volumeFraction[phase]*(neighDV.dot(neighDV))+
                                 p->volumeFraction[phase]*(neighDV.dot(neighDV)));

                eqs17+=(neighbour->mass/interpolateDensity)*((diffvolfrac)*(neighbour->vel-p->vel).dot(spikyGradient));
                eqs18+=((neighbour->mass/interpolateDensity)*(volumeFracDrift)*(spikyGradient.dot(spikyGradient)));
                //std::cout<<"test value > "<<(neighbour->volumeFraction[phase])<<"\n";
            }
            bool notANumberCheck=std::isnan(p->volumeFraction[phase]-0.001*(eqs17-eqs18));
            bool infinityCheck=std::isinf(p->volumeFraction[phase]-0.001*(eqs17-eqs18));
            if((notANumberCheck==false)&&(infinityCheck==false))
            {
                p->volumeFraction[phase]=p->volumeFraction[phase]-eqs17-eqs18;
            }
            //p->volumeFraction[phase]=p->volumeFraction[phase]-0.001*(eqs17-eqs18);
            //std::cout<<"volume fraction > "<<p->volumeFraction[phase]-0.01*(eqs17-eqs18)<<"\n";
            //std::cout<<"volume fraction > "<<p->volumeFraction[phase]<<"\n";
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
            if(p->volumeFraction[phase]<0)
            {
                p->volumeFraction[phase]=0.001f;
            }
            //adding all volume fraction values of a particle for every phase
            volumeFraction+=p->volumeFraction[phase];
        }

        //checking if values add up - NEEDS CHECKING
        if(volumeFraction!=1.0)
        {
            float scale=1.0/volumeFraction;
            for(int phase=0; phase<numberOfPhases; phase++)
            {
                p->volumeFraction[phase]=p->volumeFraction[phase]*scale;
            }
        }
    }
}

//changed intepolated density to be just the poly6, not gradient after teemu's disocvery of error in printed version of paper
float SPHSystem::poly6(vec3 position_p, vec3 position_n)
{
    vec3 diffPositions = position_p-position_n;
    double r2 = diffPositions.dot(diffPositions);
    double kernelSquared = kernel*kernel;

    return poly6_value * pow(kernelSquared-r2, 3);
}

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

    for(int i=0; i<num_particle; i++)
	{
		p=&(mem[i]);

        //original and now drift velocity is directly affecting acceleration
        p->vel.x=p->vel.x+p->acc.x*time_step/p->dens+gravity.x*time_step;
        p->vel.y=p->vel.y+p->acc.y*time_step/p->dens+gravity.y*time_step;
        p->vel.z=p->vel.z+p->acc.z*time_step/p->dens+gravity.z*time_step;

        //std::cout<<"Velocity X -> "<<p->vel.x<<"\n";
        //std::cout<<"Velocity Y -> "<<p->vel.y<<"\n";
        //std::cout<<"Velocity Z -> "<<p->vel.z<<"\n";

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
