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

#ifdef __APPLE__
#include <GLUT/glut.h>
#include <OpenGL.h>
#else
#include <GL/glut.h>
#endif

SPHSystem::SPHSystem()
{
    max_particle=30000;
	num_particle=0;

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

	time_step=0.003f;
	surf_norm=6.0f;
	surf_coe=0.1f;

	poly6_value=315.0f/(64.0f * PI * pow(kernel, 9));;
	spiky_value=-45.0f/(PI * pow(kernel, 6));
	visco_value=45.0f/(PI * pow(kernel, 6));

	grad_poly6=-945/(32 * PI * pow(kernel, 9));
	lplc_poly6=-945/(8 * PI * pow(kernel, 9));

	kernel_2=kernel*kernel;

	mem=(Particle *)malloc(sizeof(Particle)*max_particle);
	cell=(Particle **)malloc(sizeof(Particle *)*tot_cell);

	sys_running=0;

    misc=0;

    //first phase (liquid) - water - more dense, lower visc - less mass
    colour_1.x=0.2f;
    colour_1.y=0.8f;
    colour_1.z=1.0f;
    //seperate values for each phase
    individualMass_1 =0.02f;
    individualVisc_1 = 6.5f;
    self_dens_1=individualMass_1*poly6_value*pow(kernel, 6);
    self_lplc_color_1=lplc_poly6*individualMass_1*kernel_2*(0-3/4*kernel_2);
    //change rest density value individually - different liquids
    dens_1=1300;

    //second phase (liquid) - oil - less dense much higher viscocity - more mass
    colour_2.x=0.294118f;
    colour_2.y=0.0f;
    colour_2.z=0.509804f;
    //seperate values for each phase
    individualMass_2 =0.03f;
    individualVisc_2 = 350.0f;
    self_dens_2=individualMass_2*poly6_value*pow(kernel, 6);
    self_lplc_color_2=lplc_poly6*individualMass_2*kernel_2*(0-3/4*kernel_2);
    //change rest density value individually - different liquids
    dens_2=1000;

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
	advection();
}

void SPHSystem::init_system()
{
	float3 pos;
	float3 vel;

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
                add_particle(1, pos, vel);
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
                add_particle(2, pos, vel);
            }
        }
    }

	printf("Init Particle: %u\n", num_particle);
}

//blue lower density
void SPHSystem::damnScenario()
{
    float3 pos;
    float3 vel;

    vel.x=0.0f;
    vel.y=0.0f;
    vel.z=0.0f;

    volume_fraction_1=0.5;
    volume_fraction_2=0.5;

    for(pos.x=world_size.x*0.0f; pos.x<world_size.x*0.2f; pos.x+=(kernel*0.5f))
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

    for(pos.x=world_size.x*0.0f+0.5; pos.x<world_size.x*0.2f+0.5; pos.x+=(kernel*0.5f))
    {
        for(pos.y=world_size.y*0.0f; pos.y<world_size.y*0.7f; pos.y+=(kernel*0.5f))
        {
            for(pos.z=world_size.z*0.0f+0.1; pos.z<world_size.z*0.7f+0.1; pos.z+=(kernel*0.5f))
            {
                //phase, position and velocity
                add_particle(2, pos, vel);
            }
        }
    }

    printf("Init Particle: %u\n", num_particle);
}

//blue = bottom, purple = drop
void SPHSystem::dropScenario()
{
    float3 pos;
    float3 vel;

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
                add_particle(2, pos, vel);
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

    printf("Init Particle: %u\n", num_particle);
}

void SPHSystem::add_particle(int _phase, float3 pos, float3 vel)
{
	Particle *p=&(mem[num_particle]);

	p->id=num_particle;

    if(_phase==1)
    {
        p->pos=pos;
        p->vel=vel;
        p->colour=colour_1;
        p->mass=individualMass_1;
        p->visc=individualVisc_1;
        p->selfDens=self_dens_1;
        p->lplcColour=self_lplc_color_1;

        //individual rest densities for each phase
        p->restdens=dens_1;
        p->dens=p->restdens;
    }
    else if(_phase==2)
    {
        p->pos=pos;
        p->vel=vel;
        p->colour=colour_2;
        p->mass=individualMass_2;
        p->visc=individualVisc_2;
        p->selfDens=self_dens_2;
        p->lplcColour=self_lplc_color_2;

        //individual rest densities for each phase
        p->restdens=dens_2;
        p->dens=p->restdens;
    }

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

	num_particle++;
}

void SPHSystem::build_table()
{
	Particle *p;
	uint hash;

	for(uint i=0; i<tot_cell; i++)
	{
		cell[i]=NULL;
	}

	for(uint i=0; i<num_particle; i++)
	{
		p=&(mem[i]);
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

	int3 cell_pos;
	int3 near_pos;
	uint hash;

	float3 rel_pos;
	float r2;

	for(uint i=0; i<num_particle; i++)
	{
		p=&(mem[i]); 
		cell_pos=calc_cell_pos(p->pos);

		p->dens=0.0f;
		p->pres=0.0f;

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

						np=np->next;
					}
				}
			}
		}

        p->dens=p->dens+p->selfDens;
        //need to change pressure equation based on drift velocity
        //immiscable fluids pressure in equation vanishes
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

	float3 rel_pos;
	float3 rel_vel;

	float r2;
	float r;
	float kernel_r;
	float V;

	float pres_kernel;
	float visc_kernel;
	float temp_force;

	float3 grad_color;
	float lplc_color;

	for(uint i=0; i<num_particle; i++)
	{
		p=&(mem[i]); 
		cell_pos=calc_cell_pos(p->pos);

		p->acc.x=0.0f;
		p->acc.y=0.0f;
		p->acc.z=0.0f;

        p->driftVelocity.x=0.0f;
        p->driftVelocity.y=0.0f;
        p->driftVelocity.z=0.0f;

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

        //DRIFT VELOCITY
        //user constants
        float strengthFac=0.000001f;
        float diffuseConst=0.0001f;

        float pressureRelationship;
        if(misc==1)
        {
            pressureRelationship=1;
        }
        else if(misc==0)
        {
            pressureRelationship=0;
        }

        if(p->phase==1)
        {
            //volume fraction may have to be for individual particle
            float sum=(volume_fraction_1*p->dens/(volume_fraction_1*dens_1))*p->dens;
            float thirdPhase=(volume_fraction_1/volume_fraction_1)-(volume_fraction_1*p->dens/volume_fraction_1*dens_1)*(volume_fraction_1/volume_fraction_1);
            //drift velocity - second phase is where pressure relationship will affect the misc/immiscibility of the mixture
            p->driftVelocity.x=(strengthFac*(dens_1-sum)*p->acc.x)-(strengthFac*pressureRelationship*(dens_1-sum)-diffuseConst*thirdPhase);
            p->driftVelocity.y=(strengthFac*(dens_1-sum)*p->acc.x)-(strengthFac*pressureRelationship*(dens_1-sum)-diffuseConst*thirdPhase);
            p->driftVelocity.z=(strengthFac*(dens_1-sum)*p->acc.y)-(strengthFac*pressureRelationship*(dens_1-sum)-diffuseConst*thirdPhase);
        }
        else if(p->phase==2)
        {
            //volume fraction may have to be for individual particle
            float sum=(volume_fraction_2*p->dens/(volume_fraction_2*dens_2))*p->dens;
            float thirdPhase=(volume_fraction_2/volume_fraction_2)-(volume_fraction_2*p->dens/volume_fraction_2*dens_2)*(volume_fraction_2/volume_fraction_2);
            //drift velocity - second phase is where pressure relationship will affect the misc/immiscibility of the mixture
            p->driftVelocity.x=strengthFac*(dens_2-sum)*p->acc.x-(strengthFac*pressureRelationship*(dens_2-sum)-diffuseConst*thirdPhase);
            p->driftVelocity.y=strengthFac*(dens_2-sum)*p->acc.x-(strengthFac*pressureRelationship*(dens_2-sum)-diffuseConst*thirdPhase);
            p->driftVelocity.z=strengthFac*(dens_2-sum)*p->acc.y-(strengthFac*pressureRelationship*(dens_2-sum)-diffuseConst*thirdPhase);

        }

        //p->driftVelocity.x=0.007;
        //p->driftVelocity.y=0.007;
        //p->driftVelocity.z=0.007;

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

void SPHSystem::advection()
{
	Particle *p;

	for(uint i=0; i<num_particle; i++)
	{
		p=&(mem[i]);

        //potentially add drift velocity in with this step
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

int3 SPHSystem::calc_cell_pos(float3 p)
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
