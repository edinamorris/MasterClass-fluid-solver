/** File:		sph_system.h
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

#ifndef __SPHSYSTEM_H__
#define __SPHSYSTEM_H__

#include <stdio.h>
#include <stdlib.h>
#include "sph_type.h"
#include "vector"
#include "vec3.h"

class Particle
{
public:
	uint id;
    vec3 pos;
    vec3 vel;
    vec3 acc;
    vec3 ev;
    vec3 colour;
    int phase;
    vec3 driftVelocity;

    //made these individual for each particle
    float visc;
    float mass;
    float selfDens;
    float lplcColour;
    float volFrac;

	float dens;
    //individual phase density for liquids
    float restdens;
	float pres;

	float surf_norm;

    //neighbour vector
    std::vector<Particle*> neighbour;

	Particle *next;
};

class SPHSystem
{
public:
	uint max_particle;
	uint num_particle;
    uint phase1Particle;
    uint phase2Particle;

	float kernel;

    vec3 world_size;
	float cell_size;
	uint3 grid_size;
	uint tot_cell;

    vec3 gravity;
	float wall_damping;
	float gas_constant;
	float time_step;
	float surf_norm;
	float surf_coe;

	float poly6_value;
	float spiky_value;
	float visco_value;

	float grad_poly6;
    float grad_spiky;
	float lplc_poly6;

	float kernel_2;

    int misc;

    //phase 1
    vec3 colour_1;
    float individualMass_1;
    float individualVisc_1;
    float self_dens_1;
    float self_lplc_color_1;
    float dens_1;
    float volume_fraction_1;

    //phase 2
    vec3 colour_2;
    float individualMass_2;
    float individualVisc_2;
    float self_dens_2;
    float self_lplc_color_2;
    float dens_2;
    float volume_fraction_2;

	Particle *mem;
	Particle **cell;

	uint sys_running;

public:
	SPHSystem();
	~SPHSystem();
	void animation();
	void init_system();
    void damnScenario();
    void dropScenario();
    void add_particle(int _phase, vec3 pos, vec3 vel);
    //returns gradient
    vec3 poly6Grad(vec3 _iPos, vec3 _jPos);
    vec3 spikyGrad(vec3 _iPos, vec3 _jPos);

private:
	void build_table();
	void comp_dens_pres();
	void comp_force_adv();
    void driftVelocity();
	void advection();

private:
    int3 calc_cell_pos(vec3 p);
	uint calc_cell_hash(int3 cell_pos);
};

#endif
