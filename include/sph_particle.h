#ifndef __SPHPARTICLE_H__
#define __SPHPARTICLE_H__

#include <stdio.h>
#include <stdlib.h>
#include "sph_header.h"
#include "sph_type.h"
#include "vector"
#include "vec3.h"
#include "sph_type.h"

class Particle
{
public:
    Particle();

    uint id;
    vec3 pos;
    vec3 vel;
    vec3 acc;
    vec3 ev;
    vec3 colour;
    int phase;
    int misc;
    //need to be changed so not hard-coded
    //each particle must have a drift velocity calculated for each phase, so when theres two phases there will be two drift velocities
    vec3 driftVelocity [PHASES];
    float volumeFraction [PHASES];
    float prevVolumeFraction [PHASES];

    //each have individual values set
    float visc;
    float mass;
    float selfDens;
    float dens;
    float restdens;
    float pres;

    float surf_norm;
    float lplcColour;

    //neighbour vector
    std::vector<Particle*> neighbour;

    Particle *next;
};


#endif
