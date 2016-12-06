#ifndef __SPHPARTICLE_H__
#define __SPHPARTICLE_H__

#include <stdio.h>
#include <stdlib.h>
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
    //need to be changed so not hard-coded
    //each particle must have a drift velocity calculated for each phase, so when theres two phases there will be two drift velocities
    vec3 driftVelocity [2];
    float volumeFraction [2];

    //made these individual for each particle
    float visc;
    float mass;
    float selfDens;
    float lplcColour;

    float dens;
    //individual phase density for liquids
    float restdens;
    float pres;

    float surf_norm;

    //neighbour vector
    std::vector<Particle*> neighbour;

    Particle *next;
};


#endif
