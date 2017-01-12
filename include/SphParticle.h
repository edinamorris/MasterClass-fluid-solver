/// \file SphParticle.h
/// \brief all the data for a particle and it's individual attributes
/// \author Dongli Zhang modified by Edina Morris

#ifndef SPHPARTICLE_H_
#define SPHPARTICLE_H_

#include <stdio.h>
#include <stdlib.h>
#include "SphHeader.h"
#include "SphType.h"
#include "vector"
#include "Vec3.h"
#include "SphType.h"

//used by all particles in the system
class Particle
{
public:
    Particle();

    uint m_id;
    vec3 m_pos;
    vec3 m_vel;
    vec3 m_acc;
    //estimate velocity
    vec3 m_ev;
    vec3 m_colour;
    //to determine each particle's original phase
    int m_phase;
    //used to define the miscibility property, indiivdual to each particle so different phases can be miscible/immiscible with each other
    int m_misc;
    //each particle must have a drift velocity and volume fraction calculated for each phase
    vec3 m_driftVelocity [PHASES];
    float m_volumeFraction [PHASES];

    //each have individual values set
    float m_visc;
    float m_mass;
    float m_selfDens;
    float m_dens;
    float m_restdens;
    float m_pres;

    float m_surfNorm;
    float m_lplcColour;

    //neighbour vector
    std::vector<Particle*> neighbour;

    Particle *next;
};


#endif  //SPHPARTICLE_H_
