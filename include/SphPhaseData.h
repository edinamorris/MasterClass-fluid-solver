/// \file SphPhasedata.h
/// \brief Methods for assigning individual values for different phases (liquids) such as rest density, viscocity and mass
/// \author Edina Morris

#ifndef SPH_PHASEDATA_H_
#define SPH_PHASEDATA_H_

#include "Vec3.h"
#include "SphSystem.h"

//each phase has seperate default values such as density, mass and viscosity
struct Phase
{
    vec3 m_colour;
    float m_individualMass;
    float m_individualVisc;
    float m_selfDens;
    float m_selfLplcColor;
    float m_dens;
    int m_numberOfParticles;
};

//used in drift velocity and volume fraction calculations to return a phases' preliminary data such as rest density or colour
class phaseData
{
public:
    phaseData();

    SPHSystem system;
    Phase *phases = new Phase[system.getNumPhases()];

    int getNumberOfParticles(const int _phase) const {return phases[_phase].m_numberOfParticles;}
    float getDensity(const int _phase) const {return phases[_phase].m_dens;}
    float getMass(const int _phase) const {return phases[_phase].m_individualMass;}
    vec3 getColour(const int _phase) const {return phases[_phase].m_colour;}

    void addParticle(const int _phase){phases[_phase].m_numberOfParticles+=1;}

private:

};

#endif // SPH_PHASEDATA_H_

