#ifndef SPH_PHASEDATA
#define SPH_PHASEDATA

#include "vec3.h"
#include "sph_system.h"

struct Phase
{
    vec3 colour;
    float individualMass;
    float individualVisc;
    float selfDens;
    float self_lplc_color;
    float dens;
    int numberOfParticles;
};

class phaseData
{
public:
    phaseData();

    SPHSystem system;
    Phase *phases = new Phase[system.getNumPhases()];

    int getNumberOfParticles(int _phase){return phases[_phase].numberOfParticles;}
    float getDensity(int _phase){return phases[_phase].dens;}
    float getMass(int _phase){return phases[_phase].individualMass;}
    vec3 getColour(int _phase){return phases[_phase].colour;}

    void addParticle(int _phase){phases[_phase].numberOfParticles+=1;}

private:

};

#endif // SPH_PHASEDATA

