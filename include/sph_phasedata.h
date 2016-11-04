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
    float volumeFraction;
    int numberOfParticles; //in this phase
};

class phaseData
{
public:
    phaseData();

    SPHSystem system;
    Phase *phases = new Phase[system.getNumPhases()];

    int getNumberOfParticles(int _phase){return phases[_phase].numberOfParticles;}
    float getVolumeFraction(int _phase){return phases[_phase].volumeFraction;}
    float getDensity(int _phase){return phases[_phase].dens;}
    float getMass(int _phase){return phases[_phase].individualMass;}
    void addParticle(int _phase){phases[_phase].numberOfParticles+=1;}

    //need to have these updated after all particles have been added
    void setVolumeFraction(int _volFrac, int _phase){phases[_phase].volumeFraction=_volFrac;}
    void setNumberOfParticles(int _numPart, int _phase){phases[_phase].numberOfParticles=_numPart;}
private:

};

#endif // SPH_PHASEDATA

