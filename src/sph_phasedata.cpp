#include "sph_phasedata.h"

//declaring all values for the different phase liquids
phaseData::phaseData()
{
    //first liquid
    phases[0].colour.x=0.2f;
    phases[0].colour.y=0.8f;
    phases[0].colour.z=1.0f;
    phases[0].dens=1300;
    phases[0].individualMass=0.02;
    phases[0].individualVisc=6.5;
    //volume fraction of phase = number of elements in phase/total particle number = both of these need to be updated after all particles have been added
    phases[0].volumeFraction=0.0f;
    phases[0].numberOfParticles=0;

    phases[0].selfDens=phases[0].individualMass*system.poly6_value*pow(system.kernel, 6);
    phases[0].self_lplc_color=system.lplc_poly6*phases[0].individualMass*system.kernel_2*(0-3/4*system.kernel_2);

    //second liquid
    phases[1].colour.x=0.294118f;
    phases[1].colour.y=0.0f;
    phases[1].colour.z=0.509804f;
    phases[1].dens=1000;
    phases[1].individualMass=0.03;
    phases[1].individualVisc=350.0;
    //volume fraction of phase = number of elements in phase/total particle number = both of these need to be updated after all particles have been added
    phases[1].volumeFraction=0.0f;
    phases[1].numberOfParticles=0;

    phases[1].selfDens=phases[1].individualMass*system.poly6_value*pow(system.kernel, 6);
    phases[1].self_lplc_color=system.lplc_poly6*phases[1].individualMass*system.kernel_2*(0-3/4*system.kernel_2);
}
