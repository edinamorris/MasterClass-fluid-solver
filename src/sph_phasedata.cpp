#include "sph_phasedata.h"

//declaring all values for the different phase liquids
phaseData::phaseData()
{
    //higher density liquid, very low in viscocity water like substance
    //first liquid
    phases[0].colour.x=0.117647f;
    phases[0].colour.y=0.564706f;
    phases[0].colour.z=1.0f;
    phases[0].dens=1300;
    phases[0].individualMass=0.02;
    phases[0].individualVisc=6.5;

    phases[0].numberOfParticles=0;

    phases[0].selfDens=phases[0].individualMass*system.poly6_value*pow(system.kernel, 6);
    phases[0].self_lplc_color=system.lplc_poly6*phases[0].individualMass*system.kernel_2*(0-3/4*system.kernel_2);

    //lower density liquid, very high in viscoity oil like substance
    //second liquid
    phases[1].colour.x=0.862745f;
    phases[1].colour.y=0.0784314f;
    phases[1].colour.z=0.235294f;
    phases[1].dens=1000;
    phases[1].individualMass=0.03;
    phases[1].individualVisc=300.0;

    phases[1].numberOfParticles=0;

    phases[1].selfDens=phases[1].individualMass*system.poly6_value*pow(system.kernel, 6);
    phases[1].self_lplc_color=system.lplc_poly6*phases[1].individualMass*system.kernel_2*(0-3/4*system.kernel_2);
}
