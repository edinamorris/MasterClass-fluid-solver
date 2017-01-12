/// \file SphPhaseData.cpp
/// \brief Ability to change inidvidual values of the different liquids, such as colour and rest density
/// \author Edina Morris

#include "SphPhaseData.h"

//declaring all values for the different phase liquids
phaseData::phaseData()
{
    //All values that can be altered for a different phase
    //higher density liquid, very low in viscocity water like substance
    phases[0].m_colour.m_x=0.117647f;
    phases[0].m_colour.m_y=0.564706f;
    phases[0].m_colour.m_z=1.0f;
    phases[0].m_dens=1300;
    phases[0].m_individualMass=0.02;
    phases[0].m_individualVisc=6.5;

    //lower density liquid, very high in viscoity oil like substance
    phases[1].m_colour.m_x=0.862745f;
    phases[1].m_colour.m_y=0.0784314f;
    phases[1].m_colour.m_z=0.235294f;
    phases[1].m_dens=1000;
    phases[1].m_individualMass=0.03;
    phases[1].m_individualVisc=300.0;

    //Automatic assignment of rest
    for(int phase=0; phase<PHASES; phase++)
    {
        phases[phase].m_numberOfParticles=0;
        phases[phase].m_selfDens=phases[phase].m_individualMass*system.m_poly6Value*pow(system.m_kernel, 6);
        phases[phase].m_selfLplcColor=system.m_lplcPoly6*phases[phase].m_individualMass*system.m_kernel2*(0-3/4*system.m_kernel2);
    }
}
