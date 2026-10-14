#ifndef ERF_NULL_FIRE_H
#define ERF_NULL_FIRE_H

#include <ERF_FireModelBase.H>

/**
 * @file ERF_NullFire.H
 * @brief The fire model of a run without a fire (erf.fire_model = None), and of
 *        every level that does not carry the fire grid
 */
class NullFire : public FireModelBase
{
public:
    const char* Name () const override { return "None"; }
    void Init (const FireSetup& /*setup*/) override {}
    void Advance (const FireForcing& /*forcing*/) override {}
    void WritePlotfile (const std::string& /*prefix*/, amrex::Real /*time*/, int /*step*/) const override {}
};

#endif
