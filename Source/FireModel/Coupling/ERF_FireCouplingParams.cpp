#include <ERF_FireCouplingParams.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_String.H>

#include <string>
#include <vector>

FireCouplingParams
FireCouplingParams::read ()
{
    FireCouplingParams c;
    amrex::ParmParse pp("erf.fire");
    // Every message names the key and the value
    auto fail = [](const std::string& msg) { amrex::Abort("[FIRE] " + msg); };
    auto one_of = [&](const char* key, const std::string& v, const std::vector<std::string>& allowed) {
        for (const auto& a : allowed) { if (v == a) { return; } }
        std::string list;
        for (const auto& a : allowed) { list += (list.empty() ? "" : ", ") + std::string("\"") + a + "\""; }
        fail(std::string("erf.fire.") + key + " = \"" + v + "\" is not one of " + list);
    };
    auto positive = [&](const char* key, amrex::Real v) {
        if (!(v > 0.0)) { fail(std::string("erf.fire.") + key + " must be > 0, got " + std::to_string(v)); }
    };
    auto nonneg = [&](const char* key, amrex::Real v) {
        if (!(v >= 0.0)) { fail(std::string("erf.fire.") + key + " must be >= 0, got " + std::to_string(v)); }
    };

    // erf.fire.enable and erf.fire.engine selected the fire before erf.fire_model did
    for (const char* k : {"enable", "engine"}) {
        if (pp.contains(k)) {
            fail(std::string("erf.fire.") + k + " is no longer read: select the fire model with "
                 "erf.fire_model (None, or a model built into this ERF)");
        }
    }

    pp.query("fire_debug", c.fire_debug);
    pp.query("grid_ratio", c.grid_ratio);
    if (c.grid_ratio < 1) {
        fail("erf.fire.grid_ratio must be >= 1, got " + std::to_string(c.grid_ratio));
    }
    pp.query("anchor_level", c.anchor_level);
    if (c.anchor_level < -1) {
        fail("erf.fire.anchor_level must be -1 (the finest level) or a level >= 0, got "
             + std::to_string(c.anchor_level));
    }

    pp.query("heat_flux_alfg", c.heat_flux_alfg);
    pp.query("fire_atm_feedback", c.fire_atm_feedback);
    pp.query("coupling_type", c.coupling_type);
    one_of("coupling_type", c.coupling_type, {"passive", "lagged", "synchronous"});
    positive("heat_flux_alfg", c.heat_flux_alfg);
    nonneg("fire_atm_feedback", c.fire_atm_feedback);   // 0 = one-way, the documented switch
    pp.query("inject_latent", c.inject_latent);
    c.source_mode_given = pp.query("source_mode", c.source_mode);
    one_of("source_mode", c.source_mode, {"overwrite", "add"});
    pp.query("heat_open_fraction", c.heat_open_fraction);
    pp.query("heat_tendency_density", c.heat_tendency_density);
    pp.query("heat_tendency_exner", c.heat_tendency_exner);
    // The keys of the heat injection, given with passive coupling, which
    // injects no heat: warn once per key
    if (c.is_passive()) {
        for (const char* k : {"heat_flux_alfg", "fire_atm_feedback", "inject_latent", "source_mode",
                              "heat_open_fraction", "heat_tendency_density", "heat_tendency_exner"}) {
            if (pp.contains(k)) {
                amrex::Print() << "[FIRE] WARNING: erf.fire." << k << " is ignored with erf.fire.coupling_type"
                               << " = passive (no heat is injected)\n";
            }
        }
    }
    pp.query("smoke_enable", c.smoke_enable);
    pp.query("smoke_emission_factor", c.smoke_emission_factor);
    nonneg("smoke_emission_factor", c.smoke_emission_factor);
    // A smoke key without the smoke switch: a deck that sets the smoke switch
    // to false, or erf.fire_model to None, explicitly is a variant switching
    // smoke off, so warn; a deck that set neither forgot the smoke switch, so stop
    if (!c.smoke_enable && pp.contains("smoke_emission_factor")) {
        std::string fire_model;
        const bool fire_off = amrex::ParmParse("erf").query("fire_model", fire_model)
                              && amrex::toLower(fire_model) == "none";
        if (pp.contains("smoke_enable")) {
            amrex::Print() << "[FIRE] WARNING: erf.fire.smoke_emission_factor is ignored with"
                           << " erf.fire.smoke_enable = false\n";
        } else if (fire_off) {
            amrex::Print() << "[FIRE] WARNING: erf.fire.smoke_emission_factor is ignored with"
                           << " erf.fire_model = None\n";
        } else {
            fail("erf.fire.smoke_emission_factor needs erf.fire.smoke_enable = true");
        }
    }
    return c;
}
