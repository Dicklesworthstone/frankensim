//! Independent implementation of NIST SRD 107's generalized BEB Eq. (1).
//!
//! Source data content version: August 2004; parameter facts retrieved
//! 2026-10-31. This is a theoretical single-ionization model, not measured
//! collision data or the BED differential-energy model. The original Q values
//! are retained; replacing Q by one would be a further approximation.
//!
//! Equation and source constants:
//! <https://physics.nist.gov/PhysRefData/Ionization/Eqs/latex.html>,
//! <https://physics.nist.gov/PhysRefData/Ionization/intro.html>.
//! Orbital parameters:
//! <https://physics.nist.gov/cgi-bin/Ionization/ion_data.php?differential=Y&id=HI&initial=&ision=I>,
//! <https://physics.nist.gov/cgi-bin/Ionization/bun.php?Id=HeI>.
//!
//! The source forbids this nonrelativistic model above 10 keV. Direct source
//! evaluation stops there. Its mathematical BEB continuation has the declared
//! global upper envelope used to assess the omitted Maxwellian tail. That
//! envelope bounds this theoretical formula, not an unknown relativistic real
//! cross section; a resolved tail is not a physical error certificate.
//!
//! NIST SRD compilations have separate copyright/usage terms; no public-domain
//! or redistribution license is inferred from federal authorship. Only the
//! cited orbital parameter facts and independently implemented equation are
//! retained here, not the database compilation or experimental datasets.

use super::*;

// Retain the constants in NIST's published model, not a silently revised fit.
const BOHR_RADIUS_M: f64 = 0.529_18e-10;
const RYDBERG_EV: f64 = 13.6057;
const TERMS: &str = "NIST Standard Reference Database 107 compilation terms: https://www.nist.gov/open/copyright-fair-use-and-licensing-statements-srd-data-software-and-technical-series-publications ; no additional license grant asserted";

/// Ground-state neutral H single ionization, with NIST's original B=13.6057 eV.
/// This is intentionally not replaced by a newer spectroscopic ionization
/// energy; chemical calorics retain their own independently authored reference.
pub fn hydrogen(target: SpeciesId) -> ElectronCrossSectionV1 {
    model(target, 13.6057, 13.6057, 1.0, 0.5668,
        "NIST SRD 107 v3 (August 2004), generalized BEB Eq.1, neutral H 1s: https://physics.nist.gov/cgi-bin/Ionization/ion_data.php?differential=Y&id=HI&initial=&ision=I ; DOI 10.18434/T4KK5C")
}

/// Ground-state neutral He single ionization, original NIST 1s parameters.
pub fn helium(target: SpeciesId) -> ElectronCrossSectionV1 {
    model(target, 24.5870, 39.5100, 2.0, 0.8841,
        "NIST SRD 107 v3 (August 2004), generalized BEB Eq.1, neutral He 1s: https://physics.nist.gov/cgi-bin/Ionization/bun.php?Id=HeI ; DOI 10.18434/T4KK5C")
}

fn model(target: SpeciesId, b_ev: f64, u_ev: f64, occupation: f64, q: f64, citation: &str) -> ElectronCrossSectionV1 {
    let s = 4.0 * core::f64::consts::PI * BOHR_RADIUS_M * BOHR_RADIUS_M
        * occupation * (RYDBERG_EV / b_ev) * (RYDBERG_EV / b_ev);
    // For t>=1, ln(t)/t<=1/e, 0<=1-1/t-ln(t)/(t+1)<=1, and denominator>=t.
    let envelope = s * (0.5 * q / core::f64::consts::E + 2.0 - q);
    ElectronCrossSectionV1 {
        target, kind: ElectronCollisionKindV1::Ionization,
        provenance: Provenance { source: citation.into(), license: TERMS.into(), artifact: None },
        threshold: Energy::new(b_ev * ELECTRONVOLT_J),
        maximum_energy: Energy::new(10_000.0 * ELECTRONVOLT_J),
        high_energy_upper_bound: Area::new(envelope), scale: envelope,
        law: CrossSectionLaw::Beb { binding_energy: b_ev * ELECTRONVOLT_J,
            orbital_kinetic_energy: u_ev * ELECTRONVOLT_J, occupation, dipole_constant: q },
    }
}

pub(super) fn beb_value(energy: f64, binding: f64, orbital: f64, occupation: f64, q: f64) -> f64 {
    if energy <= binding { return 0.0; }
    let x = (energy - binding) / binding;
    let t = 1.0 + x;
    let log_t = crate::entropy::log_ratio(energy, binding);
    let inverse = 1.0 / t;
    let b_ev = binding / ELECTRONVOLT_J;
    let scale = 4.0 * core::f64::consts::PI * BOHR_RADIUS_M * BOHR_RADIUS_M
        * occupation * (RYDBERG_EV / b_ev) * (RYDBERG_EV / b_ev);
    // Factored threshold differences preserve positive near-threshold values.
    let bracket = 0.5 * q * log_t * (x * inverse) * (1.0 + inverse)
        + (2.0 - q) * (x * inverse - log_t / (t + 1.0));
    scale * bracket / (t + orbital / binding + 1.0)
}
