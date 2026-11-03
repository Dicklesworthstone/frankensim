//! Normal momentum carried by the original molecular mass streams.
//!
//! Gas and surface are stationary in the SAME frame. Incident molecules are
//! Maxwellian with no bulk drift; emitted molecules have a cosine-weighted
//! half-Maxwellian distribution. Their flux-weighted normal speed is
//! `sqrt(pi*R_specific*T/2)`, not the bulk mean or rms molecular speed. See
//! Nathanson (2021), section 1, <https://doi.org/10.1007/978-3-030-63963-1_27>.
//!
//! Normal points from solid into gas. Emission and removal of arriving negative
//! momentum BOTH give positive momentum to gas. The surface receives the
//! opposite load even at zero net phase mass. Reflection is a separate account
//! and must not be added again as ambient pressure on the same molecules.
//! These allocation-free scalar moments neither advance momentum nor add heat:
//! the existing molecular stream energy already includes translational energy.

use super::MolecularMassFlux;
use crate::sublimation::{SublimationError, nonnegative, positive};
use fs_math::det::sqrt;

fn represented(value: f64, nonzero: bool) -> Result<f64, SublimationError> {
    if !nonnegative(value) || (nonzero && value == 0.0) {
        Err(SublimationError::Invalid(
            "molecular normal momentum is unrepresentable",
        ))
    } else {
        Ok(value)
    }
}

/// Mean normal speed [m/s] of molecules crossing one Maxwellian half-space.
///
/// Supply the SAME specific gas constant as the original molecular/caloric
/// owner. Factoring the square roots avoids overflowing or underflowing `R*T`.
/// Finite positive inputs are required; a nonrepresentable speed refuses.
pub fn flux_weighted_normal_speed(
    specific_gas_constant_j_kg_k: f64,
    temperature_k: f64,
) -> Result<f64, SublimationError> {
    if !positive(specific_gas_constant_j_kg_k) || !positive(temperature_k) {
        return Err(SublimationError::Invalid(
            "positive finite molecular R and temperature required",
        ));
    }
    represented(
        (sqrt(specific_gas_constant_j_kg_k) * sqrt(temperature_k))
            * sqrt(core::f64::consts::FRAC_PI_2),
        true,
    )
}

/// Opposing phase beams' normal momentum rates per exposed area [Pa].
/// Both are nonnegative in the outward-solid-normal momentum balance.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MolecularPhaseMomentumFlux {
    pub outgoing_pa: f64,
    pub incoming_pa: f64,
    /// Sum of emission and sticking; excludes every reflected molecule.
    pub total_pa: f64,
}

/// Momentum of the original accepted gross phase fluxes, without resampling.
///
/// Inputs must refer to the same original species and temperature snapshot.
/// `net_kg_m2_s` is not multiplied by an assumed exhaust speed. Empty beams
/// give exact zero after validating inputs, without evaluating an unused speed.
pub fn phase_normal_momentum(
    flux: MolecularMassFlux,
    specific_gas_constant_j_kg_k: f64,
    surface_temperature_k: f64,
    gas_temperature_k: f64,
) -> Result<MolecularPhaseMomentumFlux, SublimationError> {
    if !positive(specific_gas_constant_j_kg_k)
        || !positive(surface_temperature_k)
        || !positive(gas_temperature_k)
        || !nonnegative(flux.outgoing_kg_m2_s)
        || !nonnegative(flux.incoming_kg_m2_s)
        || !flux.net_kg_m2_s.is_finite()
    {
        return Err(SublimationError::Invalid(
            "finite original molecular beams and positive R/T required",
        ));
    }
    let beam = |mass: f64, temperature: f64| {
        if mass == 0.0 {
            return Ok(0.0);
        }
        represented(
            mass * flux_weighted_normal_speed(specific_gas_constant_j_kg_k, temperature)?,
            true,
        )
    };
    let outgoing_pa = beam(flux.outgoing_kg_m2_s, surface_temperature_k)?;
    let incoming_pa = beam(flux.incoming_kg_m2_s, gas_temperature_k)?;
    let total_pa = represented(outgoing_pa + incoming_pa, false)?;
    Ok(MolecularPhaseMomentumFlux {
        outgoing_pa,
        incoming_pa,
        total_pa,
    })
}

/// Complete incident-plus-reflected load of NONSTICKING molecules [Pa].
///
/// `diffuse_reflection_fraction` explicitly declares a Maxwell scattering
/// model: this fraction is diffusely reemitted at the surface temperature and
/// the remainder reflects specularly. It is NOT inferred from a thermal-energy
/// accommodation coefficient. Fully diffuse reemission retains incident mass;
/// it is not another saturation-pressure emission stream. Internal modes are
/// equilibrated when this closure is coupled to the original reflected heat.
///
/// Carrier species use sticking=0. Specular reflection uses diffuse=0 and
/// returns `(1-sticking)*partial_pressure`. The general result is
/// `(1-sticking)*p*(1-diffuse/2 + diffuse/2*sqrt(T_surface/T_gas))`.
/// Zero pressure or sticking=1 gives exact zero. All inputs are still checked.
/// Unrepresentable nonzero intermediate moments refuse rather than disappear.
pub fn maxwell_reflected_normal_pressure(
    partial_pressure_pa: f64,
    surface_temperature_k: f64,
    gas_temperature_k: f64,
    sticking_probability: f64,
    diffuse_reflection_fraction: f64,
) -> Result<f64, SublimationError> {
    if !nonnegative(partial_pressure_pa)
        || !positive(surface_temperature_k)
        || !positive(gas_temperature_k)
        || !sticking_probability.is_finite()
        || !(0.0..=1.0).contains(&sticking_probability)
        || !diffuse_reflection_fraction.is_finite()
        || !(0.0..=1.0).contains(&diffuse_reflection_fraction)
    {
        return Err(SublimationError::Invalid(
            "finite partial pressure, positive temperatures and probabilities in [0,1] required",
        ));
    }
    if partial_pressure_pa == 0.0 || sticking_probability == 1.0 {
        return Ok(0.0);
    }
    let pressure = represented(partial_pressure_pa * (1.0 - sticking_probability), true)?;
    if diffuse_reflection_fraction == 0.0 || surface_temperature_k == gas_temperature_k {
        return Ok(pressure);
    }
    let fraction = 0.5 * diffuse_reflection_fraction;
    let diffuse = represented(pressure * fraction, true)?;
    let surface_root = sqrt(surface_temperature_k);
    let gas_root = sqrt(gas_temperature_k);
    let ratio = surface_root / gas_root;
    // A finite pressure can scale an overflowing temperature ratio back into
    // range. In this branch division by the tiny gas root precedes multiplication.
    let outgoing = if ratio.is_finite() {
        diffuse * ratio
    } else {
        (diffuse / gas_root) * surface_root
    };
    let outgoing = represented(outgoing, true)?;
    represented(pressure * (1.0 - fraction) + outgoing, true)
}

#[cfg(test)]
mod tests;
