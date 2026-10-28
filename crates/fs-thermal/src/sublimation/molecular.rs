//! Opposing molecular mass fluxes at a local solid-vapor interface.
//!
//! A net phase-change rate does not specify the energy exchanged when the
//! surface and gas have different temperatures. Retain both half-space fluxes
//! so an energy owner can use the original species caloric laws at each
//! temperature. For equilibrium internal modes in a Maxwellian ideal vapor,
//! the flux-weighted molecular energy is u_v(T) + R_specific*T/2, equivalently
//! h_v(T) - R_specific*T/2. Its datum must be aligned with solid storage by the
//! coupled material owner; no latent heat or energy-reference choice is made
//! by this scalar mass-flux operation. For the flux-weighted translational
//! energy 2RT versus the bulk 3RT/2, see Nathanson (2021), sections 1 and 3:
//! <https://doi.org/10.1007/978-3-030-63963-1_27>. Energy-dependent sticking or
//! non-Maxwellian distributions require a different boundary law.
//!
//! Assumes a local free-molecular gas, equal incoming/outgoing accommodation,
//! and specular reflection of nonsticking molecules. No gas-film resistance,
//! diffuse thermal accommodation, Knudsen-layer drift, or transport from a
//! remote pressure is inferred. Gross fluxes are nonnegative; net flux is
//! positive for sublimation and negative for deposition.

use super::{
    SublimationError, SublimationLaw, SurfaceKinetics, VaporReservoir, coefficients, finite,
    properties,
};

/// Opposing mass fluxes per unit exposed area [kg/(m² s)].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MolecularMassFlux {
    /// Molecules leaving the solid at its own surface temperature.
    pub outgoing_kg_m2_s: f64,
    /// Molecules arriving from the local vapor at its gas temperature.
    pub incoming_kg_m2_s: f64,
    /// Signed outgoing minus incoming mass, with stable equal-T evaluation.
    pub net_kg_m2_s: f64,
}

/// Local Hertz-Knudsen molecular exchange with one common accommodation.
///
/// Uses the same bounded constitutive properties and kinetic coefficients as
/// the existing signed surface-flux API. Gas pressure means volatile PARTIAL
/// pressure, not carrier-gas total pressure. A coefficient of zero seals the
/// interface. Nonfinite/out-of-domain inputs and unrepresentable gross fluxes
/// refuse; no material or gas state is changed. This constant-work scalar
/// evaluation requires no cancellation boundary of its own.
pub fn molecular_mass_flux<L: SublimationLaw + ?Sized>(
    law: &L,
    surface_temperature_k: f64,
    accommodation_coefficient: f64,
    reservoir: VaporReservoir,
) -> Result<MolecularMassFlux, SublimationError> {
    let p = properties(law, surface_temperature_k)?;
    let vapor_pressure = reservoir.partial_pressure_pa(p.molar_mass_kg_mol)?;
    let gas_temperature = reservoir.temperature_k();
    let (outgoing_coefficient, incoming_coefficient) = coefficients(
        p,
        surface_temperature_k,
        gas_temperature,
        SurfaceKinetics {
            accommodation_coefficient,
            mass_transfer_coefficient_m_s: None,
        },
    )?;
    let outgoing = finite(outgoing_coefficient * p.saturation_pressure_pa)?;
    let incoming = finite(incoming_coefficient * vapor_pressure)?;
    let net = if surface_temperature_k == gas_temperature {
        outgoing_coefficient * (p.saturation_pressure_pa - vapor_pressure)
    } else {
        outgoing - incoming
    };
    Ok(MolecularMassFlux {
        outgoing_kg_m2_s: outgoing,
        incoming_kg_m2_s: incoming,
        net_kg_m2_s: finite(net)?,
    })
}

#[cfg(test)]
mod tests;
