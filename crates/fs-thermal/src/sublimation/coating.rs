//! A solid deposit and its supporting substrate share one interface temperature.
//!
//! The deposit inventory may start at zero. The substrate retains its own mass
//! and original caloric chart; only the deposit can enter the vapor. The caller
//! must explicitly admit heterogeneous deposition on the chosen substrate.
//! No nucleation barrier, contact-angle law or arbitrary numerical seed is
//! supplied. A zero inventory on an undersaturated surface cannot evaporate.
//!
//! This is an isothermal-through-thickness surface approximation. It is useful
//! for thin deposits with negligible internal resistance, not a resolved frost
//! layer. Temperature gradients, changes of exposed area, porosity, mixed phases
//! and thermal contact resistance need their own spatial model. The prescribed
//! or finite vapor reservoir has the same thermostatted semantics as the
//! original surface operator.

use super::{
    Endpoint, SolidSurfaceState, SolverConfig, SublimationError, SublimationLaw, SurfaceKinetics,
    SurfaceStep, VaporReservoir, finish_surface_step, finite, nonnegative, positive, properties,
    solve_endpoint, transfer, validate_solver,
};

/// Original substrate storage plus the separately conserved solid deposit.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SupportedSurfaceState {
    pub substrate_mass_kg: f64,
    pub substrate_specific_enthalpy_j_kg: f64,
    /// Zero means bare substrate on which the caller admits deposition.
    pub deposit_mass_kg: f64,
    /// Frozen exposed area over this step [m²].
    pub area_m2: f64,
}

/// One accepted phase/storage transaction; no input is mutated on refusal.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SupportedSurfaceStep {
    pub state: SupportedSurfaceState,
    /// Uses the DEPOSIT caloric datum, including its carried vapor enthalpy.
    /// The energy residual here includes both substrate and deposit storage.
    pub phase: SurfaceStep,
    pub substrate_energy_change_j: f64,
    pub deposit_energy_change_j: f64,
}

/// Grow or sublime a deposit against the original substrate/deposit caloric laws.
///
/// `substrate_temperature(h_s)` must describe the substrate's continuous,
/// nondecreasing original chart. `deposit_enthalpy(T)` supplies the original
/// pure-solid deposit enthalpy, with the datum used by its sublimation latent
/// heat. The substrate's datum is independent. Both materials share T, and the
/// deposit's initial enthalpy is evaluated at the original substrate T.
///
/// For signed transfer dm (positive into vapor), conservation is solved as
///
/// M_s (h_s1-h_s0) + m_d0 (h_d1-h_d0) + dm L(T1) = 0.
///
/// The endpoint vapor carries h_d1+L; adding latent heat again would double
/// count energy. Only m_d0 bounds sublimation. The same accepted dm determines
/// reservoir mass, deposit growth/recession and every energy account. The
/// shared bounded solver checks cancellation per evaluation. This is a
/// first-order endpoint approximation; time refinement does not establish the
/// thin-layer or deposition-admission assumptions.
#[allow(clippy::too_many_arguments)]
pub fn advance_supported_surface<L: SublimationLaw + ?Sized>(
    law: &L,
    kinetics: SurfaceKinetics,
    state: SupportedSurfaceState,
    reservoir: VaporReservoir,
    duration_s: f64,
    config: SolverConfig,
    substrate_temperature: impl Fn(f64) -> Result<f64, SublimationError>,
    deposit_enthalpy: impl Fn(f64) -> Result<f64, SublimationError>,
    mut checkpoint: impl FnMut() -> Result<(), SublimationError>,
) -> Result<SupportedSurfaceStep, SublimationError> {
    if !positive(state.substrate_mass_kg)
        || !state.substrate_specific_enthalpy_j_kg.is_finite()
        || !nonnegative(state.deposit_mass_kg)
        || !nonnegative(state.area_m2)
        || !nonnegative(duration_s)
    {
        return Err(SublimationError::Invalid(
            "invalid substrate storage or deposit inventory",
        ));
    }
    validate_solver(config)?;
    checkpoint()?;
    let initial_t = substrate_temperature(state.substrate_specific_enthalpy_j_kg)?;
    properties(law, initial_t)?;
    let initial_deposit_h = finite(deposit_enthalpy(initial_t)?)?;
    let solid = SolidSurfaceState {
        mass_kg: state.deposit_mass_kg,
        specific_enthalpy_j_kg: initial_deposit_h,
        area_m2: state.area_m2,
    };
    let mut evaluations = 0;
    let mut evaluate = |h: f64| -> Result<Endpoint, SublimationError> {
        checkpoint()?;
        if evaluations >= config.max_iterations {
            return Err(SublimationError::NonConverged);
        }
        evaluations += 1;
        finite(h)?;
        let t = substrate_temperature(h)?;
        let p = properties(law, t)?;
        let deposit_h = finite(deposit_enthalpy(t)?)?;
        let dm = transfer(p, t, kinetics, solid, reservoir, duration_s)?;
        let residual = finite(
            state.substrate_mass_kg * (h - state.substrate_specific_enthalpy_j_kg)
                + state.deposit_mass_kg * (deposit_h - initial_deposit_h)
                + dm * p.latent_heat_j_kg,
        )?;
        Ok(Endpoint {
            h,
            t,
            p,
            dm,
            residual,
        })
    };
    let initial = evaluate(state.substrate_specific_enthalpy_j_kg)?;
    let endpoint = solve_endpoint(initial, state.substrate_mass_kg, config, &mut evaluate)?;
    let deposit_h = finite(deposit_enthalpy(endpoint.t)?)?;
    let phase = finish_surface_step(
        solid,
        reservoir,
        Endpoint {
            h: deposit_h,
            ..endpoint
        },
    )?;
    let substrate_energy_change_j =
        finite(state.substrate_mass_kg * (endpoint.h - state.substrate_specific_enthalpy_j_kg))?;
    let deposit_energy_change_j =
        finite(state.deposit_mass_kg * (deposit_h - initial_deposit_h) - endpoint.dm * deposit_h)?;
    checkpoint()?;
    Ok(SupportedSurfaceStep {
        state: SupportedSurfaceState {
            substrate_specific_enthalpy_j_kg: endpoint.h,
            deposit_mass_kg: phase.solid.mass_kg,
            ..state
        },
        phase,
        substrate_energy_change_j,
        deposit_energy_change_j,
    })
}

#[cfg(test)]
mod tests;
