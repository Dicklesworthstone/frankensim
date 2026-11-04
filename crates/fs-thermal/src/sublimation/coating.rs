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

/// Accepted heat storage at fixed substrate/deposit masses and exposed area.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SupportedHeatStep {
    /// Original substrate caloric coordinate at the common endpoint temperature.
    pub state: SupportedSurfaceState,
    /// Common substrate and deposit temperature [K].
    pub surface_temperature_k: f64,
    /// Signed substrate storage change [J] in its own caloric chart.
    pub substrate_energy_change_j: f64,
    /// Signed deposit storage change [J] in its independent caloric chart.
    pub deposit_energy_change_j: f64,
    /// Actual combined storage change minus the requested heat [J].
    pub energy_residual_j: f64,
}

/// Apply signed heat to a substrate and deposit at their shared temperature.
///
/// Positive `heat_j` heats the surface and negative heat cools it. The original
/// continuous, nondecreasing charts `substrate_temperature(h_s)` and
/// `deposit_enthalpy(T)` determine both heat capacities; the deposit is never
/// treated as thermally massless. With masses held fixed, the endpoint solves
///
/// M_s (h_s1-h_s0) + m_d (h_d(T1)-h_d(T0)) = heat_j.
///
/// Each material keeps its independent energy reference. Differences are taken
/// before mass weighting, rather than subtracting absolute stored energies.
/// The same solid phase law admits every trial temperature, including on a bare
/// substrate: heating outside that law's range requires a separate phase owner.
/// No phase transfer, vapor energy, latent heat or geometry change is added.
///
/// Invalid outer trials are backtracked by the shared bounded scalar solver.
/// A heat increment that cannot change the represented enthalpy is accepted
/// only if its reported residual is within the absolute/relative tolerance;
/// otherwise the operation refuses. Cancellation or failure publishes no state.
pub fn heat_supported_surface<L: SublimationLaw + ?Sized>(
    law: &L,
    state: SupportedSurfaceState,
    heat_j: f64,
    config: SolverConfig,
    substrate_temperature: impl Fn(f64) -> Result<f64, SublimationError>,
    deposit_enthalpy: impl Fn(f64) -> Result<f64, SublimationError>,
    mut checkpoint: impl FnMut() -> Result<(), SublimationError>,
) -> Result<SupportedHeatStep, SublimationError> {
    if !positive(state.substrate_mass_kg)
        || !state.substrate_specific_enthalpy_j_kg.is_finite()
        || !nonnegative(state.deposit_mass_kg)
        || !nonnegative(state.area_m2)
        || !heat_j.is_finite()
    {
        return Err(SublimationError::Invalid(
            "invalid supported heat, substrate storage or deposit inventory",
        ));
    }
    validate_solver(config)?;
    checkpoint()?;
    let initial_t = substrate_temperature(state.substrate_specific_enthalpy_j_kg)?;
    properties(law, initial_t)?;
    let initial_deposit_h = finite(deposit_enthalpy(initial_t)?)?;
    let storage = |h: f64, t: f64| -> Result<(f64, f64), SublimationError> {
        let substrate =
            finite(state.substrate_mass_kg * (h - state.substrate_specific_enthalpy_j_kg))?;
        let deposit_h = finite(deposit_enthalpy(t)?)?;
        let deposit = if state.deposit_mass_kg == 0.0 {
            0.0
        } else {
            finite(state.deposit_mass_kg * (deposit_h - initial_deposit_h))?
        };
        Ok((substrate, deposit))
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
        let (substrate, deposit) = storage(h, t)?;
        Ok(Endpoint {
            h,
            t,
            p,
            dm: 0.0,
            residual: finite(substrate + deposit - heat_j)?,
        })
    };
    let initial = evaluate(state.substrate_specific_enthalpy_j_kg)?;
    let endpoint = solve_endpoint(initial, state.substrate_mass_kg, config, &mut evaluate)?;
    let (substrate_energy_change_j, deposit_energy_change_j) = storage(endpoint.h, endpoint.t)?;
    let energy_residual_j = finite(substrate_energy_change_j + deposit_energy_change_j - heat_j)?;
    let tolerance =
        finite(config.absolute_energy_tolerance_j + config.relative_tolerance * heat_j.abs())?;
    if energy_residual_j.abs() > tolerance {
        return Err(SublimationError::NonConverged);
    }
    checkpoint()?;
    Ok(SupportedHeatStep {
        state: SupportedSurfaceState {
            substrate_specific_enthalpy_j_kg: endpoint.h,
            ..state
        },
        surface_temperature_k: endpoint.t,
        substrate_energy_change_j,
        deposit_energy_change_j,
        energy_residual_j,
    })
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
