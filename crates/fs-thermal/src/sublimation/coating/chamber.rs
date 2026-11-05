//! Supported deposits sharing one finite, thermostatted vapor chamber.
//!
//! Bare substrates can gain material without a numerical seed, and exhausted
//! deposits cannot spend substrate mass. Every local original caloric solve
//! sees the same endpoint partial pressure. The existing chamber pressure
//! closure therefore spends the vapor inventory once and permits a warmer
//! coating to supply a colder, initially bare substrate in the same interval.
//!
//! This retains the thin, isothermal coating assumptions of the parent module.
//! The caller explicitly admits heterogeneous deposition, one volatile species
//! and a common deposit/vapor enthalpy datum. Substrate energy references remain
//! independent. Gas sensible heat, nucleation, porosity, moving geometry and
//! through-coating thermal resistance are not inferred by this operation.

use super::{SupportedSurfaceState, SupportedSurfaceStep, advance_supported_surface};
use crate::sublimation::chamber::{
    ChamberExchange, ChamberSolverConfig, InitialSurface, compensated_add, solve_chamber,
};
use crate::sublimation::{
    SublimationError, SublimationLaw, SurfaceKinetics, SurfaceStep, VaporReservoir, finite,
    nonnegative, positive, properties,
};

/// One persistent substrate and its separately conserved, possibly zero deposit.
/// All phase charts describe the same volatile with the same fixed molar mass;
/// all deposit calorics use a common datum compatible with their latent heat.
#[derive(Clone, Copy)]
pub struct SupportedChamberSurface<'a> {
    pub state: SupportedSurfaceState,
    pub kinetics: SurfaceKinetics,
    pub law: &'a dyn SublimationLaw,
    pub deposit_enthalpy: &'a dyn Fn(f64) -> Result<f64, SublimationError>,
}

/// One immutable joint candidate. Surface order matches the input exactly.
#[derive(Clone, Debug, PartialEq)]
pub struct SupportedChamberStep {
    /// Each phase result carries this same final reservoir, not a copied bath.
    pub surfaces: Vec<SupportedSurfaceStep>,
    pub reservoir: VaporReservoir,
    /// Positive from deposits into vapor; negative during deposition, kg.
    pub transferred_mass_kg: f64,
    /// Signed latent energy, already included in `vapor_energy_j`, J.
    pub latent_energy_j: f64,
    /// Signed sensible-plus-latent energy delivered to the vapor thermostat, J.
    pub vapor_energy_j: f64,
    pub substrate_energy_change_j: f64,
    pub deposit_energy_change_j: f64,
    /// Actual final vapor mass minus initial mass minus accepted transfer, kg.
    pub mass_residual_kg: f64,
    /// Trial-pressure mass closure, separate from represented inventory, kg.
    pub nonlinear_mass_residual_kg: f64,
    /// Substrate and deposit storage changes plus carried vapor energy, J.
    pub energy_residual_j: f64,
    pub pressure_iterations: usize,
}

impl ChamberExchange for SupportedSurfaceStep {
    fn phase(&self) -> &SurfaceStep {
        &self.phase
    }
    fn phase_mut(&mut self) -> &mut SurfaceStep {
        &mut self.phase
    }
}

/// Advance shared vapor, original substrate calorics and supported deposits.
///
/// The existing bounded pressure root calls `advance_supported_surface` at one
/// prescribed trial pressure for each surface; final vapor mass is reconstructed
/// from those same accepted transfers. This is first-order backward-Euler mass
/// coupling, not the independent finite-bath exponential used by a single
/// supported surface. Refine physical time to resolve latent/backpressure
/// feedback. No surface borrows another surface's finite vapor copy.
///
/// `substrate_temperature(index, h)` is the original continuous nondecreasing
/// caloric chart for that substrate. Its mass and reference geometry stay fixed.
/// Source, inventory, budget or cancellation refusal returns no changed input.
/// An empty surface selection preserves the finite reservoir and zero accounts.
pub fn advance_supported_chamber(
    surfaces: &[SupportedChamberSurface<'_>],
    reservoir: VaporReservoir,
    duration_s: f64,
    config: ChamberSolverConfig,
    substrate_temperature: impl Fn(usize, f64) -> Result<f64, SublimationError>,
    mut checkpoint: impl FnMut() -> Result<(), SublimationError>,
) -> Result<SupportedChamberStep, SublimationError> {
    let solution = solve_chamber(
        surfaces.len(),
        reservoir,
        duration_s,
        config,
        |index| {
            let surface = surfaces[index];
            let state = surface.state;
            if !positive(state.substrate_mass_kg)
                || !state.substrate_specific_enthalpy_j_kg.is_finite()
                || !nonnegative(state.deposit_mass_kg)
                || !nonnegative(state.area_m2)
            {
                return Err(SublimationError::Invalid(
                    "invalid chamber substrate storage or deposit inventory",
                ));
            }
            let t = substrate_temperature(index, state.substrate_specific_enthalpy_j_kg)?;
            Ok(InitialSurface {
                temperature_k: t,
                properties: properties(surface.law, t)?,
                area_m2: state.area_m2,
                kinetics: surface.kinetics,
                has_deposit: state.deposit_mass_kg > 0.0,
            })
        },
        |index, reservoir, checkpoint| {
            let surface = surfaces[index];
            let value = advance_supported_surface(
                surface.law,
                surface.kinetics,
                surface.state,
                reservoir,
                duration_s,
                config.surface_solver,
                |h| substrate_temperature(index, h),
                surface.deposit_enthalpy,
                checkpoint,
            )?;
            let molar_mass =
                properties(surface.law, value.phase.surface_temperature_k)?.molar_mass_kg_mol;
            Ok((value, molar_mass))
        },
        &mut checkpoint,
    )?;
    let mut substrate_energy = 0.0;
    let mut substrate_correction = 0.0;
    let mut deposit_energy = 0.0;
    let mut deposit_correction = 0.0;
    for surface in &solution.surfaces {
        checkpoint()?;
        compensated_add(
            &mut substrate_energy,
            &mut substrate_correction,
            surface.substrate_energy_change_j,
        );
        compensated_add(
            &mut deposit_energy,
            &mut deposit_correction,
            surface.deposit_energy_change_j,
        );
    }
    Ok(SupportedChamberStep {
        surfaces: solution.surfaces,
        reservoir: solution.reservoir,
        transferred_mass_kg: solution.transferred_mass_kg,
        latent_energy_j: solution.latent_energy_j,
        vapor_energy_j: solution.vapor_energy_j,
        substrate_energy_change_j: finite(substrate_energy + substrate_correction)?,
        deposit_energy_change_j: finite(deposit_energy + deposit_correction)?,
        mass_residual_kg: solution.mass_residual_kg,
        nonlinear_mass_residual_kg: solution.nonlinear_mass_residual_kg,
        energy_residual_j: solution.energy_residual_j,
        pressure_iterations: solution.pressure_iterations,
    })
}

#[cfg(test)]
mod tests;
