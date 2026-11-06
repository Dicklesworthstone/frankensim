//! Original-material binding for nonlinear axisymmetric slender-jet evolution.
//!
//! The immutable resolved liquid remains attached to its evolving numerical
//! state. Density, viscosity and tension come from one admitted source point;
//! they are never reconstructed from a material name or a separate flow score.

use super::super::nonlinear::{
    NonlinearJetError, PeriodicSlenderJet, SlenderJetAdvance, SlenderJetEvolution,
    SlenderJetLimits, SlenderJetParameters, SlenderJetStepBudget,
};
use fs_exec::CancelGate;
use fs_material::liquid::surface::ResolvedLiquidSurfaceState;

/// A nonlinear periodic jet and the exact immutable liquid source it consumes.
///
/// The state is evolved with constant Newtonian coefficients. A changed
/// temperature, composition or pressure requires a newly resolved material and
/// a newly admitted initial state. This owner does not silently rebind an
/// existing trajectory or freeze a source requiring evolving rheology.
pub struct MaterialPeriodicSlenderJet<'a> {
    material: &'a ResolvedLiquidSurfaceState,
    jet: PeriodicSlenderJet,
}

impl<'a> MaterialPeriodicSlenderJet<'a> {
    /// Initialize a periodic jet from physical radii and axial velocities.
    ///
    /// Radii are sampled at cell centers, velocities at cell faces. Length,
    /// radii and velocities are SI quantities. The native owner validates
    /// geometry, representability and the declared approximation limits.
    ///
    /// # Errors
    /// Refuses cancellation, omitted surface tension, evolving non-Newtonian
    /// rheology, or any inadmissible native initial state.
    pub fn try_new(
        material: &'a ResolvedLiquidSurfaceState,
        length_m: f64,
        radii_m: &[f64],
        velocities_m_s: &[f64],
        limits: SlenderJetLimits,
        gate: &CancelGate,
    ) -> Result<Self, NonlinearJetError> {
        if gate.is_requested() {
            return Err(NonlinearJetError::Cancelled);
        }
        if material.rheology_axes().is_some() || material.rheology_pressure_axis().is_some() {
            return Err(NonlinearJetError::InvalidInput(
                "a constant Newtonian jet cannot freeze an evolving rheology source",
            ));
        }
        let surface_tension_n_m =
            material
                .surface_tension_n_m()
                .ok_or(NonlinearJetError::InvalidInput(
                    "a material capillary jet requires explicitly sourced surface tension",
                ))?;
        let parameters = SlenderJetParameters {
            density_kg_m3: material.density_kg_m3(),
            surface_tension_n_m,
            kinematic_viscosity_m2_s: material.kinematic_viscosity_m2_s(),
        };
        let jet =
            PeriodicSlenderJet::new(parameters, length_m, radii_m, velocities_m_s, limits, gate)?;
        Ok(Self { material, jet })
    }

    /// Exact source, including property selections, query coordinates and receipts.
    #[must_use]
    pub const fn material(&self) -> &'a ResolvedLiquidSurfaceState {
        self.material
    }

    /// Accepted native state; coefficients cannot be replaced through this view.
    #[must_use]
    pub const fn state(&self) -> &PeriodicSlenderJet {
        &self.jet
    }

    /// Evolve the original sourced state with the native conservative operator.
    ///
    /// Domain-limit outcomes retain the last admitted substep. Cancellation,
    /// numerical failure and exhausted work budgets preserve the entry state.
    pub fn advance(
        &mut self,
        duration_s: f64,
        budget: SlenderJetStepBudget,
        gate: &CancelGate,
    ) -> Result<SlenderJetAdvance, NonlinearJetError> {
        self.jet.advance(duration_s, budget, gate)
    }
}

impl SlenderJetEvolution for MaterialPeriodicSlenderJet<'_> {
    fn state(&self) -> &PeriodicSlenderJet {
        &self.jet
    }

    fn advance(
        &mut self,
        duration_s: f64,
        budget: SlenderJetStepBudget,
        gate: &CancelGate,
    ) -> Result<SlenderJetAdvance, NonlinearJetError> {
        self.jet.advance(duration_s, budget, gate)
    }
}
