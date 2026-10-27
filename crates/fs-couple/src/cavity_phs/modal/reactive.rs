//! Source-defined gas reactions drive the original material-derived elastic wall.
//!
//! Chemistry changes the original closed parcel at fixed volume and conserved
//! total internal energy, including formation references. Its new composition,
//! temperature and equation of state then enter the native simultaneous gas
//! pressure/work and elastic-motion solve. No fitted force or reaction-heat
//! preset is introduced, and no separate piston replaces the mesh-derived wall.
//!
//! This is a first-order reaction-then-motion split over one native modal sample.
//! Compression and external heat influence the next chemical sample. It does
//! not claim a simultaneous implicit reacting-fluid/structure solve, flame
//! propagation, chemical acoustic relaxation or structural heat transfer. The
//! structural material/geometry epoch and observation transfers remain fixed.
//!
//! The separate [`coupled`] entry solves endpoint chemistry and pressure/work
//! together while retaining this same native wall and volume chart.

pub mod coupled;

use super::{GasValveMotion, ModalGasError, ModalGasSolve, ModalGasStep};
use crate::cavity_phs::piston::reactive::{
    GasReaction, GasReactionError, GasReactionSolve, GasReactionStep,
};
use fs_exec::Cx;

/// A refusal never changes the accepted gas, anchored volume chart or modes.
#[derive(Debug)]
pub enum ReactiveModalGasError {
    Reaction(GasReactionError),
    Motion(ModalGasError),
    Invalid(&'static str),
    Cancelled,
}

impl core::fmt::Display for ReactiveModalGasError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "reactive gas/modal wall coupling: {self:?}")
    }
}
impl std::error::Error for ReactiveModalGasError {}
impl From<GasReactionError> for ReactiveModalGasError {
    fn from(error: GasReactionError) -> Self {
        match error {
            GasReactionError::Cancelled => Self::Cancelled,
            other => Self::Reaction(other),
        }
    }
}
impl From<ModalGasError> for ReactiveModalGasError {
    fn from(error: ModalGasError) -> Self {
        match error {
            ModalGasError::Cancelled => Self::Cancelled,
            other => Self::Motion(other),
        }
    }
}

/// Original chemical and mechanical stages, with their enclosing energy account.
#[derive(Clone, Debug)]
pub struct ReactiveModalGasStep {
    /// Original inventory, rates, source evaluations and isolated chemical account.
    pub reaction: GasReactionStep,
    /// Original mechanical report, whose initial gas is the reacted parcel.
    pub mechanical: ModalGasStep,
    /// Actual total gas internal-energy change from BEFORE chemistry [J].
    pub gas_energy_change_j: f64,
    /// Gas + modal storage + native modal loss + exterior work - external heat [J].
    pub total_energy_residual_j: f64,
}

/// One complete proposed state; previewing it advances neither accepted owner.
#[derive(Clone, Debug)]
pub struct ReactiveModalGasCandidate {
    motion: GasValveMotion,
    step: ReactiveModalGasStep,
}

impl ReactiveModalGasCandidate {
    pub const fn motion(&self) -> &GasValveMotion {
        &self.motion
    }
    pub const fn step(&self) -> &ReactiveModalGasStep {
        &self.step
    }
    pub fn into_parts(self) -> (GasValveMotion, ReactiveModalGasStep) {
        (self.motion, self.step)
    }
}

impl GasValveMotion {
    /// Propose one chemically reacting native modal sample without publication.
    ///
    /// The original reaction object retains its selected explicit, implicit,
    /// multichannel, subcycled and material-rate policies. It receives exactly
    /// the native motion sample duration. Chemical substeps do not advance the
    /// mechanical clock; nonlinear volume trials do not repeat chemistry.
    ///
    /// The accepted owner is cloned rather than reconstructed: its reference
    /// volume and reference modal coordinates, material identity, equilibrium,
    /// momenta and observation transfers survive composition changes. The
    /// native mechanical solve reevaluates its gas-loaded sampling gate using
    /// the reacted, midpoint and endpoint states.
    ///
    /// `heat_j` is external heat only. Formation energy is already included in
    /// gas storage, and physical modal loss leaves this enclosure; neither is
    /// credited again as gas heating. All original chemistry and mechanical
    /// admission checks apply before the whole-interval energy check below.
    pub fn reactive_step_candidate(
        &self,
        reaction: &GasReaction,
        heat_j: f64,
        exterior_pressure_pa: f64,
        chemistry: GasReactionSolve,
        controls: ModalGasSolve,
        cx: &Cx<'_>,
    ) -> Result<ReactiveModalGasCandidate, ReactiveModalGasError> {
        super::poll(cx)?;
        let duration_s = self.motion.model().sample_period_s();
        let (gas, reaction) = reaction
            .parcel_candidate(&self.gas, duration_s, chemistry, || cx.checkpoint().is_err())?
            .into_parts();
        // In particular, do not call Self::new with the reacted current volume:
        // that would reset the affine chart and permit accepted residual drift.
        let mut candidate = self.clone();
        candidate.gas = gas;
        let mechanical = candidate.step(heat_j, exterior_pressure_pa, controls, cx)?;
        let gas_energy_change_j = self.gas.mass_kg()
            * (candidate.gas.state().specific_internal_energy_j_kg()
                - self.gas.state().specific_internal_energy_j_kg());
        let total_energy_residual_j = gas_energy_change_j
            + mechanical.modal_energy_change_j
            + mechanical.motion.modal.viscous_dissipation_j
            + mechanical.exterior_pressure_work_j
            - heat_j;
        if !gas_energy_change_j.is_finite()
            || !total_energy_residual_j.is_finite()
            || total_energy_residual_j.abs() > controls.energy_tolerance_j
        {
            return Err(ReactiveModalGasError::Invalid(
                "reactive gas storage and original modal pressure work do not close whole-interval energy",
            ));
        }
        super::poll(cx)?;
        Ok(ReactiveModalGasCandidate {
            motion: candidate,
            step: ReactiveModalGasStep {
                reaction,
                mechanical,
                gas_energy_change_j,
                total_energy_residual_j,
            },
        })
    }

    /// Publish both gas and modal histories only after the complete proposal passes.
    /// A chemical, mechanical, source-domain, budget or cancellation refusal
    /// leaves the original state retryable.
    pub fn step_reactive(
        &mut self,
        reaction: &GasReaction,
        heat_j: f64,
        exterior_pressure_pa: f64,
        chemistry: GasReactionSolve,
        controls: ModalGasSolve,
        cx: &Cx<'_>,
    ) -> Result<ReactiveModalGasStep, ReactiveModalGasError> {
        let (candidate, step) = self
            .reactive_step_candidate(
                reaction,
                heat_j,
                exterior_pressure_pa,
                chemistry,
                controls,
                cx,
            )?
            .into_parts();
        super::poll(cx)?;
        *self = candidate;
        Ok(step)
    }
}
