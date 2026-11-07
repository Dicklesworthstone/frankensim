//! Fatigue carried by the original accepted generalized-alpha plate trajectory.
//!
//! The native mass, stiffness, damping and force callbacks retain ownership of
//! mechanics. Only accepted endpoints enter the shared DKT stress-life observer.
//! A fatigue/domain/cancellation refusal commits neither motion nor damage.
//! Initial acceleration is solved from the actual physical residual. Checkpoint
//! clones retain q/v/a and all unresolved rainflow histories together.
//!
//! This is a fixed-condition elastic plate estimate. Rayleigh damping describes
//! an explicitly selected macroscopic loss; it is not an additional local
//! material stress. Constitutive viscosity and stiffness rescaling refuse until
//! their local stress and calibration bindings can be supplied consistently.
//! Physical timestep refinement is necessary to resolve stress turning points.

use fs_material::fatigue::FatigueSummary;
use fs_plate::fatigue::{FatiguePlate, PlateFatigueError, PlateFatigueObservation};
use fs_time::galpha::{
    ImplicitStepTelemetry, OperatorGeneralizedAlpha, SecondOrderState, TimeSolveError,
    initialization::{InitialSolveConfig, InitialSolveError, second_order_acceleration},
    second_order_adjoint::trajectory::StructuralTrajectoryModel,
};

use crate::{PlateDynamics, PlateLoad};

#[derive(Debug)]
pub enum PlateFatigueRunError {
    Invalid(&'static str),
    Fatigue(PlateFatigueError),
    Initialization(InitialSolveError),
    Mechanics(TimeSolveError),
    Forcing(String),
}
impl core::fmt::Display for PlateFatigueRunError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native plate fatigue trajectory refused: {self:?}")
    }
}
impl std::error::Error for PlateFatigueRunError {}
impl From<PlateFatigueError> for PlateFatigueRunError {
    fn from(error: PlateFatigueError) -> Self {
        Self::Fatigue(error)
    }
}
impl From<TimeSolveError> for PlateFatigueRunError {
    fn from(error: TimeSolveError) -> Self {
        Self::Mechanics(error)
    }
}

/// Initial accepted physical state and explicit original mass-solve controls.
pub struct PlateFatigueInitial<'a> {
    pub time_s: f64,
    pub displacement: &'a [f64],
    pub velocity: &'a [f64],
    pub solve: InitialSolveConfig,
    pub maximum_workspace_components: usize,
}

/// One immutable dynamics model and its bounded, jointly accepted fatigue state.
pub struct PlateFatigueRun<'a, L: PlateLoad + ?Sized> {
    dynamics: PlateDynamics<'a, L>,
    method: OperatorGeneralizedAlpha,
    state: SecondOrderState,
    observation: PlateFatigueObservation<'a>,
}
impl<L: PlateLoad + ?Sized> Clone for PlateFatigueRun<'_, L> {
    fn clone(&self) -> Self {
        Self {
            dynamics: self.dynamics,
            method: self.method,
            state: self.state.clone(),
            observation: self.observation.clone(),
        }
    }
}

impl<'a, L: PlateLoad + ?Sized> PlateFatigueRun<'a, L> {
    /// Bind exactly the source plate's pencil and initialize native acceleration.
    /// A foreign equally sized model cannot masquerade as the gauge owner.
    pub fn new(
        plate: &'a FatiguePlate,
        dynamics: PlateDynamics<'a, L>,
        method: OperatorGeneralizedAlpha,
        initial: PlateFatigueInitial<'_>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, PlateFatigueRunError> {
        poll(cancelled)?;
        if !core::ptr::eq(plate.model(), dynamics.model())
            || dynamics.parameters().stiffness_scale != 1.0
            || dynamics.bending_viscosity().is_some()
        {
            return Err(PlateFatigueRunError::Invalid(
                "fatigue needs its original unscaled elastic plate without unobserved material viscous stress",
            ));
        }
        let n = plate.model().free;
        if !initial.time_s.is_finite()
            || initial.displacement.len() != n
            || initial.velocity.len() != n
            || initial.velocity.iter().any(|x| !x.is_finite())
        {
            return Err(PlateFatigueRunError::Invalid(
                "finite initial time and original q/v dimensions required",
            ));
        }
        // This also checks the method's stage clock before accepting a run.
        method.forcing_time(initial.time_s)?;
        let observation = plate.observe(initial.displacement, cancelled)?;
        let mut force = vec![f64::NAN; n];
        dynamics
            .forcing(initial.time_s, &mut force)
            .map_err(PlateFatigueRunError::Forcing)?;
        poll(cancelled)?;
        if force.iter().any(|f| !f.is_finite()) {
            return Err(PlateFatigueRunError::Invalid(
                "initial force must be finite and completely written",
            ));
        }
        let acceleration = second_order_acceleration(
            &dynamics,
            initial.displacement,
            initial.velocity,
            &force,
            &dynamics.mass_preconditioner(),
            initial.solve,
            initial.maximum_workspace_components,
            cancelled,
        )
        .map_err(PlateFatigueRunError::Initialization)?;
        let state = SecondOrderState::new(
            initial.time_s,
            initial.displacement,
            initial.velocity,
            &acceleration.value,
        );
        poll(cancelled)?;
        Ok(Self {
            dynamics,
            method,
            state,
            observation,
        })
    }

    pub fn state(&self) -> &SecondOrderState {
        &self.state
    }
    pub fn dynamics(&self) -> &PlateDynamics<'a, L> {
        &self.dynamics
    }
    pub fn observation(&self) -> &PlateFatigueObservation<'a> {
        &self.observation
    }
    pub fn summaries(&self) -> Vec<FatigueSummary> {
        self.observation.summaries()
    }
    pub fn preview_finish(
        &self,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<FatigueSummary>, PlateFatigueRunError> {
        Ok(self.observation.preview_finish(cancelled)?)
    }

    /// Advance with the original implicit solver and actual intermediate-time
    /// force. Only the latest telemetry and bounded mechanical state are kept.
    pub fn advance(
        &mut self,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ImplicitStepTelemetry, PlateFatigueRunError> {
        poll(cancelled)?;
        if self.observation.is_finished() {
            return Err(PlateFatigueError::Finished.into());
        }
        let time = self.method.forcing_time(self.state.t)?;
        let mut force = vec![f64::NAN; self.state.q.len()];
        self.dynamics
            .forcing(time, &mut force)
            .map_err(PlateFatigueRunError::Forcing)?;
        poll(cancelled)?;
        if force.iter().any(|f| !f.is_finite()) {
            return Err(PlateFatigueRunError::Invalid(
                "force must be finite and completely written",
            ));
        }
        let mut candidate = self.state.clone();
        let telemetry =
            self.method
                .step_controlled(&mut candidate, &self.dynamics, &force, cancelled)?;
        if candidate
            .q
            .iter()
            .chain(&candidate.v)
            .chain(&candidate.a)
            .any(|x| !x.is_finite())
        {
            return Err(PlateFatigueRunError::Invalid(
                "native endpoint must remain finite",
            ));
        }
        candidate.history.clear();
        // The observer already stages every gauge and polls before publication.
        // After it accepts, moving the mechanical candidate cannot fail. Avoid
        // cloning its bounded residues twice on every physical timestep.
        self.observation.push(&candidate.q, cancelled)?;
        self.state = candidate;
        Ok(telemetry)
    }

    /// Finalize endpoint half cycles exactly once, without changing mechanics.
    pub fn finish(
        &mut self,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<FatigueSummary>, PlateFatigueRunError> {
        Ok(self.observation.finish(cancelled)?)
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), PlateFatigueRunError> {
    if cancelled() {
        Err(PlateFatigueError::Cancelled.into())
    } else {
        Ok(())
    }
}
