//! Source heating changes an existing free body's motion without resetting it.
//!
//! Each actual thermal substep advances the original rigid-body integrator at
//! its starting inertia, then changes configuration at conserved world momenta.
//! This first-order thermal/mechanical split uses the same material inventory,
//! source expansion and principal-frame lowering as static configuration.
//! Parameter work is reported and bounded separately from caloric heat and
//! integration error. Reciprocal thermoelastic heat is not inferred.
/// Expanded no-slip spherical boundaries in the original material fluid.
pub mod creeping_flow;
/// Source-bound viscous interfaces with reciprocal finite-body heat feedback.
#[cfg(feature = "dissipative-coupling")]
pub mod damping;
/// Expanded displacement geometry and buoyancy on the same thermal clock.
pub mod hydrostatic;
/// Material-coordinate external impulse events on the same accepted trajectory.
pub mod loads;

use super::*;
use crate::principal::transfer::ConfigurationTransferControls;
use crate::principal::{PrincipalFrame, PrincipalFrameError};
use fs_ga::{Se3, Twist};
use fs_mbd::{DynamicsDiagnostics, DynamicsError, RigidBodyIntegrator, RigidBodyState, Wrench};

/// Explicit accuracy and unresolved-work allowances; no material damping knob.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HeatedMotionControls {
    /// Upper bound on each mechanical/thermal interval [s], independently of
    /// the thermal transport's own conductive CFL limit.
    pub maximum_step_s: f64,
    /// Entrywise principal-frame tensor residual, in (0, 1e-8].
    pub relative_tensor_tolerance: f64,
    /// Momentum residual per update and sum of absolute parameter work over
    /// one advance. Cooling cannot cancel prior heating to evade this bound.
    pub transfer: ConfigurationTransferControls,
    /// Sum of absolute fixed-inertia mechanical energy defects in one advance
    /// [J]. This bounds the numerical account, not physical loss or accuracy.
    pub maximum_numerical_energy_change_j: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum HeatedMotionError {
    Heating(HeatingError),
    Principal(PrincipalFrameError),
    Dynamics(DynamicsError),
    Hydrostatic(crate::hydrostatic::HydrostaticStepError),
    CreepingFlow(crate::creeping_flow::CreepingFlowError),
    Invalid(&'static str),
    Budget(&'static str),
}
impl core::fmt::Display for HeatedMotionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "heated rigid-body motion: {self:?}")
    }
}
impl std::error::Error for HeatedMotionError {}
impl From<HeatingError> for HeatedMotionError {
    fn from(e: HeatingError) -> Self {
        Self::Heating(e)
    }
}
impl From<DepositionError> for HeatedMotionError {
    fn from(e: DepositionError) -> Self {
        Self::Heating(e.into())
    }
}
impl From<PrincipalFrameError> for HeatedMotionError {
    fn from(e: PrincipalFrameError) -> Self {
        Self::Principal(e)
    }
}
impl From<DynamicsError> for HeatedMotionError {
    fn from(e: DynamicsError) -> Self {
        Self::Dynamics(e)
    }
}
impl From<crate::hydrostatic::HydrostaticError> for HeatedMotionError {
    fn from(e: crate::hydrostatic::HydrostaticError) -> Self {
        Self::Hydrostatic(e.into())
    }
}
impl From<crate::hydrostatic::HydrostaticStepError> for HeatedMotionError {
    fn from(e: crate::hydrostatic::HydrostaticStepError) -> Self {
        Self::Hydrostatic(e)
    }
}
impl From<crate::creeping_flow::CreepingFlowError> for HeatedMotionError {
    fn from(e: crate::creeping_flow::CreepingFlowError) -> Self {
        Self::CreepingFlow(e)
    }
}

/// Finite heat and actual continuing motion on one accepted time interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HeatedMotionAudit {
    pub thermal: EnergyAudit,
    pub mechanical_before: DynamicsDiagnostics,
    pub mechanical_after: DynamicsDiagnostics,
    /// Signed inertia-change work; external to the one-way caloric solve [J].
    pub parameter_work_j: f64,
    /// Sum of absolute inertia-change work [J], so opposite signs cannot hide it.
    pub absolute_parameter_work_j: f64,
    /// Signed fixed-inertia integration error, never deposited as material heat.
    pub numerical_energy_change_j: f64,
    pub absolute_numerical_energy_change_j: f64,
    /// Thermal + mechanical storage change minus heat, parameter work and
    /// numerical change. This bookkeeping closure is not reciprocal coupling.
    pub combined_residual_j: f64,
    pub maximum_transfer_momentum_residual: f64,
}

/// Immutable candidate with one source heat history and one mechanical state.
///
/// This is free flight under the supplied integrator's uniform gravity. Thermal
/// fixed points are coordinates in the moving body frame, not external world
/// supports. Shape evolution is quasi-static internally: deformation kinetic
/// energy, thermomechanical stress, contact and drag need their existing owners.
/// Caloric and thermal geometry assumptions are inherited from `HeatedBody`.
/// Clone is a complete rollback/restart point, including the physical clock.
#[derive(Clone)]
pub struct HeatedRigidBody<'a> {
    thermal: DepositionSimulation<'a>,
    motion: RigidExpansion<'a>,
}

impl<'a> HeatedRigidBody<'a> {
    /// Start motion at this body's accepted thermal time and configuration.
    /// Kinematics are specified once; later heating preserves momenta instead
    /// of rebuilding a new state at the same angular velocity.
    pub fn new(
        body: HeatedBody<'a>,
        world_from_body: Se3,
        twist_body: Twist,
        integrator: RigidBodyIntegrator,
        controls: HeatedMotionControls,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, HeatedMotionError> {
        checkpoint()?;
        controls.transfer.validate()?;
        if !controls.maximum_step_s.is_finite()
            || controls.maximum_step_s <= 0.0
            || !controls.maximum_numerical_energy_change_j.is_finite()
            || controls.maximum_numerical_energy_change_j < 0.0
        {
            return Err(HeatedMotionError::Invalid(
                "positive step and finite nonnegative energy budget required",
            ));
        }
        let principal = PrincipalFrame::from_inertia(
            body.configuration().inertia(),
            controls.relative_tensor_tolerance,
            || false,
        )?;
        let state = principal.state_from_kinematics(world_from_body, twist_body)?;
        integrator.diagnostics(state, principal.properties())?;
        checkpoint()?;
        Ok(Self {
            motion: RigidExpansion {
                time_s: body.thermal.time_s(),
                expansion: body.mechanics,
                principal,
                state,
                integrator,
                controls,
                work: MotionWork::default(),
            },
            thermal: body.thermal,
        })
    }

    pub const fn thermal(&self) -> &DepositionSimulation<'a> {
        &self.thermal
    }
    pub const fn configuration(&self) -> &ExpandedBody<'a> {
        &self.motion.expansion.configuration
    }
    pub const fn reference(&self) -> &'a ResolvedBody {
        self.motion.expansion.reference
    }
    pub const fn principal_frame(&self) -> &PrincipalFrame {
        &self.motion.principal
    }
    pub const fn state(&self) -> RigidBodyState {
        self.motion.state
    }
    pub fn world_from_body(&self) -> Result<Se3, HeatedMotionError> {
        Ok(self.motion.principal.world_from_body(self.motion.state)?)
    }
    pub fn diagnostics(&self) -> Result<DynamicsDiagnostics, HeatedMotionError> {
        Ok(self
            .motion
            .integrator
            .diagnostics(self.motion.state, self.motion.principal.properties())?)
    }

    /// Heat, evolve the existing trajectory, and publish both histories only
    /// after every thermal substep, source and numerical bound is accepted.
    pub fn advance_to(
        &self,
        target_s: f64,
        config: StepConfig,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(Self, HeatedMotionAudit), HeatedMotionError> {
        let (thermal, motion, audit) =
            self.thermal
                .advance_observed(&self.motion, target_s, config, checkpoint)?;
        let next = Self { thermal, motion };
        let audit = self.audit(&next, audit)?;
        Ok((next, audit))
    }

    /// Compose generated electrical/frictional heat and source transport on
    /// these same substeps. The original observer's heat is counted once.
    pub fn advance_with<O: DepositionObserver>(
        &self,
        observer: &O,
        target_s: f64,
        config: StepConfig,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(Self, O, HeatedMotionAudit), ObserverCouplingError<O::Error, HeatedMotionError>>
    {
        let pair = CoupledObservers::new(observer.clone(), self.motion.clone());
        let (thermal, pair, audit) = self
            .thermal
            .advance_observed(&pair, target_s, config, checkpoint)?;
        let (observer, motion) = pair.into_parts();
        let next = Self { thermal, motion };
        let audit = self
            .audit(&next, audit)
            .map_err(ObserverCouplingError::Second)?;
        Ok((next, observer, audit))
    }

    fn audit(
        &self,
        next: &Self,
        thermal: EnergyAudit,
    ) -> Result<HeatedMotionAudit, HeatedMotionError> {
        let before = self.diagnostics()?;
        let after = next.diagnostics()?;
        let work = next.motion.work;
        let combined_residual_j = finite(
            thermal.residual_j + mechanical_change(before, after)?
                - work.parameter_j
                - work.numerical_j,
        )?;
        Ok(HeatedMotionAudit {
            thermal,
            mechanical_before: before,
            mechanical_after: after,
            parameter_work_j: work.parameter_j,
            absolute_parameter_work_j: work.absolute_parameter_j,
            numerical_energy_change_j: work.numerical_j,
            absolute_numerical_energy_change_j: work.absolute_numerical_j,
            combined_residual_j,
            maximum_transfer_momentum_residual: work.maximum_momentum_residual,
        })
    }
}

#[derive(Clone, Copy, Default)]
struct MotionWork {
    parameter_j: f64,
    absolute_parameter_j: f64,
    numerical_j: f64,
    absolute_numerical_j: f64,
    maximum_momentum_residual: f64,
}

#[derive(Clone)]
struct RigidExpansion<'a> {
    expansion: ExpansionObserver<'a>,
    principal: PrincipalFrame,
    state: RigidBodyState,
    time_s: f64,
    integrator: RigidBodyIntegrator,
    controls: HeatedMotionControls,
    work: MotionWork,
}
impl DepositionObserver for RigidExpansion<'_> {
    type Error = HeatedMotionError;
    fn begin(&mut self) -> Result<(), Self::Error> {
        self.work = MotionWork::default();
        Ok(())
    }
    fn event(
        &mut self,
        event: &Event,
        thermal: &DepositionSimulation<'_>,
    ) -> Result<(), Self::Error> {
        Ok(self.expansion.event(event, thermal)?)
    }
    fn limit_timestep_s(
        &mut self,
        _: &DepositionSimulation<'_>,
        proposed_s: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<f64, Self::Error> {
        checkpoint()?;
        Ok(proposed_s.min(self.controls.maximum_step_s))
    }
    fn step(
        &mut self,
        start_s: f64,
        dt_s: f64,
        old: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        checkpoint()?;
        if start_s != self.time_s {
            return Err(HeatedMotionError::Invalid(
                "mechanical and thermal clocks disagree",
            ));
        }
        let step =
            self.integrator
                .step(self.state, self.principal.properties(), Wrench::ZERO, dt_s)?;
        let numerical_j = mechanical_change(step.diagnostics_before, step.diagnostics_after)?;
        self.finish_step(
            start_s,
            dt_s,
            old,
            thermal,
            step.state_after,
            numerical_j,
            checkpoint,
        )
    }
}

impl RigidExpansion<'_> {
    /// Accept the existing mechanics owner's candidate and refresh the SAME
    /// thermal configuration and momentum map for every supported load model.
    fn finish_step(
        &mut self,
        start_s: f64,
        dt_s: f64,
        old: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        state_after: RigidBodyState,
        numerical_j: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), HeatedMotionError> {
        self.expansion
            .step(start_s, dt_s, old, thermal, checkpoint)?;
        let principal = PrincipalFrame::from_inertia(
            self.expansion.configuration.inertia(),
            self.controls.relative_tensor_tolerance,
            || false,
        )?;
        let transfer = principal.transfer_state_from(
            &self.principal,
            state_after,
            self.controls.transfer,
            || false,
        )?;
        self.work.parameter_j = finite(self.work.parameter_j + transfer.parameter_work_j)?;
        self.work.absolute_parameter_j =
            finite(self.work.absolute_parameter_j + transfer.parameter_work_j.abs())?;
        self.work.numerical_j = finite(self.work.numerical_j + numerical_j)?;
        self.work.absolute_numerical_j =
            finite(self.work.absolute_numerical_j + numerical_j.abs())?;
        self.work.maximum_momentum_residual = self
            .work
            .maximum_momentum_residual
            .max(transfer.relative_momentum_residual);
        if self.work.absolute_parameter_j > self.controls.transfer.maximum_parameter_work_j {
            return Err(HeatedMotionError::Budget(
                "total absolute configuration work",
            ));
        }
        if self.work.absolute_numerical_j > self.controls.maximum_numerical_energy_change_j {
            return Err(HeatedMotionError::Budget(
                "total absolute fixed-inertia numerical energy change",
            ));
        }
        checkpoint()?;
        self.principal = principal;
        self.state = transfer.state_after;
        self.time_s = thermal.time_s();
        Ok(())
    }
}

fn finite(value: f64) -> Result<f64, HeatedMotionError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(HeatedMotionError::Invalid(
            "nonfinite mechanical work accounting",
        ))
    }
}
fn mechanical_change(
    before: DynamicsDiagnostics,
    after: DynamicsDiagnostics,
) -> Result<f64, HeatedMotionError> {
    finite(
        (after.translational_kinetic_energy - before.translational_kinetic_energy)
            + (after.rotational_kinetic_energy - before.rotational_kinetic_energy)
            + (after.gravitational_potential_energy - before.gravitational_potential_energy),
    )
}

/// Preserve the exact thermal owner's checkpoint error through boolean native
/// cancellation adapters, rather than replacing a budget/source error.
fn checked<T, E: Into<HeatedMotionError>>(
    checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    operation: impl FnOnce(&mut dyn FnMut() -> bool) -> Result<T, E>,
) -> Result<T, HeatedMotionError> {
    let mut failure = None;
    let result = operation(&mut || match checkpoint() {
        Ok(()) => false,
        Err(error) => {
            failure = Some(error);
            true
        }
    });
    match failure {
        Some(error) => Err(error.into()),
        None => result.map_err(Into::into),
    }
}
