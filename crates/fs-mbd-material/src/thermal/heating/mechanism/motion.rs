//! Source heat and continuing articulated motion on one accepted timeline.
//!
//! Each thermal substep advances native ABA dynamics through fs-time RK45 at
//! its starting material configuration. The endpoint's source-resolved inertias
//! then receive the same conjugate joint momenta at unchanged joint positions.
//! This is a first-order split with a stationary prescribed base. Parameter
//! work is reported separately and never deposited as caloric heat.
use super::*;
use fs_ga::{Se3, Vec3, Wrench};
use fs_mbd::articulated::{
    ArticulatedError, ArticulatedModel, BaseState, forward_dynamics, forward_kinematics,
    stationary_base_joint_momentum, stationary_base_joint_velocity,
};
use fs_time::{
    AdaptiveError, AdaptiveState, AdaptiveStatus, PiController, rk45_adaptive_controlled,
};
use std::cell::RefCell;

/// Component tolerances follow the actual generalized-coordinate ordering.
#[derive(Clone, Debug, PartialEq)]
pub struct Controls {
    /// Thermal/mechanical split interval cap [s]. Refinement is still required.
    pub maximum_step_s: f64,
    pub relative_state_tolerance: f64,
    /// Absolute position tolerances [m or rad], one per scalar joint.
    pub absolute_position_tolerance: Vec<f64>,
    /// Absolute rate tolerances [m/s or rad/s], one per scalar joint.
    pub absolute_velocity_tolerance: Vec<f64>,
    pub maximum_rk_attempts_per_substep: usize,
    pub relative_momentum_tolerance: f64,
    /// Absolute momentum tolerances [N s or N m s], one per scalar joint.
    pub absolute_momentum_tolerance: Vec<f64>,
    /// Sum of absolute thermal-configuration work over ONE advance [J].
    pub maximum_parameter_work_j: f64,
    /// Sum of absolute fixed-configuration integration defects [J].
    pub maximum_numerical_energy_change_j: f64,
}
#[derive(Clone, Debug, PartialEq)]
pub enum Error {
    Heating(HeatingError),
    Dynamics(ArticulatedError),
    Integration(AdaptiveError),
    Invalid(&'static str),
    Budget(&'static str),
}
impl core::fmt::Display for Error {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "heated articulated motion: {self:?}")
    }
}
impl std::error::Error for Error {}
impl From<HeatingError> for Error {
    fn from(e: HeatingError) -> Self {
        Self::Heating(e)
    }
}
impl From<DepositionError> for Error {
    fn from(e: DepositionError) -> Self {
        Self::Heating(e.into())
    }
}
impl From<ArticulatedError> for Error {
    fn from(e: ArticulatedError) -> Self {
        Self::Dynamics(e)
    }
}
impl From<AdaptiveError> for Error {
    fn from(e: AdaptiveError) -> Self {
        Self::Integration(e)
    }
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MechanicalEnergy {
    pub kinetic_j: f64,
    /// Datum is zero at world height g dot x = 0.
    pub gravitational_j: f64,
}
impl MechanicalEnergy {
    fn change_from(self, old: Self) -> Result<f64, Error> {
        finite((self.kinetic_j - old.kinetic_j) + (self.gravitational_j - old.gravitational_j))
    }
}
/// Accounts cover the same interval. Numerical error and thermal parameter
/// work are distinct from heat and applied joint work.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Audit {
    pub thermal: EnergyAudit,
    pub mechanical_before: MechanicalEnergy,
    pub mechanical_after: MechanicalEnergy,
    pub applied_joint_work_j: f64,
    pub parameter_work_j: f64,
    pub absolute_parameter_work_j: f64,
    pub numerical_energy_change_j: f64,
    pub absolute_numerical_energy_change_j: f64,
    /// Largest component residual divided by its absolute+relative allowance.
    pub maximum_scaled_momentum_residual: f64,
    pub rk_attempts: usize,
    pub combined_residual_j: f64,
}

/// One material mechanism, heat history, joint state and physical clock.
///
/// Constant declared joint efforts and uniform gravity act on a stationary
/// base. Joints, reference placements and thermal transport geometry remain
/// prescribed. Free thermal expansion changes each link's retained inertia;
/// deformation velocities, moving joint anchors, contacts, closed kinematic
/// loops and reciprocal thermoelastic heat require other owners. The original
/// material mass is retained; this is not a material replacement event.
#[derive(Clone)]
pub struct HeatedArticulatedMotion<'a> {
    thermal: DepositionSimulation<'a>,
    motion: Motion<'a>,
}
impl<'a> HeatedArticulatedMotion<'a> {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        heated: HeatedArticulated<'a>,
        world_from_base: Se3,
        joint_position: Vec<f64>,
        joint_velocity: Vec<f64>,
        joint_effort: Vec<f64>,
        gravity_world: Vec3,
        controls: Controls,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, Error> {
        checkpoint()?;
        let n = heated.configuration().model().dof_count();
        controls.validate(n)?;
        if joint_position.len() != n || joint_velocity.len() != n || joint_effort.len() != n {
            return Err(Error::Invalid(
                "one position, velocity and effort per scalar joint required",
            ));
        }
        let motion = Motion {
            time_s: heated.thermal.time_s(),
            expansion: heated.mechanics,
            position: joint_position,
            velocity: joint_velocity,
            effort: joint_effort,
            base: BaseState::stationary(world_from_base),
            gravity: gravity_world,
            controls,
            work: Work::default(),
        };
        forward_dynamics(
            motion.expansion.configuration.model(),
            motion.base,
            &motion.position,
            &motion.velocity,
            &motion.effort,
            motion.gravity,
            &vec![Wrench::default(); motion.expansion.configuration.model().link_count()],
        )?;
        motion.energy()?;
        checkpoint()?;
        Ok(Self {
            thermal: heated.thermal,
            motion,
        })
    }
    pub const fn thermal(&self) -> &DepositionSimulation<'a> {
        &self.thermal
    }
    pub const fn configuration(&self) -> &ExpandedArticulated<'a> {
        &self.motion.expansion.configuration
    }
    pub const fn reference(&self) -> &'a ResolvedArticulated {
        self.motion.expansion.reference
    }
    pub fn joint_position(&self) -> &[f64] {
        &self.motion.position
    }
    pub fn joint_velocity(&self) -> &[f64] {
        &self.motion.velocity
    }
    pub fn joint_momentum(&self) -> Result<Vec<f64>, Error> {
        Ok(stationary_base_joint_momentum(
            self.configuration().model(),
            self.joint_position(),
            self.joint_velocity(),
        )?)
    }
    pub fn mechanical_energy(&self) -> Result<MechanicalEnergy, Error> {
        self.motion.energy()
    }

    /// Use the existing heat scheduler. Late cancellation, material-domain,
    /// joint-limit or work refusal returns no partial state. Clone is restart.
    pub fn advance_to(
        &self,
        target_s: f64,
        config: StepConfig,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(Self, Audit), Error> {
        let (thermal, motion, audit) =
            self.thermal
                .advance_observed(&self.motion, target_s, config, checkpoint)?;
        let next = Self { thermal, motion };
        let audit = self.audit(&next, audit)?;
        Ok((next, audit))
    }
    /// Generated heat and transport share the same accepted substeps; the
    /// supplied observer's original state and errors remain attached.
    pub fn advance_with<O: DepositionObserver>(
        &self,
        observer: &O,
        target_s: f64,
        config: StepConfig,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(Self, O, Audit), ObserverCouplingError<O::Error, Error>> {
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
    fn audit(&self, next: &Self, thermal: EnergyAudit) -> Result<Audit, Error> {
        let before = self.motion.energy()?;
        let after = next.motion.energy()?;
        let w = next.motion.work;
        Ok(Audit {
            thermal,
            mechanical_before: before,
            mechanical_after: after,
            applied_joint_work_j: w.applied,
            parameter_work_j: w.parameter,
            absolute_parameter_work_j: w.absolute_parameter,
            numerical_energy_change_j: w.numerical,
            absolute_numerical_energy_change_j: w.absolute_numerical,
            maximum_scaled_momentum_residual: w.momentum_residual,
            rk_attempts: w.attempts,
            combined_residual_j: finite(
                thermal.residual_j + after.change_from(before)?
                    - w.applied
                    - w.parameter
                    - w.numerical,
            )?,
        })
    }
}
impl Controls {
    fn validate(&self, n: usize) -> Result<(), Error> {
        if n == 0
            || !positive(self.maximum_step_s)
            || !nonnegative(self.relative_state_tolerance)
            || !nonnegative(self.relative_momentum_tolerance)
            || self.maximum_rk_attempts_per_substep == 0
            || !nonnegative(self.maximum_parameter_work_j)
            || !nonnegative(self.maximum_numerical_energy_change_j)
            || [
                &self.absolute_position_tolerance,
                &self.absolute_velocity_tolerance,
                &self.absolute_momentum_tolerance,
            ]
            .iter()
            .any(|v| v.len() != n || v.iter().any(|&x| !positive(x)))
        {
            return Err(Error::Invalid(
                "nonempty motion, positive component tolerances and step/attempt bounds, finite nonnegative relative/work allowances required",
            ));
        }
        Ok(())
    }
}
#[derive(Clone, Copy, Default)]
struct Work {
    applied: f64,
    parameter: f64,
    absolute_parameter: f64,
    numerical: f64,
    absolute_numerical: f64,
    momentum_residual: f64,
    attempts: usize,
}
#[derive(Clone)]
struct Motion<'a> {
    expansion: MechanismExpansion<'a>,
    time_s: f64,
    position: Vec<f64>,
    velocity: Vec<f64>,
    effort: Vec<f64>,
    base: BaseState,
    gravity: Vec3,
    controls: Controls,
    work: Work,
}
impl Motion<'_> {
    fn energy(&self) -> Result<MechanicalEnergy, Error> {
        energy(
            self.expansion.configuration.model(),
            self.base,
            &self.position,
            &self.velocity,
            self.gravity,
        )
    }
}
impl DepositionObserver for Motion<'_> {
    type Error = Error;
    fn begin(&mut self) -> Result<(), Error> {
        self.work = Work::default();
        Ok(())
    }
    fn event(&mut self, event: &Event, thermal: &DepositionSimulation<'_>) -> Result<(), Error> {
        Ok(self.expansion.event(event, thermal)?)
    }
    fn limit_timestep_s(
        &mut self,
        _: &DepositionSimulation<'_>,
        proposed_s: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<f64, Error> {
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
    ) -> Result<(), Error> {
        checkpoint()?;
        if start_s != self.time_s || !positive(dt_s) || thermal.time_s() <= start_s {
            return Err(Error::Invalid("joint and thermal clocks disagree"));
        }
        let model = self.expansion.configuration.model();
        let before = self.energy()?;
        let n = self.position.len();
        let scales: Vec<_> = self
            .controls
            .absolute_position_tolerance
            .iter()
            .chain(&self.controls.absolute_velocity_tolerance)
            .copied()
            .collect();
        let normalized: Vec<_> = self
            .position
            .iter()
            .chain(&self.velocity)
            .zip(&scales)
            .map(|(&u, &scale)| u / scale)
            .collect();
        let mut state = AdaptiveState::new(start_s, &normalized, dt_s);
        let external = vec![Wrench::default(); model.link_count()];
        let source_error = RefCell::new(None);
        let rhs = |_: f64, u: &[f64], out: &mut [f64]| {
            let physical: Vec<_> = u.iter().zip(&scales).map(|(&x, &s)| x * s).collect();
            match forward_dynamics(
                model,
                self.base,
                &physical[..n],
                &physical[n..],
                &self.effort,
                self.gravity,
                &external,
            ) {
                Ok(response) => {
                    for i in 0..n {
                        out[i] = physical[n + i] / scales[i];
                        out[n + i] = response.generalized_acceleration[i] / scales[n + i];
                    }
                }
                Err(error) => {
                    *source_error.borrow_mut() = Some(error);
                    out.fill(f64::NAN);
                }
            }
        };
        let mut checkpoint_error = None;
        let result = rk45_adaptive_controlled(
            &mut state,
            &rhs,
            thermal.time_s(),
            self.controls.relative_state_tolerance,
            1.0,
            &PiController::default(),
            self.controls.maximum_rk_attempts_per_substep,
            &mut || match checkpoint() {
                Ok(()) => false,
                Err(error) => {
                    checkpoint_error = Some(error);
                    true
                }
            },
        );
        if let Some(error) = checkpoint_error {
            return Err(error.into());
        }
        if let Some(error) = source_error.into_inner() {
            return Err(error.into());
        }
        let report = result?;
        if report.status != AdaptiveStatus::ReachedEnd {
            return Err(Error::Budget("mechanical RK45 attempts exhausted"));
        }
        let physical: Vec<_> = state.u.iter().zip(&scales).map(|(&x, &s)| x * s).collect();
        let position = physical[..n].to_vec();
        let velocity = physical[n..].to_vec();
        let after_motion = energy(model, self.base, &position, &velocity, self.gravity)?;
        let applied = finite(
            self.effort
                .iter()
                .zip(&position)
                .zip(&self.position)
                .map(|((&force, &q), &old_q)| force * (q - old_q))
                .sum(),
        )?;
        let numerical = finite(after_motion.change_from(before)? - applied)?;
        let momentum = stationary_base_joint_momentum(model, &position, &velocity)?;
        let old_model = model.clone();
        self.expansion
            .step(start_s, dt_s, old, thermal, checkpoint)?;
        let model = self.expansion.configuration.model();
        let next_velocity = if model == &old_model {
            velocity
        } else {
            stationary_base_joint_velocity(model, &position, &momentum)?
        };
        let recovered = stationary_base_joint_momentum(model, &position, &next_velocity)?;
        let mut residual = 0.0_f64;
        for i in 0..n {
            let tolerance = finite(
                self.controls.absolute_momentum_tolerance[i]
                    + self.controls.relative_momentum_tolerance
                        * momentum[i].abs().max(recovered[i].abs()),
            )?;
            residual = residual.max(finite((momentum[i] - recovered[i]).abs() / tolerance)?);
        }
        if residual > 1.0 {
            return Err(Error::Budget("conjugate momentum transfer residual"));
        }
        let after = energy(model, self.base, &position, &next_velocity, self.gravity)?;
        let parameter = after.change_from(after_motion)?;
        self.work.applied = finite(self.work.applied + applied)?;
        self.work.parameter = finite(self.work.parameter + parameter)?;
        self.work.absolute_parameter = finite(self.work.absolute_parameter + parameter.abs())?;
        self.work.numerical = finite(self.work.numerical + numerical)?;
        self.work.absolute_numerical = finite(self.work.absolute_numerical + numerical.abs())?;
        self.work.momentum_residual = self.work.momentum_residual.max(residual);
        self.work.attempts = self
            .work
            .attempts
            .checked_add(report.attempts)
            .ok_or(Error::Budget("mechanical attempt count overflow"))?;
        if self.work.absolute_parameter > self.controls.maximum_parameter_work_j {
            return Err(Error::Budget("absolute thermal-configuration work"));
        }
        if self.work.absolute_numerical > self.controls.maximum_numerical_energy_change_j {
            return Err(Error::Budget(
                "absolute fixed-configuration numerical energy change",
            ));
        }
        checkpoint()?;
        self.position = position;
        self.velocity = next_velocity;
        self.time_s = thermal.time_s();
        Ok(())
    }
}
fn energy(
    model: &ArticulatedModel,
    base: BaseState,
    q: &[f64],
    velocity: &[f64],
    gravity: Vec3,
) -> Result<MechanicalEnergy, Error> {
    let kinematics = forward_kinematics(model, base, q, velocity)?;
    let mut kinetic = 0.0;
    let mut potential = 0.0;
    for ((link, &pose), &twist) in model
        .links()
        .iter()
        .zip(&kinematics.world_from_link)
        .zip(&kinematics.body_twist)
    {
        let inertia = link.inertia();
        kinetic = finite(kinetic + inertia.kinetic_energy(twist)?)?;
        let center = pose
            .transform_point(inertia.center_of_mass())
            .map_err(ArticulatedError::from)?;
        potential = finite(potential - inertia.mass() * gravity.dot(center))?;
    }
    Ok(MechanicalEnergy {
        kinetic_j: kinetic,
        gravitational_j: potential,
    })
}
fn finite(x: f64) -> Result<f64, Error> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(Error::Invalid(
            "nonfinite mechanical state or energy account",
        ))
    }
}
fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}
fn nonnegative(x: f64) -> bool {
    x.is_finite() && x >= 0.0
}
