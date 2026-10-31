//! Source heat and continuing articulated motion on one accepted timeline.
//!
//! Each thermal substep advances native ABA dynamics through fs-time RK45 at
//! its starting material configuration. The endpoint's source-resolved inertias
//! and material-bound joint placements then receive the same conjugate joint
//! momenta at unchanged joint positions.
//! This is a first-order split with a stationary prescribed base. Parameter
//! work is reported separately and never deposited as caloric heat.
//! Optional source-backed journal films apply resisting torques in the SAME
//! mechanical stages and deposit their integrated work on this thermal clock.
//! Ordered material dry brakes use that same work/heat path for both bodies.
//! Source-backed squeeze films share those stages through their full native
//! generalized pressure-force map, rather than an independent damping fit.
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

/// Material-driven journal loss and its receiving thermal control volume.
pub mod viscous;
use viscous::{BoundBearing, ExpandingJournalBearing, JournalBearing, JournalBearingState};
use viscous::rheology::RheologicalJournalBearing;

/// Source-bound dry sliding, generalized resisting effort and two-body heating.
pub mod dry;
use dry::{BoundBrake, DryBrakeState, DryJointBrake};
use fs_tribo::DissipationStep;

/// Material Reynolds pressure, generalized resistance and finite-wall heating.
pub mod squeeze;
use squeeze::{BoundSqueeze, SqueezeFilm, SqueezeFilmState};

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
    Material(fs_material::state_point::MaterialStatePointError),
    ThermalMaterial(fs_material::thermal::source::ThermalMaterialError),
    Interface(fs_contact::interface_binding::InterfaceBindingError),
    Tribo(fs_tribo::TriboError),
    Film(fs_tribo::resistive_film::FilmError),
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
impl From<fs_material::state_point::MaterialStatePointError> for Error {
    fn from(e: fs_material::state_point::MaterialStatePointError) -> Self {
        Self::Material(e)
    }
}
impl From<DepositionError> for Error {
    fn from(e: DepositionError) -> Self {
        Self::Heating(e.into())
    }
}
impl From<fs_material::thermal::source::ThermalMaterialError> for Error {
    fn from(e: fs_material::thermal::source::ThermalMaterialError) -> Self {
        Self::ThermalMaterial(e)
    }
}
impl From<fs_contact::interface_binding::InterfaceBindingError> for Error {
    fn from(e: fs_contact::interface_binding::InterfaceBindingError) -> Self {
        Self::Interface(e)
    }
}
impl From<fs_tribo::TriboError> for Error {
    fn from(e: fs_tribo::TriboError) -> Self {
        Self::Tribo(e)
    }
}
impl From<fs_tribo::resistive_film::FilmError> for Error {
    fn from(e: fs_tribo::resistive_film::FilmError) -> Self {
        Self::Film(e)
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
    /// Journal and squeeze-film loss deposited once in `thermal.external_heat_j` [J].
    /// Subtract this internal transfer when identifying externally supplied heat.
    pub viscous_heat_j: f64,
    /// Dry sliding loss deposited once into the two ordered finite bodies [J].
    /// This is included in thermal.external_heat_j, an internal transfer here.
    pub friction_heat_j: f64,
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
/// base. Joint types, axes and thermal transport geometry remain prescribed.
/// Free expansion changes each link's retained inertia and explicitly bound
/// material attachment positions. The original momentum transfer accounts for
/// the full changed mass matrix and gravitational potential as parameter work.
/// This is a sequence of frozen configurations: continuous expansion velocities,
/// their inertial forces, general contacts, closed kinematic loops and reciprocal
/// thermoelastic heat require other owners. The original
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
            bearings: Vec::new(),
            brakes: Vec::new(),
            squeeze: Vec::new(),
            pending: None,
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

    /// Attach up to 64 massless, isothermal journal films to this mechanism.
    /// Each explicitly named receiving cell supplies the lubricant temperature
    /// and receives all of that film's dissipated work. Film storage, pressure
    /// support and changes to the prescribed bearing geometry are not inferred.
    /// Existing accepted state is unchanged on every refusal or cancellation.
    pub fn with_journal_bearings(
        &self,
        bearings: Vec<JournalBearing<'a>>,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, Error> {
        checkpoint()?;
        if !self.motion.bearings.is_empty()
            || bearings.is_empty()
            || bearings.len() + self.motion.brakes.len() + self.motion.squeeze.len() > 64
        {
            return Err(Error::Invalid("attach 1..=64 journal films exactly once"));
        }
        let mut candidate = self.clone();
        for bearing in bearings {
            checkpoint()?;
            candidate.motion.bearings.push(BoundBearing::new(
                bearing,
                self.configuration().model(),
                &self.thermal,
                &self.motion.velocity,
            )?);
        }
        checkpoint()?;
        Ok(candidate)
    }

    /// Attach concentric films whose clearance follows the SAME source-owned
    /// region expansion as the material mechanism's inertia. Both surfaces are
    /// centered at their revolute joint; the shorter expanded cylinder sets the
    /// overlapping length. Source bindings and original reference geometry stay
    /// immutable, while actual geometry is retained in each bearing state.
    ///
    /// The existing mechanical stages freeze viscosity and geometry over one
    /// capped thermal step. The accepted endpoint refreshes both together; a
    /// closed gap, excessive geometry change or source refusal rolls back the
    /// complete coupled step. This is quasisteady circumferential shear, not
    /// squeeze-film pressure/work, end leakage or a constrained expansion solve.
    /// Radial-motion effects must be negligible in the caller's regime. Link
    /// inertia parameter work remains separate from the film's dissipated heat.
    pub fn with_expanding_journal_bearings(
        &self,
        bearings: Vec<ExpandingJournalBearing<'a>>,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, Error> {
        checkpoint()?;
        if !self.motion.bearings.is_empty()
            || bearings.is_empty()
            || bearings.len() + self.motion.brakes.len() + self.motion.squeeze.len() > 64
        {
            return Err(Error::Invalid("attach 1..=64 journal films exactly once"));
        }
        let mut candidate = self.clone();
        for bearing in bearings {
            checkpoint()?;
            candidate.motion.bearings.push(BoundBearing::new_expanding(
                bearing,
                self.configuration(),
                &self.thermal,
                &self.motion.velocity,
            )?);
        }
        checkpoint()?;
        Ok(candidate)
    }

    /// Attach source-bound, finite-gap generalized-Newtonian journals. The
    /// original T/shear line drives actual nonlinear torque at every mechanical
    /// stage, while its conjugate loss uses the SAME native RK coordinate and
    /// finite receiving cell as Newtonian journals. Source charts and optional
    /// source-owned geometry refresh together at accepted thermal endpoints.
    /// The massless isothermal/creeping-flow and 64-interface limits remain.
    pub fn with_rheological_journal_bearings(
        &self,
        bearings: Vec<RheologicalJournalBearing<'a>>,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, Error> {
        checkpoint()?;
        if !self.motion.bearings.is_empty() || bearings.is_empty()
            || bearings.len() + self.motion.brakes.len() + self.motion.squeeze.len() > 64
        {
            return Err(Error::Invalid("attach 1..=64 journal films exactly once"));
        }
        let mut candidate = self.clone();
        for bearing in bearings {
            checkpoint()?;
            candidate.motion.bearings.push(BoundBearing::new_rheological(
                bearing, self.configuration(), &self.thermal,
                &self.motion.velocity, &mut checkpoint,
            )?);
        }
        checkpoint()?;
        Ok(candidate)
    }

    /// Actual endpoint source states, resistance, torque and power in input order.
    pub fn journal_bearings(&self) -> impl ExactSizeIterator<Item = &JournalBearingState> {
        self.motion.bearings.iter().map(|b| &b.state)
    }

    /// Attach ordered child/parent dry brakes, sharing the native mechanical
    /// stages and thermal substeps with any journal films. At most 64 total
    /// dissipative interfaces are admitted. Each brake's two finite cells must
    /// be actual regions of its child and direct parent; source order matters.
    /// Normal preload/geometry/heat fractions are explicit prescribed inputs.
    pub fn with_dry_joint_brakes(
        &self,
        brakes: Vec<DryJointBrake<'a>>,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, Error> {
        checkpoint()?;
        if !self.motion.brakes.is_empty()
            || brakes.is_empty()
            || brakes.len() + self.motion.bearings.len() + self.motion.squeeze.len() > 64
        {
            return Err(Error::Invalid(
                "attach dry brakes once within the 64-interface total budget",
            ));
        }
        let mut candidate = self.clone();
        for brake in brakes {
            checkpoint()?;
            candidate.motion.brakes.push(BoundBrake::new(
                brake,
                self.configuration().model(),
                &self.motion.expansion.regions,
                &self.thermal,
                &self.motion.velocity,
                &mut checkpoint,
            )?);
        }
        checkpoint()?;
        Ok(candidate)
    }

    /// Accepted source temperatures, friction response and cumulative heat split.
    pub fn dry_joint_brakes(&self) -> impl ExactSizeIterator<Item = &DryBrakeState> {
        self.motion.brakes.iter().map(|b| &b.state)
    }

    /// Attach native Reynolds graphs in this mechanism's complete joint order.
    /// Each graph contributes its full generalized pressure resistance and
    /// integrates the conjugate power through the SAME adaptive mechanical
    /// stages. Its named finite wall receives the heat and supplies the next
    /// source temperature. Native geometry, low absolute-pressure liquid
    /// support, source and energy limits all participate in atomic acceptance.
    ///
    /// This is the controlled nonstiff ODE image, not the pressure core's usual
    /// implicit stiff-film embedding. At most 64 total journal/brake/squeeze
    /// interfaces are admitted. No film inertia, thermal store, pressure clamp,
    /// cavitation content or thermal change to the authored gap map is inferred.
    pub fn with_squeeze_films(
        &self,
        films: Vec<SqueezeFilm<'a>>,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, Error> {
        checkpoint()?;
        if !self.motion.squeeze.is_empty()
            || films.is_empty()
            || films.len() + self.motion.bearings.len() + self.motion.brakes.len() > 64
        {
            return Err(Error::Invalid("attach squeeze films once within the 64-interface total budget"));
        }
        let mut candidate = self.clone();
        for film in films {
            checkpoint()?;
            candidate.motion.squeeze.push(BoundSqueeze::new(
                film, &self.thermal, &self.motion.position, &self.motion.velocity, &mut checkpoint,
            )?);
        }
        checkpoint()?;
        Ok(candidate)
    }

    /// Actual source, pressure, generalized resistance and accepted heat history.
    pub fn squeeze_films(&self) -> impl ExactSizeIterator<Item = &SqueezeFilmState> {
        self.motion.squeeze.iter().map(|film| &film.state)
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
            viscous_heat_j: w.viscous,
            friction_heat_j: w.friction,
            parameter_work_j: w.parameter,
            absolute_parameter_work_j: w.absolute_parameter,
            numerical_energy_change_j: w.numerical,
            absolute_numerical_energy_change_j: w.absolute_numerical,
            maximum_scaled_momentum_residual: w.momentum_residual,
            rk_attempts: w.attempts,
            combined_residual_j: finite(
                thermal.residual_j + after.change_from(before)? - w.applied
                    + w.viscous
                    + w.friction
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
    viscous: f64,
    friction: f64,
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
    bearings: Vec<BoundBearing<'a>>,
    brakes: Vec<BoundBrake<'a>>,
    squeeze: Vec<BoundSqueeze<'a>>,
    pending: Option<IntegratedMotion>,
}
#[derive(Clone)]
struct IntegratedMotion {
    end_s: f64,
    position: Vec<f64>,
    velocity: Vec<f64>,
    after_motion: MechanicalEnergy,
    applied: f64,
    viscous: f64,
    friction: f64,
    dry_steps: Vec<DissipationStep>,
    squeeze_losses: Vec<f64>,
    heat_w: BTreeMap<CellId, f64>,
    numerical: f64,
    attempts: usize,
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
    fn integrate(
        &self,
        start_s: f64,
        dt_s: f64,
        end_s: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<IntegratedMotion, Error> {
        checkpoint()?;
        if start_s != self.time_s || !positive(dt_s) || end_s <= start_s {
            return Err(Error::Invalid("joint and thermal clocks disagree"));
        }
        let model = self.expansion.configuration.model();
        let before = self.energy()?;
        let n = self.position.len();
        let mut scales: Vec<_> = self
            .controls
            .absolute_position_tolerance
            .iter()
            .chain(&self.controls.absolute_velocity_tolerance)
            .copied()
            .collect();
        scales.extend(
            self.bearings
                .iter()
                .map(|b| b.definition.limits.absolute_loss_tolerance_j),
        );
        scales.extend(
            self.brakes
                .iter()
                .map(|b| b.definition.limits.absolute_loss_tolerance_j),
        );
        scales.extend(
            self.squeeze.iter().map(|film| film.definition.limits.absolute_loss_tolerance_j),
        );
        let mut normalized: Vec<_> = self
            .position
            .iter()
            .chain(&self.velocity)
            .zip(&scales)
            .map(|(&u, &scale)| u / scale)
            .collect();
        normalized.resize(scales.len(), 0.0);
        let mut state = AdaptiveState::new(start_s, &normalized, dt_s);
        let external = vec![Wrench::default(); model.link_count()];
        let source_error = RefCell::new(None);
        let rhs = |_: f64, u: &[f64], out: &mut [f64]| {
            let result = (|| -> Result<(), Error> {
                let physical: Vec<_> = u.iter().zip(&scales).map(|(&x, &s)| x * s).collect();
                let mut effort = self.effort.clone();
                for (j, bearing) in self.bearings.iter().enumerate() {
                    let (torque, power, _) = bearing.response(physical[n + bearing.dof])?;
                    effort[bearing.dof] = finite(effort[bearing.dof] + torque)?;
                    out[2 * n + j] = power / scales[2 * n + j];
                }
                for (j, brake) in self.brakes.iter().enumerate() {
                    let (force, power) = brake.response(physical[n + brake.dof])?;
                    effort[brake.dof] = finite(effort[brake.dof] + force)?;
                    let index = 2 * n + self.bearings.len() + j;
                    out[index] = power / scales[index];
                }
                for (j, film) in self.squeeze.iter().enumerate() {
                    let (resistance, _, response) = film.response(&physical[..n], &physical[n..2 * n])?;
                    for (effort, force) in effort.iter_mut().zip(resistance) {
                        *effort = finite(*effort - force)?;
                    }
                    let index = 2 * n + self.bearings.len() + self.brakes.len() + j;
                    out[index] = response.dissipated_power_w / scales[index];
                }
                let response = forward_dynamics(
                    model,
                    self.base,
                    &physical[..n],
                    &physical[n..2 * n],
                    &effort,
                    self.gravity,
                    &external,
                )?;
                for i in 0..n {
                    out[i] = physical[n + i] / scales[i];
                    out[n + i] = response.generalized_acceleration[i] / scales[n + i];
                }
                Ok(())
            })();
            if let Err(error) = result {
                *source_error.borrow_mut() = Some(error);
                out.fill(f64::NAN);
            }
        };
        let mut checkpoint_error = None;
        let result = rk45_adaptive_controlled(
            &mut state,
            &rhs,
            end_s,
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
            return Err(error);
        }
        let report = result?;
        if report.status != AdaptiveStatus::ReachedEnd {
            return Err(Error::Budget("mechanical RK45 attempts exhausted"));
        }
        let physical: Vec<_> = state.u.iter().zip(&scales).map(|(&x, &s)| x * s).collect();
        let position = physical[..n].to_vec();
        let velocity = physical[n..2 * n].to_vec();
        let after_motion = energy(model, self.base, &position, &velocity, self.gravity)?;
        let applied = finite(
            self.effort
                .iter()
                .zip(&position)
                .zip(&self.position)
                .map(|((&force, &q), &old_q)| force * (q - old_q))
                .sum(),
        )?;
        let mut viscous = 0.0;
        let mut heat_w = BTreeMap::new();
        for (bearing, &loss) in self.bearings.iter().zip(&physical[2 * n..]) {
            if !nonnegative(loss) {
                return Err(Error::Budget(
                    "journal integrated loss must be nonnegative; refine RK tolerances",
                ));
            }
            // Validate the fixed-coefficient endpoint, including the last RK stage.
            bearing.response(velocity[bearing.dof])?;
            viscous = finite(viscous + loss)?;
            let watts = heat_w.entry(bearing.definition.heat_cell).or_insert(0.0);
            *watts = finite(*watts + loss / dt_s)?;
        }
        let mut friction = 0.0;
        let mut dry_steps = Vec::with_capacity(self.brakes.len());
        for (brake, &loss) in self
            .brakes
            .iter()
            .zip(&physical[2 * n + self.bearings.len()..])
        {
            brake.response(velocity[brake.dof])?;
            let step = brake.partition_loss(loss, dt_s)?;
            friction = finite(friction + step.total_work_j())?;
            for (cell, heat) in brake
                .definition
                .cells
                .into_iter()
                .zip([step.surface_a_heat_j(), step.surface_b_heat_j()])
            {
                let watts = heat_w.entry(cell).or_insert(0.0);
                *watts = finite(*watts + heat / dt_s)?;
            }
            dry_steps.push(step);
        }
        let mut squeeze_losses = Vec::with_capacity(self.squeeze.len());
        for (film, &loss) in self.squeeze.iter().zip(
            &physical[2 * n + self.bearings.len() + self.brakes.len()..],
        ) {
            if !nonnegative(loss) {
                return Err(Error::Budget("squeeze integrated loss must be nonnegative; refine RK tolerances"));
            }
            film.response(&position, &velocity)?;
            viscous = finite(viscous + loss)?;
            let watts = heat_w.entry(film.definition.heat_cell).or_insert(0.0);
            *watts = finite(*watts + loss / dt_s)?;
            squeeze_losses.push(loss);
        }
        let numerical = finite(after_motion.change_from(before)? - applied + viscous + friction)?;
        checkpoint()?;
        Ok(IntegratedMotion {
            end_s,
            position,
            velocity,
            after_motion,
            applied,
            viscous,
            friction,
            dry_steps,
            squeeze_losses,
            heat_w,
            numerical,
            attempts: report.attempts,
        })
    }
}
impl DepositionObserver for Motion<'_> {
    type Error = Error;
    fn begin(&mut self) -> Result<(), Error> {
        self.work = Work::default();
        self.pending = None;
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
    fn requires_step_interval(&self) -> bool {
        !self.bearings.is_empty() || !self.brakes.is_empty() || !self.squeeze.is_empty()
    }
    fn step_sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        start_s: f64,
        dt_s: f64,
        end_s: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Error> {
        checkpoint()?;
        if thermal.time_s() != start_s || self.pending.is_some() {
            return Err(Error::Invalid("unmatched mechanical interval"));
        }
        // These are the exact accepted start-state material snapshots. They
        // cannot be replaced independently of this thermal/motion history.
        let candidate = self.integrate(start_s, dt_s, end_s, checkpoint)?;
        let heat = candidate.heat_w.clone();
        self.pending = Some(candidate);
        Ok(heat)
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
        let candidate = self
            .pending
            .take()
            .ok_or(Error::Invalid("missing mechanical interval"))?;
        if candidate.end_s != thermal.time_s() {
            return Err(Error::Invalid(
                "mechanical candidate endpoint differs from thermal endpoint",
            ));
        }
        let IntegratedMotion {
            position,
            velocity,
            after_motion,
            applied,
            viscous,
            friction,
            dry_steps,
            squeeze_losses,
            numerical,
            attempts,
            ..
        } = candidate;
        let model = self.expansion.configuration.model();
        let n = self.position.len();
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
        for bearing in &mut self.bearings {
            checkpoint()?;
            bearing.refresh(thermal, &next_velocity, &self.expansion.configuration, checkpoint)?;
        }
        for (brake, loss) in self.brakes.iter_mut().zip(dry_steps) {
            brake.accept(loss, thermal, &next_velocity, checkpoint)?;
        }
        for (film, loss) in self.squeeze.iter_mut().zip(squeeze_losses) {
            film.accept(loss, thermal, &position, &next_velocity, checkpoint)?;
        }
        self.work.applied = finite(self.work.applied + applied)?;
        self.work.viscous = finite(self.work.viscous + viscous)?;
        self.work.friction = finite(self.work.friction + friction)?;
        self.work.parameter = finite(self.work.parameter + parameter)?;
        self.work.absolute_parameter = finite(self.work.absolute_parameter + parameter.abs())?;
        self.work.numerical = finite(self.work.numerical + numerical)?;
        self.work.absolute_numerical = finite(self.work.absolute_numerical + numerical.abs())?;
        self.work.momentum_residual = self.work.momentum_residual.max(residual);
        self.work.attempts = self
            .work
            .attempts
            .checked_add(attempts)
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
