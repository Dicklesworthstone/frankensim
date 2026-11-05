//! Work-conjugate normal ports from the original conservative source owners.
//!
//! Signed gap is retained: a positive separation can carry JKR attachment or
//! finite-range DMT/Maugis attraction. One affine increment uses the source's
//! original mean effort and endpoint derivative, not endpoint force times dt.
//! Both wrenches use a common fixed spatial origin. Integrated infinitesimal
//! body twists must reproduce the declared gap change in that fixed frame.
//! This adapter does not discover contact, rotate a moving frame, select a snap,
//! or supply omitted Hunt--Crossley loss. Source branches and budgets remain
//! authoritative. Its own stepping state and arithmetic allocate no memory.

use super::{NormalLawSource, NormalLawSourceError};
use fs_blake3::ContentHash;
use fs_exec::CancelGate;
use fs_math::dd::Dd;

/// A fixed normal coordinate and one common spatial moment origin.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SourceNormalPortGeometry {
    /// Unit normal from the reacting body toward the acting body.
    pub normal: [f64; 3],
    /// Vector from the COMMON spatial moment origin to the contact point.
    pub moment_arm_m: [f64; 3],
}

/// Explicit numerical/operating bounds; these do not extend source validity.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SourceNormalPortControls {
    pub maximum_step_s: f64,
    pub maximum_absolute_approach_rate_m_per_s: f64,
    /// Bounds the published mean and endpoint forces, not an interior peak.
    pub maximum_force_n: f64,
    pub maximum_absolute_energy_j: f64,
    pub kinematic_absolute_tolerance_m: f64,
    pub work_absolute_tolerance_j: f64,
    pub work_relative_tolerance: f64,
    pub maximum_steps: u64,
}

/// Integrated body motion in the declared FIXED linearized work frame.
/// Rotation is an integrated infinitesimal twist, not a finite rotation angle.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct NormalPortBodyIncrement {
    /// Translational component of the integrated spatial twist at the SAME
    /// common moment origin as the wrenches; not an arbitrary COM displacement.
    pub translation_m: [f64; 3],
    pub rotation_rad: [f64; 3],
}

/// One converged candidate interval. Absolute endpoint time prevents replay.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SourceNormalPortStep {
    pub end_time_s: f64,
    /// Signed geometric separation; approach is its negative.
    pub gap_m: f64,
    pub action: NormalPortBodyIncrement,
    pub reaction: NormalPortBodyIncrement,
    pub converged: bool,
}

/// Mean action/reaction wrenches, about the SAME fixed spatial origin.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SourceNormalWrenches {
    pub action_force_n: [f64; 3],
    pub action_moment_n_m: [f64; 3],
    pub reaction_force_n: [f64; 3],
    pub reaction_moment_n_m: [f64; 3],
}

/// Accepted work and original source observations for one interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SourceNormalPortReport {
    pub source_identity: ContentHash,
    pub step: u64,
    pub from_time_s: f64,
    pub to_time_s: f64,
    pub from_approach_m: f64,
    pub to_approach_m: f64,
    pub mean_force_n: f64,
    pub endpoint_force_n: f64,
    /// Derivative of mean effort with respect to endpoint approach.
    /// Retains the source's one-sided representative at coincident piecewise
    /// corners; it is not advertised as a smooth instantaneous tangent there.
    pub mean_force_end_derivative_n_per_m: f64,
    pub wrenches: SourceNormalWrenches,
    pub stored_energy_before_j: f64,
    pub stored_energy_after_j: f64,
    /// Original mean effort times approach increment; positive into storage.
    pub source_work_j: f64,
    /// Actual work of BOTH wrenches on the supplied body increments.
    pub work_on_bodies_j: f64,
    /// Endpoint potential change minus original signed source work.
    pub source_work_residual_j: f64,
    /// Body wrench work plus source work; separately gated against cancellation.
    pub projection_work_residual_j: f64,
    pub source_work_tolerance_j: f64,
    pub projection_work_tolerance_j: f64,
}

/// Refusals leave accepted source, clock, coordinate and work count unchanged.
#[derive(Clone, Debug, PartialEq)]
pub enum SourceNormalPortError {
    Invalid(&'static str),
    Source(NormalLawSourceError),
    Cancelled,
    Nonconverged,
    StaleTime,
    StepBudget,
    Kinematics {
        residual_m: f64,
        tolerance_m: f64,
    },
    Work {
        channel: &'static str,
        residual_j: f64,
        tolerance_j: f64,
    },
}
impl core::fmt::Display for SourceNormalPortError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "source normal port refused: {self:?}")
    }
}
impl std::error::Error for SourceNormalPortError {}
impl From<NormalLawSourceError> for SourceNormalPortError {
    fn from(error: NormalLawSourceError) -> Self {
        Self::Source(error)
    }
}

/// Conservative source, fixed work geometry, and accepted endpoint only.
/// A receiver stages its candidate in the admission callback and installs it
/// infallibly after success; neither owner may publish from inside the callback.
#[derive(Clone, Debug, PartialEq)]
pub struct SourceNormalPathPort {
    source: NormalLawSource,
    geometry: SourceNormalPortGeometry,
    controls: SourceNormalPortControls,
    time_s: f64,
    approach_m: f64,
    stored_energy_j: f64,
    steps: u64,
}

impl SourceNormalPathPort {
    pub fn new(
        source: NormalLawSource,
        geometry: SourceNormalPortGeometry,
        controls: SourceNormalPortControls,
        initial_gap_m: f64,
        initial_time_s: f64,
        gate: &CancelGate,
    ) -> Result<Self, SourceNormalPortError> {
        poll(gate)?;
        validate_controls(controls)?;
        for value in geometry.normal.into_iter().chain(geometry.moment_arm_m) {
            finite(value)?;
        }
        let norm = dot(geometry.normal, geometry.normal).to_f64();
        if !norm.is_finite() || (norm - 1.0).abs() > 1e-12 {
            return Err(invalid("finite unit normal"));
        }
        if !initial_time_s.is_finite() || initial_time_s < 0.0 {
            return Err(invalid("finite nonnegative initial time"));
        }
        if let NormalLawSource::HertzPoint(point) = &source {
            if point.dissipation_s_per_m() != 0.0 {
                return Err(invalid(
                    "Hunt-Crossley loss requires its separate physical owner",
                ));
            }
        }
        let approach_m = -finite(initial_gap_m)?;
        let (stored_energy_j, force) = endpoint(&source, approach_m)?;
        bounded(
            stored_energy_j,
            controls.maximum_absolute_energy_j,
            "initial source energy",
        )?;
        bounded(force, controls.maximum_force_n, "initial source force")?;
        poll(gate)?;
        Ok(Self {
            source,
            geometry,
            controls,
            time_s: initial_time_s,
            approach_m,
            stored_energy_j,
            steps: 0,
        })
    }

    pub const fn source(&self) -> &NormalLawSource {
        &self.source
    }
    pub const fn geometry(&self) -> SourceNormalPortGeometry {
        self.geometry
    }
    pub const fn controls(&self) -> SourceNormalPortControls {
        self.controls
    }
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    pub const fn approach_m(&self) -> f64 {
        self.approach_m
    }
    pub const fn stored_energy_j(&self) -> f64 {
        self.stored_energy_j
    }
    pub const fn accepted_steps(&self) -> u64 {
        self.steps
    }

    pub fn advance(
        &mut self,
        step: SourceNormalPortStep,
        gate: &CancelGate,
    ) -> Result<SourceNormalPortReport, SourceNormalPortError> {
        self.advance_with_admission(step, gate, |report| Ok(*report))
    }

    /// Evaluate the original source and actual body work before a receiver
    /// stages its candidate. All accepted fields publish together after the
    /// final cancellation poll; receiver refusal leaves this owner untouched.
    pub fn advance_with_admission<R, E>(
        &mut self,
        step: SourceNormalPortStep,
        gate: &CancelGate,
        admit: impl FnOnce(&SourceNormalPortReport) -> Result<R, E>,
    ) -> Result<R, E>
    where
        E: From<SourceNormalPortError>,
    {
        let report = self.prepare(step, gate).map_err(E::from)?;
        let receiving = admit(&report)?;
        poll(gate).map_err(E::from)?;
        self.time_s = report.to_time_s;
        self.approach_m = report.to_approach_m;
        self.stored_energy_j = report.stored_energy_after_j;
        self.steps = report.step;
        Ok(receiving)
    }

    fn prepare(
        &self,
        step: SourceNormalPortStep,
        gate: &CancelGate,
    ) -> Result<SourceNormalPortReport, SourceNormalPortError> {
        poll(gate)?;
        if !step.converged {
            return Err(SourceNormalPortError::Nonconverged);
        }
        if !step.end_time_s.is_finite() {
            return Err(invalid("finite endpoint time"));
        }
        if step.end_time_s <= self.time_s {
            return Err(SourceNormalPortError::StaleTime);
        }
        let count = self
            .steps
            .checked_add(1)
            .ok_or(SourceNormalPortError::StepBudget)?;
        if count > self.controls.maximum_steps {
            return Err(SourceNormalPortError::StepBudget);
        }
        let dt = finite(step.end_time_s - self.time_s)?;
        if dt > self.controls.maximum_step_s {
            return Err(invalid("maximum step duration"));
        }
        let to = -finite(step.gap_m)?;
        let displacement = finite(to - self.approach_m)?;
        let rate = finite(displacement / dt)?;
        bounded(
            rate,
            self.controls.maximum_absolute_approach_rate_m_per_s,
            "maximum approach rate",
        )?;
        for body in [step.action, step.reaction] {
            for value in body.translation_m.into_iter().chain(body.rotation_rad) {
                finite(value)?;
            }
        }
        let relative_translation =
            difference(step.action.translation_m, step.reaction.translation_m);
        let relative_rotation = difference(step.action.rotation_rad, step.reaction.rotation_rad);
        let arm = self.geometry.moment_arm_m.map(Dd::from_f64);
        let rotational_motion = cross(relative_rotation, arm);
        let normal = self.geometry.normal.map(Dd::from_f64);
        let mut normal_motion = Dd::from_f64(0.0);
        for i in 0..3 {
            normal_motion =
                normal_motion + normal[i] * (relative_translation[i] + rotational_motion[i]);
        }
        let gap_residual = finite((normal_motion + Dd::from_f64(displacement)).to_f64())?;
        if gap_residual.abs() > self.controls.kinematic_absolute_tolerance_m {
            return Err(SourceNormalPortError::Kinematics {
                residual_m: gap_residual,
                tolerance_m: self.controls.kinematic_absolute_tolerance_m,
            });
        }
        // Preserve the Hertz owner's actual rate/time admission even though
        // this port is purely conservative. Signed inactive path is retained.
        if let NormalLawSource::HertzPoint(point) = &self.source {
            for approach in [self.approach_m, to] {
                point
                    .validate_kinematics(approach.max(0.0), rate, dt)
                    .map_err(NormalLawSourceError::Point)?;
            }
        }
        let (after, end_force) = endpoint(&self.source, to)?;
        bounded(
            after,
            self.controls.maximum_absolute_energy_j,
            "source endpoint energy",
        )?;
        bounded(
            end_force,
            self.controls.maximum_force_n,
            "source endpoint force",
        )?;
        let mean = finite(
            self.source
                .mean_force_between(self.approach_m, to)?
                .ok_or_else(|| invalid("original source supplies no path effort"))?,
        )?;
        bounded(mean, self.controls.maximum_force_n, "source mean force")?;
        let derivative = finite(
            self.source
                .mean_force_end_derivative(self.approach_m, to)?
                .ok_or_else(|| invalid("original source supplies no path derivative"))?,
        )?;
        poll(gate)?;
        let force = self.geometry.normal.map(|n| n * mean);
        let moment = cross(arm, force.map(Dd::from_f64)).map(Dd::to_f64);
        for value in force.into_iter().chain(moment) {
            finite(value)?;
        }
        let source_work = Dd::from_f64(mean) * Dd::from_f64(displacement);
        let potential_change = Dd::from_f64(after) - Dd::from_f64(self.stored_energy_j);
        let source_residual = finite((potential_change - source_work).to_f64())?;
        let source_work_j = finite(source_work.to_f64())?;
        let source_tolerance = self.work_gate(
            "original source path",
            source_residual,
            finite(potential_change.to_f64())?,
            source_work_j,
        )?;
        // Evaluate body work from the actual published wrenches, not by
        // assigning it the expected negative potential change.
        let mut body_work = Dd::from_f64(0.0);
        for i in 0..3 {
            body_work = body_work
                + Dd::from_f64(force[i]) * relative_translation[i]
                + Dd::from_f64(moment[i]) * relative_rotation[i];
        }
        let body_work_j = finite(body_work.to_f64())?;
        let projection_residual = finite((body_work + source_work).to_f64())?;
        let projection_tolerance = self.work_gate(
            "action/reaction projection",
            projection_residual,
            body_work_j,
            source_work_j,
        )?;
        Ok(SourceNormalPortReport {
            source_identity: self.source.identity(),
            step: count,
            from_time_s: self.time_s,
            to_time_s: step.end_time_s,
            from_approach_m: self.approach_m,
            to_approach_m: to,
            mean_force_n: mean,
            endpoint_force_n: end_force,
            mean_force_end_derivative_n_per_m: derivative,
            wrenches: SourceNormalWrenches {
                action_force_n: force,
                action_moment_n_m: moment,
                reaction_force_n: force.map(|f| -f),
                reaction_moment_n_m: moment.map(|m| -m),
            },
            stored_energy_before_j: self.stored_energy_j,
            stored_energy_after_j: after,
            source_work_j,
            work_on_bodies_j: body_work_j,
            source_work_residual_j: source_residual,
            projection_work_residual_j: projection_residual,
            source_work_tolerance_j: source_tolerance,
            projection_work_tolerance_j: projection_tolerance,
        })
    }

    fn work_gate(
        &self,
        channel: &'static str,
        residual: f64,
        a: f64,
        b: f64,
    ) -> Result<f64, SourceNormalPortError> {
        let tolerance = finite(
            self.controls.work_absolute_tolerance_j
                + self.controls.work_relative_tolerance * a.abs().max(b.abs()),
        )?;
        if residual.abs() > tolerance {
            return Err(SourceNormalPortError::Work {
                channel,
                residual_j: residual,
                tolerance_j: tolerance,
            });
        }
        Ok(tolerance)
    }
}

fn endpoint(source: &NormalLawSource, approach: f64) -> Result<(f64, f64), SourceNormalPortError> {
    let query = if source.supports_signed_approach() {
        approach
    } else {
        approach.max(0.0)
    };
    let (u, f) = source.energy_and_force(query)?;
    Ok((finite(u)?, finite(f)?))
}
fn validate_controls(c: SourceNormalPortControls) -> Result<(), SourceNormalPortError> {
    for value in [
        c.maximum_step_s,
        c.maximum_force_n,
        c.maximum_absolute_energy_j,
    ] {
        if !value.is_finite() || value <= 0.0 {
            return Err(invalid("positive finite port bound"));
        }
    }
    for value in [
        c.maximum_absolute_approach_rate_m_per_s,
        c.kinematic_absolute_tolerance_m,
        c.work_absolute_tolerance_j,
        c.work_relative_tolerance,
    ] {
        if !value.is_finite() || value < 0.0 {
            return Err(invalid("nonnegative finite port tolerance"));
        }
    }
    if c.maximum_steps == 0 {
        return Err(SourceNormalPortError::StepBudget);
    }
    Ok(())
}
fn invalid(message: &'static str) -> SourceNormalPortError {
    SourceNormalPortError::Invalid(message)
}
fn finite(value: f64) -> Result<f64, SourceNormalPortError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid("unrepresentable normal port value"))
    }
}
fn bounded(value: f64, maximum: f64, message: &'static str) -> Result<(), SourceNormalPortError> {
    if value.abs() <= maximum {
        Ok(())
    } else {
        Err(invalid(message))
    }
}
fn poll(gate: &CancelGate) -> Result<(), SourceNormalPortError> {
    if gate.is_requested() {
        Err(SourceNormalPortError::Cancelled)
    } else {
        Ok(())
    }
}
fn difference(a: [f64; 3], b: [f64; 3]) -> [Dd; 3] {
    core::array::from_fn(|i| Dd::from_f64(a[i]) - Dd::from_f64(b[i]))
}
fn dot(a: [f64; 3], b: [f64; 3]) -> Dd {
    let mut sum = Dd::from_f64(0.0);
    for i in 0..3 {
        sum = sum + Dd::from_f64(a[i]) * Dd::from_f64(b[i]);
    }
    sum
}
fn cross(a: [Dd; 3], b: [Dd; 3]) -> [Dd; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
