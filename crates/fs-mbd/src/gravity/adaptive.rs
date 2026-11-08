//! Bounded adaptive intervals composed from the original native gravity steps.

use super::{GravityError, finite, poll};
use crate::{
    DynamicsDiagnostics, DynamicsError, MassProperties, RigidBodyState, UnitQuaternion, Vec3,
};
use fs_math::dd::Dd;

/// Physical error tolerances and explicit work/step limits for one interval.
///
/// Every trial compares one native step of length h with two of length h/2.
/// The accepted state is the two-half-step state, without extrapolation. The
/// difference divided by three estimates its local error for the second-order
/// native maps; it is not a certified or accumulated global error bound.
///
/// Absolute tolerances are positive SI values. The nonnegative relative term
/// scales with each body's local change during the trial. Position uses the
/// displacement relative to its initial ballistic drift, never the absolute
/// world origin or a common velocity boost. Vector errors use Euclidean norms, orientation uses the
/// shortest SO(3) angle, and angular velocity is compared in world coordinates.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GravityAdvanceOptions {
    /// First trial duration [s]. Must lie between the declared limits.
    pub initial_step_seconds: f64,
    /// Smallest ordinary trial duration [s], before division into two halves.
    /// A shorter final remainder may succeed, but cannot be refined further.
    pub minimum_step_seconds: f64,
    /// Largest trial duration [s], before division into two halves.
    pub maximum_step_seconds: f64,
    /// Positive cap on all trials, including rejected and domain-refined trials.
    /// At most three native steps are invoked per trial. Must fit usize/3.
    pub maximum_trials: usize,
    /// Absolute COM position error tolerance [m].
    pub absolute_position_tolerance_m: f64,
    /// Absolute linear velocity error tolerance [m/s].
    pub absolute_velocity_tolerance_m_s: f64,
    /// Absolute orientation error tolerance [rad].
    pub absolute_orientation_tolerance_rad: f64,
    /// Absolute world angular velocity error tolerance [rad/s].
    pub absolute_angular_velocity_tolerance_rad_s: f64,
    /// Nonnegative relative tolerance against local trial changes.
    pub relative_tolerance: f64,
}

impl GravityAdvanceOptions {
    pub(super) fn validate(self, duration_seconds: f64) -> Result<(), GravityAdvanceError> {
        if !duration_seconds.is_finite() || duration_seconds <= 0.0 {
            return Err(DynamicsError::InvalidStepDuration.into());
        }
        if !self.minimum_step_seconds.is_finite()
            || !self.initial_step_seconds.is_finite()
            || !self.maximum_step_seconds.is_finite()
            || self.minimum_step_seconds <= 0.0
            || self.initial_step_seconds < self.minimum_step_seconds
            || self.maximum_step_seconds < self.initial_step_seconds
        {
            return Err(GravityAdvanceError::InvalidOptions(
                "step limits must satisfy 0 < minimum <= initial <= maximum",
            ));
        }
        if self.maximum_trials == 0 || self.maximum_trials > usize::MAX / 3 {
            return Err(GravityAdvanceError::InvalidOptions(
                "trial budget must lie in 1..=usize::MAX/3",
            ));
        }
        if [
            self.absolute_position_tolerance_m,
            self.absolute_velocity_tolerance_m_s,
            self.absolute_orientation_tolerance_rad,
            self.absolute_angular_velocity_tolerance_rad_s,
        ]
        .iter()
        .any(|&value| !value.is_finite() || value <= 0.0)
            || !self.relative_tolerance.is_finite()
            || self.relative_tolerance < 0.0
        {
            return Err(GravityAdvanceError::InvalidOptions(
                "absolute tolerances must be positive and relative tolerance nonnegative",
            ));
        }
        Ok(())
    }
}

/// Whole-interval refusal; no intermediate accepted state is published.
#[derive(Clone, Debug, PartialEq)]
pub enum GravityAdvanceError {
    /// Malformed tolerances, timestep limits or trial budget.
    InvalidOptions(&'static str),
    /// The original gravity/native owner refused, including cancellation.
    Gravity(GravityError),
    /// The interval did not finish within the declared trial budget.
    TrialLimit,
    /// The error estimate still exceeds tolerance at the minimum duration.
    MinimumStepReached,
    /// A forward interval or its two equal positive halves is unrepresentable.
    TimeUnderflow,
}

impl core::fmt::Display for GravityAdvanceError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(formatter, "adaptive gravity: {self:?}")
    }
}
impl std::error::Error for GravityAdvanceError {}
impl From<GravityError> for GravityAdvanceError {
    fn from(error: GravityError) -> Self {
        Self::Gravity(error)
    }
}
impl From<DynamicsError> for GravityAdvanceError {
    fn from(error: DynamicsError) -> Self {
        Self::Gravity(error.into())
    }
}

/// Work and local-error evidence for one completed interval.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct GravityAdvanceStats {
    /// All attempted full-step/two-half-step comparisons.
    pub trials: usize,
    /// Accepted trials; each contributes two native half steps.
    pub accepted_trials: usize,
    /// Error or typed expansion-domain refusals that were refined.
    pub rejected_trials: usize,
    /// Native step invocations, including discarded coarse and refused steps.
    pub native_steps: usize,
    /// Largest normalized local error estimate among accepted trials, at most 1.
    pub maximum_accepted_error_ratio: f64,
}

/// A completed interval retaining original states, diagnostics and direct work.
///
/// Mutual systems use native state vectors; a one-body field adapter may retain
/// an inline array. Variable step selection is not claimed symplectic or time
/// reversible. No energy/momentum projection or quaternion extrapolation occurs.
#[derive(Clone, Debug, PartialEq)]
pub struct GravityAdvanceReceipt<S = Vec<RigidBodyState>> {
    /// Complete requested duration [s].
    pub duration_seconds: f64,
    /// Original native body states after the entire interval.
    pub states_after: S,
    /// Original diagnostics before the first accepted half step.
    pub diagnostics_before: DynamicsDiagnostics,
    /// Original diagnostics after the last accepted half step.
    pub diagnostics_after: DynamicsDiagnostics,
    /// Compensated sum of accepted native potential increments [J]. Rejected
    /// trials contribute nothing; absolute endpoint potentials are not subtracted.
    pub potential_energy_change_j: f64,
    /// Bounded work counts and accepted local error estimates.
    pub stats: GravityAdvanceStats,
    /// Proposed initial trial duration for a subsequent call [s]. This is not
    /// a resumable controller checkpoint or a split-run equivalence guarantee.
    pub suggested_step_seconds: f64,
}

// A narrow seam for the two existing native storage shapes, not an ODE state
// encoding. The callback still owns all force, spin, domain and energy physics.
pub(super) struct NativeStep<S> {
    pub states_after: S,
    pub diagnostics_before: DynamicsDiagnostics,
    pub diagnostics_after: DynamicsDiagnostics,
    pub potential_energy_change_j: f64,
}

impl From<super::mutual::MutualGravityStepReceipt> for NativeStep<Vec<RigidBodyState>> {
    fn from(step: super::mutual::MutualGravityStepReceipt) -> Self {
        Self {
            states_after: step.states_after,
            diagnostics_before: step.diagnostics_before,
            diagnostics_after: step.diagnostics_after,
            potential_energy_change_j: step.potential_energy_change_j,
        }
    }
}

/// Adapters prevalidate their initial load before entering this private driver.
pub(super) fn advance_native<S: AsRef<[RigidBodyState]>>(
    mut states: S,
    properties: &[MassProperties],
    duration_seconds: f64,
    options: GravityAdvanceOptions,
    mut step: impl FnMut(
        &[RigidBodyState],
        f64,
        &mut dyn FnMut() -> bool,
    ) -> Result<NativeStep<S>, GravityError>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<GravityAdvanceReceipt<S>, GravityAdvanceError> {
    options.validate(duration_seconds)?;
    poll(&mut cancelled)?;
    if states.as_ref().len() != properties.len() {
        return Err(GravityError::Invalid("adaptive gravity body/property count mismatch").into());
    }
    let mut elapsed = 0.0;
    let mut proposal = options.initial_step_seconds;
    let mut stats = GravityAdvanceStats::default();
    let mut before = None;
    let mut after = None;
    let mut potential_change = Dd::ZERO;
    while elapsed < duration_seconds {
        poll(&mut cancelled)?;
        if stats.trials == options.maximum_trials {
            return Err(GravityAdvanceError::TrialLimit);
        }
        let remaining = duration_seconds - elapsed;
        let requested = proposal.min(remaining);
        let final_trial = requested == remaining;
        let mut end = if final_trial {
            duration_seconds
        } else {
            elapsed + requested
        };
        // Integrate the represented interval, never a larger nominal duration.
        if !final_trial && end - elapsed > requested {
            end = end.next_down();
        }
        let h = end - elapsed;
        let half = 0.5 * h;
        if !end.is_finite() || end <= elapsed || half <= 0.0 || half + half != h {
            return Err(GravityAdvanceError::TimeUnderflow);
        }
        if !final_trial && h < options.minimum_step_seconds {
            return Err(GravityAdvanceError::TimeUnderflow);
        }
        stats.trials += 1;
        let trial: Result<(NativeStep<S>, f64, f64), GravityError> = (|| {
            poll(&mut cancelled)?;
            stats.native_steps += 1;
            let coarse = step(states.as_ref(), h, &mut cancelled)?;
            poll(&mut cancelled)?;
            stats.native_steps += 1;
            let first = step(states.as_ref(), half, &mut cancelled)?;
            poll(&mut cancelled)?;
            stats.native_steps += 1;
            let mut fine = step(first.states_after.as_ref(), half, &mut cancelled)?;
            let error = error_ratio(
                states.as_ref(),
                coarse.states_after.as_ref(),
                fine.states_after.as_ref(),
                properties,
                h,
                options,
                &mut cancelled,
            )?;
            fine.diagnostics_before = first.diagnostics_before;
            Ok((fine, first.potential_energy_change_j, error))
        })();
        let (fine, first_change, error) = match trial {
            Ok(trial) => trial,
            Err(GravityError::OutsideExpansionDomain) if h > options.minimum_step_seconds => {
                // Initial loads were already admitted. This typed refusal can
                // be caused by the attempted path, so bounded refinement is safe.
                stats.rejected_trials += 1;
                proposal = (0.5 * h).max(options.minimum_step_seconds);
                continue;
            }
            Err(error) => return Err(error.into()),
        };
        if error > 1.0 {
            stats.rejected_trials += 1;
            if h <= options.minimum_step_seconds {
                return Err(GravityAdvanceError::MinimumStepReached);
            }
            proposal = (h * controller_factor(error).min(0.5)).max(options.minimum_step_seconds);
            continue;
        }

        potential_change = potential_change
            + Dd::from_f64(first_change)
            + Dd::from_f64(fine.potential_energy_change_j);
        finite(potential_change.to_f64())?;
        poll(&mut cancelled)?;
        before.get_or_insert(fine.diagnostics_before);
        after = Some(fine.diagnostics_after);
        states = fine.states_after;
        elapsed = end;
        stats.accepted_trials += 1;
        stats.maximum_accepted_error_ratio = stats.maximum_accepted_error_ratio.max(error);
        // A successful clipped tail must not replace the useful next proposal.
        if !final_trial || requested == proposal {
            proposal = bounded_proposal(h, controller_factor(error), options);
        }
    }
    poll(&mut cancelled)?;
    Ok(GravityAdvanceReceipt {
        duration_seconds,
        states_after: states,
        diagnostics_before: before.ok_or(GravityAdvanceError::TimeUnderflow)?,
        diagnostics_after: after.ok_or(GravityAdvanceError::TimeUnderflow)?,
        potential_energy_change_j: finite(potential_change.to_f64())?,
        stats,
        suggested_step_seconds: proposal,
    })
}

fn controller_factor(error: f64) -> f64 {
    if error == 0.0 {
        2.0
    } else if error.is_infinite() {
        0.1
    } else {
        (0.9 * fs_math::det::pow(error, -1.0 / 3.0)).clamp(0.1, 2.0)
    }
}

fn bounded_proposal(h: f64, factor: f64, options: GravityAdvanceOptions) -> f64 {
    if factor > 1.0 && h >= options.maximum_step_seconds / factor {
        options.maximum_step_seconds
    } else {
        (h * factor).clamp(options.minimum_step_seconds, options.maximum_step_seconds)
    }
}

fn error_ratio(
    initial: &[RigidBodyState],
    coarse: &[RigidBodyState],
    fine: &[RigidBodyState],
    properties: &[MassProperties],
    duration: f64,
    options: GravityAdvanceOptions,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, GravityError> {
    if coarse.len() != initial.len() || fine.len() != initial.len() {
        return Err(GravityError::Invalid(
            "adaptive gravity native state count changed",
        ));
    }
    let mut error = 0.0_f64;
    for (((&start, &coarse), &fine), &mass) in initial.iter().zip(coarse).zip(fine).zip(properties)
    {
        poll(cancelled)?;
        let cp = coarse.pose().position_world();
        let fp = fine.pose().position_world();
        let position_error = norm(fp.sub(cp))?;
        let velocity_error = norm(velocity_change(fine, coarse, mass.mass())?)?;
        let orientation_error = angle(coarse.pose().orientation(), fine.pose().orientation());
        let cw = world_angular_velocity(coarse, mass)?;
        let fw = world_angular_velocity(fine, mass)?;
        let spin_error = norm(fw.sub(cw))?;
        let mut scales = [0.0; 4];
        if options.relative_tolerance > 0.0 {
            scales[0] = deflection(coarse, start, mass.mass(), duration)?.max(deflection(
                fine,
                start,
                mass.mass(),
                duration,
            )?);
            scales[1] = norm(velocity_change(coarse, start, mass.mass())?)?
                .max(norm(velocity_change(fine, start, mass.mass())?)?);
            scales[2] = angle(start.pose().orientation(), coarse.pose().orientation())
                .max(angle(start.pose().orientation(), fine.pose().orientation()));
            let start_spin = world_angular_velocity(start, mass)?;
            scales[3] = norm(cw.sub(start_spin))?.max(norm(fw.sub(start_spin))?);
        }
        for ((difference, absolute), scale) in [
            (position_error, options.absolute_position_tolerance_m),
            (velocity_error, options.absolute_velocity_tolerance_m_s),
            (
                orientation_error,
                options.absolute_orientation_tolerance_rad,
            ),
            (
                spin_error,
                options.absolute_angular_velocity_tolerance_rad_s,
            ),
        ]
        .into_iter()
        .zip(scales)
        {
            let tolerance = finite(absolute + options.relative_tolerance * scale)?;
            // Divide by the physical tolerance first: rounding a subnormal
            // difference/3 can otherwise turn a genuine excess into acceptance.
            let ratio = (difference / tolerance) / 3.0;
            // An infinite dimensionless ratio still means an ordinary rejection;
            // a non-finite physical difference/scale already refused above.
            error = error.max(ratio);
        }
    }
    Ok(error)
}

fn norm(vector: Vec3) -> Result<f64, GravityError> {
    Ok(vector.stable_norm("adaptive gravity state difference")?)
}

fn deflection(
    after: RigidBodyState,
    before: RigidBodyState,
    mass: f64,
    duration: f64,
) -> Result<f64, GravityError> {
    let a = after.pose().position_world();
    let b = before.pose().position_world();
    let p = before.linear_momentum_world();
    let component = |a: f64, b: f64, p: f64| {
        finite(
            (Dd::from_f64(a)
                - Dd::from_f64(b)
                - (Dd::from_f64(p) / Dd::from_f64(mass)) * Dd::from_f64(duration))
            .to_f64(),
        )
    };
    norm(Vec3::new(
        component(a.x, b.x, p.x)?,
        component(a.y, b.y, p.y)?,
        component(a.z, b.z, p.z)?,
    ))
}

fn velocity_change(
    after: RigidBodyState,
    before: RigidBodyState,
    mass: f64,
) -> Result<Vec3, GravityError> {
    let a = after.linear_momentum_world();
    let b = before.linear_momentum_world();
    let component = |a: f64, b: f64| {
        let difference = a - b;
        finite(if difference.is_finite() {
            difference / mass
        } else {
            a / mass - b / mass
        })
    };
    Ok(Vec3::new(
        component(a.x, b.x)?,
        component(a.y, b.y)?,
        component(a.z, b.z)?,
    ))
}

fn world_angular_velocity(
    state: RigidBodyState,
    mass: MassProperties,
) -> Result<Vec3, GravityError> {
    Ok(state
        .pose()
        .orientation()
        .rotate_body_to_world(mass.angular_velocity_body_checked(state.angular_momentum_body())?))
}

// Quaternion sum/difference chords avoid acos(dot)'s small-angle cancellation
// and retain q/-q equivalence. Scaling also preserves subnormal chord lengths.
fn angle(a: UnitQuaternion, b: UnitQuaternion) -> f64 {
    let a = a.components();
    let b = b.components();
    let chord = |sign: f64| {
        let values = core::array::from_fn::<_, 4, _>(|i| a[i] + sign * b[i]);
        let scale = values.iter().fold(0.0_f64, |scale, v| scale.max(v.abs()));
        if scale == 0.0 {
            return 0.0;
        }
        scale
            * values
                .iter()
                .map(|v| (v / scale).powi(2))
                .sum::<f64>()
                .sqrt()
    };
    let difference = chord(-1.0);
    let sum = chord(1.0);
    4.0 * fs_math::det::atan2(difference.min(sum), difference.max(sum))
}
