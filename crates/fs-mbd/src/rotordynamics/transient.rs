//! Time-domain rotor motion at one declared constant shaft speed.
//!
//! Implicit midpoint advances the ORIGINAL full equation
//! `M v' + (C + Omega G) v + K q = f`, with `q' = v`. Native real
//! factorizations are cached: preparation costs O(n^3), subsequent steps and
//! instantaneous impulse responses cost O(n^2), with n <= 32. The step solves
//! for a velocity increment, avoiding subtraction of a large unchanged M v.
//!
//! Skew gyroscopic force changes motion without physical dissipation. Bearing
//! C and K retain all cross-coupled and active entries; their signed work is
//! not called heat. For symmetric K and C=0, midpoint preserves the quadratic
//! mechanical energy up to solve/rounding error; its phase still needs timestep
//! refinement. A caller-supplied midpoint force is second-order only when it
//! approximates the actual smooth load there. Spin ramps, nonlinear bearings,
//! shaft-model validity and event timing remain outside this constant-speed lane.

use super::{MAX_ROTOR_DOFS, RotorBearingModel, RotorError, RotorOperatingPoint, finite, poll};
use fs_la::factor::{Cholesky, FactorError, Lu, cholesky, lu};
use fs_math::dd::Dd;
use std::borrow::Cow;

pub mod drive;

/// Original generalized displacement/angle and rate coordinates.
#[derive(Debug, Clone, PartialEq)]
pub struct RotorTransientState {
    /// Original displacements/angles, using the model's coordinate units.
    pub displacement: Vec<f64>,
    /// Their time derivatives, conjugate to the supplied forces/torques.
    pub velocity: Vec<f64>,
}

/// A complete accepted finite step, including work of the original operators.
#[derive(Debug, Clone, PartialEq)]
pub struct RotorTransientStep {
    /// Reusable physical state; the borrowed input is never modified.
    pub state_after: RotorTransientState,
    /// Relative residual of the ORIGINAL impulse balance after diagonal mass
    /// scaling; zero-force rows retain the complete vector's roundoff scale.
    /// Only the part beyond the endpoint's rounding counts (see [`RotorMidpoint::step`]).
    pub relative_residual: f64,
    /// Signed generalized load work at the actual midpoint.
    pub input_work_j: f64,
    /// Signed work against C v. Active bearings can make this negative.
    pub velocity_force_work_j: f64,
    /// Signed work against K q. Only symmetric K defines elastic storage.
    pub stiffness_force_work_j: f64,
    /// Measured work against Omega G v; zero in exact arithmetic.
    pub gyroscopic_work_j: f64,
    /// Polarized kinetic change from the actual returned velocity increment.
    pub kinetic_energy_change_j: f64,
    /// Kinetic change plus original operator work minus supplied work.
    pub work_energy_residual_j: f64,
}

/// Mechanical power available to one commanded generalized-effort pattern.
///
/// These are signed NET interval-average limits for the complete pattern, not
/// per-coordinate or instantaneous limits. A source may absorb and return work
/// inside the interval. Electrical efficiency, finite source storage and peak
/// power require their own models and a resolved time grid.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RotorPowerControl {
    /// Maximum nonnegative average mechanical power supplied to the rotor [W].
    pub maximum_motoring_power_w: f64,
    /// Maximum nonnegative magnitude of average power recovered from it [W].
    pub maximum_regenerative_power_w: f64,
    /// Explicit absolute allowance [W] for final power-limit excess and the
    /// original step's independently evaluated work/energy residual divided
    /// by its duration. Zero requests exact represented closure.
    pub absolute_power_tolerance_w: f64,
}

/// Accepted motion with the original commanded effort uniformly curtailed.
#[derive(Debug, Clone, PartialEq)]
pub struct RotorPowerStep {
    /// Original full-coordinate trajectory and independent operator work.
    pub motion: RotorTransientStep,
    /// Applied force/torque is this fraction in [0,1] times the supplied command.
    /// Individual effort magnitudes can never exceed their authored command.
    pub effort_fraction: f64,
    /// Actual original source work divided by the accepted duration [W].
    pub average_input_power_w: f64,
}

/// Mechanical power limits throughout the accepted linear-velocity segment.
///
/// Unlike [`RotorPowerControl`], both ends of the segment must satisfy these
/// limits separately. They constrain the total commanded effort pattern, not
/// individual coordinates. The interpolation is the midpoint method's local
/// reconstruction; physical peaks between resolved steps require refinement.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RotorPeakPowerControl {
    /// Maximum nonnegative instantaneous mechanical input [W].
    pub maximum_motoring_power_w: f64,
    /// Maximum nonnegative magnitude of instantaneous recovered power [W].
    pub maximum_regenerative_power_w: f64,
    /// Absolute allowance for endpoint power excess and native work residual
    /// per second [W]. Zero requests exact represented closure.
    pub absolute_power_tolerance_w: f64,
}

/// Accepted peak-limited actuation with independently evaluated endpoint power.
#[derive(Debug, Clone, PartialEq)]
pub struct RotorPeakPowerStep {
    /// Original accepted motion, effort fraction and signed net average power.
    pub average: RotorPowerStep,
    /// Applied effort paired with the original velocity [W].
    pub initial_input_power_w: f64,
    /// Applied effort paired with the returned velocity [W].
    pub final_input_power_w: f64,
}

/// A failed power-controlled step publishes no motion or source debit.
#[derive(Debug, Clone, PartialEq)]
pub enum RotorPowerError {
    /// Original model, state, factorization, momentum or cancellation refusal.
    Dynamics(RotorError),
    /// Power inputs or their scalar response cannot be represented faithfully.
    Invalid(&'static str),
    /// Final original work violates the explicitly authored power allowance.
    PowerBalance {
        /// Actual signed input work divided by duration [W].
        average_input_power_w: f64,
        /// Nonnegative amount by which either signed limit is exceeded [W].
        limit_excess_w: f64,
        /// Signed native kinetic/operator/source work residual per second [W].
        work_residual_w: f64,
        /// Caller-authored absolute admission allowance [W].
        tolerance_w: f64,
    },
    /// Rounded accepted motion violates a peak power or work allowance.
    PeakPowerBalance {
        /// Applied effort paired with the initial velocity [W].
        initial_input_power_w: f64,
        /// Applied effort paired with the accepted final velocity [W].
        final_input_power_w: f64,
        /// Largest nonnegative violation of either signed endpoint limit [W].
        limit_excess_w: f64,
        /// Signed native work/energy residual per second [W].
        work_residual_w: f64,
        /// Caller-authored absolute allowance [W].
        tolerance_w: f64,
    },
}

impl From<RotorError> for RotorPowerError {
    fn from(error: RotorError) -> Self {
        Self::Dynamics(error)
    }
}

impl core::fmt::Display for RotorPowerError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "rotor power-controlled motion refused: {self:?}")
    }
}

impl std::error::Error for RotorPowerError {}

/// Instantaneous physical generalized impulse, with unchanged configuration.
#[derive(Debug, Clone, PartialEq)]
pub struct RotorImpulseResponse {
    /// Same displacements with the accepted full-inertia velocity response.
    pub state_after: RotorTransientState,
    /// Original mass-scaled momentum-increment residual, after endpoint rounding.
    pub relative_residual: f64,
    /// Physical impulse paired with the actual midpoint velocity.
    pub impulse_work_j: f64,
    /// Polarized kinetic change; includes the quadratic impulse energy.
    pub kinetic_energy_change_j: f64,
}

/// Reusable fixed-duration stepper borrowing one immutable physical model.
#[derive(Debug, Clone)]
pub struct RotorMidpoint<'a> {
    model: &'a RotorBearingModel,
    spin_rad_s: f64,
    duration_seconds: f64,
    relative_tolerance: f64,
    velocity_operator: Vec<f64>,
    stiffness_operator: Cow<'a, [f64]>,
    mass_scale: Vec<f64>,
    mass_factor: Cholesky,
    step_row_scale: Vec<f64>,
    step_factor: Lu,
}

impl<'a> RotorOperatingPoint<'a> {
    /// Prepare a midpoint trajectory and physical-impulse operator once.
    ///
    /// M is scaled by its original diagonal, and the complete nonsymmetric
    /// step matrix is row-equilibrated before native partial-pivot LU. No
    /// matrix is symmetrized, diagonalized, inverted or regularized. A singular
    /// step refuses. `relative_tolerance` in (0,1) admits the actual original
    /// momentum equation after each solve, including endpoint rounding.
    /// Native dense factorizations/solves are indivisible but bounded to 32
    /// coordinates; cancellation is polled around them and per assembly row.
    pub fn midpoint_stepper(
        &self,
        duration_seconds: f64,
        relative_tolerance: f64,
        cancelled: impl FnMut() -> bool,
    ) -> Result<RotorMidpoint<'a>, RotorError> {
        self.midpoint_stepper_with_stiffness(
            duration_seconds,
            relative_tolerance,
            Cow::Borrowed(&self.model.stiffness),
            cancelled,
        )
    }

    // The physical rigid-rotor ramp supplies its explicit Euler operator.
    // Keeping this internal prevents arbitrary operators from silently
    // replacing the declared bearing stiffness in the steady public lane.
    pub(super) fn midpoint_stepper_with_stiffness(
        &self,
        duration_seconds: f64,
        relative_tolerance: f64,
        stiffness_operator: Cow<'a, [f64]>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RotorMidpoint<'a>, RotorError> {
        poll(&mut cancelled)?;
        let half = 0.5 * duration_seconds;
        if !duration_seconds.is_finite()
            || half <= 0.0
            || !relative_tolerance.is_finite()
            || relative_tolerance <= 0.0
            || relative_tolerance >= 1.0
        {
            return Err(RotorError::Invalid(
                "positive representable duration and relative tolerance in (0,1) required",
            ));
        }
        let m = self.model;
        if stiffness_operator.len() != m.n * m.n
            || stiffness_operator.iter().any(|x| !x.is_finite())
        {
            return Err(RotorError::Invalid(
                "finite matching rotor stiffness required",
            ));
        }
        let mass_scale: Vec<_> = (0..m.n)
            .map(|i| fs_math::det::sqrt(m.mass[i * m.n + i]))
            .collect();
        let mut scaled_mass = vec![0.0; m.n * m.n];
        let mut matrix = vec![0.0; m.n * m.n];
        let mut step_row_scale = vec![0.0_f64; m.n];
        for i in 0..m.n {
            poll(&mut cancelled)?;
            for j in 0..m.n {
                let ij = i * m.n + j;
                scaled_mass[ij] = finite((m.mass[ij] / mass_scale[i]) / mass_scale[j])?;
                let entry = finite(half.mul_add(
                    half * stiffness_operator[ij],
                    half.mul_add(self.velocity_operator[ij], m.mass[ij]),
                ))?;
                matrix[ij] = finite((entry / mass_scale[i]) / mass_scale[j])?;
                step_row_scale[i] = step_row_scale[i].max(matrix[ij].abs());
            }
            if step_row_scale[i] == 0.0 {
                return Err(RotorError::TransientLinear(FactorError::Singular {
                    index: i,
                }));
            }
            for entry in &mut matrix[i * m.n..(i + 1) * m.n] {
                *entry /= step_row_scale[i];
            }
        }
        let mass_factor = cholesky(&scaled_mass, m.n).map_err(RotorError::TransientLinear)?;
        poll(&mut cancelled)?;
        let step_factor = lu(&matrix, m.n).map_err(RotorError::TransientLinear)?;
        finite(step_factor.growth)?;
        poll(&mut cancelled)?;
        Ok(RotorMidpoint {
            model: m,
            spin_rad_s: self.spin_rad_s,
            duration_seconds,
            relative_tolerance,
            velocity_operator: self.velocity_operator.clone(),
            stiffness_operator,
            mass_scale,
            mass_factor,
            step_row_scale,
            step_factor,
        })
    }
}

impl RotorMidpoint<'_> {
    /// Fixed physical duration of each accepted trajectory step.
    pub const fn duration_seconds(&self) -> f64 {
        self.duration_seconds
    }

    /// Advance the actual rotor coordinates under a force/torque sampled at
    /// the temporal midpoint. All returned work uses actual endpoint midpoints.
    /// No input is changed if the solve, original balance or cancellation refuses.
    ///
    /// The balance is admitted with its endpoint's rounding: the returned
    /// velocities are rounded, so the change the solve made is known only to
    /// their last places, and each force term only to the rounding of its
    /// products. Per row, 16 eps of `sum_j |M_ij| (|v0_j| + |v1_j|)` plus the
    /// same of each term's absolute products and of the input is the
    /// representation's, not the solve's; only what exceeds it is measured
    /// against `relative_tolerance`. A step whose change is small beside the
    /// velocity it changes is therefore not refused for that alone.
    pub fn step(
        &self,
        state: &RotorTransientState,
        midpoint_force: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RotorTransientStep, RotorError> {
        poll(&mut cancelled)?;
        self.admit_state(state)?;
        self.admit_vector(midpoint_force)?;
        let n = self.model.n;
        let h = self.duration_seconds;
        let increment = self.velocity_increment(state, midpoint_force, &mut cancelled)?;
        let velocity: Vec<_> = state
            .velocity
            .iter()
            .enumerate()
            .map(|(i, &v)| v + increment[i])
            .collect();
        let midpoint_velocity = midpoint(&state.velocity, &velocity);
        let displacement = state
            .displacement
            .iter()
            .zip(&midpoint_velocity)
            .map(|(&q, &v)| h.mul_add(v, q))
            .collect();
        let state_after = RotorTransientState {
            displacement,
            velocity,
        };
        self.admit_state(&state_after)?;
        let midpoint_position = midpoint(&state.displacement, &state_after.displacement);
        let mut delta = [0.0; MAX_ROTOR_DOFS];
        for (i, value) in delta.iter_mut().enumerate().take(n) {
            *value = state_after.velocity[i] - state.velocity[i];
        }
        let inertial = action(&self.model.mass, &delta[..n], &mut cancelled)?;
        let damping = scaled_action(
            &self.model.damping,
            &midpoint_velocity[..n],
            h,
            &mut cancelled,
        )?;
        let gyro = scaled_action(
            &self.model.gyro,
            &midpoint_velocity[..n],
            finite(h * self.spin_rad_s)?,
            &mut cancelled,
        )?;
        let stiffness = scaled_action(
            &self.stiffness_operator,
            &midpoint_position[..n],
            h,
            &mut cancelled,
        )?;
        let mut input = [0.0; MAX_ROTOR_DOFS];
        for i in 0..n {
            input[i] = h * midpoint_force[i];
        }
        self.admit_vector(&input[..n])?;
        let floor = self.endpoint_floor(
            [&state.velocity, &state_after.velocity],
            &midpoint_velocity[..n],
            &midpoint_position[..n],
            &input[..n],
            &mut cancelled,
        )?;
        let relative_residual = self.admit_balance(
            &[&inertial[..n], &damping[..n], &gyro[..n], &stiffness[..n]],
            &input[..n],
            &floor[..n],
        )?;
        let input_work_j = pairing(&input[..n], &midpoint_velocity[..n])?;
        let velocity_force_work_j = pairing(&damping[..n], &midpoint_velocity[..n])?;
        let stiffness_force_work_j = pairing(&stiffness[..n], &midpoint_velocity[..n])?;
        let gyroscopic_work_j = pairing(&gyro[..n], &midpoint_velocity[..n])?;
        let kinetic_energy_change_j = pairing(&inertial[..n], &midpoint_velocity[..n])?;
        let work_energy_residual_j = finite(
            (Dd::from_f64(kinetic_energy_change_j)
                + Dd::from_f64(velocity_force_work_j)
                + Dd::from_f64(stiffness_force_work_j)
                + Dd::from_f64(gyroscopic_work_j)
                - Dd::from_f64(input_work_j))
            .to_f64(),
        )?;
        poll(&mut cancelled)?;
        Ok(RotorTransientStep {
            state_after,
            relative_residual,
            input_work_j,
            velocity_force_work_j,
            stiffness_force_work_j,
            gyroscopic_work_j,
            kinetic_energy_change_j,
            work_energy_residual_j,
        })
    }

    /// Advance with finite commanded effort and explicit average-power limits.
    ///
    /// The applied force is `alpha * commanded_force`, with the largest
    /// admissible `alpha` in [0,1] for this original midpoint discretization.
    /// The SAME cached full nonsymmetric factorization supplies the unforced
    /// response and the force influence. Consequently the source's mean power
    /// is `alpha * (b + alpha * a)`, including the acceleration caused by the
    /// source itself. Neither initial-speed clipping nor division by a stalled
    /// shaft speed is used. The command itself supplies the finite effort cap.
    ///
    /// Both roots are considered; active/cross-coupled bearings and reversal
    /// can make source work nonmonotone in `alpha`. These are net interval
    /// limits, so an interval with zero net work can still exchange energy in
    /// both directions. Shorten the duration to resolve instantaneous limits,
    /// within-step reversal, controller changes and source-storage constraints.
    ///
    /// Two bounded response solves select the fraction, then the original
    /// [`Self::step`] independently checks the actual rounded endpoint under
    /// that effort. Its source work, signed bearing work and gyro work remain
    /// authoritative; no additional heat or energy debit is created. The final
    /// mean power and native work residual must meet `control` before return.
    /// Preparation is reused, work stays O(n²), and the extra workspace consists
    /// of fixed stack arrays; no matrix inverse or nonlinear iteration is added.
    ///
    /// # Errors
    /// Nonfinite/negative limits or allowance, unrepresentable response/root,
    /// native dynamics refusal, failed final power/energy admission, or
    /// cancellation. The original state and command are borrowed unchanged.
    pub fn step_power_limited(
        &self,
        state: &RotorTransientState,
        commanded_force: &[f64],
        control: RotorPowerControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RotorPowerStep, RotorPowerError> {
        poll(&mut cancelled)?;
        self.admit_state(state)?;
        self.admit_vector(commanded_force)?;
        if [
            control.maximum_motoring_power_w,
            control.maximum_regenerative_power_w,
            control.absolute_power_tolerance_w,
        ]
        .iter()
        .any(|x| !x.is_finite() || *x < 0.0)
        {
            return Err(RotorPowerError::Invalid(
                "finite nonnegative power limits and absolute allowance required",
            ));
        }

        let n = self.model.n;
        let zero = [0.0; MAX_ROTOR_DOFS];
        let free_increment = self.velocity_increment(state, &zero[..n], &mut cancelled)?;
        let mut force_increment = [0.0; MAX_ROTOR_DOFS];
        for i in 0..n {
            force_increment[i] = ((self.duration_seconds * commanded_force[i])
                / self.mass_scale[i])
                / self.step_row_scale[i];
        }
        self.solve_velocity_increment(&mut force_increment[..n], &mut cancelled)?;

        let mut a = Dd::ZERO;
        let mut b = Dd::ZERO;
        let half = Dd::from_f64(0.5);
        for i in 0..n {
            let force = Dd::from_f64(commanded_force[i]);
            a = a + force * (half * Dd::from_f64(force_increment[i]));
            b = b + force
                * (Dd::from_f64(state.velocity[i]) + half * Dd::from_f64(free_increment[i]));
        }
        for value in [a.hi, a.lo, b.hi, b.lo] {
            finite(value)?;
        }
        let effort_fraction = power_fraction(a, b, control)?;
        let mut applied = [0.0; MAX_ROTOR_DOFS];
        for i in 0..n {
            applied[i] = effort_fraction * commanded_force[i];
        }
        let motion = self.step(state, &applied[..n], &mut cancelled)?;
        let average_input_power_w = finite(motion.input_work_j / self.duration_seconds)?;
        let work_residual_w = finite(motion.work_energy_residual_j / self.duration_seconds)?;
        let limit_excess_w = (average_input_power_w - control.maximum_motoring_power_w)
            .max(-average_input_power_w - control.maximum_regenerative_power_w)
            .max(0.0);
        if limit_excess_w > control.absolute_power_tolerance_w
            || work_residual_w.abs() > control.absolute_power_tolerance_w
        {
            return Err(RotorPowerError::PowerBalance {
                average_input_power_w,
                limit_excess_w,
                work_residual_w,
                tolerance_w: control.absolute_power_tolerance_w,
            });
        }
        poll(&mut cancelled)?;
        Ok(RotorPowerStep {
            motion,
            effort_fraction,
            average_input_power_w,
        })
    }

    /// Curtail a held effort to satisfy both signed endpoint power limits.
    ///
    /// Source power is affine in time along the accepted linear velocity
    /// reconstruction, so its extrema occur at the endpoints even through a
    /// reversal. Initial power is linear in effort fraction and gives a single
    /// interval `[0, cap]`. Final power is quadratic in that fraction. The same
    /// stable two-root solve as [`Self::step_power_limited`] selects its largest
    /// feasible point INSIDE that interval, including active bearing models.
    /// Thus zero net work cannot conceal motoring or regenerative power spikes.
    /// Finite torque at rest is retained without dividing by shaft speed.
    ///
    /// This adds no nonlinear iteration, new factorization or heap workspace:
    /// two cached response solves and one accepted native step remain O(n²).
    /// Independent endpoint pairings and the native work residual are checked
    /// before return. This does not resolve peaks of the underlying continuous
    /// solution beyond the time interpolation, electrical efficiency, finite
    /// energy storage, per-actuator allocation, or changing shaft spin.
    ///
    /// # Errors
    /// Invalid limits or response, native refusal, failed peak/work admission,
    /// or cancellation. Borrowed state and command remain unchanged.
    pub fn step_peak_power_limited(
        &self,
        state: &RotorTransientState,
        commanded_force: &[f64],
        control: RotorPeakPowerControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RotorPeakPowerStep, RotorPowerError> {
        poll(&mut cancelled)?;
        self.admit_state(state)?;
        self.admit_vector(commanded_force)?;
        if [
            control.maximum_motoring_power_w,
            control.maximum_regenerative_power_w,
            control.absolute_power_tolerance_w,
        ]
        .iter()
        .any(|x| !x.is_finite() || *x < 0.0)
        {
            return Err(RotorPowerError::Invalid(
                "finite nonnegative peak power limits and absolute allowance required",
            ));
        }
        let n = self.model.n;
        let zero = [0.0; MAX_ROTOR_DOFS];
        let free_increment = self.velocity_increment(state, &zero[..n], &mut cancelled)?;
        let mut force_increment = [0.0; MAX_ROTOR_DOFS];
        for i in 0..n {
            force_increment[i] = ((self.duration_seconds * commanded_force[i])
                / self.mass_scale[i])
                / self.step_row_scale[i];
        }
        self.solve_velocity_increment(&mut force_increment[..n], &mut cancelled)?;

        let mut initial = Dd::ZERO;
        let mut quadratic = Dd::ZERO;
        let mut linear = Dd::ZERO;
        for i in 0..n {
            let force = Dd::from_f64(commanded_force[i]);
            initial = initial + force * Dd::from_f64(state.velocity[i]);
            quadratic = quadratic + force * Dd::from_f64(force_increment[i]);
            linear = linear
                + force * (Dd::from_f64(state.velocity[i]) + Dd::from_f64(free_increment[i]));
        }
        for value in [initial, quadratic, linear] {
            finite(value.hi)?;
            finite(value.lo)?;
        }
        let limits = RotorPowerControl {
            maximum_motoring_power_w: control.maximum_motoring_power_w,
            maximum_regenerative_power_w: control.maximum_regenerative_power_w,
            absolute_power_tolerance_w: control.absolute_power_tolerance_w,
        };
        let cap = power_fraction(Dd::ZERO, initial, limits)?;
        let cap_dd = Dd::from_f64(cap);
        let fraction = cap * power_fraction(quadratic * cap_dd * cap_dd, linear * cap_dd, limits)?;
        let mut applied = [0.0; MAX_ROTOR_DOFS];
        for i in 0..n {
            applied[i] = fraction * commanded_force[i];
        }
        let motion = self.step(state, &applied[..n], &mut cancelled)?;
        let initial_input_power_w = pairing(&applied[..n], &state.velocity)?;
        let final_input_power_w = pairing(&applied[..n], &motion.state_after.velocity)?;
        let average_input_power_w = finite(motion.input_work_j / self.duration_seconds)?;
        let work_residual_w = finite(motion.work_energy_residual_j / self.duration_seconds)?;
        let limit_excess_w = [initial_input_power_w, final_input_power_w]
            .into_iter()
            .map(|power| {
                (power - control.maximum_motoring_power_w)
                    .max(-power - control.maximum_regenerative_power_w)
            })
            .fold(0.0, f64::max);
        if limit_excess_w > control.absolute_power_tolerance_w
            || work_residual_w.abs() > control.absolute_power_tolerance_w
        {
            return Err(RotorPowerError::PeakPowerBalance {
                initial_input_power_w,
                final_input_power_w,
                limit_excess_w,
                work_residual_w,
                tolerance_w: control.absolute_power_tolerance_w,
            });
        }
        poll(&mut cancelled)?;
        Ok(RotorPeakPowerStep {
            average: RotorPowerStep {
                motion,
                effort_fraction: fraction,
                average_input_power_w,
            },
            initial_input_power_w,
            final_input_power_w,
        })
    }

    /// Apply an instantaneous generalized impulse through the ORIGINAL full M.
    /// Bearings and gyroscopic forces are finite forces and supply no impulsive
    /// correction. Positions remain unchanged. An unrepresentable momentum
    /// increment refuses the original balance instead of silently clipping it.
    pub fn apply_impulse(
        &self,
        state: &RotorTransientState,
        impulse: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RotorImpulseResponse, RotorError> {
        poll(&mut cancelled)?;
        self.admit_state(state)?;
        self.admit_vector(impulse)?;
        let n = self.model.n;
        let mut increment = [0.0; MAX_ROTOR_DOFS];
        for i in 0..n {
            increment[i] = impulse[i] / self.mass_scale[i];
        }
        self.admit_vector(&increment[..n])?;
        self.mass_factor.solve(&mut increment[..n]);
        poll(&mut cancelled)?;
        let state_after = RotorTransientState {
            displacement: state.displacement.clone(),
            velocity: state
                .velocity
                .iter()
                .enumerate()
                .map(|(i, &v)| v + increment[i] / self.mass_scale[i])
                .collect(),
        };
        self.admit_state(&state_after)?;
        let mut delta = [0.0; MAX_ROTOR_DOFS];
        for (i, value) in delta.iter_mut().enumerate().take(n) {
            *value = state_after.velocity[i] - state.velocity[i];
        }
        let momentum = action(&self.model.mass, &delta[..n], &mut cancelled)?;
        // An impulse lost to the velocity's rounding refuses (no floor).
        let relative_residual =
            self.admit_balance(&[&momentum[..n]], impulse, &[0.0; MAX_ROTOR_DOFS][..n])?;
        let midpoint_velocity = midpoint(&state.velocity, &state_after.velocity);
        let impulse_work_j = pairing(impulse, &midpoint_velocity[..n])?;
        let kinetic_energy_change_j = pairing(&momentum[..n], &midpoint_velocity[..n])?;
        poll(&mut cancelled)?;
        Ok(RotorImpulseResponse {
            state_after,
            relative_residual,
            impulse_work_j,
            kinetic_energy_change_j,
        })
    }

    // Each row of `step`'s balance is known only to ENDPOINT_ROUNDING of its
    // absolute products: the mass on both endpoint velocities, each force term
    // on its midpoint, and the input.
    fn endpoint_floor(
        &self,
        endpoints: [&[f64]; 2],
        midpoint_velocity: &[f64],
        midpoint_position: &[f64],
        input: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<[f64; MAX_ROTOR_DOFS], RotorError> {
        let n = self.model.n;
        let h = self.duration_seconds;
        let mut speeds = [0.0; MAX_ROTOR_DOFS];
        for (j, value) in speeds.iter_mut().enumerate().take(n) {
            *value = finite(endpoints[0][j].abs() + endpoints[1][j].abs())?;
        }
        let rounding = [
            absolute_action(&self.model.mass, &speeds[..n], 1.0, cancelled)?,
            absolute_action(&self.model.damping, midpoint_velocity, h, cancelled)?,
            absolute_action(
                &self.model.gyro,
                midpoint_velocity,
                finite(h * self.spin_rad_s)?,
                cancelled,
            )?,
            absolute_action(&self.stiffness_operator, midpoint_position, h, cancelled)?,
        ];
        let mut floor = [0.0; MAX_ROTOR_DOFS];
        for (i, value) in floor.iter_mut().enumerate().take(n) {
            let products = rounding.iter().map(|row| row[i]).sum::<f64>() + input[i].abs();
            *value = finite(ENDPOINT_ROUNDING * products)?;
        }
        Ok(floor)
    }

    fn velocity_increment(
        &self,
        state: &RotorTransientState,
        force: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<[f64; MAX_ROTOR_DOFS], RotorError> {
        let n = self.model.n;
        let h = self.duration_seconds;
        let mut half_position = [0.0; MAX_ROTOR_DOFS];
        for (i, value) in half_position.iter_mut().enumerate().take(n) {
            *value = (0.5 * h).mul_add(state.velocity[i], state.displacement[i]);
        }
        self.admit_vector(&half_position[..n])?;
        let velocity_force = action(&self.velocity_operator, &state.velocity, cancelled)?;
        let spring_force = action(&self.stiffness_operator, &half_position[..n], cancelled)?;
        let mut increment = [0.0; MAX_ROTOR_DOFS];
        for i in 0..n {
            increment[i] = ((h * (force[i] - velocity_force[i] - spring_force[i]))
                / self.mass_scale[i])
                / self.step_row_scale[i];
        }
        self.solve_velocity_increment(&mut increment[..n], cancelled)?;
        Ok(increment)
    }

    fn solve_velocity_increment(
        &self,
        increment: &mut [f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), RotorError> {
        self.admit_vector(increment)?;
        poll(cancelled)?;
        self.step_factor.solve(increment);
        poll(cancelled)?;
        for (value, scale) in increment.iter_mut().zip(&self.mass_scale) {
            *value = finite(*value / scale)?;
        }
        Ok(())
    }

    fn admit_state(&self, state: &RotorTransientState) -> Result<(), RotorError> {
        self.admit_vector(&state.displacement)?;
        self.admit_vector(&state.velocity)
    }

    fn admit_vector(&self, vector: &[f64]) -> Result<(), RotorError> {
        if vector.len() != self.model.n || vector.iter().any(|x| !x.is_finite()) {
            Err(RotorError::Invalid(
                "matching finite original rotor coordinates required",
            ))
        } else {
            Ok(())
        }
    }

    // `floor` is each row's admitted rounding, in force-impulse units; only the
    // residual beyond it is measured against the relative tolerance.
    fn admit_balance(
        &self,
        terms: &[&[f64]],
        rhs: &[f64],
        floor: &[f64],
    ) -> Result<f64, RotorError> {
        let mut residual_norm = 0.0_f64;
        let mut rhs_norm = 0.0_f64;
        let mut norms = [0.0_f64; 4];
        let term_norms = &mut norms[..terms.len()];
        for i in 0..self.model.n {
            let mut residual = Dd::from_f64(-rhs[i]);
            rhs_norm = rhs_norm.max(finite(rhs[i].abs() / self.mass_scale[i])?);
            for (index, term) in terms.iter().enumerate() {
                residual = residual + Dd::from_f64(term[i]);
                term_norms[index] =
                    term_norms[index].max(finite(term[i].abs() / self.mass_scale[i])?);
            }
            let beyond = (residual.to_f64().abs() - floor[i]).max(0.0);
            residual_norm = residual_norm.max(finite(beyond / self.mass_scale[i])?);
        }
        let scale = term_norms.iter().copied().fold(rhs_norm, f64::max);
        let relative = if scale == 0.0 {
            residual_norm
        } else {
            finite(
                (residual_norm / scale)
                    / (rhs_norm / scale + term_norms.iter().map(|x| x / scale).sum::<f64>()),
            )?
        };
        if relative > self.relative_tolerance {
            Err(RotorError::Residual {
                relative,
                limit: self.relative_tolerance,
            })
        } else {
            Ok(relative)
        }
    }
}

// If full command is outside the permitted interval, continuity from P(0)=0
// proves a crossing of that same violated bound. The LAST root in [0,1] is
// the largest feasible fraction, even when the quadratic first crosses the
// opposite bound or when its curvature is negative. A monotone bisection from
// zero would select the wrong branch in those cases.
fn power_fraction(a: Dd, b: Dd, control: RotorPowerControl) -> Result<f64, RotorPowerError> {
    let scale = a.abs().hi.max(b.abs().hi);
    if scale == 0.0 {
        return Ok(1.0);
    }
    let full = a / Dd::from_f64(scale) + b / Dd::from_f64(scale);
    let positive = !full.lt(Dd::ZERO);
    let limit = if positive {
        control.maximum_motoring_power_w
    } else {
        control.maximum_regenerative_power_w
    };
    // A huge allowance cannot invalidate a tiny finite response. Avoid asking
    // double-double division to represent an irrelevant infinite ratio.
    if (limit / scale).is_infinite()
        || !((Dd::from_f64(limit) / Dd::from_f64(scale)).lt(full.abs()))
    {
        return Ok(1.0);
    }
    let target = Dd::from_f64(if positive { limit } else { -limit });

    if b == Dd::ZERO {
        // In particular, start from rest without first forming target/a:
        // that ratio can underflow although its square root and the eventual
        // finite torque, velocity and work are all representable.
        let root = target.abs().sqrt() / a.abs().sqrt();
        if root.hi.is_finite() && root.lo.is_finite() && !root.lt(Dd::ZERO) && !Dd::ONE.lt(root) {
            return Ok(root.to_f64());
        }
        return Err(RotorPowerError::Invalid(
            "pure-quadratic power fraction is unrepresentable",
        ));
    }

    // Normalize by the SELECTED boundary, not a possibly enormous allowance on
    // the other power direction. Retain small roots through the stable q-form.
    let scale = Dd::from_f64(a.abs().hi.max(b.abs().hi).max(target.abs().hi));
    let quadratic = a / scale;
    let linear = b / scale;
    let constant = -target / scale;
    for (original, scaled) in [(a, quadratic), (b, linear), (-target, constant)] {
        if !scaled.hi.is_finite()
            || !scaled.lo.is_finite()
            || original != Dd::ZERO && scaled == Dd::ZERO
        {
            return Err(RotorPowerError::Invalid(
                "power-response coefficient scaling is unrepresentable",
            ));
        }
    }
    let roots = if quadratic == Dd::ZERO {
        [Some(-constant / linear), None]
    } else if constant == Dd::ZERO {
        [Some(-linear / quadratic), Some(Dd::ZERO)]
    } else {
        let discriminant = linear * linear - Dd::from_f64(4.0) * quadratic * constant;
        if discriminant.lt(Dd::ZERO) || !discriminant.hi.is_finite() || !discriminant.lo.is_finite()
        {
            return Err(RotorPowerError::Invalid(
                "no representable real power-boundary crossing",
            ));
        }
        let root = discriminant.sqrt();
        let signed_root = if linear.lt(Dd::ZERO) { -root } else { root };
        let q = Dd::from_f64(-0.5) * (linear + signed_root);
        [Some(q / quadratic), Some(constant / q)]
    };
    let mut selected = None;
    for root in roots.into_iter().flatten() {
        if root.hi.is_finite()
            && root.lo.is_finite()
            && !root.lt(Dd::ZERO)
            && !Dd::ONE.lt(root)
            && selected.is_none_or(|previous: Dd| previous.lt(root))
        {
            selected = Some(root);
        }
    }
    selected.map(Dd::to_f64).ok_or(RotorPowerError::Invalid(
        "no representable power-limited effort fraction in [0,1]",
    ))
}

fn midpoint(before: &[f64], after: &[f64]) -> [f64; MAX_ROTOR_DOFS] {
    let mut result = [0.0; MAX_ROTOR_DOFS];
    for (i, (&a, &b)) in before.iter().zip(after).enumerate() {
        result[i] = 0.5 * a + 0.5 * b;
    }
    result
}

fn action(
    matrix: &[f64],
    vector: &[f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<[f64; MAX_ROTOR_DOFS], RotorError> {
    let n = vector.len();
    let mut result = [0.0_f64; MAX_ROTOR_DOFS];
    for i in 0..n {
        poll(cancelled)?;
        for j in 0..n {
            result[i] = matrix[i * n + j].mul_add(vector[j], result[i]);
        }
        finite(result[i])?;
    }
    Ok(result)
}

fn scaled_action(
    matrix: &[f64],
    vector: &[f64],
    scale: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<[f64; MAX_ROTOR_DOFS], RotorError> {
    let mut result = action(matrix, vector, cancelled)?;
    for value in &mut result[..vector.len()] {
        *value = finite(*value * scale)?;
    }
    Ok(result)
}

/// Rounding allowed per absolute product of a rotor balance row.
const ENDPOINT_ROUNDING: f64 = 16.0 * f64::EPSILON;

// Row sums of |A_ij| |x_j| times |scale|: the size of the products an action adds.
fn absolute_action(
    matrix: &[f64],
    vector: &[f64],
    scale: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<[f64; MAX_ROTOR_DOFS], RotorError> {
    let n = vector.len();
    let mut result = [0.0_f64; MAX_ROTOR_DOFS];
    for i in 0..n {
        poll(cancelled)?;
        for j in 0..n {
            result[i] = matrix[i * n + j].abs().mul_add(vector[j].abs(), result[i]);
        }
        result[i] = finite(result[i] * scale.abs())?;
    }
    Ok(result)
}

fn pairing(force: &[f64], velocity: &[f64]) -> Result<f64, RotorError> {
    let mut sum = Dd::ZERO;
    for (&f, &v) in force.iter().zip(velocity) {
        sum = sum + Dd::from_f64(f) * Dd::from_f64(v);
    }
    finite(sum.to_f64())
}
