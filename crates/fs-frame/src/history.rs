//! Stage 3: nonlinear TIME HISTORY of a single-story frame with two
//! fiber-hinge columns — the smoke-tier concentrated-plasticity
//! idealization: story drift x maps to base-hinge curvature
//! κ = x/(h·l_p), the TRUE fiber section (Mander concrete core +
//! Menegotto–Pinto steel through fs-solid/fs-material, with all the
//! sign conventions tfz.14 pinned) returns the hinge moment, and the
//! story shear is V = 2M/h. Newmark average acceleration with Newton
//! on the section tangent; hysteresis and stiffness degradation come
//! from the fibers, not from a phenomenological spring. The
//! distributed-plasticity frame (fs-solid ForceBasedElement columns)
//! is the recorded successor.

use fs_solid::fiber::{Section, rc_section};

/// A unit-explicit, borrowed ground-acceleration record.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GroundMotion<'a> {
    /// Relative ground acceleration samples in metres per second squared.
    /// [`StoryFrame::run_checked`] applies each value at a step endpoint.
    /// [`StoryFrame::run_from_time_zero`] instead uses the first value to
    /// establish initial force equilibrium and advances through the remaining rows.
    pub acceleration_m_s2: &'a [f64],
    /// Uniform sample interval in seconds.
    pub dt_s: f64,
}

impl<'a> GroundMotion<'a> {
    /// Bind acceleration samples to their uniform sample interval.
    #[must_use]
    pub const fn new(acceleration_m_s2: &'a [f64], dt_s: f64) -> Self {
        Self {
            acceleration_m_s2,
            dt_s,
        }
    }
}

/// Explicit work and convergence limits for a checked history run.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HistoryLimits {
    /// Maximum admitted acceleration samples.
    pub max_samples: usize,
    /// Maximum Newton corrections per sample.
    pub max_newton_iterations: usize,
    /// Absolute displacement-correction tolerance in metres.
    pub displacement_tolerance_m: f64,
    /// Absolute dynamic-equilibrium residual tolerance in newtons.
    pub equilibrium_tolerance_n: f64,
}

impl HistoryLimits {
    /// Construct explicit record and nonlinear-solve limits.
    #[must_use]
    pub const fn new(
        max_samples: usize,
        max_newton_iterations: usize,
        displacement_tolerance_m: f64,
        equilibrium_tolerance_n: f64,
    ) -> Self {
        Self {
            max_samples,
            max_newton_iterations,
            displacement_tolerance_m,
            equilibrium_tolerance_n,
        }
    }
}

/// Auditable response histories from an admitted ground motion.
#[derive(Debug, Clone, PartialEq)]
pub struct HistoryResponse {
    /// Relative story displacement at every advancing endpoint, in metres.
    /// An explicitly initial-time row does not add an advancing endpoint.
    pub displacement_m: Vec<f64>,
    /// Fiber-section restoring shear at every admitted sample, in newtons.
    /// This excludes the viscous-damping contribution to a full support
    /// reaction and must not be labeled as that reaction.
    pub restoring_shear_n: Vec<f64>,
    /// Peak absolute relative story displacement, in metres.
    pub peak_abs_displacement_m: f64,
    /// Peak absolute fiber-section restoring shear, in newtons.
    pub peak_abs_restoring_shear_n: f64,
    /// Largest absolute final dynamic-equilibrium residual, in newtons.
    pub max_abs_equilibrium_residual_n: f64,
    /// Largest number of Newton corrections used by any sample.
    pub max_newton_iterations_used: usize,
}

/// Structured refusal from the bounded time-history path.
#[derive(Debug, Clone, PartialEq)]
pub enum HistoryError {
    /// A recorded-motion comparison requires at least one sample.
    EmptyRecord,
    /// An initial-time record needs at least one later sample to advance.
    MissingAdvancingSample,
    /// The sample interval was not finite and strictly positive.
    InvalidTimeStep {
        /// Rejected interval in seconds.
        dt_s: f64,
    },
    /// The record exceeds the caller's explicit work limit.
    SampleLimitExceeded {
        /// Offered sample count.
        samples: usize,
        /// Admitted sample count.
        max_samples: usize,
    },
    /// No Newton correction was permitted.
    EmptyNewtonBudget,
    /// The convergence tolerance was not finite and strictly positive.
    InvalidDisplacementTolerance {
        /// Rejected tolerance in metres.
        tolerance_m: f64,
    },
    /// The dynamic-equilibrium residual tolerance was not finite and positive.
    InvalidEquilibriumTolerance {
        /// Rejected tolerance in newtons.
        tolerance_n: f64,
    },
    /// A story parameter was outside its physical domain.
    InvalidStoryParameter {
        /// Parameter name.
        name: &'static str,
        /// Rejected value.
        value: f64,
    },
    /// A ground-acceleration sample was NaN or infinite.
    NonFiniteAcceleration {
        /// Zero-based sample index.
        sample: usize,
        /// Rejected acceleration in metres per second squared.
        acceleration_m_s2: f64,
    },
    /// Response history storage could not be reserved before integration.
    AllocationFailed {
        /// Requested number of samples in each response channel.
        samples: usize,
    },
    /// A nonlinear sample produced a non-finite dynamic quantity.
    NonFiniteState {
        /// Zero-based sample index.
        sample: usize,
    },
    /// Newton exhausted its per-sample limit without satisfying the correction tolerance.
    NewtonDidNotConverge {
        /// Zero-based sample index.
        sample: usize,
        /// Corrections attempted.
        iterations: usize,
        /// Last absolute correction in metres.
        last_correction_m: f64,
        /// Last absolute dynamic-equilibrium residual in newtons.
        last_residual_n: f64,
    },
}

impl core::fmt::Display for HistoryError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::EmptyRecord => write!(formatter, "ground-motion record is empty"),
            Self::MissingAdvancingSample => {
                write!(formatter, "initial-time ground motion needs a later sample")
            }
            Self::InvalidTimeStep { dt_s } => {
                write!(
                    formatter,
                    "ground-motion sample interval {dt_s} s is invalid"
                )
            }
            Self::SampleLimitExceeded {
                samples,
                max_samples,
            } => write!(
                formatter,
                "ground-motion record has {samples} samples, exceeding limit {max_samples}"
            ),
            Self::EmptyNewtonBudget => write!(formatter, "Newton iteration budget is zero"),
            Self::InvalidDisplacementTolerance { tolerance_m } => write!(
                formatter,
                "Newton displacement tolerance {tolerance_m} m is invalid"
            ),
            Self::InvalidEquilibriumTolerance { tolerance_n } => write!(
                formatter,
                "dynamic-equilibrium tolerance {tolerance_n} N is invalid"
            ),
            Self::InvalidStoryParameter { name, value } => {
                write!(formatter, "story parameter {name}={value} is invalid")
            }
            Self::NonFiniteAcceleration {
                sample,
                acceleration_m_s2,
            } => write!(
                formatter,
                "ground acceleration at sample {sample} is non-finite ({acceleration_m_s2} m/s^2)"
            ),
            Self::AllocationFailed { samples } => write!(
                formatter,
                "could not reserve {samples} samples for each response-history channel"
            ),
            Self::NonFiniteState { sample } => {
                write!(
                    formatter,
                    "time-history state became non-finite at sample {sample}"
                )
            }
            Self::NewtonDidNotConverge {
                sample,
                iterations,
                last_correction_m,
                last_residual_n,
            } => write!(
                formatter,
                "Newton did not converge at sample {sample} after {iterations} corrections; \
                 last |dx|={last_correction_m} m, |residual|={last_residual_n} N"
            ),
        }
    }
}

impl std::error::Error for HistoryError {}

#[derive(Debug, Clone, Copy)]
struct StepControls {
    dt_s: f64,
    damping_n_s_m: f64,
    max_newton_iterations: usize,
    displacement_tolerance_m: f64,
    equilibrium_tolerance_n: f64,
}

#[derive(Debug, Clone, Copy)]
struct StepResponse {
    displacement_m: f64,
    restoring_shear_n: f64,
    equilibrium_residual_n: f64,
    iterations: usize,
}

/// Story model parameters.
#[derive(Debug, Clone, Copy)]
pub struct StoryParams {
    /// Story height (m).
    pub h: f64,
    /// Plastic-hinge length (m).
    pub lp: f64,
    /// Story mass (kg).
    pub mass: f64,
    /// Damping ratio (Rayleigh mass-proportional at the initial
    /// period).
    pub zeta: f64,
    /// Section scale: fiber areas multiply by this (the CVaR design
    /// variable).
    pub scale: f64,
}

impl Default for StoryParams {
    fn default() -> Self {
        StoryParams {
            h: 3.0,
            lp: 0.45,
            // 280 t story mass over the two-column pair: T ≈ 0.5 s at
            // the probed k₀ ≈ 4.4e7 N/m; yield drift ratio ≈ 0.0036
            // (V_y ≈ 4.8e5 N). SI throughout — fs-material's units.
            mass: 2.8e5,
            zeta: 0.02,
            scale: 1.0,
        }
    }
}

/// The story frame state (two identical fiber-hinge columns).
#[derive(Debug, Clone)]
pub struct StoryFrame {
    /// Parameters.
    pub params: StoryParams,
    hinge: Section,
    /// Committed drift.
    pub x: f64,
    /// Committed velocity.
    pub v: f64,
    /// Committed acceleration (relative).
    pub a: f64,
}

/// Scale a section's fiber areas.
fn scaled_section(scale: f64) -> Section {
    let mut s = rc_section(0.5, 0.35, 12, 0.002);
    for f in &mut s.fibers {
        f.area *= scale;
    }
    s
}

impl StoryFrame {
    /// A story at rest.
    #[must_use]
    pub fn new(params: StoryParams) -> StoryFrame {
        StoryFrame {
            params,
            hinge: scaled_section(params.scale),
            x: 0.0,
            v: 0.0,
            a: 0.0,
        }
    }

    /// Story restoring shear and tangent stiffness at drift `x`
    /// (trial — no commit).
    #[must_use]
    pub fn restoring(&self, x: f64) -> (f64, f64) {
        let kappa = x / (self.params.h * self.params.lp);
        let st = self.hinge.respond(0.0, kappa);
        let v = 2.0 * st.m / self.params.h;
        let dv_dx = 2.0 * st.tangent[1][1] / (self.params.h * self.params.h * self.params.lp);
        (v, dv_dx)
    }

    /// Initial (elastic) story stiffness.
    #[must_use]
    pub fn initial_stiffness(&self) -> f64 {
        self.restoring(1e-9).1
    }

    fn damping_coefficient(&self) -> f64 {
        let k0 = self.initial_stiffness();
        2.0 * self.params.zeta * (k0 * self.params.mass).sqrt()
    }

    #[allow(clippy::too_many_lines)] // One auditable Newmark trial/commit transaction.
    fn advance_sample<const CHECKED: bool>(
        &mut self,
        ground_acceleration_m_s2: f64,
        sample: usize,
        controls: StepControls,
    ) -> Result<StepResponse, HistoryError> {
        let m = self.params.mass;
        let (beta, gamma) = (0.25f64, 0.5f64);
        let p_ext = -m * ground_acceleration_m_s2;
        let x0 = self.x;
        let v0 = self.v;
        let a0 = self.a;
        let mut a_new = a0;
        let h = controls.dt_s;
        let position_weight = beta * h * h;
        // Solve for acceleration itself. Recovering it from x_new - x0
        // amplifies absolute-position rounding by 1/h² and can lose a finite
        // inertial response around a preloaded or plastically displaced state.
        let kinematics = |acceleration: f64| {
            let increment = h.mul_add(v0, (h * h) * beta.mul_add(acceleration, (0.5 - beta) * a0));
            let velocity = v0 + h * gamma.mul_add(acceleration, (1.0 - gamma) * a0);
            (x0 + increment, velocity)
        };
        let mut converged = false;
        let mut last_correction_m = f64::INFINITY;
        let mut last_residual_n = f64::INFINITY;
        let mut iterations = 0usize;

        for iteration in 1..=controls.max_newton_iterations {
            let (x, v_new) = kinematics(a_new);
            let (fs, kt) = self.restoring(x);
            let residual = m * a_new + controls.damping_n_s_m * v_new + fs - p_ext;
            let dynamic_tangent = m + gamma * h * controls.damping_n_s_m + position_weight * kt;
            let correction = -residual / dynamic_tangent;
            a_new += correction;
            iterations = iteration;
            last_correction_m = (position_weight * correction).abs();
            last_residual_n = residual.abs();
            if CHECKED
                && (!a_new.is_finite()
                    || !v_new.is_finite()
                    || !fs.is_finite()
                    || !kt.is_finite()
                    || !residual.is_finite()
                    || !dynamic_tangent.is_finite()
                    || !correction.is_finite()
                    || !x.is_finite())
            {
                return Err(HistoryError::NonFiniteState { sample });
            }
            if last_correction_m < controls.displacement_tolerance_m
                && (!CHECKED || last_residual_n < controls.equilibrium_tolerance_n)
            {
                converged = true;
                break;
            }
        }

        if CHECKED && !converged {
            return Err(HistoryError::NewtonDidNotConverge {
                sample,
                iterations,
                last_correction_m,
                last_residual_n,
            });
        }

        let (x, v_new) = kinematics(a_new);
        let restoring_shear_n = if CHECKED { self.restoring(x).0 } else { 0.0 };
        let equilibrium_residual_n = if CHECKED {
            m * a_new + controls.damping_n_s_m * v_new + restoring_shear_n - p_ext
        } else {
            0.0
        };
        if CHECKED
            && (!a_new.is_finite()
                || !v_new.is_finite()
                || !restoring_shear_n.is_finite()
                || !equilibrium_residual_n.is_finite())
        {
            return Err(HistoryError::NonFiniteState { sample });
        }
        if CHECKED && equilibrium_residual_n.abs() >= controls.equilibrium_tolerance_n {
            return Err(HistoryError::NewtonDidNotConverge {
                sample,
                iterations,
                last_correction_m,
                last_residual_n: equilibrium_residual_n.abs(),
            });
        }
        let kappa = x / (self.params.h * self.params.lp);
        self.hinge.commit(0.0, kappa);
        self.x = x;
        self.v = v_new;
        self.a = a_new;
        Ok(StepResponse {
            displacement_m: x,
            restoring_shear_n,
            equilibrium_residual_n,
            iterations,
        })
    }

    /// Run the record under ground acceleration `ag` sampled at `dt`;
    /// returns the drift history. Newmark average acceleration
    /// (γ = ½, β = ¼) with Newton on the fiber tangent; the section
    /// commits once per step.
    pub fn run(&mut self, ag: &[f64], dt: f64) -> Vec<f64> {
        let damping_n_s_m = self.damping_coefficient();
        let controls = StepControls {
            dt_s: dt,
            damping_n_s_m,
            max_newton_iterations: 30,
            displacement_tolerance_m: 1e-12,
            equilibrium_tolerance_n: f64::INFINITY,
        };
        let mut drifts = Vec::with_capacity(ag.len());
        for (sample, &agi) in ag.iter().enumerate() {
            let response = self
                .advance_sample::<false>(agi, sample, controls)
                .expect("unchecked time-history step cannot return a refusal");
            drifts.push(response.displacement_m);
        }
        drifts
    }

    /// Run a unit-explicit ground-motion record under explicit work limits.
    ///
    /// Input validation and response allocation complete before integration.
    /// Integration then runs on a cloned frame and publishes both the response
    /// and the committed fiber state only after every sample converges. A
    /// refusal therefore leaves `self` unchanged.
    ///
    /// # Errors
    ///
    /// Returns [`HistoryError`] when the record, physical parameters, limits,
    /// allocation, or any nonlinear sample is inadmissible.
    pub fn run_checked(
        &mut self,
        motion: GroundMotion<'_>,
        limits: HistoryLimits,
    ) -> Result<HistoryResponse, HistoryError> {
        self.run_checked_on_grid(motion, limits, None)
    }

    /// Run samples at `t = 0, dt, ..., (N-1)dt`.
    ///
    /// The first row sets relative acceleration from the original force
    /// balance, including current restoring force and viscous damping. It does
    /// not advance or commit the material. N rows produce N-1 advancing
    /// endpoints; peak values also include the initial state. The complete
    /// source record counts against the sample budget. As in
    /// [`Self::run_checked`], every refusal leaves the frame unchanged.
    ///
    /// # Errors
    /// Original history refusals, or fewer than two source samples.
    pub fn run_from_time_zero(
        &mut self,
        motion: GroundMotion<'_>,
        limits: HistoryLimits,
    ) -> Result<HistoryResponse, HistoryError> {
        self.run_sampled_checked(motion, limits, 1)
    }

    /// Initial-time integration at a positive source-index stride, retaining
    /// the final row with a shorter last interval when needed. The fragility
    /// owner admits the common spectral band before using stride two.
    pub(crate) fn run_sampled_checked(
        &mut self,
        motion: GroundMotion<'_>,
        limits: HistoryLimits,
        stride: usize,
    ) -> Result<HistoryResponse, HistoryError> {
        assert!(stride > 0, "history sample stride must be positive");
        self.run_checked_on_grid(motion, limits, Some(stride))
    }

    #[allow(clippy::too_many_lines)] // One preflight and staged Newmark owner for both clocks.
    fn run_checked_on_grid(
        &mut self,
        motion: GroundMotion<'_>,
        limits: HistoryLimits,
        initial_time_stride: Option<usize>,
    ) -> Result<HistoryResponse, HistoryError> {
        validate_story_parameters(self.params)?;
        if motion.acceleration_m_s2.is_empty() {
            return Err(HistoryError::EmptyRecord);
        }
        if initial_time_stride.is_some() && motion.acceleration_m_s2.len() < 2 {
            return Err(HistoryError::MissingAdvancingSample);
        }
        if !motion.dt_s.is_finite() || motion.dt_s <= 0.0 {
            return Err(HistoryError::InvalidTimeStep { dt_s: motion.dt_s });
        }
        if motion.acceleration_m_s2.len() > limits.max_samples {
            return Err(HistoryError::SampleLimitExceeded {
                samples: motion.acceleration_m_s2.len(),
                max_samples: limits.max_samples,
            });
        }
        if limits.max_newton_iterations == 0 {
            return Err(HistoryError::EmptyNewtonBudget);
        }
        if !limits.displacement_tolerance_m.is_finite() || limits.displacement_tolerance_m <= 0.0 {
            return Err(HistoryError::InvalidDisplacementTolerance {
                tolerance_m: limits.displacement_tolerance_m,
            });
        }
        if !limits.equilibrium_tolerance_n.is_finite() || limits.equilibrium_tolerance_n <= 0.0 {
            return Err(HistoryError::InvalidEquilibriumTolerance {
                tolerance_n: limits.equilibrium_tolerance_n,
            });
        }
        for (sample, &acceleration_m_s2) in motion.acceleration_m_s2.iter().enumerate() {
            if !acceleration_m_s2.is_finite() {
                return Err(HistoryError::NonFiniteAcceleration {
                    sample,
                    acceleration_m_s2,
                });
            }
        }

        let initial_stiffness_n_m = self.initial_stiffness();
        if !initial_stiffness_n_m.is_finite() || initial_stiffness_n_m <= 0.0 {
            return Err(HistoryError::InvalidStoryParameter {
                name: "initial story stiffness",
                value: initial_stiffness_n_m,
            });
        }
        let damping_n_s_m =
            2.0 * self.params.zeta * (initial_stiffness_n_m * self.params.mass).sqrt();
        if !damping_n_s_m.is_finite() {
            return Err(HistoryError::InvalidStoryParameter {
                name: "derived damping coefficient",
                value: damping_n_s_m,
            });
        }

        let source_samples = motion.acceleration_m_s2.len();
        let samples = initial_time_stride.map_or(source_samples, |stride| {
            (source_samples - 1).div_ceil(stride)
        });
        let mut displacement_m = Vec::new();
        displacement_m
            .try_reserve_exact(samples)
            .map_err(|_| HistoryError::AllocationFailed { samples })?;
        let mut restoring_shear_n = Vec::new();
        restoring_shear_n
            .try_reserve_exact(samples)
            .map_err(|_| HistoryError::AllocationFailed { samples })?;

        let mut staged = self.clone();
        let mut peak_abs_displacement_m = 0.0f64;
        let mut peak_abs_restoring_shear_n = 0.0f64;
        if initial_time_stride.is_some() {
            let restoring = staged.restoring(staged.x).0;
            staged.a = -motion.acceleration_m_s2[0]
                - (damping_n_s_m * staged.v + restoring) / staged.params.mass;
            if !staged.x.is_finite()
                || !staged.v.is_finite()
                || !staged.a.is_finite()
                || !restoring.is_finite()
            {
                return Err(HistoryError::NonFiniteState { sample: 0 });
            }
            peak_abs_displacement_m = staged.x.abs();
            peak_abs_restoring_shear_n = restoring.abs();
        }
        let mut max_abs_equilibrium_residual_n = 0.0f64;
        let mut max_newton_iterations_used = 0usize;
        let mut controls = StepControls {
            dt_s: motion.dt_s,
            damping_n_s_m,
            max_newton_iterations: limits.max_newton_iterations,
            displacement_tolerance_m: limits.displacement_tolerance_m,
            equilibrium_tolerance_n: limits.equilibrium_tolerance_n,
        };
        let stride = initial_time_stride.unwrap_or(1);
        let start = initial_time_stride.unwrap_or(0);
        let last = source_samples - 1;
        let tail = (initial_time_stride.is_some() && last % stride != 0).then_some(last);
        let mut previous = 0;
        for sample in (start..source_samples).step_by(stride).chain(tail) {
            if initial_time_stride.is_some() {
                #[allow(clippy::cast_precision_loss)]
                let interval = (sample - previous) as f64;
                controls.dt_s = interval * motion.dt_s;
                if !controls.dt_s.is_finite() {
                    return Err(HistoryError::InvalidTimeStep {
                        dt_s: controls.dt_s,
                    });
                }
            }
            let response = staged.advance_sample::<true>(
                motion.acceleration_m_s2[sample],
                sample,
                controls,
            )?;
            previous = sample;
            peak_abs_displacement_m = peak_abs_displacement_m.max(response.displacement_m.abs());
            peak_abs_restoring_shear_n =
                peak_abs_restoring_shear_n.max(response.restoring_shear_n.abs());
            max_abs_equilibrium_residual_n =
                max_abs_equilibrium_residual_n.max(response.equilibrium_residual_n.abs());
            max_newton_iterations_used = max_newton_iterations_used.max(response.iterations);
            displacement_m.push(response.displacement_m);
            restoring_shear_n.push(response.restoring_shear_n);
        }

        *self = staged;
        Ok(HistoryResponse {
            displacement_m,
            restoring_shear_n,
            peak_abs_displacement_m,
            peak_abs_restoring_shear_n,
            max_abs_equilibrium_residual_n,
            max_newton_iterations_used,
        })
    }
}

fn validate_story_parameters(params: StoryParams) -> Result<(), HistoryError> {
    for (name, value) in [
        ("story height", params.h),
        ("plastic-hinge length", params.lp),
        ("story mass", params.mass),
        ("section scale", params.scale),
    ] {
        if !value.is_finite() || value <= 0.0 {
            return Err(HistoryError::InvalidStoryParameter { name, value });
        }
    }
    if !params.zeta.is_finite() || params.zeta < 0.0 {
        return Err(HistoryError::InvalidStoryParameter {
            name: "damping ratio",
            value: params.zeta,
        });
    }
    for (name, value) in [
        ("curvature denominator h*lp", params.h * params.lp),
        (
            "stiffness denominator h*h*lp",
            params.h * params.h * params.lp,
        ),
    ] {
        if !value.is_finite() || value <= 0.0 {
            return Err(HistoryError::InvalidStoryParameter { name, value });
        }
    }
    Ok(())
}

/// Peak drift RATIO (|x|max / h) of a drift history.
#[must_use]
pub fn peak_drift(drifts: &[f64], h: f64) -> f64 {
    drifts.iter().fold(0.0f64, |m, &x| m.max(x.abs())) / h
}

#[cfg(test)]
mod sampled_history_tests {
    use super::{GroundMotion, HistoryError, HistoryLimits, StoryFrame, StoryParams};

    fn limits() -> HistoryLimits {
        HistoryLimits::new(64, 30, 1e-12, 1.0)
    }

    #[test]
    fn g1_time_zero_uses_existing_force_equilibrium_and_rolls_back() {
        let mut frame = StoryFrame::new(StoryParams::default());
        frame
            .run_checked(GroundMotion::new(&[0.5; 8], 0.01), limits())
            .unwrap();
        assert!(frame.x != 0.0 && frame.v != 0.0);
        frame.a = 5_000.0; // The new initial forcing must replace this stale value.
        let mut reference = frame.clone();
        let initial_x = frame.x;
        let restoring = reference.restoring(reference.x).0;
        let damping = 2.0
            * reference.params.zeta
            * (reference.initial_stiffness() * reference.params.mass).sqrt();
        let acceleration = [0.8, -0.4, 0.2];
        reference.a =
            -acceleration[0] - (damping * reference.v + restoring) / reference.params.mass;
        let expected = reference
            .run_checked(GroundMotion::new(&acceleration[1..], 0.005), limits())
            .unwrap();
        let actual = frame
            .run_from_time_zero(GroundMotion::new(&acceleration, 0.005), limits())
            .unwrap();
        assert_eq!(actual.displacement_m, expected.displacement_m);
        assert_eq!(actual.restoring_shear_n, expected.restoring_shear_n);
        assert_eq!(actual.displacement_m.len(), acceleration.len() - 1);
        assert_eq!(
            actual.peak_abs_displacement_m,
            initial_x.abs().max(expected.peak_abs_displacement_m)
        );
        assert_eq!(
            [frame.x.to_bits(), frame.v.to_bits(), frame.a.to_bits()],
            [
                reference.x.to_bits(),
                reference.v.to_bits(),
                reference.a.to_bits()
            ]
        );

        let mut untouched = frame.clone();
        let state = [frame.x.to_bits(), frame.v.to_bits(), frame.a.to_bits()];
        assert!(matches!(
            frame.run_from_time_zero(
                GroundMotion::new(&[0.8, -0.4], 0.005),
                HistoryLimits::new(2, 1, 1e-30, 1e-30),
            ),
            Err(HistoryError::NewtonDidNotConverge { .. })
        ));
        assert_eq!(
            [frame.x.to_bits(), frame.v.to_bits(), frame.a.to_bits()],
            state
        );
        let probe = GroundMotion::new(&[0.1, -0.1], 0.005);
        assert_eq!(
            frame.run_checked(probe, limits()).unwrap(),
            untouched.run_checked(probe, limits()).unwrap()
        );
    }

    #[test]
    fn g1_strided_initial_history_retains_last_short_interval() {
        // Negligible restoring force relative to inertia: constant base
        // acceleration has the independent free-mass limit x=-a*t^2/2, v=-a*t.
        // Dyadic times also make the expected final horizon exact in f64.
        let params = StoryParams {
            mass: 1e20,
            zeta: 0.0,
            scale: 1e-6,
            ..StoryParams::default()
        };
        let dt = 0.125f64;
        for samples in 2usize..=8 {
            let acceleration = vec![1.0; samples];
            let motion = GroundMotion::new(&acceleration, dt);
            let mut fine = StoryFrame::new(params);
            let mut coarse = fine.clone();
            let fine_response = fine.run_from_time_zero(motion, limits()).unwrap();
            let coarse_response = coarse.run_sampled_checked(motion, limits(), 2).unwrap();
            assert_eq!(fine_response.displacement_m.len(), samples - 1);
            assert_eq!(
                coarse_response.displacement_m.len(),
                (samples - 1).div_ceil(2)
            );
            #[allow(clippy::cast_precision_loss)]
            let time = (samples - 1) as f64 * dt;
            for response in [fine, coarse] {
                assert!((response.x + 0.5 * time * time).abs() < 1e-14);
                assert!((response.v + time).abs() < 1e-14);
                assert!((response.a + 1.0).abs() < 1e-14);
            }
        }
        let mut frame = StoryFrame::new(params);
        assert!(matches!(
            frame.run_sampled_checked(GroundMotion::new(&[1.0, f64::NAN, 1.0], dt), limits(), 2,),
            Err(HistoryError::NonFiniteAcceleration { sample: 1, .. })
        ));
        assert!(matches!(
            frame.run_sampled_checked(
                GroundMotion::new(&[1.0; 5], dt),
                HistoryLimits::new(4, 30, 1e-12, 1.0),
                2,
            ),
            Err(HistoryError::SampleLimitExceeded {
                samples: 5,
                max_samples: 4
            })
        ));
        assert!(matches!(
            frame.run_from_time_zero(GroundMotion::new(&[1.0], dt), limits()),
            Err(HistoryError::MissingAdvancingSample)
        ));
    }
}
