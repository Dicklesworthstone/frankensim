//! Held-effort contact using the original integrated Hertz slip/spin surface.
//!
//! Radius-scaled conjugate coordinates are `x = [du_x, du_y, a*dtheta]` and
//! `w = [F_x, F_y, M/a]`. The acting wrench obeys `x = x_free + C*w`.
//! Static admission uses the joint Hertz surface, with priority over sliding.
//! Otherwise a bounded convex projection uses the original analytic support
//! function. There is no friction ellipse, velocity smoothing or elastic memory.

use core::{f64::consts::PI, fmt};

use fs_math::det;

use super::{HertzSlidingPatch, HertzSlidingPatchParameters, shape_factors};
use crate::TriboError;

/// Simultaneous original Hertz patches with full cross-contact response.
pub mod multiple;

type Vector = [f64; 3];
type Matrix = [[f64; 3]; 3];
const SPIN_RADIUS: f64 = 3.0 * PI / 16.0;

/// Original normal/interface data, frozen over this contact step.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HertzWrenchParameters {
    /// Prescribed compressive force [N], including the exact zero-load limit.
    pub normal_force_n: f64,
    /// Circular Hertz radius [m], strictly positive.
    pub contact_radius_m: f64,
    /// Uniform local static coefficient, at least `kinetic_mu`.
    pub static_mu: f64,
    /// Uniform local sliding coefficient, nonnegative.
    pub kinetic_mu: f64,
}

/// Independent numerical tolerances and physical accepted-effort budgets.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HertzWrenchConfig {
    /// Bounded bisections of the original static surface meridian.
    pub maximum_static_iterations: usize,
    /// Bounded fully corrective kinetic projection iterations.
    pub maximum_iterations: usize,
    /// Total native support evaluations, including static admission.
    pub maximum_evaluations: usize,
    /// Absolute native tangential-force residual [N].
    pub force_tolerance_n: f64,
    /// Absolute native normal-axis moment residual [N m].
    pub moment_tolerance_nm: f64,
    /// Absolute compatibility tolerance in all radius-scaled coordinates [m].
    pub displacement_tolerance_m: f64,
    /// Absolute actual-work, dual-gap and compatibility-work tolerance [J].
    pub energy_tolerance_j: f64,
    /// Relative matrix symmetry and minimum normalized Cholesky pivot tolerance.
    /// A failed pivot is refused; no regularizing compliance is added.
    pub relative_matrix_tolerance: f64,
    /// Maximum accepted resultant tangential force [N].
    pub maximum_force_n: f64,
    /// Maximum accepted normal-axis moment magnitude [N m].
    pub maximum_spin_moment_nm: f64,
}

impl HertzWrenchConfig {
    /// Admission reusable by mechanical owners before preparing a trial.
    pub fn validate(self) -> Result<(), HertzWrenchError> {
        if self.maximum_static_iterations == 0
            || self.maximum_iterations == 0
            || self.maximum_evaluations == 0
        {
            return Err(HertzWrenchError::Invalid("iteration_budget"));
        }
        for value in [
            self.force_tolerance_n,
            self.moment_tolerance_nm,
            self.displacement_tolerance_m,
            self.energy_tolerance_j,
            self.maximum_force_n,
            self.maximum_spin_moment_nm,
        ] {
            if !value.is_finite() || value < 0.0 {
                return Err(HertzWrenchError::Invalid("tolerance_or_effort_budget"));
            }
        }
        if !self.relative_matrix_tolerance.is_finite()
            || !(0.0..=1.0e-8).contains(&self.relative_matrix_tolerance)
        {
            return Err(HertzWrenchError::Invalid("relative_matrix_tolerance"));
        }
        Ok(())
    }
}

/// Released motion and the complete physical response of its mechanical owner.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HertzWrenchInput {
    /// Centre relative displacement without this contact's wrench [m].
    pub free_displacement_m: [f64; 2],
    /// Relative normal-axis rotation without this contact's wrench [rad].
    pub free_spin_angle_rad: f64,
    /// Full radius-scaled body compliance [m/N]. If the unscaled coordinates
    /// are `[du_x, du_y, dtheta]`, supply `S*C_physical*S`, `S=diag(1,1,a)`.
    /// This rung admits positive definite matrices and the exact all-zero case.
    pub scaled_compliance_m_per_n: Matrix,
}

/// Physical phase supplied by the mechanical owner, not inferred from a
/// midpoint displacement that may already include an unresolved reversal.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum HertzWrenchPhase {
    /// The original static set may support a stationary accepted contact.
    /// Also retains the original displacement-only solver's admission policy.
    StaticAdmissible,
    /// The accepted shear/spin motion is nonzero. The larger static set cannot
    /// capture this moving phase; only the original kinetic capacity is used.
    /// This does not itself locate an endpoint stopping event.
    Kinetic,
}

/// Branch selected from the same joint translational/torsional friction surface.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum HertzWrenchMode {
    /// The full static reaction is feasible, with no relative motion.
    Sticking,
    /// Nonzero relative motion selects the integrated kinetic surface.
    Sliding,
    /// The admitted phase has zero traction capacity. A kinetic phase may have
    /// zero kinetic coefficient while its original static coefficient is positive.
    Frictionless,
}

/// Candidate only: no history, mechanical state or heat owner is mutated.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HertzWrenchStep {
    /// Force acting on the side whose released relative motion was supplied [N].
    pub force_n: [f64; 2],
    /// Acting moment about the patch centre and declared normal [N m].
    pub spin_moment_nm: f64,
    /// Actual centre relative displacement [m].
    pub slip_displacement_m: [f64; 2],
    /// Actual relative normal-axis rotation [rad].
    pub spin_angle_rad: f64,
    /// Positive absorbed work `-F dot du - M*dtheta` [J].
    pub contact_work_j: f64,
    /// Original integrated local sliding dissipation [J]; unpartitioned.
    pub dissipated_energy_j: f64,
    /// `dissipated_energy_j - contact_work_j`, never added to heat [J].
    pub work_residual_j: f64,
    /// Native support optimality gap at the actual motion [J].
    pub duality_gap_j: f64,
    /// `x - x_free - C*w_acting`, using the original supplied matrix [m].
    pub compatibility_residual_m: Vector,
    /// Opposing-effort pairing with the compatibility residual [J].
    pub body_work_residual_j: f64,
    /// Norm of actual native force disagreement [N]; zero in stiction.
    pub force_residual_n: f64,
    /// Actual native moment disagreement [N m]; zero in stiction.
    pub spin_moment_residual_nm: f64,
    /// Static, sliding or zero-capacity branch.
    pub mode: HertzWrenchMode,
    /// Fully corrective kinetic iterations performed.
    pub iterations: usize,
    /// Original support evaluations, including the static boundary search.
    pub evaluations: usize,
}

/// Explicit refusal; the original contact and mechanical inputs stay unchanged.
#[derive(Clone, Debug, PartialEq)]
pub enum HertzWrenchError {
    /// Invalid input, singular/indefinite response, or nonfinite derived value.
    Invalid(&'static str),
    /// Original Hertz-patch parameter admission refused the input.
    Law(TriboError),
    /// Caller cancellation observed at a bounded evaluation/iteration boundary.
    Cancelled,
    /// The accepted physical wrench would exceed the caller's effort budget.
    EffortBudget {
        /// Required resultant force [N].
        force_n: f64,
        /// Required moment magnitude [N m].
        spin_moment_nm: f64,
    },
    /// The static boundary cannot be classified within its declared budget.
    /// This refuses instead of silently selecting the kinetic branch.
    StaticBoundaryUnresolved {
        /// Remaining force bracket width [N].
        force_bracket_n: f64,
        /// Remaining moment bracket width [N m].
        moment_bracket_nm: f64,
        /// Original support evaluations consumed.
        evaluations: usize,
    },
    /// Native support, work, or full-body compatibility failed its budget.
    Unresolved {
        /// Kinetic correction count.
        iterations: usize,
        /// Original support evaluations consumed.
        evaluations: usize,
        /// Native force disagreement [N].
        force_residual_n: f64,
        /// Native moment disagreement [N m].
        spin_moment_residual_nm: f64,
        /// Largest absolute original-matrix compatibility residual [m].
        displacement_residual_m: f64,
        /// Largest absolute work or dual-gap residual [J].
        energy_residual_j: f64,
    },
}

impl fmt::Display for HertzWrenchError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Invalid(field) => write!(f, "invalid Hertz wrench input: {field}"),
            Self::Law(error) => write!(f, "Hertz wrench law: {error}"),
            Self::Cancelled => write!(f, "Hertz wrench solve cancelled"),
            Self::EffortBudget {
                force_n,
                spin_moment_nm,
            } => write!(
                f,
                "Hertz wrench exceeds effort budget: {force_n} N, {spin_moment_nm} N m"
            ),
            Self::StaticBoundaryUnresolved {
                force_bracket_n,
                moment_bracket_nm,
                evaluations,
            } => write!(
                f,
                "Hertz static boundary unresolved after {evaluations} evaluations: {force_bracket_n} N, {moment_bracket_nm} N m brackets"
            ),
            Self::Unresolved {
                iterations,
                evaluations,
                force_residual_n,
                spin_moment_residual_nm,
                displacement_residual_m,
                energy_residual_j,
            } => write!(
                f,
                "Hertz wrench unresolved after {iterations} corrections/{evaluations} evaluations: {force_residual_n} N, {spin_moment_residual_nm} N m, {displacement_residual_m} m, {energy_residual_j} J"
            ),
        }
    }
}

impl std::error::Error for HertzWrenchError {}

/// Frozen normal/interface law. No extra material authority or storage is added.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HertzWrenchContact {
    parameters: HertzWrenchParameters,
    static_patch: HertzSlidingPatch,
    kinetic_patch: HertzSlidingPatch,
}

impl HertzWrenchContact {
    /// Admit both original Hertz support functions and static/kinetic ordering.
    pub fn new(parameters: HertzWrenchParameters) -> Result<Self, HertzWrenchError> {
        let patch = |friction_coefficient| {
            HertzSlidingPatch::new(HertzSlidingPatchParameters {
                normal_force_n: parameters.normal_force_n,
                contact_radius_m: parameters.contact_radius_m,
                friction_coefficient,
            })
            .map_err(HertzWrenchError::Law)
        };
        let static_patch = patch(parameters.static_mu)?;
        let kinetic_patch = patch(parameters.kinetic_mu)?;
        if parameters.static_mu < parameters.kinetic_mu {
            return Err(HertzWrenchError::Invalid("static_mu_below_kinetic_mu"));
        }
        Ok(Self {
            parameters,
            static_patch,
            kinetic_patch,
        })
    }

    /// Original supplied physical parameters, unchanged.
    #[must_use]
    pub const fn parameters(&self) -> HertzWrenchParameters {
        self.parameters
    }

    /// Resolve a mechanically identified phase without changing the original
    /// source parameters. A kinetic phase uses its original kinetic support for
    /// both feasibility and sliding; this prevents the larger static surface
    /// from selecting a zero-midpoint reflected velocity before physical arrest.
    /// A zero-midpoint result still requires the mechanical owner's endpoint
    /// and event checks. No stopping time or physical trajectory is inferred.
    pub fn advance_phase(
        &self,
        input: &HertzWrenchInput,
        phase: HertzWrenchPhase,
        config: HertzWrenchConfig,
        cancelled: impl FnMut() -> bool,
    ) -> Result<HertzWrenchStep, HertzWrenchError> {
        self.phase_law(phase).advance(input, config, cancelled)
    }

    fn phase_law(&self, phase: HertzWrenchPhase) -> Self {
        match phase {
            HertzWrenchPhase::StaticAdmissible => *self,
            HertzWrenchPhase::Kinetic => Self {
                parameters: HertzWrenchParameters {
                    static_mu: self.parameters.kinetic_mu,
                    ..self.parameters
                },
                static_patch: self.kinetic_patch,
                kinetic_patch: self.kinetic_patch,
            },
        }
    }

    /// Solve one held-effort candidate with all three physical cross responses.
    ///
    /// With nonzero response, zero-slip static admission precedes the convex
    /// kinetic solve. The latter minimizes `w^T C w / 2 - x_free^T w` over the
    /// original kinetic wrench set, with `w` opposing motion. Its support atoms
    /// come from the original analytic Hertz integral. At most four feasible
    /// vertices are retained; each correction solves their five-point hull.
    /// Near the objective's rounding floor, bounded convex combinations refine
    /// the native wrench residual directly. The work gap is never thermalized.
    pub fn advance(
        &self,
        input: &HertzWrenchInput,
        config: HertzWrenchConfig,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<HertzWrenchStep, HertzWrenchError> {
        config.validate()?;
        if cancelled() {
            return Err(HertzWrenchError::Cancelled);
        }
        let radius = self.parameters.contact_radius_m;
        let free = [
            input.free_displacement_m[0],
            input.free_displacement_m[1],
            radius * input.free_spin_angle_rad,
        ];
        finite_vector(free, "released_motion")?;
        if !input.free_spin_angle_rad.is_finite() {
            return Err(HertzWrenchError::Invalid("free_spin_angle_rad"));
        }
        if input.free_spin_angle_rad != 0.0 && free[2] == 0.0 {
            return Err(HertzWrenchError::Invalid("scaled_spin_underflow"));
        }
        let (matrix, factor) = admit_matrix(input.scaled_compliance_m_per_n, config)?;
        let mut problem = Problem {
            law: self,
            input,
            free,
            config,
            evaluations: 0,
            cancelled: &mut cancelled,
        };
        if self.static_patch.sliding_capacity_n() == 0.0 {
            return problem.finish([0.0; 3], free, [0.0; 3], HertzWrenchMode::Frictionless, 0);
        }
        if free == [0.0; 3] {
            return problem.finish([0.0; 3], free, [0.0; 3], HertzWrenchMode::Sticking, 0);
        }
        if let Some(factor) = factor {
            let requested = factor
                .solve(free)
                .ok_or(HertzWrenchError::Invalid("static_reaction_overflow"))?;
            if problem.static_feasible(requested)? {
                // An admitted static reaction that exceeds the effort budget is
                // refused here; a lower kinetic branch must not hide that error.
                return problem.finish(
                    requested,
                    [0.0; 3],
                    requested,
                    HertzWrenchMode::Sticking,
                    0,
                );
            }
        }
        let first = problem.support(free)?;
        if factor.is_none() {
            return problem.finish(first, free, first, HertzWrenchMode::Sliding, 0);
        }
        let mut active = Active::new(first);
        for iteration in 0..=config.maximum_iterations {
            let actual = sub(
                free,
                multiply(input.scaled_compliance_m_per_n, active.point),
            );
            finite_vector(actual, "actual_motion")?;
            let support = problem.support(actual)?;
            let step = problem.candidate(
                active.point,
                actual,
                support,
                HertzWrenchMode::Sliding,
                iteration,
            )?;
            if problem.resolved(&step) {
                return problem.accept(step);
            }
            if iteration == config.maximum_iterations
                || problem.evaluations == config.maximum_evaluations
            {
                return Err(unresolved(&step));
            }
            if (problem.cancelled)() {
                return Err(HertzWrenchError::Cancelled);
            }
            let mut next = active.correct(support, matrix, free);
            let energy_scale = dot(active.point, free).abs()
                + dot(active.point, multiply(matrix, active.point)).abs();
            if next.point == active.point
                || step.duality_gap_j.abs() <= 128.0 * f64::EPSILON * energy_scale
            {
                // At the floating-point objective floor, an O(eps) energy gap
                // can still hide O(sqrt(eps)) wrench error. Refine the native
                // residual itself. Every trial is a convex combination of two
                // feasible points, so this adds no friction regularization.
                // For r=w-support(u-Cw), its derivative is I+H*C with H>=0.
                // Thus -r decreases ||r||_C locally; bounded backtracking can
                // refine the wrench even when quadratic scores cannot rank it.
                let merit = residual_metric(sub(active.point, support), matrix);
                let mut fraction = 1.0;
                for _ in 0..32 {
                    if problem.evaluations == config.maximum_evaluations {
                        break;
                    }
                    let trial = core::array::from_fn(|i| {
                        (1.0 - fraction) * active.point[i] + fraction * support[i]
                    });
                    fraction *= 0.5;
                    if trial == active.point {
                        break;
                    }
                    let motion = sub(free, multiply(input.scaled_compliance_m_per_n, trial));
                    let native = problem.support(motion)?;
                    let refined = problem.candidate(
                        trial,
                        motion,
                        native,
                        HertzWrenchMode::Sliding,
                        iteration + 1,
                    )?;
                    if problem.resolved(&refined) {
                        return problem.accept(refined);
                    }
                    if residual_metric(sub(trial, native), matrix) < merit {
                        next = Active::new(trial);
                        break;
                    }
                }
            }
            if next.point == active.point {
                return Err(unresolved(&step));
            }
            active = next;
        }
        unreachable!("bounded correction loop returns at its last iteration")
    }
}

struct Problem<'a, F: FnMut() -> bool> {
    law: &'a HertzWrenchContact,
    input: &'a HertzWrenchInput,
    free: Vector,
    config: HertzWrenchConfig,
    evaluations: usize,
    cancelled: &'a mut F,
}

impl<F: FnMut() -> bool> Problem<'_, F> {
    fn poll(&mut self) -> Result<(), HertzWrenchError> {
        if (self.cancelled)() {
            return Err(HertzWrenchError::Cancelled);
        }
        if self.evaluations == self.config.maximum_evaluations {
            return Err(HertzWrenchError::Invalid("support_evaluation_budget"));
        }
        self.evaluations += 1;
        Ok(())
    }

    fn support(&mut self, motion: Vector) -> Result<Vector, HertzWrenchError> {
        self.poll()?;
        finite_vector(motion, "support_motion")?;
        let speed = det::hypot(motion[0], motion[1]);
        if !speed.is_finite() {
            return Err(HertzWrenchError::Invalid("support_motion_norm"));
        }
        let capacity = self.law.kinetic_patch.sliding_capacity_n();
        if motion == [0.0; 3] || capacity == 0.0 {
            return Ok([0.0; 3]);
        }
        // Original shape function, with already radius-scaled spin. Avoid
        // dividing a tiny displacement by an arbitrary time or patch radius.
        let (force, moment) = if motion[2] == 0.0 {
            (1.0, 0.0)
        } else if speed == 0.0 {
            (0.0, 1.0)
        } else {
            shape_factors(speed, motion[2].abs())
        };
        let force = capacity * force;
        let output = [
            if speed == 0.0 {
                0.0
            } else {
                force * (motion[0] / speed)
            },
            if speed == 0.0 {
                0.0
            } else {
                force * (motion[1] / speed)
            },
            (capacity * SPIN_RADIUS * moment).copysign(motion[2]),
        ];
        finite_vector(output, "native_support")?;
        Ok(output)
    }

    fn static_feasible(&mut self, requested: Vector) -> Result<bool, HertzWrenchError> {
        let force = det::hypot(requested[0], requested[1]);
        let moment = requested[2].abs();
        let capacity = self.law.static_patch.sliding_capacity_n();
        let spin_capacity = capacity * SPIN_RADIUS;
        if force > capacity || moment > spin_capacity {
            return Ok(false);
        }
        if force == 0.0 || moment == 0.0 {
            return Ok(true);
        }
        // The exact pure-motion extrema permit no other wrench component.
        // Interior shape values can round to an extremum before their other
        // component vanishes, so do not use that equality as static admission.
        if force == capacity || moment == spin_capacity {
            return Ok(false);
        }
        let mut lower = (0.0, 0.0, spin_capacity);
        let mut upper = (1.0, capacity, 0.0);
        for _ in 0..self.config.maximum_static_iterations {
            if moment <= upper.2 {
                return Ok(true);
            }
            if moment > lower.2 {
                return Ok(false);
            }
            if self.evaluations == self.config.maximum_evaluations {
                break;
            }
            self.poll()?;
            let parameter = 0.5 * (lower.0 + upper.0);
            if parameter == lower.0 || parameter == upper.0 {
                break;
            }
            let (f, m) = shape_factors(parameter, 1.0 - parameter);
            let boundary = (parameter, capacity * f, spin_capacity * m);
            if force == boundary.1 {
                return Ok(moment <= boundary.2);
            }
            if force < boundary.1 {
                upper = boundary;
            } else {
                lower = boundary;
            }
        }
        // Membership is decided only by the resolved inner/outer bracket,
        // never by silently expanding the physical static surface.
        if moment <= upper.2 {
            return Ok(true);
        }
        if moment > lower.2 {
            return Ok(false);
        }
        Err(HertzWrenchError::StaticBoundaryUnresolved {
            force_bracket_n: upper.1 - lower.1,
            moment_bracket_nm: (lower.2 - upper.2) * self.law.parameters.contact_radius_m,
            evaluations: self.evaluations,
        })
    }

    fn candidate(
        &self,
        wrench: Vector,
        actual: Vector,
        support: Vector,
        mode: HertzWrenchMode,
        iterations: usize,
    ) -> Result<HertzWrenchStep, HertzWrenchError> {
        let radius = self.law.parameters.contact_radius_m;
        let compatibility = sub(
            actual,
            sub(
                self.free,
                multiply(self.input.scaled_compliance_m_per_n, wrench),
            ),
        );
        let residual = sub(wrench, support);
        let spin_angle = actual[2] / radius;
        let moment = radius * wrench[2];
        let native_moment = radius * support[2];
        if (wrench[2] != 0.0 && moment == 0.0)
            || (support[2] != 0.0 && native_moment == 0.0)
            || (actual[2] != 0.0 && spin_angle == 0.0)
        {
            return Err(HertzWrenchError::Invalid("unrepresentable_spin_port"));
        }
        // Evaluate the physical force/displacement and moment/angle ports as
        // returned to the owner, independently of the scaled dual-gap pairing.
        let work = wrench[0] * actual[0] + wrench[1] * actual[1] + moment * spin_angle;
        let heat = support[0] * actual[0] + support[1] * actual[1] + native_moment * spin_angle;
        let gap = dot(sub(support, wrench), actual);
        let step = HertzWrenchStep {
            force_n: [-wrench[0], -wrench[1]],
            spin_moment_nm: -moment,
            slip_displacement_m: [actual[0], actual[1]],
            spin_angle_rad: spin_angle,
            contact_work_j: work,
            dissipated_energy_j: heat,
            work_residual_j: heat - work,
            duality_gap_j: gap,
            compatibility_residual_m: compatibility,
            body_work_residual_j: dot(wrench, compatibility),
            force_residual_n: det::hypot(residual[0], residual[1]),
            spin_moment_residual_nm: radius * residual[2].abs(),
            mode,
            iterations,
            evaluations: self.evaluations,
        };
        finite_vector(compatibility, "compatibility")?;
        for value in [
            step.spin_moment_nm,
            step.spin_angle_rad,
            work,
            heat,
            gap,
            step.work_residual_j,
            step.body_work_residual_j,
            step.force_residual_n,
            step.spin_moment_residual_nm,
        ] {
            if !value.is_finite() {
                return Err(HertzWrenchError::Invalid("candidate_overflow"));
            }
        }
        Ok(step)
    }

    fn resolved(&self, step: &HertzWrenchStep) -> bool {
        step.force_residual_n <= self.config.force_tolerance_n
            && step.spin_moment_residual_nm <= self.config.moment_tolerance_nm
            && step.work_residual_j.abs() <= self.config.energy_tolerance_j
            && step.duality_gap_j.abs() <= self.config.energy_tolerance_j
            && step.body_work_residual_j.abs() <= self.config.energy_tolerance_j
            && step
                .compatibility_residual_m
                .iter()
                .all(|v| v.abs() <= self.config.displacement_tolerance_m)
            && step.contact_work_j >= 0.0
            && step.dissipated_energy_j >= 0.0
    }

    fn accept(&mut self, step: HertzWrenchStep) -> Result<HertzWrenchStep, HertzWrenchError> {
        let force = det::hypot(step.force_n[0], step.force_n[1]);
        let moment = step.spin_moment_nm.abs();
        if force > self.config.maximum_force_n || moment > self.config.maximum_spin_moment_nm {
            return Err(HertzWrenchError::EffortBudget {
                force_n: force,
                spin_moment_nm: moment,
            });
        }
        if (self.cancelled)() {
            return Err(HertzWrenchError::Cancelled);
        }
        Ok(step)
    }

    fn finish(
        &mut self,
        wrench: Vector,
        actual: Vector,
        support: Vector,
        mode: HertzWrenchMode,
        iterations: usize,
    ) -> Result<HertzWrenchStep, HertzWrenchError> {
        let step = self.candidate(wrench, actual, support, mode, iterations)?;
        if !self.resolved(&step) {
            return Err(unresolved(&step));
        }
        self.accept(step)
    }
}

fn unresolved(step: &HertzWrenchStep) -> HertzWrenchError {
    HertzWrenchError::Unresolved {
        iterations: step.iterations,
        evaluations: step.evaluations,
        force_residual_n: step.force_residual_n,
        spin_moment_residual_nm: step.spin_moment_residual_nm,
        displacement_residual_m: max_abs(step.compatibility_residual_m),
        energy_residual_j: step
            .work_residual_j
            .abs()
            .max(step.duality_gap_j.abs())
            .max(step.body_work_residual_j.abs()),
    }
}

#[derive(Clone, Copy)]
struct Active {
    atoms: [Vector; 4],
    point: Vector,
    count: usize,
}

impl Active {
    fn new(point: Vector) -> Self {
        Self {
            atoms: [point; 4],
            point,
            count: 1,
        }
    }

    fn correct(self, support: Vector, matrix: Matrix, free: Vector) -> Self {
        let mut atoms = [[0.0; 3]; 5];
        atoms[..self.count].copy_from_slice(&self.atoms[..self.count]);
        atoms[self.count] = support;
        let count = self.count + 1;
        let gradient = sub(free, multiply(matrix, self.point));
        let mut best = self;
        let mut best_change = 0.0;
        // Caratheodory: four vertices suffice in three dimensions. A face has
        // at most three affine directions; all thirty faces are a fixed bound.
        for mask in 1_u32..(1_u32 << count) {
            let size = mask.count_ones() as usize;
            if size > 4 {
                continue;
            }
            let mut indices = [0; 4];
            let mut position = 0;
            for i in 0..count {
                if mask & (1 << i) != 0 {
                    indices[position] = i;
                    position += 1;
                }
            }
            let origin = atoms[indices[0]];
            let mut directions = [[0.0; 3]; 3];
            let mut scales = [0.0; 3];
            let mut independent = true;
            for i in 0..size - 1 {
                let difference = sub(atoms[indices[i + 1]], origin);
                scales[i] = max_abs(difference);
                if scales[i] == 0.0 {
                    independent = false;
                    break;
                }
                directions[i] = difference.map(|v| v / scales[i]);
            }
            if !independent {
                continue;
            }
            let mut gram = [[0.0; 3]; 3];
            let mut rhs = [0.0; 3];
            let origin_gradient = sub(free, multiply(matrix, origin));
            for i in 0..size - 1 {
                rhs[i] = dot(directions[i], origin_gradient);
                for j in 0..size - 1 {
                    gram[i][j] = dot(directions[i], multiply(matrix, directions[j]));
                }
            }
            let mut weights = [0.0; 4];
            weights[0] = 1.0;
            if size > 1 {
                let Some(factor) = Factor::new(gram, size - 1, 64.0 * f64::EPSILON) else {
                    continue;
                };
                let Some(solution) = factor.solve(rhs) else {
                    continue;
                };
                for i in 0..size - 1 {
                    weights[i + 1] = solution[i] / scales[i];
                    weights[0] -= weights[i + 1];
                }
            }
            if weights[..size]
                .iter()
                .any(|w| !w.is_finite() || *w < -64.0 * f64::EPSILON)
            {
                continue;
            }
            let sum: f64 = weights[..size].iter().map(|w| w.max(0.0)).sum();
            let mut point = [0.0; 3];
            let mut candidate = Self {
                atoms: [[0.0; 3]; 4],
                point,
                count: 0,
            };
            for i in 0..size {
                let weight = weights[i].max(0.0) / sum;
                if weight > 0.0 {
                    candidate.atoms[candidate.count] = atoms[indices[i]];
                    candidate.count += 1;
                    for j in 0..3 {
                        point[j] += weight * atoms[indices[i]][j];
                    }
                }
            }
            let delta = sub(point, self.point);
            let change = -dot(delta, gradient) + 0.5 * dot(delta, multiply(matrix, delta));
            if change.is_finite() && change < best_change {
                best_change = change;
                candidate.point = point;
                best = candidate;
            }
        }
        best
    }
}

#[derive(Clone, Copy)]
struct Factor {
    lower: Matrix,
    scale: f64,
    size: usize,
}

impl Factor {
    fn new(matrix: Matrix, size: usize, tolerance: f64) -> Option<Self> {
        let scale = matrix
            .iter()
            .take(size)
            .flat_map(|row| row.iter().take(size))
            .map(|v| v.abs())
            .fold(0.0, f64::max);
        if scale == 0.0 || !scale.is_finite() {
            return None;
        }
        let mut lower = [[0.0; 3]; 3];
        for i in 0..size {
            for j in 0..=i {
                let mut value = matrix[i][j] / scale;
                for k in 0..j {
                    value -= lower[i][k] * lower[j][k];
                }
                if i == j {
                    if !value.is_finite() || value <= tolerance {
                        return None;
                    }
                    lower[i][j] = det::sqrt(value);
                } else {
                    lower[i][j] = value / lower[j][j];
                }
            }
        }
        Some(Self { lower, scale, size })
    }

    fn solve(self, rhs: Vector) -> Option<Vector> {
        let mut x = [0.0; 3];
        for i in 0..self.size {
            let mut value = rhs[i] / self.scale;
            for j in 0..i {
                value -= self.lower[i][j] * x[j];
            }
            x[i] = value / self.lower[i][i];
        }
        for i in (0..self.size).rev() {
            let mut value = x[i];
            for j in i + 1..self.size {
                value -= self.lower[j][i] * x[j];
            }
            x[i] = value / self.lower[i][i];
        }
        x.iter().all(|v| v.is_finite()).then_some(x)
    }
}

fn admit_matrix(
    input: Matrix,
    config: HertzWrenchConfig,
) -> Result<(Matrix, Option<Factor>), HertzWrenchError> {
    for row in input {
        finite_vector(row, "scaled_compliance")?;
    }
    let scale = input.into_iter().map(max_abs).fold(0.0, f64::max);
    if scale == 0.0 {
        return Ok((input, None));
    }
    let mut matrix = input;
    for i in 0..3 {
        for j in 0..i {
            if (input[i][j] / scale - input[j][i] / scale).abs() > config.relative_matrix_tolerance
            {
                return Err(HertzWrenchError::Invalid("scaled_compliance_symmetry"));
            }
            matrix[i][j] = 0.5 * input[i][j] + 0.5 * input[j][i];
            matrix[j][i] = matrix[i][j];
        }
    }
    let factor = Factor::new(matrix, 3, config.relative_matrix_tolerance).ok_or(
        HertzWrenchError::Invalid("scaled_compliance_positive_definite"),
    )?;
    Ok((matrix, Some(factor)))
}

fn finite_vector(vector: Vector, field: &'static str) -> Result<(), HertzWrenchError> {
    if vector.into_iter().all(f64::is_finite) {
        Ok(())
    } else {
        Err(HertzWrenchError::Invalid(field))
    }
}

fn dot(a: Vector, b: Vector) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn sub(a: Vector, b: Vector) -> Vector {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn max_abs(a: Vector) -> f64 {
    a.into_iter().map(f64::abs).fold(0.0, f64::max)
}
fn multiply(a: Matrix, b: Vector) -> Vector {
    a.map(|row| dot(row, b))
}

fn residual_metric(residual: Vector, matrix: Matrix) -> f64 {
    let scale = matrix.into_iter().map(max_abs).fold(0.0, f64::max);
    dot(
        residual,
        multiply(matrix.map(|row| row.map(|v| v / scale)), residual),
    )
}
