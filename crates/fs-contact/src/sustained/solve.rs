//! Newton--Krylov solve of sustained forces against an original body response.
//!
//! Only the mechanical force-to-path response is numerically differenced.
//! The nonsmooth contact law keeps one analytic generalized derivative fixed at
//! each Newton point, including disk and normal-mode boundaries. The existing
//! fs-solver owner supplies Newton globalization and FGMRES; this module supplies
//! the contact equations and preserves the first physical source refusal.
//!
//! This is a local candidate solve, not an existence, uniqueness, contact-mode,
//! trajectory, or conservation certificate. The receiving mechanical owner must
//! independently audit the returned candidate before adopting any body state.

use super::{ContactLawError, ContactLinearization, ContactPoint, NaturalMapScales, linearize};
use fs_solver::{
    NewtonError, NewtonKrylovConfig, NewtonKrylovState, NewtonReport, NonlinearProblem,
    SolverRunProgress,
};
use std::cell::{Cell, RefCell};
use std::fmt;
use std::panic::{AssertUnwindSafe, catch_unwind};

#[cfg(test)]
mod tests;

/// Hard bound keeping contact-space work between cancellation polls small.
pub const MAX_CONTACTS: usize = 21;

/// Bounded numerical controls; no field changes the physical contact law.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SolveControls {
    /// Original Newton--Krylov policy, with residuals measured in N.
    pub newton: NewtonKrylovConfig,
    /// Maximum outer Newton attempts, in 1..=256.
    pub maximum_iterations: usize,
    /// Maximum original mechanical response calls, including final replay.
    pub maximum_response_evaluations: usize,
    /// Positive finite maximum-component force perturbation in N.
    ///
    /// This determines the numerical derivative of the original mechanical
    /// response. It must resolve that owner's integration error and force scale.
    pub difference_step_n: f64,
    /// Positive finite maximum component residual allowed after fresh replay.
    pub force_residual_tolerance_n: f64,
    /// Caller contact budget, in 1..=MAX_CONTACTS.
    pub maximum_contacts: usize,
    /// Conservative f64-equivalent entry budget for solver-owned storage.
    ///
    /// Includes Krylov bases, the configured restart's Hessenberg matrix,
    /// contact caches, Newton staging, and telemetry. The opaque original body
    /// response owns and must separately bound its own allocation and work.
    pub maximum_workspace_entries: usize,
}

impl Default for SolveControls {
    fn default() -> Self {
        Self {
            newton: NewtonKrylovConfig::default(),
            maximum_iterations: 64,
            maximum_response_evaluations: 16_384,
            difference_step_n: 1.0e-4,
            force_residual_tolerance_n: 1.0e-8,
            maximum_contacts: MAX_CONTACTS,
            maximum_workspace_entries: 262_144,
        }
    }
}

/// A numerical candidate and the final independently recomputed path.
#[derive(Clone, Debug, PartialEq)]
pub struct SolveReport {
    /// Applied held force per contact, in [normal, tangent 1, tangent 2] N.
    pub forces_n: Vec<[f64; 3]>,
    /// Final original response: [endpoint gap, tangent slip 1, tangent slip 2] m.
    pub paths_m: Vec<[f64; 3]>,
    /// Unmodified report of the original nonlinear solver.
    pub nonlinear: NewtonReport,
    /// Maximum absolute component of the freshly evaluated contact residual.
    pub maximum_residual_n: f64,
    /// Original mechanical response calls actually spent, including replay.
    pub response_evaluations: usize,
}

/// First source refusal or numerical failure; no Clone bound on source errors.
#[derive(Debug)]
pub enum SolveError<E> {
    /// An invalid control or input domain.
    Invalid(&'static str),
    /// Nonempty contact count exceeds the caller or hard budget.
    ContactCount {
        /// Supplied contacts.
        supplied: usize,
        /// Effective maximum.
        maximum: usize,
    },
    /// Contact-indexed inputs have different lengths.
    Dimensions {
        /// Friction entries.
        contacts: usize,
        /// Scale entries.
        scales: usize,
        /// Force triples.
        forces: usize,
    },
    /// Solver-owned storage would exceed its admitted budget.
    WorkspaceBudget {
        /// Conservative required f64-equivalent entries.
        required_entries: usize,
        /// Caller budget.
        maximum_entries: usize,
    },
    /// The next original response would exceed the caller budget.
    ResponseBudget {
        /// Calls already spent.
        evaluations: usize,
        /// Caller budget.
        maximum: usize,
    },
    /// Cancellation observed before publishing a candidate.
    Cancelled,
    /// Exact original mechanical source error.
    Response(E),
    /// Original response unwound; available only with panic=unwind.
    ResponsePanicked,
    /// Original response did not return one path triple per contact.
    ResponseShape {
        /// Contact count.
        expected: usize,
        /// Returned path triples.
        actual: usize,
    },
    /// A response, difference stencil or action is unrepresentable.
    NonFinite(&'static str),
    /// Original contact-law refusal at a specific contact.
    Law {
        /// Contact index.
        contact: usize,
        /// Exact law error.
        error: ContactLawError,
    },
    /// Exact original Newton constructor/configuration error.
    Newton(NewtonError),
    /// Original controlled solve did not finish its invocation.
    Interrupted {
        /// Native interruption receipt.
        progress: SolverRunProgress,
        /// Native last-committed-state report.
        nonlinear: NewtonReport,
    },
    /// Original solver stopped without numerical convergence.
    NotConverged(NewtonReport),
    /// Fresh final response failed the caller's force-residual gate.
    Residual {
        /// Recomputed maximum absolute residual in N.
        maximum_residual_n: f64,
        /// Caller tolerance in N.
        tolerance_n: f64,
        /// Native solver report, retained even if it claimed convergence.
        nonlinear: NewtonReport,
    },
}

impl<E: fmt::Display> fmt::Display for SolveError<E> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Invalid(reason) => write!(f, "invalid sustained-force solve: {reason}"),
            Self::ContactCount { supplied, maximum } => {
                write!(
                    f,
                    "sustained-force contact count {supplied} outside 1..={maximum}"
                )
            }
            Self::Dimensions {
                contacts,
                scales,
                forces,
            } => write!(
                f,
                "sustained-force dimensions differ: {contacts} contacts, {scales} scales, {forces} forces"
            ),
            Self::WorkspaceBudget {
                required_entries,
                maximum_entries,
            } => write!(
                f,
                "sustained-force workspace requires {required_entries} entries, budget {maximum_entries}"
            ),
            Self::ResponseBudget {
                evaluations,
                maximum,
            } => write!(
                f,
                "sustained-force response budget exhausted after {evaluations}/{maximum} calls"
            ),
            Self::Cancelled => write!(f, "sustained-force solve cancelled"),
            Self::Response(error) => write!(f, "original mechanical response: {error}"),
            Self::ResponsePanicked => write!(f, "original mechanical response panicked"),
            Self::ResponseShape { expected, actual } => write!(
                f,
                "original mechanical response returned {actual} paths for {expected} contacts"
            ),
            Self::NonFinite(operation) => {
                write!(f, "non-finite sustained-force operation: {operation}")
            }
            Self::Law { contact, error } => write!(f, "contact {contact}: {error}"),
            Self::Newton(error) => write!(f, "sustained-force Newton solver: {error}"),
            Self::Interrupted { progress, .. } => {
                write!(f, "sustained-force solver interrupted: {progress:?}")
            }
            Self::NotConverged(report) => write!(
                f,
                "sustained-force solver did not converge: {:?}",
                report.diagnosis
            ),
            Self::Residual {
                maximum_residual_n,
                tolerance_n,
                ..
            } => write!(
                f,
                "sustained-force replay residual {maximum_residual_n} N exceeds {tolerance_n} N"
            ),
        }
    }
}

impl<E: std::error::Error + 'static> std::error::Error for SolveError<E> {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Response(error) => Some(error),
            Self::Law { error, .. } => Some(error),
            Self::Newton(error) => Some(error),
            _ => None,
        }
    }
}

struct Linearization {
    point: Vec<f64>,
    paths_m: Vec<[f64; 3]>,
    blocks: Vec<ContactLinearization>,
}

struct Problem<'a, E, F, C> {
    friction_coefficients: &'a [f64],
    scales: &'a [NaturalMapScales],
    controls: SolveControls,
    response: RefCell<&'a mut F>,
    cancelled: RefCell<&'a mut C>,
    evaluations: Cell<usize>,
    failed: Cell<bool>,
    failure: RefCell<Option<SolveError<E>>>,
    cached: RefCell<Option<Linearization>>,
}

impl<E, F, C> Problem<'_, E, F, C>
where
    F: FnMut(&[[f64; 3]], &mut dyn FnMut() -> bool) -> Result<Vec<[f64; 3]>, E>,
    C: FnMut() -> bool,
{
    fn fail(&self, error: SolveError<E>) {
        if !self.failed.replace(true) {
            *self.failure.borrow_mut() = Some(error);
        }
    }

    fn take_failure(&self) -> Option<SolveError<E>> {
        self.failure.borrow_mut().take()
    }

    fn poll(&self) -> bool {
        if self.failed.get() {
            return false;
        }
        let stopped = (self.cancelled.borrow_mut())();
        if stopped {
            self.fail(SolveError::Cancelled);
            false
        } else {
            true
        }
    }

    fn checked_point(&self, x: &[f64]) -> bool {
        if x.len() != self.dimension() || x.iter().any(|value| !value.is_finite()) {
            self.fail(SolveError::NonFinite(
                "incomplete or non-finite force coordinates",
            ));
            false
        } else {
            true
        }
    }

    /// Every actual invocation is charged, including refusal, unwind and replay.
    fn evaluate_paths(&self, x: &[f64]) -> Option<Vec<[f64; 3]>> {
        if !self.poll() || !self.checked_point(x) {
            return None;
        }
        let evaluations = self.evaluations.get();
        if evaluations >= self.controls.maximum_response_evaluations {
            self.fail(SolveError::ResponseBudget {
                evaluations,
                maximum: self.controls.maximum_response_evaluations,
            });
            return None;
        }
        self.evaluations.set(evaluations + 1);
        let forces = force_triples(x);
        let returned = catch_unwind(AssertUnwindSafe(|| {
            (self.response.borrow_mut())(&forces, &mut || !self.poll())
        }));
        if self.failed.get() {
            return None;
        }
        let paths = match returned {
            Ok(Ok(paths)) => paths,
            Ok(Err(error)) => {
                self.fail(SolveError::Response(error));
                return None;
            }
            Err(_) => {
                self.fail(SolveError::ResponsePanicked);
                return None;
            }
        };
        if paths.len() != self.friction_coefficients.len() {
            self.fail(SolveError::ResponseShape {
                expected: self.friction_coefficients.len(),
                actual: paths.len(),
            });
            return None;
        }
        if paths.iter().flatten().any(|value| !value.is_finite()) {
            self.fail(SolveError::NonFinite("original mechanical path"));
            return None;
        }
        self.poll().then_some(paths)
    }

    fn blocks_at(&self, x: &[f64], paths: &[[f64; 3]]) -> Option<Vec<ContactLinearization>> {
        let mut blocks = Vec::with_capacity(self.friction_coefficients.len());
        for (contact, path) in paths.iter().enumerate() {
            if !self.poll() {
                return None;
            }
            let point = ContactPoint {
                force_n: [x[3 * contact], x[3 * contact + 1], x[3 * contact + 2]],
                gap_m: path[0],
                slip_increment_m: [path[1], path[2]],
                friction_coefficient: self.friction_coefficients[contact],
            };
            match linearize(point, self.scales[contact]) {
                Ok(block) => blocks.push(block),
                Err(error) => {
                    self.fail(SolveError::Law { contact, error });
                    return None;
                }
            }
        }
        Some(blocks)
    }

    fn ensure_linearization(&self, x: &[f64]) -> bool {
        if !self.poll() || !self.checked_point(x) {
            return false;
        }
        if self.cached.borrow().as_ref().is_some_and(|cached| {
            cached
                .point
                .iter()
                .zip(x)
                .all(|(a, b)| a.to_bits() == b.to_bits())
        }) {
            return true;
        }
        let Some(paths_m) = self.evaluate_paths(x) else {
            return false;
        };
        let Some(blocks) = self.blocks_at(x, &paths_m) else {
            return false;
        };
        *self.cached.borrow_mut() = Some(Linearization {
            point: x.to_vec(),
            paths_m,
            blocks,
        });
        true
    }
}

impl<E, F, C> NonlinearProblem for Problem<'_, E, F, C>
where
    F: FnMut(&[[f64; 3]], &mut dyn FnMut() -> bool) -> Result<Vec<[f64; 3]>, E>,
    C: FnMut() -> bool,
{
    fn dimension(&self) -> usize {
        3 * self.friction_coefficients.len()
    }

    fn residual(&self, x: &[f64], output: &mut [f64]) {
        output.fill(f64::NAN);
        if output.len() != self.dimension() {
            self.fail(SolveError::Invalid("residual output dimension"));
            return;
        }
        if !self.ensure_linearization(x) {
            return;
        }
        let cached = self.cached.borrow();
        let base = cached
            .as_ref()
            .expect("successful linearization retains its base");
        for (triple, block) in output.chunks_exact_mut(3).zip(&base.blocks) {
            triple.copy_from_slice(&block.residual_n);
        }
    }

    fn jacobian_apply(&self, x: &[f64], direction: &[f64], output: &mut [f64]) {
        output.fill(f64::NAN);
        if !self.poll() {
            return;
        }
        if direction.len() != self.dimension()
            || output.len() != self.dimension()
            || direction.iter().any(|value| !value.is_finite())
        {
            self.fail(SolveError::Invalid(
                "finite complete Jacobian directions required",
            ));
            return;
        }
        if !self.ensure_linearization(x) {
            return;
        }
        let scale = direction.iter().copied().map(f64::abs).fold(0.0, f64::max);
        if scale == 0.0 {
            output.fill(0.0);
            return;
        }
        let delta = self.controls.difference_step_n / scale;
        if !delta.is_finite() || delta <= 0.0 {
            self.fail(SolveError::NonFinite(
                "mechanical response difference scale",
            ));
            return;
        }
        let shifted: Vec<f64> = x
            .iter()
            .zip(direction)
            .map(|(&value, &change)| delta.mul_add(change, value))
            .collect();
        if shifted.iter().any(|value| !value.is_finite()) || shifted.as_slice() == x {
            self.fail(SolveError::NonFinite(
                "unchanged or non-finite force stencil",
            ));
            return;
        }
        // Difference ONLY the original mechanical path. The cache continues to
        // contain the exact base point and the same fixed contact-law blocks.
        let Some(shifted_paths) = self.evaluate_paths(&shifted) else {
            return;
        };
        let cached = self.cached.borrow();
        let base = cached
            .as_ref()
            .expect("successful linearization retains its base");
        for (contact, shifted_path) in shifted_paths.iter().enumerate() {
            if !self.poll() {
                output.fill(f64::NAN);
                return;
            }
            let force_direction = [
                direction[3 * contact],
                direction[3 * contact + 1],
                direction[3 * contact + 2],
            ];
            let mut path_direction = [0.0; 3];
            for component in 0..3 {
                path_direction[component] =
                    (shifted_path[component] - base.paths_m[contact][component]) / delta;
            }
            match base.blocks[contact].action(force_direction, path_direction) {
                Ok(value) => output[3 * contact..3 * contact + 3].copy_from_slice(&value),
                Err(error) => {
                    self.fail(SolveError::Law { contact, error });
                    output.fill(f64::NAN);
                    return;
                }
            }
        }
    }
}

fn force_triples(flat: &[f64]) -> Vec<[f64; 3]> {
    flat.chunks_exact(3).map(|v| [v[0], v[1], v[2]]).collect()
}

fn admit<E>(
    friction_coefficients: &[f64],
    scales: &[NaturalMapScales],
    initial_forces_n: &[[f64; 3]],
    controls: SolveControls,
) -> Result<(), SolveError<E>> {
    if controls.maximum_contacts == 0
        || controls.maximum_contacts > MAX_CONTACTS
        || controls.maximum_iterations == 0
        || controls.maximum_iterations > 256
        || controls.maximum_response_evaluations == 0
        || controls.newton.linear_restart > 64
        || controls.newton.max_linear_cycles > 128
        || !controls.difference_step_n.is_finite()
        || controls.difference_step_n <= 0.0
        || !controls.force_residual_tolerance_n.is_finite()
        || controls.force_residual_tolerance_n <= 0.0
    {
        return Err(SolveError::Invalid(
            "positive finite scales and bounded contact, iteration, restart and cycle budgets required",
        ));
    }
    controls.newton.validate().map_err(SolveError::Newton)?;
    let contacts = friction_coefficients.len();
    if contacts == 0 || contacts > controls.maximum_contacts {
        return Err(SolveError::ContactCount {
            supplied: contacts,
            maximum: controls.maximum_contacts,
        });
    }
    if scales.len() != contacts || initial_forces_n.len() != contacts {
        return Err(SolveError::Dimensions {
            contacts,
            scales: scales.len(),
            forces: initial_forces_n.len(),
        });
    }
    if friction_coefficients
        .iter()
        .any(|mu| !mu.is_finite() || *mu < 0.0)
        || initial_forces_n
            .iter()
            .flatten()
            .any(|force| !force.is_finite())
    {
        return Err(SolveError::Invalid(
            "finite initial forces and nonnegative finite friction required",
        ));
    }
    for (contact, scale) in scales.iter().enumerate() {
        scale
            .validate()
            .map_err(|error| SolveError::Law { contact, error })?;
    }
    // Conservative f64-equivalent entries, including vector headers, growth
    // reserve and simultaneously live Newton checkpoints/reports. FGMRES uses
    // the configured restart m even when m exceeds the contact dimension.
    // Counts are safe after the hard contact/iteration/restart/cycle admission.
    let n = 3 * contacts;
    let m = controls.newton.linear_restart;
    let required_entries = (m + 1) * m
        + 2 * (m + 1) * n
        + 64 * n
        + 16 * m
        + 64 * controls.maximum_iterations
        + 4 * (m * controls.newton.max_linear_cycles + 1)
        + 64 * contacts
        + 256;
    if required_entries > controls.maximum_workspace_entries {
        return Err(SolveError::WorkspaceBudget {
            required_entries,
            maximum_entries: controls.maximum_workspace_entries,
        });
    }
    Ok(())
}

/// Solve normal configuration complementarity and Coulomb forces jointly.
///
/// `response` must evaluate the supplied held forces against the same accepted
/// original mechanical state every time. It returns one [endpoint normal gap,
/// accumulated tangent slip 1, accumulated tangent slip 2] triple in m per
/// contact, expressed in the same fixed frame as that contact's force in N.
/// It must not mutate accepted physical history. Native material evolution,
/// inertia and all interactions belong inside this original response.
///
/// The callback receives a cancellation predicate so nested original solvers
/// can cooperate without borrowing the outer predicate twice. A source refusal
/// latches permanently and stops further original-response evaluations.
/// Allocation/work inside the callback remains its owner's responsibility.
///
/// A converged iterate is always replayed through a fresh original response,
/// then checked against `force_residual_tolerance_n`. This replay is the last
/// successful response call: a caller may retain its complete original body
/// candidates, verify exact force/path identity, and independently audit them.
/// Numerical convergence alone does not authorize physical state adoption.
///
/// # Errors
/// Refuses invalid controls/input, resource exhaustion, source errors/unwinds,
/// cancellation, non-finite derivatives, native solver interruption/stall, or a
/// final force residual outside tolerance. The first source error is preserved.
pub fn solve<E>(
    friction_coefficients: &[f64],
    scales: &[NaturalMapScales],
    initial_forces_n: &[[f64; 3]],
    controls: SolveControls,
    response: &mut impl FnMut(&[[f64; 3]], &mut dyn FnMut() -> bool) -> Result<Vec<[f64; 3]>, E>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<SolveReport, SolveError<E>> {
    admit(friction_coefficients, scales, initial_forces_n, controls)?;
    let problem = Problem {
        friction_coefficients,
        scales,
        controls,
        response: RefCell::new(response),
        cancelled: RefCell::new(cancelled),
        evaluations: Cell::new(0),
        failed: Cell::new(false),
        failure: RefCell::new(None),
        cached: RefCell::new(None),
    };
    if !problem.poll() {
        return Err(problem
            .take_failure()
            .expect("failed poll retains its cause"));
    }
    let initial: Vec<f64> = initial_forces_n.iter().flatten().copied().collect();
    let initialized = NewtonKrylovState::new(&problem, initial, controls.newton);
    if let Some(error) = problem.take_failure() {
        return Err(error);
    }
    let mut nonlinear = initialized.map_err(SolveError::Newton)?;
    let run = nonlinear.run_controlled(&problem, controls.maximum_iterations, &mut || {
        !problem.poll()
    });
    if let Some(error) = problem.take_failure() {
        return Err(error);
    }
    if run.progress != SolverRunProgress::Complete {
        return Err(SolveError::Interrupted {
            progress: run.progress,
            nonlinear: run.report,
        });
    }
    if !run.report.converged {
        return Err(SolveError::NotConverged(run.report));
    }
    let paths_m = problem.evaluate_paths(&nonlinear.x);
    if let Some(error) = problem.take_failure() {
        return Err(error);
    }
    let paths_m = paths_m.expect("successful response retains its paths");
    let blocks = problem.blocks_at(&nonlinear.x, &paths_m);
    if let Some(error) = problem.take_failure() {
        return Err(error);
    }
    let maximum_residual_n = blocks
        .expect("successful law evaluation retains its blocks")
        .iter()
        .flat_map(|block| &block.residual_n)
        .copied()
        .map(f64::abs)
        .fold(0.0, f64::max);
    if maximum_residual_n > controls.force_residual_tolerance_n {
        return Err(SolveError::Residual {
            maximum_residual_n,
            tolerance_n: controls.force_residual_tolerance_n,
            nonlinear: run.report,
        });
    }
    if !problem.poll() {
        return Err(problem
            .take_failure()
            .expect("failed poll retains its cause"));
    }
    Ok(SolveReport {
        forces_n: force_triples(&nonlinear.x),
        paths_m,
        nonlinear: run.report,
        maximum_residual_n,
        response_evaluations: problem.evaluations.get(),
    })
}
