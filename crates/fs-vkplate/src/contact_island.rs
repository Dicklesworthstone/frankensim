//! Simultaneous normal contact through the frozen body compliance operator.
//!
//! With signed contact footprints `N`, the caller supplies `C = N^T A^-1 N`.
//! All laws see the same final penetration: `d + C phi(d) = d_pred`. The
//! constitutive secants and body inverse remain owned by the caller. In
//! particular, an initially open contact must remain in the system: another
//! contact can close it through an off-diagonal response.

use std::cell::{Cell, RefCell};

use fs_solver::{CheckedCgConfig, CheckedCgError, LinearOp, checked_cg};
use fs_sparse::precond::Precond;

/// A candidate the rounded-pose audit refuses is solved again, one eighth
/// tighter each time. The retries change nothing about what the audit accepts,
/// and every retry spends the same Newton budget.
pub(crate) const ROUNDED_POSE_RETRIES: usize = 3;
pub(crate) const RETRY_TIGHTENING: f64 = 0.125;

/// Which independent check refused a candidate. The core audits the response
/// it solved; the body adapter repeats that audit at the rounded nodal pose.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RefusalSite {
    CoreAudit,
    RoundedPose,
}

/// Accuracy and bounded work for one simultaneous contact solve.
#[derive(Clone, Copy, Debug)]
pub struct Controls {
    /// Relative tolerance for both penetration and constitutive force.
    pub relative_tolerance: f64,
    /// Absolute penetration tolerance, in metres.
    pub displacement_tolerance: f64,
    /// Absolute force tolerance, in newtons.
    pub force_tolerance: f64,
    /// Coupled Newton steps per response component. Independent scalar laws may use the
    /// existing fixed 60-iteration safeguarded solve at the default budget.
    /// Reducing this budget disables that scalar fast path.
    pub max_newton_iterations: usize,
    /// CG iteration budget per Newton step, including defect corrections.
    pub max_cg_iterations: usize,
    pub max_backtracks: usize,
}

impl Default for Controls {
    fn default() -> Self {
        Self {
            relative_tolerance: 1e-10,
            displacement_tolerance: 1e-14,
            force_tolerance: 1e-12,
            max_newton_iterations: 32,
            max_cg_iterations: 256,
            max_backtracks: 24,
        }
    }
}

impl Controls {
    pub(crate) fn validate(self) -> Result<(), SolveError> {
        if !(self.relative_tolerance.is_finite()
            && self.relative_tolerance > 0.0
            && self.relative_tolerance < 1.0
            && self.displacement_tolerance.is_finite()
            && self.displacement_tolerance > 0.0
            && self.force_tolerance.is_finite()
            && self.force_tolerance > 0.0)
            || self.max_cg_iterations == 0
        {
            return Err(SolveError::InvalidInput(
                "invalid accuracy or CG work limit",
            ));
        }
        Ok(())
    }
}

/// No approximate force is published when the solve is refused.
#[derive(Clone, Debug, PartialEq)]
pub enum SolveError {
    /// Failure of the joint compliant normal/Coulomb stage. No contact state
    /// or nodal correction from that stage has been published.
    Friction(crate::friction_island::SolveError),
    InvalidInput(&'static str),
    InvalidLaw {
        contact: usize,
    },
    InvalidResponse,
    LinearSolve {
        iterations: usize,
        residual: f64,
    },
    NoConvergence {
        site: RefusalSite,
        newton_iterations: usize,
        displacement_residual: f64,
        force_residual: f64,
    },
}

impl SolveError {
    /// Count the Newton steps of refused candidates in a later refusal.
    pub(crate) fn after_newton(mut self, spent: usize) -> Self {
        if let Self::NoConvergence {
            newton_iterations, ..
        } = &mut self
        {
            *newton_iterations += spent;
        }
        self
    }
}

impl std::fmt::Display for SolveError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Friction(error) => error.fmt(f),
            Self::InvalidInput(reason) => write!(f, "normal contact: {reason}"),
            Self::InvalidLaw { contact } => write!(
                f,
                "normal contact {contact}: force and slope must be finite and nonnegative"
            ),
            Self::InvalidResponse => write!(f, "normal contact: non-finite body response"),
            Self::LinearSolve {
                iterations,
                residual,
            } => write!(
                f,
                "normal contact: CG failed after {iterations} iterations (residual {residual:e})"
            ),
            Self::NoConvergence {
                site,
                newton_iterations,
                displacement_residual,
                force_residual,
            } => write!(
                f,
                "normal contact: refused at {site:?} after {newton_iterations} Newton steps (scaled penetration {displacement_residual:e}, force {force_residual:e})"
            ),
        }
    }
}

impl std::error::Error for SolveError {}

/// Residuals are maxima of component errors divided by their declared
/// absolute-plus-relative tolerances. Both must be at most one to succeed.
#[derive(Clone, Copy, Debug, Default)]
pub struct SolveReport {
    /// Coupled Newton steps summed over components; zero when scalar laws suffice.
    pub newton_iterations: usize,
    pub cg_iterations: usize,
    pub compliance_applications: usize,
    pub displacement_residual: f64,
    pub force_residual: f64,
    /// Candidates the core accepted and the rounded-pose audit then refused
    /// before this one was published. Their Newton work is in the counts above.
    pub rounded_pose_retries: usize,
}

/// Linear storage in the number of contacts. Vectors are retained between
/// samples; the shared checked-CG driver still owns its Krylov allocations.
#[derive(Clone, Default)]
pub(crate) struct Workspace {
    accepted: bool,
    d: Vec<f64>,
    force: Vec<f64>,
    slope: Vec<f64>,
    residual: Vec<f64>,
    response: Vec<f64>,
    scale: Vec<f64>,
    step: Vec<f64>,
    root: Vec<f64>,
    rhs: Vec<f64>,
    inverse_diagonal: Vec<f64>,
    op_input: Vec<f64>,
    op_response: Vec<f64>,
    trial_d: Vec<f64>,
    trial_force: Vec<f64>,
    trial_slope: Vec<f64>,
    trial_response: Vec<f64>,
}

struct Compliance<F> {
    apply: RefCell<F>,
    applications: Cell<usize>,
    failed: Cell<bool>,
}

impl<F: FnMut(&[f64], &mut [f64])> Compliance<F> {
    fn apply(&self, x: &[f64], y: &mut [f64]) -> Result<(), SolveError> {
        if !x.iter().all(|v| v.is_finite()) {
            self.failed.set(true);
            return Err(SolveError::InvalidResponse);
        }
        if x.iter().all(|&v| v == 0.0) {
            y.fill(0.0);
            return Ok(());
        }
        (self.apply.borrow_mut())(x, y);
        self.applications.set(self.applications.get() + 1);
        if !y.iter().all(|v| v.is_finite()) {
            self.failed.set(true);
            return Err(SolveError::InvalidResponse);
        }
        Ok(())
    }
}

// I + sqrt(P) C sqrt(P) is SPD even with zero slopes, redundant contact
// footprints, or dark body modes. Never divide by a constitutive slope or
// invert the (possibly singular) contact compliance.
struct NewtonOp<'a, F> {
    compliance: &'a Compliance<F>,
    root: &'a [f64],
    scratch: RefCell<(&'a mut [f64], &'a mut [f64])>,
}

impl<F: FnMut(&[f64], &mut [f64])> LinearOp for NewtonOp<'_, F> {
    fn n(&self) -> usize {
        self.root.len()
    }

    fn apply(&self, x: &[f64], y: &mut [f64]) {
        let mut scratch = self.scratch.borrow_mut();
        let (input, response) = &mut *scratch;
        for i in 0..x.len() {
            input[i] = self.root[i] * x[i];
        }
        if self.compliance.apply(input, response).is_err() {
            y.fill(f64::NAN);
            return;
        }
        for i in 0..x.len() {
            y[i] = self.root[i].mul_add(response[i], x[i]);
            if !y[i].is_finite() {
                self.compliance.failed.set(true);
            }
        }
    }

    // I + sqrt(P) C sqrt(P) with a symmetric compliance C: symmetric.
    fn apply_transpose(&self, x: &[f64], y: &mut [f64]) {
        self.apply(x, y);
    }
}

struct Diagonal<'a>(&'a [f64]);

impl Precond for Diagonal<'_> {
    fn apply(&self, r: &[f64], z: &mut [f64]) {
        for i in 0..r.len() {
            z[i] = self.0[i] * r[i];
        }
    }
}

fn law_value(
    law: &mut impl FnMut(usize, f64) -> (f64, f64),
    contact: usize,
    d: f64,
) -> Result<(f64, f64), SolveError> {
    if !d.is_finite() {
        return Err(SolveError::InvalidResponse);
    }
    let (force, slope) = law(contact, d);
    if !force.is_finite() || force < 0.0 || !slope.is_finite() || slope < 0.0 {
        return Err(SolveError::InvalidLaw { contact });
    }
    Ok((force, slope))
}

impl Workspace {
    /// Empty after a refusal; valid only following a successful `solve`.
    pub(crate) fn forces(&self) -> &[f64] {
        if self.accepted { &self.force } else { &[] }
    }

    /// Actual `d_pred - C f`, independently checked against every force law.
    pub(crate) fn penetrations(&self) -> &[f64] {
        if self.accepted { &self.response } else { &[] }
    }

    /// Solve all supplied contacts simultaneously. `apply_compliance` must be
    /// one fixed symmetric positive-semidefinite operator for the whole call.
    /// `evaluate_law` must be a frozen, monotone force law and its derivative;
    /// history variables must not advance during trial evaluations.
    ///
    /// The diagonal bound is a nonnegative approximation to `diag(C)` used
    /// only for preconditioning. A dry-body upper bound is sufficient. Signed
    /// repeated footprint entries must be combined before forming that bound.
    /// This method changes scratch only; the caller publishes the forces and
    /// corrected body state together after success.
    /// `acceptance` is the largest scaled residual, at most one, this call may
    /// return. It tightens the iteration only; `Controls` still defines the
    /// laws' tolerance, which the rounded-pose audit applies at one.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn solve(
        &mut self,
        predicted: &[f64],
        previous_force: &[f64],
        compliance_diagonal_bound: &[f64],
        controls: Controls,
        acceptance: f64,
        apply_compliance: impl FnMut(&[f64], &mut [f64]),
        mut evaluate_law: impl FnMut(usize, f64) -> (f64, f64),
    ) -> Result<SolveReport, SolveError> {
        self.accepted = false;
        let n = predicted.len();
        if previous_force.len() != n || compliance_diagonal_bound.len() != n {
            return Err(SolveError::InvalidInput("contact array lengths differ"));
        }
        if !predicted.iter().all(|x| x.is_finite())
            || !previous_force
                .iter()
                .chain(compliance_diagonal_bound)
                .all(|x| x.is_finite() && *x >= 0.0)
        {
            return Err(SolveError::InvalidInput(
                "invalid prediction, warm force, or compliance bound",
            ));
        }
        controls.validate()?;
        if !(acceptance > 0.0 && acceptance <= 1.0) {
            return Err(SolveError::InvalidInput("acceptance must lie in (0, 1]"));
        }
        for v in [
            &mut self.d,
            &mut self.force,
            &mut self.slope,
            &mut self.residual,
            &mut self.response,
            &mut self.scale,
            &mut self.step,
            &mut self.root,
            &mut self.rhs,
            &mut self.inverse_diagonal,
            &mut self.op_input,
            &mut self.op_response,
            &mut self.trial_d,
            &mut self.trial_force,
            &mut self.trial_slope,
            &mut self.trial_response,
        ] {
            v.resize(n, 0.0);
        }
        if n == 0 {
            self.accepted = true;
            return Ok(SolveReport::default());
        }
        let compliance = Compliance {
            apply: RefCell::new(apply_compliance),
            applications: Cell::new(0),
            failed: Cell::new(false),
        };
        // Preserve the useful previous force, while including all off-diagonal
        // displacements in its warm-start pose.
        compliance.apply(previous_force, &mut self.response)?;
        for i in 0..n {
            self.d[i] = predicted[i] - self.response[i];
            (self.force[i], self.slope[i]) = law_value(&mut evaluate_law, i, self.d[i])?;
        }
        compliance.apply(&self.force, &mut self.response)?;
        let mut report = SolveReport::default();
        loop {
            report.displacement_residual = 0.0;
            report.force_residual = 0.0;
            for i in 0..n {
                let actual = predicted[i] - self.response[i];
                self.residual[i] = self.d[i] - actual;
                self.scale[i] = controls.displacement_tolerance
                    + controls.relative_tolerance
                        * predicted[i]
                            .abs()
                            .max(self.d[i].abs())
                            .max(self.response[i].abs());
                let (actual_force, _) = law_value(&mut evaluate_law, i, actual)?;
                let force_scale = controls.force_tolerance
                    + controls.relative_tolerance * self.force[i].abs().max(actual_force.abs());
                if !self.residual[i].is_finite()
                    || !self.scale[i].is_finite()
                    || !force_scale.is_finite()
                {
                    return Err(SolveError::InvalidResponse);
                }
                report.displacement_residual = report
                    .displacement_residual
                    .max(self.residual[i].abs() / self.scale[i]);
                report.force_residual = report
                    .force_residual
                    .max((self.force[i] - actual_force).abs() / force_scale);
            }
            report.compliance_applications = compliance.applications.get();
            if report.displacement_residual <= acceptance && report.force_residual <= acceptance {
                for i in 0..n {
                    self.response[i] = predicted[i] - self.response[i];
                }
                self.accepted = true;
                return Ok(report);
            }
            let exhausted = || SolveError::NoConvergence {
                site: RefusalSite::CoreAudit,
                newton_iterations: report.newton_iterations,
                displacement_residual: report.displacement_residual,
                force_residual: report.force_residual,
            };
            if report.newton_iterations == controls.max_newton_iterations {
                return Err(exhausted());
            }
            for i in 0..n {
                self.root[i] = self.slope[i].sqrt();
                self.rhs[i] = self.root[i] * self.residual[i];
                let diagonal = self.slope[i].mul_add(compliance_diagonal_bound[i], 1.0);
                if !diagonal.is_finite() || !self.rhs[i].is_finite() {
                    return Err(SolveError::InvalidResponse);
                }
                self.inverse_diagonal[i] = diagonal.recip();
            }
            let op = NewtonOp {
                compliance: &compliance,
                root: &self.root,
                scratch: RefCell::new((&mut self.op_input, &mut self.op_response)),
            };
            let linear = checked_cg(
                &op,
                &Diagonal(&self.inverse_diagonal),
                &self.rhs,
                CheckedCgConfig {
                    tolerance: (0.1 * controls.relative_tolerance).clamp(32.0 * f64::EPSILON, 1e-3),
                    max_iterations: controls.max_cg_iterations,
                    max_corrections: 2,
                },
                |_| {
                    if compliance.failed.get() {
                        Err(SolveError::InvalidResponse)
                    } else {
                        Ok(())
                    }
                },
            )
            .map_err(|error| match error {
                CheckedCgError::InvalidInput(reason) => SolveError::InvalidInput(reason),
                CheckedCgError::Interrupted(error) => error,
            })?;
            report.cg_iterations += linear.report.iters;
            if !linear.report.converged_euclidean() {
                return Err(SolveError::LinearSolve {
                    iterations: linear.report.iters,
                    residual: linear.report.rel_residual,
                });
            }
            for i in 0..n {
                self.op_input[i] = self.root[i] * linear.x[i];
            }
            compliance.apply(&self.op_input, &mut self.step)?;
            for i in 0..n {
                self.step[i] -= self.residual[i];
            }
            let mut length = 1.0_f64;
            let mut improved = false;
            for _ in 0..=controls.max_backtracks {
                for i in 0..n {
                    self.trial_d[i] = length.mul_add(self.step[i], self.d[i]);
                    (self.trial_force[i], self.trial_slope[i]) =
                        law_value(&mut evaluate_law, i, self.trial_d[i])?;
                }
                compliance.apply(&self.trial_force, &mut self.trial_response)?;
                // Hold the scaling fixed throughout this line search.
                let mut merit = 0.0_f64;
                for i in 0..n {
                    let r = self.trial_d[i] - (predicted[i] - self.trial_response[i]);
                    merit = merit.max(r.abs() / self.scale[i]);
                }
                if merit <= report.displacement_residual * (1.0 - 1e-4 * length) {
                    improved = true;
                    break;
                }
                length *= 0.5;
            }
            if !improved {
                return Err(exhausted());
            }
            std::mem::swap(&mut self.d, &mut self.trial_d);
            std::mem::swap(&mut self.force, &mut self.trial_force);
            std::mem::swap(&mut self.slope, &mut self.trial_slope);
            std::mem::swap(&mut self.response, &mut self.trial_response);
            report.newton_iterations += 1;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn apply_dense(c: [[f64; 2]; 2], x: &[f64], y: &mut [f64]) {
        for i in 0..2 {
            y[i] = c[i][0] * x[0] + c[i][1] * x[1];
        }
    }

    #[test]
    fn stiff_coupled_laws_match_independent_dense_elimination() {
        let c = [[2e-6, -1e-6], [-1e-6, 2e-6]];
        let stiffness = [1e10, 3e10];
        let predicted = [0.1, 0.3];
        let mut workspace = Workspace::default();
        let report = workspace
            .solve(
                &predicted,
                &[0.0; 2],
                &[2e-6; 2],
                Controls::default(),
                1.0,
                |x, y| apply_dense(c, x, y),
                |i, d| {
                    if d > 0.0 {
                        (stiffness[i] * d, stiffness[i])
                    } else {
                        (0.0, 0.0)
                    }
                },
            )
            .unwrap();
        // Eliminate the two force unknowns in (diag(1/k) + C) f = d_pred.
        let a = c[0][0] + stiffness[0].recip();
        let b = c[1][1] + stiffness[1].recip();
        let determinant = a * b - c[0][1] * c[1][0];
        let expected = [
            (b * predicted[0] - c[0][1] * predicted[1]) / determinant,
            (a * predicted[1] - c[1][0] * predicted[0]) / determinant,
        ];
        for i in 0..2 {
            assert!((workspace.forces()[i] / expected[i] - 1.0).abs() < 1e-10);
            assert!(
                (stiffness[i] * workspace.penetrations()[i] / workspace.forces()[i] - 1.0).abs()
                    < 1e-10
            );
        }
        assert!(report.force_residual <= 1.0 && report.displacement_residual <= 1.0);
    }

    #[test]
    fn another_contact_can_close_an_initially_open_gap() {
        let c = [[1.0, -0.5], [-0.5, 1.0]];
        let mut workspace = Workspace::default();
        workspace
            .solve(
                &[1.0, -0.2],
                &[0.0; 2],
                &[1.0; 2],
                Controls::default(),
                1.0,
                |x, y| apply_dense(c, x, y),
                |_, d| if d > 0.0 { (d, 1.0) } else { (0.0, 0.0) },
            )
            .unwrap();
        assert!((workspace.forces()[0] - 38.0 / 75.0).abs() < 1e-11);
        assert!((workspace.forces()[1] - 2.0 / 75.0).abs() < 1e-11);
        assert!(workspace.penetrations()[1] > 0.0);
    }

    #[test]
    fn singular_compliance_is_admitted_but_work_exhaustion_publishes_nothing() {
        let c = [[1.0; 2]; 2];
        let mut workspace = Workspace::default();
        workspace
            .solve(
                &[1.0; 2],
                &[0.2; 2],
                &[1.0; 2],
                Controls::default(),
                1.0,
                |x, y| apply_dense(c, x, y),
                |_, d| (d.max(0.0), if d > 0.0 { 1.0 } else { 0.0 }),
            )
            .unwrap();
        for &force in workspace.forces() {
            assert!((force - 1.0 / 3.0).abs() < 1e-11);
        }
        let refused = workspace.solve(
            &[1.0; 2],
            &[0.0; 2],
            &[1.0; 2],
            Controls {
                max_newton_iterations: 0,
                ..Controls::default()
            },
            1.0,
            |x, y| apply_dense(c, x, y),
            |_, d| (d.max(0.0), if d > 0.0 { 1.0 } else { 0.0 }),
        );
        assert!(matches!(
            refused,
            Err(SolveError::NoConvergence {
                site: RefusalSite::CoreAudit,
                ..
            })
        ));
        assert!(workspace.forces().is_empty() && workspace.penetrations().is_empty());
    }
}
