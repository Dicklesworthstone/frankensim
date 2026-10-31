//! Pseudo-arclength continuation through limit points: path-following
//! of equilibrium branches where plain load control fails
//! (snap-through/snap-back). The algorithm is GENERIC over
//! [`PathResidual`] — the conformance battery validates it against the
//! closed-form von Mises truss before trusting it on continuum
//! fixtures, and [`crate::HyperProblem`] plugs in for the real thing.
//!
//! Keller's bordered corrector: at arc state (u, λ) with predictor
//! tangent (u̇, λ̇), each Newton step solves the full bordered system
//! `[K, -F_ext; u̇ᵀ, λ̇]`. Unlike a Schur reduction using K⁻¹, this stays
//! regular at a simple fold even when K is singular. Small systems use
//! pivoted LU; larger systems use matrix-free FGMRES on the generally
//! nonsymmetric border. Step control reduces on failed correction and grows on
//! fast success. Events (limit points via λ̇ sign change, branch
//! points via stiffness-definiteness change without a λ̇ flip) are
//! recorded on the path; branch switching perturbs along the pencil's
//! null direction.
//!
//! The trajectory IS the state: [`PathState`] is a plain cloneable
//! value — checkpoint by clone, resume by handing it back. The battery
//! proves split runs are bit-identical to straight runs.

use crate::SolidError;
use fs_solver::nonlinear::FgmresState;
use fs_solver::op::LinearOp;
use fs_sparse::Csr;
use fs_sparse::precond::IdentityPrecond;
use std::fmt::Write as _;

/// A parameterized equilibrium residual `R(u, λ) = R_int(u) − λ F_ext`.
pub trait PathResidual {
    /// DOF count.
    fn ndof(&self) -> usize;
    /// The residual at (u, λ).
    ///
    /// # Errors
    /// Material/state refusals as [`SolidError`].
    fn residual(&self, u: &[f64], lambda: f64) -> Result<Vec<f64>, SolidError>;
    /// The tangent stiffness at (u, λ).
    ///
    /// # Errors
    /// Material/state refusals as [`SolidError`].
    fn tangent(&self, u: &[f64], lambda: f64) -> Result<Csr, SolidError>;
    /// `−∂R/∂λ` (the reference load vector, zero on pinned DOFs).
    fn load_vector(&self) -> Vec<f64>;
}

/// Continuation controls.
#[derive(Debug, Clone, Copy)]
pub struct ArcSettings {
    /// Initial arclength step.
    pub ds: f64,
    /// Smallest step before declaring a stall.
    pub ds_min: f64,
    /// Largest step.
    pub ds_max: f64,
    /// Newton iterations per corrector.
    pub max_corrector_iters: usize,
    /// Absolute corrector gate, applied separately to the force-residual
    /// Euclidean norm and the tangent-plane arclength defect. The caller
    /// must use a consistent displacement/load-factor coordinate scaling.
    pub tol: f64,
}

impl Default for ArcSettings {
    fn default() -> Self {
        ArcSettings {
            ds: 0.05,
            ds_min: 1e-6,
            ds_max: 0.5,
            max_corrector_iters: 12,
            tol: 1e-9,
        }
    }
}

/// A path event (evidence rows).
#[derive(Debug, Clone, PartialEq)]
pub enum PathEvent {
    /// λ̇ changed sign: a limit point (snap-through/back) between the
    /// two recorded steps.
    LimitPoint {
        /// Load at detection.
        lambda: f64,
        /// Path step index.
        step: usize,
    },
    /// The tangent stiffness lost positive definiteness without a λ̇
    /// flip: a bifurcation (branch point) candidate.
    BranchPoint {
        /// Load at detection.
        lambda: f64,
        /// Path step index.
        step: usize,
    },
}

/// The checkpointable trajectory: clone to checkpoint, hand back to
/// resume. Plain data — no hidden solver state.
#[derive(Debug, Clone)]
pub struct PathState {
    /// Current displacement.
    pub u: Vec<f64>,
    /// Current load factor.
    pub lambda: f64,
    /// Previous predictor tangent (u̇, λ̇), unit arc norm.
    pub tangent: (Vec<f64>, f64),
    /// Current arclength step.
    pub ds: f64,
    /// Steps taken so far.
    pub step: usize,
    /// Recorded events.
    pub events: Vec<PathEvent>,
    /// (λ, ‖u‖∞ proxy) trace rows.
    pub trace: Vec<(f64, f64)>,
    /// Sign of the stiffness-definiteness indicator at the last step.
    pub last_definite: bool,
    /// Set by [`switch_branch`]: the next predictor honors the stored
    /// tangent (the bifurcation null direction) instead of recomputing
    /// the load-direction tangent.
    pub pending_switch: bool,
}

impl PathState {
    /// A fresh path at the unloaded state.
    #[must_use]
    pub fn start(ndof: usize, settings: &ArcSettings) -> PathState {
        PathState {
            u: vec![0.0; ndof],
            lambda: 0.0,
            tangent: (vec![0.0; ndof], 1.0),
            ds: settings.ds,
            step: 0,
            events: Vec::new(),
            trace: vec![(0.0, 0.0)],
            last_definite: true,
            pending_switch: false,
        }
    }

    /// Evidence rows (ledger-style JSON) for the trace and events.
    #[must_use]
    pub fn trace_json(&self) -> String {
        let mut s = String::from("[");
        for (i, (l, d)) in self.trace.iter().enumerate() {
            if i > 0 {
                s.push(',');
            }
            let _ = write!(s, "{{\"lambda\":{l:.6e},\"defl\":{d:.6e}}}");
        }
        s.push(']');
        s
    }
}

fn norm(x: &[f64]) -> f64 {
    x.iter().fold(0.0_f64, |value, &x| value.hypot(x))
}

fn invalid(what: &str) -> SolidError {
    SolidError::InvalidInput { what: what.into() }
}

/// Row equilibration changes neither the equilibrium equation nor the arc
/// constraint. It avoids hiding the arc row behind a large force unit scale.
struct Bordered<'a> {
    k: &'a Csr,
    load: &'a [f64],
    tangent: &'a [f64],
    lambda_tangent: f64,
    scales: Vec<f64>,
}

impl LinearOp for Bordered<'_> {
    fn n(&self) -> usize {
        self.load.len() + 1
    }

    fn apply(&self, x: &[f64], y: &mut [f64]) {
        let n = self.load.len();
        for (i, value) in y[..n].iter_mut().enumerate() {
            let (columns, values) = self.k.row(i);
            *value = columns
                .iter()
                .zip(values)
                .map(|(&j, &entry)| (entry / self.scales[i]) * x[j])
                .sum::<f64>()
                - (self.load[i] / self.scales[i]) * x[n];
        }
        y[n] = (self.tangent.iter().zip(x).map(|(t, v)| t * v).sum::<f64>()
            + self.lambda_tangent * x[n])
            / self.scales[n];
    }

    fn apply_transpose(&self, x: &[f64], y: &mut [f64]) {
        let n = self.load.len();
        y.fill(0.0);
        for (i, value) in x[..n].iter().enumerate() {
            let (columns, values) = self.k.row(i);
            for (&j, &entry) in columns.iter().zip(values) {
                y[j] += (entry / self.scales[i]) * value;
            }
            y[n] -= (self.load[i] / self.scales[i]) * value;
        }
        let last = x[n] / self.scales[n];
        for (value, t) in y[..n].iter_mut().zip(self.tangent) {
            *value += t * last;
        }
        y[n] += self.lambda_tangent * last;
    }
}

fn bordered_solve(
    k: &Csr,
    load: &[f64],
    tangent: &[f64],
    lambda_tangent: f64,
    rhs: &[f64],
    arc_rhs: f64,
) -> Result<Vec<f64>, SolidError> {
    let n = load.len();
    if k.nrows() != n || k.ncols() != n || tangent.len() != n || rhs.len() != n {
        return Err(invalid("continuation tangent/residual dimensions disagree"));
    }
    let mut scales = Vec::with_capacity(n + 1);
    for (i, force) in load.iter().enumerate() {
        let (_, values) = k.row(i);
        if values.iter().any(|x| !x.is_finite()) {
            return Err(invalid("nonfinite continuation tangent"));
        }
        scales.push(
            values
                .iter()
                .fold(force.abs(), |s, x| s.max(x.abs()))
                .max(f64::MIN_POSITIVE),
        );
    }
    scales.push(
        tangent
            .iter()
            .fold(lambda_tangent.abs(), |s, x| s.max(x.abs()))
            .max(f64::MIN_POSITIVE),
    );
    let op = Bordered {
        k,
        load,
        tangent,
        lambda_tangent,
        scales,
    };
    let b: Vec<_> = rhs
        .iter()
        .copied()
        .chain([arc_rhs])
        .zip(&op.scales)
        .map(|(value, scale)| value / scale)
        .collect();
    if b.iter().any(|x| !x.is_finite()) {
        return Err(invalid(
            "nonfinite equilibrated continuation right hand side",
        ));
    }
    let relative_residual = |x: &[f64]| {
        if x.iter().any(|value| !value.is_finite()) {
            return f64::INFINITY;
        }
        let mut residual = vec![0.0; n + 1];
        op.apply(x, &mut residual);
        for (value, rhs) in residual.iter_mut().zip(&b) {
            *value -= rhs;
        }
        norm(&residual) / norm(&b).max(f64::MIN_POSITIVE)
    };
    if n < 1024 {
        // The dimension cap applies to the entire bordered dense matrix.
        let size = n + 1;
        let mut dense = vec![0.0; size * size];
        for i in 0..n {
            let (columns, values) = k.row(i);
            for (&j, &value) in columns.iter().zip(values) {
                dense[i * size + j] = value / op.scales[i];
            }
            dense[i * size + n] = -load[i] / op.scales[i];
            dense[n * size + i] = tangent[i] / op.scales[n];
        }
        dense[n * size + n] = lambda_tangent / op.scales[n];
        if let Ok(factor) = fs_la::factor::lu(&dense, size) {
            let mut x = b.clone();
            factor.solve(&mut x);
            if relative_residual(&x) <= 1e-10 {
                return Ok(x);
            }
        }
    }
    let mut state = FgmresState::new(&b, (n + 1).min(64));
    let _ = state.run(&op, &IdentityPrecond, &b, 1e-11, 320);
    let rel_residual = relative_residual(&state.x);
    if !rel_residual.is_finite() || rel_residual > 1e-10 {
        return Err(SolidError::SolveFailed {
            iters: state.iters,
            rel_residual,
        });
    }
    Ok(state.x)
}

/// Definiteness probe: a Cholesky attempt on the (fixture-scale)
/// tangent stiffness. A path tangent Rayleigh quotient CANNOT see the
/// transverse singularity of a symmetric bifurcation — factorization
/// failure can. Gated to n ≤ 1024 (above that, detection is skipped
/// and branch events are the caller's pencil-monitoring job).
fn is_positive_definite(k: &Csr) -> bool {
    let n = k.nrows();
    if n > 1024 {
        return true;
    }
    fs_la::factor::cholesky(&k.to_dense(), n).is_ok()
}

/// Advance the path by `steps` pseudo-arclength steps.
///
/// # Errors
/// [`SolidError::NewtonStalled`] when the step controller bottoms out;
/// solver/material refusals propagate. A failed step leaves the complete
/// checkpoint unchanged; earlier accepted steps in the call remain committed.
pub fn advance<P: PathResidual>(
    problem: &P,
    state: &mut PathState,
    settings: &ArcSettings,
    steps: usize,
) -> Result<(), SolidError> {
    let n = problem.ndof();
    let tangent_length = norm(&state.tangent.0).hypot(state.tangent.1);
    if n == 0
        || state.u.len() != n
        || state.tangent.0.len() != n
        || state
            .u
            .iter()
            .chain(&state.tangent.0)
            .any(|x| !x.is_finite())
        || !state.lambda.is_finite()
        || !state.tangent.1.is_finite()
        || !tangent_length.is_finite()
        || tangent_length <= 0.0
        || !settings.ds_min.is_finite()
        || settings.ds_min <= 0.0
        || !settings.ds_max.is_finite()
        || settings.ds_max < settings.ds_min
        || !state.ds.is_finite()
        || state.ds < settings.ds_min
        || state.ds > settings.ds_max
        || !settings.tol.is_finite()
        || settings.tol <= 0.0
        || settings.max_corrector_iters == 0
        || state.step.checked_add(steps).is_none()
    {
        return Err(invalid(
            "finite dimension-matched path and positive ordered continuation controls required",
        ));
    }
    let fext = problem.load_vector();
    if fext.len() != n || fext.iter().any(|x| !x.is_finite()) {
        return Err(invalid(
            "finite dimension-matched continuation load required",
        ));
    }
    for _ in 0..steps {
        let k = problem.tangent(&state.u, state.lambda)?;
        if k.nrows() != n
            || k.ncols() != n
            || (0..n).any(|row| k.row(row).1.iter().any(|x| !x.is_finite()))
        {
            return Err(invalid(
                "finite dimension-matched continuation tangent required",
            ));
        }
        let (u_dot, lam_dot): (Vec<f64>, f64) = if state.pending_switch {
            // Honor the branch-switch null direction: the arc
            // constraint then PINS progress along the mode, so the
            // corrector cannot fall back to the fundamental branch.
            let length = norm(&state.tangent.0).hypot(state.tangent.1);
            (
                state.tangent.0.iter().map(|x| x / length).collect(),
                state.tangent.1 / length,
            )
        } else {
            // The last row orients the new tangent along the previous one.
            // No inverse of K is needed at a simple fold.
            let v = bordered_solve(
                &k,
                &fext,
                &state.tangent.0,
                state.tangent.1,
                &vec![0.0; n],
                1.0,
            )?;
            let length = norm(&v);
            if !length.is_finite() || length == 0.0 {
                return Err(invalid("unresolved continuation predictor direction"));
            }
            (v[..n].iter().map(|x| x / length).collect(), v[n] / length)
        };
        // Events: λ̇ sign change = limit point; definiteness change
        // without a λ̇ flip = branch-point candidate.
        let limit_point = state.step > 0 && lam_dot * state.tangent.1 < 0.0;
        let definite = is_positive_definite(&k);
        let branch_point =
            state.step > 0 && definite != state.last_definite && lam_dot * state.tangent.1 > 0.0;
        // Corrector with step halving.
        let mut ds = state.ds;
        loop {
            match correct(
                problem,
                &fext,
                &state.u,
                state.lambda,
                &u_dot,
                lam_dot,
                ds,
                settings,
            ) {
                Ok((u1, l1, iters)) => {
                    if limit_point {
                        state.events.push(PathEvent::LimitPoint {
                            lambda: state.lambda,
                            step: state.step,
                        });
                    }
                    if branch_point {
                        state.events.push(PathEvent::BranchPoint {
                            lambda: state.lambda,
                            step: state.step,
                        });
                    }
                    state.last_definite = definite;
                    state.pending_switch = false;
                    state.u = u1;
                    state.lambda = l1;
                    state.step += 1;
                    let defl = state.u.iter().fold(0.0f64, |m, &x| m.max(x.abs()));
                    state.trace.push((l1, defl));
                    state.tangent = (u_dot, lam_dot);
                    // Grow on fast success, up to the cap.
                    state.ds = if iters <= 3 {
                        (ds * 1.5).min(settings.ds_max)
                    } else {
                        ds
                    };
                    break;
                }
                Err(_) if ds > settings.ds_min => {
                    ds = (ds * 0.25).max(settings.ds_min);
                }
                Err(e) => return Err(e),
            }
        }
    }
    Ok(())
}

/// One bordered-Newton correction from a predictor point.
#[allow(clippy::too_many_arguments)]
fn correct<P: PathResidual>(
    problem: &P,
    fext: &[f64],
    u0: &[f64],
    l0: f64,
    u_dot: &[f64],
    lam_dot: f64,
    ds: f64,
    settings: &ArcSettings,
) -> Result<(Vec<f64>, f64, usize), SolidError> {
    let mut u: Vec<f64> = u0.iter().zip(u_dot).map(|(a, b)| a + ds * b).collect();
    let mut lam = l0 + ds * lam_dot;
    let mut history = Vec::new();
    for it in 0..=settings.max_corrector_iters {
        if u.iter().any(|x| !x.is_finite()) || !lam.is_finite() {
            return Err(invalid("nonfinite continuation trial"));
        }
        let r = problem.residual(&u, lam)?;
        if r.len() != u.len() || r.iter().any(|x| !x.is_finite()) {
            return Err(invalid(
                "finite dimension-matched continuation residual required",
            ));
        }
        let rn = norm(&r);
        let arc = u
            .iter()
            .zip(u0)
            .zip(u_dot)
            .map(|((u, u0), t)| (u - u0) * t)
            .sum::<f64>()
            + (lam - l0) * lam_dot
            - ds;
        history.push(rn);
        if rn <= settings.tol && arc.abs() <= settings.tol {
            return Ok((u, lam, it));
        }
        if it == settings.max_corrector_iters {
            break;
        }
        let k = problem.tangent(&u, lam)?;
        // Solve the complete Newton equations, including any accumulated arc
        // defect. At a simple fold this matrix is regular although K is not.
        let neg_r: Vec<f64> = r.iter().map(|x| -x).collect();
        let delta = bordered_solve(&k, fext, u_dot, lam_dot, &neg_r, -arc)?;
        for (ui, du) in u.iter_mut().zip(&delta) {
            *ui += du;
        }
        lam += delta[u.len()];
    }
    Err(SolidError::NewtonStalled { history })
}

/// Branch switching: seed a small perturbation along the buckling
/// mode and make the mode the NEXT PREDICTOR direction
/// (`pending_switch`): the arc constraint then pins the first step's
/// progress along the null direction, and the bordered corrector
/// lands on the bifurcated branch — a state-only perturbation would
/// simply relax back to the (still-existing, unstable) fundamental
/// branch, and a basin-scale jump would invert elements.
pub fn switch_branch(state: &mut PathState, mode: &[[f64; 2]], amplitude: f64) {
    let mmax = mode
        .iter()
        .fold(0.0f64, |m, v| m.max(v[0].abs()).max(v[1].abs()))
        .max(1e-30);
    let mut dir = vec![0.0f64; state.u.len()];
    for (node, m) in mode.iter().enumerate() {
        if 2 * node + 1 < state.u.len() {
            state.u[2 * node] += amplitude * m[0] / mmax;
            state.u[2 * node + 1] += amplitude * m[1] / mmax;
            dir[2 * node] = m[0] / mmax;
            dir[2 * node + 1] = m[1] / mmax;
        }
    }
    let n2: f64 = dir.iter().map(|x| x * x).sum::<f64>().sqrt().max(1e-30);
    for x in &mut dir {
        *x /= n2;
    }
    if amplitude < 0.0 {
        for x in &mut dir {
            *x = -*x;
        }
    }
    state.tangent = (dir, 0.0);
    state.pending_switch = true;
}
