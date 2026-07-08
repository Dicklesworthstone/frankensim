//! Pseudo-arclength continuation through limit points: path-following
//! of equilibrium branches where plain load control fails
//! (snap-through/snap-back). The algorithm is GENERIC over
//! [`PathResidual`] — the conformance battery validates it against the
//! closed-form von Mises truss before trusting it on continuum
//! fixtures, and [`crate::HyperProblem`] plugs in for the real thing.
//!
//! Keller's bordered corrector: at arc state (u, λ) with predictor
//! tangent (u̇, λ̇), each Newton step solves the bordered system via two
//! tangent solves (K a = −R, K b = F_ext) — MINRES underneath, so
//! passing a limit point (singular K in the load direction) is
//! routine. Step control halves on failed correction and grows on
//! fast success. Events (limit points via λ̇ sign change, branch
//! points via stiffness-definiteness change without a λ̇ flip) are
//! recorded on the path; branch switching perturbs along the pencil's
//! null direction.
//!
//! The trajectory IS the state: [`PathState`] is a plain cloneable
//! value — checkpoint by clone, resume by handing it back. The battery
//! proves split runs are bit-identical to straight runs.

use crate::SolidError;
use fs_solver::krylov::MinresState;
use fs_solver::op::CsrOp;
use fs_sparse::Csr;
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
    /// Corrector residual gate.
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

fn minres_solve(k: &Csr, b: &[f64]) -> Result<Vec<f64>, SolidError> {
    let n = k.nrows();
    if n <= 1024 {
        // Fixture scale: one pivoted dense LU per solve —
        // indefinite-capable (limit points!), deterministic, and an
        // order faster than unpreconditioned Krylov on thin-structure
        // tangents. Production scale falls through to MINRES.
        if let Ok(f) = fs_la::factor::lu(&k.to_dense(), n) {
            let mut x = b.to_vec();
            f.solve(&mut x);
            // Residual gate (LU on a singular tangent can be junk).
            let mut r = vec![0.0f64; n];
            k.spmv(&x, &mut r);
            let bn: f64 = b.iter().map(|v| v * v).sum::<f64>().sqrt();
            let rn: f64 = r
                .iter()
                .zip(b)
                .map(|(a, v)| (a - v) * (a - v))
                .sum::<f64>()
                .sqrt();
            if rn <= 1e-6 * bn.max(1e-300) {
                return Ok(x);
            }
        }
    }
    let op = CsrOp::symmetric(k.clone());
    let mut st = MinresState::new(&op, b);
    let _ = st.run(&op, 1e-11, 20_000);
    if !st.rel_residual().is_finite() || st.rel_residual() > 1e-6 {
        return Err(SolidError::SolveFailed {
            iters: st.iters,
            rel_residual: st.rel_residual(),
        });
    }
    Ok(st.x)
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
/// solver/material refusals propagate.
pub fn advance<P: PathResidual>(
    problem: &P,
    state: &mut PathState,
    settings: &ArcSettings,
    steps: usize,
) -> Result<(), SolidError> {
    let n = problem.ndof();
    let fext = problem.load_vector();
    for _ in 0..steps {
        let k = problem.tangent(&state.u, state.lambda)?;
        let (u_dot, lam_dot) = if state.pending_switch {
            // Honor the branch-switch null direction: the arc
            // constraint then PINS progress along the mode, so the
            // corrector cannot fall back to the fundamental branch.
            state.pending_switch = false;
            (state.tangent.0.clone(), state.tangent.1)
        } else {
            // Predictor tangent: K u̇ = λ̇ F_ext with arc normalization.
            let v = minres_solve(&k, &fext)?;
            let norm2: f64 = v.iter().map(|x| x * x).sum::<f64>() + 1.0;
            let mut lam_dot = 1.0 / norm2.sqrt();
            let mut u_dot: Vec<f64> = v.iter().map(|x| x * lam_dot).collect();
            // Orient along the previous tangent (path continuity).
            let dot: f64 = u_dot
                .iter()
                .zip(&state.tangent.0)
                .map(|(a, b)| a * b)
                .sum::<f64>()
                + lam_dot * state.tangent.1;
            if dot < 0.0 {
                for x in &mut u_dot {
                    *x = -*x;
                }
                lam_dot = -lam_dot;
            }
            (u_dot, lam_dot)
        };
        // Events: λ̇ sign change = limit point; definiteness change
        // without a λ̇ flip = branch-point candidate.
        if state.step > 0 && lam_dot * state.tangent.1 < 0.0 {
            state.events.push(PathEvent::LimitPoint {
                lambda: state.lambda,
                step: state.step,
            });
        }
        let definite = is_positive_definite(&k);
        if state.step > 0 && definite != state.last_definite && lam_dot * state.tangent.1 > 0.0 {
            state.events.push(PathEvent::BranchPoint {
                lambda: state.lambda,
                step: state.step,
            });
        }
        state.last_definite = definite;
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
                    state.u = u1;
                    state.lambda = l1;
                    state.step += 1;
                    let defl = state
                        .u
                        .iter()
                        .fold(0.0f64, |m, &x| m.max(x.abs()));
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
                    ds *= 0.25;
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
    for it in 0..settings.max_corrector_iters {
        let r = problem.residual(&u, lam)?;
        let rn: f64 = r.iter().map(|x| x * x).sum::<f64>().sqrt();
        if rn < settings.tol {
            return Ok((u, lam, it));
        }
        let k = problem.tangent(&u, lam)?;
        // Bordered solve (Keller): K a = −R, K b = F_ext; the arc
        // constraint u̇ᵀΔu + λ̇Δλ = 0 pins the mixture.
        let neg_r: Vec<f64> = r.iter().map(|x| -x).collect();
        let a = minres_solve(&k, &neg_r)?;
        let b = minres_solve(&k, fext)?;
        let num: f64 = u_dot.iter().zip(&a).map(|(t, x)| t * x).sum();
        let den: f64 = u_dot.iter().zip(&b).map(|(t, x)| t * x).sum::<f64>() + lam_dot;
        let dlam = -num / den;
        for ((ui, ai), bi) in u.iter_mut().zip(&a).zip(&b) {
            *ui += ai + dlam * bi;
        }
        lam += dlam;
    }
    Err(SolidError::NewtonStalled {
        history: vec![f64::NAN],
    })
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
