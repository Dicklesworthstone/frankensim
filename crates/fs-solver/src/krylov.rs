//! Resumable Krylov methods: CG (SPD), MINRES (symmetric indefinite),
//! GMRES(m) (general — and the transposed-solve workhorse). Every
//! state is plain data: `clone()` is a checkpoint, and interrupting
//! between `run` calls loses nothing (CG/MINRES resume per ITERATION,
//! GMRES per RESTART CYCLE — its Arnoldi basis is deliberately not
//! checkpointed mid-cycle). Inner products are the crate's
//! deterministic fixed-shape reduction.

use crate::op::LinearOp;
use crate::{dot, norm2};
use fs_sparse::precond::Precond;

/// Why a solve stopped short (the structured alternative to a
/// timeout mystery).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StallDiagnosis {
    /// Residual plateaued: it stayed within ±5% of the window start over
    /// the last 50 iterations (flat, not diverging) — preconditioner
    /// quality is the first suspect.
    Plateau,
    /// Residual still falling when the budget ran out — raise the
    /// budget or improve the preconditioner.
    BudgetExhausted,
    /// A non-finite quantity appeared (breakdown, or an indefinite
    /// operator handed to CG).
    Breakdown,
}

/// WHAT a reported relative residual is, and in which norm — the
/// provenance `SolveReport::rel_residual` cannot carry on its own.
///
/// The five solvers in this crate do not all report the same quantity,
/// so `converged == true` does not by itself mean
/// `‖b − Ax‖₂/‖b‖₂ < tol`. Ask the producing state for its claim
/// ([`CgState::residual_claim`] and siblings) before treating a report
/// as a Euclidean correctness statement.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ResidualClaim {
    /// `‖b − Ax‖₂/‖b‖₂`, recomputed from the current iterate by an
    /// explicit operator application. The strongest claim: it is the
    /// Euclidean relative residual, not an estimate of it.
    /// (GMRES(m), FGMRES.)
    TrueEuclidean(f64),
    /// A recurrence's estimate of `‖b − Ax‖₂/‖b‖₂`, propagated rather
    /// than recomputed (CG's updated `r`, MINRES's Givens `η`). Equal
    /// to the true residual in exact arithmetic; in floating point it
    /// DRIFTS, and the drift is one-sided in the dangerous direction —
    /// the estimate keeps falling after the true residual stagnates.
    RecursiveEstimate(f64),
    /// A recurrence's estimate of the residual in the PRECONDITIONER's
    /// M-norm, `‖r‖_M/‖b‖_M`, not the Euclidean norm (P-MINRES). For an
    /// ill-conditioned M the two differ by orders of magnitude: with
    /// `M = diag(1, 1e-12)`, `r = (0,1)` against `b = (1,0)` gives
    /// `‖r‖_M/‖b‖_M = 1e-6` while `‖r‖₂/‖b‖₂ = 1`.
    MNormEstimate(f64),
}

impl ResidualClaim {
    /// The reported number, whatever it means.
    #[must_use]
    pub fn value(&self) -> f64 {
        match self {
            ResidualClaim::TrueEuclidean(v)
            | ResidualClaim::RecursiveEstimate(v)
            | ResidualClaim::MNormEstimate(v) => *v,
        }
    }

    /// True only for a recomputed Euclidean relative residual.
    #[must_use]
    pub fn is_true_euclidean(&self) -> bool {
        matches!(self, ResidualClaim::TrueEuclidean(_))
    }

    /// A stable one-line description for receipts and reports.
    #[must_use]
    pub fn provenance(&self) -> &'static str {
        match self {
            ResidualClaim::TrueEuclidean(_) => "recomputed Euclidean relative residual",
            ResidualClaim::RecursiveEstimate(_) => {
                "recursively propagated estimate of the Euclidean relative residual"
            }
            ResidualClaim::MNormEstimate(_) => {
                "recurrence estimate of the relative residual in the preconditioner's M-norm"
            }
        }
    }
}

/// Solve outcome with the full residual history (error transparency).
///
/// `rel_residual` is a bare number: its NORM and PROVENANCE depend on
/// the producing solver (see [`ResidualClaim`]), and this type cannot
/// carry them. A driver holding only a `SolveReport` therefore may not
/// read `converged` as `‖b − Ax‖₂/‖b‖₂ < tol`; it may read it only as
/// "the residual this solver reports met `tol`". Use the producing
/// state's `residual_claim()` when the distinction matters.
#[derive(Debug, Clone)]
pub struct SolveReport {
    /// Iterations performed (cumulative across resumes).
    pub iters: usize,
    /// Final relative residual as the producing solver measures it:
    /// recomputed `‖b − Ax‖₂/‖b‖₂` for GMRES/FGMRES, the recursive
    /// estimate of it for CG and MINRES, and the M-NORM estimate
    /// `‖r‖_M/‖b‖_M` for P-MINRES. See [`ResidualClaim`].
    pub rel_residual: f64,
    /// `rel_residual < tol` — the reported residual claim met the
    /// tolerance. NOT, on its own, a Euclidean correctness statement.
    pub converged: bool,
    /// `rel_residual` after each iteration, in the same measure.
    pub history: Vec<f64>,
    /// Present iff not converged.
    pub diagnosis: Option<StallDiagnosis>,
}

fn diagnose(history: &[f64], tol: f64) -> Option<StallDiagnosis> {
    let last = *history.last()?;
    if !last.is_finite() {
        return Some(StallDiagnosis::Breakdown);
    }
    if last < tol {
        return None;
    }
    // Plateau = no MATERIAL progress (< 5% relative) over the last
    // 50-iteration window; anything shorter-lived reads as budget.
    let window = 50.min(history.len());
    if history.len() >= 20 {
        let prev = history[history.len() - window];
        // The residual must be roughly FLAT — within ±5% of the window start.
        // The lower bound alone (`last > prev*0.95`) also matches a DIVERGING
        // tail (last ≫ prev), which the CONTRACT says must never read as
        // Plateau (singular-system CG diverges to a large finite value). The
        // upper bound routes divergence to BudgetExhausted instead.
        if prev.is_finite() && last > prev * 0.95 && last < prev * 1.05 {
            return Some(StallDiagnosis::Plateau);
        }
    }
    Some(StallDiagnosis::BudgetExhausted)
}

fn report(iters: usize, history: &[f64], bnorm_rel: f64, tol: f64) -> SolveReport {
    SolveReport {
        iters,
        rel_residual: bnorm_rel,
        converged: bnorm_rel < tol,
        history: history.to_vec(),
        diagnosis: if bnorm_rel < tol {
            None
        } else {
            diagnose(history, tol)
        },
    }
}

// ------------------------------------------------------------------------ CG

/// Resumable preconditioned conjugate gradients (SPD operators).
#[derive(Debug, Clone)]
pub struct CgState {
    /// Current iterate.
    pub x: Vec<f64>,
    r: Vec<f64>,
    p: Vec<f64>,
    rz: f64,
    bnorm: f64,
    /// Iterations performed so far (across resumes).
    pub iters: usize,
    /// Relative residual after each iteration.
    pub history: Vec<f64>,
}

impl CgState {
    /// Start a solve of A·x = b from x₀ = 0.
    #[must_use]
    pub fn new<A: LinearOp, P: Precond>(a: &A, m: &P, b: &[f64]) -> CgState {
        let n = a.n();
        assert_eq!(b.len(), n, "rhs length mismatch");
        let r = b.to_vec();
        let mut z = vec![0.0f64; n];
        m.apply(&r, &mut z);
        let rz = dot(&r, &z);
        CgState {
            x: vec![0.0f64; n],
            r,
            p: z,
            rz,
            bnorm: norm2(b).max(f64::MIN_POSITIVE),
            iters: 0,
            history: Vec::new(),
        }
    }

    /// Current relative residual ESTIMATE: `r` is the recursively
    /// updated residual, never recomputed as `b − Ax`, so this drifts
    /// from the true Euclidean residual as the iteration proceeds.
    #[must_use]
    pub fn rel_residual(&self) -> f64 {
        norm2(&self.r) / self.bnorm
    }

    /// The typed residual claim: CG reports a recursive ESTIMATE of the
    /// Euclidean relative residual.
    #[must_use]
    pub fn residual_claim(&self) -> ResidualClaim {
        ResidualClaim::RecursiveEstimate(self.rel_residual())
    }

    /// Run until `tol` or `max_iters` ADDITIONAL iterations; call
    /// again to continue bitwise-identically to a straight run.
    pub fn run<A: LinearOp, P: Precond>(
        &mut self,
        a: &A,
        m: &P,
        tol: f64,
        max_iters: usize,
    ) -> SolveReport {
        let n = a.n();
        let mut ap = vec![0.0f64; n];
        let mut z = vec![0.0f64; n];
        for _ in 0..max_iters {
            if self.rel_residual() < tol {
                break;
            }
            a.apply(&self.p, &mut ap);
            let pap = dot(&self.p, &ap);
            let alpha = self.rz / pap;
            for (i, api) in ap.iter().enumerate().take(n) {
                self.x[i] = alpha.mul_add(self.p[i], self.x[i]);
                self.r[i] = alpha.mul_add(-api, self.r[i]);
            }
            m.apply(&self.r, &mut z);
            let rz_new = dot(&self.r, &z);
            let beta = rz_new / self.rz;
            self.rz = rz_new;
            for (pi, zi) in self.p.iter_mut().zip(&z) {
                *pi = beta.mul_add(*pi, *zi);
            }
            self.iters += 1;
            self.history.push(self.rel_residual());
        }
        report(self.iters, &self.history, self.rel_residual(), tol)
    }
}

// -------------------------------------------------------------------- MINRES

/// Resumable MINRES (symmetric, possibly indefinite operators;
/// Paige–Saunders Lanczos + Givens). Unpreconditioned v1 —
/// preconditioned MINRES needs an SPD split preconditioner (recorded
/// no-claim).
#[derive(Debug, Clone)]
pub struct MinresState {
    /// Current iterate.
    pub x: Vec<f64>,
    v_prev: Vec<f64>,
    v: Vec<f64>,
    w_prev2: Vec<f64>,
    w_prev1: Vec<f64>,
    beta: f64,
    c_km1: f64,
    c_km2: f64,
    s_km1: f64,
    s_km2: f64,
    eta: f64,
    bnorm: f64,
    /// Iterations performed so far.
    pub iters: usize,
    /// Relative residual estimate after each iteration.
    pub history: Vec<f64>,
}

impl MinresState {
    /// Start a solve of A·x = b from x₀ = 0.
    #[must_use]
    pub fn new<A: LinearOp>(a: &A, b: &[f64]) -> MinresState {
        let n = a.n();
        assert_eq!(b.len(), n, "rhs length mismatch");
        let bnorm = norm2(b).max(f64::MIN_POSITIVE);
        let beta = norm2(b);
        let v: Vec<f64> = b.iter().map(|&bi| bi / beta).collect();
        MinresState {
            x: vec![0.0f64; n],
            v_prev: vec![0.0f64; n],
            v,
            w_prev2: vec![0.0f64; n],
            w_prev1: vec![0.0f64; n],
            beta,
            c_km1: 1.0,
            c_km2: 1.0,
            s_km1: 0.0,
            s_km2: 0.0,
            eta: beta,
            bnorm,
            iters: 0,
            history: Vec::new(),
        }
    }

    /// Current relative residual ESTIMATE (|η|/‖b‖ — exact in exact
    /// arithmetic; the battery cross-checks against the true
    /// residual).
    #[must_use]
    pub fn rel_residual(&self) -> f64 {
        self.eta.abs() / self.bnorm
    }

    /// The typed residual claim: MINRES reports a Givens-recurrence
    /// ESTIMATE of the Euclidean relative residual.
    #[must_use]
    pub fn residual_claim(&self) -> ResidualClaim {
        ResidualClaim::RecursiveEstimate(self.rel_residual())
    }

    /// Run until `tol` or `max_iters` additional iterations.
    pub fn run<A: LinearOp>(&mut self, a: &A, tol: f64, max_iters: usize) -> SolveReport {
        let n = a.n();
        let mut p = vec![0.0f64; n];
        for _ in 0..max_iters {
            if self.rel_residual() < tol {
                break;
            }
            // Lanczos step.
            a.apply(&self.v, &mut p);
            let alpha = dot(&self.v, &p);
            for (i, pi) in p.iter_mut().enumerate().take(n) {
                *pi = alpha.mul_add(-self.v[i], self.beta.mul_add(-self.v_prev[i], *pi));
            }
            let beta_next = norm2(&p);
            // Givens: eliminate the subdiagonal of the tridiagonal.
            let delta = self
                .c_km1
                .mul_add(alpha, -(self.c_km2 * self.s_km1 * self.beta));
            let rho1 = fs_math::det::sqrt(delta.mul_add(delta, beta_next * beta_next));
            let rho2 = self
                .s_km1
                .mul_add(alpha, self.c_km2 * self.c_km1 * self.beta);
            let rho3 = self.s_km2 * self.beta;
            let c_k = delta / rho1;
            let s_k = beta_next / rho1;
            // Direction update and iterate step.
            for i in 0..n {
                let w_k = (self.v[i] - rho3 * self.w_prev2[i] - rho2 * self.w_prev1[i]) / rho1;
                self.x[i] = (c_k * self.eta).mul_add(w_k, self.x[i]);
                self.w_prev2[i] = self.w_prev1[i];
                self.w_prev1[i] = w_k;
            }
            self.eta *= -s_k;
            // Roll the Lanczos pair and Givens memory.
            for (i, pi) in p.iter().enumerate().take(n) {
                let v_next = pi / beta_next;
                self.v_prev[i] = self.v[i];
                self.v[i] = v_next;
            }
            self.beta = beta_next;
            self.c_km2 = self.c_km1;
            self.c_km1 = c_k;
            self.s_km2 = self.s_km1;
            self.s_km1 = s_k;
            self.iters += 1;
            self.history.push(self.rel_residual());
        }
        report(self.iters, &self.history, self.rel_residual(), tol)
    }
}

// -------------------------------------------------------------------- PMINRES

/// Resumable PRECONDITIONED MINRES for symmetric (possibly indefinite)
/// operators with an SPD preconditioner M ≈ |A|⁻¹: Lanczos in the
/// M-inner product (Choi/Paige–Saunders shape). The saddle-point
/// driver (Stokes block preconditioners, bead avuw) is the consumer.
/// `rel_residual` is the residual estimate in the M-norm — exact in
/// exact arithmetic; the battery cross-checks the TRUE residual.
#[derive(Debug, Clone)]
pub struct PminresState {
    /// Current iterate.
    pub x: Vec<f64>,
    v_prev: Vec<f64>,
    v: Vec<f64>,
    z: Vec<f64>,
    w_prev2: Vec<f64>,
    w_prev1: Vec<f64>,
    beta: f64,
    c_km1: f64,
    c_km2: f64,
    s_km1: f64,
    s_km2: f64,
    eta: f64,
    bnorm: f64,
    /// Iterations performed so far.
    pub iters: usize,
    /// Relative (M-norm) residual estimate after each iteration.
    pub history: Vec<f64>,
}

impl PminresState {
    /// Start a solve of A·x = b from x₀ = 0 with SPD preconditioner `m`.
    ///
    /// # Panics
    /// If `m` is not positive on b (⟨b, M·b⟩ ≤ 0 — an SPD violation).
    #[must_use]
    pub fn new<A: LinearOp, P: Precond>(a: &A, m: &P, b: &[f64]) -> PminresState {
        let n = a.n();
        assert_eq!(b.len(), n, "rhs length mismatch");
        let mut z = vec![0.0f64; n];
        m.apply(b, &mut z);
        let bz = dot(b, &z);
        // ⟨b, Mb⟩ == 0 on a NONZERO b means M is singular on the very
        // vector the Lanczos basis starts from: β collapses to the
        // smallest denormal, v = b/β overflows, and every subsequent
        // claim is noise. Admitting it (the old `bz >= 0.0`) is
        // fail-open exactly where positivity is lost.
        assert!(
            bz > 0.0 || b.iter().all(|value| *value == 0.0),
            "preconditioner must be SPD: ⟨b, Mb⟩ = {bz} is not positive on a nonzero rhs"
        );
        let beta = fs_math::det::sqrt(bz).max(f64::MIN_POSITIVE);
        let v: Vec<f64> = b.iter().map(|&bi| bi / beta).collect();
        for zi in &mut z {
            *zi /= beta;
        }
        PminresState {
            x: vec![0.0f64; n],
            v_prev: vec![0.0f64; n],
            v,
            z,
            w_prev2: vec![0.0f64; n],
            w_prev1: vec![0.0f64; n],
            beta,
            c_km1: 1.0,
            c_km2: 1.0,
            s_km1: 0.0,
            s_km2: 0.0,
            eta: beta,
            bnorm: beta,
            iters: 0,
            history: Vec::new(),
        }
    }

    /// Current relative residual estimate IN THE M-NORM
    /// (`‖r‖_M/‖b‖_M`), not the Euclidean norm — for an ill-conditioned
    /// preconditioner the two are unrelated in magnitude.
    #[must_use]
    pub fn rel_residual(&self) -> f64 {
        self.eta.abs() / self.bnorm
    }

    /// The typed residual claim: P-MINRES reports an M-NORM estimate.
    /// A driver that needs a Euclidean statement must recompute
    /// `b − Ax` itself.
    #[must_use]
    pub fn residual_claim(&self) -> ResidualClaim {
        ResidualClaim::MNormEstimate(self.rel_residual())
    }

    /// Run until `tol` or `max_iters` additional iterations.
    ///
    /// # Panics
    /// If the preconditioner loses positivity mid-flight — i.e. the next
    /// Lanczos direction `p` is nonzero but `⟨p, Mp⟩` does not clear the
    /// rounding floor of that inner product. That is a REFUSAL, not a
    /// convergence: see the comment at the guard.
    pub fn run<A: LinearOp, P: Precond>(
        &mut self,
        a: &A,
        m: &P,
        tol: f64,
        max_iters: usize,
    ) -> SolveReport {
        let n = a.n();
        let mut p = vec![0.0f64; n];
        let mut z_next = vec![0.0f64; n];
        for _ in 0..max_iters {
            if self.rel_residual() < tol {
                break;
            }
            // Lanczos step in the M-inner product: operate on z = M·v.
            a.apply(&self.z, &mut p);
            let alpha = dot(&self.z, &p);
            for (i, pi) in p.iter_mut().enumerate().take(n) {
                *pi = alpha.mul_add(-self.v[i], self.beta.mul_add(-self.v_prev[i], *pi));
            }
            m.apply(&p, &mut z_next);
            let vz = dot(&p, &z_next);
            // ⟨p, Mp⟩ is the SQUARED M-norm of the next Lanczos
            // direction. Two very different states share `vz ≈ 0` and
            // must not be conflated:
            //
            //   * p == 0 — Lanczos HAPPY breakdown. The Krylov space is
            //     invariant, the final Givens step with β = 0 makes the
            //     iterate the exact minimizer, and η → 0 is a true
            //     residual claim.
            //   * p != 0 with ⟨p, Mp⟩ at or below the rounding floor of
            //     that inner product — M is NOT positive definite on p.
            //
            // The old guard (`vz >= -1e-30`) admitted the second case:
            // β_next collapsed to the smallest denormal, s_k ≈ 0, η ≈ 0,
            // and the next tolerance check reported `converged` with
            // `rel_residual ≈ 1e-308` for an ARBITRARY iterate. A guard
            // whose message is "preconditioner lost positivity" must not
            // pass exactly where positivity is lost.
            let p_norm = norm2(&p);
            let happy_breakdown = p_norm == 0.0;
            let beta_next = if happy_breakdown {
                0.0
            } else {
                let rounding_floor = 4.0 * f64::EPSILON * p_norm * norm2(&z_next);
                assert!(
                    vz > rounding_floor,
                    "preconditioner lost positivity: ⟨p, Mp⟩ = {vz} does not clear the rounding \
                     floor {rounding_floor} at ‖p‖ = {p_norm}; M is not positive definite on the \
                     Lanczos direction, so the M-norm residual estimate carries no claim"
                );
                fs_math::det::sqrt(vz).max(f64::MIN_POSITIVE)
            };
            // Givens: identical bookkeeping to plain MINRES.
            let delta = self
                .c_km1
                .mul_add(alpha, -(self.c_km2 * self.s_km1 * self.beta));
            let rho1 = fs_math::det::sqrt(delta.mul_add(delta, beta_next * beta_next));
            if rho1 == 0.0 {
                // Degenerate step (α and β both zero): there is no
                // rotation to apply and no residual claim to make. Stop
                // without touching the iterate — an unresolved solve, not
                // a converged one.
                break;
            }
            let rho2 = self
                .s_km1
                .mul_add(alpha, self.c_km2 * self.c_km1 * self.beta);
            let rho3 = self.s_km2 * self.beta;
            let c_k = delta / rho1;
            let s_k = beta_next / rho1;
            // Direction update uses the PRECONDITIONED vector z.
            for i in 0..n {
                let w_k = (self.z[i] - rho3 * self.w_prev2[i] - rho2 * self.w_prev1[i]) / rho1;
                self.x[i] = (c_k * self.eta).mul_add(w_k, self.x[i]);
                self.w_prev2[i] = self.w_prev1[i];
                self.w_prev1[i] = w_k;
            }
            self.eta *= -s_k;
            self.beta = beta_next;
            self.c_km2 = self.c_km1;
            self.c_km1 = c_k;
            self.s_km2 = self.s_km1;
            self.s_km1 = s_k;
            self.iters += 1;
            self.history.push(self.rel_residual());
            if happy_breakdown {
                // β_next = 0: the Lanczos pair cannot be rolled (there is
                // no next direction), and the step just applied is the
                // exact one. Stopping here keeps NaNs out of the state.
                break;
            }
            // Roll the Lanczos pair.
            for i in 0..n {
                self.v_prev[i] = self.v[i];
                self.v[i] = p[i] / beta_next;
                self.z[i] = z_next[i] / beta_next;
            }
        }
        report(self.iters, &self.history, self.rel_residual(), tol)
    }
}

// --------------------------------------------------------------------- GMRES

/// Resumable restarted GMRES(m) for general operators. Resume
/// granularity is the RESTART CYCLE: state between cycles is just the
/// iterate (plus counters), so a split run at cycle boundaries is
/// bitwise-equal to a straight run.
#[derive(Debug, Clone)]
pub struct GmresState {
    /// Current iterate.
    pub x: Vec<f64>,
    /// Restart length m.
    pub restart: usize,
    bnorm: f64,
    rel: f64,
    /// Inner iterations performed so far.
    pub iters: usize,
    /// Relative residual after each completed cycle.
    pub history: Vec<f64>,
}

impl GmresState {
    /// Start a solve of A·x = b from x₀ = 0 with restart length m.
    #[must_use]
    pub fn new(b: &[f64], restart: usize) -> GmresState {
        assert!(restart >= 1, "restart length must be >= 1");
        let bnorm = norm2(b).max(f64::MIN_POSITIVE);
        GmresState {
            x: vec![0.0f64; b.len()],
            restart,
            bnorm,
            rel: 1.0,
            iters: 0,
            history: Vec::new(),
        }
    }

    /// Current relative residual (from the last completed cycle) —
    /// recomputed as `‖b − Ax‖₂/‖b‖₂`, not a recurrence estimate.
    #[must_use]
    pub fn rel_residual(&self) -> f64 {
        self.rel
    }

    /// The typed residual claim: GMRES recomputes the TRUE Euclidean
    /// relative residual at every cycle end.
    #[must_use]
    pub fn residual_claim(&self) -> ResidualClaim {
        ResidualClaim::TrueEuclidean(self.rel)
    }

    /// Run up to `max_cycles` ADDITIONAL restart cycles (transposed
    /// solves pass `transposed = true` and use `apply_transpose` —
    /// same machinery, same preconditioner slot).
    pub fn run<A: LinearOp>(
        &mut self,
        a: &A,
        b: &[f64],
        tol: f64,
        max_cycles: usize,
        transposed: bool,
    ) -> SolveReport {
        let n = a.n();
        let m = self.restart;
        let apply = |x: &[f64], y: &mut [f64]| {
            if transposed {
                a.apply_transpose(x, y);
            } else {
                a.apply(x, y);
            }
        };
        let mut scratch = vec![0.0f64; n];
        for _ in 0..max_cycles {
            // Residual.
            apply(&self.x, &mut scratch);
            let mut r: Vec<f64> = b.iter().zip(&scratch).map(|(bi, ai)| bi - ai).collect();
            let beta = norm2(&r);
            self.rel = beta / self.bnorm;
            if self.rel < tol {
                break;
            }
            for ri in &mut r {
                *ri /= beta;
            }
            // Arnoldi with modified Gram–Schmidt; Givens-triangularized
            // Hessenberg solved at cycle end.
            let mut basis: Vec<Vec<f64>> = vec![r];
            let mut h = vec![0.0f64; (m + 1) * m];
            let mut cs = vec![0.0f64; m];
            let mut sn = vec![0.0f64; m];
            let mut g = vec![0.0f64; m + 1];
            g[0] = beta;
            let mut cols = 0usize;
            for j in 0..m {
                let mut w = vec![0.0f64; n];
                apply(&basis[j], &mut w);
                for (i, vi) in basis.iter().enumerate() {
                    let hij = dot(vi, &w);
                    h[i * m + j] = hij;
                    for (wk, vk) in w.iter_mut().zip(vi) {
                        *wk = hij.mul_add(-vk, *wk);
                    }
                }
                let hj1 = norm2(&w);
                h[(j + 1) * m + j] = hj1;
                // Apply accumulated Givens rotations to column j.
                for i in 0..j {
                    let t = cs[i].mul_add(h[i * m + j], sn[i] * h[(i + 1) * m + j]);
                    h[(i + 1) * m + j] = (-sn[i]).mul_add(h[i * m + j], cs[i] * h[(i + 1) * m + j]);
                    h[i * m + j] = t;
                }
                // New rotation killing h[j+1][j].
                let denom = fs_math::det::sqrt(h[j * m + j].mul_add(h[j * m + j], hj1 * hj1));
                cs[j] = h[j * m + j] / denom;
                sn[j] = hj1 / denom;
                h[j * m + j] = denom;
                h[(j + 1) * m + j] = 0.0;
                g[j + 1] = -sn[j] * g[j];
                g[j] *= cs[j];
                cols = j + 1;
                self.iters += 1;
                if hj1 == 0.0 || (g[j + 1].abs() / self.bnorm) < tol {
                    break;
                }
                for wk in &mut w {
                    *wk /= hj1;
                }
                basis.push(w);
            }
            // Back-substitute y and update x.
            let mut y = vec![0.0f64; cols];
            for i in (0..cols).rev() {
                let mut acc = g[i];
                for j in (i + 1)..cols {
                    acc = h[i * m + j].mul_add(-y[j], acc);
                }
                y[i] = acc / h[i * m + i];
            }
            for (yj, bj) in y.iter().zip(&basis) {
                for (xi, bji) in self.x.iter_mut().zip(bj) {
                    *xi = yj.mul_add(*bji, *xi);
                }
            }
            // True residual for the cycle-end history entry.
            apply(&self.x, &mut scratch);
            let rtrue: f64 = {
                let diff: Vec<f64> = b.iter().zip(&scratch).map(|(bi, ai)| bi - ai).collect();
                norm2(&diff)
            };
            self.rel = rtrue / self.bnorm;
            self.history.push(self.rel);
            if self.rel < tol {
                break;
            }
        }
        report(self.iters, &self.history, self.rel, tol)
    }
}

#[cfg(test)]
mod tests {
    use super::{StallDiagnosis, diagnose};

    #[test]
    fn diverging_tail_is_not_diagnosed_as_plateau() {
        // Regression (CONTRACT: "singular-system CG diverges — must never read
        // Plateau"). A monotonically GROWING residual satisfies the old lower
        // bound `last > prev*0.95` and was mislabeled Plateau. With the flat
        // ±5% band it must fall through to BudgetExhausted instead.
        let diverging: Vec<f64> = (0..60).map(|i| fs_math::det::powi(1.05, i)).collect();
        assert_eq!(
            diagnose(&diverging, 1e-8),
            Some(StallDiagnosis::BudgetExhausted),
            "a diverging residual must not read as Plateau"
        );
    }

    #[test]
    fn genuinely_flat_residual_is_plateau() {
        // A residual pinned near a fixed value across the window is the real
        // plateau (e.g. restarted-GMRES stagnation) and must still be labeled.
        let flat = vec![2.0_f64; 60];
        assert_eq!(diagnose(&flat, 1e-8), Some(StallDiagnosis::Plateau));
        // Slow-but-real progress (< 5% over the window) is a plateau too.
        let creeping: Vec<f64> = (0..60).map(|i| 2.0 - 0.0005 * f64::from(i)).collect();
        assert_eq!(diagnose(&creeping, 1e-8), Some(StallDiagnosis::Plateau));
    }

    #[test]
    fn converged_and_broken_are_classified() {
        assert_eq!(diagnose(&[1.0, 1e-12], 1e-8), None); // below tol
        assert_eq!(
            diagnose(&[1.0, f64::NAN], 1e-8),
            Some(StallDiagnosis::Breakdown)
        );
    }
}
