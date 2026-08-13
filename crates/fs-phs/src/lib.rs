//! # fs-phs — port-Hamiltonian systems
//!
//! Structural passivity: `dx/dt = (J - R) grad H(x) + G u`,
//! `y = G^T grad H(x)` with `J` skew-symmetric and `R` symmetric
//! positive semidefinite. Passivity is then a property of the
//! FORMULATION — `dH/dt = -gradH^T R gradH + u^T y <= u^T y` — not of
//! an audit after the fact, and ANY power-conjugate interconnection of
//! such systems is again passive by construction.
//!
//! Memoryless dissipative ports ([`bernoulli_volume_flow`],
//! [`quasistatic_aperture_opening`]) are first-class: a valve, a
//! vocal fold, and a beating reed are the same jet plus a 1-DOF
//! [`mass_spring_damper`] or its quasistatic reduction.
//!
//! Three pillars:
//! - **Admission**: `J`/`R` structure is verified at construction
//!   (skewness bound, symmetric eigenvalue check for `R >= 0`) and
//!   re-verified after every composition. No false passivity: refuse,
//!   never repair silently.
//! - **Gonzalez discrete gradients** ([`step`]): the midpoint discrete
//!   gradient satisfies `dg . (x1 - x0) = H(x1) - H(x0)` IDENTICALLY,
//!   so the discrete energy balance
//!   `H1 - H0 = -dt dg^T R dg + dt u^T y` holds exactly (to solver
//!   tolerance) for ANY storage — quadratic or not, damped or driven.
//!   Order 2 (the Gonzalez midpoint choice is pinned; see
//!   [`discrete_gradient`]).
//! - **Structure preservation**: power-preserving interconnection
//!   ([`interconnect`]) and Galerkin reduction ([`reduce_galerkin`])
//!   return systems in the same admitted form.
//!
//! The independent runtime check is the SUPPLY-RATE audit
//! ([`StepRecord::supply_defect`]): `H1 - H0 - dt u^T y <= 0` must
//! hold for a passive system; a symmetrized `J` or sign-flipped `R`
//! smuggled past admission (via [`PortHamiltonian::from_raw_parts`])
//! violates it observably — that is the mutation battery's job.
//!
//! DEFERRED (recorded, v1 is ODE-form): constrained/DAE Dirac
//! structures (rigid interconnection, Kirchhoff laws) wait for their
//! first consumer; fs-time strategy wiring lives with fs-time (L3,
//! above this crate).

// Dense row-major matrix kernels index by (row, col) throughout;
// iterator rewrites of these loops obscure the algebra (same call the
// other numerical crates made — fs-spectral, fs-sos, fs-tropical).
#![allow(clippy::needless_range_loop)]

use fs_la::factor::lu;
use fs_math::det;

/// Energy storage: the Hamiltonian and its exact gradient.
///
/// Implementations supply ANALYTIC gradients; the discrete-gradient
/// integrator turns them into exact discrete energy accounting.
pub trait Storage {
    /// `H(x)` — total stored energy, joules.
    fn hamiltonian(&self, x: &[f64]) -> f64;
    /// `grad H(x)` into `out` (efforts).
    fn gradient(&self, x: &[f64], out: &mut [f64]);
}

/// Quadratic storage `H = x^T Q x / 2` with `Q` symmetric PSD
/// (row-major, verified at construction).
#[derive(Debug, Clone, PartialEq)]
pub struct QuadraticStorage {
    n: usize,
    q: Vec<f64>,
}

impl QuadraticStorage {
    /// Admit a symmetric PSD `Q` (row-major `n x n`).
    ///
    /// # Errors
    /// [`PhsError::NotSymmetric`] / [`PhsError::NotPsd`].
    pub fn new(q: Vec<f64>, n: usize) -> Result<Self, PhsError> {
        require_symmetric(&q, n, "Q")?;
        require_psd(&q, n, "Q")?;
        Ok(QuadraticStorage { n, q })
    }

    /// The stiffness/compliance matrix (row-major).
    #[must_use]
    pub fn matrix(&self) -> &[f64] {
        &self.q
    }
}

impl Storage for QuadraticStorage {
    fn hamiltonian(&self, x: &[f64]) -> f64 {
        let mut acc = 0.0;
        for i in 0..self.n {
            let mut row = 0.0;
            for j in 0..self.n {
                row += self.q[i * self.n + j] * x[j];
            }
            acc += x[i] * row;
        }
        0.5 * acc
    }

    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        for i in 0..self.n {
            let mut row = 0.0;
            for j in 0..self.n {
                row += self.q[i * self.n + j] * x[j];
            }
            out[i] = row;
        }
    }
}

/// Separable non-quadratic storage: `H = sum_i h_i(x_i)` with per-state
/// scalar laws — the shape nonlinear springs and nonlinear
/// capacitors take.
pub struct SeparableStorage {
    laws: Vec<ScalarLaw>,
}

/// One scalar energy law `h(x)` with analytic derivative.
pub struct ScalarLaw {
    /// Energy `h(x)`.
    pub h: fn(f64) -> f64,
    /// Derivative `h'(x)`.
    pub dh: fn(f64) -> f64,
}

impl SeparableStorage {
    /// Build from per-state laws.
    #[must_use]
    pub fn new(laws: Vec<ScalarLaw>) -> Self {
        SeparableStorage { laws }
    }
}

impl Storage for SeparableStorage {
    fn hamiltonian(&self, x: &[f64]) -> f64 {
        self.laws.iter().zip(x).map(|(l, &xi)| (l.h)(xi)).sum()
    }

    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        for ((l, &xi), o) in self.laws.iter().zip(x).zip(out.iter_mut()) {
            *o = (l.dh)(xi);
        }
    }
}

/// Direct sum of two storages over a state split (interconnection and
/// reduction compose through this).
pub struct SumStorage {
    /// Left storage and its state dimension.
    pub a: (Box<dyn Storage>, usize),
    /// Right storage and its state dimension.
    pub b: (Box<dyn Storage>, usize),
}

impl Storage for SumStorage {
    fn hamiltonian(&self, x: &[f64]) -> f64 {
        let na = self.a.1;
        self.a.0.hamiltonian(&x[..na]) + self.b.0.hamiltonian(&x[na..])
    }

    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        let na = self.a.1;
        self.a.0.gradient(&x[..na], &mut out[..na]);
        self.b.0.gradient(&x[na..], &mut out[na..]);
    }
}

/// Typed refusal.
#[derive(Debug, Clone, PartialEq)]
pub enum PhsError {
    /// A structure matrix is not (skew-)symmetric within tolerance.
    NotSymmetric {
        /// Which matrix.
        what: &'static str,
    },
    /// `R` (or `Q`) has a negative eigenvalue beyond tolerance.
    NotPsd {
        /// Which matrix.
        what: &'static str,
    },
    /// Dimension mismatch between parts.
    Dimension {
        /// What disagreed.
        what: &'static str,
    },
    /// The implicit step's Newton iteration exhausted its budget.
    NewtonStalled {
        /// Residual norm at exhaustion.
        residual: f64,
    },
    /// Port pairing indices out of range or duplicated.
    BadPortPairing,
    /// The eigenvalue admission solve failed.
    Eigen,
}

impl core::fmt::Display for PhsError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            PhsError::NotSymmetric { what } => write!(f, "{what} violates its symmetry class"),
            PhsError::NotPsd { what } => write!(f, "{what} is not positive semidefinite"),
            PhsError::Dimension { what } => write!(f, "dimension mismatch: {what}"),
            PhsError::NewtonStalled { residual } => {
                write!(
                    f,
                    "discrete-gradient Newton stalled at residual {residual:.3e}"
                )
            }
            PhsError::BadPortPairing => write!(f, "port pairing indices invalid"),
            PhsError::Eigen => write!(f, "eigenvalue admission solve failed"),
        }
    }
}

impl std::error::Error for PhsError {}

const SYM_TOL: f64 = 1.0e-10;
const PSD_TOL: f64 = 1.0e-10;

fn matrix_scale(m: &[f64]) -> f64 {
    // Relative floor, not 1.0: an absolute floor silently loosens
    // admission for small-magnitude (SI-unit) matrices (review
    // finding: R entries ~1e-6 with 100% relative asymmetry admitted).
    m.iter()
        .fold(0.0f64, |a, &v| a.max(v.abs()))
        .max(f64::MIN_POSITIVE)
}

fn require_skew(j: &[f64], n: usize, what: &'static str) -> Result<(), PhsError> {
    let scale = matrix_scale(j);
    for i in 0..n {
        for k in 0..=i {
            if (j[i * n + k] + j[k * n + i]).abs() > SYM_TOL * scale {
                return Err(PhsError::NotSymmetric { what });
            }
        }
    }
    Ok(())
}

fn require_symmetric(r: &[f64], n: usize, what: &'static str) -> Result<(), PhsError> {
    let scale = matrix_scale(r);
    for i in 0..n {
        for k in 0..i {
            if (r[i * n + k] - r[k * n + i]).abs() > SYM_TOL * scale {
                return Err(PhsError::NotSymmetric { what });
            }
        }
    }
    Ok(())
}

fn require_psd(r: &[f64], n: usize, what: &'static str) -> Result<(), PhsError> {
    if n == 0 {
        return Ok(());
    }
    let (values, _vectors) = fs_la::eigen::jacobi_eigh(r, n);
    let scale = matrix_scale(r);
    for &lam in &values {
        if lam < -PSD_TOL * scale {
            return Err(PhsError::NotPsd { what });
        }
    }
    Ok(())
}

/// An admitted port-Hamiltonian system.
pub struct PortHamiltonian {
    n: usize,
    m: usize,
    j: Vec<f64>,
    r: Vec<f64>,
    g: Vec<f64>,
    storage: Box<dyn Storage>,
}

impl PortHamiltonian {
    /// Admit a system: `J` (`n x n`, skew), `R` (`n x n`, symmetric
    /// PSD), `G` (`n x m`), all row-major.
    ///
    /// # Errors
    /// [`PhsError`] naming the violated structure property.
    pub fn new(
        n: usize,
        m: usize,
        j: Vec<f64>,
        r: Vec<f64>,
        g: Vec<f64>,
        storage: Box<dyn Storage>,
    ) -> Result<Self, PhsError> {
        if j.len() != n * n || r.len() != n * n || g.len() != n * m {
            return Err(PhsError::Dimension {
                what: "J/R/G sizes vs (n, m)",
            });
        }
        require_skew(&j, n, "J")?;
        require_symmetric(&r, n, "R")?;
        require_psd(&r, n, "R")?;
        Ok(PortHamiltonian {
            n,
            m,
            j,
            r,
            g,
            storage,
        })
    }

    /// Bypass admission — for mutation batteries and trusted callers
    /// ONLY. The supply-rate audit exists precisely to catch a caller
    /// who violates this trust; nothing else in the crate will.
    #[must_use]
    pub fn from_raw_parts(
        n: usize,
        m: usize,
        j: Vec<f64>,
        r: Vec<f64>,
        g: Vec<f64>,
        storage: Box<dyn Storage>,
    ) -> Self {
        PortHamiltonian {
            n,
            m,
            j,
            r,
            g,
            storage,
        }
    }

    /// State dimension.
    #[must_use]
    pub fn state_dim(&self) -> usize {
        self.n
    }

    /// Port count.
    #[must_use]
    pub fn port_dim(&self) -> usize {
        self.m
    }

    /// Stored energy at `x`.
    #[must_use]
    pub fn hamiltonian(&self, x: &[f64]) -> f64 {
        self.storage.hamiltonian(x)
    }

    /// Port output `y = G^T grad H(x)`.
    #[must_use]
    pub fn output(&self, x: &[f64]) -> Vec<f64> {
        let mut e = vec![0.0; self.n];
        self.storage.gradient(x, &mut e);
        self.output_from_effort(&e)
    }

    fn output_from_effort(&self, e: &[f64]) -> Vec<f64> {
        let mut y = vec![0.0; self.m];
        for p in 0..self.m {
            let mut acc = 0.0;
            for i in 0..self.n {
                acc += self.g[i * self.m + p] * e[i];
            }
            y[p] = acc;
        }
        y
    }

    /// Structure matrices (row-major J, R, G) — read-only views.
    #[must_use]
    pub fn structure(&self) -> (&[f64], &[f64], &[f64]) {
        (&self.j, &self.r, &self.g)
    }
}

/// Gonzalez midpoint discrete gradient:
/// `dg(a, b) = gradH(m) + [(H(b) - H(a) - gradH(m).(b-a)) / |b-a|^2]
/// (b - a)` with `m = (a+b)/2`; for `b == a` it is `gradH(a)`.
///
/// Satisfies `dg . (b - a) = H(b) - H(a)` EXACTLY by construction —
/// the property the discrete energy balance rests on. Order 2 in time
/// when used in the midpoint step (the pinned formula choice).
#[must_use]
pub fn discrete_gradient(storage: &dyn Storage, a: &[f64], b: &[f64]) -> Vec<f64> {
    let n = a.len();
    let mid: Vec<f64> = a
        .iter()
        .zip(b)
        .map(|(&p, &q)| f64::midpoint(p, q))
        .collect();
    let mut dg = vec![0.0; n];
    storage.gradient(&mid, &mut dg);
    let mut dx_norm_sq = 0.0;
    let mut mid_dot = 0.0;
    for i in 0..n {
        let dx = b[i] - a[i];
        dx_norm_sq += dx * dx;
        mid_dot += dg[i] * dx;
    }
    // Cancellation guard (review finding): for |dx| below the roundoff
    // scale the correction quotient amplifies eps*|H|/|dx| without
    // bound; grad(mid) is already exact to available precision there.
    let state_scale = a
        .iter()
        .chain(b.iter())
        .fold(f64::MIN_POSITIVE, |acc, &v| acc.max(v.abs()));
    let dx_floor = 1.0e-14 * state_scale;
    if dx_norm_sq > dx_floor * dx_floor {
        let corr = (storage.hamiltonian(b) - storage.hamiltonian(a) - mid_dot) / dx_norm_sq;
        for i in 0..n {
            dg[i] += corr * (b[i] - a[i]);
        }
    }
    dg
}

/// One completed discrete-gradient step with its exact energy ledger.
#[derive(Debug, Clone)]
pub struct StepRecord {
    /// New state.
    pub x: Vec<f64>,
    /// Port output at the discrete gradient (power-conjugate to the
    /// held input over the step).
    pub y: Vec<f64>,
    /// `H(x1) - H(x0)`.
    pub delta_h: f64,
    /// `dt * dg^T R dg` — energy dissipated this step (>= 0 for an
    /// admitted `R`).
    pub dissipated: f64,
    /// `dt * u^T y` — energy supplied through the ports this step.
    pub supplied: f64,
    /// Newton iterations used.
    pub newton_iters: usize,
    /// Residual norm of the ACCEPTED implicit-solve iterate. The
    /// supply audit is blind below `~n * |dg|_inf * solver_residual`
    /// (the ledger inherits this much slack), so audit thresholds must
    /// scale with it — the crate discloses rather than hides the band.
    pub solver_residual: f64,
}

impl StepRecord {
    /// The discrete balance residual
    /// `delta_h + dissipated - supplied` — zero to solver tolerance
    /// for an ADMITTED (skew-J) system, where it restates the update
    /// equation: a SOLVER diagnostic, not a structure audit. Under a
    /// symmetrized J it equals the energy-pumping term `dt dg^T J_sym
    /// dg` instead.
    #[must_use]
    pub fn balance_residual(&self) -> f64 {
        self.delta_h + self.dissipated - self.supplied
    }

    /// The INDEPENDENT passivity audit: `delta_h - supplied` must be
    /// `<= tol` for a passive system (energy never exceeds supply).
    /// This is what catches a symmetrized `J` or sign-flipped `R`
    /// smuggled past admission.
    #[must_use]
    pub fn supply_defect(&self) -> f64 {
        self.delta_h - self.supplied
    }
}

/// Newton tolerance (relative to state scale) and budget for the
/// implicit discrete-gradient solve.
const NEWTON_TOL: f64 = 1.0e-13;
const NEWTON_MAX: usize = 50;

/// One Gonzalez discrete-gradient step of size `dt` under held input
/// `u`: solves `x1 = x0 + dt (J - R) dg(x0, x1) + dt G u` by Newton
/// with a finite-difference Jacobian, then reports the exact ledger.
///
/// # Errors
/// [`PhsError::NewtonStalled`] if the implicit solve exhausts its
/// budget; [`PhsError::Dimension`] on input-size mismatch.
pub fn step(sys: &PortHamiltonian, x0: &[f64], u: &[f64], dt: f64) -> Result<StepRecord, PhsError> {
    let n = sys.n;
    if x0.len() != n || u.len() != sys.m {
        return Err(PhsError::Dimension {
            what: "state/input length",
        });
    }
    if !dt.is_finite() || dt < 0.0 {
        return Err(PhsError::Dimension {
            what: "dt must be finite and non-negative",
        });
    }
    // Forcing term dt*G*u (held over the step).
    let mut gu = vec![0.0; n];
    for i in 0..n {
        for p in 0..sys.m {
            gu[i] += sys.g[i * sys.m + p] * u[p];
        }
    }
    let residual = |x1: &[f64], dg: &[f64]| -> Vec<f64> {
        let mut res = vec![0.0; n];
        for i in 0..n {
            let mut flow = 0.0;
            for k in 0..n {
                flow += (sys.j[i * n + k] - sys.r[i * n + k]) * dg[k];
            }
            res[i] = x1[i] - x0[i] - dt * (flow + gu[i]);
        }
        res
    };
    // State-RELATIVE stopping tolerance: an absolute floor here leaks
    // O(tol) state residual per step straight into the energy ledger
    // (executed: 3e-8 relative H drift on a lossless ladder with
    // millivolt-scale states under an absolute 1e-12). The input term
    // enters as the STATE INCREMENT dt*|G u| (review finding: bare
    // |G u| is a rate and inflates the tolerance by 1/dt).
    let scale = x0
        .iter()
        .map(|v| v.abs())
        .chain(gu.iter().map(|v| (dt * v).abs()))
        .fold(1.0e-30f64, f64::max);
    let mut x1 = x0.to_vec();
    let mut iters = 0usize;
    let mut stagnant = 0usize;
    let mut best: Option<(Vec<f64>, f64)> = None;
    let mut r_init: Option<f64> = None;
    loop {
        let dg = discrete_gradient(sys.storage.as_ref(), x0, &x1);
        let res = residual(&x1, &dg);
        let rnorm = res.iter().fold(0.0f64, |a, &v| a.max(v.abs()));
        if r_init.is_none() {
            r_init = Some(rnorm);
        }
        let improved = best.as_ref().is_none_or(|(_, b)| rnorm < *b);
        if improved {
            best = Some((x1.clone(), rnorm));
        }
        if rnorm <= NEWTON_TOL * scale {
            break;
        }
        // Stagnation with an approximate (finite-difference) Jacobian:
        // Newton residuals are not monotone far from the solution, so
        // allow a few non-improving iterates before accepting the best
        // one (review finding: first-non-improving aborts spuriously).
        if improved {
            stagnant = 0;
        } else {
            stagnant += 1;
        }
        if stagnant >= 3 || iters >= NEWTON_MAX {
            let (bx, brnorm) = best.expect("at least one iterate");
            // Acceptance scale includes the ITERATE's own magnitude:
            // the entry scale is built from x0 and the input increment
            // and can be far below where the solution actually lands
            // (executed: quasi-static reed, x0 ~ 0, solution ~ 2.5e-4,
            // FD-Jacobian noise floor ~ 4e-14 — a refusal at a
            // perfectly converged iterate). The accepted residual is
            // DISCLOSED in StepRecord::solver_residual either way.
            // ONE acceptance criterion (replacing a chain of
            // scale-chasing special cases the executed tests kept
            // falsifying): the best iterate is accepted iff its
            // residual sits >= 6 orders below BOTH candidate scales —
            // the iterate's own magnitude and the initial residual. A
            // genuinely stuck solve leaves the residual comparable to
            // one of them; a solve stagnating at the FD-Jacobian noise
            // floor clears both. The achieved residual is DISCLOSED in
            // StepRecord::solver_residual for audit scaling either
            // way.
            let sc_acc = bx.iter().fold(scale, |acc, &v| acc.max(v.abs()));
            let ref_scale = sc_acc.max(r_init.unwrap_or(sc_acc));
            if brnorm <= 1.0e-6 * ref_scale {
                x1 = bx;
                break;
            }
            return Err(PhsError::NewtonStalled { residual: brnorm });
        }
        let jac = fd_jacobian(sys, x0, &x1, scale, &residual);
        let fact = lu(&jac, n).map_err(|_| PhsError::NewtonStalled { residual: rnorm })?;
        let mut delta = res;
        fact.solve(&mut delta);
        for i in 0..n {
            x1[i] -= delta[i];
        }
        iters += 1;
    }
    // Ledger at the converged discrete gradient.
    let solver_residual = {
        let dg = discrete_gradient(sys.storage.as_ref(), x0, &x1);
        let res = residual(&x1, &dg);
        res.iter().fold(0.0f64, |a, &v| a.max(v.abs()))
    };
    Ok(ledger(sys, x0, x1, u, dt, iters, solver_residual))
}

/// Assemble the exact per-step energy ledger at the accepted iterate.
fn ledger(
    sys: &PortHamiltonian,
    x0: &[f64],
    x1: Vec<f64>,
    u: &[f64],
    dt: f64,
    newton_iters: usize,
    solver_residual: f64,
) -> StepRecord {
    let n = sys.n;
    let dg = discrete_gradient(sys.storage.as_ref(), x0, &x1);
    let mut dissipated = 0.0;
    for i in 0..n {
        let mut row = 0.0;
        for k in 0..n {
            row += sys.r[i * n + k] * dg[k];
        }
        dissipated += dg[i] * row;
    }
    dissipated *= dt;
    let y = sys.output_from_effort(&dg);
    let supplied = dt * u.iter().zip(&y).map(|(&a, &b)| a * b).sum::<f64>();
    let delta_h = sys.storage.hamiltonian(&x1) - sys.storage.hamiltonian(x0);
    StepRecord {
        x: x1,
        y,
        delta_h,
        dissipated,
        supplied,
        newton_iters,
        solver_residual,
    }
}

/// The implicit residual closure shared between the Newton loop and
/// its Jacobian.
type ResidualFn<'a> = &'a dyn Fn(&[f64], &[f64]) -> Vec<f64>;

/// Central-difference Jacobian of the implicit residual in `x1` (h^2
/// accuracy keeps the approximate-Newton linear-convergence tail
/// short).
fn fd_jacobian(
    sys: &PortHamiltonian,
    x0: &[f64],
    x1: &[f64],
    scale: f64,
    residual: ResidualFn<'_>,
) -> Vec<f64> {
    let n = x0.len();
    let mut jac = vec![0.0; n * n];
    for col in 0..n {
        let h = 1.0e-6 * (scale + x1[col].abs());
        let mut xp = x1.to_vec();
        let mut xm = x1.to_vec();
        xp[col] += h;
        xm[col] -= h;
        let dgp = discrete_gradient(sys.storage.as_ref(), x0, &xp);
        let dgm = discrete_gradient(sys.storage.as_ref(), x0, &xm);
        let rp = residual(&xp, &dgp);
        let rm = residual(&xm, &dgm);
        for row in 0..n {
            jac[row * n + col] = (rp[row] - rm[row]) / (2.0 * h);
        }
    }
    jac
}

/// Power-preserving interconnection: pair port `pa.0` of `a` with port
/// `pa.1` of `b` (for each pair) through the canonical skew coupling
/// `u_a = -y_b`, `u_b = +y_a`. Unpaired ports of both systems remain
/// external (a's first, then b's). The composite is re-admitted: its
/// `J` must verify skew and its `R` PSD, or the composition refuses.
///
/// # Errors
/// [`PhsError::BadPortPairing`] on invalid indices; admission errors
/// if the composite violates structure (cannot happen for admitted
/// inputs — the check is the no-false-passivity discipline).
pub fn interconnect(
    a: PortHamiltonian,
    b: PortHamiltonian,
    pairs: &[(usize, usize)],
) -> Result<PortHamiltonian, PhsError> {
    let (na, nb) = (a.n, b.n);
    let n = na + nb;
    // Validate pairing.
    let mut used_a = vec![false; a.m];
    let mut used_b = vec![false; b.m];
    for &(pa, pb) in pairs {
        if pa >= a.m || pb >= b.m || used_a[pa] || used_b[pb] {
            return Err(PhsError::BadPortPairing);
        }
        used_a[pa] = true;
        used_b[pb] = true;
    }
    // Composite J: block-diag(J_a, J_b) plus the skew coupling
    // -Ga_p Gb_p^T / +Gb_p Ga_p^T over the paired columns.
    let mut j = vec![0.0; n * n];
    let mut r = vec![0.0; n * n];
    for i in 0..na {
        for k in 0..na {
            j[i * n + k] = a.j[i * na + k];
            r[i * n + k] = a.r[i * na + k];
        }
    }
    for i in 0..nb {
        for k in 0..nb {
            j[(na + i) * n + na + k] = b.j[i * nb + k];
            r[(na + i) * n + na + k] = b.r[i * nb + k];
        }
    }
    for &(pa, pb) in pairs {
        for i in 0..na {
            for k in 0..nb {
                let c = a.g[i * a.m + pa] * b.g[k * b.m + pb];
                j[i * n + na + k] -= c;
                j[(na + k) * n + i] += c;
            }
        }
    }
    // External ports: unpaired of a then unpaired of b.
    let ext_a: Vec<usize> = (0..a.m).filter(|&p| !used_a[p]).collect();
    let ext_b: Vec<usize> = (0..b.m).filter(|&p| !used_b[p]).collect();
    let m = ext_a.len() + ext_b.len();
    let mut g = vec![0.0; n * m];
    for (col, &p) in ext_a.iter().enumerate() {
        for i in 0..na {
            g[i * m + col] = a.g[i * a.m + p];
        }
    }
    for (col, &p) in ext_b.iter().enumerate() {
        for i in 0..nb {
            g[(na + i) * m + ext_a.len() + col] = b.g[i * b.m + p];
        }
    }
    let storage = Box::new(SumStorage {
        a: (a.storage, na),
        b: (b.storage, nb),
    });
    PortHamiltonian::new(n, m, j, r, g, storage)
}

/// Reduced storage: `H_r(xr) = H(V xr)`, gradient `V^T gradH(V xr)`.
struct ReducedStorage {
    v: Vec<f64>,
    n: usize,
    k: usize,
    inner: Box<dyn Storage>,
}

impl ReducedStorage {
    fn lift(&self, xr: &[f64]) -> Vec<f64> {
        let mut x = vec![0.0; self.n];
        for l in 0..self.n {
            let mut acc = 0.0;
            for c in 0..self.k {
                acc += self.v[l * self.k + c] * xr[c];
            }
            x[l] = acc;
        }
        x
    }
}

impl Storage for ReducedStorage {
    fn hamiltonian(&self, xr: &[f64]) -> f64 {
        self.inner.hamiltonian(&self.lift(xr))
    }
    fn gradient(&self, xr: &[f64], out: &mut [f64]) {
        let x = self.lift(xr);
        let mut e = vec![0.0; self.n];
        self.inner.gradient(&x, &mut e);
        for c in 0..self.k {
            let mut acc = 0.0;
            for l in 0..self.n {
                acc += self.v[l * self.k + c] * e[l];
            }
            out[c] = acc;
        }
    }
}

/// Structure-preserving Galerkin reduction: project onto the columns
/// of `v` (`n x k`, row-major; caller supplies an orthonormal or at
/// least full-rank basis). `J_r = V^T J V` is skew and `R_r = V^T R V`
/// is PSD AUTOMATICALLY — that is the point of Galerkin on a pHS — and
/// the reduced system is re-admitted anyway (no false passivity).
///
/// The reduced storage evaluates `H(V x_r)`: the reduced energy is the
/// full energy of the reconstructed state, so `H_r(x_r(0)) <= H(x(0))`
/// deficit at t = 0 is exactly the energy outside the basis — the
/// honest, certified part of the reduction error story. A-priori
/// trajectory error bounds are a recorded no-claim.
///
/// # Errors
/// [`PhsError::Dimension`]; re-admission failures.
pub fn reduce_galerkin(
    sys: &PortHamiltonian,
    v: &[f64],
    k: usize,
) -> Result<PortHamiltonian, PhsError> {
    let n = sys.n;
    if v.len() != n * k || k > n {
        return Err(PhsError::Dimension {
            what: "basis size vs (n, k) with k <= n",
        });
    }
    let vt_m_v = |mat: &[f64]| -> Vec<f64> {
        // V^T M V for row-major n x n.
        let mut mv = vec![0.0; n * k];
        for i in 0..n {
            for c in 0..k {
                let mut acc = 0.0;
                for l in 0..n {
                    acc += mat[i * n + l] * v[l * k + c];
                }
                mv[i * k + c] = acc;
            }
        }
        let mut out = vec![0.0; k * k];
        for rrow in 0..k {
            for c in 0..k {
                let mut acc = 0.0;
                for l in 0..n {
                    acc += v[l * k + rrow] * mv[l * k + c];
                }
                out[rrow * k + c] = acc;
            }
        }
        out
    };
    let jr = vt_m_v(&sys.j);
    let rr = vt_m_v(&sys.r);
    let mut gr = vec![0.0; k * sys.m];
    for c in 0..k {
        for p in 0..sys.m {
            let mut acc = 0.0;
            for l in 0..n {
                acc += v[l * k + c] * sys.g[l * sys.m + p];
            }
            gr[c * sys.m + p] = acc;
        }
    }
    // The reduction needs the full system's storage; move semantics
    // would consume `sys`, so the reduced system re-uses it through a
    // shared quadratic copy when available. v1 keeps it simple: the
    // caller passes the system by reference and the reduced storage
    // clones the QUADRATIC case; non-quadratic reduction is deferred
    // with the DAE work (first-consumer trigger).
    let inner: Box<dyn Storage> = {
        let mut q = vec![0.0; n * n];
        let mut ok = true;
        for col in 0..n {
            let mut basis = vec![0.0; n];
            basis[col] = 1.0;
            let mut grad = vec![0.0; n];
            sys.storage.gradient(&basis, &mut grad);
            for row in 0..n {
                q[row * n + col] = grad[row];
            }
        }
        // Verify linearity (quadratic H): gradient at a probe point
        // must equal Q times the probe. The tolerance is PER
        // COMPONENT and relative to that row's own magnitudes — an
        // absolute floor here is the same hazard `matrix_scale`
        // documents (a nonlinear storage whose gradients sit below
        // the floor silently passes as quadratic and reduces to a
        // fully wrong surrogate), and even a GLOBAL relative scale
        // leaks the mixed-scale case where one large component (a
        // momentum row) masks a tiny nonlinear one (executed: a
        // 1e-9-stiffness Duffing q-row hid behind an O(1) p-row).
        // The row's absolute-sum dot product bounds its roundoff, so
        // a legitimately cancelling quadratic row still passes.
        let probe: Vec<f64> = (0..n).map(|i| 0.1f64.mul_add(i as f64, 0.3)).collect();
        let mut grad = vec![0.0; n];
        sys.storage.gradient(&probe, &mut grad);
        for i in 0..n {
            let mut acc = 0.0;
            let mut acc_abs = 0.0;
            for l in 0..n {
                let term = q[i * n + l] * probe[l];
                acc += term;
                acc_abs += term.abs();
            }
            if (acc - grad[i]).abs() > 1.0e-8 * (grad[i].abs() + acc_abs) {
                ok = false;
            }
        }
        if !ok {
            return Err(PhsError::Dimension {
                what: "reduction requires quadratic storage in v1 (non-quadratic deferred)",
            });
        }
        Box::new(QuadraticStorage::new(q, n)?)
    };
    let storage = Box::new(ReducedStorage {
        v: v.to_vec(),
        n,
        k,
        inner,
    });
    PortHamiltonian::new(k, sys.m, jr, rr, gr, storage)
}

// ---------------------------------------------------------------------
// The standard zoo
// ---------------------------------------------------------------------

/// Mass-spring-damper as a pHS: state `x = [q, p]`,
/// `H = p^2/(2m) + k q^2/2`, force input on `p`.
///
/// # Errors
/// Admission errors on non-physical parameters (`m <= 0`, `k < 0`,
/// `c < 0`).
pub fn mass_spring_damper(m: f64, k: f64, c: f64) -> Result<PortHamiltonian, PhsError> {
    if m <= 0.0 || k < 0.0 || c < 0.0 {
        return Err(PhsError::NotPsd {
            what: "msd parameters",
        });
    }
    let q = vec![k, 0.0, 0.0, 1.0 / m];
    let storage = Box::new(QuadraticStorage::new(q, 2)?);
    let j = vec![0.0, 1.0, -1.0, 0.0];
    let r = vec![0.0, 0.0, 0.0, c];
    let g = vec![0.0, 1.0];
    PortHamiltonian::new(2, 1, j, r, g, storage)
}

/// Lumped Helmholtz resonator as a 1-DOF acoustic pHS.
///
/// Acoustic mass `ρ L_eff / S` and compliance `V / (ρ c²)` with the
/// unflanged end correction `L_eff = L + 2 (8/3π) a`. State is
/// volume displacement and its momentum; the port is pressure
/// (effort) × volume velocity (flow). A bottle, a vented enclosure,
/// and a Helmholtz cavity are the same object.
///
/// # Errors
/// Admission errors on non-physical geometry or gas.
pub fn helmholtz_resonator(
    volume: f64,
    neck_radius: f64,
    neck_length: f64,
    density: f64,
    sound_speed: f64,
    resistance: f64,
) -> Result<PortHamiltonian, PhsError> {
    let (m_ac, stiffness, resistance) =
        helmholtz_parts(volume, neck_radius, neck_length, density, sound_speed, resistance)?;
    mass_spring_damper(m_ac, stiffness, resistance)
}

/// Flow-driven dual of [`helmholtz_resonator`]: `u` is injected
/// volume velocity, `y` is cavity pressure. A plate monopole dumps
/// `A v` into this port and reads `p` back.
///
/// # Errors
/// Same admission as [`helmholtz_resonator`].
pub fn helmholtz_resonator_flow(
    volume: f64,
    neck_radius: f64,
    neck_length: f64,
    density: f64,
    sound_speed: f64,
    resistance: f64,
) -> Result<PortHamiltonian, PhsError> {
    let (m_ac, stiffness, resistance) =
        helmholtz_parts(volume, neck_radius, neck_length, density, sound_speed, resistance)?;
    let q = vec![stiffness, 0.0, 0.0, 1.0 / m_ac];
    let storage = Box::new(QuadraticStorage::new(q, 2)?);
    let j = vec![0.0, 1.0, -1.0, 0.0];
    let r = vec![0.0, 0.0, 0.0, resistance];
    let g = vec![1.0, 0.0];
    PortHamiltonian::new(2, 1, j, r, g, storage)
}

fn helmholtz_parts(
    volume: f64,
    neck_radius: f64,
    neck_length: f64,
    density: f64,
    sound_speed: f64,
    resistance: f64,
) -> Result<(f64, f64, f64), PhsError> {
    if !(volume > 0.0
        && neck_radius > 0.0
        && neck_length >= 0.0
        && density > 0.0
        && sound_speed > 0.0
        && resistance >= 0.0
        && volume.is_finite()
        && neck_radius.is_finite()
        && neck_length.is_finite()
        && density.is_finite()
        && sound_speed.is_finite()
        && resistance.is_finite())
    {
        return Err(PhsError::NotPsd {
            what: "helmholtz resonator parameters",
        });
    }
    let pi = core::f64::consts::PI;
    let area = pi * neck_radius * neck_radius;
    let l_eff = neck_length + 2.0 * (8.0 / (3.0 * pi)) * neck_radius;
    let m_ac = density * l_eff / area;
    let c_ac = volume / (density * sound_speed * sound_speed);
    Ok((m_ac, 1.0 / c_ac, resistance))
}

/// Lossless LC ladder (discrete transmission line) with `cells` LC
/// cells: states alternate `[q_1, phi_1, q_2, phi_2, ...]` with
/// `H = sum q^2/(2C) + phi^2/(2L)`; the port drives the first cell.
///
/// # Errors
/// Admission errors on non-physical parameters.
pub fn lc_ladder(
    cells: usize,
    inductance: f64,
    capacitance: f64,
) -> Result<PortHamiltonian, PhsError> {
    if cells == 0 || inductance <= 0.0 || capacitance <= 0.0 {
        return Err(PhsError::Dimension {
            what: "lc ladder parameters",
        });
    }
    let n = 2 * cells;
    let mut q = vec![0.0; n * n];
    for cell in 0..cells {
        q[(2 * cell) * n + 2 * cell] = 1.0 / capacitance;
        q[(2 * cell + 1) * n + 2 * cell + 1] = 1.0 / inductance;
    }
    let mut j = vec![0.0; n * n];
    for cell in 0..cells {
        let (qi, pi) = (2 * cell, 2 * cell + 1);
        j[qi * n + pi] = 1.0;
        j[pi * n + qi] = -1.0;
        if cell + 1 < cells {
            let qn = 2 * (cell + 1);
            j[pi * n + qn] = -1.0;
            j[qn * n + pi] = 1.0;
        }
    }
    let r = vec![0.0; n * n];
    let mut g = vec![0.0; n];
    g[0] = 1.0;
    let storage = Box::new(QuadraticStorage::new(q, n)?);
    PortHamiltonian::new(n, 1, j, r, g, storage)
}

/// Modal bank from mass-normalized modes — the first-class bridge from
/// the eig/plate beads: per mode `i`, `H_i = (p_i^2 + w_i^2 q_i^2)/2`,
/// damping `R = diag(0, 2 zeta_i w_i)` per mode, and the drive port
/// weights `phi_i` (mode shape at the drive point, mass-normalized).
///
/// # Errors
/// Admission errors on negative damping or mismatched lengths.
pub fn modal_bank(
    omegas: &[f64],
    zetas: &[f64],
    drive: &[f64],
) -> Result<PortHamiltonian, PhsError> {
    let nm = omegas.len();
    if zetas.len() != nm || drive.len() != nm {
        return Err(PhsError::Dimension {
            what: "modal bank lengths",
        });
    }
    let n = 2 * nm;
    let mut q = vec![0.0; n * n];
    let mut j = vec![0.0; n * n];
    let mut r = vec![0.0; n * n];
    let mut g = vec![0.0; n];
    for i in 0..nm {
        if omegas[i] <= 0.0 || zetas[i] < 0.0 {
            return Err(PhsError::NotPsd {
                what: "modal parameters",
            });
        }
        let (qi, pi) = (2 * i, 2 * i + 1);
        q[qi * n + qi] = omegas[i] * omegas[i];
        q[pi * n + pi] = 1.0;
        j[qi * n + pi] = 1.0;
        j[pi * n + qi] = -1.0;
        r[pi * n + pi] = 2.0 * zetas[i] * omegas[i];
        g[pi] = drive[i];
    }
    let storage = Box::new(QuadraticStorage::new(q, n)?);
    PortHamiltonian::new(n, 1, j, r, g, storage)
}

/// Duffing storage: `H = p^2/(2m) + k q^2/2 + k3 q^4/4` with analytic
/// gradient (parameters baked in; separable laws need 'static fns).
struct DuffingStorage {
    m: f64,
    k: f64,
    k3: f64,
}

impl Storage for DuffingStorage {
    fn hamiltonian(&self, x: &[f64]) -> f64 {
        let (q, p) = (x[0], x[1]);
        p * p / (2.0 * self.m) + 0.5 * self.k * q * q + 0.25 * self.k3 * det::powi(q, 4)
    }
    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        let (q, p) = (x[0], x[1]);
        out[0] = self.k * q + self.k3 * det::powi(q, 3);
        out[1] = p / self.m;
    }
}

/// Nonlinear (Duffing-type) spring-mass: `H = p^2/(2m) + k q^2/2 +
/// k3 q^4/4` via a bespoke storage (parameters baked in) — the
/// non-quadratic exercise for the Gonzalez gradient.
///
/// # Errors
/// Admission errors on non-physical parameters.
pub fn duffing_oscillator(m: f64, k: f64, k3: f64, c: f64) -> Result<PortHamiltonian, PhsError> {
    if m <= 0.0 || k < 0.0 || k3 < 0.0 || c < 0.0 {
        return Err(PhsError::NotPsd {
            what: "duffing parameters",
        });
    }
    let j = vec![0.0, 1.0, -1.0, 0.0];
    let r = vec![0.0, 0.0, 0.0, c];
    let g = vec![0.0, 1.0];
    PortHamiltonian::new(2, 1, j, r, g, Box::new(DuffingStorage { m, k, k3 }))
}

// ---------------------------------------------------------------------
// Memoryless dissipative ports
// ---------------------------------------------------------------------

/// Two-sided Bernoulli jet through a slit.
///
/// `U = w h sgn(Δp) √(2|Δp|/ρ)` for `h > 0` and `ρ > 0`, else 0.
/// Dissipativity `Δp · U ≥ 0` is a property of the law: backflow
/// reverses with the pressure difference. This is the same port for a
/// relief valve, a vocal fold, a leaflet, or a beating reed — music
/// is not a special case.
#[must_use]
pub fn bernoulli_volume_flow(width_m: f64, opening_m: f64, dp_pa: f64, density: f64) -> f64 {
    if !(opening_m > 0.0)
        || !(width_m > 0.0)
        || !(density > 0.0)
        || !width_m.is_finite()
        || !opening_m.is_finite()
        || !dp_pa.is_finite()
        || !density.is_finite()
        || dp_pa.abs() < 1.0e-12
    {
        return 0.0;
    }
    let mag = width_m * opening_m * det::sqrt(2.0 * dp_pa.abs() / density);
    if dp_pa < 0.0 { -mag } else { mag }
}

/// Quasistatic linearly restoring aperture: `h = H max(0, 1 − Δp/P_c)`.
///
/// This is the infinite-stiffness / zero-mass reduction of a 1-DOF
/// [`mass_spring_damper`] whose restoring force balances a
/// `P_c`-scaled pressure on the face. Finite-mass leaflets keep the
/// pHS and feed its opening into [`bernoulli_volume_flow`].
#[must_use]
pub fn quasistatic_aperture_opening(
    rest_opening_m: f64,
    closing_pressure_pa: f64,
    dp_pa: f64,
) -> f64 {
    if !(rest_opening_m > 0.0)
        || !(closing_pressure_pa > 0.0)
        || !rest_opening_m.is_finite()
        || !closing_pressure_pa.is_finite()
        || !dp_pa.is_finite()
    {
        return 0.0;
    }
    let open = (1.0 - dp_pa / closing_pressure_pa).clamp(0.0, 1.0);
    rest_opening_m * open
}

#[cfg(test)]
mod valve_ports {
    use super::{bernoulli_volume_flow, quasistatic_aperture_opening};

    #[test]
    fn bernoulli_jet_is_odd_dissipative_and_scales_as_sqrt_dp() {
        let u = bernoulli_volume_flow(0.01, 4.0e-4, 100.0, 1.2);
        let back = bernoulli_volume_flow(0.01, 4.0e-4, -100.0, 1.2);
        assert!((u + back).abs() < 1.0e-16);
        assert!(100.0 * u >= 0.0);
        assert!((-100.0) * back >= 0.0);
        let u4 = bernoulli_volume_flow(0.01, 4.0e-4, 400.0, 1.2);
        assert!((u4 / u - 2.0).abs() < 1.0e-12);
        assert_eq!(bernoulli_volume_flow(0.01, 0.0, 100.0, 1.2), 0.0);
    }

    #[test]
    fn quasistatic_aperture_closes_at_the_named_pressure() {
        assert_eq!(quasistatic_aperture_opening(4.0e-4, 1_000.0, 1_000.0), 0.0);
        assert!((quasistatic_aperture_opening(4.0e-4, 1_000.0, 0.0) - 4.0e-4).abs() < 1.0e-16);
        assert!((quasistatic_aperture_opening(4.0e-4, 1_000.0, 400.0) - 2.4e-4).abs() < 1.0e-16);
    }
}
