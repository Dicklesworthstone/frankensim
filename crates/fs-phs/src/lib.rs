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

/// Descriptor (implicit) port-Hamiltonian system
/// `E ẋ = (J − R) grad H + G u` with possibly singular `E`.
///
/// The last `n - n_diff` coordinates are algebraic. The composite
/// `J` of a 0-junction is a Dirac structure on the extended space
/// `[x_a, x_b, λ]`.
pub struct DescriptorPortHamiltonian {
    n: usize,
    n_diff: usize,
    m: usize,
    j: Vec<f64>,
    r: Vec<f64>,
    g: Vec<f64>,
    storage: Box<dyn Storage>,
}

impl DescriptorPortHamiltonian {
    /// Differential-state count (the leading block of `E = I`).
    #[must_use]
    pub fn differential_dim(&self) -> usize {
        self.n_diff
    }

    /// Full extended-state count, including multipliers.
    #[must_use]
    pub fn state_dim(&self) -> usize {
        self.n
    }

    /// Port count.
    #[must_use]
    pub fn port_dim(&self) -> usize {
        self.m
    }

    /// Composite Dirac `J` (row-major).
    #[must_use]
    pub fn dirac_j(&self) -> &[f64] {
        &self.j
    }

    /// Stored energy (multipliers do not enter `H`).
    #[must_use]
    pub fn hamiltonian(&self, x: &[f64]) -> f64 {
        self.storage.hamiltonian(x)
    }

    /// Port output `y = G^T grad H`.
    #[must_use]
    pub fn output(&self, x: &[f64]) -> Vec<f64> {
        let mut e = vec![0.0; self.n];
        self.storage.gradient(x, &mut e);
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
}

/// Isolated 0-junction of two 1-port systems: `p_a = p_b`, `U_a + U_b = U_ext`.
///
/// This is [`common_effort_star`] of two members. When both ports are
/// admittance-causal (`u` = force, `y` = velocity) the same Dirac is
/// the mechanical 1-junction (common `v`, forces split).
///
/// # Errors
/// [`PhsError::BadPortPairing`] unless both are 1-port; admission
/// errors on the composite `J`/`R`.
pub fn common_effort_dirac(
    a: PortHamiltonian,
    b: PortHamiltonian,
) -> Result<DescriptorPortHamiltonian, PhsError> {
    common_effort_star(vec![a, b])
}

/// Dual reading of [`common_effort_dirac`]: admittance ports share
/// flow and split effort. A string and a plate at one attachment
/// are this junction.
///
/// # Errors
/// As [`common_effort_dirac`].
pub fn common_flow_dirac(
    a: PortHamiltonian,
    b: PortHamiltonian,
) -> Result<DescriptorPortHamiltonian, PhsError> {
    common_effort_dirac(a, b)
}

/// 1-junction on one named port pair; leftover ports stay external.
///
/// `y_a[port_a] = y_b[port_b]` and `u_a[port_a] + u_b[port_b] = 0`.
/// Every other column of `G` is kept, in order (remaining `a`, then
/// remaining `b`). A bow on a string–plate join, a blow on a
/// string–plate–duct, and a side load on a bus are this object.
/// Two 1-port members with no leftover are a closed 1-junction
/// (`m = 0`); [`common_flow_dirac`] is the same join with an
/// external force on the junction.
///
/// # Errors
/// [`PhsError::BadPortPairing`] on a bad index; admission errors
/// on the composite.
pub fn join_port(
    a: PortHamiltonian,
    b: PortHamiltonian,
    port_a: usize,
    port_b: usize,
) -> Result<DescriptorPortHamiltonian, PhsError> {
    if port_a >= a.m || port_b >= b.m {
        return Err(PhsError::BadPortPairing);
    }
    let (na, nb) = (a.n, b.n);
    let n_diff = na + nb;
    let n = n_diff + 1;
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
    let lam = n_diff;
    for i in 0..na {
        let ga = a.g[i * a.m + port_a];
        j[i * n + lam] = ga;
        j[lam * n + i] = -ga;
    }
    for i in 0..nb {
        let gb = b.g[i * b.m + port_b];
        j[(na + i) * n + lam] = -gb;
        j[lam * n + na + i] = gb;
    }
    let ext_a: Vec<usize> = (0..a.m).filter(|&p| p != port_a).collect();
    let ext_b: Vec<usize> = (0..b.m).filter(|&p| p != port_b).collect();
    let m = ext_a.len() + ext_b.len();
    let g = if m == 0 {
        Vec::new()
    } else {
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
        g
    };
    require_skew(&j, n, "J")?;
    require_symmetric(&r, n, "R")?;
    require_psd(&r, n, "R")?;
    let inner = Box::new(SumStorage {
        a: (a.storage, na),
        b: (b.storage, nb),
    });
    let lambda = Box::new(QuadraticStorage::new(vec![0.0], 1)?);
    let storage = Box::new(SumStorage {
        a: (inner, n_diff),
        b: (lambda, 1),
    });
    Ok(DescriptorPortHamiltonian {
        n,
        n_diff,
        m,
        j,
        r,
        g,
        storage,
    })
}

/// Kirchhoff star of `N ≥ 2` one-port systems: every output `y` is
/// equal and the inputs sum to `u_ext`.
///
/// Extended state `[x_1, …, x_N, λ_1, …, λ_{N-1}]` with
/// `U_i = λ_i` for `i < N`, `U_N = u_ext − Σ λ`, and algebraic
/// rows `y_i = y_N`. A three-pipe wye, a string–plate–cavity
/// pressure node, and a bus of reservoirs are the same object.
///
/// # Errors
/// [`PhsError::BadPortPairing`] unless every member is 1-port and
/// `N ≥ 2`; admission errors on the composite.
pub fn common_effort_star(
    systems: Vec<PortHamiltonian>,
) -> Result<DescriptorPortHamiltonian, PhsError> {
    if systems.len() < 2 || systems.iter().any(|s| s.m != 1) {
        return Err(PhsError::BadPortPairing);
    }
    let n_sys = systems.len();
    let n_lambda = n_sys - 1;
    let mut off = Vec::with_capacity(n_sys);
    let mut n_diff = 0usize;
    for s in &systems {
        off.push(n_diff);
        n_diff += s.n;
    }
    let n = n_diff + n_lambda;
    let mut j = vec![0.0; n * n];
    let mut r = vec![0.0; n * n];
    for (s, &o) in systems.iter().zip(&off) {
        for i in 0..s.n {
            for k in 0..s.n {
                j[(o + i) * n + o + k] = s.j[i * s.n + k];
                r[(o + i) * n + o + k] = s.r[i * s.n + k];
            }
        }
    }
    let last = n_sys - 1;
    let o_last = off[last];
    for ell in 0..n_lambda {
        let lam = n_diff + ell;
        let o = off[ell];
        for i in 0..systems[ell].n {
            j[(o + i) * n + lam] = systems[ell].g[i];
            j[lam * n + o + i] = -systems[ell].g[i];
        }
        for i in 0..systems[last].n {
            j[(o_last + i) * n + lam] = -systems[last].g[i];
            j[lam * n + o_last + i] = systems[last].g[i];
        }
    }
    let mut g = vec![0.0; n];
    for i in 0..systems[last].n {
        g[o_last + i] = systems[last].g[i];
    }
    require_skew(&j, n, "J")?;
    require_symmetric(&r, n, "R")?;
    require_psd(&r, n, "R")?;
    let mut acc: Option<(Box<dyn Storage>, usize)> = None;
    for s in systems {
        let piece = (s.storage, s.n);
        acc = Some(match acc {
            None => piece,
            Some(prev) => {
                let dim = prev.1 + piece.1;
                (Box::new(SumStorage { a: prev, b: piece }), dim)
            }
        });
    }
    let (inner, inner_n) = acc.expect("N >= 2");
    let lambda = Box::new(QuadraticStorage::new(
        vec![0.0; n_lambda * n_lambda],
        n_lambda,
    )?);
    let storage = Box::new(SumStorage {
        a: (inner, inner_n),
        b: (lambda, n_lambda),
    });
    Ok(DescriptorPortHamiltonian {
        n,
        n_diff,
        m: 1,
        j,
        r,
        g,
        storage,
    })
}

/// Implicit descriptor step: Gonzalez on the differential block,
/// algebraic residual `0 = ((J−R) e + G u)_alg` at the new state
/// with effort `e = [∇H(x), λ]` (multipliers are not stored energy).
///
/// # Errors
/// Dimension / `dt` refusals; [`PhsError::NewtonStalled`].
pub fn step_descriptor(
    sys: &DescriptorPortHamiltonian,
    x0: &[f64],
    u: &[f64],
    dt: f64,
) -> Result<StepRecord, PhsError> {
    let n = sys.n;
    let n_d = sys.n_diff;
    if x0.len() != n || u.len() != sys.m {
        return Err(PhsError::Dimension {
            what: "descriptor state/input length",
        });
    }
    if !dt.is_finite() || dt < 0.0 {
        return Err(PhsError::Dimension {
            what: "dt must be finite and non-negative",
        });
    }
    let mut gu = vec![0.0; n];
    for i in 0..n {
        for p in 0..sys.m {
            gu[i] += sys.g[i * sys.m + p] * u[p];
        }
    }
    // Multipliers are efforts, not stored coordinates: e = [∇H(x), λ].
    // Gonzalez still owns the differential block; λ enters J as itself
    // so Ga λ / −Gb λ actually drive the ports.
    let effort = |x1: &[f64]| -> Vec<f64> {
        let mut e = discrete_gradient(sys.storage.as_ref(), x0, x1);
        for i in n_d..n {
            e[i] = f64::midpoint(x0[i], x1[i]);
        }
        e
    };
    let residual = |x1: &[f64], dg: &[f64]| -> Vec<f64> {
        let mut res = vec![0.0; n];
        for i in 0..n {
            let mut flow = 0.0;
            for k in 0..n {
                flow += (sys.j[i * n + k] - sys.r[i * n + k]) * dg[k];
            }
            if i < n_d {
                res[i] = x1[i] - x0[i] - dt * (flow + gu[i]);
            } else {
                res[i] = flow + gu[i];
            }
        }
        res
    };
    let scale = x0
        .iter()
        .map(|v| v.abs())
        .chain(gu.iter().map(|v| (dt * v).abs()))
        .chain(u.iter().map(|v| v.abs()))
        .fold(1.0e-30f64, f64::max);
    let mut x1 = x0.to_vec();
    let mut iters = 0usize;
    let mut stagnant = 0usize;
    let mut best: Option<(Vec<f64>, f64)> = None;
    let mut r_init: Option<f64> = None;
    loop {
        let dg = effort(&x1);
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
        if improved {
            stagnant = 0;
        } else {
            stagnant += 1;
        }
        if stagnant >= 3 || iters >= NEWTON_MAX {
            let (bx, brnorm) = best.expect("at least one iterate");
            let sc_acc = bx.iter().fold(scale, |acc, &v| acc.max(v.abs()));
            let ref_scale = sc_acc.max(r_init.unwrap_or(sc_acc));
            if brnorm <= 1.0e-6 * ref_scale {
                x1 = bx;
                break;
            }
            return Err(PhsError::NewtonStalled { residual: brnorm });
        }
        let jac = {
            let mut jac = vec![0.0; n * n];
            for col in 0..n {
                let h = 1.0e-6 * (1.0 + scale + x1[col].abs());
                let mut xp = x1.clone();
                let mut xm = x1.clone();
                xp[col] += h;
                xm[col] -= h;
                let dgp = effort(&xp);
                let dgm = effort(&xm);
                let rp = residual(&xp, &dgp);
                let rm = residual(&xm, &dgm);
                for row in 0..n {
                    jac[row * n + col] = (rp[row] - rm[row]) / (2.0 * h);
                }
            }
            jac
        };
        let Ok(fact) = lu(&jac, n) else {
            iters += 1;
            continue;
        };
        let mut dx = res;
        fact.solve(&mut dx);
        for i in 0..n {
            x1[i] -= dx[i];
        }
        iters += 1;
    }
    let dg = effort(&x1);
    let mut y = vec![0.0; sys.m];
    for p in 0..sys.m {
        let mut acc = 0.0;
        for i in 0..n {
            acc += sys.g[i * sys.m + p] * dg[i];
        }
        y[p] = acc;
    }
    let h0 = sys.storage.hamiltonian(x0);
    let h1 = sys.storage.hamiltonian(&x1);
    let mut dissipated = 0.0;
    for i in 0..n {
        let mut rd = 0.0;
        for k in 0..n {
            rd += sys.r[i * n + k] * dg[k];
        }
        dissipated += dg[i] * rd;
    }
    dissipated *= dt;
    let mut supplied = 0.0;
    for p in 0..sys.m {
        supplied += u[p] * y[p];
    }
    supplied *= dt;
    let res = residual(&x1, &dg);
    let solver_residual = res.iter().fold(0.0f64, |a, &v| a.max(v.abs()));
    Ok(StepRecord {
        x: x1,
        y,
        delta_h: h1 - h0,
        dissipated,
        supplied,
        newton_iters: iters,
        solver_residual,
    })
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

/// Transformer of ratio `n`: `u_a = n y_b`, `u_b = −n y_a`.
///
/// Power is identically zero (`n y_b y_a − n y_a y_b = 0`). A
/// plate area, a hydraulic ram, and a lever are this object:
/// force = `n` × pressure, volume-flow = `n` × velocity.
///
/// # Errors
/// [`PhsError::BadPortPairing`] on a bad port index; [`PhsError::NotPsd`]
/// on a non-finite ratio; admission errors on the composite.
pub fn transformer(
    a: PortHamiltonian,
    b: PortHamiltonian,
    port_a: usize,
    port_b: usize,
    ratio: f64,
) -> Result<PortHamiltonian, PhsError> {
    if port_a >= a.m || port_b >= b.m {
        return Err(PhsError::BadPortPairing);
    }
    if !ratio.is_finite() {
        return Err(PhsError::NotPsd {
            what: "transformer ratio",
        });
    }
    let (na, nb) = (a.n, b.n);
    let n = na + nb;
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
    for i in 0..na {
        for k in 0..nb {
            let c = ratio * a.g[i * a.m + port_a] * b.g[k * b.m + port_b];
            // Opposite gyrator sign: u_a = n y_b, u_b = −n y_a.
            j[i * n + na + k] += c;
            j[(na + k) * n + i] -= c;
        }
    }
    let ext_a: Vec<usize> = (0..a.m).filter(|&p| p != port_a).collect();
    let ext_b: Vec<usize> = (0..b.m).filter(|&p| p != port_b).collect();
    let m = ext_a.len() + ext_b.len();
    let mut g = vec![0.0; n * m.max(1)];
    if m > 0 {
        g.truncate(n * m);
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
    } else {
        g = Vec::new();
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
    let (m_ac, stiffness, resistance) = helmholtz_parts(
        volume,
        neck_radius,
        neck_length,
        density,
        sound_speed,
        resistance,
    )?;
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
    let (m_ac, stiffness, resistance) = helmholtz_parts(
        volume,
        neck_radius,
        neck_length,
        density,
        sound_speed,
        resistance,
    )?;
    let q = vec![stiffness, 0.0, 0.0, 1.0 / m_ac];
    let storage = Box::new(QuadraticStorage::new(q, 2)?);
    let j = vec![0.0, 1.0, -1.0, 0.0];
    let r = vec![0.0, 0.0, 0.0, resistance];
    let g = vec![1.0, 0.0];
    PortHamiltonian::new(2, 1, j, r, g, storage)
}

/// Mouth baffle for the low-`ka` radiation load.
///
/// Coefficients match `fs_duct::Termination` (Levine–Schwinger
/// unflanged, flanged 0.8216). This crate stays L2 and does not
/// depend on the duct TMM; the numbers are the same physical fit.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MouthFlange {
    /// Unflanged open pipe: `R / Z0 = (ka)²/4`, end correction `0.6133 a`.
    Unflanged,
    /// Infinite-baffle / flanged: `R / Z0 = (ka)²/2`, end correction `0.8216 a`.
    Flanged,
}

/// Compact-piston radiation impedance `(R, X)` [Pa s / m³] under the
/// workspace `e^{-iωt}` convention (mass-like `X < 0`).
///
/// Refused above `ka = 1` — the same ceiling as the duct TMM. This is
/// the resistive *and* reactive load of a compact opening, not a
/// far-field observer model.
///
/// # Errors
/// [`PhsError::NotPsd`] on non-physical inputs or `ka > 1`.
pub fn compact_radiation_impedance(
    density: f64,
    sound_speed: f64,
    radius: f64,
    omega: f64,
    flange: MouthFlange,
) -> Result<(f64, f64), PhsError> {
    if !(density > 0.0
        && sound_speed > 0.0
        && radius > 0.0
        && omega > 0.0
        && density.is_finite()
        && sound_speed.is_finite()
        && radius.is_finite()
        && omega.is_finite())
    {
        return Err(PhsError::NotPsd {
            what: "compact radiation parameters",
        });
    }
    let ka = omega / sound_speed * radius;
    if ka > 1.0 {
        return Err(PhsError::NotPsd {
            what: "compact radiation ka > 1 (low-ka fit ceiling)",
        });
    }
    let area = core::f64::consts::PI * radius * radius;
    let z0 = density * sound_speed / area;
    let (r_over, x_over) = match flange {
        MouthFlange::Unflanged => (0.25 * ka * ka, -0.6133 * ka),
        MouthFlange::Flanged => (0.50 * ka * ka, -0.8216 * ka),
    };
    Ok((r_over * z0, x_over * z0))
}

/// [`helmholtz_resonator`] whose damper is the compact-mouth radiation
/// resistance evaluated at the lossless natural frequency.
///
/// Added mass stays the neck end correction already inside
/// [`helmholtz_parts`]; this only turns `Re Z_rad(ω₀)` into the pHS
/// `R` entry. Frequency-dependent `R(ω)` is the tabulated-impedance
/// path, not this lumped zoo member.
///
/// # Errors
/// Admission errors from the resonator or the radiation fit.
pub fn helmholtz_resonator_radiating(
    volume: f64,
    neck_radius: f64,
    neck_length: f64,
    density: f64,
    sound_speed: f64,
    flange: MouthFlange,
) -> Result<PortHamiltonian, PhsError> {
    let (m_ac, stiffness, _) =
        helmholtz_parts(volume, neck_radius, neck_length, density, sound_speed, 0.0)?;
    let omega0 = det::sqrt(stiffness / m_ac);
    let (r_rad, _) =
        compact_radiation_impedance(density, sound_speed, neck_radius, omega0, flange)?;
    mass_spring_damper(m_ac, stiffness, r_rad)
}

/// Series connection of two impedance-causal 1-ports: same `u`,
/// `y = y_a + y_b`.
///
/// This is the ODE-form of stacked impedances (`Z = Z_a + Z_b`). It is
/// **not** Kirchhoff common-effort (same `p`, opposite `U`), which
/// remains a deferred DAE Dirac structure. Both systems must expose
/// exactly one port.
///
/// # Errors
/// [`PhsError::BadPortPairing`] unless both are 1-port; admission
/// errors on the composite.
pub fn series_impedance_ports(
    a: PortHamiltonian,
    b: PortHamiltonian,
) -> Result<PortHamiltonian, PhsError> {
    if a.m != 1 || b.m != 1 {
        return Err(PhsError::BadPortPairing);
    }
    let (na, nb) = (a.n, b.n);
    let n = na + nb;
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
    let mut g = vec![0.0; n];
    for i in 0..na {
        g[i] = a.g[i];
    }
    for i in 0..nb {
        g[na + i] = b.g[i];
    }
    let storage = Box::new(SumStorage {
        a: (a.storage, na),
        b: (b.storage, nb),
    });
    PortHamiltonian::new(n, 1, j, r, g, storage)
}

/// Two-port common-effort capacitor: both ports see the same pressure
/// `p = q / C`, and `q̇ = U₁ + U₂`.
///
/// This is the ODE image of a Kirchhoff effort junction with storage.
/// Identifying two *impedance-causal* ports at the same `p` without
/// this capacitor remains a DAE Dirac structure.
///
/// # Errors
/// [`PhsError::NotPsd`] if `compliance` is not positive and finite.
pub fn common_effort_capacitor(compliance: f64) -> Result<PortHamiltonian, PhsError> {
    if !(compliance > 0.0 && compliance.is_finite()) {
        return Err(PhsError::NotPsd {
            what: "common-effort compliance",
        });
    }
    let storage = Box::new(QuadraticStorage::new(vec![1.0 / compliance], 1)?);
    PortHamiltonian::new(1, 2, vec![0.0], vec![0.0], vec![1.0, 1.0], storage)
}

/// One implicit step of two 1-port systems sharing effort.
///
/// Finds `U_a` such that `p_a(U_a) = p_b(U_ext − U_a)` by Newton,
/// then commits those Gonzalez steps. This is the index-1 Kirchhoff
/// current law for impedance-causal ports (same `p`, flows sum to
/// `u_ext`). It is **not** a Dirac structure on the composite `J`
/// (that remains DAE). The capacitor [`common_effort_capacitor`] is
/// the ODE regularization of the same junction.
///
/// # Errors
/// [`PhsError::BadPortPairing`] unless both are 1-port;
/// [`PhsError::NewtonStalled`] if the algebraic split does not
/// close; step errors from either system.
pub fn kirchhoff_parallel_step(
    a: &PortHamiltonian,
    xa: &[f64],
    b: &PortHamiltonian,
    xb: &[f64],
    u_ext: f64,
    dt: f64,
) -> Result<(StepRecord, StepRecord), PhsError> {
    if a.m != 1 || b.m != 1 {
        return Err(PhsError::BadPortPairing);
    }
    if !u_ext.is_finite() {
        return Err(PhsError::Dimension {
            what: "external flow must be finite",
        });
    }
    let mut ua = 0.5 * u_ext;
    let mut rec_a = step(a, xa, &[ua], dt)?;
    let mut rec_b = step(b, xb, &[u_ext - ua], dt)?;
    for _ in 0..12 {
        let residual = rec_a.y[0] - rec_b.y[0];
        let scale = 1.0 + rec_a.y[0].abs() + rec_b.y[0].abs();
        if residual.abs() <= 1.0e-8 * scale {
            return Ok((rec_a, rec_b));
        }
        let h = 1.0e-6 * (1.0 + ua.abs());
        let pa_plus = step(a, xa, &[ua + h], dt)?.y[0];
        let pb_minus = step(b, xb, &[u_ext - ua - h], dt)?.y[0];
        let deriv = (pa_plus - pb_minus - residual) / h;
        if deriv.abs() < 1.0e-18 {
            break;
        }
        ua -= residual / deriv;
        rec_a = step(a, xa, &[ua], dt)?;
        rec_b = step(b, xb, &[u_ext - ua], dt)?;
    }
    let residual = rec_a.y[0] - rec_b.y[0];
    let scale = 1.0 + rec_a.y[0].abs() + rec_b.y[0].abs();
    if residual.abs() <= 1.0e-6 * scale {
        return Ok((rec_a, rec_b));
    }
    Err(PhsError::NewtonStalled {
        residual: residual.abs(),
    })
}

/// Regularized Coulomb traction [N]: `F = −μ N tanh(v / v_reg)`.
///
/// Dissipative by construction (`F v ≤ 0`). This is the memoryless
/// friction port — a bow, a brake, and a fault are the same law. The
/// stick reaction at exact rest is the `v_reg → 0` limit, not a
/// complementary constraint.
#[must_use]
pub fn regularized_coulomb(mu: f64, normal_n: f64, velocity: f64, v_reg: f64) -> f64 {
    if !(mu >= 0.0 && v_reg > 0.0 && mu.is_finite() && v_reg.is_finite())
        || !normal_n.is_finite()
        || !velocity.is_finite()
    {
        return 0.0;
    }
    -mu * normal_n.abs() * det::tanh(velocity / v_reg)
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

/// [`lc_ladder`] with a resistive termination on the last cell's
/// flux — the lumped image of a compact radiation (or any passive
/// load) at the far end of a discrete waveguide. `r_load = 0` is
/// the lossless ladder.
///
/// # Errors
/// Admission errors on non-physical parameters.
pub fn lc_ladder_terminated(
    cells: usize,
    inductance: f64,
    capacitance: f64,
    r_load: f64,
) -> Result<PortHamiltonian, PhsError> {
    if r_load < 0.0 || !r_load.is_finite() {
        return Err(PhsError::NotPsd {
            what: "lc ladder termination resistance",
        });
    }
    let mut sys = lc_ladder(cells, inductance, capacitance)?;
    if r_load == 0.0 {
        return Ok(sys);
    }
    let n = sys.n;
    let last_p = n - 1;
    sys.r[last_p * n + last_p] = r_load;
    PortHamiltonian::new(n, sys.m, sys.j, sys.r, sys.g, sys.storage)
}

/// Inviscid acoustic cylinder as an LC ladder.
///
/// Per cell: inertance `ρ dx / A`, compliance `A dx /(ρ c²)`.
/// The inlet port is flow-in / effort-out (`u` = volume velocity,
/// `y` = pressure). `inlets = 2` puts two identical columns on the
/// first charge so a blow and a transformer body can share the
/// mouth. An open far end loads the last flux with compact-mouth
/// `Re Z_rad` at the quarter-wave pin; a closed end is lossless.
/// A stepped bore is [`acoustic_chain`]. Wall losses are the
/// optional [`ViscothermalPin`] on that chain — not Stokes–
/// Kirchhoff bulk `α` and not ISO 9613.
///
/// # Errors
/// Non-physical geometry, `cells < 2`, `inlets` not in `{1, 2}`,
/// or a radiation-fit refusal.
pub fn acoustic_cylinder(
    length: f64,
    radius: f64,
    density: f64,
    sound_speed: f64,
    cells: usize,
    open: bool,
    inlets: usize,
) -> Result<PortHamiltonian, PhsError> {
    acoustic_waveguide(
        length,
        radius,
        density,
        sound_speed,
        cells,
        open,
        inlets,
        &[],
    )
}

/// Open side branch on an [`acoustic_chain`]: a neck inertance
/// shunted to atmosphere at a station along the line.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AcousticTap {
    /// Station along the total length in `[0, 1]`.
    pub station: f64,
    /// Chimney / neck length [m] (end correction is added).
    pub neck_length: f64,
    /// Neck radius [m].
    pub neck_radius: f64,
}

/// Wide-tube viscothermal pin for an [`acoustic_chain`].
///
/// First-order Zwikker–Kosten at one frequency: shear number
/// `r_v = a √(ω ρ / μ)`, series wall resistance
/// `R = ω L √2 / r_v`, thermal shunt
/// `G = ω C (γ−1) √2 / (r_v √Pr)` when `r_v ≥ 10`. Below that
/// shear number the pin is the same Poiseuille + isothermal-tending
/// shunt as `fs_duct::LossModel::AllRegime` (series `8 μ / a²` on
/// `L`, thermal `G` on `C`, inertance `4/3` of the inviscid value).
/// Evaluated at the quarter-wave pin. Zero viscosity is the
/// lossless mutation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ViscothermalPin {
    /// Dynamic viscosity `μ` [Pa s].
    pub dynamic_viscosity: f64,
    /// Heat-capacity ratio `γ`.
    pub gamma: f64,
    /// Prandtl number `μ c_p / κ`.
    pub prandtl: f64,
}

/// One cylindrical run in an [`acoustic_chain`].
///
/// A muffler chamber, a constriction, and one slice of a cone
/// approximated by cylinders are this object. `cells` is the LC
/// resolution of this run only.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AcousticSection {
    /// Axial length [m].
    pub length: f64,
    /// Inner radius [m].
    pub radius: f64,
    /// LC cells in this run (`>= 1`; the chain as a whole needs `>= 2`).
    pub cells: usize,
}

/// Inviscid uniform waveguide with optional open side branches.
///
/// Each open tap is a neck inertance `ρ (ℓ + 0.6 a)/A` shunted to
/// `p = 0` with compact-mouth `Re Z_rad` — a tone hole, a relief
/// port, and a side vent are this object. Closed pads add no states.
/// A stepped bore is [`acoustic_chain`].
///
/// # Errors
/// As [`acoustic_cylinder`], plus a bad tap station or neck.
pub fn acoustic_waveguide(
    length: f64,
    radius: f64,
    density: f64,
    sound_speed: f64,
    cells: usize,
    open: bool,
    inlets: usize,
    taps: &[AcousticTap],
) -> Result<PortHamiltonian, PhsError> {
    acoustic_chain(
        &[AcousticSection {
            length,
            radius,
            cells,
        }],
        density,
        sound_speed,
        open,
        inlets,
        taps,
        None,
    )
}

/// Concatenated LC waveguide, optionally wide-tube lossy.
///
/// Each section has its own `L = ρ dx / A`, `C = A dx /(ρ c²)`.
/// Adjacent runs share the same series `J` as cells inside a run:
/// the interface is the area jump, not a special junction. Open
/// taps sit at a station of the total length. An open far end
/// loads the last flux with compact-mouth `Re Z_rad` at the
/// quarter-wave pin of the total length, using the last radius.
/// A [`ViscothermalPin`] adds first-order Zwikker–Kosten `R` on
/// each inertance and thermal `G` on each compliance at that
/// same pin. Chimney taps stay inviscid (the duct TMM's no-claim).
/// Full-spectrum TMM remains the all-regime fallback.
///
/// # Errors
/// Empty chain, non-physical geometry, `inlets` not in `{1, 2}`,
/// fewer than two cells in total, a bad tap, a bad pin, or a
/// radiation-fit refusal.
pub fn acoustic_chain(
    sections: &[AcousticSection],
    density: f64,
    sound_speed: f64,
    open: bool,
    inlets: usize,
    taps: &[AcousticTap],
    viscothermal: Option<&ViscothermalPin>,
) -> Result<PortHamiltonian, PhsError> {
    let mut n_cells = 0usize;
    let mut last_radius = 0.0;
    for s in sections {
        if !(s.length > 0.0 && s.radius > 0.0) || s.cells == 0 {
            return Err(PhsError::NotPsd {
                what: "acoustic section parameters",
            });
        }
        n_cells += s.cells;
        last_radius = s.radius;
    }
    if sections.is_empty()
        || n_cells < 2
        || !(density > 0.0 && sound_speed > 0.0)
        || (inlets != 1 && inlets != 2)
    {
        return Err(PhsError::NotPsd {
            what: "acoustic chain parameters",
        });
    }
    if let Some(pin) = viscothermal {
        if !(pin.dynamic_viscosity >= 0.0
            && pin.dynamic_viscosity.is_finite()
            && pin.gamma > 1.0
            && pin.gamma.is_finite()
            && pin.prandtl > 0.0
            && pin.prandtl.is_finite())
        {
            return Err(PhsError::NotPsd {
                what: "viscothermal pin parameters",
            });
        }
    }
    let n_line = 2 * n_cells;
    let n_tap = taps.len();
    let n = n_line + n_tap;
    let mut inertances = Vec::with_capacity(n_cells);
    let mut compliances = Vec::with_capacity(n_cells);
    let mut radii = Vec::with_capacity(n_cells);
    let mut cell_x0 = Vec::with_capacity(n_cells);
    let mut x_acc = 0.0;
    for s in sections {
        let area = core::f64::consts::PI * s.radius * s.radius;
        let dx = s.length / s.cells as f64;
        let inertance = density * dx / area;
        let compliance = area * dx / (density * sound_speed * sound_speed);
        for _ in 0..s.cells {
            cell_x0.push(x_acc);
            inertances.push(inertance);
            compliances.push(compliance);
            radii.push(s.radius);
            x_acc += dx;
        }
    }
    let total_length = x_acc;
    let mut q = vec![0.0; n * n];
    let mut j = vec![0.0; n * n];
    let mut r = vec![0.0; n * n];
    for cell in 0..n_cells {
        let (qi, pi) = (2 * cell, 2 * cell + 1);
        q[qi * n + qi] = 1.0 / compliances[cell];
        q[pi * n + pi] = 1.0 / inertances[cell];
        j[qi * n + pi] = 1.0;
        j[pi * n + qi] = -1.0;
        if cell + 1 < n_cells {
            let qn = 2 * (cell + 1);
            j[pi * n + qn] = -1.0;
            j[qn * n + pi] = 1.0;
        }
    }
    let omega_pin = core::f64::consts::PI * sound_speed / (2.0 * total_length);
    if let Some(pin) = viscothermal {
        for cell in 0..n_cells {
            let (r_series, g_shunt, l_scale) = all_regime_series_and_shunt(
                inertances[cell],
                compliances[cell],
                radii[cell],
                density,
                omega_pin,
                pin,
            );
            let (qi, pi) = (2 * cell, 2 * cell + 1);
            q[pi * n + pi] = 1.0 / (inertances[cell] * l_scale);
            r[qi * n + qi] += g_shunt;
            r[pi * n + pi] += r_series;
        }
    }
    if open {
        let r_load = compact_radiation_impedance(
            density,
            sound_speed,
            last_radius,
            omega_pin,
            MouthFlange::Unflanged,
        )
        .map(|(rr, _)| rr)?;
        let last_p = n_line - 1;
        r[last_p * n + last_p] = r_load;
    }
    for (t, tap) in taps.iter().enumerate() {
        if !(tap.station >= 0.0
            && tap.station <= 1.0
            && tap.neck_length >= 0.0
            && tap.neck_radius > 0.0)
        {
            return Err(PhsError::NotPsd {
                what: "acoustic tap station and neck",
            });
        }
        let pos = tap.station * total_length;
        let mut cell = 0usize;
        for (i, &x0) in cell_x0.iter().enumerate() {
            if x0 <= pos {
                cell = i;
            }
        }
        let tap_q = 2 * cell;
        let phi = n_line + t;
        let a_h = core::f64::consts::PI * tap.neck_radius * tap.neck_radius;
        let l_eff = tap.neck_length + 0.6 * tap.neck_radius;
        let l_h = density * l_eff.max(1.0e-6) / a_h;
        q[phi * n + phi] = 1.0 / l_h;
        j[tap_q * n + phi] = -1.0;
        j[phi * n + tap_q] = 1.0;
        let r_ac = compact_radiation_impedance(
            density,
            sound_speed,
            tap.neck_radius,
            omega_pin,
            MouthFlange::Unflanged,
        )
        .map(|(rr, _)| rr)
        .unwrap_or(0.0);
        r[phi * n + phi] = r_ac;
    }
    let mut g = vec![0.0; n * inlets];
    g[0] = 1.0;
    if inlets == 2 {
        g[1] = 1.0;
    }
    PortHamiltonian::new(n, inlets, j, r, g, Box::new(QuadraticStorage::new(q, n)?))
}

/// All-regime wall law at one `ω`: wide-tube ZK for `r_v ≥ 10`,
/// Poiseuille + isothermal-tending shunt below. Returns
/// `(R_series, G_shunt, L_scale)`.
fn all_regime_series_and_shunt(
    inertance: f64,
    compliance: f64,
    radius: f64,
    density: f64,
    omega: f64,
    pin: &ViscothermalPin,
) -> (f64, f64, f64) {
    if !(pin.dynamic_viscosity > 0.0 && omega > 0.0 && radius > 0.0 && density > 0.0) {
        return (0.0, 0.0, 1.0);
    }
    let rv = radius * det::sqrt(omega * density / pin.dynamic_viscosity);
    if !(rv > 0.0 && rv.is_finite()) {
        return (0.0, 0.0, 1.0);
    }
    const WIDE_TUBE_SHEAR: f64 = 10.0;
    if rv >= WIDE_TUBE_SHEAR {
        let eps = core::f64::consts::SQRT_2 / rv;
        let r_series = omega * inertance * eps;
        let g_shunt = omega * compliance * (pin.gamma - 1.0) * eps / det::sqrt(pin.prandtl);
        return (r_series.max(0.0), g_shunt.max(0.0), 1.0);
    }
    let r_series = 8.0 * pin.dynamic_viscosity / (density * radius * radius) * inertance;
    let rt = rv * det::sqrt(pin.prandtl);
    let g_shunt = (pin.gamma - 1.0) * compliance * omega * (rt * rt / 16.0).min(0.5);
    (r_series.max(0.0), g_shunt.max(0.0), 4.0 / 3.0)
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
    modal_bank_ports(omegas, zetas, &[drive])
}

/// Mass-normalized modal bank with one drive column per port.
///
/// A plate that sees both an attachment force and a face pressure
/// is this object: two columns of `φ` at two points.
///
/// # Errors
/// Admission errors on negative damping, empty ports, or mismatched
/// lengths.
pub fn modal_bank_ports(
    omegas: &[f64],
    zetas: &[f64],
    drives: &[&[f64]],
) -> Result<PortHamiltonian, PhsError> {
    let nm = omegas.len();
    let m = drives.len();
    if zetas.len() != nm || m == 0 || drives.iter().any(|d| d.len() != nm) {
        return Err(PhsError::Dimension {
            what: "modal bank lengths",
        });
    }
    let n = 2 * nm;
    let mut q = vec![0.0; n * n];
    let mut j = vec![0.0; n * n];
    let mut r = vec![0.0; n * n];
    let mut g = vec![0.0; n * m];
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
        for (p, drive) in drives.iter().enumerate() {
            g[pi * m + p] = drive[i];
        }
    }
    let storage = Box::new(QuadraticStorage::new(q, n)?);
    PortHamiltonian::new(n, m, j, r, g, storage)
}

/// Moving-end taut waveguide: free-fixed eigenfunctions
/// `φ_k(x) = √(2/μL) cos((k−½)π x/L)`, 1-port at the free end.
///
/// `φ_k(0) ≠ 0`, so the port output is the attachment velocity.
/// Fixed-fixed sines have `φ(0)=0` and cannot Dirac-join a body.
/// A cable, a stay, and a string on a moving support are this
/// object.
///
/// # Errors
/// Non-physical geometry or damping.
pub fn moving_end_waveguide(
    n_modes: usize,
    length: f64,
    tension: f64,
    lin_density: f64,
    zetas: &[f64],
) -> Result<PortHamiltonian, PhsError> {
    if n_modes == 0
        || zetas.len() != n_modes
        || !(length > 0.0 && tension > 0.0 && lin_density > 0.0)
    {
        return Err(PhsError::NotPsd {
            what: "moving-end waveguide parameters",
        });
    }
    let c = det::sqrt(tension / lin_density);
    let phi0 = det::sqrt(2.0 / (lin_density * length));
    let pi = core::f64::consts::PI;
    let omegas: Vec<f64> = (0..n_modes)
        .map(|k| (k as f64 + 0.5) * pi * c / length)
        .collect();
    let drive = vec![phi0; n_modes];
    modal_bank(&omegas, zetas, &drive)
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
