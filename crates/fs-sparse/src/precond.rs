//! Preconditioner components (plan §6.2): Chebyshev polynomial smoothers,
//! zero-fill incomplete factorizations, and smoothed-aggregation AMG —
//! the toolkit the FLUX solver stack composes.
//!
//! Everything exposes the OPERATOR interface ([`Precond::apply`]); nothing
//! here requires assembling what can be applied. Deterministic by
//! construction: fixed iteration orders, index-order aggregation with
//! lowest-index tie-breaks (P2 applies to SETUP, not just solves), power
//! iterations from fixed start vectors, and NO platform libm anywhere in
//! solver state (contract rule since the eigensolver divergence).
//!
//! v1 scope notes: ILU/IC are sequential (level scheduling is a recorded
//! perf refinement; the bead itself labels them the bandwidth-hostile
//! fallback, not the default). The AMG coarsest level solves with
//! IC(0)-preconditioned CG in-crate (a dense direct coarse solve joins
//! solver-stack integration). Supernodal Cholesky is deferred per its own
//! "not a load-bearing wall" scope cap.

use crate::Csr;
use crate::ops::{spgemm, transpose};

/// The operator interface: z ← M⁻¹·r (approximately).
pub trait Precond {
    /// Apply the preconditioner to `r`, writing `z` (same length).
    fn apply(&self, r: &[f64], z: &mut [f64]);
}

/// Identity (the "no preconditioner" baseline).
#[derive(Debug, Clone, Copy)]
pub struct IdentityPrecond;

impl Precond for IdentityPrecond {
    fn apply(&self, r: &[f64], z: &mut [f64]) {
        z.copy_from_slice(r);
    }
}

// ---------------------------------------------------------------------------
// Spectral estimates and conservative bounds
// ---------------------------------------------------------------------------

/// Estimate λ_max of the SPD operator by power iteration from a FIXED
/// start vector, multiplied by `safety` (≥ 1).
///
/// This is a heuristic, not an upper bound: the starting vector can be
/// orthogonal, or nearly orthogonal, to a dominant eigenspace. No fixed
/// safety factor repairs that omission. Use [`spectral_radius_upper_bound`]
/// when a smoother or relaxation step needs an enclosing spectral bound.
#[must_use]
pub fn lambda_max_estimate(a: &Csr, iters: usize, safety: f64) -> f64 {
    let n = a.nrows();
    // Fixed deterministic start: varying pattern, no libm, no RNG.
    let mut v: Vec<f64> = (0..n)
        .map(|i| 1.0 + 0.3 * (((i * 2_654_435_761) % 97) as f64) / 97.0)
        .collect();
    let mut w = vec![0.0f64; n];
    let mut lam = 0.0f64;
    for _ in 0..iters {
        let nrm = v.iter().map(|x| x * x).sum::<f64>().sqrt();
        for x in &mut v {
            *x /= nrm;
        }
        a.spmv(&v, &mut w);
        lam = v.iter().zip(&w).map(|(x, y)| x * y).sum::<f64>();
        std::mem::swap(&mut v, &mut w);
    }
    lam * safety
}

/// Conservative upper bound on the spectral radius of a square matrix.
///
/// Every eigenvalue satisfies `|lambda| <= max_i sum_j |a_ij|`, including
/// for a nonsymmetric row-scaled operator such as `D^-1 A`. Each nontrivial
/// sum is rounded upward, so summation cannot lower the bound below the
/// exact row sum of the supplied binary64 coefficients. This does not prove
/// symmetry or positive definiteness. An empty/zero matrix returns zero;
/// an unrepresentable bound returns infinity.
///
/// # Panics
/// Non-square input or a non-finite matrix coefficient.
#[must_use]
pub fn spectral_radius_upper_bound(a: &Csr) -> f64 {
    assert_eq!(
        a.nrows(),
        a.ncols(),
        "spectral bound requires a square matrix"
    );
    let mut bound = 0.0_f64;
    for row in 0..a.nrows() {
        let mut sum = 0.0;
        for &value in a.row(row).1 {
            assert!(
                value.is_finite(),
                "spectral bound requires finite coefficients"
            );
            let term = value.abs();
            if sum == 0.0 {
                sum = term;
            } else if term != 0.0 {
                sum = (sum + term).next_up();
            }
        }
        bound = bound.max(sum);
    }
    bound
}

// ---------------------------------------------------------------------------
// Chebyshev polynomial smoother
// ---------------------------------------------------------------------------

/// Degree-k Chebyshev smoother targeting the band [λ_max/α, λ_max]:
/// pure SpMV chains, no sequential dependencies — the many-core-honest
/// smoother. Damps the targeted band by the known Chebyshev factor.
#[derive(Debug, Clone)]
pub struct Chebyshev {
    a: Csr,
    degree: usize,
    lo: f64,
    hi: f64,
}

impl Chebyshev {
    /// Build with a conservative absolute-row-sum upper bound and the
    /// standard band divisor α (30.0 is the common default). The historical
    /// power estimate with safety 1.1 is retained when it gives a wider band;
    /// it can never lower the enclosing upper endpoint.
    ///
    /// # Panics
    /// Invalid degree/divisor, non-square/non-finite input, or a zero,
    /// infinite, or unrepresentable positive spectral band. The caller still
    /// owns the SPD assumption.
    #[must_use]
    pub fn new(a: &Csr, degree: usize, alpha: f64) -> Chebyshev {
        assert!(degree >= 1, "Chebyshev degree must be >= 1");
        assert!(
            alpha.is_finite() && alpha > 1.0,
            "band divisor must exceed 1 and be finite"
        );
        let bound = spectral_radius_upper_bound(a);
        assert!(
            bound.is_finite() && bound > 0.0,
            "Chebyshev needs a finite positive spectral bound"
        );
        let estimate = lambda_max_estimate(a, 30, 1.1);
        let hi = if estimate.is_finite() {
            estimate.max(bound)
        } else {
            bound
        };
        let lo = hi / alpha;
        assert!(
            lo > 0.0 && lo < hi,
            "Chebyshev band endpoints must be positive and distinct"
        );
        Chebyshev {
            a: a.clone(),
            degree,
            lo,
            hi,
        }
    }

    /// The band this smoother targets (evidence for ledgering).
    #[must_use]
    pub fn band(&self) -> (f64, f64) {
        (self.lo, self.hi)
    }
}

impl Precond for Chebyshev {
    /// z ≈ A⁻¹·r restricted to the band: the classic three-term Chebyshev
    /// iteration with zero initial guess.
    fn apply(&self, r: &[f64], z: &mut [f64]) {
        let n = r.len();
        let theta = f64::midpoint(self.hi, self.lo);
        let delta = f64::midpoint(self.hi, -self.lo);
        // First step: z1 = r/θ.
        let mut zk: Vec<f64> = r.iter().map(|&v| v / theta).collect();
        let mut zk_prev = vec![0.0f64; n];
        let mut rho_prev = delta / theta;
        let mut az = vec![0.0f64; n];
        for _ in 1..self.degree {
            self.a.spmv(&zk, &mut az);
            // residual of current iterate: rk = r − A·zk.
            let rho = 1.0 / (2.0 * theta / delta - rho_prev);
            let c1 = 2.0 * rho / delta;
            let c2 = rho * rho_prev;
            for i in 0..n {
                let rk = r[i] - az[i];
                let znew = c1.mul_add(rk, zk[i]) + c2 * (zk[i] - zk_prev[i]);
                zk_prev[i] = zk[i];
                zk[i] = znew;
            }
            rho_prev = rho;
        }
        z.copy_from_slice(&zk);
    }
}

// ---------------------------------------------------------------------------
// ILU(0) / IC(0)
// ---------------------------------------------------------------------------

#[path = "precond/ilu0.rs"]
mod incomplete;
pub use incomplete::{Ilu0, IluBreakdown, IluControlledError};

/// Factor. # Errors: [`IluBreakdown`] with the failing row.
pub fn ilu0(a: &Csr) -> Result<Ilu0, IluBreakdown> {
    match ilu0_controlled(a, || Ok::<(), core::convert::Infallible>(())) {
        Ok(factor) => Ok(factor),
        Err(IluControlledError::Breakdown(error)) => Err(error),
        Err(IluControlledError::Capacity) => panic!("ilu0 sparse allocation failed"),
        Err(IluControlledError::Checkpoint(never)) => match never {},
    }
}

/// The same zero-fill factorization with bounded setup checkpoints.
///
/// Finite positive or negative pivots are admitted without shifting. Missing
/// or zero pivots and nonfinite factor arithmetic return typed breakdown;
/// checkpoint refusal publishes no factor and preserves the original error.
/// Storage is O(rows + stored entries). Factor work depends on the original
/// sparse row intersections; no fill or dense matrix is introduced.
pub fn ilu0_controlled<E>(
    a: &Csr,
    checkpoint: impl FnMut() -> Result<(), E>,
) -> Result<Ilu0, IluControlledError<E>> {
    assert_eq!(a.nrows(), a.ncols(), "ilu0 requires a square matrix");
    incomplete::factor(&a.row_ptr, &a.col_idx, &a.vals, checkpoint)
}

impl Precond for Ilu0 {
    /// z = U⁻¹·L⁻¹·r (forward then backward substitution on the pattern).
    fn apply(&self, r: &[f64], z: &mut [f64]) {
        self.apply_controlled(r, z, || Ok::<(), core::convert::Infallible>(()))
            .expect("ilu0 non-finite inverse action");
    }
}

// ---------------------------------------------------------------------------
// PCG reference driver
// ---------------------------------------------------------------------------

/// Result of a PCG solve: the evidence object.
#[derive(Debug, Clone)]
pub struct PcgReport {
    /// Finite iterate updates completed before convergence, cap, or breakdown.
    pub iters: usize,
    /// Recomputed ‖b − A·x‖₂/‖b‖₂ for the returned iterate, never the
    /// recursive CG residual. For a zero RHS this is zero only for a zero
    /// residual, and infinity otherwise. Non-finite residuals give infinity.
    pub rel_residual: f64,
    /// Whether the recomputed residual met the tolerance within the cap.
    pub converged: bool,
}

/// A nonnegative scalar times an exact power of two. Dot products and norms
/// need not themselves fit f64 when their ratios (the CG coefficients) do.
#[derive(Clone, Copy)]
struct PcgMagnitude {
    value: f64,
    exponent: i32,
}

fn pcg_scale(maximum: f64) -> (f64, i32) {
    let bits = (maximum.to_bits() >> 52).max(1);
    (f64::from_bits(bits << 52), bits as i32 - 1023)
}

fn pcg_norm(values: &[f64]) -> Option<PcgMagnitude> {
    let mut squared = 0.0;
    let mut maximum = 0.0_f64;
    for &value in values {
        if !value.is_finite() {
            return None;
        }
        squared += value * value;
        maximum = maximum.max(value.abs());
    }
    if maximum == 0.0 || squared.is_normal() {
        return Some(PcgMagnitude {
            value: squared.sqrt(),
            exponent: 0,
        });
    }
    let (scale, exponent) = pcg_scale(maximum);
    let inverse = scale.recip();
    let squared: f64 = values.iter().map(|v| (v * inverse).powi(2)).sum();
    Some(PcgMagnitude {
        value: squared.sqrt(),
        exponent,
    })
}

fn pcg_positive_dot(a: &[f64], b: &[f64]) -> Option<PcgMagnitude> {
    let mut sum = 0.0;
    let (mut max_a, mut max_b) = (0.0_f64, 0.0_f64);
    for (&a, &b) in a.iter().zip(b) {
        if !a.is_finite() || !b.is_finite() {
            return None;
        }
        sum += a * b;
        max_a = max_a.max(a.abs());
        max_b = max_b.max(b.abs());
    }
    // Preserve the established accumulation and ratio bits at ordinary
    // scales. Subnormal dot products also need scaling: their quantization
    // can otherwise turn a perfectly conditioned solve into a wrong step.
    if sum.is_normal() && sum > 0.0 {
        return Some(PcgMagnitude {
            value: sum,
            exponent: 0,
        });
    }
    if max_a == 0.0 || max_b == 0.0 {
        return None;
    }
    let (scale_a, exponent_a) = pcg_scale(max_a);
    let (scale_b, exponent_b) = pcg_scale(max_b);
    let (inverse_a, inverse_b) = (scale_a.recip(), scale_b.recip());
    let sum: f64 = a
        .iter()
        .zip(b)
        .map(|(&a, &b)| (a * inverse_a) * (b * inverse_b))
        .sum();
    (sum.is_finite() && sum > 0.0).then_some(PcgMagnitude {
        value: sum,
        exponent: exponent_a + exponent_b,
    })
}

fn pcg_normalize(value: f64) -> (f64, i32) {
    let (value, correction) = if value < f64::MIN_POSITIVE {
        (value * f64::MIN_POSITIVE.recip(), -1022)
    } else {
        (value, 0)
    };
    let exponent = ((value.to_bits() >> 52) as i32) - 1023 + correction;
    let mantissa = f64::from_bits((value.to_bits() & ((1_u64 << 52) - 1)) | (1023_u64 << 52));
    (mantissa, exponent)
}

fn pcg_ratio(numerator: PcgMagnitude, denominator: PcgMagnitude) -> f64 {
    if numerator.value == 0.0 {
        return 0.0;
    }
    if denominator.value == 0.0 {
        return f64::INFINITY;
    }
    if numerator.exponent == denominator.exponent {
        return numerator.value / denominator.value;
    }
    let (numerator_value, numerator_exponent) = pcg_normalize(numerator.value);
    let (denominator_value, denominator_exponent) = pcg_normalize(denominator.value);
    let mut value = numerator_value / denominator_value;
    let mut exponent =
        numerator.exponent - denominator.exponent + numerator_exponent - denominator_exponent;
    if value < 1.0 {
        value *= 2.0;
        exponent -= 1;
    }
    // Exact normal power-of-two factors defer subnormal rounding until the
    // final scale. An overflowing/underflowing final ratio stays explicit.
    while exponent > 1023 {
        value *= f64::from_bits(2046_u64 << 52);
        exponent -= 1023;
    }
    while exponent < -1022 {
        value *= f64::MIN_POSITIVE;
        exponent += 1022;
    }
    value * f64::from_bits(((exponent + 1023) as u64) << 52)
}

fn pcg_relative_residual(r: &[f64], bnorm: PcgMagnitude) -> f64 {
    pcg_norm(r).map_or(f64::INFINITY, |norm| {
        let ratio = pcg_ratio(norm, bnorm);
        if ratio == 0.0 && norm.value != 0.0 {
            f64::from_bits(1)
        } else {
            ratio
        }
    })
}

fn pcg_report(
    a: &Csr,
    b: &[f64],
    x: &[f64],
    residual: &mut [f64],
    bnorm: PcgMagnitude,
    iters: usize,
    tol: f64,
) -> PcgReport {
    a.spmv(x, residual);
    for (r, &b) in residual.iter_mut().zip(b) {
        *r = b - *r;
    }
    let rel_residual = pcg_relative_residual(residual, bnorm);
    let converged = if tol == 0.0 {
        residual.iter().all(|&value| value == 0.0)
    } else {
        rel_residual <= tol
    };
    PcgReport {
        iters,
        rel_residual,
        converged,
    }
}

/// Preconditioned conjugate gradients (reference driver — the solver
/// stack supersedes this; kept minimal and deterministic).
///
/// Positive, finite preconditioned curvature is required. Numerical breakdown
/// returns `converged: false` with the last finite iterate and its recomputed
/// residual. Iterate updates are staged before publication. A recursive
/// convergence candidate is checked against `b - A*x`; residual drift restarts
/// the search from that actual residual without resetting the iteration cap.
/// Norms and dot-product ratios use power-of-two scaling when ordinary
/// accumulation would overflow, underflow, or become subnormal. This does not
/// make unrepresentable matrix products or solution values representable.
///
/// # Panics
/// Non-square/mismatched shapes, non-finite RHS/initial iterate, or a negative
/// or non-finite tolerance. Matrix and preconditioner SPD are caller obligations.
pub fn pcg<P: Precond>(
    a: &Csr,
    b: &[f64],
    x: &mut [f64],
    m: &P,
    tol: f64,
    max_iters: usize,
) -> PcgReport {
    let n = b.len();
    assert_eq!(a.nrows(), n, "pcg: rows must match RHS");
    assert_eq!(a.ncols(), n, "pcg requires a square matrix");
    assert_eq!(x.len(), n, "pcg: iterate must match RHS");
    assert!(
        tol.is_finite() && tol >= 0.0,
        "pcg: tolerance must be finite and nonnegative"
    );
    assert!(
        x.iter().all(|v| v.is_finite()),
        "pcg: initial iterate must be finite"
    );
    let bnorm = pcg_norm(b).expect("pcg: RHS must be finite");
    let mut r = vec![0.0f64; n];
    let entry = pcg_report(a, b, x, &mut r, bnorm, 0, tol);
    if entry.converged || max_iters == 0 || !r.iter().all(|v| v.is_finite()) {
        return entry;
    }
    let mut z = vec![0.0f64; n];
    m.apply(&r, &mut z);
    let mut p = z.clone();
    let mut rz = pcg_positive_dot(&r, &z);
    let mut ap = vec![0.0f64; n];
    let mut next_x = vec![0.0f64; n];
    let mut completed = 0;
    for _ in 0..max_iters {
        let Some(current_rz) = rz else {
            break;
        };
        a.spmv(&p, &mut ap);
        let Some(pap) = pcg_positive_dot(&p, &ap) else {
            break;
        };
        let alpha = pcg_ratio(current_rz, pap);
        if !alpha.is_finite() || alpha <= 0.0 {
            break;
        }
        for i in 0..n {
            next_x[i] = alpha.mul_add(p[i], x[i]);
            ap[i] = (-alpha).mul_add(ap[i], r[i]);
        }
        if !next_x.iter().chain(&ap).all(|v| v.is_finite()) {
            break;
        }
        x.copy_from_slice(&next_x);
        std::mem::swap(&mut r, &mut ap);
        completed += 1;
        if pcg_relative_residual(&r, bnorm) <= tol || completed == max_iters {
            let report = pcg_report(a, b, x, &mut ap, bnorm, completed, tol);
            if report.converged || completed == max_iters {
                return report;
            }
            // The recurrence lost the true residual. Restart from the
            // recomputed vector while retaining both x and the spent budget.
            std::mem::swap(&mut r, &mut ap);
            m.apply(&r, &mut z);
            rz = pcg_positive_dot(&r, &z);
            p.copy_from_slice(&z);
            continue;
        }
        m.apply(&r, &mut z);
        let Some(rz_new) = pcg_positive_dot(&r, &z) else {
            break;
        };
        let beta = pcg_ratio(rz_new, current_rz);
        if !beta.is_finite() || beta < 0.0 {
            break;
        }
        for i in 0..n {
            ap[i] = beta.mul_add(p[i], z[i]);
        }
        if !ap.iter().all(|v| v.is_finite()) {
            break;
        }
        std::mem::swap(&mut p, &mut ap);
        rz = Some(rz_new);
    }
    pcg_report(a, b, x, &mut r, bnorm, completed, tol)
}

// ---------------------------------------------------------------------------
// Smoothed-aggregation AMG
// ---------------------------------------------------------------------------

/// One AMG level: the operator, its smoother, and the prolongator down
/// from the finer level (absent on the finest).
struct AmgLevel {
    a: Csr,
    smoother: Chebyshev,
    /// Prolongator from THIS (coarser) level up to the finer one.
    p_from_coarse: Option<Csr>,
}

/// Smoothed-aggregation AMG V-cycle preconditioner.
pub struct SaAmg {
    levels: Vec<AmgLevel>,
    coarse_ilu: Option<Ilu0>,
    /// Per-level unknown counts (evidence: operator complexity).
    pub level_sizes: Vec<usize>,
}

/// Deterministic greedy aggregation on the strength graph: index-order
/// roots, lowest-index attachment (P2 on setup). Returns (aggregate id
/// per node, aggregate count).
fn aggregate(a: &Csr, theta: f64) -> (Vec<usize>, usize) {
    const UNASSIGNED: usize = usize::MAX;
    let n = a.nrows();
    let mut agg = vec![UNASSIGNED; n];
    let strong = |i: usize, j: usize, aij: f64| -> bool {
        let aii = a.get(i, i).abs();
        let ajj = a.get(j, j).abs();
        aij.abs() > theta * (aii * ajj).sqrt()
    };
    let mut count = 0;
    // Pass 1: roots in index order over fully-unaggregated neighborhoods.
    for i in 0..n {
        if agg[i] != UNASSIGNED {
            continue;
        }
        let (cols, vals) = a.row(i);
        let neighborhood_free = cols
            .iter()
            .zip(vals)
            .all(|(&j, &v)| j == i || !strong(i, j, v) || agg[j] == UNASSIGNED);
        if neighborhood_free {
            agg[i] = count;
            for (&j, &v) in cols.iter().zip(vals) {
                if j != i && strong(i, j, v) {
                    agg[j] = count;
                }
            }
            count += 1;
        }
    }
    // Pass 2: attach leftovers to the lowest-id neighboring aggregate.
    for i in 0..n {
        if agg[i] != UNASSIGNED {
            continue;
        }
        let (cols, vals) = a.row(i);
        let mut best = UNASSIGNED;
        for (&j, &v) in cols.iter().zip(vals) {
            if j != i && strong(i, j, v) && agg[j] != UNASSIGNED {
                best = best.min(agg[j]);
            }
        }
        if best == UNASSIGNED {
            // Isolated node: its own aggregate (keeps the prolongator
            // full-rank; happens on disconnected/weak rows).
            agg[i] = count;
            count += 1;
        } else {
            agg[i] = best;
        }
    }
    (agg, count)
}

/// Build the Jacobi-smoothed prolongator P = (I − ω·D⁻¹·A)·P₀ where P₀ is
/// the piecewise-constant tentative prolongator over the aggregates.
fn smoothed_prolongator(a: &Csr, agg: &[usize], n_agg: usize) -> Csr {
    let n = a.nrows();
    // P0 as CSR: one entry per row.
    let mut row_ptr = vec![0usize; n + 1];
    let mut col_idx = Vec::with_capacity(n);
    let mut vals = Vec::with_capacity(n);
    for (i, &g) in agg.iter().enumerate() {
        col_idx.push(g);
        vals.push(1.0);
        row_ptr[i + 1] = i + 1;
    }
    let p0 = Csr::from_parts(n, n_agg, row_ptr, col_idx, vals);
    // The power estimate normally gives the historical Jacobi weight, but
    // it can miss a dominant mode. Admit that weight only when the outward
    // row-sum bound proves omega*rho(D^-1 A) < 2. Otherwise use the standard
    // 4/3 target with the conservative bound. D^-1 A need not be symmetric,
    // so its power-iteration Rayleigh quotient is not a bound either.
    let dinv_a = scale_rows_by_dinv(a);
    let lam = lambda_max_estimate(&dinv_a, 20, 1.0);
    let bound = spectral_radius_upper_bound(&dinv_a);
    assert!(
        bound.is_finite() && bound > 0.0,
        "AMG needs a finite positive scaled spectral bound"
    );
    let estimate = 4.0 / (3.0 * lam);
    let omega = if estimate > 0.0 && (estimate * bound).next_up() < 2.0 {
        estimate
    } else {
        ((4.0 / 3.0) / bound).next_down()
    };
    // P = P0 − ω·(D⁻¹A)·P0.
    let ap0 = spgemm(&dinv_a, &p0);
    // P = P0 + (−ω)·AP0: merge the two patterns.
    add_scaled(&p0, &ap0, -omega)
}

/// D⁻¹·A (rows scaled by the inverse diagonal).
fn scale_rows_by_dinv(a: &Csr) -> Csr {
    let n = a.nrows();
    let mut row_ptr = vec![0usize; n + 1];
    let mut col_idx = Vec::new();
    let mut vals = Vec::new();
    for i in 0..n {
        let d = a.get(i, i);
        let inv = if d == 0.0 { 0.0 } else { 1.0 / d };
        let (cols, vs) = a.row(i);
        for (&c, &v) in cols.iter().zip(vs) {
            col_idx.push(c);
            vals.push(v * inv);
        }
        row_ptr[i + 1] = col_idx.len();
    }
    Csr::from_parts(n, a.ncols(), row_ptr, col_idx, vals)
}

/// x + s·y over merged CSR patterns.
fn add_scaled(x: &Csr, y: &Csr, s: f64) -> Csr {
    let n = x.nrows();
    let mut row_ptr = vec![0usize; n + 1];
    let mut col_idx = Vec::new();
    let mut vals = Vec::new();
    for i in 0..n {
        let (cx, vx) = x.row(i);
        let (cy, vy) = y.row(i);
        let (mut a, mut b) = (0usize, 0usize);
        while a < cx.len() || b < cy.len() {
            let ca = cx.get(a).copied().unwrap_or(usize::MAX);
            let cb = cy.get(b).copied().unwrap_or(usize::MAX);
            let c = ca.min(cb);
            let mut v = 0.0f64;
            if ca == c {
                v += vx[a];
                a += 1;
            }
            if cb == c {
                v = s.mul_add(vy[b], v);
                b += 1;
            }
            col_idx.push(c);
            vals.push(v);
            row_ptr[i + 1] = col_idx.len();
        }
        row_ptr[i + 1] = col_idx.len();
    }
    Csr::from_parts(n, x.ncols(), row_ptr, col_idx, vals)
}

impl SaAmg {
    /// Set up the hierarchy: strength θ, Chebyshev smoother degree, and
    /// the coarsest-size cutoff. Deterministic across reruns (tested).
    #[must_use]
    pub fn new(a: &Csr, theta: f64, smoother_degree: usize) -> SaAmg {
        let mut levels = Vec::new();
        let mut level_sizes = Vec::new();
        let mut current = a.clone();
        let mut p_down: Option<Csr> = None;
        loop {
            let n = current.nrows();
            level_sizes.push(n);
            let smoother = Chebyshev::new(&current, smoother_degree, 30.0);
            let done = n <= 64;
            let (agg, n_agg) = if done {
                (Vec::new(), 0)
            } else {
                aggregate(&current, theta)
            };
            // Stagnation guard: coarsening must actually shrink.
            let stagnated = !done && n_agg * 10 >= n * 9;
            levels.push(AmgLevel {
                a: current.clone(),
                smoother,
                p_from_coarse: p_down.take(),
            });
            if done || stagnated {
                break;
            }
            let p = smoothed_prolongator(&current, &agg, n_agg);
            let pt = transpose(&p);
            let ap = spgemm(&current, &p);
            current = spgemm(&pt, &ap);
            p_down = Some(p);
        }
        // Coarsest-level ILU(0) for the PCG coarse solve (shift-retry on
        // breakdown: SPD Galerkin products rarely need it, but be safe).
        let coarse = &levels.last().unwrap().a;
        let coarse_ilu = ilu0(coarse).ok();
        SaAmg {
            levels,
            coarse_ilu,
            level_sizes,
        }
    }

    /// Operator complexity: Σ nnz(A_ℓ)/nnz(A₀) — the memory-honesty metric.
    #[must_use]
    pub fn operator_complexity(&self) -> f64 {
        let base = self.levels[0].a.nnz().max(1) as f64;
        self.levels.iter().map(|l| l.a.nnz() as f64).sum::<f64>() / base
    }

    fn vcycle(&self, level: usize, r: &[f64], z: &mut [f64]) {
        let lvl = &self.levels[level];
        let n = lvl.a.nrows();
        if level + 1 == self.levels.len() {
            // Coarsest: IC-preconditioned CG to tight tolerance.
            z.fill(0.0);
            if let Some(ilu) = &self.coarse_ilu {
                pcg(&lvl.a, r, z, ilu, 1e-12, 4 * n.max(16));
            } else {
                pcg(&lvl.a, r, z, &IdentityPrecond, 1e-12, 4 * n.max(16));
            }
            return;
        }
        // Pre-smooth.
        lvl.smoother.apply(r, z);
        // Residual: rr = r − A z.
        let mut az = vec![0.0f64; n];
        lvl.a.spmv(z, &mut az);
        let rr: Vec<f64> = r.iter().zip(&az).map(|(a, b)| a - b).collect();
        // Restrict with Pᵀ (the coarser level's stored prolongator).
        let p = self.levels[level + 1]
            .p_from_coarse
            .as_ref()
            .expect("non-coarsest level must carry a prolongator");
        let nc = p.ncols();
        let pt = transpose(p);
        let mut rc = vec![0.0f64; nc];
        pt.spmv(&rr, &mut rc);
        // Coarse correction.
        let mut zc = vec![0.0f64; nc];
        self.vcycle(level + 1, &rc, &mut zc);
        // Prolongate and add.
        let mut up = vec![0.0f64; n];
        p.spmv(&zc, &mut up);
        for (zi, ui) in z.iter_mut().zip(&up) {
            *zi += ui;
        }
        // Post-smooth on the corrected residual.
        lvl.a.spmv(z, &mut az);
        let r2: Vec<f64> = r.iter().zip(&az).map(|(a, b)| a - b).collect();
        let mut dz = vec![0.0f64; n];
        lvl.smoother.apply(&r2, &mut dz);
        for (zi, di) in z.iter_mut().zip(&dz) {
            *zi += di;
        }
    }
}

impl Precond for SaAmg {
    fn apply(&self, r: &[f64], z: &mut [f64]) {
        self.vcycle(0, r, z);
    }
}

#[cfg(test)]
mod spectral_safety_tests {
    use super::{
        Chebyshev, Precond, lambda_max_estimate, smoothed_prolongator, spectral_radius_upper_bound,
    };
    use crate::{Coo, Csr};

    #[test]
    fn g0_row_sum_bound_rounds_outward_and_handles_zero() {
        let empty = Csr::from_parts(0, 0, vec![0], vec![], vec![]);
        assert_eq!(spectral_radius_upper_bound(&empty), 0.0);
        let zero = Csr::from_parts(2, 2, vec![0, 0, 0], vec![], vec![]);
        assert_eq!(spectral_radius_upper_bound(&zero), 0.0);
        let mut coo = Coo::new(2, 2);
        coo.push(0, 0, 1.0);
        coo.push(0, 1, 2.0_f64.powi(-53));
        coo.push(1, 1, 1.0);
        assert!(spectral_radius_upper_bound(&coo.assemble()) >= 1.0_f64.next_up());
    }

    #[test]
    fn g0_chebyshev_damps_a_dominant_mode_missed_by_power_iteration() {
        // The fixed power seed is (1,b), the eigenvector for lambda=1.
        // Its perpendicular (-b,1) has lambda=2. Thirty iterations cannot
        // recover from the near-zero dominant projection introduced by f64
        // rounding; multiplying the estimate by 1.1 is still unsafe.
        let b = 1.0 + 0.3 * ((2_654_435_761_usize % 97) as f64) / 97.0;
        let d = 1.0 + b * b;
        let mut coo = Coo::new(2, 2);
        coo.push(0, 0, (1.0 + 2.0 * b * b) / d);
        coo.push(0, 1, -b / d);
        coo.push(1, 0, -b / d);
        coo.push(1, 1, (b * b + 2.0) / d);
        let a = coo.assemble();
        assert!(lambda_max_estimate(&a, 30, 1.1) < 1.2);
        let smoother = Chebyshev::new(&a, 8, 30.0);
        assert!(smoother.band().1 >= 2.0);
        let rhs = [-b, 1.0];
        let mut result = [0.0; 2];
        smoother.apply(&rhs, &mut result);
        let mut action = [0.0; 2];
        a.spmv(&result, &mut action);
        let residual_squared = rhs
            .iter()
            .zip(action)
            .map(|(&r, ax)| (r - ax).powi(2))
            .sum::<f64>();
        // The old band magnifies this residual by over 27,000, despite SPD.
        assert!(residual_squared < 0.04 * d, "{result:?}, {action:?}");
    }

    #[test]
    fn g3_chebyshev_uses_the_bound_when_power_normalization_overflows_or_underflows() {
        for scale in [1e-300, 1e300] {
            let a = Csr::from_parts(1, 1, vec![0, 1], vec![0], vec![scale]);
            let smoother = Chebyshev::new(&a, 8, 30.0);
            assert!(smoother.band().1.is_finite());
            assert!(smoother.band().1 >= scale);
            let mut result = [0.0];
            smoother.apply(&[scale], &mut result);
            assert!((result[0] - 1.0).abs() < 0.2, "scale {scale}, {result:?}");
        }
    }

    #[test]
    fn g0_amg_prolongator_damps_a_mode_missed_by_its_power_seed() {
        // The seed repeats after 97 entries. The antisymmetric (0,97) mode
        // of this SPD block has lambda=1.99 and nearly zero initial weight.
        let n = 98;
        let mut coo = Coo::new(n, n);
        for i in 0..n {
            coo.push(i, i, 1.0);
        }
        coo.push(0, n - 1, -0.99);
        coo.push(n - 1, 0, -0.99);
        let a = coo.assemble();
        assert!(lambda_max_estimate(&a, 20, 1.0) < 1.1);
        let aggregates: Vec<_> = (0..n).collect();
        let p = smoothed_prolongator(&a, &aggregates, n);
        let mut mode = vec![0.0; n];
        mode[0] = 1.0;
        mode[n - 1] = -1.0;
        let mut smoothed = vec![0.0; n];
        p.spmv(&mode, &mut smoothed);
        assert!(smoothed.iter().map(|x| x * x).sum::<f64>() < 2.0);
    }
}
