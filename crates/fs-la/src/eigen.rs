//! Symmetric eigensolvers (plan §6.1): matrix-free Lanczos and LOBPCG for
//! extremal eigenpairs, plus a dense cyclic-Jacobi symmetric eigensolver
//! (the small-problem workhorse and the Rayleigh–Ritz kernel).
//!
//! MATRIX-FREE: operators are closures `apply(x, y)` writing y = A·x —
//! fs-sparse SpMV, stencils, and composed operators all plug in without
//! coupling fs-la to any matrix format.
//!
//! RESUMABLE: both iterative solvers expose plain-data state objects;
//! `checkpoint = clone`, and resuming N+M steps is BITWISE identical to
//! running N+M straight (tested) — the P2 property extended to
//! eigeniteration.
//!
//! Honesty: every returned [`EigenPair`] carries the numerically recomputed
//! residual norm ‖A·v − λ·v‖₂ rather than an internal recurrence estimate.
//! Ordinary `f64` arithmetic and a caller-supplied operator do not by
//! themselves confer certificate authority; higher layers must admit the
//! operator proposition and evidence before making a scientific claim.
//!
//! v1 notes: Lanczos uses FULL reorthogonalization (the conservative
//! special case of selective reorth; the ω-recurrence optimization is a
//! recorded refinement). LOBPCG uses the X−X_prev conjugate direction and
//! identity preconditioning by default. Eigenvector adjoints are the
//! recorded follow-up; dλ/dp = vᵀ(∂A/∂p)v composes caller-side with
//! fs-ad duals.

use crate::factor::qr;

/// Conservative cap for temporary `f64` elements admitted by the dense
/// eigensolver kernels. This is an OVERFLOW/OOM guard, not a promise that an
/// admitted problem will fit the caller's memory budget: higher layers still
/// have to account for retained state, allocator overhead, and concurrent
/// work. The cap is intentionally expressed in elements so every shape
/// product is checked before bytes or allocations are formed.
pub const MAX_EIGEN_WORK_ELEMENTS: usize = 64 * 1024 * 1024;

/// Schema version for [`JacobiEighAdmission`].
pub const JACOBI_EIGH_ADMISSION_SCHEMA_VERSION: u32 = 1;

/// Conservative cap for scalar work performed between cooperative polls by
/// the service wrapper. It bounds in-house algebra only; a caller-supplied
/// operator must separately honor its own bounded-call contract.
pub const MAX_EIGEN_UNPOLLED_SCALAR_WORK: usize = 128 * 1024 * 1024;

/// Sealed evidence that one dense Jacobi eigendecomposition fits the checked
/// shape arithmetic and practical aggregate-work cap.
///
/// This receipt is allocation-free structural evidence. It does not reserve
/// memory, prove that the allocator can satisfy the request, or make a
/// bounded-latency claim for the cubic sweep work.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct JacobiEighAdmission {
    schema_version: u32,
    dimension: usize,
    matrix_entries: usize,
    aggregate_work_elements: usize,
    work_element_cap: usize,
}

impl JacobiEighAdmission {
    /// Admission schema version.
    #[must_use]
    pub const fn schema_version(self) -> u32 {
        self.schema_version
    }

    /// Admitted square-matrix dimension.
    #[must_use]
    pub const fn dimension(self) -> usize {
        self.dimension
    }

    /// Checked `n * n` entry count used by every dense square allocation.
    #[must_use]
    pub const fn matrix_entries(self) -> usize {
        self.matrix_entries
    }

    /// Conservative aggregate element envelope `4 * n * n + 3 * n`.
    #[must_use]
    pub const fn aggregate_work_elements(self) -> usize {
        self.aggregate_work_elements
    }

    /// Practical aggregate element cap applied by this admission schema.
    #[must_use]
    pub const fn work_element_cap(self) -> usize {
        self.work_element_cap
    }
}

/// Allocation-free refusal from [`admit_jacobi_eigh`].
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum JacobiEighAdmissionError {
    /// The square matrix shape `n * n` left the target's `usize` domain.
    MatrixShapeOverflow {
        /// Requested matrix dimension.
        dimension: usize,
    },
    /// The aggregate formula `4 * n * n + 3 * n` left `usize`.
    AggregateWorkOverflow {
        /// Requested matrix dimension.
        dimension: usize,
    },
    /// The representable aggregate formula exceeds the practical work cap.
    WorkCapExceeded {
        /// Requested matrix dimension.
        dimension: usize,
        /// Exact representable aggregate element requirement.
        required_elements: usize,
        /// Aggregate element cap applied by this schema.
        cap_elements: usize,
    },
}

impl core::fmt::Display for JacobiEighAdmissionError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match *self {
            Self::MatrixShapeOverflow { dimension } => {
                write!(
                    formatter,
                    "Jacobi matrix shape {dimension} * {dimension} overflowed"
                )
            }
            Self::AggregateWorkOverflow { dimension } => write!(
                formatter,
                "Jacobi aggregate work overflowed for dimension {dimension}"
            ),
            Self::WorkCapExceeded {
                dimension,
                required_elements,
                cap_elements,
            } => write!(
                formatter,
                "Jacobi dimension {dimension} needs {required_elements} aggregate elements but the practical cap is {cap_elements}"
            ),
        }
    }
}

impl std::error::Error for JacobiEighAdmissionError {}

/// Check the complete dense Jacobi shape and aggregate-work envelope without
/// allocating.
///
/// Arithmetic refusal precedes the practical cap so callers can distinguish
/// an unrepresentable proposition from a representable but inadmissible one.
/// Dimension zero retains the low-level kernel's established empty-matrix
/// semantics.
///
/// # Errors
/// Returns [`JacobiEighAdmissionError::MatrixShapeOverflow`] when `n * n`
/// cannot be represented, [`JacobiEighAdmissionError::AggregateWorkOverflow`]
/// when `4 * n * n + 3 * n` cannot be represented, or
/// [`JacobiEighAdmissionError::WorkCapExceeded`] when that aggregate exceeds
/// [`MAX_EIGEN_WORK_ELEMENTS`].
#[must_use]
pub fn admit_jacobi_eigh(n: usize) -> Result<JacobiEighAdmission, JacobiEighAdmissionError> {
    let matrix_entries = n
        .checked_mul(n)
        .ok_or(JacobiEighAdmissionError::MatrixShapeOverflow { dimension: n })?;
    let dense_work = matrix_entries
        .checked_mul(4)
        .ok_or(JacobiEighAdmissionError::AggregateWorkOverflow { dimension: n })?;
    let linear_work = n
        .checked_mul(3)
        .ok_or(JacobiEighAdmissionError::AggregateWorkOverflow { dimension: n })?;
    let aggregate_work_elements = dense_work
        .checked_add(linear_work)
        .ok_or(JacobiEighAdmissionError::AggregateWorkOverflow { dimension: n })?;
    if aggregate_work_elements > MAX_EIGEN_WORK_ELEMENTS {
        return Err(JacobiEighAdmissionError::WorkCapExceeded {
            dimension: n,
            required_elements: aggregate_work_elements,
            cap_elements: MAX_EIGEN_WORK_ELEMENTS,
        });
    }
    Ok(JacobiEighAdmission {
        schema_version: JACOBI_EIGH_ADMISSION_SCHEMA_VERSION,
        dimension: n,
        matrix_entries,
        aggregate_work_elements,
        work_element_cap: MAX_EIGEN_WORK_ELEMENTS,
    })
}

fn checked_work_len(rows: usize, cols: usize) -> Option<usize> {
    rows.checked_mul(cols)
        .filter(|&len| len <= MAX_EIGEN_WORK_ELEMENTS)
}

fn checked_lanczos_storage(n: usize, m: usize, k_want: usize) -> Option<()> {
    let pairs = k_want.min(m);
    let basis_elements = checked_work_len(n, m)?;
    let pair_elements = checked_work_len(n, pairs)?;
    let dense_elements = checked_work_len(m, m)?;
    // Count both accepted state and the service's rollback snapshot, plus
    // current/output Ritz vectors, dense extraction work, and linear scratch.
    let aggregate_elements = basis_elements
        .checked_mul(2)?
        .checked_add(pair_elements.checked_mul(2)?)?
        .checked_add(dense_elements.checked_mul(4)?)?
        .checked_add(n.checked_mul(4)?)?;
    if aggregate_elements > MAX_EIGEN_WORK_ELEMENTS {
        return None;
    }
    Some(())
}

fn checked_lanczos_work(n: usize, m: usize, k_want: usize) -> Option<()> {
    checked_lanczos_storage(n, m, k_want)?;
    let pairs = k_want.min(m);
    let basis_elements = checked_work_len(n, m)?;
    let pair_elements = checked_work_len(n, pairs)?;
    let dense_elements = checked_work_len(m, m)?;

    // Bound one service step at the largest dimension reached in the tick:
    // two-pass full reorthogonalization, Ritz-vector reconstruction, residual
    // evaluation, and all 60 cyclic-Jacobi sweeps. The 1024 multiplier is a
    // conservative scalar-operation allowance per Jacobi matrix row.
    let basis_factor = pairs.checked_mul(2)?.checked_add(12)?;
    let vector_work = basis_elements
        .checked_mul(basis_factor)?
        .checked_add(pair_elements.checked_mul(8)?)?;
    let dense_work = dense_elements.checked_mul(m)?.checked_mul(1024)?;
    let scalar_work = vector_work.checked_add(dense_work)?;
    (scalar_work <= MAX_EIGEN_UNPOLLED_SCALAR_WORK).then_some(())
}

/// One converged (or best-effort) eigenpair with its recomputed residual.
#[derive(Debug, Clone)]
pub struct EigenPair {
    /// Ritz value.
    pub value: f64,
    /// Ritz vector (unit 2-norm).
    pub vector: Vec<f64>,
    /// Numerically recomputed residual ‖A·v − λ·v‖₂.
    pub residual: f64,
}

// ---------------------------------------------------------------------------
// Dense cyclic Jacobi (symmetric)
// ---------------------------------------------------------------------------

/// Eigendecomposition of a dense SYMMETRIC matrix (row-major n×n) by
/// cyclic Jacobi: returns (eigenvalues ascending, eigenvectors as columns
/// of a row-major n×n matrix). Deterministic sweep order; ascending sort
/// with lowest-index tie-break (P2).
///
/// This trusted low-level kernel caps aggregate workspace but has no `Cx` and
/// makes no bounded-latency claim for its cubic scalar work. Resumable service
/// callers apply a separate checked cubic-work admission before reaching it.
///
/// # Panics
/// The trusted low-level kernel panics if its shape or aggregate-work
/// arithmetic overflows, its aggregate dense workspace exceeds the practical
/// cap, or `a.len()` does not equal the admitted `n * n` shape.
#[must_use]
pub fn jacobi_eigh(a: &[f64], n: usize) -> (Vec<f64>, Vec<f64>) {
    let admission =
        admit_jacobi_eigh(n).expect("Jacobi aggregate workspace must fit the practical work cap");
    let square = admission.matrix_entries();
    assert_eq!(a.len(), square, "a must be n*n = {square}");
    let mut m = a.to_vec();
    let mut v = vec![0.0f64; square];
    for i in 0..n {
        v[i * n + i] = 1.0;
    }
    for _sweep in 0..60 {
        let mut rotated = false;
        for p in 0..n {
            for q in p + 1..n {
                let apq = m[p * n + q];
                let scale = (m[p * n + p].abs() + m[q * n + q].abs()).max(f64::MIN_POSITIVE);
                if apq.abs() <= 1e-16 * scale {
                    continue;
                }
                rotated = true;
                let theta = (m[q * n + q] - m[p * n + p]) / (2.0 * apq);
                let t = if theta >= 0.0 {
                    1.0 / (theta + (1.0 + theta * theta).sqrt())
                } else {
                    -1.0 / (-theta + (1.0 + theta * theta).sqrt())
                };
                let c = 1.0 / (1.0 + t * t).sqrt();
                let s = c * t;
                // Update rows/cols p and q of the symmetric matrix.
                for k in 0..n {
                    let (mkp, mkq) = (m[k * n + p], m[k * n + q]);
                    m[k * n + p] = c.mul_add(mkp, -(s * mkq));
                    m[k * n + q] = s.mul_add(mkp, c * mkq);
                }
                for k in 0..n {
                    let (mpk, mqk) = (m[p * n + k], m[q * n + k]);
                    m[p * n + k] = c.mul_add(mpk, -(s * mqk));
                    m[q * n + k] = s.mul_add(mpk, c * mqk);
                }
                for k in 0..n {
                    let (vkp, vkq) = (v[k * n + p], v[k * n + q]);
                    v[k * n + p] = c.mul_add(vkp, -(s * vkq));
                    v[k * n + q] = s.mul_add(vkp, c * vkq);
                }
            }
        }
        if !rotated {
            break;
        }
    }
    // Sort ascending, lowest original index first on ties.
    let mut order: Vec<usize> = (0..n).collect();
    let vals: Vec<f64> = (0..n).map(|i| m[i * n + i]).collect();
    order.sort_by(|&x, &y| vals[x].total_cmp(&vals[y]).then(x.cmp(&y)));
    let mut evals = vec![0.0f64; n];
    let mut evecs = vec![0.0f64; square];
    for (slot, &src) in order.iter().enumerate() {
        evals[slot] = vals[src];
        for i in 0..n {
            evecs[i * n + slot] = v[i * n + src];
        }
    }
    (evals, evecs)
}

// ---------------------------------------------------------------------------
// Lanczos (full reorthogonalization)
// ---------------------------------------------------------------------------

/// Resumable Lanczos state: plain data; `clone()` IS a checkpoint.
#[derive(Debug, Clone)]
pub struct LanczosState {
    n: usize,
    basis: Vec<Vec<f64>>,
    alphas: Vec<f64>,
    betas: Vec<f64>,
    /// Residual direction for the next expansion step.
    next: Vec<f64>,
    /// Sticky breakdown marker: the Krylov space is invariant and MUST
    /// NOT be re-expanded (re-entering with the zeroed residual would
    /// push a zero basis vector and corrupt the Ritz extraction).
    exhausted: bool,
}

impl LanczosState {
    /// Allocation-free admission decision for a fresh state's first bounded
    /// advance and Ritz extraction at dimension `n`.
    #[must_use]
    pub fn initial_work_is_admitted(n: usize, additional_steps: usize, k_want: usize) -> bool {
        n > 0
            && additional_steps > 0
            && n <= MAX_EIGEN_WORK_ELEMENTS
            && checked_lanczos_work(n, additional_steps.min(n), k_want).is_some()
    }

    /// Fresh state with the DETERMINISTIC start vector (fixed pattern,
    /// normalized — no RNG, replay-stable by construction).
    ///
    /// # Panics
    /// Use [`Self::try_new`] for untrusted dimensions; this trusted constructor
    /// panics when its checked practical-work preflight refuses.
    #[must_use]
    pub fn new(n: usize) -> LanczosState {
        Self::try_new(n).expect("Lanczos dimension must fit the practical work cap")
    }

    /// Checked counterpart to [`Self::new`]. Refuses a zero or impractically
    /// large dimension before allocating the start vector. Admission is only
    /// a shape guard; it is not a complete memory-budget guarantee.
    #[must_use]
    pub fn try_new(n: usize) -> Option<LanczosState> {
        if n == 0 || n > MAX_EIGEN_WORK_ELEMENTS {
            return None;
        }
        // fs-math strict sin, NOT std sin: the start vector must be
        // cross-ISA bit-identical (platform libm here cost us a golden-hash
        // divergence on x86 — caught by the trj verification pipeline).
        let mut v0: Vec<f64> = (0..n)
            .map(|i| 1.0 + 0.5 * fs_math::det::sin((i as f64) * 0.7))
            .collect();
        let nrm = norm2(&v0);
        for x in &mut v0 {
            *x /= nrm;
        }
        Some(LanczosState {
            n,
            basis: Vec::new(),
            alphas: Vec::new(),
            betas: Vec::new(),
            next: v0,
            exhausted: false,
        })
    }

    /// Fresh state seeded with a CALLER-SUPPLIED start direction
    /// (normalized here) — the warm-start hook for parameter
    /// continuation: seed with the previous parameter point's Ritz
    /// vector. Returns `None` for an empty, non-finite, or
    /// numerically zero seed instead of guessing.
    #[must_use]
    pub fn with_start(v0: &[f64]) -> Option<LanczosState> {
        let n = v0.len();
        if n == 0 || n > MAX_EIGEN_WORK_ELEMENTS || v0.iter().any(|x| !x.is_finite()) {
            return None;
        }
        let start = normalized_copy(v0)?;
        Some(LanczosState {
            n,
            basis: Vec::new(),
            alphas: Vec::new(),
            betas: Vec::new(),
            next: start,
            exhausted: false,
        })
    }

    /// Krylov dimension built so far.
    #[must_use]
    pub fn dim(&self) -> usize {
        self.basis.len()
    }

    /// Whether the Krylov space is exhausted (invariant subspace
    /// found); further expansion requests are no-ops.
    #[must_use]
    pub fn exhausted(&self) -> bool {
        self.exhausted
    }

    /// Whether advancing by at most `additional_steps` and extracting
    /// `k_want` pairs stays within the aggregate-memory and unpolled scalar
    /// work caps. The estimate is conservative and allocation-free.
    #[must_use]
    pub fn work_is_admitted(&self, additional_steps: usize, k_want: usize) -> bool {
        let target = if self.exhausted {
            self.basis.len()
        } else {
            match self.basis.len().checked_add(additional_steps) {
                Some(target) => target.min(self.n),
                None => return false,
            }
        };
        checked_lanczos_work(self.n, target, k_want).is_some()
    }
}

fn norm2(v: &[f64]) -> f64 {
    // Preserve the established operation bits on the ordinary finite path.
    // Only fall back to scaled accumulation when naive squaring overflowed or
    // flushed a genuinely nonzero vector to zero.
    let naive_squared = v.iter().map(|x| x * x).sum::<f64>();
    if naive_squared.is_finite() && (naive_squared > 0.0 || v.iter().all(|x| *x == 0.0)) {
        return naive_squared.sqrt();
    }
    scaled_norm2(v.iter().copied())
}

fn norm2_values(values: impl IntoIterator<Item = f64>) -> f64 {
    // Residual loops historically used mul_add accumulation. Retain those
    // ordinary-path bits while simultaneously preparing an xLASSQ fallback
    // for overflow or underflow.
    let mut naive_squared = 0.0f64;
    let mut any_nonzero = false;
    let mut scale = 0.0f64;
    let mut ssq = 1.0f64;
    for x in values {
        if x.is_nan() {
            return f64::NAN;
        }
        if x.is_infinite() {
            return f64::INFINITY;
        }
        naive_squared = x.mul_add(x, naive_squared);
        let ax = x.abs();
        if ax == 0.0 {
            continue;
        }
        any_nonzero = true;
        if scale < ax {
            let ratio = scale / ax;
            ssq = 1.0 + ssq * ratio * ratio;
            scale = ax;
        } else {
            let ratio = ax / scale;
            ssq += ratio * ratio;
        }
    }
    if naive_squared.is_finite() && (naive_squared > 0.0 || !any_nonzero) {
        naive_squared.sqrt()
    } else if scale == 0.0 {
        0.0
    } else {
        scale * ssq.sqrt()
    }
}

fn scaled_norm2(values: impl IntoIterator<Item = f64>) -> f64 {
    // LAPACK xLASSQ-style scaled sum of squares: finite input does not
    // overflow/underflow merely because its representation is very large or
    // very small. The final norm can still be infinite when the mathematical
    // norm itself is outside f64's range.
    let mut scale = 0.0f64;
    let mut ssq = 1.0f64;
    for x in values {
        if x.is_nan() {
            return f64::NAN;
        }
        if x.is_infinite() {
            return f64::INFINITY;
        }
        let ax = x.abs();
        if ax == 0.0 {
            continue;
        }
        if scale < ax {
            let ratio = scale / ax;
            ssq = 1.0 + ssq * ratio * ratio;
            scale = ax;
        } else {
            let ratio = ax / scale;
            ssq += ratio * ratio;
        }
    }
    if scale == 0.0 {
        0.0
    } else {
        scale * ssq.sqrt()
    }
}

fn max_abs(v: &[f64]) -> f64 {
    v.iter().fold(0.0f64, |scale, x| scale.max(x.abs()))
}

fn normalized_copy(v: &[f64]) -> Option<Vec<f64>> {
    let nrm = norm2(v);
    if nrm.is_finite() && nrm > 0.0 {
        return Some(v.iter().map(|x| x / nrm).collect());
    }
    let scale = max_abs(v);
    if !(scale.is_finite() && scale > 0.0) {
        return None;
    }
    let mut out: Vec<f64> = v.iter().map(|x| x / scale).collect();
    let nrm = norm2(&out);
    if !(nrm.is_finite() && nrm > 0.0) {
        return None;
    }
    for x in &mut out {
        *x /= nrm;
    }
    Some(out)
}

fn norm_relative_to(v: &[f64], scale: f64) -> f64 {
    if scale == 0.0 {
        return norm2(v);
    }
    v.iter()
        .map(|x| {
            let relative = x / scale;
            relative * relative
        })
        .sum::<f64>()
        .sqrt()
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    let mut acc = 0.0f64;
    for (x, y) in a.iter().zip(b) {
        acc = x.mul_add(*y, acc);
    }
    acc
}

/// Advance Lanczos by `steps` expansions (full reorthogonalization —
/// deterministic, conservative), then extract the `k_want` extremal Ritz
/// pairs (`largest` picks which end). Residuals in the returned pairs are
/// recomputed through the supplied operator rather than inferred from the
/// recurrence.
///
/// # Panics
/// This trusted low-level entry point panics when the requested advance would
/// exceed its checked aggregate-memory cap. It has no cancellation context;
/// callers that need bounded cooperative work should use the spectral service,
/// which applies [`LanczosState::work_is_admitted`] as a typed preflight.
pub fn lanczos_run<Op>(
    op: &Op,
    state: &mut LanczosState,
    steps: usize,
    k_want: usize,
    largest: bool,
) -> Vec<EigenPair>
where
    Op: Fn(&[f64], &mut [f64]),
{
    assert!(
        checked_lanczos_storage(
            state.n,
            state.basis.len().saturating_add(steps).min(state.n),
            k_want,
        )
        .is_some(),
        "Lanczos advance exceeds aggregate memory cap"
    );
    let n = state.n;
    let mut w = vec![0.0f64; n];
    for _ in 0..steps {
        if state.exhausted {
            break;
        }
        let v = state.next.clone();
        op(&v, &mut w);
        let action_scale = max_abs(&w);
        let action_norm_scaled = norm_relative_to(&w, action_scale);
        let alpha = dot(&v, &w);
        for (wi, &vi) in w.iter_mut().zip(&v) {
            *wi = (-alpha).mul_add(vi, *wi);
        }
        if let Some(prev) = state.basis.last() {
            let beta_prev = *state.betas.last().unwrap_or(&0.0);
            for (wi, &pi) in w.iter_mut().zip(prev) {
                *wi = (-beta_prev).mul_add(pi, *wi);
            }
        }
        state.basis.push(v);
        state.alphas.push(alpha);
        // FULL reorthogonalization (two passes — "twice is enough").
        for _ in 0..2 {
            for b in &state.basis {
                let c = dot(&w, b);
                for (wi, &bi) in w.iter_mut().zip(b) {
                    *wi = (-c).mul_add(bi, *wi);
                }
            }
        }
        let beta = norm2(&w);
        let relative_beta = if action_scale == 0.0 || action_norm_scaled == 0.0 {
            beta
        } else {
            norm_relative_to(&w, action_scale) / action_norm_scaled
        };
        if beta == 0.0 || relative_beta <= 64.0 * f64::EPSILON {
            // Invariant subspace found: mark the state STICKILY
            // exhausted so a later call cannot re-enter and push the
            // zeroed residual as a basis vector.
            state.betas.push(0.0);
            state.next = vec![0.0; n];
            state.exhausted = true;
            break;
        }
        state.betas.push(beta);
        state.next = normalized_copy(&w).unwrap_or_else(|| vec![f64::NAN; n]);
    }
    ritz_pairs(op, state, k_want, largest)
}

/// Rayleigh–Ritz on the current Krylov basis.
fn ritz_pairs<Op>(op: &Op, state: &LanczosState, k_want: usize, largest: bool) -> Vec<EigenPair>
where
    Op: Fn(&[f64], &mut [f64]),
{
    let m = state.basis.len();
    if m == 0 {
        return Vec::new();
    }
    // Dense tridiagonal T.
    let square =
        checked_work_len(m, m).expect("Lanczos Ritz shape must fit the practical dense-work cap");
    let mut t = vec![0.0f64; square];
    for i in 0..m {
        t[i * m + i] = state.alphas[i];
        if i + 1 < m {
            t[i * m + i + 1] = state.betas[i];
            t[(i + 1) * m + i] = state.betas[i];
        }
    }
    let (evals, evecs) = jacobi_eigh(&t, m);
    let n = state.n;
    let take = k_want.min(m);
    let idx: Vec<usize> = if largest {
        (m - take..m).rev().collect()
    } else {
        (0..take).collect()
    };
    let mut out = Vec::with_capacity(take);
    let mut av = vec![0.0f64; n];
    for &j in &idx {
        let mut y = vec![0.0f64; n];
        for (i, b) in state.basis.iter().enumerate() {
            let c = evecs[i * m + j];
            for (yk, &bk) in y.iter_mut().zip(b) {
                *yk = c.mul_add(bk, *yk);
            }
        }
        let nrm = norm2(&y);
        for v in &mut y {
            *v /= nrm;
        }
        op(&y, &mut av);
        let lam = evals[j];
        let residual = norm2_values(av.iter().zip(&y).map(|(a, &v)| (-lam).mul_add(v, *a)));
        out.push(EigenPair {
            value: lam,
            vector: y,
            residual,
        });
    }
    out
}

// ---------------------------------------------------------------------------
// LOBPCG
// ---------------------------------------------------------------------------

/// Resumable LOBPCG state (plain data; `clone()` IS a checkpoint).
#[derive(Debug, Clone)]
pub struct LobpcgState {
    n: usize,
    b: usize,
    /// Current block X, row-major n×b, orthonormal columns.
    x: Vec<f64>,
    /// Previous block (for the conjugate direction); empty on iteration 0.
    x_prev: Vec<f64>,
    /// Iterations completed.
    pub iters: usize,
}

#[derive(Debug, Clone, Copy)]
struct LobpcgShape {
    two_b: usize,
    three_b: usize,
    n_b: usize,
    n_two_b: usize,
    n_three_b: usize,
}

fn checked_lobpcg_shape(n: usize, b: usize) -> Option<LobpcgShape> {
    if b == 0 {
        return None;
    }
    let two_b = b.checked_mul(2)?;
    let three_b = b.checked_mul(3)?;
    if three_b > n {
        return None;
    }
    let n_b = checked_work_len(n, b)?;
    let n_two_b = checked_work_len(n, two_b)?;
    let n_three_b = checked_work_len(n, three_b)?;
    let reduced_square = checked_work_len(three_b, three_b)?;

    // Conservative simultaneous-work upper bound: retained/current blocks
    // plus the service rollback snapshot, residual/output blocks, Z plus its
    // rank-scratch/Q/AZ copies, and the reduced Rayleigh--Ritz matrix. This
    // deliberately overcounts iteration one. It guards shape arithmetic and
    // impractical work; it does not replace the caller's operation-memory
    // budget.
    let aggregate = n_b
        .checked_mul(4)?
        .checked_add(n_two_b)?
        .checked_add(n_three_b.checked_mul(6)?)?
        .checked_add(reduced_square.checked_mul(3)?)?;
    if aggregate > MAX_EIGEN_WORK_ELEMENTS {
        return None;
    }
    let tall_skinny_work = n_three_b.checked_mul(three_b)?.checked_mul(32)?;
    // Cyclic Jacobi performs at most 60 sweeps, each visiting O((3b)^2)
    // rotations whose row/column/eigenvector updates are O(3b). Use a
    // deliberately loose constant so this remains a real unpolled-work cap.
    let reduced_work = reduced_square.checked_mul(three_b)?.checked_mul(1024)?;
    let scalar_work = tall_skinny_work.checked_add(reduced_work)?;
    if scalar_work > MAX_EIGEN_UNPOLLED_SCALAR_WORK {
        return None;
    }
    Some(LobpcgShape {
        two_b,
        three_b,
        n_b,
        n_two_b,
        n_three_b,
    })
}

impl LobpcgState {
    /// Whether all LOBPCG search-space products and the conservative aggregate
    /// work estimate pass the constructor preflight. This performs no
    /// allocation and makes the same decision as [`Self::try_new`] and
    /// [`Self::with_block`].
    #[must_use]
    pub fn shape_is_admitted(n: usize, b: usize) -> bool {
        checked_lobpcg_shape(n, b).is_some()
    }

    /// Deterministic orthonormal start block (fixed pattern, QR-cleaned).
    ///
    /// # Panics
    /// Use [`Self::try_new`] for untrusted dimensions; this trusted constructor
    /// panics when any checked shape or aggregate-work preflight refuses.
    #[must_use]
    pub fn new(n: usize, b: usize) -> LobpcgState {
        Self::try_new(n, b)
            .expect("LOBPCG shape must satisfy 1 <= b, 3b <= n, and the practical work cap")
    }

    /// Checked counterpart to [`Self::new`]. Every `b`, `2b`, `3b`, `n*b`,
    /// `n*(2b)`, `n*(3b)`, and `(3b)^2` shape is preflighted before the first
    /// allocation. A conservative `n*(3b)^2`/reduced-solve scalar-work bound
    /// also caps time between service polls. These practical caps are refusal
    /// guards, not a complete memory-budget or convergence claim.
    #[must_use]
    pub fn try_new(n: usize, b: usize) -> Option<LobpcgState> {
        let shape = checked_lobpcg_shape(n, b)?;
        let mut x0 = vec![0.0f64; shape.n_b];
        for j in 0..b {
            for i in 0..n {
                // Strict sin (see LanczosState::new — same libm hazard).
                x0[i * b + j] = 1.0
                    + fs_math::det::sin((i as f64) * 0.31 + (j as f64) * 1.7)
                    + if i % (j + 2) == 0 { 0.5 } else { 0.0 };
            }
        }
        Some(LobpcgState {
            n,
            b,
            x: orthonormalize(&x0, n, b),
            x_prev: Vec::new(),
            iters: 0,
        })
    }

    /// Fresh state seeded with a CALLER-SUPPLIED row-major n×b block
    /// (QR-cleaned here) — the warm-start hook for parameter
    /// continuation: seed with the previous parameter point's Ritz
    /// block. Returns `None` for a wrong-size, non-finite block, a
    /// block size violating `1 <= b`, `3b <= n`, or hostile sizes
    /// whose products overflow `usize` (checked arithmetic; refusal,
    /// never a wrap or panic).
    #[must_use]
    pub fn with_block(n: usize, b: usize, x0: &[f64]) -> Option<LobpcgState> {
        let shape = checked_lobpcg_shape(n, b)?;
        if x0.len() != shape.n_b || x0.iter().any(|x| !x.is_finite()) {
            return None;
        }
        // Refusal semantics: a rank-deficient seed must refuse rather
        // than be silently completed with arbitrary orthogonal
        // directions.
        if independent_columns(x0, n, b).len() < b {
            return None;
        }
        Some(LobpcgState {
            n,
            b,
            x: orthonormalize(x0, n, b),
            x_prev: Vec::new(),
            iters: 0,
        })
    }
}

/// Indices of numerically independent columns of a row-major n×k
/// block, detected by two-pass modified Gram–Schmidt on a column-scaled
/// SCRATCH copy (relative tolerance 1e-12). Scaling before the norm and dot
/// products makes the decision invariant to finite power-of-two-sized input
/// ranges instead of overflowing or underflowing the rank test.
fn independent_columns(block: &[f64], n: usize, k: usize) -> Vec<usize> {
    let mut kept: Vec<usize> = Vec::new();
    let mut qcols: Vec<Vec<f64>> = Vec::new();
    for j in 0..k {
        let mut col: Vec<f64> = (0..n).map(|i| block[i * k + j]).collect();
        let scale = max_abs(&col);
        if !(scale.is_finite() && scale > 0.0) {
            continue;
        }
        for x in &mut col {
            *x /= scale;
        }
        let orig = norm2(&col);
        for _ in 0..2 {
            for q in &qcols {
                let c = dot(&col, q);
                for (x, &qi) in col.iter_mut().zip(q) {
                    *x = (-c).mul_add(qi, *x);
                }
            }
        }
        let nrm = norm2(&col);
        if nrm > 0.0 && nrm > 1e-12 * orig {
            for x in &mut col {
                *x /= nrm;
            }
            qcols.push(col);
            kept.push(j);
        }
    }
    kept
}

fn scale_columns_for_qr(block: &mut [f64], n: usize, k: usize) -> Option<()> {
    if block.len() != checked_work_len(n, k)? {
        return None;
    }
    for j in 0..k {
        let mut scale = 0.0f64;
        for i in 0..n {
            scale = scale.max(block[i * k + j].abs());
        }
        if !(scale.is_finite() && scale > 0.0) {
            return None;
        }
        // Preserve the established operation bits for ordinary, QR-safe
        // columns. Rescale only when the factorization's raw sum of squares
        // could overflow or flush every square to zero.
        let safe_high = (f64::MAX / n.max(1) as f64).sqrt() * 0.5;
        let safe_low = f64::MIN_POSITIVE.sqrt();
        if scale > safe_high || scale < safe_low {
            for i in 0..n {
                block[i * k + j] /= scale;
            }
        }
    }
    Some(())
}

/// Orthonormalize the columns of a row-major n×k block via QR (returns
/// Q's first k columns explicitly).
fn orthonormalize(block: &[f64], n: usize, k: usize) -> Vec<f64> {
    let len = checked_work_len(n, k).expect("QR shape was preflighted by the constructor");
    assert_eq!(block.len(), len, "QR block must be n*k");
    let mut scaled = block.to_vec();
    scale_columns_for_qr(&mut scaled, n, k).expect("QR columns must be finite and nonzero");
    let f = qr(&scaled, n, k);
    let mut q = vec![0.0f64; len];
    for j in 0..k {
        let mut e = vec![0.0f64; n];
        e[j] = 1.0;
        f.apply_q(&mut e);
        for i in 0..n {
            q[i * k + j] = e[i];
        }
    }
    q
}

/// Run LOBPCG for `steps` iterations toward the `largest` or smallest
/// end of the spectrum, with an optional preconditioner `prec` (maps a
/// residual block column to a preconditioned direction; identity =
/// `|r, out| out.copy_from_slice(r)`). Returns the current Ritz pairs
/// with numerically recomputed residuals.
pub fn lobpcg_run<Op, Pr>(
    op: &Op,
    state: &mut LobpcgState,
    steps: usize,
    largest: bool,
    prec: &Pr,
) -> Vec<EigenPair>
where
    Op: Fn(&[f64], &mut [f64]),
    Pr: Fn(&[f64], &mut [f64]),
{
    let (n, b) = (state.n, state.b);
    let shape =
        checked_lobpcg_shape(n, b).expect("LOBPCG state shape was preflighted by its constructor");
    for _ in 0..steps {
        // Ritz values for the current X.
        let ax = apply_block(op, &state.x, n, b);
        let lam = block_rayleigh(&state.x, &ax, n, b);
        // Residual block W = AX − X·diag(λ), preconditioned column-wise.
        let mut w = vec![0.0f64; shape.n_b];
        for j in 0..b {
            let mut col = vec![0.0f64; n];
            for i in 0..n {
                col[i] = (-lam[j]).mul_add(state.x[i * b + j], ax[i * b + j]);
            }
            let mut pcol = vec![0.0f64; n];
            prec(&col, &mut pcol);
            for i in 0..n {
                w[i * b + j] = pcol[i];
            }
        }
        // Search space Z = [X | W | P], QR-orthonormalized.
        let p_cols = if state.x_prev.is_empty() { 0 } else { b };
        let (full_zc, z_len) = if p_cols == 0 {
            (shape.two_b, shape.n_two_b)
        } else {
            (shape.three_b, shape.n_three_b)
        };
        let mut z = vec![0.0f64; z_len];
        for i in 0..n {
            for j in 0..b {
                z[i * full_zc + j] = state.x[i * b + j];
                z[i * full_zc + b + j] = w[i * b + j];
            }
            for j in 0..p_cols {
                // P = X − X_prev (the classic conjugate direction).
                z[i * full_zc + shape.two_b + j] = state.x[i * b + j] - state.x_prev[i * b + j];
            }
        }
        // Drop numerically dependent Z columns BEFORE QR: at
        // convergence W ~ 0 and Householder completion would otherwise
        // inject arbitrary directions outside span(Z) into the Ritz
        // space. Scratch-column scaling makes the rank decision safe across
        // finite magnitudes while ordinary QR-safe full-rank paths retain
        // their established operation bits.
        let kept = independent_columns(&z, n, full_zc);
        let zc = kept.len();
        if zc < b {
            // Search space collapsed to the converged block: stop.
            break;
        }
        let z = if zc == full_zc {
            z
        } else {
            let compact_len = checked_work_len(n, zc)
                .expect("rank-truncated LOBPCG shape fits the constructor preflight");
            let mut zk = vec![0.0f64; compact_len];
            for i in 0..n {
                for (slot, &j) in kept.iter().enumerate() {
                    zk[i * zc + slot] = z[i * full_zc + j];
                }
            }
            zk
        };
        let zq = orthonormalize(&z, n, zc);
        // Rayleigh–Ritz on the reduced space.
        let az = apply_block(op, &zq, n, zc);
        let reduced_square =
            checked_work_len(zc, zc).expect("reduced LOBPCG shape fits the constructor preflight");
        let mut s = vec![0.0f64; reduced_square];
        for p in 0..zc {
            for q in 0..zc {
                let mut acc = 0.0f64;
                for i in 0..n {
                    acc = zq[i * zc + p].mul_add(az[i * zc + q], acc);
                }
                s[p * zc + q] = acc;
            }
        }
        // Exact symmetrization (Jacobi assumes it; roundoff breaks it).
        for p in 0..zc {
            for q in p + 1..zc {
                let avg = f64::midpoint(s[p * zc + q], s[q * zc + p]);
                s[p * zc + q] = avg;
                s[q * zc + p] = avg;
            }
        }
        let (evals, evecs) = jacobi_eigh(&s, zc);
        // Take b extremal Ritz vectors.
        let cols: Vec<usize> = if largest {
            (zc - b..zc).rev().collect()
        } else {
            (0..b).collect()
        };
        let mut x_new = vec![0.0f64; shape.n_b];
        for (slot, &c) in cols.iter().enumerate() {
            for i in 0..n {
                let mut acc = 0.0f64;
                for p in 0..zc {
                    acc = zq[i * zc + p].mul_add(evecs[p * zc + c], acc);
                }
                x_new[i * b + slot] = acc;
            }
        }
        let _ = evals;
        state.x_prev = std::mem::take(&mut state.x);
        state.x = orthonormalize(&x_new, n, b);
        state.iters += 1;
    }
    // Final Ritz pairs with recomputed residuals.
    let ax = apply_block(op, &state.x, n, b);
    let lam = block_rayleigh(&state.x, &ax, n, b);
    let mut out = Vec::with_capacity(b);
    for j in 0..b {
        let mut v = vec![0.0f64; n];
        for i in 0..n {
            v[i] = state.x[i * b + j];
        }
        let residual = norm2_values((0..n).map(|i| (-lam[j]).mul_add(v[i], ax[i * b + j])));
        out.push(EigenPair {
            value: lam[j],
            vector: v,
            residual,
        });
    }
    // Deterministic presentation: ascending (or descending for largest).
    out.sort_by(|a, c| {
        if largest {
            c.value.total_cmp(&a.value)
        } else {
            a.value.total_cmp(&c.value)
        }
    });
    out
}

fn apply_block<Op>(op: &Op, x: &[f64], n: usize, k: usize) -> Vec<f64>
where
    Op: Fn(&[f64], &mut [f64]),
{
    let len = checked_work_len(n, k).expect("operator block shape was preflighted");
    assert_eq!(x.len(), len, "operator input block must be n*k");
    let mut out = vec![0.0f64; len];
    let mut col = vec![0.0f64; n];
    let mut res = vec![0.0f64; n];
    for j in 0..k {
        for i in 0..n {
            col[i] = x[i * k + j];
        }
        op(&col, &mut res);
        for i in 0..n {
            out[i * k + j] = res[i];
        }
    }
    out
}

fn block_rayleigh(x: &[f64], ax: &[f64], n: usize, b: usize) -> Vec<f64> {
    let mut lam = vec![0.0f64; b];
    for (j, l) in lam.iter_mut().enumerate() {
        let mut acc = 0.0f64;
        for i in 0..n {
            acc = x[i * b + j].mul_add(ax[i * b + j], acc);
        }
        *l = acc;
    }
    lam
}
