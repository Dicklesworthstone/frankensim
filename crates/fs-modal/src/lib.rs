//! fs-modal — the first-class vibration eigenproblem facility (bead
//! frankensim-fsim-vibration-eig-jw6yq, musical-acoustics program):
//! K·φ = λ·M·φ for symmetric K and symmetric positive-definite M, with
//! residual-derived CERTIFIED INTERVALS per mode and INERTIA-CERTIFIED
//! spectrum slicing ("no missed modes" as evidence, not hope).
//!
//! WHY — everything modal (bars, plates, soundboards, bells, resonators)
//! needs eigenpairs whose completeness can be trusted: a missed body mode is
//! an audible hole. The two load-bearing certificates:
//!
//! - **Count certificate**: for a window (a, b], the number of eigenvalues
//!   inside is EXACTLY `neg(K − b·M) − neg(K − a·M)` by Sylvester's law,
//!   where `neg` is the negative-inertia count of the sparse LDLᵀ
//!   factorization (`fs_sparse::direct`). [`slice_window`] refuses
//!   (`WindowUnresolved`) unless it converges exactly that many modes.
//! - **Value certificate**: with `φᵀMφ = 1`, the residual `r = Kφ − λ̂Mφ`
//!   bounds the distance to a TRUE eigenvalue: `min_i |λ_i − λ̂| ≤
//!   ‖r‖_{M⁻¹}` (standard bound for the symmetric-definite pencil; proof:
//!   substitute z = M^{1/2}φ to reduce to the symmetric standard problem).
//!   Every [`ModePair`] carries that bound as `interval`.
//!
//! Strategy selection is explicit, not magical: [`eigh_gen_dense`] for
//! small/dense pencils (Cholesky reduction + cyclic Jacobi, both fs-la);
//! [`slice_window`] for sparse pencils (shift-invert Lanczos in the M-inner
//! product over `fs_sparse::direct` factorizations, with full
//! reorthogonalization and RESTART-WITH-DEFLATION so degenerate clusters —
//! which single-vector Lanczos provably cannot resolve in one Krylov space —
//! are still recovered and certified); [`shift_invert_modes_mfree`] as the
//! matrix-free core for callers that bring their own shifted inverse.
//! Mode shapes are mass-normalized (`φᵀMφ = 1`) everywhere.
//!
//! Determinism: deterministic start vectors (integer LCG, no platform libm),
//! fixed reorthogonalization order, index-ordered tie-breaks; repeat runs are
//! bitwise identical (tested).
//!
//! The [`time`] module is the time-domain realization of mass-normalized damped modes (exact zero-order-hold stepping and
//! observer pressure); it moved here from `fs-couple` (IRP M1.1-W1a) so that strings, solids and couplers share one acyclic owner.

use fs_la::eigen::jacobi_eigh;
use fs_la::eigen_complex::eig as complex_eig;
use fs_la::factor::cholesky;
use fs_math::c64::C64;
use fs_sparse::{Coo, Csr, DirectOrdering, LdltError, LdltFactor, LdltOptions, SymbolicLdlt};

/// Exact zero-order-hold realization of mass-normalized damped modes (moved from fs-couple, M1.1-W1a).
pub mod time;

/// Physical displacement/momentum transfer between source modal bases.
pub mod transfer;

mod windows;
pub use windows::{MultiSliceReport, SliceTimingReceipt, WindowPlan, slice_windows};

mod damping;
pub use damping::{
    CertifiedModalDamping, FarCouplingBound, MAX_DAMPING_CLUSTER, ModalDampingBlock,
    ModalDampingBudget, ROUCHE_SAFETY, certify_modal_damping, certify_modal_damping_windowed,
    ResidualFarBounds, certify_modal_damping_owned_diagonal,
    certify_modal_damping_windowed_far_rows, damping_window_neighbours, residual_far_row_bounds,
};

/// Crate version, re-exported for provenance stamping.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
/// Arithmetic contract of the certified serial slice, also used by windows.
/// Existing AMD cross-ISA digest is retained unchanged at version1.
pub const MODAL_SLICE_BIT_SEMANTICS_VERSION: u32 = 1;

/// Typed refusals. Display strings carry stable `FS-MODAL-*` codes.
#[derive(Debug, Clone, PartialEq)]
pub enum ModalError {
    /// K and M shapes disagree or are not square.
    DimensionMismatch {
        /// (rows, cols) of K.
        k_shape: (usize, usize),
        /// (rows, cols) of M.
        m_shape: (usize, usize),
    },
    /// The mass matrix is not positive definite (its LDLᵀ inertia reported
    /// negative eigenvalues, or the dense Cholesky failed). The pencil
    /// theory this crate certifies requires SPD M; a semidefinite or
    /// indefinite mass needs a different formulation and is refused.
    MassNotSpd {
        /// Number of negative mass eigenvalues (0 when reported by the
        /// dense path, which stops at the first failing pivot).
        negative: usize,
    },
    /// A requested mass-diagonal coordinate change cannot be represented
    /// with finite positive factors or finite transformed matrix entries.
    Equilibration {
        /// First mass row or transformed matrix row that cannot be represented.
        index: usize,
    },
    /// A sparse factorization stage failed.
    Factor {
        /// The shift at which (K − σM) was being factored (NaN for the mass
        /// factorization).
        shift: f64,
        /// Underlying direct-solver refusal.
        source: LdltError,
    },
    /// The window is malformed (low ≥ high, or non-finite endpoints).
    InvalidWindow {
        /// Lower endpoint.
        low: f64,
        /// Upper endpoint.
        high: f64,
    },
    /// The iteration budget ended before the inertia-certified count of
    /// modes converged: the slice is INCOMPLETE and no partial claim is
    /// made. This is also the mutation gate: a harvester that skipped a
    /// cluster cannot return success.
    WindowUnresolved {
        /// Inertia-certified number of eigenvalues in the window.
        expected: usize,
        /// Modes actually converged within budget.
        converged: usize,
        /// The window.
        window: (f64, f64),
    },
    /// The dense complex eigensolver failed on the linearized quadratic
    /// pencil.
    QuadraticEig {
        /// Diagnostic from fs-la.
        detail: String,
    },
    /// No partial multi-window report is returned after cancellation/drain.
    Cancelled,
    /// Invalid window plan or resource limit.
    WindowPlan {
        field: &'static str,
    },
    WindowBudget {
        requested: usize,
        cap: usize,
    },
    MemoryBudget {
        estimated_bytes: usize,
        cap: usize,
    },
    Executor {
        detail: String,
    },
    /// Certified modal damping input refused: budget, frequency, shape,
    /// non-finite, asymmetric or negative entry, window, or column bound.
    DampingInput {
        /// What was refused.
        field: &'static str,
        /// Mode or row-major entry index.
        index: usize,
    },
    /// A diagonal estimate or a block has real poles; the certificate is
    /// for lightly damped (underdamped) modes only.
    DampingNotUnderdamped {
        /// Offending modes.
        modes: Vec<usize>,
    },
    /// The smallest certified cluster exceeds the budget's cap.
    DampingClusterTooLarge {
        /// The cluster.
        modes: Vec<usize>,
        /// The cap.
        cap: usize,
    },
    /// A required block needs coupling entries outside the supplied window.
    DampingBlockSpansWindow {
        /// The block.
        modes: Vec<usize>,
        /// The window [rad/s].
        window: f64,
    },
}

impl core::fmt::Display for ModalError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            ModalError::DimensionMismatch { k_shape, m_shape } => write!(
                f,
                "FS-MODAL-DIMENSION-MISMATCH: K is {}x{}, M is {}x{}",
                k_shape.0, k_shape.1, m_shape.0, m_shape.1
            ),
            ModalError::MassNotSpd { negative } => write!(
                f,
                "FS-MODAL-MASS-NOT-SPD: mass matrix has {negative} negative eigenvalue(s)"
            ),
            ModalError::Equilibration { index } => write!(
                f,
                "FS-MODAL-EQUILIBRATION: nonpositive or unrepresentable mass-coordinate scale at row {index}"
            ),
            ModalError::Factor { shift, source } => {
                write!(f, "FS-MODAL-FACTOR at shift {shift}: {source}")
            }
            ModalError::InvalidWindow { low, high } => {
                write!(f, "FS-MODAL-INVALID-WINDOW: ({low}, {high}]")
            }
            ModalError::WindowUnresolved {
                expected,
                converged,
                window,
            } => write!(
                f,
                "FS-MODAL-WINDOW-UNRESOLVED: inertia certifies {expected} mode(s) in \
                 ({}, {}] but only {converged} converged within budget",
                window.0, window.1
            ),
            ModalError::QuadraticEig { detail } => {
                write!(f, "FS-MODAL-QUADRATIC-EIG: {detail}")
            }
            ModalError::Cancelled => write!(
                f,
                "FS-MODAL-CANCELLED: all windows drained; no partial report"
            ),
            ModalError::WindowPlan { field } => write!(f, "FS-MODAL-WINDOW-PLAN: {field}"),
            ModalError::WindowBudget { requested, cap } => write!(
                f,
                "FS-MODAL-WINDOW-BUDGET: {requested} windows exceeds cap {cap}"
            ),
            ModalError::MemoryBudget {
                estimated_bytes,
                cap,
            } => write!(
                f,
                "FS-MODAL-MEMORY-BUDGET: payload {estimated_bytes} exceeds cap {cap}"
            ),
            ModalError::Executor { detail } => write!(f, "FS-MODAL-EXECUTOR: {detail}"),
            ModalError::DampingInput { field, index } => {
                write!(f, "FS-MODAL-DAMPING-INPUT: {field} (index {index})")
            }
            ModalError::DampingNotUnderdamped { modes } => write!(
                f,
                "FS-MODAL-DAMPING-NOT-UNDERDAMPED: modes {modes:?} have real poles"
            ),
            ModalError::DampingClusterTooLarge { modes, cap } => write!(
                f,
                "FS-MODAL-DAMPING-CLUSTER: certified cluster of {} modes {modes:?} exceeds cap {cap}",
                modes.len()
            ),
            ModalError::DampingBlockSpansWindow { modes, window } => write!(
                f,
                "FS-MODAL-DAMPING-WINDOW: block {modes:?} needs coupling beyond window {window} rad/s"
            ),
        }
    }
}

impl std::error::Error for ModalError {}

/// One certified eigenpair of the pencil (K, M).
#[derive(Debug, Clone)]
pub struct ModePair {
    /// Ritz value λ̂.
    pub lambda: f64,
    /// Mass-normalized mode shape: `φᵀMφ = 1` (enforced numerically).
    pub phi: Vec<f64>,
    /// `‖Kφ − λ̂Mφ‖_{M⁻¹}` — the certified distance bound to a true
    /// eigenvalue of the pencil.
    pub residual: f64,
    /// `[λ̂ − residual, λ̂ + residual]`: contains at least one TRUE
    /// eigenvalue of the pencil (SPD-M theory; see the crate docs).
    pub interval: (f64, f64),
}

/// Work accounting for one [`slice_window`] call.
#[derive(Debug, Clone, Copy)]
pub struct SliceStats {
    /// The shift-invert pole actually used.
    pub shift: f64,
    /// Number of (K − σM) factorization attempts (window endpoints and any
    /// re-picked interior shifts included).
    pub factorizations: usize,
    /// Total Lanczos steps across restarts.
    pub lanczos_iters: usize,
    /// Deflation restarts taken (each recovers eigenvector directions a
    /// single Krylov space cannot represent, e.g. degenerate clusters).
    pub restarts: usize,
    /// nnz(L) of the interior-shift factorization.
    pub factor_nnz_l: usize,
    /// Peak frontal bytes of the interior-shift factorization.
    pub factor_peak_bytes: usize,
    /// Delayed pivots in the interior-shift factorization.
    pub pivots_delayed: usize,
}

/// Result of an inertia-certified window slice.
#[derive(Debug, Clone)]
pub struct SliceReport {
    /// The window (low, high] that was sliced.
    pub window: (f64, f64),
    /// Eigenvalue count strictly below the LOW endpoint (Sylvester inertia
    /// of K − low·M).
    pub below_low: usize,
    /// Eigenvalue count strictly below the HIGH endpoint.
    pub below_high: usize,
    /// `below_high − below_low`: the certified in-window count.
    pub expected: usize,
    /// The converged modes, ascending in λ. `modes.len() == expected` is an
    /// invariant of a returned report (otherwise [`ModalError::WindowUnresolved`]).
    pub modes: Vec<ModePair>,
    /// Work accounting.
    pub stats: SliceStats,
}

/// Options for [`slice_window`].
#[derive(Debug, Clone, Copy)]
pub struct SliceOptions {
    /// Total Lanczos-step budget across restarts. `0` selects
    /// `min(n, 20 + 6·expected)`.
    pub max_lanczos: usize,
    /// Maximum deflation restarts. `0` selects `2 + expected`.
    pub max_restarts: usize,
    /// Relative convergence tolerance on the shift-invert Ritz residual
    /// estimate before a mode is accepted for explicit certification.
    pub ritz_tol: f64,
    /// Options for the sparse LDLᵀ factorizations.
    pub ldlt: LdltOptions,
    /// Fill-reducing ordering.
    pub ordering: DirectOrdering,
    /// Solve in the congruent coordinates x = D y, D_ii = 1/sqrt(M_ii).
    /// This changes conditioning but neither the eigenvalues nor the inertia
    /// count; certificates are recomputed in the original coordinates.
    pub mass_diagonal_equilibration: bool,
}

impl Default for SliceOptions {
    fn default() -> SliceOptions {
        SliceOptions {
            max_lanczos: 0,
            max_restarts: 0,
            ritz_tol: 1e-10,
            ldlt: LdltOptions::default(),
            ordering: DirectOrdering::Amd,
            mass_diagonal_equilibration: false,
        }
    }
}

// ---------------------------------------------------------------------------
// Pencil assembly helpers
// ---------------------------------------------------------------------------

fn check_shapes(k: &Csr, m: &Csr) -> Result<usize, ModalError> {
    let n = k.nrows();
    if k.ncols() != n || m.nrows() != n || m.ncols() != n {
        return Err(ModalError::DimensionMismatch {
            k_shape: (k.nrows(), k.ncols()),
            m_shape: (m.nrows(), m.ncols()),
        });
    }
    Ok(n)
}

/// K − σ·M with the UNION pattern of K and M regardless of σ (stored zeros
/// are kept), so one symbolic analysis serves every shift.
fn shifted_pencil(k: &Csr, m: &Csr, sigma: f64) -> Csr {
    let n = k.nrows();
    let mut coo = Coo::new(n, n);
    for r in 0..n {
        let (cols, vals) = k.row(r);
        for (&c, &v) in cols.iter().zip(vals) {
            coo.push(r, c, v);
        }
        let (mc, mv) = m.row(r);
        for (&c, &v) in mc.iter().zip(mv) {
            coo.push(r, c, -sigma * v);
        }
    }
    coo.assemble()
}

fn mass_coordinate_scales(m: &Csr) -> Result<Vec<f64>, ModalError> {
    let mut scales = Vec::with_capacity(m.nrows());
    for r in 0..m.nrows() {
        let (cols, vals) = m.row(r);
        let diagonal = cols.binary_search(&r).ok().map(|i| vals[i]);
        let Some(diagonal) = diagonal.filter(|v| v.is_finite() && *v > 0.0) else {
            return Err(ModalError::Equilibration { index: r });
        };
        let scale = diagonal.sqrt().recip();
        if !scale.is_finite() || scale <= 0.0 {
            return Err(ModalError::Equilibration { index: r });
        }
        scales.push(scale);
    }
    Ok(scales)
}

fn congruence_scaled(a: &Csr, scales: &[f64]) -> Result<Csr, ModalError> {
    let n = a.nrows();
    let mut coo = Coo::new(n, n);
    for r in 0..n {
        let (cols, vals) = a.row(r);
        for (&c, &v) in cols.iter().zip(vals) {
            let value = (v * scales[r]) * scales[c];
            if !value.is_finite() || (v != 0.0 && value == 0.0) {
                return Err(ModalError::Equilibration { index: r });
            }
            coo.push(r, c, value);
        }
    }
    Ok(coo.assemble())
}

/// Deterministic pseudo-random fill for start vectors: integer LCG only —
/// no platform libm feeds solver state (workspace doctrine).
fn lcg_fill(v: &mut [f64], seed: u64) {
    let mut s = seed.wrapping_mul(0x9E37_79B9_7F4A_7C15).wrapping_add(1);
    for slot in v.iter_mut() {
        s = s
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        *slot = ((s >> 11) as f64) / (1u64 << 53) as f64 - 0.5;
    }
}

// ---------------------------------------------------------------------------
// M-inner-product shift-invert Lanczos core
// ---------------------------------------------------------------------------

/// Internal result of one deflated Lanczos sweep.
struct SweepOutcome {
    /// (θ shift-invert Ritz value, Ritz vector in problem space).
    accepted: Vec<(f64, Vec<f64>)>,
    iters: usize,
}

/// One Lanczos sweep in the M-inner product on `op = solve∘M`, with full
/// reorthogonalization against both the growing basis and the already
/// `converged` (deflated) directions. Accepts Ritz pairs whose shift-invert
/// residual estimate `|β_j · y_last|` is below `tol·|θ|` and whose mapped
/// eigenvalue λ = σ + 1/θ lies in `window` (when given).
#[allow(clippy::too_many_arguments)] // internal core shared by two public fronts
#[allow(clippy::too_many_lines)] // one coherent Krylov sweep
fn lanczos_sweep(
    n: usize,
    sigma: f64,
    window: Option<(f64, f64)>,
    budget: usize,
    need: usize,
    tol: f64,
    seed: u64,
    apply_m: &mut dyn FnMut(&[f64], &mut [f64]),
    solve_shifted: &mut dyn FnMut(&[f64], &mut [f64]),
    converged: &[(f64, Vec<f64>)],
) -> SweepOutcome {
    let mut basis: Vec<Vec<f64>> = Vec::new(); // M-orthonormal v_j
    let mut mbasis: Vec<Vec<f64>> = Vec::new(); // M·v_j (cached)
    let mut alphas: Vec<f64> = Vec::new();
    let mut betas: Vec<f64> = Vec::new(); // β_j links v_j → v_{j+1}
    let mut scratch_m = vec![0.0f64; n];
    let mut operator_scale = 0.0_f64;

    // Converged directions in M-image form for deflation.
    let deflate_m: Vec<Vec<f64>> = converged
        .iter()
        .map(|(_, phi)| {
            let mut mp = vec![0.0f64; n];
            apply_m(phi, &mut mp);
            mp
        })
        .collect();

    // Deterministic start vector, deflated and M-normalized.
    let mut v = vec![0.0f64; n];
    lcg_fill(&mut v, seed);
    for (idx, (_, phi)) in converged.iter().enumerate() {
        let c: f64 = v.iter().zip(&deflate_m[idx]).map(|(a, b)| a * b).sum();
        for (slot, p) in v.iter_mut().zip(phi) {
            *slot -= c * p;
        }
    }
    apply_m(&v, &mut scratch_m);
    let nrm: f64 = v
        .iter()
        .zip(&scratch_m)
        .map(|(a, b)| a * b)
        .sum::<f64>()
        .max(0.0)
        .sqrt();
    if nrm <= 0.0 {
        return SweepOutcome {
            accepted: Vec::new(),
            iters: 0,
        };
    }
    for slot in &mut v {
        *slot /= nrm;
    }

    let mut iters = 0usize;
    let mut accepted: Vec<(f64, Vec<f64>)> = Vec::new();
    while iters < budget && basis.len() < n {
        // Cache M·v_j.
        let mut mv = vec![0.0f64; n];
        apply_m(&v, &mut mv);
        basis.push(v.clone());
        mbasis.push(mv);
        let j = basis.len() - 1;

        // w = op v_j = (K − σM)⁻¹ M v_j.
        let mut w = vec![0.0f64; n];
        solve_shifted(&mbasis[j], &mut w);

        // Three-term coefficients (computed via explicit M-inner products;
        // full reorthogonalization below subsumes the recurrence anyway).
        let alpha: f64 = w.iter().zip(&mbasis[j]).map(|(a, b)| a * b).sum();
        alphas.push(alpha);
        // Full reorthogonalization against the basis (two passes, classic
        // Kahan twice-is-enough) and the deflated converged directions.
        for _pass in 0..2 {
            for (bi, mb) in basis.iter().zip(&mbasis) {
                let c: f64 = w.iter().zip(mb).map(|(a, b)| a * b).sum();
                for (slot, b) in w.iter_mut().zip(bi) {
                    *slot -= c * b;
                }
            }
            for (idx, (_, phi)) in converged.iter().enumerate() {
                let c: f64 = w.iter().zip(&deflate_m[idx]).map(|(a, b)| a * b).sum();
                for (slot, p) in w.iter_mut().zip(phi) {
                    *slot -= c * p;
                }
            }
        }
        apply_m(&w, &mut scratch_m);
        let beta = w
            .iter()
            .zip(&scratch_m)
            .map(|(a, b)| a * b)
            .sum::<f64>()
            .max(0.0)
            .sqrt();
        iters += 1;
        // The shifted inverse has units of inverse eigenvalue. An absolute
        // beta cutoff (formerly 1e-14) mistakes tightly clustered inverse
        // values in wide acoustic bands for an invariant Krylov subspace.
        // Scale to the observed tridiagonal coefficients, not max(scale, 1):
        // that floor would reintroduce dependence on the physical units.
        operator_scale = operator_scale.max(alpha.abs()).max(beta);
        let breakdown = beta == 0.0
            || (operator_scale.is_finite() && beta <= (16.0 * f64::EPSILON) * operator_scale);

        // Ritz extraction on the tridiagonal T_j (small dense; j ≤ budget).
        let jn = basis.len();
        let mut t = vec![0.0f64; jn * jn];
        for (i, &a) in alphas.iter().enumerate() {
            t[i * jn + i] = a;
        }
        for i in 0..jn.saturating_sub(1) {
            t[i * jn + (i + 1)] = betas[i];
            t[(i + 1) * jn + i] = betas[i];
        }
        let (theta, y) = jacobi_eigh(&t, jn);

        // Which Ritz pairs are converged in the shift-invert metric?
        let mut ready: Vec<(f64, usize)> = Vec::new();
        for (ri, &th) in theta.iter().enumerate() {
            if th == 0.0 {
                continue;
            }
            let lambda = sigma + 1.0 / th;
            if let Some((lo, hi)) = window
                && !(lambda > lo && lambda <= hi)
            {
                continue;
            }
            let est = beta * y[(jn - 1) * jn + ri].abs();
            if est <= tol * th.abs().max(f64::MIN_POSITIVE) {
                ready.push((th, ri));
            }
        }
        // Stop on breakdown, exhaustion, or as soon as this sweep has
        // produced everything the caller still needs — leaving the rest of
        // the total budget for deflation restarts (a single Krylov space
        // provably yields at most one Ritz direction per distinct
        // eigenvalue, so degenerate clusters REQUIRE later sweeps).
        let stop = breakdown || basis.len() == n || iters == budget || ready.len() >= need;
        if stop {
            // Harvest all ready Ritz pairs into problem space.
            for &(_th, ri) in &ready {
                let mut phi = vec![0.0f64; n];
                for (bi, b) in basis.iter().enumerate() {
                    let c = y[bi * jn + ri];
                    for (slot, bv) in phi.iter_mut().zip(b) {
                        *slot = c.mul_add(*bv, *slot);
                    }
                }
                let th = theta[ri];
                accepted.push((sigma + 1.0 / th, phi));
            }
            break;
        }
        // Advance.
        betas.push(beta);
        for (slot, wv) in v.iter_mut().zip(&w) {
            *slot = *wv / beta;
        }
    }
    SweepOutcome { accepted, iters }
}

// ---------------------------------------------------------------------------
// Public sparse path: inertia-certified spectrum slicing
// ---------------------------------------------------------------------------

/// Slice the window `(low, high]` of the pencil (K, M): certify the exact
/// in-window eigenvalue count by Sylvester inertia at both endpoints, then
/// harvest exactly that many eigenpairs by shift-invert Lanczos in the
/// M-inner product (full reorthogonalization + deflation restarts), each
/// with a certified interval. Returns [`ModalError::WindowUnresolved`] —
/// never a silently short list — if the budget ends first.
///
/// Endpoints that graze an eigenvalue make (K − σM) singular; the
/// factorization then refuses ([`ModalError::Factor`]) and the caller picks
/// a nearby endpoint — refusal beats a fabricated count.
///
/// # Errors
/// See [`ModalError`]. Non-finite matrix entries surface as
/// [`ModalError::Factor`] wrapping the direct solver's named refusal.
#[allow(clippy::too_many_lines)] // one coherent certified pipeline
pub fn slice_window(
    k: &Csr,
    m: &Csr,
    window: (f64, f64),
    opts: &SliceOptions,
) -> Result<SliceReport, ModalError> {
    let (lo, hi) = window;
    if !(lo.is_finite() && hi.is_finite() && lo < hi) {
        return Err(ModalError::InvalidWindow { low: lo, high: hi });
    }
    let prepared = PreparedPencil::new(k, m, lo, opts, None, &mut || Ok(()))?;
    prepared.slice(window, opts, None, &mut || Ok(()))
}

/// Single owner for mass admission, scaling and the shared shifted pattern.
/// Single-window and parallel-window calls use the identical numeric body.
struct PreparedPencil<'a> {
    k: &'a Csr,
    m: &'a Csr,
    m_fact: LdltFactor,
    scales: Option<Vec<f64>>,
    scaled: Option<(Csr, Csr)>,
    p_sym: SymbolicLdlt,
}
impl<'a> PreparedPencil<'a> {
    fn new(
        k: &'a Csr,
        m: &'a Csr,
        lo: f64,
        opts: &SliceOptions,
        coordinates: Option<&[[f64; 3]]>,
        checkpoint: &mut dyn FnMut() -> Result<(), ModalError>,
    ) -> Result<Self, ModalError> {
        check_shapes(k, m)?;
        checkpoint()?;

        // Mass factorization: SPD gate + the M⁻¹ norm used by every certificate.
        let m_sym = SymbolicLdlt::analyze_with_coordinates(m, opts.ordering, coordinates).map_err(
            |source| ModalError::Factor {
                shift: f64::NAN,
                source,
            },
        )?;
        checkpoint()?;
        let m_fact = m_sym
            .factor(m, &opts.ldlt)
            .map_err(|source| ModalError::Factor {
                shift: f64::NAN,
                source,
            })?;
        if m_fact.inertia().negative > 0 {
            return Err(ModalError::MassNotSpd {
                negative: m_fact.inertia().negative,
            });
        }

        let scales = if opts.mass_diagonal_equilibration {
            Some(mass_coordinate_scales(m)?)
        } else {
            None
        };
        let scaled = scales
            .as_ref()
            .map(|scale| Ok((congruence_scaled(k, scale)?, congruence_scaled(m, scale)?)))
            .transpose()?;
        let (work_k, work_m) = scaled.as_ref().map_or((k, m), |(ks, ms)| (ks, ms));

        // One symbolic analysis serves every shift (union pattern).
        let a_lo = shifted_pencil(work_k, work_m, lo);
        checkpoint()?;
        let p_sym = SymbolicLdlt::analyze_with_coordinates(&a_lo, opts.ordering, coordinates)
            .map_err(|source| ModalError::Factor { shift: lo, source })?;
        Ok(Self {
            k,
            m,
            m_fact,
            scales,
            scaled,
            p_sym,
        })
    }

    fn work(&self) -> (&Csr, &Csr) {
        self.scaled
            .as_ref()
            .map_or((self.k, self.m), |(k, m)| (k, m))
    }

    fn slice(
        &self,
        window: (f64, f64),
        opts: &SliceOptions,
        counts: Option<(usize, usize)>,
        checkpoint: &mut dyn FnMut() -> Result<(), ModalError>,
    ) -> Result<SliceReport, ModalError> {
        let (k, m) = (self.k, self.m);
        let n = k.nrows();
        let (lo, hi) = window;
        let (work_k, work_m) = self.work();
        let (p_sym, m_fact, scales) = (&self.p_sym, &self.m_fact, &self.scales);
        let mut factorizations = 0usize;
        let mut factor_at =
            |sigma: f64, mat: &Csr, count: &mut usize| -> Result<LdltFactor, ModalError> {
                checkpoint()?;
                *count += 1;
                p_sym
                    .factor(mat, &opts.ldlt)
                    .map_err(|source| ModalError::Factor {
                        shift: sigma,
                        source,
                    })
            };

        let (below_low, below_high) = if let Some(counts) = counts {
            counts
        } else {
            let a_lo = shifted_pencil(work_k, work_m, lo);
            let f_lo = factor_at(lo, &a_lo, &mut factorizations)?;
            let a_hi = shifted_pencil(work_k, work_m, hi);
            let f_hi = factor_at(hi, &a_hi, &mut factorizations)?;
            let below_low = f_lo.inertia().negative;
            let below_high = f_hi.inertia().negative;
            drop(f_lo);
            drop(f_hi);
            (below_low, below_high)
        };
        let expected = below_high - below_low;

        if expected == 0 {
            return Ok(SliceReport {
                window,
                below_low,
                below_high,
                expected: 0,
                modes: Vec::new(),
                stats: SliceStats {
                    shift: f64::NAN,
                    factorizations,
                    lanczos_iters: 0,
                    restarts: 0,
                    factor_nnz_l: 0,
                    factor_peak_bytes: 0,
                    pivots_delayed: 0,
                },
            });
        }

        // Interior shift: midpoint first, deterministic golden-ratio fallbacks
        // if the factorization refuses (shift grazing an eigenvalue).
        let fractions: [f64; 4] = [0.5, 0.381_966_011_250_105, 0.618_033_988_749_895, 0.25];
        let mut chosen: Option<(f64, LdltFactor)> = None;
        let mut last_err = None;
        for &fr in &fractions {
            let sigma = fr.mul_add(hi - lo, lo);
            let a_mid = shifted_pencil(work_k, work_m, sigma);
            match factor_at(sigma, &a_mid, &mut factorizations) {
                Ok(f) => {
                    chosen = Some((sigma, f));
                    break;
                }
                Err(
                    e @ ModalError::Factor {
                        source: LdltError::PivotBreakdown { .. },
                        ..
                    },
                ) => last_err = Some(e),
                Err(e) => return Err(e),
            }
        }
        let Some((sigma, f_mid)) = chosen else {
            return Err(last_err.expect("at least one shift attempt was made"));
        };
        let factor_nnz_l = f_mid.stats().nnz_l;
        let factor_peak_bytes = f_mid.stats().peak_front_bytes;
        let pivots_delayed = f_mid.stats().pivots_delayed;

        let max_restarts = if opts.max_restarts == 0 {
            2 + expected
        } else {
            opts.max_restarts
        };
        // Total budget spans all restarts; each sweep gets its own slice sized
        // to what it still needs, so an early sweep cannot starve the deflation
        // restarts that degenerate clusters require.
        let budget_total = if opts.max_lanczos == 0 {
            n.min((20 + 6 * expected) * (1 + max_restarts))
        } else {
            opts.max_lanczos
        };

        let mut apply_m = |x: &[f64], y: &mut [f64]| work_m.spmv(x, y);
        let mut solve_shifted = |b: &[f64], out: &mut [f64]| {
            let x = f_mid.solve(b);
            out.copy_from_slice(&x);
        };

        // Restart-with-deflation until the certified count is met.
        let mut converged: Vec<(f64, Vec<f64>)> = Vec::new();
        let mut iters_used = 0usize;
        let mut restarts = 0usize;
        while converged.len() < expected {
            checkpoint()?;
            if restarts > max_restarts || iters_used >= budget_total {
                return Err(ModalError::WindowUnresolved {
                    expected,
                    converged: converged.len(),
                    window,
                });
            }
            let need = expected - converged.len();
            let sweep = lanczos_sweep(
                n,
                sigma,
                Some(window),
                (20 + 6 * need).min(budget_total - iters_used),
                need,
                opts.ritz_tol,
                0xF5_0DA1 + restarts as u64,
                &mut apply_m,
                &mut solve_shifted,
                &converged,
            );
            iters_used += sweep.iters;
            let found_any = !sweep.accepted.is_empty();
            for (lambda, phi) in sweep.accepted {
                if converged.len() < expected {
                    converged.push((lambda, phi));
                }
            }
            restarts += 1;
            if !found_any && iters_used >= budget_total {
                return Err(ModalError::WindowUnresolved {
                    expected,
                    converged: converged.len(),
                    window,
                });
            }
        }

        // Explicit certification: mass-normalize, compute the true pencil
        // residual, and derive the M⁻¹-norm interval per mode.
        let mut modes: Vec<ModePair> = Vec::with_capacity(expected);
        let mut scratch = vec![0.0f64; n];
        for (lambda, mut phi) in converged {
            checkpoint()?;
            if let Some(scales) = scales {
                for (coordinate, scale) in phi.iter_mut().zip(scales) {
                    *coordinate *= scale;
                }
            }
            m.spmv(&phi, &mut scratch);
            let mnorm: f64 = phi
                .iter()
                .zip(&scratch)
                .map(|(a, b)| a * b)
                .sum::<f64>()
                .max(f64::MIN_POSITIVE)
                .sqrt();
            for slot in &mut phi {
                *slot /= mnorm;
            }
            let mut kphi = vec![0.0f64; n];
            k.spmv(&phi, &mut kphi);
            m.spmv(&phi, &mut scratch);
            let r: Vec<f64> = kphi
                .iter()
                .zip(&scratch)
                .map(|(kk, mm)| lambda.mul_add(-mm, *kk))
                .collect();
            let minv_r = m_fact.solve(&r);
            let bound = r
                .iter()
                .zip(&minv_r)
                .map(|(a, b)| a * b)
                .sum::<f64>()
                .max(0.0)
                .sqrt();
            modes.push(ModePair {
                lambda,
                phi,
                residual: bound,
                interval: (lambda - bound, lambda + bound),
            });
        }
        modes.sort_by(|a, b| a.lambda.total_cmp(&b.lambda));

        Ok(SliceReport {
            window,
            below_low,
            below_high,
            expected,
            modes,
            stats: SliceStats {
                shift: sigma,
                factorizations,
                lanczos_iters: iters_used,
                restarts: restarts.saturating_sub(1),
                factor_nnz_l,
                factor_peak_bytes,
                pivots_delayed,
            },
        })
    }
}

// ---------------------------------------------------------------------------
// Matrix-free core
// ---------------------------------------------------------------------------

/// Matrix-free shift-invert harvest: the caller supplies `apply_m` (y = M·x)
/// and `solve_shifted` (y = (K − σM)⁻¹·x) closures; this runs the same
/// M-inner-product Lanczos core as [`slice_window`] and returns up to `want`
/// eigenpairs nearest σ (mass-normalized), ascending in λ.
///
/// NO-CLAIM: without an M-solve there is no M⁻¹-norm certificate, so this
/// path returns raw Ritz pairs WITHOUT certified intervals and performs no
/// inertia count certification — the assembled [`slice_window`] path is the
/// certified front door.
///
/// # Errors
/// [`ModalError::WindowUnresolved`] (with `expected = want`) if the budget
/// ends before `want` pairs converge.
pub fn shift_invert_modes_mfree(
    n: usize,
    sigma: f64,
    want: usize,
    budget: usize,
    mut apply_m: impl FnMut(&[f64], &mut [f64]),
    mut solve_shifted: impl FnMut(&[f64], &mut [f64]),
) -> Result<Vec<(f64, Vec<f64>)>, ModalError> {
    let mut converged: Vec<(f64, Vec<f64>)> = Vec::new();
    let mut iters = 0usize;
    let mut restarts = 0usize;
    while converged.len() < want {
        if iters >= budget || restarts > want + 2 {
            return Err(ModalError::WindowUnresolved {
                expected: want,
                converged: converged.len(),
                window: (f64::NEG_INFINITY, f64::INFINITY),
            });
        }
        let need = want - converged.len();
        let sweep = lanczos_sweep(
            n,
            sigma,
            None,
            (20 + 6 * need).min(budget - iters),
            need,
            1e-10,
            0xF5_0DA1 + restarts as u64,
            &mut apply_m,
            &mut solve_shifted,
            &converged,
        );
        iters += sweep.iters;
        // Keep pairs nearest the pole first (largest |θ| ⇔ nearest σ).
        let mut got = sweep.accepted;
        got.sort_by(|a, b| {
            let da = (a.0 - sigma).abs();
            let db = (b.0 - sigma).abs();
            da.total_cmp(&db)
        });
        let found_any = !got.is_empty();
        for pair in got {
            if converged.len() < want {
                converged.push(pair);
            }
        }
        restarts += 1;
        if !found_any && iters >= budget {
            return Err(ModalError::WindowUnresolved {
                expected: want,
                converged: converged.len(),
                window: (f64::NEG_INFINITY, f64::INFINITY),
            });
        }
    }
    converged.sort_by(|a, b| a.0.total_cmp(&b.0));
    Ok(converged)
}

// ---------------------------------------------------------------------------
// Dense path
// ---------------------------------------------------------------------------

/// Dense generalized symmetric-definite eigensolve: Cholesky-reduce
/// (M = LLᵀ, C = L⁻¹KL⁻ᵀ), cyclic-Jacobi the reduced matrix, back-transform,
/// mass-normalize, and attach the same M⁻¹-norm certified interval as the
/// sparse path. Inputs are row-major n×n; K must be numerically symmetric
/// and M positive definite. Returns ALL n modes ascending.
///
/// # Errors
/// [`ModalError::MassNotSpd`] if the Cholesky factorization of M fails.
///
/// # Panics
/// On shape mismatches between the slices and `n` (programmer error).
#[allow(clippy::too_many_lines)] // one coherent dense pipeline
pub fn eigh_gen_dense(k: &[f64], m: &[f64], n: usize) -> Result<Vec<ModePair>, ModalError> {
    assert_eq!(k.len(), n * n, "k must be n*n");
    assert_eq!(m.len(), n * n, "m must be n*n");
    if n == 0 {
        return Ok(Vec::new());
    }
    let chol = cholesky(m, n).map_err(|_| ModalError::MassNotSpd { negative: 0 })?;
    // Extract L once (dense lower triangle).
    let mut l = vec![0.0f64; n * n];
    for i in 0..n {
        for j in 0..=i {
            l[i * n + j] = chol.l(i, j);
        }
    }
    // Forward substitution: solve L·X = B for one column vector b in place.
    let fwd = |l: &[f64], b: &mut [f64]| {
        for i in 0..n {
            let mut acc = b[i];
            for j in 0..i {
                acc = l[i * n + j].mul_add(-b[j], acc);
            }
            b[i] = acc / l[i * n + i];
        }
    };
    // Back substitution: solve Lᵀ·x = b in place.
    let bwd = |l: &[f64], b: &mut [f64]| {
        for i in (0..n).rev() {
            let mut acc = b[i];
            for j in i + 1..n {
                acc = l[j * n + i].mul_add(-b[j], acc);
            }
            b[i] = acc / l[i * n + i];
        }
    };
    // Y = L⁻¹·K (column by column), then C = L⁻¹·Yᵀ, symmetrized.
    let mut y = vec![0.0f64; n * n]; // column-major scratch: y[.., col]
    let mut col = vec![0.0f64; n];
    for c in 0..n {
        for (i, slot) in col.iter_mut().enumerate() {
            *slot = k[i * n + c];
        }
        fwd(&l, &mut col);
        for i in 0..n {
            y[i + c * n] = col[i];
        }
    }
    let mut cmat = vec![0.0f64; n * n];
    for c in 0..n {
        // column c of C = L⁻¹ · (row c of Y) = L⁻¹ · Yᵀ e_c.
        for (i, slot) in col.iter_mut().enumerate() {
            *slot = y[c + i * n];
        }
        fwd(&l, &mut col);
        for i in 0..n {
            cmat[i * n + c] = col[i];
        }
    }
    for i in 0..n {
        for j in 0..i {
            let avg = f64::midpoint(cmat[i * n + j], cmat[j * n + i]);
            cmat[i * n + j] = avg;
            cmat[j * n + i] = avg;
        }
    }
    let (vals, vecs) = jacobi_eigh(&cmat, n);

    let mut modes = Vec::with_capacity(n);
    let mut scratch = vec![0.0f64; n];
    for (idx, &lambda) in vals.iter().enumerate() {
        // φ = L⁻ᵀ y (columns of `vecs` are eigenvectors).
        let mut phi: Vec<f64> = (0..n).map(|i| vecs[i * n + idx]).collect();
        bwd(&l, &mut phi);
        // Mass-normalize.
        for (i, slot) in scratch.iter_mut().enumerate() {
            let mut acc = 0.0;
            for (j, p) in phi.iter().enumerate() {
                acc = m[i * n + j].mul_add(*p, acc);
            }
            *slot = acc;
        }
        let mnorm: f64 = phi
            .iter()
            .zip(&scratch)
            .map(|(a, b)| a * b)
            .sum::<f64>()
            .max(f64::MIN_POSITIVE)
            .sqrt();
        for slot in &mut phi {
            *slot /= mnorm;
        }
        // Residual and certificate.
        let mut r = vec![0.0f64; n];
        for (i, slot) in r.iter_mut().enumerate() {
            let mut kk = 0.0;
            let mut mm = 0.0;
            for (j, p) in phi.iter().enumerate() {
                kk = k[i * n + j].mul_add(*p, kk);
                mm = m[i * n + j].mul_add(*p, mm);
            }
            *slot = lambda.mul_add(-mm, kk);
        }
        let mut minv_r = r.clone();
        chol.solve(&mut minv_r);
        let bound = r
            .iter()
            .zip(&minv_r)
            .map(|(a, b)| a * b)
            .sum::<f64>()
            .max(0.0)
            .sqrt();
        modes.push(ModePair {
            lambda,
            phi,
            residual: bound,
            interval: (lambda - bound, lambda + bound),
        });
    }
    Ok(modes)
}

// ---------------------------------------------------------------------------
// Quadratic (damped) eigenvalues
// ---------------------------------------------------------------------------

/// All 2n eigenvalues of the quadratic pencil `(λ²M + λC + K)·φ = 0` (dense
/// row-major inputs, SPD M). Exact power-of-two coordinate and time scales
/// precede a Cholesky mass whitening and dimensionless companion linearization.
/// All damping entries, including gyroscopic terms, reach the fs-la complex QR
/// solver; no approximate diagonalization or modal decoupling is performed.
/// Exactly diagonal pencils use the scaled analytic scalar quadratic and
/// product-of-roots recovery, retaining each slow overdamped relaxation pole.
/// Negative stiffness and active damping retain their unstable poles.
/// Returns eigenvalues sorted canonically by (re, im); eigenvectors are not
/// computed (the fs-la complex path is values-only — recorded boundary).
///
/// # Errors
/// [`ModalError::MassNotSpd`] if M fails Cholesky; [`ModalError::QuadraticEig`]
/// if the complex QR iteration fails to converge or coefficients, transformed
/// entries or resulting poles are not finite and representable. Scaling does
/// not certify accuracy of a highly separated or defective matrix spectrum.
///
/// # Panics
/// On shape mismatches (programmer error).
pub fn quadratic_eigenvalues(
    m: &[f64],
    c: &[f64],
    k: &[f64],
    n: usize,
) -> Result<Vec<C64>, ModalError> {
    assert_eq!(m.len(), n * n, "m must be n*n");
    assert_eq!(c.len(), n * n, "c must be n*n");
    assert_eq!(k.len(), n * n, "k must be n*n");
    if n == 0 {
        return Ok(Vec::new());
    }
    if n == 1 {
        return quadratic_scalar_eigenvalues(m[0], c[0], k[0]);
    }
    let invalid = || ModalError::QuadraticEig {
        detail: "quadratic pencil scaling requires finite, representable coefficients and poles".into(),
    };
    if m.iter().any(|v| !v.is_finite()) || (0..n).any(|i| m[i*n+i] <= 0.0) {
        return Err(ModalError::MassNotSpd { negative: 0 });
    }
    if c.iter().chain(k).any(|v| !v.is_finite()) { return Err(invalid()); }
    if (0..n).all(|i| (0..n).all(|j| i == j
        || (m[i*n+j] == 0.0 && c[i*n+j] == 0.0 && k[i*n+j] == 0.0))) {
        let mut poles = Vec::with_capacity(2*n);
        for i in 0..n {
            poles.extend(quadratic_scalar_eigenvalues(m[i*n+i],c[i*n+i],k[i*n+i])?);
        }
        poles.sort_by(|a,b| a.re.total_cmp(&b.re).then_with(|| a.im.total_cmp(&b.im)));
        return Ok(poles);
    }
    // x_i = 2^(-e_i) y_i equilibrates the mass diagonal without rounding a
    // square root into every coefficient. lambda = 2^t mu keeps BOTH rate
    // and squared-rate blocks representable before any division by mass.
    let exponents: Vec<_> = (0..n).map(|i| quadratic_fraction(m[i*n+i]).1.div_euclid(2)).collect();
    let mut time_exponent = i32::MIN;
    for i in 0..n { for j in 0..n {
        let coordinate = exponents[i] + exponents[j];
        if c[i*n+j] != 0.0 {
            time_exponent = time_exponent.max(quadratic_fraction(c[i*n+j]).1 - coordinate);
        }
        if k[i*n+j] != 0.0 {
            time_exponent = time_exponent.max((quadratic_fraction(k[i*n+j]).1 - coordinate + 1).div_euclid(2));
        }
    }}
    if time_exponent == i32::MIN { time_exponent = 0; }
    let transform = |a: &[f64], order: i32| -> Result<Vec<f64>,ModalError> {
        a.iter().enumerate().map(|(index,&value)| {
            let exponent = -exponents[index/n] - exponents[index%n] - order*time_exponent;
            let scaled = quadratic_power_scale(value,exponent);
            if !scaled.is_finite() || (value != 0.0 && scaled == 0.0) { Err(invalid()) }
            else { Ok(scaled) }
        }).collect()
    };
    let mass = transform(m,0)?;
    let chol = cholesky(&mass,n).map_err(|_| ModalError::MassNotSpd { negative: 0 })?;
    let mut damping = transform(c,1)?;
    let mut stiffness = transform(k,2)?;
    quadratic_mass_whiten(&mut damping,n,&chol);
    quadratic_mass_whiten(&mut stiffness,n,&chol);
    if damping.iter().chain(&stiffness).any(|v| !v.is_finite()) { return Err(invalid()); }
    // z = L^T y gives mu^2 I + mu L^-1 C L^-T + L^-1 K L^-T.
    // Unlike M^-1 K, the whitened symmetric blocks stay symmetric; an
    // original antisymmetric gyroscopic block is retained with its sign.
    let dim = 2 * n;
    let mut a = vec![C64::ZERO; dim * dim];
    for i in 0..n {
        a[i * dim + (n + i)] = C64::new(1.0, 0.0);
        for j in 0..n {
            a[(n+i)*dim+j] = C64::from_re(-stiffness[i*n+j]);
            a[(n+i)*dim+n+j] = C64::from_re(-damping[i*n+j]);
        }
    }
    let mut poles = complex_eig(&a, dim).map_err(|e| ModalError::QuadraticEig {
        detail: format!("{e:?}"),
    })?;
    for pole in &mut poles {
        let original = *pole;
        pole.re = quadratic_power_scale(pole.re,time_exponent);
        pole.im = quadratic_power_scale(pole.im,time_exponent);
        if !pole.re.is_finite() || !pole.im.is_finite()
            || (original != C64::ZERO && *pole == C64::ZERO) { return Err(invalid()); }
    }
    poles.sort_by(|a,b| a.re.total_cmp(&b.re).then_with(|| a.im.total_cmp(&b.im)));
    Ok(poles)
}

// In-place congruence L^-1 A L^-T, retaining nonsymmetric damping too.
fn quadratic_mass_whiten(a: &mut [f64], n: usize, l: &fs_la::factor::Cholesky) {
    let mut column = vec![0.0;n];
    let forward = |column: &mut [f64]| {
        for i in 0..n {
            for j in 0..i { column[i] = (-l.l(i,j)).mul_add(column[j],column[i]); }
            column[i] /= l.l(i,i);
        }
    };
    for j in 0..n {
        for i in 0..n { column[i] = a[i*n+j]; }
        forward(&mut column);
        for i in 0..n { a[i*n+j] = column[i]; }
    }
    for i in 0..n {
        column.copy_from_slice(&a[i*n..(i+1)*n]);
        forward(&mut column);
        a[i*n..(i+1)*n].copy_from_slice(&column);
    }
}

// Exact binary rescaling except at a final subnormal/overflow boundary.
// Chunking avoids constructing 2^exponent when that factor alone overflows.
fn quadratic_power_scale(mut value: f64, mut exponent: i32) -> f64 {
    while exponent != 0 {
        let step = exponent.clamp(-512,512);
        value *= f64::from_bits(((1023+step) as u64) << 52);
        exponent -= step;
    }
    value
}

fn quadratic_fraction(value: f64) -> (f64,i32) {
    let shift = if value.abs() < f64::MIN_POSITIVE { 52 } else { 0 };
    let bits = (value.abs() * if shift == 0 { 1.0 } else { 4_503_599_627_370_496.0 }).to_bits();
    let exponent = ((bits >> 52) & 0x7ff) as i32 - 1023 - shift;
    let mantissa = f64::from_bits((bits & 0x000f_ffff_ffff_ffff) | (1023_u64 << 52));
    (mantissa,exponent)
}

// Solve m lambda^2 + c lambda + k = 0 without squaring an unscaled rate or
// subtracting the nearly equal roots of a strongly overdamped oscillator.
fn quadratic_scalar_eigenvalues(m: f64, c: f64, k: f64) -> Result<Vec<C64>, ModalError> {
    if !m.is_finite() || m <= 0.0 {
        return Err(ModalError::MassNotSpd { negative: 0 });
    }
    let invalid = || ModalError::QuadraticEig {
        detail: "scalar pencil requires finite coefficients and representable poles".into(),
    };
    if !c.is_finite() || !k.is_finite() { return Err(invalid()); }
    let damping = c / m;
    if k == 0.0 && c != 0.0 && damping == 0.0 { return Err(invalid()); }
    let mut poles = if k == 0.0 {
        vec![C64::from_re(-damping), C64::ZERO]
    } else {
        // sqrt(|k|)/sqrt(m) can be finite even when k/m over/underflows.
        let omega = k.abs().sqrt() / m.sqrt();
        let decay = if damping.is_finite() { 0.5 * damping } else { (0.5 * c) / m };
        if !omega.is_finite() || !decay.is_finite() { return Err(invalid()); }
        let scale = decay.abs().max(omega);
        let d = decay / scale;
        let w = omega / scale;
        let discriminant = if k > 0.0 {
            quadratic_scalar_discriminant(m, c, k)
        } else {
            d * d + w * w
        };
        if discriminant < 0.0 {
            let imaginary = scale * (-discriminant).sqrt();
            vec![C64::new(-decay, -imaginary), C64::new(-decay, imaginary)]
        } else if discriminant == 0.0 {
            vec![C64::from_re(-decay); 2]
        } else {
            let far = -scale * (d + discriminant.sqrt().copysign(d));
            // far has magnitude >= sqrt(|k|/m). Form the product quotient
            // in two factors, avoiding both k/m and omega^2 overflow.
            let near = k.signum() * omega * (omega / far);
            vec![C64::from_re(far), C64::from_re(near)]
        }
    };
    if poles.iter().any(|pole| !pole.re.is_finite() || !pole.im.is_finite()
        || (k != 0.0 && pole.re == 0.0 && pole.im == 0.0)) {
        return Err(invalid());
    }
    poles.sort_by(|a, b| a.re.total_cmp(&b.re).then_with(|| a.im.total_cmp(&b.im)));
    Ok(poles)
}

// (c^2 - 4 m k) / max(c^2, 4 m k), with positive m,k. Normalize the
// ORIGINAL coefficients by powers of two, then compensate the difference of
// products. Comparing independently rounded sqrt(k)/sqrt(m) and c/(2m)
// instead would spuriously split exactly critical roots (e.g. 175,770,847).
fn quadratic_scalar_discriminant(m: f64, c: f64, k: f64) -> f64 {
    if c == 0.0 { return -1.0; }
    let (mm, me) = quadratic_fraction(m);
    let (cm, ce) = quadratic_fraction(c);
    let (km, ke) = quadratic_fraction(k);
    let c_exponent = 2 * ce;
    let mk_exponent = me + ke + 2;
    let exponent = c_exponent.max(mk_exponent);
    let scaled_c = cm * 2.0_f64.powi(c_exponent - exponent);
    let scaled_m = mm * 2.0_f64.powi(mk_exponent - exponent);
    let product = scaled_m * km;
    let product_error = (-scaled_m).mul_add(km, product);
    let difference = scaled_c.mul_add(cm, -product) + product_error;
    difference / (scaled_c * cm).max(product)
}

/// Largest M-Gram defect `max |φ_iᵀMφ_j − δ_ij|` a cluster may carry and
/// still be tightened by [`quadratic_radii`].
pub const QUADRATIC_GRAM_DEFECT_MAX: f64 = 1.0e-12;

/// Quadratic (Kato–Temple) eigenvalue radii for an inertia-certified slice.
///
/// The value certificate `|λ − λ̂| ≤ r` has a floating-point floor in `r` on
/// fine meshes. A certified slice also knows its spectral gaps, because
/// `modes.len() == expected` counts every eigenvalue in the window. That
/// admits the quadratic residual bound `|λ_j − θ_j| ≤ ‖R‖² / δ`:
/// Temple (1928) and Kato (1949) for one vector; Mathias, SIAM J. Matrix
/// Anal. Appl. 19(2), 1998 for a cluster of M-orthonormal Ritz vectors.
/// `R` is the cluster's `M⁻¹`-norm residual block and `δ` the distance from
/// its Ritz values to every eigenvalue outside it.
///
/// How the radii are built:
/// 1. With `F² = Σ r_i²` over the whole report, `F ≥ ‖R‖₂`. By Kahan's
///    theorem, the `n` modes match `n` distinct eigenvalues, each within `F`
///    of its Ritz value. If every `[θ_i − F, θ_i + F]` lies inside the window,
///    those are exactly the window's eigenvalues.
/// 2. Overlapping widened intervals form clusters.
/// 3. `δ` is bounded below by the neighbouring clusters' widened intervals
///    and by the window edges. The window counts `[low, high)`; every other
///    eigenvalue is `< low` or `≥ high`.
/// 4. A cluster of several modes is tightened only if its measured M-Gram
///    defect `η` is at most [`QUADRATIC_GRAM_DEFECT_MAX`]. The cluster bound
///    is stated for Rayleigh–Ritz values. The diagonal Ritz values `θ` differ
///    from those by at most `e = ‖QᵀKQ − diag θ‖_F` (Weyl), and the
///    Rayleigh–Ritz residual is no larger than `R`. So the gap becomes
///    `δ − e`, and `e + 2η·max|θ|` is added.
/// 5. Each radius is
///    `min(r_i, Σ_c r² / (δ_c − e) + e + 2η·max|θ| + 4ε·max|θ|)`. The last
///    term covers the rounding of the computed Ritz values themselves.
///
/// Modes that fail any step keep their linear residual. The report's
/// `interval`s, and every digest over them, are unchanged.
///
/// # Panics
/// If a mode shape's length disagrees with `k`/`m` (programmer error).
#[must_use]
pub fn quadratic_radii(report: &SliceReport, k: &Csr, m: &Csr) -> Vec<f64> {
    let modes = &report.modes;
    let mut radii: Vec<f64> = modes.iter().map(|mode| mode.residual).collect();
    let n = modes.len();
    if n == 0 || n != report.expected {
        return radii;
    }
    let (low, high) = report.window;
    let frobenius = modes
        .iter()
        .map(|mode| mode.residual * mode.residual)
        .sum::<f64>()
        .sqrt();
    if !frobenius.is_finite()
        || modes
            .iter()
            .any(|mode| mode.lambda - frobenius < low || mode.lambda + frobenius >= high)
    {
        return radii;
    }
    let mut clusters: Vec<(usize, usize)> = Vec::new();
    let mut start = 0;
    for i in 1..=n {
        if i == n || modes[i].lambda - frobenius > modes[i - 1].lambda + frobenius {
            clusters.push((start, i));
            start = i;
        }
    }
    for (c, &(first, end)) in clusters.iter().enumerate() {
        let below = if c == 0 {
            low
        } else {
            modes[clusters[c - 1].1 - 1].lambda + frobenius
        };
        let above = if c + 1 == clusters.len() {
            high
        } else {
            modes[clusters[c + 1].0].lambda - frobenius
        };
        let gap = (modes[first].lambda - below).min(above - modes[end - 1].lambda);
        if !(gap > 0.0 && gap.is_finite()) {
            continue;
        }
        let mut defect = 0.0_f64;
        let mut off_diagonal_sq = 0.0_f64;
        if end - first > 1 {
            let mut stiff_phi = vec![0.0; k.nrows()];
            let mut mass_phi = vec![0.0; m.nrows()];
            for i in first..end {
                assert_eq!(modes[i].phi.len(), m.nrows(), "mode/mass size mismatch");
                assert_eq!(
                    modes[i].phi.len(),
                    k.nrows(),
                    "mode/stiffness size mismatch"
                );
                k.spmv(&modes[i].phi, &mut stiff_phi);
                m.spmv(&modes[i].phi, &mut mass_phi);
                for (j, other) in modes.iter().enumerate().take(end).skip(first) {
                    let gram: f64 = other.phi.iter().zip(&mass_phi).map(|(a, b)| a * b).sum();
                    let rayleigh: f64 = other.phi.iter().zip(&stiff_phi).map(|(a, b)| a * b).sum();
                    let (identity, diagonal) = if i == j {
                        (1.0, modes[i].lambda)
                    } else {
                        (0.0, 0.0)
                    };
                    defect = defect.max((gram - identity).abs());
                    off_diagonal_sq += (rayleigh - diagonal).powi(2);
                }
            }
            if !(defect <= QUADRATIC_GRAM_DEFECT_MAX) {
                continue;
            }
        }
        let galerkin = off_diagonal_sq.sqrt();
        let effective_gap = gap - galerkin;
        if !(effective_gap > 0.0) {
            continue;
        }
        let block: f64 = modes[first..end]
            .iter()
            .map(|mode| mode.residual * mode.residual)
            .sum();
        let scale = modes[first..end]
            .iter()
            .map(|mode| mode.lambda.abs())
            .fold(0.0_f64, f64::max);
        // The computed Rayleigh quotient itself carries rounding of a few ulps.
        let quadratic =
            block / effective_gap + galerkin + 2.0 * defect * scale + 4.0 * f64::EPSILON * scale;
        for radius in &mut radii[first..end] {
            *radius = radius.min(quadratic);
        }
    }
    radii
}

/// Rayleigh-recentred, quadratically certified eigenvalues of a COMPLETE
/// dense solve, such as [`eigh_gen_dense`].
///
/// `pairs` must hold every mode of the `n×n` pencil `(K, M)`. Each λ is
/// replaced by the Rayleigh quotient of its own mass-normalized vector in the
/// original pencil, which for that vector minimizes the residual, and its
/// `M⁻¹` residual and interval are recomputed. The second item is the
/// per-mode [`quadratic_radii`]. With all `n` modes present the eigenvalue
/// count is exact by dimension, so the whole real line is the certified
/// window.
///
/// Why: an f64 eigenvector of an ill-conditioned pencil keeps a linear
/// residual near `ε·λ_max·‖φ‖` however it is computed, because rounding φ
/// alone excites the stiffest modes. The quadratic radius bounds the
/// eigenvalue error itself. [`eigh_gen_dense`]'s own output is unchanged.
///
/// # Errors
/// [`ModalError::DimensionMismatch`], or a refusal of the mass
/// factorization ([`ModalError::Factor`], [`ModalError::MassNotSpd`]).
///
/// # Panics
/// If `pairs` does not hold one mode per row, or a mode's length disagrees
/// with `k` (programmer error).
pub fn recentred_dense_radii(
    pairs: &[ModePair],
    k: &Csr,
    m: &Csr,
) -> Result<(Vec<ModePair>, Vec<f64>), ModalError> {
    check_shapes(k, m)?;
    let n = k.nrows();
    assert_eq!(
        pairs.len(),
        n,
        "a complete dense solve has one mode per row"
    );
    let opts = SliceOptions::default();
    let m_fact = SymbolicLdlt::analyze_with_coordinates(m, opts.ordering, None)
        .and_then(|symbolic| symbolic.factor(m, &opts.ldlt))
        .map_err(|source| ModalError::Factor {
            shift: f64::NAN,
            source,
        })?;
    if m_fact.inertia().negative > 0 {
        return Err(ModalError::MassNotSpd {
            negative: m_fact.inertia().negative,
        });
    }
    let (mut stiff_phi, mut mass_phi) = (vec![0.0; n], vec![0.0; n]);
    let mut modes = Vec::with_capacity(n);
    for pair in pairs {
        assert_eq!(pair.phi.len(), n, "mode/pencil size mismatch");
        let mut phi = pair.phi.clone();
        m.spmv(&phi, &mut mass_phi);
        let norm = phi
            .iter()
            .zip(&mass_phi)
            .map(|(a, b)| a * b)
            .sum::<f64>()
            .max(f64::MIN_POSITIVE)
            .sqrt();
        for slot in &mut phi {
            *slot /= norm;
        }
        k.spmv(&phi, &mut stiff_phi);
        m.spmv(&phi, &mut mass_phi);
        let lambda = phi.iter().zip(&stiff_phi).map(|(a, b)| a * b).sum::<f64>()
            / phi.iter().zip(&mass_phi).map(|(a, b)| a * b).sum::<f64>();
        let r: Vec<f64> = stiff_phi
            .iter()
            .zip(&mass_phi)
            .map(|(kk, mm)| lambda.mul_add(-mm, *kk))
            .collect();
        let minv_r = m_fact.solve(&r);
        let bound = r
            .iter()
            .zip(&minv_r)
            .map(|(a, b)| a * b)
            .sum::<f64>()
            .max(0.0)
            .sqrt();
        modes.push(ModePair {
            lambda,
            phi,
            residual: bound,
            interval: (lambda - bound, lambda + bound),
        });
    }
    modes.sort_by(|a, b| a.lambda.total_cmp(&b.lambda));
    let report = SliceReport {
        window: (f64::NEG_INFINITY, f64::INFINITY),
        below_low: 0,
        below_high: n,
        expected: n,
        modes,
        stats: SliceStats {
            shift: f64::NAN,
            factorizations: 1,
            lanczos_iters: 0,
            restarts: 0,
            factor_nnz_l: 0,
            factor_peak_bytes: 0,
            pivots_delayed: 0,
        },
    };
    let radii = quadratic_radii(&report, k, m);
    Ok((report.modes, radii))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// tridiag(-1, 2, -1): eigenvalues 2 − 2cos(iπ/(n+1)).
    fn tridiag_k(n: usize) -> Csr {
        let mut coo = Coo::new(n, n);
        for i in 0..n {
            coo.push(i, i, 2.0);
            if i > 0 {
                coo.push(i, i - 1, -1.0);
            }
            if i + 1 < n {
                coo.push(i, i + 1, -1.0);
            }
        }
        coo.assemble()
    }

    fn diag(n: usize, value: f64) -> Csr {
        let mut coo = Coo::new(n, n);
        for i in 0..n {
            coo.push(i, i, value);
        }
        coo.assemble()
    }

    fn analytic_tridiag(n: usize) -> Vec<f64> {
        (1..=n)
            .map(|i| 2.0 - 2.0 * (i as f64 * std::f64::consts::PI / (n as f64 + 1.0)).cos())
            .collect()
    }

    #[test]
    fn mass_coordinate_equilibration_preserves_the_original_pencil_certificates() {
        // This is one generalized pencil written in coordinates with very
        // different units. Its exact eigenvalues remain the tridiagonal ones.
        let n = 12;
        let coordinates = [0.2, 0.5, 1.0, 2.0, 5.0];
        let scale: Vec<f64> = (0..n).map(|i| coordinates[i % coordinates.len()]).collect();
        let mut k_coo = Coo::new(n, n);
        let mut m_coo = Coo::new(n, n);
        for i in 0..n {
            k_coo.push(i, i, 2.0 / scale[i].powi(2));
            m_coo.push(i, i, 1.0 / scale[i].powi(2));
            if i + 1 < n {
                let coupling = -1.0 / (scale[i] * scale[i + 1]);
                k_coo.push(i, i + 1, coupling);
                k_coo.push(i + 1, i, coupling);
            }
        }
        let k = k_coo.assemble();
        let m = m_coo.assemble();
        let truth = analytic_tridiag(n);
        let window = (
            f64::midpoint(truth[0], truth[1]),
            f64::midpoint(truth[5], truth[6]),
        );
        for mass_diagonal_equilibration in [false, true] {
            let options = SliceOptions {
                mass_diagonal_equilibration,
                ..SliceOptions::default()
            };
            let report = slice_window(&k, &m, window, &options).expect("certified slice");
            assert_eq!(report.expected, 5);
            assert_eq!(report.below_low, 1);
            for (mode, expected) in report.modes.iter().zip(&truth[1..6]) {
                // The analytic cosine is itself evaluated in f64; allow its
                // rounding interval to overlap the computed residual interval.
                let oracle_roundoff = 4.0 * f64::EPSILON * expected.abs().max(1.0);
                assert!(
                    mode.interval.0 <= *expected + oracle_roundoff
                        && *expected - oracle_roundoff <= mode.interval.1,
                    "equilibrated={mass_diagonal_equilibration}, lambda={}, interval={:?}, analytic={expected}",
                    mode.lambda,
                    mode.interval
                );
                assert!((mode.lambda - expected).abs() < 1e-12);
                let mut mass_phi = vec![0.0; n];
                m.spmv(&mode.phi, &mut mass_phi);
                let mass_norm: f64 = mode.phi.iter().zip(&mass_phi).map(|(a, b)| a * b).sum();
                assert!((mass_norm - 1.0).abs() < 1e-10);
            }
        }
    }

    #[test]
    fn spring_mass_chain_window_is_certified_and_intervals_contain_truth() {
        let n = 40;
        let mass = 0.5;
        let k = tridiag_k(n);
        let m = diag(n, mass);
        // Pencil eigenvalues = analytic/mass. Slice modes 3..=7 (1-based).
        let truth: Vec<f64> = analytic_tridiag(n).iter().map(|l| l / mass).collect();
        let lo = f64::midpoint(truth[1], truth[2]);
        let hi = f64::midpoint(truth[6], truth[7]);
        let rep = slice_window(&k, &m, (lo, hi), &SliceOptions::default()).expect("slice");
        assert_eq!(rep.expected, 5);
        assert_eq!(rep.below_low, 2);
        assert_eq!(rep.modes.len(), 5);
        for (mode, want) in rep.modes.iter().zip(&truth[2..7]) {
            assert!(
                mode.interval.0 <= *want && *want <= mode.interval.1,
                "certified interval [{}, {}] must contain analytic {want}",
                mode.interval.0,
                mode.interval.1
            );
            // Mass-normalization is enforced.
            let mut mv = vec![0.0; n];
            m.spmv(&mode.phi, &mut mv);
            let mnorm: f64 = mode.phi.iter().zip(&mv).map(|(a, b)| a * b).sum();
            assert!((mnorm - 1.0).abs() < 1e-10, "phi'Mphi = {mnorm}");
        }
        println!(
            "{{\"suite\":\"fs-modal\",\"case\":\"spring-mass-window\",\"expected\":5,\"iters\":{},\"verdict\":\"pass\"}}",
            rep.stats.lanczos_iters
        );
    }

    #[test]
    fn plate_five_point_squared_identity() {
        // Simply-supported plate via the squared-Laplacian identity:
        // K = L², so pencil eigenvalues are λ_L² — an independent analytic
        // oracle exercising a genuinely 2D, degenerate-cluster spectrum.
        let s = 8;
        let n = s * s;
        let mut coo = Coo::new(n, n);
        for i in 0..s {
            for j in 0..s {
                let u = i * s + j;
                coo.push(u, u, 4.0);
                if i > 0 {
                    coo.push(u, u - s, -1.0);
                }
                if i + 1 < s {
                    coo.push(u, u + s, -1.0);
                }
                if j > 0 {
                    coo.push(u, u - 1, -1.0);
                }
                if j + 1 < s {
                    coo.push(u, u + 1, -1.0);
                }
            }
        }
        let lap = coo.assemble();
        let k2 = fs_sparse::ops::spgemm(&lap, &lap);
        let m = diag(n, 1.0);
        let mut truth: Vec<f64> = Vec::new();
        for i in 1..=s {
            for j in 1..=s {
                let li = 2.0 - 2.0 * (i as f64 * std::f64::consts::PI / (s as f64 + 1.0)).cos();
                let lj = 2.0 - 2.0 * (j as f64 * std::f64::consts::PI / (s as f64 + 1.0)).cos();
                let l = li + lj;
                truth.push(l * l);
            }
        }
        truth.sort_by(f64::total_cmp);
        // Window around modes 2..=6 — includes TWO DEGENERATE pairs
        // (λ_ij = λ_ji), forcing the deflation-restart path. Endpoints must
        // split DISTINCT values: truth[4] == truth[5] is itself a degenerate
        // pair, so a midpoint there would BE an eigenvalue and the endpoint
        // factorization would (correctly) refuse as singular.
        assert!(truth[0] < truth[1] && truth[5] < truth[6], "fixture layout");
        let lo = f64::midpoint(truth[0], truth[1]);
        let hi = f64::midpoint(truth[5], truth[6]);
        let rep = slice_window(&k2, &m, (lo, hi), &SliceOptions::default()).expect("slice");
        assert_eq!(rep.expected, 5);
        for (mode, want) in rep.modes.iter().zip(&truth[1..6]) {
            assert!(
                mode.interval.0 <= *want && *want <= mode.interval.1,
                "interval must contain analytic {want}, got [{}, {}]",
                mode.interval.0,
                mode.interval.1
            );
        }
        // NOTE: full reorthogonalization often recovers degenerate copies
        // through rounding without a restart — the restart path is proven
        // deterministically by `exact_degeneracy_recovered_by_deflation`.
        println!(
            "{{\"suite\":\"fs-modal\",\"case\":\"plate-5pt-squared\",\"expected\":5,\"restarts\":{},\"verdict\":\"pass\"}}",
            rep.stats.restarts
        );
    }

    #[test]
    fn dense_path_matches_analytic_and_slice() {
        let n = 12;
        let k = tridiag_k(n);
        let m = diag(n, 2.0);
        let modes = eigh_gen_dense(&k.to_dense(), &m.to_dense(), n).expect("dense");
        let truth: Vec<f64> = analytic_tridiag(n).iter().map(|l| l / 2.0).collect();
        assert_eq!(modes.len(), n);
        for (mode, want) in modes.iter().zip(&truth) {
            assert!(
                mode.interval.0 <= *want && *want <= mode.interval.1,
                "dense interval must contain {want}"
            );
            assert!(mode.residual < 1e-10);
        }
        // Cross-strategy agreement on a window.
        let rep = slice_window(
            &k,
            &m,
            (truth[2] - 0.01, truth[5] + 0.01),
            &SliceOptions::default(),
        )
        .expect("slice");
        // The half-open window (truth[2] − 0.01, truth[5] + 0.01] contains
        // exactly the four eigenvalues with indices 2..=5.
        assert_eq!(rep.expected, 4);
        for (sparse_mode, dense_mode) in rep.modes.iter().zip(&modes[2..6]) {
            assert!((sparse_mode.lambda - dense_mode.lambda).abs() < 1e-9);
        }
    }

    #[test]
    fn metamorphic_shift_invariance_and_scaling_covariance() {
        let n = 30;
        let k = tridiag_k(n);
        let m = diag(n, 1.0);
        let truth = analytic_tridiag(n);
        let (lo, hi) = (
            f64::midpoint(truth[3], truth[4]),
            f64::midpoint(truth[8], truth[9]),
        );
        let base = slice_window(&k, &m, (lo, hi), &SliceOptions::default()).expect("base");

        // eig(K + cM, M) = eig(K, M) + c.
        let c = 3.25;
        let shifted_k = {
            let mut coo = Coo::new(n, n);
            for r in 0..n {
                let (cols, vals) = k.row(r);
                for (&cc, &v) in cols.iter().zip(vals) {
                    coo.push(r, cc, v);
                }
                coo.push(r, r, c);
            }
            coo.assemble()
        };
        let rep = slice_window(&shifted_k, &m, (lo + c, hi + c), &SliceOptions::default())
            .expect("shifted");
        assert_eq!(rep.expected, base.expected);
        for (a, b) in rep.modes.iter().zip(&base.modes) {
            assert!(((a.lambda - c) - b.lambda).abs() < 1e-8, "shift invariance");
        }

        // eig(sK, M) = s·eig(K, M).
        let s = 4.0;
        let scaled_k = {
            let mut coo = Coo::new(n, n);
            for r in 0..n {
                let (cols, vals) = k.row(r);
                for (&cc, &v) in cols.iter().zip(vals) {
                    coo.push(r, cc, s * v);
                }
            }
            coo.assemble()
        };
        let rep = slice_window(&scaled_k, &m, (s * lo, s * hi), &SliceOptions::default())
            .expect("scaled");
        assert_eq!(rep.expected, base.expected);
        for (a, b) in rep.modes.iter().zip(&base.modes) {
            assert!((a.lambda / s - b.lambda).abs() < 1e-8, "scaling covariance");
        }
    }

    #[test]
    fn exact_degeneracy_recovered_by_deflation() {
        // K = diag(1, 1): the whole space is one degenerate eigenplane, so
        // the first Lanczos step yields w = 0 EXACTLY (no other eigenvalue
        // leaves a rounding residue for later iterations to amplify) and
        // the sweep breaks down holding ONE copy. The second copy is
        // reachable only through a deflated restart — the deterministic
        // proof that the restart path is load-bearing. (With mixed spectra,
        // full-reorth rounding often surfaces copies without a restart;
        // that path is exercised by the plate fixture.)
        let n = 2;
        let mut coo = Coo::new(n, n);
        coo.push(0, 0, 1.0);
        coo.push(1, 1, 1.0);
        let k = coo.assemble();
        let m = diag(n, 1.0);
        let rep = slice_window(&k, &m, (0.5, 1.5), &SliceOptions::default()).expect("slice");
        assert_eq!(rep.expected, 2, "inertia certifies the double eigenvalue");
        assert_eq!(rep.modes.len(), 2);
        assert!(
            rep.stats.restarts >= 1,
            "the second copy is unreachable without a deflated restart"
        );
        for mode in &rep.modes {
            assert!((mode.lambda - 1.0).abs() < 1e-10);
        }
        // The two recovered shapes must be genuinely independent
        // (M-orthogonal, here M = I).
        let dot: f64 = rep.modes[0]
            .phi
            .iter()
            .zip(&rep.modes[1].phi)
            .map(|(a, b)| a * b)
            .sum();
        assert!(
            dot.abs() < 1e-8,
            "deflated copies must be orthogonal, dot = {dot}"
        );
        println!(
            "{{\"suite\":\"fs-modal\",\"case\":\"exact-degeneracy-deflation\",\"restarts\":{},\"verdict\":\"pass\"}}",
            rep.stats.restarts
        );
    }

    #[test]
    fn mass_not_spd_refuses() {
        let n = 6;
        let k = tridiag_k(n);
        let mut coo = Coo::new(n, n);
        for i in 0..n {
            coo.push(i, i, if i == 3 { -1.0 } else { 1.0 });
        }
        let m = coo.assemble();
        let err = slice_window(&k, &m, (0.1, 1.0), &SliceOptions::default()).unwrap_err();
        assert!(matches!(err, ModalError::MassNotSpd { negative: 1 }));
        assert!(err.to_string().contains("FS-MODAL-MASS-NOT-SPD"));
    }

    #[test]
    fn starved_budget_fires_the_count_certificate() {
        // MUTATION GATE: a harvester that cannot deliver the inertia-
        // certified count must refuse loudly, never return a short list.
        let n = 40;
        let k = tridiag_k(n);
        let m = diag(n, 1.0);
        let truth = analytic_tridiag(n);
        let opts = SliceOptions {
            max_lanczos: 3,
            max_restarts: 1,
            ..SliceOptions::default()
        };
        let err = slice_window(
            &k,
            &m,
            (
                f64::midpoint(truth[1], truth[2]),
                f64::midpoint(truth[9], truth[10]),
            ),
            &opts,
        )
        .unwrap_err();
        match err {
            ModalError::WindowUnresolved {
                expected,
                converged,
                ..
            } => {
                assert_eq!(expected, 8);
                assert!(converged < 8);
            }
            other => panic!("expected WindowUnresolved, got {other}"),
        }
    }

    #[test]
    fn invalid_window_and_empty_window() {
        let k = tridiag_k(8);
        let m = diag(8, 1.0);
        assert!(matches!(
            slice_window(&k, &m, (2.0, 1.0), &SliceOptions::default()),
            Err(ModalError::InvalidWindow { .. })
        ));
        // A window with no eigenvalues returns an EMPTY certified report.
        let rep = slice_window(&k, &m, (-1.0, -0.5), &SliceOptions::default()).expect("empty");
        assert_eq!(rep.expected, 0);
        assert!(rep.modes.is_empty());
    }

    #[test]
    fn matrix_free_core_matches_assembled_path() {
        let n = 25;
        let k = tridiag_k(n);
        let m = diag(n, 1.0);
        let truth = analytic_tridiag(n);
        let sigma = f64::midpoint(truth[4], truth[5]);
        let a = shifted_pencil(&k, &m, sigma);
        let sym = SymbolicLdlt::analyze(&a, DirectOrdering::Amd).expect("analyze");
        let fac = sym.factor(&a, &LdltOptions::default()).expect("factor");
        let got = shift_invert_modes_mfree(
            n,
            sigma,
            4,
            200,
            |x, y| m.spmv(x, y),
            |b, out| out.copy_from_slice(&fac.solve(b)),
        )
        .expect("mfree");
        assert_eq!(got.len(), 4);
        // The four nearest analytic eigenvalues to sigma, ascending.
        let mut nearest: Vec<f64> = truth.clone();
        nearest.sort_by(|x, y| (x - sigma).abs().total_cmp(&(y - sigma).abs()));
        let mut want: Vec<f64> = nearest[..4].to_vec();
        want.sort_by(f64::total_cmp);
        for (pair, w) in got.iter().zip(&want) {
            assert!(
                (pair.0 - w).abs() < 1e-8,
                "mfree {} vs analytic {w}",
                pair.0
            );
        }
    }

    #[test]
    fn quadratic_scalar_strong_damping_preserves_slow_and_unstable_poles() {
        // The unbalanced 2x2 companion previously returned a zero slow pole:
        // sqrt((2e10)^2 - 4) rounds to 2e10. Both signs of C and K matter.
        for scale in [1.0e-250, 1.0, 1.0e250] {
            for (c, k, expected) in [
                (4.0e10, 4.0, [-4.0e10, -1.0e-10]),
                (-4.0e10, 4.0, [1.0e-10, 4.0e10]),
                (4.0e10, -4.0, [-4.0e10, 1.0e-10]),
                (-4.0e10, -4.0, [-1.0e-10, 4.0e10]),
            ] {
                let poles = quadratic_scalar_eigenvalues(scale, scale * c, scale * k).unwrap();
                for (pole, expected) in poles.iter().zip(expected) {
                    assert_eq!(pole.im, 0.0);
                    assert!((pole.re / expected - 1.0).abs() < 2.0e-14,
                        "m={scale} c={c} k={k}: {pole:?} expected {expected}");
                }
            }
        }
    }

    #[test]
    fn quadratic_scalar_preserves_oscillatory_critical_free_and_buckling_roots() {
        for (m, c, k, expected) in [
            (2.0, 4.0, 10.0, [(-1.0, -2.0), (-1.0, 2.0)]),
            (2.0, -4.0, 10.0, [(1.0, -2.0), (1.0, 2.0)]),
            (4.0, 8.0, 4.0, [(-1.0, 0.0), (-1.0, 0.0)]),
            (175.0, 770.0, 847.0, [(-2.2, 0.0), (-2.2, 0.0)]),
            (1.0, 0.0, -9.0, [(-3.0, 0.0), (3.0, 0.0)]),
            (2.0, 4.0, 0.0, [(-2.0, 0.0), (0.0, 0.0)]),
            (2.0, -4.0, 0.0, [(0.0, 0.0), (2.0, 0.0)]),
            (2.0, 0.0, 0.0, [(0.0, 0.0), (0.0, 0.0)]),
        ] {
            let poles = quadratic_scalar_eigenvalues(m, c, k).unwrap();
            for (pole, (re, im)) in poles.iter().zip(expected) {
                assert!((pole.re - re).abs() < 2.0e-14 && (pole.im - im).abs() < 2.0e-14,
                    "m={m} c={c} k={k}: {pole:?} expected ({re},{im})");
            }
        }
        // One representable step on either side must retain its actual regime;
        // there is no tolerance band snapped to critical damping.
        let below = f64::from_bits(770.0_f64.to_bits() - 1);
        let above = f64::from_bits(770.0_f64.to_bits() + 1);
        let underdamped = quadratic_scalar_eigenvalues(175.0, below, 847.0).unwrap();
        let overdamped = quadratic_scalar_eigenvalues(175.0, above, 847.0).unwrap();
        assert!(underdamped[0].im < 0.0 && underdamped[1].im > 0.0);
        assert!(overdamped.iter().all(|pole| pole.im == 0.0));
        assert!(overdamped[0].re < overdamped[1].re);
        for (m, k, frequency) in [(1.0e-200, 1.0e200, 1.0e200),
            (1.0e200, 1.0e-200, 1.0e-200)] {
            let poles = quadratic_scalar_eigenvalues(m, 0.0, k).unwrap();
            assert_eq!(poles[0].re, 0.0);
            assert_eq!(poles[1].re, 0.0);
            assert!((poles[0].im / -frequency - 1.0).abs() < 2.0e-14);
            assert!((poles[1].im / frequency - 1.0).abs() < 2.0e-14);
        }
    }

    #[test]
    fn quadratic_scalar_refuses_invalid_or_unrepresentable_poles() {
        for mass in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(matches!(quadratic_scalar_eigenvalues(mass, 1.0, 1.0),
                Err(ModalError::MassNotSpd { .. })));
        }
        for (m, c, k) in [(1.0, f64::NAN, 1.0), (1.0, 1.0, f64::INFINITY),
            (1.0e-300, 1.0e300, 0.0), (1.0, 1.0e308, 1.0e-308),
            (1.0e308, 1.0e-308, 0.0)] {
            assert!(matches!(quadratic_scalar_eigenvalues(m, c, k),
                Err(ModalError::QuadraticEig { .. })));
        }
        // The public entry point dispatches to the same scalar implementation.
        assert_eq!(quadratic_eigenvalues(&[1.0], &[4.0e10], &[4.0], 1).unwrap(),
            quadratic_scalar_eigenvalues(1.0, 4.0e10, 4.0).unwrap());
    }

    #[test]
    fn quadratic_matrix_scaling_preserves_nonproportional_damped_resonances() {
        // G1: in the original coordinates M=I, K=diag(4,9), and
        // C=[[8,sqrt(46)],[sqrt(46),12]]/5 is positive definite and
        // nonproportional. Its determinant factors INDEPENDENTLY as
        // (lambda^2 + lambda + 6)(lambda^2 + 3 lambda + 6).
        // Apply the nonorthogonal physical chart T=[[1,1/2],[0,2]] before
        // changing the time unit. Neither damping coupling may be discarded.
        let b = 46.0_f64.sqrt() / 5.0;
        let mass = [1.0, 0.5, 0.5, 4.25];
        let damping = [1.6, 0.8 + 2.0*b, 0.8 + 2.0*b, 10.0 + 2.0*b];
        let stiffness = [4.0, 2.0, 2.0, 37.0];
        let expected = [(-1.5, -3.75_f64.sqrt()), (-1.5, 3.75_f64.sqrt()),
            (-0.5, -5.75_f64.sqrt()), (-0.5, 5.75_f64.sqrt())];
        for (amplitude, rate, chart) in [
            (1.0, 1.0, [1.0, 1.0]),
            (1.0e-200, 1.0e200, [1.0, 1.0]),
            (1.0e200, 1.0e-200, [1.0, 1.0]),
            (1.0e-250, 1.0, [1.0, 1.0]),
            (1.0e250, 1.0, [1.0, 1.0]),
            (1.0, 1.0, [1.0e100, -1.0e-100]),
        ] {
            let transform = |a: &[f64], factor: f64| -> Vec<f64> {
                (0..4).map(|i| a[i] * factor * chart[i/2] * chart[i%2]).collect()
            };
            let m = transform(&mass, amplitude);
            let c = transform(&damping, amplitude * rate);
            let k = transform(&stiffness, (amplitude * rate) * rate);
            let poles = quadratic_eigenvalues(&m, &c, &k, 2).unwrap();
            let mut normalized: Vec<_> = poles.iter().map(|p| C64::new(p.re/rate,p.im/rate)).collect();
            // Pair by the exact damping class and conjugate sign, independent
            // of a last-bit real tie between members of the same pair.
            normalized.sort_by(|a,b| a.im.total_cmp(&b.im));
            let mut expected = expected;
            expected.sort_by(|a,b| a.1.total_cmp(&b.1));
            for (actual,(re,im)) in normalized.iter().zip(expected) {
                assert!((actual.re-re).abs()<2e-12 && (actual.im-im).abs()<2e-12,
                    "amplitude={amplitude} rate={rate} chart={chart:?}: {actual:?} expected ({re},{im})");
            }
        }
    }

    #[test]
    fn quadratic_diagonal_bank_preserves_slow_critical_and_oscillatory_modes() {
        let m = [1.0,0.0,0.0, 0.0,2.0,0.0, 0.0,0.0,3.0];
        let c = [4e10,0.0,0.0, 0.0,0.0,0.0, 0.0,0.0,6.0];
        let k = [4.0,0.0,0.0, 0.0,18.0,0.0, 0.0,0.0,3.0];
        let roots = quadratic_eigenvalues(&m,&c,&k,3).unwrap();
        let expected: [(f64,f64);6] = [(-4e10,0.0),(-1.0,0.0),(-1.0,0.0),(-1e-10,0.0),(0.0,-3.0),(0.0,3.0)];
        for (root,(re,im)) in roots.iter().zip(expected) {
            let scale = re.abs().max(im.abs());
            assert!((root.re-re).abs()<=2e-14*scale && (root.im-im).abs()<=2e-14*scale,
                "{root:?} expected ({re},{im})");
        }
    }

    #[test]
    fn quadratic_mass_whitening_retains_gyroscopic_frequency_splitting() {
        // G1: M=I, K=diag(4,9), C=[[0,3],[-3,0]] is conservative
        // gyroscopic coupling. det Q=lambda^4+22 lambda^2+36, so
        // omega^2=11 +/- sqrt(85). Symmetrizing C would erase the split.
        // The same nonorthogonal chart as above makes M genuinely coupled.
        let frequencies = [(11.0-85.0_f64.sqrt()).sqrt(),(11.0+85.0_f64.sqrt()).sqrt()];
        for rate in [1.0,1e150,1e-150] {
            let mass: Vec<_> = [1.0,0.5,0.5,4.25].iter().map(|m|m/rate).collect();
            let stiffness: Vec<_> = [4.0,2.0,2.0,37.0].iter().map(|k|k*rate).collect();
            let poles = quadratic_eigenvalues(&mass,&[0.0,6.0,-6.0,0.0],&stiffness,2).unwrap();
            let mut positive: Vec<_> = poles.iter().filter(|p|p.im>0.0).map(|p| {
                assert!((p.re/rate).abs()<2e-12);
                p.im/rate
            }).collect();
            positive.sort_by(f64::total_cmp);
            assert_eq!(positive.len(),2);
            for (actual,expected) in positive.iter().zip(frequencies) {
                assert!((actual-expected).abs()<2e-12,"{actual} expected {expected}");
            }
        }
    }

    #[test]
    fn quadratic_light_damping_matches_per_mode_roots() {
        // Proportional damping C = αK: each undamped ω² = λ maps to the
        // exact quadratic roots μ = (−αλ ± sqrt(α²λ² − 4λ))/2.
        let n = 3;
        let k = tridiag_k(n);
        let m = diag(n, 1.0);
        let alpha = 0.02;
        let kd = k.to_dense();
        let md = m.to_dense();
        let cd: Vec<f64> = kd.iter().map(|v| alpha * v).collect();
        let mut roots = quadratic_eigenvalues(&md, &cd, &kd, n).expect("quadratic");
        let mut want: Vec<(f64, f64)> = Vec::new();
        for lambda in analytic_tridiag(n) {
            let disc = alpha * alpha * lambda * lambda - 4.0 * lambda;
            assert!(disc < 0.0, "light damping keeps modes underdamped");
            let re = -alpha * lambda / 2.0;
            let im = (-disc).sqrt() / 2.0;
            want.push((re, -im));
            want.push((re, im));
        }
        // Sort by imaginary part (all magnitudes distinct here): the ± pair
        // real parts agree only to roundoff, so a (re, im) sort would
        // misalign conjugates by last-ulp ties.
        want.sort_by(|a, b| a.1.total_cmp(&b.1));
        roots.sort_by(|a, b| a.im.total_cmp(&b.im));
        assert_eq!(roots.len(), want.len());
        for (r, w) in roots.iter().zip(&want) {
            assert!(
                (r.re - w.0).abs() < 1e-9 && (r.im - w.1).abs() < 1e-9,
                "quadratic root ({}, {}) vs exact ({}, {})",
                r.re,
                r.im,
                w.0,
                w.1
            );
        }
    }

    #[test]
    fn slicing_is_bitwise_deterministic() {
        let n = 30;
        let k = tridiag_k(n);
        let m = diag(n, 1.0);
        let truth = analytic_tridiag(n);
        let win = (
            f64::midpoint(truth[2], truth[3]),
            f64::midpoint(truth[6], truth[7]),
        );
        let a = slice_window(&k, &m, win, &SliceOptions::default()).expect("a");
        let b = slice_window(&k, &m, win, &SliceOptions::default()).expect("b");
        for (x, y) in a.modes.iter().zip(&b.modes) {
            assert_eq!(x.lambda.to_bits(), y.lambda.to_bits());
            assert!(
                x.phi
                    .iter()
                    .zip(&y.phi)
                    .all(|(p, q)| p.to_bits() == q.to_bits()),
                "mode shapes must be bitwise identical across reruns"
            );
        }
        println!(
            "{{\"suite\":\"fs-modal\",\"case\":\"bitwise-determinism\",\"verdict\":\"pass\"}}"
        );
    }

    /// Re-certify a perturbed shape under `M = I`: Rayleigh quotient,
    /// Euclidean residual, linear interval.
    fn recertify(k: &Csr, phi: &[f64]) -> ModePair {
        let norm = phi.iter().map(|v| v * v).sum::<f64>().sqrt();
        let phi: Vec<f64> = phi.iter().map(|v| v / norm).collect();
        let mut kphi = vec![0.0; phi.len()];
        k.spmv(&phi, &mut kphi);
        let lambda: f64 = phi.iter().zip(&kphi).map(|(a, b)| a * b).sum();
        let residual = kphi
            .iter()
            .zip(&phi)
            .map(|(kv, v)| (kv - lambda * v).powi(2))
            .sum::<f64>()
            .sqrt();
        ModePair {
            lambda,
            phi,
            residual,
            interval: (lambda - residual, lambda + residual),
        }
    }

    #[test]
    fn quadratic_radii_contain_the_truth_far_below_the_linear_bound() {
        // Perturbed Ritz vectors of tridiag(-1, 2, -1) carry residuals near
        // 1e-4; the quadratic radius r^2/gap must still contain the analytic
        // eigenvalue while sitting orders of magnitude below r.
        let n = 40;
        let k = tridiag_k(n);
        let m = diag(n, 1.0);
        let truth = analytic_tridiag(n);
        let window = (
            f64::midpoint(truth[0], truth[1]),
            f64::midpoint(truth[5], truth[6]),
        );
        let mut report = slice_window(&k, &m, window, &SliceOptions::default()).expect("slice");
        assert_eq!(report.modes.len(), 5);
        for mode in &mut report.modes {
            let shifted: Vec<f64> = mode
                .phi
                .iter()
                .enumerate()
                .map(|(j, v)| v + 1.0e-4 * (0.37 * j as f64).sin())
                .collect();
            *mode = recertify(&k, &shifted);
        }
        let radii = quadratic_radii(&report, &k, &m);
        for ((mode, &radius), &exact) in report.modes.iter().zip(&radii).zip(&truth[1..6]) {
            assert!(mode.residual > 1.0e-6, "the perturbation must be live");
            // r^2/gap with r ~ 1e-4 and gaps ~ 1e-2 tightens by r/gap ~ 1e-2.
            assert!(
                radius < 5.0e-2 * mode.residual,
                "{radius} vs {}",
                mode.residual
            );
            assert!(
                (mode.lambda - exact).abs() <= radius,
                "{} vs {exact} at radius {radius}",
                mode.lambda
            );
        }
    }

    #[test]
    fn quadratic_radii_tighten_an_orthonormal_cluster_and_refuse_a_defective_one() {
        // K = diag(1, 1 + 1e-6, 3, 4): a near-double pair inside (0.5, 2.0).
        // Hand-built orthonormal Ritz pair (1, 0, e, 0) and (0, 1, 0, e) has
        // residuals ~2e-3 >> the pair's 1e-6 split, so the widened intervals
        // overlap and the two modes form one cluster.
        let mut coo = Coo::new(4, 4);
        for (i, v) in [1.0, 1.0 + 1.0e-6, 3.0, 4.0].into_iter().enumerate() {
            coo.push(i, i, v);
        }
        let k = coo.assemble();
        let m = diag(4, 1.0);
        // A real report supplies the stats; its window and modes are replaced.
        let mut report = slice_window(
            &tridiag_k(8),
            &diag(8, 1.0),
            (0.05, 0.5),
            &SliceOptions::default(),
        )
        .expect("donor report");
        let e = 1.0e-3;
        report.window = (0.5, 2.0);
        report.expected = 2;
        report.modes = vec![
            recertify(&k, &[1.0, 0.0, e, 0.0]),
            recertify(&k, &[0.0, 1.0, 0.0, e]),
        ];
        let radii = quadratic_radii(&report, &k, &m);
        for ((mode, &radius), exact) in report.modes.iter().zip(&radii).zip([1.0, 1.0 + 1.0e-6]) {
            assert!(mode.residual > 1.0e-3, "the perturbation must be live");
            assert!(
                radius < 2.0e-2 * mode.residual,
                "{radius} vs {}",
                mode.residual
            );
            assert!(
                (mode.lambda - exact).abs() <= radius,
                "{} vs {exact} at radius {radius}",
                mode.lambda
            );
        }
        // Two copies of one shape are not an M-orthonormal cluster: no tightening.
        let mut defective = report.clone();
        defective.modes = vec![
            recertify(&k, &[1.0, 0.0, e, 0.0]),
            recertify(&k, &[1.0, 0.0, e, 0.0]),
        ];
        let kept = quadratic_radii(&defective, &k, &m);
        for (mode, &radius) in defective.modes.iter().zip(&kept) {
            assert_eq!(radius.to_bits(), mode.residual.to_bits());
        }
    }

    #[test]
    fn quadratic_radii_keep_the_linear_bound_when_a_mode_reaches_the_window_edge() {
        let n = 30;
        let k = tridiag_k(n);
        let m = diag(n, 1.0);
        let truth = analytic_tridiag(n);
        let mut report = slice_window(
            &k,
            &m,
            (
                f64::midpoint(truth[0], truth[1]),
                f64::midpoint(truth[3], truth[4]),
            ),
            &SliceOptions::default(),
        )
        .expect("slice");
        let shifted: Vec<f64> = report.modes[0]
            .phi
            .iter()
            .enumerate()
            .map(|(j, v)| v + 1.0e-3 * (0.37 * j as f64).sin())
            .collect();
        report.modes[0] = recertify(&k, &shifted);
        // Move the low edge inside the first mode's widened interval.
        report.window.0 = report.modes[0].lambda - 0.5 * report.modes[0].residual;
        let radii = quadratic_radii(&report, &k, &m);
        for (mode, &radius) in report.modes.iter().zip(&radii) {
            assert_eq!(radius.to_bits(), mode.residual.to_bits());
        }
    }
}
