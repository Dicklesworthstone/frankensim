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

use fs_la::eigen::jacobi_eigh;
use fs_la::eigen_complex::eig as complex_eig;
use fs_la::factor::cholesky;
use fs_math::c64::C64;
use fs_sparse::{Coo, Csr, DirectOrdering, LdltError, LdltFactor, LdltOptions, SymbolicLdlt};

/// Crate version, re-exported for provenance stamping.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

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
}

impl Default for SliceOptions {
    fn default() -> SliceOptions {
        SliceOptions {
            max_lanczos: 0,
            max_restarts: 0,
            ritz_tol: 1e-10,
            ldlt: LdltOptions::default(),
            ordering: DirectOrdering::Amd,
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
            if let Some((lo, hi)) = window {
                if !(lambda > lo && lambda <= hi) {
                    continue;
                }
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
        let stop = beta <= 1e-14 || basis.len() == n || iters == budget || ready.len() >= need;
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
    let n = check_shapes(k, m)?;
    let (lo, hi) = window;
    if !(lo.is_finite() && hi.is_finite() && lo < hi) {
        return Err(ModalError::InvalidWindow { low: lo, high: hi });
    }

    // Mass factorization: SPD gate + the M⁻¹ norm used by every certificate.
    let m_sym = SymbolicLdlt::analyze(m, opts.ordering).map_err(|source| ModalError::Factor {
        shift: f64::NAN,
        source,
    })?;
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

    // One symbolic analysis serves every shift (union pattern).
    let a_lo = shifted_pencil(k, m, lo);
    let p_sym = SymbolicLdlt::analyze(&a_lo, opts.ordering)
        .map_err(|source| ModalError::Factor { shift: lo, source })?;
    let mut factorizations = 0usize;
    let factor_at = |sigma: f64, mat: &Csr, count: &mut usize| -> Result<LdltFactor, ModalError> {
        *count += 1;
        p_sym
            .factor(mat, &opts.ldlt)
            .map_err(|source| ModalError::Factor {
                shift: sigma,
                source,
            })
    };

    let f_lo = factor_at(lo, &a_lo, &mut factorizations)?;
    let a_hi = shifted_pencil(k, m, hi);
    let f_hi = factor_at(hi, &a_hi, &mut factorizations)?;
    let below_low = f_lo.inertia().negative;
    let below_high = f_hi.inertia().negative;
    let expected = below_high - below_low;
    drop(f_lo);
    drop(f_hi);

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
        let a_mid = shifted_pencil(k, m, sigma);
        match factor_at(sigma, &a_mid, &mut factorizations) {
            Ok(f) => {
                chosen = Some((sigma, f));
                break;
            }
            Err(e) => last_err = Some(e),
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

    let mut apply_m = |x: &[f64], y: &mut [f64]| m.spmv(x, y);
    let mut solve_shifted = |b: &[f64], out: &mut [f64]| {
        let x = f_mid.solve(b);
        out.copy_from_slice(&x);
    };

    // Restart-with-deflation until the certified count is met.
    let mut converged: Vec<(f64, Vec<f64>)> = Vec::new();
    let mut iters_used = 0usize;
    let mut restarts = 0usize;
    while converged.len() < expected {
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
/// row-major inputs, SPD M) via the companion linearization
/// `[[0, I], [−M⁻¹K, −M⁻¹C]]` and the fs-la dense complex eigensolver.
/// Returns eigenvalues sorted canonically by (re, im); eigenvectors are not
/// computed (the fs-la complex path is values-only — recorded boundary).
///
/// # Errors
/// [`ModalError::MassNotSpd`] if M fails Cholesky; [`ModalError::QuadraticEig`]
/// if the complex QR iteration fails to converge.
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
    let chol = cholesky(m, n).map_err(|_| ModalError::MassNotSpd { negative: 0 })?;
    let minv_times = |b: &[f64]| -> Vec<f64> {
        let mut x = b.to_vec();
        chol.solve(&mut x);
        x
    };
    // Companion A (2n×2n complex, row-major).
    let dim = 2 * n;
    let mut a = vec![C64::ZERO; dim * dim];
    for i in 0..n {
        a[i * dim + (n + i)] = C64::new(1.0, 0.0);
    }
    // Rows n..2n: [−M⁻¹K, −M⁻¹C] built column-block by column-block.
    let mut colbuf = vec![0.0f64; n];
    for cidx in 0..n {
        for (i, slot) in colbuf.iter_mut().enumerate() {
            *slot = k[i * n + cidx];
        }
        let x = minv_times(&colbuf);
        for (i, xv) in x.iter().enumerate() {
            a[(n + i) * dim + cidx] = C64::new(-xv, 0.0);
        }
        for (i, slot) in colbuf.iter_mut().enumerate() {
            *slot = c[i * n + cidx];
        }
        let x = minv_times(&colbuf);
        for (i, xv) in x.iter().enumerate() {
            a[(n + i) * dim + (n + cidx)] = C64::new(-xv, 0.0);
        }
    }
    complex_eig(&a, dim).map_err(|e| ModalError::QuadraticEig {
        detail: format!("{e:?}"),
    })
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
}
