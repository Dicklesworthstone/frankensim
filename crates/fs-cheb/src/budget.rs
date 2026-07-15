//! Typed budgets, admission preflight, and cancellable entry points
//! (bead frankensim-sj31i.55, slice 1).
//!
//! The classic fs-cheb entry points describe work as bounded merely
//! because a scalar maximum exists; a caller can select enormous
//! degrees or matrix dimensions and drive unadmitted allocation and
//! O(n³) work with no `Cx`, no memory admission, and no cancellation.
//! This module adds the RESOURCE CONTRACT: a [`ChebBudget`] declares
//! typed caps, [`ChebAdmission`] derives the WORST-CASE samples,
//! coefficients, work operations, and peak temporary bytes with
//! CHECKED `u128` formulas BEFORE any allocation (a saturating or
//! overflowing size refuses instead of iterating), and the budgeted
//! entry points thread an explicit [`Cx`] with cancellation polls at
//! bounded round/sweep boundaries. Terminal states are EXPLICIT:
//! `Complete` carries a [`WorkReceipt`]; `Cancelled` carries the spent
//! receipt plus a RESUME point where scientifically meaningful;
//! refusals are typed [`ChebError`]s — never a panic from size
//! arithmetic.
//!
//! Slice-1 coverage: the adaptive [`Cheb1`] constructor (resumable),
//! the Dirichlet collocation eigensolve (partial eigenvalues on
//! cancellation), and the fixed-grid root scan (no partial claim — an
//! incomplete scan is not a root set). The classic panicking APIs are
//! unchanged this slice; `cheb2`/`colleague`/`fourier`/
//! `orr_sommerfeld` budgeting is recorded follow-up scope in
//! CONTRACT.md and the bead.

use crate::{Cheb1, PLATEAU_REL, affine_from_reference, diff_matrix, fma};
use fs_exec::Cx;

/// Version of the budget/admission schema: bump when a cap default or
/// a worst-case formula changes meaning.
pub const CHEB_BUDGET_SCHEMA_VERSION: u32 = 2;

const F64_BYTES: u128 = core::mem::size_of::<f64>() as u128;
const SAMPLE_POLL_STRIDE: usize = 64;
const ROOT_SCAN_FACTOR: u128 = 8;
const ROOT_EVALS_PER_CELL: u128 = 90;
const ROOT_OPS_PER_COEFFICIENT_EVAL: u128 = 16;
const ADAPTIVE_TEMP_BYTES_PER_GRID_POINT: u128 = 64;
const ADAPTIVE_OPS_PER_GRID_STAGE: u128 = 4096;
const FD_SURROGATE_DIM: u128 = 64;
// fs-la's current blocked LU allocates 163_840 f64 packing scalars per
// GEMM call. Keep more than 3x headroom in this schema so the admission
// remains conservative across tuning-only panel changes.
const LU_PACK_WORKSPACE_SCALARS: u128 = 1 << 19;

/// Typed caps for fs-cheb construction and transform work. Construct
/// via [`ChebBudget::default`] and override fields; non-exhaustive so
/// new axes can join without breaking callers.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ChebBudget {
    /// Maximum retained coefficients in one function object.
    pub max_coefficients: usize,
    /// Maximum TOTAL samples across all adaptive rounds.
    pub max_samples: usize,
    /// Maximum collocation dimension (`n + 1` for the Lobatto grid).
    pub max_eigen_dim: usize,
    /// Maximum abstract work operations (checked complexity formulas).
    pub max_work_ops: u64,
    /// Maximum peak temporary bytes (checked size formulas).
    pub max_temp_bytes: u64,
}

impl ChebBudget {
    /// The v1 cap schedule.
    pub const V1: ChebBudget = ChebBudget {
        max_coefficients: 1 << 20,
        max_samples: 1 << 22,
        max_eigen_dim: 4096,
        max_work_ops: 1 << 42,
        max_temp_bytes: 1 << 32,
    };
}

impl Default for ChebBudget {
    fn default() -> Self {
        ChebBudget::V1
    }
}

/// Typed refusals and terminal diagnoses. Every size/complexity
/// refusal happens BEFORE allocation or function evaluation.
#[derive(Debug, Clone, PartialEq)]
pub enum ChebError {
    /// The domain is not finite with `a < b`.
    Domain {
        /// Offending endpoints (bit patterns preserved via Display).
        a: f64,
        /// Right endpoint.
        b: f64,
    },
    /// The problem shape is structurally invalid (e.g. eigensolve needs
    /// at least one interior point).
    Shape {
        /// What is wrong.
        what: &'static str,
    },
    /// A checked worst-case formula exceeds its declared cap.
    CapExceeded {
        /// Which cap.
        what: &'static str,
        /// Worst-case need (exact, no saturation).
        need: u128,
        /// The declared cap.
        cap: u128,
    },
    /// A checked size/complexity formula leaves the representable
    /// domain — refused rather than saturated and iterated.
    Overflow {
        /// Which formula.
        what: &'static str,
    },
    /// The plateau was not reached within the admitted degree cap
    /// (non-smooth or too oscillatory input — the classic API panics
    /// here; the budgeted API refuses).
    Unresolved {
        /// The degree cap that failed to resolve the function.
        max_degree: usize,
    },
    /// A sample or transform coefficient became non-finite.
    NonFinite {
        /// Where.
        what: &'static str,
    },
    /// A numerical precondition failed inside admitted work (e.g. a
    /// singular shifted operator).
    Numerical {
        /// What failed.
        what: &'static str,
    },
    /// Cancellation drained at a bounded boundary and this operation
    /// has no acceptance-capable partial result to return.
    Cancelled,
}

impl core::fmt::Display for ChebError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            ChebError::Domain { a, b } => write!(
                f,
                "domain must be finite with a < b (got [{a}, {b}]; bits {:016X}/{:016X})",
                a.to_bits(),
                b.to_bits()
            ),
            ChebError::Shape { what } => write!(f, "invalid problem shape: {what}"),
            ChebError::CapExceeded { what, need, cap } => write!(
                f,
                "budget refused before allocation: {what} needs {need}, cap {cap}; \
                 raise the explicit ChebBudget or shrink the request"
            ),
            ChebError::Overflow { what } => write!(
                f,
                "checked size formula `{what}` leaves the representable domain; \
                 the request is impossible, not merely expensive"
            ),
            ChebError::Unresolved { max_degree } => write!(
                f,
                "function not resolved at degree {max_degree} (non-smooth or too \
                 oscillatory; raise max_degree or split the domain)"
            ),
            ChebError::NonFinite { what } => {
                write!(f, "{what} must be representable as finite f64")
            }
            ChebError::Numerical { what } => write!(f, "numerical precondition failed: {what}"),
            ChebError::Cancelled => write!(
                f,
                "cancelled at a bounded boundary; no acceptance-capable partial \
                 result exists for this operation"
            ),
        }
    }
}

impl std::error::Error for ChebError {}

/// Deterministic diagnostics for one admitted run (fixed traversal
/// order, no time source). This is not authenticated acceptance
/// evidence; its public fields deliberately remain ordinary data.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct WorkReceipt {
    /// Budget schema in force.
    pub schema_version: u32,
    /// Adaptive rounds / eigen shifts / scan chunks completed.
    pub rounds_completed: u32,
    /// Function samples or inverse-iteration sweeps actually spent.
    pub samples_spent: usize,
    /// Worst-case operations ADMITTED for the run (the preflight bound).
    pub ops_admitted: u64,
}

/// Terminal state of a budgeted adaptive construction.
#[derive(Debug, Clone, PartialEq)]
pub enum BuildRun {
    /// The plateau was reached; the function is complete.
    Complete {
        /// The constructed function.
        function: Cheb1,
        /// What the run spent.
        receipt: WorkReceipt,
    },
    /// Cancellation drained at a round boundary. `resume_from` is the
    /// grid size the next call should start at
    /// ([`try_build_budgeted`]'s `start_degree`) — resumption is
    /// deterministic and bitwise-equivalent to an uncancelled run.
    Cancelled {
        /// Grid size to resume at.
        resume_from: usize,
        /// What the run spent before draining.
        receipt: WorkReceipt,
    },
}

/// Terminal state of a budgeted eigensolve.
#[derive(Debug, Clone, PartialEq)]
pub enum EigsRun {
    /// All requested fixed-sweep eigenvalue estimates were produced.
    /// This state does not certify convergence or residual quality.
    Complete {
        /// The eigenvalues, smallest-first.
        eigs: Vec<f64>,
        /// What the run spent.
        receipt: WorkReceipt,
    },
    /// Cancellation drained at a shift/sweep boundary. The retained
    /// prefix contains completed fixed-sweep estimates only; it is
    /// diagnostic output, not a convergence certificate.
    Cancelled {
        /// Fixed-sweep estimates completed before the drain (may be empty).
        partial_eigs: Vec<f64>,
        /// What the run spent before draining.
        receipt: WorkReceipt,
    },
}

/// Evidence that a request's worst case fits the budget: the checked
/// preflight numbers, derived BEFORE any allocation. Sealed fields —
/// holding one means the formulas actually ran.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ChebAdmission {
    schema_version: u32,
    samples_admitted: usize,
    coefficients_admitted: usize,
    ops_admitted: u64,
    temp_bytes_admitted: u64,
}

impl ChebAdmission {
    /// Budget schema the preflight ran under.
    #[must_use]
    pub fn schema_version(&self) -> u32 {
        self.schema_version
    }

    /// Worst-case total samples admitted.
    #[must_use]
    pub fn samples_admitted(&self) -> usize {
        self.samples_admitted
    }

    /// Worst-case retained coefficients admitted.
    #[must_use]
    pub fn coefficients_admitted(&self) -> usize {
        self.coefficients_admitted
    }

    /// Worst-case abstract operations admitted.
    #[must_use]
    pub fn ops_admitted(&self) -> u64 {
        self.ops_admitted
    }

    /// Worst-case peak temporary bytes admitted.
    #[must_use]
    pub fn temp_bytes_admitted(&self) -> u64 {
        self.temp_bytes_admitted
    }
}

fn cap_check(what: &'static str, need: u128, cap: u128) -> Result<(), ChebError> {
    if need > cap {
        return Err(ChebError::CapExceeded { what, need, cap });
    }
    Ok(())
}

fn checked_add(what: &'static str, left: u128, right: u128) -> Result<u128, ChebError> {
    left.checked_add(right).ok_or(ChebError::Overflow { what })
}

fn checked_mul(what: &'static str, left: u128, right: u128) -> Result<u128, ChebError> {
    left.checked_mul(right).ok_or(ChebError::Overflow { what })
}

fn checked_sum(what: &'static str, terms: &[u128]) -> Result<u128, ChebError> {
    terms
        .iter()
        .try_fold(0u128, |sum, &term| checked_add(what, sum, term))
}

fn admitted_u64(what: &'static str, need: u128) -> Result<u64, ChebError> {
    u64::try_from(need).map_err(|_| ChebError::Overflow { what })
}

fn domain_check(a: f64, b: f64) -> Result<(), ChebError> {
    if a.is_finite() && b.is_finite() && a < b {
        Ok(())
    } else {
        Err(ChebError::Domain { a, b })
    }
}

fn checkpoint(cx: &Cx<'_>) -> Result<(), ChebError> {
    cx.checkpoint().map_err(|_| ChebError::Cancelled)
}

/// Worst-case preflight for the adaptive constructor. Grid sizes
/// double from `start` through the largest power of two not exceeding
/// `max_degree`. The temporary envelope covers the sampled `f64`
/// values plus DCT-II complex data/scratch, FFT twiddles, output, and
/// six-step sub-plan workspace. All arithmetic is checked `u128`.
///
/// # Errors
/// [`ChebError`] — domain, overflow, or cap refusals.
pub fn admit_adaptive_build(
    a: f64,
    b: f64,
    max_degree: usize,
    start: usize,
    budget: &ChebBudget,
) -> Result<ChebAdmission, ChebError> {
    domain_check(a, b)?;
    let degree_cap = max_degree.max(16);
    let final_grid = 1usize
        .checked_shl(degree_cap.ilog2())
        .ok_or(ChebError::Overflow {
            what: "largest adaptive grid",
        })?;
    if start < 16 || !start.is_power_of_two() {
        return Err(ChebError::Shape {
            what: "resume grid must be a power of two and at least 16",
        });
    }
    if start > final_grid {
        return Err(ChebError::Shape {
            what: "resume grid exceeds the admitted degree cap",
        });
    }
    cap_check(
        "retained coefficients",
        final_grid as u128,
        budget.max_coefficients as u128,
    )?;
    // Exact geometric sum start + 2*start + ... + final_grid.
    let twice_final = checked_mul("adaptive sample sum", 2, final_grid as u128)?;
    let total_samples = twice_final
        .checked_sub(start as u128)
        .ok_or(ChebError::Overflow {
            what: "adaptive sample sum",
        })?;
    cap_check(
        "adaptive samples",
        total_samples,
        budget.max_samples as u128,
    )?;
    let samples_admitted = usize::try_from(total_samples).map_err(|_| ChebError::Overflow {
        what: "total adaptive samples",
    })?;
    // A deliberately conservative abstract-work envelope: every grid
    // point receives a fixed allowance for mapping/sampling, strict trig,
    // FFT-plan twiddles, radix stages, scaling, and the plateau scan. The
    // geometric factor covers all earlier rounds.
    let log2 = u128::from(final_grid.ilog2().max(1));
    let stage_count = checked_add("adaptive work", log2, 1)?;
    let ops = checked_mul(
        "adaptive work",
        checked_mul("adaptive work", twice_final, stage_count)?,
        ADAPTIVE_OPS_PER_GRID_STAGE,
    )?;
    cap_check("adaptive work", ops, u128::from(budget.max_work_ops))?;
    let temp_bytes = checked_mul(
        "adaptive temporary bytes",
        final_grid as u128,
        ADAPTIVE_TEMP_BYTES_PER_GRID_POINT,
    )?;
    cap_check(
        "adaptive temporary bytes",
        temp_bytes,
        u128::from(budget.max_temp_bytes),
    )?;
    Ok(ChebAdmission {
        schema_version: CHEB_BUDGET_SCHEMA_VERSION,
        samples_admitted,
        coefficients_admitted: final_grid,
        ops_admitted: admitted_u64("adaptive work", ops)?,
        temp_bytes_admitted: admitted_u64("adaptive temporary bytes", temp_bytes)?,
    })
}

/// Worst-case preflight for the Dirichlet collocation eigensolve:
/// dense (n+1)² differentiation matrices, an O(m³) matrix square, and
/// per-shift LU + 100 inverse-power sweeps on the (n−1)² interior
/// block. Exact `u128` formulas; `usize::MAX`-shaped inputs refuse
/// before any allocation.
///
/// # Errors
/// [`ChebError`] — shape, overflow, or cap refusals.
#[allow(clippy::too_many_lines)] // every downstream allocation/work term stays explicit and auditable
pub fn admit_dirichlet_eigs(
    n: usize,
    k: usize,
    budget: &ChebBudget,
) -> Result<ChebAdmission, ChebError> {
    if n < 2 {
        return Err(ChebError::Shape {
            what: "the Dirichlet eigensolve needs n >= 2 (at least one interior point)",
        });
    }
    if k == 0 {
        return Err(ChebError::Shape {
            what: "requesting zero eigenvalues is a caller bug, not work",
        });
    }
    let m = n.checked_add(1).ok_or(ChebError::Overflow {
        what: "collocation dimension n + 1",
    })?;
    cap_check(
        "collocation dimension",
        m as u128,
        budget.max_eigen_dim as u128,
    )?;
    if k > n - 1 {
        return Err(ChebError::Shape {
            what: "cannot request more eigenvalues than interior points",
        });
    }
    if k > 64 {
        return Err(ChebError::Shape {
            what: "the fixed 64-point FD surrogate supplies at most 64 shifts \
                   (the classic API silently shorts here; the budgeted API refuses)",
        });
    }
    let m_u128 = m as u128;
    let m2 = checked_mul("collocation dimension squared", m_u128, m_u128)?;
    let ni = (n - 1) as u128;
    let ni2 = checked_mul("interior dimension squared", ni, ni)?;
    let ni3 = checked_mul("interior dimension cubed", ni2, ni)?;

    // Persistent D/D², interior/shifted matrices, LU storage and its
    // worst blocked-update matrix, fixed Jacobi input/copies/eigenvectors,
    // iteration vectors, pivots, result storage, and fs-la GEMM packs.
    // Several lifetimes do not overlap; summing them remains deliberately
    // conservative and stable under harmless scope changes downstream.
    let two_m2 = checked_mul("eigensolve temporary scalars", 2, m2)?;
    let five_ni2 = checked_mul("eigensolve temporary scalars", 5, ni2)?;
    let linear_ni = checked_mul("eigensolve temporary scalars", 128, ni)?;
    let fixed_jacobi = checked_mul(
        "eigensolve temporary scalars",
        8,
        checked_mul(
            "eigensolve temporary scalars",
            FD_SURROGATE_DIM,
            FD_SURROGATE_DIM,
        )?,
    )?;
    let lu_pack_workspace = if ni > 32 {
        LU_PACK_WORKSPACE_SCALARS
    } else {
        0
    };
    let result_scalars = checked_mul("eigensolve temporary scalars", 4, k as u128)?;
    let temp_scalars = checked_sum(
        "eigensolve temporary scalars",
        &[
            two_m2,
            five_ni2,
            linear_ni,
            fixed_jacobi,
            lu_pack_workspace,
            result_scalars,
            1024,
        ],
    )?;
    let temp_bytes = checked_mul("eigensolve temporary bytes", temp_scalars, F64_BYTES)?;
    cap_check(
        "eigensolve temporary bytes",
        temp_bytes,
        u128::from(budget.max_temp_bytes),
    )?;
    // Conservative scalar-operation envelope: matrix construction/square,
    // all 60 cyclic-Jacobi sweeps (including three dense row/column/vector
    // updates per rotation), blocked LU, 100 two-triangular solves with
    // normalization, and the final Rayleigh product for every shift.
    let m3 = checked_mul("collocation dimension cubed", m2, m_u128)?;
    let matrix_work = checked_sum(
        "eigensolve work",
        &[
            checked_mul("eigensolve work", 4, m3)?,
            checked_mul("eigensolve work", 32, m2)?,
        ],
    )?;
    let fd_rotations = FD_SURROGATE_DIM * (FD_SURROGATE_DIM - 1) / 2;
    let jacobi_work = checked_mul(
        "eigensolve work",
        checked_mul(
            "eigensolve work",
            checked_mul("eigensolve work", 60, fd_rotations)?,
            FD_SURROGATE_DIM,
        )?,
        64,
    )?;
    let per_shift_work = checked_sum(
        "eigensolve work",
        &[
            checked_mul("eigensolve work", 8, ni3)?,
            checked_mul("eigensolve work", 512, ni2)?,
            checked_mul("eigensolve work", 128, ni)?,
        ],
    )?;
    let ops = checked_sum(
        "eigensolve work",
        &[
            matrix_work,
            jacobi_work,
            checked_mul("eigensolve work", k as u128, per_shift_work)?,
        ],
    )?;
    cap_check("eigensolve work", ops, u128::from(budget.max_work_ops))?;
    Ok(ChebAdmission {
        schema_version: CHEB_BUDGET_SCHEMA_VERSION,
        samples_admitted: 0,
        coefficients_admitted: 0,
        ops_admitted: admitted_u64("eigensolve work", ops)?,
        temp_bytes_admitted: admitted_u64("eigensolve temporary bytes", temp_bytes)?,
    })
}

/// Worst-case preflight for the fixed-grid root scan: `8·len` scan
/// cells (minimum 64), up to 89 safeguarded-refinement evaluations in
/// every cell, checked normalization/differentiation, and retained root
/// candidates. The bound does not assume exact arithmetic limits the
/// number of observed floating-point sign changes.
///
/// # Errors
/// [`ChebError`] — overflow or cap refusals.
pub fn admit_root_scan(coeff_len: usize, budget: &ChebBudget) -> Result<ChebAdmission, ChebError> {
    cap_check(
        "root-scan coefficients",
        coeff_len as u128,
        budget.max_coefficients as u128,
    )?;
    let samples_u128 = checked_mul(
        "root-scan sample count",
        coeff_len as u128,
        ROOT_SCAN_FACTOR,
    )?
    .max(64);
    let samples = usize::try_from(samples_u128).map_err(|_| ChebError::Overflow {
        what: "root-scan sample count",
    })?;
    cap_check(
        "root-scan samples",
        samples as u128,
        budget.max_samples as u128,
    )?;
    let evals = checked_add(
        "root-scan evaluation count",
        checked_mul(
            "root-scan evaluation count",
            samples_u128,
            ROOT_EVALS_PER_CELL,
        )?,
        2,
    )?;
    let eval_work = checked_mul(
        "root-scan work",
        checked_mul("root-scan work", evals, coeff_len as u128)?,
        ROOT_OPS_PER_COEFFICIENT_EVAL,
    )?;
    let linear_work = checked_mul(
        "root-scan work",
        128,
        checked_add("root-scan work", coeff_len as u128, samples_u128)?,
    )?;
    let ops = checked_add("root-scan work", eval_work, linear_work)?;
    cap_check("root-scan work", ops, u128::from(budget.max_work_ops))?;
    // Reference coefficients + derivative recurrence/output + retained
    // candidates: three coefficient-length f64 buffers cover every peak.
    let temp_bytes = checked_mul(
        "root-scan temporary bytes",
        checked_mul("root-scan temporary bytes", 3, coeff_len as u128)?,
        F64_BYTES,
    )?;
    cap_check(
        "root-scan temporary bytes",
        temp_bytes,
        u128::from(budget.max_temp_bytes),
    )?;
    Ok(ChebAdmission {
        schema_version: CHEB_BUDGET_SCHEMA_VERSION,
        samples_admitted: samples,
        coefficients_admitted: coeff_len,
        ops_admitted: admitted_u64("root-scan work", ops)?,
        temp_bytes_admitted: admitted_u64("root-scan temporary bytes", temp_bytes)?,
    })
}

/// Fallible mirror of the private sampler: identical sample sequence
/// (bitwise), refusing a non-finite sample instead of panicking.
fn sample_first_kind_checked<F: Fn(f64) -> f64>(
    f: &F,
    a: f64,
    b: f64,
    n: usize,
    cx: &Cx<'_>,
    samples_spent: &mut usize,
) -> Result<Vec<f64>, ChebError> {
    checkpoint(cx)?;
    let mut vals = Vec::with_capacity(n);
    for k in 0..n {
        checkpoint(cx)?;
        let theta = std::f64::consts::PI * (k as f64 + 0.5) / (n as f64);
        let t = fs_math::det::cos(theta);
        let x = affine_from_reference(t, a, b);
        let y = f(x);
        *samples_spent = samples_spent.checked_add(1).ok_or(ChebError::Overflow {
            what: "adaptive samples spent",
        })?;
        checkpoint(cx)?;
        if !y.is_finite() {
            return Err(ChebError::NonFinite {
                what: "Cheb1 sample",
            });
        }
        vals.push(y);
    }
    Ok(vals)
}

fn coeffs_at_checked<F: Fn(f64) -> f64>(
    f: &F,
    a: f64,
    b: f64,
    n: usize,
    cx: &Cx<'_>,
    samples_spent: &mut usize,
) -> Result<Vec<f64>, ChebError> {
    let vals = sample_first_kind_checked(f, a, b, n, cx, samples_spent)?;
    checkpoint(cx)?;
    let mut c = fs_fft::dct2(&vals);
    checkpoint(cx)?;
    let scale = 2.0 / n as f64;
    for (index, v) in c.iter_mut().enumerate() {
        if index.is_multiple_of(SAMPLE_POLL_STRIDE) {
            checkpoint(cx)?;
        }
        *v *= scale;
    }
    checkpoint(cx)?;
    if !c.iter().all(|coefficient| coefficient.is_finite()) {
        return Err(ChebError::NonFinite {
            what: "Chebyshev transform coefficient",
        });
    }
    Ok(c)
}

/// Budgeted, cancellable, RESUMABLE adaptive construction. Semantics
/// are bitwise-identical to [`Cheb1::build`] on the happy path (same
/// sample sequence, same transform, same plateau rule); refusals are
/// typed instead of panics; cancellation drains at a round boundary
/// and returns the resume grid size. Pass a prior run's `resume_from`
/// as `start_degree` to continue deterministically.
///
/// # Errors
/// [`ChebError`] — domain/overflow/cap refusals before any sampling,
/// [`ChebError::NonFinite`] on unrepresentable samples,
/// [`ChebError::Unresolved`] when the admitted degree cap cannot
/// resolve the function.
pub fn try_build_budgeted<F: Fn(f64) -> f64>(
    f: &F,
    a: f64,
    b: f64,
    max_degree: usize,
    start_degree: Option<usize>,
    budget: &ChebBudget,
    cx: &Cx<'_>,
) -> Result<BuildRun, ChebError> {
    domain_check(a, b)?;
    let start = start_degree
        .unwrap_or(16)
        .max(16)
        .checked_next_power_of_two()
        .ok_or(ChebError::Overflow {
            what: "resume grid next power of two",
        })?;
    let admission = admit_adaptive_build(a, b, max_degree, start, budget)?;
    let degree_cap = max_degree.max(16);
    let mut n = start;
    let mut rounds_completed = 0u32;
    let mut samples_spent = 0usize;
    let cancelled =
        |resume_from: usize, rounds_completed: u32, samples_spent: usize| BuildRun::Cancelled {
            resume_from,
            receipt: WorkReceipt {
                schema_version: CHEB_BUDGET_SCHEMA_VERSION,
                rounds_completed,
                samples_spent,
                ops_admitted: admission.ops_admitted(),
            },
        };
    loop {
        // Bounded tile boundary: one adaptive round.
        if cx.checkpoint().is_err() {
            return Ok(cancelled(n, rounds_completed, samples_spent));
        }
        let mut coeffs = match coeffs_at_checked(f, a, b, n, cx, &mut samples_spent) {
            Ok(coeffs) => coeffs,
            Err(ChebError::Cancelled) => {
                return Ok(cancelled(n, rounds_completed, samples_spent));
            }
            Err(error) => return Err(error),
        };
        rounds_completed += 1;
        let maxc = coeffs
            .iter()
            .fold(0.0f64, |m, &c| m.max(c.abs()))
            .max(f64::MIN_POSITIVE);
        let tail = &coeffs[3 * n / 4..];
        if tail.iter().all(|&c| c.abs() <= PLATEAU_REL * maxc) {
            let keep = coeffs
                .iter()
                .rposition(|&c| c.abs() > PLATEAU_REL * maxc)
                .map_or(1, |p| p + 1);
            if cx.checkpoint().is_err() {
                return Ok(cancelled(n, rounds_completed, samples_spent));
            }
            coeffs.truncate(keep);
            if cx.checkpoint().is_err() {
                return Ok(cancelled(n, rounds_completed, samples_spent));
            }
            return Ok(BuildRun::Complete {
                function: Cheb1 { a, b, coeffs },
                receipt: WorkReceipt {
                    schema_version: CHEB_BUDGET_SCHEMA_VERSION,
                    rounds_completed,
                    samples_spent,
                    ops_admitted: admission.ops_admitted(),
                },
            });
        }
        n = n.checked_mul(2).ok_or(ChebError::Overflow {
            what: "adaptive grid doubling",
        })?;
        if n > degree_cap {
            return Err(ChebError::Unresolved {
                max_degree: degree_cap,
            });
        }
    }
}

/// Budgeted, cancellable Dirichlet collocation eigensolve. Semantics
/// match [`crate::dirichlet_laplace_eigs`] bitwise on the happy path;
/// impossible shapes refuse before allocation; cancellation drains at
/// every public allocation/opaque-kernel boundary and every ten inverse
/// sweeps. A cancelled run retains only fully completed fixed-sweep
/// estimates; those values do not carry convergence authority.
///
/// # Errors
/// [`ChebError`] — shape/overflow/cap refusals, or
/// [`ChebError::Numerical`] on a singular shifted operator.
#[allow(clippy::too_many_lines)] // mirrors the classic solver with polls + receipts inline
pub fn dirichlet_laplace_eigs_budgeted(
    n: usize,
    k: usize,
    budget: &ChebBudget,
    cx: &Cx<'_>,
) -> Result<EigsRun, ChebError> {
    let admission = admit_dirichlet_eigs(n, k, budget)?;
    let cancelled = |partial_eigs: Vec<f64>, rounds: u32, sweeps: usize| EigsRun::Cancelled {
        partial_eigs,
        receipt: WorkReceipt {
            schema_version: CHEB_BUDGET_SCHEMA_VERSION,
            rounds_completed: rounds,
            samples_spent: sweeps,
            ops_admitted: admission.ops_admitted(),
        },
    };
    if cx.checkpoint().is_err() {
        return Ok(cancelled(Vec::new(), 0, 0));
    }
    let m = n + 1;
    if cx.checkpoint().is_err() {
        return Ok(cancelled(Vec::new(), 0, 0));
    }
    let d = diff_matrix(n);
    if cx.checkpoint().is_err() {
        return Ok(cancelled(Vec::new(), 0, 0));
    }
    let mut d2 = vec![0.0f64; m * m];
    if cx.checkpoint().is_err() {
        return Ok(cancelled(Vec::new(), 0, 0));
    }
    fma::dsq_into_dispatch(&d, m, &mut d2);
    if cx.checkpoint().is_err() {
        return Ok(cancelled(Vec::new(), 0, 0));
    }
    let ni = n - 1;
    let mut a = vec![0.0f64; ni * ni];
    for i in 0..ni {
        if i.is_multiple_of(8) && cx.checkpoint().is_err() {
            return Ok(cancelled(Vec::new(), 0, 0));
        }
        for j in 0..ni {
            a[i * ni + j] = -d2[(i + 1) * m + (j + 1)];
        }
    }
    let nf = 64usize;
    let h = 2.0 / (nf as f64 + 1.0);
    if cx.checkpoint().is_err() {
        return Ok(cancelled(Vec::new(), 0, 0));
    }
    let mut fd = vec![0.0f64; nf * nf];
    for i in 0..nf {
        if i.is_multiple_of(8) && cx.checkpoint().is_err() {
            return Ok(cancelled(Vec::new(), 0, 0));
        }
        fd[i * nf + i] = 2.0 / (h * h);
        if i + 1 < nf {
            fd[i * nf + i + 1] = -1.0 / (h * h);
            fd[(i + 1) * nf + i] = -1.0 / (h * h);
        }
    }
    if cx.checkpoint().is_err() {
        return Ok(cancelled(Vec::new(), 0, 0));
    }
    let (fd_eigs, _) = fs_la::eigen::jacobi_eigh(&fd, nf);
    if cx.checkpoint().is_err() {
        return Ok(cancelled(Vec::new(), 0, 0));
    }
    let mut eigs = Vec::with_capacity(k);
    if cx.checkpoint().is_err() {
        return Ok(cancelled(eigs, 0, 0));
    }
    let mut shifted = vec![0.0f64; a.len()];
    let mut sweeps_total = 0usize;
    let receipt = |rounds: u32, sweeps: usize| WorkReceipt {
        schema_version: CHEB_BUDGET_SCHEMA_VERSION,
        rounds_completed: rounds,
        samples_spent: sweeps,
        ops_admitted: admission.ops_admitted(),
    };
    for (shift_index, &fd_est) in fd_eigs.iter().take(k).enumerate() {
        // Shift boundary poll: only completed fixed-sweep estimates survive.
        if cx.checkpoint().is_err() {
            return Ok(cancelled(eigs, shift_index as u32, sweeps_total));
        }
        let mu = fd_est * 0.95;
        shifted.copy_from_slice(&a);
        for i in 0..ni {
            shifted[i * ni + i] -= mu;
        }
        if cx.checkpoint().is_err() {
            return Ok(cancelled(eigs, shift_index as u32, sweeps_total));
        }
        let lu = fs_la::factor::lu(&shifted, ni).map_err(|_| ChebError::Numerical {
            what: "shifted collocation operator is singular",
        })?;
        if cx.checkpoint().is_err() {
            return Ok(cancelled(eigs, shift_index as u32, sweeps_total));
        }
        let mut v: Vec<f64> = (0..ni)
            .map(|i| 1.0 + 0.25 * (((i * 7 + 3) % 11) as f64))
            .collect();
        for sweep in 0..100 {
            if sweep % 10 == 0 && cx.checkpoint().is_err() {
                return Ok(cancelled(eigs, shift_index as u32, sweeps_total));
            }
            let nrm = v.iter().map(|x| x * x).sum::<f64>().sqrt();
            if !nrm.is_finite() || nrm == 0.0 {
                return Err(ChebError::Numerical {
                    what: "inverse-iteration vector has no finite non-zero norm",
                });
            }
            for x in &mut v {
                *x /= nrm;
            }
            lu.solve(&mut v);
            if !v.iter().all(|value| value.is_finite()) {
                return Err(ChebError::Numerical {
                    what: "inverse iteration produced a non-finite vector",
                });
            }
            sweeps_total += 1;
            if (sweep + 1) % 10 == 0 && cx.checkpoint().is_err() {
                return Ok(cancelled(eigs, shift_index as u32, sweeps_total));
            }
        }
        let nrm2: f64 = v.iter().map(|x| x * x).sum();
        if !nrm2.is_finite() || nrm2 == 0.0 {
            return Err(ChebError::Numerical {
                what: "Rayleigh vector has no finite non-zero squared norm",
            });
        }
        if cx.checkpoint().is_err() {
            return Ok(cancelled(eigs, shift_index as u32, sweeps_total));
        }
        let mut av = vec![0.0f64; ni];
        fma::matvec_into_dispatch(&a, &v, ni, &mut av);
        let estimate = v.iter().zip(&av).map(|(x, y)| x * y).sum::<f64>() / nrm2;
        if !estimate.is_finite() {
            return Err(ChebError::Numerical {
                what: "Rayleigh estimate is non-finite",
            });
        }
        if cx.checkpoint().is_err() {
            return Ok(cancelled(eigs, shift_index as u32, sweeps_total));
        }
        eigs.push(estimate);
    }
    if cx.checkpoint().is_err() {
        return Ok(cancelled(eigs, k as u32, sweeps_total));
    }
    Ok(EigsRun::Complete {
        eigs,
        receipt: receipt(k as u32, sweeps_total),
    })
}

fn normalize_root_coefficients(coeffs: &mut [f64], cx: &Cx<'_>) -> Result<(), ChebError> {
    checkpoint(cx)?;
    let cmax = coeffs
        .iter()
        .fold(0.0f64, |scale, coefficient| scale.max(coefficient.abs()));
    if !cmax.is_finite() || cmax == 0.0 {
        return Err(ChebError::Numerical {
            what: "root normalization requires finite non-zero coefficients",
        });
    }
    let scale = crate::normalization_power_of_two(cmax);
    for (index, coefficient) in coeffs.iter_mut().enumerate() {
        if index.is_multiple_of(SAMPLE_POLL_STRIDE) {
            checkpoint(cx)?;
        }
        let original = *coefficient;
        let normalized = original / scale;
        if !normalized.is_finite() || normalized * scale != original {
            return Err(ChebError::Numerical {
                what: "root normalization would lose coefficient information",
            });
        }
        *coefficient = normalized;
    }
    checkpoint(cx)
}

fn reference_derivative_checked(reference: &Cheb1, cx: &Cx<'_>) -> Result<Cheb1, ChebError> {
    checkpoint(cx)?;
    let n = reference.coeffs.len();
    if n == 1 {
        return Ok(Cheb1 {
            a: -1.0,
            b: 1.0,
            coeffs: vec![0.0],
        });
    }
    let mut recurrence = vec![0.0f64; n];
    for k in (1..n).rev() {
        if k.is_multiple_of(SAMPLE_POLL_STRIDE) {
            checkpoint(cx)?;
        }
        let above = if k + 2 < n { recurrence[k + 1] } else { 0.0 };
        recurrence[k - 1] = (2.0 * k as f64).mul_add(reference.coeffs[k], above);
        if !recurrence[k - 1].is_finite() {
            return Err(ChebError::Numerical {
                what: "reference derivative coefficient is non-finite",
            });
        }
    }
    checkpoint(cx)?;
    let coefficients: Vec<f64> = recurrence[..n - 1].to_vec();
    checkpoint(cx)?;
    Ok(Cheb1 {
        a: -1.0,
        b: 1.0,
        coeffs: coefficients,
    })
}

fn finite_root_eval(function: &Cheb1, t: f64) -> Result<f64, ChebError> {
    let value = function.eval_reference(t);
    if value.is_finite() {
        Ok(value)
    } else {
        Err(ChebError::NonFinite {
            what: "root scan evaluation",
        })
    }
}

fn ensure_resolvable_root(reference: &Cheb1, derivative: &Cheb1, t: f64) -> Result<(), ChebError> {
    let slope = finite_root_eval(derivative, t)?;
    let degree_scale = reference.degree().max(1) as f64;
    let slope_floor = 64.0 * 1.490_116_119_384_765_6e-8 * degree_scale;
    if slope.abs() > slope_floor {
        Ok(())
    } else {
        Err(ChebError::Numerical {
            what: "fixed-grid root scan cannot resolve a multiple or ill-conditioned root; use colleague/certified root evidence",
        })
    }
}

fn refine_root_checked(
    reference: &Cheb1,
    derivative: &Cheb1,
    mut lo: f64,
    mut hi: f64,
    cx: &Cx<'_>,
) -> Result<f64, ChebError> {
    for _ in 0..40 {
        checkpoint(cx)?;
        let mid = f64::midpoint(lo, hi);
        let value = finite_root_eval(reference, mid)?;
        if value == 0.0 {
            ensure_resolvable_root(reference, derivative, mid)?;
            return Ok(mid);
        }
        let lo_value = finite_root_eval(reference, lo)?;
        if lo_value != 0.0 && lo_value.is_sign_negative() != value.is_sign_negative() {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    let mut root = f64::midpoint(lo, hi);
    for _ in 0..4 {
        checkpoint(cx)?;
        let value = finite_root_eval(reference, root)?;
        let slope = finite_root_eval(derivative, root)?;
        if slope == 0.0 {
            break;
        }
        let step = value / slope;
        if !step.is_finite() {
            break;
        }
        let candidate = root - step;
        if !candidate.is_finite() || candidate < lo || candidate > hi {
            break;
        }
        root = candidate;
    }
    checkpoint(cx)?;
    ensure_resolvable_root(reference, derivative, root)?;
    Ok(root)
}

impl Cheb1 {
    /// Budgeted, cancellable root scan. The scan sequence, refinement,
    /// and conditioning rule match [`Cheb1::roots`] bitwise on the
    /// happy path; admission preflights the scan/refinement work;
    /// cancellation polls before allocation/evaluation, throughout
    /// normalization and refinement, every 64 scan cells, and before
    /// finalization. It refuses WITHOUT a partial result — an incomplete
    /// scan is not a root-set claim.
    ///
    /// # Errors
    /// [`ChebError`] — cap/overflow refusals before evaluation;
    /// [`ChebError::Numerical`] for the identically-zero polynomial,
    /// non-finite evaluations, or an unresolvable multiple root;
    /// [`ChebError::Cancelled`] on a drained cancellation.
    pub fn roots_budgeted(&self, budget: &ChebBudget, cx: &Cx<'_>) -> Result<Vec<f64>, ChebError> {
        let admission = admit_root_scan(self.coeffs.len(), budget)?;
        checkpoint(cx)?;
        if !self.coeffs.iter().any(|&coefficient| coefficient != 0.0) {
            return Err(ChebError::Numerical {
                what: "the identically zero polynomial has a continuum of roots",
            });
        }
        checkpoint(cx)?;
        let mut reference_coeffs = self.coeffs.clone();
        normalize_root_coefficients(&mut reference_coeffs, cx)?;
        let reference = Cheb1 {
            a: -1.0,
            b: 1.0,
            coeffs: reference_coeffs,
        };
        let derivative = reference_derivative_checked(&reference, cx)?;
        checkpoint(cx)?;
        let mut roots_t = Vec::with_capacity(reference.degree());
        let samples = admission.samples_admitted();
        let mut prev_t = -1.0;
        let mut prev_v = finite_root_eval(&reference, prev_t)?;
        for k in 1..=samples {
            // Bounded tile boundary: one poll per 64 scan cells.
            if k % 64 == 1 && cx.checkpoint().is_err() {
                return Err(ChebError::Cancelled);
            }
            let t = 2.0 * (k as f64) / (samples as f64) - 1.0;
            let value = finite_root_eval(&reference, t)?;
            if prev_v == 0.0 {
                ensure_resolvable_root(&reference, &derivative, prev_t)?;
                if roots_t.len() >= reference.degree() {
                    return Err(ChebError::Numerical {
                        what: "floating-point root candidates exceed the polynomial degree bound",
                    });
                }
                roots_t.push(prev_t);
            } else if value != 0.0 && prev_v.is_sign_negative() != value.is_sign_negative() {
                let root = refine_root_checked(&reference, &derivative, prev_t, t, cx)?;
                if roots_t.len() >= reference.degree() {
                    return Err(ChebError::Numerical {
                        what: "floating-point root candidates exceed the polynomial degree bound",
                    });
                }
                roots_t.push(root);
            }
            prev_t = t;
            prev_v = value;
        }
        if prev_v == 0.0 {
            ensure_resolvable_root(&reference, &derivative, prev_t)?;
            if roots_t.len() >= reference.degree() {
                return Err(ChebError::Numerical {
                    what: "floating-point root candidates exceed the polynomial degree bound",
                });
            }
            roots_t.push(prev_t);
        }
        checkpoint(cx)?;
        for (index, root) in roots_t.iter_mut().enumerate() {
            if index.is_multiple_of(SAMPLE_POLL_STRIDE) {
                checkpoint(cx)?;
            }
            *root = affine_from_reference(*root, self.a, self.b);
        }
        checkpoint(cx)?;
        Ok(roots_t)
    }
}
