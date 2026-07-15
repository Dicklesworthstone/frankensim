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
pub const CHEB_BUDGET_SCHEMA_VERSION: u32 = 1;

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

/// What one admitted, budgeted run actually spent — deterministic and
/// replayable (fixed traversal order, no time source).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct WorkReceipt {
    /// Budget schema in force.
    pub schema_version: u32,
    /// Adaptive rounds / eigen shifts / scan chunks completed.
    pub rounds_completed: u32,
    /// Samples or matrix evaluations actually spent.
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
    /// All requested eigenvalues converged.
    Complete {
        /// The eigenvalues, smallest-first.
        eigs: Vec<f64>,
        /// What the run spent.
        receipt: WorkReceipt,
    },
    /// Cancellation drained at a shift/sweep boundary. The prefix of
    /// converged eigenvalues IS scientifically meaningful (each shift
    /// converges independently) and is retained; the remainder was
    /// never computed.
    Cancelled {
        /// Eigenvalues converged before the drain (may be empty).
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

/// Worst-case preflight for the adaptive constructor: grid sizes
/// double from `start` to the power-of-two degree cap, so total
/// samples are bounded by twice the final grid; per-round transform
/// work is O(n log n); peak temporaries are two f64 buffers of the
/// final grid. All formulas are exact `u128` arithmetic — an
/// unrepresentable request refuses as [`ChebError::Overflow`].
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
    let final_grid = degree_cap
        .checked_next_power_of_two()
        .ok_or(ChebError::Overflow {
            what: "degree cap next power of two",
        })?;
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
    // 16 + 32 + ... + final_grid < 2 * final_grid.
    let total_samples = 2u128 * final_grid as u128;
    cap_check(
        "adaptive samples",
        total_samples,
        budget.max_samples as u128,
    )?;
    let samples_admitted = usize::try_from(total_samples).map_err(|_| ChebError::Overflow {
        what: "total adaptive samples",
    })?;
    // Per round: one DCT-II (~5 n log2 n) plus the plateau scan (n).
    let log2 = u128::from(final_grid.ilog2().max(1));
    let ops = 2u128 * (5 * final_grid as u128 * log2 + final_grid as u128);
    cap_check("adaptive work", ops, u128::from(budget.max_work_ops))?;
    // Peak temporaries: samples + coefficients at the final grid.
    let temp_bytes = 2u128 * final_grid as u128 * 8;
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
    let m2 = (m as u128) * (m as u128);
    let ni = (n - 1) as u128;
    let ni2 = ni * ni;
    // d + d2 (m² each), a + shifted + LU copy (ni² each), the fixed
    // 64² FD surrogate, and two ni-length iteration vectors.
    let temp_bytes = (2 * m2 + 3 * ni2 + 64 * 64 + 2 * ni) * 8;
    cap_check(
        "eigensolve temporary bytes",
        temp_bytes,
        u128::from(budget.max_temp_bytes),
    )?;
    // dsq O(m³) + surrogate Jacobi (fixed 64³ class) + per shift:
    // LU ~ni³/3 + 100 solve/normalize sweeps ~ni² each + Rayleigh ~ni².
    let ops = m2 * (m as u128) + 64u128 * 64 * 64 * 8 + (k as u128) * (ni2 * ni / 3 + 101 * ni2);
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
/// samples (each an O(len) Clenshaw evaluation) plus up to 44
/// refinement evaluations per detected cell.
///
/// # Errors
/// [`ChebError`] — overflow or cap refusals.
pub fn admit_root_scan(coeff_len: usize, budget: &ChebBudget) -> Result<ChebAdmission, ChebError> {
    let samples = coeff_len
        .checked_mul(8)
        .ok_or(ChebError::Overflow {
            what: "root-scan sample count",
        })?
        .max(64);
    cap_check(
        "root-scan samples",
        samples as u128,
        budget.max_samples as u128,
    )?;
    let evals = samples as u128 * 45; // scan + worst-case refinement
    let ops = evals * (coeff_len as u128);
    cap_check("root-scan work", ops, u128::from(budget.max_work_ops))?;
    let temp_bytes = 2u128 * (coeff_len as u128) * 8; // reference + derivative copies
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
) -> Result<Vec<f64>, ChebError> {
    let mut vals = Vec::with_capacity(n);
    for k in 0..n {
        let theta = std::f64::consts::PI * (k as f64 + 0.5) / (n as f64);
        let t = fs_math::det::cos(theta);
        let x = affine_from_reference(t, a, b);
        let y = f(x);
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
) -> Result<Vec<f64>, ChebError> {
    let vals = sample_first_kind_checked(f, a, b, n)?;
    let mut c = fs_fft::dct2(&vals);
    let scale = 2.0 / n as f64;
    for v in &mut c {
        *v *= scale;
    }
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
    loop {
        // Bounded tile boundary: one adaptive round.
        if cx.checkpoint().is_err() {
            return Ok(BuildRun::Cancelled {
                resume_from: n,
                receipt: WorkReceipt {
                    schema_version: CHEB_BUDGET_SCHEMA_VERSION,
                    rounds_completed,
                    samples_spent,
                    ops_admitted: admission.ops_admitted(),
                },
            });
        }
        let coeffs = coeffs_at_checked(f, a, b, n)?;
        samples_spent += n;
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
            return Ok(BuildRun::Complete {
                function: Cheb1 {
                    a,
                    b,
                    coeffs: coeffs[..keep].to_vec(),
                },
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
/// shift and 10-sweep boundaries, retaining the independently
/// converged eigenvalue prefix.
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
    let m = n + 1;
    let d = diff_matrix(n);
    let mut d2 = vec![0.0f64; m * m];
    fma::dsq_into_dispatch(&d, m, &mut d2);
    let ni = n - 1;
    let mut a = vec![0.0f64; ni * ni];
    for i in 0..ni {
        for j in 0..ni {
            a[i * ni + j] = -d2[(i + 1) * m + (j + 1)];
        }
    }
    let nf = 64usize;
    let h = 2.0 / (nf as f64 + 1.0);
    let mut fd = vec![0.0f64; nf * nf];
    for i in 0..nf {
        fd[i * nf + i] = 2.0 / (h * h);
        if i + 1 < nf {
            fd[i * nf + i + 1] = -1.0 / (h * h);
            fd[(i + 1) * nf + i] = -1.0 / (h * h);
        }
    }
    let (fd_eigs, _) = fs_la::eigen::jacobi_eigh(&fd, nf);
    let mut eigs = Vec::with_capacity(k);
    let mut shifted = vec![0.0f64; a.len()];
    let mut sweeps_total = 0usize;
    let receipt = |rounds: u32, sweeps: usize| WorkReceipt {
        schema_version: CHEB_BUDGET_SCHEMA_VERSION,
        rounds_completed: rounds,
        samples_spent: sweeps,
        ops_admitted: admission.ops_admitted(),
    };
    for (shift_index, &fd_est) in fd_eigs.iter().take(k).enumerate() {
        // Shift boundary poll: the converged prefix stays valid.
        if cx.checkpoint().is_err() {
            return Ok(EigsRun::Cancelled {
                partial_eigs: eigs,
                receipt: receipt(shift_index as u32, sweeps_total),
            });
        }
        let mu = fd_est * 0.95;
        shifted.copy_from_slice(&a);
        for i in 0..ni {
            shifted[i * ni + i] -= mu;
        }
        let lu = fs_la::factor::lu(&shifted, ni).map_err(|_| ChebError::Numerical {
            what: "shifted collocation operator is singular",
        })?;
        let mut v: Vec<f64> = (0..ni)
            .map(|i| 1.0 + 0.25 * (((i * 7 + 3) % 11) as f64))
            .collect();
        for sweep in 0..100 {
            if sweep % 10 == 0 && cx.checkpoint().is_err() {
                return Ok(EigsRun::Cancelled {
                    partial_eigs: eigs,
                    receipt: receipt(shift_index as u32, sweeps_total),
                });
            }
            let nrm = v.iter().map(|x| x * x).sum::<f64>().sqrt();
            for x in &mut v {
                *x /= nrm;
            }
            lu.solve(&mut v);
            sweeps_total += 1;
        }
        let nrm2: f64 = v.iter().map(|x| x * x).sum();
        let mut av = vec![0.0f64; ni];
        fma::matvec_into_dispatch(&a, &v, ni, &mut av);
        eigs.push(v.iter().zip(&av).map(|(x, y)| x * y).sum::<f64>() / nrm2);
    }
    Ok(EigsRun::Complete {
        eigs,
        receipt: receipt(k as u32, sweeps_total),
    })
}

impl Cheb1 {
    /// Budgeted, cancellable root scan. The scan sequence, refinement,
    /// and conditioning rule match [`Cheb1::roots`] bitwise on the
    /// happy path; admission preflights the scan/refinement work;
    /// cancellation polls every 64 scan cells and refuses WITHOUT a
    /// partial result — an incomplete scan is not a root-set claim, so
    /// there is deliberately nothing acceptance-capable to return.
    ///
    /// # Errors
    /// [`ChebError`] — cap/overflow refusals before evaluation;
    /// [`ChebError::Numerical`] for the identically-zero polynomial,
    /// non-finite evaluations, or an unresolvable multiple root;
    /// [`ChebError::Cancelled`] on a drained cancellation.
    pub fn roots_budgeted(&self, budget: &ChebBudget, cx: &Cx<'_>) -> Result<Vec<f64>, ChebError> {
        let _admission = admit_root_scan(self.coeffs.len(), budget)?;
        if !self.coeffs.iter().any(|&coefficient| coefficient != 0.0) {
            return Err(ChebError::Numerical {
                what: "the identically zero polynomial has a continuum of roots",
            });
        }
        let mut reference_coeffs = self.coeffs.clone();
        crate::normalize_coefficients_exact(&mut reference_coeffs, "Chebyshev root scan");
        let reference = Cheb1 {
            a: -1.0,
            b: 1.0,
            coeffs: reference_coeffs,
        };
        let derivative = reference.differentiate();
        let resolvable = |t: f64| -> Result<(), ChebError> {
            let slope = derivative.eval_reference(t);
            let degree_scale = reference.degree().max(1) as f64;
            let slope_floor = 64.0 * 1.490_116_119_384_765_6e-8 * degree_scale;
            if slope.is_finite() && slope.abs() > slope_floor {
                Ok(())
            } else {
                Err(ChebError::Numerical {
                    what: "fixed-grid root scan cannot resolve a multiple or \
                           ill-conditioned root; use colleague/certified root evidence",
                })
            }
        };
        let finite = |v: f64| -> Result<f64, ChebError> {
            if v.is_finite() {
                Ok(v)
            } else {
                Err(ChebError::NonFinite {
                    what: "root scan evaluation",
                })
            }
        };

        let mut roots_t = Vec::new();
        let samples = self.coeffs.len().saturating_mul(8).max(64);
        let mut prev_t = -1.0;
        let mut prev_v = finite(reference.eval_reference(prev_t))?;
        for k in 1..=samples {
            // Bounded tile boundary: one poll per 64 scan cells.
            if k % 64 == 1 && cx.checkpoint().is_err() {
                return Err(ChebError::Cancelled);
            }
            let t = 2.0 * (k as f64) / (samples as f64) - 1.0;
            let v = finite(reference.eval_reference(t))?;
            if prev_v == 0.0 {
                resolvable(prev_t)?;
                roots_t.push(prev_t);
            } else if v != 0.0 && prev_v.is_sign_negative() != v.is_sign_negative() {
                let root = reference.bisect_newton_reference(&derivative, prev_t, t);
                roots_t.push(root);
            }
            prev_t = t;
            prev_v = v;
        }
        if prev_v == 0.0 {
            resolvable(prev_t)?;
            roots_t.push(prev_t);
        }
        Ok(roots_t
            .into_iter()
            .map(|t| affine_from_reference(t, self.a, self.b))
            .collect())
    }
}
