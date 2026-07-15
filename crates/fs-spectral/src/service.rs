//! Generic resumable symmetric eigensolver service (bead
//! `frankensim-ext-spectral-eigensolver-service-bfid`).
//!
//! Ownership rule (normative): generic operator spectra live HERE, at
//! L1; domain crates assemble operators, adapt them DOWNWARD to
//! [`SymmetricOp`], and interpret results (an fs-solver `LinearOp`
//! adapter is an L3 shim, deliberately not defined here). This module
//! wraps fs-la's deterministic Lanczos / LOBPCG backends behind one
//! service with the house contracts: resumable plain-data state
//! (clone = checkpoint, split runs bitwise-equal to straight runs),
//! bounded cancellable ticks under an explicit `Cx`, typed refusals
//! that never corrupt accepted state, and warm-start hooks for
//! parameter continuation. The existing dense monitor path
//! (`spectral_gap`, `GapHealthMonitor`, `propagate`) remains the
//! interpretation layer above this backend.
//!
//! Honesty boundaries: a converged Ritz pair's residual `r` certifies
//! (Weyl) that SOME eigenvalue lies within `r` of the Ritz value —
//! per-pair intervals certify existence, not distinctness. Cluster
//! reports are therefore "multiplicity at the achieved resolution",
//! never exact multiplicity claims; the certified gap lower bound is
//! between cluster interval hulls and can be zero when clusters touch.

use fs_exec::Cx;
use fs_la::eigen::{EigenPair, LanczosState, LobpcgState, lanczos_run, lobpcg_run};
use std::cell::Cell;

/// Admission cap on backend steps per tick: one tick is a bounded,
/// cancellation-polled unit of work, never an unbounded burst.
pub const MAX_STEPS_PER_TICK: usize = 1024;

/// Typed refusals for the eigensolver service. (No `Eq`: the
/// unconverged variant carries an f64 residual.)
#[derive(Debug, Clone, PartialEq)]
pub enum ServiceError {
    /// A dimension, count, or block size is unusable.
    InvalidQuery {
        /// The rejected condition.
        what: &'static str,
    },
    /// The operator's dimension does not match the service state.
    DimensionMismatch {
        /// Service dimension.
        expected: usize,
        /// Operator dimension.
        got: usize,
    },
    /// The operator produced a non-finite value (the tick's state
    /// mutation was rolled back; the service remains usable).
    NonFiniteOperator,
    /// A warm-start seed was rejected (wrong size, non-finite, zero,
    /// or rank-deficient).
    InvalidSeed,
    /// The Krylov space is exhausted (invariant subspace) with fewer
    /// pairs than requested; more ticks cannot help.
    SubspaceExhausted {
        /// Pairs available from the exhausted subspace.
        available: usize,
    },
    /// The tick budget ran out before the tolerance was met. The state
    /// inside is still valid and resumable.
    Unconverged {
        /// Ticks consumed.
        ticks: usize,
        /// Worst residual among the wanted pairs at the last tick.
        worst_residual: f64,
    },
    /// Cooperative cancellation observed; the in-flight tick was
    /// rolled back and the service is resumable.
    Cancelled,
}

impl std::fmt::Display for ServiceError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ServiceError::InvalidQuery { what } => write!(f, "invalid eigen query: {what}"),
            ServiceError::DimensionMismatch { expected, got } => write!(
                f,
                "operator dimension {got} does not match service dimension {expected}"
            ),
            ServiceError::NonFiniteOperator => {
                write!(f, "operator application produced a non-finite value")
            }
            ServiceError::InvalidSeed => write!(f, "warm-start seed rejected"),
            ServiceError::SubspaceExhausted { available } => write!(
                f,
                "Krylov space exhausted with only {available} pair(s) available"
            ),
            ServiceError::Unconverged {
                ticks,
                worst_residual,
            } => write!(
                f,
                "eigen service unconverged after {ticks} ticks (worst residual \
                 {worst_residual:e}); state remains resumable"
            ),
            ServiceError::Cancelled => {
                write!(
                    f,
                    "cancelled at a tick boundary; state rolled back and resumable"
                )
            }
        }
    }
}

impl std::error::Error for ServiceError {}

/// The L1-safe symmetric operator abstraction. Domain crates adapt
/// their operator types (e.g. fs-solver `LinearOp`s) downward to this
/// trait; fs-spectral never imports upward.
pub trait SymmetricOp {
    /// Operator dimension (square).
    fn dim(&self) -> usize;
    /// `y ← A·x` (symmetric A; the service does not verify symmetry —
    /// a nonsymmetric operator voids every claim here).
    fn apply(&self, x: &[f64], y: &mut [f64]);
}

/// Row-major dense symmetric operator (small systems and tests).
#[derive(Debug, Clone)]
pub struct DenseSymOp {
    n: usize,
    a: Vec<f64>,
}

impl DenseSymOp {
    /// Wrap a row-major n×n matrix. Refuses wrong sizes (checked
    /// arithmetic — hostile n cannot overflow) or non-finite entries;
    /// exact symmetry is checked because it is cheap and load-bearing
    /// for every downstream claim.
    pub fn new(n: usize, a: Vec<f64>) -> Result<DenseSymOp, ServiceError> {
        let len = n.checked_mul(n).ok_or(ServiceError::InvalidQuery {
            what: "dense operator dimension overflows",
        })?;
        if n == 0 || a.len() != len {
            return Err(ServiceError::InvalidQuery {
                what: "dense operator must be square and non-empty",
            });
        }
        if a.iter().any(|x| !x.is_finite()) {
            return Err(ServiceError::NonFiniteOperator);
        }
        for i in 0..n {
            for j in (i + 1)..n {
                if a[i * n + j] != a[j * n + i] {
                    return Err(ServiceError::InvalidQuery {
                        what: "dense operator is not exactly symmetric",
                    });
                }
            }
        }
        Ok(DenseSymOp { n, a })
    }
}

impl SymmetricOp for DenseSymOp {
    fn dim(&self) -> usize {
        self.n
    }

    fn apply(&self, x: &[f64], y: &mut [f64]) {
        for i in 0..self.n {
            let row = &self.a[i * self.n..(i + 1) * self.n];
            let mut acc = 0.0f64;
            for (aij, xj) in row.iter().zip(x) {
                acc = aij.mul_add(*xj, acc);
            }
            y[i] = acc;
        }
    }
}

/// Which fs-la backend drives the service.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EigenBackend {
    /// Krylov tridiagonalization with full reorthogonalization —
    /// extremal pairs of large sparse operators.
    Lanczos,
    /// Blocked preconditioned iteration — clustered/multiple extremal
    /// pairs (block size = the query's `k`).
    Lobpcg,
}

/// What the caller wants.
#[derive(Debug, Clone, Copy)]
pub struct EigenQuery {
    /// Number of extremal pairs.
    pub k: usize,
    /// Which end of the spectrum.
    pub largest: bool,
    /// Convergence tolerance on TRUE residual norms.
    pub tol: f64,
    /// Backend steps per `tick` call (resume + cancellation
    /// granularity; admission-capped by [`MAX_STEPS_PER_TICK`]).
    pub steps_per_tick: usize,
}

/// Plain-data resumable state: `clone()` IS the checkpoint.
#[derive(Debug, Clone)]
enum BackendState {
    Lanczos(LanczosState),
    Lobpcg(LobpcgState),
}

/// One converged (or in-progress) eigenvalue with its certificate.
#[derive(Debug, Clone)]
pub struct CertifiedEigenvalue {
    /// The Ritz value.
    pub value: f64,
    /// TRUE operator residual of the Ritz pair.
    pub residual: f64,
    /// Weyl containment: some eigenvalue of the operator lies in
    /// `[value − residual, value + residual]`. Existence, not
    /// distinctness.
    pub interval: (f64, f64),
    /// The Ritz vector (unit norm) — warm-start fodder.
    pub vector: Vec<f64>,
}

/// A cluster of eigenvalue intervals that overlap at the achieved
/// resolution.
#[derive(Debug, Clone)]
pub struct EigenCluster {
    /// Hull of the member intervals.
    pub hull: (f64, f64),
    /// Members at this resolution — NOT an exact multiplicity claim.
    pub count: usize,
}

/// Multiplicity-aware gap report over the converged pairs.
#[derive(Debug, Clone)]
pub struct GapReport {
    /// Clusters in ascending hull order.
    pub clusters: Vec<EigenCluster>,
    /// Certified lower bound of the gap between the first two cluster
    /// hulls (0 when they touch or only one cluster exists).
    pub leading_gap_lower_bound: f64,
}

/// Progress after one tick.
#[derive(Debug, Clone)]
pub struct EigenProgress {
    /// Current pairs (ascending by value).
    pub pairs: Vec<CertifiedEigenvalue>,
    /// Whether every wanted pair meets the tolerance.
    pub converged: bool,
    /// The Krylov space is exhausted: no further tick can enlarge it
    /// (Lanczos invariant-subspace breakdown; sticky).
    pub subspace_exhausted: bool,
    /// Ticks consumed so far.
    pub ticks: usize,
}

/// The resumable eigensolver service.
#[derive(Debug, Clone)]
pub struct EigenService {
    n: usize,
    query: EigenQuery,
    state: BackendState,
    ticks: usize,
}

fn validate_query(n: usize, query: &EigenQuery, backend: EigenBackend) -> Result<(), ServiceError> {
    if query.k == 0 || query.k > n {
        return Err(ServiceError::InvalidQuery {
            what: "k must satisfy 1 <= k <= n",
        });
    }
    if !(query.tol > 0.0 && query.tol.is_finite()) {
        return Err(ServiceError::InvalidQuery {
            what: "tolerance must be positive and finite",
        });
    }
    if query.steps_per_tick == 0 || query.steps_per_tick > MAX_STEPS_PER_TICK {
        return Err(ServiceError::InvalidQuery {
            what: "steps_per_tick must be in 1..=MAX_STEPS_PER_TICK",
        });
    }
    if backend == EigenBackend::Lobpcg {
        let three_k = query.k.checked_mul(3).ok_or(ServiceError::InvalidQuery {
            what: "LOBPCG block size overflows",
        })?;
        if three_k > n {
            return Err(ServiceError::InvalidQuery {
                what: "LOBPCG block needs 3k <= n",
            });
        }
    }
    Ok(())
}

impl EigenService {
    /// Cold start with the backend's deterministic seed.
    pub fn new(
        backend: EigenBackend,
        n: usize,
        query: EigenQuery,
    ) -> Result<EigenService, ServiceError> {
        validate_query(n, &query, backend)?;
        let state = match backend {
            EigenBackend::Lanczos => BackendState::Lanczos(LanczosState::new(n)),
            EigenBackend::Lobpcg => BackendState::Lobpcg(LobpcgState::new(n, query.k)),
        };
        Ok(EigenService {
            n,
            query,
            state,
            ticks: 0,
        })
    }

    /// Warm start from previous Ritz vectors (parameter continuation).
    /// Lanczos seeds with the FIRST vector; LOBPCG seeds the whole
    /// block. Rejected seeds (wrong size, non-finite, zero,
    /// rank-deficient, overflowing dimensions) refuse rather than
    /// silently cold-start.
    pub fn warm(
        backend: EigenBackend,
        n: usize,
        query: EigenQuery,
        seed_vectors: &[Vec<f64>],
    ) -> Result<EigenService, ServiceError> {
        validate_query(n, &query, backend)?;
        if seed_vectors.is_empty() || seed_vectors.iter().any(|v| v.len() != n) {
            return Err(ServiceError::InvalidSeed);
        }
        let state = match backend {
            EigenBackend::Lanczos => BackendState::Lanczos(
                LanczosState::with_start(&seed_vectors[0]).ok_or(ServiceError::InvalidSeed)?,
            ),
            EigenBackend::Lobpcg => {
                if seed_vectors.len() < query.k {
                    return Err(ServiceError::InvalidSeed);
                }
                // Checked size arithmetic: hostile n·k must refuse,
                // never wrap or abort in the allocator.
                let len = n.checked_mul(query.k).ok_or(ServiceError::InvalidSeed)?;
                let mut block = vec![0.0f64; len];
                for (j, v) in seed_vectors.iter().take(query.k).enumerate() {
                    for i in 0..n {
                        block[i * query.k + j] = v[i];
                    }
                }
                BackendState::Lobpcg(
                    LobpcgState::with_block(n, query.k, &block).ok_or(ServiceError::InvalidSeed)?,
                )
            }
        };
        Ok(EigenService {
            n,
            query,
            state,
            ticks: 0,
        })
    }

    /// Ticks consumed.
    #[must_use]
    pub fn ticks(&self) -> usize {
        self.ticks
    }

    /// The service dimension.
    #[must_use]
    pub fn dim(&self) -> usize {
        self.n
    }

    /// Advance one bounded tick (`steps_per_tick` backend steps with a
    /// cancellation poll between steps) and report. Error paths
    /// (dimension mismatch, non-finite operator output, cancellation)
    /// ROLL BACK the in-flight tick, so accepted state is never
    /// corrupted; `clone()` before or after any tick is a valid
    /// checkpoint and split runs replay bitwise-identically.
    pub fn tick(
        &mut self,
        op: &dyn SymmetricOp,
        cx: &Cx<'_>,
    ) -> Result<EigenProgress, ServiceError> {
        if op.dim() != self.n {
            return Err(ServiceError::DimensionMismatch {
                expected: self.n,
                got: op.dim(),
            });
        }
        let snapshot = self.state.clone();
        let nonfinite = Cell::new(false);
        let apply = |x: &[f64], y: &mut [f64]| {
            op.apply(x, y);
            if y.iter().any(|v| !v.is_finite()) {
                nonfinite.set(true);
            }
        };
        let mut pairs: Vec<EigenPair> = Vec::new();
        for _ in 0..self.query.steps_per_tick {
            if cx.checkpoint().is_err() {
                self.state = snapshot;
                return Err(ServiceError::Cancelled);
            }
            pairs = match &mut self.state {
                BackendState::Lanczos(state) => {
                    lanczos_run(&apply, state, 1, self.query.k, self.query.largest)
                }
                BackendState::Lobpcg(state) => lobpcg_run(
                    &apply,
                    state,
                    1,
                    self.query.largest,
                    &|r: &[f64], out: &mut [f64]| out.copy_from_slice(r),
                ),
            };
            if nonfinite.get() {
                self.state = snapshot;
                return Err(ServiceError::NonFiniteOperator);
            }
        }
        self.ticks += 1;
        if pairs
            .iter()
            .any(|p| !p.value.is_finite() || !p.residual.is_finite())
        {
            self.state = snapshot;
            self.ticks -= 1;
            return Err(ServiceError::NonFiniteOperator);
        }
        let mut certified = Vec::with_capacity(pairs.len());
        for p in pairs {
            let lo = p.value - p.residual;
            let hi = p.value + p.residual;
            if !(lo.is_finite() && hi.is_finite()) {
                self.state = snapshot;
                self.ticks -= 1;
                return Err(ServiceError::NonFiniteOperator);
            }
            certified.push(CertifiedEigenvalue {
                value: p.value,
                residual: p.residual,
                interval: (lo, hi),
                vector: p.vector,
            });
        }
        certified.sort_by(|a, b| a.value.total_cmp(&b.value));
        // The backends return exactly the wanted extremal pairs, so
        // convergence is: enough of them, and every one within
        // tolerance.
        let converged = certified.len() >= self.query.k
            && certified.iter().all(|p| p.residual <= self.query.tol);
        let subspace_exhausted = match &self.state {
            BackendState::Lanczos(state) => state.exhausted(),
            BackendState::Lobpcg(_) => false,
        };
        Ok(EigenProgress {
            pairs: certified,
            converged,
            subspace_exhausted,
            ticks: self.ticks,
        })
    }

    /// Drive `tick` until convergence, subspace exhaustion, or the
    /// tick budget. Budget exhaustion and subspace exhaustion are
    /// typed errors; the service remains valid and resumable.
    pub fn run_to_tolerance(
        &mut self,
        op: &dyn SymmetricOp,
        cx: &Cx<'_>,
        max_ticks: usize,
    ) -> Result<EigenProgress, ServiceError> {
        if max_ticks == 0 {
            return Err(ServiceError::InvalidQuery {
                what: "max_ticks must be positive",
            });
        }
        let mut last: Option<EigenProgress> = None;
        for _ in 0..max_ticks {
            let progress = self.tick(op, cx)?;
            if progress.converged {
                return Ok(progress);
            }
            if progress.subspace_exhausted {
                if progress.pairs.len() >= self.query.k
                    && progress.pairs.iter().all(|p| p.residual <= self.query.tol)
                {
                    return Ok(progress);
                }
                return Err(ServiceError::SubspaceExhausted {
                    available: progress.pairs.len(),
                });
            }
            last = Some(progress);
        }
        let worst = last
            .as_ref()
            .map(|p| p.pairs.iter().map(|c| c.residual).fold(0.0f64, f64::max))
            .unwrap_or(f64::INFINITY);
        Err(ServiceError::Unconverged {
            ticks: self.ticks,
            worst_residual: worst,
        })
    }
}

/// Cluster pairs by interval overlap and report the leading gap's
/// certified lower bound. Sorting is by interval LOWER bound (then
/// upper, then value — deterministic total order), and merging
/// extends BOTH hull endpoints, so transitively bridging intervals
/// (e.g. `[0,0]`, `[1,1]`, `[-1,5]`) collapse into one cluster
/// instead of reporting a false gap. Non-finite intervals are
/// refused by `tick` before they can reach a report; any that arrive
/// here through hand-built pairs are excluded from clustering.
#[must_use]
pub fn gap_report(pairs: &[CertifiedEigenvalue]) -> GapReport {
    let mut sorted: Vec<&CertifiedEigenvalue> = pairs
        .iter()
        .filter(|p| p.interval.0.is_finite() && p.interval.1.is_finite())
        .collect();
    sorted.sort_by(|a, b| {
        a.interval
            .0
            .total_cmp(&b.interval.0)
            .then(a.interval.1.total_cmp(&b.interval.1))
            .then(a.value.total_cmp(&b.value))
    });
    let mut clusters: Vec<EigenCluster> = Vec::new();
    for p in sorted {
        match clusters.last_mut() {
            Some(c) if p.interval.0 <= c.hull.1 => {
                c.hull.1 = c.hull.1.max(p.interval.1);
                c.count += 1;
            }
            _ => clusters.push(EigenCluster {
                hull: p.interval,
                count: 1,
            }),
        }
    }
    let leading_gap_lower_bound = if clusters.len() >= 2 {
        (clusters[1].hull.0 - clusters[0].hull.1).max(0.0)
    } else {
        0.0
    };
    GapReport {
        clusters,
        leading_gap_lower_bound,
    }
}
