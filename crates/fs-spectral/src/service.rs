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
//! Honesty boundaries: for an actually symmetric operator, a converged Ritz
//! pair's residual `r` gives the usual numerical Weyl enclosure that SOME
//! eigenvalue lies within `r` of the Ritz value. [`SymmetricOp`] is a caller
//! precondition, not an admitted symmetry witness, and ordinary `f64`
//! residual arithmetic is not theorem-grade interval evidence. These values
//! are therefore DRAFT numerical estimates until a higher truth/admission
//! layer binds the operator proposition and validates the evidence. Per-pair
//! intervals establish neither distinctness nor exact multiplicity.

use fs_exec::Cx;
use fs_la::eigen::{EigenPair, LanczosState, LobpcgState, lanczos_run, lobpcg_run};
use std::cell::Cell;

/// Admission cap on backend steps per tick. Backend constructors separately
/// cap aggregate memory and in-house scalar work between polls; a custom
/// [`SymmetricOp`] remains responsible for bounding one `apply` call.
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
    /// An operator method panicked while panic unwinding was enabled. The
    /// tick's state mutation was rolled back; `panic = "abort"` builds cannot
    /// contain an operator panic.
    OperatorFault,
    /// A warm-start seed was rejected (wrong size, non-finite, zero, or a
    /// rank-deficient selected LOBPCG block).
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
    /// The tick budget ended before the backend had produced the requested
    /// number of pairs. This is distinct from convergence failure because no
    /// meaningful worst-of-k residual exists yet.
    Incomplete {
        /// Cumulative service ticks consumed.
        ticks: usize,
        /// Pairs currently available.
        available: usize,
        /// Pairs requested by the query.
        requested: usize,
    },
    /// The next bounded tick would exceed the admitted aggregate-memory or
    /// unpolled scalar-work cap. State is unchanged; the caller must choose a
    /// smaller problem or an explicit restart policy.
    WorkLimit {
        /// Rejected backend lane.
        backend: &'static str,
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
            ServiceError::OperatorFault => {
                write!(f, "operator method panicked; state rolled back")
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
            ServiceError::Incomplete {
                ticks,
                available,
                requested,
            } => write!(
                f,
                "eigen service incomplete after {ticks} ticks ({available}/{requested} pairs); \
                 state remains resumable"
            ),
            ServiceError::WorkLimit { backend } => write!(
                f,
                "{backend} tick exceeds aggregate memory or unpolled-work cap; state unchanged"
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
    /// Convergence tolerance on numerically recomputed residual norms.
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

/// One converged (or in-progress) Ritz estimate with a derived numerical
/// residual interval. Fields are sealed so callers cannot inject a reversed
/// or unrelated interval. Despite the historical type name, this is not an
/// admitted scientific-authority token; see the module honesty boundary.
#[derive(Debug, Clone)]
pub struct CertifiedEigenvalue {
    /// The Ritz value.
    value: f64,
    /// Reported operator residual of the Ritz pair. Service-produced values
    /// are recomputed by the backend; the public draft constructor is not an
    /// authority verifier.
    residual: f64,
    /// Derived numerical interval, never caller-supplied.
    interval: (f64, f64),
    /// Candidate Ritz vector — warm-start fodder. Service-produced vectors are
    /// normalized by the backend; draft construction only checks finiteness.
    vector: Vec<f64>,
}

impl CertifiedEigenvalue {
    /// Build a draft Ritz estimate from a finite value, nonnegative finite
    /// residual, and finite nonempty vector. Interval endpoints are derived
    /// here and widened by one representable value when the residual is
    /// nonzero. This prevents malformed interval injection; it does not turn
    /// caller-supplied residuals into admitted scientific evidence.
    pub fn from_residual(
        value: f64,
        residual: f64,
        vector: Vec<f64>,
    ) -> Result<Self, ServiceError> {
        if !value.is_finite()
            || !residual.is_finite()
            || residual < 0.0
            || vector.is_empty()
            || vector.iter().any(|x| !x.is_finite())
        {
            return Err(ServiceError::InvalidQuery {
                what: "Ritz value/residual/vector must be finite and residual nonnegative",
            });
        }
        let raw_lo = value - residual;
        let raw_hi = value + residual;
        if !(raw_lo.is_finite() && raw_hi.is_finite() && raw_lo <= raw_hi) {
            return Err(ServiceError::InvalidQuery {
                what: "derived Ritz interval must have finite ordered endpoints",
            });
        }
        let interval = if residual > 0.0 {
            (raw_lo.next_down(), raw_hi.next_up())
        } else {
            (value, value)
        };
        if !(interval.0.is_finite() && interval.1.is_finite() && interval.0 <= interval.1) {
            return Err(ServiceError::InvalidQuery {
                what: "outward-widened Ritz interval must remain finite and ordered",
            });
        }
        Ok(Self {
            value,
            residual,
            interval,
            vector,
        })
    }

    /// Ritz value.
    #[must_use]
    pub fn value(&self) -> f64 {
        self.value
    }

    /// Reported residual norm.
    #[must_use]
    pub fn residual(&self) -> f64 {
        self.residual
    }

    /// Derived finite, ordered numerical interval.
    #[must_use]
    pub fn interval(&self) -> (f64, f64) {
        self.interval
    }

    /// Ritz vector for warm starts.
    #[must_use]
    pub fn vector(&self) -> &[f64] {
        &self.vector
    }
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
    /// Numerical lower bound between the first two cluster hulls (0 when they
    /// touch or only one cluster exists). This is not an authority token.
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
        if backend == EigenBackend::Lanczos
            && !LanczosState::initial_work_is_admitted(n, query.steps_per_tick, query.k)
        {
            return Err(ServiceError::InvalidQuery {
                what: "Lanczos first tick exceeds aggregate memory or unpolled-work cap",
            });
        }
        let state = match backend {
            EigenBackend::Lanczos => BackendState::Lanczos(LanczosState::try_new(n).ok_or(
                ServiceError::InvalidQuery {
                    what: "Lanczos dimension exceeds the practical work cap",
                },
            )?),
            EigenBackend::Lobpcg => BackendState::Lobpcg(LobpcgState::try_new(n, query.k).ok_or(
                ServiceError::InvalidQuery {
                    what: "LOBPCG aggregate workspace exceeds the practical work cap",
                },
            )?),
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
        if seed_vectors.is_empty()
            || seed_vectors.iter().any(|v| {
                v.len() != n
                    || v.iter().any(|value| !value.is_finite())
                    || v.iter().all(|value| *value == 0.0)
            })
        {
            return Err(ServiceError::InvalidSeed);
        }
        if backend == EigenBackend::Lanczos
            && !LanczosState::initial_work_is_admitted(n, query.steps_per_tick, query.k)
        {
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
                if !LobpcgState::shape_is_admitted(n, query.k) {
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
    /// (dimension mismatch, non-finite operator output, operator unwind,
    /// cancellation)
    /// ROLL BACK the in-flight tick, so accepted state is never
    /// corrupted; `clone()` before or after any tick is a valid
    /// checkpoint and split runs replay bitwise-identically.
    pub fn tick(
        &mut self,
        op: &dyn SymmetricOp,
        cx: &Cx<'_>,
    ) -> Result<EigenProgress, ServiceError> {
        let op_dim = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| op.dim()))
            .map_err(|_| ServiceError::OperatorFault)?;
        if op_dim != self.n {
            return Err(ServiceError::DimensionMismatch {
                expected: self.n,
                got: op_dim,
            });
        }
        if cx.checkpoint().is_err() {
            return Err(ServiceError::Cancelled);
        }
        if matches!(
            &self.state,
            BackendState::Lanczos(state)
                if !state.work_is_admitted(self.query.steps_per_tick, self.query.k)
        ) {
            return Err(ServiceError::WorkLimit { backend: "Lanczos" });
        }
        let snapshot = self.state.clone();
        let nonfinite = Cell::new(false);
        let operator_fault = Cell::new(false);
        let cancelled = Cell::new(false);
        let apply = |x: &[f64], y: &mut [f64]| {
            if cx.checkpoint().is_err() {
                cancelled.set(true);
                y.fill(0.0);
                return;
            }
            if let Err(payload) =
                std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| op.apply(x, y)))
            {
                operator_fault.set(true);
                std::panic::resume_unwind(payload);
            }
            if y.iter().any(|v| !v.is_finite()) {
                nonfinite.set(true);
            }
            if cx.checkpoint().is_err() {
                cancelled.set(true);
            }
        };
        let mut pairs: Vec<EigenPair> = Vec::new();
        for _ in 0..self.query.steps_per_tick {
            if cx.checkpoint().is_err() {
                self.state = snapshot;
                return Err(ServiceError::Cancelled);
            }
            let step =
                std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| match &mut self.state {
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
                }));
            pairs = match step {
                Ok(pairs) => pairs,
                Err(_) if operator_fault.get() => {
                    self.state = snapshot;
                    return Err(ServiceError::OperatorFault);
                }
                Err(payload) => {
                    self.state = snapshot;
                    std::panic::resume_unwind(payload)
                }
            };
            if cancelled.get() {
                self.state = snapshot;
                return Err(ServiceError::Cancelled);
            }
            if nonfinite.get() {
                self.state = snapshot;
                return Err(ServiceError::NonFiniteOperator);
            }
        }
        if pairs.iter().any(|p| {
            !p.value.is_finite()
                || !p.residual.is_finite()
                || p.residual < 0.0
                || p.vector.is_empty()
                || p.vector.iter().any(|x| !x.is_finite())
        }) {
            self.state = snapshot;
            return Err(ServiceError::NonFiniteOperator);
        }
        let mut estimates = Vec::with_capacity(pairs.len());
        for p in pairs {
            let estimate = match CertifiedEigenvalue::from_residual(p.value, p.residual, p.vector) {
                Ok(estimate) => estimate,
                Err(_) => {
                    self.state = snapshot;
                    return Err(ServiceError::NonFiniteOperator);
                }
            };
            estimates.push(estimate);
        }
        estimates.sort_by(|a, b| a.value.total_cmp(&b.value));
        // A request can arrive during the final backend application or Ritz
        // extraction. Poll once more immediately before committing the state;
        // otherwise a one-step tick could acknowledge success after a request.
        if cx.checkpoint().is_err() {
            self.state = snapshot;
            return Err(ServiceError::Cancelled);
        }
        self.ticks = match self.ticks.checked_add(1) {
            Some(ticks) => ticks,
            None => {
                self.state = snapshot;
                return Err(ServiceError::InvalidQuery {
                    what: "service tick counter overflow",
                });
            }
        };
        // The backends return exactly the wanted extremal pairs, so
        // convergence is: enough of them, and every one within
        // tolerance.
        let converged = estimates.len() >= self.query.k
            && estimates.iter().all(|p| p.residual <= self.query.tol);
        let subspace_exhausted = match &self.state {
            BackendState::Lanczos(state) => state.exhausted(),
            BackendState::Lobpcg(_) => false,
        };
        Ok(EigenProgress {
            pairs: estimates,
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
        let last = last.expect("positive max_ticks always produces progress or returns early");
        if last.pairs.len() < self.query.k {
            return Err(ServiceError::Incomplete {
                ticks: self.ticks,
                available: last.pairs.len(),
                requested: self.query.k,
            });
        }
        let worst = last
            .pairs
            .iter()
            .take(self.query.k)
            .map(|c| c.residual)
            .reduce(f64::max)
            .expect("validated k and pair-count gate guarantee a residual");
        debug_assert!(worst.is_finite() && worst >= 0.0);
        Err(ServiceError::Unconverged {
            ticks: self.ticks,
            worst_residual: worst,
        })
    }
}

/// Cluster pairs by interval overlap and report the leading numerical gap
/// lower bound. Sorting is by interval LOWER bound (then
/// upper, then value — deterministic total order), and merging
/// extends BOTH hull endpoints, so transitively bridging intervals
/// (e.g. `[0,0]`, `[1,1]`, `[-1,5]`) collapse into one cluster
/// instead of reporting a false gap. Non-finite intervals are
/// refused by `tick` before they can reach a report; any that arrive
/// here through malformed internal data are excluded from clustering.
#[must_use]
pub fn gap_report(pairs: &[CertifiedEigenvalue]) -> GapReport {
    let mut sorted: Vec<&CertifiedEigenvalue> = pairs
        .iter()
        .filter(|p| {
            p.interval.0.is_finite() && p.interval.1.is_finite() && p.interval.0 <= p.interval.1
        })
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
        let raw_gap = clusters[1].hull.0 - clusters[0].hull.1;
        if raw_gap.is_finite() {
            raw_gap.max(0.0)
        } else {
            // Both hull endpoints are finite and ordered, so positive
            // overflow means the real gap exceeds f64::MAX. Clamp to the
            // largest finite conservative lower bound rather than emitting an
            // infinite numerical artifact.
            f64::MAX
        }
    } else {
        0.0
    };
    GapReport {
        clusters,
        leading_gap_lower_bound,
    }
}
