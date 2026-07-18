//! CMA-ES as natural-gradient IGO (plan §9.3, Bet 6): weighted
//! recombination with log-rank weights, rank-µ + rank-1 covariance
//! updates, cumulative step-size adaptation — the standard Hansen
//! couplings, which ARE the natural-gradient couplings on the Gaussian
//! statistical manifold. Rank-based selection makes the evolution
//! invariant to monotone transforms of the objective BY CONSTRUCTION
//! (property-tested bitwise, not cited).
//!
//! Determinism: sampling from a keyed Philox stream, `total_cmp` ranking
//! with lowest-index tie-breaks, fixed eigendecomposition cadence via
//! the landed cyclic Jacobi — the trajectory is a pure function of the
//! seed.

use fs_la::eigen::jacobi_eigh;
use fs_rand::StreamKey;

/// Kernel id for CMA sampling streams (stable registry).
const K_CMA: u32 = 0xD1F0;

/// Domain stride used by the versioned checked BIPOP seed derivation.
/// The checked entry point refuses before any callback when the complete
/// conservative restart range could wrap this coordinate.
const BIPOP_RESTART_SEED_STRIDE: u64 = 0x9E37_79B9;

/// Largest zero-based population-doubling rung launched by BIPOP.
const BIPOP_LARGE_RUN_CAP: u32 = 8;

/// Local generation envelope assigned to one restart.
const BIPOP_GENERATIONS_PER_RESTART: usize = 250;

/// Sum of the large-population scales that can finance a later small restart.
/// Production terminates immediately after rung eight, so only rungs `0..=7`
/// contribute: `2^0 + ... + 2^7 = 2^8 - 1`.
const BIPOP_PRE_TERMINAL_LARGE_RUN_SCALE_SUM: usize = (1usize << BIPOP_LARGE_RUN_CAP as usize) - 1;

/// A full `250*lambda` local cap admits the initial callback and at most 249
/// complete populations: `1 + 249*lambda` callbacks in total.
const BIPOP_FULL_ENVELOPE_GENERATIONS: usize = BIPOP_GENERATIONS_PER_RESTART - 1;

/// One `(seed, kernel, tile)` Philox stream has exactly `2^64` distinct
/// counter positions. Consuming all of them once is valid; requesting one more
/// would reuse position zero because `Stream` advances its `u64` index modulo
/// this cardinality.
const FS_RAND_STREAM_COUNTER_CARDINALITY: u128 = 1u128 << 64;

/// Schema version for [`BipopAdmission`].
pub const BIPOP_ADMISSION_SCHEMA_VERSION: u32 = 3;

/// Schema version for [`BipopRestartRecord`].
pub const BIPOP_RESTART_SCHEMA_VERSION: u32 = 1;

/// Sealed result of callback-free BIPOP input and arithmetic admission.
///
/// The restart bound combines two independent proofs. One objective evaluation
/// is the minimum spend of a launched restart, so `total_budget - 1` is a
/// budget bound. Independently, production admits nine large rungs (`0..=8`)
/// and terminates immediately after publishing rung eight. Only rungs `0..=7`
/// can therefore finance later small restarts. A full local envelope on rung
/// `r` spends at most `1 + 249*base_lambda*2^r` callbacks. Before the `j`th
/// small restart can launch, its predecessor small spend still trails the
/// cumulative large spend. While more than `base_lambda` aggregate callbacks
/// remain, every continuing small restart completes at least one population
/// and therefore spends at least `base_lambda + 1` callbacks. At most
/// `ceil((8 + 249*base_lambda*(2^8-1))/(base_lambda+1))` such restarts can
/// launch. Once no more than `base_lambda` callbacks remain, at most that many
/// further launches are possible: every continuing launch spends at least one
/// callback, while a zero-callback generated-start refusal terminates the run.
/// Adding the eight preceding large ordinals gives the scheduler bound below.
/// Admission uses the minimum of this theorem and the budget bound.
///
/// The receipt also proves that neither the shared restart-perturbation stream
/// nor any per-restart CMA stream can reuse a Philox counter coordinate. Holding
/// it proves only that these formulas were representable; it is not an
/// authenticated identity for the start, target, sigma, seed, or callback.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BipopAdmission {
    schema_version: u32,
    stream_semantics_version: u32,
    jacobi_admission: Option<fs_la::eigen::JacobiEighAdmission>,
    dimension: usize,
    total_budget: usize,
    base_lambda: usize,
    max_large_lambda: usize,
    max_local_budget: usize,
    max_restart_ordinal: u64,
    max_matrix_entries: usize,
    max_population_entries: usize,
    max_restart_stream_blocks: u128,
    max_cma_stream_blocks: u128,
}

impl BipopAdmission {
    /// Admission schema used for the receipt.
    #[must_use]
    pub fn schema_version(&self) -> u32 {
        self.schema_version
    }

    /// fs-rand stream semantics whose exact two-block normal consumption is
    /// assumed by both counter-range proofs.
    #[must_use]
    pub fn stream_semantics_version(&self) -> u32 {
        self.stream_semantics_version
    }

    /// Exact fs-la dense-Jacobi capability receipt when the aggregate budget
    /// can reach at least one complete CMA generation. Initial-evaluation-only
    /// schedules return `None` because they never invoke that dependency.
    #[must_use]
    pub fn jacobi_admission(&self) -> Option<fs_la::eigen::JacobiEighAdmission> {
        self.jacobi_admission
    }

    /// Decision-vector dimension admitted for every restart.
    #[must_use]
    pub fn dimension(&self) -> usize {
        self.dimension
    }

    /// Hard aggregate callback budget.
    #[must_use]
    pub fn total_budget(&self) -> usize {
        self.total_budget
    }

    /// Standard population at the first large and every small restart.
    #[must_use]
    pub fn base_lambda(&self) -> usize {
        self.base_lambda
    }

    /// Conservative population representability bound for the admitted ladder.
    #[must_use]
    pub fn max_large_lambda(&self) -> usize {
        self.max_large_lambda
    }

    /// Largest pre-minimum local budget formula (`lambda * 250`).
    #[must_use]
    pub fn max_local_budget(&self) -> usize {
        self.max_local_budget
    }

    /// Conservative largest restart ordinal under both budget and ladder caps.
    #[must_use]
    pub fn max_restart_ordinal(&self) -> u64 {
        self.max_restart_ordinal
    }

    /// Largest dense square-matrix element count needed by one CMA run.
    #[must_use]
    pub fn max_matrix_entries(&self) -> usize {
        self.max_matrix_entries
    }

    /// Largest population-coordinate element count needed by one CMA run.
    #[must_use]
    pub fn max_population_entries(&self) -> usize {
        self.max_population_entries
    }

    /// Most Philox counter blocks consumed by the shared restart stream.
    #[must_use]
    pub fn max_restart_stream_blocks(&self) -> u128 {
        self.max_restart_stream_blocks
    }

    /// Most Philox counter blocks consumed by any one CMA restart stream.
    #[must_use]
    pub fn max_cma_stream_blocks(&self) -> u128 {
        self.max_cma_stream_blocks
    }
}

/// Structured refusal from [`admit_bipop`] or [`try_bipop_cmaes`].
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BipopError {
    /// A zero-dimensional decision cannot enter CMA.
    EmptyStart,
    /// One supplied coordinate is NaN or infinite.
    NonFiniteStart {
        /// Coordinate index.
        component: usize,
        /// Exact IEEE-754 payload.
        bits: u64,
    },
    /// Sigma must be finite and strictly positive (including refusing ±0).
    InvalidSigma {
        /// Exact IEEE-754 payload.
        bits: u64,
    },
    /// `Some(target)` must be finite; `None` is the typed no-target spelling.
    NonFiniteTarget {
        /// Exact IEEE-754 payload.
        bits: u64,
    },
    /// No restart can be launched under a zero callback budget.
    ZeroBudget,
    /// A dimension-derived storage formula is not representable.
    DimensionOverflow {
        /// Stable formula label.
        what: &'static str,
    },
    /// The budget/ladder restart envelope cannot be represented as a u64 ordinal.
    RestartOrdinalOverflow {
        /// Requested aggregate budget.
        total_budget: usize,
    },
    /// The conservative derived-seed range would wrap u64.
    SeedRangeOverflow {
        /// Root seed.
        seed: u64,
        /// Largest conservatively reachable restart ordinal.
        max_restart_ordinal: u64,
    },
    /// The large-population ladder cannot be represented.
    PopulationOverflow {
        /// Zero-based large-run rung.
        large_run: u32,
    },
    /// `lambda * 250` cannot be represented.
    LocalBudgetOverflow {
        /// Population whose local budget overflowed.
        lambda: usize,
    },
    /// A Philox stream would need more than its `2^64` distinct positions.
    RandomCounterRangeOverflow {
        /// Stable stream-domain label.
        stream: &'static str,
        /// Exact requested block count, or `None` if even the `u128` product
        /// was not representable.
        required_blocks: Option<u128>,
    },
    /// The shared fs-la dense-Jacobi authority refused this decision dimension.
    EigensolverAdmission {
        /// Exact dependency refusal; no duplicated capability formula is used.
        error: fs_la::eigen::JacobiEighAdmissionError,
    },
    /// Checked scheduler accounting failed at a named restart boundary.
    ArithmeticOverflow {
        /// Restart ordinal being prepared or finalized.
        restart: u64,
        /// Stable formula label.
        what: &'static str,
    },
    /// A finite admitted root plus a finite perturbation produced a non-finite
    /// restart coordinate; the affected restart was not invoked.
    GeneratedStartNonFinite {
        /// Restart ordinal.
        restart: u64,
        /// Coordinate index.
        component: usize,
        /// Exact IEEE-754 payload.
        bits: u64,
    },
    /// CMA generated a non-finite decision coordinate. Earlier callbacks or
    /// restarts may already have completed, but the affected candidate was not
    /// passed to the objective.
    GeneratedCandidateNonFinite {
        /// Restart ordinal.
        restart: u64,
        /// One-based CMA generation.
        generation: usize,
        /// Zero-based candidate within the population.
        candidate: usize,
        /// Coordinate index.
        component: usize,
        /// Exact IEEE-754 payload.
        bits: u64,
    },
    /// A CMA run violated the local hard budget assumed by the scheduler.
    InternalBudgetViolation {
        /// Restart ordinal.
        restart: u64,
        /// Callbacks reported by CMA.
        spent: usize,
        /// Local cap supplied to CMA.
        allocated: usize,
    },
    /// Aggregate trace accounting exceeded the admitted hard budget.
    InternalAggregateBudgetViolation {
        /// Restart ordinal whose report crossed the boundary.
        restart: u64,
        /// Aggregate trace end after the restart.
        total_spent: usize,
        /// Admitted aggregate cap.
        total_budget: usize,
    },
    /// An admitted scheduler reached a state ruled out by its preflight and
    /// loop invariants.
    InternalInvariant {
        /// Stable invariant label.
        what: &'static str,
    },
    /// Generated evidence failed its own structural validator.
    GeneratedLedgerInvalid {
        /// Restart index when the invariant is local.
        restart: Option<usize>,
        /// Stable validator invariant.
        invariant: &'static str,
    },
}

impl core::fmt::Display for BipopError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match *self {
            Self::EmptyStart => {
                formatter.write_str("BIPOP start must contain at least one coordinate")
            }
            Self::NonFiniteStart { component, bits } => write!(
                formatter,
                "BIPOP start component {component} is non-finite (bits 0x{bits:016x})"
            ),
            Self::InvalidSigma { bits } => write!(
                formatter,
                "BIPOP sigma must be finite and strictly positive (bits 0x{bits:016x})"
            ),
            Self::NonFiniteTarget { bits } => write!(
                formatter,
                "BIPOP target must be finite when present (bits 0x{bits:016x}); use None for no target"
            ),
            Self::ZeroBudget => formatter.write_str("BIPOP total callback budget must be positive"),
            Self::DimensionOverflow { what } => {
                write!(formatter, "BIPOP dimension formula `{what}` overflowed")
            }
            Self::RestartOrdinalOverflow { total_budget } => write!(
                formatter,
                "BIPOP restart envelope for budget {total_budget} cannot be represented by the u64 ordinal"
            ),
            Self::SeedRangeOverflow {
                seed,
                max_restart_ordinal,
            } => write!(
                formatter,
                "BIPOP seed range from root {seed} through restart {max_restart_ordinal} would wrap"
            ),
            Self::PopulationOverflow { large_run } => write!(
                formatter,
                "BIPOP population ladder overflowed at large run {large_run}"
            ),
            Self::LocalBudgetOverflow { lambda } => write!(
                formatter,
                "BIPOP local budget formula overflowed for population {lambda}"
            ),
            Self::RandomCounterRangeOverflow {
                stream,
                required_blocks,
            } => match required_blocks {
                Some(required_blocks) => write!(
                    formatter,
                    "BIPOP {stream} stream requires {required_blocks} Philox blocks, exceeding the 2^64-position counter domain"
                ),
                None => write!(
                    formatter,
                    "BIPOP {stream} stream block count overflowed the u128 admission accumulator"
                ),
            },
            Self::EigensolverAdmission { error } => {
                write!(formatter, "BIPOP eigensolver admission refused: {error}")
            }
            Self::ArithmeticOverflow { restart, what } => write!(
                formatter,
                "BIPOP scheduler formula `{what}` overflowed at restart {restart}"
            ),
            Self::GeneratedStartNonFinite {
                restart,
                component,
                bits,
            } => write!(
                formatter,
                "BIPOP restart {restart} start component {component} became non-finite (bits 0x{bits:016x})"
            ),
            Self::GeneratedCandidateNonFinite {
                restart,
                generation,
                candidate,
                component,
                bits,
            } => write!(
                formatter,
                "BIPOP restart {restart} generation {generation} candidate {candidate} component {component} became non-finite (bits 0x{bits:016x})"
            ),
            Self::InternalBudgetViolation {
                restart,
                spent,
                allocated,
            } => write!(
                formatter,
                "BIPOP restart {restart} spent {spent} callbacks under local cap {allocated}"
            ),
            Self::InternalAggregateBudgetViolation {
                restart,
                total_spent,
                total_budget,
            } => write!(
                formatter,
                "BIPOP restart {restart} advanced aggregate callbacks to {total_spent} under hard cap {total_budget}"
            ),
            Self::InternalInvariant { what } => {
                write!(formatter, "BIPOP internal invariant failed: {what}")
            }
            Self::GeneratedLedgerInvalid { restart, invariant } => match restart {
                Some(restart) => write!(
                    formatter,
                    "generated BIPOP restart {restart} violates {invariant}"
                ),
                None => write!(formatter, "generated BIPOP ledger violates {invariant}"),
            },
        }
    }
}

impl std::error::Error for BipopError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::EigensolverAdmission { error } => Some(error),
            _ => None,
        }
    }
}

fn checked_random_counter_blocks(
    stream: &'static str,
    factors: &[u128],
) -> Result<u128, BipopError> {
    let required_blocks = factors
        .iter()
        .try_fold(1u128, |product, factor| product.checked_mul(*factor));
    let Some(required_blocks) = required_blocks else {
        return Err(BipopError::RandomCounterRangeOverflow {
            stream,
            required_blocks: None,
        });
    };
    if required_blocks > FS_RAND_STREAM_COUNTER_CARDINALITY {
        return Err(BipopError::RandomCounterRangeOverflow {
            stream,
            required_blocks: Some(required_blocks),
        });
    }
    Ok(required_blocks)
}

fn restart_stream_block_bound(
    dimension: usize,
    max_restart_ordinal: u64,
) -> Result<u128, BipopError> {
    checked_random_counter_blocks(
        "restart-perturbation",
        &[2, dimension as u128, u128::from(max_restart_ordinal)],
    )
}

fn cma_stream_block_bound(
    dimension: usize,
    lambda: usize,
    max_generations: usize,
) -> Result<u128, BipopError> {
    checked_random_counter_blocks(
        "CMA-candidate",
        &[
            2,
            dimension as u128,
            lambda as u128,
            max_generations as u128,
        ],
    )
}

fn checked_dense_matrix_allocation(max_matrix_entries: usize) -> Result<usize, BipopError> {
    let bytes = max_matrix_entries
        .checked_mul(core::mem::size_of::<f64>())
        .ok_or(BipopError::DimensionOverflow {
            what: "dense covariance bytes",
        })?;
    if bytes > isize::MAX as usize {
        return Err(BipopError::DimensionOverflow {
            what: "dense covariance address space",
        });
    }
    Ok(bytes)
}

fn scheduler_max_restart_ordinal(
    base_lambda: usize,
    total_budget: usize,
) -> Result<u64, BipopError> {
    let preterminal_large_spend = (base_lambda as u128)
        .checked_mul(BIPOP_FULL_ENVELOPE_GENERATIONS as u128)
        .and_then(|spend| spend.checked_mul(BIPOP_PRE_TERMINAL_LARGE_RUN_SCALE_SUM as u128))
        .and_then(|spend| spend.checked_add(u128::from(BIPOP_LARGE_RUN_CAP)))
        .ok_or(BipopError::RestartOrdinalOverflow { total_budget })?;
    let continuing_small_spend = (base_lambda as u128)
        .checked_add(1)
        .ok_or(BipopError::RestartOrdinalOverflow { total_budget })?;
    let max_full_generation_small_count = preterminal_large_spend
        .checked_add(continuing_small_spend - 1)
        .map(|numerator| numerator / continuing_small_spend)
        .ok_or(BipopError::RestartOrdinalOverflow { total_budget })?;
    let max_small_restart_count = max_full_generation_small_count
        .checked_add(base_lambda as u128)
        .ok_or(BipopError::RestartOrdinalOverflow { total_budget })?;
    let scheduler_max_restart_ordinal = u128::from(BIPOP_LARGE_RUN_CAP)
        .checked_add(max_small_restart_count)
        .ok_or(BipopError::RestartOrdinalOverflow { total_budget })?;
    let budget_max_restart_ordinal = (total_budget - 1) as u128;
    u64::try_from(budget_max_restart_ordinal.min(scheduler_max_restart_ordinal))
        .map_err(|_| BipopError::RestartOrdinalOverflow { total_budget })
}

/// Validate BIPOP inputs and every conservative arithmetic/storage envelope
/// before allocating scheduler state or invoking the objective.
///
/// `target = None` disables target stopping. This avoids using a non-finite
/// floating-point sentinel on the checked API.
///
/// # Errors
/// Returns [`BipopError`] for malformed input or an unrepresentable envelope.
pub fn admit_bipop(
    x0: &[f64],
    sigma0: f64,
    total_budget: usize,
    target: Option<f64>,
    seed: u64,
) -> Result<BipopAdmission, BipopError> {
    if x0.is_empty() {
        return Err(BipopError::EmptyStart);
    }
    for (component, value) in x0.iter().enumerate() {
        if !value.is_finite() {
            return Err(BipopError::NonFiniteStart {
                component,
                bits: value.to_bits(),
            });
        }
    }
    if !sigma0.is_finite() || sigma0 <= 0.0 {
        return Err(BipopError::InvalidSigma {
            bits: sigma0.to_bits(),
        });
    }
    if let Some(target) = target {
        if !target.is_finite() {
            return Err(BipopError::NonFiniteTarget {
                bits: target.to_bits(),
            });
        }
    }
    if total_budget == 0 {
        return Err(BipopError::ZeroBudget);
    }

    let dimension = x0.len();
    let lambda_offset = (3.0 * fs_math::det::ln(dimension as f64)).floor();
    if !lambda_offset.is_finite() || lambda_offset < 0.0 || lambda_offset > usize::MAX as f64 {
        return Err(BipopError::DimensionOverflow {
            what: "base population",
        });
    }
    let base_lambda =
        4usize
            .checked_add(lambda_offset as usize)
            .ok_or(BipopError::DimensionOverflow {
                what: "base population",
            })?;
    let jacobi_admission = if total_budget > base_lambda {
        Some(
            fs_la::eigen::admit_jacobi_eigh(dimension)
                .map_err(|error| BipopError::EigensolverAdmission { error })?,
        )
    } else {
        None
    };

    let max_restart_ordinal = scheduler_max_restart_ordinal(base_lambda, total_budget)?;
    let seed_delta = max_restart_ordinal
        .checked_mul(BIPOP_RESTART_SEED_STRIDE)
        .ok_or(BipopError::SeedRangeOverflow {
            seed,
            max_restart_ordinal,
        })?;
    seed.checked_add(seed_delta)
        .ok_or(BipopError::SeedRangeOverflow {
            seed,
            max_restart_ordinal,
        })?;

    let budget_ladder_cap = u32::try_from(max_restart_ordinal)
        .unwrap_or(u32::MAX)
        .min(BIPOP_LARGE_RUN_CAP);
    let scale = 1usize
        .checked_shl(budget_ladder_cap)
        .ok_or(BipopError::PopulationOverflow {
            large_run: budget_ladder_cap,
        })?;
    let max_large_lambda =
        base_lambda
            .checked_mul(scale)
            .ok_or(BipopError::PopulationOverflow {
                large_run: budget_ladder_cap,
            })?;
    let max_local_budget = max_large_lambda
        .checked_mul(BIPOP_GENERATIONS_PER_RESTART)
        .ok_or(BipopError::LocalBudgetOverflow {
            lambda: max_large_lambda,
        })?;
    let max_matrix_entries =
        dimension
            .checked_mul(dimension)
            .ok_or(BipopError::DimensionOverflow {
                what: "dense covariance entries",
            })?;
    let max_population_entries =
        dimension
            .checked_mul(max_large_lambda)
            .ok_or(BipopError::DimensionOverflow {
                what: "population coordinate entries",
            })?;
    checked_dense_matrix_allocation(max_matrix_entries)?;
    max_population_entries
        .checked_mul(core::mem::size_of::<f64>())
        .ok_or(BipopError::DimensionOverflow {
            what: "population coordinate bytes",
        })?;

    let max_restart_stream_blocks = restart_stream_block_bound(dimension, max_restart_ordinal)?;
    let mut max_cma_stream_blocks = 0u128;
    for large_run in 0..=budget_ladder_cap {
        let scale = 1usize
            .checked_shl(large_run)
            .ok_or(BipopError::PopulationOverflow { large_run })?;
        let lambda = base_lambda
            .checked_mul(scale)
            .ok_or(BipopError::PopulationOverflow { large_run })?;
        let local_envelope = lambda
            .checked_mul(BIPOP_GENERATIONS_PER_RESTART)
            .ok_or(BipopError::LocalBudgetOverflow { lambda })?;
        let allocated_budget = local_envelope.min(total_budget);
        // CMA spends one callback at the start, then admits only complete
        // populations: `1 + generations * lambda <= allocated_budget`.
        let max_generations = (allocated_budget - 1) / lambda;
        let blocks = cma_stream_block_bound(dimension, lambda, max_generations)?;
        max_cma_stream_blocks = max_cma_stream_blocks.max(blocks);
    }

    Ok(BipopAdmission {
        schema_version: BIPOP_ADMISSION_SCHEMA_VERSION,
        stream_semantics_version: fs_rand::STREAM_SEMANTICS_VERSION,
        jacobi_admission,
        dimension,
        total_budget,
        base_lambda,
        max_large_lambda,
        max_local_budget,
        max_restart_ordinal,
        max_matrix_entries,
        max_population_entries,
        max_restart_stream_blocks,
        max_cma_stream_blocks,
    })
}

/// Tunables (defaults follow Hansen's standard settings).
#[derive(Debug, Clone)]
pub struct CmaParams {
    /// Population size λ (default 4 + ⌊3·ln n⌋).
    pub lambda: usize,
    /// Initial step size σ₀.
    pub sigma0: f64,
    /// Evaluation budget.
    pub max_evals: usize,
    /// Target objective value (stop when reached).
    pub f_target: f64,
    /// Generations between eigendecompositions (SPD refresh cadence).
    pub eigen_interval: usize,
}

impl CmaParams {
    /// Standard defaults for dimension `n`.
    #[must_use]
    pub fn standard(n: usize, sigma0: f64, max_evals: usize, f_target: f64) -> CmaParams {
        let lambda = 4 + (3.0 * fs_math::det::ln(n as f64)).floor() as usize;
        CmaParams {
            lambda,
            sigma0,
            max_evals,
            f_target,
            eigen_interval: 1,
        }
    }
}

/// Why one CMA run stopped.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CmaStopReason {
    /// The requested objective target was reached.
    TargetReached,
    /// The local budget could not admit another complete population.
    BudgetExhausted,
    /// TolX/TolFun stopped a run.
    Stagnated,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct CmaGeneratedPointError {
    generation: usize,
    candidate: usize,
    component: usize,
    bits: u64,
}

/// Run evidence.
#[derive(Debug, Clone)]
pub struct CmaReport {
    /// Best point found.
    pub x_best: Vec<f64>,
    /// Best objective value.
    pub f_best: f64,
    /// Objective evaluations consumed.
    pub evals: usize,
    /// Generations run.
    pub generations: usize,
    /// Whether `f_target` was reached.
    pub converged: bool,
    /// Final step size (diagnostic).
    pub sigma: f64,
}

/// Full-covariance CMA-ES from `x0`. Deterministic per `seed`.
///
/// # Panics
/// Panics before the affected callback if the initial or a generated decision
/// point contains a non-finite coordinate. The fallible BIPOP surface projects
/// the same guard into [`BipopError::GeneratedCandidateNonFinite`].
#[must_use]
pub fn cmaes<F: FnMut(&[f64]) -> f64>(
    f: &mut F,
    x0: &[f64],
    params: &CmaParams,
    seed: u64,
) -> CmaReport {
    cmaes_with_stop(f, x0, params, seed).0
}

#[allow(clippy::too_many_lines)] // the algorithm is one coherent loop
fn cmaes_with_stop<F: FnMut(&[f64]) -> f64>(
    f: &mut F,
    x0: &[f64],
    params: &CmaParams,
    seed: u64,
) -> (CmaReport, CmaStopReason) {
    cmaes_with_stop_target(f, x0, params, seed, Some(params.f_target)).unwrap_or_else(|error| {
        panic!(
            "CMA generated a non-finite query at generation {} candidate {} component {} (bits 0x{:016x})",
            error.generation, error.candidate, error.component, error.bits
        )
    })
}

#[allow(clippy::too_many_lines)] // the algorithm is one coherent loop
fn cmaes_with_stop_target<F: FnMut(&[f64]) -> f64>(
    f: &mut F,
    x0: &[f64],
    params: &CmaParams,
    seed: u64,
    target: Option<f64>,
) -> Result<(CmaReport, CmaStopReason), CmaGeneratedPointError> {
    let n = x0.len();
    assert!(n >= 1, "dimension must be positive");
    if let Some((component, value)) = x0.iter().enumerate().find(|(_, value)| !value.is_finite()) {
        return Err(CmaGeneratedPointError {
            generation: 0,
            candidate: 0,
            component,
            bits: value.to_bits(),
        });
    }
    let lambda = params.lambda.max(4);
    let mu = lambda / 2;
    // Log-rank recombination weights (Hansen standard).
    let raw: Vec<f64> = (0..mu)
        .map(|i| {
            fs_math::det::ln(f64::midpoint(lambda as f64, 1.0)) - fs_math::det::ln(i as f64 + 1.0)
        })
        .collect();
    let wsum: f64 = raw.iter().sum();
    let weights: Vec<f64> = raw.iter().map(|w| w / wsum).collect();
    let mu_eff = 1.0 / weights.iter().map(|w| w * w).sum::<f64>();
    let nf = n as f64;
    // Standard strategy parameters (the IGO/natural-gradient couplings).
    let cc = (4.0 + mu_eff / nf) / (nf + 4.0 + 2.0 * mu_eff / nf);
    let cs = (mu_eff + 2.0) / (nf + mu_eff + 5.0);
    let c1 = 2.0 / ((nf + 1.3) * (nf + 1.3) + mu_eff);
    let cmu =
        (1.0 - c1).min(2.0 * (mu_eff - 2.0 + 1.0 / mu_eff) / ((nf + 2.0) * (nf + 2.0) + mu_eff));
    let damps = 1.0 + 2.0 * (fs_math::det::sqrt((mu_eff - 1.0) / (nf + 1.0)) - 1.0).max(0.0) + cs;
    // E‖N(0,I)‖ (Hansen's approximation).
    let chi_n = fs_math::det::sqrt(nf) * (1.0 - 1.0 / (4.0 * nf) + 1.0 / (21.0 * nf * nf));

    let mut mean = x0.to_vec();
    let mut sigma = params.sigma0;
    let mut cov = vec![0.0f64; n * n];
    for i in 0..n {
        cov[i * n + i] = 1.0;
    }
    let mut p_c = vec![0.0f64; n];
    let mut p_s = vec![0.0f64; n];
    // Eigen state: C = B·diag(d²)·Bᵀ; sqrt factors refreshed on cadence.
    let mut b_mat = cov.clone();
    let mut d_sqrt = vec![1.0f64; n];
    let mut stream = StreamKey {
        seed,
        kernel: K_CMA,
        tile: 0,
    }
    .stream();

    let mut x_best = mean.clone();
    let mut f_best = f(&mean);
    let mut evals = 1usize;
    let mut generations = 0usize;
    if target.is_some_and(|target| f_best <= target) {
        return Ok((
            CmaReport {
                x_best,
                f_best,
                evals,
                generations,
                converged: true,
                sigma,
            },
            CmaStopReason::TargetReached,
        ));
    }
    let mut stop_reason = CmaStopReason::BudgetExhausted;
    // TolFun stagnation: generations since a meaningful f_best improvement.
    let mut gens_since_improve = 0usize;

    let mut zs: Vec<Vec<f64>> = vec![vec![0.0; n]; lambda];
    let mut ys: Vec<Vec<f64>> = vec![vec![0.0; n]; lambda];
    let mut fitness: Vec<f64> = vec![0.0; lambda];

    while evals
        .checked_add(lambda)
        .is_some_and(|next_generation| next_generation <= params.max_evals)
    {
        generations += 1;
        // Refresh eigendecomposition on cadence (SPD maintenance).
        if generations % params.eigen_interval.max(1) == 1 || params.eigen_interval <= 1 {
            // Symmetrize (roundoff hygiene) then eigh; floor eigenvalues.
            for i in 0..n {
                for j in i + 1..n {
                    let avg = f64::midpoint(cov[i * n + j], cov[j * n + i]);
                    cov[i * n + j] = avg;
                    cov[j * n + i] = avg;
                }
            }
            let (vals, vecs) = jacobi_eigh(&cov, n);
            let vmax = vals.last().copied().unwrap_or(1.0).max(f64::MIN_POSITIVE);
            for (k, &v) in vals.iter().enumerate() {
                d_sqrt[k] = fs_math::det::sqrt(v.max(1e-14 * vmax));
            }
            b_mat.copy_from_slice(&vecs);
        }
        // Sample λ candidates: x = m + σ·B·diag(d)·z.
        for (k, z) in zs.iter_mut().enumerate() {
            for zi in z.iter_mut() {
                *zi = stream.next_normal();
            }
            let y = &mut ys[k];
            for i in 0..n {
                let mut acc = 0.0f64;
                for j in 0..n {
                    acc = (b_mat[i * n + j] * d_sqrt[j]).mul_add(z[j], acc);
                }
                y[i] = acc;
            }
            let x: Vec<f64> = mean
                .iter()
                .zip(y.iter())
                .map(|(m, yi)| sigma.mul_add(*yi, *m))
                .collect();
            if let Some((component, value)) =
                x.iter().enumerate().find(|(_, value)| !value.is_finite())
            {
                return Err(CmaGeneratedPointError {
                    generation: generations,
                    candidate: k,
                    component,
                    bits: value.to_bits(),
                });
            }
            fitness[k] = f(&x);
            evals += 1;
            if fitness[k] < f_best {
                if f_best - fitness[k] > 1e-12 * (1.0 + f_best.abs()) {
                    gens_since_improve = 0;
                }
                f_best = fitness[k];
                x_best = x;
            }
        }
        gens_since_improve += 1;
        if target.is_some_and(|target| f_best <= target) {
            return Ok((
                CmaReport {
                    x_best,
                    f_best,
                    evals,
                    generations,
                    converged: true,
                    sigma,
                },
                CmaStopReason::TargetReached,
            ));
        }
        // Rank (total_cmp, lowest index on ties — P2).
        let mut order: Vec<usize> = (0..lambda).collect();
        order.sort_by(|&a, &b| fitness[a].total_cmp(&fitness[b]).then(a.cmp(&b)));
        // Weighted recombination in y-space.
        let mut y_w = vec![0.0f64; n];
        for (w, &idx) in weights.iter().zip(&order) {
            for i in 0..n {
                y_w[i] = w.mul_add(ys[idx][i], y_w[i]);
            }
        }
        // Mean update.
        for i in 0..n {
            mean[i] = sigma.mul_add(y_w[i], mean[i]);
        }
        // CSA path: p_s ← (1−cs)p_s + √(cs(2−cs)µeff)·C^{−1/2}·y_w,
        // with C^{−1/2} = B·diag(1/d)·Bᵀ.
        let mut c_inv_half_yw = vec![0.0f64; n];
        for i in 0..n {
            // t = Bᵀ y_w
            let mut acc = 0.0f64;
            for j in 0..n {
                acc = b_mat[j * n + i].mul_add(y_w[j], acc);
            }
            c_inv_half_yw[i] = acc / d_sqrt[i];
        }
        let mut tmp = vec![0.0f64; n];
        for i in 0..n {
            let mut acc = 0.0f64;
            for j in 0..n {
                acc = b_mat[i * n + j].mul_add(c_inv_half_yw[j], acc);
            }
            tmp[i] = acc;
        }
        let csn = fs_math::det::sqrt(cs * (2.0 - cs) * mu_eff);
        for i in 0..n {
            p_s[i] = (1.0 - cs).mul_add(p_s[i], csn * tmp[i]);
        }
        let ps_norm = fs_math::det::sqrt(p_s.iter().map(|t| t * t).sum::<f64>());
        // Step-size update (the natural-gradient-consistent coupling).
        sigma *= fs_math::det::exp((cs / damps) * (ps_norm / chi_n - 1.0));
        // STAGNATION STOP: once the search distribution has collapsed
        // (σ·√λmax(C) negligible vs σ₀) the run is dead — keep sampling
        // and it just burns budget polishing whatever basin it's in.
        // BIPOP's restart ladder DEPENDS on dead runs terminating
        // (measured during bring-up: without this, failed runs consumed
        // their entire 120k budget at f ≈ 1 on rastrigin).
        let spread = sigma * d_sqrt.iter().fold(0.0f64, |m, &d| m.max(d));
        if spread < 1e-12 * params.sigma0 || gens_since_improve > 120 {
            // TolX OR TolFun: σ-collapse alone fires too slowly inside a
            // per-run budget (measured: a λ=150 local-basin run burned
            // 120k evals with f stalled for hundreds of generations) —
            // the f-stall criterion is what actually frees the budget.
            stop_reason = CmaStopReason::Stagnated;
            break;
        }
        // Rank-1 path with stall indicator h_σ.
        let h_sig = ps_norm
            / fs_math::det::sqrt(
                1.0 - fs_math::det::powi(
                    1.0 - cs,
                    2 * i32::try_from(generations.min(100_000)).expect("generation count"),
                ),
            )
            < (1.4 + 2.0 / (nf + 1.0)) * chi_n;
        let ccn = fs_math::det::sqrt(cc * (2.0 - cc) * mu_eff);
        for i in 0..n {
            let h = if h_sig { ccn * y_w[i] } else { 0.0 };
            p_c[i] = (1.0 - cc).mul_add(p_c[i], h);
        }
        // Covariance update: rank-1 + rank-µ.
        let delta_h = if h_sig { 0.0 } else { cc * (2.0 - cc) };
        for i in 0..n {
            for j in 0..n {
                let mut rank_mu = 0.0f64;
                for (w, &idx) in weights.iter().zip(&order) {
                    rank_mu = (w * ys[idx][i]).mul_add(ys[idx][j], rank_mu);
                }
                let rank1 = p_c[i] * p_c[j];
                cov[i * n + j] = (1.0 - c1 - cmu).mul_add(
                    cov[i * n + j],
                    c1.mul_add(rank1 + delta_h * cov[i * n + j], cmu * rank_mu),
                );
            }
        }
    }
    Ok((
        CmaReport {
            x_best,
            f_best,
            evals,
            generations,
            converged: false,
            sigma,
        },
        stop_reason,
    ))
}

/// Which BIPOP budget lane launched a restart.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BipopLane {
    /// The doubling population ladder.
    Large,
    /// The base-population interleave.
    Small,
}

/// One immutable, versioned BIPOP restart receipt.
///
/// Point and objective values retain their exact `f64` bits. The aggregate
/// trace interval is half-open and counts objective invocations even though
/// this first ledger tranche does not retain the objective payloads themselves.
#[derive(Debug, Clone)]
pub struct BipopRestartRecord {
    schema_version: u32,
    ordinal: u64,
    lane: BipopLane,
    lambda: usize,
    allocated_budget: usize,
    seed: u64,
    start: Vec<f64>,
    trace_start: usize,
    trace_end: usize,
    stop_reason: CmaStopReason,
    report: CmaReport,
}

impl BipopRestartRecord {
    /// Restart-record schema version.
    #[must_use]
    pub fn schema_version(&self) -> u32 {
        self.schema_version
    }

    /// Zero-based restart ordinal.
    #[must_use]
    pub fn ordinal(&self) -> u64 {
        self.ordinal
    }

    /// Large or small BIPOP budget lane.
    #[must_use]
    pub fn lane(&self) -> BipopLane {
        self.lane
    }

    /// Population size used by this restart.
    #[must_use]
    pub fn lambda(&self) -> usize {
        self.lambda
    }

    /// Local evaluation cap assigned to this restart.
    #[must_use]
    pub fn allocated_budget(&self) -> usize {
        self.allocated_budget
    }

    /// CMA stream seed derived for this restart.
    #[must_use]
    pub fn seed(&self) -> u64 {
        self.seed
    }

    /// Exact start point supplied to this restart.
    #[must_use]
    pub fn start(&self) -> &[f64] {
        &self.start
    }

    /// Start of this restart's half-open aggregate evaluation interval.
    #[must_use]
    pub fn trace_start(&self) -> usize {
        self.trace_start
    }

    /// End of this restart's half-open aggregate evaluation interval.
    #[must_use]
    pub fn trace_end(&self) -> usize {
        self.trace_end
    }

    /// Causal terminal classification retained from the CMA run.
    #[must_use]
    pub fn stop_reason(&self) -> CmaStopReason {
        self.stop_reason
    }

    /// Complete CMA result for this restart.
    #[must_use]
    pub fn report(&self) -> &CmaReport {
        &self.report
    }
}

/// Structured refusal from [`BipopReport::validate_ledger`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BipopLedgerError {
    restart: Option<usize>,
    invariant: &'static str,
}

impl BipopLedgerError {
    fn global(invariant: &'static str) -> Self {
        Self {
            restart: None,
            invariant,
        }
    }

    fn at(restart: usize, invariant: &'static str) -> Self {
        Self {
            restart: Some(restart),
            invariant,
        }
    }

    /// Restart index associated with the refusal, if it is local.
    #[must_use]
    pub fn restart(&self) -> Option<usize> {
        self.restart
    }

    /// Stable invariant name suitable for structured diagnostics.
    #[must_use]
    pub fn invariant(&self) -> &'static str {
        self.invariant
    }
}

impl core::fmt::Display for BipopLedgerError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self.restart {
            Some(restart) => write!(
                formatter,
                "BIPOP restart {restart} violates {}",
                self.invariant
            ),
            None => write!(formatter, "BIPOP ledger violates {}", self.invariant),
        }
    }
}

impl std::error::Error for BipopLedgerError {}

fn cma_reports_match_bits(left: &CmaReport, right: &CmaReport) -> bool {
    left.f_best.to_bits() == right.f_best.to_bits()
        && left.evals == right.evals
        && left.generations == right.generations
        && left.converged == right.converged
        && left.sigma.to_bits() == right.sigma.to_bits()
        && left.x_best.len() == right.x_best.len()
        && left
            .x_best
            .iter()
            .zip(&right.x_best)
            .all(|(left, right)| left.to_bits() == right.to_bits())
}

/// BIPOP restart evidence.
#[derive(Debug, Clone)]
pub struct BipopReport {
    /// Compatibility projection of [`Self::best_record`].
    pub best: CmaReport,
    /// Compatibility projection of every restart's population size.
    pub schedule: Vec<usize>,
    /// Compatibility projection of the terminal aggregate trace offset.
    pub total_evals: usize,
    total_budget: usize,
    records: Vec<BipopRestartRecord>,
    best_restart: usize,
}

impl BipopReport {
    /// Ordered immutable restart ledger.
    #[must_use]
    pub fn records(&self) -> &[BipopRestartRecord] {
        &self.records
    }

    /// Admitted hard aggregate callback budget that governs every record's
    /// local allocation and the ledger's terminal completeness.
    #[must_use]
    pub fn total_budget(&self) -> usize {
        self.total_budget
    }

    /// Index of the earliest restart attaining the best objective under
    /// `f64::total_cmp`.
    #[must_use]
    pub fn best_restart(&self) -> usize {
        self.best_restart
    }

    /// Named record from which [`Self::best`] is projected.
    #[must_use]
    pub fn best_record(&self) -> Option<&BipopRestartRecord> {
        self.records.get(self.best_restart)
    }

    /// Recheck the ordered ledger and every compatibility projection.
    ///
    /// This is a structural validator over retained evidence. It does not yet
    /// authenticate the first start, root seed, sigma, target, or callback
    /// semantics against an external input identity.
    ///
    /// # Errors
    /// Returns a [`BipopLedgerError`] naming the first deterministic invariant
    /// violation.
    #[allow(clippy::too_many_lines)] // one ordered pass mirrors the versioned record schema
    pub fn validate_ledger(&self) -> Result<(), BipopLedgerError> {
        let first = self
            .records
            .first()
            .ok_or_else(|| BipopLedgerError::global("nonempty"))?;
        if self.schedule.len() != self.records.len() {
            return Err(BipopLedgerError::global("schedule-length"));
        }
        if self.total_budget == 0 {
            return Err(BipopLedgerError::global("positive-total-budget"));
        }
        let base_lambda = first.lambda;
        let base_seed = first.seed;
        let point_dim = first.start.len();
        if point_dim == 0 {
            return Err(BipopLedgerError::at(0, "positive-point-dimension"));
        }
        let matrix_entries = point_dim
            .checked_mul(point_dim)
            .ok_or_else(|| BipopLedgerError::at(0, "dense-matrix-admission"))?;
        checked_dense_matrix_allocation(matrix_entries)
            .map_err(|_| BipopLedgerError::at(0, "dense-matrix-admission"))?;
        checked_random_counter_blocks(
            "restart-perturbation",
            &[2, point_dim as u128, (self.records.len() - 1) as u128],
        )
        .map_err(|_| BipopLedgerError::global("restart-counter-range"))?;
        let expected_base_lambda = 4 + (3.0 * fs_math::det::ln(point_dim as f64)).floor() as usize;
        if base_lambda != expected_base_lambda {
            return Err(BipopLedgerError::at(0, "base-population"));
        }
        // Mirror callback-free admission, not merely the realized trace. An
        // early target hit cannot retroactively authenticate a problem whose
        // hard budget could have entered the refused Jacobi dependency.
        if self.total_budget > expected_base_lambda {
            fs_la::eigen::admit_jacobi_eigh(point_dim)
                .map_err(|_| BipopLedgerError::at(0, "eigensolver-admission"))?;
        }
        let admitted_max_ordinal = scheduler_max_restart_ordinal(base_lambda, self.total_budget)
            .map_err(|_| BipopLedgerError::global("restart-envelope"))?;
        let last_ordinal = u64::try_from(self.records.len() - 1)
            .map_err(|_| BipopLedgerError::global("restart-envelope"))?;
        if last_ordinal > admitted_max_ordinal {
            return Err(BipopLedgerError::global("restart-envelope"));
        }
        let admitted_seed_delta = admitted_max_ordinal
            .checked_mul(BIPOP_RESTART_SEED_STRIDE)
            .ok_or_else(|| BipopLedgerError::global("admission-seed-range"))?;
        base_seed
            .checked_add(admitted_seed_delta)
            .ok_or_else(|| BipopLedgerError::global("admission-seed-range"))?;
        let mut cursor = 0usize;
        let mut large_budget_used = 0usize;
        let mut small_budget_used = 0usize;
        let mut large_runs = 0u32;

        for (index, record) in self.records.iter().enumerate() {
            // Production stops immediately after publishing rung eight, so no
            // tenth-large or trailing-small record can be authentic even if
            // its local lane arithmetic is otherwise self-consistent.
            if large_runs > BIPOP_LARGE_RUN_CAP {
                return Err(BipopLedgerError::at(index, "large-run-cap"));
            }
            if record.schema_version != BIPOP_RESTART_SCHEMA_VERSION {
                return Err(BipopLedgerError::at(index, "schema-version"));
            }
            let expected_ordinal = u64::try_from(index)
                .map_err(|_| BipopLedgerError::at(index, "ordinal-overflow"))?;
            if record.ordinal != expected_ordinal {
                return Err(BipopLedgerError::at(index, "ordinal"));
            }
            let expected_seed = expected_ordinal
                .checked_mul(BIPOP_RESTART_SEED_STRIDE)
                .and_then(|delta| base_seed.checked_add(delta))
                .ok_or_else(|| BipopLedgerError::at(index, "derived-seed-overflow"))?;
            if record.seed != expected_seed {
                return Err(BipopLedgerError::at(index, "derived-seed"));
            }
            if record.start.len() != point_dim || record.report.x_best.len() != point_dim {
                return Err(BipopLedgerError::at(index, "point-dimension"));
            }
            if record.start.iter().any(|value| !value.is_finite()) {
                return Err(BipopLedgerError::at(index, "finite-start"));
            }
            if record.report.x_best.iter().any(|value| !value.is_finite()) {
                return Err(BipopLedgerError::at(index, "finite-best-point"));
            }
            if !record.report.sigma.is_finite() || record.report.sigma.is_sign_negative() {
                return Err(BipopLedgerError::at(index, "finite-nonnegative-sigma"));
            }

            let expected_lane = if large_budget_used <= small_budget_used {
                BipopLane::Large
            } else {
                BipopLane::Small
            };
            if record.lane != expected_lane {
                return Err(BipopLedgerError::at(index, "lane-selection"));
            }
            let expected_lambda = match expected_lane {
                BipopLane::Large => 1usize
                    .checked_shl(large_runs)
                    .and_then(|scale| base_lambda.checked_mul(scale))
                    .ok_or_else(|| BipopLedgerError::at(index, "population-overflow"))?,
                BipopLane::Small => base_lambda,
            };
            if record.lambda != expected_lambda || self.schedule[index] != record.lambda {
                return Err(BipopLedgerError::at(index, "population-schedule"));
            }
            if record.trace_start != cursor {
                return Err(BipopLedgerError::at(index, "trace-start"));
            }
            let expected_end = cursor
                .checked_add(record.report.evals)
                .ok_or_else(|| BipopLedgerError::at(index, "trace-overflow"))?;
            if record.trace_end != expected_end {
                return Err(BipopLedgerError::at(index, "trace-end"));
            }
            if record.report.evals > record.allocated_budget {
                return Err(BipopLedgerError::at(index, "local-budget"));
            }
            if record.allocated_budget == 0 {
                return Err(BipopLedgerError::at(index, "positive-local-budget"));
            }
            let accounted_evals = record
                .report
                .generations
                .checked_mul(record.lambda)
                .and_then(|samples| samples.checked_add(1))
                .ok_or_else(|| BipopLedgerError::at(index, "evaluation-overflow"))?;
            if record.report.evals != accounted_evals {
                return Err(BipopLedgerError::at(index, "generation-accounting"));
            }
            match record.stop_reason {
                CmaStopReason::TargetReached => {
                    if !record.report.converged {
                        return Err(BipopLedgerError::at(index, "terminal-reason"));
                    }
                }
                CmaStopReason::BudgetExhausted => {
                    let next_generation = record
                        .report
                        .evals
                        .checked_add(record.lambda)
                        .ok_or_else(|| BipopLedgerError::at(index, "evaluation-overflow"))?;
                    if record.report.converged || next_generation <= record.allocated_budget {
                        return Err(BipopLedgerError::at(index, "terminal-reason"));
                    }
                }
                CmaStopReason::Stagnated => {
                    if record.report.converged || record.report.generations == 0 {
                        return Err(BipopLedgerError::at(index, "terminal-reason"));
                    }
                }
            }
            if index + 1 < self.records.len() && record.report.converged {
                return Err(BipopLedgerError::at(index, "continued-after-convergence"));
            }
            // Validate the report's intrinsic generation/terminal arithmetic
            // first so an overflowing forged report cannot hide behind a second
            // scheduler violation. Then bind the local cap to the production
            // `lambda * 250` envelope used by the admission theorem.
            let local_envelope = record
                .lambda
                .checked_mul(BIPOP_GENERATIONS_PER_RESTART)
                .ok_or_else(|| BipopLedgerError::at(index, "local-budget-overflow"))?;
            let remaining = self
                .total_budget
                .checked_sub(cursor)
                .ok_or_else(|| BipopLedgerError::at(index, "aggregate-budget"))?;
            let expected_allocated_budget = local_envelope.min(remaining);
            if record.allocated_budget != expected_allocated_budget {
                return Err(BipopLedgerError::at(index, "allocated-budget"));
            }
            cma_stream_block_bound(point_dim, record.lambda, record.report.generations)
                .map_err(|_| BipopLedgerError::at(index, "candidate-counter-range"))?;

            cursor = expected_end;
            match record.lane {
                BipopLane::Large => {
                    large_budget_used = large_budget_used
                        .checked_add(record.report.evals)
                        .ok_or_else(|| BipopLedgerError::at(index, "lane-budget-overflow"))?;
                    large_runs = large_runs
                        .checked_add(1)
                        .ok_or_else(|| BipopLedgerError::at(index, "large-run-overflow"))?;
                }
                BipopLane::Small => {
                    small_budget_used = small_budget_used
                        .checked_add(record.report.evals)
                        .ok_or_else(|| BipopLedgerError::at(index, "lane-budget-overflow"))?;
                }
            }
        }

        if cursor != self.total_evals {
            return Err(BipopLedgerError::global("total-evaluations"));
        }
        let last = self
            .records
            .last()
            .ok_or_else(|| BipopLedgerError::global("nonempty"))?;
        if cursor < self.total_budget && !last.report.converged && large_runs <= BIPOP_LARGE_RUN_CAP
        {
            return Err(BipopLedgerError::global("nonterminal-prefix"));
        }
        let mut expected_best = 0usize;
        for index in 1..self.records.len() {
            if self.records[index]
                .report
                .f_best
                .total_cmp(&self.records[expected_best].report.f_best)
                .is_lt()
            {
                expected_best = index;
            }
        }
        if self.best_restart != expected_best {
            return Err(BipopLedgerError::global("best-restart"));
        }
        if !cma_reports_match_bits(&self.best, &self.records[expected_best].report) {
            return Err(BipopLedgerError::global("best-projection"));
        }
        Ok(())
    }
}

/// Fallible BIPOP-CMA-ES under one admitted hard callback budget.
///
/// `target = None` disables target stopping. Every raw-input, conservative
/// seed-range, dimension, population, and local-budget formula is checked
/// before the first callback. Per-restart start generation and accounting are
/// checked again before that restart becomes visible to the callback.
///
/// # Errors
/// Returns [`BipopError`] without invoking `f` for raw-input or preflight
/// refusal. A later generated-start or generated-candidate refusal can follow
/// completed callbacks or earlier restarts, but the affected decision point is
/// never passed to the objective or published.
#[allow(clippy::too_many_lines)] // scheduler and record publication are one atomic state machine
pub fn try_bipop_cmaes<F: FnMut(&[f64]) -> f64>(
    f: &mut F,
    x0: &[f64],
    sigma0: f64,
    total_budget: usize,
    target: Option<f64>,
    seed: u64,
) -> Result<BipopReport, BipopError> {
    let admission = admit_bipop(x0, sigma0, total_budget, target, seed)?;
    run_admitted_bipop(f, x0, sigma0, target, seed, admission)
}

/// Legacy BIPOP spelling retained as a checked compatibility projection.
///
/// Finite targets map to `Some(target)`. Historical `-∞` means no target and
/// maps to `None`; all other malformed inputs now refuse before callbacks and
/// panic at this legacy boundary. New callers should use [`try_bipop_cmaes`]
/// for typed refusal.
#[must_use]
pub fn bipop_cmaes<F: FnMut(&[f64]) -> f64>(
    f: &mut F,
    x0: &[f64],
    sigma0: f64,
    total_budget: usize,
    f_target: f64,
    seed: u64,
) -> BipopReport {
    let target = if f_target.to_bits() == f64::NEG_INFINITY.to_bits() {
        None
    } else {
        Some(f_target)
    };
    try_bipop_cmaes(f, x0, sigma0, total_budget, target, seed)
        .unwrap_or_else(|error| panic!("BIPOP input or scheduler refused: {error}"))
}

#[allow(clippy::too_many_lines)] // scheduler and record publication are one atomic state machine
fn run_admitted_bipop<F: FnMut(&[f64]) -> f64>(
    f: &mut F,
    x0: &[f64],
    sigma0: f64,
    target: Option<f64>,
    seed: u64,
    admission: BipopAdmission,
) -> Result<BipopReport, BipopError> {
    let base_lambda = admission.base_lambda;
    let total_budget = admission.total_budget;
    let dimension = admission.dimension;
    let mut records: Vec<BipopRestartRecord> = Vec::new();
    let mut total_evals = 0usize;
    let mut best_restart: Option<usize> = None;
    let mut large_runs = 0u32;
    let mut restart = 0u64;
    let mut small_budget_used = 0usize;
    let mut large_budget_used = 0usize;
    let mut restart_stream_blocks_used = 0u128;
    // Deterministic restart-start perturbations (a restart from the SAME
    // point with a tiny sigma is just a polish run and cannot escape a
    // local basin — measured during bring-up on rastrigin).
    let mut restart_stream = StreamKey {
        seed,
        kernel: K_CMA,
        tile: 1,
    }
    .stream();
    while total_evals < total_budget {
        if restart > admission.max_restart_ordinal {
            return Err(BipopError::InternalInvariant {
                what: "restart ordinal exceeded its admitted scheduler bound",
            });
        }
        // BIPOP rule: run LARGE next if its cumulative budget lags.
        let run_large = large_budget_used <= small_budget_used;
        let lambda = if run_large {
            let scale = 1usize
                .checked_shl(large_runs)
                .ok_or(BipopError::PopulationOverflow {
                    large_run: large_runs,
                })?;
            base_lambda
                .checked_mul(scale)
                .ok_or(BipopError::PopulationOverflow {
                    large_run: large_runs,
                })?
        } else {
            base_lambda
        };
        // Per-run budget scales with the population (≈250 generations):
        // handing a small-λ run half the TOTAL budget just polishes one
        // local minimum expensively — the doubling ladder must be reached
        // (measured during bring-up on rastrigin).
        let local_envelope = lambda
            .checked_mul(BIPOP_GENERATIONS_PER_RESTART)
            .ok_or(BipopError::LocalBudgetOverflow { lambda })?;
        if lambda > admission.max_large_lambda || local_envelope > admission.max_local_budget {
            return Err(BipopError::InternalInvariant {
                what: "runtime population or local budget exceeded admission",
            });
        }
        let remaining =
            total_budget
                .checked_sub(total_evals)
                .ok_or(BipopError::ArithmeticOverflow {
                    restart,
                    what: "remaining aggregate budget",
                })?;
        let budget = local_envelope.min(remaining);
        if budget == 0 {
            return Err(BipopError::InternalInvariant {
                what: "a launched restart must receive at least one callback",
            });
        }
        let max_generations = (budget - 1) / lambda;
        let planned_cma_stream_blocks = cma_stream_block_bound(dimension, lambda, max_generations)?;
        if planned_cma_stream_blocks > admission.max_cma_stream_blocks {
            return Err(BipopError::InternalInvariant {
                what: "runtime CMA stream exceeded its admitted counter bound",
            });
        }
        // Reserve every scheduler addition against the full local cap before
        // this restart becomes visible to the callback. The post-run checks
        // below then narrow these envelopes to the actual callback count.
        let trace_start = total_evals;
        let trace_cap = trace_start
            .checked_add(budget)
            .ok_or(BipopError::ArithmeticOverflow {
                restart,
                what: "aggregate trace envelope",
            })?;
        if trace_cap > total_budget {
            return Err(BipopError::InternalAggregateBudgetViolation {
                restart,
                total_spent: trace_cap,
                total_budget,
            });
        }
        let (large_budget_cap, small_budget_cap) = if run_large {
            (
                large_budget_used
                    .checked_add(budget)
                    .ok_or(BipopError::ArithmeticOverflow {
                        restart,
                        what: "large-lane budget envelope",
                    })?,
                small_budget_used,
            )
        } else {
            (
                large_budget_used,
                small_budget_used
                    .checked_add(budget)
                    .ok_or(BipopError::ArithmeticOverflow {
                        restart,
                        what: "small-lane budget envelope",
                    })?,
            )
        };
        let next_restart = restart
            .checked_add(1)
            .ok_or(BipopError::ArithmeticOverflow {
                restart,
                what: "restart ordinal",
            })?;
        let large_runs_after = if run_large {
            large_runs
                .checked_add(1)
                .ok_or(BipopError::ArithmeticOverflow {
                    restart,
                    what: "large-run count",
                })?
        } else {
            large_runs
        };
        let params = CmaParams {
            lambda,
            sigma0,
            max_evals: budget,
            f_target: target.unwrap_or(f64::NEG_INFINITY),
            eigen_interval: 1,
        };
        // Restarts after the first launch from a perturbed start point.
        let start: Vec<f64> = if restart == 0 {
            x0.to_vec()
        } else {
            let blocks_for_start = restart_stream_block_bound(dimension, 1)?;
            let next_restart_stream_blocks = restart_stream_blocks_used
                .checked_add(blocks_for_start)
                .ok_or(BipopError::RandomCounterRangeOverflow {
                    stream: "restart-perturbation",
                    required_blocks: None,
                })?;
            if next_restart_stream_blocks > admission.max_restart_stream_blocks
                || next_restart_stream_blocks > FS_RAND_STREAM_COUNTER_CARDINALITY
            {
                return Err(BipopError::InternalInvariant {
                    what: "runtime restart stream exceeded its admitted counter bound",
                });
            }
            let start = x0
                .iter()
                .map(|&v| sigma0.mul_add(restart_stream.next_normal(), v))
                .collect();
            restart_stream_blocks_used = next_restart_stream_blocks;
            start
        };
        for (component, value) in start.iter().enumerate() {
            if !value.is_finite() {
                return Err(BipopError::GeneratedStartNonFinite {
                    restart,
                    component,
                    bits: value.to_bits(),
                });
            }
        }
        let derived_seed = restart
            .checked_mul(BIPOP_RESTART_SEED_STRIDE)
            .and_then(|delta| seed.checked_add(delta))
            .ok_or(BipopError::SeedRangeOverflow {
                seed,
                max_restart_ordinal: restart,
            })?;
        let (rep, stop_reason) = cmaes_with_stop_target(f, &start, &params, derived_seed, target)
            .map_err(|error| BipopError::GeneratedCandidateNonFinite {
            restart,
            generation: error.generation,
            candidate: error.candidate,
            component: error.component,
            bits: error.bits,
        })?;
        let actual_cma_stream_blocks = cma_stream_block_bound(dimension, lambda, rep.generations)?;
        if actual_cma_stream_blocks > planned_cma_stream_blocks {
            return Err(BipopError::InternalInvariant {
                what: "actual CMA stream consumption exceeded its preflight",
            });
        }
        if rep.evals > budget {
            return Err(BipopError::InternalBudgetViolation {
                restart,
                spent: rep.evals,
                allocated: budget,
            });
        }
        let trace_end =
            trace_start
                .checked_add(rep.evals)
                .ok_or(BipopError::ArithmeticOverflow {
                    restart,
                    what: "aggregate trace end",
                })?;
        if trace_end > total_budget {
            return Err(BipopError::InternalAggregateBudgetViolation {
                restart,
                total_spent: trace_end,
                total_budget,
            });
        }
        if trace_end > trace_cap {
            return Err(BipopError::InternalInvariant {
                what: "actual trace exceeded its admitted envelope",
            });
        }
        let record_index = records.len();
        let is_better = best_restart.is_none_or(|best_index| {
            rep.f_best
                .total_cmp(&records[best_index].report.f_best)
                .is_lt()
        });
        records.push(BipopRestartRecord {
            schema_version: BIPOP_RESTART_SCHEMA_VERSION,
            ordinal: restart,
            lane: if run_large {
                BipopLane::Large
            } else {
                BipopLane::Small
            },
            lambda,
            allocated_budget: budget,
            seed: derived_seed,
            start,
            trace_start,
            trace_end,
            stop_reason,
            report: rep,
        });
        if is_better {
            best_restart = Some(record_index);
        }
        total_evals = trace_end;
        if run_large {
            large_budget_used = large_budget_used
                .checked_add(records[record_index].report.evals)
                .ok_or(BipopError::ArithmeticOverflow {
                    restart,
                    what: "large-lane evaluation total",
                })?;
            if large_budget_used > large_budget_cap {
                return Err(BipopError::InternalInvariant {
                    what: "actual large-lane total exceeded its admitted envelope",
                });
            }
        } else {
            small_budget_used = small_budget_used
                .checked_add(records[record_index].report.evals)
                .ok_or(BipopError::ArithmeticOverflow {
                    restart,
                    what: "small-lane evaluation total",
                })?;
            if small_budget_used > small_budget_cap {
                return Err(BipopError::InternalInvariant {
                    what: "actual small-lane total exceeded its admitted envelope",
                });
            }
        }
        large_runs = large_runs_after;
        if records[record_index].report.converged {
            break;
        }
        restart = next_restart;
        if large_runs > BIPOP_LARGE_RUN_CAP {
            // Cap the LADDER, not total restarts: small runs are cheap
            // and interleave freely; counting them against the cap
            // stalled the ladder at λ ≈ 64 (measured during bring-up).
            break;
        }
    }
    let best_restart = best_restart.ok_or(BipopError::InternalInvariant {
        what: "positive admitted budget must launch one restart",
    })?;
    let schedule = records.iter().map(BipopRestartRecord::lambda).collect();
    let best = records[best_restart].report.clone();
    let report = BipopReport {
        best,
        schedule,
        total_evals,
        total_budget,
        records,
        best_restart,
    };
    report
        .validate_ledger()
        .map_err(|error| BipopError::GeneratedLedgerInvalid {
            restart: error.restart,
            invariant: error.invariant,
        })?;
    Ok(report)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// G0: both random domains count every semantic axis, accept exactly the
    /// full `2^64` Philox coordinate set once, and refuse the first reuse. The
    /// small exact KATs are mutation guards for the Box-Muller factor, decision
    /// dimension, population/generation product, and restart ordinal.
    #[test]
    fn random_counter_bounds_count_every_axis_and_refuse_reuse() {
        assert_eq!(
            restart_stream_block_bound(3, 5).expect("small restart bound"),
            2 * 3 * 5
        );
        assert_eq!(
            cma_stream_block_bound(3, 7, 11).expect("small CMA bound"),
            2 * 3 * 7 * 11
        );
        assert_eq!(
            checked_random_counter_blocks(
                "exact-cardinality",
                &[2, FS_RAND_STREAM_COUNTER_CARDINALITY / 2],
            )
            .expect("every counter coordinate may be consumed once"),
            FS_RAND_STREAM_COUNTER_CARDINALITY
        );

        let required_blocks = FS_RAND_STREAM_COUNTER_CARDINALITY + 2;
        assert_eq!(
            checked_random_counter_blocks(
                "first-reuse",
                &[2, FS_RAND_STREAM_COUNTER_CARDINALITY / 2 + 1],
            )
            .expect_err("the first coordinate reuse must refuse"),
            BipopError::RandomCounterRangeOverflow {
                stream: "first-reuse",
                required_blocks: Some(required_blocks),
            }
        );
        assert_eq!(
            checked_random_counter_blocks("u128-overflow", &[u128::MAX, 2])
                .expect_err("the proof accumulator must also fail closed"),
            BipopError::RandomCounterRangeOverflow {
                stream: "u128-overflow",
                required_blocks: None,
            }
        );
    }

    /// G0: checked multiplication alone is not sufficient for a Rust `Vec`.
    /// One dense covariance allocation must also fit the target's `isize`
    /// address-difference domain before any objective callback is possible.
    #[test]
    fn dense_covariance_allocation_refuses_first_unaddressable_size() {
        let element_bytes = core::mem::size_of::<f64>();
        let max_entries = (isize::MAX as usize) / element_bytes;
        let exact_bytes = checked_dense_matrix_allocation(max_entries)
            .expect("the last whole f64 allocation below isize::MAX admits");
        assert!(exact_bytes <= isize::MAX as usize);

        let first_unaddressable = max_entries
            .checked_add(1)
            .expect("one entry beyond the f64 address-space boundary fits usize");
        assert_eq!(
            checked_dense_matrix_allocation(first_unaddressable),
            Err(BipopError::DimensionOverflow {
                what: "dense covariance address space",
            })
        );
        assert_eq!(
            checked_dense_matrix_allocation(usize::MAX),
            Err(BipopError::DimensionOverflow {
                what: "dense covariance bytes",
            }),
            "the checked byte-product branch is independently retained"
        );
    }

    /// G0/G4: the infallible direct-CMA compatibility surface remains a panic
    /// projection, but it shares the finite-query guard with fallible BIPOP and
    /// never exposes the overflowing candidate to user code. This retained seed
    /// has a positive first normal, so `MAX + MAX*z` deterministically overflows
    /// at generation one, candidate zero, component zero.
    #[test]
    fn legacy_cma_panics_before_nonfinite_candidate_callback() {
        let calls = std::cell::Cell::new(0usize);
        let mut objective = |point: &[f64]| {
            calls.set(calls.get() + 1);
            assert!(point.iter().all(|value| value.is_finite()));
            0.0
        };
        let params = CmaParams::standard(1, f64::MAX, 5, f64::NEG_INFINITY);
        let panic = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            cmaes(&mut objective, &[f64::MAX], &params, 0xB1_90_00_02)
        }));

        assert!(
            panic.is_err(),
            "legacy CMA must project typed refusal as panic"
        );
        assert_eq!(
            calls.get(),
            1,
            "only the finite initial point may reach the callback"
        );
    }

    /// G0/G3: production stops immediately after its ninth large run (rungs
    /// zero through eight). A locally self-consistent trailing small record is
    /// therefore a forged history and must fail the independent validator.
    #[test]
    fn ledger_refuses_any_record_after_ninth_large_run() {
        let base_lambda = 4usize;
        let total_budget = 258_056usize;
        let mut records = Vec::new();
        let mut schedule = Vec::new();
        let mut cursor = 0usize;
        let mut push_record =
            |lane: BipopLane, lambda: usize, allocated_budget: usize, generations: usize| {
                let index = records.len();
                let ordinal = u64::try_from(index).expect("fixture ordinal fits");
                let evals = generations
                    .checked_mul(lambda)
                    .and_then(|samples| samples.checked_add(1))
                    .expect("fixture evaluation count fits");
                let trace_end = cursor
                    .checked_add(evals)
                    .expect("fixture trace offset fits");
                let seed = ordinal
                    .checked_mul(BIPOP_RESTART_SEED_STRIDE)
                    .expect("fixture seed fits");
                let report = CmaReport {
                    x_best: vec![0.0],
                    f_best: 0.0,
                    evals,
                    generations,
                    converged: false,
                    sigma: 1.0,
                };
                records.push(BipopRestartRecord {
                    schema_version: BIPOP_RESTART_SCHEMA_VERSION,
                    ordinal,
                    lane,
                    lambda,
                    allocated_budget,
                    seed,
                    start: vec![0.0],
                    trace_start: cursor,
                    trace_end,
                    stop_reason: CmaStopReason::Stagnated,
                    report,
                });
                schedule.push(lambda);
                cursor = trace_end;
            };

        // Each first-eight large rung is paired with enough small-lane
        // generations to make the cumulative spends exactly equal. Every
        // allocated cap is the production local envelope under this large
        // retained aggregate budget.
        for large_run in 0..BIPOP_LARGE_RUN_CAP {
            let scale = 1usize << large_run;
            let large_lambda = base_lambda * scale;
            push_record(
                BipopLane::Large,
                large_lambda,
                large_lambda * BIPOP_GENERATIONS_PER_RESTART,
                1,
            );
            push_record(
                BipopLane::Small,
                base_lambda,
                base_lambda * BIPOP_GENERATIONS_PER_RESTART,
                scale,
            );
        }

        // Rung eight consumes exactly the remaining aggregate allocation and
        // is the ninth/terminal large run. The following locally plausible
        // small record is forbidden before any of its fields can self-justify.
        let final_large_lambda = base_lambda * (1usize << BIPOP_LARGE_RUN_CAP);
        push_record(
            BipopLane::Large,
            final_large_lambda,
            final_large_lambda * BIPOP_GENERATIONS_PER_RESTART,
            1,
        );
        push_record(
            BipopLane::Small,
            base_lambda,
            base_lambda * BIPOP_GENERATIONS_PER_RESTART,
            1,
        );
        assert_eq!(records.len(), 18);
        assert_eq!(cursor, 3_086);
        let report = BipopReport {
            best: records[0].report.clone(),
            schedule,
            total_evals: cursor,
            total_budget,
            records,
            best_restart: 0,
        };

        let error = report
            .validate_ledger()
            .expect_err("post-cap history must fail closed");
        assert_eq!(error.restart(), Some(17));
        assert_eq!(error.invariant(), "large-run-cap");
    }

    /// G0/G3: a forged record cannot mint more large-lane spend than the
    /// production `lambda * 250` envelope. Without this independent check, an
    /// oversized first large record could justify arbitrarily many small
    /// records and invalidate the scheduler theorem used by admission.
    #[test]
    fn ledger_refuses_allocated_budget_above_local_envelope() {
        let lambda = 4usize;
        let allocated_budget = lambda * BIPOP_GENERATIONS_PER_RESTART + 1;
        let generations = (allocated_budget - 1) / lambda;
        let run = CmaReport {
            x_best: vec![0.0],
            f_best: 0.0,
            evals: allocated_budget,
            generations,
            converged: false,
            sigma: 1.0,
        };
        let record = BipopRestartRecord {
            schema_version: BIPOP_RESTART_SCHEMA_VERSION,
            ordinal: 0,
            lane: BipopLane::Large,
            lambda,
            allocated_budget,
            seed: 0,
            start: vec![0.0],
            trace_start: 0,
            trace_end: allocated_budget,
            stop_reason: CmaStopReason::BudgetExhausted,
            report: run.clone(),
        };
        let report = BipopReport {
            best: run.clone(),
            schedule: vec![lambda],
            total_evals: allocated_budget,
            total_budget: allocated_budget,
            records: vec![record],
            best_restart: 0,
        };

        let error = report
            .validate_ledger()
            .expect_err("an oversized local cap is not production-authentic");
        assert_eq!(error.restart(), Some(0));
        assert_eq!(error.invariant(), "allocated-budget");
    }

    /// G0/G3: a prefix can be locally self-consistent yet impossible as a
    /// terminal report under its retained aggregate budget. The validator must
    /// not accept truncation merely because every surviving interval closes.
    #[test]
    fn ledger_refuses_nonterminal_prefix_under_retained_budget() {
        let mut objective = |_point: &[f64]| 0.0;
        let mut report = try_bipop_cmaes(&mut objective, &[0.0], 0.5, 20, None, 0)
            .expect("production fixture admits");
        assert!(report.records.len() > 1, "fixture must contain a suffix");
        report.records.pop();
        report.schedule.pop();
        report.total_evals = report
            .records
            .last()
            .expect("the truncated prefix remains nonempty")
            .trace_end;

        let error = report
            .validate_ledger()
            .expect_err("a nonterminal prefix is not a complete run");
        assert_eq!(error.restart(), None);
        assert_eq!(error.invariant(), "nonterminal-prefix");
    }

    /// G0/G3: target convergence is an immediate scheduler terminal. A later
    /// record cannot be authenticated by giving that record another locally
    /// plausible converged report and relying on the final-record check.
    #[test]
    fn ledger_refuses_any_record_after_target_convergence() {
        let lambda = 4usize;
        let total_budget = 1_001usize;
        let run = CmaReport {
            x_best: vec![0.0],
            f_best: 0.0,
            evals: 1,
            generations: 0,
            converged: true,
            sigma: 1.0,
        };
        let make_record = |index: usize, lane: BipopLane| BipopRestartRecord {
            schema_version: BIPOP_RESTART_SCHEMA_VERSION,
            ordinal: u64::try_from(index).expect("fixture ordinal fits"),
            lane,
            lambda,
            allocated_budget: 1_000,
            seed: u64::try_from(index)
                .expect("fixture ordinal fits")
                .checked_mul(BIPOP_RESTART_SEED_STRIDE)
                .expect("fixture seed fits"),
            start: vec![0.0],
            trace_start: index,
            trace_end: index + 1,
            stop_reason: CmaStopReason::TargetReached,
            report: run.clone(),
        };
        let report = BipopReport {
            best: run.clone(),
            schedule: vec![lambda, lambda],
            total_evals: 2,
            total_budget,
            records: vec![
                make_record(0, BipopLane::Large),
                make_record(1, BipopLane::Small),
            ],
            best_restart: 0,
        };

        let error = report
            .validate_ledger()
            .expect_err("production stops immediately at the first target hit");
        assert_eq!(error.restart(), Some(0));
        assert_eq!(error.invariant(), "continued-after-convergence");
    }

    /// G0/G3: a lucky target hit on the initial callback cannot bypass the
    /// dependency authority that preflight required because the admitted hard
    /// budget could have entered a complete generation.
    #[test]
    fn ledger_early_target_cannot_bypass_reachable_jacobi_authority() {
        let dimension = 4_096usize;
        let lambda = 28usize;
        let total_budget = 29usize;
        let run = CmaReport {
            x_best: vec![0.0; dimension],
            f_best: 0.0,
            evals: 1,
            generations: 0,
            converged: true,
            sigma: 1.0,
        };
        let record = BipopRestartRecord {
            schema_version: BIPOP_RESTART_SCHEMA_VERSION,
            ordinal: 0,
            lane: BipopLane::Large,
            lambda,
            allocated_budget: total_budget,
            seed: 0,
            start: vec![0.0; dimension],
            trace_start: 0,
            trace_end: 1,
            stop_reason: CmaStopReason::TargetReached,
            report: run.clone(),
        };
        let report = BipopReport {
            best: run,
            schedule: vec![lambda],
            total_evals: 1,
            total_budget,
            records: vec![record],
            best_restart: 0,
        };

        let error = report
            .validate_ledger()
            .expect_err("early convergence cannot weaken callback-free admission");
        assert_eq!(error.restart(), Some(0));
        assert_eq!(error.invariant(), "eigensolver-admission");
    }

    /// G0/G3: retained decision points must preserve the finite-query guards
    /// enforced by production. Objective outputs remain intentionally outside
    /// this check.
    #[test]
    fn ledger_refuses_nonfinite_start_and_best_points() {
        let make_report = || {
            let mut objective = |_point: &[f64]| 0.0;
            try_bipop_cmaes(&mut objective, &[0.0], 0.5, 1, None, 0)
                .expect("one-callback fixture admits")
        };

        let mut bad_start = make_report();
        bad_start.records[0].start[0] = f64::NAN;
        let error = bad_start
            .validate_ledger()
            .expect_err("non-finite retained start must refuse");
        assert_eq!(error.restart(), Some(0));
        assert_eq!(error.invariant(), "finite-start");

        let mut bad_best = make_report();
        bad_best.records[0].report.x_best[0] = f64::INFINITY;
        let error = bad_best
            .validate_ledger()
            .expect_err("non-finite retained best point must refuse");
        assert_eq!(error.restart(), Some(0));
        assert_eq!(error.invariant(), "finite-best-point");

        let mut bad_sigma = make_report();
        bad_sigma.records[0].report.sigma = f64::INFINITY;
        let error = bad_sigma
            .validate_ledger()
            .expect_err("non-finite derived step size must refuse");
        assert_eq!(error.restart(), Some(0));
        assert_eq!(error.invariant(), "finite-nonnegative-sigma");

        let mut negative_sigma = make_report();
        negative_sigma.records[0].report.sigma = -0.0;
        let error = negative_sigma
            .validate_ledger()
            .expect_err("a forged negative-sign step size must refuse");
        assert_eq!(error.restart(), Some(0));
        assert_eq!(error.invariant(), "finite-nonnegative-sigma");
    }

    /// G0: an overflowing hypothetical next generation is not evidence that a
    /// run exhausted its local budget. The validator must refuse the arithmetic
    /// boundary instead of treating `checked_add(None)` as "no generation fits."
    #[test]
    fn ledger_refuses_wrapped_next_generation_accounting() {
        let lambda = 6usize;
        let evals = usize::MAX - 2;
        let generations = (evals - 1) / lambda;
        assert_eq!(generations * lambda + 1, evals);

        let run = CmaReport {
            x_best: vec![0.0, 0.0],
            f_best: 0.0,
            evals,
            generations,
            converged: false,
            sigma: 1.0,
        };
        let record = BipopRestartRecord {
            schema_version: BIPOP_RESTART_SCHEMA_VERSION,
            ordinal: 0,
            lane: BipopLane::Large,
            lambda,
            allocated_budget: usize::MAX,
            seed: 7,
            start: vec![0.0, 0.0],
            trace_start: 0,
            trace_end: evals,
            stop_reason: CmaStopReason::BudgetExhausted,
            report: run.clone(),
        };
        let report = BipopReport {
            best: run,
            schedule: vec![lambda],
            total_evals: evals,
            total_budget: usize::MAX,
            records: vec![record],
            best_restart: 0,
        };

        let error = report
            .validate_ledger()
            .expect_err("overflowing next-generation accounting must fail closed");
        assert_eq!(error.restart(), Some(0));
        assert_eq!(error.invariant(), "evaluation-overflow");
    }
}
