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

/// Schema version for [`BipopAdmission`].
pub const BIPOP_ADMISSION_SCHEMA_VERSION: u32 = 1;

/// Schema version for [`BipopRestartRecord`].
pub const BIPOP_RESTART_SCHEMA_VERSION: u32 = 1;

/// Sealed result of callback-free BIPOP input and arithmetic admission.
///
/// The maxima are deliberately conservative. In particular, one objective
/// evaluation is the minimum spend of a launched restart, so
/// `total_budget - 1` bounds every reachable restart ordinal without assuming
/// convergence or stagnation behavior. Holding this receipt proves only that
/// the checked formulas below were representable; it is not an authenticated
/// identity for the start, target, sigma, seed, or callback.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BipopAdmission {
    schema_version: u32,
    dimension: usize,
    total_budget: usize,
    base_lambda: usize,
    max_large_lambda: usize,
    max_local_budget: usize,
    max_restart_ordinal: u64,
    max_matrix_entries: usize,
    max_population_entries: usize,
}

impl BipopAdmission {
    /// Admission schema used for the receipt.
    #[must_use]
    pub fn schema_version(&self) -> u32 {
        self.schema_version
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

    /// Conservative largest restart ordinal (`total_budget - 1`).
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
    /// The aggregate budget cannot be represented as a u64 restart ordinal.
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
                "BIPOP budget {total_budget} cannot be represented by the u64 restart ordinal"
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

impl std::error::Error for BipopError {}

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

    let max_restart_ordinal_usize = total_budget - 1;
    let max_restart_ordinal = u64::try_from(max_restart_ordinal_usize)
        .map_err(|_| BipopError::RestartOrdinalOverflow { total_budget })?;
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

    let budget_ladder_cap = u32::try_from(max_restart_ordinal_usize)
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
    max_matrix_entries
        .checked_mul(core::mem::size_of::<f64>())
        .ok_or(BipopError::DimensionOverflow {
            what: "dense covariance bytes",
        })?;
    max_population_entries
        .checked_mul(core::mem::size_of::<f64>())
        .ok_or(BipopError::DimensionOverflow {
            what: "population coordinate bytes",
        })?;

    Ok(BipopAdmission {
        schema_version: BIPOP_ADMISSION_SCHEMA_VERSION,
        dimension,
        total_budget,
        base_lambda,
        max_large_lambda,
        max_local_budget,
        max_restart_ordinal,
        max_matrix_entries,
        max_population_entries,
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
    cmaes_with_stop_target(f, x0, params, seed, Some(params.f_target))
}

#[allow(clippy::too_many_lines)] // the algorithm is one coherent loop
fn cmaes_with_stop_target<F: FnMut(&[f64]) -> f64>(
    f: &mut F,
    x0: &[f64],
    params: &CmaParams,
    seed: u64,
    target: Option<f64>,
) -> (CmaReport, CmaStopReason) {
    let n = x0.len();
    assert!(n >= 1, "dimension must be positive");
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
        return (
            CmaReport {
                x_best,
                f_best,
                evals,
                generations,
                converged: true,
                sigma,
            },
            CmaStopReason::TargetReached,
        );
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
            return (
                CmaReport {
                    x_best,
                    f_best,
                    evals,
                    generations,
                    converged: true,
                    sigma,
                },
                CmaStopReason::TargetReached,
            );
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
    (
        CmaReport {
            x_best,
            f_best,
            evals,
            generations,
            converged: false,
            sigma,
        },
        stop_reason,
    )
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
    records: Vec<BipopRestartRecord>,
    best_restart: usize,
}

impl BipopReport {
    /// Ordered immutable restart ledger.
    #[must_use]
    pub fn records(&self) -> &[BipopRestartRecord] {
        &self.records
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
        let base_lambda = first.lambda;
        let base_seed = first.seed;
        let point_dim = first.start.len();
        if point_dim == 0 {
            return Err(BipopLedgerError::at(0, "positive-point-dimension"));
        }
        let expected_base_lambda = 4 + (3.0 * fs_math::det::ln(point_dim as f64)).floor() as usize;
        if base_lambda != expected_base_lambda {
            return Err(BipopLedgerError::at(0, "base-population"));
        }
        let mut cursor = 0usize;
        let mut large_budget_used = 0usize;
        let mut small_budget_used = 0usize;
        let mut large_runs = 0u32;

        for (index, record) in self.records.iter().enumerate() {
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
/// refusal. A later generated-start refusal can follow completed earlier
/// restarts, but the refused restart itself is never invoked or published.
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
    let mut records: Vec<BipopRestartRecord> = Vec::new();
    let mut total_evals = 0usize;
    let mut best_restart: Option<usize> = None;
    let mut large_runs = 0u32;
    let mut restart = 0u64;
    let mut small_budget_used = 0usize;
    let mut large_budget_used = 0usize;
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
            x0.iter()
                .map(|&v| sigma0.mul_add(restart_stream.next_normal(), v))
                .collect()
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
        let (rep, stop_reason) = cmaes_with_stop_target(f, &start, &params, derived_seed, target);
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
