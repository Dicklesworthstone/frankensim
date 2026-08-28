//! Generation-oriented CMA family implementations.
//!
//! This module keeps four materially different covariance models behind one
//! deterministic ask/tell contract: active full CMA-ES, separable CMA-ES,
//! LM-CMA-ES, and LM-MA-ES.  Sampling uses one keyed Philox stream and ranking
//! always breaks equal objective values by candidate index.

#![allow(clippy::needless_range_loop)] // explicit indices are numeric-kernel coordinates

use core::fmt;

use fs_la::eigen::{JacobiEighAdmissionError, admit_jacobi_eigh, jacobi_eigh};
use fs_rand::{Stream, StreamKey};

/// Stable fs-rand kernel identity for generation-oriented CMA sampling.
pub const CMA_FAMILY_STREAM_KERNEL: u32 = 0xD1F1;
const EIGENVALUE_FLOOR: f64 = 1.0e-14;

/// Covariance representation used by a [`CmaOptimizer`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CmaFamily {
    /// Dense covariance with Hansen active (negative-weight) adaptation.
    Full,
    /// Diagonal covariance using sep-CMA rates and active recombination.
    Separable,
    /// Loshchilov's limited-memory Cholesky-factor CMA-ES.
    LmCma,
    /// Loshchilov, Glasmachers, and Beyer's limited-memory matrix adaptation.
    LmMa,
}

/// Coarse asymptotic order for one part of a CMA implementation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CmaComplexityOrder {
    /// O(n).
    Linear,
    /// O(mn), where `m` is the configured limited-memory capacity.
    MemoryLinear,
    /// O(m²n), where `m` is the configured limited-memory capacity.
    MemoryQuadratic,
    /// O(n²).
    Quadratic,
    /// O(n³).
    Cubic,
}

/// Honest allocation and work metadata for an admitted optimizer.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CmaComplexity {
    /// Work to transform one sampled standard-normal vector.
    pub sampling_per_candidate: CmaComplexityOrder,
    /// Dominant work in one completed generation update.
    pub update_per_generation: CmaComplexityOrder,
    /// Maximum persistent floating-point scalar count, including best-so-far.
    pub persistent_scalars: usize,
    /// Floating-point scalars held while an ask batch is outstanding.
    pub pending_generation_scalars: usize,
    /// Conservative peak floating-point workspace during transactional `tell`.
    pub update_workspace_scalars: usize,
    /// Dense square-matrix entries retained by the strategy.
    pub dense_matrix_entries: usize,
    /// Configured limited-memory direction capacity, or zero.
    pub memory_capacity: usize,
}

/// Configuration for the generation-oriented CMA family API.
#[derive(Debug, Clone)]
pub struct CmaConfig {
    /// Covariance representation and update family.
    pub family: CmaFamily,
    /// Initial distribution mean.
    pub mean: Vec<f64>,
    /// Initial global step size.
    pub sigma: f64,
    /// Hard objective-evaluation budget. Only complete populations are spent.
    pub max_evaluations: usize,
    /// Deterministic Philox study seed.
    pub seed: u64,
    /// Optional population size; the default is `4 + floor(3 ln(n))`.
    pub population_size: Option<usize>,
    /// Optional bounded-memory capacity for LM-CMA and LM-MA.
    pub memory: Option<usize>,
}

impl CmaConfig {
    /// Construct the reference-default configuration for one family.
    #[must_use]
    pub fn standard(
        family: CmaFamily,
        mean: Vec<f64>,
        sigma: f64,
        max_evaluations: usize,
        seed: u64,
    ) -> Self {
        Self {
            family,
            mean,
            sigma,
            max_evaluations,
            seed,
            population_size: None,
            memory: None,
        }
    }
}

/// Allocation-free result of validating a [`CmaConfig`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CmaAdmission {
    /// Exact fs-rand draw-semantics version used by admitted sampling.
    pub stream_semantics_version: u32,
    /// Stable fs-rand kernel identity used by admitted sampling.
    pub stream_kernel: u32,
    /// Decision-space dimension.
    pub dimension: usize,
    /// Population size used by every admitted generation.
    pub population_size: usize,
    /// Number of selected parents.
    pub parent_count: usize,
    /// Complete generations admitted by the budget.
    pub max_generations: usize,
    /// Exact maximum number of evaluations the optimizer can spend.
    pub admitted_evaluations: usize,
    /// Maximum Philox blocks consumed by standard-normal sampling.
    pub normal_stream_blocks: u64,
    /// Work and storage metadata.
    pub complexity: CmaComplexity,
}

/// Typed validation, ask/tell-contract, and numerical failures.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CmaFamilyError {
    /// An optimizer needs at least one decision variable.
    EmptyMean,
    /// One initial-mean coordinate is not finite.
    NonFiniteMean {
        /// Zero-based bad coordinate.
        coordinate: usize,
        /// Exact IEEE-754 payload.
        bits: u64,
    },
    /// The initial step size must be finite and strictly positive.
    InvalidSigma {
        /// Exact IEEE-754 payload.
        bits: u64,
    },
    /// Population sizes below four do not define this implementation.
    InvalidPopulation {
        /// Rejected population size.
        population_size: usize,
    },
    /// No complete population fits the objective budget.
    BudgetTooSmall {
        /// Requested evaluation budget.
        max_evaluations: usize,
        /// Complete-population cost.
        population_size: usize,
    },
    /// A limited-memory capacity must be nonzero.
    InvalidMemory {
        /// Rejected memory capacity.
        memory: usize,
    },
    /// Memory capacity was supplied for a dense or diagonal strategy.
    MemoryNotApplicable {
        /// Family that does not consume limited-memory configuration.
        family: CmaFamily,
    },
    /// Checked shape or storage arithmetic overflowed.
    ShapeOverflow {
        /// Checked product or aggregate that overflowed.
        context: &'static str,
    },
    /// The exact Philox draw envelope exceeds the stream counter domain.
    RandomCounterOverflow,
    /// Dense Jacobi workspace admission refused the dimension.
    DenseEigensolver(JacobiEighAdmissionError),
    /// `ask` was called before the outstanding batch was told.
    AskAlreadyPending {
        /// Outstanding generation.
        generation: u64,
    },
    /// No complete population remains in the admitted budget.
    BudgetExhausted {
        /// Evaluations left in the admitted complete-generation envelope.
        remaining: usize,
        /// Population cost of another generation.
        required: usize,
    },
    /// `tell` was called without an outstanding ask batch.
    NoPendingAsk,
    /// The batch belongs to a different generation.
    GenerationMismatch {
        /// Outstanding generation.
        expected: u64,
        /// Supplied batch generation.
        actual: u64,
    },
    /// The batch does not match the optimizer's outstanding batch.
    BatchMismatch,
    /// Objective count differs from the population size.
    ObjectiveCount {
        /// Population size.
        expected: usize,
        /// Supplied objective count.
        actual: usize,
    },
    /// Objective values must be finite.
    NonFiniteObjective {
        /// Candidate index.
        candidate: usize,
        /// Exact IEEE-754 payload.
        bits: u64,
    },
    /// An internal numerical result left the finite floating-point domain.
    NumericalFailure {
        /// Numerical operation being checked.
        stage: &'static str,
        /// Coordinate or flattened entry index.
        coordinate: usize,
        /// Exact IEEE-754 payload.
        bits: u64,
    },
}

impl fmt::Display for CmaFamilyError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match *self {
            Self::EmptyMean => formatter.write_str("CMA mean must not be empty"),
            Self::NonFiniteMean { coordinate, bits } => write!(
                formatter,
                "CMA mean coordinate {coordinate} is non-finite (bits 0x{bits:016x})"
            ),
            Self::InvalidSigma { bits } => write!(
                formatter,
                "CMA sigma must be finite and positive (bits 0x{bits:016x})"
            ),
            Self::InvalidPopulation { population_size } => {
                write!(formatter, "CMA population {population_size} is below four")
            }
            Self::BudgetTooSmall {
                max_evaluations,
                population_size,
            } => write!(
                formatter,
                "CMA budget {max_evaluations} cannot admit population {population_size}"
            ),
            Self::InvalidMemory { memory } => {
                write!(formatter, "CMA limited memory {memory} must be nonzero")
            }
            Self::MemoryNotApplicable { family } => {
                write!(formatter, "memory capacity does not apply to {family:?}")
            }
            Self::ShapeOverflow { context } => {
                write!(formatter, "CMA checked shape overflowed in {context}")
            }
            Self::RandomCounterOverflow => {
                formatter.write_str("CMA normal draws exceed the Philox counter domain")
            }
            Self::DenseEigensolver(error) => write!(formatter, "{error}"),
            Self::AskAlreadyPending { generation } => {
                write!(formatter, "CMA generation {generation} is still pending")
            }
            Self::BudgetExhausted {
                remaining,
                required,
            } => write!(
                formatter,
                "CMA budget has {remaining} evaluations left; {required} are required"
            ),
            Self::NoPendingAsk => formatter.write_str("CMA tell has no pending ask batch"),
            Self::GenerationMismatch { expected, actual } => write!(
                formatter,
                "CMA tell expected generation {expected}, received {actual}"
            ),
            Self::BatchMismatch => {
                formatter.write_str("CMA tell batch does not match the pending ask")
            }
            Self::ObjectiveCount { expected, actual } => write!(
                formatter,
                "CMA tell expected {expected} objectives, received {actual}"
            ),
            Self::NonFiniteObjective { candidate, bits } => write!(
                formatter,
                "CMA objective {candidate} is non-finite (bits 0x{bits:016x})"
            ),
            Self::NumericalFailure {
                stage,
                coordinate,
                bits,
            } => write!(
                formatter,
                "CMA numerical failure in {stage} at coordinate {coordinate} (bits 0x{bits:016x})"
            ),
        }
    }
}

impl std::error::Error for CmaFamilyError {}

impl From<JacobiEighAdmissionError> for CmaFamilyError {
    fn from(error: JacobiEighAdmissionError) -> Self {
        Self::DenseEigensolver(error)
    }
}

/// One immutable candidate population returned by [`CmaOptimizer::ask`].
#[derive(Debug, Clone)]
pub struct CmaAsk {
    generation: u64,
    signature: u64,
    candidates: Vec<Vec<f64>>,
    isotropic_steps: Vec<Vec<f64>>,
    distribution_steps: Vec<Vec<f64>>,
}

impl CmaAsk {
    /// Zero-based generation identifier.
    #[must_use]
    pub const fn generation(&self) -> u64 {
        self.generation
    }

    /// Candidate points in deterministic population order.
    #[must_use]
    pub fn candidates(&self) -> &[Vec<f64>] {
        &self.candidates
    }

    /// Number of candidates in this complete generation.
    #[must_use]
    pub fn len(&self) -> usize {
        self.candidates.len()
    }

    /// Whether this batch contains no candidates.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.candidates.is_empty()
    }
}

/// Earliest best-so-far observation under deterministic total ordering.
#[derive(Debug, Clone, PartialEq)]
pub struct CmaBest {
    /// Decision point.
    pub point: Vec<f64>,
    /// Objective value.
    pub objective: f64,
    /// Zero-based generation containing the observation.
    pub generation: u64,
    /// Candidate index within that generation.
    pub candidate: usize,
}

/// Shape diagnostics whose payload matches the strategy's real storage.
#[derive(Debug, Clone, PartialEq)]
pub enum CmaShapeSnapshot {
    /// Dense covariance diagnostics. No dense matrix is copied into snapshots.
    Full {
        /// Covariance diagonal.
        diagonal: Vec<f64>,
        /// Smallest eigenvalue after the SPD floor.
        min_eigenvalue: f64,
        /// Largest eigenvalue after the SPD floor.
        max_eigenvalue: f64,
        /// Count of strictly negative active recombination weights.
        negative_weight_count: usize,
    },
    /// Separable covariance diagnostics.
    Diagonal {
        /// Exact diagonal variances.
        variances: Vec<f64>,
        /// Count of strictly negative active recombination weights.
        negative_weight_count: usize,
    },
    /// Bounded-memory diagnostics without copying the stored directions.
    LimitedMemory {
        /// Current number of retained direction vectors.
        vectors: usize,
        /// Configured vector capacity.
        capacity: usize,
        /// Euclidean norm of each retained update direction.
        direction_norms: Vec<f64>,
    },
}

/// Browser-sized diagnostics after a completed generation.
#[derive(Debug, Clone, PartialEq)]
pub struct CmaSnapshot {
    /// Strategy family.
    pub family: CmaFamily,
    /// Number of completed generations.
    pub generation: u64,
    /// Exact objective evaluations consumed.
    pub evaluations: usize,
    /// Current search-distribution mean.
    pub mean: Vec<f64>,
    /// Current global step size.
    pub sigma: f64,
    /// Best finite observation, if any generation has completed.
    pub best: Option<CmaBest>,
    /// Representation-specific shape diagnostics.
    pub shape: CmaShapeSnapshot,
    /// Stable complexity and storage metadata.
    pub complexity: CmaComplexity,
}

/// Validate shape, budget, random-counter, and dense-work requirements.
#[allow(clippy::too_many_lines)] // one receipt keeps all coupled admission bounds together
pub fn admit_cma(config: &CmaConfig) -> Result<CmaAdmission, CmaFamilyError> {
    let n = config.mean.len();
    if n == 0 {
        return Err(CmaFamilyError::EmptyMean);
    }
    for (coordinate, &value) in config.mean.iter().enumerate() {
        if !value.is_finite() {
            return Err(CmaFamilyError::NonFiniteMean {
                coordinate,
                bits: value.to_bits(),
            });
        }
    }
    if !config.sigma.is_finite() || config.sigma <= 0.0 {
        return Err(CmaFamilyError::InvalidSigma {
            bits: config.sigma.to_bits(),
        });
    }
    let default_population = 4usize
        .checked_add((3.0 * fs_math::det::ln(n as f64)).floor() as usize)
        .ok_or(CmaFamilyError::ShapeOverflow {
            context: "default population",
        })?;
    let lambda = config.population_size.unwrap_or(default_population);
    if lambda < 4 {
        return Err(CmaFamilyError::InvalidPopulation {
            population_size: lambda,
        });
    }
    let max_generations = config.max_evaluations / lambda;
    if max_generations == 0 {
        return Err(CmaFamilyError::BudgetTooSmall {
            max_evaluations: config.max_evaluations,
            population_size: lambda,
        });
    }
    let memory = match config.family {
        CmaFamily::LmCma | CmaFamily::LmMa => config.memory.unwrap_or(default_population),
        family => {
            if config.memory.is_some() {
                return Err(CmaFamilyError::MemoryNotApplicable { family });
            }
            0
        }
    };
    if matches!(config.family, CmaFamily::LmCma | CmaFamily::LmMa) && memory == 0 {
        return Err(CmaFamilyError::InvalidMemory { memory });
    }
    let mu = lambda / 2;
    let square = if config.family == CmaFamily::Full {
        admit_jacobi_eigh(n)?.matrix_entries()
    } else {
        0
    };
    let pending = lambda
        .checked_mul(n)
        .and_then(|value| value.checked_mul(3))
        .ok_or(CmaFamilyError::ShapeOverflow {
            context: "pending generation",
        })?;
    let common = n
        .checked_mul(2)
        .and_then(|value| value.checked_add(mu))
        .and_then(|value| value.checked_add(4))
        .ok_or(CmaFamilyError::ShapeOverflow {
            context: "common state",
        })?;
    let (sampling, update, dense_entries, strategy_scalars, update_workspace) = match config.family
    {
        CmaFamily::Full => (
            CmaComplexityOrder::Quadratic,
            CmaComplexityOrder::Cubic,
            square.checked_mul(2).ok_or(CmaFamilyError::ShapeOverflow {
                context: "full dense state",
            })?,
            square
                .checked_mul(2)
                .and_then(|value| value.checked_add(n.checked_mul(3)?))
                .and_then(|value| value.checked_add(lambda))
                .and_then(|value| value.checked_add(5))
                .ok_or(CmaFamilyError::ShapeOverflow {
                    context: "full strategy state",
                })?,
            square
                .checked_mul(7)
                .and_then(|value| value.checked_add(n.checked_mul(9)?))
                .and_then(|value| value.checked_add(lambda.checked_mul(2)?))
                .and_then(|value| value.checked_add(6))
                .ok_or(CmaFamilyError::ShapeOverflow {
                    context: "full update workspace",
                })?,
        ),
        CmaFamily::Separable => (
            CmaComplexityOrder::Linear,
            CmaComplexityOrder::Linear,
            0,
            n.checked_mul(3)
                .and_then(|value| value.checked_add(lambda))
                .and_then(|value| value.checked_add(5))
                .ok_or(CmaFamilyError::ShapeOverflow {
                    context: "separable strategy state",
                })?,
            n.checked_mul(6)
                .and_then(|value| value.checked_add(lambda.checked_mul(2)?))
                .and_then(|value| value.checked_add(6))
                .ok_or(CmaFamilyError::ShapeOverflow {
                    context: "separable update workspace",
                })?,
        ),
        CmaFamily::LmCma => (
            CmaComplexityOrder::MemoryLinear,
            CmaComplexityOrder::MemoryQuadratic,
            0,
            memory
                .checked_mul(n.checked_mul(2).and_then(|v| v.checked_add(2)).ok_or(
                    CmaFamilyError::ShapeOverflow {
                        context: "LM-CMA record",
                    },
                )?)
                .and_then(|value| value.checked_add(n))
                .and_then(|value| value.checked_add(lambda))
                .and_then(|value| value.checked_add(mu))
                .and_then(|value| value.checked_add(6))
                .ok_or(CmaFamilyError::ShapeOverflow {
                    context: "LM-CMA strategy state",
                })?,
            memory
                .checked_mul(n.checked_mul(2).and_then(|v| v.checked_add(2)).ok_or(
                    CmaFamilyError::ShapeOverflow {
                        context: "LM-CMA cloned record",
                    },
                )?)
                .and_then(|value| value.checked_add(n.checked_mul(4)?))
                .and_then(|value| value.checked_add(lambda.checked_mul(4)?))
                .and_then(|value| value.checked_add(mu))
                .and_then(|value| value.checked_add(9))
                .ok_or(CmaFamilyError::ShapeOverflow {
                    context: "LM-CMA update workspace",
                })?,
        ),
        CmaFamily::LmMa => (
            CmaComplexityOrder::MemoryLinear,
            CmaComplexityOrder::MemoryLinear,
            0,
            memory
                .checked_mul(n.checked_add(2).ok_or(CmaFamilyError::ShapeOverflow {
                    context: "LM-MA path",
                })?)
                .and_then(|value| value.checked_add(n))
                .and_then(|value| value.checked_add(1))
                .ok_or(CmaFamilyError::ShapeOverflow {
                    context: "LM-MA strategy state",
                })?,
            memory
                .checked_mul(n.checked_add(2).ok_or(CmaFamilyError::ShapeOverflow {
                    context: "LM-MA cloned path",
                })?)
                .and_then(|value| value.checked_add(n.checked_mul(3)?))
                .and_then(|value| value.checked_add(2))
                .ok_or(CmaFamilyError::ShapeOverflow {
                    context: "LM-MA update workspace",
                })?,
        ),
    };
    let persistent_scalars =
        common
            .checked_add(strategy_scalars)
            .ok_or(CmaFamilyError::ShapeOverflow {
                context: "persistent state",
            })?;
    let blocks = (max_generations as u128)
        .checked_mul(lambda as u128)
        .and_then(|value| value.checked_mul(n as u128))
        .and_then(|value| value.checked_mul(2))
        .ok_or(CmaFamilyError::RandomCounterOverflow)?;
    let normal_stream_blocks =
        u64::try_from(blocks).map_err(|_| CmaFamilyError::RandomCounterOverflow)?;
    Ok(CmaAdmission {
        stream_semantics_version: fs_rand::STREAM_SEMANTICS_VERSION,
        stream_kernel: CMA_FAMILY_STREAM_KERNEL,
        dimension: n,
        population_size: lambda,
        parent_count: mu,
        max_generations,
        admitted_evaluations: max_generations * lambda,
        normal_stream_blocks,
        complexity: CmaComplexity {
            sampling_per_candidate: sampling,
            update_per_generation: update,
            persistent_scalars,
            pending_generation_scalars: pending,
            update_workspace_scalars: update_workspace,
            dense_matrix_entries: dense_entries,
            memory_capacity: memory,
        },
    })
}

#[derive(Debug, Clone)]
struct SharedParameters {
    n: usize,
    lambda: usize,
    weights: Vec<f64>,
    mu_eff: f64,
    chi_n: f64,
}

impl SharedParameters {
    fn new(admission: CmaAdmission) -> Self {
        let lambda = admission.population_size;
        let mu = admission.parent_count;
        let midpoint = f64::midpoint(lambda as f64, 1.0);
        let mut weights: Vec<f64> = (1..=mu)
            .map(|rank| fs_math::det::ln(midpoint) - fs_math::det::ln(rank as f64))
            .collect();
        normalize(&mut weights);
        let mu_eff = 1.0 / weights.iter().map(|weight| weight * weight).sum::<f64>();
        let nf = admission.dimension as f64;
        let chi_n = fs_math::det::sqrt(nf) * (1.0 - 1.0 / (4.0 * nf) + 1.0 / (21.0 * nf * nf));
        Self {
            n: admission.dimension,
            lambda,
            weights,
            mu_eff,
            chi_n,
        }
    }
}

#[derive(Debug, Clone)]
struct FullState {
    covariance: Vec<f64>,
    basis: Vec<f64>,
    axis_scales: Vec<f64>,
    p_c: Vec<f64>,
    p_sigma: Vec<f64>,
    all_weights: Vec<f64>,
    c_c: f64,
    c_sigma: f64,
    c1: f64,
    c_mu: f64,
    damping: f64,
}

impl FullState {
    fn new(shared: &SharedParameters) -> Self {
        let n = shared.n;
        let nf = n as f64;
        let c_c = (4.0 + shared.mu_eff / nf) / (nf + 4.0 + 2.0 * shared.mu_eff / nf);
        let c_sigma = (shared.mu_eff + 2.0) / (nf + shared.mu_eff + 5.0);
        let c1 = (shared.lambda as f64 / 6.0).min(1.0) * 2.0 / ((nf + 1.3).powi(2) + shared.mu_eff);
        let c_mu = (1.0 - c1).min(
            2.0 * (0.25 + shared.mu_eff + 1.0 / shared.mu_eff - 2.0)
                / ((nf + 2.0).powi(2) + shared.mu_eff),
        );
        let damping = 1.0
            + 2.0 * (fs_math::det::sqrt((shared.mu_eff - 1.0) / (nf + 1.0)) - 1.0).max(0.0)
            + c_sigma;
        let all_weights = active_weights(shared.lambda, shared.mu_eff, nf, c1, c_mu);
        let mut covariance = vec![0.0; n * n];
        let mut basis = vec![0.0; n * n];
        for index in 0..n {
            covariance[index * n + index] = 1.0;
            basis[index * n + index] = 1.0;
        }
        Self {
            covariance,
            basis,
            axis_scales: vec![1.0; n],
            p_c: vec![0.0; n],
            p_sigma: vec![0.0; n],
            all_weights,
            c_c,
            c_sigma,
            c1,
            c_mu,
            damping,
        }
    }

    fn sample(&self, z: &[f64], output: &mut [f64]) {
        let n = z.len();
        for row in 0..n {
            let mut value = 0.0;
            for column in 0..n {
                value = (self.basis[row * n + column] * self.axis_scales[column])
                    .mul_add(z[column], value);
            }
            output[row] = value;
        }
    }

    fn update(
        &mut self,
        shared: &SharedParameters,
        generation: u64,
        order: &[usize],
        zs: &[Vec<f64>],
        ys: &[Vec<f64>],
        sigma: &mut f64,
    ) -> Result<Vec<f64>, CmaFamilyError> {
        let n = shared.n;
        let y_w = weighted_sum(&shared.weights, order, ys, n);
        let z_w = weighted_sum(&shared.weights, order, zs, n);
        // C^{-1/2} y_w = B z_w for C = B D² Bᵀ.  The eigenbasis rotation is
        // essential once B is no longer identity; using raw z_w here silently
        // breaks full CMA's coordinate-rotation equivariance.
        let mut whitened_step = vec![0.0; n];
        for row in 0..n {
            for column in 0..n {
                whitened_step[row] =
                    self.basis[row * n + column].mul_add(z_w[column], whitened_step[row]);
            }
        }
        let path_scale = fs_math::det::sqrt(self.c_sigma * (2.0 - self.c_sigma) * shared.mu_eff);
        for coordinate in 0..n {
            self.p_sigma[coordinate] = (1.0 - self.c_sigma).mul_add(
                self.p_sigma[coordinate],
                path_scale * whitened_step[coordinate],
            );
        }
        let path_norm = norm(&self.p_sigma);
        *sigma *=
            fs_math::det::exp((self.c_sigma / self.damping) * (path_norm / shared.chi_n - 1.0));
        ensure_positive(*sigma, "full step size", 0)?;
        let decay_power = i32::try_from((generation + 1).min(100_000)).unwrap_or(100_000) * 2;
        let h_sigma = path_norm
            / fs_math::det::sqrt(1.0 - fs_math::det::powi(1.0 - self.c_sigma, decay_power))
            < (1.4 + 2.0 / (n as f64 + 1.0)) * shared.chi_n;
        let covariance_path_scale = fs_math::det::sqrt(self.c_c * (2.0 - self.c_c) * shared.mu_eff);
        for coordinate in 0..n {
            self.p_c[coordinate] = (1.0 - self.c_c).mul_add(
                self.p_c[coordinate],
                if h_sigma {
                    covariance_path_scale * y_w[coordinate]
                } else {
                    0.0
                },
            );
        }
        ensure_vector_finite(&self.p_c, "full covariance path")?;
        let delta_h = if h_sigma {
            0.0
        } else {
            self.c_c * (2.0 - self.c_c)
        };
        let weight_sum: f64 = self.all_weights.iter().sum();
        let old_factor = 1.0 + self.c1 * delta_h - self.c1 - self.c_mu * weight_sum;
        let covariance_weights = active_covariance_weights(&self.all_weights, order, zs, n);
        let old_covariance = self.covariance.clone();
        for row in 0..n {
            for column in row..n {
                let mut rank_mu = 0.0;
                for (&effective, &candidate) in covariance_weights.iter().zip(order) {
                    rank_mu =
                        (effective * ys[candidate][row]).mul_add(ys[candidate][column], rank_mu);
                }
                let index = row * n + column;
                let next = old_factor.mul_add(
                    old_covariance[index],
                    self.c1 * self.p_c[row] * self.p_c[column] + self.c_mu * rank_mu,
                );
                ensure_finite(next, "full covariance", index)?;
                self.covariance[index] = next;
                self.covariance[column * n + row] = next;
            }
        }
        self.refresh_factor()?;
        Ok(y_w)
    }

    fn refresh_factor(&mut self) -> Result<(), CmaFamilyError> {
        let n = self.axis_scales.len();
        let (eigenvalues, eigenvectors) = jacobi_eigh(&self.covariance, n);
        ensure_vector_finite(&eigenvalues, "full eigendecomposition")?;
        let maximum = eigenvalues
            .last()
            .copied()
            .unwrap_or(1.0)
            .max(f64::MIN_POSITIVE);
        for (index, &value) in eigenvalues.iter().enumerate() {
            let floored = value.max(EIGENVALUE_FLOOR * maximum);
            ensure_positive(floored, "full eigenvalue", index)?;
            self.axis_scales[index] = fs_math::det::sqrt(floored);
        }
        self.basis = eigenvectors;
        for row in 0..n {
            for column in row..n {
                let mut value = 0.0;
                for axis in 0..n {
                    value = (self.basis[row * n + axis]
                        * self.axis_scales[axis]
                        * self.axis_scales[axis])
                        .mul_add(self.basis[column * n + axis], value);
                }
                self.covariance[row * n + column] = value;
                self.covariance[column * n + row] = value;
            }
        }
        Ok(())
    }
}

#[derive(Debug, Clone)]
struct SeparableState {
    variances: Vec<f64>,
    p_c: Vec<f64>,
    p_sigma: Vec<f64>,
    all_weights: Vec<f64>,
    c_c: f64,
    c_sigma: f64,
    c1: f64,
    c_mu: f64,
    damping: f64,
}

impl SeparableState {
    fn new(shared: &SharedParameters) -> Self {
        let nf = shared.n as f64;
        let sqrt_n = fs_math::det::sqrt(nf);
        let c_c =
            (1.0 + 1.0 / nf + shared.mu_eff / nf) / (sqrt_n + 1.0 / nf + 2.0 * shared.mu_eff / nf);
        let c_sigma = (shared.mu_eff + 2.0) / (nf + shared.mu_eff + 5.0);
        let c1 = 1.0 / (nf + 2.0 * sqrt_n + shared.mu_eff / nf);
        let c_mu = (1.0 - c1).min(
            (0.25 + shared.mu_eff + 1.0 / shared.mu_eff - 2.0)
                / (nf + 4.0 * sqrt_n + shared.mu_eff / 2.0),
        );
        let damping = 1.0
            + 2.0 * (fs_math::det::sqrt((shared.mu_eff - 1.0) / (nf + 1.0)) - 1.0).max(0.0)
            + c_sigma;
        let all_weights = active_weights(shared.lambda, shared.mu_eff, nf, c1, c_mu);
        Self {
            variances: vec![1.0; shared.n],
            p_c: vec![0.0; shared.n],
            p_sigma: vec![0.0; shared.n],
            all_weights,
            c_c,
            c_sigma,
            c1,
            c_mu,
            damping,
        }
    }

    fn sample(&self, z: &[f64], output: &mut [f64]) {
        for ((out, &normal), &variance) in output.iter_mut().zip(z).zip(&self.variances) {
            *out = fs_math::det::sqrt(variance) * normal;
        }
    }

    fn update(
        &mut self,
        shared: &SharedParameters,
        generation: u64,
        order: &[usize],
        zs: &[Vec<f64>],
        ys: &[Vec<f64>],
        sigma: &mut f64,
    ) -> Result<Vec<f64>, CmaFamilyError> {
        let y_w = weighted_sum(&shared.weights, order, ys, shared.n);
        let z_w = weighted_sum(&shared.weights, order, zs, shared.n);
        let path_scale = fs_math::det::sqrt(self.c_sigma * (2.0 - self.c_sigma) * shared.mu_eff);
        for coordinate in 0..shared.n {
            self.p_sigma[coordinate] = (1.0 - self.c_sigma)
                .mul_add(self.p_sigma[coordinate], path_scale * z_w[coordinate]);
        }
        let path_norm = norm(&self.p_sigma);
        *sigma *=
            fs_math::det::exp((self.c_sigma / self.damping) * (path_norm / shared.chi_n - 1.0));
        ensure_positive(*sigma, "separable step size", 0)?;
        let decay_power = i32::try_from((generation + 1).min(100_000)).unwrap_or(100_000) * 2;
        let h_sigma = path_norm
            / fs_math::det::sqrt(1.0 - fs_math::det::powi(1.0 - self.c_sigma, decay_power))
            < (1.4 + 2.0 / (shared.n as f64 + 1.0)) * shared.chi_n;
        let covariance_path_scale = fs_math::det::sqrt(self.c_c * (2.0 - self.c_c) * shared.mu_eff);
        let variance_floor = EIGENVALUE_FLOOR
            * self
                .variances
                .iter()
                .copied()
                .fold(f64::MIN_POSITIVE, f64::max);
        let weight_sum: f64 = self.all_weights.iter().sum();
        let covariance_weights = active_covariance_weights(&self.all_weights, order, zs, shared.n);
        let mut rank_mu = vec![0.0; shared.n];
        for (&effective, &candidate) in covariance_weights.iter().zip(order) {
            for coordinate in 0..shared.n {
                rank_mu[coordinate] = (effective * ys[candidate][coordinate])
                    .mul_add(ys[candidate][coordinate], rank_mu[coordinate]);
            }
        }
        let delta_h = if h_sigma {
            0.0
        } else {
            self.c_c * (2.0 - self.c_c)
        };
        let old_factor = 1.0 + self.c1 * delta_h - self.c1 - self.c_mu * weight_sum;
        for coordinate in 0..shared.n {
            self.p_c[coordinate] = (1.0 - self.c_c).mul_add(
                self.p_c[coordinate],
                if h_sigma {
                    covariance_path_scale * y_w[coordinate]
                } else {
                    0.0
                },
            );
            let next = old_factor.mul_add(
                self.variances[coordinate],
                self.c1 * self.p_c[coordinate] * self.p_c[coordinate]
                    + self.c_mu * rank_mu[coordinate],
            );
            ensure_finite(next, "separable variance", coordinate)?;
            self.variances[coordinate] = next.max(variance_floor);
            ensure_positive(self.variances[coordinate], "separable variance", coordinate)?;
        }
        Ok(y_w)
    }
}

#[derive(Debug, Clone)]
struct LmCmaRecord {
    generation: u64,
    p: Vec<f64>,
    v: Vec<f64>,
    b: f64,
    d: f64,
}

#[derive(Debug, Clone)]
struct LmCmaState {
    p_c: Vec<f64>,
    records: Vec<LmCmaRecord>,
    capacity: usize,
    n_steps: u64,
    a: f64,
    inverse_a: f64,
    c_c: f64,
    c1: f64,
    weights: Vec<f64>,
    mu_eff: f64,
    success_path: f64,
    previous_objectives: Option<Vec<f64>>,
}

impl LmCmaState {
    fn new(shared: &SharedParameters, capacity: usize) -> Self {
        let c1 = 1.0 / (10.0 * fs_math::det::ln(shared.n as f64 + 1.0));
        let a = fs_math::det::sqrt(1.0 - c1);
        // Loshchilov's reference uses ln(mu + 0.5) - ln(rank), which is the
        // same recombination rule already owned by SharedParameters for the
        // default mu = floor(lambda / 2). Reusing it prevents the limited-
        // memory implementation from silently drifting to ln(mu + 1).
        let weights = shared.weights.clone();
        let mu_eff = 1.0 / weights.iter().map(|weight| weight * weight).sum::<f64>();
        Self {
            p_c: vec![0.0; shared.n],
            records: Vec::with_capacity(capacity),
            capacity,
            n_steps: u64::try_from(capacity).unwrap_or(u64::MAX),
            a,
            inverse_a: 1.0 / a,
            c_c: 1.0 / capacity as f64,
            c1,
            weights,
            mu_eff,
            success_path: 0.0,
            previous_objectives: None,
        }
    }

    fn sample(&self, z: &[f64], output: &mut [f64]) {
        output.copy_from_slice(z);
        for record in &self.records {
            // A_{i+1} z = a A_i z + b_i p_i (v_i^T z): every rank-one
            // projection uses the original isotropic sample, while only the
            // accumulated A_i z is progressively transformed.
            let projection = dot(&record.v, z);
            for (out, &path) in output.iter_mut().zip(&record.p) {
                *out = self.a.mul_add(*out, record.b * path * projection);
            }
        }
    }

    #[cfg(test)]
    fn inverse_transform(&self, input: &[f64]) -> Vec<f64> {
        let mut output = input.to_vec();
        for record in &self.records {
            let projection = dot(&record.v, &output);
            for (out, &direction) in output.iter_mut().zip(&record.v) {
                *out = self
                    .inverse_a
                    .mul_add(*out, -record.d * direction * projection);
            }
        }
        output
    }

    fn update(
        &mut self,
        shared: &SharedParameters,
        generation: u64,
        order: &[usize],
        ys: &[Vec<f64>],
        objectives: &[f64],
        sigma: &mut f64,
    ) -> Result<Vec<f64>, CmaFamilyError> {
        let y_w = weighted_sum(&self.weights, order, ys, shared.n);
        let path_scale = fs_math::det::sqrt(self.c_c * (2.0 - self.c_c) * self.mu_eff);
        for coordinate in 0..shared.n {
            self.p_c[coordinate] =
                (1.0 - self.c_c).mul_add(self.p_c[coordinate], path_scale * y_w[coordinate]);
        }
        ensure_vector_finite(&self.p_c, "LM-CMA path")?;
        let record = LmCmaRecord {
            generation,
            p: self.p_c.clone(),
            v: vec![0.0; shared.n],
            b: 0.0,
            d: 0.0,
        };
        let recompute_from = self.insert_record(record);
        self.recompute_inverse_directions(recompute_from)?;
        if let Some(previous) = &self.previous_objectives {
            let success = population_success_rule(objectives, previous);
            self.success_path = 0.7f64.mul_add(self.success_path, 0.3 * success);
            *sigma *= fs_math::det::exp(self.success_path);
            ensure_positive(*sigma, "LM-CMA step size", 0)?;
        }
        self.previous_objectives = Some(objectives.to_vec());
        Ok(y_w)
    }

    fn insert_record(&mut self, record: LmCmaRecord) -> usize {
        if self.records.len() < self.capacity {
            let inserted = self.records.len();
            self.records.push(record);
            return inserted;
        }
        let mut closest = 0usize;
        let mut closest_gap = u64::MAX;
        for index in 0..self.records.len().saturating_sub(1) {
            let gap = self.records[index + 1]
                .generation
                .saturating_sub(self.records[index].generation);
            if gap < closest_gap {
                closest_gap = gap;
                closest = index;
            }
        }
        let remove = if closest_gap < self.n_steps {
            closest + 1
        } else {
            0
        };
        self.records.remove(remove);
        self.records.push(record);
        remove
    }

    /// Rebuild every inverse direction whose prefix transform changed.
    ///
    /// The corrected July 2014 LM-CMA reference recomputes `v = A^-1 p`
    /// after its temporal-memory replacement. Retaining the old `v` values
    /// reproduces the corruption in the original publication: once a stored
    /// predecessor is removed, every later inverse direction still encodes
    /// that deleted factor and `A`/`A^-1` cease to be inverses. In a 5,040-D
    /// plateaued objective this made sampled coordinates grow past 1e170.
    /// Records before `start` have an unchanged prefix; rebuilding only the
    /// suffix is both exact and the minimal corrected work.
    fn recompute_inverse_directions(&mut self, start: usize) -> Result<(), CmaFamilyError> {
        let ratio = self.c1 / (1.0 - self.c1);
        for index in start..self.records.len() {
            let (prefix, suffix) = self.records.split_at_mut(index);
            let record = &mut suffix[0];
            let mut inverse_direction = record.p.clone();
            for prior in prefix {
                let projection = dot(&prior.v, &inverse_direction);
                for (value, &direction) in inverse_direction.iter_mut().zip(&prior.v) {
                    *value = self
                        .inverse_a
                        .mul_add(*value, -prior.d * direction * projection);
                }
            }
            ensure_vector_finite(&inverse_direction, "LM-CMA inverse transform")?;
            let inverse_norm_squared = norm_squared(&inverse_direction);
            let root = fs_math::det::sqrt(ratio.mul_add(inverse_norm_squared, 1.0));
            let b = self.a * ratio / (root + 1.0);
            let d = ratio / (self.a * root * (root + 1.0));
            ensure_finite(b, "LM-CMA factor coefficient", index)?;
            ensure_finite(d, "LM-CMA inverse coefficient", index)?;
            record.v = inverse_direction;
            record.b = b;
            record.d = d;
        }
        Ok(())
    }
}

#[derive(Debug, Clone)]
struct LmMaState {
    p_sigma: Vec<f64>,
    paths: Vec<Vec<f64>>,
    c_d: Vec<f64>,
    c_c: Vec<f64>,
    c_sigma: f64,
}

impl LmMaState {
    fn new(shared: &SharedParameters, capacity: usize) -> Self {
        let mut c_d = Vec::with_capacity(capacity);
        let mut c_c = Vec::with_capacity(capacity);
        let mut d_denominator = shared.n as f64;
        let mut c_denominator = shared.n as f64;
        for _ in 0..capacity {
            c_d.push((1.0 / d_denominator).min(1.0));
            c_c.push((shared.lambda as f64 / c_denominator).min(1.0));
            d_denominator *= 1.5;
            c_denominator *= 4.0;
        }
        Self {
            p_sigma: vec![0.0; shared.n],
            paths: vec![vec![0.0; shared.n]; capacity],
            c_d,
            c_c,
            c_sigma: (2.0 * shared.lambda as f64 / shared.n as f64).min(1.0),
        }
    }

    fn sample(&self, generation: u64, z: &[f64], output: &mut [f64]) {
        output.copy_from_slice(z);
        let active = usize::try_from(generation)
            .unwrap_or(usize::MAX)
            .min(self.paths.len());
        for path_index in 0..active {
            let projection = dot(&self.paths[path_index], output);
            let rate = self.c_d[path_index];
            for (out, &path) in output.iter_mut().zip(&self.paths[path_index]) {
                *out = (1.0 - rate).mul_add(*out, rate * path * projection);
            }
        }
    }

    fn update(
        &mut self,
        shared: &SharedParameters,
        order: &[usize],
        zs: &[Vec<f64>],
        ys: &[Vec<f64>],
        sigma: &mut f64,
    ) -> Result<Vec<f64>, CmaFamilyError> {
        let y_w = weighted_sum(&shared.weights, order, ys, shared.n);
        let z_w = weighted_sum(&shared.weights, order, zs, shared.n);
        let path_scale = fs_math::det::sqrt(self.c_sigma * (2.0 - self.c_sigma) * shared.mu_eff);
        for coordinate in 0..shared.n {
            self.p_sigma[coordinate] = (1.0 - self.c_sigma)
                .mul_add(self.p_sigma[coordinate], path_scale * z_w[coordinate]);
        }
        let exponent = self.c_sigma * 0.5 * (norm_squared(&self.p_sigma) / shared.n as f64 - 1.0);
        *sigma *= fs_math::det::exp(exponent);
        ensure_positive(*sigma, "LM-MA step size", 0)?;
        for (path_index, path) in self.paths.iter_mut().enumerate() {
            let rate = self.c_c[path_index];
            let scale = fs_math::det::sqrt(shared.mu_eff * rate * (2.0 - rate));
            for coordinate in 0..shared.n {
                path[coordinate] = (1.0 - rate).mul_add(path[coordinate], scale * z_w[coordinate]);
            }
            ensure_vector_finite(path, "LM-MA direction path")?;
        }
        Ok(y_w)
    }
}

#[derive(Debug, Clone)]
enum Strategy {
    Full(FullState),
    Separable(SeparableState),
    LmCma(LmCmaState),
    LmMa(LmMaState),
}

impl Strategy {
    fn sample(&self, generation: u64, z: &[f64], output: &mut [f64]) {
        match self {
            Self::Full(state) => state.sample(z, output),
            Self::Separable(state) => state.sample(z, output),
            Self::LmCma(state) => state.sample(z, output),
            Self::LmMa(state) => state.sample(generation, z, output),
        }
    }

    #[allow(clippy::too_many_arguments)] // dispatch mirrors the four published recurrences
    fn update(
        &mut self,
        shared: &SharedParameters,
        generation: u64,
        order: &[usize],
        zs: &[Vec<f64>],
        ys: &[Vec<f64>],
        objectives: &[f64],
        sigma: &mut f64,
    ) -> Result<Vec<f64>, CmaFamilyError> {
        match self {
            Self::Full(state) => state.update(shared, generation, order, zs, ys, sigma),
            Self::Separable(state) => state.update(shared, generation, order, zs, ys, sigma),
            Self::LmCma(state) => state.update(shared, generation, order, ys, objectives, sigma),
            Self::LmMa(state) => state.update(shared, order, zs, ys, sigma),
        }
    }
}

/// Stateful, deterministic optimizer implementing the unified ask/tell API.
#[derive(Debug, Clone)]
pub struct CmaOptimizer {
    family: CmaFamily,
    admission: CmaAdmission,
    shared: SharedParameters,
    mean: Vec<f64>,
    sigma: f64,
    stream: Stream,
    strategy: Strategy,
    generation: u64,
    evaluations: usize,
    pending: Option<(u64, u64)>,
    best: Option<CmaBest>,
}

impl CmaOptimizer {
    /// Validate and allocate a new optimizer.
    pub fn new(config: CmaConfig) -> Result<Self, CmaFamilyError> {
        let admission = admit_cma(&config)?;
        let shared = SharedParameters::new(admission);
        let strategy = match config.family {
            CmaFamily::Full => Strategy::Full(FullState::new(&shared)),
            CmaFamily::Separable => Strategy::Separable(SeparableState::new(&shared)),
            CmaFamily::LmCma => Strategy::LmCma(LmCmaState::new(
                &shared,
                admission.complexity.memory_capacity,
            )),
            CmaFamily::LmMa => Strategy::LmMa(LmMaState::new(
                &shared,
                admission.complexity.memory_capacity,
            )),
        };
        Ok(Self {
            family: config.family,
            admission,
            shared,
            mean: config.mean,
            sigma: config.sigma,
            stream: StreamKey {
                seed: config.seed,
                kernel: CMA_FAMILY_STREAM_KERNEL,
                tile: 0,
            }
            .stream(),
            strategy,
            generation: 0,
            evaluations: 0,
            pending: None,
            best: None,
        })
    }

    /// Return the immutable admission receipt used to allocate this state.
    #[must_use]
    pub const fn admission(&self) -> CmaAdmission {
        self.admission
    }

    /// Generate exactly one complete candidate population.
    pub fn ask(&mut self) -> Result<CmaAsk, CmaFamilyError> {
        if let Some((generation, _)) = self.pending {
            return Err(CmaFamilyError::AskAlreadyPending { generation });
        }
        let remaining = self
            .admission
            .admitted_evaluations
            .saturating_sub(self.evaluations);
        if remaining < self.shared.lambda {
            return Err(CmaFamilyError::BudgetExhausted {
                remaining,
                required: self.shared.lambda,
            });
        }
        let mut stream = self.stream;
        let mut candidates = Vec::with_capacity(self.shared.lambda);
        let mut isotropic_steps = Vec::with_capacity(self.shared.lambda);
        let mut distribution_steps = Vec::with_capacity(self.shared.lambda);
        for candidate in 0..self.shared.lambda {
            let mut z = vec![0.0; self.shared.n];
            for value in &mut z {
                *value = stream.next_normal();
            }
            let mut y = vec![0.0; self.shared.n];
            self.strategy.sample(self.generation, &z, &mut y);
            let mut point = vec![0.0; self.shared.n];
            for coordinate in 0..self.shared.n {
                point[coordinate] = self.sigma.mul_add(y[coordinate], self.mean[coordinate]);
                ensure_finite(point[coordinate], "candidate generation", coordinate).map_err(
                    |error| match error {
                        CmaFamilyError::NumericalFailure { bits, .. } => {
                            CmaFamilyError::NumericalFailure {
                                stage: "candidate generation",
                                coordinate: candidate
                                    .checked_mul(self.shared.n)
                                    .and_then(|base| base.checked_add(coordinate))
                                    .unwrap_or(usize::MAX),
                                bits,
                            }
                        }
                        other => other,
                    },
                )?;
            }
            isotropic_steps.push(z);
            distribution_steps.push(y);
            candidates.push(point);
        }
        let signature = batch_signature(self.generation, &candidates);
        self.stream = stream;
        self.pending = Some((self.generation, signature));
        Ok(CmaAsk {
            generation: self.generation,
            signature,
            candidates,
            isotropic_steps,
            distribution_steps,
        })
    }

    /// Complete one outstanding generation using finite objective values.
    pub fn tell(
        &mut self,
        batch: &CmaAsk,
        objectives: &[f64],
    ) -> Result<CmaSnapshot, CmaFamilyError> {
        let Some((expected_generation, expected_signature)) = self.pending else {
            return Err(CmaFamilyError::NoPendingAsk);
        };
        if batch.generation != expected_generation {
            return Err(CmaFamilyError::GenerationMismatch {
                expected: expected_generation,
                actual: batch.generation,
            });
        }
        if batch.signature != expected_signature
            || batch_signature(batch.generation, &batch.candidates) != expected_signature
        {
            return Err(CmaFamilyError::BatchMismatch);
        }
        if objectives.len() != self.shared.lambda {
            return Err(CmaFamilyError::ObjectiveCount {
                expected: self.shared.lambda,
                actual: objectives.len(),
            });
        }
        for (candidate, &objective) in objectives.iter().enumerate() {
            if !objective.is_finite() {
                return Err(CmaFamilyError::NonFiniteObjective {
                    candidate,
                    bits: objective.to_bits(),
                });
            }
        }
        let mut order: Vec<usize> = (0..self.shared.lambda).collect();
        order.sort_by(|&left, &right| {
            objectives[left]
                .total_cmp(&objectives[right])
                .then(left.cmp(&right))
        });
        let mut strategy = self.strategy.clone();
        let mut sigma = self.sigma;
        let y_w = strategy.update(
            &self.shared,
            self.generation,
            &order,
            &batch.isotropic_steps,
            &batch.distribution_steps,
            objectives,
            &mut sigma,
        )?;
        let mut mean = self.mean.clone();
        for coordinate in 0..self.shared.n {
            mean[coordinate] = self.sigma.mul_add(y_w[coordinate], mean[coordinate]);
            ensure_finite(mean[coordinate], "mean update", coordinate)?;
        }
        let winner = order[0];
        let next_best = match &self.best {
            Some(best) if !objectives[winner].total_cmp(&best.objective).is_lt() => {
                Some(best.clone())
            }
            _ => Some(CmaBest {
                point: batch.candidates[winner].clone(),
                objective: objectives[winner],
                generation: self.generation,
                candidate: winner,
            }),
        };
        self.strategy = strategy;
        self.sigma = sigma;
        self.mean = mean;
        self.best = next_best;
        self.evaluations += self.shared.lambda;
        self.generation += 1;
        self.pending = None;
        Ok(self.snapshot())
    }

    /// Snapshot current diagnostics without exposing dense internal matrices.
    #[must_use]
    pub fn snapshot(&self) -> CmaSnapshot {
        let shape = match &self.strategy {
            Strategy::Full(state) => CmaShapeSnapshot::Full {
                diagonal: (0..self.shared.n)
                    .map(|index| state.covariance[index * self.shared.n + index])
                    .collect(),
                min_eigenvalue: state.axis_scales.first().map_or(1.0, |value| value * value),
                max_eigenvalue: state.axis_scales.last().map_or(1.0, |value| value * value),
                negative_weight_count: state
                    .all_weights
                    .iter()
                    .filter(|&&weight| weight < 0.0)
                    .count(),
            },
            Strategy::Separable(state) => CmaShapeSnapshot::Diagonal {
                variances: state.variances.clone(),
                negative_weight_count: state
                    .all_weights
                    .iter()
                    .filter(|&&weight| weight < 0.0)
                    .count(),
            },
            Strategy::LmCma(state) => CmaShapeSnapshot::LimitedMemory {
                vectors: state.records.len(),
                capacity: state.capacity,
                direction_norms: state.records.iter().map(|record| norm(&record.p)).collect(),
            },
            Strategy::LmMa(state) => CmaShapeSnapshot::LimitedMemory {
                vectors: self.generation.min(state.paths.len() as u64) as usize,
                capacity: state.paths.len(),
                direction_norms: state
                    .paths
                    .iter()
                    .take(self.generation.min(state.paths.len() as u64) as usize)
                    .map(|path| norm(path))
                    .collect(),
            },
        };
        CmaSnapshot {
            family: self.family,
            generation: self.generation,
            evaluations: self.evaluations,
            mean: self.mean.clone(),
            sigma: self.sigma,
            best: self.best.clone(),
            shape,
            complexity: self.admission.complexity,
        }
    }
}

fn normalize(weights: &mut [f64]) {
    let sum: f64 = weights.iter().sum();
    for weight in weights {
        *weight /= sum;
    }
}

fn active_weights(lambda: usize, mu_eff: f64, n: f64, c1: f64, c_mu: f64) -> Vec<f64> {
    let midpoint = f64::midpoint(lambda as f64, 1.0);
    let mut weights: Vec<f64> = (1..=lambda)
        .map(|rank| fs_math::det::ln(midpoint) - fs_math::det::ln(rank as f64))
        .collect();
    let positive_sum: f64 = weights.iter().copied().filter(|weight| *weight > 0.0).sum();
    for weight in weights.iter_mut().filter(|weight| **weight > 0.0) {
        *weight /= positive_sum;
    }
    let negative_sum: f64 = weights
        .iter()
        .copied()
        .filter(|weight| *weight < 0.0)
        .map(f64::abs)
        .sum();
    let negative_square_sum: f64 = weights
        .iter()
        .copied()
        .filter(|weight| *weight < 0.0)
        .map(|weight| weight * weight)
        .sum();
    if negative_sum > 0.0 && c_mu > 0.0 {
        let negative_mu_eff = negative_sum * negative_sum / negative_square_sum;
        let alpha = (1.0 + c1 / c_mu)
            .min(1.0 + 2.0 * negative_mu_eff / (mu_eff + 2.0))
            .min((1.0 - c1 - c_mu) / (n * c_mu));
        for weight in weights.iter_mut().filter(|weight| **weight < 0.0) {
            *weight *= alpha / negative_sum;
        }
    }
    weights
}

fn active_covariance_weights(
    weights: &[f64],
    order: &[usize],
    isotropic_steps: &[Vec<f64>],
    dimension: usize,
) -> Vec<f64> {
    weights
        .iter()
        .zip(order)
        .map(|(&weight, &candidate)| {
            if weight >= 0.0 {
                return weight;
            }
            let norm_squared = norm_squared(&isotropic_steps[candidate]);
            if norm_squared > f64::MIN_POSITIVE {
                weight * dimension as f64 / norm_squared
            } else {
                0.0
            }
        })
        .collect()
}

fn weighted_sum(
    weights: &[f64],
    order: &[usize],
    vectors: &[Vec<f64>],
    dimension: usize,
) -> Vec<f64> {
    let mut output = vec![0.0; dimension];
    for (&weight, &candidate) in weights.iter().zip(order) {
        for (out, &value) in output.iter_mut().zip(&vectors[candidate]) {
            *out = weight.mul_add(value, *out);
        }
    }
    output
}

fn norm_squared(vector: &[f64]) -> f64 {
    vector
        .iter()
        .fold(0.0, |sum, value| value.mul_add(*value, sum))
}

fn norm(vector: &[f64]) -> f64 {
    fs_math::det::sqrt(norm_squared(vector))
}

fn dot(left: &[f64], right: &[f64]) -> f64 {
    left.iter()
        .zip(right)
        .fold(0.0, |sum, (&a, &b)| a.mul_add(b, sum))
}

fn population_success_rule(current: &[f64], previous: &[f64]) -> f64 {
    let lambda = current.len();
    let mut combined: Vec<(f64, bool, usize)> = previous
        .iter()
        .copied()
        .enumerate()
        .map(|(index, value)| (value, false, index))
        .chain(
            current
                .iter()
                .copied()
                .enumerate()
                .map(|(index, value)| (value, true, index)),
        )
        .collect();
    combined.sort_by(|left, right| {
        left.0
            .total_cmp(&right.0)
            .then(left.1.cmp(&right.1))
            .then(left.2.cmp(&right.2))
    });
    let mut current_rank_sum = 0usize;
    let mut previous_rank_sum = 0usize;
    for (position, (_, is_current, _)) in combined.iter().enumerate() {
        let descending_rank = 2 * lambda - position;
        if *is_current {
            current_rank_sum += descending_rank;
        } else {
            previous_rank_sum += descending_rank;
        }
    }
    (current_rank_sum as f64 - previous_rank_sum as f64) / (lambda * lambda) as f64 - 0.25
}

fn batch_signature(generation: u64, candidates: &[Vec<f64>]) -> u64 {
    let mut hash = 0xcbf2_9ce4_8422_2325u64 ^ generation;
    for point in candidates {
        for &value in point {
            hash ^= value.to_bits();
            hash = hash.wrapping_mul(0x0000_0100_0000_01b3);
        }
    }
    hash
}

fn ensure_finite(value: f64, stage: &'static str, coordinate: usize) -> Result<(), CmaFamilyError> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(CmaFamilyError::NumericalFailure {
            stage,
            coordinate,
            bits: value.to_bits(),
        })
    }
}

fn ensure_positive(
    value: f64,
    stage: &'static str,
    coordinate: usize,
) -> Result<(), CmaFamilyError> {
    if value > 0.0 {
        ensure_finite(value, stage, coordinate)
    } else {
        Err(CmaFamilyError::NumericalFailure {
            stage,
            coordinate,
            bits: value.to_bits(),
        })
    }
}

fn ensure_vector_finite(vector: &[f64], stage: &'static str) -> Result<(), CmaFamilyError> {
    for (coordinate, &value) in vector.iter().enumerate() {
        ensure_finite(value, stage, coordinate)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sphere(point: &[f64]) -> f64 {
        norm_squared(point)
    }

    fn cigar(point: &[f64]) -> f64 {
        point[0].mul_add(
            point[0],
            1.0e6 * point[1..].iter().map(|value| value * value).sum::<f64>(),
        )
    }

    fn rotated_quadratic(point: &[f64]) -> f64 {
        assert_eq!(point.len(), 8);
        let mut transformed = point.to_vec();
        let scale = 1.0 / fs_math::det::sqrt(8.0);
        let mut width = 1usize;
        while width < transformed.len() {
            for base in (0..transformed.len()).step_by(width * 2) {
                for offset in 0..width {
                    let left = transformed[base + offset];
                    let right = transformed[base + width + offset];
                    transformed[base + offset] = left + right;
                    transformed[base + width + offset] = left - right;
                }
            }
            width *= 2;
        }
        transformed
            .iter()
            .enumerate()
            .map(|(index, value)| {
                let exponent = index as f64 / 7.0;
                fs_math::det::exp(fs_math::det::ln(1.0e4) * exponent) * (scale * value).powi(2)
            })
            .sum()
    }

    fn complete_generation(
        optimizer: &mut CmaOptimizer,
        objective: fn(&[f64]) -> f64,
    ) -> CmaSnapshot {
        let batch = optimizer.ask().expect("generation must be admitted");
        let values: Vec<f64> = batch
            .candidates()
            .iter()
            .map(|point| objective(point))
            .collect();
        optimizer
            .tell(&batch, &values)
            .expect("finite generation must update")
    }

    fn assert_same_bits(left: &CmaSnapshot, right: &CmaSnapshot) {
        assert_eq!(left.family, right.family);
        assert_eq!(left.generation, right.generation);
        assert_eq!(left.evaluations, right.evaluations);
        assert_eq!(left.sigma.to_bits(), right.sigma.to_bits());
        for (&a, &b) in left.mean.iter().zip(&right.mean) {
            assert_eq!(a.to_bits(), b.to_bits());
        }
        match (&left.best, &right.best) {
            (Some(a), Some(b)) => {
                assert_eq!(a.objective.to_bits(), b.objective.to_bits());
                assert_eq!(a.generation, b.generation);
                assert_eq!(a.candidate, b.candidate);
                for (&x, &y) in a.point.iter().zip(&b.point) {
                    assert_eq!(x.to_bits(), y.to_bits());
                }
            }
            (None, None) => {}
            _ => panic!("replay best-state presence differs"),
        }
        assert_eq!(left.shape, right.shape);
        assert_eq!(left.complexity, right.complexity);
    }

    #[test]
    fn active_weight_fixture_matches_hansen_equations() {
        let weights = active_weights(7, 2.6, 5.0, 0.081, 0.28);
        let expected = [
            0.585_645_106_509_765_1,
            0.292_822_553_254_882_54,
            0.121_532_340_235_352_46,
            0.0,
            -0.085_715_365_120_049_72,
            -0.155_749_917_845_457_3,
            -0.214_963_288_463_064_4,
        ];
        for (&actual, expected) in weights.iter().zip(expected) {
            assert!((actual - expected).abs() < 2.0e-15);
        }
        assert!((weights[..3].iter().sum::<f64>() - 1.0).abs() < 2.0e-15);
        assert_eq!(weights.iter().filter(|&&weight| weight < 0.0).count(), 3);
    }

    #[test]
    fn full_small_population_c1_matches_pycma_reference() {
        let mut config = CmaConfig::standard(CmaFamily::Full, vec![0.0], 1.0, 4, 1);
        config.population_size = Some(4);
        let shared = SharedParameters::new(admit_cma(&config).expect("admission"));
        let full = FullState::new(&shared);
        let expected = (4.0 / 6.0) * 2.0 / ((1.0_f64 + 1.3).powi(2) + shared.mu_eff);
        assert!((full.c1 - expected).abs() < f64::EPSILON);
    }

    #[test]
    fn separable_active_update_matches_diagonal_hansen_equation() {
        let mut config = CmaConfig::standard(CmaFamily::Separable, vec![0.0; 3], 1.0, 6, 1);
        config.population_size = Some(6);
        let shared = SharedParameters::new(admit_cma(&config).expect("admission"));
        let mut state = SeparableState::new(&shared);
        assert_eq!(
            state
                .all_weights
                .iter()
                .filter(|&&weight| weight < 0.0)
                .count(),
            3
        );
        let zs = vec![
            vec![0.5, -1.0, 0.25],
            vec![1.5, 0.25, -0.75],
            vec![-0.75, 0.5, 1.25],
            vec![0.125, -0.25, 0.5],
            vec![-1.25, 0.75, -0.5],
            vec![0.875, 1.125, -1.5],
        ];
        let ys = zs.clone();
        let order = [0, 1, 2, 3, 4, 5];

        let y_w = weighted_sum(&shared.weights, &order, &ys, shared.n);
        let z_w = weighted_sum(&shared.weights, &order, &zs, shared.n);
        let path_scale = fs_math::det::sqrt(state.c_sigma * (2.0 - state.c_sigma) * shared.mu_eff);
        let expected_p_sigma: Vec<f64> = z_w.iter().map(|value| path_scale * value).collect();
        let path_norm = norm(&expected_p_sigma);
        let h_sigma = path_norm
            / fs_math::det::sqrt(1.0 - fs_math::det::powi(1.0 - state.c_sigma, 2))
            < (1.4 + 2.0 / (shared.n as f64 + 1.0)) * shared.chi_n;
        let covariance_path_scale =
            fs_math::det::sqrt(state.c_c * (2.0 - state.c_c) * shared.mu_eff);
        let expected_p_c: Vec<f64> = y_w
            .iter()
            .map(|value| {
                if h_sigma {
                    covariance_path_scale * value
                } else {
                    0.0
                }
            })
            .collect();
        let delta_h = if h_sigma {
            0.0
        } else {
            state.c_c * (2.0 - state.c_c)
        };
        let old_factor = 1.0 + state.c1 * delta_h
            - state.c1
            - state.c_mu * state.all_weights.iter().sum::<f64>();
        let expected_variances: Vec<f64> = (0..shared.n)
            .map(|coordinate| {
                let rank_mu = order
                    .iter()
                    .enumerate()
                    .map(|(rank, &candidate)| {
                        let weight = state.all_weights[rank];
                        let effective = if weight < 0.0 {
                            weight * shared.n as f64 / norm_squared(&zs[candidate])
                        } else {
                            weight
                        };
                        effective * ys[candidate][coordinate].powi(2)
                    })
                    .sum::<f64>();
                old_factor + state.c1 * expected_p_c[coordinate].powi(2) + state.c_mu * rank_mu
            })
            .collect();

        let mut sigma = 1.0;
        state
            .update(&shared, 0, &order, &zs, &ys, &mut sigma)
            .expect("active diagonal update");
        for ((&actual, &expected), (&actual_path, &expected_path)) in state
            .variances
            .iter()
            .zip(&expected_variances)
            .zip(state.p_c.iter().zip(&expected_p_c))
        {
            assert!((actual - expected).abs() < 2.0e-15);
            assert!((actual_path - expected_path).abs() < 2.0e-15);
        }
    }

    #[test]
    fn full_csa_whitening_rotates_with_a_nonidentity_eigenbasis() {
        let mut config = CmaConfig::standard(CmaFamily::Full, vec![0.0; 2], 1.0, 4, 1);
        config.population_size = Some(4);
        let shared = SharedParameters::new(admit_cma(&config).expect("admission"));
        let mut identity = FullState::new(&shared);
        let mut rotated = identity.clone();
        // Exact 90-degree rotation: Q [x,y] = [-y,x]. Covariance remains I,
        // while the square-root factor is a genuinely nonidentity basis.
        rotated.basis = vec![0.0, -1.0, 1.0, 0.0];
        let zs = vec![
            vec![0.5, -1.0],
            vec![1.5, 0.25],
            vec![-0.75, 0.5],
            vec![0.125, -0.25],
        ];
        let identity_steps = zs.clone();
        let rotated_steps: Vec<Vec<f64>> = zs.iter().map(|z| vec![-z[1], z[0]]).collect();
        let order = [0, 1, 2, 3];
        let mut identity_sigma = 1.0;
        let mut rotated_sigma = 1.0;
        identity
            .update(
                &shared,
                0,
                &order,
                &zs,
                &identity_steps,
                &mut identity_sigma,
            )
            .expect("identity update");
        rotated
            .update(&shared, 0, &order, &zs, &rotated_steps, &mut rotated_sigma)
            .expect("rotated update");
        assert!((rotated.p_sigma[0] + identity.p_sigma[1]).abs() < 2.0e-15);
        assert!((rotated.p_sigma[1] - identity.p_sigma[0]).abs() < 2.0e-15);
        assert!((rotated_sigma - identity_sigma).abs() < 2.0e-15);
    }

    #[test]
    fn lm_ma_two_path_transform_matches_algorithm_one_progressively() {
        let state = LmMaState {
            p_sigma: vec![0.0; 2],
            paths: vec![vec![1.0, 2.0], vec![-1.0, 0.5]],
            c_d: vec![0.2, 0.3],
            c_c: vec![0.0, 0.0],
            c_sigma: 0.1,
        };
        let z = [0.4, -0.7];

        let mut generation_zero = [0.0; 2];
        state.sample(0, &z, &mut generation_zero);
        assert_eq!(generation_zero.map(f64::to_bits), z.map(f64::to_bits));

        let mut generation_one = [0.0; 2];
        state.sample(1, &z, &mut generation_one);
        // d1 = .8 z + .2 p1 (p1^T z), with p1^T z = -1.
        let first = [0.12, -0.96];
        for (&actual, expected) in generation_one.iter().zip(first) {
            assert!((actual - expected).abs() < 2.0e-15);
        }

        let mut generation_two = [0.0; 2];
        state.sample(2, &z, &mut generation_two);
        // Algorithm 1 line 10 reuses the progressively transformed d:
        // p2^T d1 = -.6; d2 = .7 d1 + .3 p2 (-.6).
        let second = [0.264, -0.762];
        for (&actual, expected) in generation_two.iter().zip(second) {
            assert!((actual - expected).abs() < 2.0e-15);
        }
    }

    #[test]
    fn reference_parameter_fixtures_match_sep_lmcma_and_lmma() {
        let config = CmaConfig::standard(CmaFamily::LmMa, vec![0.0; 100], 1.0, 16, 1);
        let admission = admit_cma(&CmaConfig {
            population_size: Some(16),
            memory: Some(16),
            ..config
        })
        .expect("fixture admission");
        let shared = SharedParameters::new(admission);
        let separable = SeparableState::new(&shared);
        let expected_cc =
            (1.0 + 0.01 + shared.mu_eff / 100.0) / (10.0 + 0.01 + 2.0 * shared.mu_eff / 100.0);
        assert!((separable.c_c - expected_cc).abs() < 1.0e-15);
        assert!((separable.c1 - 1.0 / (120.0 + shared.mu_eff / 100.0)).abs() < 1.0e-15);

        let lm_cma = LmCmaState::new(&shared, 16);
        assert!((lm_cma.c_c - 0.0625).abs() < f64::EPSILON);
        assert!((lm_cma.c1 - 1.0 / (10.0 * fs_math::det::ln(101.0))).abs() < f64::EPSILON);
        assert_eq!(
            lm_cma.weights.iter().map(|value| value.to_bits()).collect::<Vec<_>>(),
            shared.weights.iter().map(|value| value.to_bits()).collect::<Vec<_>>()
        );

        let lm_ma = LmMaState::new(&shared, 4);
        assert!((lm_ma.c_sigma - 0.32).abs() < f64::EPSILON);
        assert!((lm_ma.c_d[0] - 0.01).abs() < f64::EPSILON);
        assert!((lm_ma.c_d[1] - 1.0 / 150.0).abs() < f64::EPSILON);
        assert!((lm_ma.c_c[0] - 0.16).abs() < f64::EPSILON);
        assert!((lm_ma.c_c[1] - 0.04).abs() < f64::EPSILON);
    }

    #[test]
    fn reference_memory_default_is_dimension_based_when_population_is_overridden() {
        for family in [CmaFamily::LmCma, CmaFamily::LmMa] {
            let mut config = CmaConfig::standard(family, vec![0.0; 100], 1.0, 16, 1);
            config.population_size = Some(16);
            let admission = admit_cma(&config).expect("limited-memory admission");
            assert_eq!(admission.complexity.memory_capacity, 17);
        }
    }

    #[test]
    fn complexity_receipts_count_actual_float_storage_and_workspace() {
        let receipt = |family, memory| {
            let mut config = CmaConfig::standard(family, vec![0.0; 3], 1.0, 6, 1);
            config.population_size = Some(6);
            config.memory = memory;
            admit_cma(&config).expect("complexity admission").complexity
        };

        let full = receipt(CmaFamily::Full, None);
        assert_eq!(full.persistent_scalars, 51);
        assert_eq!(full.pending_generation_scalars, 54);
        assert_eq!(full.update_workspace_scalars, 108);
        assert_eq!(full.dense_matrix_entries, 18);

        let separable = receipt(CmaFamily::Separable, None);
        assert_eq!(separable.persistent_scalars, 33);
        assert_eq!(separable.pending_generation_scalars, 54);
        assert_eq!(separable.update_workspace_scalars, 36);
        assert_eq!(separable.dense_matrix_entries, 0);

        let lm_cma = receipt(CmaFamily::LmCma, Some(2));
        assert_eq!(lm_cma.persistent_scalars, 47);
        assert_eq!(lm_cma.update_workspace_scalars, 64);
        assert_eq!(lm_cma.dense_matrix_entries, 0);

        let lm_ma = receipt(CmaFamily::LmMa, Some(2));
        assert_eq!(lm_ma.persistent_scalars, 27);
        assert_eq!(lm_ma.update_workspace_scalars, 21);
        assert_eq!(lm_ma.dense_matrix_entries, 0);
    }

    #[test]
    fn ask_tell_contract_and_exact_budget_are_enforced() {
        let mut config = CmaConfig::standard(CmaFamily::Separable, vec![1.0; 5], 0.5, 31, 9);
        config.population_size = Some(10);
        let mut optimizer = CmaOptimizer::new(config).expect("valid config");
        assert_eq!(optimizer.admission().admitted_evaluations, 30);

        let first = optimizer.ask().expect("first ask");
        assert!(matches!(
            optimizer.ask(),
            Err(CmaFamilyError::AskAlreadyPending { generation: 0 })
        ));
        assert!(matches!(
            optimizer.tell(&first, &[0.0; 9]),
            Err(CmaFamilyError::ObjectiveCount {
                expected: 10,
                actual: 9
            })
        ));
        let mut modified = first.clone();
        modified.candidates[0][0] += 1.0;
        assert!(matches!(
            optimizer.tell(&modified, &[0.0; 10]),
            Err(CmaFamilyError::BatchMismatch)
        ));
        optimizer
            .tell(&first, &[0.0; 10])
            .expect("valid retry after refusal");
        assert!(matches!(
            optimizer.tell(&first, &[0.0; 10]),
            Err(CmaFamilyError::NoPendingAsk)
        ));
        for _ in 0..2 {
            let batch = optimizer.ask().expect("budgeted ask");
            optimizer.tell(&batch, &[0.0; 10]).expect("budgeted tell");
        }
        assert_eq!(optimizer.snapshot().evaluations, 30);
        assert!(matches!(
            optimizer.ask(),
            Err(CmaFamilyError::BudgetExhausted {
                remaining: 0,
                required: 10
            })
        ));
    }

    #[test]
    fn validation_and_finite_refusals_are_typed_and_retryable() {
        assert!(matches!(
            CmaOptimizer::new(CmaConfig::standard(CmaFamily::Full, Vec::new(), 1.0, 10, 1)),
            Err(CmaFamilyError::EmptyMean)
        ));
        assert!(matches!(
            CmaOptimizer::new(CmaConfig::standard(
                CmaFamily::Separable,
                vec![f64::NAN],
                1.0,
                10,
                1
            )),
            Err(CmaFamilyError::NonFiniteMean { coordinate: 0, .. })
        ));
        let mut inapplicable = CmaConfig::standard(CmaFamily::Full, vec![0.0; 3], 1.0, 10, 1);
        inapplicable.memory = Some(2);
        assert!(matches!(
            admit_cma(&inapplicable),
            Err(CmaFamilyError::MemoryNotApplicable {
                family: CmaFamily::Full
            })
        ));

        let mut optimizer = CmaOptimizer::new(CmaConfig::standard(
            CmaFamily::LmMa,
            vec![0.0; 16],
            1.0,
            20,
            5,
        ))
        .expect("finite optimizer");
        let batch = optimizer.ask().expect("ask");
        let mut objectives = vec![0.0; batch.len()];
        objectives[2] = f64::INFINITY;
        assert!(matches!(
            optimizer.tell(&batch, &objectives),
            Err(CmaFamilyError::NonFiniteObjective { candidate: 2, .. })
        ));
        objectives[2] = 0.0;
        optimizer
            .tell(&batch, &objectives)
            .expect("finite retry consumes the still-pending batch");
    }

    #[test]
    fn objective_ties_select_lowest_candidate_indices() {
        let mut config = CmaConfig::standard(CmaFamily::Separable, vec![1.0; 6], 0.5, 8, 77);
        config.population_size = Some(8);
        let mut optimizer = CmaOptimizer::new(config).expect("optimizer");
        let batch = optimizer.ask().expect("ask");
        let expected_step = weighted_sum(
            &optimizer.shared.weights,
            &(0..optimizer.shared.lambda).collect::<Vec<_>>(),
            &batch.distribution_steps,
            optimizer.shared.n,
        );
        let expected_mean: Vec<f64> = expected_step
            .iter()
            .map(|step| 0.5f64.mul_add(*step, 1.0))
            .collect();
        let snapshot = optimizer
            .tell(&batch, &[1.0; 8])
            .expect("tied objectives are valid");
        for (&actual, &expected) in snapshot.mean.iter().zip(&expected_mean) {
            assert_eq!(actual.to_bits(), expected.to_bits());
        }
        assert_eq!(snapshot.best.expect("best").candidate, 0);
    }

    #[test]
    fn all_families_replay_bit_for_bit_and_respect_rank_transforms() {
        for family in [
            CmaFamily::Full,
            CmaFamily::Separable,
            CmaFamily::LmCma,
            CmaFamily::LmMa,
        ] {
            let config = CmaConfig::standard(family, vec![2.0; 12], 0.8, 400, 0x0A11_CE55);
            let mut raw = CmaOptimizer::new(config.clone()).expect("raw optimizer");
            let mut transformed = CmaOptimizer::new(config.clone()).expect("transformed optimizer");
            for _ in 0..20 {
                let raw_batch = raw.ask().expect("raw ask");
                let transformed_batch = transformed.ask().expect("transformed ask");
                for (left, right) in raw_batch
                    .candidates
                    .iter()
                    .zip(&transformed_batch.candidates)
                {
                    for (&a, &b) in left.iter().zip(right) {
                        assert_eq!(a.to_bits(), b.to_bits());
                    }
                }
                let values: Vec<f64> = raw_batch
                    .candidates
                    .iter()
                    .map(|point| sphere(point))
                    .collect();
                let mapped: Vec<f64> = values
                    .iter()
                    .map(|value| fs_math::det::ln(1.0 + value))
                    .collect();
                let raw_snapshot = raw.tell(&raw_batch, &values).expect("raw tell");
                let mapped_snapshot = transformed
                    .tell(&transformed_batch, &mapped)
                    .expect("mapped tell");
                for (&a, &b) in raw_snapshot.mean.iter().zip(&mapped_snapshot.mean) {
                    assert_eq!(a.to_bits(), b.to_bits());
                }
                assert_eq!(
                    raw_snapshot.sigma.to_bits(),
                    mapped_snapshot.sigma.to_bits()
                );
            }
            let mut replay = CmaOptimizer::new(config).expect("replay optimizer");
            let mut snapshot = replay.snapshot();
            for _ in 0..20 {
                snapshot = complete_generation(&mut replay, sphere);
            }
            assert_same_bits(&raw.snapshot(), &snapshot);
        }
    }

    #[test]
    fn covariance_representations_remain_positive() {
        for family in [CmaFamily::Full, CmaFamily::Separable] {
            let initial = cigar(&vec![3.0; 10]);
            let mut optimizer =
                CmaOptimizer::new(CmaConfig::standard(family, vec![3.0; 10], 1.0, 1_000, 17))
                    .expect("optimizer");
            for _ in 0..60 {
                complete_generation(&mut optimizer, cigar);
            }
            assert!(
                optimizer
                    .snapshot()
                    .best
                    .as_ref()
                    .is_some_and(|best| best.objective < initial * 0.1),
                "{family:?} must make positive progress on the axis-aligned cigar"
            );
            match optimizer.snapshot().shape {
                CmaShapeSnapshot::Full {
                    diagonal,
                    min_eigenvalue,
                    max_eigenvalue,
                    negative_weight_count,
                } => {
                    assert!(negative_weight_count > 0);
                    assert!(min_eigenvalue > 0.0 && max_eigenvalue >= min_eigenvalue);
                    assert!(
                        diagonal
                            .iter()
                            .all(|value| value.is_finite() && *value > 0.0)
                    );
                }
                CmaShapeSnapshot::Diagonal {
                    variances,
                    negative_weight_count,
                } => {
                    assert!(negative_weight_count > 0);
                    assert!(
                        variances
                            .iter()
                            .all(|value| value.is_finite() && *value > 0.0)
                    );
                }
                CmaShapeSnapshot::LimitedMemory { .. } => panic!("wrong covariance snapshot"),
            }
        }
    }

    #[test]
    fn limited_memory_is_capped_and_never_reports_dense_state() {
        for family in [CmaFamily::LmCma, CmaFamily::LmMa] {
            let mut config = CmaConfig::standard(family, vec![2.0; 32], 1.0, 500, 23);
            config.memory = Some(3);
            let mut optimizer = CmaOptimizer::new(config).expect("limited-memory optimizer");
            for _ in 0..12 {
                complete_generation(&mut optimizer, sphere);
            }
            let snapshot = optimizer.snapshot();
            assert_eq!(snapshot.complexity.dense_matrix_entries, 0);
            match snapshot.shape {
                CmaShapeSnapshot::LimitedMemory {
                    vectors,
                    capacity,
                    direction_norms,
                } => {
                    assert_eq!((vectors, capacity, direction_norms.len()), (3, 3, 3));
                    assert!(direction_norms.iter().all(|value| value.is_finite()));
                }
                _ => panic!("wrong limited-memory snapshot"),
            }
        }
    }

    #[test]
    fn lm_cma_inverse_tracks_temporal_record_replacement() {
        let mut config = CmaConfig::standard(CmaFamily::LmCma, vec![2.0; 8], 0.8, 80, 41);
        config.population_size = Some(4);
        config.memory = Some(3);
        let mut optimizer = CmaOptimizer::new(config).expect("LM-CMA optimizer");
        for _ in 0..12 {
            complete_generation(&mut optimizer, rotated_quadratic);
        }
        let Strategy::LmCma(state) = &optimizer.strategy else {
            panic!("wrong strategy");
        };
        let isotropic = [0.25, -0.5, 0.75, -1.0, 1.25, -1.5, 1.75, -2.0];
        let mut transformed = [0.0; 8];
        state.sample(&isotropic, &mut transformed);
        let recovered = state.inverse_transform(&transformed);
        for (&actual, expected) in recovered.iter().zip(isotropic) {
            assert!(
                (actual - expected).abs() < 2.0e-12,
                "A^-1(Az) must recover z after bounded-memory replacement"
            );
        }
    }

    #[test]
    fn lm_cma_5040d_plateau_stays_finite_through_memory_replacement() {
        const N: usize = 5_040;
        const POPULATION: usize = 16;
        const GENERATIONS: usize = 40;
        let mut config = CmaConfig::standard(
            CmaFamily::LmCma,
            vec![0.0; N],
            0.01,
            POPULATION * GENERATIONS,
            0x4731_5040,
        );
        config.population_size = Some(POPULATION);
        config.memory = Some(12);
        let mut optimizer = CmaOptimizer::new(config).expect("5,040-D LM-CMA optimizer");
        for _ in 0..GENERATIONS {
            let batch = optimizer.ask().expect("plateau generation must be admitted");
            let maximum_coordinate = batch
                .candidates()
                .iter()
                .flatten()
                .map(|value| value.abs())
                .fold(0.0, f64::max);
            assert!(
                maximum_coordinate.is_finite() && maximum_coordinate < 1.0e6,
                "LM-CMA transform escaped its finite search scale: {maximum_coordinate:e}"
            );
            let snapshot = optimizer
                .tell(&batch, &vec![1.0; POPULATION])
                .expect("finite plateau generation must update");
            assert!(snapshot.sigma.is_finite());
            assert!(snapshot.mean.iter().all(|value| value.is_finite()));
            let CmaShapeSnapshot::LimitedMemory {
                direction_norms, ..
            } = snapshot.shape
            else {
                panic!("wrong shape receipt");
            };
            assert!(direction_norms.iter().all(|value| value.is_finite()));
        }
    }

    #[test]
    fn dimension_5040_admits_and_runs_without_dense_quadratic_storage() {
        const N: usize = 5_040;
        for family in [CmaFamily::Separable, CmaFamily::LmCma, CmaFamily::LmMa] {
            let mut config = CmaConfig::standard(family, vec![0.25; N], 0.1, 64, 29);
            if matches!(family, CmaFamily::LmCma | CmaFamily::LmMa) {
                config.memory = Some(4);
            }
            let admission = admit_cma(&config).expect("large linear-state admission");
            assert_eq!(admission.complexity.dense_matrix_entries, 0);
            assert!(admission.complexity.persistent_scalars < 20 * N);
            let mut optimizer = CmaOptimizer::new(config).expect("large linear-state optimizer");
            let snapshot = complete_generation(&mut optimizer, sphere);
            assert_eq!(snapshot.evaluations, admission.population_size);
        }
    }

    #[test]
    fn all_families_make_real_progress_on_sphere() {
        for family in [
            CmaFamily::Full,
            CmaFamily::Separable,
            CmaFamily::LmCma,
            CmaFamily::LmMa,
        ] {
            let mut optimizer =
                CmaOptimizer::new(CmaConfig::standard(family, vec![4.0; 16], 1.5, 4_000, 101))
                    .expect("optimizer");
            for _ in 0..200 {
                complete_generation(&mut optimizer, sphere);
            }
            let best = optimizer
                .snapshot()
                .best
                .expect("best observation")
                .objective;
            assert!(best < 2.0, "{family:?} sphere best {best}");
        }
    }

    #[test]
    fn full_cma_handles_rotated_ill_conditioned_quadratic() {
        let initial = rotated_quadratic(&[3.0; 8]);
        let mut optimizer = CmaOptimizer::new(CmaConfig::standard(
            CmaFamily::Full,
            vec![3.0; 8],
            1.0,
            5_000,
            0xBEEF,
        ))
        .expect("full optimizer");
        for _ in 0..350 {
            complete_generation(&mut optimizer, rotated_quadratic);
        }
        let best = optimizer
            .snapshot()
            .best
            .expect("best observation")
            .objective;
        assert!(
            best < initial * 1.0e-4,
            "rotated best {best}, initial {initial}"
        );
    }
}
