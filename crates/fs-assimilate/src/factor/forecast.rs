//! Linear Gaussian time prediction through square-root right rotations.
//!
//! For `x_next = F x + b + w`, independent zero-mean `w ~ N(0,Q)`,
//! triangularize `[F U sqrt(D), U_q sqrt(D_q)]` by orthogonal right Givens
//! rotations. Its row Gram matrix is `F P F^T + Q`, but that matrix is never
//! formed or differenced. Positive pivots define a new unit-upper UD factor.
//! Exact zero rows retain zero pivots; no tolerance chooses a numerical rank.
//! These are numerical Gaussian predictions, not validated-flow enclosures.

use fs_blake3::{ContentHash, DomainHasher};

mod nonlinear;
pub use nonlinear::NonlinearDynamics;

use super::{FactorBelief, compensated_dot, packed_upper_index, packed_upper_len, try_f64_vec};
use crate::{
    AssimError, Cx, WorkPlan, WorkProgress, canonical_f64_bits, canonicalize_zero,
    validate_leaf_identity, validate_state_dimension,
};

/// One checked discrete-time affine transition and independent process noise.
///
/// The matrix and offset describe one complete step, not a derivative to be
/// multiplied by time again. The caller owns units, discretization accuracy,
/// independence of process noise from the prior, and the noise model's truth.
#[derive(Debug, Clone, PartialEq)]
pub struct LinearDynamics {
    transition: Vec<Vec<f64>>,
    offset: Vec<f64>,
    noise: FactorBelief,
    identity: String,
    content_id: ContentHash,
}

impl LinearDynamics {
    /// Admit a square finite transition, finite affine offset, and zero-mean
    /// factor-form process noise. Singular/zero noise is allowed without
    /// adding jitter. All three dimensions must agree and fit the factor cap.
    pub fn new(
        mut transition: Vec<Vec<f64>>,
        mut offset: Vec<f64>,
        noise: FactorBelief,
        identity: impl Into<String>,
        cx: &Cx<'_>,
    ) -> Result<Self, AssimError> {
        let n = noise.dim();
        validate_state_dimension(n)?;
        if transition.len() != n || offset.len() != n || transition.iter().any(|row| row.len() != n)
        {
            return Err(invalid(
                "transition, offset, and process-noise dimensions differ",
            ));
        }
        let nu = n as u128;
        let mut progress = WorkProgress::new(
            cx,
            WorkPlan::checked(
                0,
                0,
                0,
                0,
                8 * nu * nu + 8 * nu,
                0,
                32 * nu * nu + 32 * nu + 512,
            )?,
        )?;
        progress.checkpoint("linear dynamics admission")?;
        for value in transition.iter_mut().flatten().chain(offset.iter_mut()) {
            if !value.is_finite() {
                return Err(invalid("non-finite transition or offset"));
            }
            *value = canonicalize_zero(*value);
            progress.scalar("linear dynamics admission", 1)?;
        }
        if noise.mean().iter().any(|value| *value != 0.0) {
            return Err(invalid(
                "process noise must have zero mean; use the affine offset",
            ));
        }
        let identity = identity.into();
        validate_leaf_identity("linear dynamics", &identity)?;
        let mut hash = DomainHasher::new("org.frankensim.fs-assimilate.linear-dynamics.v1");
        hash.update(&(identity.len() as u64).to_le_bytes());
        hash.update(identity.as_bytes());
        hash.update(&(n as u64).to_le_bytes());
        progress.hash_bytes("linear dynamics identity", identity.len() as u128 + 16)?;
        for value in transition
            .iter()
            .flatten()
            .chain(&offset)
            .chain(&noise.upper)
            .chain(&noise.diag)
        {
            hash.update(&canonical_f64_bits(*value).to_le_bytes());
            progress.hash_bytes("linear dynamics identity", 8)?;
        }
        progress.checkpoint("linear dynamics publication")?;
        Ok(Self {
            transition,
            offset,
            noise,
            identity,
            content_id: hash.finalize(),
        })
    }

    /// Number of physical state variables evolved by this transition.
    #[must_use]
    pub fn dimension(&self) -> usize {
        self.offset.len()
    }

    /// Caller-declared dynamics identity, not authentication of a model.
    #[must_use]
    pub fn identity(&self) -> &str {
        &self.identity
    }

    /// Content hash binding the declaration and every transition/noise value.
    #[must_use]
    pub const fn content_id(&self) -> ContentHash {
        self.content_id
    }

    /// Predict a factor belief, retaining covariance as `U D U^T`.
    ///
    /// `O(n^3)` work, `O(n^2)` fallible flat scratch, one cumulative work/poll
    /// accountant. Cancellation or numeric refusal leaves the prior unchanged.
    /// No covariance-contraction claim applies: process noise can increase it.
    pub fn predict(&self, prior: &FactorBelief, cx: &Cx<'_>) -> Result<FactorBelief, AssimError> {
        if prior.dim() != self.dimension() {
            return Err(AssimError::DimMismatch {
                state: prior.dim(),
                operator: self.dimension(),
            });
        }
        self.predict_prefix(prior, cx)
    }

    // Evolve a physical prefix and retain a static, possibly correlated suffix.
    // With U=[U_xx U_xb; 0 U_bb], diag(F,I) U has cross block F U_xb.
    // Refactor only F U_xx D_x U_xx^T F^T+Q. Nuisance columns/pivots stay
    // unchanged, so their learned cross-covariance is never marginalized away.
    pub(crate) fn predict_prefix(
        &self,
        prior: &FactorBelief,
        cx: &Cx<'_>,
    ) -> Result<FactorBelief, AssimError> {
        let n = self.dimension();
        let d = prior.dim();
        if d < n {
            return Err(AssimError::DimMismatch {
                state: d,
                operator: n,
            });
        }
        validate_state_dimension(d)?;
        // Dimensions are already capped at 256; these products fit u128/usize.
        let mut progress = WorkProgress::new(
            cx,
            WorkPlan::checked(
                0,
                0,
                0,
                0,
                prediction_work(n, d),
                0,
                0,
            )?,
        )?;
        predict_prefix_factors(
            &self.transition, PredictionMean::Affine(&self.offset), &self.noise,
            prior, &mut progress,
        )
    }
}

fn prediction_work(n: usize, d: usize) -> u128 {
    let (n, d) = (n as u128, d as u128);
    32 * n * n * n + 16 * n * n * d + 32 * d * d + 128
}

#[derive(Clone, Copy)]
enum PredictionMean<'a> {
    Affine(&'a [f64]),
    // Use f(mean) directly: subtracting J*mean to form an affine offset and
    // adding it back can erase a perfectly representable nonlinear prediction.
    Evaluated(&'a [f64]),
}

fn predict_prefix_factors(
    transition: &[Vec<f64>],
    prediction_mean: PredictionMean<'_>,
    noise: &FactorBelief,
    prior: &FactorBelief,
    progress: &mut WorkProgress<'_, '_>,
) -> Result<FactorBelief, AssimError> {
    let n = transition.len();
    let d = prior.dim();
    let du = d as u128;
    progress.checkpoint("linear prediction admission")?;
    let mut mean = try_f64_vec(d, "linear prediction mean", progress)?;
    let mut upper = try_f64_vec(
        packed_upper_len(d)?,
        "linear prediction factor",
        progress,
    )?;
    let mut diag = try_f64_vec(d, "linear prediction pivots", progress)?;
    let mut roots = try_f64_vec(2 * n * n, "linear prediction roots", progress)?;
    let mut column = try_f64_vec(n, "linear prediction column", progress)?;
    let width = 2 * n;
    for row in 0..n {
        mean[row] = match prediction_mean {
            PredictionMean::Affine(offset) => {
                let (value, _) = compensated_dot(
                    &transition[row], &prior.mean[..n], "linear prediction mean",
                    "linear prediction mean", progress,
                )?;
                finite(value + offset[row], "linear prediction mean")?
            }
            PredictionMean::Evaluated(value) => finite(value[row], "nonlinear prediction mean")?,
        };
        progress.scalar("prediction mean", 1)?;
    }
    for j in 0..n {
        let scale = prior.diag[j].sqrt();
        for (k, value) in column.iter_mut().enumerate() {
            *value = if k > j || scale == 0.0 {
                0.0
            } else if k == j {
                scale
            } else {
                finite(
                    prior.upper[packed_upper_index(d, k, j)] * scale,
                    "prior square root",
                )?
            };
            progress.scalar("prior square root", 2)?;
        }
        let noise_scale = noise.diag[j].sqrt();
        for row in 0..n {
            let (value, _) = compensated_dot(
                &transition[row],
                &column,
                "propagated square root",
                "propagated square root",
                progress,
            )?;
            roots[row * width + j] = value;
            roots[row * width + n + j] = if row > j || noise_scale == 0.0 {
                0.0
            } else if row == j {
                noise_scale
            } else {
                finite(
                    noise.upper[packed_upper_index(n, row, j)] * noise_scale,
                    "process square root",
                )?
            };
            progress.scalar("process square root", 2)?;
        }
    }
    // Keep static nuisance factors exactly, and propagate their cross block.
    for j in n..d {
        mean[j] = prior.mean[j];
        diag[j] = prior.diag[j];
        for (k, value) in column.iter_mut().enumerate() {
            *value = prior.upper[packed_upper_index(d, k, j)];
        }
        for row in 0..j {
            upper[packed_upper_index(d, row, j)] = if row < n {
                compensated_dot(
                    &transition[row],
                    &column,
                    "forecast cross covariance",
                    "forecast cross covariance",
                    progress,
                )?
                .0
            } else {
                prior.upper[packed_upper_index(d, row, j)]
            };
            progress.scalar("forecast nuisance retention", 1)?;
        }
    }
    let mut active_columns = width;
    for row in (0..n).rev() {
        progress.checkpoint("linear prediction orthogonalization")?;
        let pivot_column = active_columns - 1;
        for k in 0..pivot_column {
            let a = roots[row * width + pivot_column];
            let b = roots[row * width + k];
            if b == 0.0 {
                continue;
            }
            let radius = finite(a.hypot(b), "forecast root norm")?;
            let cosine = a / radius;
            let sine = b / radius;
            for i in 0..row {
                let old_pivot = roots[i * width + pivot_column];
                let old_other = roots[i * width + k];
                roots[i * width + pivot_column] = finite(
                    cosine.mul_add(old_pivot, sine * old_other),
                    "forecast rotation",
                )?;
                roots[i * width + k] = finite(
                    cosine.mul_add(old_other, -sine * old_pivot),
                    "forecast rotation",
                )?;
                progress.scalar("forecast rotation", 8)?;
            }
            roots[row * width + pivot_column] = radius;
            roots[row * width + k] = 0.0;
            progress.scalar("forecast rotation", 8)?;
        }
        let pivot = roots[row * width + pivot_column];
        if pivot == 0.0 {
            // A zero row consumes no root column: earlier rows may still
            // carry variance in the candidate pivot column.
            continue;
        }
        diag[row] = finite(pivot * pivot, "forecast variance")?;
        if diag[row] == 0.0 {
            return Err(AssimError::CovarianceCertificationUnresolved);
        }
        for i in 0..row {
            upper[packed_upper_index(d, i, row)] = finite(
                roots[i * width + pivot_column] / pivot,
                "forecast unit-upper factor",
            )?;
            progress.scalar("forecast factor extraction", 1)?;
        }
        active_columns -= 1;
    }
    let prediction = FactorBelief { mean, upper, diag };
    // Finite UD entries alone do not imply a representable covariance:
    // a large off-diagonal unit factor can overflow a component variance.
    // Every published component must remain usable by the factor API.
    for row in 0..d {
        finite(
            prediction.variance(row).expect("bounded factor row"),
            "forecast marginal variance",
        )?;
        progress.scalar("forecast marginal variance", 3 * du)?;
    }
    progress.checkpoint("linear prediction publication")?;
    Ok(prediction)
}

pub(crate) const fn invalid(reason: &'static str) -> AssimError {
    AssimError::InvalidForecast { reason }
}

fn finite(value: f64, stage: &'static str) -> Result<f64, AssimError> {
    if value.is_finite() {
        Ok(canonicalize_zero(value))
    } else {
        Err(AssimError::NonFiniteComputation { stage })
    }
}
