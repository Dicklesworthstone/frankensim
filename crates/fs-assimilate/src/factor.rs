//! Stable factor-form scalar assimilation substrate (bead sj31i.37).
//!
//! The dense Joseph path in `lib.rs` is `O(n^3)` per scalar observation and
//! can only revalidate the posterior as PSD; it cannot express the Loewner
//! contraction claim `P_prior - P_posterior >= 0` with executable evidence.
//! This module carries the belief in `U D U^T` factor form and applies the
//! Thornton-Bierman scalar UD measurement update, which is `O(n^2)`,
//! square-root-free, and PSD by construction. Every update returns a
//! [`ContractionReceipt`] with an executable `Certified` / `Refuted` /
//! `Unresolved` verdict:
//!
//! - the alpha chain `a_j = a_{j-1} + f_j g_j` (`a_0 = R`,
//!   `a_n = h P h^T + R`) is evaluated simultaneously in binary64 and in
//!   outward-rounded [`Interval`] arithmetic;
//! - each factor pivot contracts as `d'_j = d_j * a_{j-1} / a_j <= d_j` in
//!   exact arithmetic; the interval ratios make that claim independently
//!   checkable instead of prose;
//! - the exact scalar identity for the measurement direction,
//!   `h P' h^T = (h P h^T) * R / a_n`, is checked against an interval
//!   enclosure computed from the same chain;
//! - a verdict is `Certified` only when every enclosure is decisive and the
//!   point arithmetic agrees, `Refuted` when a point check or an enclosure
//!   is decisively violated, and `Unresolved` otherwise. Unresolved never
//!   advertises a contraction.
//!
//! Derivation note (the U column update): with `f = U^T h`, `g = D f`, the
//! posterior factors satisfy `U' D' U'^T = U [D - g g^T / a_n] U^T`. The
//! inner rank-one downdate factors as `U~ D' U~^T` with
//! `u~_{ij} = -g_i f_j / a_{j-1}` (verified symbolically for n = 2, 3
//! against the direct expansion), so `U' = U U~` gives
//! `u'_{ij} = u_{ij} - (f_j / a_{j-1}) * sum_{k=i}^{j-1} u_{ik} g_k`.
//! Processing columns ascending with one auxiliary vector `b` (each
//! element read before it is overwritten) reproduces the original-factor
//! sums in place.
//!
//! The dense Joseph updater remains as the independent oracle used by
//! [`verify_factor_assimilation`] and the conformance tests.

use fs_ivl::Interval;

use crate::{
    AssimError, Belief, Cx, Observation, WorkPlan, WorkProgress, assimilate, checked_work_add,
    checked_work_mul,
};

/// Contraction authority attached to one factor-form scalar update.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ContractionState {
    /// Every enclosure is decisive: pivots contract, the innovation
    /// enclosure is wholly positive, and the measurement-direction identity
    /// holds. The contraction claim carries executable evidence.
    Certified,
    /// A point check or an enclosure is decisively violated: the computed
    /// update expanded a direction or broke the exact identity. The claim
    /// is affirmatively rejected, not merely unproven.
    Refuted,
    /// Enclosures straddle or are indecisive: floating arithmetic cannot
    /// resolve the contraction claim. No contraction is advertised.
    Unresolved,
}

impl ContractionState {
    /// Stable wire name for receipts and structured logs.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::Certified => "certified",
            Self::Refuted => "refuted",
            Self::Unresolved => "unresolved",
        }
    }
}

/// Checked diagnostic for the scalar misfit monotonicity law: in exact
/// arithmetic the weighted scalar misfit cannot increase through the
/// optimal update (`after = before * (R / a_n)^2`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MisfitVerdict {
    /// Computed misfit did not increase beyond the computed arithmetic
    /// error envelope.
    NonIncreasing,
    /// Computed misfit increased beyond the envelope: the diagnostic
    /// rejects the monotonicity law for this update.
    Violated,
    /// The error envelope cannot decide the comparison.
    Indeterminate,
}

impl MisfitVerdict {
    /// Stable wire name for receipts and structured logs.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::NonIncreasing => "non-increasing",
            Self::Violated => "violated",
            Self::Indeterminate => "indeterminate",
        }
    }
}

/// Executable evidence for one scalar covariance contraction.
///
/// The receipt is data, not a comment: the measurement variances, the
/// identity enclosure, and the verdict are recomputed independently by
/// [`verify_factor_assimilation`], and the whole record is bound into
/// `identity`.
#[derive(Debug, Clone, PartialEq)]
pub struct ContractionReceipt {
    state: ContractionState,
    innovation_variance: f64,
    innovation_enclosure: (f64, f64),
    max_pivot_ratio: f64,
    first_expanding_pivot: Option<usize>,
    measurement_variance_prior: f64,
    measurement_variance_posterior: f64,
    measurement_identity_enclosure: (f64, f64),
    misfit_before: f64,
    misfit_after: f64,
    misfit_verdict: MisfitVerdict,
    identity: String,
}

impl ContractionReceipt {
    /// The executable verdict.
    #[must_use]
    pub const fn state(&self) -> ContractionState {
        self.state
    }

    /// Computed innovation variance `h P h^T + R`.
    #[must_use]
    pub const fn innovation_variance(&self) -> f64 {
        self.innovation_variance
    }

    /// Outward-rounded enclosure of the innovation variance `(lo, hi)`.
    #[must_use]
    pub const fn innovation_enclosure(&self) -> (f64, f64) {
        self.innovation_enclosure
    }

    /// Largest computed pivot ratio `d'_j / d_j` (zero pivots excluded).
    #[must_use]
    pub const fn max_pivot_ratio(&self) -> f64 {
        self.max_pivot_ratio
    }

    /// First pivot whose computed diagonal definitely expanded or turned
    /// negative, if any.
    #[must_use]
    pub const fn first_expanding_pivot(&self) -> Option<usize> {
        self.first_expanding_pivot
    }

    /// Computed prior measurement-direction variance `h P h^T`.
    #[must_use]
    pub const fn measurement_variance_prior(&self) -> f64 {
        self.measurement_variance_prior
    }

    /// Computed posterior measurement-direction variance `h P' h^T`.
    #[must_use]
    pub const fn measurement_variance_posterior(&self) -> f64 {
        self.measurement_variance_posterior
    }

    /// Interval enclosure of the exact identity value
    /// `(h P h^T) * R / a_n` for the posterior measurement variance.
    #[must_use]
    pub const fn measurement_identity_enclosure(&self) -> (f64, f64) {
        self.measurement_identity_enclosure
    }

    /// Weighted scalar misfit before the update.
    #[must_use]
    pub const fn misfit_before(&self) -> f64 {
        self.misfit_before
    }

    /// Weighted scalar misfit after the update.
    #[must_use]
    pub const fn misfit_after(&self) -> f64 {
        self.misfit_after
    }

    /// Checked misfit monotonicity diagnostic.
    #[must_use]
    pub const fn misfit_verdict(&self) -> MisfitVerdict {
        self.misfit_verdict
    }

    /// Domain-separated receipt identity
    /// `scalar-contraction:v1:<64 lowercase hex>`.
    #[must_use]
    pub fn identity(&self) -> &str {
        &self.identity
    }
}

/// A Gaussian belief carried in `U D U^T` factor form.
///
/// `U` is unit upper triangular stored packed by rows (strict triangle),
/// `D` is an exactly non-negative diagonal. The factor form keeps every
/// posterior PSD by construction and each scalar update `O(n^2)`.
#[derive(Debug, Clone, PartialEq)]
pub struct FactorBelief {
    mean: Vec<f64>,
    upper: Vec<f64>,
    diag: Vec<f64>,
}

fn packed_upper_len(n: usize) -> Result<usize, AssimError> {
    let units = checked_work_mul(n as u128, (n as u128).saturating_sub(1), "factor shape")?;
    usize::try_from(units / 2).map_err(|_| AssimError::WorkPlanOverflow {
        phase: "factor shape",
    })
}

fn packed_upper_index(n: usize, row: usize, column: usize) -> usize {
    debug_assert!(row < column && column < n);
    // row*(2n - row - 1)/2 + (column - row - 1); inputs are
    // dimension-checked by construction and n <= MAX_DENSE_STATE_DIM.
    row * (2 * n - row - 1) / 2 + (column - row - 1)
}

fn try_f64_vec(
    len: usize,
    stage: &'static str,
    progress: &mut WorkProgress<'_, '_>,
) -> Result<Vec<f64>, AssimError> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(len)
        .map_err(|_| AssimError::AllocationRefused {
            stage,
            requested_bytes: (len as u128).saturating_mul(8),
        })?;
    values.resize(len, 0.0);
    progress.scalar(stage, len as u128)?;
    Ok(values)
}

impl FactorBelief {
    /// State dimension.
    #[must_use]
    pub fn dim(&self) -> usize {
        self.mean.len()
    }

    /// Read-only mean.
    #[must_use]
    pub fn mean(&self) -> &[f64] {
        &self.mean
    }

    /// Factor diagonal entry `d_j` (exactly non-negative).
    #[must_use]
    pub fn diag(&self, index: usize) -> Option<f64> {
        self.diag.get(index).copied()
    }

    /// Strict upper-triangle factor entry `u_{i,j}` for `i < j`.
    #[must_use]
    pub fn upper_entry(&self, row: usize, column: usize) -> Option<f64> {
        if row >= column || column >= self.dim() {
            return None;
        }
        Some(self.upper[packed_upper_index(self.dim(), row, column)])
    }

    /// Component variance `P_{i,i}` reconstructed from the factor.
    #[must_use]
    pub fn variance(&self, component: usize) -> Option<f64> {
        let n = self.dim();
        if component >= n {
            return None;
        }
        let mut total = self.diag[component];
        for column in (component + 1)..n {
            let u = self.upper[packed_upper_index(n, component, column)];
            total = (u * u * self.diag[column]) + total;
        }
        Some(total)
    }

    /// A diagonal factor belief (`U = I`); variances must be finite and
    /// non-negative, and signed-zero is canonicalized to `+0.0`.
    ///
    /// # Errors
    /// Returns [`AssimError`] for mismatched lengths, an empty state, a
    /// non-finite mean, or a negative or non-finite variance.
    pub fn diagonal(mean: Vec<f64>, variances: Vec<f64>, cx: &Cx<'_>) -> Result<Self, AssimError> {
        if mean.is_empty() {
            return Err(AssimError::EmptyBelief);
        }
        if mean.len() != variances.len() {
            return Err(AssimError::DiagonalDimensionMismatch {
                means: mean.len(),
                variances: variances.len(),
            });
        }
        let n = mean.len();
        let plan = factor_construction_work_plan(n)?;
        let mut progress = WorkProgress::new(cx, plan)?;
        progress.checkpoint("initial")?;
        for (index, entry) in mean.iter().enumerate() {
            if !entry.is_finite() {
                return Err(AssimError::NonFiniteMean { index });
            }
            progress.scalar("factor construction", 1)?;
        }
        let mut diag = try_f64_vec(n, "factor diagonal", &mut progress)?;
        for (index, (target, variance)) in diag.iter_mut().zip(&variances).enumerate() {
            if !variance.is_finite() || *variance < 0.0 {
                return Err(AssimError::NegativeVariance { index });
            }
            *target = if *variance == 0.0 { 0.0 } else { *variance };
            progress.scalar("factor construction", 1)?;
        }
        let upper = try_f64_vec(packed_upper_len(n)?, "factor upper", &mut progress)?;
        Ok(Self { mean, upper, diag })
    }

    /// Factor a validated dense belief covariance into `U D U^T` form with
    /// exact zero-pivot handling: a validated PSD input yields exact zero
    /// diagonals for singular directions, whose factor columns stay exact
    /// zeros. A negative pivot from rounding pathology refuses rather than
    /// clamping.
    ///
    /// # Errors
    /// Returns [`AssimError`] for a non-finite intermediate or a pivot that
    /// proves the admitted covariance is not representable as a factor pair.
    pub fn from_belief(belief: &Belief, cx: &Cx<'_>) -> Result<Self, AssimError> {
        let n = belief.dim();
        let plan = factor_construction_work_plan(n)?;
        let mut progress = WorkProgress::new(cx, plan)?;
        progress.checkpoint("initial")?;
        let mean = belief.mean().to_vec();
        let mut upper = try_f64_vec(packed_upper_len(n)?, "factor upper", &mut progress)?;
        let mut diag = try_f64_vec(n, "factor diagonal", &mut progress)?;
        let cov = belief.covariance();
        // Unit-UPPER U D U^T decomposition: entry (i, j) with i <= j is
        // sum_{k >= j} u_{i,k} d_k u_{j,k}, so pivots resolve BACKWARD
        // (column n-1 down to 0). A forward pass belongs to the unit-LOWER
        // L D L^T form and silently misfactors non-diagonal inputs.
        for j in (0..n).rev() {
            let mut pivot = cov[j][j];
            for (k, diag_k) in diag.iter().enumerate().skip(j + 1) {
                let u_jk = upper[packed_upper_index(n, j, k)];
                pivot -= u_jk * u_jk * *diag_k;
                progress.scalar("factor construction", 2)?;
            }
            if !pivot.is_finite() {
                return Err(AssimError::NonFiniteComputation {
                    stage: "factor pivot",
                });
            }
            if pivot < 0.0 {
                return Err(AssimError::CovarianceCertificationUnresolved);
            }
            diag[j] = if pivot == 0.0 { 0.0 } else { pivot };
            if diag[j] == 0.0 {
                continue;
            }
            for i in 0..j {
                let mut entry = cov[i][j];
                for (k, diag_k) in diag.iter().enumerate().skip(j + 1) {
                    let u_ik = upper[packed_upper_index(n, i, k)];
                    let u_jk = upper[packed_upper_index(n, j, k)];
                    entry -= u_ik * *diag_k * u_jk;
                    progress.scalar("factor construction", 2)?;
                }
                let scaled = entry / pivot;
                progress.scalar("factor construction", 1)?;
                if !scaled.is_finite() {
                    return Err(AssimError::NonFiniteComputation {
                        stage: "factor upper entry",
                    });
                }
                let index = packed_upper_index(n, i, j);
                upper[index] = scaled;
            }
        }
        Ok(Self { mean, upper, diag })
    }

    /// Reconstruct the dense covariance `U D U^T` (independent-oracle lane;
    /// `O(n^3)`, not part of the production update path).
    ///
    /// Entry `(i, j)` with `i <= j` is `sum_{k >= j} u_{i,k} d_k u_{j,k}`
    /// because `U` is upper triangular; the result is mirrored exactly.
    #[must_use]
    pub fn to_dense_covariance(&self) -> Vec<Vec<f64>> {
        let n = self.dim();
        let mut cov = vec![vec![0.0; n]; n];
        for i in 0..n {
            for j in i..n {
                let mut entry = 0.0;
                for (k, diag_k) in diag_iter_from(&self.diag, j) {
                    let u_ik = if k == i {
                        1.0
                    } else {
                        self.upper[packed_upper_index(n, i, k)]
                    };
                    let u_jk = if k == j {
                        1.0
                    } else {
                        self.upper[packed_upper_index(n, j, k)]
                    };
                    entry += u_ik * *diag_k * u_jk;
                }
                cov[i][j] = entry;
                cov[j][i] = entry;
            }
        }
        cov
    }
}

fn diag_iter_from(diag: &[f64], from: usize) -> impl Iterator<Item = (usize, &f64)> {
    diag.iter().enumerate().skip(from)
}

/// Work totals for the factor substrate, exposed so a parent workflow can
/// seal admission before constructing inputs (checked arithmetic).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ScalarFactorWorkEstimate {
    /// Factor construction from a dense belief (`O(n^3)` one-time).
    pub construction: u128,
    /// One stable scalar update with its contraction receipt (`O(n^2)`).
    pub update_with_receipt: u128,
    /// Independent checker lane (`O(n^3)` oracle comparison).
    pub independent_check: u128,
}

/// Checked work envelopes for the factor substrate at dimension `dim`.
///
/// # Errors
/// Returns [`AssimError::WorkPlanOverflow`] when the shape cannot be
/// represented.
pub fn scalar_factor_work_estimate(dim: usize) -> Result<ScalarFactorWorkEstimate, AssimError> {
    let n = dim as u128;
    let n2 = checked_work_mul(n, n, "factor estimate")?;
    let n3 = checked_work_mul(n2, n, "factor estimate")?;
    let construction = checked_work_add(
        checked_work_mul(2, n3, "factor estimate")?,
        checked_work_mul(8, n2, "factor estimate")?,
        "factor estimate",
    )?;
    let update_with_receipt = checked_work_add(
        checked_work_mul(24, n2, "factor estimate")?,
        checked_work_mul(64, n, "factor estimate")?,
        "factor estimate",
    )?;
    let independent_check = checked_work_add(
        checked_work_mul(4, n3, "factor estimate")?,
        checked_work_mul(24, n2, "factor estimate")?,
        "factor estimate",
    )?;
    Ok(ScalarFactorWorkEstimate {
        construction,
        update_with_receipt,
        independent_check,
    })
}

fn factor_construction_work_plan(dim: usize) -> Result<WorkPlan, AssimError> {
    let estimate = scalar_factor_work_estimate(dim)?;
    WorkPlan::checked(0, 0, 0, 0, estimate.construction, 0, 64)
}

fn factor_update_work_plan(dim: usize) -> Result<WorkPlan, AssimError> {
    let estimate = scalar_factor_work_estimate(dim)?;
    WorkPlan::checked(0, 0, 0, 8, estimate.update_with_receipt, 0, 256)
}

fn factor_checker_work_plan(dim: usize) -> Result<WorkPlan, AssimError> {
    let estimate = scalar_factor_work_estimate(dim)?;
    WorkPlan::checked(
        0,
        0,
        0,
        8,
        estimate.independent_check,
        estimate.construction,
        256,
    )
}

/// Compensated (Neumaier) dot product; returns the sum and a rounding
/// envelope `2 * eps * sum|terms|` for certificate comparisons.
fn compensated_dot(
    a: &[f64],
    b: &[f64],
    stage: &'static str,
    phase: &'static str,
    progress: &mut WorkProgress<'_, '_>,
) -> Result<(f64, f64), AssimError> {
    debug_assert_eq!(a.len(), b.len());
    let mut sum = 0.0_f64;
    let mut err = 0.0_f64;
    let mut scale = 0.0_f64;
    for (left, right) in a.iter().zip(b) {
        let product = *left * *right;
        let next = sum + product;
        if sum.abs() >= product.abs() {
            err += (sum - next) + product;
        } else {
            err += (product - next) + sum;
        }
        sum = next;
        scale += product.abs();
        progress.scalar(phase, 3)?;
        if !sum.is_finite() {
            return Err(AssimError::NonFiniteComputation { stage });
        }
    }
    let total = sum + err;
    if !total.is_finite() {
        return Err(AssimError::NonFiniteComputation { stage });
    }
    Ok((total, 2.0_f64 * f64::EPSILON * scale))
}

/// The result of one stable scalar update: the posterior factor belief and
/// its executable contraction receipt.
#[derive(Debug, Clone, PartialEq)]
pub struct FactorAssimilation {
    belief: FactorBelief,
    receipt: ContractionReceipt,
}

impl FactorAssimilation {
    /// The factor-form posterior.
    #[must_use]
    pub const fn belief(&self) -> &FactorBelief {
        &self.belief
    }

    /// The executable contraction receipt.
    #[must_use]
    pub const fn receipt(&self) -> &ContractionReceipt {
        &self.receipt
    }
}

const RECEIPT_DOMAIN: &str = "org.frankensim.fs-assimilate.scalar-contraction.v1";
const RECEIPT_ID_PREFIX: &str = "scalar-contraction:v1:";

fn receipt_identity(
    prior: &FactorBelief,
    obs: &Observation,
    innovation_variance: f64,
    state: ContractionState,
) -> String {
    let mut hasher = fs_blake3::DomainHasher::new(RECEIPT_DOMAIN);
    hasher.update(&(prior.mean.len() as u64).to_le_bytes());
    for entry in &prior.mean {
        hasher.update(&entry.to_le_bytes());
    }
    hasher.update(&(prior.upper.len() as u64).to_le_bytes());
    for entry in &prior.upper {
        hasher.update(&entry.to_le_bytes());
    }
    hasher.update(&(prior.diag.len() as u64).to_le_bytes());
    for entry in &prior.diag {
        hasher.update(&entry.to_le_bytes());
    }
    hasher.update(&(obs.operator().len() as u64).to_le_bytes());
    for entry in obs.operator() {
        hasher.update(&entry.to_le_bytes());
    }
    hasher.update(&obs.value().to_le_bytes());
    hasher.update(&obs.noise_var().to_le_bytes());
    hasher.update(&(obs.instrument().len() as u64).to_le_bytes());
    hasher.update(obs.instrument().as_bytes());
    hasher.update(&innovation_variance.to_le_bytes());
    hasher.update(state.name().as_bytes());
    format!("{RECEIPT_ID_PREFIX}{}", hasher.finalize())
}

/// One stable scalar factor update with its executable contraction receipt.
///
/// The covariance path is the Thornton-Bierman UD measurement update:
/// `f = U^T h`, `g = D f`, `a_j = a_{j-1} + f_j g_j` with `a_0 = R`,
/// `d'_j = d_j a_{j-1}/a_j`, and an in-place column update of `U` with one
/// auxiliary vector — `O(n^2)` total with `4n` flat scratch. The dense
/// Joseph updater is not consulted; the oracle comparison lives in
/// [`verify_factor_assimilation`].
///
/// # Errors
/// Returns [`AssimError`] for a dimension mismatch, a non-finite or
/// non-positive innovation variance, non-finite arithmetic, an exhausted
/// budget, or observed cancellation. A numerical contraction verdict is
/// never an error: it is reported as [`ContractionState`] data.
pub fn assimilate_scalar(
    prior: &FactorBelief,
    obs: &Observation,
    cx: &Cx<'_>,
) -> Result<FactorAssimilation, AssimError> {
    let n = prior.dim();
    if obs.operator().len() != n {
        return Err(AssimError::DimMismatch {
            state: n,
            operator: obs.operator().len(),
        });
    }
    let plan = factor_update_work_plan(n)?;
    let mut progress = WorkProgress::new(cx, plan)?;
    progress.checkpoint("initial")?;

    let h = obs.operator();
    let noise = obs.noise_var();

    // f = U^T h (f_i = sum_{j<=i} u_{j,i} h_j, u_{i,i} = 1); g = D f.
    let mut f = try_f64_vec(n, "factor sensitivity", &mut progress)?;
    for i in 0..n {
        let mut total = h[i];
        for j in 0..i {
            total = prior.upper[packed_upper_index(n, j, i)].mul_add(h[j], total);
            progress.scalar("factor-update", 1)?;
        }
        if !total.is_finite() {
            return Err(AssimError::NonFiniteComputation {
                stage: "factor sensitivity",
            });
        }
        f[i] = total;
    }
    let mut g = try_f64_vec(n, "factor scaled sensitivity", &mut progress)?;
    for (slot, (f_i, d_i)) in g.iter_mut().zip(f.iter().zip(&prior.diag)) {
        *slot = f_i * d_i;
        progress.scalar("factor-update", 1)?;
        if !slot.is_finite() {
            return Err(AssimError::NonFiniteComputation {
                stage: "factor scaled sensitivity",
            });
        }
    }

    // Alpha chain in binary64 and outward-rounded interval arithmetic.
    // a_j is nondecreasing in exact arithmetic because f_j g_j = d_j f_j^2
    // is non-negative; a non-finite or non-positive point value refuses,
    // while the interval chain drives the receipt verdict only.
    let mut alpha = try_f64_vec(n + 1, "factor alpha chain", &mut progress)?;
    let mut alpha_ivl = vec![Interval::point(noise); n + 1];
    alpha[0] = noise;
    for j in 0..n {
        let next = f[j].mul_add(g[j], alpha[j]);
        progress.scalar("factor-update", 2)?;
        if !next.is_finite() {
            return Err(AssimError::NonFiniteComputation {
                stage: "factor alpha chain",
            });
        }
        alpha[j + 1] = next;
        alpha_ivl[j + 1] = alpha_ivl[j] + Interval::point(f[j]) * Interval::point(g[j]);
    }
    let innovation_variance = alpha[n];
    if innovation_variance <= 0.0 {
        return Err(AssimError::SingularInnovation);
    }
    let chain = alpha_ivl[n];
    let chain_decisive = chain.lo().is_finite() && chain.hi().is_finite();

    // Kalman gain FIRST: K = U g / a_n with the PRIOR factor,
    // (U g)_i = sum_{k>=i} u_{i,k} g_k. The downdate below rewrites `upper`
    // in place, so any gain computed after it would silently use U'.
    let mut gain = try_f64_vec(n, "factor gain", &mut progress)?;
    for i in 0..n {
        let mut total = g[i];
        for (k, g_k) in g.iter().enumerate().skip(i + 1) {
            total = prior.upper[packed_upper_index(n, i, k)].mul_add(*g_k, total);
            progress.scalar("factor-update", 1)?;
        }
        let gain_entry = total / innovation_variance;
        progress.scalar("factor-update", 1)?;
        if !gain_entry.is_finite() {
            return Err(AssimError::NonFiniteComputation {
                stage: "factor gain",
            });
        }
        gain[i] = gain_entry;
    }

    // Mean update with compensated prediction.
    let (predicted, _prediction_envelope) = compensated_dot(
        h,
        &prior.mean,
        "observation prediction",
        "factor-update",
        &mut progress,
    )?;
    let innovation = obs.value() - predicted;
    progress.scalar("factor-update", 1)?;
    if !innovation.is_finite() {
        return Err(AssimError::NonFiniteComputation {
            stage: "observation innovation",
        });
    }
    let mut mean = prior.mean.clone();
    for (mean_i, gain_i) in mean.iter_mut().zip(&gain) {
        *mean_i = gain_i.mul_add(innovation, *mean_i);
        progress.scalar("factor-update", 2)?;
        if !mean_i.is_finite() {
            return Err(AssimError::NonFiniteComputation {
                stage: "factor posterior mean",
            });
        }
    }

    // Factor downdate with per-pivot interval ratios.
    let mut diag = prior.diag.clone();
    let mut upper = prior.upper.clone();
    let mut b = try_f64_vec(n, "factor update scratch", &mut progress)?;
    b.copy_from_slice(&g);
    let mut max_pivot_ratio = 0.0_f64;
    let mut first_expanding_pivot = None;
    let mut any_refuted = false;
    for j in 0..n {
        let ratio = alpha[j] / alpha[j + 1];
        progress.scalar("factor-update", 1)?;
        if !ratio.is_finite() {
            return Err(AssimError::NonFiniteComputation {
                stage: "factor pivot ratio",
            });
        }
        if diag[j] != 0.0 {
            let updated = diag[j] * ratio;
            progress.scalar("factor-update", 1)?;
            if !updated.is_finite() {
                return Err(AssimError::NonFiniteComputation {
                    stage: "factor pivot update",
                });
            }
            if updated < 0.0 || updated > diag[j] {
                if first_expanding_pivot.is_none() {
                    first_expanding_pivot = Some(j);
                }
            }
            let pivot_ratio = updated / diag[j];
            if pivot_ratio > max_pivot_ratio {
                max_pivot_ratio = pivot_ratio;
            }
            diag[j] = updated;
        }
        let ratio_ivl = alpha_ivl[j] / alpha_ivl[j + 1];
        if ratio_ivl.lo() > 1.0 {
            // The computed chain itself decisively expanded this pivot:
            // reachable only when the enclosures have gone haywire, since
            // f_j g_j = d_j f_j^2 is non-negative by construction.
            any_refuted = true;
        }
        // Column update of U (in place; b advances with the original entry).
        if j >= 1 {
            let lambda = -f[j] / alpha[j];
            progress.scalar("factor-update", 1)?;
            if !lambda.is_finite() {
                return Err(AssimError::NonFiniteComputation {
                    stage: "factor column scale",
                });
            }
            for i in 0..j {
                let index = packed_upper_index(n, i, j);
                let original = upper[index];
                let updated_entry = lambda.mul_add(b[i], original);
                b[i] = original.mul_add(g[j], b[i]);
                progress.scalar("factor-update", 4)?;
                if !updated_entry.is_finite() || !b[i].is_finite() {
                    return Err(AssimError::NonFiniteComputation {
                        stage: "factor column update",
                    });
                }
                upper[index] = updated_entry;
            }
        }
    }

    // Measurement-direction identity: q = h P h^T = f . g and
    // q' = h P' h^T via the updated factor; exact law q' = q R / a_n.
    let (q_prior, q_prior_envelope) = compensated_dot(
        &f,
        &g,
        "prior measurement variance",
        "factor-receipt",
        &mut progress,
    )?;
    let mut f_post = try_f64_vec(n, "posterior sensitivity", &mut progress)?;
    for i in 0..n {
        let mut total = h[i];
        for j in 0..i {
            total = upper[packed_upper_index(n, j, i)].mul_add(h[j], total);
            progress.scalar("factor-receipt", 1)?;
        }
        if !total.is_finite() {
            return Err(AssimError::NonFiniteComputation {
                stage: "posterior sensitivity",
            });
        }
        f_post[i] = total;
    }
    let mut scaled_post = try_f64_vec(n, "posterior scaled sensitivity", &mut progress)?;
    for (slot, (f_i, d_i)) in scaled_post.iter_mut().zip(f_post.iter().zip(&diag)) {
        *slot = f_i * d_i;
        progress.scalar("factor-receipt", 1)?;
    }
    let (q_post, q_post_envelope) = compensated_dot(
        &f_post,
        &scaled_post,
        "posterior measurement variance",
        "factor-receipt",
        &mut progress,
    )?;

    let identity_enclosure = (chain - Interval::point(noise)) * Interval::point(noise) / chain;
    let identity_decisive =
        identity_enclosure.lo().is_finite() && identity_enclosure.hi().is_finite();
    let identity_holds = identity_decisive
        && identity_enclosure.lo() - q_post_envelope <= q_post
        && q_post <= identity_enclosure.hi() + q_post_envelope;
    let identity_violated = identity_decisive
        && (q_post + q_post_envelope < identity_enclosure.lo()
            || q_post - q_post_envelope > identity_enclosure.hi());

    // Misfit monotonicity diagnostic (checked, not prose).
    let misfit_before = innovation * innovation / noise;
    let (posterior_prediction, _post_pred_envelope) = compensated_dot(
        h,
        &mean,
        "posterior prediction",
        "factor-receipt",
        &mut progress,
    )?;
    let posterior_residual = obs.value() - posterior_prediction;
    let misfit_after = posterior_residual * posterior_residual / noise;
    progress.scalar("factor-receipt", 4)?;
    if !misfit_before.is_finite() || !misfit_after.is_finite() {
        return Err(AssimError::NonFiniteComputation {
            stage: "misfit diagnostic",
        });
    }
    let misfit_envelope = 2.0_f64
        * f64::EPSILON
        * (misfit_before.abs() + misfit_after.abs() + q_prior_envelope + q_post_envelope);
    let misfit_verdict = if misfit_after <= misfit_before + misfit_envelope {
        MisfitVerdict::NonIncreasing
    } else if misfit_after > misfit_before + 16.0 * misfit_envelope {
        MisfitVerdict::Violated
    } else {
        MisfitVerdict::Indeterminate
    };

    // Verdict algebra. The pivot contraction d'_j = d_j a_{j-1}/a_j <= d_j
    // is an exact-arithmetic theorem enforced by the point checks above;
    // the interval chain certifies that the floating computation stayed
    // consistent with that semantics (decisive, wholly positive), and the
    // measurement-direction identity enclosure must contain the computed
    // posterior value. A decisive breach refutes, an indecisive enclosure
    // leaves the claim unresolved, and a fully decisive clean bill
    // certifies. Outward-rounded interval ratios are never used to certify:
    // a 1-ULP outward nudge would make even exact-zero steps straddle 1.
    let state = if first_expanding_pivot.is_some() || any_refuted || identity_violated {
        ContractionState::Refuted
    } else if !chain_decisive || !identity_decisive || !identity_holds {
        ContractionState::Unresolved
    } else {
        ContractionState::Certified
    };

    let receipt = ContractionReceipt {
        state,
        innovation_variance,
        innovation_enclosure: (chain.lo(), chain.hi()),
        max_pivot_ratio,
        first_expanding_pivot,
        measurement_variance_prior: q_prior,
        measurement_variance_posterior: q_post,
        measurement_identity_enclosure: (identity_enclosure.lo(), identity_enclosure.hi()),
        misfit_before,
        misfit_after,
        misfit_verdict,
        identity: receipt_identity(prior, obs, innovation_variance, state),
    };
    progress.checkpoint("finalize")?;
    Ok(FactorAssimilation {
        belief: FactorBelief { mean, upper, diag },
        receipt,
    })
}

/// Verdict of the independent checker lane.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CheckVerdict {
    /// The factor posterior agrees with the dense Joseph oracle within the
    /// computed scale tolerance and the receipt re-verifies independently.
    Verified,
    /// A decisive disagreement between the factor result and the
    /// independent oracle, or an independently recomputed value that
    /// contradicts the reported receipt.
    Discrepancy,
}

impl CheckVerdict {
    /// Stable wire name for receipts and structured logs.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::Verified => "verified",
            Self::Discrepancy => "discrepancy",
        }
    }
}

/// Independent check of one factor assimilation against the dense Joseph
/// oracle and an independently recomputed contraction identity.
#[derive(Debug, Clone, PartialEq)]
pub struct IndependentCheck {
    verdict: CheckVerdict,
    max_abs_diff: f64,
    tolerance: f64,
    recomputed_innovation_variance: f64,
    identity_consistent: bool,
}

impl IndependentCheck {
    /// The checker verdict.
    #[must_use]
    pub const fn verdict(&self) -> CheckVerdict {
        self.verdict
    }

    /// Largest absolute entry difference between the reconstructed factor
    /// covariance and the dense Joseph oracle covariance.
    #[must_use]
    pub const fn max_abs_diff(&self) -> f64 {
        self.max_abs_diff
    }

    /// The computed scale tolerance used for the comparison.
    #[must_use]
    pub const fn tolerance(&self) -> f64 {
        self.tolerance
    }

    /// Innovation variance recomputed independently from the dense prior.
    #[must_use]
    pub const fn recomputed_innovation_variance(&self) -> f64 {
        self.recomputed_innovation_variance
    }

    /// True when the independently recomputed identity contains the
    /// reported posterior measurement variance.
    #[must_use]
    pub const fn identity_consistent(&self) -> bool {
        self.identity_consistent
    }
}

/// Independently verify one factor assimilation: reconstruct the factor
/// covariance densely, compare it against the dense Joseph oracle, and
/// recompute the innovation variance and the measurement-direction identity
/// from the dense inputs.
///
/// This lane is deliberately `O(n^3)`: it is the oracle, not the production
/// path. A `Discrepancy` is data about the update, not an error.
///
/// # Errors
/// Returns [`AssimError`] when the oracle update itself fails or budgets or
/// cancellation refuse.
pub fn verify_factor_assimilation(
    prior: &Belief,
    obs: &Observation,
    result: &FactorAssimilation,
    cx: &Cx<'_>,
) -> Result<IndependentCheck, AssimError> {
    let n = prior.dim();
    let plan = factor_checker_work_plan(n)?;
    let mut progress = WorkProgress::new(cx, plan)?;
    progress.checkpoint("initial")?;

    let oracle = assimilate(prior, obs, cx)?;
    let reconstructed = result.belief.to_dense_covariance();
    let oracle_cov = oracle.covariance();

    let mut max_abs_diff = 0.0_f64;
    let mut scale = 0.0_f64;
    for (i, row) in reconstructed.iter().enumerate() {
        for (j, entry) in row.iter().enumerate() {
            let diff = (entry - oracle_cov[i][j]).abs();
            progress.scalar("factor-checker", 1)?;
            if diff > max_abs_diff {
                max_abs_diff = diff;
            }
            let magnitude = oracle_cov[i][j].abs();
            if magnitude > scale {
                scale = magnitude;
            }
        }
    }
    let tolerance = 64.0_f64 * f64::EPSILON * scale.max(1.0) * (n as f64);
    let agrees = max_abs_diff <= tolerance;

    // Independent recomputation from dense inputs.
    let h = obs.operator();
    let mut ph = Vec::with_capacity(n);
    for row in prior.covariance() {
        ph.push(
            compensated_dot(
                row,
                h,
                "checker covariance-times-operator",
                "factor-checker",
                &mut progress,
            )?
            .0,
        );
    }
    let (measurement_variance, _env) = compensated_dot(
        h,
        &ph,
        "checker measurement variance",
        "factor-checker",
        &mut progress,
    )?;
    let recomputed_innovation = measurement_variance + obs.noise_var();
    let recomputed_identity = measurement_variance * obs.noise_var() / recomputed_innovation;
    let reported_post = result.receipt.measurement_variance_posterior();
    let identity_width = (result.receipt.measurement_identity_enclosure().1
        - result.receipt.measurement_identity_enclosure().0)
        .abs();
    let identity_consistent = (reported_post - recomputed_identity).abs()
        <= identity_width.max(64.0_f64 * f64::EPSILON * recomputed_identity.abs().max(1.0));
    let receipt_consistent = (result.receipt.innovation_variance() - recomputed_innovation).abs()
        <= 64.0_f64 * f64::EPSILON * recomputed_innovation.abs().max(1.0);

    let verdict = if agrees && identity_consistent && receipt_consistent {
        CheckVerdict::Verified
    } else {
        CheckVerdict::Discrepancy
    };
    progress.checkpoint("finalize")?;
    Ok(IndependentCheck {
        verdict,
        max_abs_diff,
        tolerance,
        recomputed_innovation_variance: recomputed_innovation,
        identity_consistent,
    })
}

/// The complete no-mock scalar path: dense prior to factor belief to
/// stable update to independent checker.
#[derive(Debug, Clone, PartialEq)]
pub struct CheckedAssimilation {
    posterior: FactorBelief,
    receipt: ContractionReceipt,
    check: IndependentCheck,
}

impl CheckedAssimilation {
    /// The factor-form posterior belief.
    #[must_use]
    pub const fn posterior(&self) -> &FactorBelief {
        &self.posterior
    }

    /// The executable contraction receipt.
    #[must_use]
    pub const fn receipt(&self) -> &ContractionReceipt {
        &self.receipt
    }

    /// The independent checker result.
    #[must_use]
    pub const fn check(&self) -> &IndependentCheck {
        &self.check
    }
}

/// Compose the production scalar path end to end: factor the dense prior,
/// apply the stable scalar update, and independently check the result
/// against the dense oracle. This is the no-mock lane bead sj31i.37
/// demands: every advertised scalar contraction carries executable evidence
/// and an independent checker verdict.
///
/// # Errors
/// Returns [`AssimError`] from any composed stage; a numerical verdict is
/// carried as data, never as an error.
pub fn assimilate_belief_scalar_checked(
    prior: &Belief,
    obs: &Observation,
    cx: &Cx<'_>,
) -> Result<CheckedAssimilation, AssimError> {
    let factor_prior = FactorBelief::from_belief(prior, cx)?;
    let result = assimilate_scalar(&factor_prior, obs, cx)?;
    let check = verify_factor_assimilation(prior, obs, &result, cx)?;
    Ok(CheckedAssimilation {
        posterior: result.belief,
        receipt: result.receipt,
        check,
    })
}

/// Structured log suite identity for the scalar factor lane.
pub const SCALAR_FACTOR_LOG_SUITE: &str = "fs-assimilate/scalar-factor";

/// Emit the bounded structured log for one checked scalar assimilation.
///
/// The detail line carries the bead-mandated fields without unbounded
/// payloads: state dimension, method identity, conditioning indicators
/// (innovation enclosure width, max pivot ratio), residual/contraction/
/// misfit bounds (identity enclosure, misfit verdict), contraction and
/// checker dispositions, refusal/no-claim reason when the verdict is not
/// certified, and the content roots that retain prior/observation/verdict.
/// The event is emitted through `fs_obs` and the line is validated against
/// the wire schema; a failure event additionally passes the failure-record
/// lint.
///
/// # Errors
/// Returns [`fs_obs::SchemaError`] when the constructed event or line
/// violates the observability wire schema.
pub fn emit_checked_assimilation_log(
    checked: &CheckedAssimilation,
    dim: usize,
    planned_work: u128,
    emitter: &mut fs_obs::Emitter,
) -> Result<fs_obs::Event, fs_obs::SchemaError> {
    let receipt = checked.receipt();
    let check = checked.check();
    let pass =
        receipt.state() == ContractionState::Certified && check.verdict() == CheckVerdict::Verified;
    let enclosure_width = (receipt.measurement_identity_enclosure().1
        - receipt.measurement_identity_enclosure().0)
        .abs();
    let detail = format!(
        "{{\"dim\":{dim},\"method\":\"bierman-ud/v1\",\
         \"contraction\":\"{}\",\"misfit\":\"{}\",\"checker\":\"{}\",\
         \"innovation_variance\":{:e},\"identity_enclosure_width\":{:e},\
         \"max_pivot_ratio\":{:e},\"misfit_before\":{:e},\"misfit_after\":{:e},\
         \"oracle_max_abs_diff\":{:e},\"oracle_tolerance\":{:e},\
         \"planned_work\":{planned_work},\"receipt_root\":\"{}\",\
         \"no_claim\":\"{}\"}}",
        receipt.state().name(),
        receipt.misfit_verdict().name(),
        check.verdict().name(),
        receipt.innovation_variance(),
        enclosure_width,
        receipt.max_pivot_ratio(),
        receipt.misfit_before(),
        receipt.misfit_after(),
        check.max_abs_diff(),
        check.tolerance(),
        receipt.identity(),
        if receipt.state() == ContractionState::Certified {
            "none"
        } else {
            "contraction-not-advertised"
        },
    );
    let event = emitter.emit(
        if pass {
            fs_obs::Severity::Info
        } else {
            fs_obs::Severity::Error
        },
        fs_obs::EventKind::ConformanceCase {
            suite: SCALAR_FACTOR_LOG_SUITE.to_string(),
            case: "sensor-update-checked".to_string(),
            pass,
            detail,
            seed: 0,
        },
        None,
    );
    if !pass {
        fs_obs::lint_failure_record(&event)?;
    }
    let line = event.to_jsonl();
    fs_obs::validate_line(&line)?;
    Ok(event)
}
