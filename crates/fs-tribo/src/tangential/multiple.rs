//! Simultaneous original tangential laws on one full mechanical load line.
//!
//! Each contact keeps two tangent axes. Off-diagonal contact blocks are physical
//! recoil and are retained in `du = du_free - C F`. Block Gauss--Seidel calls the
//! existing two-axis matrix adapter, always from the original accepted history.
//! No sweep advances constitutive time or deposits heat a second time.
//!
//! The prepared matrix is immutable. Its pivoted PSD factor is used only for
//! admission and a supported minimum-norm zero-slip reaction, which seeds the
//! original laws before kinetic alternatives. Physical matrix entries are never
//! replaced by that factor or a regularized response. Contact order is the
//! deterministic block/root-selection order; nonunique friction problems do
//! not acquire a uniqueness or convergence claim from this numerical embedding.

use core::fmt;

use fs_math::det;

use super::matrix::{
    MatrixTangentialConfig, MatrixTangentialError, advance_with_matrix_compliance,
};
use super::{TangentialLaw, TangentialStep, TangentialStepInput};

/// Explicit preparation, nonlinear-work, local and collective closure budgets.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MultipleTangentialConfig {
    /// Contact ceiling, 1 through 128; each contact has exactly two axes.
    pub maximum_contacts: usize,
    /// Ceiling on supplied row-major matrix entries, 4 through 65536.
    /// Preparation retains at most three such dense arrays plus O(contacts).
    pub maximum_matrix_terms: usize,
    /// Maximum complete deterministic block sweeps, 1 through 4096.
    pub maximum_sweeps: usize,
    /// Total original-law evaluations across every contact and sweep.
    /// Between 1 and 1048576; a remaining global budget caps each local solve.
    pub maximum_law_evaluations: usize,
    /// Original two-axis solve budget; its displacement and work tolerances also
    /// gate each contact against the complete global matrix after every sweep.
    pub local: MatrixTangentialConfig,
    /// Absolute tolerance on summed storage plus heat minus actual work [J].
    pub total_energy_tolerance_j: f64,
}

impl MultipleTangentialConfig {
    /// Validate controls without allocating or evaluating a constitutive law.
    pub fn validate(self) -> Result<(), MultipleTangentialError> {
        if !(1..=128).contains(&self.maximum_contacts)
            || !(4..=65536).contains(&self.maximum_matrix_terms)
            || !(1..=4096).contains(&self.maximum_sweeps)
            || !(1..=1048576).contains(&self.maximum_law_evaluations)
            || !self.total_energy_tolerance_j.is_finite()
            || self.total_energy_tolerance_j <= 0.0
            || self.local.validate().is_err()
        {
            return Err(MultipleTangentialError::Invalid(
                "multiple tangential budgets",
            ));
        }
        Ok(())
    }
}

/// One borrowed law/history and one original released-motion/environment input.
pub struct TangentialContactInput<'a, L: TangentialLaw> {
    /// Each site may use a different instance/parameters of the same law type.
    pub law: &'a L,
    /// Original accepted history, unchanged by every trial and refusal.
    pub state: &'a L::State,
    /// Zero scalar compliance; the prepared matrix owns all mechanical recoil.
    /// Every site must use the same positive held interval.
    pub input: TangentialStepInput,
}

/// All candidates publish together only after original global work closes.
#[derive(Clone, Debug, PartialEq)]
pub struct MultipleTangentialStep<S> {
    /// Original law candidates in caller contact order. Compatibility residuals
    /// use the complete physical matrix, including every cross-contact block.
    pub contacts: Vec<TangentialStep<S>>,
    /// Local `Delta H_i + D_i - F_i dot (du_free - C F)_i` [J].
    pub contact_work_residuals_j: Vec<f64>,
    /// Complete block sweeps performed.
    pub sweeps: usize,
    /// Original-law calls spent by all local solves and their backtracking.
    pub law_evaluations: usize,
    /// Largest contact's Euclidean full-matrix displacement residual [m].
    pub maximum_compatibility_residual_m: f64,
    /// Sum of original reversible contact-energy changes [J].
    pub stored_change_j: f64,
    /// Sum of original physical dissipated work [J], never an iteration deficit.
    pub dissipated_j: f64,
    /// Summed storage plus heat minus complete actual interface work [J].
    pub body_work_residual_j: f64,
}

/// A failure publishes no contact history and no partial mechanical update.
#[derive(Clone, Debug, PartialEq)]
pub enum MultipleTangentialError {
    /// Invalid shape, budgets, nonpassive matrix or inconsistent physical input.
    Invalid(&'static str),
    /// Original local law/matrix refusal, retaining its contact index and cause.
    Contact {
        /// Original caller contact index.
        index: usize,
        /// Original local constitutive or matrix-solve refusal.
        error: MatrixTangentialError,
    },
    /// Caller cancellation observed at a bounded preparation/solve boundary.
    Cancelled,
    /// Global compatibility/work did not close under the caller's work budget.
    Unresolved {
        /// Complete block sweeps performed before refusal.
        sweeps: usize,
        /// Original-law calls spent by completed local solves.
        law_evaluations: usize,
        /// Largest local displacement or work-equivalent displacement [m].
        residual_m: f64,
        /// Signed total actual work residual [J].
        work_residual_j: f64,
    },
}

impl fmt::Display for MultipleTangentialError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "simultaneous tangential contacts: {self:?}")
    }
}
impl std::error::Error for MultipleTangentialError {}

/// Immutable full mechanical response admitted and factored once for reuse.
#[derive(Clone, Debug, PartialEq)]
pub struct PreparedTangentialContacts {
    count: usize,
    compliance: Vec<f64>,
    diagonal: Vec<[[f64; 2]; 2]>,
    factor: PositiveFactor,
    config: MultipleTangentialConfig,
}

impl PreparedTangentialContacts {
    /// Admit a full row-major `(2 * count)^2` compliance matrix [m/N].
    ///
    /// Contact order and tangent axes remain fixed across reuse. A changed body
    /// response requires a new owner; preparation does not change stored shear.
    /// Empty contact sets belong to the mechanics owner's no-contact path.
    pub fn new(
        count: usize,
        compliance_m_per_n: Vec<f64>,
        config: MultipleTangentialConfig,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, MultipleTangentialError> {
        poll(&mut cancelled)?;
        config.validate()?;
        if count == 0 || count > config.maximum_contacts {
            return Err(MultipleTangentialError::Invalid("contact count"));
        }
        let dimension = 2 * count;
        let terms = dimension * dimension;
        if terms > config.maximum_matrix_terms || compliance_m_per_n.len() != terms {
            return Err(MultipleTangentialError::Invalid(
                "full compliance shape/term budget",
            ));
        }
        let factor = PositiveFactor::new(
            &compliance_m_per_n,
            dimension,
            config.local.matrix_relative_tolerance,
            &mut cancelled,
        )?;
        let diagonal = (0..count)
            .map(|i| {
                let row = 2 * i;
                [
                    [
                        compliance_m_per_n[row * dimension + row],
                        compliance_m_per_n[row * dimension + row + 1],
                    ],
                    [
                        compliance_m_per_n[(row + 1) * dimension + row],
                        compliance_m_per_n[(row + 1) * dimension + row + 1],
                    ],
                ]
            })
            .collect();
        poll(&mut cancelled)?;
        Ok(Self {
            count,
            compliance: compliance_m_per_n,
            diagonal,
            factor,
            config,
        })
    }

    /// Number of two-axis contact sites.
    #[must_use]
    pub const fn contact_count(&self) -> usize {
        self.count
    }
    /// Complete original physical matrix, in unchanged caller coordinate order.
    #[must_use]
    pub fn compliance_m_per_n(&self) -> &[f64] {
        &self.compliance
    }
    /// Move the original matrix to a mechanical result without a dense copy.
    #[must_use]
    pub fn into_compliance_m_per_n(self) -> Vec<f64> {
        self.compliance
    }
    /// Exact explicit preparation and solve controls.
    #[must_use]
    pub const fn config(&self) -> MultipleTangentialConfig {
        self.config
    }

    /// Resolve all original laws against all cross-contact recoil together.
    ///
    /// Each native call borrows the same original accepted history. The local
    /// released motion subtracts other contacts' latest forces, while the
    /// existing matrix adapter owns the whole diagonal 2-by-2 response. After
    /// each sweep the complete matrix independently reconstructs every actual
    /// displacement and local work, followed by one collective energy gate.
    ///
    /// The supported global zero-slip reaction is tried first. When forces are
    /// nonunique, its minimum-norm seed and original caller block order specify
    /// the numerical selection; they do not solve every friction-cone feasibility
    /// problem. A source/force-budget refusal never becomes a favorable mode.
    pub fn advance<L: TangentialLaw>(
        &self,
        sites: &[TangentialContactInput<'_, L>],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MultipleTangentialStep<L::State>, MultipleTangentialError> {
        poll(&mut cancelled)?;
        if sites.len() != self.count {
            return Err(MultipleTangentialError::Invalid("contact input count"));
        }
        let dt = sites[0].input.dt_s;
        let mut released = Vec::with_capacity(2 * self.count);
        let mut energy_before = Vec::with_capacity(self.count);
        for (index, site) in sites.iter().enumerate() {
            poll(&mut cancelled)?;
            if !dt.is_finite()
                || dt <= 0.0
                || site.input.dt_s != dt
                || site.input.compliance_m_per_n != 0.0
                || site.input.slip_increment_m.iter().any(|x| !x.is_finite())
            {
                return Err(MultipleTangentialError::Invalid(
                    "shared interval and matrix-owned released motion",
                ));
            }
            site.input
                .solve
                .validate()
                .map_err(|error| MultipleTangentialError::Contact {
                    index,
                    error: MatrixTangentialError::Law(error),
                })?;
            let energy = number(site.law.stored_energy_j(site.state))?;
            if energy < 0.0 {
                return Err(MultipleTangentialError::Invalid("negative contact storage"));
            }
            energy_before.push(energy);
            released.extend_from_slice(&site.input.slip_increment_m);
        }
        let dimension = released.len();
        let mut forces = self.factor.reaction(&released, &mut cancelled)?;
        let supported = (0..self.count).all(|i| {
            let r = [
                released[2 * i] - row_product(&self.compliance, dimension, 2 * i, &forces),
                released[2 * i + 1] - row_product(&self.compliance, dimension, 2 * i + 1, &forces),
            ];
            det::hypot(r[0], r[1]) <= self.config.local.displacement_tolerance_m
        });
        if !supported {
            forces.fill(0.0);
        }
        let mut candidates: Vec<Option<TangentialStep<L::State>>> =
            (0..self.count).map(|_| None).collect();
        let mut work_residuals = vec![0.0; self.count];
        let mut evaluations = 0;
        let mut last_residual = f64::INFINITY;
        let mut last_work = f64::INFINITY;
        for sweep in 1..=self.config.maximum_sweeps {
            for (index, site) in sites.iter().enumerate() {
                poll(&mut cancelled)?;
                let remaining = self.config.maximum_law_evaluations - evaluations;
                if remaining == 0 {
                    return Err(unresolved(sweep - 1, evaluations, last_residual, last_work));
                }
                let mut input = site.input;
                for axis in 0..2 {
                    let row = 2 * index + axis;
                    let mut recoil = 0.0;
                    for (column, &force) in forces.iter().enumerate() {
                        if column / 2 != index {
                            recoil =
                                number(recoil + self.compliance[row * dimension + column] * force)?;
                        }
                    }
                    input.slip_increment_m[axis] = number(released[row] - recoil)?;
                }
                let local = MatrixTangentialConfig {
                    maximum_evaluations: remaining.min(self.config.local.maximum_evaluations),
                    ..self.config.local
                };
                let step = advance_with_matrix_compliance(
                    site.law,
                    site.state,
                    &input,
                    self.diagonal[index],
                    local,
                    &mut cancelled,
                )
                .map_err(|error| match error {
                    MatrixTangentialError::Cancelled => MultipleTangentialError::Cancelled,
                    error => MultipleTangentialError::Contact { index, error },
                })?;
                evaluations += step.law_evaluations;
                forces[2 * index..2 * index + 2].copy_from_slice(&step.contact.traction_n);
                candidates[index] = Some(step.contact);
            }
            let mut stored = 0.0;
            let mut heat = 0.0;
            let mut work = 0.0;
            let mut maximum_residual = 0.0_f64;
            let mut maximum_work = 0.0_f64;
            for (index, site) in sites.iter().enumerate() {
                poll(&mut cancelled)?;
                let candidate = candidates[index].as_mut().expect("complete contact sweep");
                let mut actual = [0.0; 2];
                for axis in 0..2 {
                    let row = 2 * index + axis;
                    actual[axis] = number(
                        released[row] - row_product(&self.compliance, dimension, row, &forces),
                    )?;
                    candidate.compatibility_residual_m[axis] =
                        number(candidate.slip_increment_m[axis] - actual[axis])?;
                }
                let change =
                    number(site.law.stored_energy_j(&candidate.state) - energy_before[index])?;
                if number(change - candidate.stored_change_j)?.abs()
                    > self.config.local.energy_tolerance_j
                    || !candidate.dissipated_j.is_finite()
                    || candidate.dissipated_j < 0.0
                {
                    return Err(MultipleTangentialError::Invalid(
                        "original contact storage/heat",
                    ));
                }
                let local_work = number(dot(candidate.traction_n, actual))?;
                work_residuals[index] = number(change + candidate.dissipated_j - local_work)?;
                maximum_residual = maximum_residual.max(number(det::hypot(
                    candidate.compatibility_residual_m[0],
                    candidate.compatibility_residual_m[1],
                ))?);
                maximum_work = maximum_work.max(work_residuals[index].abs());
                stored = number(stored + change)?;
                heat = number(heat + candidate.dissipated_j)?;
                work = number(work + local_work)?;
            }
            last_work = number(stored + heat - work)?;
            last_residual = maximum_residual
                .max(
                    self.config.local.displacement_tolerance_m
                        * (maximum_work / self.config.local.energy_tolerance_j),
                )
                .max(
                    self.config.local.displacement_tolerance_m
                        * (last_work.abs() / self.config.total_energy_tolerance_j),
                );
            if maximum_residual <= self.config.local.displacement_tolerance_m
                && maximum_work <= self.config.local.energy_tolerance_j
                && last_work.abs() <= self.config.total_energy_tolerance_j
            {
                poll(&mut cancelled)?;
                return Ok(MultipleTangentialStep {
                    contacts: candidates
                        .into_iter()
                        .map(|candidate| candidate.expect("complete contact sweep"))
                        .collect(),
                    contact_work_residuals_j: work_residuals,
                    sweeps: sweep,
                    law_evaluations: evaluations,
                    maximum_compatibility_residual_m: maximum_residual,
                    stored_change_j: stored,
                    dissipated_j: heat,
                    body_work_residual_j: last_work,
                });
            }
        }
        Err(unresolved(
            self.config.maximum_sweeps,
            evaluations,
            last_residual,
            last_work,
        ))
    }
}

// Complete-pivot Cholesky admission, normalized by the original largest entry.
// A discarded tail has infinity norm at most the explicit arithmetic tolerance;
// the original matrix, not this factor, is used in every physical residual.
#[derive(Clone, Debug, PartialEq)]
struct PositiveFactor {
    dimension: usize,
    scale: f64,
    rank: usize,
    permutation: Vec<usize>,
    lower: Vec<f64>,
    gram_lower: Vec<f64>,
}

impl PositiveFactor {
    fn new(
        matrix: &[f64],
        n: usize,
        tolerance: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, MultipleTangentialError> {
        let mut scale = 0.0_f64;
        for row in matrix.chunks(n) {
            poll(cancelled)?;
            for &value in row {
                scale = scale.max(number(value)?.abs());
            }
        }
        let mut result = Self {
            dimension: n,
            scale,
            rank: 0,
            permutation: (0..n).collect(),
            lower: vec![0.0; n * n],
            gram_lower: Vec::new(),
        };
        if scale == 0.0 {
            return Ok(result);
        }
        let mut work = vec![0.0; n * n];
        for i in 0..n {
            poll(cancelled)?;
            for j in 0..n {
                let a = matrix[i * n + j] / scale;
                let b = matrix[j * n + i] / scale;
                if (a - b).abs() > tolerance {
                    return Err(MultipleTangentialError::Invalid(
                        "nonsymmetric full compliance",
                    ));
                }
                work[i * n + j] = 0.5 * a + 0.5 * b;
            }
        }
        for k in 0..n {
            let mut tail_norm = 0.0_f64;
            let mut pivot = k;
            for i in k..n {
                poll(cancelled)?;
                let mut sum = 0.0;
                for j in k..n {
                    sum = number(sum + work[i * n + j].abs())?;
                }
                tail_norm = tail_norm.max(sum);
                if work[i * n + i] > work[pivot * n + pivot] {
                    pivot = i;
                }
            }
            if tail_norm <= tolerance {
                break;
            }
            if work[pivot * n + pivot] <= 0.0 {
                return Err(MultipleTangentialError::Invalid(
                    "nonpositive-semidefinite full compliance",
                ));
            }
            if pivot != k {
                for j in 0..n {
                    work.swap(k * n + j, pivot * n + j);
                }
                for i in 0..n {
                    work.swap(i * n + k, i * n + pivot);
                }
                for j in 0..k {
                    result.lower.swap(k * n + j, pivot * n + j);
                }
                result.permutation.swap(k, pivot);
            }
            let diagonal = number(det::sqrt(work[k * n + k]))?;
            result.lower[k * n + k] = diagonal;
            for i in k + 1..n {
                result.lower[i * n + k] = number(work[i * n + k] / diagonal)?;
            }
            for i in k + 1..n {
                poll(cancelled)?;
                for j in k + 1..=i {
                    let value = number(
                        work[i * n + j] - result.lower[i * n + k] * result.lower[j * n + k],
                    )?;
                    work[i * n + j] = value;
                    work[j * n + i] = value;
                }
            }
            result.rank += 1;
        }
        // The Schur scratch is no longer needed. Release it before allocating
        // the range Gram factor so preparation stays within three dense arrays
        // including the caller's retained physical compliance.
        drop(work);
        let r = result.rank;
        if r > 0 && r < n {
            result.gram_lower = vec![0.0; r * r];
            for i in 0..r {
                for j in 0..=i {
                    poll(cancelled)?;
                    let mut value = 0.0;
                    for row in 0..n {
                        value =
                            number(value + result.lower[row * n + i] * result.lower[row * n + j])?;
                    }
                    for k in 0..j {
                        value = number(
                            value - result.gram_lower[i * r + k] * result.gram_lower[j * r + k],
                        )?;
                    }
                    result.gram_lower[i * r + j] = if i == j {
                        if value <= 0.0 {
                            return Err(MultipleTangentialError::Invalid(
                                "unresolved supported compliance factor",
                            ));
                        }
                        det::sqrt(value)
                    } else {
                        number(value / result.gram_lower[j * r + j])?
                    };
                }
            }
        }
        poll(cancelled)?;
        Ok(result)
    }

    fn reaction(
        &self,
        released: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<f64>, MultipleTangentialError> {
        let n = self.dimension;
        let r = self.rank;
        let mut result = vec![0.0; n];
        if r == 0 {
            return Ok(result);
        }
        let mut y = vec![0.0; r];
        for i in 0..r {
            poll(cancelled)?;
            let mut value = released[self.permutation[i]];
            for j in 0..i {
                value -= self.lower[i * n + j] * y[j];
            }
            y[i] = value / self.lower[i * n + i];
        }
        if r == n {
            for i in (0..r).rev() {
                poll(cancelled)?;
                let mut value = y[i];
                for j in i + 1..r {
                    value -= self.lower[j * n + i] * y[j];
                }
                y[i] = value / self.lower[i * n + i];
            }
            for i in 0..n {
                result[self.permutation[i]] = y[i] / self.scale;
            }
        } else {
            // L y = u gives the supported displacement. Minimum-norm force is
            // L (L^T L)^-1 y / scale; no arbitrary independent rows are forced.
            for i in 0..r {
                poll(cancelled)?;
                let mut value = y[i];
                for j in 0..i {
                    value -= self.gram_lower[i * r + j] * y[j];
                }
                y[i] = value / self.gram_lower[i * r + i];
            }
            for i in (0..r).rev() {
                poll(cancelled)?;
                let mut value = y[i];
                for j in i + 1..r {
                    value -= self.gram_lower[j * r + i] * y[j];
                }
                y[i] = value / self.gram_lower[i * r + i];
            }
            for row in 0..n {
                poll(cancelled)?;
                let mut value = 0.0;
                for j in 0..r {
                    value += self.lower[row * n + j] * y[j];
                }
                result[self.permutation[row]] = value / self.scale;
            }
        }
        // An unrepresentable numerical seed need not rule out a bounded kinetic
        // solution. It never changes the physical matrix or accepted history.
        if result.iter().any(|x| !x.is_finite()) {
            result.fill(0.0);
        }
        Ok(result)
    }
}

fn row_product(matrix: &[f64], n: usize, row: usize, vector: &[f64]) -> f64 {
    matrix[row * n..(row + 1) * n]
        .iter()
        .zip(vector)
        .map(|(&a, &b)| a * b)
        .sum()
}
fn dot(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0] * b[0] + a[1] * b[1]
}
fn number(value: f64) -> Result<f64, MultipleTangentialError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(MultipleTangentialError::Invalid(
            "nonfinite multiple-contact arithmetic",
        ))
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), MultipleTangentialError> {
    if cancelled() {
        Err(MultipleTangentialError::Cancelled)
    } else {
        Ok(())
    }
}
fn unresolved(
    sweeps: usize,
    law_evaluations: usize,
    residual_m: f64,
    work_residual_j: f64,
) -> MultipleTangentialError {
    MultipleTangentialError::Unresolved {
        sweeps,
        law_evaluations,
        residual_m,
        work_residual_j,
    }
}
