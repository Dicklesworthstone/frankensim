//! Bounded interface quasi-Newton acceleration with inverse least squares.
//!
//! For an interface map `g = G(x)`, keep differences of the *vector* residual
//! `r = g - x` and of `g`. Solve `min_c ||V c + r||_2`, then propose `g + W c`.
//! This is IQN-ILS (Anderson acceleration with unit mixing after startup), not
//! a scalar secant on a signed mean: opposing interface modes cannot cancel
//! before entering the least-squares problem.
//!
//! A twice-orthogonalized modified Gram-Schmidt factorization drops dependent
//! columns using a relative test. Normal equations are never formed. Newest
//! secants have priority, memory is bounded by the configured history, and an
//! empty usable history falls back to the caller's relaxation. All coordinates
//! must use a compatible scale; physical bounds, convergence, energy balance,
//! and cancellation remain the owning driver's responsibility.

#![forbid(unsafe_code)]

use core::fmt;

/// Hard bound on the small dense least-squares problem.
pub const MAX_HISTORY: usize = 64;

/// Memory and rank-filter policy for an interface accelerator.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IqnIlsConfig {
    /// Maximum retained secants, in `1..=MAX_HISTORY`.
    pub max_history: usize,
    /// Discard a column when its orthogonal norm is at most this fraction of
    /// its original norm. Must be finite and strictly between zero and one.
    pub relative_rank_tolerance: f64,
}

impl Default for IqnIlsConfig {
    fn default() -> Self {
        Self {
            max_history: 8,
            relative_rank_tolerance: 1.0e-10,
        }
    }
}

/// A rejected input or an unrepresentable intermediate. Failed steps do not
/// change either the history or the previous accepted sample.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum IqnIlsError {
    /// An interface must have at least one coordinate.
    EmptyInterface,
    /// The history limit is zero or exceeds [`MAX_HISTORY`].
    InvalidHistoryLimit(usize),
    /// An invalid rank threshold, represented without a floating-point Eq.
    InvalidRankTolerance(u64),
    /// Fallback relaxation must be finite; its magnitude is driver policy.
    InvalidRelaxation(u64),
    /// Each input must have the dimension declared at construction.
    DimensionMismatch {
        /// The declared dimension.
        expected: usize,
        /// The observed input length.
        found: usize,
    },
    /// A non-finite input or intermediate, with its first responsible stage.
    NonFinite {
        /// Stable arithmetic-stage name.
        stage: &'static str,
        /// IEEE-754 representation of the rejected value.
        value_bits: u64,
    },
}

impl fmt::Display for IqnIlsError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "IQN-ILS: {self:?}")
    }
}

impl std::error::Error for IqnIlsError {}

/// An interface update, not a convergence or conservation certificate.
#[derive(Debug, Clone, PartialEq)]
pub struct IqnIlsStep {
    /// The proposed next interface state.
    pub values: Vec<f64>,
    /// Number of independent residual-difference columns used.
    pub used_columns: usize,
    /// Fallback omega when no secant is usable; otherwise one, the weight of
    /// `G(x)` in `G(x) + W c`. The latter update is not scalar relaxation.
    pub relaxation_omega: f64,
}

#[derive(Debug, Clone, PartialEq)]
struct Sample {
    residual: Vec<f64>,
    image: Vec<f64>,
}

#[derive(Debug, Clone, PartialEq)]
struct Secant {
    residual: Vec<f64>,
    image: Vec<f64>,
}

/// Accelerator state for one fixed interface ordering and one fixed map.
///
/// Clone the complete value to retain acceleration history for replay. Call
/// [`Self::reset`] when the map, coordinate ordering, or scaling changes; a
/// reference vector alone is not an IQN-ILS checkpoint.
#[derive(Debug, Clone, PartialEq)]
pub struct IqnIls {
    dimension: usize,
    config: IqnIlsConfig,
    previous: Option<Sample>,
    history: Vec<Secant>,
}

impl IqnIls {
    /// Construct an empty bounded history. No interface-sized allocation is
    /// made until the first step.
    ///
    /// # Errors
    /// Refuses an empty interface, invalid history limit, or invalid threshold.
    pub fn new(dimension: usize, config: IqnIlsConfig) -> Result<Self, IqnIlsError> {
        if dimension == 0 {
            return Err(IqnIlsError::EmptyInterface);
        }
        if !(1..=MAX_HISTORY).contains(&config.max_history) {
            return Err(IqnIlsError::InvalidHistoryLimit(config.max_history));
        }
        let tolerance = config.relative_rank_tolerance;
        if !(tolerance.is_finite() && tolerance > 0.0 && tolerance < 1.0) {
            return Err(IqnIlsError::InvalidRankTolerance(tolerance.to_bits()));
        }
        Ok(Self {
            dimension,
            config,
            previous: None,
            history: Vec::new(),
        })
    }

    /// Discard all history before changing the interface map or coordinates.
    pub fn reset(&mut self) {
        self.previous = None;
        self.history.clear();
    }

    /// Number of retained secants (some may be filtered at the next step).
    #[must_use]
    pub fn history_len(&self) -> usize {
        self.history.len()
    }

    /// Propose an update from the actual map evaluation `image = G(current)`.
    ///
    /// `fallback_omega` is used only when no independent secant is available.
    /// Any finite factor is accepted, including signed Aitken factors and zero
    /// (an explicit no-op fallback). The driver owns relaxation bounds and the
    /// stopping rule. Exact zero residual returns the map without extrapolation.
    ///
    /// # Errors
    /// Refuses wrong lengths, invalid fallback relaxation, non-finite inputs,
    /// or unrepresentable differences, factorization, or updates. No failure
    /// changes the accepted sample/history, so retry is deterministic.
    pub fn step(
        &mut self,
        current: &[f64],
        image: &[f64],
        fallback_omega: f64,
    ) -> Result<IqnIlsStep, IqnIlsError> {
        self.check_shape(current)?;
        self.check_shape(image)?;
        if !fallback_omega.is_finite() {
            return Err(IqnIlsError::InvalidRelaxation(fallback_omega.to_bits()));
        }
        for &value in current {
            finite("current interface", value)?;
        }
        for &value in image {
            finite("interface map", value)?;
        }
        let residual = difference(image, current, "interface residual")?;
        let candidate = self
            .previous
            .as_ref()
            .map(|previous| -> Result<Secant, IqnIlsError> {
                Ok(Secant {
                    residual: difference(&residual, &previous.residual, "residual difference")?,
                    image: difference(image, &previous.image, "map difference")?,
                })
            })
            .transpose()?
            .filter(|secant: &Secant| max_abs(&secant.residual) > 0.0);

        // Borrow old columns while building the proposal. Commit state only
        // after *all* arithmetic succeeds; no rollback copies are needed.
        let columns: Vec<&Secant> = candidate
            .iter()
            .chain(self.history.iter().rev())
            .take(self.config.max_history)
            .collect();
        let proposal = self.propose(current, image, &residual, fallback_omega, &columns)?;

        if let Some(secant) = candidate {
            if self.history.len() == self.config.max_history {
                self.history.remove(0);
            }
            self.history.push(secant);
        }
        self.previous = Some(Sample {
            residual,
            image: image.to_vec(),
        });
        Ok(proposal)
    }

    fn check_shape(&self, values: &[f64]) -> Result<(), IqnIlsError> {
        if values.len() != self.dimension {
            return Err(IqnIlsError::DimensionMismatch {
                expected: self.dimension,
                found: values.len(),
            });
        }
        Ok(())
    }

    fn propose(
        &self,
        current: &[f64],
        image: &[f64],
        residual: &[f64],
        fallback_omega: f64,
        columns: &[&Secant],
    ) -> Result<IqnIlsStep, IqnIlsError> {
        let residual_scale = max_abs(residual);
        if residual_scale == 0.0 {
            return Ok(IqnIlsStep {
                values: image.to_vec(),
                used_columns: 0,
                relaxation_omega: fallback_omega,
            });
        }
        let qr = FilteredQr::factor(columns, self.config.relative_rank_tolerance)?;
        if qr.q.is_empty() {
            let values = current
                .iter()
                .zip(residual)
                .map(|(&x, &r)| finite("fallback update", x + fallback_omega * r))
                .collect::<Result<Vec<_>, _>>()?;
            return Ok(IqnIlsStep {
                values,
                used_columns: 0,
                relaxation_omega: fallback_omega,
            });
        }

        // Scale the RHS independently to avoid overflowing its dot products.
        // The column and RHS scalings are undone in the map correction.
        let rhs: Vec<f64> = residual.iter().map(|&r| -r / residual_scale).collect();
        let coefficients = qr.solve(&rhs)?;
        let mut values = Vec::with_capacity(self.dimension);
        for (row, &mapped) in image.iter().enumerate() {
            let mut correction = 0.0;
            for (column, &coefficient) in qr.w.iter().zip(&coefficients) {
                correction = finite("map correction", correction + column[row] * coefficient)?;
            }
            values.push(finite("accelerated update", mapped + correction * residual_scale)?);
        }
        Ok(IqnIlsStep {
            values,
            used_columns: qr.q.len(),
            relaxation_omega: 1.0,
        })
    }
}

struct FilteredQr {
    q: Vec<Vec<f64>>,
    // Upper triangular R, stored by columns.
    r: Vec<Vec<f64>>,
    // Map differences scaled by the same factor as the accepted V columns.
    w: Vec<Vec<f64>>,
}

impl FilteredQr {
    fn factor(columns: &[&Secant], tolerance: f64) -> Result<Self, IqnIlsError> {
        let mut qr = Self {
            q: Vec::new(),
            r: Vec::new(),
            w: Vec::new(),
        };
        for column in columns {
            // No interface can have more independent modes than coordinates,
            // even if a caller selects a threshold below round-off.
            if qr.q.len() == column.residual.len() {
                break;
            }
            let scale = max_abs(&column.residual);
            if scale == 0.0 {
                continue;
            }
            let mut v: Vec<f64> = column.residual.iter().map(|&value| value / scale).collect();
            let original_norm = norm(&v)?;
            let mut triangular = vec![0.0; qr.q.len()];
            // Reorthogonalization is necessary near dependent interface modes.
            for _ in 0..2 {
                for (basis, coefficient) in qr.q.iter().zip(&mut triangular) {
                    let projection = dot(basis, &v)?;
                    *coefficient = finite("QR coefficient", *coefficient + projection)?;
                    for (value, &q) in v.iter_mut().zip(basis) {
                        *value = finite("QR orthogonalization", *value - projection * q)?;
                    }
                }
            }
            let orthogonal_norm = norm(&v)?;
            if orthogonal_norm <= tolerance * original_norm {
                continue;
            }
            for value in &mut v {
                *value /= orthogonal_norm;
            }
            let w = column
                .image
                .iter()
                .map(|&value| finite("scaled map difference", value / scale))
                .collect::<Result<Vec<_>, _>>()?;
            triangular.push(orthogonal_norm);
            qr.q.push(v);
            qr.r.push(triangular);
            qr.w.push(w);
        }
        Ok(qr)
    }

    fn solve(&self, rhs: &[f64]) -> Result<Vec<f64>, IqnIlsError> {
        let mut coefficients = self
            .q
            .iter()
            .map(|basis| dot(basis, rhs))
            .collect::<Result<Vec<_>, _>>()?;
        for row in (0..coefficients.len()).rev() {
            let mut value = coefficients[row];
            for (column, &coefficient) in coefficients.iter().enumerate().skip(row + 1) {
                value = finite("QR back substitution", value - self.r[column][row] * coefficient)?;
            }
            coefficients[row] = finite("QR solution", value / self.r[row][row])?;
        }
        Ok(coefficients)
    }
}

fn finite(stage: &'static str, value: f64) -> Result<f64, IqnIlsError> {
    if !value.is_finite() {
        return Err(IqnIlsError::NonFinite {
            stage,
            value_bits: value.to_bits(),
        });
    }
    Ok(value)
}

fn difference(a: &[f64], b: &[f64], stage: &'static str) -> Result<Vec<f64>, IqnIlsError> {
    a.iter().zip(b).map(|(&a, &b)| finite(stage, a - b)).collect()
}

fn max_abs(values: &[f64]) -> f64 {
    values.iter().fold(0.0_f64, |maximum, &value| maximum.max(value.abs()))
}

fn dot(a: &[f64], b: &[f64]) -> Result<f64, IqnIlsError> {
    a.iter().zip(b).try_fold(0.0, |sum, (&a, &b)| finite("QR dot product", sum + a * b))
}

fn norm(values: &[f64]) -> Result<f64, IqnIlsError> {
    finite("QR norm", fs_math::det::sqrt(dot(values, values)?))
}
