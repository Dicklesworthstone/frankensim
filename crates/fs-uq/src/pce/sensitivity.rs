//! Variance-based global sensitivity of an orthonormal Hermite surrogate.
use std::collections::BTreeMap;

use super::PceModel;

/// One nonempty ANOVA support and its fraction of surrogate variance.
#[derive(Clone, Debug, PartialEq)]
pub struct SobolComponent {
    /// Zero-based germ indices, in strictly increasing order. Polynomial
    /// degree is not interaction order: h_5(x_i) still has support [i].
    pub variables: Vec<usize>,
    /// Sum of squared coefficients with exactly this support / variance.
    pub index: f64,
}

/// Global sensitivity under INDEPENDENT standard-normal germs.
///
/// These are algebraic indices of the supplied surrogate, not confidence
/// bounds or a certificate of agreement with the underlying simulator.
/// Correlated physical inputs must not be identified with independent germs:
/// after a KL/whitening map these indices describe the germs, not individual
/// physical coordinates. Components are lexicographically ordered by support;
/// only supports present with nonzero coefficients are stored.
#[derive(Clone, Debug, PartialEq)]
pub struct SobolIndices {
    /// S_i = Var(E[Y | xi_i]) / Var(Y). Includes all univariate degrees.
    pub first_order: Vec<f64>,
    /// S_Ti = sum of all ANOVA components containing germ i.
    pub total_order: Vec<f64>,
    /// Exact-support components (including singletons). They sum to one up
    /// to floating-point rounding. Total-order indices need not sum to one.
    pub components: Vec<SobolComponent>,
}

/// Invalid public PCE data, or an undefined variance normalization.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum PceSensitivityError {
    /// No basis terms were supplied.
    EmptyModel,
    /// Basis and coefficient arrays have different lengths.
    CoefficientCount {
        /// Number of supplied multi-indices.
        terms: usize,
        /// Number of supplied coefficients.
        coefficients: usize,
    },
    /// A basis multi-index does not match the declared germ dimension.
    TermDimension {
        /// Offending term's position in the supplied model.
        term: usize,
        /// Declared germ dimension.
        expected: usize,
        /// Multi-index length.
        actual: usize,
    },
    /// A coefficient (including the constant) is NaN or infinite.
    NonFiniteCoefficient {
        /// Offending coefficient's position.
        term: usize,
    },
    /// Duplicate basis terms cannot be treated as orthogonal contributions.
    DuplicateBasis {
        /// First occurrence in the supplied model.
        first: usize,
        /// Repeated occurrence in the supplied model.
        second: usize,
    },
    /// Every nonconstant coefficient is zero; normalization is undefined.
    ZeroVariance,
}

impl std::fmt::Display for PceSensitivityError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::EmptyModel => write!(f, "PCE sensitivity requires basis terms"),
            Self::CoefficientCount { terms, coefficients } => write!(
                f, "PCE sensitivity has {terms} terms but {coefficients} coefficients"
            ),
            Self::TermDimension { term, expected, actual } => write!(
                f, "PCE sensitivity term {term} has dimension {actual}, expected {expected}"
            ),
            Self::NonFiniteCoefficient { term } => write!(
                f, "PCE sensitivity coefficient {term} is non-finite"
            ),
            Self::DuplicateBasis { first, second } => write!(
                f, "PCE sensitivity terms {first} and {second} have the same basis index"
            ),
            Self::ZeroVariance => write!(
                f, "Sobol indices are undefined for a constant surrogate"
            ),
        }
    }
}

impl std::error::Error for PceSensitivityError {}

#[derive(Clone, Copy, Default)]
struct EnergySum {
    sum: f64,
    correction: f64,
}

impl EnergySum {
    fn add(&mut self, value: f64) {
        let adjusted = value - self.correction;
        let next = self.sum + adjusted;
        self.correction = (next - self.sum) - adjusted;
        self.sum = next;
    }
}

impl PceModel {
    /// Decompose the surrogate's variance into main effects and interactions.
    ///
    /// Orthogonality gives V_u = sum(c_alpha^2 : support(alpha) = u).
    /// No simulator calls, fitting, quadrature or new linear solves are needed.
    /// The coefficients are scaled before squaring so a change of output units
    /// does not overflow/underflow the variance normalization. A very small
    /// component may still round to zero relative to the dominant component.
    ///
    /// Public model fields are checked before use. Unique sparse bases and any
    /// term ordering are accepted; a missing constant means zero mean for this
    /// calculation. Summation follows canonical basis/support order, not the
    /// supplied term order. A constant model returns `ZeroVariance`, not a
    /// misleading all-zero sensitivity report.
    ///
    /// # Errors
    /// Returns `PceSensitivityError` for mismatched dimensions/counts,
    /// duplicate basis terms, non-finite coefficients, or zero variance.
    pub fn sobol_indices(&self) -> Result<SobolIndices, PceSensitivityError> {
        if self.indices.len() != self.coefficients.len() {
            return Err(PceSensitivityError::CoefficientCount {
                terms: self.indices.len(),
                coefficients: self.coefficients.len(),
            });
        }
        if self.indices.is_empty() {
            return Err(PceSensitivityError::EmptyModel);
        }
        let mut ordered = BTreeMap::<&[usize], usize>::new();
        let mut scale = 0.0f64;
        for (term, (alpha, coefficient)) in self.indices.iter()
            .zip(&self.coefficients).enumerate()
        {
            if alpha.len() != self.dim {
                return Err(PceSensitivityError::TermDimension {
                    term, expected: self.dim, actual: alpha.len(),
                });
            }
            if !coefficient.is_finite() {
                return Err(PceSensitivityError::NonFiniteCoefficient { term });
            }
            if let Some(first) = ordered.insert(alpha.as_slice(), term) {
                return Err(PceSensitivityError::DuplicateBasis { first, second: term });
            }
            if alpha.iter().any(|&degree| degree != 0) {
                scale = scale.max(coefficient.abs());
            }
        }
        if scale == 0.0 {
            return Err(PceSensitivityError::ZeroVariance);
        }
        let mut groups = BTreeMap::<Vec<usize>, EnergySum>::new();
        for (alpha, term) in ordered {
            let coefficient = self.coefficients[term];
            if coefficient == 0.0 {
                continue;
            }
            let support: Vec<_> = alpha.iter().enumerate()
                .filter_map(|(i, &degree)| (degree != 0).then_some(i)).collect();
            if support.is_empty() {
                continue;
            }
            let scaled = coefficient / scale;
            groups.entry(support).or_default().add(scaled * scaled);
        }
        // At least one scaled coefficient is exactly +/-1, so this total is
        // positive. Each squared term is <= 1; no representable in-memory
        // model can make their sum overflow f64.
        let mut variance = EnergySum::default();
        for energy in groups.values() {
            variance.add(energy.sum);
        }
        let mut first_order = vec![0.0; self.dim];
        let mut totals = vec![EnergySum::default(); self.dim];
        let mut components = Vec::with_capacity(groups.len());
        for (variables, energy) in groups {
            let index = energy.sum / variance.sum;
            if variables.len() == 1 {
                first_order[variables[0]] = index;
            }
            for &variable in &variables {
                totals[variable].add(energy.sum);
            }
            components.push(SobolComponent { variables, index });
        }
        Ok(SobolIndices {
            first_order,
            total_order: totals.into_iter().map(|energy| energy.sum / variance.sum).collect(),
            components,
        })
    }
}
