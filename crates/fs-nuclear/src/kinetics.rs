//! Frozen-coefficient point-neutron kinetics, distinct from decay heat.
//!
//! The neutron amplitude `n` and effective precursor coordinates `C_i` obey
//! `n' = (rho-beta)n/Lambda + sum(lambda_i C_i) + S` and
//! `C_i' = beta_i n/Lambda - lambda_i C_i`. Their common normalization is
//! supplied by the caller. Neither an amplitude nor its time integral is a
//! fission count without a separately admitted flux/rate normalization.
//!
//! A positive shifted Taylor exponential followed by scaling and squaring
//! integrates each constant-reactivity/source interval. Two augmented states
//! integrate the external source and neutron amplitude in that same solve;
//! the returned integral is not an endpoint-power rectangle. Small delayed
//! group sets use O(g^2) storage and O(g^3 (terms + squarings)) arithmetic.
//! Stiffness increases the scaling count logarithmically, not through explicit
//! prompt-lifetime substeps. No benchmark or cross-ISA replay claim is made.
//!
//! Series truncation and floating-point operation estimates are budgeted.
//! These are numerical estimates, not interval certificates, data uncertainty
//! bounds, or proof of a spatial/thermal reactor model. See Argonne
//! SAS4A/SASSYS-1, sections 4.2 and 4.3, for the point-kinetics equations.

use std::sync::Arc;

/// One effective delayed-neutron group, not a delayed-energy reservoir.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DelayedNeutronGroup {
    /// Dimensionless effective fraction beta_i, in [0, 1).
    pub delayed_fraction: f64,
    /// Positive precursor decay constant [1/s].
    pub decay_constant_s_inv: f64,
}

/// Immutable caller-identified parameters for a fixed spatial/spectral shape.
#[derive(Clone, Debug, PartialEq)]
pub struct PointKineticsModel {
    data_id: String,
    generation_time_s: f64,
    groups: Vec<DelayedNeutronGroup>,
    delayed_fraction: f64,
}

impl PointKineticsModel {
    /// Admit explicit effective neutron data. No isotope constants, initial
    /// precursor equilibrium, or criticality calculation are supplied here.
    pub fn new(
        data_id: String,
        generation_time_s: f64,
        groups: Vec<DelayedNeutronGroup>,
    ) -> Result<Self, KineticsError> {
        if data_id.trim().is_empty() || !positive(generation_time_s) {
            return Err(KineticsError::InvalidInput(
                "data identity and positive generation time required",
            ));
        }
        if groups.iter().any(|group| {
            !nonnegative(group.delayed_fraction) || !positive(group.decay_constant_s_inv)
        }) {
            return Err(KineticsError::InvalidInput("invalid delayed-neutron group"));
        }
        let delayed_fraction = checked_sum(groups.iter().map(|group| group.delayed_fraction))?;
        if delayed_fraction >= 1.0 {
            return Err(KineticsError::InvalidInput(
                "total delayed fraction must be less than one",
            ));
        }
        Ok(Self {
            data_id,
            generation_time_s,
            groups,
            delayed_fraction,
        })
    }

    /// Original data/revision/applicability label; this label is not validation.
    #[must_use]
    pub fn data_id(&self) -> &str {
        &self.data_id
    }

    /// Effective prompt neutron generation time [s].
    #[must_use]
    pub const fn generation_time_s(&self) -> f64 {
        self.generation_time_s
    }

    /// Ordered, immutable neutron precursor data.
    #[must_use]
    pub fn groups(&self) -> &[DelayedNeutronGroup] {
        &self.groups
    }

    /// Sum of the effective delayed fractions.
    #[must_use]
    pub const fn delayed_fraction(&self) -> f64 {
        self.delayed_fraction
    }

    /// Explicitly construct the precursor coordinates for constant `n`.
    /// Calling this is an authored equilibrium-history assumption; `new` on
    /// the state always requires its own complete initial precursor vector.
    pub fn equilibrium_precursors(
        &self,
        neutron_population: f64,
    ) -> Result<Vec<f64>, KineticsError> {
        if !nonnegative(neutron_population) {
            return Err(KineticsError::InvalidInput("invalid neutron population"));
        }
        self.groups
            .iter()
            .map(|group| {
                finite(
                    (group.delayed_fraction / self.generation_time_s) / group.decay_constant_s_inv
                        * neutron_population,
                )
            })
            .collect()
    }
}

/// Work limits and absolute numerical tolerances for one accepted interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct KineticsBudget {
    /// Maximum number of delayed-neutron groups; dense dimension is groups + 3.
    pub max_delayed_groups: usize,
    /// Maximum matrix squarings after scaling the positive operator.
    pub max_squarings: u32,
    /// Maximum nonconstant terms in the initial Taylor exponential.
    pub max_taylor_terms: usize,
    /// Allowed estimated L1 error of the population/precursor vector and
    /// allowed absolute population-balance residual, in the caller's units.
    pub population_tolerance: f64,
    /// Allowed estimated error in integral(n dt), in population units times s.
    pub integral_tolerance_population_s: f64,
}

impl KineticsBudget {
    fn validate(self) -> Result<(), KineticsError> {
        if !nonnegative(self.population_tolerance)
            || !nonnegative(self.integral_tolerance_population_s)
        {
            return Err(KineticsError::InvalidInput(
                "kinetics tolerances must be finite and nonnegative",
            ));
        }
        Ok(())
    }
}

/// Accepted neutron amplitude and effective precursor history.
/// Clone is an independent checkpoint; only immutable model data are shared.
#[derive(Clone, Debug, PartialEq)]
pub struct PointKineticsState {
    model: Arc<PointKineticsModel>,
    neutron_population: f64,
    precursor_populations: Vec<f64>,
}

impl PointKineticsState {
    /// Supply the actual initial history, in one consistent normalization.
    pub fn new(
        model: Arc<PointKineticsModel>,
        neutron_population: f64,
        precursor_populations: Vec<f64>,
    ) -> Result<Self, KineticsError> {
        if !nonnegative(neutron_population)
            || precursor_populations.len() != model.groups.len()
            || !precursor_populations.iter().copied().all(nonnegative)
        {
            return Err(KineticsError::InvalidInput(
                "invalid neutron or precursor inventory",
            ));
        }
        checked_sum(
            std::iter::once(neutron_population).chain(precursor_populations.iter().copied()),
        )?;
        Ok(Self {
            model,
            neutron_population,
            precursor_populations,
        })
    }

    /// Immutable source parameters.
    #[must_use]
    pub fn model(&self) -> &PointKineticsModel {
        &self.model
    }

    /// Current neutron amplitude, in the caller's fixed normalization.
    #[must_use]
    pub const fn neutron_population(&self) -> f64 {
        self.neutron_population
    }

    /// Current effective precursor coordinates in the original group order.
    #[must_use]
    pub fn precursor_populations(&self) -> &[f64] {
        &self.precursor_populations
    }

    /// Prepare a complete independent candidate. Reactivity rho is finite and
    /// less than one; source is nonnegative population units per second.
    /// Coefficients are frozen over this interval. Coupled feedback and time
    /// variation require explicitly partitioned, independently refined steps.
    pub fn trial(
        &self,
        duration_s: f64,
        reactivity: f64,
        external_source_per_s: f64,
        budget: KineticsBudget,
        mut checkpoint: impl FnMut() -> Result<(), KineticsError>,
    ) -> Result<(Self, KineticsStep), KineticsError> {
        checkpoint()?;
        budget.validate()?;
        if !nonnegative(duration_s)
            || !reactivity.is_finite()
            || reactivity >= 1.0
            || !nonnegative(external_source_per_s)
        {
            return Err(KineticsError::InvalidInput(
                "invalid duration, reactivity, or neutron source",
            ));
        }
        if self.model.groups.len() > budget.max_delayed_groups {
            return Err(KineticsError::Budget("delayed-neutron groups"));
        }
        let initial_population = checked_sum(
            std::iter::once(self.neutron_population)
                .chain(self.precursor_populations.iter().copied()),
        )?;
        let external_neutrons = finite(external_source_per_s * duration_s)?;
        if duration_s > 0.0 && external_source_per_s > 0.0 && external_neutrons == 0.0 {
            return Err(KineticsError::Arithmetic);
        }
        if duration_s == 0.0 || (initial_population == 0.0 && external_neutrons == 0.0) {
            checkpoint()?;
            return Ok((
                self.clone(),
                KineticsStep {
                    duration_s,
                    ..KineticsStep::default()
                },
            ));
        }
        if budget.max_taylor_terms == 0 {
            return Err(KineticsError::Budget("Taylor terms"));
        }

        // In normalized time u=t/dt the source coordinate is S*dt, and the
        // last coordinate is integral(n du). Thus neither source magnitude
        // nor the physical units of the requested integral inflate scaling.
        let dim = self
            .model
            .groups
            .len()
            .checked_add(3)
            .ok_or(KineticsError::Arithmetic)?;
        let source_index = dim - 2;
        let integral_index = dim - 1;
        let mut matrix = zero_matrix(dim)?;
        let prompt_rate =
            finite((reactivity - self.model.delayed_fraction) / self.model.generation_time_s)?;
        matrix[0] = finite(prompt_rate * duration_s)?;
        let mut shift = (-matrix[0]).max(0.0);
        for (index, group) in self.model.groups.iter().enumerate() {
            checkpoint()?;
            let row = index + 1;
            let decay = finite(group.decay_constant_s_inv * duration_s)?;
            let production =
                finite((group.delayed_fraction / self.model.generation_time_s) * duration_s)?;
            matrix[row] = decay;
            matrix[row * dim] = production;
            matrix[row * dim + row] = -decay;
            shift = shift.max(decay);
        }
        matrix[source_index] = 1.0;
        matrix[integral_index * dim] = 1.0;
        let physical_norm = matrix_one_norm(&matrix, dim)?;
        for index in 0..dim {
            matrix[index * dim + index] = finite(matrix[index * dim + index] + shift)?;
        }
        let shifted_norm = matrix_one_norm(&matrix, dim)?;
        let mut scaled_norm = shifted_norm;
        let mut scale = 1.0;
        let mut squarings = 0;
        while scaled_norm > 0.5 {
            checkpoint()?;
            if squarings >= budget.max_squarings {
                return Err(KineticsError::Budget("matrix squarings"));
            }
            scaled_norm *= 0.5;
            scale *= 0.5;
            squarings += 1;
        }
        if scale == 0.0 {
            return Err(KineticsError::Arithmetic);
        }
        for entry in &mut matrix {
            *entry *= scale;
        }
        let shifted_scale = finite(shift * scale)?;
        let unshift = (-shifted_scale).exp();

        // The Taylor matrix is nonnegative. After term k, successive norm
        // ratios are at most ||B||/(k+1), giving a geometric series-tail bound.
        let mut term = identity_matrix(dim)?;
        let mut exponential = term.clone();
        let mut tail = f64::INFINITY;
        let mut term_roundoff = 0.0;
        let mut sum_roundoff = 0.0;
        let mut terms = 0;
        for order in 1..=budget.max_taylor_terms {
            checkpoint()?;
            let previous_norm = matrix_one_norm(&term, dim)?;
            term = matrix_product(&term, &matrix, dim, &mut checkpoint)?;
            let reciprocal = 1.0 / order as f64;
            for entry in &mut term {
                *entry *= reciprocal;
            }
            let term_norm = matrix_one_norm(&term, dim)?;
            term_roundoff = finite(
                term_roundoff * scaled_norm * reciprocal
                    + gamma(2 * dim + 2) * previous_norm * scaled_norm * reciprocal,
            )?;
            for (value, addition) in exponential.iter_mut().zip(&term) {
                *value = finite(*value + addition)?;
            }
            let exponential_norm = matrix_one_norm(&exponential, dim)?;
            sum_roundoff = finite(sum_roundoff + term_roundoff + f64::EPSILON * exponential_norm)?;
            let ratio = scaled_norm / (order as f64 + 1.0);
            tail = finite(term_norm * ratio / (1.0 - ratio))?;
            terms = order;
            if tail <= f64::EPSILON * exponential_norm * 0.125 {
                break;
            }
        }
        // Include coefficient formation, summation and exp/scale roundoff as
        // an estimate. libm is not interval-enclosed, so this is not a proof.
        let mut matrix_error = finite(
            (tail + sum_roundoff) * unshift
                + 8.0 * f64::EPSILON * (physical_norm + shifted_norm) * scale,
        )?;
        for entry in &mut exponential {
            *entry = finite(*entry * unshift)?;
        }
        matrix_error =
            finite(matrix_error + 2.0 * f64::EPSILON * matrix_one_norm(&exponential, dim)?)?;
        for _ in 0..squarings {
            checkpoint()?;
            let norm = matrix_one_norm(&exponential, dim)?;
            let squared = matrix_product(&exponential, &exponential, dim, &mut checkpoint)?;
            matrix_error =
                finite(matrix_error * (2.0 * norm + matrix_error) + gamma(2 * dim) * norm * norm)?;
            exponential = squared;
        }

        let mut initial = Vec::with_capacity(dim);
        initial.push(self.neutron_population);
        initial.extend_from_slice(&self.precursor_populations);
        initial.push(external_neutrons);
        initial.push(0.0);
        let initial_norm = finite(initial_population + external_neutrons)?;
        let operator_norm = matrix_one_norm(&exponential, dim)?;
        let mut final_values = vec![0.0; dim];
        for row in 0..dim {
            checkpoint()?;
            final_values[row] = checked_sum(
                (0..dim).map(|column| exponential[row * dim + column] * initial[column]),
            )?;
        }
        let integral = finite(final_values[integral_index] * duration_s)?;
        let population_error = finite(
            matrix_error * initial_norm
                + gamma(2 * dim) * operator_norm * initial_norm
                + f64::EPSILON * external_neutrons * operator_norm
                + f64::MIN_POSITIVE * (dim * dim) as f64,
        )?;
        let integral_error = finite(population_error * duration_s + f64::EPSILON * integral)?;
        let reactivity_production = finite((reactivity / self.model.generation_time_s) * integral)?;
        let balance_residual = checked_sum(
            std::iter::once(final_values[0] - self.neutron_population)
                .chain(
                    final_values[1..source_index]
                        .iter()
                        .zip(&self.precursor_populations)
                        .map(|(&end, &start)| end - start),
                )
                .chain([-reactivity_production, -external_neutrons]),
        )?;
        if population_error > budget.population_tolerance
            || balance_residual.abs() > budget.population_tolerance
        {
            return Err(KineticsError::Budget("population accuracy"));
        }
        if integral_error > budget.integral_tolerance_population_s {
            return Err(KineticsError::Budget("integrated population accuracy"));
        }
        let next = Self::new(
            Arc::clone(&self.model),
            final_values[0],
            final_values[1..source_index].to_vec(),
        )?;
        checkpoint()?;
        Ok((
            next,
            KineticsStep {
                duration_s,
                integrated_population_s: integral,
                external_neutrons,
                population_balance_residual: balance_residual,
                population_error_estimate: population_error,
                integral_error_estimate_population_s: integral_error,
                taylor_terms: terms,
                squarings,
            },
        ))
    }

    /// Publish only a fully accepted candidate; any refusal leaves self intact.
    pub fn advance(
        &mut self,
        duration_s: f64,
        reactivity: f64,
        external_source_per_s: f64,
        budget: KineticsBudget,
        checkpoint: impl FnMut() -> Result<(), KineticsError>,
    ) -> Result<KineticsStep, KineticsError> {
        let (next, report) = self.trial(
            duration_s,
            reactivity,
            external_source_per_s,
            budget,
            checkpoint,
        )?;
        *self = next;
        Ok(report)
    }
}

/// Integrated source response and numerical work of an accepted interval.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct KineticsStep {
    /// Accepted interval duration [s].
    pub duration_s: f64,
    /// Integral of neutron amplitude, in population units times seconds.
    pub integrated_population_s: f64,
    /// Authored external source integrated over the interval, S*dt.
    pub external_neutrons: f64,
    /// delta(n + sum C_i) - rho/Lambda * integral(n dt) - S*dt.
    pub population_balance_residual: f64,
    /// Absolute population-vector L1 numerical error estimate, not a certificate.
    pub population_error_estimate: f64,
    /// Absolute numerical error estimate of the integral, in population*s.
    pub integral_error_estimate_population_s: f64,
    /// Terms used in the small scaled exponential, excluding the identity.
    pub taylor_terms: usize,
    /// Matrix squarings, bounded independently of the Taylor terms.
    pub squarings: u32,
}

/// Explicit temperature-feedback fit; coefficients carry no automatic isotope
/// or reactor authority. Logarithmic fits can represent a supplied Doppler
/// model, while linear fits represent a locally calibrated derivative.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum TemperatureFeedbackLaw {
    /// rho(T)-rho(Tref) = coefficient * (T-Tref).
    Linear { coefficient_per_kelvin: f64 },
    /// rho(T)-rho(Tref) = coefficient * ln(T/Tref).
    Logarithmic { coefficient: f64 },
}

/// One identified feedback law with an admitted absolute-temperature domain.
#[derive(Clone, Debug, PartialEq)]
pub struct TemperatureFeedback {
    data_id: String,
    reference_temperature_k: f64,
    minimum_temperature_k: f64,
    maximum_temperature_k: f64,
    law: TemperatureFeedbackLaw,
}

impl TemperatureFeedback {
    /// Admit a fit and its domain. The reference temperature must be in-domain.
    pub fn new(
        data_id: String,
        reference_temperature_k: f64,
        minimum_temperature_k: f64,
        maximum_temperature_k: f64,
        law: TemperatureFeedbackLaw,
    ) -> Result<Self, KineticsError> {
        let coefficient = match law {
            TemperatureFeedbackLaw::Linear {
                coefficient_per_kelvin,
            } => coefficient_per_kelvin,
            TemperatureFeedbackLaw::Logarithmic { coefficient } => coefficient,
        };
        if data_id.trim().is_empty()
            || !coefficient.is_finite()
            || !positive(minimum_temperature_k)
            || !positive(maximum_temperature_k)
            || !reference_temperature_k.is_finite()
            || reference_temperature_k < minimum_temperature_k
            || reference_temperature_k > maximum_temperature_k
        {
            return Err(KineticsError::InvalidInput(
                "invalid temperature-feedback fit or domain",
            ));
        }
        Ok(Self {
            data_id,
            reference_temperature_k,
            minimum_temperature_k,
            maximum_temperature_k,
            law,
        })
    }

    /// Caller-provided identity of this feedback fit.
    #[must_use]
    pub fn data_id(&self) -> &str {
        &self.data_id
    }

    /// Evaluate only within the admitted domain. Add the returned shift to
    /// independently authored control/reference reactivity before `trial`.
    /// This evaluator does not itself close a thermal feedback iteration.
    pub fn reactivity_shift(&self, temperature_k: f64) -> Result<f64, KineticsError> {
        if !temperature_k.is_finite()
            || temperature_k < self.minimum_temperature_k
            || temperature_k > self.maximum_temperature_k
        {
            return Err(KineticsError::InvalidInput(
                "temperature is outside feedback-fit domain",
            ));
        }
        match self.law {
            TemperatureFeedbackLaw::Linear {
                coefficient_per_kelvin,
            } => finite(coefficient_per_kelvin * (temperature_k - self.reference_temperature_k)),
            TemperatureFeedbackLaw::Logarithmic { coefficient } => {
                let relative =
                    (temperature_k - self.reference_temperature_k) / self.reference_temperature_k;
                let logarithm = if relative.is_finite() && relative > -1.0 {
                    relative.ln_1p()
                } else {
                    temperature_k.ln() - self.reference_temperature_k.ln()
                };
                finite(coefficient * logarithm)
            }
        }
    }
}

/// Refusal without publication of a neutron history.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum KineticsError {
    /// Invalid model, state, feedback, interval or accuracy input.
    InvalidInput(&'static str),
    /// Caller-supplied work or numerical accuracy limit exceeded.
    Budget(&'static str),
    /// Required arithmetic or storage cannot be represented.
    Arithmetic,
    /// Caller checkpoint requested cancellation.
    Cancelled,
}

impl core::fmt::Display for KineticsError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidInput(detail) => write!(formatter, "invalid kinetics input: {detail}"),
            Self::Budget(detail) => write!(formatter, "kinetics budget exceeded: {detail}"),
            Self::Arithmetic => {
                formatter.write_str("kinetics arithmetic or storage is not representable")
            }
            Self::Cancelled => formatter.write_str("kinetics advance cancelled"),
        }
    }
}
impl std::error::Error for KineticsError {}

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn nonnegative(value: f64) -> bool {
    value.is_finite() && value >= 0.0
}
fn finite(value: f64) -> Result<f64, KineticsError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(KineticsError::Arithmetic)
    }
}

fn checked_sum(values: impl IntoIterator<Item = f64>) -> Result<f64, KineticsError> {
    let (mut total, mut correction) = (0.0_f64, 0.0_f64);
    for value in values {
        let next = finite(total + value)?;
        correction += if total.abs() >= value.abs() {
            (total - next) + value
        } else {
            (value - next) + total
        };
        total = next;
    }
    finite(total + correction)
}

fn zero_matrix(dim: usize) -> Result<Vec<f64>, KineticsError> {
    let count = dim.checked_mul(dim).ok_or(KineticsError::Arithmetic)?;
    let mut matrix = Vec::new();
    matrix
        .try_reserve_exact(count)
        .map_err(|_| KineticsError::Budget("matrix memory"))?;
    matrix.resize(count, 0.0);
    Ok(matrix)
}

fn identity_matrix(dim: usize) -> Result<Vec<f64>, KineticsError> {
    let mut matrix = zero_matrix(dim)?;
    for index in 0..dim {
        matrix[index * dim + index] = 1.0;
    }
    Ok(matrix)
}

fn matrix_one_norm(matrix: &[f64], dim: usize) -> Result<f64, KineticsError> {
    let mut norm = 0.0_f64;
    for column in 0..dim {
        norm = norm.max(checked_sum(
            (0..dim).map(|row| matrix[row * dim + column].abs()),
        )?);
    }
    Ok(norm)
}

fn matrix_product(
    left: &[f64],
    right: &[f64],
    dim: usize,
    checkpoint: &mut impl FnMut() -> Result<(), KineticsError>,
) -> Result<Vec<f64>, KineticsError> {
    let mut product = zero_matrix(dim)?;
    for row in 0..dim {
        checkpoint()?;
        // All entries are nonnegative: no cancellation, pivoting or post-hoc
        // clipping can turn a negative population into an accepted candidate.
        for inner in 0..dim {
            let value = left[row * dim + inner];
            if value == 0.0 {
                continue;
            }
            for column in 0..dim {
                product[row * dim + column] += value * right[inner * dim + column];
            }
        }
        if !product[row * dim..(row + 1) * dim]
            .iter()
            .copied()
            .all(nonnegative)
        {
            return Err(KineticsError::Arithmetic);
        }
    }
    Ok(product)
}

fn gamma(operations: usize) -> f64 {
    let scaled = operations as f64 * f64::EPSILON;
    scaled / (1.0 - scaled)
}

#[cfg(test)]
mod tests;
