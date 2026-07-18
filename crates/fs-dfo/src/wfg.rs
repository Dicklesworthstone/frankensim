//! Typed kernels for the Walking Fish Group benchmark family.
//!
//! The production kernel currently implements normalized WFG1 through WFG4
//! compositions for arbitrary valid objective, position-parameter, and
//! distance-parameter counts.  Its transformations and shapes follow the
//! corrected WFG toolkit as represented by jMetal revision
//! `ea7e882f6b8f94b99535921674e62cda7986f20e`.  Inputs are already normalized
//! to `[0, 1]`; accepting the heterogeneous canonical bounds
//! `z_i in [0, 2(i + 1)]` is deliberately left to a later adapter.
//!
//! Determinism is structural within the evaluator: reductions have a fixed
//! left-to-right order and transcendental calls use [`fs_math::det`].  This
//! module does not yet claim the complete WFG5-WFG9 suite, executable parity
//! with an external oracle, optimizer convergence, cancellation coverage,
//! cross-ISA bit stability, or performance evidence.

#![deny(unsafe_code)]

const CORRECTION_EPSILON: f64 = 1.0e-10;
const S_MULTI_A: f64 = 30.0;
const S_MULTI_B: f64 = 10.0;
const S_MULTI_CENTER: f64 = 0.35;

/// Structured refusal from a typed normalized WFG evaluator.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WfgError {
    /// A multiobjective WFG problem needs at least two objectives.
    TooFewObjectives {
        /// Supplied objective count.
        objectives: usize,
    },
    /// Each of the `M - 1` position groups must contain a parameter.
    TooFewPositionParameters {
        /// Supplied position-parameter count.
        position_parameters: usize,
        /// Supplied objective count.
        objectives: usize,
    },
    /// Position parameters cannot be divided into equal objective groups.
    PositionParametersNotDivisible {
        /// Supplied position-parameter count.
        position_parameters: usize,
        /// Required group count, equal to `objectives - 1`.
        groups: usize,
    },
    /// The WFG distance block must not be empty.
    NoDistanceParameters,
    /// WFG2 and WFG3 reduce distance coordinates in adjacent pairs.
    DistanceParametersNotEven {
        /// Supplied distance-parameter count.
        distance_parameters: usize,
    },
    /// The total decision dimension could not be represented by `usize`.
    DimensionOverflow {
        /// Supplied position-parameter count.
        position_parameters: usize,
        /// Supplied distance-parameter count.
        distance_parameters: usize,
    },
    /// The normalized decision vector has the wrong dimension.
    WrongInputLength {
        /// Dimension admitted by the problem specification.
        expected: usize,
        /// Supplied decision-vector length.
        actual: usize,
    },
    /// A normalized coordinate is NaN or infinite.
    NonFiniteInput {
        /// Zero-based coordinate index.
        component: usize,
        /// Exact IEEE-754 payload.
        bits: u64,
    },
    /// A finite normalized coordinate lies outside `[0, 1]`.
    InputOutOfRange {
        /// Zero-based coordinate index.
        component: usize,
        /// Exact IEEE-754 payload.
        bits: u64,
    },
    /// Storage for a validated evaluation could not be reserved.
    AllocationFailed {
        /// Stable name of the requested vector.
        what: &'static str,
        /// Number of `f64` elements requested.
        elements: usize,
    },
}

impl core::fmt::Display for WfgError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match *self {
            Self::TooFewObjectives { objectives } => write!(
                formatter,
                "WFG requires at least two objectives, received {objectives}"
            ),
            Self::TooFewPositionParameters {
                position_parameters,
                objectives,
            } => write!(
                formatter,
                "WFG with {objectives} objectives requires at least {} position parameters, received {position_parameters}",
                objectives.saturating_sub(1)
            ),
            Self::PositionParametersNotDivisible {
                position_parameters,
                groups,
            } => write!(
                formatter,
                "WFG position count {position_parameters} is not divisible by its {groups} objective groups"
            ),
            Self::NoDistanceParameters => {
                formatter.write_str("WFG requires at least one distance parameter")
            }
            Self::DistanceParametersNotEven {
                distance_parameters,
            } => write!(
                formatter,
                "WFG2 and WFG3 require an even distance count, received {distance_parameters}"
            ),
            Self::DimensionOverflow {
                position_parameters,
                distance_parameters,
            } => write!(
                formatter,
                "WFG decision dimension {position_parameters} + {distance_parameters} overflowed usize"
            ),
            Self::WrongInputLength { expected, actual } => write!(
                formatter,
                "WFG expected {expected} normalized coordinates, received {actual}"
            ),
            Self::NonFiniteInput { component, bits } => write!(
                formatter,
                "WFG normalized coordinate {component} is non-finite (bits 0x{bits:016x})"
            ),
            Self::InputOutOfRange { component, bits } => write!(
                formatter,
                "WFG normalized coordinate {component} lies outside [0, 1] (bits 0x{bits:016x})"
            ),
            Self::AllocationFailed { what, elements } => write!(
                formatter,
                "WFG could not reserve {elements} elements for {what}"
            ),
        }
    }
}

impl std::error::Error for WfgError {}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct WfgSpec {
    objectives: usize,
    position_parameters: usize,
    distance_parameters: usize,
    dimension: usize,
    position_group_size: usize,
}

impl WfgSpec {
    fn new(
        objectives: usize,
        position_parameters: usize,
        distance_parameters: usize,
    ) -> Result<Self, WfgError> {
        if objectives < 2 {
            return Err(WfgError::TooFewObjectives { objectives });
        }
        let groups = objectives - 1;
        if position_parameters < groups {
            return Err(WfgError::TooFewPositionParameters {
                position_parameters,
                objectives,
            });
        }
        if !position_parameters.is_multiple_of(groups) {
            return Err(WfgError::PositionParametersNotDivisible {
                position_parameters,
                groups,
            });
        }
        if distance_parameters == 0 {
            return Err(WfgError::NoDistanceParameters);
        }
        let Some(dimension) = position_parameters.checked_add(distance_parameters) else {
            return Err(WfgError::DimensionOverflow {
                position_parameters,
                distance_parameters,
            });
        };

        Ok(Self {
            objectives,
            position_parameters,
            distance_parameters,
            dimension,
            position_group_size: position_parameters / groups,
        })
    }

    fn validate_input(self, input: &[f64]) -> Result<(), WfgError> {
        if input.len() != self.dimension {
            return Err(WfgError::WrongInputLength {
                expected: self.dimension,
                actual: input.len(),
            });
        }
        for (component, &value) in input.iter().enumerate() {
            if !value.is_finite() {
                return Err(WfgError::NonFiniteInput {
                    component,
                    bits: value.to_bits(),
                });
            }
            if !(0.0..=1.0).contains(&value) {
                return Err(WfgError::InputOutOfRange {
                    component,
                    bits: value.to_bits(),
                });
            }
        }
        Ok(())
    }
}

macro_rules! wfg_accessors {
    () => {
        /// Number of objective values produced by each evaluation.
        #[must_use]
        pub const fn objectives(self) -> usize {
            self.spec.objectives
        }

        /// Number of position-related decision parameters.
        #[must_use]
        pub const fn position_parameters(self) -> usize {
            self.spec.position_parameters
        }

        /// Number of distance-related decision parameters.
        #[must_use]
        pub const fn distance_parameters(self) -> usize {
            self.spec.distance_parameters
        }

        /// Total normalized decision dimension, `k + l`.
        #[must_use]
        pub const fn dimension(self) -> usize {
            self.spec.dimension
        }
    };
}

/// A validated normalized WFG1 problem definition.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Wfg1 {
    spec: WfgSpec,
}

impl Wfg1 {
    /// Validate a normalized WFG1 problem definition.
    ///
    /// # Errors
    ///
    /// Returns a structured refusal for an invalid objective count, position
    /// partition, distance count, or checked total dimension.
    pub fn new(
        objectives: usize,
        position_parameters: usize,
        distance_parameters: usize,
    ) -> Result<Self, WfgError> {
        Ok(Self {
            spec: WfgSpec::new(objectives, position_parameters, distance_parameters)?,
        })
    }

    wfg_accessors!();

    /// Evaluate one normalized decision vector.
    ///
    /// # Errors
    ///
    /// Returns a structured refusal for invalid input or failed intermediate
    /// storage admission.  No transformation runs before input validation.
    pub fn evaluate_normalized(&self, input: &[f64]) -> Result<WfgEvaluation, WfgError> {
        self.spec.validate_input(input)?;
        let transformed = wfg1_transform(input, self.spec.position_parameters)?;
        let reduced = wfg1_reduce(&transformed, self.spec)?;
        let positioned = identity_positioned(&reduced)?;
        let shape = wfg1_shape(&positioned)?;
        finish_evaluation(transformed, reduced, positioned, shape)
    }
}

/// A validated normalized WFG2 problem definition.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Wfg2 {
    spec: WfgSpec,
}

impl Wfg2 {
    /// Validate a normalized WFG2 problem definition.
    ///
    /// # Errors
    ///
    /// Returns a structured refusal for invalid shared dimensions or an odd
    /// distance count, which cannot be reduced in canonical adjacent pairs.
    pub fn new(
        objectives: usize,
        position_parameters: usize,
        distance_parameters: usize,
    ) -> Result<Self, WfgError> {
        Ok(Self {
            spec: even_distance_spec(objectives, position_parameters, distance_parameters)?,
        })
    }

    wfg_accessors!();

    /// Evaluate one normalized decision vector.
    ///
    /// # Errors
    ///
    /// Returns a structured refusal for invalid input or failed intermediate
    /// storage admission.  No transformation runs before input validation.
    pub fn evaluate_normalized(&self, input: &[f64]) -> Result<WfgEvaluation, WfgError> {
        self.spec.validate_input(input)?;
        let transformed = wfg23_transform(input, self.spec)?;
        let reduced = equal_group_reduce(&transformed, self.spec)?;
        let positioned = identity_positioned(&reduced)?;
        let shape = wfg2_shape(&positioned)?;
        finish_evaluation(transformed, reduced, positioned, shape)
    }
}

/// A validated normalized WFG3 problem definition.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Wfg3 {
    spec: WfgSpec,
}

impl Wfg3 {
    /// Validate a normalized WFG3 problem definition.
    ///
    /// # Errors
    ///
    /// Returns a structured refusal for invalid shared dimensions or an odd
    /// distance count, which cannot be reduced in canonical adjacent pairs.
    pub fn new(
        objectives: usize,
        position_parameters: usize,
        distance_parameters: usize,
    ) -> Result<Self, WfgError> {
        Ok(Self {
            spec: even_distance_spec(objectives, position_parameters, distance_parameters)?,
        })
    }

    wfg_accessors!();

    /// Evaluate one normalized decision vector.
    ///
    /// # Errors
    ///
    /// Returns a structured refusal for invalid input or failed intermediate
    /// storage admission.  No transformation runs before input validation.
    pub fn evaluate_normalized(&self, input: &[f64]) -> Result<WfgEvaluation, WfgError> {
        self.spec.validate_input(input)?;
        let transformed = wfg23_transform(input, self.spec)?;
        let reduced = equal_group_reduce(&transformed, self.spec)?;
        let positioned = wfg3_positioned(&reduced)?;
        let shape = linear_shape(&positioned)?;
        finish_evaluation(transformed, reduced, positioned, shape)
    }
}

/// A validated normalized WFG4 problem definition.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Wfg4 {
    spec: WfgSpec,
}

impl Wfg4 {
    /// Validate a normalized WFG4 problem definition.
    ///
    /// # Errors
    ///
    /// Returns a structured refusal for an invalid objective count, position
    /// partition, distance count, or checked total dimension.
    pub fn new(
        objectives: usize,
        position_parameters: usize,
        distance_parameters: usize,
    ) -> Result<Self, WfgError> {
        Ok(Self {
            spec: WfgSpec::new(objectives, position_parameters, distance_parameters)?,
        })
    }

    wfg_accessors!();

    /// Evaluate one normalized decision vector.
    ///
    /// The returned intermediate vectors make benchmark receipts and
    /// independent recomputation possible without duplicating private logic.
    ///
    /// # Errors
    ///
    /// Returns a structured refusal for invalid input or failed intermediate
    /// storage admission.  No transformation runs before input validation.
    pub fn evaluate_normalized(&self, input: &[f64]) -> Result<WfgEvaluation, WfgError> {
        self.spec.validate_input(input)?;
        let transformed = wfg4_transform(input)?;
        let reduced = equal_group_reduce(&transformed, self.spec)?;
        let positioned = identity_positioned(&reduced)?;
        let shape = concave_shape(&positioned)?;
        finish_evaluation(transformed, reduced, positioned, shape)
    }
}

/// One normalized WFG evaluation with replay-relevant intermediates.
#[derive(Debug, Clone, PartialEq)]
pub struct WfgEvaluation {
    transformed: Vec<f64>,
    reduced: Vec<f64>,
    positioned: Vec<f64>,
    shape: Vec<f64>,
    objectives: Vec<f64>,
}

impl WfgEvaluation {
    /// Final transformation vector immediately before objective reduction.
    #[must_use]
    pub fn transformed(&self) -> &[f64] {
        &self.transformed
    }

    /// The `M` fixed-order objective reductions, conventionally named `t`.
    #[must_use]
    pub fn reduced(&self) -> &[f64] {
        &self.reduced
    }

    /// The `M` post-degeneracy coordinates, conventionally named `x`.
    #[must_use]
    pub fn positioned(&self) -> &[f64] {
        &self.positioned
    }

    /// The `M` WFG shape values before scale and distance are applied.
    #[must_use]
    pub fn shape(&self) -> &[f64] {
        &self.shape
    }

    /// The scaled WFG objective vector.
    #[must_use]
    pub fn objectives(&self) -> &[f64] {
        &self.objectives
    }

    /// Consume the evaluation and return its objective vector.
    #[must_use]
    pub fn into_objectives(self) -> Vec<f64> {
        self.objectives
    }
}

fn reserved_vec(what: &'static str, elements: usize) -> Result<Vec<f64>, WfgError> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(elements)
        .map_err(|_| WfgError::AllocationFailed { what, elements })?;
    Ok(values)
}

fn even_distance_spec(
    objectives: usize,
    position_parameters: usize,
    distance_parameters: usize,
) -> Result<WfgSpec, WfgError> {
    let spec = WfgSpec::new(objectives, position_parameters, distance_parameters)?;
    if !distance_parameters.is_multiple_of(2) {
        return Err(WfgError::DistanceParametersNotEven {
            distance_parameters,
        });
    }
    Ok(spec)
}

fn finish_evaluation(
    transformed: Vec<f64>,
    reduced: Vec<f64>,
    positioned: Vec<f64>,
    shape: Vec<f64>,
) -> Result<WfgEvaluation, WfgError> {
    debug_assert_eq!(reduced.len(), positioned.len());
    debug_assert_eq!(positioned.len(), shape.len());
    let objective_count = positioned.len();
    let distance = positioned[objective_count - 1];
    let mut objectives = reserved_vec("objectives", objective_count)?;
    for (index, &shape_value) in shape.iter().enumerate() {
        let scale = 2.0 * (index + 1) as f64;
        objectives.push(scale.mul_add(shape_value, distance));
    }
    Ok(WfgEvaluation {
        transformed,
        reduced,
        positioned,
        shape,
        objectives,
    })
}

fn identity_positioned(reduced: &[f64]) -> Result<Vec<f64>, WfgError> {
    let mut positioned = reserved_vec("positioned coordinates", reduced.len())?;
    positioned.extend_from_slice(reduced);
    Ok(positioned)
}

fn wfg3_positioned(reduced: &[f64]) -> Result<Vec<f64>, WfgError> {
    let objective_count = reduced.len();
    let distance = reduced[objective_count - 1];
    let mut positioned = reserved_vec("positioned coordinates", objective_count)?;
    for (index, &coordinate) in reduced[..objective_count - 1].iter().enumerate() {
        if index == 0 {
            positioned.push(coordinate);
        } else {
            positioned.push(distance.mul_add(coordinate - 0.5, 0.5));
        }
    }
    positioned.push(distance);
    Ok(positioned)
}

fn wfg1_transform(input: &[f64], position_parameters: usize) -> Result<Vec<f64>, WfgError> {
    let mut transformed = reserved_vec("WFG1 transformed coordinates", input.len())?;
    for (index, &value) in input.iter().enumerate() {
        let shifted = if index < position_parameters {
            value
        } else {
            b_flat(s_linear(value, 0.35), 0.8, 0.75, 0.85)
        };
        transformed.push(b_poly(shifted, 0.02));
    }
    Ok(transformed)
}

fn wfg23_transform(input: &[f64], spec: WfgSpec) -> Result<Vec<f64>, WfgError> {
    let compressed_distance = spec.distance_parameters / 2;
    let transformed_len = spec.position_parameters + compressed_distance;
    let mut transformed = reserved_vec("WFG2/WFG3 transformed coordinates", transformed_len)?;
    transformed.extend_from_slice(&input[..spec.position_parameters]);
    for pair in input[spec.position_parameters..].chunks_exact(2) {
        let shifted = [s_linear(pair[0], 0.35), s_linear(pair[1], 0.35)];
        transformed.push(r_nonsep(&shifted, 2));
    }
    Ok(transformed)
}

fn wfg4_transform(input: &[f64]) -> Result<Vec<f64>, WfgError> {
    let mut transformed = reserved_vec("WFG4 transformed coordinates", input.len())?;
    for &value in input {
        transformed.push(s_multi(value));
    }
    Ok(transformed)
}

fn wfg1_reduce(transformed: &[f64], spec: WfgSpec) -> Result<Vec<f64>, WfgError> {
    let mut reduced = reserved_vec("WFG1 reduced coordinates", spec.objectives)?;
    for group in 0..spec.objectives - 1 {
        let start = group * spec.position_group_size;
        let end = start + spec.position_group_size;
        reduced.push(linearly_weighted_reduction(&transformed[start..end], start));
    }
    reduced.push(linearly_weighted_reduction(
        &transformed[spec.position_parameters..],
        spec.position_parameters,
    ));
    Ok(reduced)
}

fn equal_group_reduce(transformed: &[f64], spec: WfgSpec) -> Result<Vec<f64>, WfgError> {
    let mut reduced = reserved_vec("reduced coordinates", spec.objectives)?;
    for group in 0..spec.objectives - 1 {
        let start = group * spec.position_group_size;
        let end = start + spec.position_group_size;
        reduced.push(equal_weight_reduction(&transformed[start..end]));
    }
    reduced.push(equal_weight_reduction(
        &transformed[spec.position_parameters..],
    ));
    Ok(reduced)
}

fn wfg1_shape(positioned: &[f64]) -> Result<Vec<f64>, WfgError> {
    let mut shape = convex_shape(positioned)?;
    let last = shape.len() - 1;
    shape[last] = mixed_shape(positioned, 5, 1.0);
    Ok(shape)
}

fn wfg2_shape(positioned: &[f64]) -> Result<Vec<f64>, WfgError> {
    let mut shape = convex_shape(positioned)?;
    let last = shape.len() - 1;
    shape[last] = disc_shape(positioned, 5, 1.0, 1.0);
    Ok(shape)
}

/// Snap only roundoff-sized excursions at the WFG unit-interval boundaries.
/// This is not a general clamp.
fn correct_to_01(value: f64) -> f64 {
    if (-CORRECTION_EPSILON..=0.0).contains(&value) {
        0.0
    } else if (1.0..=1.0 + CORRECTION_EPSILON).contains(&value) {
        1.0
    } else {
        value
    }
}

fn b_poly(value: f64, alpha: f64) -> f64 {
    debug_assert!(alpha > 0.0);
    correct_to_01(fs_math::det::pow(value, alpha))
}

fn b_flat(value: f64, height: f64, lower: f64, upper: f64) -> f64 {
    debug_assert!((0.0..=1.0).contains(&height));
    debug_assert!(0.0 < lower && lower < upper && upper < 1.0);
    let below = if value < lower {
        -height * (lower - value) / lower
    } else {
        0.0
    };
    let above = if value > upper {
        -(1.0 - height) * (value - upper) / (1.0 - upper)
    } else {
        0.0
    };
    correct_to_01(height + below - above)
}

fn s_linear(value: f64, zero: f64) -> f64 {
    debug_assert!(0.0 < zero && zero < 1.0);
    let denominator = if value <= zero { zero } else { 1.0 - zero };
    correct_to_01((value - zero).abs() / denominator)
}

fn r_nonsep(values: &[f64], subproblem_size: usize) -> f64 {
    debug_assert!(!values.is_empty());
    debug_assert!(subproblem_size > 0 && subproblem_size <= values.len());
    let mut numerator = 0.0;
    for (index, &value) in values.iter().enumerate() {
        numerator += value;
        for offset in 1..subproblem_size {
            numerator += (value - values[(index + offset) % values.len()]).abs();
        }
    }
    let half_ceiling = subproblem_size.div_ceil(2) as f64;
    let size = subproblem_size as f64;
    let denominator =
        values.len() as f64 * half_ceiling * (1.0 + 2.0 * size - 2.0 * half_ceiling) / size;
    correct_to_01(numerator / denominator)
}

/// Canonical WFG `s_multi(y, A=30, B=10, C=0.35)` transformation.
fn s_multi(value: f64) -> f64 {
    let denominator = if value <= S_MULTI_CENTER {
        2.0 * S_MULTI_CENTER
    } else {
        2.0 * (S_MULTI_CENTER - 1.0)
    };
    let ratio = (value - S_MULTI_CENTER).abs() / denominator;
    let phase = (4.0 * S_MULTI_A + 2.0) * core::f64::consts::PI * (0.5 - ratio);
    let quadratic = 4.0 * S_MULTI_B * ratio * ratio;
    correct_to_01((1.0 + fs_math::det::cos(phase) + quadratic) / (S_MULTI_B + 2.0))
}

fn equal_weight_reduction(values: &[f64]) -> f64 {
    debug_assert!(!values.is_empty());
    correct_to_01(values.iter().sum::<f64>() / values.len() as f64)
}

fn linearly_weighted_reduction(values: &[f64], global_start: usize) -> f64 {
    debug_assert!(!values.is_empty());
    let mut numerator = 0.0;
    let mut denominator = 0.0;
    for (offset, &value) in values.iter().enumerate() {
        let weight = 2.0 * (global_start + offset + 1) as f64;
        numerator = value.mul_add(weight, numerator);
        denominator += weight;
    }
    correct_to_01(numerator / denominator)
}

fn linear_shape(positioned: &[f64]) -> Result<Vec<f64>, WfgError> {
    let objectives = positioned.len();
    let mut prefix = reserved_vec("linear shape prefix", objectives)?;
    let mut product = 1.0;
    prefix.push(product);
    for &coordinate in &positioned[..objectives - 1] {
        product *= coordinate;
        prefix.push(product);
    }

    let mut shape = reserved_vec("linear shape", objectives)?;
    for objective in 1..=objectives {
        let prefix_index = objectives - objective;
        let value = if objective == 1 {
            prefix[prefix_index]
        } else {
            prefix[prefix_index] * (1.0 - positioned[prefix_index])
        };
        shape.push(value);
    }
    Ok(shape)
}

fn convex_shape(positioned: &[f64]) -> Result<Vec<f64>, WfgError> {
    let objectives = positioned.len();
    let mut prefix = reserved_vec("convex shape prefix", objectives)?;
    let mut product = 1.0;
    prefix.push(product);
    for &coordinate in &positioned[..objectives - 1] {
        product *= 1.0 - fs_math::det::cos(coordinate * core::f64::consts::FRAC_PI_2);
        prefix.push(product);
    }

    let mut shape = reserved_vec("convex shape", objectives)?;
    for objective in 1..=objectives {
        let prefix_index = objectives - objective;
        let value = if objective == 1 {
            prefix[prefix_index]
        } else {
            prefix[prefix_index]
                * (1.0 - fs_math::det::sin(positioned[prefix_index] * core::f64::consts::FRAC_PI_2))
        };
        shape.push(value);
    }
    Ok(shape)
}

fn mixed_shape(positioned: &[f64], waves: usize, alpha: f64) -> f64 {
    debug_assert!(waves > 0 && alpha > 0.0);
    let waves = waves as f64;
    let denominator = 2.0 * waves * core::f64::consts::PI;
    let ripple =
        fs_math::det::cos(denominator.mul_add(positioned[0], core::f64::consts::FRAC_PI_2))
            / denominator;
    fs_math::det::pow(1.0 - positioned[0] - ripple, alpha)
}

fn disc_shape(positioned: &[f64], waves: usize, alpha: f64, beta: f64) -> f64 {
    debug_assert!(waves > 0 && alpha > 0.0 && beta > 0.0);
    let powered = fs_math::det::pow(positioned[0], beta);
    let ripple = fs_math::det::cos(waves as f64 * powered * core::f64::consts::PI);
    1.0 - fs_math::det::pow(positioned[0], alpha) * ripple * ripple
}

fn concave_shape(reduced: &[f64]) -> Result<Vec<f64>, WfgError> {
    let objectives = reduced.len();
    let mut sine_prefix = reserved_vec("concave sine prefix", objectives)?;
    let mut product = 1.0;
    sine_prefix.push(product);
    for &coordinate in &reduced[..objectives - 1] {
        product *= fs_math::det::sin(coordinate * core::f64::consts::FRAC_PI_2);
        sine_prefix.push(product);
    }

    let mut shape = reserved_vec("concave shape", objectives)?;
    for objective in 1..=objectives {
        let prefix_index = objectives - objective;
        let value = if objective == 1 {
            sine_prefix[prefix_index]
        } else {
            sine_prefix[prefix_index]
                * fs_math::det::cos(reduced[prefix_index] * core::f64::consts::FRAC_PI_2)
        };
        shape.push(value);
    }
    Ok(shape)
}

#[cfg(test)]
mod tests {
    use super::*;

    const TOLERANCE: f64 = 1.0e-11;

    fn assert_close(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() <= TOLERANCE,
            "actual={actual:.17e}, expected={expected:.17e}"
        );
    }

    fn assert_slice_close(actual: &[f64], expected: &[f64]) {
        assert_eq!(actual.len(), expected.len());
        for (&actual, &expected) in actual.iter().zip(expected) {
            assert_close(actual, expected);
        }
    }

    fn assert_slice_bits_eq(actual: &[f64], expected: &[f64]) {
        assert_eq!(actual.len(), expected.len());
        for (index, (&actual, &expected)) in actual.iter().zip(expected).enumerate() {
            assert_eq!(
                actual.to_bits(),
                expected.to_bits(),
                "bit mismatch at index {index}: actual={actual:.17e}, expected={expected:.17e}"
            );
        }
    }

    fn assert_evaluation_bits_eq(actual: &WfgEvaluation, expected: &WfgEvaluation) {
        assert_slice_bits_eq(actual.transformed(), expected.transformed());
        assert_slice_bits_eq(actual.reduced(), expected.reduced());
        assert_slice_bits_eq(actual.positioned(), expected.positioned());
        assert_slice_bits_eq(actual.shape(), expected.shape());
        assert_slice_bits_eq(actual.objectives(), expected.objectives());
    }

    fn assert_common_input_refusals(evaluate: impl Fn(&[f64]) -> Result<WfgEvaluation, WfgError>) {
        assert_eq!(
            evaluate(&[0.0; 5]).unwrap_err(),
            WfgError::WrongInputLength {
                expected: 6,
                actual: 5,
            }
        );

        for (value, expected) in [
            (
                f64::NAN,
                WfgError::NonFiniteInput {
                    component: 2,
                    bits: f64::NAN.to_bits(),
                },
            ),
            (
                f64::INFINITY,
                WfgError::NonFiniteInput {
                    component: 2,
                    bits: f64::INFINITY.to_bits(),
                },
            ),
            (
                -f64::EPSILON,
                WfgError::InputOutOfRange {
                    component: 2,
                    bits: (-f64::EPSILON).to_bits(),
                },
            ),
            (
                1.0 + f64::EPSILON,
                WfgError::InputOutOfRange {
                    component: 2,
                    bits: (1.0 + f64::EPSILON).to_bits(),
                },
            ),
        ] {
            let mut input = [0.0; 6];
            input[2] = value;
            assert_eq!(evaluate(&input).unwrap_err(), expected);
        }
    }

    #[test]
    fn specification_refuses_malformed_dimensions() {
        for error in [
            Wfg1::new(1, 1, 2).unwrap_err(),
            Wfg2::new(1, 1, 2).unwrap_err(),
            Wfg3::new(1, 1, 2).unwrap_err(),
        ] {
            assert_eq!(error, WfgError::TooFewObjectives { objectives: 1 });
        }
        for error in [
            Wfg1::new(4, 4, 2).unwrap_err(),
            Wfg2::new(4, 4, 2).unwrap_err(),
            Wfg3::new(4, 4, 2).unwrap_err(),
        ] {
            assert_eq!(
                error,
                WfgError::PositionParametersNotDivisible {
                    position_parameters: 4,
                    groups: 3,
                }
            );
        }
        for error in [
            Wfg1::new(3, 4, 0).unwrap_err(),
            Wfg2::new(3, 4, 0).unwrap_err(),
            Wfg3::new(3, 4, 0).unwrap_err(),
        ] {
            assert_eq!(error, WfgError::NoDistanceParameters);
        }

        assert_eq!(
            Wfg4::new(1, 1, 1).unwrap_err(),
            WfgError::TooFewObjectives { objectives: 1 }
        );
        assert_eq!(
            Wfg4::new(4, 2, 1).unwrap_err(),
            WfgError::TooFewPositionParameters {
                position_parameters: 2,
                objectives: 4,
            }
        );
        assert_eq!(
            Wfg4::new(4, 4, 1).unwrap_err(),
            WfgError::PositionParametersNotDivisible {
                position_parameters: 4,
                groups: 3,
            }
        );
        assert_eq!(
            Wfg4::new(3, 4, 0).unwrap_err(),
            WfgError::NoDistanceParameters
        );
        assert_eq!(
            Wfg4::new(2, usize::MAX, 1).unwrap_err(),
            WfgError::DimensionOverflow {
                position_parameters: usize::MAX,
                distance_parameters: 1,
            }
        );
        assert_eq!(
            Wfg2::new(3, 4, 3).unwrap_err(),
            WfgError::DistanceParametersNotEven {
                distance_parameters: 3,
            }
        );
        assert_eq!(
            Wfg3::new(3, 4, 3).unwrap_err(),
            WfgError::DistanceParametersNotEven {
                distance_parameters: 3,
            }
        );
        assert!(Wfg1::new(3, 4, 3).is_ok());
    }

    #[test]
    fn all_variant_specs_expose_their_admitted_dimensions() {
        let wfg1 = Wfg1::new(3, 4, 3).unwrap();
        assert_eq!(wfg1.objectives(), 3);
        assert_eq!(wfg1.position_parameters(), 4);
        assert_eq!(wfg1.distance_parameters(), 3);
        assert_eq!(wfg1.dimension(), 7);

        let wfg2 = Wfg2::new(4, 6, 2).unwrap();
        assert_eq!(wfg2.objectives(), 4);
        assert_eq!(wfg2.position_parameters(), 6);
        assert_eq!(wfg2.distance_parameters(), 2);
        assert_eq!(wfg2.dimension(), 8);

        let wfg3 = Wfg3::new(5, 8, 4).unwrap();
        assert_eq!(wfg3.objectives(), 5);
        assert_eq!(wfg3.position_parameters(), 8);
        assert_eq!(wfg3.distance_parameters(), 4);
        assert_eq!(wfg3.dimension(), 12);
    }

    #[test]
    fn input_admission_is_exact_and_structured() {
        let wfg1 = Wfg1::new(3, 4, 2).unwrap();
        let wfg2 = Wfg2::new(3, 4, 2).unwrap();
        let wfg3 = Wfg3::new(3, 4, 2).unwrap();
        let wfg4 = Wfg4::new(3, 4, 2).unwrap();

        assert_common_input_refusals(|input| wfg1.evaluate_normalized(input));
        assert_common_input_refusals(|input| wfg2.evaluate_normalized(input));
        assert_common_input_refusals(|input| wfg3.evaluate_normalized(input));
        assert_common_input_refusals(|input| wfg4.evaluate_normalized(input));
    }

    #[test]
    fn correction_snaps_only_boundary_roundoff() {
        assert_eq!(correct_to_01(-CORRECTION_EPSILON), 0.0);
        assert_eq!(correct_to_01(-0.0).to_bits(), 0.0f64.to_bits());
        assert_eq!(correct_to_01(0.0), 0.0);
        assert_eq!(correct_to_01(1.0), 1.0);
        assert_eq!(correct_to_01(1.0 + CORRECTION_EPSILON), 1.0);
        assert_eq!(
            correct_to_01(-2.0 * CORRECTION_EPSILON),
            -2.0 * CORRECTION_EPSILON
        );
        assert_eq!(
            correct_to_01(1.0 + 2.0 * CORRECTION_EPSILON),
            1.0 + 2.0 * CORRECTION_EPSILON
        );
    }

    #[test]
    fn canonical_s_multi_anchors_match_corrected_toolkit() {
        assert_close(s_multi(0.0), 1.0);
        assert_close(s_multi(S_MULTI_CENTER), 0.0);
        assert_close(s_multi(1.0), 1.0);
    }

    #[test]
    fn canonical_wfg1_through_wfg3_primitive_anchors_match_the_toolkit() {
        assert_close(b_poly(0.25, 2.0), 0.0625);
        assert_close(b_flat(0.5, 0.8, 0.75, 0.85), 0.533_333_333_333_333_3);
        assert_close(b_flat(0.8, 0.8, 0.75, 0.85), 0.8);
        assert_close(b_flat(0.9, 0.8, 0.75, 0.85), 0.866_666_666_666_666_7);
        assert_close(s_linear(0.0, 0.35), 1.0);
        assert_close(s_linear(0.35, 0.35), 0.0);
        assert_close(s_linear(1.0, 0.35), 1.0);
        assert_close(r_nonsep(&[0.2, 0.8], 2), 0.733_333_333_333_333_4);

        let positioned = [0.2, 0.4, 0.6];
        assert_slice_close(&linear_shape(&positioned).unwrap(), &[0.08, 0.12, 0.8]);
        assert_slice_close(
            &convex_shape(&positioned).unwrap(),
            &[
                0.009_347_373_623_712_36,
                0.020_175_225_787_320_738,
                0.690_983_005_625_052_5,
            ],
        );
        assert_close(mixed_shape(&positioned, 5, 1.0), 0.8);
        assert_close(disc_shape(&positioned, 5, 1.0, 1.0), 0.8);
    }

    #[test]
    fn pinned_reference_probe_matches_wfg1_full_pipeline() {
        // Frozen from an independent direct f64 port of the corrected
        // normalized equations at the pinned jMetal revision.  M=4, k=6,
        // and l=4 expose global weights and every b_flat region.
        let problem = Wfg1::new(4, 6, 4).unwrap();
        let evaluation = problem
            .evaluate_normalized(&[0.11, 0.23, 0.37, 0.49, 0.58, 0.72, 0.61, 0.85, 0.89, 0.97])
            .unwrap();

        assert_slice_close(
            evaluation.transformed(),
            &[
                0.956_814_732_462_447_4,
                0.971_034_268_446_919_2,
                0.980_311_358_064_300_5,
                0.985_834_293_575_054_4,
                0.989_164_587_101_819,
                0.993_451_454_455_177_4,
                0.983_109_231_740_204,
                0.995_547_072_784_466_3,
                0.995_547_072_784_466_3,
                0.998_730_538_334_589_8,
            ],
        );
        assert_slice_close(
            evaluation.reduced(),
            &[
                0.966_294_423_118_761_9,
                0.983_467_321_213_302_8,
                0.991_502_878_385_469,
                0.993_922_654_201_860_4,
            ],
        );
        assert_slice_close(evaluation.positioned(), evaluation.reduced());
        assert_slice_close(
            evaluation.shape(),
            &[
                0.910_175_423_133_181_4,
                0.000_082_168_919_713_291_96,
                0.000_319_343_833_051_351_03,
                0.005_954_899_528_811_435,
            ],
        );
        assert_slice_close(
            evaluation.objectives(),
            &[
                2.814_273_500_468_223_3,
                0.994_251_329_880_713_6,
                0.995_838_717_200_168_6,
                1.041_561_850_432_351_8,
            ],
        );
    }

    #[test]
    fn pinned_reference_probe_matches_wfg2_full_pipeline() {
        // The two unequal distance pairs make adjacent-pair compression and
        // the post-compression distance slice independently observable.
        let problem = Wfg2::new(4, 6, 4).unwrap();
        let evaluation = problem
            .evaluate_normalized(&[0.11, 0.23, 0.37, 0.49, 0.58, 0.72, 0.61, 0.85, 0.89, 0.97])
            .unwrap();

        assert_slice_close(
            evaluation.transformed(),
            &[
                0.11,
                0.23,
                0.37,
                0.49,
                0.58,
                0.72,
                0.635_897_435_897_435_9,
                0.676_923_076_923_076_8,
            ],
        );
        assert_slice_close(
            evaluation.reduced(),
            &[0.17, 0.43, 0.65, 0.656_410_256_410_256_3],
        );
        assert_slice_close(evaluation.positioned(), evaluation.reduced());
        assert_slice_close(
            evaluation.shape(),
            &[
                0.003_715_970_218_770_366,
                0.001_146_770_920_965_606,
                0.013_282_367_711_360_762,
                0.865_038_253_555_139_7,
            ],
        );
        assert_slice_close(
            evaluation.objectives(),
            &[
                0.663_842_196_847_797,
                0.660_997_340_094_118_7,
                0.736_104_462_678_420_9,
                7.576_716_284_851_374,
            ],
        );
    }

    #[test]
    fn pinned_reference_probe_matches_wfg3_full_pipeline() {
        // The interior nonzero distance exercises both degenerate position
        // coordinates from A=[1, 0, 0] before the linear shape is applied.
        let problem = Wfg3::new(4, 6, 4).unwrap();
        let evaluation = problem
            .evaluate_normalized(&[0.11, 0.23, 0.37, 0.49, 0.58, 0.72, 0.61, 0.85, 0.89, 0.97])
            .unwrap();

        assert_slice_close(
            evaluation.transformed(),
            &[
                0.11,
                0.23,
                0.37,
                0.49,
                0.58,
                0.72,
                0.635_897_435_897_435_9,
                0.676_923_076_923_076_8,
            ],
        );
        assert_slice_close(
            evaluation.reduced(),
            &[0.17, 0.43, 0.65, 0.656_410_256_410_256_3],
        );
        assert_slice_close(
            evaluation.positioned(),
            &[
                0.17,
                0.454_051_282_051_282_03,
                0.598_461_538_461_538_4,
                0.656_410_256_410_256_3,
            ],
        );
        assert_slice_close(
            evaluation.shape(),
            &[
                0.046_194_478_895_463_506,
                0.030_994_239_053_254_446,
                0.092_811_282_051_282_06,
                0.83,
            ],
        );
        assert_slice_close(
            evaluation.objectives(),
            &[
                0.748_799_214_201_183_3,
                0.780_387_212_623_274,
                1.213_277_948_717_948_6,
                7.296_410_256_410_256,
            ],
        );
    }

    #[test]
    fn wfg2_pair_reduction_is_swap_invariant_but_not_regrouping_invariant() {
        let problem = Wfg2::new(4, 6, 4).unwrap();
        let original = problem
            .evaluate_normalized(&[0.11, 0.23, 0.37, 0.49, 0.58, 0.72, 0.61, 0.85, 0.89, 0.97])
            .unwrap();
        let pair_swapped = problem
            .evaluate_normalized(&[0.11, 0.23, 0.37, 0.49, 0.58, 0.72, 0.85, 0.61, 0.97, 0.89])
            .unwrap();
        let cross_regrouped = problem
            .evaluate_normalized(&[0.11, 0.23, 0.37, 0.49, 0.58, 0.72, 0.61, 0.89, 0.85, 0.97])
            .unwrap();

        assert_slice_close(original.transformed(), pair_swapped.transformed());
        assert_slice_close(original.reduced(), pair_swapped.reduced());
        assert_slice_close(original.shape(), pair_swapped.shape());
        assert_slice_close(original.objectives(), pair_swapped.objectives());
        assert!((original.transformed()[6] - cross_regrouped.transformed()[6]).abs() > 0.05);
        assert!((original.transformed()[7] - cross_regrouped.transformed()[7]).abs() > 0.02);
    }

    #[test]
    fn wfg3_zero_distance_collapses_only_degenerate_position_coordinates() {
        let evaluation = Wfg3::new(4, 6, 4)
            .unwrap()
            .evaluate_normalized(&[0.1, 0.3, 0.2, 0.6, 0.4, 0.8, 0.35, 0.35, 0.35, 0.35])
            .unwrap();

        assert_slice_close(evaluation.reduced(), &[0.2, 0.4, 0.6, 0.0]);
        assert_slice_close(evaluation.positioned(), &[0.2, 0.5, 0.5, 0.0]);
        assert_slice_close(evaluation.shape(), &[0.05, 0.05, 0.1, 0.8]);
        assert_slice_close(evaluation.objectives(), &[0.1, 0.2, 0.6, 6.4]);
    }

    #[test]
    fn wfg1_through_wfg3_repeated_evaluations_are_bitwise_identical() {
        let input = [0.11, 0.23, 0.37, 0.49, 0.58, 0.72, 0.61, 0.85, 0.89, 0.97];

        let problem = Wfg1::new(4, 6, 4).unwrap();
        let first = problem.evaluate_normalized(&input).unwrap();
        let second = problem.evaluate_normalized(&input).unwrap();
        assert_evaluation_bits_eq(&first, &second);
        assert_slice_bits_eq(&first.clone().into_objectives(), second.objectives());

        let problem = Wfg2::new(4, 6, 4).unwrap();
        let first = problem.evaluate_normalized(&input).unwrap();
        let second = problem.evaluate_normalized(&input).unwrap();
        assert_evaluation_bits_eq(&first, &second);
        assert_slice_bits_eq(&first.clone().into_objectives(), second.objectives());

        let problem = Wfg3::new(4, 6, 4).unwrap();
        let first = problem.evaluate_normalized(&input).unwrap();
        let second = problem.evaluate_normalized(&input).unwrap();
        assert_evaluation_bits_eq(&first, &second);
        assert_slice_bits_eq(&first.clone().into_objectives(), second.objectives());
    }

    #[test]
    fn canonical_m3_k4_l20_extreme_is_recovered() {
        let problem = Wfg4::new(3, 4, 20).unwrap();
        let evaluation = problem
            .evaluate_normalized(&vec![S_MULTI_CENTER; problem.dimension()])
            .unwrap();

        assert_slice_close(evaluation.transformed(), &[0.0; 24]);
        assert_slice_close(evaluation.reduced(), &[0.0, 0.0, 0.0]);
        assert_slice_close(evaluation.positioned(), evaluation.reduced());
        assert_slice_close(evaluation.shape(), &[0.0, 0.0, 1.0]);
        assert_slice_close(evaluation.objectives(), &[0.0, 0.0, 6.0]);
    }

    #[test]
    fn interior_half_reduction_catches_shape_index_and_scale_mutants() {
        let problem = Wfg4::new(3, 4, 2).unwrap();
        let evaluation = problem
            .evaluate_normalized(&[0.35, 0.0, 0.35, 0.0, 0.35, 0.0])
            .unwrap();

        assert_slice_close(evaluation.reduced(), &[0.5, 0.5, 0.5]);
        assert_slice_close(evaluation.positioned(), evaluation.reduced());
        assert_slice_close(evaluation.shape(), &[0.5, 0.5, 0.5_f64.sqrt()]);
        assert_slice_close(
            evaluation.objectives(),
            &[1.5, 2.5, 0.5 + 3.0 * 2.0_f64.sqrt()],
        );
    }

    #[test]
    fn pinned_reference_interior_probe_matches_full_pipeline() {
        // Frozen from an independent f64 port of the corrected equations at
        // jMetal revision ea7e882f6b8f94b99535921674e62cda7986f20e.
        // Non-anchor coordinates make A/B/phase transformation mutants
        // observable; asymmetric position and distance blocks make slicing,
        // shape, scale, and distance mutants observable in the same KAT.
        let problem = Wfg4::new(3, 4, 20).unwrap();
        let input = [
            0.00, 0.37, 0.74, 0.10, 0.47, 0.84, 0.20, 0.57, 0.94, 0.30, 0.67, 0.03, 0.40, 0.77,
            0.13, 0.50, 0.87, 0.23, 0.60, 0.97, 0.33, 0.70, 0.06, 0.43,
        ];
        let evaluation = problem.evaluate_normalized(&input).unwrap();

        for (index, expected) in [
            (0, 1.0),
            (1, 0.006_941_067_221_959_935),
            (2, 0.409_084_749_531_246_85),
            (3, 0.489_959_990_197_517_96),
            (4, 0.168_487_021_762_804_94),
            (23, 0.093_942_962_058_725_45),
        ] {
            assert_close(evaluation.transformed()[index], expected);
        }
        assert_slice_close(
            evaluation.reduced(),
            &[
                0.503_470_533_610_979_9,
                0.449_522_369_864_382_43,
                0.354_308_186_433_816_06,
            ],
        );
        assert_slice_close(evaluation.positioned(), evaluation.reduced());
        assert_slice_close(
            evaluation.shape(),
            &[
                0.461_320_042_089_157_9,
                0.540_957_680_606_666,
                0.703_241_499_457_699_8,
            ],
        );
        assert_slice_close(
            evaluation.objectives(),
            &[
                1.276_948_270_612_131_8,
                2.518_138_908_860_48,
                4.573_757_183_180_015,
            ],
        );
    }

    #[test]
    fn concave_anchors_generalize_across_objective_counts() {
        for objectives in [2, 3, 5] {
            let position_parameters = 2 * (objectives - 1);
            let problem = Wfg4::new(objectives, position_parameters, 2).unwrap();

            let center = problem
                .evaluate_normalized(&vec![S_MULTI_CENTER; problem.dimension()])
                .unwrap();
            let mut expected_center = vec![0.0; objectives];
            expected_center[objectives - 1] = 2.0 * objectives as f64;
            assert_slice_close(center.objectives(), &expected_center);

            let boundary = problem
                .evaluate_normalized(&vec![0.0; problem.dimension()])
                .unwrap();
            let mut expected_boundary = vec![1.0; objectives];
            expected_boundary[0] = 3.0;
            assert_slice_close(boundary.transformed(), &vec![1.0; problem.dimension()]);
            assert_slice_close(boundary.objectives(), &expected_boundary);
        }
    }

    #[test]
    fn asymmetric_groups_use_the_exact_declared_slices() {
        let problem = Wfg4::new(4, 6, 3).unwrap();
        let input = [0.05, 0.15, 0.25, 0.45, 0.65, 0.75, 0.85, 0.95, 0.35];
        let evaluation = problem.evaluate_normalized(&input).unwrap();
        let transformed: Vec<f64> = input.into_iter().map(s_multi).collect();
        let expected = [
            (transformed[0] + transformed[1]) / 2.0,
            (transformed[2] + transformed[3]) / 2.0,
            (transformed[4] + transformed[5]) / 2.0,
            (transformed[6] + transformed[7] + transformed[8]) / 3.0,
        ];
        assert_slice_close(evaluation.reduced(), &expected);
    }

    #[test]
    fn two_element_reduction_groups_are_bitwise_pair_swap_invariant() {
        let problem = Wfg4::new(3, 4, 2).unwrap();
        let input = [0.07, 0.21, 0.46, 0.72, 0.18, 0.91];
        let permuted = [0.21, 0.07, 0.72, 0.46, 0.91, 0.18];
        let first = problem.evaluate_normalized(&input).unwrap();
        let second = problem.evaluate_normalized(&permuted).unwrap();

        assert_eq!(first.reduced(), second.reduced());
        assert_eq!(first.positioned(), second.positioned());
        assert_eq!(first.shape(), second.shape());
        assert_eq!(first.objectives(), second.objectives());
    }

    #[test]
    fn wrong_center_and_same_objective_shape_mutants_are_observable() {
        let problem = Wfg4::new(3, 4, 2).unwrap();
        let center_input = [S_MULTI_CENTER; 6];
        let canonical_center = problem.evaluate_normalized(&center_input).unwrap();

        let wrong_center_distance = center_input[4..]
            .iter()
            .map(|&value| {
                let denominator = if value <= 0.5 { 1.0 } else { -1.0 };
                let ratio = (value - 0.5).abs() / denominator;
                let phase = (4.0 * S_MULTI_A + 2.0) * core::f64::consts::PI * (0.5 - ratio);
                let quadratic = 4.0 * S_MULTI_B * ratio * ratio;
                correct_to_01((1.0 + fs_math::det::cos(phase) + quadratic) / (S_MULTI_B + 2.0))
            })
            .sum::<f64>()
            / 2.0;
        assert!(wrong_center_distance > 0.05);

        let interior = problem
            .evaluate_normalized(&[0.07, 0.21, 0.46, 0.72, 0.18, 0.91])
            .unwrap();
        let convex_mutant_first = (1.0
            - fs_math::det::cos(interior.reduced()[0] * core::f64::consts::FRAC_PI_2))
            * (1.0 - fs_math::det::cos(interior.reduced()[1] * core::f64::consts::FRAC_PI_2));
        assert_ne!(convex_mutant_first.to_bits(), interior.shape()[0].to_bits());
        assert_close(canonical_center.reduced()[2], 0.0);
    }

    #[test]
    fn repeated_evaluation_is_bitwise_identical() {
        let problem = Wfg4::new(5, 8, 7).unwrap();
        let input: Vec<f64> = (0..problem.dimension())
            .map(|index| ((index * 37 + 19) % 101) as f64 / 100.0)
            .collect();
        let first = problem.evaluate_normalized(&input).unwrap();
        let second = problem.evaluate_normalized(&input).unwrap();

        assert_evaluation_bits_eq(&first, &second);
        assert_slice_bits_eq(&first.clone().into_objectives(), second.objectives());
    }
}
