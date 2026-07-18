//! Typed kernels for the Walking Fish Group benchmark family.
//!
//! The production kernel implements the complete normalized WFG1 through WFG9
//! compositions for arbitrary valid objective, position-parameter, and
//! distance-parameter counts.  Its transformations and shapes follow the
//! corrected WFG toolkit as represented by jMetal revision
//! `ea7e882f6b8f94b99535921674e62cda7986f20e`. Callers may evaluate either
//! normalized coordinates in `[0, 1]` or canonical heterogeneous coordinates
//! `z_i in [0, 2(i + 1)]`; both paths enter the same normalized kernel after
//! structured admission.
//!
//! Determinism is structural within the evaluator: reductions have a fixed
//! left-to-right order and transcendental calls use [`fs_math::det`].  This
//! module does not yet claim executable parity with an external oracle,
//! optimizer convergence, cancellation coverage, cross-ISA bit stability, or
//! performance evidence. Direct WFG6/WFG9 nonseparable reductions and
//! WFG7/WFG8 conditioning retain canonical fixed-order arithmetic; no
//! subquadratic complexity claim is made.

#![deny(unsafe_code)]

const CORRECTION_EPSILON: f64 = 1.0e-10;
const S_MULTI_A: f64 = 30.0;
const S_MULTI_B: f64 = 10.0;
const S_MULTI_WFG9_B: f64 = 95.0;
const S_MULTI_CENTER: f64 = 0.35;
const S_DECEPT_CENTER: f64 = 0.35;
const S_DECEPT_RADIUS: f64 = 0.001;
const S_DECEPT_FLOOR: f64 = 0.05;
const B_PARAM_A: f64 = 0.98 / 49.98;
const B_PARAM_B: f64 = 0.02;
const B_PARAM_C: f64 = 50.0;

/// Structured refusal from a typed WFG evaluator.
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
    /// The decision vector has the wrong dimension.
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
    /// A canonical-domain coordinate is NaN or infinite.
    NonFiniteCanonicalInput {
        /// Zero-based coordinate index.
        component: usize,
        /// Exact IEEE-754 payload.
        bits: u64,
    },
    /// A finite canonical coordinate lies outside `[0, 2(i + 1)]`.
    CanonicalInputOutOfRange {
        /// Zero-based coordinate index.
        component: usize,
        /// Exact IEEE-754 payload of the rejected coordinate.
        bits: u64,
        /// Exact IEEE-754 payload of the coordinate's inclusive upper bound.
        upper_bound_bits: u64,
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
                "WFG expected {expected} decision coordinates, received {actual}"
            ),
            Self::NonFiniteInput { component, bits } => write!(
                formatter,
                "WFG normalized coordinate {component} is non-finite (bits 0x{bits:016x})"
            ),
            Self::InputOutOfRange { component, bits } => write!(
                formatter,
                "WFG normalized coordinate {component} lies outside [0, 1] (bits 0x{bits:016x})"
            ),
            Self::NonFiniteCanonicalInput { component, bits } => write!(
                formatter,
                "WFG canonical coordinate {component} is non-finite (bits 0x{bits:016x})"
            ),
            Self::CanonicalInputOutOfRange {
                component,
                bits,
                upper_bound_bits,
            } => write!(
                formatter,
                "WFG canonical coordinate {component} lies outside its inclusive bound (value bits 0x{bits:016x}, upper-bound bits 0x{upper_bound_bits:016x})"
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

    fn validate_length(self, input: &[f64]) -> Result<(), WfgError> {
        if input.len() != self.dimension {
            return Err(WfgError::WrongInputLength {
                expected: self.dimension,
                actual: input.len(),
            });
        }
        Ok(())
    }

    fn validate_input(self, input: &[f64]) -> Result<(), WfgError> {
        self.validate_length(input)?;
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

    fn validate_canonical_input(self, input: &[f64]) -> Result<(), WfgError> {
        self.validate_length(input)?;
        for (component, &value) in input.iter().enumerate() {
            if !value.is_finite() {
                return Err(WfgError::NonFiniteCanonicalInput {
                    component,
                    bits: value.to_bits(),
                });
            }
            let upper_bound = canonical_upper_bound(component);
            if !(0.0..=upper_bound).contains(&value) {
                return Err(WfgError::CanonicalInputOutOfRange {
                    component,
                    bits: value.to_bits(),
                    upper_bound_bits: upper_bound.to_bits(),
                });
            }
        }
        Ok(())
    }

    fn normalize_canonical(self, input: &[f64]) -> Result<Vec<f64>, WfgError> {
        self.normalize_canonical_with(input, reserved_vec)
    }

    fn normalize_canonical_with(
        self,
        input: &[f64],
        reserve: impl FnOnce(&'static str, usize) -> Result<Vec<f64>, WfgError>,
    ) -> Result<Vec<f64>, WfgError> {
        self.validate_canonical_input(input)?;
        let mut normalized = reserve("canonical normalized coordinates", self.dimension)?;
        for (component, &value) in input.iter().enumerate() {
            normalized.push(value / canonical_upper_bound(component));
        }
        Ok(normalized)
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

        /// Total decision dimension, `k + l`.
        #[must_use]
        pub const fn dimension(self) -> usize {
            self.spec.dimension
        }

        /// Evaluate one decision vector in the canonical heterogeneous WFG
        /// domain, where coordinate `i` is bounded by `[0, 2(i + 1)]`.
        ///
        /// # Errors
        ///
        /// Returns a structured refusal for a wrong dimension, a non-finite
        /// coordinate, a coordinate outside its indexed canonical bound, or
        /// failed normalization/intermediate storage admission. No WFG
        /// transformation runs before the complete canonical vector is
        /// admitted and normalized.
        pub fn evaluate_canonical(&self, input: &[f64]) -> Result<WfgEvaluation, WfgError> {
            let normalized = self.spec.normalize_canonical(input)?;
            self.evaluate_normalized(&normalized)
        }
    };
}

/// A validated WFG1 problem definition.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Wfg1 {
    spec: WfgSpec,
}

impl Wfg1 {
    /// Validate a WFG1 problem definition.
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

/// A validated WFG2 problem definition.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Wfg2 {
    spec: WfgSpec,
}

impl Wfg2 {
    /// Validate a WFG2 problem definition.
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

/// A validated WFG3 problem definition.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Wfg3 {
    spec: WfgSpec,
}

impl Wfg3 {
    /// Validate a WFG3 problem definition.
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

/// A validated WFG4 problem definition.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Wfg4 {
    spec: WfgSpec,
}

impl Wfg4 {
    /// Validate a WFG4 problem definition.
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

/// A validated WFG5 problem definition.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Wfg5 {
    spec: WfgSpec,
}

impl Wfg5 {
    /// Validate a WFG5 problem definition.
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
    /// storage admission. No transformation runs before input validation.
    pub fn evaluate_normalized(&self, input: &[f64]) -> Result<WfgEvaluation, WfgError> {
        self.spec.validate_input(input)?;
        let transformed = wfg5_transform(input)?;
        let reduced = equal_group_reduce(&transformed, self.spec)?;
        let positioned = identity_positioned(&reduced)?;
        let shape = concave_shape(&positioned)?;
        finish_evaluation(transformed, reduced, positioned, shape)
    }
}

/// A validated WFG6 problem definition.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Wfg6 {
    spec: WfgSpec,
}

impl Wfg6 {
    /// Validate a WFG6 problem definition.
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
    /// storage admission. No transformation runs before input validation.
    pub fn evaluate_normalized(&self, input: &[f64]) -> Result<WfgEvaluation, WfgError> {
        self.spec.validate_input(input)?;
        let transformed = wfg6_transform(input, self.spec.position_parameters)?;
        let reduced =
            nonseparable_group_reduce(&transformed, self.spec, "WFG6 reduced coordinates")?;
        let positioned = identity_positioned(&reduced)?;
        let shape = concave_shape(&positioned)?;
        finish_evaluation(transformed, reduced, positioned, shape)
    }
}

/// A validated WFG7 problem definition.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Wfg7 {
    spec: WfgSpec,
}

impl Wfg7 {
    /// Validate a WFG7 problem definition.
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
    /// storage admission. No transformation runs before input validation.
    pub fn evaluate_normalized(&self, input: &[f64]) -> Result<WfgEvaluation, WfgError> {
        self.spec.validate_input(input)?;
        let transformed = wfg7_transform(input, self.spec.position_parameters)?;
        let reduced = equal_group_reduce(&transformed, self.spec)?;
        let positioned = identity_positioned(&reduced)?;
        let shape = concave_shape(&positioned)?;
        finish_evaluation(transformed, reduced, positioned, shape)
    }
}

/// A validated WFG8 problem definition.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Wfg8 {
    spec: WfgSpec,
}

impl Wfg8 {
    /// Validate a WFG8 problem definition.
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
    /// storage admission. No transformation runs before input validation.
    pub fn evaluate_normalized(&self, input: &[f64]) -> Result<WfgEvaluation, WfgError> {
        self.spec.validate_input(input)?;
        let transformed = wfg8_transform(input, self.spec.position_parameters)?;
        let reduced = equal_group_reduce(&transformed, self.spec)?;
        let positioned = identity_positioned(&reduced)?;
        let shape = concave_shape(&positioned)?;
        finish_evaluation(transformed, reduced, positioned, shape)
    }
}

/// A validated WFG9 problem definition.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Wfg9 {
    spec: WfgSpec,
}

impl Wfg9 {
    /// Validate a WFG9 problem definition.
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
    /// storage admission. No transformation runs before input validation.
    pub fn evaluate_normalized(&self, input: &[f64]) -> Result<WfgEvaluation, WfgError> {
        self.spec.validate_input(input)?;
        let transformed = wfg9_transform(input, self.spec.position_parameters)?;
        let reduced =
            nonseparable_group_reduce(&transformed, self.spec, "WFG9 reduced coordinates")?;
        let positioned = identity_positioned(&reduced)?;
        let shape = concave_shape(&positioned)?;
        finish_evaluation(transformed, reduced, positioned, shape)
    }
}

/// One WFG evaluation with replay-relevant normalized intermediates.
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

fn canonical_upper_bound(component: usize) -> f64 {
    2.0 * (component as f64 + 1.0)
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

fn wfg5_transform(input: &[f64]) -> Result<Vec<f64>, WfgError> {
    let mut transformed = reserved_vec("WFG5 transformed coordinates", input.len())?;
    for &value in input {
        transformed.push(s_decept(
            value,
            S_DECEPT_CENTER,
            S_DECEPT_RADIUS,
            S_DECEPT_FLOOR,
        ));
    }
    Ok(transformed)
}

fn wfg6_transform(input: &[f64], position_parameters: usize) -> Result<Vec<f64>, WfgError> {
    let mut transformed = reserved_vec("WFG6 transformed coordinates", input.len())?;
    transformed.extend_from_slice(&input[..position_parameters]);
    for &value in &input[position_parameters..] {
        transformed.push(s_linear(value, 0.35));
    }
    Ok(transformed)
}

fn wfg7_transform(input: &[f64], position_parameters: usize) -> Result<Vec<f64>, WfgError> {
    let mut transformed = reserved_vec("WFG7 transformed coordinates", input.len())?;
    for (index, &value) in input[..position_parameters].iter().enumerate() {
        let conditioning = equal_weight_reduction(&input[index + 1..]);
        transformed.push(b_param(
            value,
            conditioning,
            B_PARAM_A,
            B_PARAM_B,
            B_PARAM_C,
        ));
    }
    for &value in &input[position_parameters..] {
        transformed.push(s_linear(value, 0.35));
    }
    Ok(transformed)
}

fn wfg8_transform(input: &[f64], position_parameters: usize) -> Result<Vec<f64>, WfgError> {
    let mut transformed = reserved_vec("WFG8 transformed coordinates", input.len())?;
    transformed.extend_from_slice(&input[..position_parameters]);
    for (offset, &value) in input[position_parameters..].iter().enumerate() {
        let index = position_parameters + offset;
        let conditioning = equal_weight_reduction(&input[..index]);
        let biased = b_param(value, conditioning, B_PARAM_A, B_PARAM_B, B_PARAM_C);
        transformed.push(s_linear(biased, 0.35));
    }
    Ok(transformed)
}

fn wfg9_transform(input: &[f64], position_parameters: usize) -> Result<Vec<f64>, WfgError> {
    let mut transformed = reserved_vec("WFG9 transformed coordinates", input.len())?;
    for (index, &value) in input.iter().enumerate() {
        let biased = if index + 1 < input.len() {
            let conditioning = equal_weight_reduction(&input[index + 1..]);
            b_param(value, conditioning, B_PARAM_A, B_PARAM_B, B_PARAM_C)
        } else {
            value
        };
        let final_value = if index < position_parameters {
            s_decept(biased, S_DECEPT_CENTER, S_DECEPT_RADIUS, S_DECEPT_FLOOR)
        } else {
            s_multi_with_parameters(biased, S_MULTI_A, S_MULTI_WFG9_B, S_MULTI_CENTER)
        };
        transformed.push(final_value);
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

fn nonseparable_group_reduce(
    transformed: &[f64],
    spec: WfgSpec,
    allocation_name: &'static str,
) -> Result<Vec<f64>, WfgError> {
    let mut reduced = reserved_vec(allocation_name, spec.objectives)?;
    for group in 0..spec.objectives - 1 {
        let start = group * spec.position_group_size;
        let end = start + spec.position_group_size;
        reduced.push(r_nonsep(&transformed[start..end], spec.position_group_size));
    }
    reduced.push(r_nonsep(
        &transformed[spec.position_parameters..],
        spec.distance_parameters,
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

fn b_param(value: f64, conditioning: f64, a: f64, b: f64, c: f64) -> f64 {
    debug_assert!((0.0..=1.0).contains(&value));
    debug_assert!((0.0..=1.0).contains(&conditioning));
    debug_assert!(0.0 < a && a < 1.0);
    debug_assert!(0.0 < b && b < c);
    let exponent_coordinate =
        a - (1.0 - 2.0 * conditioning) * ((0.5 - conditioning).floor() + a).abs();
    let exponent = (c - b).mul_add(exponent_coordinate, b);
    correct_to_01(fs_math::det::pow(value, exponent))
}

fn s_decept(value: f64, center: f64, radius: f64, floor_value: f64) -> f64 {
    debug_assert!((0.0..=1.0).contains(&value));
    debug_assert!(0.0 < center - radius && center + radius < 1.0);
    debug_assert!((0.0..=1.0).contains(&floor_value));
    let below = (value - center + radius).floor()
        * (1.0 - floor_value + (center - radius) / radius)
        / (center - radius);
    let above = (center + radius - value).floor()
        * (1.0 - floor_value + (1.0 - center - radius) / radius)
        / (1.0 - center - radius);
    correct_to_01(1.0 + ((value - center).abs() - radius) * (below + above + 1.0 / radius))
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
    s_multi_with_parameters(value, S_MULTI_A, S_MULTI_B, S_MULTI_CENTER)
}

fn s_multi_with_parameters(value: f64, modes: f64, bias: f64, center: f64) -> f64 {
    debug_assert!(modes > 0.0 && bias > 0.0);
    debug_assert!(0.0 < center && center < 1.0);
    let denominator = if value <= center {
        2.0 * center
    } else {
        2.0 * (center - 1.0)
    };
    let ratio = (value - center).abs() / denominator;
    let phase = (4.0 * modes + 2.0) * core::f64::consts::PI * (0.5 - ratio);
    let quadratic = 4.0 * bias * ratio * ratio;
    correct_to_01((1.0 + fs_math::det::cos(phase) + quadratic) / (bias + 2.0))
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
    use core::cell::Cell;

    const TOLERANCE: f64 = 1.0e-11;

    fn assert_close(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() <= TOLERANCE,
            "actual={actual:.17e}, expected={expected:.17e}"
        );
    }

    fn assert_tiny_relative_close(actual: f64, expected: f64) {
        assert!(
            actual.is_finite() && expected.is_finite() && expected != 0.0,
            "tiny relative oracle requires finite values and a nonzero reference: actual={actual:.17e}, expected={expected:.17e}"
        );
        let relative_error = ((actual - expected) / expected).abs();
        assert!(
            relative_error <= 1.0e-12,
            "actual={actual:.17e}, expected={expected:.17e}, relative_error={relative_error:.17e}"
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

    fn assert_canonical_refusal_precedes_reservation(
        spec: WfgSpec,
        input: &[f64],
        expected: WfgError,
    ) {
        let reservation_calls = Cell::new(0);
        let actual = spec
            .normalize_canonical_with(input, |_, _| {
                reservation_calls.set(reservation_calls.get() + 1);
                Ok(Vec::new())
            })
            .unwrap_err();
        assert_eq!(actual, expected);
        assert_eq!(reservation_calls.get(), 0);
    }

    #[test]
    fn canonical_normalization_preserves_allocation_failure_and_admission_order() {
        let spec = WfgSpec::new(4, 6, 6).unwrap();
        let reservation_calls = Cell::new(0);
        let actual = spec
            .normalize_canonical_with(&[0.0; 12], |what, elements| {
                reservation_calls.set(reservation_calls.get() + 1);
                Err(WfgError::AllocationFailed { what, elements })
            })
            .unwrap_err();
        assert_eq!(
            actual,
            WfgError::AllocationFailed {
                what: "canonical normalized coordinates",
                elements: 12,
            }
        );
        assert_eq!(reservation_calls.get(), 1);

        assert_canonical_refusal_precedes_reservation(
            spec,
            &[f64::NAN; 11],
            WfgError::WrongInputLength {
                expected: 12,
                actual: 11,
            },
        );

        let mut nonfinite = [0.0; 12];
        nonfinite[3] = f64::from_bits(0x7ff8_1234_5678_9abc);
        assert_canonical_refusal_precedes_reservation(
            spec,
            &nonfinite,
            WfgError::NonFiniteCanonicalInput {
                component: 3,
                bits: nonfinite[3].to_bits(),
            },
        );

        let mut out_of_range = [0.0; 12];
        let upper_bound = canonical_upper_bound(7);
        out_of_range[7] = f64::from_bits(upper_bound.to_bits() + 1);
        assert_canonical_refusal_precedes_reservation(
            spec,
            &out_of_range,
            WfgError::CanonicalInputOutOfRange {
                component: 7,
                bits: out_of_range[7].to_bits(),
                upper_bound_bits: upper_bound.to_bits(),
            },
        );
    }

    fn assert_common_input_refusals(evaluate: impl Fn(&[f64]) -> Result<WfgEvaluation, WfgError>) {
        assert_eq!(
            evaluate(&[0.0; 5]).unwrap_err(),
            WfgError::WrongInputLength {
                expected: 6,
                actual: 5,
            }
        );
        assert_eq!(
            evaluate(&[0.0; 7]).unwrap_err(),
            WfgError::WrongInputLength {
                expected: 6,
                actual: 7,
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
                f64::NEG_INFINITY,
                WfgError::NonFiniteInput {
                    component: 2,
                    bits: f64::NEG_INFINITY.to_bits(),
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
            Wfg5::new(1, 1, 2).unwrap_err(),
            Wfg6::new(1, 1, 2).unwrap_err(),
            Wfg7::new(1, 1, 2).unwrap_err(),
            Wfg8::new(1, 1, 2).unwrap_err(),
            Wfg9::new(1, 1, 2).unwrap_err(),
        ] {
            assert_eq!(error, WfgError::TooFewObjectives { objectives: 1 });
        }
        for error in [
            Wfg1::new(4, 4, 2).unwrap_err(),
            Wfg2::new(4, 4, 2).unwrap_err(),
            Wfg3::new(4, 4, 2).unwrap_err(),
            Wfg5::new(4, 4, 2).unwrap_err(),
            Wfg6::new(4, 4, 2).unwrap_err(),
            Wfg7::new(4, 4, 2).unwrap_err(),
            Wfg8::new(4, 4, 2).unwrap_err(),
            Wfg9::new(4, 4, 2).unwrap_err(),
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
            Wfg5::new(3, 4, 0).unwrap_err(),
            Wfg6::new(3, 4, 0).unwrap_err(),
            Wfg7::new(3, 4, 0).unwrap_err(),
            Wfg8::new(3, 4, 0).unwrap_err(),
            Wfg9::new(3, 4, 0).unwrap_err(),
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
        assert!(Wfg5::new(3, 4, 3).is_ok());
        assert!(Wfg6::new(3, 4, 3).is_ok());
        assert!(Wfg7::new(3, 4, 3).is_ok());
        assert!(Wfg8::new(3, 4, 3).is_ok());
        assert!(Wfg9::new(3, 4, 3).is_ok());
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

        let wfg5 = Wfg5::new(3, 4, 3).unwrap();
        assert_eq!(wfg5.objectives(), 3);
        assert_eq!(wfg5.position_parameters(), 4);
        assert_eq!(wfg5.distance_parameters(), 3);
        assert_eq!(wfg5.dimension(), 7);

        let wfg6 = Wfg6::new(4, 6, 2).unwrap();
        assert_eq!(wfg6.objectives(), 4);
        assert_eq!(wfg6.position_parameters(), 6);
        assert_eq!(wfg6.distance_parameters(), 2);
        assert_eq!(wfg6.dimension(), 8);

        let wfg7 = Wfg7::new(5, 8, 4).unwrap();
        assert_eq!(wfg7.objectives(), 5);
        assert_eq!(wfg7.position_parameters(), 8);
        assert_eq!(wfg7.distance_parameters(), 4);
        assert_eq!(wfg7.dimension(), 12);

        let wfg8 = Wfg8::new(3, 4, 3).unwrap();
        assert_eq!(wfg8.objectives(), 3);
        assert_eq!(wfg8.position_parameters(), 4);
        assert_eq!(wfg8.distance_parameters(), 3);
        assert_eq!(wfg8.dimension(), 7);

        let wfg9 = Wfg9::new(4, 6, 5).unwrap();
        assert_eq!(wfg9.objectives(), 4);
        assert_eq!(wfg9.position_parameters(), 6);
        assert_eq!(wfg9.distance_parameters(), 5);
        assert_eq!(wfg9.dimension(), 11);
    }

    #[test]
    fn input_admission_is_exact_and_structured() {
        let wfg1 = Wfg1::new(3, 4, 2).unwrap();
        let wfg2 = Wfg2::new(3, 4, 2).unwrap();
        let wfg3 = Wfg3::new(3, 4, 2).unwrap();
        let wfg4 = Wfg4::new(3, 4, 2).unwrap();
        let wfg5 = Wfg5::new(3, 4, 2).unwrap();
        let wfg6 = Wfg6::new(3, 4, 2).unwrap();
        let wfg7 = Wfg7::new(3, 4, 2).unwrap();
        let wfg8 = Wfg8::new(3, 4, 2).unwrap();
        let wfg9 = Wfg9::new(3, 4, 2).unwrap();

        assert_common_input_refusals(|input| wfg1.evaluate_normalized(input));
        assert_common_input_refusals(|input| wfg2.evaluate_normalized(input));
        assert_common_input_refusals(|input| wfg3.evaluate_normalized(input));
        assert_common_input_refusals(|input| wfg4.evaluate_normalized(input));
        assert_common_input_refusals(|input| wfg5.evaluate_normalized(input));
        assert_common_input_refusals(|input| wfg6.evaluate_normalized(input));
        assert_common_input_refusals(|input| wfg7.evaluate_normalized(input));
        assert_common_input_refusals(|input| wfg8.evaluate_normalized(input));
        assert_common_input_refusals(|input| wfg9.evaluate_normalized(input));
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
    fn canonical_wfg5_through_wfg7_bias_anchors_match_the_toolkit() {
        assert_close(
            s_decept(
                S_DECEPT_CENTER,
                S_DECEPT_CENTER,
                S_DECEPT_RADIUS,
                S_DECEPT_FLOOR,
            ),
            0.0,
        );
        assert_close(
            s_decept(
                S_DECEPT_CENTER - S_DECEPT_RADIUS,
                S_DECEPT_CENTER,
                S_DECEPT_RADIUS,
                S_DECEPT_FLOOR,
            ),
            1.0,
        );
        assert_close(
            s_decept(
                S_DECEPT_CENTER + S_DECEPT_RADIUS,
                S_DECEPT_CENTER,
                S_DECEPT_RADIUS,
                S_DECEPT_FLOOR,
            ),
            1.0,
        );
        assert_close(
            s_decept(
                S_DECEPT_CENTER - S_DECEPT_RADIUS / 2.0,
                S_DECEPT_CENTER,
                S_DECEPT_RADIUS,
                S_DECEPT_FLOOR,
            ),
            0.5,
        );
        assert_close(
            s_decept(
                S_DECEPT_CENTER + S_DECEPT_RADIUS / 2.0,
                S_DECEPT_CENTER,
                S_DECEPT_RADIUS,
                S_DECEPT_FLOOR,
            ),
            0.5,
        );
        assert_close(s_decept(0.0, 0.35, 0.001, 0.05), 0.050_000_000_000_029_354);
        assert_close(s_decept(1.0, 0.35, 0.001, 0.05), 0.049_999_999_999_958_63);

        assert_close(
            b_param(0.25, 0.0, B_PARAM_A, B_PARAM_B, B_PARAM_C),
            0.972_654_947_412_285_5,
        );
        assert_close(b_param(0.25, 0.5, B_PARAM_A, B_PARAM_B, B_PARAM_C), 0.25);
        assert_tiny_relative_close(
            b_param(0.25, 1.0, B_PARAM_A, B_PARAM_B, B_PARAM_C),
            7.888_609_052_210_118e-31,
        );
        assert_eq!(
            b_param(0.0, 1.0, B_PARAM_A, B_PARAM_B, B_PARAM_C).to_bits(),
            0.0_f64.to_bits(),
        );
        assert_close(
            b_param(0.99, 1.0, B_PARAM_A, B_PARAM_B, B_PARAM_C),
            0.605_006_067_137_536_4,
        );
        assert_eq!(
            b_param(1.0, 1.0, B_PARAM_A, B_PARAM_B, B_PARAM_C).to_bits(),
            1.0_f64.to_bits(),
        );
        assert_close(r_nonsep(&[0.37], 1), 0.37);
        assert_close(r_nonsep(&[0.1, 0.4, 0.9], 3), 0.766_666_666_666_666_7);
        assert_close(
            r_nonsep(&[0.1, 0.2, 0.4, 0.7, 0.9], 5),
            0.713_333_333_333_333_3,
        );
    }

    #[test]
    fn canonical_wfg9_multimodal_bias_95_anchors_match_the_toolkit() {
        assert_close(s_multi_with_parameters(0.0, 30.0, 95.0, 0.35), 1.0);
        assert_close(s_multi_with_parameters(0.35, 30.0, 95.0, 0.35), 0.0);
        assert_close(s_multi_with_parameters(1.0, 30.0, 95.0, 0.35), 1.0);
        assert_close(
            s_multi_with_parameters(0.17, 30.0, 95.0, 0.35),
            0.273_397_480_864_936_37,
        );
        assert_ne!(
            s_multi_with_parameters(0.17, 30.0, 95.0, 0.35).to_bits(),
            s_multi(0.17).to_bits(),
        );
        assert_close(
            s_multi_with_parameters(0.73, 30.0, 95.0, 0.35),
            0.340_027_377_732_080_7,
        );
        assert_ne!(
            s_multi_with_parameters(0.73, 30.0, 95.0, 0.35).to_bits(),
            s_multi(0.73).to_bits(),
        );
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
    fn pinned_reference_probe_matches_wfg5_full_pipeline() {
        // Frozen from an independent direct f64 port of the corrected
        // normalized equations. Every coordinate is an interior deceptive
        // probe and every objective group is asymmetric.
        let evaluation = Wfg5::new(4, 6, 4)
            .unwrap()
            .evaluate_normalized(&[0.12, 0.25, 0.38, 0.51, 0.64, 0.77, 0.09, 0.31, 0.53, 0.75])
            .unwrap();

        assert_slice_close(
            evaluation.transformed(),
            &[
                0.376_647_564_469_933_3,
                0.730_515_759_312_329_3,
                0.957_550_077_041_600_6,
                0.767_257_318_952_224,
                0.576_964_560_862_847_5,
                0.386_671_802_773_470_9,
                0.294_985_673_352_457_3,
                0.893_839_541_547_281_3,
                0.737_981_510_015_396_8,
                0.415_947_611_710_298_15,
            ],
        );
        assert_slice_close(
            evaluation.reduced(),
            &[
                0.553_581_661_891_131_3,
                0.862_403_697_996_912_3,
                0.481_818_181_818_159_2,
                0.585_688_584_156_358_4,
            ],
        );
        assert_slice_close(evaluation.positioned(), evaluation.reduced());
        assert_slice_close(
            evaluation.shape(),
            &[
                0.512_409_246_208_475,
                0.542_546_918_388_73,
                0.163_855_316_365_618_98,
                0.645_159_701_969_735_5,
            ],
        );
        assert_slice_close(
            evaluation.objectives(),
            &[
                1.610_507_076_573_308_5,
                2.755_876_257_711_278_8,
                1.568_820_482_350_072_2,
                5.746_966_199_914_242,
            ],
        );
    }

    #[test]
    fn pinned_reference_probe_matches_wfg6_full_pipeline() {
        // M=4, k=9, l=5 forces A=3 in each position reduction and odd A=5
        // in the distinct distance reduction, exposing pair/even/fixed-A
        // mutants in the same frozen trace.
        let evaluation = Wfg6::new(4, 9, 5)
            .unwrap()
            .evaluate_normalized(&[
                0.07, 0.19, 0.31, 0.46, 0.58, 0.73, 0.14, 0.67, 0.92, 0.02, 0.27, 0.44, 0.69, 0.96,
            ])
            .unwrap();

        assert_slice_close(
            evaluation.transformed(),
            &[
                0.07,
                0.19,
                0.31,
                0.46,
                0.58,
                0.73,
                0.14,
                0.67,
                0.92,
                0.942_857_142_857_142_8,
                0.228_571_428_571_428_48,
                0.138_461_538_461_538_5,
                0.523_076_923_076_923,
                0.938_461_538_461_538_5,
            ],
        );
        assert_slice_close(
            evaluation.reduced(),
            &[
                0.254_999_999_999_999_95,
                0.474_999_999_999_999_9,
                0.808_333_333_333_333_5,
                0.803_076_923_076_923_1,
            ],
        );
        assert_slice_close(evaluation.positioned(), evaluation.reduced());
        assert_slice_close(
            evaluation.shape(),
            &[
                0.252_777_739_904_556_76,
                0.078_489_574_532_045_62,
                0.286_332_678_195_076_85,
                0.920_845_480_141_026_3,
            ],
        );
        assert_slice_close(
            evaluation.objectives(),
            &[
                1.308_632_402_886_036_8,
                1.117_035_221_205_105_6,
                2.521_072_992_247_384_5,
                8.169_840_764_205_134,
            ],
        );
    }

    #[test]
    fn pinned_reference_probe_matches_wfg7_full_pipeline() {
        // Suffix means straddle 0.5 so both b_param exponent branches remain
        // interior and every conditioned position coordinate is observable.
        let evaluation = Wfg7::new(4, 6, 4)
            .unwrap()
            .evaluate_normalized(&[0.12, 0.25, 0.38, 0.51, 0.64, 0.77, 0.09, 0.31, 0.53, 0.75])
            .unwrap();

        assert_slice_close(
            evaluation.transformed(),
            &[
                0.135_933_152_438_588_94,
                0.251_703_991_513_018_6,
                0.098_057_458_291_359_09,
                0.189_538_890_886_922_36,
                0.645_622_781_437_027_8,
                0.802_211_677_331_793_9,
                0.742_857_142_857_142_9,
                0.114_285_714_285_714_24,
                0.276_923_076_923_077,
                0.615_384_615_384_615_4,
            ],
        );
        assert_slice_close(
            evaluation.reduced(),
            &[
                0.193_818_571_975_803_77,
                0.143_798_174_589_140_74,
                0.723_917_229_384_410_8,
                0.437_362_637_362_637_36,
            ],
        );
        assert_slice_close(evaluation.positioned(), evaluation.reduced());
        assert_slice_close(
            evaluation.shape(),
            &[
                0.060_921_741_434_431_45,
                0.028_211_043_903_512_478,
                0.292_153_307_271_704_35,
                0.954_012_119_143_502,
            ],
        );
        assert_slice_close(
            evaluation.objectives(),
            &[
                0.559_206_120_231_500_3,
                0.550_206_812_976_687_3,
                2.190_282_480_992_863,
                8.069_459_590_510_654,
            ],
        );
    }

    #[test]
    fn wfg7_conditioning_uses_original_strict_suffixes() {
        let problem = Wfg7::new(4, 6, 4).unwrap();
        let input = [0.12, 0.25, 0.38, 0.51, 0.64, 0.77, 0.09, 0.31, 0.53, 0.75];
        let baseline = problem.evaluate_normalized(&input).unwrap();

        let mut second_changed = input;
        second_changed[1] = 0.35;
        let second_changed = problem.evaluate_normalized(&second_changed).unwrap();
        assert!(
            second_changed.transformed()[0] < baseline.transformed()[0],
            "the first conditioner must include the changed second coordinate"
        );
        assert!(
            second_changed.transformed()[1] > baseline.transformed()[1],
            "the changed base must retain its original strict-suffix conditioner"
        );
        assert_slice_bits_eq(
            &baseline.transformed()[2..],
            &second_changed.transformed()[2..],
        );

        let mut last_changed = input;
        last_changed[9] = 0.65;
        let last_changed = problem.evaluate_normalized(&last_changed).unwrap();
        for (index, (&baseline_value, &changed_value)) in baseline.transformed()
            [..problem.position_parameters()]
            .iter()
            .zip(&last_changed.transformed()[..problem.position_parameters()])
            .enumerate()
        {
            assert_ne!(
                baseline_value.to_bits(),
                changed_value.to_bits(),
                "position coordinate {index} ignored the final suffix member"
            );
        }
        assert_slice_bits_eq(
            &baseline.transformed()[6..9],
            &last_changed.transformed()[6..9],
        );
        assert_ne!(
            baseline.transformed()[9].to_bits(),
            last_changed.transformed()[9].to_bits()
        );
    }

    #[test]
    fn wfg5_through_wfg7_admit_the_smallest_valid_dimension() {
        let input = [0.2, 0.7];
        let evaluations = [
            Wfg5::new(2, 1, 1)
                .unwrap()
                .evaluate_normalized(&input)
                .unwrap(),
            Wfg6::new(2, 1, 1)
                .unwrap()
                .evaluate_normalized(&input)
                .unwrap(),
            Wfg7::new(2, 1, 1)
                .unwrap()
                .evaluate_normalized(&input)
                .unwrap(),
        ];

        for evaluation in evaluations {
            assert_eq!(evaluation.transformed().len(), 2);
            assert_eq!(evaluation.reduced().len(), 2);
            assert_slice_close(evaluation.positioned(), evaluation.reduced());
            assert_eq!(evaluation.shape().len(), 2);
            assert_eq!(evaluation.objectives().len(), 2);
            assert!(
                evaluation
                    .objectives()
                    .iter()
                    .all(|value| value.is_finite())
            );
        }
    }

    #[test]
    fn wfg5_through_wfg7_repeated_evaluations_are_bitwise_identical() {
        let input = [0.12, 0.25, 0.38, 0.51, 0.64, 0.77, 0.09, 0.31, 0.53, 0.75];

        let problem = Wfg5::new(4, 6, 4).unwrap();
        let first = problem.evaluate_normalized(&input).unwrap();
        let second = problem.evaluate_normalized(&input).unwrap();
        assert_evaluation_bits_eq(&first, &second);
        assert_slice_bits_eq(&first.clone().into_objectives(), second.objectives());

        let problem = Wfg6::new(4, 6, 4).unwrap();
        let first = problem.evaluate_normalized(&input).unwrap();
        let second = problem.evaluate_normalized(&input).unwrap();
        assert_evaluation_bits_eq(&first, &second);
        assert_slice_bits_eq(&first.clone().into_objectives(), second.objectives());

        let problem = Wfg7::new(4, 6, 4).unwrap();
        let first = problem.evaluate_normalized(&input).unwrap();
        let second = problem.evaluate_normalized(&input).unwrap();
        assert_evaluation_bits_eq(&first, &second);
        assert_slice_bits_eq(&first.clone().into_objectives(), second.objectives());
    }

    #[test]
    fn pinned_reference_probe_matches_wfg8_full_pipeline() {
        // Frozen from an independent direct f64 port of the corrected
        // normalized equations at the pinned jMetal revision. Prefix means
        // [.5, .44, .49875, .54, .505] cross both b_param branches.
        let evaluation = Wfg8::new(4, 6, 5)
            .unwrap()
            .evaluate_normalized(&[
                0.12, 0.24, 0.36, 0.48, 0.84, 0.96, 0.08, 0.91, 0.87, 0.19, 0.73,
            ])
            .unwrap();

        assert_slice_close(
            evaluation.transformed(),
            &[
                0.12,
                0.24,
                0.36,
                0.48,
                0.84,
                0.96,
                0.771_428_571_428_571_4,
                0.877_152_197_225_494_5,
                0.800_456_750_371_557,
                0.999_192_021_745_912_9,
                0.424_120_392_160_531_7,
            ],
        );
        assert_slice_close(
            evaluation.reduced(),
            &[0.18, 0.42, 0.9, 0.774_469_986_586_413_5],
        );
        assert_slice_close(evaluation.positioned(), evaluation.reduced());
        assert_slice_close(
            evaluation.shape(),
            &[
                0.168_890_377_004_149_13,
                0.026_749_607_838_002_145,
                0.220_446_220_845_134_77,
                0.960_293_685_676_943_1,
            ],
        );
        assert_slice_close(
            evaluation.objectives(),
            &[
                1.112_250_740_594_711_7,
                0.881_468_417_938_422_1,
                2.097_147_311_657_222_3,
                8.456_819_472_001_959,
            ],
        );
    }

    #[test]
    fn pinned_reference_probe_matches_wfg9_full_pipeline() {
        // M=4, k=9, l=5 combines strict-suffix bias on all but the final
        // coordinate, deceptive positions, B=95 multimodal distance, and
        // nonseparable A=3/A=5 reductions in one independent frozen trace.
        let evaluation = Wfg9::new(4, 9, 5)
            .unwrap()
            .evaluate_normalized(&[
                0.07, 0.19, 0.31, 0.46, 0.58, 0.73, 0.14, 0.67, 0.92, 0.02, 0.27, 0.44, 0.69, 0.96,
            ])
            .unwrap();

        assert_slice_close(
            evaluation.transformed(),
            &[
                0.249_936_027_976_360_4,
                0.089_311_334_172_644_36,
                0.066_006_191_427_145_03,
                0.101_235_059_537_476_36,
                0.260_128_051_086_625_3,
                0.814_545_522_483_265_2,
                0.050_000_917_378_806_41,
                0.306_294_634_436_505_95,
                0.161_810_791_293_283_24,
                1.0,
                0.999_999_999_983_413_2,
                0.999_999_999_989_152_6,
                0.999_999_790_065_218_1,
                0.880_237_676_103_566,
            ],
        );
        assert_slice_close(
            evaluation.reduced(),
            &[
                0.190_162_149_962_168_51,
                0.671_525_080_815_087,
                0.257_213_535_223_232_26,
                0.389_222_459_800_570_6,
            ],
        );
        assert_slice_close(evaluation.positioned(), evaluation.reduced());
        assert_slice_close(
            evaluation.shape(),
            &[
                0.100_629_811_296_251_55,
                0.235_362_840_713_073_62,
                0.145_192_651_346_676,
                0.955_718_090_382_763_3,
            ],
        );
        assert_slice_close(
            evaluation.objectives(),
            &[
                0.590_482_082_393_073_6,
                1.330_673_822_652_865,
                1.260_378_367_880_626_6,
                8.034_967_182_862_678,
            ],
        );
    }

    #[test]
    fn wfg8_conditioning_uses_original_strict_prefixes() {
        let problem = Wfg8::new(4, 6, 5).unwrap();
        let input = [
            0.12, 0.24, 0.36, 0.48, 0.84, 0.96, 0.08, 0.91, 0.87, 0.19, 0.73,
        ];
        let baseline = problem.evaluate_normalized(&input).unwrap();

        let mut first_changed = input;
        first_changed[0] = 0.22;
        let first_changed = problem.evaluate_normalized(&first_changed).unwrap();
        assert_ne!(
            baseline.transformed()[0].to_bits(),
            first_changed.transformed()[0].to_bits()
        );
        assert_slice_bits_eq(
            &baseline.transformed()[1..6],
            &first_changed.transformed()[1..6],
        );
        for (index, (&baseline_value, &changed_value)) in baseline.transformed()[6..]
            .iter()
            .zip(&first_changed.transformed()[6..])
            .enumerate()
        {
            assert_ne!(
                baseline_value.to_bits(),
                changed_value.to_bits(),
                "distance coordinate {} ignored the first prefix member",
                index + 6
            );
        }

        let mut middle_changed = input;
        middle_changed[7] = 0.71;
        let middle_changed = problem.evaluate_normalized(&middle_changed).unwrap();
        assert_slice_bits_eq(
            &baseline.transformed()[..7],
            &middle_changed.transformed()[..7],
        );
        for (index, (&baseline_value, &changed_value)) in baseline.transformed()[7..]
            .iter()
            .zip(&middle_changed.transformed()[7..])
            .enumerate()
        {
            assert_ne!(
                baseline_value.to_bits(),
                changed_value.to_bits(),
                "coordinate {} escaped its forward prefix cone",
                index + 7
            );
        }

        let mut last_changed = input;
        last_changed[10] = 0.63;
        let last_changed = problem.evaluate_normalized(&last_changed).unwrap();
        assert_slice_bits_eq(
            &baseline.transformed()[..10],
            &last_changed.transformed()[..10],
        );
        assert_ne!(
            baseline.transformed()[10].to_bits(),
            last_changed.transformed()[10].to_bits()
        );

        let mut sum_preserved = input;
        sum_preserved[0] += 0.05;
        sum_preserved[1] -= 0.05;
        let sum_preserved = problem.evaluate_normalized(&sum_preserved).unwrap();
        assert_ne!(
            baseline.transformed()[0].to_bits(),
            sum_preserved.transformed()[0].to_bits()
        );
        assert_ne!(
            baseline.transformed()[1].to_bits(),
            sum_preserved.transformed()[1].to_bits()
        );
        assert_slice_bits_eq(
            &baseline.transformed()[2..6],
            &sum_preserved.transformed()[2..6],
        );
        assert_slice_close(
            &baseline.transformed()[6..],
            &sum_preserved.transformed()[6..],
        );
    }

    #[test]
    fn wfg9_conditioning_uses_original_strict_suffixes() {
        let problem = Wfg9::new(4, 9, 5).unwrap();
        let input = [
            0.07, 0.19, 0.31, 0.46, 0.58, 0.73, 0.14, 0.67, 0.92, 0.02, 0.27, 0.44, 0.69, 0.96,
        ];
        let baseline = problem.evaluate_normalized(&input).unwrap();

        let mut second_changed = input;
        second_changed[1] = 0.29;
        let second_changed = problem.evaluate_normalized(&second_changed).unwrap();
        assert_ne!(
            baseline.transformed()[0].to_bits(),
            second_changed.transformed()[0].to_bits()
        );
        assert_ne!(
            baseline.transformed()[1].to_bits(),
            second_changed.transformed()[1].to_bits()
        );
        assert_slice_bits_eq(
            &baseline.transformed()[2..],
            &second_changed.transformed()[2..],
        );

        let mut final_changed = input;
        final_changed[13] = 0.86;
        let final_changed = problem.evaluate_normalized(&final_changed).unwrap();
        for (index, (&baseline_value, &changed_value)) in baseline
            .transformed()
            .iter()
            .zip(final_changed.transformed())
            .enumerate()
        {
            assert_ne!(
                baseline_value.to_bits(),
                changed_value.to_bits(),
                "coordinate {index} ignored the final suffix member"
            );
        }

        let mut sum_preserved = input;
        sum_preserved[12] += 0.05;
        sum_preserved[13] -= 0.05;
        let sum_preserved = problem.evaluate_normalized(&sum_preserved).unwrap();
        assert_slice_close(
            &baseline.transformed()[..12],
            &sum_preserved.transformed()[..12],
        );
        assert_ne!(
            baseline.transformed()[12].to_bits(),
            sum_preserved.transformed()[12].to_bits()
        );
        assert_ne!(
            baseline.transformed()[13].to_bits(),
            sum_preserved.transformed()[13].to_bits()
        );
        assert_close(
            baseline.transformed()[13],
            s_multi_with_parameters(input[13], 30.0, 95.0, 0.35),
        );
    }

    #[test]
    fn pinned_minimum_dimension_traces_cover_wfg8_and_wfg9() {
        let input = [0.2, 0.7];
        let wfg8 = Wfg8::new(2, 1, 1)
            .unwrap()
            .evaluate_normalized(&input)
            .unwrap();
        assert_slice_close(wfg8.transformed(), &[0.2, 0.789_749_348_925_109_4]);
        assert_slice_close(wfg8.reduced(), wfg8.transformed());
        assert_slice_close(wfg8.positioned(), wfg8.reduced());
        assert_slice_close(
            wfg8.shape(),
            &[0.309_016_994_374_947_4, 0.951_056_516_295_153_5],
        );
        assert_slice_close(
            wfg8.objectives(),
            &[1.407_783_337_675_004_2, 4.593_975_414_105_723],
        );

        let wfg9 = Wfg9::new(2, 1, 1)
            .unwrap()
            .evaluate_normalized(&input)
            .unwrap();
        assert_slice_close(
            wfg9.transformed(),
            &[0.050_000_000_000_040_234, 0.303_400_357_978_124_2],
        );
        assert_slice_close(wfg9.reduced(), wfg9.transformed());
        assert_slice_close(wfg9.positioned(), wfg9.reduced());
        assert_slice_close(
            wfg9.shape(),
            &[0.078_459_095_727_907_95, 0.996_917_333_733_123],
        );
        assert_slice_close(
            wfg9.objectives(),
            &[0.460_318_549_433_940_1, 4.291_069_692_910_616],
        );
    }

    #[test]
    fn wfg8_and_wfg9_repeated_evaluations_are_bitwise_identical() {
        let wfg8_input = [
            0.12, 0.24, 0.36, 0.48, 0.84, 0.96, 0.08, 0.91, 0.87, 0.19, 0.73,
        ];
        let problem = Wfg8::new(4, 6, 5).unwrap();
        let first = problem.evaluate_normalized(&wfg8_input).unwrap();
        let second = problem.evaluate_normalized(&wfg8_input).unwrap();
        assert_evaluation_bits_eq(&first, &second);
        assert_slice_bits_eq(&first.clone().into_objectives(), second.objectives());

        let wfg9_input = [
            0.07, 0.19, 0.31, 0.46, 0.58, 0.73, 0.14, 0.67, 0.92, 0.02, 0.27, 0.44, 0.69, 0.96,
        ];
        let problem = Wfg9::new(4, 9, 5).unwrap();
        let first = problem.evaluate_normalized(&wfg9_input).unwrap();
        let second = problem.evaluate_normalized(&wfg9_input).unwrap();
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
