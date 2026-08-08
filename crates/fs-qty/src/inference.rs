//! Shared inference dimensional core (bead sj31i.7.1): the one
//! non-confusable dimensional algebra for as-built, assimilation, and OED
//! boundaries.
//!
//! The six-base [`Dims`] foundation and the semantic/affine scalar carriers
//! live in `lib.rs` and [`crate::semantic`]; this module freezes the ROLE
//! algebra on top of them: which dimensions a state slot, covariance entry,
//! information entry, observation, operator column, residual, cost, or
//! utility carry, and where affine absolute temperature is trapped instead
//! of silently linearized. Schemas carry canonical bytes and a
//! domain-separated content identity so two runs can compare structure by
//! hash instead of by faith.
//!
//! This core defines shared types only; it cannot mint domain evidence.

use crate::semantic::{QuantitySpec, SemanticType};
use crate::{Dims, semantic};

/// Maximum admitted state slots for one schema: matches the dense
/// assimilation envelope so a schema cannot name states no consumer can
/// update.
pub const MAX_STATE_SLOTS: usize = 256;

/// Canonical schema version for the inference dimensional core.
pub const INFERENCE_CORE_SCHEMA_VERSION: u8 = 1;

/// Domain for the content identity of one inference schema.
const SCHEMA_IDENTITY_DOMAIN: &str = "org.frankensim.fs-qty.inference-core.v1";

/// Wire prefix of the schema identity.
pub const SCHEMA_IDENTITY_PREFIX: &str = "fs-qty-inference-core:v1:";

/// Structured refusals of the inference dimensional core. Every mismatch
/// names both sides; nothing fails by default policy.
#[derive(Debug, Clone, PartialEq)]
pub enum InferenceError {
    /// A state schema must carry at least one slot.
    EmptyStateSchema,
    /// A state schema exceeds the admitted dense envelope.
    StateSlotLimit {
        /// Requested slots.
        slots: usize,
        /// Maximum admitted slots.
        max: usize,
    },
    /// A slot index outside the schema.
    SlotOutOfRange {
        /// Requested slot.
        slot: usize,
        /// Schema slot count.
        slots: usize,
    },
    /// Dimension exponent arithmetic overflowed the checked range.
    DimensionOverflow {
        /// Stable algebra stage.
        stage: &'static str,
    },
    /// A dimension admission compared incompatible quantities.
    DimensionMismatch {
        /// Role under admission.
        role: &'static str,
        /// Expected dimension rendering.
        expected: String,
        /// Supplied dimension rendering.
        actual: String,
    },
    /// A linear operator row cannot consume an affine absolute-temperature
    /// slot: affine points support differences, not scaling or summation.
    AffineSlotThroughLinearOperator {
        /// Offending state slot.
        slot: usize,
    },
    /// A decision measure mixed cost and utility semantics, or mixed a
    /// decision measure with a physical quantity.
    DecisionMeasureMismatch {
        /// Left operand rendering.
        left: &'static str,
        /// Right operand rendering.
        right: &'static str,
    },
    /// A decision measure received a non-finite value.
    NonFiniteDecisionValue,
    /// A canonical decoding carried a stale or unknown schema version.
    StaleSchemaVersion {
        /// Version this build decodes.
        expected: u8,
        /// Version found in the bytes.
        found: u8,
    },
    /// A canonical decoding had the wrong byte length or shape.
    MalformedSchemaBytes {
        /// Stable decode stage.
        stage: &'static str,
    },
}

impl core::fmt::Display for InferenceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::EmptyStateSchema => write!(f, "state schema must carry at least one slot"),
            Self::StateSlotLimit { slots, max } => {
                write!(f, "state schema slots {slots} exceed the admitted envelope {max}")
            }
            Self::SlotOutOfRange { slot, slots } => {
                write!(f, "slot {slot} is outside a {slots}-slot schema")
            }
            Self::DimensionOverflow { stage } => {
                write!(f, "dimension exponent overflow during {stage}")
            }
            Self::DimensionMismatch {
                role,
                expected,
                actual,
            } => write!(
                f,
                "dimension mismatch for {role}: expected {expected}, got {actual}"
            ),
            Self::AffineSlotThroughLinearOperator { slot } => write!(
                f,
                "state slot {slot} is affine absolute temperature; a linear operator cannot scale or sum affine points"
            ),
            Self::DecisionMeasureMismatch { left, right } => {
                write!(f, "cannot combine {left} with {right}")
            }
            Self::NonFiniteDecisionValue => {
                write!(f, "decision measure values must be finite")
            }
            Self::StaleSchemaVersion { expected, found } => {
                write!(f, "schema version {found} is not decodable by version {expected}")
            }
            Self::MalformedSchemaBytes { stage } => {
                write!(f, "malformed canonical schema bytes during {stage}")
            }
        }
    }
}

impl std::error::Error for InferenceError {}

/// One state-slot schema: an exact [`QuantitySpec`] (dimensional or
/// semantic, including the affine temperature kinds). The newtype exists so
/// role APIs cannot be called with a bare spec by accident.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct SlotSchema {
    spec: QuantitySpec,
}

impl SlotSchema {
    /// Admit one slot schema.
    #[must_use]
    pub const fn new(spec: QuantitySpec) -> Self {
        Self { spec }
    }

    /// The admitted quantity schema.
    #[must_use]
    pub const fn spec(self) -> QuantitySpec {
        self.spec
    }

    /// Six-base dimensions of this slot.
    #[must_use]
    pub fn dims(self) -> Dims {
        match self.spec {
            QuantitySpec::Dimensional(dims) => dims,
            QuantitySpec::Semantic(semantic) => semantic.expected_dims(),
        }
    }

    /// True when the slot is affine absolute temperature, which a linear
    /// operator row must refuse.
    #[must_use]
    pub fn is_affine_absolute_temperature(self) -> bool {
        matches!(
            self.spec,
            QuantitySpec::Semantic(semantic)
                if semantic.kind() == semantic::QuantityKind::AbsoluteTemperature
        )
    }

    /// The same quantity as a temperature difference, for operator
    /// admission after an explicit reference subtraction. Linear slots map
    /// to themselves; absolute temperature maps to its difference kind.
    ///
    /// # Errors
    /// This conversion is total; `Result` keeps the admission boundary
    /// uniform for callers that chain fallible steps.
    pub fn as_difference(self) -> Result<Self, InferenceError> {
        let mapped = match self.spec {
            QuantitySpec::Semantic(semantic_type)
                if semantic_type.kind() == semantic::QuantityKind::AbsoluteTemperature =>
            {
                QuantitySpec::Semantic(SemanticType::new(
                    semantic::QuantityKind::TemperatureDifference,
                    semantic_type.form(),
                ))
            }
            other => other,
        };
        Ok(Self { spec: mapped })
    }
}

/// An ordered, checked state-slot schema: the dimensional identity of one
/// state vector.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StateSchema {
    slots: Vec<SlotSchema>,
}

impl StateSchema {
    /// Admit an ordered slot list. Order is semantic: slot `i` owns
    /// covariance row/column `i` and operator column `i`.
    ///
    /// # Errors
    /// Returns [`InferenceError`] for an empty or oversized slot list.
    pub fn try_new(slots: Vec<SlotSchema>) -> Result<Self, InferenceError> {
        if slots.is_empty() {
            return Err(InferenceError::EmptyStateSchema);
        }
        if slots.len() > MAX_STATE_SLOTS {
            return Err(InferenceError::StateSlotLimit {
                slots: slots.len(),
                max: MAX_STATE_SLOTS,
            });
        }
        Ok(Self { slots })
    }

    /// Slot count.
    #[must_use]
    pub fn len(&self) -> usize {
        self.slots.len()
    }

    /// True when the schema carries no slots (never constructible through
    /// [`StateSchema::try_new`]).
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.slots.is_empty()
    }

    /// One admitted slot.
    ///
    /// # Errors
    /// Returns [`InferenceError::SlotOutOfRange`] outside the schema.
    pub fn slot(&self, slot: usize) -> Result<SlotSchema, InferenceError> {
        self.slots
            .get(slot)
            .copied()
            .ok_or(InferenceError::SlotOutOfRange {
                slot,
                slots: self.slots.len(),
            })
    }

    /// Read-only slot slice in declaration order.
    #[must_use]
    pub fn slots(&self) -> &[SlotSchema] {
        &self.slots
    }

    /// Admit a value schema into a state slot, naming both sides on
    /// refusal.
    ///
    /// # Errors
    /// Returns [`InferenceError::DimensionMismatch`] when the supplied
    /// schema's dimensions differ from the slot's.
    pub fn admit_value(
        &self,
        slot: usize,
        value: QuantitySpec,
    ) -> Result<(), InferenceError> {
        let slot_schema = self.slot(slot)?;
        let expected = slot_schema.dims();
        let actual = match value {
            QuantitySpec::Dimensional(dims) => dims,
            QuantitySpec::Semantic(semantic) => semantic.expected_dims(),
        };
        if expected != actual {
            return Err(InferenceError::DimensionMismatch {
                role: "state value",
                expected: expected.unit_string(),
                actual: actual.unit_string(),
            });
        }
        Ok(())
    }

    /// Canonical bytes: version byte, slot count, then each slot's
    /// [`QuantitySpec::canonical_bytes`] in declaration order.
    #[must_use]
    pub fn canonical_bytes(&self) -> Vec<u8> {
        let mut bytes = Vec::with_capacity(5 + self.slots.len() * 12);
        bytes.push(INFERENCE_CORE_SCHEMA_VERSION);
        bytes.extend_from_slice(&(self.slots.len() as u32).to_le_bytes());
        for slot in &self.slots {
            bytes.extend_from_slice(&slot.spec.canonical_bytes());
        }
        bytes
    }

    /// Decode canonical bytes, refusing stale versions and malformed
    /// shapes.
    ///
    /// # Errors
    /// Returns [`InferenceError::StaleSchemaVersion`] for any version other
    /// than [`INFERENCE_CORE_SCHEMA_VERSION`] and
    /// [`InferenceError::MalformedSchemaBytes`] for a truncated or ragged
    /// payload.
    pub fn from_canonical_bytes(bytes: &[u8]) -> Result<Self, InferenceError> {
        let (&version, rest) = bytes.split_first().ok_or(InferenceError::MalformedSchemaBytes {
            stage: "schema version",
        })?;
        if version != INFERENCE_CORE_SCHEMA_VERSION {
            return Err(InferenceError::StaleSchemaVersion {
                expected: INFERENCE_CORE_SCHEMA_VERSION,
                found: version,
            });
        }
        let (count_bytes, rest) = rest.split_first_chunk::<4>().ok_or(
            InferenceError::MalformedSchemaBytes {
                stage: "slot count",
            },
        )?;
        let count = u32::from_le_bytes(*count_bytes) as usize;
        if count == 0 || count > MAX_STATE_SLOTS {
            return Err(InferenceError::StateSlotLimit {
                slots: count,
                max: MAX_STATE_SLOTS,
            });
        }
        let (chunks, remainder) = rest.as_chunks::<12>();
        if !remainder.is_empty() || chunks.len() != count {
            return Err(InferenceError::MalformedSchemaBytes {
                stage: "slot payload",
            });
        }
        let mut slots = Vec::with_capacity(count);
        for chunk in chunks {
            let spec = QuantitySpec::from_canonical_bytes(chunk).map_err(|_| {
                InferenceError::MalformedSchemaBytes {
                    stage: "slot descriptor",
                }
            })?;
            slots.push(SlotSchema::new(spec));
        }
        Ok(Self { slots })
    }

    /// Domain-separated content identity
    /// `fs-qty-inference-core:v1:<64 lowercase hex>`.
    #[must_use]
    pub fn identity(&self) -> String {
        let mut hasher = fs_blake3::DomainHasher::new(SCHEMA_IDENTITY_DOMAIN);
        hasher.update(&self.canonical_bytes());
        format!("{SCHEMA_IDENTITY_PREFIX}{}", hasher.finalize())
    }
}

/// Covariance dimensional algebra over one state schema: entry `(i, j)`
/// carries `dims(i) * dims(j)`. Covariances of affine quantities are second
/// moments of deviations and carry no origin, so slot forms do not reach
/// entry dimensions.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CovarianceSchema {
    state: StateSchema,
}

impl CovarianceSchema {
    /// Derive the covariance schema for a state schema.
    #[must_use]
    pub const fn over(state: StateSchema) -> Self {
        Self { state }
    }

    /// The state schema this covariance varies over.
    #[must_use]
    pub const fn state(&self) -> &StateSchema {
        &self.state
    }

    /// Dimensions of covariance entry `(i, j)`: the checked product of the
    /// two slot dimensions.
    ///
    /// # Errors
    /// Returns [`InferenceError::SlotOutOfRange`] outside the schema and
    /// [`InferenceError::DimensionOverflow`] on exponent overflow.
    pub fn entry_dims(&self, row: usize, column: usize) -> Result<Dims, InferenceError> {
        let row_dims = self.state.slot(row)?.dims();
        let column_dims = self.state.slot(column)?.dims();
        row_dims
            .checked_plus(column_dims)
            .ok_or(InferenceError::DimensionOverflow {
                stage: "covariance entry",
            })
    }

    /// Dimensions of information-matrix (inverse covariance) entry
    /// `(i, j)`: the checked inverse of the covariance entry dimensions.
    ///
    /// # Errors
    /// See [`Self::entry_dims`].
    pub fn information_entry_dims(
        &self,
        row: usize,
        column: usize,
    ) -> Result<Dims, InferenceError> {
        self.entry_dims(row, column)?
            .checked_times(-1)
            .ok_or(InferenceError::DimensionOverflow {
                stage: "information entry",
            })
    }
}

/// Observation dimensional schema: the reading's quantity and the derived
/// noise-variance and residual dimensions.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct ObservationSchema {
    value: QuantitySpec,
}

impl ObservationSchema {
    /// Admit an observation value schema.
    #[must_use]
    pub const fn new(value: QuantitySpec) -> Self {
        Self { value }
    }

    /// The reading's quantity schema.
    #[must_use]
    pub const fn value(self) -> QuantitySpec {
        self.value
    }

    /// Reading dimensions.
    #[must_use]
    pub fn value_dims(self) -> Dims {
        match self.value {
            QuantitySpec::Dimensional(dims) => dims,
            QuantitySpec::Semantic(semantic) => semantic.expected_dims(),
        }
    }

    /// Noise-variance dimensions: the square of the reading dimensions.
    ///
    /// # Errors
    /// Returns [`InferenceError::DimensionOverflow`] on exponent overflow.
    pub fn noise_variance_dims(self) -> Result<Dims, InferenceError> {
        self.value_dims()
            .checked_times(2)
            .ok_or(InferenceError::DimensionOverflow {
                stage: "noise variance",
            })
    }

    /// Residual dimensions equal the reading dimensions.
    #[must_use]
    pub fn residual_dims(self) -> Dims {
        self.value_dims()
    }

    /// Admit a declared noise variance, naming both sides on refusal.
    ///
    /// # Errors
    /// Returns [`InferenceError::DimensionMismatch`] when the declared
    /// variance is not the squared reading dimension.
    pub fn admit_noise_variance(&self, variance: Dims) -> Result<(), InferenceError> {
        let expected = self.noise_variance_dims()?;
        if expected != variance {
            return Err(InferenceError::DimensionMismatch {
                role: "observation noise variance",
                expected: expected.unit_string(),
                actual: variance.unit_string(),
            });
        }
        Ok(())
    }
}

/// Observation-operator dimensional schema: one linear row mapping a state
/// schema to an observation. Column `i` carries
/// `output_dims - state_dims(i)`, so `sum_i h_i x_i` is dimensionally
/// homogeneous in the reading. Affine absolute-temperature slots are
/// refused: a linear map cannot scale or sum affine points.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OperatorSchema {
    output: ObservationSchema,
    state: StateSchema,
}

impl OperatorSchema {
    /// Admit an operator row from an observation schema over a state
    /// schema.
    ///
    /// # Errors
    /// Returns [`InferenceError::AffineSlotThroughLinearOperator`] when any
    /// state slot is affine absolute temperature.
    pub fn try_new(
        output: ObservationSchema,
        state: StateSchema,
    ) -> Result<Self, InferenceError> {
        for (slot, slot_schema) in state.slots().iter().enumerate() {
            if slot_schema.is_affine_absolute_temperature() {
                return Err(InferenceError::AffineSlotThroughLinearOperator { slot });
            }
        }
        Ok(Self { output, state })
    }

    /// The observation schema this row produces.
    #[must_use]
    pub const fn output(&self) -> ObservationSchema {
        self.output
    }

    /// The state schema this row consumes.
    #[must_use]
    pub const fn state(&self) -> &StateSchema {
        &self.state
    }

    /// Dimensions of operator column `i`.
    ///
    /// # Errors
    /// Returns [`InferenceError::SlotOutOfRange`] outside the schema and
    /// [`InferenceError::DimensionOverflow`] on exponent overflow.
    pub fn column_dims(&self, slot: usize) -> Result<Dims, InferenceError> {
        let output_dims = self.output.value_dims();
        let state_dims = self.state.slot(slot)?.dims();
        output_dims
            .checked_minus(state_dims)
            .ok_or(InferenceError::DimensionOverflow {
                stage: "operator column",
            })
    }

    /// Admit a declared operator column coefficient, naming both sides on
    /// refusal.
    ///
    /// # Errors
    /// Returns [`InferenceError::DimensionMismatch`] when the declared
    /// coefficient does not carry `output - state(slot)` dimensions.
    pub fn admit_column(&self, slot: usize, coefficient: Dims) -> Result<(), InferenceError> {
        let expected = self.column_dims(slot)?;
        if expected != coefficient {
            return Err(InferenceError::DimensionMismatch {
                role: "operator column coefficient",
                expected: expected.unit_string(),
                actual: coefficient.unit_string(),
            });
        }
        Ok(())
    }
}

/// Decision-measure semantics: cost and utility are decision quantities,
/// not physical quantities. Same-kind algebra is checked; cross-kind
/// combination and any physical/dimensional mixing refuse by type.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum DecisionMeasure {
    /// An objective cost (lower is better).
    Cost(f64),
    /// A decision utility (higher is better).
    Utility(f64),
}

impl DecisionMeasure {
    /// Admit a cost.
    ///
    /// # Errors
    /// Returns [`InferenceError::NonFiniteDecisionValue`] for NaN or
    /// infinite values.
    pub fn cost(value: f64) -> Result<Self, InferenceError> {
        if !value.is_finite() {
            return Err(InferenceError::NonFiniteDecisionValue);
        }
        Ok(Self::Cost(value))
    }

    /// Admit a utility.
    ///
    /// # Errors
    /// Returns [`InferenceError::NonFiniteDecisionValue`] for NaN or
    /// infinite values.
    pub fn utility(value: f64) -> Result<Self, InferenceError> {
        if !value.is_finite() {
            return Err(InferenceError::NonFiniteDecisionValue);
        }
        Ok(Self::Utility(value))
    }

    /// The measure kind name for diagnostics.
    #[must_use]
    pub const fn kind_name(self) -> &'static str {
        match self {
            Self::Cost(_) => "cost",
            Self::Utility(_) => "utility",
        }
    }

    /// The carried value.
    #[must_use]
    pub const fn value(self) -> f64 {
        match self {
            Self::Cost(value) | Self::Utility(value) => value,
        }
    }

    /// Same-kind addition; cross-kind combination refuses, naming both
    /// sides. There is no favorable default.
    ///
    /// # Errors
    /// Returns [`InferenceError::DecisionMeasureMismatch`] for mixed
    /// cost/utility or non-finite results.
    pub fn checked_add(self, other: Self) -> Result<Self, InferenceError> {
        match (self, other) {
            (Self::Cost(left), Self::Cost(right)) => Self::cost(left + right),
            (Self::Utility(left), Self::Utility(right)) => Self::utility(left + right),
            (left, right) => Err(InferenceError::DecisionMeasureMismatch {
                left: left.kind_name(),
                right: right.kind_name(),
            }),
        }
    }

    /// Same-kind subtraction; cross-kind combination refuses.
    ///
    /// # Errors
    /// Returns [`InferenceError::DecisionMeasureMismatch`] for mixed
    /// cost/utility or non-finite results.
    pub fn checked_sub(self, other: Self) -> Result<Self, InferenceError> {
        match (self, other) {
            (Self::Cost(left), Self::Cost(right)) => Self::cost(left - right),
            (Self::Utility(left), Self::Utility(right)) => Self::utility(left - right),
            (left, right) => Err(InferenceError::DecisionMeasureMismatch {
                left: left.kind_name(),
                right: right.kind_name(),
            }),
        }
    }

    /// Scaling by a finite scalar preserves the measure kind.
    ///
    /// # Errors
    /// Returns [`InferenceError::NonFiniteDecisionValue`] for non-finite
    /// factors or results.
    pub fn checked_scale(self, factor: f64) -> Result<Self, InferenceError> {
        if !factor.is_finite() {
            return Err(InferenceError::NonFiniteDecisionValue);
        }
        match self {
            Self::Cost(value) => Self::cost(value * factor),
            Self::Utility(value) => Self::utility(value * factor),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::semantic::{QuantityKind, SemanticType, ValueForm};
    use crate::units;

    const LENGTH: Dims = Dims([1, 0, 0, 0, 0, 0]);
    const TIME: Dims = Dims([0, 0, 1, 0, 0, 0]);
    const TEMPERATURE: Dims = Dims([0, 0, 0, 1, 0, 0]);

    fn slot(dims: Dims) -> SlotSchema {
        SlotSchema::new(QuantitySpec::dimensional(dims))
    }

    fn absolute_temperature_slot() -> SlotSchema {
        SlotSchema::new(QuantitySpec::Semantic(SemanticType::new(
            QuantityKind::AbsoluteTemperature,
            ValueForm::Static,
        )))
    }

    #[test]
    fn g0_state_schema_admission_and_bounds() {
        assert!(StateSchema::try_new(Vec::new()).is_err());
        let one = StateSchema::try_new(vec![slot(LENGTH)]).expect("one slot");
        assert_eq!(one.len(), 1);
        let max = StateSchema::try_new(vec![slot(LENGTH); MAX_STATE_SLOTS]).expect("max slots");
        assert_eq!(max.len(), MAX_STATE_SLOTS);
        let over = StateSchema::try_new(vec![slot(LENGTH); MAX_STATE_SLOTS + 1]);
        assert!(matches!(
            over,
            Err(InferenceError::StateSlotLimit { .. })
        ));
        assert!(matches!(
            one.slot(1),
            Err(InferenceError::SlotOutOfRange { slot: 1, slots: 1 })
        ));
    }

    #[test]
    fn g0_covariance_and_information_entry_algebra() {
        let state = StateSchema::try_new(vec![slot(LENGTH), slot(TIME), slot(Dims::NONE)])
            .expect("state");
        let covariance = CovarianceSchema::over(state);
        assert_eq!(
            covariance.entry_dims(0, 0).expect("entry"),
            Dims([2, 0, 0, 0, 0, 0])
        );
        assert_eq!(
            covariance.entry_dims(0, 1).expect("entry"),
            Dims([1, 0, 1, 0, 0, 0])
        );
        assert_eq!(
            covariance.entry_dims(0, 2).expect("entry"),
            LENGTH,
            "dimensionless slot leaves the other slot's dimensions"
        );
        assert_eq!(
            covariance.information_entry_dims(0, 1).expect("info"),
            Dims([-1, 0, -1, 0, 0, 0])
        );
        assert!(covariance.entry_dims(3, 0).is_err());
    }

    #[test]
    fn g0_observation_noise_and_residual_dimensions() {
        let velocity = ObservationSchema::new(QuantitySpec::dimensional(Dims([1, 0, -1, 0, 0, 0])));
        assert_eq!(
            velocity.noise_variance_dims().expect("variance dims"),
            Dims([2, 0, -2, 0, 0, 0])
        );
        assert_eq!(velocity.residual_dims(), Dims([1, 0, -1, 0, 0, 0]));
        assert!(velocity.admit_noise_variance(Dims([2, 0, -2, 0, 0, 0])).is_ok());
        let refusal = velocity
            .admit_noise_variance(Dims([2, 0, 0, 0, 0, 0]))
            .expect_err("wrong variance dims refuse");
        let rendered = refusal.to_string();
        assert!(rendered.contains("expected m^2·s^-2"), "{rendered}");
        assert!(rendered.contains("got m^2"), "{rendered}");
    }

    #[test]
    fn g0_operator_column_algebra_and_affine_trap() {
        let state = StateSchema::try_new(vec![slot(LENGTH), slot(TIME)]).expect("state");
        let output = ObservationSchema::new(QuantitySpec::dimensional(Dims([1, 0, -1, 0, 0, 0])));
        let operator = OperatorSchema::try_new(output, state).expect("operator");
        assert_eq!(operator.column_dims(0).expect("column"), Dims([0, 0, -1, 0, 0, 0]));
        assert_eq!(operator.column_dims(1).expect("column"), Dims([1, 0, -2, 0, 0, 0]));
        assert!(operator.admit_column(0, Dims([0, 0, -1, 0, 0, 0])).is_ok());
        assert!(operator.admit_column(0, LENGTH).is_err());

        let affine_state =
            StateSchema::try_new(vec![absolute_temperature_slot()]).expect("affine state");
        let refusal = OperatorSchema::try_new(output, affine_state)
            .expect_err("affine slot cannot pass a linear operator");
        assert!(matches!(
            refusal,
            InferenceError::AffineSlotThroughLinearOperator { slot: 0 }
        ));
    }

    #[test]
    fn g0_affine_difference_conversion_enables_operator_admission() {
        let affine_state =
            StateSchema::try_new(vec![absolute_temperature_slot()]).expect("affine state");
        let as_difference = affine_state
            .slot(0)
            .expect("slot")
            .as_difference()
            .expect("difference conversion");
        assert!(!as_difference.is_affine_absolute_temperature());
        assert_eq!(as_difference.dims(), TEMPERATURE);
        let difference_state = StateSchema::try_new(vec![as_difference]).expect("delta state");
        let output = ObservationSchema::new(QuantitySpec::Semantic(SemanticType::new(
            QuantityKind::TemperatureDifference,
            ValueForm::Static,
        )));
        assert!(OperatorSchema::try_new(output, difference_state).is_ok());
    }

    #[test]
    fn g0_decision_measures_refuse_cross_kind_mixing() {
        let cost_a = DecisionMeasure::cost(1.5).expect("cost");
        let cost_b = DecisionMeasure::cost(2.5).expect("cost");
        let utility = DecisionMeasure::utility(3.0).expect("utility");
        assert_eq!(
            cost_a.checked_add(cost_b).expect("same kind").value().to_bits(),
            4.0_f64.to_bits()
        );
        let refusal = cost_a
            .checked_add(utility)
            .expect_err("cost+utility refuses");
        let rendered = refusal.to_string();
        assert!(rendered.contains("cost"), "{rendered}");
        assert!(rendered.contains("utility"), "{rendered}");
        assert!(DecisionMeasure::cost(f64::NAN).is_err());
        assert!(cost_a.checked_scale(f64::INFINITY).is_err());
        assert_eq!(
            utility.checked_sub(DecisionMeasure::utility(1.0).expect("u")).expect("same kind").value().to_bits(),
            2.0_f64.to_bits()
        );
    }

    #[test]
    fn g0_canonical_round_trip_and_stale_version_refusal() {
        let state = StateSchema::try_new(vec![
            slot(LENGTH),
            absolute_temperature_slot(),
            slot(Dims::NONE),
        ])
        .expect("state");
        let bytes = state.canonical_bytes();
        let decoded = StateSchema::from_canonical_bytes(&bytes).expect("round trip");
        assert_eq!(decoded, state);

        let mut stale = bytes.clone();
        stale[0] = INFERENCE_CORE_SCHEMA_VERSION + 1;
        assert!(matches!(
            StateSchema::from_canonical_bytes(&stale),
            Err(InferenceError::StaleSchemaVersion { .. })
        ));
        assert!(matches!(
            StateSchema::from_canonical_bytes(&bytes[..bytes.len() - 1]),
            Err(InferenceError::MalformedSchemaBytes { .. })
        ));
    }

    #[test]
    fn g5_identity_is_deterministic_and_schema_sensitive() {
        let left = StateSchema::try_new(vec![slot(LENGTH), slot(TIME)]).expect("left");
        let right = StateSchema::try_new(vec![slot(TIME), slot(LENGTH)]).expect("right");
        assert_eq!(left.identity(), left.identity());
        assert_ne!(left.identity(), right.identity());
        assert!(left.identity().starts_with(SCHEMA_IDENTITY_PREFIX));
        let mutated = StateSchema::try_new(vec![slot(Dims([2, 0, 0, 0, 0, 0])), slot(TIME)])
            .expect("mutated");
        assert_ne!(left.identity(), mutated.identity());
    }

    #[test]
    fn g3_slot_permutation_changes_identity_but_preserves_algebra() {
        let forward = StateSchema::try_new(vec![slot(LENGTH), slot(TIME)]).expect("forward");
        let reversed = StateSchema::try_new(vec![slot(TIME), slot(LENGTH)]).expect("reversed");
        assert_ne!(forward.identity(), reversed.identity());
        let forward_covariance = CovarianceSchema::over(forward.clone());
        let reversed_covariance = CovarianceSchema::over(reversed);
        assert_eq!(
            forward_covariance.entry_dims(0, 1).expect("entry"),
            reversed_covariance.entry_dims(1, 0).expect("entry")
        );
    }

    #[test]
    fn g3_unit_rescaling_preserves_dimensions() {
        // mm and m are the same dimension with different scales; the core
        // compares dimensions, so rescaled values admit identically.
        let meters = units::meters(2.0);
        let millimeters = units::millimeters(2000.0);
        assert_eq!(meters.value().to_bits(), millimeters.value().to_bits());
        let state = StateSchema::try_new(vec![slot(LENGTH)]).expect("state");
        assert!(state.admit_value(0, QuantitySpec::dimensional(LENGTH)).is_ok());
        let refusal = state
            .admit_value(0, QuantitySpec::dimensional(TIME))
            .expect_err("time is not length");
        let rendered = refusal.to_string();
        assert!(rendered.contains("expected m, got s"), "{rendered}");
    }
}
