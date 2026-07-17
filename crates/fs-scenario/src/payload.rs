//! Typed, versioned scenario payload values.
//!
//! This module is deliberately a data boundary, not a solver API.  Every
//! numeric sample carries either an exact six-base [`Dims`] contract or a
//! validated [`SemanticType`]; basis, frame, orientation parity, and
//! continuity/reset semantics travel with that contract.  Constructors reject
//! malformed shapes and non-finite or dimensionally inconsistent data before
//! it can enter scenario IR.
//!
//! [`canonical_payload_bytes`] and [`decode_payload`] define a compact,
//! deterministic V1 envelope.  The decoder is bounded before every allocation
//! and reconstructs values through the same validating constructors.  This is
//! a payload codec only: embedding it in the scenario or FrankenScript IR is a
//! one-way integration owned by those higher-level modules; [`crate::ir`]
//! embeds these exact bytes in canonical scenario v2 artifacts.

use crate::frame::FrameId;
use core::fmt;
use fs_qty::chemistry::{ChemistryError, SpeciesId};
use fs_qty::semantic::{
    AngleDomain, Composition, CompositionBasis, PhasorAmplitude, PhasorQty, QuantityKind,
    SemanticError, SemanticQty, SemanticType, StrainBasis, StrainComponent, ValueForm,
};
use fs_qty::{Dims, QtyAny};

/// Canonical payload-envelope version written by this build.
pub const PAYLOAD_WIRE_VERSION: u16 = 1;

/// Maximum aggregate count accepted by a validated in-memory payload.
pub const MAX_PAYLOAD_ITEMS: usize = 1_048_576;

/// Maximum byte length of any identifier in the payload algebra.
pub const MAX_PAYLOAD_ID_BYTES: usize = 128;

/// Conservative hard envelope bound implied by the admitted item/id ceilings.
///
/// The closed V1 codec uses fewer than `MAX_PAYLOAD_ID_BYTES + 32` bytes per
/// charged aggregate item, including length fields and fixed-width numeric
/// records. This bound keeps default decoding closed over every payload
/// admitted by the constructors while remaining explicit.
pub const MAX_PAYLOAD_WIRE_BYTES: usize = MAX_PAYLOAD_ITEMS * (MAX_PAYLOAD_ID_BYTES + 32) + 4_096;

const PAYLOAD_MAGIC: &[u8; 8] = b"FSPAYLD\0";
const TIME_DIMS: Dims = Dims([0, 0, 1, 0, 0, 0]);

/// A canonical ASCII identifier used for bases, events, artifacts, fields,
/// components, and ports.
///
/// Identifiers begin with an ASCII alphanumeric byte.  Remaining bytes may
/// additionally contain `_`, `-`, `.`, `:`, or `/`.  The intentionally small
/// grammar avoids Unicode-normalization and whitespace aliases in hashes.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct PayloadId(String);

impl PayloadId {
    /// Validate and retain an identifier.
    ///
    /// # Errors
    /// Returns [`PayloadError::InvalidIdentifier`] for empty, overlong, or
    /// non-canonical text.
    pub fn new(value: impl Into<String>) -> Result<Self, PayloadError> {
        let value = value.into();
        if value.is_empty() {
            return Err(PayloadError::InvalidIdentifier {
                value,
                reason: "identifier must not be empty",
            });
        }
        if value.len() > MAX_PAYLOAD_ID_BYTES {
            return Err(PayloadError::InvalidIdentifier {
                value,
                reason: "identifier exceeds the 128-byte limit",
            });
        }
        let bytes = value.as_bytes();
        if !bytes[0].is_ascii_alphanumeric() {
            return Err(PayloadError::InvalidIdentifier {
                value,
                reason: "identifier must begin with an ASCII letter or digit",
            });
        }
        if !bytes.iter().all(|byte| {
            byte.is_ascii_alphanumeric() || matches!(byte, b'_' | b'-' | b'.' | b':' | b'/')
        }) {
            return Err(PayloadError::InvalidIdentifier {
                value,
                reason: "identifier contains whitespace, non-ASCII, or unsupported punctuation",
            });
        }
        Ok(Self(value))
    }

    /// Borrow the canonical text.
    #[must_use]
    pub fn as_str(&self) -> &str {
        &self.0
    }
}

impl fmt::Display for PayloadId {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(&self.0)
    }
}

/// Contract applied to every numeric component of one payload.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum QuantityContract {
    /// Exact six-base SI exponent vector, without a stronger semantic claim.
    Dimensions(Dims),
    /// Dimension plus the sealed `fs-qty` physical kind and value form.
    Semantic(SemanticType),
    /// Components carry distinct contracts.  Admitted only by
    /// [`CharacteristicState`], whose named descriptors hold those contracts.
    Heterogeneous,
}

impl QuantityContract {
    /// Homogeneous dimension vector, or `None` for a heterogeneous carrier.
    #[must_use]
    pub const fn dims(self) -> Option<Dims> {
        match self {
            Self::Dimensions(dims) => Some(dims),
            Self::Semantic(semantic) => Some(semantic.expected_dims()),
            Self::Heterogeneous => None,
        }
    }
}

/// Behavior under an orientation-reversing frame transform.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum OrientationParity {
    /// Even parity: the value does not change sign under reflection.
    Even,
    /// Odd parity: the value changes sign under reflection.
    Odd,
}

/// Temporal reference-origin semantics carried by a payload.
///
/// This tag describes whether the reference origin is reset; it does not claim
/// mathematical continuity of the sample path. In particular, a `StepLeft`
/// table may contain jumps while retaining one uninterrupted reference origin.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum ReferenceSemantics {
    /// One uninterrupted reference origin is used for the complete history.
    Continuous,
    /// The reference origin is reset when the named scenario event occurs.
    ResetAtEvent(PayloadId),
}

/// Metadata shared by every payload variant.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct PayloadMeta {
    contract: QuantityContract,
    basis_key: PayloadId,
    frame: FrameId,
    orientation: OrientationParity,
    reference: ReferenceSemantics,
}

impl PayloadMeta {
    /// Construct complete payload metadata and validate semantic kind/form.
    pub fn new(
        contract: QuantityContract,
        basis_key: PayloadId,
        frame: FrameId,
        orientation: OrientationParity,
        reference: ReferenceSemantics,
    ) -> Result<Self, PayloadError> {
        validate_contract(contract, "payload metadata")?;
        Ok(Self {
            contract,
            basis_key,
            frame,
            orientation,
            reference,
        })
    }

    /// Quantity contract.
    #[must_use]
    pub const fn contract(&self) -> QuantityContract {
        self.contract
    }

    /// Required canonical basis name.
    #[must_use]
    pub const fn basis_key(&self) -> &PayloadId {
        &self.basis_key
    }

    /// Frame in which components are expressed.
    #[must_use]
    pub const fn frame(&self) -> FrameId {
        self.frame
    }

    /// Orientation parity.
    #[must_use]
    pub const fn orientation(&self) -> OrientationParity {
        self.orientation
    }

    /// Continuity/reset behavior.
    #[must_use]
    pub const fn reference(&self) -> &ReferenceSemantics {
        &self.reference
    }
}

/// Supported distribution families.
///
/// `Normal` and `Uniform` carry exactly two shape-compatible parameter values;
/// `Empirical` carries one or more samples.  This layer proves shape, units,
/// and finiteness only.  Correlation and sampling algorithms belong to the
/// stochastic owner and are deliberately not claimed here.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum DistributionFamily {
    /// Mean and standard deviation.
    Normal,
    /// Inclusive lower and upper endpoints.
    Uniform,
    /// Canonically ordered empirical support samples.
    Empirical,
}

/// Interpolation used between admitted table samples.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum TableInterpolation {
    /// Hold the left sample until the next sample time.
    StepLeft,
    /// Component-wise linear interpolation.
    Linear,
}

/// Explicit behavior when evaluating outside a table's time domain.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum OutsideDomainPolicy {
    /// Refuse evaluation outside the closed sample interval.
    Refuse,
    /// Clamp to the nearest endpoint sample.
    Clamp,
    /// Repeat with the table span as period.  Requires at least two samples.
    Periodic,
}

/// Fixed value, time table, or distribution parameters for one payload shape.
#[derive(Debug, Clone, PartialEq)]
pub enum SampleSource<T> {
    /// One fixed value.
    Fixed(T),
    /// Strictly increasing coherent-SI seconds and one value per time.
    Table {
        /// Sample times, each with dimensions of time.
        times: Vec<QtyAny>,
        /// Values paired positionally with `times`.
        values: Vec<T>,
        /// Between-sample interpolation.
        interpolation: TableInterpolation,
        /// Outside-domain behavior.
        outside: OutsideDomainPolicy,
    },
    /// Closed family and its shape-compatible parameters/support.
    Distribution {
        /// Distribution family.
        family: DistributionFamily,
        /// Normal/uniform parameters or empirical support values.
        parameters: Vec<T>,
    },
}

impl<T> SampleSource<T> {
    /// Construct a fixed source.
    #[must_use]
    pub const fn fixed(value: T) -> Self {
        Self::Fixed(value)
    }

    /// Construct a time table after checking its structural time axis.
    ///
    /// Payload constructors repeat these checks so manually assembled enum
    /// variants cannot bypass admission.
    pub fn table(
        times: Vec<QtyAny>,
        values: Vec<T>,
        interpolation: TableInterpolation,
        outside: OutsideDomainPolicy,
    ) -> Result<Self, PayloadError> {
        validate_table_axis(&times, values.len(), outside)?;
        Ok(Self::Table {
            times,
            values,
            interpolation,
            outside,
        })
    }

    /// Construct distribution parameters/support after checking family arity.
    pub fn distribution(
        family: DistributionFamily,
        parameters: Vec<T>,
    ) -> Result<Self, PayloadError> {
        validate_distribution_arity(family, parameters.len())?;
        Ok(Self::Distribution { family, parameters })
    }

    /// Time axis for a table source.
    #[must_use]
    pub fn times(&self) -> Option<&[QtyAny]> {
        match self {
            Self::Table { times, .. } => Some(times),
            Self::Fixed(_) | Self::Distribution { .. } => None,
        }
    }
}

/// Scalar payload.
#[derive(Debug, Clone, PartialEq)]
pub struct ScalarPayload {
    meta: PayloadMeta,
    source: SampleSource<QtyAny>,
}

impl ScalarPayload {
    /// Construct and validate every scalar sample.
    pub fn new(meta: PayloadMeta, source: SampleSource<QtyAny>) -> Result<Self, PayloadError> {
        require_homogeneous(&meta, "scalar payload")?;
        validate_source(&source, |quantity, index| {
            validate_quantity(*quantity, meta.contract, "scalar payload", index)
        })?;
        validate_scalar_distribution(&source, "scalar distribution")?;
        enforce_aggregate_limit(
            aggregate_source_items(&source, |_| 0, aggregate_meta_items(&meta)?).ok_or(
                PayloadError::CountOverflow {
                    context: "scalar payload items",
                },
            )?,
            MAX_PAYLOAD_ITEMS,
        )?;
        Ok(Self { meta, source })
    }

    /// Payload metadata.
    #[must_use]
    pub const fn meta(&self) -> &PayloadMeta {
        &self.meta
    }

    /// Value source.
    #[must_use]
    pub const fn source(&self) -> &SampleSource<QtyAny> {
        &self.source
    }
}

/// Vector payload with a stable component count across all samples.
#[derive(Debug, Clone, PartialEq)]
pub struct VectorPayload {
    meta: PayloadMeta,
    source: SampleSource<Vec<QtyAny>>,
    components: usize,
}

impl VectorPayload {
    /// Construct a nonempty vector source and prove shape/contract agreement.
    pub fn new(meta: PayloadMeta, source: SampleSource<Vec<QtyAny>>) -> Result<Self, PayloadError> {
        Self::new_with_item_limit(meta, source, MAX_PAYLOAD_ITEMS)
    }

    /// Construct with a smaller aggregate-item ceiling.
    ///
    /// This is useful to admission planners that reserve part of a scenario's
    /// global budget before accepting one nested payload.
    pub fn new_with_item_limit(
        meta: PayloadMeta,
        source: SampleSource<Vec<QtyAny>>,
        max_items: usize,
    ) -> Result<Self, PayloadError> {
        require_homogeneous(&meta, "vector payload")?;
        let mut components = None;
        validate_source(&source, |sample, sample_index| {
            if sample.is_empty() {
                return Err(PayloadError::Empty {
                    context: "vector components",
                });
            }
            check_item_count("vector components", sample.len())?;
            match components {
                None => components = Some(sample.len()),
                Some(expected) if expected != sample.len() => {
                    return Err(PayloadError::ShapeMismatch {
                        context: "vector sample",
                        expected,
                        actual: sample.len(),
                    });
                }
                Some(_) => {}
            }
            for (component, quantity) in sample.iter().enumerate() {
                validate_quantity(
                    *quantity,
                    meta.contract,
                    "vector component",
                    flattened_index(sample_index, sample.len(), component),
                )?;
            }
            Ok(())
        })?;
        validate_vector_distribution(&source, "vector distribution")?;
        enforce_aggregate_limit(
            aggregate_source_items(&source, |sample| sample.len(), aggregate_meta_items(&meta)?)
                .ok_or(PayloadError::CountOverflow {
                    context: "vector payload items",
                })?,
            max_items.min(MAX_PAYLOAD_ITEMS),
        )?;
        let components = components.ok_or(PayloadError::Empty {
            context: "vector components",
        })?;
        Ok(Self {
            meta,
            source,
            components,
        })
    }

    /// Payload metadata.
    #[must_use]
    pub const fn meta(&self) -> &PayloadMeta {
        &self.meta
    }

    /// Stable component count.
    #[must_use]
    pub const fn components(&self) -> usize {
        self.components
    }

    /// Value source.
    #[must_use]
    pub const fn source(&self) -> &SampleSource<Vec<QtyAny>> {
        &self.source
    }
}

/// Row-major tensor payload with a stable rectangular shape.
#[derive(Debug, Clone, PartialEq)]
pub struct TensorPayload {
    meta: PayloadMeta,
    rows: usize,
    columns: usize,
    source: SampleSource<Vec<QtyAny>>,
}

impl TensorPayload {
    /// Construct a tensor and validate every row-major sample.
    pub fn new(
        meta: PayloadMeta,
        rows: usize,
        columns: usize,
        source: SampleSource<Vec<QtyAny>>,
    ) -> Result<Self, PayloadError> {
        require_homogeneous(&meta, "tensor payload")?;
        if rows == 0 || columns == 0 {
            return Err(PayloadError::Empty {
                context: "tensor shape",
            });
        }
        let expected = rows
            .checked_mul(columns)
            .ok_or(PayloadError::CountOverflow {
                context: "tensor shape",
            })?;
        check_item_count("tensor components", expected)?;
        validate_source(&source, |sample, sample_index| {
            if sample.len() != expected {
                return Err(PayloadError::ShapeMismatch {
                    context: "tensor sample",
                    expected,
                    actual: sample.len(),
                });
            }
            for (component, quantity) in sample.iter().enumerate() {
                validate_quantity(
                    *quantity,
                    meta.contract,
                    "tensor component",
                    flattened_index(sample_index, expected, component),
                )?;
            }
            Ok(())
        })?;
        validate_vector_distribution(&source, "tensor distribution")?;
        enforce_aggregate_limit(
            aggregate_source_items(&source, |sample| sample.len(), aggregate_meta_items(&meta)?)
                .ok_or(PayloadError::CountOverflow {
                    context: "tensor payload items",
                })?,
            MAX_PAYLOAD_ITEMS,
        )?;
        Ok(Self {
            meta,
            rows,
            columns,
            source,
        })
    }

    /// Payload metadata.
    #[must_use]
    pub const fn meta(&self) -> &PayloadMeta {
        &self.meta
    }

    /// Row count.
    #[must_use]
    pub const fn rows(&self) -> usize {
        self.rows
    }

    /// Column count.
    #[must_use]
    pub const fn columns(&self) -> usize {
        self.columns
    }

    /// Row-major value source.
    #[must_use]
    pub const fn source(&self) -> &SampleSource<Vec<QtyAny>> {
        &self.source
    }
}

/// Complex phasor payload. The required semantic Peak/RMS contract and paired
/// `fs-qty` carrier prevent samples or real/imaginary components from drifting
/// in physical kind or amplitude convention.
#[derive(Debug, Clone, PartialEq)]
pub struct ComplexPhasorPayload {
    meta: PayloadMeta,
    source: SampleSource<PhasorQty>,
}

impl ComplexPhasorPayload {
    /// Construct and validate a phasor source against the metadata contract.
    pub fn new(meta: PayloadMeta, source: SampleSource<PhasorQty>) -> Result<Self, PayloadError> {
        require_homogeneous(&meta, "complex phasor payload")?;
        let QuantityContract::Semantic(semantic) = meta.contract else {
            return Err(PayloadError::PhasorSemanticContractRequired);
        };
        if !matches!(semantic.form(), ValueForm::Peak | ValueForm::Rms)
            || !semantic.kind().admits_phasor()
        {
            return Err(PayloadError::PhasorSemanticContractRequired);
        }
        validate_source(&source, |phasor, index| {
            validate_phasor(*phasor, meta.contract, index)
        })?;
        validate_phasor_distribution(&source)?;
        enforce_aggregate_limit(
            aggregate_source_items(&source, |_| 0, aggregate_meta_items(&meta)?).ok_or(
                PayloadError::CountOverflow {
                    context: "phasor payload items",
                },
            )?,
            MAX_PAYLOAD_ITEMS,
        )?;
        Ok(Self { meta, source })
    }

    /// Payload metadata.
    #[must_use]
    pub const fn meta(&self) -> &PayloadMeta {
        &self.meta
    }

    /// Phasor source.
    #[must_use]
    pub const fn source(&self) -> &SampleSource<PhasorQty> {
        &self.source
    }
}

/// One species quantity in a canonical species-axis bundle.
#[derive(Debug, Clone, PartialEq)]
pub struct SpeciesValue {
    species: SpeciesId,
    quantity: QtyAny,
}

impl SpeciesValue {
    /// Pair a validated `fs-qty` species identifier with a six-base quantity.
    #[must_use]
    pub const fn new(species: SpeciesId, quantity: QtyAny) -> Self {
        Self { species, quantity }
    }

    /// Species identifier.
    #[must_use]
    pub const fn species(&self) -> &SpeciesId {
        &self.species
    }

    /// Runtime six-base quantity.
    #[must_use]
    pub const fn quantity(&self) -> QtyAny {
        self.quantity
    }
}

/// Species bundle whose axis is nonempty, unique, sorted, and stable across
/// every table sample or distribution parameter.
#[derive(Debug, Clone, PartialEq)]
pub struct SpeciesBundle {
    meta: PayloadMeta,
    source: SampleSource<Vec<SpeciesValue>>,
    species: Vec<SpeciesId>,
}

impl SpeciesBundle {
    /// Construct a species bundle with a canonical, shape-stable axis.
    pub fn new(
        meta: PayloadMeta,
        source: SampleSource<Vec<SpeciesValue>>,
    ) -> Result<Self, PayloadError> {
        Self::new_with_item_limit(meta, source, MAX_PAYLOAD_ITEMS)
    }

    /// Construct while charging the retained canonical-axis cache to a caller
    /// supplied aggregate-item ceiling.
    pub fn new_with_item_limit(
        meta: PayloadMeta,
        source: SampleSource<Vec<SpeciesValue>>,
        max_items: usize,
    ) -> Result<Self, PayloadError> {
        require_homogeneous(&meta, "species bundle")?;
        let axis_items = first_source_value(&source)
            .map(Vec::len)
            .ok_or(PayloadError::Empty {
                context: "species bundle",
            })?;
        let axis_cache_items = axis_items
            .checked_mul(2)
            .ok_or(PayloadError::CountOverflow {
                context: "species-axis cache items",
            })?;
        let retained_items = aggregate_source_items(
            &source,
            |sample| sample.len().saturating_mul(2),
            aggregate_meta_items(&meta)?,
        )
        .and_then(|count| count.checked_add(axis_cache_items))
        .ok_or(PayloadError::CountOverflow {
            context: "species payload items",
        })?;
        enforce_aggregate_limit(retained_items, max_items.min(MAX_PAYLOAD_ITEMS))?;
        let mut species: Option<Vec<SpeciesId>> = None;
        validate_source(&source, |sample, sample_index| {
            if sample.is_empty() {
                return Err(PayloadError::Empty {
                    context: "species bundle",
                });
            }
            check_item_count("species bundle", sample.len())?;
            for entry in sample {
                let length = entry.species.as_str().len();
                if length > MAX_PAYLOAD_ID_BYTES {
                    return Err(PayloadError::IdentifierByteLimit {
                        limit: MAX_PAYLOAD_ID_BYTES,
                        actual: length,
                    });
                }
            }
            for pair in sample.windows(2) {
                if pair[0].species >= pair[1].species {
                    return Err(PayloadError::NonCanonicalSpeciesAxis {
                        previous: try_copy_text(
                            pair[0].species.as_str(),
                            "species-axis diagnostic",
                        )?,
                        next: try_copy_text(pair[1].species.as_str(), "species-axis diagnostic")?,
                    });
                }
            }
            match &species {
                None => {
                    let mut axis = Vec::new();
                    axis.try_reserve_exact(sample.len()).map_err(|_| {
                        PayloadError::AllocationRefused {
                            context: "species axis",
                            count: sample.len(),
                        }
                    })?;
                    for entry in sample {
                        axis.push(SpeciesId::new(try_copy_text(
                            entry.species.as_str(),
                            "species-axis identity",
                        )?)?);
                    }
                    species = Some(axis);
                }
                Some(expected)
                    if expected.len() != sample.len()
                        || expected
                            .iter()
                            .zip(sample)
                            .any(|(expected, entry)| expected != &entry.species) =>
                {
                    return Err(PayloadError::SpeciesAxisMismatch {
                        sample: sample_index,
                    });
                }
                Some(_) => {}
            }
            for (entry_index, entry) in sample.iter().enumerate() {
                validate_quantity(
                    entry.quantity,
                    meta.contract,
                    "species quantity",
                    flattened_index(sample_index, sample.len(), entry_index),
                )?;
            }
            Ok(())
        })?;
        validate_species_distribution(&meta, &source)?;
        validate_composition_samples(&meta, &source)?;
        let species = species.ok_or(PayloadError::Empty {
            context: "species bundle",
        })?;
        Ok(Self {
            meta,
            source,
            species,
        })
    }

    /// Payload metadata.
    #[must_use]
    pub const fn meta(&self) -> &PayloadMeta {
        &self.meta
    }

    /// Canonical species axis.
    #[must_use]
    pub fn species(&self) -> &[SpeciesId] {
        &self.species
    }

    /// Bundle source.
    #[must_use]
    pub const fn source(&self) -> &SampleSource<Vec<SpeciesValue>> {
        &self.source
    }
}

/// Direction of a characteristic component relative to a boundary.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum CharacteristicDirection {
    /// Characteristic enters the modeled domain.
    Incoming,
    /// Characteristic leaves the modeled domain.
    Outgoing,
    /// Zero-speed or tangential characteristic.
    Stationary,
}

/// Named characteristic component with its own physical contract.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct CharacteristicComponent {
    name: PayloadId,
    direction: CharacteristicDirection,
    contract: QuantityContract,
}

impl CharacteristicComponent {
    /// Construct one homogeneous characteristic component.
    pub fn new(
        name: PayloadId,
        direction: CharacteristicDirection,
        contract: QuantityContract,
    ) -> Result<Self, PayloadError> {
        validate_contract(contract, "characteristic component")?;
        if contract == QuantityContract::Heterogeneous {
            return Err(PayloadError::HeterogeneousContract {
                context: "characteristic component",
            });
        }
        Ok(Self {
            name,
            direction,
            contract,
        })
    }

    /// Stable component name.
    #[must_use]
    pub const fn name(&self) -> &PayloadId {
        &self.name
    }

    /// Boundary-relative direction.
    #[must_use]
    pub const fn direction(&self) -> CharacteristicDirection {
        self.direction
    }

    /// Component-specific quantity contract.
    #[must_use]
    pub const fn contract(&self) -> QuantityContract {
        self.contract
    }
}

/// Characteristic variables with an explicit direction for each stable
/// component.
#[derive(Debug, Clone, PartialEq)]
pub struct CharacteristicState {
    meta: PayloadMeta,
    components: Vec<CharacteristicComponent>,
    source: SampleSource<Vec<QtyAny>>,
}

impl CharacteristicState {
    /// Construct a nonempty characteristic state and validate all samples.
    pub fn new(
        meta: PayloadMeta,
        components: Vec<CharacteristicComponent>,
        source: SampleSource<Vec<QtyAny>>,
    ) -> Result<Self, PayloadError> {
        require_heterogeneous(&meta, "characteristic state")?;
        if components.is_empty() {
            return Err(PayloadError::Empty {
                context: "characteristic components",
            });
        }
        check_item_count("characteristic components", components.len())?;
        for pair in components.windows(2) {
            if pair[0].name >= pair[1].name {
                return Err(PayloadError::NonCanonicalComponentAxis {
                    previous: try_copy_text(
                        pair[0].name.as_str(),
                        "characteristic-axis diagnostic",
                    )?,
                    next: try_copy_text(pair[1].name.as_str(), "characteristic-axis diagnostic")?,
                });
            }
        }
        validate_source(&source, |sample, sample_index| {
            if sample.len() != components.len() {
                return Err(PayloadError::ShapeMismatch {
                    context: "characteristic sample",
                    expected: components.len(),
                    actual: sample.len(),
                });
            }
            for (component_index, (component, quantity)) in
                components.iter().zip(sample).enumerate()
            {
                validate_quantity(
                    *quantity,
                    component.contract,
                    "characteristic component",
                    flattened_index(sample_index, components.len(), component_index),
                )?;
            }
            Ok(())
        })?;
        validate_vector_distribution(&source, "characteristic distribution")?;
        enforce_aggregate_limit(
            aggregate_meta_items(&meta)?
                .checked_add(components.len().saturating_mul(2))
                .and_then(|count| aggregate_source_items(&source, |sample| sample.len(), count))
                .ok_or(PayloadError::CountOverflow {
                    context: "characteristic payload items",
                })?,
            MAX_PAYLOAD_ITEMS,
        )?;
        Ok(Self {
            meta,
            components,
            source,
        })
    }

    /// Payload metadata.
    #[must_use]
    pub const fn meta(&self) -> &PayloadMeta {
        &self.meta
    }

    /// Named component descriptors.
    #[must_use]
    pub fn components(&self) -> &[CharacteristicComponent] {
        &self.components
    }

    /// Characteristic source.
    #[must_use]
    pub const fn source(&self) -> &SampleSource<Vec<QtyAny>> {
        &self.source
    }
}

/// Reference to a field trace stored in a content-addressed artifact.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FieldTraceRef {
    meta: PayloadMeta,
    artifact: PayloadId,
    field: PayloadId,
}

impl FieldTraceRef {
    /// Construct a field-trace reference from canonical identifiers.
    pub fn new(
        meta: PayloadMeta,
        artifact: PayloadId,
        field: PayloadId,
    ) -> Result<Self, PayloadError> {
        require_homogeneous(&meta, "field trace reference")?;
        enforce_aggregate_limit(
            aggregate_meta_items(&meta)?
                .checked_add(2)
                .ok_or(PayloadError::CountOverflow {
                    context: "field trace identity items",
                })?,
            MAX_PAYLOAD_ITEMS,
        )?;
        Ok(Self {
            meta,
            artifact,
            field,
        })
    }

    /// Payload metadata.
    #[must_use]
    pub const fn meta(&self) -> &PayloadMeta {
        &self.meta
    }

    /// Content-addressed artifact key or stable ledger alias.
    #[must_use]
    pub const fn artifact(&self) -> &PayloadId {
        &self.artifact
    }

    /// Field name within the artifact.
    #[must_use]
    pub const fn field(&self) -> &PayloadId {
        &self.field
    }
}

/// Reference to a named component port.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PortRef {
    meta: PayloadMeta,
    component: PayloadId,
    port: PayloadId,
}

impl PortRef {
    /// Construct a port reference from canonical identifiers.
    pub fn new(
        meta: PayloadMeta,
        component: PayloadId,
        port: PayloadId,
    ) -> Result<Self, PayloadError> {
        require_homogeneous(&meta, "port reference")?;
        enforce_aggregate_limit(
            aggregate_meta_items(&meta)?
                .checked_add(2)
                .ok_or(PayloadError::CountOverflow {
                    context: "port identity items",
                })?,
            MAX_PAYLOAD_ITEMS,
        )?;
        Ok(Self {
            meta,
            component,
            port,
        })
    }

    /// Payload metadata.
    #[must_use]
    pub const fn meta(&self) -> &PayloadMeta {
        &self.meta
    }

    /// Component identifier.
    #[must_use]
    pub const fn component(&self) -> &PayloadId {
        &self.component
    }

    /// Port identifier.
    #[must_use]
    pub const fn port(&self) -> &PayloadId {
        &self.port
    }
}

/// Closed set of typed scenario payloads admitted by V1.
#[derive(Debug, Clone, PartialEq)]
pub enum Payload {
    /// Scalar quantity.
    Scalar(ScalarPayload),
    /// Fixed-width vector quantity.
    Vector(VectorPayload),
    /// Rectangular row-major tensor quantity.
    Tensor(TensorPayload),
    /// Complete paired complex phasor.
    ComplexPhasor(ComplexPhasorPayload),
    /// Canonically keyed chemical species quantities.
    SpeciesBundle(SpeciesBundle),
    /// Directed characteristic variables.
    CharacteristicState(CharacteristicState),
    /// External field-trace reference.
    FieldTraceRef(FieldTraceRef),
    /// Component-port reference.
    PortRef(PortRef),
}

/// Stable structural tag for [`Payload`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum PayloadKind {
    /// [`Payload::Scalar`].
    Scalar,
    /// [`Payload::Vector`].
    Vector,
    /// [`Payload::Tensor`].
    Tensor,
    /// [`Payload::ComplexPhasor`].
    ComplexPhasor,
    /// [`Payload::SpeciesBundle`].
    SpeciesBundle,
    /// [`Payload::CharacteristicState`].
    CharacteristicState,
    /// [`Payload::FieldTraceRef`].
    FieldTraceRef,
    /// [`Payload::PortRef`].
    PortRef,
}

/// Scenario-preflight footprint produced by one checkpointed payload walk.
pub(crate) struct PayloadValidationStats {
    /// Numeric scalar slots retained by the source representation.
    pub(crate) dynamic_scalars: usize,
    /// Aggregate bytes retained by canonical identities.
    pub(crate) identity_bytes: usize,
    /// Longest retained identity component.
    pub(crate) max_identity_component_bytes: usize,
    /// Number of identity components actually visited by the accounting walk.
    pub(crate) identity_components: usize,
}

/// Payload-accounting refusal with caller-owned cancellation propagation.
pub(crate) enum PayloadValidationStatsError<E> {
    /// A checked payload count overflowed or exceeded its closed bound.
    Payload(PayloadError),
    /// The caller refused continuation at a payload tile boundary.
    Checkpoint(E),
}

impl Payload {
    /// Closed structural payload kind.
    #[must_use]
    pub const fn kind(&self) -> PayloadKind {
        match self {
            Self::Scalar(_) => PayloadKind::Scalar,
            Self::Vector(_) => PayloadKind::Vector,
            Self::Tensor(_) => PayloadKind::Tensor,
            Self::ComplexPhasor(_) => PayloadKind::ComplexPhasor,
            Self::SpeciesBundle(_) => PayloadKind::SpeciesBundle,
            Self::CharacteristicState(_) => PayloadKind::CharacteristicState,
            Self::FieldTraceRef(_) => PayloadKind::FieldTraceRef,
            Self::PortRef(_) => PayloadKind::PortRef,
        }
    }

    /// Metadata shared by every variant.
    #[must_use]
    pub const fn meta(&self) -> &PayloadMeta {
        match self {
            Self::Scalar(value) => value.meta(),
            Self::Vector(value) => value.meta(),
            Self::Tensor(value) => value.meta(),
            Self::ComplexPhasor(value) => value.meta(),
            Self::SpeciesBundle(value) => value.meta(),
            Self::CharacteristicState(value) => value.meta(),
            Self::FieldTraceRef(value) => value.meta(),
            Self::PortRef(value) => value.meta(),
        }
    }

    /// Homogeneous six-base dimensions, absent for characteristic states.
    #[must_use]
    pub const fn homogeneous_dims(&self) -> Option<Dims> {
        self.meta().contract().dims()
    }

    /// Checked number of dynamic scalar values retained by this payload,
    /// including table coordinates as well as table values.
    pub fn bounded_dynamic_scalar_count(&self) -> Result<usize, PayloadError> {
        payload_dynamic_scalar_count(self)
    }

    /// Checked total bytes in every basis/event/component/species/reference id.
    pub fn identity_bytes(&self) -> Result<usize, PayloadError> {
        self.identity_stats().map(|(total, _)| total)
    }

    /// Longest individual identity component in bytes.
    pub fn max_identity_component_bytes(&self) -> Result<usize, PayloadError> {
        self.identity_stats().map(|(_, maximum)| maximum)
    }

    /// Checked aggregate and maximum identity byte counts in one traversal.
    pub fn identity_stats(&self) -> Result<(usize, usize), PayloadError> {
        payload_identity_stats(self)
    }

    /// Compute the complete scenario-validation footprint while polling at
    /// every source sample and retained identity component actually traversed.
    pub(crate) fn validation_stats_with_checkpoint<E>(
        &self,
        checkpoint: &mut impl FnMut(&'static str) -> Result<(), E>,
    ) -> Result<PayloadValidationStats, PayloadValidationStatsError<E>> {
        let dynamic_scalars = payload_dynamic_scalar_count_with_checkpoint(self, checkpoint)?;
        let (identity_bytes, max_identity_component_bytes, identity_components) =
            payload_identity_stats_with_checkpoint(self, checkpoint)?;
        Ok(PayloadValidationStats {
            dynamic_scalars,
            identity_bytes,
            max_identity_component_bytes,
            identity_components,
        })
    }
}

fn payload_dynamic_scalar_count(payload: &Payload) -> Result<usize, PayloadError> {
    let mut checkpoint = |_: &'static str| Ok::<(), core::convert::Infallible>(());
    match payload_dynamic_scalar_count_with_checkpoint(payload, &mut checkpoint) {
        Ok(count) => Ok(count),
        Err(PayloadValidationStatsError::Payload(error)) => Err(error),
        Err(PayloadValidationStatsError::Checkpoint(never)) => match never {},
    }
}

fn payload_dynamic_scalar_count_with_checkpoint<E>(
    payload: &Payload,
    checkpoint: &mut impl FnMut(&'static str) -> Result<(), E>,
) -> Result<usize, PayloadValidationStatsError<E>> {
    let count = match payload {
        Payload::Scalar(value) => sum_source_with_checkpoint(&value.source, |_| 1, checkpoint)?,
        Payload::Vector(value) => sum_source_with_checkpoint(&value.source, Vec::len, checkpoint)?,
        Payload::Tensor(value) => sum_source_with_checkpoint(&value.source, Vec::len, checkpoint)?,
        Payload::ComplexPhasor(value) => {
            sum_source_with_checkpoint(&value.source, |_| 2, checkpoint)?
        }
        Payload::SpeciesBundle(value) => {
            sum_source_with_checkpoint(&value.source, Vec::len, checkpoint)?
        }
        Payload::CharacteristicState(value) => {
            sum_source_with_checkpoint(&value.source, Vec::len, checkpoint)?
        }
        Payload::FieldTraceRef(_) | Payload::PortRef(_) => 0,
    };
    enforce_aggregate_limit(count, MAX_PAYLOAD_ITEMS)
        .map_err(PayloadValidationStatsError::Payload)?;
    Ok(count)
}

fn sum_source_with_checkpoint<T, E>(
    source: &SampleSource<T>,
    mut count_value: impl FnMut(&T) -> usize,
    checkpoint: &mut impl FnMut(&'static str) -> Result<(), E>,
) -> Result<usize, PayloadValidationStatsError<E>> {
    let mut total = 0_usize;
    match source {
        SampleSource::Fixed(value) => add_source_count(&mut total, count_value(value), checkpoint)?,
        SampleSource::Table { times, values, .. } => {
            total = total
                .checked_add(times.len())
                .ok_or(PayloadValidationStatsError::Payload(
                    PayloadError::CountOverflow {
                        context: "payload dynamic scalar count",
                    },
                ))?;
            for value in values {
                add_source_count(&mut total, count_value(value), checkpoint)?;
            }
        }
        SampleSource::Distribution { parameters, .. } => {
            for value in parameters {
                add_source_count(&mut total, count_value(value), checkpoint)?;
            }
        }
    }
    Ok(total)
}

fn add_source_count<E>(
    total: &mut usize,
    amount: usize,
    checkpoint: &mut impl FnMut(&'static str) -> Result<(), E>,
) -> Result<(), PayloadValidationStatsError<E>> {
    checkpoint("typed payload scalar accounting")
        .map_err(PayloadValidationStatsError::Checkpoint)?;
    *total = (*total)
        .checked_add(amount)
        .ok_or(PayloadValidationStatsError::Payload(
            PayloadError::CountOverflow {
                context: "payload dynamic scalar count",
            },
        ))?;
    Ok(())
}

fn payload_identity_stats(payload: &Payload) -> Result<(usize, usize), PayloadError> {
    let mut checkpoint = |_: &'static str| Ok::<(), core::convert::Infallible>(());
    match payload_identity_stats_with_checkpoint(payload, &mut checkpoint) {
        Ok((total, maximum, _)) => Ok((total, maximum)),
        Err(PayloadValidationStatsError::Payload(error)) => Err(error),
        Err(PayloadValidationStatsError::Checkpoint(never)) => match never {},
    }
}

fn payload_identity_stats_with_checkpoint<E>(
    payload: &Payload,
    checkpoint: &mut impl FnMut(&'static str) -> Result<(), E>,
) -> Result<(usize, usize, usize), PayloadValidationStatsError<E>> {
    let mut total = 0_usize;
    let mut maximum = 0_usize;
    let mut components = 0_usize;
    let mut add = |text: &str| -> Result<(), PayloadValidationStatsError<E>> {
        checkpoint("typed payload identity accounting")
            .map_err(PayloadValidationStatsError::Checkpoint)?;
        total = total
            .checked_add(text.len())
            .ok_or(PayloadValidationStatsError::Payload(
                PayloadError::CountOverflow {
                    context: "payload identity bytes",
                },
            ))?;
        maximum = maximum.max(text.len());
        components = components
            .checked_add(1)
            .ok_or(PayloadValidationStatsError::Payload(
                PayloadError::CountOverflow {
                    context: "payload identity component count",
                },
            ))?;
        Ok(())
    };
    add(payload.meta().basis_key.as_str())?;
    if let ReferenceSemantics::ResetAtEvent(event) = &payload.meta().reference {
        add(event.as_str())?;
    }
    match payload {
        Payload::SpeciesBundle(value) => {
            // The bundle retains a canonical species-axis cache in addition to
            // the source samples. Admission accounting must charge both owned
            // copies, just as the aggregate-item bound does.
            for species in &value.species {
                add(species.as_str())?;
            }
            match &value.source {
                SampleSource::Fixed(sample) => {
                    for entry in sample {
                        add(entry.species.as_str())?;
                    }
                }
                SampleSource::Table { values, .. } => {
                    for sample in values {
                        for entry in sample {
                            add(entry.species.as_str())?;
                        }
                    }
                }
                SampleSource::Distribution { parameters, .. } => {
                    for sample in parameters {
                        for entry in sample {
                            add(entry.species.as_str())?;
                        }
                    }
                }
            }
        }
        Payload::CharacteristicState(value) => {
            for component in &value.components {
                add(component.name.as_str())?;
            }
        }
        Payload::FieldTraceRef(value) => {
            add(value.artifact.as_str())?;
            add(value.field.as_str())?;
        }
        Payload::PortRef(value) => {
            add(value.component.as_str())?;
            add(value.port.as_str())?;
        }
        Payload::Scalar(_)
        | Payload::Vector(_)
        | Payload::Tensor(_)
        | Payload::ComplexPhasor(_) => {}
    }
    Ok((total, maximum, components))
}

/// Decode-allocation limits.  Counts are aggregate across nested vectors, not
/// merely per-field maxima.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PayloadDecodeLimits {
    /// Maximum complete envelope byte count.
    pub max_bytes: usize,
    /// Maximum aggregate strings/vector/table items.
    pub max_items: usize,
    /// Maximum bytes in one identifier.
    pub max_identifier_bytes: usize,
}

impl PayloadDecodeLimits {
    /// Conservative default for an admitted scenario payload.
    pub const DEFAULT: Self = Self {
        max_bytes: MAX_PAYLOAD_WIRE_BYTES,
        max_items: MAX_PAYLOAD_ITEMS,
        max_identifier_bytes: MAX_PAYLOAD_ID_BYTES,
    };
}

impl Default for PayloadDecodeLimits {
    fn default() -> Self {
        Self::DEFAULT
    }
}

/// Structured payload validation or canonical-codec refusal.
#[derive(Debug, Clone, PartialEq)]
pub enum PayloadError {
    /// Canonical identifier grammar was violated.
    InvalidIdentifier {
        /// Rejected text.
        value: String,
        /// Stable diagnostic.
        reason: &'static str,
    },
    /// A required aggregate was empty.
    Empty {
        /// Named aggregate.
        context: &'static str,
    },
    /// A count exceeded the admitted in-memory limit.
    TooManyItems {
        /// Named aggregate.
        context: &'static str,
        /// Admitted maximum.
        limit: usize,
        /// Supplied count.
        actual: usize,
    },
    /// Checked shape arithmetic overflowed.
    CountOverflow {
        /// Named aggregate.
        context: &'static str,
    },
    /// Parallel vectors or samples disagreed in shape.
    ShapeMismatch {
        /// Named aggregate.
        context: &'static str,
        /// Required count.
        expected: usize,
        /// Supplied count.
        actual: usize,
    },
    /// A table time was invalid.
    InvalidTableTime {
        /// Time index.
        index: usize,
        /// Stable diagnostic.
        reason: &'static str,
    },
    /// Distribution parameter/support count violated its family contract.
    DistributionArity {
        /// Distribution family.
        family: DistributionFamily,
        /// Required count description.
        expected: &'static str,
        /// Supplied count.
        actual: usize,
    },
    /// A semantic kind/form pair is not admitted by `fs-qty`.
    InvalidSemanticContract {
        /// Boundary validating the declaration.
        context: &'static str,
        /// Rejected descriptor.
        semantic: SemanticType,
    },
    /// A heterogeneous contract appeared outside its one admitted carrier.
    HeterogeneousContract {
        /// Boundary that requires a homogeneous quantity.
        context: &'static str,
    },
    /// A distribution family has no honest interpretation for this carrier.
    UnsupportedDistribution {
        /// Named carrier.
        context: &'static str,
        /// Rejected family.
        family: DistributionFamily,
    },
    /// Normal scale or uniform/empirical order was invalid.
    DistributionOrder {
        /// Named carrier.
        context: &'static str,
        /// Stable parameter/support index.
        index: usize,
        /// Stable diagnostic.
        reason: &'static str,
    },
    /// Numeric sample was not finite.
    NonFinite {
        /// Named value.
        context: &'static str,
        /// Stable flattened sample index.
        index: usize,
        /// Rejected value.
        value: f64,
    },
    /// Runtime dimensions disagreed with the declared contract.
    DimensionMismatch {
        /// Named value.
        context: &'static str,
        /// Stable flattened sample index.
        index: usize,
        /// Required dimensions.
        expected: Dims,
        /// Supplied dimensions.
        actual: Dims,
    },
    /// A semantic quantity constructor refused a sample.
    Semantic(SemanticError),
    /// A chemistry identifier constructor refused decoded bytes.
    Chemistry(ChemistryError),
    /// Complex phasors require a stable semantic kind and Peak/RMS convention.
    PhasorSemanticContractRequired,
    /// Phasor kind or peak/RMS convention disagreed with metadata.
    PhasorContractMismatch {
        /// Stable flattened sample index.
        index: usize,
        /// Declared semantic type when present.
        expected: Option<SemanticType>,
        /// Actual physical kind.
        actual_kind: QuantityKind,
        /// Actual peak/RMS convention.
        actual_amplitude: PhasorAmplitude,
    },
    /// Species axes were not strictly increasing and unique.
    NonCanonicalSpeciesAxis {
        /// Previous identifier.
        previous: String,
        /// Following non-increasing identifier.
        next: String,
    },
    /// A later species sample changed the canonical axis.
    SpeciesAxisMismatch {
        /// Stable sample index.
        sample: usize,
    },
    /// Characteristic descriptors were not strictly name-sorted and unique.
    NonCanonicalComponentAxis {
        /// Previous name.
        previous: String,
        /// Following non-increasing name.
        next: String,
    },
    /// Complete wire input exceeded its byte budget.
    ByteLimit {
        /// Admitted maximum.
        limit: usize,
        /// Supplied count.
        actual: usize,
    },
    /// Aggregate decoded items exceeded their cumulative budget.
    ItemLimit {
        /// Admitted maximum.
        limit: usize,
        /// Count requested after this field.
        actual: usize,
    },
    /// One decoded identifier exceeded its byte budget.
    IdentifierByteLimit {
        /// Admitted maximum.
        limit: usize,
        /// Supplied count.
        actual: usize,
    },
    /// Input ended before a complete primitive could be read.
    Truncated {
        /// Byte offset at the failed read.
        at: usize,
        /// Bytes required by the primitive.
        needed: usize,
    },
    /// A decoded allocation could not be reserved safely.
    AllocationRefused {
        /// Aggregate being allocated.
        context: &'static str,
        /// Requested element count or string bytes.
        count: usize,
    },
    /// Wire contained a finite but noncanonical floating-point representation.
    NonCanonicalFloat {
        /// Byte offset of the encoded `f64`.
        at: usize,
        /// Rejected IEEE-754 bits.
        bits: u64,
    },
    /// Envelope magic did not identify this codec.
    InvalidMagic,
    /// Wire version is not implemented by this decoder.
    UnsupportedVersion {
        /// Version found on the wire.
        found: u16,
    },
    /// A closed enum discriminant was unknown.
    InvalidTag {
        /// Byte offset of the tag.
        at: usize,
        /// Enum being decoded.
        context: &'static str,
        /// Unknown byte.
        tag: u8,
    },
    /// Identifier bytes were not UTF-8.
    InvalidUtf8 {
        /// Byte offset of the string payload.
        at: usize,
    },
    /// Bytes followed the one complete payload.
    TrailingBytes {
        /// First trailing byte.
        at: usize,
    },
}

impl fmt::Display for PayloadError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidIdentifier { value, reason } => {
                write!(formatter, "invalid payload identifier {value:?}: {reason}")
            }
            Self::Empty { context } => write!(formatter, "{context} must not be empty"),
            Self::TooManyItems {
                context,
                limit,
                actual,
            } => write!(formatter, "{context} has {actual} items; limit is {limit}"),
            Self::CountOverflow { context } => write!(formatter, "{context} count overflowed"),
            Self::ShapeMismatch {
                context,
                expected,
                actual,
            } => write!(
                formatter,
                "{context} shape mismatch: expected {expected}, got {actual}"
            ),
            Self::InvalidTableTime { index, reason } => {
                write!(formatter, "invalid table time at {index}: {reason}")
            }
            Self::DistributionArity {
                family,
                expected,
                actual,
            } => write!(
                formatter,
                "{family:?} distribution expects {expected}; got {actual} values"
            ),
            Self::InvalidSemanticContract { context, semantic } => write!(
                formatter,
                "{context} rejects unsupported semantic kind/form {semantic:?}"
            ),
            Self::HeterogeneousContract { context } => {
                write!(
                    formatter,
                    "{context} requires a homogeneous quantity contract"
                )
            }
            Self::UnsupportedDistribution { context, family } => {
                write!(
                    formatter,
                    "{context} does not admit {family:?} distribution data"
                )
            }
            Self::DistributionOrder {
                context,
                index,
                reason,
            } => write!(
                formatter,
                "invalid {context} at parameter/support index {index}: {reason}"
            ),
            Self::NonFinite {
                context,
                index,
                value,
            } => write!(formatter, "{context}[{index}] is non-finite: {value}"),
            Self::DimensionMismatch {
                context,
                index,
                expected,
                actual,
            } => write!(
                formatter,
                "{context}[{index}] dimensions {actual:?} do not match {expected:?}"
            ),
            Self::Semantic(error) => write!(formatter, "semantic payload refusal: {error}"),
            Self::Chemistry(error) => write!(formatter, "chemistry payload refusal: {error}"),
            Self::PhasorSemanticContractRequired => formatter.write_str(
                "complex phasor payloads require a phasor-capable semantic Peak or RMS contract",
            ),
            Self::PhasorContractMismatch {
                index,
                expected,
                actual_kind,
                actual_amplitude,
            } => write!(
                formatter,
                "phasor[{index}] {actual_kind:?}/{actual_amplitude:?} does not match {expected:?}"
            ),
            Self::NonCanonicalSpeciesAxis { previous, next } => write!(
                formatter,
                "species axis is not strictly canonical: {previous:?} before {next:?}"
            ),
            Self::SpeciesAxisMismatch { sample } => {
                write!(formatter, "species axis changed at sample {sample}")
            }
            Self::NonCanonicalComponentAxis { previous, next } => write!(
                formatter,
                "characteristic axis is not strictly canonical: {previous:?} before {next:?}"
            ),
            Self::ByteLimit { limit, actual } => {
                write!(
                    formatter,
                    "payload has {actual} bytes; decode limit is {limit}"
                )
            }
            Self::ItemLimit { limit, actual } => write!(
                formatter,
                "payload requests {actual} aggregate items; decode limit is {limit}"
            ),
            Self::IdentifierByteLimit { limit, actual } => write!(
                formatter,
                "identifier has {actual} bytes; decode limit is {limit}"
            ),
            Self::Truncated { at, needed } => write!(
                formatter,
                "payload truncated at byte {at}; primitive needs {needed} bytes"
            ),
            Self::AllocationRefused { context, count } => {
                write!(
                    formatter,
                    "allocation refused for {context} with {count} items"
                )
            }
            Self::NonCanonicalFloat { at, bits } => write!(
                formatter,
                "noncanonical floating-point bits 0x{bits:016x} at byte {at}"
            ),
            Self::InvalidMagic => formatter.write_str("invalid payload envelope magic"),
            Self::UnsupportedVersion { found } => {
                write!(formatter, "unsupported payload wire version {found}")
            }
            Self::InvalidTag { at, context, tag } => {
                write!(formatter, "invalid {context} tag {tag} at byte {at}")
            }
            Self::InvalidUtf8 { at } => {
                write!(formatter, "invalid UTF-8 identifier at byte {at}")
            }
            Self::TrailingBytes { at } => {
                write!(formatter, "trailing payload bytes begin at {at}")
            }
        }
    }
}

impl std::error::Error for PayloadError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Semantic(error) => Some(error),
            Self::Chemistry(error) => Some(error),
            _ => None,
        }
    }
}

impl From<SemanticError> for PayloadError {
    fn from(error: SemanticError) -> Self {
        Self::Semantic(error)
    }
}

impl From<ChemistryError> for PayloadError {
    fn from(error: ChemistryError) -> Self {
        Self::Chemistry(error)
    }
}

fn check_item_count(context: &'static str, count: usize) -> Result<(), PayloadError> {
    if count > MAX_PAYLOAD_ITEMS {
        Err(PayloadError::TooManyItems {
            context,
            limit: MAX_PAYLOAD_ITEMS,
            actual: count,
        })
    } else {
        Ok(())
    }
}

fn validate_contract(
    contract: QuantityContract,
    context: &'static str,
) -> Result<(), PayloadError> {
    if let QuantityContract::Semantic(semantic) = contract
        && !semantic.kind().admits_scalar_form(semantic.form())
    {
        return Err(PayloadError::InvalidSemanticContract { context, semantic });
    }
    Ok(())
}

fn require_homogeneous(meta: &PayloadMeta, context: &'static str) -> Result<(), PayloadError> {
    validate_contract(meta.contract, context)?;
    if meta.contract == QuantityContract::Heterogeneous {
        Err(PayloadError::HeterogeneousContract { context })
    } else {
        Ok(())
    }
}

fn require_heterogeneous(meta: &PayloadMeta, context: &'static str) -> Result<(), PayloadError> {
    if meta.contract == QuantityContract::Heterogeneous {
        Ok(())
    } else {
        Err(PayloadError::ShapeMismatch {
            context,
            expected: 0,
            actual: 1,
        })
    }
}

fn validate_table_axis(
    times: &[QtyAny],
    value_count: usize,
    outside: OutsideDomainPolicy,
) -> Result<(), PayloadError> {
    if times.is_empty() {
        return Err(PayloadError::Empty {
            context: "table samples",
        });
    }
    check_item_count("table samples", times.len())?;
    if times.len() != value_count {
        return Err(PayloadError::ShapeMismatch {
            context: "table times and values",
            expected: times.len(),
            actual: value_count,
        });
    }
    if outside == OutsideDomainPolicy::Periodic && times.len() < 2 {
        return Err(PayloadError::InvalidTableTime {
            index: 0,
            reason: "periodic tables require at least two samples",
        });
    }
    let mut previous = None;
    for (index, time) in times.iter().enumerate() {
        if time.dims != TIME_DIMS {
            return Err(PayloadError::InvalidTableTime {
                index,
                reason: "time must use the six-base seconds dimension",
            });
        }
        if !time.value.is_finite() {
            return Err(PayloadError::InvalidTableTime {
                index,
                reason: "time must be finite",
            });
        }
        if previous.is_some_and(|prior| time.value <= prior) {
            return Err(PayloadError::InvalidTableTime {
                index,
                reason: "times must be strictly increasing",
            });
        }
        previous = Some(time.value);
    }
    Ok(())
}

fn validate_distribution_arity(
    family: DistributionFamily,
    count: usize,
) -> Result<(), PayloadError> {
    check_item_count("distribution values", count)?;
    let valid = match family {
        DistributionFamily::Normal | DistributionFamily::Uniform => count == 2,
        DistributionFamily::Empirical => count >= 1,
    };
    if valid {
        Ok(())
    } else {
        Err(PayloadError::DistributionArity {
            family,
            expected: match family {
                DistributionFamily::Normal => "mean and standard deviation",
                DistributionFamily::Uniform => "lower and upper endpoints",
                DistributionFamily::Empirical => "one or more support samples",
            },
            actual: count,
        })
    }
}

fn validate_source<T>(
    source: &SampleSource<T>,
    mut validate_value: impl FnMut(&T, usize) -> Result<(), PayloadError>,
) -> Result<(), PayloadError> {
    match source {
        SampleSource::Fixed(value) => validate_value(value, 0),
        SampleSource::Table {
            times,
            values,
            outside,
            ..
        } => {
            validate_table_axis(times, values.len(), *outside)?;
            for (index, value) in values.iter().enumerate() {
                validate_value(value, index)?;
            }
            Ok(())
        }
        SampleSource::Distribution { family, parameters } => {
            validate_distribution_arity(*family, parameters.len())?;
            for (index, value) in parameters.iter().enumerate() {
                validate_value(value, index)?;
            }
            Ok(())
        }
    }
}

fn first_source_value<T>(source: &SampleSource<T>) -> Option<&T> {
    match source {
        SampleSource::Fixed(value) => Some(value),
        SampleSource::Table { values, .. } => values.first(),
        SampleSource::Distribution { parameters, .. } => parameters.first(),
    }
}

fn validate_quantity(
    quantity: QtyAny,
    contract: QuantityContract,
    context: &'static str,
    index: usize,
) -> Result<(), PayloadError> {
    if !quantity.value.is_finite() {
        return Err(PayloadError::NonFinite {
            context,
            index,
            value: quantity.value,
        });
    }
    let expected = contract
        .dims()
        .ok_or(PayloadError::HeterogeneousContract { context })?;
    if quantity.dims != expected {
        return Err(PayloadError::DimensionMismatch {
            context,
            index,
            expected,
            actual: quantity.dims,
        });
    }
    if let QuantityContract::Semantic(semantic) = contract {
        SemanticQty::new(quantity, semantic)?;
    }
    Ok(())
}

const fn flattened_index(sample: usize, width: usize, component: usize) -> usize {
    sample.saturating_mul(width).saturating_add(component)
}

fn validate_phasor(
    phasor: PhasorQty,
    contract: QuantityContract,
    index: usize,
) -> Result<(), PayloadError> {
    let real = phasor.real();
    let imaginary = phasor.imaginary();
    for (component_index, quantity) in [real, imaginary].into_iter().enumerate() {
        if !quantity.value.is_finite() {
            return Err(PayloadError::NonFinite {
                context: "phasor component",
                index: flattened_index(index, 2, component_index),
                value: quantity.value,
            });
        }
        let expected_dims = contract.dims().ok_or(PayloadError::HeterogeneousContract {
            context: "phasor payload",
        })?;
        if quantity.dims != expected_dims {
            return Err(PayloadError::DimensionMismatch {
                context: "phasor component",
                index: flattened_index(index, 2, component_index),
                expected: expected_dims,
                actual: quantity.dims,
            });
        }
    }
    if let QuantityContract::Semantic(semantic) = contract {
        let expected_amplitude = match semantic.form() {
            ValueForm::Peak => Some(PhasorAmplitude::Peak),
            ValueForm::Rms => Some(PhasorAmplitude::Rms),
            ValueForm::Static | ValueForm::Instantaneous => None,
        };
        if semantic.kind() != phasor.kind() || expected_amplitude != Some(phasor.amplitude()) {
            return Err(PayloadError::PhasorContractMismatch {
                index,
                expected: Some(semantic),
                actual_kind: phasor.kind(),
                actual_amplitude: phasor.amplitude(),
            });
        }
    }
    Ok(())
}

fn aggregate_meta_items(meta: &PayloadMeta) -> Result<usize, PayloadError> {
    let reset = usize::from(matches!(
        &meta.reference,
        ReferenceSemantics::ResetAtEvent(_)
    ));
    1_usize
        .checked_add(reset)
        .ok_or(PayloadError::CountOverflow {
            context: "payload metadata items",
        })
}

fn aggregate_source_items<T>(
    source: &SampleSource<T>,
    nested_items: impl Fn(&T) -> usize,
    base: usize,
) -> Option<usize> {
    match source {
        SampleSource::Fixed(value) => base.checked_add(nested_items(value)),
        SampleSource::Table { times, values, .. } => {
            let outer = times.len().checked_add(values.len())?;
            values
                .iter()
                .try_fold(base.checked_add(outer)?, |count, value| {
                    count.checked_add(nested_items(value))
                })
        }
        SampleSource::Distribution { parameters, .. } => parameters
            .iter()
            .try_fold(base.checked_add(parameters.len())?, |count, value| {
                count.checked_add(nested_items(value))
            }),
    }
}

fn enforce_aggregate_limit(count: usize, limit: usize) -> Result<(), PayloadError> {
    if count > limit {
        Err(PayloadError::TooManyItems {
            context: "payload aggregate items",
            limit,
            actual: count,
        })
    } else {
        Ok(())
    }
}

fn validate_scalar_distribution(
    source: &SampleSource<QtyAny>,
    context: &'static str,
) -> Result<(), PayloadError> {
    let SampleSource::Distribution { family, parameters } = source else {
        return Ok(());
    };
    match family {
        DistributionFamily::Normal if parameters[1].value < 0.0 => {
            Err(PayloadError::DistributionOrder {
                context,
                index: 1,
                reason: "normal scale must be nonnegative",
            })
        }
        DistributionFamily::Uniform if parameters[0].value > parameters[1].value => {
            Err(PayloadError::DistributionOrder {
                context,
                index: 1,
                reason: "uniform lower endpoint exceeds upper endpoint",
            })
        }
        DistributionFamily::Empirical => ensure_strict_order(
            parameters,
            |left, right| canonical_value_cmp(left.value, right.value),
            context,
        ),
        DistributionFamily::Normal | DistributionFamily::Uniform => Ok(()),
    }
}

fn validate_vector_distribution(
    source: &SampleSource<Vec<QtyAny>>,
    context: &'static str,
) -> Result<(), PayloadError> {
    let SampleSource::Distribution { family, parameters } = source else {
        return Ok(());
    };
    match family {
        DistributionFamily::Normal => {
            if parameters[1].iter().any(|value| value.value < 0.0) {
                Err(PayloadError::DistributionOrder {
                    context,
                    index: 1,
                    reason: "every normal scale component must be nonnegative",
                })
            } else {
                Ok(())
            }
        }
        DistributionFamily::Uniform => {
            if parameters[0]
                .iter()
                .zip(&parameters[1])
                .any(|(lower, upper)| lower.value > upper.value)
            {
                Err(PayloadError::DistributionOrder {
                    context,
                    index: 1,
                    reason: "a uniform lower component exceeds its upper component",
                })
            } else {
                Ok(())
            }
        }
        DistributionFamily::Empirical => ensure_strict_order(
            parameters,
            |left, right| qty_slice_cmp(left, right),
            context,
        ),
    }
}

fn validate_phasor_distribution(source: &SampleSource<PhasorQty>) -> Result<(), PayloadError> {
    let SampleSource::Distribution { family, parameters } = source else {
        return Ok(());
    };
    match family {
        DistributionFamily::Normal | DistributionFamily::Uniform => {
            Err(PayloadError::UnsupportedDistribution {
                context: "complex phasor payload",
                family: *family,
            })
        }
        DistributionFamily::Empirical => ensure_strict_order(
            parameters,
            |left, right| {
                canonical_value_cmp(left.real().value, right.real().value).then_with(|| {
                    canonical_value_cmp(left.imaginary().value, right.imaginary().value)
                })
            },
            "phasor empirical support",
        ),
    }
}

fn validate_species_distribution(
    meta: &PayloadMeta,
    source: &SampleSource<Vec<SpeciesValue>>,
) -> Result<(), PayloadError> {
    let SampleSource::Distribution { family, parameters } = source else {
        return Ok(());
    };
    if matches!(
        meta.contract,
        QuantityContract::Semantic(semantic)
            if matches!(semantic.kind(), QuantityKind::Composition(_))
    ) && !matches!(family, DistributionFamily::Empirical)
    {
        return Err(PayloadError::UnsupportedDistribution {
            context: "semantic composition species bundle",
            family: *family,
        });
    }
    match family {
        DistributionFamily::Normal => {
            if parameters[1].iter().any(|entry| entry.quantity.value < 0.0) {
                Err(PayloadError::DistributionOrder {
                    context: "species normal distribution",
                    index: 1,
                    reason: "every normal scale component must be nonnegative",
                })
            } else {
                Ok(())
            }
        }
        DistributionFamily::Uniform => {
            if parameters[0]
                .iter()
                .zip(&parameters[1])
                .any(|(lower, upper)| lower.quantity.value > upper.quantity.value)
            {
                Err(PayloadError::DistributionOrder {
                    context: "species uniform distribution",
                    index: 1,
                    reason: "a uniform lower component exceeds its upper component",
                })
            } else {
                Ok(())
            }
        }
        DistributionFamily::Empirical => ensure_strict_order(
            parameters,
            |left, right| {
                left.iter()
                    .zip(right)
                    .map(|(left, right)| {
                        canonical_value_cmp(left.quantity.value, right.quantity.value)
                    })
                    .find(|order| !order.is_eq())
                    .unwrap_or_else(|| left.len().cmp(&right.len()))
            },
            "species empirical support",
        ),
    }
}

fn validate_composition_samples(
    meta: &PayloadMeta,
    source: &SampleSource<Vec<SpeciesValue>>,
) -> Result<(), PayloadError> {
    let QuantityContract::Semantic(semantic) = meta.contract else {
        return Ok(());
    };
    let QuantityKind::Composition(basis) = semantic.kind() else {
        return Ok(());
    };
    validate_source(source, |sample, _| {
        // Type pin (batch triage): with the receiver uninferred here, both
        // `Extend::extend` and `PayloadByteSink::extend for Vec<u8>` are
        // method candidates (E0034); Composition::new fixes f64 anyway.
        let mut fractions: Vec<f64> = Vec::new();
        fractions
            .try_reserve_exact(sample.len())
            .map_err(|_| PayloadError::AllocationRefused {
                context: "composition validation scratch",
                count: sample.len(),
            })?;
        fractions.extend(sample.iter().map(|entry| entry.quantity.value));
        Composition::new(basis, fractions)?;
        Ok(())
    })
}

fn try_copy_text(text: &str, context: &'static str) -> Result<String, PayloadError> {
    let mut copy = String::new();
    copy.try_reserve_exact(text.len())
        .map_err(|_| PayloadError::AllocationRefused {
            context,
            count: text.len(),
        })?;
    copy.push_str(text);
    Ok(copy)
}

fn qty_slice_cmp(left: &[QtyAny], right: &[QtyAny]) -> core::cmp::Ordering {
    left.iter()
        .zip(right)
        .map(|(left, right)| canonical_value_cmp(left.value, right.value))
        .find(|order| !order.is_eq())
        .unwrap_or_else(|| left.len().cmp(&right.len()))
}

fn canonical_value_cmp(left: f64, right: f64) -> core::cmp::Ordering {
    let left = if left == 0.0 { 0.0 } else { left };
    let right = if right == 0.0 { 0.0 } else { right };
    left.total_cmp(&right)
}

fn ensure_strict_order<T>(
    values: &[T],
    mut compare: impl FnMut(&T, &T) -> core::cmp::Ordering,
    context: &'static str,
) -> Result<(), PayloadError> {
    for (index, pair) in values.windows(2).enumerate() {
        if !compare(&pair[0], &pair[1]).is_lt() {
            return Err(PayloadError::DistributionOrder {
                context,
                index: index + 1,
                reason: "empirical support must be strictly canonical and duplicate-free",
            });
        }
    }
    Ok(())
}

trait PayloadByteSink {
    fn extend(&mut self, bytes: &[u8]) -> Result<(), PayloadError>;

    fn push(&mut self, byte: u8) -> Result<(), PayloadError> {
        self.extend(core::slice::from_ref(&byte))
    }
}

impl PayloadByteSink for Vec<u8> {
    fn extend(&mut self, bytes: &[u8]) -> Result<(), PayloadError> {
        self.extend_from_slice(bytes);
        Ok(())
    }
}

#[derive(Default)]
struct PayloadByteCounter {
    len: usize,
}

impl PayloadByteSink for PayloadByteCounter {
    fn extend(&mut self, bytes: &[u8]) -> Result<(), PayloadError> {
        self.len = self
            .len
            .checked_add(bytes.len())
            .ok_or(PayloadError::CountOverflow {
                context: "canonical payload byte length",
            })?;
        Ok(())
    }
}

fn encode_payload(out: &mut impl PayloadByteSink, payload: &Payload) -> Result<(), PayloadError> {
    out.extend(PAYLOAD_MAGIC)?;
    put_u16(out, PAYLOAD_WIRE_VERSION)?;
    match payload {
        Payload::Scalar(value) => {
            put_u8(out, 0)?;
            write_meta(out, &value.meta)?;
            write_source(out, &value.source, write_qty)?;
        }
        Payload::Vector(value) => {
            put_u8(out, 1)?;
            write_meta(out, &value.meta)?;
            write_source(out, &value.source, |out, values| {
                write_qty_slice(out, values)
            })?;
        }
        Payload::Tensor(value) => {
            put_u8(out, 2)?;
            write_meta(out, &value.meta)?;
            put_len(out, value.rows)?;
            put_len(out, value.columns)?;
            write_source(out, &value.source, |out, values| {
                write_qty_slice(out, values)
            })?;
        }
        Payload::ComplexPhasor(value) => {
            put_u8(out, 3)?;
            write_meta(out, &value.meta)?;
            write_source(out, &value.source, write_phasor)?;
        }
        Payload::SpeciesBundle(value) => {
            put_u8(out, 4)?;
            write_meta(out, &value.meta)?;
            write_source(out, &value.source, |out, values| {
                write_species_sample(out, values)
            })?;
        }
        Payload::CharacteristicState(value) => {
            put_u8(out, 5)?;
            write_meta(out, &value.meta)?;
            put_len(out, value.components.len())?;
            for component in &value.components {
                write_id(out, &component.name)?;
                put_u8(
                    out,
                    match component.direction {
                        CharacteristicDirection::Incoming => 0,
                        CharacteristicDirection::Outgoing => 1,
                        CharacteristicDirection::Stationary => 2,
                    },
                )?;
                write_contract(out, component.contract)?;
            }
            write_source(out, &value.source, |out, values| {
                write_qty_slice(out, values)
            })?;
        }
        Payload::FieldTraceRef(value) => {
            put_u8(out, 6)?;
            write_meta(out, &value.meta)?;
            write_id(out, &value.artifact)?;
            write_id(out, &value.field)?;
        }
        Payload::PortRef(value) => {
            put_u8(out, 7)?;
            write_meta(out, &value.meta)?;
            write_id(out, &value.component)?;
            write_id(out, &value.port)?;
        }
    }
    Ok(())
}

/// Exact canonical V1 byte length without materializing the envelope.
///
/// This traverses the same encoder grammar as [`try_canonical_payload_bytes`],
/// so scenario-IR emission can admit its complete output allocation before
/// requesting heap storage.
pub fn canonical_payload_byte_len(payload: &Payload) -> Result<usize, PayloadError> {
    let mut counter = PayloadByteCounter::default();
    encode_payload(&mut counter, payload)?;
    if counter.len > MAX_PAYLOAD_WIRE_BYTES {
        return Err(PayloadError::ByteLimit {
            limit: MAX_PAYLOAD_WIRE_BYTES,
            actual: counter.len,
        });
    }
    Ok(counter.len)
}

/// Encode one admitted payload after fallibly reserving its exact V1 length.
///
/// Floating-point zero is normalized to positive zero; all other finite values
/// retain their exact IEEE-754 bit pattern. No partially encoded buffer is
/// published when the exact reservation is refused.
pub fn try_canonical_payload_bytes(payload: &Payload) -> Result<Vec<u8>, PayloadError> {
    let len = canonical_payload_byte_len(payload)?;
    let mut out = Vec::new();
    out.try_reserve_exact(len)
        .map_err(|_| PayloadError::AllocationRefused {
            context: "canonical payload bytes",
            count: len,
        })?;
    encode_payload(&mut out, payload)?;
    debug_assert_eq!(out.len(), len);
    Ok(out)
}

/// Encode one admitted payload as canonical V1 bytes.
///
/// This compatibility convenience delegates to the fallible exact-reservation
/// path. Resource-authoritative callers use [`try_canonical_payload_bytes`].
#[must_use]
pub fn canonical_payload_bytes(payload: &Payload) -> Vec<u8> {
    try_canonical_payload_bytes(payload)
        .expect("validated payload must fit the closed canonical V1 envelope")
}

/// Decode one canonical payload under [`PayloadDecodeLimits::DEFAULT`].
pub fn decode_payload(bytes: &[u8]) -> Result<Payload, PayloadError> {
    decode_payload_with_limits(bytes, PayloadDecodeLimits::DEFAULT)
}

/// Decode one canonical payload with caller-supplied allocation limits.
///
/// The byte limit is checked before parsing.  Every vector/string length is
/// charged to one aggregate item budget before allocation, and no trailing
/// extension bytes are accepted in V1.
pub fn decode_payload_with_limits(
    bytes: &[u8],
    limits: PayloadDecodeLimits,
) -> Result<Payload, PayloadError> {
    if bytes.len() > limits.max_bytes {
        return Err(PayloadError::ByteLimit {
            limit: limits.max_bytes,
            actual: bytes.len(),
        });
    }
    let mut decoder = Decoder::new(bytes, limits);
    if decoder.take(PAYLOAD_MAGIC.len())? != PAYLOAD_MAGIC {
        return Err(PayloadError::InvalidMagic);
    }
    let version = decoder.u16()?;
    if version != PAYLOAD_WIRE_VERSION {
        return Err(PayloadError::UnsupportedVersion { found: version });
    }
    let variant_at = decoder.at;
    let variant = decoder.u8()?;
    if variant > 7 {
        return Err(PayloadError::InvalidTag {
            at: variant_at,
            context: "payload variant",
            tag: variant,
        });
    }
    let meta = read_meta(&mut decoder)?;
    let payload = match variant {
        0 => Payload::Scalar(ScalarPayload::new(
            meta,
            read_source(&mut decoder, 14, |decoder| decoder.qty())?,
        )?),
        1 => Payload::Vector(VectorPayload::new(
            meta,
            read_source(&mut decoder, 4, read_qty_vec)?,
        )?),
        2 => {
            let rows = decoder.raw_len()?;
            let columns = decoder.raw_len()?;
            Payload::Tensor(TensorPayload::new(
                meta,
                rows,
                columns,
                read_source(&mut decoder, 4, read_qty_vec)?,
            )?)
        }
        3 => Payload::ComplexPhasor(ComplexPhasorPayload::new(
            meta,
            read_source(&mut decoder, 30, read_phasor)?,
        )?),
        4 => {
            let source = read_source(&mut decoder, 4, read_species_sample)?;
            let axis_items = first_source_value(&source).map_or(0, Vec::len);
            let axis_cache_items =
                axis_items
                    .checked_mul(2)
                    .ok_or(PayloadError::CountOverflow {
                        context: "species-axis cache items",
                    })?;
            decoder.charge_items(axis_cache_items)?;
            Payload::SpeciesBundle(SpeciesBundle::new_with_item_limit(
                meta,
                source,
                decoder.limits.max_items,
            )?)
        }
        5 => {
            let count = decoder.len("characteristic components")?;
            decoder.ensure_can_charge(count)?;
            let mut components = decoder.allocate_vec(count, "characteristic components", 7)?;
            for _ in 0..count {
                let name = decoder.id()?;
                let at = decoder.at;
                let direction = match decoder.u8()? {
                    0 => CharacteristicDirection::Incoming,
                    1 => CharacteristicDirection::Outgoing,
                    2 => CharacteristicDirection::Stationary,
                    tag => {
                        return Err(PayloadError::InvalidTag {
                            at,
                            context: "characteristic direction",
                            tag,
                        });
                    }
                };
                let contract = read_contract(&mut decoder)?;
                components.push(CharacteristicComponent::new(name, direction, contract)?);
            }
            Payload::CharacteristicState(CharacteristicState::new(
                meta,
                components,
                read_source(&mut decoder, 4, read_qty_vec)?,
            )?)
        }
        6 => Payload::FieldTraceRef(FieldTraceRef::new(meta, decoder.id()?, decoder.id()?)?),
        7 => Payload::PortRef(PortRef::new(meta, decoder.id()?, decoder.id()?)?),
        _ => {
            return Err(PayloadError::InvalidTag {
                at: variant_at,
                context: "payload variant",
                tag: variant,
            });
        }
    };
    if decoder.at != bytes.len() {
        return Err(PayloadError::TrailingBytes { at: decoder.at });
    }
    Ok(payload)
}

fn put_u8(out: &mut impl PayloadByteSink, value: u8) -> Result<(), PayloadError> {
    out.push(value)
}

fn put_u16(out: &mut impl PayloadByteSink, value: u16) -> Result<(), PayloadError> {
    out.extend(&value.to_le_bytes())
}

fn put_u32(out: &mut impl PayloadByteSink, value: u32) -> Result<(), PayloadError> {
    out.extend(&value.to_le_bytes())
}

fn put_u64(out: &mut impl PayloadByteSink, value: u64) -> Result<(), PayloadError> {
    out.extend(&value.to_le_bytes())
}

fn put_len(out: &mut impl PayloadByteSink, value: usize) -> Result<(), PayloadError> {
    let value = u32::try_from(value).map_err(|_| PayloadError::CountOverflow {
        context: "canonical payload collection length",
    })?;
    put_u32(out, value)
}

fn canonical_bits(value: f64) -> u64 {
    if value == 0.0 { 0 } else { value.to_bits() }
}

fn write_qty(out: &mut impl PayloadByteSink, quantity: &QtyAny) -> Result<(), PayloadError> {
    put_u64(out, canonical_bits(quantity.value))?;
    for exponent in quantity.dims.0 {
        put_u8(out, exponent as u8)?;
    }
    Ok(())
}

fn write_dims(out: &mut impl PayloadByteSink, dims: Dims) -> Result<(), PayloadError> {
    for exponent in dims.0 {
        put_u8(out, exponent as u8)?;
    }
    Ok(())
}

fn write_id(out: &mut impl PayloadByteSink, id: &PayloadId) -> Result<(), PayloadError> {
    put_len(out, id.0.len())?;
    out.extend(id.0.as_bytes())
}

fn write_species_id(out: &mut impl PayloadByteSink, id: &SpeciesId) -> Result<(), PayloadError> {
    put_len(out, id.as_str().len())?;
    out.extend(id.as_str().as_bytes())
}

fn write_meta(out: &mut impl PayloadByteSink, meta: &PayloadMeta) -> Result<(), PayloadError> {
    write_contract(out, meta.contract)?;
    write_id(out, &meta.basis_key)?;
    put_u32(out, meta.frame.0)?;
    put_u8(
        out,
        match meta.orientation {
            OrientationParity::Even => 0,
            OrientationParity::Odd => 1,
        },
    )?;
    match &meta.reference {
        ReferenceSemantics::Continuous => put_u8(out, 0)?,
        ReferenceSemantics::ResetAtEvent(event) => {
            put_u8(out, 1)?;
            write_id(out, event)?;
        }
    }
    Ok(())
}

fn write_contract(
    out: &mut impl PayloadByteSink,
    contract: QuantityContract,
) -> Result<(), PayloadError> {
    match contract {
        QuantityContract::Dimensions(dims) => {
            put_u8(out, 0)?;
            write_dims(out, dims)?;
        }
        QuantityContract::Semantic(semantic) => {
            put_u8(out, 1)?;
            write_semantic_type(out, semantic)?;
        }
        QuantityContract::Heterogeneous => put_u8(out, 2)?,
    }
    Ok(())
}

fn write_source<T, S: PayloadByteSink>(
    out: &mut S,
    source: &SampleSource<T>,
    mut write_value: impl FnMut(&mut S, &T) -> Result<(), PayloadError>,
) -> Result<(), PayloadError> {
    match source {
        SampleSource::Fixed(value) => {
            put_u8(out, 0)?;
            write_value(out, value)?;
        }
        SampleSource::Table {
            times,
            values,
            interpolation,
            outside,
        } => {
            put_u8(out, 1)?;
            put_u8(
                out,
                match interpolation {
                    TableInterpolation::StepLeft => 0,
                    TableInterpolation::Linear => 1,
                },
            )?;
            put_u8(
                out,
                match outside {
                    OutsideDomainPolicy::Refuse => 0,
                    OutsideDomainPolicy::Clamp => 1,
                    OutsideDomainPolicy::Periodic => 2,
                },
            )?;
            put_len(out, times.len())?;
            for time in times {
                write_qty(out, time)?;
            }
            for value in values {
                write_value(out, value)?;
            }
        }
        SampleSource::Distribution { family, parameters } => {
            put_u8(out, 2)?;
            put_u8(
                out,
                match family {
                    DistributionFamily::Normal => 0,
                    DistributionFamily::Uniform => 1,
                    DistributionFamily::Empirical => 2,
                },
            )?;
            put_len(out, parameters.len())?;
            for parameter in parameters {
                write_value(out, parameter)?;
            }
        }
    }
    Ok(())
}

fn write_qty_slice(out: &mut impl PayloadByteSink, values: &[QtyAny]) -> Result<(), PayloadError> {
    put_len(out, values.len())?;
    for value in values {
        write_qty(out, value)?;
    }
    Ok(())
}

fn write_phasor(out: &mut impl PayloadByteSink, phasor: &PhasorQty) -> Result<(), PayloadError> {
    write_qty(out, &phasor.real())?;
    write_qty(out, &phasor.imaginary())?;
    write_quantity_kind(out, phasor.kind())?;
    put_u8(
        out,
        match phasor.amplitude() {
            PhasorAmplitude::Peak => 0,
            PhasorAmplitude::Rms => 1,
        },
    )
}

fn write_species_sample(
    out: &mut impl PayloadByteSink,
    values: &[SpeciesValue],
) -> Result<(), PayloadError> {
    put_len(out, values.len())?;
    for value in values {
        write_species_id(out, &value.species)?;
        write_qty(out, &value.quantity)?;
    }
    Ok(())
}

fn write_semantic_type(
    out: &mut impl PayloadByteSink,
    semantic: SemanticType,
) -> Result<(), PayloadError> {
    write_quantity_kind(out, semantic.kind())?;
    put_u8(
        out,
        match semantic.form() {
            ValueForm::Static => 0,
            ValueForm::Instantaneous => 1,
            ValueForm::Peak => 2,
            ValueForm::Rms => 3,
        },
    )
}

fn write_quantity_kind(
    out: &mut impl PayloadByteSink,
    kind: QuantityKind,
) -> Result<(), PayloadError> {
    match kind {
        QuantityKind::AbsoluteTemperature => put_u8(out, 0)?,
        QuantityKind::TemperatureDifference => put_u8(out, 1)?,
        QuantityKind::Angle(domain) => {
            put_u8(out, 2)?;
            write_angle_domain(out, domain)?;
        }
        QuantityKind::AngularVelocity(domain) => {
            put_u8(out, 3)?;
            write_angle_domain(out, domain)?;
        }
        QuantityKind::Torque => put_u8(out, 4)?,
        QuantityKind::Energy => put_u8(out, 5)?,
        QuantityKind::Pressure => put_u8(out, 6)?,
        QuantityKind::Stress => put_u8(out, 7)?,
        QuantityKind::Strain { basis, component } => {
            put_u8(out, 8)?;
            put_u8(
                out,
                match basis {
                    StrainBasis::Tensor => 0,
                    StrainBasis::Engineering => 1,
                },
            )?;
            put_u8(
                out,
                match component {
                    StrainComponent::Normal => 0,
                    StrainComponent::Shear => 1,
                },
            )?;
        }
        QuantityKind::Composition(basis) => {
            put_u8(out, 9)?;
            put_u8(
                out,
                match basis {
                    CompositionBasis::MassFraction => 0,
                    CompositionBasis::MoleFraction => 1,
                    CompositionBasis::VolumeFraction => 2,
                },
            )?;
        }
        QuantityKind::Mass => put_u8(out, 10)?,
        QuantityKind::Amount => put_u8(out, 11)?,
        QuantityKind::MolarMass => put_u8(out, 12)?,
        QuantityKind::MassConcentration => put_u8(out, 13)?,
        QuantityKind::AmountConcentration => put_u8(out, 14)?,
        QuantityKind::Entropy => put_u8(out, 15)?,
        QuantityKind::HeatCapacity => put_u8(out, 16)?,
        QuantityKind::AcousticPressure => put_u8(out, 17)?,
        QuantityKind::AcousticPower => put_u8(out, 18)?,
    }
    Ok(())
}

fn write_angle_domain(
    out: &mut impl PayloadByteSink,
    domain: AngleDomain,
) -> Result<(), PayloadError> {
    put_u8(
        out,
        match domain {
            AngleDomain::Mechanical => 0,
            AngleDomain::Electrical => 1,
        },
    )
}

struct Decoder<'a> {
    bytes: &'a [u8],
    at: usize,
    limits: PayloadDecodeLimits,
    items: usize,
}

impl<'a> Decoder<'a> {
    fn new(bytes: &'a [u8], limits: PayloadDecodeLimits) -> Self {
        let limits = PayloadDecodeLimits {
            max_bytes: limits.max_bytes,
            max_items: limits.max_items.min(MAX_PAYLOAD_ITEMS),
            max_identifier_bytes: limits.max_identifier_bytes.min(MAX_PAYLOAD_ID_BYTES),
        };
        Self {
            bytes,
            at: 0,
            limits,
            items: 0,
        }
    }

    fn take(&mut self, count: usize) -> Result<&'a [u8], PayloadError> {
        let end = self.at.checked_add(count).ok_or(PayloadError::Truncated {
            at: self.at,
            needed: count,
        })?;
        let value = self
            .bytes
            .get(self.at..end)
            .ok_or(PayloadError::Truncated {
                at: self.at,
                needed: count,
            })?;
        self.at = end;
        Ok(value)
    }

    fn u8(&mut self) -> Result<u8, PayloadError> {
        Ok(self.take(1)?[0])
    }

    fn u16(&mut self) -> Result<u16, PayloadError> {
        let mut bytes = [0; 2];
        bytes.copy_from_slice(self.take(2)?);
        Ok(u16::from_le_bytes(bytes))
    }

    fn u32(&mut self) -> Result<u32, PayloadError> {
        let mut bytes = [0; 4];
        bytes.copy_from_slice(self.take(4)?);
        Ok(u32::from_le_bytes(bytes))
    }

    fn u64(&mut self) -> Result<u64, PayloadError> {
        let mut bytes = [0; 8];
        bytes.copy_from_slice(self.take(8)?);
        Ok(u64::from_le_bytes(bytes))
    }

    fn len(&mut self, _context: &'static str) -> Result<usize, PayloadError> {
        let count = self.raw_len()?;
        self.charge_items(count)?;
        Ok(count)
    }

    fn raw_len(&mut self) -> Result<usize, PayloadError> {
        Ok(self.u32()? as usize)
    }

    fn charge_items(&mut self, count: usize) -> Result<(), PayloadError> {
        let actual = self.ensure_can_charge(count)?;
        self.items = actual;
        Ok(())
    }

    fn ensure_can_charge(&self, count: usize) -> Result<usize, PayloadError> {
        let actual = self
            .items
            .checked_add(count)
            .ok_or(PayloadError::ItemLimit {
                limit: self.limits.max_items,
                actual: usize::MAX,
            })?;
        if actual > self.limits.max_items {
            return Err(PayloadError::ItemLimit {
                limit: self.limits.max_items,
                actual,
            });
        }
        Ok(actual)
    }

    fn dims(&mut self) -> Result<Dims, PayloadError> {
        let mut exponents = [0_i8; 6];
        for exponent in &mut exponents {
            *exponent = self.u8()? as i8;
        }
        Ok(Dims(exponents))
    }

    fn qty(&mut self) -> Result<QtyAny, PayloadError> {
        let at = self.at;
        let bits = self.u64()?;
        if bits == (-0.0_f64).to_bits() {
            return Err(PayloadError::NonCanonicalFloat { at, bits });
        }
        let value = f64::from_bits(bits);
        let dims = self.dims()?;
        Ok(QtyAny::new(value, dims))
    }

    fn remaining(&self) -> usize {
        self.bytes.len().saturating_sub(self.at)
    }

    fn allocate_vec<T>(
        &self,
        count: usize,
        context: &'static str,
        minimum_wire_bytes: usize,
    ) -> Result<Vec<T>, PayloadError> {
        let needed = count
            .checked_mul(minimum_wire_bytes)
            .ok_or(PayloadError::AllocationRefused { context, count })?;
        if self.remaining() < needed {
            return Err(PayloadError::Truncated {
                at: self.at,
                needed,
            });
        }
        let mut values = Vec::new();
        values
            .try_reserve_exact(count)
            .map_err(|_| PayloadError::AllocationRefused { context, count })?;
        Ok(values)
    }

    fn string(&mut self) -> Result<String, PayloadError> {
        let length = self.u32()? as usize;
        if length > self.limits.max_identifier_bytes {
            return Err(PayloadError::IdentifierByteLimit {
                limit: self.limits.max_identifier_bytes,
                actual: length,
            });
        }
        self.charge_items(1)?;
        let at = self.at;
        let bytes = self.take(length)?;
        let text = core::str::from_utf8(bytes).map_err(|_| PayloadError::InvalidUtf8 { at })?;
        let mut owned = String::new();
        owned
            .try_reserve_exact(length)
            .map_err(|_| PayloadError::AllocationRefused {
                context: "payload identifier",
                count: length,
            })?;
        owned.push_str(text);
        Ok(owned)
    }

    fn id(&mut self) -> Result<PayloadId, PayloadError> {
        PayloadId::new(self.string()?)
    }

    fn species_id(&mut self) -> Result<SpeciesId, PayloadError> {
        Ok(SpeciesId::new(self.string()?)?)
    }
}

fn read_meta(decoder: &mut Decoder<'_>) -> Result<PayloadMeta, PayloadError> {
    let contract = read_contract(decoder)?;
    let basis_key = decoder.id()?;
    let frame = FrameId(decoder.u32()?);
    let orientation_at = decoder.at;
    let orientation = match decoder.u8()? {
        0 => OrientationParity::Even,
        1 => OrientationParity::Odd,
        tag => {
            return Err(PayloadError::InvalidTag {
                at: orientation_at,
                context: "orientation parity",
                tag,
            });
        }
    };
    let reference_at = decoder.at;
    let reference = match decoder.u8()? {
        0 => ReferenceSemantics::Continuous,
        1 => ReferenceSemantics::ResetAtEvent(decoder.id()?),
        tag => {
            return Err(PayloadError::InvalidTag {
                at: reference_at,
                context: "reference semantics",
                tag,
            });
        }
    };
    PayloadMeta::new(contract, basis_key, frame, orientation, reference)
}

fn read_contract(decoder: &mut Decoder<'_>) -> Result<QuantityContract, PayloadError> {
    let contract_at = decoder.at;
    let contract = match decoder.u8()? {
        0 => QuantityContract::Dimensions(decoder.dims()?),
        1 => QuantityContract::Semantic(read_semantic_type(decoder)?),
        2 => QuantityContract::Heterogeneous,
        tag => {
            return Err(PayloadError::InvalidTag {
                at: contract_at,
                context: "quantity contract",
                tag,
            });
        }
    };
    Ok(contract)
}

fn read_source<T>(
    decoder: &mut Decoder<'_>,
    minimum_value_bytes: usize,
    mut read_value: impl FnMut(&mut Decoder<'_>) -> Result<T, PayloadError>,
) -> Result<SampleSource<T>, PayloadError> {
    let source_at = decoder.at;
    match decoder.u8()? {
        0 => Ok(SampleSource::Fixed(read_value(decoder)?)),
        1 => {
            let interpolation_at = decoder.at;
            let interpolation = match decoder.u8()? {
                0 => TableInterpolation::StepLeft,
                1 => TableInterpolation::Linear,
                tag => {
                    return Err(PayloadError::InvalidTag {
                        at: interpolation_at,
                        context: "table interpolation",
                        tag,
                    });
                }
            };
            let outside_at = decoder.at;
            let outside = match decoder.u8()? {
                0 => OutsideDomainPolicy::Refuse,
                1 => OutsideDomainPolicy::Clamp,
                2 => OutsideDomainPolicy::Periodic,
                tag => {
                    return Err(PayloadError::InvalidTag {
                        at: outside_at,
                        context: "outside-domain policy",
                        tag,
                    });
                }
            };
            let count = decoder.len("table samples")?;
            decoder.charge_items(count)?;
            let mut times = decoder.allocate_vec(count, "table times", 14)?;
            for _ in 0..count {
                times.push(decoder.qty()?);
            }
            let mut values = decoder.allocate_vec(count, "table values", minimum_value_bytes)?;
            for _ in 0..count {
                values.push(read_value(decoder)?);
            }
            SampleSource::table(times, values, interpolation, outside)
        }
        2 => {
            let family_at = decoder.at;
            let family = match decoder.u8()? {
                0 => DistributionFamily::Normal,
                1 => DistributionFamily::Uniform,
                2 => DistributionFamily::Empirical,
                tag => {
                    return Err(PayloadError::InvalidTag {
                        at: family_at,
                        context: "distribution family",
                        tag,
                    });
                }
            };
            let count = decoder.len("distribution values")?;
            let mut parameters =
                decoder.allocate_vec(count, "distribution values", minimum_value_bytes)?;
            for _ in 0..count {
                parameters.push(read_value(decoder)?);
            }
            SampleSource::distribution(family, parameters)
        }
        tag => Err(PayloadError::InvalidTag {
            at: source_at,
            context: "sample source",
            tag,
        }),
    }
}

fn read_qty_vec(decoder: &mut Decoder<'_>) -> Result<Vec<QtyAny>, PayloadError> {
    let count = decoder.len("quantity vector")?;
    let mut values = decoder.allocate_vec(count, "quantity vector", 14)?;
    for _ in 0..count {
        values.push(decoder.qty()?);
    }
    Ok(values)
}

fn read_phasor(decoder: &mut Decoder<'_>) -> Result<PhasorQty, PayloadError> {
    let real = decoder.qty()?;
    let imaginary = decoder.qty()?;
    let kind = read_quantity_kind(decoder)?;
    let amplitude_at = decoder.at;
    let amplitude = match decoder.u8()? {
        0 => PhasorAmplitude::Peak,
        1 => PhasorAmplitude::Rms,
        tag => {
            return Err(PayloadError::InvalidTag {
                at: amplitude_at,
                context: "phasor amplitude",
                tag,
            });
        }
    };
    Ok(PhasorQty::new(real, imaginary, kind, amplitude)?)
}

fn read_species_sample(decoder: &mut Decoder<'_>) -> Result<Vec<SpeciesValue>, PayloadError> {
    let count = decoder.len("species values")?;
    decoder.ensure_can_charge(count)?;
    let mut values = decoder.allocate_vec(count, "species values", 19)?;
    for _ in 0..count {
        values.push(SpeciesValue::new(decoder.species_id()?, decoder.qty()?));
    }
    Ok(values)
}

fn read_semantic_type(decoder: &mut Decoder<'_>) -> Result<SemanticType, PayloadError> {
    let kind = read_quantity_kind(decoder)?;
    let form_at = decoder.at;
    let form = match decoder.u8()? {
        0 => ValueForm::Static,
        1 => ValueForm::Instantaneous,
        2 => ValueForm::Peak,
        3 => ValueForm::Rms,
        tag => {
            return Err(PayloadError::InvalidTag {
                at: form_at,
                context: "semantic value form",
                tag,
            });
        }
    };
    Ok(SemanticType::new(kind, form))
}

fn read_quantity_kind(decoder: &mut Decoder<'_>) -> Result<QuantityKind, PayloadError> {
    let at = decoder.at;
    Ok(match decoder.u8()? {
        0 => QuantityKind::AbsoluteTemperature,
        1 => QuantityKind::TemperatureDifference,
        2 => QuantityKind::Angle(read_angle_domain(decoder)?),
        3 => QuantityKind::AngularVelocity(read_angle_domain(decoder)?),
        4 => QuantityKind::Torque,
        5 => QuantityKind::Energy,
        6 => QuantityKind::Pressure,
        7 => QuantityKind::Stress,
        8 => {
            let basis_at = decoder.at;
            let basis = match decoder.u8()? {
                0 => StrainBasis::Tensor,
                1 => StrainBasis::Engineering,
                tag => {
                    return Err(PayloadError::InvalidTag {
                        at: basis_at,
                        context: "strain basis",
                        tag,
                    });
                }
            };
            let component_at = decoder.at;
            let component = match decoder.u8()? {
                0 => StrainComponent::Normal,
                1 => StrainComponent::Shear,
                tag => {
                    return Err(PayloadError::InvalidTag {
                        at: component_at,
                        context: "strain component",
                        tag,
                    });
                }
            };
            QuantityKind::Strain { basis, component }
        }
        9 => {
            let basis_at = decoder.at;
            let basis = match decoder.u8()? {
                0 => CompositionBasis::MassFraction,
                1 => CompositionBasis::MoleFraction,
                2 => CompositionBasis::VolumeFraction,
                tag => {
                    return Err(PayloadError::InvalidTag {
                        at: basis_at,
                        context: "composition basis",
                        tag,
                    });
                }
            };
            QuantityKind::Composition(basis)
        }
        10 => QuantityKind::Mass,
        11 => QuantityKind::Amount,
        12 => QuantityKind::MolarMass,
        13 => QuantityKind::MassConcentration,
        14 => QuantityKind::AmountConcentration,
        15 => QuantityKind::Entropy,
        16 => QuantityKind::HeatCapacity,
        17 => QuantityKind::AcousticPressure,
        18 => QuantityKind::AcousticPower,
        tag => {
            return Err(PayloadError::InvalidTag {
                at,
                context: "semantic quantity kind",
                tag,
            });
        }
    })
}

fn read_angle_domain(decoder: &mut Decoder<'_>) -> Result<AngleDomain, PayloadError> {
    let at = decoder.at;
    match decoder.u8()? {
        0 => Ok(AngleDomain::Mechanical),
        1 => Ok(AngleDomain::Electrical),
        tag => Err(PayloadError::InvalidTag {
            at,
            context: "angle domain",
            tag,
        }),
    }
}
