//! Entity-bound sensor declarations and owner-neutral observation operators.
//!
//! A sensor is part of the scenario, not an anonymous measurement vector. The
//! declaration retains where the instrument is mounted, what it measures, how
//! mounting changes the reading, its placement uncertainty, its dynamic model,
//! and the calibration authority (or the fact that it is virtual).
//!
//! Compilation is deliberately geometry-neutral. The caller supplies an
//! admitted point component or patch-average restriction row; this module
//! checks and retains that row, applies the declared affine mount model, and
//! propagates explicitly supplied placement sensitivity into reading variance.
//! Deriving the row or sensitivity from a mesh/field remains the owning
//! geometry/solver layer's responsibility.

use core::fmt;

use fs_blake3::{ContentHash, DomainHasher, hash_bytes};
use fs_exec::Cx;

use crate::entity::{
    EntityCatalog, EntityId, EntityKind, EntityRef, EvidenceTier, KindExpectation, ResolutionFault,
};
use crate::ir::write_ir;
use crate::scenario::Scenario;

/// Schema version bound into every scenario-sensor identity.
pub const SENSOR_SCHEMA_VERSION: u32 = 1;

/// Domain-separated identity namespace for scenario sensors.
pub const SENSOR_IDENTITY_DOMAIN: &str = "org.frankensim.scenario.sensor.v1";

/// Schema version bound into every compiled sensor-set identity.
pub const SENSOR_SET_SCHEMA_VERSION: u32 = 1;

/// Domain-separated identity namespace for catalog-checked sensor sets.
pub const SENSOR_SET_IDENTITY_DOMAIN: &str = "org.frankensim.scenario.sensor-set.v1";

/// Version of the canonical scenario-bound sensor-extension envelope.
pub const SENSOR_EXTENSION_SCHEMA_VERSION: u16 = 1;

/// Domain-separated identity namespace for scenario-bound sensor extensions.
pub const SENSOR_EXTENSION_IDENTITY_DOMAIN: &str = "org.frankensim.scenario.sensor-extension.v1";

/// Generic Design Ledger artifact kind for canonical sensor extensions.
pub const SENSOR_EXTENSION_ARTIFACT_KIND: &str = "scenario-sensor-extension-v1";

/// Maximum admitted bytes in one canonical sensor-extension envelope.
pub const MAX_SENSOR_EXTENSION_WIRE_BYTES: usize = 64 * 1024 * 1024;

const SENSOR_EXTENSION_MAGIC: &[u8; 8] = b"FSSENX1\0";

/// Maximum dense state dimension accepted by the v1 operator compiler.
pub const MAX_SENSOR_STATE_DIMENSION: usize = 65_536;

/// Maximum support terms in one point/patch restriction row.
pub const MAX_SENSOR_SUPPORT_TERMS: usize = 4_096;

/// Maximum bytes in one retained sensor identity or authority string.
pub const MAX_SENSOR_TEXT_BYTES: usize = 4_096;

/// Conservative default admission limits for one catalog-checked sensor set.
pub const DEFAULT_SENSOR_SET_BUDGET: SensorSetBudget = SensorSetBudget {
    max_sensors: 4_096,
    max_work: 16_777_216,
};

const PATCH_WEIGHT_TOLERANCE: f64 = 64.0 * f64::EPSILON;

/// The closed v1 set of scenario sensor families.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SensorKind {
    /// Contact temperature sensor.
    Thermocouple,
    /// Resistance temperature detector.
    Rtd,
    /// Volumetric-flow instrument.
    FlowMeter,
    /// Static pressure tap.
    PressureTap,
    /// Area/patch temperature observation.
    IrCameraRegion,
}

impl SensorKind {
    const fn tag(self) -> u8 {
        match self {
            Self::Thermocouple => 1,
            Self::Rtd => 2,
            Self::FlowMeter => 3,
            Self::PressureTap => 4,
            Self::IrCameraRegion => 5,
        }
    }

    /// Stable human-readable label.
    #[must_use]
    pub const fn label(self) -> &'static str {
        match self {
            Self::Thermocouple => "thermocouple",
            Self::Rtd => "rtd",
            Self::FlowMeter => "flow-meter",
            Self::PressureTap => "pressure-tap",
            Self::IrCameraRegion => "ir-camera-region",
        }
    }

    /// Physical quantity measured by this family.
    #[must_use]
    pub const fn quantity(self) -> SensorQuantity {
        match self {
            Self::Thermocouple | Self::Rtd | Self::IrCameraRegion => SensorQuantity::Temperature,
            Self::FlowMeter => SensorQuantity::VolumetricFlow,
            Self::PressureTap => SensorQuantity::Pressure,
        }
    }

    const fn admits_entity(self, kind: EntityKind) -> bool {
        match self {
            Self::Thermocouple | Self::Rtd | Self::PressureTap => {
                matches!(kind, EntityKind::Region | EntityKind::Surface)
            }
            Self::FlowMeter => matches!(kind, EntityKind::Surface | EntityKind::Interface),
            Self::IrCameraRegion => matches!(kind, EntityKind::Surface),
        }
    }
}

/// The output quantity and SI dimensions of a sensor row.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SensorQuantity {
    /// Absolute or relative temperature, K.
    Temperature,
    /// Volume per time, m³/s.
    VolumetricFlow,
    /// Pressure, Pa.
    Pressure,
}

impl SensorQuantity {
    const fn tag(self) -> u8 {
        match self {
            Self::Temperature => 1,
            Self::VolumetricFlow => 2,
            Self::Pressure => 3,
        }
    }

    /// SI base exponents `[m, kg, s, K, A, mol]`.
    #[must_use]
    pub const fn dims(self) -> [i8; 6] {
        match self {
            Self::Temperature => [0, 0, 0, 1, 0, 0],
            Self::VolumetricFlow => [3, 0, -1, 0, 0, 0],
            Self::Pressure => [-1, 1, -2, 0, 0, 0],
        }
    }
}

/// Explicit placement uncertainty and its declared reading sensitivity.
///
/// `standard_uncertainty_m` is an axis-aligned one-standard-uncertainty
/// declaration in the entity's local coordinates. `reading_sensitivity_per_m`
/// is the caller-supplied local derivative of the final sensor reading with
/// respect to placement. This module propagates their diagonal first-order
/// variance; it does not derive either input from geometry.
#[derive(Debug, Clone, PartialEq)]
pub struct PlacementUncertainty {
    kind: PlacementUncertaintyKind,
}

#[derive(Debug, Clone, PartialEq)]
enum PlacementUncertaintyKind {
    DeclaredExact {
        source: String,
    },
    AxisAligned {
        standard_uncertainty_m: [f64; 3],
        reading_sensitivity_per_m: [f64; 3],
        source: String,
    },
}

impl PlacementUncertainty {
    /// Declare exact placement explicitly.
    ///
    /// # Errors
    ///
    /// Refuses a blank or oversized source identity.
    pub fn declared_exact(source: impl Into<String>) -> Result<Self, SensorError> {
        let source = source.into();
        validate_text("placement.source", &source)?;
        Ok(Self {
            kind: PlacementUncertaintyKind::DeclaredExact { source },
        })
    }

    /// Declare axis-aligned standard uncertainty and reading sensitivity.
    ///
    /// # Errors
    ///
    /// Refuses non-finite/negative standard uncertainties, non-finite
    /// sensitivities, an all-zero uncertainty declaration, or a malformed
    /// source identity.
    pub fn axis_aligned(
        standard_uncertainty_m: [f64; 3],
        reading_sensitivity_per_m: [f64; 3],
        source: impl Into<String>,
    ) -> Result<Self, SensorError> {
        let source = source.into();
        validate_text("placement.source", &source)?;
        if standard_uncertainty_m
            .iter()
            .any(|value| !value.is_finite() || *value < 0.0)
        {
            return Err(SensorError::InvalidPlacementUncertainty {
                what: "standard placement uncertainties must be finite and non-negative",
            });
        }
        if standard_uncertainty_m.iter().all(|value| *value == 0.0) {
            return Err(SensorError::InvalidPlacementUncertainty {
                what: "an all-zero declaration must use declared_exact so the authority is explicit",
            });
        }
        if reading_sensitivity_per_m
            .iter()
            .any(|value| !value.is_finite())
        {
            return Err(SensorError::InvalidPlacementUncertainty {
                what: "placement reading sensitivities must be finite",
            });
        }
        Ok(Self {
            kind: PlacementUncertaintyKind::AxisAligned {
                standard_uncertainty_m,
                reading_sensitivity_per_m,
                source,
            },
        })
    }

    /// Whether placement is explicitly declared exact for this model.
    #[must_use]
    pub const fn is_declared_exact(&self) -> bool {
        matches!(&self.kind, PlacementUncertaintyKind::DeclaredExact { .. })
    }

    /// Axis-aligned standard placement uncertainty in metres, when declared.
    #[must_use]
    pub const fn standard_uncertainty_m(&self) -> Option<[f64; 3]> {
        match &self.kind {
            PlacementUncertaintyKind::DeclaredExact { .. } => None,
            PlacementUncertaintyKind::AxisAligned {
                standard_uncertainty_m,
                ..
            } => Some(*standard_uncertainty_m),
        }
    }

    /// Declared final-reading sensitivity per metre, when applicable.
    #[must_use]
    pub const fn reading_sensitivity_per_m(&self) -> Option<[f64; 3]> {
        match &self.kind {
            PlacementUncertaintyKind::DeclaredExact { .. } => None,
            PlacementUncertaintyKind::AxisAligned {
                reading_sensitivity_per_m,
                ..
            } => Some(*reading_sensitivity_per_m),
        }
    }

    /// Named uncertainty authority.
    #[must_use]
    pub fn source(&self) -> &str {
        match &self.kind {
            PlacementUncertaintyKind::DeclaredExact { source }
            | PlacementUncertaintyKind::AxisAligned { source, .. } => source,
        }
    }

    /// Propagated first-order reading variance.
    ///
    /// # Errors
    ///
    /// Refuses finite-input overflow.
    pub fn propagated_variance(&self) -> Result<f64, SensorError> {
        match &self.kind {
            PlacementUncertaintyKind::DeclaredExact { .. } => Ok(0.0),
            PlacementUncertaintyKind::AxisAligned {
                standard_uncertainty_m,
                reading_sensitivity_per_m,
                ..
            } => {
                let mut variance = 0.0;
                for (sigma, sensitivity) in
                    standard_uncertainty_m.iter().zip(reading_sensitivity_per_m)
                {
                    let contribution = sigma * sensitivity;
                    variance = contribution.mul_add(contribution, variance);
                    if !variance.is_finite() {
                        return Err(SensorError::NonFiniteComputation {
                            operation: "placement variance propagation",
                        });
                    }
                }
                Ok(variance)
            }
        }
    }
}

/// Entity-local sensor support.
#[derive(Debug, Clone, PartialEq)]
pub struct SensorLocation {
    entity: EntityRef,
    local_position_m: [f64; 3],
    placement: PlacementUncertainty,
}

impl SensorLocation {
    /// Construct one entity-bound local support.
    ///
    /// # Errors
    ///
    /// Refuses a non-finite local coordinate or a reference whose declared
    /// kind expectation does not admit its identity's embedded kind.
    pub fn new(
        entity: EntityRef,
        local_position_m: [f64; 3],
        placement: PlacementUncertainty,
    ) -> Result<Self, SensorError> {
        if local_position_m.iter().any(|value| !value.is_finite()) {
            return Err(SensorError::NonFiniteField {
                field: "location.local_position_m",
            });
        }
        if !entity.expect().admits(entity.target().kind()) {
            return Err(SensorError::EntityExpectationMismatch {
                target: entity.target(),
                expected: entity.expect(),
            });
        }
        Ok(Self {
            entity,
            local_position_m: canonicalize_zero3(local_position_m),
            placement,
        })
    }

    /// Persistent entity reference.
    #[must_use]
    pub const fn entity(&self) -> EntityRef {
        self.entity
    }

    /// Local coordinates in metres.
    #[must_use]
    pub const fn local_position_m(&self) -> [f64; 3] {
        self.local_position_m
    }

    /// Placement uncertainty declaration.
    #[must_use]
    pub const fn placement(&self) -> &PlacementUncertainty {
        &self.placement
    }
}

/// One sparse term in a point/patch restriction row.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ObservationTerm {
    component: usize,
    weight: f64,
}

impl ObservationTerm {
    /// Construct one weighted state component.
    #[must_use]
    pub const fn new(component: usize, weight: f64) -> Self {
        Self { component, weight }
    }

    /// State component.
    #[must_use]
    pub const fn component(self) -> usize {
        self.component
    }

    /// Restriction weight.
    #[must_use]
    pub const fn weight(self) -> f64 {
        self.weight
    }
}

/// Declared restriction from a solver state to one raw field sample.
#[derive(Debug, Clone, PartialEq)]
pub struct ObservationSupport {
    state_dimension: usize,
    kind: ObservationSupportKind,
}

#[derive(Debug, Clone, PartialEq)]
enum ObservationSupportKind {
    Point { component: usize },
    PatchAverage { terms: Vec<ObservationTerm> },
}

impl ObservationSupport {
    /// Construct a point selector.
    ///
    /// # Errors
    ///
    /// Refuses a zero/oversized dimension or out-of-range component.
    pub fn point(state_dimension: usize, component: usize) -> Result<Self, SensorError> {
        validate_state_dimension(state_dimension)?;
        if component >= state_dimension {
            return Err(SensorError::ComponentOutOfRange {
                component,
                state_dimension,
            });
        }
        Ok(Self {
            state_dimension,
            kind: ObservationSupportKind::Point { component },
        })
    }

    /// Construct a checked patch-average row.
    ///
    /// # Errors
    ///
    /// Refuses empty/oversized support, out-of-range or duplicate components,
    /// non-finite/non-positive weights, or weights that do not sum to one
    /// within the documented floating admission tolerance.
    pub fn patch_average(
        state_dimension: usize,
        mut terms: Vec<ObservationTerm>,
    ) -> Result<Self, SensorError> {
        validate_state_dimension(state_dimension)?;
        if terms.is_empty() || terms.len() > MAX_SENSOR_SUPPORT_TERMS {
            return Err(SensorError::SupportSize {
                actual: terms.len(),
                limit: MAX_SENSOR_SUPPORT_TERMS,
            });
        }
        for term in &terms {
            if term.component >= state_dimension {
                return Err(SensorError::ComponentOutOfRange {
                    component: term.component,
                    state_dimension,
                });
            }
            if !term.weight.is_finite() || term.weight <= 0.0 {
                return Err(SensorError::InvalidSupportWeight {
                    component: term.component,
                });
            }
        }
        terms.sort_unstable_by_key(|term| term.component);
        if let Some(component) = terms.windows(2).find_map(|pair| {
            if pair[0].component == pair[1].component {
                Some(pair[0].component)
            } else {
                None
            }
        }) {
            return Err(SensorError::DuplicateComponent { component });
        }
        let mut weight_sum = 0.0;
        for term in &terms {
            weight_sum += term.weight;
            if !weight_sum.is_finite() {
                return Err(SensorError::NonFiniteComputation {
                    operation: "patch weight accumulation",
                });
            }
        }
        if (weight_sum - 1.0).abs() > PATCH_WEIGHT_TOLERANCE {
            return Err(SensorError::PatchWeightsDoNotSumToOne { actual: weight_sum });
        }
        Ok(Self {
            state_dimension,
            kind: ObservationSupportKind::PatchAverage { terms },
        })
    }

    /// Dense state dimension.
    #[must_use]
    pub const fn state_dimension(&self) -> usize {
        self.state_dimension
    }

    /// Whether this support is a non-negative, unit-sum patch average.
    #[must_use]
    pub const fn is_patch_average(&self) -> bool {
        matches!(&self.kind, ObservationSupportKind::PatchAverage { .. })
    }

    /// Sparse declaration terms.
    #[must_use]
    pub fn terms(&self) -> ObservationTerms<'_> {
        match &self.kind {
            ObservationSupportKind::Point { component } => ObservationTerms::Point {
                component: *component,
                emitted: false,
            },
            ObservationSupportKind::PatchAverage { terms } => ObservationTerms::Patch(terms.iter()),
        }
    }
}

/// Borrowing iterator over declared observation terms.
pub enum ObservationTerms<'a> {
    /// One point term.
    Point {
        /// Selected component.
        component: usize,
        /// Whether the term was emitted.
        emitted: bool,
    },
    /// Patch terms.
    Patch(core::slice::Iter<'a, ObservationTerm>),
}

impl Iterator for ObservationTerms<'_> {
    type Item = ObservationTerm;

    fn next(&mut self) -> Option<Self::Item> {
        match self {
            Self::Point { component, emitted } => {
                if *emitted {
                    None
                } else {
                    *emitted = true;
                    Some(ObservationTerm::new(*component, 1.0))
                }
            }
            Self::Patch(terms) => terms.next().copied(),
        }
    }
}

/// How mounting maps the raw field sample to the sensor's steady reading.
#[derive(Debug, Clone, PartialEq)]
pub struct SensorMount {
    kind: SensorMountKind,
}

#[derive(Debug, Clone, PartialEq)]
enum SensorMountKind {
    DeclaredIdeal {
        source: String,
    },
    Affine {
        gain: f64,
        offset: f64,
        source: String,
    },
}

impl SensorMount {
    /// Declare an ideal mount explicitly.
    ///
    /// # Errors
    ///
    /// Refuses a malformed source identity.
    pub fn declared_ideal(source: impl Into<String>) -> Result<Self, SensorError> {
        let source = source.into();
        validate_text("mount.source", &source)?;
        Ok(Self {
            kind: SensorMountKind::DeclaredIdeal { source },
        })
    }

    /// Declare an affine mount correction.
    ///
    /// # Errors
    ///
    /// Refuses a non-finite/non-positive gain, non-finite offset, or malformed
    /// source identity.
    pub fn affine(gain: f64, offset: f64, source: impl Into<String>) -> Result<Self, SensorError> {
        let source = source.into();
        validate_text("mount.source", &source)?;
        if !gain.is_finite() || gain <= 0.0 {
            return Err(SensorError::InvalidMount {
                what: "affine mount gain must be finite and positive",
            });
        }
        if !offset.is_finite() {
            return Err(SensorError::InvalidMount {
                what: "affine mount offset must be finite",
            });
        }
        Ok(Self {
            kind: SensorMountKind::Affine {
                gain,
                offset: canonicalize_zero(offset),
                source,
            },
        })
    }

    /// Whether mounting is explicitly idealized as reading-preserving.
    #[must_use]
    pub const fn is_declared_ideal(&self) -> bool {
        matches!(&self.kind, SensorMountKind::DeclaredIdeal { .. })
    }

    /// `(gain, offset)` used by the steady operator.
    #[must_use]
    pub const fn affine_parts(&self) -> (f64, f64) {
        match &self.kind {
            SensorMountKind::DeclaredIdeal { .. } => (1.0, 0.0),
            SensorMountKind::Affine { gain, offset, .. } => (*gain, *offset),
        }
    }

    /// Named mount authority.
    #[must_use]
    pub fn source(&self) -> &str {
        match &self.kind {
            SensorMountKind::DeclaredIdeal { source } | SensorMountKind::Affine { source, .. } => {
                source
            }
        }
    }
}

/// Declared sensor dynamics.
#[derive(Debug, Clone, PartialEq)]
pub struct SensorDynamics {
    kind: SensorDynamicsKind,
}

#[derive(Debug, Clone, PartialEq)]
enum SensorDynamicsKind {
    DeclaredInstantaneous {
        source: String,
    },
    FirstOrder {
        time_constant_s: f64,
        source: String,
    },
}

impl SensorDynamics {
    /// Declare instantaneous response explicitly.
    ///
    /// # Errors
    ///
    /// Refuses a malformed source identity.
    pub fn declared_instantaneous(source: impl Into<String>) -> Result<Self, SensorError> {
        let source = source.into();
        validate_text("dynamics.source", &source)?;
        Ok(Self {
            kind: SensorDynamicsKind::DeclaredInstantaneous { source },
        })
    }

    /// Declare first-order response.
    ///
    /// # Errors
    ///
    /// Refuses a non-finite/non-positive time constant or malformed source.
    pub fn first_order(
        time_constant_s: f64,
        source: impl Into<String>,
    ) -> Result<Self, SensorError> {
        let source = source.into();
        validate_text("dynamics.source", &source)?;
        if !time_constant_s.is_finite() || time_constant_s <= 0.0 {
            return Err(SensorError::InvalidDynamics {
                what: "first-order time constant must be finite and positive",
            });
        }
        Ok(Self {
            kind: SensorDynamicsKind::FirstOrder {
                time_constant_s,
                source,
            },
        })
    }

    /// Whether dynamics are explicitly idealized as instantaneous.
    #[must_use]
    pub const fn is_declared_instantaneous(&self) -> bool {
        matches!(&self.kind, SensorDynamicsKind::DeclaredInstantaneous { .. })
    }

    /// Time constant, or `None` for an explicitly instantaneous declaration.
    #[must_use]
    pub const fn time_constant_s(&self) -> Option<f64> {
        match &self.kind {
            SensorDynamicsKind::DeclaredInstantaneous { .. } => None,
            SensorDynamicsKind::FirstOrder {
                time_constant_s, ..
            } => Some(*time_constant_s),
        }
    }

    /// Named dynamics authority.
    #[must_use]
    pub fn source(&self) -> &str {
        match &self.kind {
            SensorDynamicsKind::DeclaredInstantaneous { source }
            | SensorDynamicsKind::FirstOrder { source, .. } => source,
        }
    }
}

/// Calibration authority for a physical or virtual sensor.
#[derive(Debug, Clone, PartialEq)]
pub struct SensorCalibration {
    kind: SensorCalibrationKind,
}

#[derive(Debug, Clone, PartialEq)]
enum SensorCalibrationKind {
    Physical {
        certificate_ref: String,
        date: String,
        source: String,
        instrument_variance: f64,
    },
    Virtual {
        definition_ref: String,
    },
}

impl SensorCalibration {
    /// Construct a physical calibration record.
    ///
    /// # Errors
    ///
    /// Refuses malformed identities, an invalid calendar date, or non-finite/
    /// non-positive instrument variance.
    pub fn physical(
        certificate_ref: impl Into<String>,
        date: impl Into<String>,
        source: impl Into<String>,
        instrument_variance: f64,
    ) -> Result<Self, SensorError> {
        let certificate_ref = certificate_ref.into();
        let date = date.into();
        let source = source.into();
        validate_text("calibration.certificate_ref", &certificate_ref)?;
        validate_text("calibration.source", &source)?;
        validate_date(&date)?;
        if !instrument_variance.is_finite() || instrument_variance <= 0.0 {
            return Err(SensorError::InvalidInstrumentVariance);
        }
        Ok(Self {
            kind: SensorCalibrationKind::Physical {
                certificate_ref,
                date,
                source,
                instrument_variance,
            },
        })
    }

    /// Construct a virtual sensor record.
    ///
    /// # Errors
    ///
    /// Refuses a malformed definition identity.
    pub fn virtual_sensor(definition_ref: impl Into<String>) -> Result<Self, SensorError> {
        let definition_ref = definition_ref.into();
        validate_text("calibration.definition_ref", &definition_ref)?;
        Ok(Self {
            kind: SensorCalibrationKind::Virtual { definition_ref },
        })
    }

    /// Whether this is a virtual probe.
    #[must_use]
    pub const fn is_virtual(&self) -> bool {
        matches!(&self.kind, SensorCalibrationKind::Virtual { .. })
    }

    /// Positive instrument variance for a physical sensor.
    #[must_use]
    pub const fn instrument_variance(&self) -> Option<f64> {
        match &self.kind {
            SensorCalibrationKind::Physical {
                instrument_variance,
                ..
            } => Some(*instrument_variance),
            SensorCalibrationKind::Virtual { .. } => None,
        }
    }

    /// Physical calibration certificate identity, when applicable.
    #[must_use]
    pub fn certificate_ref(&self) -> Option<&str> {
        match &self.kind {
            SensorCalibrationKind::Physical {
                certificate_ref, ..
            } => Some(certificate_ref),
            SensorCalibrationKind::Virtual { .. } => None,
        }
    }

    /// Physical calibration date (`YYYY-MM-DD`), when applicable.
    #[must_use]
    pub fn date(&self) -> Option<&str> {
        match &self.kind {
            SensorCalibrationKind::Physical { date, .. } => Some(date),
            SensorCalibrationKind::Virtual { .. } => None,
        }
    }

    /// Calibration/provider source, when applicable.
    #[must_use]
    pub fn source(&self) -> Option<&str> {
        match &self.kind {
            SensorCalibrationKind::Physical { source, .. } => Some(source),
            SensorCalibrationKind::Virtual { .. } => None,
        }
    }

    /// Virtual probe definition/model reference, when applicable.
    #[must_use]
    pub fn definition_ref(&self) -> Option<&str> {
        match &self.kind {
            SensorCalibrationKind::Physical { .. } => None,
            SensorCalibrationKind::Virtual { definition_ref } => Some(definition_ref),
        }
    }
}

/// One complete scenario sensor declaration.
#[derive(Debug, Clone, PartialEq)]
pub struct ScenarioSensor {
    name: String,
    kind: SensorKind,
    location: SensorLocation,
    support: ObservationSupport,
    mount: SensorMount,
    dynamics: SensorDynamics,
    calibration: SensorCalibration,
    placement_candidate: bool,
}

impl ScenarioSensor {
    /// Admit one complete sensor declaration.
    ///
    /// # Errors
    ///
    /// Refuses a malformed name or a sensor/entity-kind combination that the
    /// closed v1 family table does not support.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        name: impl Into<String>,
        kind: SensorKind,
        location: SensorLocation,
        support: ObservationSupport,
        mount: SensorMount,
        dynamics: SensorDynamics,
        calibration: SensorCalibration,
        placement_candidate: bool,
    ) -> Result<Self, SensorError> {
        let name = name.into();
        validate_text("sensor.name", &name)?;
        let entity_kind = location.entity.target().kind();
        if !kind.admits_entity(entity_kind) {
            return Err(SensorError::UnsupportedEntityKind { kind, entity_kind });
        }
        Ok(Self {
            name,
            kind,
            location,
            support,
            mount,
            dynamics,
            calibration,
            placement_candidate,
        })
    }

    /// Declared sensor name.
    #[must_use]
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Sensor family.
    #[must_use]
    pub const fn kind(&self) -> SensorKind {
        self.kind
    }

    /// Measured quantity.
    #[must_use]
    pub const fn quantity(&self) -> SensorQuantity {
        self.kind.quantity()
    }

    /// Entity-local support.
    #[must_use]
    pub const fn location(&self) -> &SensorLocation {
        &self.location
    }

    /// Declared point/patch restriction.
    #[must_use]
    pub const fn support(&self) -> &ObservationSupport {
        &self.support
    }

    /// Mount model.
    #[must_use]
    pub const fn mount(&self) -> &SensorMount {
        &self.mount
    }

    /// Dynamics declaration.
    #[must_use]
    pub const fn dynamics(&self) -> &SensorDynamics {
        &self.dynamics
    }

    /// Calibration declaration.
    #[must_use]
    pub const fn calibration(&self) -> &SensorCalibration {
        &self.calibration
    }

    /// Whether this sensor is a placement candidate.
    #[must_use]
    pub const fn is_placement_candidate(&self) -> bool {
        self.placement_candidate
    }

    /// Deterministic identity over every semantic declaration field.
    #[must_use]
    pub fn identity(&self) -> ContentHash {
        let mut hasher = DomainHasher::new(SENSOR_IDENTITY_DOMAIN);
        absorb_u64(&mut hasher, u64::from(SENSOR_SCHEMA_VERSION));
        absorb_bytes(&mut hasher, self.name.as_bytes());
        hasher.update(&[self.kind.tag(), self.kind.quantity().tag()]);
        absorb_entity_ref(&mut hasher, self.location.entity);
        absorb_f64s(&mut hasher, &self.location.local_position_m);
        absorb_placement(&mut hasher, &self.location.placement);
        absorb_support(&mut hasher, &self.support);
        absorb_mount(&mut hasher, &self.mount);
        absorb_dynamics(&mut hasher, &self.dynamics);
        absorb_calibration(&mut hasher, &self.calibration);
        hasher.update(&[u8::from(self.placement_candidate)]);
        hasher.finalize()
    }

    /// Compile the declaration into a dense steady observation row.
    ///
    /// # Errors
    ///
    /// Refuses allocation failure, finite-input overflow while applying the
    /// mount gain, or placement-variance overflow.
    pub fn compile(&self) -> Result<CompiledSensorOperator, SensorError> {
        let state_dimension = self.support.state_dimension();
        let mut operator = Vec::new();
        operator.try_reserve_exact(state_dimension).map_err(|_| {
            SensorError::AllocationRefused {
                resource: "compiled dense sensor operator",
            }
        })?;
        operator.resize(state_dimension, 0.0);
        let (gain, offset) = self.mount.affine_parts();
        for term in self.support.terms() {
            let weight = gain * term.weight;
            if !weight.is_finite() {
                return Err(SensorError::NonFiniteComputation {
                    operation: "mount gain application",
                });
            }
            operator[term.component] = canonicalize_zero(weight);
        }
        let placement_variance = self.location.placement.propagated_variance()?;
        Ok(CompiledSensorOperator {
            sensor_identity: self.identity(),
            name: self.name.clone(),
            quantity: self.quantity(),
            entity: self.location.entity.target(),
            local_position_m: self.location.local_position_m,
            operator,
            offset,
            placement_variance,
            instrument_variance: self.calibration.instrument_variance(),
            virtual_sensor: self.calibration.is_virtual(),
            placement_candidate: self.placement_candidate,
        })
    }
}

/// Canonical sidecar binding checked sensor declarations to one exact scenario.
///
/// The sidecar deliberately preserves scenario IR v2. Its scenario root is the
/// BLAKE3 content hash of [`write_ir`] output, and its own semantic identity
/// binds that root plus every ordered [`ScenarioSensor::identity`]. Decoding
/// reconstructs every declaration through the same checked constructors used
/// at authoring time.
#[derive(Debug, Clone, PartialEq)]
pub struct ScenarioSensorExtension {
    scenario_ir_hash: ContentHash,
    sensors: Vec<ScenarioSensor>,
    identity: ContentHash,
}

impl ScenarioSensorExtension {
    /// Bind an ordered sensor collection to the canonical IR of `scenario`.
    ///
    /// # Errors
    ///
    /// Refuses more than [`DEFAULT_SENSOR_SET_BUDGET`] sensors or duplicate
    /// exact names. Caller order is semantic.
    pub fn new(
        scenario: &Scenario,
        sensors: Vec<ScenarioSensor>,
    ) -> Result<Self, SensorExtensionError> {
        let canonical_ir = write_ir(scenario);
        Self::from_hash(hash_bytes(canonical_ir.as_bytes()), sensors)
    }

    fn from_hash(
        scenario_ir_hash: ContentHash,
        sensors: Vec<ScenarioSensor>,
    ) -> Result<Self, SensorExtensionError> {
        if sensors.len() > DEFAULT_SENSOR_SET_BUDGET.max_sensors {
            return Err(SensorExtensionError::LimitExceeded {
                resource: "sensor declarations",
                requested: sensors.len(),
                limit: DEFAULT_SENSOR_SET_BUDGET.max_sensors,
            });
        }
        for second in 1..sensors.len() {
            if let Some(first) = sensors[..second]
                .iter()
                .position(|sensor| sensor.name() == sensors[second].name())
            {
                return Err(SensorExtensionError::DuplicateName { first, second });
            }
        }
        let identity = sensor_extension_identity(scenario_ir_hash, &sensors);
        Ok(Self {
            scenario_ir_hash,
            sensors,
            identity,
        })
    }

    /// Hash of the exact canonical scenario IR this extension accompanies.
    #[must_use]
    pub const fn scenario_ir_hash(&self) -> ContentHash {
        self.scenario_ir_hash
    }

    /// Ordered checked sensor declarations.
    #[must_use]
    pub fn sensors(&self) -> &[ScenarioSensor] {
        &self.sensors
    }

    /// Domain-separated identity of the scenario root and ordered sensors.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    /// Whether `scenario` emits the exact canonical IR root bound here.
    #[must_use]
    pub fn verifies_scenario(&self, scenario: &Scenario) -> bool {
        hash_bytes(write_ir(scenario).as_bytes()) == self.scenario_ir_hash
    }

    /// Materialize the canonical bounded wire envelope.
    ///
    /// # Errors
    ///
    /// Refuses checked length overflow, the wire-byte ceiling, or exact output
    /// reservation failure before publishing bytes.
    pub fn try_encode(&self) -> Result<Vec<u8>, SensorExtensionError> {
        let encoded_len = sensor_extension_wire_len(self)?;
        let mut out = Vec::new();
        out.try_reserve_exact(encoded_len).map_err(|_| {
            SensorExtensionError::AllocationRefused {
                resource: "sensor extension wire bytes",
                requested: encoded_len,
            }
        })?;
        out.extend_from_slice(SENSOR_EXTENSION_MAGIC);
        out.extend_from_slice(&SENSOR_EXTENSION_SCHEMA_VERSION.to_le_bytes());
        out.extend_from_slice(self.scenario_ir_hash.as_bytes());
        push_wire_u32(&mut out, self.sensors.len())?;
        for sensor in &self.sensors {
            encode_sensor(&mut out, sensor)?;
            out.extend_from_slice(sensor.identity().as_bytes());
        }
        out.extend_from_slice(self.identity.as_bytes());
        debug_assert_eq!(out.len(), encoded_len);
        Ok(out)
    }

    /// Decode one complete canonical bounded wire envelope.
    ///
    /// Every declaration is reconstructed through its checked constructors.
    /// Embedded per-sensor and extension identities are recomputed, and the
    /// final byte-for-byte re-encoding check rejects alternate spellings such
    /// as negative zero.
    ///
    /// # Errors
    ///
    /// Refuses oversized, truncated, future-version, malformed-tag,
    /// constructor-invalid, identity-mismatched, noncanonical, or trailing
    /// input without publishing a partial extension.
    pub fn decode(bytes: &[u8]) -> Result<Self, SensorExtensionError> {
        if bytes.len() > MAX_SENSOR_EXTENSION_WIRE_BYTES {
            return Err(SensorExtensionError::LimitExceeded {
                resource: "sensor extension wire bytes",
                requested: bytes.len(),
                limit: MAX_SENSOR_EXTENSION_WIRE_BYTES,
            });
        }
        let mut reader = SensorWireReader::new(bytes);
        if reader.read_array::<8>("magic")? != *SENSOR_EXTENSION_MAGIC {
            return Err(SensorExtensionError::BadMagic);
        }
        let observed_version = reader.read_u16("schema version")?;
        if observed_version != SENSOR_EXTENSION_SCHEMA_VERSION {
            return Err(SensorExtensionError::UnsupportedVersion {
                observed: observed_version,
                supported: SENSOR_EXTENSION_SCHEMA_VERSION,
            });
        }
        let scenario_ir_hash = reader.read_hash("scenario IR hash")?;
        let sensor_count =
            reader.read_count("sensor declarations", DEFAULT_SENSOR_SET_BUDGET.max_sensors)?;
        let mut sensors = Vec::new();
        sensors.try_reserve_exact(sensor_count).map_err(|_| {
            SensorExtensionError::AllocationRefused {
                resource: "decoded sensor declarations",
                requested: sensor_count,
            }
        })?;
        for row in 0..sensor_count {
            let sensor = decode_sensor(&mut reader, row)?;
            let embedded_identity = reader.read_hash("sensor identity")?;
            if embedded_identity != sensor.identity() {
                return Err(SensorExtensionError::SensorIdentityMismatch { row });
            }
            sensors.push(sensor);
        }
        let embedded_identity = reader.read_hash("extension identity")?;
        if reader.remaining() != 0 {
            return Err(SensorExtensionError::TrailingBytes {
                count: reader.remaining(),
            });
        }
        let extension = Self::from_hash(scenario_ir_hash, sensors)?;
        if extension.identity != embedded_identity {
            return Err(SensorExtensionError::ExtensionIdentityMismatch);
        }
        if extension.try_encode()?.as_slice() != bytes {
            return Err(SensorExtensionError::NonCanonical);
        }
        Ok(extension)
    }
}

/// Canonical sensor-extension admission or wire refusal.
#[derive(Debug, Clone, PartialEq)]
pub enum SensorExtensionError {
    /// A named count or byte limit was exceeded.
    LimitExceeded {
        /// Stable resource name.
        resource: &'static str,
        /// Observed or requested units.
        requested: usize,
        /// Admitted units.
        limit: usize,
    },
    /// Checked wire-length arithmetic overflowed.
    LengthOverflow,
    /// A fallible exact reservation was refused.
    AllocationRefused {
        /// Stable resource name.
        resource: &'static str,
        /// Requested elements or bytes.
        requested: usize,
    },
    /// The fixed envelope magic did not match.
    BadMagic,
    /// The wire declares a schema this build does not implement.
    UnsupportedVersion {
        /// Version found on the wire.
        observed: u16,
        /// Only version admitted by this build.
        supported: u16,
    },
    /// Input ended before one complete field was available.
    Truncated {
        /// Stable field name.
        field: &'static str,
    },
    /// A closed wire tag used an unknown value.
    InvalidTag {
        /// Stable tagged field name.
        field: &'static str,
        /// Unknown tag.
        observed: u8,
    },
    /// A length-prefixed text field was not UTF-8.
    InvalidUtf8 {
        /// Stable text field name.
        field: &'static str,
    },
    /// One reconstructed declaration failed its checked constructor.
    InvalidSensor {
        /// Sensor row.
        row: usize,
        /// Constructor refusal.
        source: SensorError,
    },
    /// Two ordered declarations repeat one exact name.
    DuplicateName {
        /// First row carrying the name.
        first: usize,
        /// Later row carrying the same name.
        second: usize,
    },
    /// A stored per-sensor identity does not match the decoded declaration.
    SensorIdentityMismatch {
        /// Sensor row.
        row: usize,
    },
    /// The stored extension identity does not match its semantic fields.
    ExtensionIdentityMismatch,
    /// Input decoded semantically but was not the one canonical byte spelling.
    NonCanonical,
    /// Complete canonical input was followed by extra bytes.
    TrailingBytes {
        /// Unconsumed byte count.
        count: usize,
    },
}

impl fmt::Display for SensorExtensionError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::LimitExceeded {
                resource,
                requested,
                limit,
            } => write!(
                formatter,
                "sensor extension {resource} request {requested} exceeds limit {limit}"
            ),
            Self::LengthOverflow => formatter.write_str("sensor extension wire length overflowed"),
            Self::AllocationRefused {
                resource,
                requested,
            } => write!(
                formatter,
                "sensor extension allocation for {requested} {resource} units was refused"
            ),
            Self::BadMagic => formatter.write_str("invalid sensor extension magic"),
            Self::UnsupportedVersion {
                observed,
                supported,
            } => write!(
                formatter,
                "unsupported sensor extension schema {observed}; supported schema is {supported}"
            ),
            Self::Truncated { field } => {
                write!(formatter, "sensor extension ended inside {field}")
            }
            Self::InvalidTag { field, observed } => {
                write!(
                    formatter,
                    "sensor extension {field} has unknown tag {observed}"
                )
            }
            Self::InvalidUtf8 { field } => {
                write!(formatter, "sensor extension {field} is not UTF-8")
            }
            Self::InvalidSensor { row, source } => {
                write!(formatter, "sensor extension row {row} is invalid: {source}")
            }
            Self::DuplicateName { first, second } => write!(
                formatter,
                "sensor extension rows {first} and {second} repeat the same exact name"
            ),
            Self::SensorIdentityMismatch { row } => {
                write!(
                    formatter,
                    "sensor extension row {row} identity does not verify"
                )
            }
            Self::ExtensionIdentityMismatch => {
                formatter.write_str("sensor extension identity does not verify")
            }
            Self::NonCanonical => {
                formatter.write_str("sensor extension has a noncanonical byte spelling")
            }
            Self::TrailingBytes { count } => {
                write!(formatter, "sensor extension has {count} trailing bytes")
            }
        }
    }
}

impl core::error::Error for SensorExtensionError {}

fn sensor_extension_identity(
    scenario_ir_hash: ContentHash,
    sensors: &[ScenarioSensor],
) -> ContentHash {
    let mut hasher = DomainHasher::new(SENSOR_EXTENSION_IDENTITY_DOMAIN);
    absorb_u64(&mut hasher, u64::from(SENSOR_EXTENSION_SCHEMA_VERSION));
    hasher.update(scenario_ir_hash.as_bytes());
    absorb_u64(&mut hasher, sensors.len() as u64);
    for sensor in sensors {
        hasher.update(sensor.identity().as_bytes());
    }
    hasher.finalize()
}

fn sensor_extension_wire_len(
    extension: &ScenarioSensorExtension,
) -> Result<usize, SensorExtensionError> {
    let mut length = SENSOR_EXTENSION_MAGIC.len()
        + core::mem::size_of::<u16>()
        + 32
        + core::mem::size_of::<u32>()
        + 32;
    for sensor in &extension.sensors {
        checked_wire_add(&mut length, sensor_wire_len(sensor)?)?;
        checked_wire_add(&mut length, 32)?;
    }
    if length > MAX_SENSOR_EXTENSION_WIRE_BYTES {
        return Err(SensorExtensionError::LimitExceeded {
            resource: "sensor extension wire bytes",
            requested: length,
            limit: MAX_SENSOR_EXTENSION_WIRE_BYTES,
        });
    }
    Ok(length)
}

fn sensor_wire_len(sensor: &ScenarioSensor) -> Result<usize, SensorExtensionError> {
    let mut length = wire_text_len(&sensor.name)?;
    checked_wire_add(&mut length, 1 + 1 + 32 + 1 + 3 * 8)?;
    checked_wire_add(
        &mut length,
        match &sensor.location.placement.kind {
            PlacementUncertaintyKind::DeclaredExact { source } => 1 + wire_text_len(source)?,
            PlacementUncertaintyKind::AxisAligned { source, .. } => {
                1 + 6 * 8 + wire_text_len(source)?
            }
        },
    )?;
    checked_wire_add(
        &mut length,
        match &sensor.support.kind {
            ObservationSupportKind::Point { .. } => 1 + 2 * 4,
            ObservationSupportKind::PatchAverage { terms } => 1 + 2 * 4 + terms.len() * (4 + 8),
        },
    )?;
    checked_wire_add(
        &mut length,
        match &sensor.mount.kind {
            SensorMountKind::DeclaredIdeal { source } => 1 + wire_text_len(source)?,
            SensorMountKind::Affine { source, .. } => 1 + 2 * 8 + wire_text_len(source)?,
        },
    )?;
    checked_wire_add(
        &mut length,
        match &sensor.dynamics.kind {
            SensorDynamicsKind::DeclaredInstantaneous { source } => 1 + wire_text_len(source)?,
            SensorDynamicsKind::FirstOrder { source, .. } => 1 + 8 + wire_text_len(source)?,
        },
    )?;
    checked_wire_add(
        &mut length,
        match &sensor.calibration.kind {
            SensorCalibrationKind::Physical {
                certificate_ref,
                date,
                source,
                ..
            } => {
                1 + wire_text_len(certificate_ref)?
                    + wire_text_len(date)?
                    + wire_text_len(source)?
                    + 8
            }
            SensorCalibrationKind::Virtual { definition_ref } => 1 + wire_text_len(definition_ref)?,
        },
    )?;
    checked_wire_add(&mut length, 1)?;
    Ok(length)
}

fn wire_text_len(text: &str) -> Result<usize, SensorExtensionError> {
    core::mem::size_of::<u32>()
        .checked_add(text.len())
        .ok_or(SensorExtensionError::LengthOverflow)
}

fn checked_wire_add(total: &mut usize, value: usize) -> Result<(), SensorExtensionError> {
    *total = total
        .checked_add(value)
        .ok_or(SensorExtensionError::LengthOverflow)?;
    Ok(())
}

fn encode_sensor(out: &mut Vec<u8>, sensor: &ScenarioSensor) -> Result<(), SensorExtensionError> {
    push_wire_text(out, &sensor.name)?;
    out.push(sensor.kind.tag());
    encode_entity_ref(out, sensor.location.entity);
    for value in sensor.location.local_position_m {
        push_wire_f64(out, value);
    }
    match &sensor.location.placement.kind {
        PlacementUncertaintyKind::DeclaredExact { source } => {
            out.push(1);
            push_wire_text(out, source)?;
        }
        PlacementUncertaintyKind::AxisAligned {
            standard_uncertainty_m,
            reading_sensitivity_per_m,
            source,
        } => {
            out.push(2);
            for value in standard_uncertainty_m {
                push_wire_f64(out, *value);
            }
            for value in reading_sensitivity_per_m {
                push_wire_f64(out, *value);
            }
            push_wire_text(out, source)?;
        }
    }
    match &sensor.support.kind {
        ObservationSupportKind::Point { component } => {
            out.push(1);
            push_wire_u32(out, sensor.support.state_dimension)?;
            push_wire_u32(out, *component)?;
        }
        ObservationSupportKind::PatchAverage { terms } => {
            out.push(2);
            push_wire_u32(out, sensor.support.state_dimension)?;
            push_wire_u32(out, terms.len())?;
            for term in terms {
                push_wire_u32(out, term.component)?;
                push_wire_f64(out, term.weight);
            }
        }
    }
    match &sensor.mount.kind {
        SensorMountKind::DeclaredIdeal { source } => {
            out.push(1);
            push_wire_text(out, source)?;
        }
        SensorMountKind::Affine {
            gain,
            offset,
            source,
        } => {
            out.push(2);
            push_wire_f64(out, *gain);
            push_wire_f64(out, *offset);
            push_wire_text(out, source)?;
        }
    }
    match &sensor.dynamics.kind {
        SensorDynamicsKind::DeclaredInstantaneous { source } => {
            out.push(1);
            push_wire_text(out, source)?;
        }
        SensorDynamicsKind::FirstOrder {
            time_constant_s,
            source,
        } => {
            out.push(2);
            push_wire_f64(out, *time_constant_s);
            push_wire_text(out, source)?;
        }
    }
    match &sensor.calibration.kind {
        SensorCalibrationKind::Physical {
            certificate_ref,
            date,
            source,
            instrument_variance,
        } => {
            out.push(1);
            push_wire_text(out, certificate_ref)?;
            push_wire_text(out, date)?;
            push_wire_text(out, source)?;
            push_wire_f64(out, *instrument_variance);
        }
        SensorCalibrationKind::Virtual { definition_ref } => {
            out.push(2);
            push_wire_text(out, definition_ref)?;
        }
    }
    out.push(u8::from(sensor.placement_candidate));
    Ok(())
}

#[allow(clippy::too_many_lines)]
fn decode_sensor(
    reader: &mut SensorWireReader<'_>,
    row: usize,
) -> Result<ScenarioSensor, SensorExtensionError> {
    let invalid = |source| SensorExtensionError::InvalidSensor { row, source };
    let name = reader.read_text("sensor name")?;
    let kind = decode_sensor_kind(reader.read_u8("sensor kind")?)?;
    let entity = decode_entity_ref(reader)?;
    let mut local_position_m = [0.0; 3];
    for value in &mut local_position_m {
        *value = reader.read_f64("sensor local position")?;
    }
    let placement = match reader.read_u8("placement tag")? {
        1 => PlacementUncertainty::declared_exact(reader.read_text("placement source")?),
        2 => {
            let mut standard_uncertainty_m = [0.0; 3];
            for value in &mut standard_uncertainty_m {
                *value = reader.read_f64("placement uncertainty")?;
            }
            let mut reading_sensitivity_per_m = [0.0; 3];
            for value in &mut reading_sensitivity_per_m {
                *value = reader.read_f64("placement sensitivity")?;
            }
            PlacementUncertainty::axis_aligned(
                standard_uncertainty_m,
                reading_sensitivity_per_m,
                reader.read_text("placement source")?,
            )
        }
        observed => {
            return Err(SensorExtensionError::InvalidTag {
                field: "placement",
                observed,
            });
        }
    }
    .map_err(invalid)?;
    let support = match reader.read_u8("support tag")? {
        1 => ObservationSupport::point(
            reader.read_count("sensor state dimension", MAX_SENSOR_STATE_DIMENSION)?,
            reader.read_count("sensor component", MAX_SENSOR_STATE_DIMENSION)?,
        ),
        2 => {
            let state_dimension =
                reader.read_count("sensor state dimension", MAX_SENSOR_STATE_DIMENSION)?;
            let term_count = reader.read_count("sensor support terms", MAX_SENSOR_SUPPORT_TERMS)?;
            let mut terms = Vec::new();
            terms.try_reserve_exact(term_count).map_err(|_| {
                SensorExtensionError::AllocationRefused {
                    resource: "decoded sensor support terms",
                    requested: term_count,
                }
            })?;
            for _ in 0..term_count {
                terms.push(ObservationTerm::new(
                    reader.read_count("sensor component", MAX_SENSOR_STATE_DIMENSION)?,
                    reader.read_f64("sensor support weight")?,
                ));
            }
            ObservationSupport::patch_average(state_dimension, terms)
        }
        observed => {
            return Err(SensorExtensionError::InvalidTag {
                field: "support",
                observed,
            });
        }
    }
    .map_err(invalid)?;
    let mount = match reader.read_u8("mount tag")? {
        1 => SensorMount::declared_ideal(reader.read_text("mount source")?),
        2 => SensorMount::affine(
            reader.read_f64("mount gain")?,
            reader.read_f64("mount offset")?,
            reader.read_text("mount source")?,
        ),
        observed => {
            return Err(SensorExtensionError::InvalidTag {
                field: "mount",
                observed,
            });
        }
    }
    .map_err(invalid)?;
    let dynamics = match reader.read_u8("dynamics tag")? {
        1 => SensorDynamics::declared_instantaneous(reader.read_text("dynamics source")?),
        2 => SensorDynamics::first_order(
            reader.read_f64("dynamics time constant")?,
            reader.read_text("dynamics source")?,
        ),
        observed => {
            return Err(SensorExtensionError::InvalidTag {
                field: "dynamics",
                observed,
            });
        }
    }
    .map_err(invalid)?;
    let calibration = match reader.read_u8("calibration tag")? {
        1 => SensorCalibration::physical(
            reader.read_text("calibration certificate")?,
            reader.read_text("calibration date")?,
            reader.read_text("calibration source")?,
            reader.read_f64("instrument variance")?,
        ),
        2 => SensorCalibration::virtual_sensor(reader.read_text("virtual definition")?),
        observed => {
            return Err(SensorExtensionError::InvalidTag {
                field: "calibration",
                observed,
            });
        }
    }
    .map_err(invalid)?;
    let placement_candidate = match reader.read_u8("placement candidate")? {
        0 => false,
        1 => true,
        observed => {
            return Err(SensorExtensionError::InvalidTag {
                field: "placement candidate",
                observed,
            });
        }
    };
    let location = SensorLocation::new(entity, local_position_m, placement).map_err(invalid)?;
    ScenarioSensor::new(
        name,
        kind,
        location,
        support,
        mount,
        dynamics,
        calibration,
        placement_candidate,
    )
    .map_err(invalid)
}

fn decode_sensor_kind(tag: u8) -> Result<SensorKind, SensorExtensionError> {
    match tag {
        1 => Ok(SensorKind::Thermocouple),
        2 => Ok(SensorKind::Rtd),
        3 => Ok(SensorKind::FlowMeter),
        4 => Ok(SensorKind::PressureTap),
        5 => Ok(SensorKind::IrCameraRegion),
        observed => Err(SensorExtensionError::InvalidTag {
            field: "sensor kind",
            observed,
        }),
    }
}

fn encode_entity_ref(out: &mut Vec<u8>, reference: EntityRef) {
    out.push(entity_kind_tag(reference.target().kind()));
    out.extend_from_slice(reference.target().digest().as_bytes());
    out.push(expectation_tag(reference.expect()));
}

fn decode_entity_ref(reader: &mut SensorWireReader<'_>) -> Result<EntityRef, SensorExtensionError> {
    let kind = decode_entity_kind(reader.read_u8("entity kind")?)?;
    let digest = reader.read_hash("entity digest")?;
    let expectation = decode_expectation(reader.read_u8("entity expectation")?)?;
    Ok(EntityRef::new(
        EntityId::from_wire(kind, digest),
        expectation,
    ))
}

fn decode_entity_kind(tag: u8) -> Result<EntityKind, SensorExtensionError> {
    match tag {
        1 => Ok(EntityKind::Assembly),
        2 => Ok(EntityKind::Part),
        3 => Ok(EntityKind::Region),
        4 => Ok(EntityKind::Surface),
        5 => Ok(EntityKind::Interface),
        observed => Err(SensorExtensionError::InvalidTag {
            field: "entity kind",
            observed,
        }),
    }
}

fn expectation_tag(expectation: KindExpectation) -> u8 {
    match expectation {
        KindExpectation::Exact(EntityKind::Assembly) => 1,
        KindExpectation::Exact(EntityKind::Part) => 2,
        KindExpectation::Exact(EntityKind::Region) => 3,
        KindExpectation::Exact(EntityKind::Surface) => 4,
        KindExpectation::Exact(EntityKind::Interface) => 5,
        KindExpectation::Domain => 6,
        KindExpectation::Boundary => 7,
        KindExpectation::Any => 8,
    }
}

fn decode_expectation(tag: u8) -> Result<KindExpectation, SensorExtensionError> {
    match tag {
        1 => Ok(KindExpectation::Exact(EntityKind::Assembly)),
        2 => Ok(KindExpectation::Exact(EntityKind::Part)),
        3 => Ok(KindExpectation::Exact(EntityKind::Region)),
        4 => Ok(KindExpectation::Exact(EntityKind::Surface)),
        5 => Ok(KindExpectation::Exact(EntityKind::Interface)),
        6 => Ok(KindExpectation::Domain),
        7 => Ok(KindExpectation::Boundary),
        8 => Ok(KindExpectation::Any),
        observed => Err(SensorExtensionError::InvalidTag {
            field: "entity expectation",
            observed,
        }),
    }
}

fn push_wire_u32(out: &mut Vec<u8>, value: usize) -> Result<(), SensorExtensionError> {
    let value = u32::try_from(value).map_err(|_| SensorExtensionError::LengthOverflow)?;
    out.extend_from_slice(&value.to_le_bytes());
    Ok(())
}

fn push_wire_f64(out: &mut Vec<u8>, value: f64) {
    out.extend_from_slice(&canonicalize_zero(value).to_bits().to_le_bytes());
}

fn push_wire_text(out: &mut Vec<u8>, text: &str) -> Result<(), SensorExtensionError> {
    push_wire_u32(out, text.len())?;
    out.extend_from_slice(text.as_bytes());
    Ok(())
}

struct SensorWireReader<'a> {
    bytes: &'a [u8],
    position: usize,
}

impl<'a> SensorWireReader<'a> {
    const fn new(bytes: &'a [u8]) -> Self {
        Self { bytes, position: 0 }
    }

    const fn remaining(&self) -> usize {
        self.bytes.len().saturating_sub(self.position)
    }

    fn read_array<const N: usize>(
        &mut self,
        field: &'static str,
    ) -> Result<[u8; N], SensorExtensionError> {
        let end = self
            .position
            .checked_add(N)
            .ok_or(SensorExtensionError::Truncated { field })?;
        let source = self
            .bytes
            .get(self.position..end)
            .ok_or(SensorExtensionError::Truncated { field })?;
        let mut value = [0u8; N];
        value.copy_from_slice(source);
        self.position = end;
        Ok(value)
    }

    fn read_u8(&mut self, field: &'static str) -> Result<u8, SensorExtensionError> {
        Ok(self.read_array::<1>(field)?[0])
    }

    fn read_u16(&mut self, field: &'static str) -> Result<u16, SensorExtensionError> {
        Ok(u16::from_le_bytes(self.read_array(field)?))
    }

    fn read_u32(&mut self, field: &'static str) -> Result<u32, SensorExtensionError> {
        Ok(u32::from_le_bytes(self.read_array(field)?))
    }

    fn read_f64(&mut self, field: &'static str) -> Result<f64, SensorExtensionError> {
        Ok(f64::from_bits(u64::from_le_bytes(self.read_array(field)?)))
    }

    fn read_hash(&mut self, field: &'static str) -> Result<ContentHash, SensorExtensionError> {
        Ok(ContentHash(self.read_array(field)?))
    }

    fn read_count(
        &mut self,
        resource: &'static str,
        limit: usize,
    ) -> Result<usize, SensorExtensionError> {
        let requested = self.read_u32(resource)? as usize;
        if requested > limit {
            return Err(SensorExtensionError::LimitExceeded {
                resource,
                requested,
                limit,
            });
        }
        Ok(requested)
    }

    fn read_text(&mut self, field: &'static str) -> Result<String, SensorExtensionError> {
        let length = self.read_count(field, MAX_SENSOR_TEXT_BYTES)?;
        let end = self
            .position
            .checked_add(length)
            .ok_or(SensorExtensionError::Truncated { field })?;
        let source = self
            .bytes
            .get(self.position..end)
            .ok_or(SensorExtensionError::Truncated { field })?;
        let text = core::str::from_utf8(source)
            .map_err(|_| SensorExtensionError::InvalidUtf8 { field })?;
        let mut owned = String::new();
        owned
            .try_reserve_exact(length)
            .map_err(|_| SensorExtensionError::AllocationRefused {
                resource: field,
                requested: length,
            })?;
        owned.push_str(text);
        self.position = end;
        Ok(owned)
    }
}

/// Explicit limits for one catalog-checked sensor-set compilation.
///
/// `max_work` counts machine-independent checkpoints: pre-publication
/// boundaries, sensor identity derivations, exact-name comparisons, catalog
/// resolution/compilation boundaries, and ordered receipt rows.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SensorSetBudget {
    /// Maximum number of sensor declarations.
    pub max_sensors: usize,
    /// Maximum deterministic logical work units.
    pub max_work: u128,
}

impl Default for SensorSetBudget {
    fn default() -> Self {
        DEFAULT_SENSOR_SET_BUDGET
    }
}

/// Exact preflight shape for one catalog-checked sensor set.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SensorSetPlan {
    /// Sensor declarations in caller order.
    pub sensors: usize,
    /// Exact pairwise name comparisons needed to prove uniqueness.
    pub duplicate_comparisons: u128,
    /// Deterministic logical work units.
    pub planned_work: u128,
}

/// Preflight one sensor set without allocating or consulting the catalog.
///
/// # Errors
///
/// Refuses a sensor-count limit, checked work arithmetic overflow, or a work
/// plan larger than `budget.max_work`.
pub fn plan_sensor_set(
    sensors: &[ScenarioSensor],
    budget: SensorSetBudget,
) -> Result<SensorSetPlan, SensorSetError> {
    if sensors.len() > budget.max_sensors {
        return Err(SensorSetError::LimitExceeded {
            resource: "sensor declarations",
            requested: sensors.len(),
            limit: budget.max_sensors,
        });
    }
    let count = sensors.len() as u128;
    let duplicate_comparisons = count
        .checked_mul(count.saturating_sub(1))
        .and_then(|value| value.checked_div(2))
        .ok_or(SensorSetError::WorkPlanOverflow {
            phase: "duplicate-name comparisons",
        })?;
    let sensor_work = count
        .checked_mul(4)
        .ok_or(SensorSetError::WorkPlanOverflow {
            phase: "per-sensor work",
        })?;
    let planned_work = 2u128
        .checked_add(sensor_work)
        .and_then(|value| value.checked_add(duplicate_comparisons))
        .ok_or(SensorSetError::WorkPlanOverflow {
            phase: "total sensor-set work",
        })?;
    if planned_work > budget.max_work {
        return Err(SensorSetError::WorkExceeded {
            requested: planned_work,
            limit: budget.max_work,
        });
    }
    Ok(SensorSetPlan {
        sensors: sensors.len(),
        duplicate_comparisons,
        planned_work,
    })
}

/// One catalog resolution and compiled operator in stable caller row order.
#[derive(Debug, Clone, PartialEq)]
pub struct CompiledSensorBinding {
    row: usize,
    requested_entity: EntityId,
    current_entity: EntityId,
    supersession_hops: usize,
    evidence_tier: EvidenceTier,
    operator: CompiledSensorOperator,
}

impl CompiledSensorBinding {
    /// Zero-based caller row retained by the set receipt.
    #[must_use]
    pub const fn row(&self) -> usize {
        self.row
    }

    /// Entity identity authored by the sensor declaration.
    #[must_use]
    pub const fn requested_entity(&self) -> EntityId {
        self.requested_entity
    }

    /// Current entity identity selected through the catalog.
    #[must_use]
    pub const fn current_entity(&self) -> EntityId {
        self.current_entity
    }

    /// Number of catalog supersession hops followed.
    #[must_use]
    pub const fn supersession_hops(&self) -> usize {
        self.supersession_hops
    }

    /// Weakest evidence tier on the catalog resolution path.
    #[must_use]
    pub const fn evidence_tier(&self) -> EvidenceTier {
        self.evidence_tier
    }

    /// Operator whose entity accessor names the resolved current entity.
    #[must_use]
    pub const fn operator(&self) -> &CompiledSensorOperator {
        &self.operator
    }
}

/// All-or-nothing catalog-checked compilation of an ordered sensor set.
#[derive(Debug, Clone, PartialEq)]
pub struct CompiledSensorSet {
    identity: ContentHash,
    catalog_receipt_root: ContentHash,
    bindings: Vec<CompiledSensorBinding>,
}

impl CompiledSensorSet {
    /// Domain-separated identity of the exact catalog snapshot and ordered
    /// compiled bindings.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    /// Exact catalog receipt root bound into [`Self::identity`].
    #[must_use]
    pub const fn catalog_receipt_root(&self) -> ContentHash {
        self.catalog_receipt_root
    }

    /// Compiled bindings in caller row order.
    #[must_use]
    pub fn bindings(&self) -> &[CompiledSensorBinding] {
        &self.bindings
    }
}

/// Compile an ordered sensor set under [`DEFAULT_SENSOR_SET_BUDGET`].
///
/// # Errors
///
/// Returns a typed resource, cancellation, duplicate-name, catalog-resolution,
/// resolved-kind, or operator-compilation refusal. No partial set is returned.
pub fn compile_sensor_set(
    sensors: &[ScenarioSensor],
    catalog: &EntityCatalog,
    cx: &Cx<'_>,
) -> Result<CompiledSensorSet, SensorSetError> {
    compile_sensor_set_with_budget(sensors, catalog, SensorSetBudget::default(), cx)
}

/// Compile an ordered sensor set under an explicit budget.
///
/// The exact catalog receipt root is conservatively semantic: even an
/// unrelated catalog receipt changes the set identity. The root must be
/// externally pinned if the caller needs tail-truncation detection.
///
/// Cancellation is checked before preflight and at every planned set boundary.
/// Individual bounded entity-resolution walks and dense operator compilation
/// are not internally interruptible in schema v1.
///
/// # Errors
///
/// Returns a typed resource, cancellation, duplicate-name, catalog-resolution,
/// resolved-kind, or operator-compilation refusal. No partial set is returned.
pub fn compile_sensor_set_with_budget(
    sensors: &[ScenarioSensor],
    catalog: &EntityCatalog,
    budget: SensorSetBudget,
    cx: &Cx<'_>,
) -> Result<CompiledSensorSet, SensorSetError> {
    cx.checkpoint().map_err(|_| SensorSetError::Cancelled {
        phase: "initial",
        completed: 0,
        planned: 0,
    })?;
    let plan = plan_sensor_set(sensors, budget)?;
    let mut completed = 0u128;
    sensor_set_checkpoint(cx, "post-preflight", &mut completed, plan.planned_work)?;

    let identities = sensor_identities(sensors, cx, &mut completed, plan.planned_work)?;
    check_unique_sensor_names(sensors, cx, &mut completed, plan.planned_work)?;

    let mut bindings = Vec::new();
    bindings
        .try_reserve_exact(plan.sensors)
        .map_err(|_| SensorSetError::AllocationRefused {
            resource: "compiled sensor bindings",
            requested: plan.sensors,
        })?;
    for (row, (sensor, sensor_identity)) in sensors.iter().zip(identities).enumerate() {
        sensor_set_checkpoint(
            cx,
            "entity-resolution boundary",
            &mut completed,
            plan.planned_work,
        )?;
        let resolution = catalog
            .resolve(sensor.location().entity())
            .map_err(|fault| SensorSetError::Resolution { row, fault })?;
        let actual = resolution.current().kind();
        if !sensor.kind().admits_entity(actual) {
            return Err(SensorSetError::ResolvedEntityKind {
                row,
                current: resolution.current(),
                actual,
                sensor_kind: sensor.kind(),
            });
        }
        let mut operator = sensor
            .compile()
            .map_err(|source| SensorSetError::Compilation { row, source })?;
        operator.entity = resolution.current();
        debug_assert_eq!(operator.sensor_identity(), sensor_identity);
        sensor_set_checkpoint(
            cx,
            "operator compilation",
            &mut completed,
            plan.planned_work,
        )?;
        bindings.push(CompiledSensorBinding {
            row,
            requested_entity: resolution.requested(),
            current_entity: resolution.current(),
            supersession_hops: resolution.hops(),
            evidence_tier: resolution.tier(),
            operator,
        });
    }

    let catalog_receipt_root = catalog.receipt_root();
    let identity = sensor_set_identity(
        catalog_receipt_root,
        &bindings,
        cx,
        &mut completed,
        plan.planned_work,
    )?;
    sensor_set_checkpoint(cx, "pre-publication", &mut completed, plan.planned_work)?;
    debug_assert_eq!(completed, plan.planned_work);
    Ok(CompiledSensorSet {
        identity,
        catalog_receipt_root,
        bindings,
    })
}

/// Dense, owner-neutral observation operator compiled from a scenario sensor.
#[derive(Debug, Clone, PartialEq)]
pub struct CompiledSensorOperator {
    sensor_identity: ContentHash,
    name: String,
    quantity: SensorQuantity,
    entity: EntityId,
    local_position_m: [f64; 3],
    operator: Vec<f64>,
    offset: f64,
    placement_variance: f64,
    instrument_variance: Option<f64>,
    virtual_sensor: bool,
    placement_candidate: bool,
}

impl CompiledSensorOperator {
    /// Scenario-sensor content identity.
    #[must_use]
    pub const fn sensor_identity(&self) -> ContentHash {
        self.sensor_identity
    }

    /// Stable instrument identity accepted by bounded leaf-identity consumers.
    #[must_use]
    pub fn instrument_identity(&self) -> String {
        self.sensor_identity.to_hex()
    }

    /// Declared sensor name.
    #[must_use]
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Measured quantity.
    #[must_use]
    pub const fn quantity(&self) -> SensorQuantity {
        self.quantity
    }

    /// Bound entity identity.
    #[must_use]
    pub const fn entity(&self) -> EntityId {
        self.entity
    }

    /// Entity-local coordinates in metres.
    #[must_use]
    pub const fn local_position_m(&self) -> [f64; 3] {
        self.local_position_m
    }

    /// Dense row after applying mount gain.
    #[must_use]
    pub fn operator(&self) -> &[f64] {
        &self.operator
    }

    /// Affine mount offset in output units.
    #[must_use]
    pub const fn offset(&self) -> f64 {
        self.offset
    }

    /// Placement contribution to reading variance.
    #[must_use]
    pub const fn placement_variance(&self) -> f64 {
        self.placement_variance
    }

    /// Calibration/instrument contribution to reading variance.
    #[must_use]
    pub const fn instrument_variance(&self) -> Option<f64> {
        self.instrument_variance
    }

    /// Whether this operator belongs to a virtual sensor.
    #[must_use]
    pub const fn is_virtual(&self) -> bool {
        self.virtual_sensor
    }

    /// Whether this operator is a placement candidate.
    #[must_use]
    pub const fn is_placement_candidate(&self) -> bool {
        self.placement_candidate
    }

    /// Evaluate the same compiled affine operator used for measurement
    /// handoff and corpus comparison.
    ///
    /// # Errors
    ///
    /// Refuses dimension mismatch, non-finite state values, or finite-input
    /// overflow.
    pub fn predict(&self, state: &[f64]) -> Result<f64, SensorError> {
        if state.len() != self.operator.len() {
            return Err(SensorError::StateDimensionMismatch {
                expected: self.operator.len(),
                actual: state.len(),
            });
        }
        let mut value = self.offset;
        for (weight, component) in self.operator.iter().zip(state) {
            if !component.is_finite() {
                return Err(SensorError::NonFiniteField {
                    field: "sensor state",
                });
            }
            value = weight.mul_add(*component, value);
            if !value.is_finite() {
                return Err(SensorError::NonFiniteComputation {
                    operation: "sensor prediction",
                });
            }
        }
        Ok(value)
    }

    /// Compare one measured value against the compiled prediction.
    ///
    /// The residual uses exactly the same operator returned for assimilation,
    /// avoiding a second probe definition.
    ///
    /// # Errors
    ///
    /// Refuses a non-finite measurement or any prediction refusal.
    pub fn compare(&self, measured: f64, state: &[f64]) -> Result<SensorComparison, SensorError> {
        if !measured.is_finite() {
            return Err(SensorError::NonFiniteField {
                field: "measured reading",
            });
        }
        let predicted = self.predict(state)?;
        let residual = measured - predicted;
        if !residual.is_finite() {
            return Err(SensorError::NonFiniteComputation {
                operation: "sensor comparison residual",
            });
        }
        Ok(SensorComparison {
            predicted,
            measured,
            residual: canonicalize_zero(residual),
        })
    }

    /// Produce the exact parts for a linear-Gaussian observation consumer.
    ///
    /// The affine mount equation is converted from
    /// `measured = gain * H * state + offset + noise` to the standard
    /// `adjusted_value = operator * state + noise` form by subtracting the
    /// retained offset. The returned operator is cloned so its owner cannot
    /// mutate this compiled artifact.
    ///
    /// # Errors
    ///
    /// Virtual sensors refuse because they carry no physical measurement-noise
    /// authority. Physical sensors refuse non-finite measurement, variance
    /// overflow, or allocation failure.
    pub fn observation_parts(&self, measured: f64) -> Result<SensorObservationParts, SensorError> {
        if self.virtual_sensor {
            return Err(SensorError::VirtualSensorHasNoNoiseAuthority);
        }
        if !measured.is_finite() {
            return Err(SensorError::NonFiniteField {
                field: "measured reading",
            });
        }
        let adjusted_value = measured - self.offset;
        if !adjusted_value.is_finite() {
            return Err(SensorError::NonFiniteComputation {
                operation: "mount offset removal",
            });
        }
        let instrument_variance = self
            .instrument_variance
            .ok_or(SensorError::VirtualSensorHasNoNoiseAuthority)?;
        let noise_variance = instrument_variance + self.placement_variance;
        if !noise_variance.is_finite() || noise_variance <= 0.0 {
            return Err(SensorError::NonFiniteComputation {
                operation: "sensor noise aggregation",
            });
        }
        let mut operator = Vec::new();
        operator
            .try_reserve_exact(self.operator.len())
            .map_err(|_| SensorError::AllocationRefused {
                resource: "observation operator handoff",
            })?;
        operator.extend_from_slice(&self.operator);
        Ok(SensorObservationParts {
            operator,
            adjusted_value: canonicalize_zero(adjusted_value),
            noise_variance,
            instrument_identity: self.instrument_identity(),
        })
    }
}

/// Exact owner-neutral parts of one linear-Gaussian observation.
#[derive(Debug, Clone, PartialEq)]
pub struct SensorObservationParts {
    operator: Vec<f64>,
    adjusted_value: f64,
    noise_variance: f64,
    instrument_identity: String,
}

impl SensorObservationParts {
    /// Dense observation row.
    #[must_use]
    pub fn operator(&self) -> &[f64] {
        &self.operator
    }

    /// Measured value after removing the affine mount offset.
    #[must_use]
    pub const fn adjusted_value(&self) -> f64 {
        self.adjusted_value
    }

    /// Instrument plus propagated placement variance.
    #[must_use]
    pub const fn noise_variance(&self) -> f64 {
        self.noise_variance
    }

    /// Stable instrument identity.
    #[must_use]
    pub fn instrument_identity(&self) -> &str {
        &self.instrument_identity
    }
}

/// One predicted-versus-measured comparison using a compiled sensor.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SensorComparison {
    predicted: f64,
    measured: f64,
    residual: f64,
}

impl SensorComparison {
    /// Simulated sensor reading.
    #[must_use]
    pub const fn predicted(self) -> f64 {
        self.predicted
    }

    /// Supplied measured reading.
    #[must_use]
    pub const fn measured(self) -> f64 {
        self.measured
    }

    /// `measured - predicted`.
    #[must_use]
    pub const fn residual(self) -> f64 {
        self.residual
    }
}

/// Typed scenario-sensor refusal.
#[derive(Debug, Clone, PartialEq)]
pub enum SensorError {
    /// A required retained string is empty.
    EmptyText {
        /// Field name.
        field: &'static str,
    },
    /// A retained string exceeds the public bound.
    TextTooLong {
        /// Field name.
        field: &'static str,
        /// Actual bytes.
        actual: usize,
        /// Maximum bytes.
        limit: usize,
    },
    /// A public numeric field is non-finite.
    NonFiniteField {
        /// Field name.
        field: &'static str,
    },
    /// The reference expectation disagrees with the target identity kind.
    EntityExpectationMismatch {
        /// Target identity.
        target: EntityId,
        /// Declared expectation.
        expected: KindExpectation,
    },
    /// A sensor family cannot bind to this entity kind.
    UnsupportedEntityKind {
        /// Sensor family.
        kind: SensorKind,
        /// Entity kind.
        entity_kind: EntityKind,
    },
    /// Placement uncertainty is malformed.
    InvalidPlacementUncertainty {
        /// Diagnosis.
        what: &'static str,
    },
    /// State dimension is outside the v1 bound.
    InvalidStateDimension {
        /// Requested dimension.
        actual: usize,
        /// Maximum admitted dimension.
        limit: usize,
    },
    /// One sparse component is outside the state.
    ComponentOutOfRange {
        /// Component index.
        component: usize,
        /// State dimension.
        state_dimension: usize,
    },
    /// Point/patch support has invalid cardinality.
    SupportSize {
        /// Actual terms.
        actual: usize,
        /// Maximum terms.
        limit: usize,
    },
    /// A support weight is non-finite or non-positive.
    InvalidSupportWeight {
        /// Component carrying the weight.
        component: usize,
    },
    /// A sparse support repeats a component.
    DuplicateComponent {
        /// Repeated component.
        component: usize,
    },
    /// Patch weights are not a declared average.
    PatchWeightsDoNotSumToOne {
        /// Observed sum.
        actual: f64,
    },
    /// Mount model is malformed.
    InvalidMount {
        /// Diagnosis.
        what: &'static str,
    },
    /// Dynamics model is malformed.
    InvalidDynamics {
        /// Diagnosis.
        what: &'static str,
    },
    /// Calibration date is not a valid `YYYY-MM-DD` date.
    InvalidCalibrationDate,
    /// Physical instrument variance is not finite and strictly positive.
    InvalidInstrumentVariance,
    /// Dense state length differs from the compiled operator.
    StateDimensionMismatch {
        /// Expected length.
        expected: usize,
        /// Actual length.
        actual: usize,
    },
    /// A finite-input computation overflowed.
    NonFiniteComputation {
        /// Operation.
        operation: &'static str,
    },
    /// A virtual sensor cannot become a physical observation without a noise
    /// authority.
    VirtualSensorHasNoNoiseAuthority,
    /// Fallible allocation refused before publication.
    AllocationRefused {
        /// Resource being allocated.
        resource: &'static str,
    },
}

/// Typed refusal from catalog-checked sensor-set compilation.
#[derive(Debug, Clone, PartialEq)]
pub enum SensorSetError {
    /// The sensor collection exceeds its explicit cap.
    LimitExceeded {
        /// Stable resource name.
        resource: &'static str,
        /// Requested elements.
        requested: usize,
        /// Admitted elements.
        limit: usize,
    },
    /// Checked work-plan arithmetic overflowed.
    WorkPlanOverflow {
        /// Phase whose arithmetic overflowed.
        phase: &'static str,
    },
    /// The exact work plan exceeds its explicit cap.
    WorkExceeded {
        /// Requested logical units.
        requested: u128,
        /// Admitted logical units.
        limit: u128,
    },
    /// A preflighted allocation was refused.
    AllocationRefused {
        /// Stable resource name.
        resource: &'static str,
        /// Requested elements.
        requested: usize,
    },
    /// Two caller rows use the same exact sensor name.
    DuplicateName {
        /// First row carrying the name.
        first: usize,
        /// Later row carrying the same name.
        second: usize,
    },
    /// The authored entity reference did not resolve.
    Resolution {
        /// Sensor row.
        row: usize,
        /// Structured catalog refusal.
        fault: ResolutionFault,
    },
    /// A resolved entity is not admitted by the sensor family.
    ResolvedEntityKind {
        /// Sensor row.
        row: usize,
        /// Resolved current identity.
        current: EntityId,
        /// Resolved current kind.
        actual: EntityKind,
        /// Sensor family imposing the narrower contract.
        sensor_kind: SensorKind,
    },
    /// One admitted declaration did not compile.
    Compilation {
        /// Sensor row.
        row: usize,
        /// Operator refusal.
        source: SensorError,
    },
    /// Cancellation was observed before publication.
    Cancelled {
        /// Stable checkpoint phase.
        phase: &'static str,
        /// Fully completed logical work units.
        completed: u128,
        /// Exact preflighted units, or zero before preflight.
        planned: u128,
    },
}

impl fmt::Display for SensorSetError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::LimitExceeded {
                resource,
                requested,
                limit,
            } => write!(
                formatter,
                "sensor-set {resource} request {requested} exceeds limit {limit}"
            ),
            Self::WorkPlanOverflow { phase } => {
                write!(formatter, "sensor-set work plan overflowed during {phase}")
            }
            Self::WorkExceeded { requested, limit } => write!(
                formatter,
                "sensor-set work request {requested} exceeds limit {limit}"
            ),
            Self::AllocationRefused {
                resource,
                requested,
            } => write!(
                formatter,
                "sensor-set allocation for {requested} {resource} elements was refused"
            ),
            Self::DuplicateName { first, second } => write!(
                formatter,
                "sensor rows {first} and {second} repeat the same exact name"
            ),
            Self::Resolution { row, fault } => {
                write!(
                    formatter,
                    "sensor row {row} entity resolution refused with {}: {fault}",
                    fault.code()
                )
            }
            Self::ResolvedEntityKind {
                row,
                current,
                actual,
                sensor_kind,
            } => write!(
                formatter,
                "sensor row {row} resolves to {current} ({actual}), which is not admitted by {}",
                sensor_kind.label()
            ),
            Self::Compilation { row, source } => {
                write!(formatter, "sensor row {row} did not compile: {source}")
            }
            Self::Cancelled {
                phase,
                completed,
                planned,
            } => write!(
                formatter,
                "sensor-set compilation cancelled during {phase} after {completed}/{planned} work units"
            ),
        }
    }
}

impl core::error::Error for SensorSetError {}

impl fmt::Display for SensorError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyText { field } => write!(formatter, "{field} must not be empty"),
            Self::TextTooLong {
                field,
                actual,
                limit,
            } => write!(
                formatter,
                "{field} is {actual} bytes; the admitted maximum is {limit}"
            ),
            Self::NonFiniteField { field } => write!(formatter, "{field} must be finite"),
            Self::EntityExpectationMismatch { target, expected } => write!(
                formatter,
                "sensor entity reference {target} declares expectation {expected}, which does not admit the target identity kind"
            ),
            Self::UnsupportedEntityKind { kind, entity_kind } => write!(
                formatter,
                "{} sensors cannot bind to {entity_kind} entities in schema v{SENSOR_SCHEMA_VERSION}",
                kind.label()
            ),
            Self::InvalidPlacementUncertainty { what }
            | Self::InvalidMount { what }
            | Self::InvalidDynamics { what } => formatter.write_str(what),
            Self::InvalidStateDimension { actual, limit } => write!(
                formatter,
                "sensor state dimension {actual} is outside 1..={limit}"
            ),
            Self::ComponentOutOfRange {
                component,
                state_dimension,
            } => write!(
                formatter,
                "sensor component {component} is outside state dimension {state_dimension}"
            ),
            Self::SupportSize { actual, limit } => write!(
                formatter,
                "sensor support has {actual} terms; the admitted range is 1..={limit}"
            ),
            Self::InvalidSupportWeight { component } => write!(
                formatter,
                "sensor support weight for component {component} must be finite and positive"
            ),
            Self::DuplicateComponent { component } => {
                write!(formatter, "sensor support repeats component {component}")
            }
            Self::PatchWeightsDoNotSumToOne { actual } => write!(
                formatter,
                "patch-average weights sum to {actual:e}, not one within {PATCH_WEIGHT_TOLERANCE:e}"
            ),
            Self::InvalidCalibrationDate => {
                formatter.write_str("calibration date must be a valid YYYY-MM-DD date")
            }
            Self::InvalidInstrumentVariance => formatter.write_str(
                "physical instrument variance must be finite and strictly positive",
            ),
            Self::StateDimensionMismatch { expected, actual } => write!(
                formatter,
                "sensor operator expects state dimension {expected}, got {actual}"
            ),
            Self::NonFiniteComputation { operation } => {
                write!(formatter, "{operation} produced a non-finite result")
            }
            Self::VirtualSensorHasNoNoiseAuthority => formatter.write_str(
                "virtual sensor has no physical measurement-noise authority; use predict/compare for simulation or attach a physical calibration",
            ),
            Self::AllocationRefused { resource } => {
                write!(formatter, "{resource} allocation was refused before publication")
            }
        }
    }
}

impl core::error::Error for SensorError {}

fn sensor_identities(
    sensors: &[ScenarioSensor],
    cx: &Cx<'_>,
    completed: &mut u128,
    planned: u128,
) -> Result<Vec<ContentHash>, SensorSetError> {
    let mut identities = Vec::new();
    identities
        .try_reserve_exact(sensors.len())
        .map_err(|_| SensorSetError::AllocationRefused {
            resource: "sensor identities",
            requested: sensors.len(),
        })?;
    for sensor in sensors {
        identities.push(sensor.identity());
        sensor_set_checkpoint(cx, "sensor identity", completed, planned)?;
    }
    Ok(identities)
}

fn check_unique_sensor_names(
    sensors: &[ScenarioSensor],
    cx: &Cx<'_>,
    completed: &mut u128,
    planned: u128,
) -> Result<(), SensorSetError> {
    for second in 1..sensors.len() {
        for first in 0..second {
            let duplicate = sensors[first].name() == sensors[second].name();
            sensor_set_checkpoint(cx, "duplicate-name comparison", completed, planned)?;
            if duplicate {
                return Err(SensorSetError::DuplicateName { first, second });
            }
        }
    }
    Ok(())
}

fn sensor_set_identity(
    catalog_receipt_root: ContentHash,
    bindings: &[CompiledSensorBinding],
    cx: &Cx<'_>,
    completed: &mut u128,
    planned: u128,
) -> Result<ContentHash, SensorSetError> {
    let mut hasher = DomainHasher::new(SENSOR_SET_IDENTITY_DOMAIN);
    absorb_u64(&mut hasher, u64::from(SENSOR_SET_SCHEMA_VERSION));
    hasher.update(catalog_receipt_root.as_bytes());
    absorb_u64(&mut hasher, bindings.len() as u64);
    for binding in bindings {
        absorb_u64(&mut hasher, binding.row as u64);
        hasher.update(binding.operator.sensor_identity().as_bytes());
        absorb_entity_id(&mut hasher, binding.requested_entity);
        absorb_entity_id(&mut hasher, binding.current_entity);
        absorb_u64(&mut hasher, binding.supersession_hops as u64);
        absorb_bytes(&mut hasher, binding.evidence_tier.label().as_bytes());
        sensor_set_checkpoint(cx, "sensor-set identity row", completed, planned)?;
    }
    Ok(hasher.finalize())
}

fn sensor_set_checkpoint(
    cx: &Cx<'_>,
    phase: &'static str,
    completed: &mut u128,
    planned: u128,
) -> Result<(), SensorSetError> {
    cx.checkpoint().map_err(|_| SensorSetError::Cancelled {
        phase,
        completed: *completed,
        planned,
    })?;
    *completed = completed
        .checked_add(1)
        .ok_or(SensorSetError::WorkPlanOverflow {
            phase: "completed-work accounting",
        })?;
    Ok(())
}

fn validate_text(field: &'static str, value: &str) -> Result<(), SensorError> {
    if value.is_empty() {
        return Err(SensorError::EmptyText { field });
    }
    if value.len() > MAX_SENSOR_TEXT_BYTES {
        return Err(SensorError::TextTooLong {
            field,
            actual: value.len(),
            limit: MAX_SENSOR_TEXT_BYTES,
        });
    }
    Ok(())
}

fn validate_state_dimension(state_dimension: usize) -> Result<(), SensorError> {
    if state_dimension == 0 || state_dimension > MAX_SENSOR_STATE_DIMENSION {
        return Err(SensorError::InvalidStateDimension {
            actual: state_dimension,
            limit: MAX_SENSOR_STATE_DIMENSION,
        });
    }
    Ok(())
}

fn validate_date(date: &str) -> Result<(), SensorError> {
    let bytes = date.as_bytes();
    if bytes.len() != 10 || bytes[4] != b'-' || bytes[7] != b'-' {
        return Err(SensorError::InvalidCalibrationDate);
    }
    let year = parse_decimal(&bytes[0..4]).ok_or(SensorError::InvalidCalibrationDate)?;
    let month = parse_decimal(&bytes[5..7]).ok_or(SensorError::InvalidCalibrationDate)?;
    let day = parse_decimal(&bytes[8..10]).ok_or(SensorError::InvalidCalibrationDate)?;
    let leap = year % 4 == 0 && (year % 100 != 0 || year % 400 == 0);
    let days = match month {
        1 | 3 | 5 | 7 | 8 | 10 | 12 => 31,
        4 | 6 | 9 | 11 => 30,
        2 if leap => 29,
        2 => 28,
        _ => return Err(SensorError::InvalidCalibrationDate),
    };
    if day == 0 || day > days {
        return Err(SensorError::InvalidCalibrationDate);
    }
    Ok(())
}

fn parse_decimal(bytes: &[u8]) -> Option<u32> {
    let mut value = 0u32;
    for byte in bytes {
        if !byte.is_ascii_digit() {
            return None;
        }
        value = value
            .checked_mul(10)?
            .checked_add(u32::from(*byte - b'0'))?;
    }
    Some(value)
}

fn canonicalize_zero(value: f64) -> f64 {
    if value == 0.0 { 0.0 } else { value }
}

fn canonicalize_zero3(mut value: [f64; 3]) -> [f64; 3] {
    for component in &mut value {
        *component = canonicalize_zero(*component);
    }
    value
}

fn absorb_u64(hasher: &mut DomainHasher, value: u64) {
    hasher.update(&value.to_le_bytes());
}

fn absorb_bytes(hasher: &mut DomainHasher, bytes: &[u8]) {
    absorb_u64(hasher, bytes.len() as u64);
    hasher.update(bytes);
}

fn absorb_f64s(hasher: &mut DomainHasher, values: &[f64]) {
    absorb_u64(hasher, values.len() as u64);
    for value in values {
        hasher.update(&canonicalize_zero(*value).to_bits().to_le_bytes());
    }
}

fn absorb_entity_id(hasher: &mut DomainHasher, entity: EntityId) {
    hasher.update(&[entity_kind_tag(entity.kind())]);
    hasher.update(entity.digest().as_bytes());
}

fn absorb_entity_ref(hasher: &mut DomainHasher, entity: EntityRef) {
    absorb_entity_id(hasher, entity.target());
    hasher.update(&[match entity.expect() {
        KindExpectation::Exact(EntityKind::Assembly) => 1,
        KindExpectation::Exact(EntityKind::Part) => 2,
        KindExpectation::Exact(EntityKind::Region) => 3,
        KindExpectation::Exact(EntityKind::Surface) => 4,
        KindExpectation::Exact(EntityKind::Interface) => 5,
        KindExpectation::Domain => 6,
        KindExpectation::Boundary => 7,
        KindExpectation::Any => 8,
    }]);
}

const fn entity_kind_tag(kind: EntityKind) -> u8 {
    match kind {
        EntityKind::Assembly => 1,
        EntityKind::Part => 2,
        EntityKind::Region => 3,
        EntityKind::Surface => 4,
        EntityKind::Interface => 5,
    }
}

fn absorb_placement(hasher: &mut DomainHasher, placement: &PlacementUncertainty) {
    match &placement.kind {
        PlacementUncertaintyKind::DeclaredExact { source } => {
            hasher.update(&[1]);
            absorb_bytes(hasher, source.as_bytes());
        }
        PlacementUncertaintyKind::AxisAligned {
            standard_uncertainty_m,
            reading_sensitivity_per_m,
            source,
        } => {
            hasher.update(&[2]);
            absorb_f64s(hasher, standard_uncertainty_m);
            absorb_f64s(hasher, reading_sensitivity_per_m);
            absorb_bytes(hasher, source.as_bytes());
        }
    }
}

fn absorb_support(hasher: &mut DomainHasher, support: &ObservationSupport) {
    match &support.kind {
        ObservationSupportKind::Point { component } => {
            hasher.update(&[1]);
            absorb_u64(hasher, support.state_dimension as u64);
            absorb_u64(hasher, *component as u64);
        }
        ObservationSupportKind::PatchAverage { terms } => {
            hasher.update(&[2]);
            absorb_u64(hasher, support.state_dimension as u64);
            absorb_u64(hasher, terms.len() as u64);
            for term in terms {
                absorb_u64(hasher, term.component as u64);
                hasher.update(&canonicalize_zero(term.weight).to_bits().to_le_bytes());
            }
        }
    }
}

fn absorb_mount(hasher: &mut DomainHasher, mount: &SensorMount) {
    match &mount.kind {
        SensorMountKind::DeclaredIdeal { source } => {
            hasher.update(&[1]);
            absorb_bytes(hasher, source.as_bytes());
        }
        SensorMountKind::Affine {
            gain,
            offset,
            source,
        } => {
            hasher.update(&[2]);
            absorb_f64s(hasher, &[*gain, *offset]);
            absorb_bytes(hasher, source.as_bytes());
        }
    }
}

fn absorb_dynamics(hasher: &mut DomainHasher, dynamics: &SensorDynamics) {
    match &dynamics.kind {
        SensorDynamicsKind::DeclaredInstantaneous { source } => {
            hasher.update(&[1]);
            absorb_bytes(hasher, source.as_bytes());
        }
        SensorDynamicsKind::FirstOrder {
            time_constant_s,
            source,
        } => {
            hasher.update(&[2]);
            absorb_f64s(hasher, &[*time_constant_s]);
            absorb_bytes(hasher, source.as_bytes());
        }
    }
}

fn absorb_calibration(hasher: &mut DomainHasher, calibration: &SensorCalibration) {
    match &calibration.kind {
        SensorCalibrationKind::Physical {
            certificate_ref,
            date,
            source,
            instrument_variance,
        } => {
            hasher.update(&[1]);
            absorb_bytes(hasher, certificate_ref.as_bytes());
            absorb_bytes(hasher, date.as_bytes());
            absorb_bytes(hasher, source.as_bytes());
            absorb_f64s(hasher, &[*instrument_variance]);
        }
        SensorCalibrationKind::Virtual { definition_ref } => {
            hasher.update(&[2]);
            absorb_bytes(hasher, definition_ref.as_bytes());
        }
    }
}
