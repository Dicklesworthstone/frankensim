//! Level-E instrumented rig: the specification as DATA, and a fail-closed
//! ingest that turns one measured run into an admitted corpus record.
//!
//! Owned hardware is the only source that removes the two chronic blockers of
//! published validation data — missing raw records and unlicensable detail —
//! so the corpus needs a path for it. Physical procurement has calendar
//! constraints code does not, so the software half lands first and is
//! exercised against synthetic runs; nothing here waits on a rig existing.
//!
//! # What this module refuses
//!
//! Ingest is a gate, not a converter. It refuses a run whose channels do not
//! match the declared instrument set, whose sensors lack calibration
//! certificates, or whose measured energy does not balance. Non-blind
//! refusals name the sensor or the numbers, because the point of ingesting
//! through a gate is to learn what is wrong with the data before it becomes
//! evidence. Blind refusals preserve the same typed cause while redacting
//! measured values and interval endpoints.
//!
//! # The energy-balance gate
//!
//! A thermal rig run has a physical closure condition: the electrical power
//! into the heater must leave as enthalpy in the coolant stream.
//!
//! ```text
//!   P_in  ≈  rho · c_p · Q · (T_out − T_in)
//! ```
//!
//! This is a check on the DATA, not on any model — no solver is involved, and
//! it catches the mundane failures that quietly poison a corpus: a
//! mis-scaled flow channel, swapped inlet/outlet thermocouples, an unlogged
//! secondary heat path, a decimal error in declared power.
//!
//! It is adjudicated against a finite interval enclosure built from the
//! instruments' OWN bounded half-widths. Product cross terms are retained and
//! every arithmetic endpoint is widened outward. That has a deliberate
//! consequence: a balance channel declaring
//! [`MeasurementUncertainty::Unstated`] or only a covariance cannot pass,
//! because neither supplies a finite support band. That is a feature. A
//! Level-E dataset exists precisely to carry metrology, and admitting an
//! unquantified run would produce exactly the unfalsifiable record this corpus
//! is meant to exclude.
//!
//! # Blind holdout
//!
//! A [`DatasetPartition::BlindHoldout`] run is SEALED after ingest: the gate
//! still runs, so a blind run is known to be physically sane, but its
//! measured values do not surface through ordinary verdict/readings accessors.
//! The exact raw and metrology bytes are nevertheless retained through an
//! explicit persistence surface. Sealing after checking is the ordering that
//! gives both properties; checking after sealing would give neither.

use crate::corpus::{
    Availability, DatasetPartition, MAX_CORPUS_TEXT_BYTES, MeasurementUncertainty,
    SENSOR_RECORD_SCHEMA_VERSION, SensorRecord, decode_sensor_record, encode_sensor_record,
    valid_date, valid_slug, validate_sensor_record,
};
use fs_blake3::ContentHash;
use fs_qty::Dims;

/// SI exponents, `[m, kg, s, K, …]`.
const TEMPERATURE_DIMS: Dims = Dims([0, 0, 0, 1, 0, 0]);
const POWER_DIMS: Dims = Dims([2, 1, -3, 0, 0, 0]);
const VOLUME_FLOW_DIMS: Dims = Dims([3, 0, -1, 0, 0, 0]);
const PRESSURE_DIMS: Dims = Dims([-1, 1, -2, 0, 0, 0]);

const RIG_MAGIC: &[u8; 8] = b"FSVVRIG\0";
const RIG_IDENTITY_DOMAIN: &str = "org.frankensim.fs-vvreg.rig-run.v2";
const RAW_READINGS_DOMAIN: &str = "org.frankensim.fs-vvreg.rig-readings.v1";
const BALANCE_MODEL_TAG: u8 = 1;
const ABSOLUTE_INTERVAL_POLICY_TAG: u8 = 1;
const OUTWARD_ROUNDING_POLICY_TAG: u8 = 1;

/// Schema version of the retained rig-run record.
pub const RIG_SCHEMA_VERSION: u32 = 2;

/// Maximum instruments admitted on one rig.
pub const MAX_INSTRUMENTS: usize = 1_024;
/// Maximum canonical bytes retained for one rig run.
pub const MAX_RETAINED_RIG_BYTES: usize = 16 * 1024 * 1024;

/// What one instrument channel measures, and how the balance uses it.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum ChannelRole {
    /// Electrical power into the heated plate, W.
    HeaterPower,
    /// Coolant temperature entering the enclosure, K.
    InletTemperature,
    /// Coolant temperature leaving the enclosure, K.
    OutletTemperature,
    /// Volumetric coolant flow, m³/s.
    VolumeFlow,
    /// Ambient temperature, K. Recorded, not part of the closure.
    Ambient,
    /// Any other measured channel: recorded and calibration-checked, but the
    /// balance does not consume it.
    Auxiliary,
}

impl ChannelRole {
    /// Stable lowercase role name.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::HeaterPower => "heater-power",
            Self::InletTemperature => "inlet-temperature",
            Self::OutletTemperature => "outlet-temperature",
            Self::VolumeFlow => "volume-flow",
            Self::Ambient => "ambient",
            Self::Auxiliary => "auxiliary",
        }
    }

    /// The dimensions a channel in this role must declare, if the role fixes
    /// them. `Auxiliary` deliberately does not.
    #[must_use]
    pub const fn expected_dims(self) -> Option<Dims> {
        match self {
            Self::HeaterPower => Some(POWER_DIMS),
            Self::InletTemperature | Self::OutletTemperature | Self::Ambient => {
                Some(TEMPERATURE_DIMS)
            }
            Self::VolumeFlow => Some(VOLUME_FLOW_DIMS),
            Self::Auxiliary => None,
        }
    }

    /// Roles the energy-balance closure requires exactly once.
    #[must_use]
    pub const fn balance_roles() -> [Self; 4] {
        [
            Self::HeaterPower,
            Self::InletTemperature,
            Self::OutletTemperature,
            Self::VolumeFlow,
        ]
    }
}

/// A structured rig specification or ingest failure.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RigError {
    /// A required identifier was blank.
    BlankIdentifier {
        /// Which identifier.
        field: &'static str,
    },
    /// A declared scalar was not finite, or violated its sign requirement.
    InvalidScalar {
        /// Offending field.
        field: &'static str,
        /// Violated requirement.
        requirement: &'static str,
    },
    /// A required identifier violated the bounded canonical slug grammar.
    InvalidIdentifier {
        /// Offending identifier.
        field: &'static str,
        /// Stable grammar requirement.
        requirement: &'static str,
    },
    /// A corpus sensor row failed the canonical sensor admission boundary.
    InvalidSensor {
        /// Sensor identifier, if one was supplied.
        sensor_id: String,
        /// Canonical corpus refusal.
        detail: String,
    },
    /// One run reading is malformed.
    InvalidReading {
        /// Sensor identifier carried by the reading.
        sensor_id: String,
        /// Offending reading field.
        field: &'static str,
        /// Stable requirement.
        requirement: &'static str,
    },
    /// An acquisition date is not a real canonical ISO date.
    InvalidDate {
        /// Offending date field.
        field: &'static str,
        /// Supplied value.
        value: String,
    },
    /// A run lies outside one sensor certificate's declared validity window.
    CalibrationOutOfWindow {
        /// Sensor identifier.
        sensor_id: String,
        /// Run acquisition date.
        acquired_on: String,
        /// Certificate issue date.
        issued_on: String,
        /// Optional last-valid date.
        valid_through: Option<String>,
    },
    /// More instruments than one rig admits.
    TooManyInstruments {
        /// Supplied count.
        have: usize,
        /// Admitted maximum.
        max: usize,
    },
    /// More readings than one run admits.
    TooManyReadings {
        /// Supplied count.
        have: usize,
        /// Admitted maximum.
        max: usize,
    },
    /// A channel's declared dimensions do not match its role.
    RoleDimensions {
        /// Sensor identifier.
        sensor_id: String,
        /// Declared role.
        role: ChannelRole,
    },
    /// Two instruments share one sensor identifier.
    DuplicateSensor {
        /// Repeated identifier.
        sensor_id: String,
    },
    /// A balance role is missing from the specification.
    MissingRole {
        /// The absent role.
        role: ChannelRole,
    },
    /// A balance role is declared more than once.
    DuplicateRole {
        /// The repeated role.
        role: ChannelRole,
    },
    /// A sensor has no calibration certificate.
    ///
    /// Level E exists to own its metrology; an uncalibrated channel is not a
    /// weaker measurement, it is an unknown one.
    Uncalibrated {
        /// Sensor identifier.
        sensor_id: String,
        /// The declared reason the certificate is absent.
        reason: String,
    },
    /// A run supplied a reading for a sensor the rig does not declare.
    UnknownReading {
        /// Sensor identifier in the run.
        sensor_id: String,
    },
    /// A declared instrument has no reading in the run.
    MissingReading {
        /// Sensor identifier in the specification.
        sensor_id: String,
    },
    /// A run supplied two readings for one sensor.
    DuplicateReading {
        /// Repeated sensor identifier.
        sensor_id: String,
    },
    /// The measured energy does not close within the instruments' own bands.
    EnergyBalanceViolated {
        /// Rendered comparison, including every contributing channel.
        detail: String,
    },
    /// The balance cannot be adjudicated, so it cannot be passed.
    EnergyBalanceUnadjudicated {
        /// Why no band exists.
        reason: String,
    },
    /// A measured support lies wholly in the resolved wrong direction.
    PhysicalDirectionViolated {
        /// Physical interval whose sign is resolved incorrectly.
        field: &'static str,
    },
    /// A retained rig record exceeds its fixed resource cap.
    RetainedRecordTooLarge {
        /// Supplied byte count.
        have: usize,
        /// Admitted maximum.
        max: usize,
    },
    /// Retained bytes use an unsupported rig schema.
    UnsupportedSchema {
        /// Observed schema version.
        observed: u32,
    },
    /// Retained bytes cannot be decoded as a complete rig record.
    MalformedRetainedRecord {
        /// Stable refusal reason.
        reason: &'static str,
    },
    /// Decoded semantics do not re-encode to the supplied bytes.
    NonCanonicalRetainedRecord,
    /// A capacity reservation failed.
    AllocationFailed,
}

impl core::fmt::Display for RigError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::BlankIdentifier { field } => write!(f, "{field} is blank"),
            Self::InvalidScalar { field, requirement }
            | Self::InvalidIdentifier { field, requirement } => {
                write!(f, "{field} must be {requirement}")
            }
            Self::InvalidSensor { sensor_id, detail } => {
                write!(f, "sensor {sensor_id:?} is invalid: {detail}")
            }
            Self::InvalidReading {
                sensor_id,
                field,
                requirement,
            } => write!(
                f,
                "reading for sensor {sensor_id:?} has invalid {field}: must be {requirement}"
            ),
            Self::InvalidDate { field, value } => {
                write!(f, "{field} {value:?} is not a real YYYY-MM-DD date")
            }
            Self::CalibrationOutOfWindow {
                sensor_id,
                acquired_on,
                issued_on,
                valid_through,
            } => write!(
                f,
                "sensor {sensor_id:?} calibration [{issued_on}, {}] does not cover acquisition {acquired_on}",
                valid_through.as_deref().unwrap_or("open")
            ),
            Self::TooManyInstruments { have, max } => {
                write!(f, "{have} instruments exceeds the admitted maximum {max}")
            }
            Self::TooManyReadings { have, max } => {
                write!(f, "{have} readings exceeds the admitted maximum {max}")
            }
            Self::RoleDimensions { sensor_id, role } => write!(
                f,
                "sensor {sensor_id:?} declares dimensions incompatible with role {}",
                role.name()
            ),
            Self::DuplicateSensor { sensor_id } => {
                write!(f, "sensor identifier {sensor_id:?} is declared twice")
            }
            Self::MissingRole { role } => write!(
                f,
                "the rig declares no {} channel, so the energy balance cannot close",
                role.name()
            ),
            Self::DuplicateRole { role } => {
                write!(f, "the rig declares more than one {} channel", role.name())
            }
            Self::Uncalibrated { sensor_id, reason } => write!(
                f,
                "sensor {sensor_id:?} has no calibration certificate ({reason})"
            ),
            Self::UnknownReading { sensor_id } => write!(
                f,
                "the run reports sensor {sensor_id:?}, which the rig does not declare"
            ),
            Self::MissingReading { sensor_id } => {
                write!(
                    f,
                    "the run has no reading for declared sensor {sensor_id:?}"
                )
            }
            Self::DuplicateReading { sensor_id } => {
                write!(f, "the run reports sensor {sensor_id:?} twice")
            }
            Self::EnergyBalanceViolated { detail } => {
                write!(f, "the measured energy does not close:\n{detail}")
            }
            Self::EnergyBalanceUnadjudicated { reason } => {
                write!(f, "the energy balance cannot be adjudicated: {reason}")
            }
            Self::PhysicalDirectionViolated { field } => {
                write!(
                    f,
                    "{field} lies wholly outside the required positive direction"
                )
            }
            Self::RetainedRecordTooLarge { have, max } => {
                write!(
                    f,
                    "{have} retained rig bytes exceeds the admitted maximum {max}"
                )
            }
            Self::UnsupportedSchema { observed } => {
                write!(f, "unsupported retained rig schema {observed}")
            }
            Self::MalformedRetainedRecord { reason } => {
                write!(f, "malformed retained rig record: {reason}")
            }
            Self::NonCanonicalRetainedRecord => {
                f.write_str("retained rig bytes are semantically valid but not canonical")
            }
            Self::AllocationFailed => write!(f, "a rig allocation failed"),
        }
    }
}

impl std::error::Error for RigError {}

/// Declared coolant properties used by the closure.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CoolantProperties {
    density_kg_per_m3: f64,
    specific_heat_j_per_kg_k: f64,
    relative_half_width: f64,
}

impl CoolantProperties {
    /// Declare coolant density and specific heat with a shared relative
    /// half-width covering both.
    ///
    /// The band is declared rather than looked up: property tables carry
    /// their own uncertainty, and folding it in here keeps the balance
    /// adjudication honest about every term it uses, not only the measured
    /// ones.
    ///
    /// # Errors
    /// A non-positive density or specific heat, an unrepresentable property
    /// product, or a relative half-width outside `[0, 1)`.
    pub fn new(
        density_kg_per_m3: f64,
        specific_heat_j_per_kg_k: f64,
        relative_half_width: f64,
    ) -> Result<Self, RigError> {
        for (value, field) in [
            (density_kg_per_m3, "coolant density"),
            (specific_heat_j_per_kg_k, "coolant specific heat"),
        ] {
            if !value.is_finite() || value <= 0.0 {
                return Err(RigError::InvalidScalar {
                    field,
                    requirement: "finite and positive",
                });
            }
        }
        let property_product = density_kg_per_m3 * specific_heat_j_per_kg_k;
        if !property_product.is_finite() || property_product <= 0.0 {
            return Err(RigError::InvalidScalar {
                field: "coolant density-specific-heat product",
                requirement: "finite and positive",
            });
        }
        if !relative_half_width.is_finite()
            || relative_half_width < 0.0
            || relative_half_width >= 1.0
        {
            return Err(RigError::InvalidScalar {
                field: "coolant relative half-width",
                requirement: "finite and in [0, 1)",
            });
        }
        Ok(Self {
            density_kg_per_m3,
            specific_heat_j_per_kg_k,
            relative_half_width,
        })
    }

    /// Declared density, kg/m³.
    #[must_use]
    pub const fn density_kg_per_m3(self) -> f64 {
        self.density_kg_per_m3
    }

    /// Declared specific heat, J/(kg·K).
    #[must_use]
    pub const fn specific_heat_j_per_kg_k(self) -> f64 {
        self.specific_heat_j_per_kg_k
    }

    /// Declared relative half-width on the property product.
    #[must_use]
    pub const fn relative_half_width(self) -> f64 {
        self.relative_half_width
    }
}

/// One instrumented channel: its corpus sensor record plus its rig role.
#[derive(Debug, Clone, PartialEq)]
pub struct InstrumentSpec {
    sensor: SensorRecord,
    role: ChannelRole,
}

impl InstrumentSpec {
    /// Bind a corpus sensor record to a rig role.
    ///
    /// Reusing [`SensorRecord`] is deliberate: calibration, placement and
    /// uncertainty then have ONE representation shared with the corpus, so an
    /// ingested run cannot describe its metrology differently from the
    /// dataset it becomes.
    ///
    /// # Errors
    /// Any canonical corpus sensor validation failure, or dimensions
    /// incompatible with the role.
    pub fn new(sensor: SensorRecord, role: ChannelRole) -> Result<Self, RigError> {
        validate_sensor_record(&sensor).map_err(|error| RigError::InvalidSensor {
            sensor_id: diagnostic_sensor_id(&sensor.id),
            detail: error.to_string(),
        })?;
        if let Some(expected) = role.expected_dims()
            && sensor.quantity_dims != expected
        {
            return Err(RigError::RoleDimensions {
                sensor_id: sensor.id.clone(),
                role,
            });
        }
        Ok(Self { sensor, role })
    }

    /// The corpus sensor record.
    #[must_use]
    pub const fn sensor(&self) -> &SensorRecord {
        &self.sensor
    }

    /// The rig role.
    #[must_use]
    pub const fn role(&self) -> ChannelRole {
        self.role
    }

    /// Stated finite-support half-width in SI units.
    ///
    /// A covariance is not silently converted into a support interval. The
    /// balance gate has no identity-bound coverage multiplier, so covariance
    /// and unstated cases remain unadjudicable.
    fn bounded_half_width_si(&self) -> Result<f64, RigError> {
        match &self.sensor.uncertainty {
            MeasurementUncertainty::Bounded { half_width } => Ok(half_width.value),
            MeasurementUncertainty::CovarianceDiagonal { .. } => {
                Err(RigError::EnergyBalanceUnadjudicated {
                    reason: format!(
                        "balance sensor {:?} declares only CovarianceDiagonal; one sigma is not a finite support band",
                        self.sensor.id
                    ),
                })
            }
            MeasurementUncertainty::Unstated => Err(RigError::EnergyBalanceUnadjudicated {
                reason: format!(
                    "balance sensor {:?} declares MeasurementUncertainty::Unstated",
                    self.sensor.id
                ),
            }),
        }
    }
}

/// A rig specification: the instrument set and the coolant it moves.
#[derive(Debug, Clone, PartialEq)]
pub struct RigSpec {
    id: String,
    coolant: CoolantProperties,
    instruments: Vec<InstrumentSpec>,
}

impl RigSpec {
    /// Admit a rig specification.
    ///
    /// Every instrument must carry a calibration certificate, sensor
    /// identifiers must be unique, and each of
    /// [`ChannelRole::balance_roles`] must appear exactly once — the closure
    /// is not optional equipment.
    ///
    /// # Errors
    /// [`RigError`] for a blank rig id, an empty or oversized instrument set,
    /// a duplicated sensor id, an uncalibrated sensor, or a missing or
    /// duplicated balance role.
    pub fn new(
        id: impl Into<String>,
        coolant: CoolantProperties,
        instruments: Vec<InstrumentSpec>,
    ) -> Result<Self, RigError> {
        let id = id.into();
        if id.trim().is_empty() {
            return Err(RigError::BlankIdentifier { field: "rig id" });
        }
        if !valid_slug(&id) {
            return Err(RigError::InvalidIdentifier {
                field: "rig id",
                requirement: "a bounded lowercase ASCII slug",
            });
        }
        if instruments.is_empty() {
            return Err(RigError::MissingRole {
                role: ChannelRole::HeaterPower,
            });
        }
        if instruments.len() > MAX_INSTRUMENTS {
            return Err(RigError::TooManyInstruments {
                have: instruments.len(),
                max: MAX_INSTRUMENTS,
            });
        }

        let mut instruments = instruments;
        instruments.sort_by(|a, b| a.sensor.id.cmp(&b.sensor.id));
        for pair in instruments.windows(2) {
            if pair[0].sensor.id == pair[1].sensor.id {
                return Err(RigError::DuplicateSensor {
                    sensor_id: pair[0].sensor.id.clone(),
                });
            }
        }
        for instrument in &instruments {
            if let Availability::Unavailable { reason } = &instrument.sensor.calibration {
                return Err(RigError::Uncalibrated {
                    sensor_id: instrument.sensor.id.clone(),
                    reason: reason.clone(),
                });
            }
        }
        for role in ChannelRole::balance_roles() {
            let count = instruments
                .iter()
                .filter(|instrument| instrument.role == role)
                .count();
            match count {
                1 => {}
                0 => return Err(RigError::MissingRole { role }),
                _ => return Err(RigError::DuplicateRole { role }),
            }
        }

        Ok(Self {
            id,
            coolant,
            instruments,
        })
    }

    /// Stable rig identifier.
    #[must_use]
    pub fn id(&self) -> &str {
        &self.id
    }

    /// Declared coolant properties.
    #[must_use]
    pub const fn coolant(&self) -> CoolantProperties {
        self.coolant
    }

    /// Instruments in deterministic sensor-id order.
    #[must_use]
    pub fn instruments(&self) -> &[InstrumentSpec] {
        &self.instruments
    }

    /// The single instrument holding one balance role.
    #[must_use]
    fn role_instrument(&self, role: ChannelRole) -> Option<&InstrumentSpec> {
        self.instruments
            .iter()
            .find(|instrument| instrument.role == role)
    }
}

/// One channel reading, in coherent SI units.
#[derive(Debug, Clone, PartialEq)]
pub struct Reading {
    /// Sensor identifier, matching a declared instrument.
    pub sensor_id: String,
    /// Measured value in coherent SI units.
    pub value_si: f64,
}

/// One acquisition run, before ingest.
#[derive(Debug, Clone, PartialEq)]
pub struct RigRun {
    /// Stable run identifier.
    pub run_id: String,
    /// Acquisition date used to check every calibration validity window.
    pub acquired_on: String,
    /// Declared partition. Declared AT INGEST, never inferred afterwards.
    pub partition: DatasetPartition,
    /// Channel readings.
    pub readings: Vec<Reading>,
}

/// One finite closed interval in coherent SI units.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ClosedInterval {
    lo: f64,
    hi: f64,
}

impl ClosedInterval {
    /// Inclusive lower endpoint.
    #[must_use]
    pub const fn lo(self) -> f64 {
        self.lo
    }

    /// Inclusive upper endpoint.
    #[must_use]
    pub const fn hi(self) -> f64 {
        self.hi
    }

    /// Whether this interval contains a finite scalar.
    #[must_use]
    pub fn contains(self, value: f64) -> bool {
        value.is_finite() && self.lo <= value && value <= self.hi
    }

    fn around(center: f64, half_width: f64, field: &'static str) -> Result<Self, RigError> {
        let lo = next_down(center - half_width);
        let hi = next_up(center + half_width);
        Self::new_finite(lo, hi, field)
    }

    fn subtract(self, rhs: Self, field: &'static str) -> Result<Self, RigError> {
        let lo = next_down(self.lo - rhs.hi);
        let hi = next_up(self.hi - rhs.lo);
        Self::new_finite(lo, hi, field)
    }

    fn multiply(self, rhs: Self, field: &'static str) -> Result<Self, RigError> {
        let products = [
            (self.lo, rhs.lo),
            (self.lo, rhs.hi),
            (self.hi, rhs.lo),
            (self.hi, rhs.hi),
        ];
        let mut lo = f64::INFINITY;
        let mut hi = f64::NEG_INFINITY;
        for (lhs, rhs) in products {
            lo = lo.min(next_down(lhs * rhs));
            hi = hi.max(next_up(lhs * rhs));
        }
        Self::new_finite(lo, hi, field)
    }

    fn new_finite(lo: f64, hi: f64, field: &'static str) -> Result<Self, RigError> {
        if !lo.is_finite() || !hi.is_finite() || lo > hi {
            return Err(RigError::InvalidScalar {
                field,
                requirement: "a finite ordered outward enclosure",
            });
        }
        Ok(Self { lo, hi })
    }
}

/// The adjudicated energy-closure verdict retained with an ingested run.
#[derive(Debug, Clone, PartialEq)]
pub enum EnergyBalance {
    /// The closure holds within the instruments' combined stated band.
    Balanced {
        /// Electrical power in, W.
        power_in_w: f64,
        /// Outward enclosure of electrical input power, W.
        power_interval_w: ClosedInterval,
        /// Enthalpy rise carried by the coolant, W.
        heat_out_w: f64,
        /// Outward enclosure of coolant enthalpy rise, W.
        heat_interval_w: ClosedInterval,
        /// Signed `power_in − heat_out`, W.
        residual_w: f64,
        /// Exact outward residual enclosure `power − heat`, W.
        residual_interval_w: ClosedInterval,
    },
    /// The run passed the gate, but the numbers are withheld because the run
    /// is a sealed blind holdout.
    Sealed,
}

/// An admitted rig run.
#[derive(Clone, PartialEq)]
pub struct IngestedRun {
    run_id: String,
    rig_id: String,
    acquired_on: String,
    partition: DatasetPartition,
    sealed: bool,
    balance: EnergyBalance,
    calibrated_channels: usize,
    readings: Vec<Reading>,
    instruments: Vec<InstrumentSpec>,
    raw_readings_identity: ContentHash,
    retained_bytes: Box<[u8]>,
    identity: ContentHash,
}

#[derive(Clone, Copy)]
struct BalanceComputation {
    power_in: f64,
    power_interval: ClosedInterval,
    heat_out: f64,
    heat_interval: ClosedInterval,
    residual: f64,
    residual_interval: ClosedInterval,
    flow: f64,
    delta_t: f64,
}

impl core::fmt::Debug for IngestedRun {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        let mut debug = formatter.debug_struct("IngestedRun");
        debug
            .field("run_id", &self.run_id)
            .field("rig_id", &self.rig_id)
            .field("acquired_on", &self.acquired_on)
            .field("partition", &self.partition)
            .field("sealed", &self.sealed)
            .field("balance", &self.balance)
            .field("calibrated_channels", &self.calibrated_channels)
            .field("instruments", &self.instruments)
            .field("raw_readings_identity", &self.raw_readings_identity)
            .field("retained_byte_len", &self.retained_bytes.len())
            .field("identity", &self.identity);
        if self.sealed {
            debug.field("readings", &"<sealed>");
        } else {
            debug.field("readings", &self.readings);
        }
        debug.finish()
    }
}

impl IngestedRun {
    /// Stable run identifier.
    #[must_use]
    pub fn run_id(&self) -> &str {
        &self.run_id
    }

    /// The rig that produced it.
    #[must_use]
    pub fn rig_id(&self) -> &str {
        &self.rig_id
    }

    /// Acquisition date against which calibration windows were checked.
    #[must_use]
    pub fn acquired_on(&self) -> &str {
        &self.acquired_on
    }

    /// Declared partition.
    #[must_use]
    pub const fn partition(&self) -> DatasetPartition {
        self.partition
    }

    /// Whether the measured values are sealed from ordinary use.
    #[must_use]
    pub const fn sealed(&self) -> bool {
        self.sealed
    }

    /// The retained closure verdict. [`EnergyBalance::Sealed`] for a blind
    /// holdout: the gate ran and passed, and the numbers do not leak.
    #[must_use]
    pub const fn balance(&self) -> &EnergyBalance {
        &self.balance
    }

    /// Channels whose calibration certificate was verified present.
    #[must_use]
    pub const fn calibrated_channels(&self) -> usize {
        self.calibrated_channels
    }

    /// Canonical readings in sensor-id order for an unsealed run.
    ///
    /// Blind values stay off this ordinary-use surface. The explicit
    /// [`Self::retained_bytes`] persistence surface still carries them.
    #[must_use]
    pub fn readings(&self) -> Option<&[Reading]> {
        (!self.sealed).then_some(self.readings.as_slice())
    }

    /// Exact admitted metrology snapshot in sensor-id order.
    #[must_use]
    pub fn instruments(&self) -> &[InstrumentSpec] {
        &self.instruments
    }

    /// Content identity of the canonical raw-reading rows.
    #[must_use]
    pub const fn raw_readings_identity(&self) -> ContentHash {
        self.raw_readings_identity
    }

    /// Exact canonical retained evidence bytes.
    ///
    /// This explicit persistence surface includes blind raw values. It is not
    /// an ordinary validation accessor and supplies no secrecy or access
    /// control; callers storing blind bytes must protect the artifact.
    #[must_use]
    pub fn retained_bytes(&self) -> &[u8] {
        &self.retained_bytes
    }

    /// Decode, revalidate, re-ingest, and canonical-byte-check a retained run.
    ///
    /// # Errors
    /// [`RigError`] when the bytes are malformed, noncanonical, out of bounds,
    /// or no longer satisfy the fail-closed ingest gate.
    pub fn from_retained_bytes(bytes: &[u8]) -> Result<Self, RigError> {
        decode_retained_run(bytes)
    }

    /// Domain-separated identity of the admitted run. An integrity address,
    /// not authentication.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
}

/// Ingest one run against its rig, gating on channel match, calibration, and
/// the measured energy balance.
///
/// # Errors
/// [`RigError`] for an unknown, missing, duplicated, or non-finite reading;
/// an uncalibrated channel; or an energy balance that is violated or cannot
/// be adjudicated.
pub fn ingest(spec: &RigSpec, run: &RigRun) -> Result<IngestedRun, RigError> {
    let sorted = canonical_readings(spec, run)?;
    let computed = compute_balance(spec, &sorted)?;

    if !computed.residual_interval.contains(0.0) {
        let detail = if run.partition == DatasetPartition::BlindHoldout {
            render_blind_balance(spec)
        } else {
            render_balance(spec, &computed)
        };
        return Err(RigError::EnergyBalanceViolated { detail });
    }

    let mut retained_readings = Vec::new();
    retained_readings
        .try_reserve_exact(sorted.len())
        .map_err(|_| RigError::AllocationFailed)?;
    retained_readings.extend(sorted.into_iter().cloned());
    let retained_bytes = encode_retained_run(spec, run, &retained_readings)?;
    let raw_reading_bytes = encode_readings(&retained_readings)?;
    let raw_readings_identity = fs_blake3::hash_domain(RAW_READINGS_DOMAIN, &raw_reading_bytes);
    let identity = fs_blake3::hash_domain(RIG_IDENTITY_DOMAIN, &retained_bytes);

    let sealed = run.partition == DatasetPartition::BlindHoldout;
    let balance = if sealed {
        EnergyBalance::Sealed
    } else {
        EnergyBalance::Balanced {
            power_in_w: computed.power_in,
            power_interval_w: computed.power_interval,
            heat_out_w: computed.heat_out,
            heat_interval_w: computed.heat_interval,
            residual_w: computed.residual,
            residual_interval_w: computed.residual_interval,
        }
    };

    Ok(IngestedRun {
        run_id: run.run_id.clone(),
        rig_id: spec.id.clone(),
        acquired_on: run.acquired_on.clone(),
        partition: run.partition,
        sealed,
        balance,
        calibrated_channels: spec.instruments.len(),
        readings: retained_readings,
        instruments: spec.instruments.clone(),
        raw_readings_identity,
        retained_bytes: retained_bytes.into_boxed_slice(),
        identity,
    })
}

fn canonical_readings<'a>(spec: &RigSpec, run: &'a RigRun) -> Result<Vec<&'a Reading>, RigError> {
    if run.run_id.trim().is_empty() {
        return Err(RigError::BlankIdentifier { field: "run id" });
    }
    if !valid_slug(&run.run_id) {
        return Err(RigError::InvalidIdentifier {
            field: "run id",
            requirement: "a bounded lowercase ASCII slug",
        });
    }
    if !valid_date(&run.acquired_on) {
        return Err(RigError::InvalidDate {
            field: "run acquisition date",
            value: run.acquired_on.clone(),
        });
    }
    if run.readings.len() > MAX_INSTRUMENTS {
        return Err(RigError::TooManyReadings {
            have: run.readings.len(),
            max: MAX_INSTRUMENTS,
        });
    }

    let mut sorted = Vec::new();
    sorted
        .try_reserve_exact(run.readings.len())
        .map_err(|_| RigError::AllocationFailed)?;
    for reading in &run.readings {
        if !valid_slug(&reading.sensor_id) {
            return Err(RigError::InvalidReading {
                sensor_id: diagnostic_sensor_id(&reading.sensor_id),
                field: "sensor id",
                requirement: "a bounded lowercase ASCII slug",
            });
        }
        if !reading.value_si.is_finite() {
            return Err(RigError::InvalidReading {
                sensor_id: reading.sensor_id.clone(),
                field: "value_si",
                requirement: "finite",
            });
        }
        sorted.push(reading);
    }
    sorted.sort_by(|a, b| a.sensor_id.cmp(&b.sensor_id));
    for pair in sorted.windows(2) {
        if pair[0].sensor_id == pair[1].sensor_id {
            return Err(RigError::DuplicateReading {
                sensor_id: pair[0].sensor_id.clone(),
            });
        }
    }
    for reading in &sorted {
        if spec
            .instruments
            .binary_search_by(|instrument| instrument.sensor.id.cmp(&reading.sensor_id))
            .is_err()
        {
            return Err(RigError::UnknownReading {
                sensor_id: reading.sensor_id.clone(),
            });
        }
    }
    validate_instruments_for_run(spec, run, &sorted)?;
    Ok(sorted)
}

fn validate_instruments_for_run(
    spec: &RigSpec,
    run: &RigRun,
    sorted: &[&Reading],
) -> Result<(), RigError> {
    for instrument in &spec.instruments {
        if sorted
            .binary_search_by(|reading| reading.sensor_id.cmp(&instrument.sensor.id))
            .is_err()
        {
            return Err(RigError::MissingReading {
                sensor_id: instrument.sensor.id.clone(),
            });
        }
        validate_sensor_record(&instrument.sensor).map_err(|error| RigError::InvalidSensor {
            sensor_id: instrument.sensor.id.clone(),
            detail: error.to_string(),
        })?;
        // RigSpec::new already refuses an uncalibrated instrument; re-checking
        // here keeps ingest self-standing rather than relying on a
        // precondition a future caller could route around.
        match &instrument.sensor.calibration {
            Availability::Unavailable { reason } => {
                return Err(RigError::Uncalibrated {
                    sensor_id: instrument.sensor.id.clone(),
                    reason: reason.clone(),
                });
            }
            Availability::Available(calibration)
                if run.acquired_on.as_str() < calibration.issued_on.as_str()
                    || calibration
                        .valid_through
                        .as_ref()
                        .is_some_and(|last| run.acquired_on.as_str() > last.as_str()) =>
            {
                return Err(RigError::CalibrationOutOfWindow {
                    sensor_id: instrument.sensor.id.clone(),
                    acquired_on: run.acquired_on.clone(),
                    issued_on: calibration.issued_on.clone(),
                    valid_through: calibration.valid_through.clone(),
                });
            }
            Availability::Available(_) => {}
        }
    }
    Ok(())
}

fn compute_balance(spec: &RigSpec, sorted: &[&Reading]) -> Result<BalanceComputation, RigError> {
    let value_of = |role: ChannelRole| -> f64 {
        let instrument = spec
            .role_instrument(role)
            .expect("RigSpec::new admitted every balance role");
        sorted
            .iter()
            .find(|reading| reading.sensor_id == instrument.sensor.id)
            .expect("every declared instrument has a reading")
            .value_si
    };
    let band_of = |role: ChannelRole| -> Result<f64, RigError> {
        spec.role_instrument(role)
            .expect("RigSpec::new admitted every balance role")
            .bounded_half_width_si()
    };

    let power_in = value_of(ChannelRole::HeaterPower);
    let inlet = value_of(ChannelRole::InletTemperature);
    let outlet = value_of(ChannelRole::OutletTemperature);
    let flow = value_of(ChannelRole::VolumeFlow);

    let power_interval = ClosedInterval::around(
        power_in,
        band_of(ChannelRole::HeaterPower)?,
        "power interval",
    )?;
    let inlet_interval = ClosedInterval::around(
        inlet,
        band_of(ChannelRole::InletTemperature)?,
        "inlet-temperature interval",
    )?;
    let outlet_interval = ClosedInterval::around(
        outlet,
        band_of(ChannelRole::OutletTemperature)?,
        "outlet-temperature interval",
    )?;
    let flow_interval =
        ClosedInterval::around(flow, band_of(ChannelRole::VolumeFlow)?, "flow interval")?;
    let delta_interval = outlet_interval.subtract(inlet_interval, "coolant-rise interval")?;

    require_strictly_positive(power_interval, "heater-power interval")?;
    require_strictly_positive(flow_interval, "volume-flow interval")?;
    require_strictly_positive(delta_interval, "coolant-rise interval")?;

    let coolant = spec.coolant;
    let property_product = coolant.density_kg_per_m3 * coolant.specific_heat_j_per_kg_k;
    let property_interval = coolant_property_interval(coolant)?;
    require_strictly_positive(property_interval, "coolant-property interval")?;

    let heat_interval = property_interval
        .multiply(flow_interval, "coolant-property-flow interval")?
        .multiply(delta_interval, "heat-output interval")?;
    let delta_t = outlet - inlet;
    let heat_out = property_product * flow * delta_t;
    let residual = power_in - heat_out;
    if !heat_out.is_finite() || !residual.is_finite() {
        return Err(RigError::InvalidScalar {
            field: "nominal energy balance",
            requirement: "finite",
        });
    }
    let residual_interval = power_interval.subtract(heat_interval, "energy-residual interval")?;

    Ok(BalanceComputation {
        power_in,
        power_interval,
        heat_out,
        heat_interval,
        residual,
        residual_interval,
        flow,
        delta_t,
    })
}

fn render_balance(spec: &RigSpec, computed: &BalanceComputation) -> String {
    use core::fmt::Write as _;
    let mut out = String::new();
    let _ = writeln!(
        out,
        "  rig {}: power in {} W, heat out {} W",
        spec.id, computed.power_in, computed.heat_out
    );
    let _ = writeln!(
        out,
        "  power interval [{}, {}] W and heat interval [{}, {}] W do not overlap",
        computed.power_interval.lo,
        computed.power_interval.hi,
        computed.heat_interval.lo,
        computed.heat_interval.hi
    );
    let _ = writeln!(
        out,
        "  residual {} W has outward enclosure [{}, {}] W",
        computed.residual, computed.residual_interval.lo, computed.residual_interval.hi
    );
    let _ = writeln!(
        out,
        "  flow {} m^3/s, coolant rise {} K, rho {} kg/m^3, c_p {} J/(kg K)",
        computed.flow,
        computed.delta_t,
        spec.coolant.density_kg_per_m3,
        spec.coolant.specific_heat_j_per_kg_k
    );
    for instrument in &spec.instruments {
        let band = match instrument.sensor.uncertainty {
            MeasurementUncertainty::Bounded { half_width } => {
                format!("{}", half_width.value)
            }
            MeasurementUncertainty::CovarianceDiagonal { .. } => "covariance-only".to_string(),
            MeasurementUncertainty::Unstated => "unstated".to_string(),
        };
        let _ = writeln!(
            out,
            "  channel {} [{}]: half-width {band}",
            instrument.sensor.id,
            instrument.role.name()
        );
    }
    out
}

fn render_blind_balance(spec: &RigSpec) -> String {
    use core::fmt::Write as _;

    let mut out = String::new();
    let _ = writeln!(
        out,
        "rig {}: blind-holdout energy-balance intervals do not overlap",
        spec.id
    );
    let _ = writeln!(out, "  measured values and interval endpoints: REDACTED");
    for instrument in &spec.instruments {
        let _ = writeln!(
            out,
            "  channel {} [{}]",
            instrument.sensor.id,
            instrument.role.name()
        );
    }
    out
}

fn require_strictly_positive(
    interval: ClosedInterval,
    field: &'static str,
) -> Result<(), RigError> {
    if interval.hi <= 0.0 {
        return Err(RigError::PhysicalDirectionViolated { field });
    }
    if interval.lo <= 0.0 {
        return Err(RigError::EnergyBalanceUnadjudicated {
            reason: format!(
                "{field} touches or spans zero, so its positive direction is unresolved"
            ),
        });
    }
    Ok(())
}

fn coolant_property_interval(coolant: CoolantProperties) -> Result<ClosedInterval, RigError> {
    // `rho * cp` must itself be enclosed before the declared relative band is
    // applied. Rounding the product to one f64 first and then widening
    // `product +/- product*r` can remain one ULP inside the exact real
    // expression for adversarial neighboring inputs.
    let density = ClosedInterval::new_finite(
        coolant.density_kg_per_m3,
        coolant.density_kg_per_m3,
        "coolant-density point interval",
    )?;
    let specific_heat = ClosedInterval::new_finite(
        coolant.specific_heat_j_per_kg_k,
        coolant.specific_heat_j_per_kg_k,
        "coolant-specific-heat point interval",
    )?;
    let nominal_product = density.multiply(specific_heat, "coolant-property product interval")?;
    let relative_factor = ClosedInterval::around(
        1.0,
        coolant.relative_half_width,
        "coolant-property relative-factor interval",
    )?;
    nominal_product.multiply(relative_factor, "coolant-property interval")
}

fn encode_retained_run(
    spec: &RigSpec,
    run: &RigRun,
    readings: &[Reading],
) -> Result<Vec<u8>, RigError> {
    let mut sensor_rows = Vec::new();
    sensor_rows
        .try_reserve_exact(spec.instruments.len())
        .map_err(|_| RigError::AllocationFailed)?;
    for instrument in &spec.instruments {
        let row =
            encode_sensor_record(&instrument.sensor).map_err(|error| RigError::InvalidSensor {
                sensor_id: instrument.sensor.id.clone(),
                detail: error.to_string(),
            })?;
        sensor_rows.push((role_tag(instrument.role), row));
    }

    let mut estimated = RIG_MAGIC.len() + 4 + 3 + 4 + 4 + 24 + 4 + 4 + 4 + 1 + 4;
    estimated = checked_encoded_add(estimated, spec.id.len())?;
    estimated = checked_encoded_add(estimated, run.run_id.len())?;
    estimated = checked_encoded_add(estimated, run.acquired_on.len())?;
    for (_, row) in &sensor_rows {
        estimated = checked_encoded_add(estimated, 1 + 4 + row.len())?;
    }
    for reading in readings {
        estimated = checked_encoded_add(estimated, 4 + reading.sensor_id.len() + 8)?;
    }
    if estimated > MAX_RETAINED_RIG_BYTES {
        return Err(RigError::RetainedRecordTooLarge {
            have: estimated,
            max: MAX_RETAINED_RIG_BYTES,
        });
    }

    let mut out = Vec::new();
    out.try_reserve_exact(estimated)
        .map_err(|_| RigError::AllocationFailed)?;
    out.extend_from_slice(RIG_MAGIC);
    out.extend_from_slice(&RIG_SCHEMA_VERSION.to_le_bytes());
    out.push(BALANCE_MODEL_TAG);
    out.push(ABSOLUTE_INTERVAL_POLICY_TAG);
    out.push(OUTWARD_ROUNDING_POLICY_TAG);
    out.extend_from_slice(&SENSOR_RECORD_SCHEMA_VERSION.to_le_bytes());
    push_text(&mut out, &spec.id);
    push_f64(&mut out, spec.coolant.density_kg_per_m3);
    push_f64(&mut out, spec.coolant.specific_heat_j_per_kg_k);
    push_f64(&mut out, spec.coolant.relative_half_width);
    push_count(&mut out, sensor_rows.len());
    for (role, row) in sensor_rows {
        out.push(role);
        push_bytes(&mut out, &row);
    }
    push_text(&mut out, &run.run_id);
    push_text(&mut out, &run.acquired_on);
    out.push(partition_tag(run.partition));
    push_count(&mut out, readings.len());
    for reading in readings {
        push_text(&mut out, &reading.sensor_id);
        push_f64(&mut out, reading.value_si);
    }
    debug_assert_eq!(out.len(), estimated);
    Ok(out)
}

fn encode_readings(readings: &[Reading]) -> Result<Vec<u8>, RigError> {
    let mut estimated = 4_usize;
    for reading in readings {
        estimated = checked_encoded_add(estimated, 4 + reading.sensor_id.len() + 8)?;
    }
    if estimated > MAX_RETAINED_RIG_BYTES {
        return Err(RigError::RetainedRecordTooLarge {
            have: estimated,
            max: MAX_RETAINED_RIG_BYTES,
        });
    }
    let mut out = Vec::new();
    out.try_reserve_exact(estimated)
        .map_err(|_| RigError::AllocationFailed)?;
    push_count(&mut out, readings.len());
    for reading in readings {
        push_text(&mut out, &reading.sensor_id);
        push_f64(&mut out, reading.value_si);
    }
    Ok(out)
}

fn decode_retained_run(bytes: &[u8]) -> Result<IngestedRun, RigError> {
    if bytes.len() > MAX_RETAINED_RIG_BYTES {
        return Err(RigError::RetainedRecordTooLarge {
            have: bytes.len(),
            max: MAX_RETAINED_RIG_BYTES,
        });
    }
    let mut reader = RigReader::new(bytes);
    if reader.take(RIG_MAGIC.len())? != RIG_MAGIC {
        return Err(RigError::MalformedRetainedRecord {
            reason: "bad FSVVRIG magic",
        });
    }
    let schema = reader.u32()?;
    if schema != RIG_SCHEMA_VERSION {
        return Err(RigError::UnsupportedSchema { observed: schema });
    }
    for (observed, expected, reason) in [
        (
            reader.u8()?,
            BALANCE_MODEL_TAG,
            "unsupported balance model tag",
        ),
        (
            reader.u8()?,
            ABSOLUTE_INTERVAL_POLICY_TAG,
            "unsupported interval policy tag",
        ),
        (
            reader.u8()?,
            OUTWARD_ROUNDING_POLICY_TAG,
            "unsupported rounding policy tag",
        ),
    ] {
        if observed != expected {
            return Err(RigError::MalformedRetainedRecord { reason });
        }
    }
    if reader.u32()? != SENSOR_RECORD_SCHEMA_VERSION {
        return Err(RigError::MalformedRetainedRecord {
            reason: "unsupported nested sensor schema",
        });
    }

    let rig_id = reader.text()?;
    let coolant = CoolantProperties::new(reader.f64()?, reader.f64()?, reader.f64()?)?;
    let instrument_count = reader.count(MAX_INSTRUMENTS)?;
    let mut instruments = Vec::new();
    instruments
        .try_reserve_exact(instrument_count)
        .map_err(|_| RigError::AllocationFailed)?;
    for _ in 0..instrument_count {
        let role = parse_role(reader.u8()?)?;
        let sensor_bytes = reader.bytes()?;
        let sensor =
            decode_sensor_record(sensor_bytes).map_err(|error| RigError::InvalidSensor {
                sensor_id: "<retained-sensor>".to_string(),
                detail: error.to_string(),
            })?;
        instruments.push(InstrumentSpec::new(sensor, role)?);
    }
    let spec = RigSpec::new(rig_id, coolant, instruments)?;

    let run_id = reader.text()?;
    let acquired_on = reader.text()?;
    let partition = parse_partition(reader.u8()?)?;
    let reading_count = reader.count(MAX_INSTRUMENTS)?;
    let mut readings = Vec::new();
    readings
        .try_reserve_exact(reading_count)
        .map_err(|_| RigError::AllocationFailed)?;
    for _ in 0..reading_count {
        readings.push(Reading {
            sensor_id: reader.text()?,
            value_si: reader.f64()?,
        });
    }
    if reader.remaining() != 0 {
        return Err(RigError::MalformedRetainedRecord {
            reason: "trailing bytes",
        });
    }

    let run = RigRun {
        run_id,
        acquired_on,
        partition,
        readings,
    };
    let ingested = ingest(&spec, &run)?;
    if ingested.retained_bytes() != bytes {
        return Err(RigError::NonCanonicalRetainedRecord);
    }
    Ok(ingested)
}

fn checked_encoded_add(current: usize, additional: usize) -> Result<usize, RigError> {
    current
        .checked_add(additional)
        .ok_or(RigError::RetainedRecordTooLarge {
            have: usize::MAX,
            max: MAX_RETAINED_RIG_BYTES,
        })
}

fn push_count(out: &mut Vec<u8>, value: usize) {
    out.extend_from_slice(&(value as u32).to_le_bytes());
}

fn push_text(out: &mut Vec<u8>, value: &str) {
    push_count(out, value.len());
    out.extend_from_slice(value.as_bytes());
}

fn push_bytes(out: &mut Vec<u8>, value: &[u8]) {
    push_count(out, value.len());
    out.extend_from_slice(value);
}

fn push_f64(out: &mut Vec<u8>, value: f64) {
    out.extend_from_slice(&value.to_bits().to_le_bytes());
}

const fn role_tag(role: ChannelRole) -> u8 {
    match role {
        ChannelRole::HeaterPower => 1,
        ChannelRole::InletTemperature => 2,
        ChannelRole::OutletTemperature => 3,
        ChannelRole::VolumeFlow => 4,
        ChannelRole::Ambient => 5,
        ChannelRole::Auxiliary => 6,
    }
}

fn parse_role(tag: u8) -> Result<ChannelRole, RigError> {
    match tag {
        1 => Ok(ChannelRole::HeaterPower),
        2 => Ok(ChannelRole::InletTemperature),
        3 => Ok(ChannelRole::OutletTemperature),
        4 => Ok(ChannelRole::VolumeFlow),
        5 => Ok(ChannelRole::Ambient),
        6 => Ok(ChannelRole::Auxiliary),
        _ => Err(RigError::MalformedRetainedRecord {
            reason: "invalid channel-role tag",
        }),
    }
}

const fn partition_tag(partition: DatasetPartition) -> u8 {
    match partition {
        DatasetPartition::Training => 1,
        DatasetPartition::Calibration => 2,
        DatasetPartition::Validation => 3,
        DatasetPartition::BlindHoldout => 4,
    }
}

fn parse_partition(tag: u8) -> Result<DatasetPartition, RigError> {
    match tag {
        1 => Ok(DatasetPartition::Training),
        2 => Ok(DatasetPartition::Calibration),
        3 => Ok(DatasetPartition::Validation),
        4 => Ok(DatasetPartition::BlindHoldout),
        _ => Err(RigError::MalformedRetainedRecord {
            reason: "invalid dataset-partition tag",
        }),
    }
}

fn diagnostic_sensor_id(value: &str) -> String {
    if value.len() <= 128 {
        return value.to_string();
    }
    let prefix: String = value.chars().take(64).collect();
    format!("{prefix}...({} bytes)", value.len())
}

fn next_up(value: f64) -> f64 {
    if value.is_nan() || value == f64::INFINITY {
        return value;
    }
    if value == 0.0 {
        return f64::from_bits(1);
    }
    let bits = value.to_bits();
    f64::from_bits(if value > 0.0 { bits + 1 } else { bits - 1 })
}

fn next_down(value: f64) -> f64 {
    if value.is_nan() || value == f64::NEG_INFINITY {
        return value;
    }
    if value == 0.0 {
        return -f64::from_bits(1);
    }
    let bits = value.to_bits();
    f64::from_bits(if value > 0.0 { bits - 1 } else { bits + 1 })
}

struct RigReader<'a> {
    bytes: &'a [u8],
    cursor: usize,
}

impl<'a> RigReader<'a> {
    const fn new(bytes: &'a [u8]) -> Self {
        Self { bytes, cursor: 0 }
    }

    fn remaining(&self) -> usize {
        self.bytes.len().saturating_sub(self.cursor)
    }

    fn take(&mut self, count: usize) -> Result<&'a [u8], RigError> {
        let end = self
            .cursor
            .checked_add(count)
            .ok_or(RigError::MalformedRetainedRecord {
                reason: "length overflow",
            })?;
        let value = self
            .bytes
            .get(self.cursor..end)
            .ok_or(RigError::MalformedRetainedRecord {
                reason: "truncated bytes",
            })?;
        self.cursor = end;
        Ok(value)
    }

    fn u8(&mut self) -> Result<u8, RigError> {
        Ok(self.take(1)?[0])
    }

    fn u32(&mut self) -> Result<u32, RigError> {
        let bytes: [u8; 4] =
            self.take(4)?
                .try_into()
                .map_err(|_| RigError::MalformedRetainedRecord {
                    reason: "truncated u32",
                })?;
        Ok(u32::from_le_bytes(bytes))
    }

    fn f64(&mut self) -> Result<f64, RigError> {
        let bytes: [u8; 8] =
            self.take(8)?
                .try_into()
                .map_err(|_| RigError::MalformedRetainedRecord {
                    reason: "truncated f64",
                })?;
        Ok(f64::from_bits(u64::from_le_bytes(bytes)))
    }

    fn count(&mut self, max: usize) -> Result<usize, RigError> {
        let count = self.u32()? as usize;
        if count > max {
            return Err(RigError::MalformedRetainedRecord {
                reason: "collection count exceeds its cap",
            });
        }
        Ok(count)
    }

    fn text(&mut self) -> Result<String, RigError> {
        let count = self.count(MAX_CORPUS_TEXT_BYTES)?;
        let bytes = self.take(count)?;
        let text = core::str::from_utf8(bytes).map_err(|_| RigError::MalformedRetainedRecord {
            reason: "invalid UTF-8",
        })?;
        Ok(text.to_string())
    }

    fn bytes(&mut self) -> Result<&'a [u8], RigError> {
        let count = self.u32()? as usize;
        if count > MAX_RETAINED_RIG_BYTES {
            return Err(RigError::MalformedRetainedRecord {
                reason: "nested byte row exceeds its cap",
            });
        }
        self.take(count)
    }
}

/// Dimensions a pressure channel must declare, for callers assembling an
/// auxiliary instrument set.
#[must_use]
pub const fn pressure_dims() -> Dims {
    PRESSURE_DIMS
}

#[cfg(test)]
mod tests {
    use super::{
        ChannelRole, CoolantProperties, InstrumentSpec, POWER_DIMS, Reading, RigError, RigRun,
        RigSpec, TEMPERATURE_DIMS, VOLUME_FLOW_DIMS, coolant_property_interval,
        decode_retained_run, encode_retained_run,
    };
    use crate::corpus::{
        Availability, CalibrationRecord, DatasetPartition, MeasurementUncertainty, SensorRecord,
    };
    use fs_qty::QtyAny;

    #[test]
    fn coolant_product_band_encloses_exact_pre_rounding_upper_endpoint() {
        // Exact-rational falsifier for the tempting `fl(rho*cp) +/-
        // fl(fl(rho*cp)*r)` construction. Its old upper endpoint was
        // 0x3ffe_6666_6666_666d, but the exact real rho*cp*(1+r) lies above
        // that value. 0x...66e is the least binary64 upper bound.
        let one = 1.0_f64.to_bits();
        let coolant = CoolantProperties::new(f64::from_bits(one - 1), f64::from_bits(one + 4), 0.9)
            .expect("counterexample declarations are finite and positive");
        let interval = coolant_property_interval(coolant).expect("outward interval exists");
        let exact_upper_ceil = f64::from_bits(0x3ffe_6666_6666_666e);
        assert!(
            interval.hi() >= exact_upper_ceil,
            "upper endpoint {:?} under-encloses exact rho*cp*(1+r) {:?}",
            interval.hi(),
            exact_upper_ceil
        );
    }

    #[test]
    fn retained_decoder_refuses_semantically_valid_noncanonical_reading_order() {
        fn sensor(id: &str, dims: fs_qty::Dims) -> SensorRecord {
            SensorRecord {
                id: id.to_string(),
                instrument_id: Availability::Available(format!("instrument-{id}")),
                raw_channel: format!("ch_{id}"),
                quantity_dims: dims,
                calibration: Availability::Available(CalibrationRecord {
                    certificate_id: format!("certificate-{id}"),
                    certificate_hash: fs_blake3::hash_domain(
                        "test.rig.noncanonical.calibration",
                        id.as_bytes(),
                    ),
                    issued_on: "2026-01-01".to_string(),
                    valid_through: None,
                }),
                placement: Availability::Unavailable {
                    reason: "synthetic noncanonical-codec fixture".to_string(),
                },
                uncertainty: MeasurementUncertainty::Bounded {
                    half_width: QtyAny::new(0.0, dims),
                },
            }
        }

        let instruments = [
            ("power", POWER_DIMS, ChannelRole::HeaterPower),
            ("t-in", TEMPERATURE_DIMS, ChannelRole::InletTemperature),
            ("t-out", TEMPERATURE_DIMS, ChannelRole::OutletTemperature),
            ("flow", VOLUME_FLOW_DIMS, ChannelRole::VolumeFlow),
        ]
        .into_iter()
        .map(|(id, dims, role)| {
            InstrumentSpec::new(sensor(id, dims), role).expect("valid synthetic instrument")
        })
        .collect();
        let spec = RigSpec::new(
            "codec-rig",
            CoolantProperties::new(1.0, 1.0, 0.0).expect("valid coolant"),
            instruments,
        )
        .expect("valid rig");
        let run = RigRun {
            run_id: "codec-run".to_string(),
            acquired_on: "2026-06-01".to_string(),
            partition: DatasetPartition::Validation,
            // Deliberately valid but not canonical sensor-id order.
            readings: vec![
                Reading {
                    sensor_id: "power".to_string(),
                    value_si: 1.0,
                },
                Reading {
                    sensor_id: "t-in".to_string(),
                    value_si: 1.0,
                },
                Reading {
                    sensor_id: "t-out".to_string(),
                    value_si: 2.0,
                },
                Reading {
                    sensor_id: "flow".to_string(),
                    value_si: 1.0,
                },
            ],
        };
        let noncanonical =
            encode_retained_run(&spec, &run, &run.readings).expect("structurally valid frame");
        assert_eq!(
            decode_retained_run(&noncanonical)
                .expect_err("semantic re-ingest must expose noncanonical order"),
            RigError::NonCanonicalRetainedRecord
        );
    }
}
