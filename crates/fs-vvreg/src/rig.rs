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
//! certificates, or whose measured energy does not balance. Each refusal
//! names the sensor or the numbers, because the point of ingesting through a
//! gate is to learn what is wrong with the data before it becomes evidence.
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
//! It is adjudicated against the instruments' OWN stated uncertainty, summed
//! conservatively. That has a deliberate consequence: a run whose channels
//! declare [`MeasurementUncertainty::Unstated`] cannot pass, because there is
//! no band to adjudicate against. That is a feature. A Level-E dataset exists
//! precisely to carry metrology, and admitting an unquantified run would
//! produce exactly the unfalsifiable record this corpus is meant to exclude.
//!
//! # Blind holdout
//!
//! A [`DatasetPartition::BlindHoldout`] run is SEALED after ingest: the gate
//! still runs, so a blind run is known to be physically sane, but its
//! measured values do not survive into the retained verdict. Sealing after
//! checking is the ordering that gives both properties; checking after
//! sealing would give neither.

use crate::corpus::{Availability, DatasetPartition, MeasurementUncertainty, SensorRecord};
use fs_blake3::ContentHash;
use fs_qty::Dims;

/// SI exponents, `[m, kg, s, K, …]`.
const TEMPERATURE_DIMS: Dims = Dims([0, 0, 0, 1, 0, 0]);
const POWER_DIMS: Dims = Dims([2, 1, -3, 0, 0, 0]);
const VOLUME_FLOW_DIMS: Dims = Dims([3, 0, -1, 0, 0, 0]);
const PRESSURE_DIMS: Dims = Dims([-1, 1, -2, 0, 0, 0]);

const RIG_IDENTITY_DOMAIN: &str = "org.frankensim.fs-vvreg.rig-run.v1";

/// Schema version of the retained rig-run record.
pub const RIG_SCHEMA_VERSION: u32 = 1;

/// Maximum instruments admitted on one rig.
pub const MAX_INSTRUMENTS: usize = 1_024;

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
    /// More instruments than one rig admits.
    TooManyInstruments {
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
    /// A capacity reservation failed.
    AllocationFailed,
}

impl core::fmt::Display for RigError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::BlankIdentifier { field } => write!(f, "{field} is blank"),
            Self::InvalidScalar { field, requirement } => {
                write!(f, "{field} must be {requirement}")
            }
            Self::TooManyInstruments { have, max } => {
                write!(f, "{have} instruments exceeds the admitted maximum {max}")
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
    /// A non-positive density or specific heat, or a negative/non-finite
    /// relative half-width.
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
        if !relative_half_width.is_finite() || relative_half_width < 0.0 {
            return Err(RigError::InvalidScalar {
                field: "coolant relative half-width",
                requirement: "finite and non-negative",
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
    /// A blank identifier, or declared dimensions incompatible with the role.
    pub fn new(sensor: SensorRecord, role: ChannelRole) -> Result<Self, RigError> {
        if sensor.id.trim().is_empty() {
            return Err(RigError::BlankIdentifier { field: "sensor id" });
        }
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

    /// Stated one-sided half-width in SI units, when the channel declares a
    /// bounded one.
    ///
    /// `CovarianceDiagonal` is converted by taking the square root of the
    /// variance — a one-sigma half-width, stated as such rather than silently
    /// treated as a coverage interval.
    #[must_use]
    pub fn half_width_si(&self) -> Option<f64> {
        match &self.sensor.uncertainty {
            MeasurementUncertainty::Bounded { half_width } => Some(half_width.value.abs()),
            MeasurementUncertainty::CovarianceDiagonal { variance } => {
                Some(variance.value.abs().sqrt())
            }
            MeasurementUncertainty::Unstated => None,
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
    /// Declared partition. Declared AT INGEST, never inferred afterwards.
    pub partition: DatasetPartition,
    /// Channel readings.
    pub readings: Vec<Reading>,
}

/// The adjudicated energy-closure verdict retained with an ingested run.
#[derive(Debug, Clone, PartialEq)]
pub enum EnergyBalance {
    /// The closure holds within the instruments' combined stated band.
    Balanced {
        /// Electrical power in, W.
        power_in_w: f64,
        /// Enthalpy rise carried by the coolant, W.
        heat_out_w: f64,
        /// Signed `power_in − heat_out`, W.
        residual_w: f64,
        /// Combined conservative half-width the residual was judged against.
        allowed_w: f64,
    },
    /// The run passed the gate, but the numbers are withheld because the run
    /// is a sealed blind holdout.
    Sealed,
}

/// An admitted rig run.
#[derive(Debug, Clone, PartialEq)]
pub struct IngestedRun {
    run_id: String,
    rig_id: String,
    partition: DatasetPartition,
    sealed: bool,
    balance: EnergyBalance,
    calibrated_channels: usize,
    identity: ContentHash,
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
    if run.run_id.trim().is_empty() {
        return Err(RigError::BlankIdentifier { field: "run id" });
    }

    let mut sorted: Vec<&Reading> = run.readings.iter().collect();
    sorted.sort_by(|a, b| a.sensor_id.cmp(&b.sensor_id));
    for pair in sorted.windows(2) {
        if pair[0].sensor_id == pair[1].sensor_id {
            return Err(RigError::DuplicateReading {
                sensor_id: pair[0].sensor_id.clone(),
            });
        }
    }
    for reading in &sorted {
        if !reading.value_si.is_finite() {
            return Err(RigError::InvalidScalar {
                field: "reading value",
                requirement: "finite",
            });
        }
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
    for instrument in &spec.instruments {
        if sorted
            .binary_search_by(|reading| reading.sensor_id.cmp(&instrument.sensor.id))
            .is_err()
        {
            return Err(RigError::MissingReading {
                sensor_id: instrument.sensor.id.clone(),
            });
        }
        // RigSpec::new already refuses an uncalibrated instrument; re-checking
        // here keeps ingest self-standing rather than relying on a
        // precondition a future caller could route around.
        if let Availability::Unavailable { reason } = &instrument.sensor.calibration {
            return Err(RigError::Uncalibrated {
                sensor_id: instrument.sensor.id.clone(),
                reason: reason.clone(),
            });
        }
    }

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
    let band_of = |role: ChannelRole| -> Option<f64> {
        spec.role_instrument(role)
            .expect("RigSpec::new admitted every balance role")
            .half_width_si()
    };

    let power_in = value_of(ChannelRole::HeaterPower);
    let inlet = value_of(ChannelRole::InletTemperature);
    let outlet = value_of(ChannelRole::OutletTemperature);
    let flow = value_of(ChannelRole::VolumeFlow);

    let coolant = spec.coolant;
    let delta_t = outlet - inlet;
    let heat_out = coolant.density_kg_per_m3 * coolant.specific_heat_j_per_kg_k * flow * delta_t;

    // Conservative first-order propagation: relative half-widths of a product
    // add, and the temperature difference's band is the sum of its endpoints'.
    // Deliberately not a quadrature sum — an rss here would report a tighter
    // band than the declarations support.
    let (Some(power_band), Some(inlet_band), Some(outlet_band), Some(flow_band)) = (
        band_of(ChannelRole::HeaterPower),
        band_of(ChannelRole::InletTemperature),
        band_of(ChannelRole::OutletTemperature),
        band_of(ChannelRole::VolumeFlow),
    ) else {
        return Err(RigError::EnergyBalanceUnadjudicated {
            reason: "a balance channel declares MeasurementUncertainty::Unstated, so there is no band to judge the closure against; Level E exists to carry metrology".to_string(),
        });
    };

    let delta_band = inlet_band + outlet_band;
    if !(delta_t.abs() > 0.0) {
        return Err(RigError::EnergyBalanceUnadjudicated {
            reason:
                "the measured coolant temperature rise is zero, so no enthalpy flow can be resolved"
                    .to_string(),
        });
    }
    let relative =
        flow_band / flow.abs() + delta_band / delta_t.abs() + coolant.relative_half_width;
    let heat_band = heat_out.abs() * relative;
    let allowed = power_band + heat_band;
    let residual = power_in - heat_out;

    for (value, field) in [
        (heat_out, "heat out"),
        (allowed, "allowed band"),
        (residual, "residual"),
    ] {
        if !value.is_finite() {
            return Err(RigError::InvalidScalar {
                field: match field {
                    "heat out" => "heat out",
                    "allowed band" => "allowed band",
                    _ => "residual",
                },
                requirement: "finite",
            });
        }
    }

    if residual.abs() > allowed {
        return Err(RigError::EnergyBalanceViolated {
            detail: render_balance(spec, power_in, heat_out, residual, allowed, flow, delta_t),
        });
    }

    let sealed = run.partition == DatasetPartition::BlindHoldout;
    let balance = if sealed {
        EnergyBalance::Sealed
    } else {
        EnergyBalance::Balanced {
            power_in_w: power_in,
            heat_out_w: heat_out,
            residual_w: residual,
            allowed_w: allowed,
        }
    };
    let identity = run_identity(spec, run, power_in, heat_out, residual, allowed);

    Ok(IngestedRun {
        run_id: run.run_id.clone(),
        rig_id: spec.id.clone(),
        partition: run.partition,
        sealed,
        balance,
        calibrated_channels: spec.instruments.len(),
        identity,
    })
}

fn render_balance(
    spec: &RigSpec,
    power_in: f64,
    heat_out: f64,
    residual: f64,
    allowed: f64,
    flow: f64,
    delta_t: f64,
) -> String {
    use core::fmt::Write as _;
    let mut out = String::new();
    let _ = writeln!(
        out,
        "  rig {}: power in {power_in} W, heat out {heat_out} W",
        spec.id
    );
    let _ = writeln!(
        out,
        "  residual {residual} W exceeds the combined stated band {allowed} W"
    );
    let _ = writeln!(
        out,
        "  flow {flow} m^3/s, coolant rise {delta_t} K, rho {} kg/m^3, c_p {} J/(kg K)",
        spec.coolant.density_kg_per_m3, spec.coolant.specific_heat_j_per_kg_k
    );
    for instrument in &spec.instruments {
        let band = instrument
            .half_width_si()
            .map_or_else(|| "unstated".to_string(), |value| format!("{value}"));
        let _ = writeln!(
            out,
            "  channel {} [{}]: half-width {band}",
            instrument.sensor.id,
            instrument.role.name()
        );
    }
    out
}

fn run_identity(
    spec: &RigSpec,
    run: &RigRun,
    power_in: f64,
    heat_out: f64,
    residual: f64,
    allowed: f64,
) -> ContentHash {
    let mut hasher = fs_blake3::Blake3::new();
    let mut field = |bytes: &[u8]| {
        hasher.update(&(bytes.len() as u64).to_le_bytes());
        hasher.update(bytes);
    };
    field(RIG_IDENTITY_DOMAIN.as_bytes());
    field(&RIG_SCHEMA_VERSION.to_le_bytes());
    field(spec.id.as_bytes());
    field(run.run_id.as_bytes());
    field(run.partition.name().as_bytes());
    for instrument in &spec.instruments {
        field(instrument.sensor.id.as_bytes());
        field(instrument.role.name().as_bytes());
    }
    for value in [power_in, heat_out, residual, allowed] {
        let canonical = if value == 0.0 { 0.0 } else { value };
        field(&canonical.to_bits().to_le_bytes());
    }
    let preimage = hasher.finalize();
    fs_blake3::hash_domain(RIG_IDENTITY_DOMAIN, preimage.as_bytes())
}

/// Dimensions a pressure channel must declare, for callers assembling an
/// auxiliary instrument set.
#[must_use]
pub const fn pressure_dims() -> Dims {
    PRESSURE_DIMS
}
