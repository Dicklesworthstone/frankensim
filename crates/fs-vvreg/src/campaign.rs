//! Level-E physical cooling campaign and operating matrix schema
//! (bead `frankensim-extreal-program-f85xj.4.5.4`).
//!
//! Models the physical / synthetic Level-E campaign execution:
//! - Operating matrix with adversarial vent, contact (TIM), and fan variations
//! - Mandatory safety checklist and overtemperature boundaries
//! - Retention and partition declarations (with at least one candidate withheld/blind holdout)
//! - Explicit [`CampaignDisposition::NoData`] until physical acquisition occurs
//! - Deterministic content-addressed manifest hashing

use crate::corpus::{
    Availability, DatasetPartition, MAX_CORPUS_TEXT_BYTES, valid_date, valid_slug,
};
use crate::rig::{RigError, RigRun, RigSpec};
use fs_blake3::{ContentHash, hash_bytes};

/// Schema version for the Level-E campaign manifest.
pub const LEVEL_E_CAMPAIGN_SCHEMA_VERSION: u32 = 1;

/// Domain separator for Level-E campaign manifest content hashing.
pub const CAMPAIGN_MANIFEST_DOMAIN: &str = "org.frankensim.fs-vvreg.level-e-campaign.v1";

/// Maximum points in one operating matrix.
pub const MAX_MATRIX_POINTS: usize = 256;

/// Controllable enclosure vent configuration.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum VentConfig {
    /// Fully unobstructed vent path.
    NominalOpen,
    /// 50% flow restriction on outlet.
    Restricted50,
    /// Fully blocked / closed vent (adversarial overheating condition).
    Blocked100,
}

impl VentConfig {
    /// Canonical slug representation.
    #[must_use]
    pub const fn slug(self) -> &'static str {
        match self {
            Self::NominalOpen => "nominal-open",
            Self::Restricted50 => "restricted-50",
            Self::Blocked100 => "blocked-100",
        }
    }
}

/// Controllable fan operating state.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum FanState {
    /// Full rated RPM.
    Nominal100,
    /// Reduced voltage / 50% speed.
    Reduced50,
    /// Locked rotor / stalled fan (adversarial loss of forced convection).
    Stalled,
}

impl FanState {
    /// Canonical slug representation.
    #[must_use]
    pub const fn slug(self) -> &'static str {
        match self {
            Self::Nominal100 => "nominal-100",
            Self::Reduced50 => "reduced-50",
            Self::Stalled => "stalled",
        }
    }
}

/// Thermal Interface Material (TIM) contact condition between heater and plate.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum TimCondition {
    /// Standard thermal paste (low thermal resistance).
    StandardGrease,
    /// High thermal resistance pad (degraded contact).
    DegradedPad,
    /// Dry metal-to-metal contact with air gap (adversarial poor contact).
    DryContact,
}

impl TimCondition {
    /// Canonical slug representation.
    #[must_use]
    pub const fn slug(self) -> &'static str {
        match self {
            Self::StandardGrease => "standard-grease",
            Self::DegradedPad => "degraded-pad",
            Self::DryContact => "dry-contact",
        }
    }
}

/// One planned operating condition in the campaign matrix.
#[derive(Debug, Clone, PartialEq)]
pub struct OperatingMatrixPoint {
    /// Unique point identifier (canonical slug).
    pub point_id: String,
    /// Target heater power in Watts (must be positive).
    pub heater_power_w: f64,
    /// Target volumetric flow in m³/s (must be positive).
    pub flow_rate_m3_s: f64,
    /// Vent configuration.
    pub vent: VentConfig,
    /// Fan state.
    pub fan: FanState,
    /// Thermal interface condition.
    pub tim: TimCondition,
    /// Target dataset partition.
    pub partition: DatasetPartition,
}

/// Pre-flight safety and readiness checklist signed by human operator.
#[derive(Debug, Clone, PartialEq)]
pub struct SafetyChecklist {
    /// Operator identifier (canonical slug).
    pub operator_id: String,
    /// Overtemperature safety cutoff in Kelvin (e.g. 363.15 K = 90 °C).
    pub overtemp_cutoff_k: f64,
    /// Hardware emergency power interrupt verified.
    pub emergency_stop_verified: bool,
    /// Sensor wiring and polarity verified.
    pub wiring_verified: bool,
    /// All calibration certificates in date.
    pub calibrations_verified: bool,
    /// Checklist sign-off date (canonical ISO YYYY-MM-DD).
    pub signoff_date: String,
}

/// Disposition of the physical campaign.
#[derive(Debug, Clone, PartialEq)]
pub enum CampaignDisposition {
    /// Explicit placeholder prior to physical hardware execution.
    NoData {
        /// Documented procurement / assembly reason.
        reason: String,
    },
    /// Synthetic trial run executed to validate software ingestion.
    SyntheticTrial {
        /// Admitted synthetic runs.
        retained_runs: Vec<RigRun>,
    },
    /// Actual physical execution with owned raw data.
    PhysicalAcquired {
        /// Admitted physical runs.
        retained_runs: Vec<RigRun>,
    },
}

/// Comprehensive Level-E campaign manifest.
#[derive(Debug, Clone, PartialEq)]
pub struct LevelECampaignManifest {
    /// Campaign identifier (canonical slug).
    pub campaign_id: String,
    /// Rig specification identifier.
    pub rig_id: String,
    /// Planned operating matrix points.
    pub matrix: Vec<OperatingMatrixPoint>,
    /// Signed pre-flight safety checklist.
    pub safety: SafetyChecklist,
    /// Campaign execution disposition.
    pub disposition: CampaignDisposition,
}

impl LevelECampaignManifest {
    /// Validate all structural and semantic invariants of the campaign manifest.
    ///
    /// # Errors
    /// Returns [`RigError`] on invalid slug, non-positive scalar, empty matrix,
    /// unverified safety check, or missing blind-holdout partition.
    pub fn validate(&self) -> Result<(), RigError> {
        if !valid_slug(&self.campaign_id) {
            return Err(RigError::InvalidIdentifier {
                field: "campaign_id",
                requirement: "campaign_id must be a canonical lowercase ASCII slug (<= 64 bytes)",
            });
        }
        if !valid_slug(&self.rig_id) {
            return Err(RigError::InvalidIdentifier {
                field: "rig_id",
                requirement: "rig_id must be a canonical lowercase ASCII slug (<= 64 bytes)",
            });
        }
        if self.matrix.is_empty() {
            return Err(RigError::InvalidScalar {
                field: "matrix",
                requirement: "operating matrix must not be empty",
            });
        }
        if self.matrix.len() > MAX_MATRIX_POINTS {
            return Err(RigError::TooManyInstruments {
                have: self.matrix.len(),
                max: MAX_MATRIX_POINTS,
            });
        }

        let mut has_blind_holdout = false;
        for pt in &self.matrix {
            if !valid_slug(&pt.point_id) {
                return Err(RigError::InvalidIdentifier {
                    field: "matrix.point_id",
                    requirement: "point_id must be a canonical lowercase ASCII slug",
                });
            }
            if !pt.heater_power_w.is_finite() || pt.heater_power_w <= 0.0 {
                return Err(RigError::InvalidScalar {
                    field: "matrix.heater_power_w",
                    requirement: "heater power target must be strictly positive and finite",
                });
            }
            if !pt.flow_rate_m3_s.is_finite() || pt.flow_rate_m3_s <= 0.0 {
                return Err(RigError::InvalidScalar {
                    field: "matrix.flow_rate_m3_s",
                    requirement: "flow rate target must be strictly positive and finite",
                });
            }
            if pt.partition == DatasetPartition::BlindHoldout {
                has_blind_holdout = true;
            }
        }

        if !has_blind_holdout {
            return Err(RigError::InvalidScalar {
                field: "matrix.partition",
                requirement: "operating matrix must declare at least one BlindHoldout partition",
            });
        }

        if !valid_slug(&self.safety.operator_id) {
            return Err(RigError::InvalidIdentifier {
                field: "safety.operator_id",
                requirement: "operator_id must be a canonical lowercase ASCII slug",
            });
        }
        if !self.safety.overtemp_cutoff_k.is_finite() || self.safety.overtemp_cutoff_k <= 273.15 {
            return Err(RigError::InvalidScalar {
                field: "safety.overtemp_cutoff_k",
                requirement: "overtemperature cutoff must be finite and above 273.15 K (0 °C)",
            });
        }
        if !valid_date(&self.safety.signoff_date) {
            return Err(RigError::InvalidDate {
                field: "safety.signoff_date",
                value: self.safety.signoff_date.clone(),
            });
        }
        if !self.safety.emergency_stop_verified
            || !self.safety.wiring_verified
            || !self.safety.calibrations_verified
        {
            return Err(RigError::InvalidScalar {
                field: "safety",
                requirement: "all safety checklist items must be verified before execution",
            });
        }

        Ok(())
    }

    /// Compute the deterministic content hash for this campaign manifest.
    #[must_use]
    pub fn content_hash(&self) -> ContentHash {
        let mut buf = Vec::new();
        buf.extend_from_slice(CAMPAIGN_MANIFEST_DOMAIN.as_bytes());
        buf.extend_from_slice(&LEVEL_E_CAMPAIGN_SCHEMA_VERSION.to_le_bytes());
        buf.extend_from_slice(self.campaign_id.as_bytes());
        buf.extend_from_slice(self.rig_id.as_bytes());
        for pt in &self.matrix {
            buf.extend_from_slice(pt.point_id.as_bytes());
            buf.extend_from_slice(&pt.heater_power_w.to_bits().to_le_bytes());
            buf.extend_from_slice(&pt.flow_rate_m3_s.to_bits().to_le_bytes());
            buf.extend_from_slice(pt.vent.slug().as_bytes());
            buf.extend_from_slice(pt.fan.slug().as_bytes());
            buf.extend_from_slice(pt.tim.slug().as_bytes());
            buf.extend_from_slice(pt.partition.name().as_bytes());
        }
        buf.extend_from_slice(self.safety.operator_id.as_bytes());
        buf.extend_from_slice(&self.safety.overtemp_cutoff_k.to_bits().to_le_bytes());
        buf.push(u8::from(self.safety.emergency_stop_verified));
        buf.push(u8::from(self.safety.wiring_verified));
        buf.push(u8::from(self.safety.calibrations_verified));
        buf.extend_from_slice(self.safety.signoff_date.as_bytes());
        match &self.disposition {
            CampaignDisposition::NoData { reason } => {
                buf.push(0);
                buf.extend_from_slice(reason.as_bytes());
            }
            CampaignDisposition::SyntheticTrial { retained_runs } => {
                buf.push(1);
                buf.extend_from_slice(&(retained_runs.len() as u32).to_le_bytes());
                for r in retained_runs {
                    buf.extend_from_slice(r.run_id.as_bytes());
                }
            }
            CampaignDisposition::PhysicalAcquired { retained_runs } => {
                buf.push(2);
                buf.extend_from_slice(&(retained_runs.len() as u32).to_le_bytes());
                for r in retained_runs {
                    buf.extend_from_slice(r.run_id.as_bytes());
                }
            }
        }
        hash_bytes(&buf)
    }
}
