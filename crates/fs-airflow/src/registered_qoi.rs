//! Registered production extraction of thermal Quantities of Interest (QoIs)
//! from real solved fields (bead `frankensim-s2l9v.1`).
//!
//! Provides:
//! - Canonical QoI semantic identifiers, query routing, and supported output mapping;
//! - Strict distinction between scalar QoIs, fields, and reports (non-scalar kinds fail closed);
//! - Admitted work and memory budget pre-flight planning;
//! - Deterministic traversal, canonical ordering, and stable tie-breaking for extrema;
//! - Strict SI unit binding and physical bounds enforcement (non-negative absolute temperatures);
//! - Cooperative cancellation polling at bounded work tiles;
//! - Immutable candidate contribution rows with content-addressed provenance.

use core::fmt;
use std::collections::BTreeSet;

use fs_blake3::{ContentHash, hash_domain};
use fs_conduction::{ConductionMesh, ConductionSolution};
use fs_evidence::uncertainty::EngineeringUncertaintyBudget;
use fs_evidence::{ModelEvidence, ProvenanceHash};
use fs_exec::Cx;

use crate::OperatingPoint;
use crate::qoi::{QoiError, ThermalQoiDeclarations, ThermalQoiKind, extract_thermal_qois};

const CANDIDATE_ROW_DOMAIN: &str = "org.frankensim.fs-airflow.candidate-qoi-row.v1";

/// Canonical semantic family of a registered thermal Quantity of Interest.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum QoiSemanticId {
    /// Maximum temperature across a designated component/junction region.
    JunctionMaximum,
    /// Area-weighted mean temperature across a surface boundary.
    SurfaceMeanTemperature,
    /// Temperature spread (max - min) across a surface boundary.
    SurfaceTemperatureSpread,
    /// Surface face temperature standard deviation.
    SurfaceTemperatureStdDev,
    /// Total enclosure pressure drop across the airflow path.
    PressureDrop,
    /// Required electrical fan power from airflow operating point.
    FanPower,
    /// Thermal margin against requirement (limit - max_temp).
    ThermalMargin,
}

impl QoiSemanticId {
    /// All canonical semantic identifiers in deterministic registry order.
    pub const ALL: [Self; 7] = [
        Self::JunctionMaximum,
        Self::SurfaceMeanTemperature,
        Self::SurfaceTemperatureSpread,
        Self::SurfaceTemperatureStdDev,
        Self::PressureDrop,
        Self::FanPower,
        Self::ThermalMargin,
    ];

    /// Canonical string identifier.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::JunctionMaximum => "thermal.junction_maximum",
            Self::SurfaceMeanTemperature => "thermal.surface_mean",
            Self::SurfaceTemperatureSpread => "thermal.surface_spread",
            Self::SurfaceTemperatureStdDev => "thermal.surface_std_dev",
            Self::PressureDrop => "airflow.pressure_drop",
            Self::FanPower => "airflow.fan_power",
            Self::ThermalMargin => "thermal.thermal_margin",
        }
    }

    /// Primary SI units for this QoI family.
    #[must_use]
    pub const fn units(self) -> &'static str {
        match self {
            Self::JunctionMaximum | Self::SurfaceMeanTemperature => "kelvin",
            Self::SurfaceTemperatureSpread | Self::SurfaceTemperatureStdDev => "kelvin",
            Self::PressureDrop => "pascal",
            Self::FanPower => "watt",
            Self::ThermalMargin => "kelvin",
        }
    }

    /// Physical QoI classification kind.
    #[must_use]
    pub const fn qoi_kind(self) -> ThermalQoiKind {
        match self {
            Self::JunctionMaximum | Self::SurfaceMeanTemperature => {
                ThermalQoiKind::AbsoluteTemperature
            }
            Self::SurfaceTemperatureSpread
            | Self::SurfaceTemperatureStdDev
            | Self::ThermalMargin => ThermalQoiKind::TemperatureDifference,
            Self::PressureDrop | Self::FanPower => ThermalQoiKind::TemperatureDifference,
        }
    }

    /// Parse from user-facing query name (accepting canonical or common aliases).
    pub fn parse(name: &str) -> Option<Self> {
        let normalized = name.trim().to_lowercase();
        match normalized.as_str() {
            "thermal.junction_maximum"
            | "junction_maximum"
            | "junction_temp"
            | "junction_max"
            | "t_j_max"
            | "max_temperature" => Some(Self::JunctionMaximum),

            "thermal.surface_mean"
            | "surface_mean"
            | "surface_mean_temp"
            | "case_mean_temp"
            | "t_surface_mean" => Some(Self::SurfaceMeanTemperature),

            "thermal.surface_spread"
            | "surface_spread"
            | "surface_temp_spread"
            | "case_spread"
            | "delta_t_surface" => Some(Self::SurfaceTemperatureSpread),

            "thermal.surface_std_dev" | "surface_std_dev" | "surface_temperature_std_dev" => {
                Some(Self::SurfaceTemperatureStdDev)
            }

            "airflow.pressure_drop" | "pressure_drop" | "delta_p" | "system_resistance" => {
                Some(Self::PressureDrop)
            }

            "airflow.fan_power" | "fan_power" | "fan_input_power" | "p_fan" => {
                Some(Self::FanPower)
            }

            "thermal.thermal_margin" | "thermal_margin" | "margin" | "junction_margin" => {
                Some(Self::ThermalMargin)
            }

            _ => None,
        }
    }
}

impl fmt::Display for QoiSemanticId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.as_str())
    }
}

/// Declared output kind from project configuration.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum OutputKind {
    /// Scalar QoI (supported by this extractor).
    Scalar,
    /// Spatially distributed field (produced by solver/field stages, not scalar QoI).
    Field,
    /// Structured or tabular report (produced by report stage).
    Report,
    /// Unrecognized output kind.
    Other(String),
}

impl OutputKind {
    /// Parse from wire string.
    pub fn parse(s: &str) -> Self {
        match s.trim().to_lowercase().as_str() {
            "scalar" | "qoi" => Self::Scalar,
            "field" | "spatial_field" => Self::Field,
            "report" | "summary" => Self::Report,
            other => Self::Other(other.to_string()),
        }
    }
}

/// One requested output to extract.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OutputQuery {
    /// Raw output name from project configuration.
    pub name: String,
    /// Wire kind.
    pub kind: OutputKind,
    /// Optional region name constraint.
    pub region: Option<String>,
}

impl OutputQuery {
    /// Create a scalar query.
    #[must_use]
    pub fn scalar(name: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            kind: OutputKind::Scalar,
            region: None,
        }
    }

    /// Create a scalar query scoped to a region.
    #[must_use]
    pub fn scalar_with_region(name: impl Into<String>, region: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            kind: OutputKind::Scalar,
            region: Some(region.into()),
        }
    }
}

/// Work and memory limits for registered QoI extraction.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct QoiExecutionLimits {
    /// Maximum allowed mesh elements (default 10,000,000).
    pub max_elements: usize,
    /// Maximum allowed mesh vertices (default 5,000,000).
    pub max_vertices: usize,
    /// Maximum allowed output queries in one batch (default 1,000).
    pub max_queries: usize,
    /// Maximum allowed working memory in bytes (default 512 MiB).
    pub max_memory_bytes: usize,
}

impl Default for QoiExecutionLimits {
    fn default() -> Self {
        Self {
            max_elements: 10_000_000,
            max_vertices: 5_000_000,
            max_queries: 1_000,
            max_memory_bytes: 512 * 1024 * 1024,
        }
    }
}

/// One immutable candidate QoI row emitted by registered extraction.
#[derive(Debug, Clone, PartialEq)]
pub struct CandidateQoiRow {
    /// Canonical semantic identifier.
    pub semantic_id: QoiSemanticId,
    /// User query name that requested this QoI.
    pub query_name: String,
    /// Physical kind (absolute temperature vs difference vs other).
    pub kind: ThermalQoiKind,
    /// SI base unit label.
    pub units: &'static str,
    /// Nominal scalar value.
    pub value: f64,
    /// Eight-source engineering uncertainty budget.
    pub uncertainty: EngineeringUncertaintyBudget,
    /// Associated model evidence.
    pub evidence: ModelEvidence,
    /// Content-addressed identity hash of this candidate row.
    pub identity_hash: ContentHash,
    /// Lineage parent digests.
    pub source_lineage: Vec<ContentHash>,
    /// Associated region name, if any.
    pub region_name: Option<String>,
    /// Tied extremum vertex index, if applicable.
    pub tie_witness_vertex: Option<usize>,
}

/// Receipt proving complete registered thermal extraction over a set of queries.
#[derive(Debug, Clone, PartialEq)]
pub struct RegisteredQoiExtractionReceipt {
    /// Extracted candidate rows in canonical sorted order.
    pub rows: Vec<CandidateQoiRow>,
    /// Number of requested queries evaluated.
    pub requested_query_count: usize,
    /// Number of emitted candidate rows.
    pub emitted_qoi_count: usize,
    /// Bounded work items executed.
    pub executed_work_items: usize,
    /// Content-addressed provenance hash.
    pub provenance: ProvenanceHash,
}

/// Structured error from registered thermal QoI extraction.
#[derive(Debug, Clone, PartialEq)]
pub enum RegisteredQoiError {
    /// Requested output name is not recognized in the thermal registry.
    UnsupportedOutputName {
        /// Requested name.
        name: String,
        /// Canonical supported alternatives.
        suggested: Vec<&'static str>,
    },
    /// Requested output is non-scalar (e.g. field or report).
    NonScalarOutputKind {
        /// Output name.
        name: String,
        /// Declared kind.
        kind: String,
        /// Downstream owner responsible for this kind.
        downstream_owner: &'static str,
    },
    /// Duplicate query name in the batch.
    DuplicateQuery {
        /// Duplicate name.
        name: String,
    },
    /// Requested region was not found in declarations.
    RegionNotFound {
        /// Requested name.
        requested: String,
        /// Available region names.
        available: Vec<String>,
    },
    /// Input work exceeds execution limits.
    WorkLimitExceeded {
        /// Quantity name.
        field: &'static str,
        /// Actual items.
        actual: usize,
        /// Maximum limit.
        limit: usize,
    },
    /// Input memory plan exceeds execution limits.
    MemoryLimitExceeded {
        /// Estimated bytes.
        estimated_bytes: usize,
        /// Maximum limit.
        limit_bytes: usize,
    },
    /// Produced scalar value is non-finite.
    NonFiniteScalar {
        /// Semantic ID.
        semantic_id: QoiSemanticId,
        /// Non-finite value.
        value: f64,
    },
    /// Produced absolute temperature is below physical absolute zero.
    NegativeAbsoluteTemperature {
        /// Value in kelvin.
        kelvin: f64,
    },
    /// Computation was cancelled.
    Cancelled,
    /// Underlying thermal QoI error.
    UnderlyingQoi(QoiError),
}

impl fmt::Display for RegisteredQoiError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnsupportedOutputName { name, suggested } => {
                write!(
                    f,
                    "unsupported output name `{name}`; supported names: {}",
                    suggested.join(", ")
                )
            }
            Self::NonScalarOutputKind {
                name,
                kind,
                downstream_owner,
            } => {
                write!(
                    f,
                    "output `{name}` has non-scalar kind `{kind}`; scalar QoI extractor only evaluates scalar quantities (defer to {downstream_owner})"
                )
            }
            Self::DuplicateQuery { name } => {
                write!(f, "duplicate output query name `{name}` in request batch")
            }
            Self::RegionNotFound {
                requested,
                available,
            } => {
                write!(
                    f,
                    "requested region `{requested}` not found; available: {}",
                    available.join(", ")
                )
            }
            Self::WorkLimitExceeded {
                field,
                actual,
                limit,
            } => {
                write!(
                    f,
                    "work limit exceeded for {field}: actual {actual} > limit {limit}"
                )
            }
            Self::MemoryLimitExceeded {
                estimated_bytes,
                limit_bytes,
            } => {
                write!(
                    f,
                    "memory limit exceeded: estimated {estimated_bytes} B > limit {limit_bytes} B"
                )
            }
            Self::NonFiniteScalar { semantic_id, value } => {
                write!(
                    f,
                    "non-finite scalar {value} produced for QoI `{}`",
                    semantic_id.as_str()
                )
            }
            Self::NegativeAbsoluteTemperature { kelvin } => {
                write!(
                    f,
                    "absolute temperature {kelvin} K is below physical absolute zero (0 K)"
                )
            }
            Self::Cancelled => write!(f, "registered thermal QoI extraction was cancelled"),
            Self::UnderlyingQoi(error) => write!(f, "thermal QoI error: {error}"),
        }
    }
}

impl std::error::Error for RegisteredQoiError {}

impl From<QoiError> for RegisteredQoiError {
    fn from(error: QoiError) -> Self {
        Self::UnderlyingQoi(error)
    }
}

/// Execute registered thermal QoI extraction over a set of output queries.
///
/// Pre-flights limits and memory plans, validates query kinds, extracts
/// underlying QoIs via [`extract_thermal_qois`], and binds candidate rows with
/// canonical identities and deterministic sorting.
///
/// # Errors
/// Refuses non-scalar output kinds, unknown query names, duplicate queries,
/// memory/work limit breaches, non-finite values, or uncancelled interruptions.
pub fn extract_registered_qois(
    queries: &[OutputQuery],
    mesh: &ConductionMesh,
    solution: &ConductionSolution,
    operating_point: &OperatingPoint,
    declarations: &ThermalQoiDeclarations<'_>,
    limits: QoiExecutionLimits,
    cx: &Cx<'_>,
) -> Result<RegisteredQoiExtractionReceipt, RegisteredQoiError> {
    if cx.checkpoint().is_err() {
        return Err(RegisteredQoiError::Cancelled);
    }

    // 1. Pre-flight execution limits
    if mesh.element_count() > limits.max_elements {
        return Err(RegisteredQoiError::WorkLimitExceeded {
            field: "element count",
            actual: mesh.element_count(),
            limit: limits.max_elements,
        });
    }
    if mesh.vertex_count() > limits.max_vertices {
        return Err(RegisteredQoiError::WorkLimitExceeded {
            field: "vertex count",
            actual: mesh.vertex_count(),
            limit: limits.max_vertices,
        });
    }
    if queries.len() > limits.max_queries {
        return Err(RegisteredQoiError::WorkLimitExceeded {
            field: "query count",
            actual: queries.len(),
            limit: limits.max_queries,
        });
    }

    let estimated_memory = mesh.vertex_count() * core::mem::size_of::<f64>()
        + mesh.element_count() * 4 * core::mem::size_of::<usize>()
        + queries.len() * 1024;
    if estimated_memory > limits.max_memory_bytes {
        return Err(RegisteredQoiError::MemoryLimitExceeded {
            estimated_bytes: estimated_memory,
            limit_bytes: limits.max_memory_bytes,
        });
    }

    // 2. Validate queries: duplicate check, kind check, name resolution
    let mut seen_names = BTreeSet::new();
    let mut parsed_queries = Vec::with_capacity(queries.len());

    let suggested_names: Vec<&'static str> =
        QoiSemanticId::ALL.iter().map(|id| id.as_str()).collect();

    for query in queries {
        if !seen_names.insert(query.name.clone()) {
            return Err(RegisteredQoiError::DuplicateQuery {
                name: query.name.clone(),
            });
        }

        match &query.kind {
            OutputKind::Field => {
                return Err(RegisteredQoiError::NonScalarOutputKind {
                    name: query.name.clone(),
                    kind: "field".to_string(),
                    downstream_owner: "conduction / field export stage",
                });
            }
            OutputKind::Report => {
                return Err(RegisteredQoiError::NonScalarOutputKind {
                    name: query.name.clone(),
                    kind: "report".to_string(),
                    downstream_owner: "fs-report stage",
                });
            }
            OutputKind::Other(k) => {
                return Err(RegisteredQoiError::NonScalarOutputKind {
                    name: query.name.clone(),
                    kind: k.clone(),
                    downstream_owner: "unknown non-scalar producer",
                });
            }
            OutputKind::Scalar => {}
        }

        let semantic_id = QoiSemanticId::parse(&query.name).ok_or_else(|| {
            RegisteredQoiError::UnsupportedOutputName {
                name: query.name.clone(),
                suggested: suggested_names.clone(),
            }
        })?;

        // If region constraint is specified, check against declared regions
        if let Some(region_name) = &query.region {
            let available = vec![
                declarations.junction_region.name().to_string(),
                declarations.surface_region.name().to_string(),
            ];
            let matched = declarations.junction_region.name() == region_name
                || declarations.surface_region.name() == region_name;
            if !matched {
                return Err(RegisteredQoiError::RegionNotFound {
                    requested: region_name.clone(),
                    available,
                });
            }
        }

        parsed_queries.push((query, semantic_id));
    }

    if cx.checkpoint().is_err() {
        return Err(RegisteredQoiError::Cancelled);
    }

    // 3. Extract core QoI set using canonical extractor
    let qoi_set = extract_thermal_qois(mesh, solution, operating_point, declarations)?;

    if cx.checkpoint().is_err() {
        return Err(RegisteredQoiError::Cancelled);
    }

    // 4. Map queries to candidate rows
    let mut rows = Vec::with_capacity(parsed_queries.len());

    for (query, semantic_id) in parsed_queries {
        let (value, uncertainty, evidence, region_name, tie_vertex) = match semantic_id {
            QoiSemanticId::JunctionMaximum => {
                let jm = &qoi_set.junction_maximum;
                let val = jm.qoi.evidence.value.value();
                let unc = jm.qoi.uncertainty.clone();
                let ev = jm.qoi.evidence.model.clone();
                let reg = Some(declarations.junction_region.name().to_string());
                let v_idx = jm.vertex;
                (val, unc, ev, reg, Some(v_idx))
            }
            QoiSemanticId::SurfaceMeanTemperature => {
                let sm = &qoi_set.uniformity.mean_temperature;
                let val = sm.evidence.value.value();
                let unc = sm.uncertainty.clone();
                let ev = sm.evidence.model.clone();
                let reg = Some(declarations.surface_region.name().to_string());
                (val, unc, ev, reg, None)
            }
            QoiSemanticId::SurfaceTemperatureSpread => {
                let ss = &qoi_set.uniformity.spread;
                let val = ss.evidence.value.value();
                let unc = ss.uncertainty.clone();
                let ev = ss.evidence.model.clone();
                let reg = Some(declarations.surface_region.name().to_string());
                (val, unc, ev, reg, None)
            }
            QoiSemanticId::SurfaceTemperatureStdDev => {
                let sd = &qoi_set.uniformity.face_mean_standard_deviation;
                let val = sd.evidence.value.value();
                let unc = sd.uncertainty.clone();
                let ev = sd.evidence.model.clone();
                let reg = Some(declarations.surface_region.name().to_string());
                (val, unc, ev, reg, None)
            }
            QoiSemanticId::PressureDrop => {
                let pd = &qoi_set.pressure_drop;
                let val = pd.evidence.value.value();
                let unc = pd.uncertainty.clone();
                let ev = pd.evidence.model.clone();
                (val, unc, ev, None, None)
            }
            QoiSemanticId::FanPower => {
                let fp = &qoi_set.fan_power;
                let val = fp.evidence.value.value();
                let unc = fp.uncertainty.clone();
                let ev = fp.evidence.model.clone();
                (val, unc, ev, None, None)
            }
            QoiSemanticId::ThermalMargin => {
                let tm = &qoi_set.thermal_margin;
                let val = tm.evidence.value.value();
                let unc = tm.uncertainty.clone();
                let ev = tm.evidence.model.clone();
                (val, unc, ev, None, None)
            }
        };

        // Finiteness validation
        if !value.is_finite() {
            return Err(RegisteredQoiError::NonFiniteScalar {
                semantic_id,
                value,
            });
        }

        // Absolute zero validation
        if semantic_id.qoi_kind() == ThermalQoiKind::AbsoluteTemperature && value < 0.0 {
            return Err(RegisteredQoiError::NegativeAbsoluteTemperature { kelvin: value });
        }

        // Canonical candidate row identity
        let mut payload = Vec::new();
        payload.extend_from_slice(semantic_id.as_str().as_bytes());
        payload.extend_from_slice(query.name.as_bytes());
        payload.extend_from_slice(semantic_id.units().as_bytes());
        payload.extend_from_slice(&value.to_bits().to_le_bytes());
        payload.extend_from_slice(uncertainty.qoi().as_bytes());
        if let Some(r) = &region_name {
            payload.extend_from_slice(r.as_bytes());
        }
        if let Some(v) = tie_vertex {
            payload.extend_from_slice(&v.to_le_bytes());
        }
        let identity_hash = hash_domain(CANDIDATE_ROW_DOMAIN, &payload);

        let source_lineage = vec![identity_hash];

        rows.push(CandidateQoiRow {
            semantic_id,
            query_name: query.name.clone(),
            kind: semantic_id.qoi_kind(),
            units: semantic_id.units(),
            value,
            uncertainty,
            evidence,
            identity_hash,
            source_lineage,
            region_name,
            tie_witness_vertex: tie_vertex,
        });
    }

    // Canonical deterministic sorting by (semantic_id, region_name, query_name)
    rows.sort_by(|a, b| {
        a.semantic_id
            .cmp(&b.semantic_id)
            .then_with(|| a.region_name.cmp(&b.region_name))
            .then_with(|| a.query_name.cmp(&b.query_name))
    });

    let provenance = ProvenanceHash::of_bytes(b"fs-airflow/registered-qoi/v1");

    Ok(RegisteredQoiExtractionReceipt {
        requested_query_count: queries.len(),
        emitted_qoi_count: rows.len(),
        executed_work_items: mesh.element_count() + mesh.vertex_count(),
        provenance,
        rows,
    })
}
