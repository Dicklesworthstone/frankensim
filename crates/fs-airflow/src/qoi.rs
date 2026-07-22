//! Evidence-bearing quantities of interest for the thermal vertical.
//!
//! This module is the narrow consumer seam between the steady conduction
//! field and the enclosure-airflow operating point. It deliberately does not
//! infer temperature error from a residual, treat an estimated airflow band
//! as a physical certificate, or turn a missing receipt into zero. Every QoI
//! carries an exactly-eight-source [`EngineeringUncertaintyBudget`]; sources
//! without a valid propagation path are explicit [`TermValue::Unknown`].

use core::fmt;
use std::collections::BTreeSet;

use fs_blake3::{ContentHash, hash_domain};
use fs_conduction::{ConductionMesh, ConductionSolution};
use fs_evidence::uncertainty::{
    BudgetTotal, EngineeringUncertaintyBudget, EngineeringUncertaintyKind,
    EngineeringUncertaintyTerm, TermValue, UncertaintyArtifactRef, UncertaintyError,
};
use fs_evidence::{
    Evidence, ModelEvidence, NumericalCertificate, ProvenanceHash, SensitivitySummary,
    StatisticalCertificate, ValidityDomain,
};
use fs_qty::{Power, Pressure, Temperature};

use crate::{OperatingPoint, SourceProvenance};

const QOI_IDENTITY_DOMAIN: &str = "org.frankensim.fs-airflow.thermal-qoi.v1";
const QOI_TERM_DOMAIN: &str = "org.frankensim.fs-airflow.thermal-qoi.term.v1";
const MAX_REGION_ENTRIES: usize = 1_000_000;

/// A declared set of conduction vertices representing one component region.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct JunctionRegion {
    name: String,
    vertices: Vec<usize>,
}

impl JunctionRegion {
    /// Admit a non-empty component region in canonical vertex order.
    pub fn try_new(name: impl Into<String>, mut vertices: Vec<usize>) -> Result<Self, QoiError> {
        let name = admit_name("junction region", name.into())?;
        canonicalize_indices("junction vertices", &mut vertices)?;
        Ok(Self { name, vertices })
    }

    /// Stable region name.
    #[must_use]
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Canonically ordered vertex indices.
    #[must_use]
    pub fn vertices(&self) -> &[usize] {
        &self.vertices
    }
}

/// A declared set of boundary faces representing one reporting surface.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SurfaceRegion {
    name: String,
    boundary_faces: Vec<usize>,
}

impl SurfaceRegion {
    /// Admit a non-empty surface region in canonical boundary-face order.
    pub fn try_new(
        name: impl Into<String>,
        mut boundary_faces: Vec<usize>,
    ) -> Result<Self, QoiError> {
        let name = admit_name("surface region", name.into())?;
        canonicalize_indices("surface boundary faces", &mut boundary_faces)?;
        Ok(Self {
            name,
            boundary_faces,
        })
    }

    /// Stable region name.
    #[must_use]
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Canonically ordered boundary-face slots.
    #[must_use]
    pub fn boundary_faces(&self) -> &[usize] {
        &self.boundary_faces
    }
}

/// Total fan-efficiency input needed to turn air power into fan input power.
#[derive(Debug, Clone, PartialEq)]
pub struct FanPowerSpec {
    /// Nominal total efficiency in `(0, 1]`.
    total_efficiency: f64,
    /// Declared absolute efficiency half-width. The full interval must remain
    /// inside `(0, 1]`.
    efficiency_half_width: f64,
    /// Authority for the efficiency value and its allowance.
    source: SourceProvenance,
}

impl FanPowerSpec {
    /// Validate a fan-power conversion input.
    pub fn try_new(
        total_efficiency: f64,
        efficiency_half_width: f64,
        source: SourceProvenance,
    ) -> Result<Self, QoiError> {
        if source.citation.trim().is_empty() || source.identifier.trim().is_empty() {
            return Err(QoiError::invalid(
                "fan efficiency source",
                "citation and stable identifier must both be non-empty",
            ));
        }
        let low = total_efficiency - efficiency_half_width;
        let high = total_efficiency + efficiency_half_width;
        if !(total_efficiency.is_finite()
            && efficiency_half_width.is_finite()
            && efficiency_half_width >= 0.0
            && low > 0.0
            && high <= 1.0)
        {
            return Err(QoiError::invalid(
                "fan efficiency",
                "nominal +/- half-width must be a finite interval inside (0, 1]",
            ));
        }
        Ok(Self {
            total_efficiency,
            efficiency_half_width,
            source,
        })
    }

    /// Nominal total efficiency.
    #[must_use]
    pub const fn total_efficiency(&self) -> f64 {
        self.total_efficiency
    }

    /// Declared absolute efficiency half-width.
    #[must_use]
    pub const fn efficiency_half_width(&self) -> f64 {
        self.efficiency_half_width
    }

    /// Authority for the efficiency value and allowance.
    #[must_use]
    pub const fn source(&self) -> &SourceProvenance {
        &self.source
    }
}

/// Maximum admissible component temperature from a cited requirement.
#[derive(Debug, Clone, PartialEq)]
pub struct ThermalRequirement {
    /// Maximum permitted temperature in kelvin.
    maximum_temperature: Temperature,
    /// Requirement/standard authority.
    source: SourceProvenance,
}

impl ThermalRequirement {
    /// Admit a finite positive absolute-temperature requirement.
    pub fn try_new(
        maximum_temperature: Temperature,
        source: SourceProvenance,
    ) -> Result<Self, QoiError> {
        if !(maximum_temperature.value().is_finite() && maximum_temperature.value() > 0.0) {
            return Err(QoiError::invalid(
                "thermal requirement",
                "maximum temperature must be finite and positive in kelvin",
            ));
        }
        if source.citation.trim().is_empty() || source.identifier.trim().is_empty() {
            return Err(QoiError::invalid(
                "thermal requirement source",
                "citation and stable identifier must both be non-empty",
            ));
        }
        Ok(Self {
            maximum_temperature,
            source,
        })
    }

    /// Maximum permitted temperature.
    #[must_use]
    pub const fn maximum_temperature(&self) -> Temperature {
        self.maximum_temperature
    }

    /// Requirement/standard authority.
    #[must_use]
    pub const fn source(&self) -> &SourceProvenance {
        &self.source
    }
}

/// One typed QoI with its legacy evidence carrier and rich engineering budget.
#[derive(Debug, Clone, PartialEq)]
pub struct ThermalQoi<T> {
    /// Typed value plus the existing numerical/statistical/model evidence.
    pub evidence: Evidence<T>,
    /// Exactly-eight-source engineering uncertainty budget.
    pub uncertainty: EngineeringUncertaintyBudget,
}

/// Deterministically selected maximum and its canonical vertex witness.
#[derive(Debug, Clone, PartialEq)]
pub struct JunctionMaximum {
    /// Evidence-bearing maximum temperature.
    pub qoi: ThermalQoi<Temperature>,
    /// Lowest canonical vertex index among equal maxima.
    pub vertex: usize,
}

/// Area-weighted surface temperature and two uniformity diagnostics.
#[derive(Debug, Clone, PartialEq)]
pub struct SurfaceUniformity {
    /// Exact P1 face-integral mean divided by selected area.
    pub mean_temperature: ThermalQoi<Temperature>,
    /// Maximum minus minimum selected surface-vertex temperature.
    pub spread: ThermalQoi<Temperature>,
    /// Area-weighted standard deviation of face-mean temperatures.
    pub face_mean_standard_deviation: ThermalQoi<Temperature>,
}

/// The five E05.10 QoI families emitted together from one source snapshot.
#[derive(Debug, Clone, PartialEq)]
pub struct ThermalQoiSet {
    /// Maximum temperature over the declared component region.
    pub junction_maximum: JunctionMaximum,
    /// Enclosure pressure drop at the solved operating point.
    pub pressure_drop: ThermalQoi<Pressure>,
    /// Fan input power `Delta p * Q / eta_total`.
    pub fan_power: ThermalQoi<Power>,
    /// Surface mean/spread/std family.
    pub uniformity: SurfaceUniformity,
    /// `requirement maximum - junction maximum`.
    pub thermal_margin: ThermalQoi<Temperature>,
}

impl ThermalQoiSet {
    /// Every emitted budget, including all three uniformity records.
    #[must_use]
    pub fn budgets(&self) -> [&EngineeringUncertaintyBudget; 7] {
        [
            &self.junction_maximum.qoi.uncertainty,
            &self.pressure_drop.uncertainty,
            &self.fan_power.uncertainty,
            &self.uniformity.mean_temperature.uncertainty,
            &self.uniformity.spread.uncertainty,
            &self.uniformity.face_mean_standard_deviation.uncertainty,
            &self.thermal_margin.uncertainty,
        ]
    }

    /// True only when every QoI retains an explicit unknown source. This is a
    /// useful integration assertion while the vertical lacks complete DWR,
    /// model-validation, and measurement propagation.
    #[must_use]
    pub fn all_totals_are_honestly_unknown(&self) -> bool {
        self.budgets()
            .iter()
            .all(|budget| matches!(budget.total(), BudgetTotal::Unknown { .. }))
    }
}

/// Deterministic refusal from thermal QoI extraction.
#[derive(Debug, Clone, PartialEq)]
pub enum QoiError {
    /// A caller declaration or upstream record is malformed.
    InvalidInput {
        /// Stable input category.
        field: &'static str,
        /// Actionable diagnosis.
        detail: String,
    },
    /// Margin extraction has no cited requirement.
    MissingRequirement,
    /// The eight-term evidence layer refused a malformed budget.
    Uncertainty(UncertaintyError),
}

impl QoiError {
    fn invalid(field: &'static str, detail: impl Into<String>) -> Self {
        Self::InvalidInput {
            field,
            detail: detail.into(),
        }
    }
}

impl fmt::Display for QoiError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidInput { field, detail } => write!(f, "invalid {field}: {detail}"),
            Self::MissingRequirement => write!(
                f,
                "thermal margin requires a cited maximum-temperature requirement; no default is admissible"
            ),
            Self::Uncertainty(error) => write!(f, "uncertainty budget refused: {error}"),
        }
    }
}

impl std::error::Error for QoiError {}

impl From<UncertaintyError> for QoiError {
    fn from(error: UncertaintyError) -> Self {
        Self::Uncertainty(error)
    }
}

/// Extract all E05.10 QoI families from one conduction/airflow snapshot.
///
/// The field values are useful estimates, but temperature extrema and surface
/// metrics receive `NumericalKind::NoClaim` until a DWR/refinement path binds
/// an actual error enclosure. Pressure and fan power retain the upstream
/// estimate bands. Every budget separately retains an unknown model-form term,
/// so conditional parameter/BC envelopes cannot become product authority.
///
/// # Errors
/// Refuses malformed regions, mismatched/non-finite fields, invalid operating
/// evidence, invalid fan efficiency, or an absent requirement.
pub fn extract_thermal_qois(
    mesh: &ConductionMesh,
    solution: &ConductionSolution,
    operating_point: &OperatingPoint,
    junction_region: &JunctionRegion,
    surface_region: &SurfaceRegion,
    fan_power: &FanPowerSpec,
    requirement: Option<&ThermalRequirement>,
) -> Result<ThermalQoiSet, QoiError> {
    let requirement = requirement.ok_or(QoiError::MissingRequirement)?;
    validate_solution(mesh, solution)?;
    validate_operating_point(operating_point)?;
    validate_region_indices(mesh, junction_region, surface_region)?;

    let solution_id = solution_identity(mesh, solution);
    let operating_id = operating_identity(operating_point);
    let fan_id = fan_power_identity(fan_power);
    let requirement_id = requirement_identity(requirement);

    let (maximum, maximum_vertex) = junction_maximum(solution, junction_region);
    let temperature_model = conduction_model(solution);
    let maximum_identity = qoi_identity(
        "junction-maximum",
        &[solution_id],
        &region_identity(junction_region.name(), junction_region.vertices()),
    );
    let maximum_budget = unknown_budget("thermal-junction-maximum", "K", maximum_identity, &[])?;
    let maximum_evidence =
        no_claim_temperature(maximum, temperature_model.clone(), maximum_identity);
    let junction_maximum = JunctionMaximum {
        qoi: ThermalQoi {
            evidence: maximum_evidence,
            uncertainty: maximum_budget,
        },
        vertex: maximum_vertex,
    };

    let surface = surface_statistics(mesh, solution, surface_region)?;
    let surface_region_id = region_identity(surface_region.name(), surface_region.boundary_faces());
    let mean_identity = qoi_identity("surface-mean", &[solution_id], &surface_region_id);
    let spread_identity = qoi_identity("surface-spread", &[solution_id], &surface_region_id);
    let std_identity = qoi_identity("surface-face-mean-std", &[solution_id], &surface_region_id);
    let uniformity = SurfaceUniformity {
        mean_temperature: ThermalQoi {
            evidence: no_claim_temperature(
                Temperature::new(surface.mean),
                temperature_model.clone(),
                mean_identity,
            ),
            uncertainty: unknown_budget("thermal-surface-mean", "K", mean_identity, &[])?,
        },
        spread: ThermalQoi {
            evidence: no_claim_temperature(
                Temperature::new(surface.spread),
                temperature_model.clone(),
                spread_identity,
            ),
            uncertainty: unknown_budget("thermal-surface-spread", "K", spread_identity, &[])?,
        },
        face_mean_standard_deviation: ThermalQoi {
            evidence: no_claim_temperature(
                Temperature::new(surface.standard_deviation),
                temperature_model.clone(),
                std_identity,
            ),
            uncertainty: unknown_budget("thermal-surface-face-mean-std", "K", std_identity, &[])?,
        },
    };

    let pressure_value = operating_point.pressure.value.value();
    let pressure_half_width = interval_half_width(
        pressure_value,
        operating_point.pressure.numerical.lo,
        operating_point.pressure.numerical.hi,
    )?;
    let pressure_identity = qoi_identity("pressure-drop", &[operating_id], b"operating-pressure");
    let pressure_term = TermValue::interval(0.0, pressure_half_width)?;
    let pressure_drop = ThermalQoi {
        evidence: operating_point.pressure.clone(),
        uncertainty: unknown_budget(
            "thermal-pressure-drop",
            "Pa",
            pressure_identity,
            &[(
                EngineeringUncertaintyKind::BoundaryConditions,
                pressure_term,
                "thermal-qoi-operating-envelope",
                operating_id,
            )],
        )?,
    };

    let fan_power_qoi = fan_power_qoi(operating_point, fan_power, operating_id, fan_id)?;

    let margin = requirement.maximum_temperature.value() - maximum.value();
    if !margin.is_finite() {
        return Err(QoiError::invalid(
            "thermal margin",
            "requirement minus junction maximum was non-finite",
        ));
    }
    let margin_identity = qoi_identity(
        "thermal-margin",
        &[maximum_identity, requirement_id],
        requirement.source.identifier.as_bytes(),
    );
    let thermal_margin = ThermalQoi {
        evidence: no_claim_temperature(
            Temperature::new(margin),
            temperature_model,
            margin_identity,
        ),
        uncertainty: unknown_budget("thermal-margin", "K", margin_identity, &[])?,
    };

    Ok(ThermalQoiSet {
        junction_maximum,
        pressure_drop,
        fan_power: fan_power_qoi,
        uniformity,
        thermal_margin,
    })
}

fn fan_power_qoi(
    operating_point: &OperatingPoint,
    spec: &FanPowerSpec,
    operating_id: ContentHash,
    fan_id: ContentHash,
) -> Result<ThermalQoi<Power>, QoiError> {
    let flow = &operating_point.flow;
    let pressure = &operating_point.pressure;
    let eta = spec.total_efficiency;
    let eta_low = eta - spec.efficiency_half_width;
    let eta_high = eta + spec.efficiency_half_width;
    let nominal_air_power = pressure.value.value() * flow.value.value();
    let nominal = nominal_air_power / eta;
    let low = pressure.numerical.lo.max(0.0) * flow.numerical.lo.max(0.0) / eta_high;
    let high = pressure.numerical.hi * flow.numerical.hi / eta_low;
    if !(nominal.is_finite()
        && low.is_finite()
        && high.is_finite()
        && 0.0 <= low
        && low <= nominal
        && nominal <= high)
    {
        return Err(QoiError::invalid(
            "fan power",
            "pressure/flow/efficiency envelope did not produce an ordered finite non-negative power interval",
        ));
    }

    let boundary_low = pressure.numerical.lo.max(0.0) * flow.numerical.lo.max(0.0) / eta;
    let boundary_high = pressure.numerical.hi * flow.numerical.hi / eta;
    let boundary_half_width = interval_half_width(nominal, boundary_low, boundary_high)?;
    let efficiency_low_power = nominal_air_power / eta_high;
    let efficiency_high_power = nominal_air_power / eta_low;
    let efficiency_half_width =
        interval_half_width(nominal, efficiency_low_power, efficiency_high_power)?;
    let identity = qoi_identity(
        "fan-power",
        &[operating_id, fan_id],
        b"pressure-flow-efficiency",
    );

    let mut model = pressure.model.clone();
    model
        .cards
        .push(format!("fan-efficiency:{}", spec.source.identifier));
    model.cards.sort_unstable();
    model.cards.dedup();
    model.assumptions.push(
        "fan input power uses Delta-p times volume-flow divided by declared total efficiency"
            .to_string(),
    );
    model.assumptions.sort_unstable();
    model.assumptions.dedup();
    let leaf = ProvenanceHash::of_bytes(identity.as_bytes());
    let provenance = ProvenanceHash::chain(
        "thermal-qoi-fan-power-v1",
        &[pressure.provenance, flow.provenance, leaf],
    );
    let evidence = Evidence {
        value: Power::new(nominal),
        qoi: nominal,
        numerical: NumericalCertificate::estimate(low, high),
        statistical: StatisticalCertificate::None,
        model,
        sensitivity: SensitivitySummary::default(),
        provenance,
        adjoint_ref: None,
    };
    let uncertainty = unknown_budget(
        "thermal-fan-input-power",
        "W",
        identity,
        &[
            (
                EngineeringUncertaintyKind::Parameters,
                TermValue::interval(0.0, efficiency_half_width)?,
                "thermal-qoi-efficiency-envelope",
                fan_id,
            ),
            (
                EngineeringUncertaintyKind::BoundaryConditions,
                TermValue::interval(0.0, boundary_half_width)?,
                "thermal-qoi-operating-envelope",
                operating_id,
            ),
        ],
    )?;
    Ok(ThermalQoi {
        evidence,
        uncertainty,
    })
}

fn no_claim_temperature(
    value: Temperature,
    model: ModelEvidence,
    identity: ContentHash,
) -> Evidence<Temperature> {
    Evidence {
        value,
        qoi: value.value(),
        numerical: NumericalCertificate::no_claim(),
        statistical: StatisticalCertificate::None,
        model,
        sensitivity: SensitivitySummary::default(),
        provenance: ProvenanceHash::of_bytes(identity.as_bytes()),
        adjoint_ref: None,
    }
}

fn conduction_model(solution: &ConductionSolution) -> ModelEvidence {
    ModelEvidence {
        cards: vec!["fs-conduction:steady-p1".to_string()],
        assumptions: vec![format!(
            "temperature QoI is a direct functional of a steady P1 field with {} material receipt(s); no DWR-to-QoI bound is attached",
            solution.report.material_receipts
        )],
        validity: ValidityDomain::unconstrained(),
        discrepancy_rel: f64::INFINITY,
        in_domain: true,
    }
}

fn unknown_budget(
    qoi: &str,
    unit: &str,
    qoi_identity: ContentHash,
    known: &[(
        EngineeringUncertaintyKind,
        TermValue,
        &'static str,
        ContentHash,
    )],
) -> Result<EngineeringUncertaintyBudget, QoiError> {
    let mut terms = Vec::with_capacity(EngineeringUncertaintyKind::ALL.len());
    for kind in EngineeringUncertaintyKind::ALL {
        let (value, role, digest) = known
            .iter()
            .find(|(candidate, _, _, _)| *candidate == kind)
            .map(|(_, value, role, digest)| (value.clone(), *role, *digest))
            .unwrap_or((
                TermValue::unknown(unknown_reason(kind))?,
                "thermal-qoi-evidence-gap",
                term_identity(qoi_identity, kind),
            ));
        let provenance = UncertaintyArtifactRef::new(role, digest)?;
        terms.push(EngineeringUncertaintyTerm::try_new(
            kind, value, provenance,
        )?);
    }
    EngineeringUncertaintyBudget::try_new(qoi, unit, terms).map_err(Into::into)
}

fn unknown_reason(kind: EngineeringUncertaintyKind) -> &'static str {
    match kind {
        EngineeringUncertaintyKind::Roundoff => {
            "no retained outward-roundoff propagation receipt reaches this QoI"
        }
        EngineeringUncertaintyKind::SolverAlgebraic => {
            "the recomputed residual has no admitted inverse-stability map into this QoI unit"
        }
        EngineeringUncertaintyKind::Discretization => {
            "no retained DWR or refinement-ladder bound is attached to this QoI"
        }
        EngineeringUncertaintyKind::Geometry => {
            "no import, meshing, registration, or as-built geometry uncertainty is propagated"
        }
        EngineeringUncertaintyKind::Parameters => {
            "no complete material and manufacturing parameter propagation is attached"
        }
        EngineeringUncertaintyKind::BoundaryConditions => {
            "no complete boundary and operating-condition propagation is attached"
        }
        EngineeringUncertaintyKind::ModelForm => {
            "no externally validated model-form discrepancy envelope is attached"
        }
        EngineeringUncertaintyKind::Measurement => {
            "no calibrated observation or comparison-data uncertainty is attached"
        }
        _ => "no authoritative propagation is attached for this uncertainty family",
    }
}

fn validate_solution(mesh: &ConductionMesh, solution: &ConductionSolution) -> Result<(), QoiError> {
    if solution.temperature.len() != mesh.vertex_count() {
        return Err(QoiError::invalid(
            "temperature field",
            format!(
                "{} nodal values for {} mesh vertices",
                solution.temperature.len(),
                mesh.vertex_count()
            ),
        ));
    }
    if solution.temperature.iter().any(|value| !value.is_finite()) {
        return Err(QoiError::invalid(
            "temperature field",
            "every nodal temperature must be finite",
        ));
    }
    Ok(())
}

fn validate_operating_point(point: &OperatingPoint) -> Result<(), QoiError> {
    validate_estimate(
        "operating pressure",
        point.pressure.value.value(),
        point.pressure.numerical.lo,
        point.pressure.numerical.hi,
    )?;
    validate_estimate(
        "operating flow",
        point.flow.value.value(),
        point.flow.numerical.lo,
        point.flow.numerical.hi,
    )?;
    if point.pressure.value.value() < 0.0
        || point.pressure.numerical.lo < 0.0
        || point.flow.value.value() < 0.0
        || point.flow.numerical.lo < 0.0
    {
        return Err(QoiError::invalid(
            "operating point",
            "pressure and volume-flow estimates must be non-negative",
        ));
    }
    Ok(())
}

fn validate_estimate(field: &'static str, value: f64, lo: f64, hi: f64) -> Result<(), QoiError> {
    if !(value.is_finite() && lo.is_finite() && hi.is_finite() && lo <= value && value <= hi) {
        return Err(QoiError::invalid(
            field,
            "value and bounds must be finite, ordered, and self-containing",
        ));
    }
    Ok(())
}

fn validate_region_indices(
    mesh: &ConductionMesh,
    junction: &JunctionRegion,
    surface: &SurfaceRegion,
) -> Result<(), QoiError> {
    if let Some(&vertex) = junction
        .vertices()
        .iter()
        .find(|&&vertex| vertex >= mesh.vertex_count())
    {
        return Err(QoiError::invalid(
            "junction region",
            format!(
                "vertex {vertex} is outside {} vertices",
                mesh.vertex_count()
            ),
        ));
    }
    if let Some(&face) = surface
        .boundary_faces()
        .iter()
        .find(|&&face| face >= mesh.boundary().len())
    {
        return Err(QoiError::invalid(
            "surface region",
            format!(
                "boundary-face slot {face} is outside {} boundary faces",
                mesh.boundary().len()
            ),
        ));
    }
    Ok(())
}

fn junction_maximum(
    solution: &ConductionSolution,
    region: &JunctionRegion,
) -> (Temperature, usize) {
    let mut vertex = region.vertices()[0];
    let mut maximum = solution.temperature[vertex];
    for &candidate in &region.vertices()[1..] {
        let value = solution.temperature[candidate];
        if value > maximum || (value.to_bits() == maximum.to_bits() && candidate < vertex) {
            maximum = value;
            vertex = candidate;
        }
    }
    (Temperature::new(maximum), vertex)
}

struct SurfaceStatistics {
    mean: f64,
    spread: f64,
    standard_deviation: f64,
}

fn surface_statistics(
    mesh: &ConductionMesh,
    solution: &ConductionSolution,
    region: &SurfaceRegion,
) -> Result<SurfaceStatistics, QoiError> {
    let mut area = 0.0;
    let mut integral = 0.0;
    let mut face_means = Vec::with_capacity(region.boundary_faces().len());
    let mut vertices = BTreeSet::new();
    for &slot in region.boundary_faces() {
        let face = &mesh.boundary()[slot];
        let mean = face
            .vertices
            .iter()
            .map(|&vertex| solution.temperature[vertex as usize])
            .sum::<f64>()
            / 3.0;
        area += face.area;
        integral = face.area.mul_add(mean, integral);
        face_means.push((face.area, mean));
        vertices.extend(face.vertices.iter().map(|&vertex| vertex as usize));
    }
    if !(area.is_finite() && area > 0.0 && integral.is_finite()) {
        return Err(QoiError::invalid(
            "surface region",
            "selected area and temperature integral must be finite with positive area",
        ));
    }
    let mean = integral / area;
    let mut variance_integral = 0.0;
    for (face_area, face_mean) in face_means {
        let delta = face_mean - mean;
        variance_integral = face_area.mul_add(delta * delta, variance_integral);
    }
    let variance = (variance_integral / area).max(0.0);
    let mut minimum = f64::INFINITY;
    let mut maximum = f64::NEG_INFINITY;
    for vertex in vertices {
        let temperature = solution.temperature[vertex];
        minimum = minimum.min(temperature);
        maximum = maximum.max(temperature);
    }
    let spread = maximum - minimum;
    let standard_deviation = fs_math::det::sqrt(variance);
    if !(mean.is_finite() && spread.is_finite() && standard_deviation.is_finite()) {
        return Err(QoiError::invalid(
            "surface statistics",
            "mean, spread, or standard deviation was non-finite",
        ));
    }
    Ok(SurfaceStatistics {
        mean,
        spread,
        standard_deviation,
    })
}

fn interval_half_width(value: f64, lo: f64, hi: f64) -> Result<f64, QoiError> {
    validate_estimate("QoI envelope", value, lo, hi)?;
    let half_width = (value - lo).abs().max((hi - value).abs());
    if half_width.is_finite() {
        Ok(half_width)
    } else {
        Err(QoiError::invalid(
            "QoI envelope",
            "half-width arithmetic overflowed",
        ))
    }
}

fn admit_name(field: &'static str, value: String) -> Result<String, QoiError> {
    if value.trim().is_empty() || value.len() > 1024 {
        Err(QoiError::invalid(field, "name must contain 1..=1024 bytes"))
    } else {
        Ok(value)
    }
}

fn canonicalize_indices(field: &'static str, values: &mut Vec<usize>) -> Result<(), QoiError> {
    if values.is_empty() || values.len() > MAX_REGION_ENTRIES {
        return Err(QoiError::invalid(
            field,
            format!("entry count must lie in 1..={MAX_REGION_ENTRIES}"),
        ));
    }
    values.sort_unstable();
    if values.windows(2).any(|pair| pair[0] == pair[1]) {
        return Err(QoiError::invalid(
            field,
            "duplicate indices are not admissible",
        ));
    }
    Ok(())
}

fn solution_identity(mesh: &ConductionMesh, solution: &ConductionSolution) -> ContentHash {
    let mut bytes = Vec::new();
    push_usize(&mut bytes, mesh.vertex_count());
    push_usize(&mut bytes, mesh.element_count());
    for position in mesh.positions() {
        for coordinate in position {
            bytes.extend_from_slice(&coordinate.to_bits().to_le_bytes());
        }
    }
    for tet in &mesh.complex().tets {
        for vertex in tet {
            bytes.extend_from_slice(&vertex.to_le_bytes());
        }
    }
    for value in &solution.temperature {
        bytes.extend_from_slice(&value.to_bits().to_le_bytes());
    }
    bytes.extend_from_slice(&solution.report.final_residual.to_bits().to_le_bytes());
    bytes.extend_from_slice(&solution.report.residual_threshold.to_bits().to_le_bytes());
    bytes.extend_from_slice(&solution.report.energy.closure_w.to_bits().to_le_bytes());
    bytes.extend_from_slice(solution.report.material_provenance.tag().as_bytes());
    push_usize(&mut bytes, solution.report.material_receipts);
    hash_domain(QOI_IDENTITY_DOMAIN, &bytes)
}

fn operating_identity(point: &OperatingPoint) -> ContentHash {
    let mut bytes = Vec::new();
    for value in [
        point.flow.value.value(),
        point.flow.numerical.lo,
        point.flow.numerical.hi,
        point.pressure.value.value(),
        point.pressure.numerical.lo,
        point.pressure.numerical.hi,
        point.nominal_root.flow.lo(),
        point.nominal_root.flow.hi(),
        point.leakage_fraction,
    ] {
        bytes.extend_from_slice(&value.to_bits().to_le_bytes());
    }
    bytes.extend_from_slice(&point.flow.provenance.0.to_le_bytes());
    bytes.extend_from_slice(&point.pressure.provenance.0.to_le_bytes());
    for card in &point.pressure.model.cards {
        push_string(&mut bytes, card);
    }
    hash_domain(QOI_IDENTITY_DOMAIN, &bytes)
}

fn fan_power_identity(spec: &FanPowerSpec) -> ContentHash {
    let mut bytes = Vec::new();
    bytes.extend_from_slice(&spec.total_efficiency.to_bits().to_le_bytes());
    bytes.extend_from_slice(&spec.efficiency_half_width.to_bits().to_le_bytes());
    push_string(&mut bytes, &spec.source.citation);
    push_string(&mut bytes, &spec.source.identifier);
    hash_domain(QOI_IDENTITY_DOMAIN, &bytes)
}

fn requirement_identity(requirement: &ThermalRequirement) -> ContentHash {
    let mut bytes = Vec::new();
    bytes.extend_from_slice(
        &requirement
            .maximum_temperature
            .value()
            .to_bits()
            .to_le_bytes(),
    );
    push_string(&mut bytes, &requirement.source.citation);
    push_string(&mut bytes, &requirement.source.identifier);
    hash_domain(QOI_IDENTITY_DOMAIN, &bytes)
}

fn region_identity(name: &str, indices: &[usize]) -> Vec<u8> {
    let mut bytes = Vec::new();
    push_string(&mut bytes, name);
    for &index in indices {
        push_usize(&mut bytes, index);
    }
    bytes
}

fn qoi_identity(label: &str, parents: &[ContentHash], context: &[u8]) -> ContentHash {
    let mut bytes = Vec::new();
    push_string(&mut bytes, label);
    for parent in parents {
        bytes.extend_from_slice(parent.as_bytes());
    }
    bytes.extend_from_slice(context);
    hash_domain(QOI_IDENTITY_DOMAIN, &bytes)
}

fn term_identity(qoi: ContentHash, kind: EngineeringUncertaintyKind) -> ContentHash {
    let mut bytes = Vec::new();
    bytes.extend_from_slice(qoi.as_bytes());
    bytes.extend_from_slice(kind.name().as_bytes());
    hash_domain(QOI_TERM_DOMAIN, &bytes)
}

fn push_string(bytes: &mut Vec<u8>, value: &str) {
    push_usize(bytes, value.len());
    bytes.extend_from_slice(value.as_bytes());
}

fn push_usize(bytes: &mut Vec<u8>, value: usize) {
    bytes.extend_from_slice(&(value as u64).to_le_bytes());
}
