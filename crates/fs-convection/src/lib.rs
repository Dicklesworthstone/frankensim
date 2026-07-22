//! Regime-bounded convection correlations whose evidence travels with the
//! predicted heat-transfer coefficient.
//!
//! This crate implements the inexpensive correlation rung below resolved
//! airflow. Every formula has one [`CorrelationCard`], uses the shared
//! [`fs_evidence::ValidityDomain`], refuses missing or out-of-domain groups,
//! and returns `Evidence<HeatTransferCoefficient>`. A naked, unitless `h`
//! cannot cross this API.
//!
//! Formula implementation can be numerically verified; physical prediction
//! quality cannot. Empirical-card discrepancy allowances remain model-form
//! evidence, and analytic limiting rows name their idealizing assumptions.

use core::fmt;
use std::collections::BTreeMap;

use fs_evidence::{
    Ambition, Evidence, ModelCard, ModelEvidence, NumericalCertificate, ProvenanceHash,
    SensitivitySummary, StatisticalCertificate, ValidityDomain,
};
use fs_math::det;
use fs_qty::{Length, Qty, Temperature};

/// Thermal conductivity in coherent SI W/(m K).
pub type ThermalConductivity = Qty<1, 1, -3, -1, 0, 0>;
/// Convective heat-transfer coefficient in coherent SI W/(m² K).
pub type HeatTransferCoefficient = Qty<0, 1, -3, -1, 0, 0>;

/// Stable identifier for one implemented heat-transfer relation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum CorrelationId {
    /// Circular duct, fully developed laminar, constant wall temperature.
    CircularDuctLaminarCwt,
    /// Circular duct, fully developed laminar, constant wall heat flux.
    CircularDuctLaminarChf,
    /// Circular duct, thermally developing Hausen relation.
    CircularDuctHausen,
    /// Rectangular duct, fully developed laminar, constant wall temperature.
    RectangularDuctLaminarCwt,
    /// Rectangular duct, fully developed laminar, constant wall heat flux.
    RectangularDuctLaminarChf,
    /// Smooth turbulent tube, Dittus-Boelter relation.
    DittusBoelter,
    /// Smooth turbulent tube, Gnielinski relation.
    Gnielinski,
    /// Isothermal laminar flat plate, average coefficient.
    FlatPlateLaminarAverage,
    /// Isothermal turbulent flat plate with leading-edge correction.
    FlatPlateTurbulentAverage,
    /// Circular cylinder in crossflow, Churchill-Bernstein relation.
    ChurchillBernsteinCylinder,
    /// Isothermal vertical plate in natural convection, Churchill-Chu relation.
    ChurchillChuVerticalPlate,
}

impl CorrelationId {
    /// Complete catalog order. This is stable audit order, not a selection
    /// preference.
    pub const ALL: [CorrelationId; 11] = [
        Self::CircularDuctLaminarCwt,
        Self::CircularDuctLaminarChf,
        Self::CircularDuctHausen,
        Self::RectangularDuctLaminarCwt,
        Self::RectangularDuctLaminarChf,
        Self::DittusBoelter,
        Self::Gnielinski,
        Self::FlatPlateLaminarAverage,
        Self::FlatPlateTurbulentAverage,
        Self::ChurchillBernsteinCylinder,
        Self::ChurchillChuVerticalPlate,
    ];

    /// Stable model-card identifier.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::CircularDuctLaminarCwt => "convection.circular-duct-laminar-cwt",
            Self::CircularDuctLaminarChf => "convection.circular-duct-laminar-chf",
            Self::CircularDuctHausen => "convection.circular-duct-hausen-developing",
            Self::RectangularDuctLaminarCwt => "convection.rectangular-duct-laminar-cwt",
            Self::RectangularDuctLaminarChf => "convection.rectangular-duct-laminar-chf",
            Self::DittusBoelter => "convection.dittus-boelter",
            Self::Gnielinski => "convection.gnielinski",
            Self::FlatPlateLaminarAverage => "convection.flat-plate-laminar-average",
            Self::FlatPlateTurbulentAverage => "convection.flat-plate-turbulent-average",
            Self::ChurchillBernsteinCylinder => "convection.churchill-bernstein-cylinder",
            Self::ChurchillChuVerticalPlate => "convection.churchill-chu-vertical-plate",
        }
    }

    /// Stable formula description retained with reports and failures.
    #[must_use]
    pub const fn formula(self) -> &'static str {
        match self {
            Self::CircularDuctLaminarCwt => "Nu = 3.66",
            Self::CircularDuctLaminarChf => "Nu = 4.36",
            Self::CircularDuctHausen => {
                "Nu = 3.66 + 0.0668 Gz / (1 + 0.04 Gz^(2/3)); Gz = Re Pr / (L/Dh)"
            }
            Self::RectangularDuctLaminarCwt => {
                "Nu = 7.541(1 - 2.610a + 4.970a2 - 5.119a3 + 2.702a4 - 0.548a5)"
            }
            Self::RectangularDuctLaminarChf => {
                "Nu = 8.235(1 - 2.0421a + 3.0853a2 - 2.4765a3 + 1.0578a4 - 0.1861a5)"
            }
            Self::DittusBoelter => "Nu = 0.023 Re^0.8 Pr^n; n=0.4 heating, 0.3 cooling",
            Self::Gnielinski => {
                "Nu = (f/8)(Re-1000)Pr / (1 + 12.7 sqrt(f/8)(Pr^(2/3)-1)); f=(0.79 ln Re-1.64)^-2"
            }
            Self::FlatPlateLaminarAverage => "Nu_L = 0.664 Re_L^(1/2) Pr^(1/3)",
            Self::FlatPlateTurbulentAverage => "Nu_L = (0.037 Re_L^0.8 - 871) Pr^(1/3)",
            Self::ChurchillBernsteinCylinder => {
                "Nu_D = 0.3 + 0.62 Re^(1/2) Pr^(1/3)/(1+(0.4/Pr)^(2/3))^(1/4) * (1+(Re/282000)^(5/8))^(4/5)"
            }
            Self::ChurchillChuVerticalPlate => {
                "Nu_L = (0.825 + 0.387 Ra^(1/6)/(1+(0.492/Pr)^(9/16))^(8/27))^2"
            }
        }
    }
}

/// Whether the fluid is heated or cooled in the Dittus-Boelter convention.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ThermalDirection {
    /// The wall heats the bulk fluid (`n = 0.4`).
    #[default]
    HeatingFluid,
    /// The wall cools the bulk fluid (`n = 0.3`).
    CoolingFluid,
}

/// Why a model-form discrepancy number appears on a card.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DiscrepancyBasis {
    /// Closed-form ideal limit; zero denotes no empirical fit residual, not
    /// zero error for a real apparatus.
    AnalyticIdealLimit,
    /// Conservative engineering allowance used until a retained validation
    /// dataset replaces it. It is not represented as source-published data.
    EngineeringAllowance,
}

/// Bibliographic authority retained by a correlation card.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SourceProvenance {
    /// Human-readable source citation including edition or journal location.
    pub citation: &'static str,
    /// Stable DOI, ISBN, report identifier, or catalog identity.
    pub identifier: &'static str,
}

/// One immutable correlation card plus the source semantics not carried by
/// the generic evidence model card.
#[derive(Debug, Clone, PartialEq)]
pub struct CorrelationCard {
    /// Stable implemented relation.
    pub id: CorrelationId,
    /// Shared FrankenSim model card and validity box.
    pub model: ModelCard,
    /// Formula source.
    pub source: SourceProvenance,
    /// Whether the discrepancy is an ideal-limit zero or an explicit
    /// engineering allowance.
    pub discrepancy_basis: DiscrepancyBasis,
}

/// Dimensionless inputs used to evaluate one card.
///
/// Optional fields stay explicit because different cards require different
/// group sets. Evaluation reports a missing axis; no constructor guesses one.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CorrelationInputs {
    reynolds: Option<f64>,
    prandtl: Option<f64>,
    length_over_hydraulic_diameter: Option<f64>,
    aspect_ratio: Option<f64>,
    rayleigh: Option<f64>,
    direction: ThermalDirection,
}

impl CorrelationInputs {
    /// Forced-convection point with Reynolds and Prandtl numbers.
    #[must_use]
    pub const fn forced(reynolds: f64, prandtl: f64) -> Self {
        Self {
            reynolds: Some(reynolds),
            prandtl: Some(prandtl),
            length_over_hydraulic_diameter: None,
            aspect_ratio: None,
            rayleigh: None,
            direction: ThermalDirection::HeatingFluid,
        }
    }

    /// Natural-convection point with Rayleigh and Prandtl numbers.
    #[must_use]
    pub const fn natural(rayleigh: f64, prandtl: f64) -> Self {
        Self {
            reynolds: None,
            prandtl: Some(prandtl),
            length_over_hydraulic_diameter: None,
            aspect_ratio: None,
            rayleigh: Some(rayleigh),
            direction: ThermalDirection::HeatingFluid,
        }
    }

    /// Declare the duct length-to-hydraulic-diameter ratio.
    #[must_use]
    pub const fn with_length_ratio(mut self, ratio: f64) -> Self {
        self.length_over_hydraulic_diameter = Some(ratio);
        self
    }

    /// Declare the rectangular-duct minor/major side ratio in `(0, 1]`.
    #[must_use]
    pub const fn with_aspect_ratio(mut self, ratio: f64) -> Self {
        self.aspect_ratio = Some(ratio);
        self
    }

    /// Declare the Dittus-Boelter heat-flow direction.
    #[must_use]
    pub const fn with_direction(mut self, direction: ThermalDirection) -> Self {
        self.direction = direction;
        self
    }

    fn groups(self) -> Result<BTreeMap<String, f64>, CorrelationError> {
        let mut groups = BTreeMap::new();
        if let Some(value) = self.reynolds {
            insert_group(&mut groups, "Re", value)?;
        }
        if let Some(value) = self.prandtl {
            insert_group(&mut groups, "Pr", value)?;
        }
        if let Some(value) = self.length_over_hydraulic_diameter {
            insert_group(&mut groups, "L_over_Dh", value)?;
        }
        if let Some(value) = self.aspect_ratio {
            insert_group(&mut groups, "aspect_ratio", value)?;
        }
        if let Some(value) = self.rayleigh {
            insert_group(&mut groups, "Ra", value)?;
        }
        if let (Some(re), Some(pr)) = (self.reynolds, self.prandtl) {
            insert_group(&mut groups, "Pe", re * pr)?;
        }
        Ok(groups)
    }
}

fn insert_group(
    groups: &mut BTreeMap<String, f64>,
    name: &'static str,
    value: f64,
) -> Result<(), CorrelationError> {
    if !(value.is_finite() && value > 0.0) {
        return Err(CorrelationError::InvalidGroup {
            axis: name,
            value_bits: value.to_bits(),
        });
    }
    groups.insert(name.to_string(), value);
    Ok(())
}

/// One failed validity-axis check.
#[derive(Debug, Clone, PartialEq)]
pub struct DomainViolation {
    /// Stable dimensionless-group axis.
    pub axis: String,
    /// Supplied value, or `None` when the required group was absent.
    pub value: Option<f64>,
    /// Inclusive lower validity bound.
    pub low: f64,
    /// Inclusive upper validity bound.
    pub high: f64,
}

/// Structured correlation refusal.
#[derive(Debug, Clone, PartialEq)]
pub enum CorrelationError {
    /// A supplied group is non-finite or not strictly positive.
    InvalidGroup {
        /// Axis name.
        axis: &'static str,
        /// Exact rejected bits.
        value_bits: u64,
    },
    /// Required groups are missing or outside the card validity box.
    OutOfDomain {
        /// Refused correlation.
        correlation: CorrelationId,
        /// Deterministically axis-sorted failures.
        violations: Vec<DomainViolation>,
    },
    /// Conductivity or characteristic length cannot produce a physical `h`.
    InvalidDimensionalInput {
        /// Input name.
        field: &'static str,
        /// Exact rejected bits.
        value_bits: u64,
        /// Required physical condition.
        requirement: &'static str,
    },
    /// Formula arithmetic produced no positive finite Nusselt number or `h`.
    NonFiniteResult {
        /// Stage that failed.
        stage: &'static str,
        /// Exact result bits.
        value_bits: u64,
    },
}

impl fmt::Display for CorrelationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidGroup { axis, value_bits } => write!(
                f,
                "dimensionless group {axis} must be finite and positive; got bits 0x{value_bits:016x}"
            ),
            Self::OutOfDomain {
                correlation,
                violations,
            } => {
                write!(
                    f,
                    "{} refused outside its validity domain:",
                    correlation.name()
                )?;
                for violation in violations {
                    match violation.value {
                        Some(value) => write!(
                            f,
                            " {}={value} not in [{}, {}];",
                            violation.axis, violation.low, violation.high
                        )?,
                        None => write!(
                            f,
                            " {} missing, required in [{}, {}];",
                            violation.axis, violation.low, violation.high
                        )?,
                    }
                }
                Ok(())
            }
            Self::InvalidDimensionalInput {
                field,
                value_bits,
                requirement,
            } => write!(
                f,
                "{field} refused: bits 0x{value_bits:016x}; {requirement}"
            ),
            Self::NonFiniteResult { stage, value_bits } => write!(
                f,
                "{stage} produced a non-positive or non-finite value with bits 0x{value_bits:016x}"
            ),
        }
    }
}

impl core::error::Error for CorrelationError {}

/// A Nusselt prediction whose model evidence and source card cannot be
/// detached from the scalar.
#[derive(Debug, Clone, PartialEq)]
pub struct NusseltEvaluation {
    card: CorrelationCard,
    evidence: Evidence<f64>,
    groups: BTreeMap<String, f64>,
}

/// A conduction Robin row paired with the exact correlation evidence that
/// produced its coefficient.
///
/// Fields are private so the boundary condition and evidence cannot drift
/// after construction.
#[derive(Debug, Clone, PartialEq)]
pub struct CorrelationRobinBoundary {
    coefficient: Evidence<HeatTransferCoefficient>,
    reference_temperature: Temperature,
    boundary: fs_conduction::ThermalBc,
}

impl CorrelationRobinBoundary {
    /// Evidence-bearing coefficient used by the lowered boundary row.
    #[must_use]
    pub const fn coefficient(&self) -> &Evidence<HeatTransferCoefficient> {
        &self.coefficient
    }

    /// Declared ambient/reference temperature.
    #[must_use]
    pub const fn reference_temperature(&self) -> Temperature {
        self.reference_temperature
    }

    /// Conduction boundary row whose raw scalar is pinned to
    /// [`Self::coefficient`].
    #[must_use]
    pub const fn boundary_condition(&self) -> &fs_conduction::ThermalBc {
        &self.boundary
    }
}

impl NusseltEvaluation {
    /// Card used for this evaluation.
    #[must_use]
    pub const fn card(&self) -> &CorrelationCard {
        &self.card
    }

    /// Evidence-bearing dimensionless Nusselt number.
    #[must_use]
    pub const fn evidence(&self) -> &Evidence<f64> {
        &self.evidence
    }

    /// Exact evaluated group point in sorted order.
    #[must_use]
    pub const fn groups(&self) -> &BTreeMap<String, f64> {
        &self.groups
    }

    /// Convert `Nu` to `h = Nu k / L` while retaining the same model card,
    /// validity verdict, discrepancy allowance, and operation provenance.
    ///
    /// # Errors
    /// Refuses non-positive/non-finite conductivity, length, or result.
    pub fn heat_transfer_coefficient(
        &self,
        conductivity: ThermalConductivity,
        characteristic_length: Length,
    ) -> Result<Evidence<HeatTransferCoefficient>, CorrelationError> {
        let k = conductivity.value();
        let length = characteristic_length.value();
        if !(k.is_finite() && k > 0.0) {
            return Err(CorrelationError::InvalidDimensionalInput {
                field: "fluid thermal conductivity",
                value_bits: k.to_bits(),
                requirement: "supply a finite positive value in W/(m K)",
            });
        }
        if !(length.is_finite() && length > 0.0) {
            return Err(CorrelationError::InvalidDimensionalInput {
                field: "characteristic length",
                value_bits: length.to_bits(),
                requirement: "supply a finite positive value in m",
            });
        }
        let h = self.evidence.value * k / length;
        if !(h.is_finite() && h > 0.0) {
            return Err(CorrelationError::NonFiniteResult {
                stage: "Nu-to-h conversion",
                value_bits: h.to_bits(),
            });
        }
        Ok(Evidence {
            value: HeatTransferCoefficient::new(h),
            qoi: h,
            // Formula arithmetic is deterministic but is not a rigorous
            // forward-error enclosure, so the numerical slice stays Estimate.
            numerical: NumericalCertificate::estimate(h, h),
            statistical: StatisticalCertificate::None,
            model: self.evidence.model.clone(),
            sensitivity: SensitivitySummary::default(),
            provenance: ProvenanceHash::chain("convection-nu-to-h-v1", &[self.evidence.provenance]),
            adjoint_ref: None,
        })
    }

    /// Lower this evaluation into a conduction Robin row while retaining the
    /// exact evidence-bearing coefficient beside it.
    ///
    /// # Errors
    /// Refuses invalid conductivity, characteristic length, reference
    /// temperature, or any failure of the conduction boundary admission gate.
    pub fn robin_boundary(
        &self,
        conductivity: ThermalConductivity,
        characteristic_length: Length,
        reference_temperature: Temperature,
    ) -> Result<CorrelationRobinBoundary, CorrelationError> {
        let temperature = reference_temperature.value();
        if !temperature.is_finite() {
            return Err(CorrelationError::InvalidDimensionalInput {
                field: "Robin reference temperature",
                value_bits: temperature.to_bits(),
                requirement: "supply a finite coherent-SI absolute temperature",
            });
        }
        let coefficient = self.heat_transfer_coefficient(conductivity, characteristic_length)?;
        let boundary = fs_conduction::ThermalBc::robin(coefficient.value.value(), temperature)
            .map_err(|_| CorrelationError::NonFiniteResult {
                stage: "conduction Robin boundary admission",
                value_bits: coefficient.value.value().to_bits(),
            })?;
        Ok(CorrelationRobinBoundary {
            coefficient,
            reference_temperature,
            boundary,
        })
    }
}

/// Complete correlation catalog in stable audit order.
#[must_use]
pub fn correlation_catalog() -> Vec<CorrelationCard> {
    CorrelationId::ALL.into_iter().map(card_for).collect()
}

/// Evaluate one relation at a dimensionless point.
///
/// # Errors
/// [`CorrelationError::InvalidGroup`], [`CorrelationError::OutOfDomain`],
/// or [`CorrelationError::NonFiniteResult`].
pub fn evaluate(
    id: CorrelationId,
    inputs: CorrelationInputs,
) -> Result<NusseltEvaluation, CorrelationError> {
    let card = card_for(id);
    let groups = inputs.groups()?;
    let violations = domain_violations(&card.model.validity, &groups);
    if !violations.is_empty() {
        return Err(CorrelationError::OutOfDomain {
            correlation: id,
            violations,
        });
    }
    let nu = evaluate_formula(id, inputs, &groups);
    if !(nu.is_finite() && nu > 0.0) {
        return Err(CorrelationError::NonFiniteResult {
            stage: "Nusselt correlation",
            value_bits: nu.to_bits(),
        });
    }
    let model = ModelEvidence::from_card(&card.model, &groups);
    debug_assert!(model.in_domain);
    let provenance = evaluation_provenance(id, &groups, inputs.direction);
    Ok(NusseltEvaluation {
        evidence: Evidence {
            value: nu,
            qoi: nu,
            numerical: NumericalCertificate::estimate(nu, nu),
            statistical: StatisticalCertificate::None,
            model,
            sensitivity: SensitivitySummary::default(),
            provenance,
            adjoint_ref: None,
        },
        card,
        groups,
    })
}

fn domain_violations(
    domain: &ValidityDomain,
    groups: &BTreeMap<String, f64>,
) -> Vec<DomainViolation> {
    domain
        .bounds()
        .iter()
        .filter_map(|(axis, &(low, high))| {
            let value = groups.get(axis).copied();
            if value.is_some_and(|v| v.is_finite() && v >= low && v <= high) {
                None
            } else {
                Some(DomainViolation {
                    axis: axis.clone(),
                    value,
                    low,
                    high,
                })
            }
        })
        .collect()
}

fn evaluation_provenance(
    id: CorrelationId,
    groups: &BTreeMap<String, f64>,
    direction: ThermalDirection,
) -> ProvenanceHash {
    let mut bytes = Vec::new();
    bytes.extend_from_slice(b"org.frankensim.fs-convection.evaluation.v1\0");
    bytes.extend_from_slice(id.name().as_bytes());
    bytes.push(match direction {
        ThermalDirection::HeatingFluid => 0,
        ThermalDirection::CoolingFluid => 1,
    });
    for (name, value) in groups {
        bytes.extend_from_slice(&(name.len() as u64).to_le_bytes());
        bytes.extend_from_slice(name.as_bytes());
        bytes.extend_from_slice(&value.to_bits().to_le_bytes());
    }
    ProvenanceHash::of_bytes(&bytes)
}

fn evaluate_formula(
    id: CorrelationId,
    inputs: CorrelationInputs,
    groups: &BTreeMap<String, f64>,
) -> f64 {
    let group = |name: &str| groups[name];
    match id {
        CorrelationId::CircularDuctLaminarCwt => 3.66,
        CorrelationId::CircularDuctLaminarChf => 4.36,
        CorrelationId::CircularDuctHausen => {
            let graetz = group("Re") * group("Pr") / group("L_over_Dh");
            3.66 + 0.0668 * graetz / (1.0 + 0.04 * det::pow(graetz, 2.0 / 3.0))
        }
        CorrelationId::RectangularDuctLaminarCwt => {
            let a = group("aspect_ratio");
            let shape = (-0.548f64)
                .mul_add(a, 2.702)
                .mul_add(a, -5.119)
                .mul_add(a, 4.970)
                .mul_add(a, -2.610)
                .mul_add(a, 1.0);
            7.541 * shape
        }
        CorrelationId::RectangularDuctLaminarChf => {
            let a = group("aspect_ratio");
            let shape = (-0.1861f64)
                .mul_add(a, 1.0578)
                .mul_add(a, -2.4765)
                .mul_add(a, 3.0853)
                .mul_add(a, -2.0421)
                .mul_add(a, 1.0);
            8.235 * shape
        }
        CorrelationId::DittusBoelter => {
            let exponent = match inputs.direction {
                ThermalDirection::HeatingFluid => 0.4,
                ThermalDirection::CoolingFluid => 0.3,
            };
            0.023 * det::pow(group("Re"), 0.8) * det::pow(group("Pr"), exponent)
        }
        CorrelationId::Gnielinski => {
            let re = group("Re");
            let pr = group("Pr");
            let friction = 1.0 / det::pow(0.79f64.mul_add(det::ln(re), -1.64), 2.0);
            let numerator = (friction / 8.0) * (re - 1_000.0) * pr;
            let denominator =
                1.0 + 12.7 * det::sqrt(friction / 8.0) * (det::pow(pr, 2.0 / 3.0) - 1.0);
            numerator / denominator
        }
        CorrelationId::FlatPlateLaminarAverage => {
            0.664 * det::sqrt(group("Re")) * det::pow(group("Pr"), 1.0 / 3.0)
        }
        CorrelationId::FlatPlateTurbulentAverage => {
            (0.037 * det::pow(group("Re"), 0.8) - 871.0) * det::pow(group("Pr"), 1.0 / 3.0)
        }
        CorrelationId::ChurchillBernsteinCylinder => {
            let re = group("Re");
            let pr = group("Pr");
            let numerator = 0.62 * det::sqrt(re) * det::pow(pr, 1.0 / 3.0);
            let prandtl = det::pow(1.0 + det::pow(0.4 / pr, 2.0 / 3.0), 1.0 / 4.0);
            let reynolds = det::pow(1.0 + det::pow(re / 282_000.0, 5.0 / 8.0), 4.0 / 5.0);
            0.3 + numerator / prandtl * reynolds
        }
        CorrelationId::ChurchillChuVerticalPlate => {
            let ra = group("Ra");
            let pr = group("Pr");
            let denominator = det::pow(1.0 + det::pow(0.492 / pr, 9.0 / 16.0), 8.0 / 27.0);
            let term = 0.825 + 0.387 * det::pow(ra, 1.0 / 6.0) / denominator;
            term * term
        }
    }
}

fn card_for(id: CorrelationId) -> CorrelationCard {
    const SHAH_LONDON: SourceProvenance = SourceProvenance {
        citation: "R. K. Shah and A. L. London, Laminar Flow Forced Convection in Ducts, Academic Press, 1978",
        identifier: "ISBN 978-0-12-020051-1",
    };
    const DITTUS_BOELTER: SourceProvenance = SourceProvenance {
        citation: "F. W. Dittus and L. M. K. Boelter, Heat Transfer in Automobile Radiators of the Tubular Type, University of California Publications in Engineering 2(13), 1930, 443-461",
        identifier: "OCLC 8882310",
    };
    const GNIELINSKI: SourceProvenance = SourceProvenance {
        citation: "V. Gnielinski, New Equations for Heat and Mass Transfer in Turbulent Pipe and Channel Flow, International Chemical Engineering 16(2), 1976, 359-368",
        identifier: "ISSN 0020-6318; volume 16 issue 2 pages 359-368",
    };
    const FLAT_PLATE: SourceProvenance = SourceProvenance {
        citation: "E. Pohlhausen, Der Waermeaustausch zwischen festen Koerpern und Fluessigkeiten mit kleiner Reibung und kleiner Waermeleitung, ZAMM 1, 1921, 115-121",
        identifier: "ZAMM 1 (1921) 115-121",
    };
    const CHURCHILL_BERNSTEIN: SourceProvenance = SourceProvenance {
        citation: "S. W. Churchill and M. Bernstein, A Correlating Equation for Forced Convection from Gases and Liquids to a Circular Cylinder in Crossflow, Journal of Heat Transfer 99(2), 1977, 300-306",
        identifier: "doi:10.1115/1.3450685",
    };
    const CHURCHILL_CHU: SourceProvenance = SourceProvenance {
        citation: "S. W. Churchill and H. H. S. Chu, Correlating Equations for Laminar and Turbulent Free Convection from a Vertical Plate, International Journal of Heat and Mass Transfer 18, 1975, 1323-1329",
        identifier: "doi:10.1016/0017-9310(75)90243-4",
    };

    let fully_developed = || {
        ValidityDomain::unconstrained()
            .with("L_over_Dh", 50.0, 1.0e6)
            .with("Pr", 0.6, 1_000.0)
            .with("Re", 1.0, 2_300.0)
    };
    let (source, validity, assumptions, failures, discrepancy_rel, discrepancy_basis) = match id {
        CorrelationId::CircularDuctLaminarCwt => (
            SHAH_LONDON,
            fully_developed(),
            vec![
                "circular smooth duct",
                "fully developed laminar flow",
                "constant wall temperature",
                "constant properties",
            ],
            vec![
                "entrance effects",
                "transition or turbulence",
                "non-circular geometry",
            ],
            0.0,
            DiscrepancyBasis::AnalyticIdealLimit,
        ),
        CorrelationId::CircularDuctLaminarChf => (
            SHAH_LONDON,
            fully_developed(),
            vec![
                "circular smooth duct",
                "fully developed laminar flow",
                "uniform wall heat flux",
                "constant properties",
            ],
            vec![
                "entrance effects",
                "transition or turbulence",
                "non-circular geometry",
            ],
            0.0,
            DiscrepancyBasis::AnalyticIdealLimit,
        ),
        CorrelationId::CircularDuctHausen => (
            SHAH_LONDON,
            ValidityDomain::unconstrained()
                .with("L_over_Dh", 0.05, 1_000.0)
                .with("Pr", 0.6, 1_000.0)
                .with("Re", 1.0, 2_300.0),
            vec![
                "circular smooth duct",
                "hydrodynamically developed laminar flow",
                "thermally developing flow",
                "constant wall temperature",
                "constant properties",
            ],
            vec![
                "simultaneously developing velocity field",
                "transition or turbulence",
                "strong property variation",
            ],
            0.15,
            DiscrepancyBasis::EngineeringAllowance,
        ),
        CorrelationId::RectangularDuctLaminarCwt => (
            SHAH_LONDON,
            fully_developed().with("aspect_ratio", 0.001, 1.0),
            vec![
                "rectangular smooth duct",
                "fully developed laminar flow",
                "constant wall temperature on all four walls",
                "aspect ratio is minor side over major side",
            ],
            vec![
                "unequal wall temperatures",
                "developing flow",
                "rough or interrupted fins",
            ],
            0.0,
            DiscrepancyBasis::AnalyticIdealLimit,
        ),
        CorrelationId::RectangularDuctLaminarChf => (
            SHAH_LONDON,
            fully_developed().with("aspect_ratio", 0.001, 1.0),
            vec![
                "rectangular smooth duct",
                "fully developed laminar flow",
                "uniform heat flux on all four walls",
                "aspect ratio is minor side over major side",
            ],
            vec![
                "unequal wall heat flux",
                "developing flow",
                "rough or interrupted fins",
            ],
            0.0,
            DiscrepancyBasis::AnalyticIdealLimit,
        ),
        CorrelationId::DittusBoelter => (
            DITTUS_BOELTER,
            ValidityDomain::unconstrained()
                .with("L_over_Dh", 10.0, 1.0e6)
                .with("Pr", 0.7, 120.0)
                .with("Re", 10_000.0, 120_000.0),
            vec![
                "smooth circular tube",
                "fully developed turbulent flow",
                "modest wall-to-bulk property variation",
                "heating/cooling exponent declared",
            ],
            vec![
                "transition",
                "roughness",
                "strong property variation",
                "non-circular secondary flow",
            ],
            0.25,
            DiscrepancyBasis::EngineeringAllowance,
        ),
        CorrelationId::Gnielinski => (
            GNIELINSKI,
            ValidityDomain::unconstrained()
                .with("L_over_Dh", 10.0, 1.0e6)
                .with("Pr", 0.5, 2_000.0)
                .with("Re", 3_000.0, 5.0e6),
            vec![
                "smooth circular tube",
                "fully developed forced convection",
                "Darcy friction factor relation embedded",
                "constant properties",
            ],
            vec![
                "roughness",
                "swirl or strong secondary flow",
                "strong property variation",
            ],
            0.10,
            DiscrepancyBasis::EngineeringAllowance,
        ),
        CorrelationId::FlatPlateLaminarAverage => (
            FLAT_PLATE,
            ValidityDomain::unconstrained()
                .with("Pr", 0.6, 50.0)
                .with("Re", 1.0, 5.0e5),
            vec![
                "isothermal flat plate",
                "zero pressure gradient",
                "laminar attached boundary layer",
                "average from leading edge",
            ],
            vec![
                "transition",
                "separation",
                "finite-width edge effects",
                "surface roughness",
            ],
            0.10,
            DiscrepancyBasis::EngineeringAllowance,
        ),
        CorrelationId::FlatPlateTurbulentAverage => (
            FLAT_PLATE,
            ValidityDomain::unconstrained()
                .with("Pr", 0.6, 60.0)
                .with("Re", 5.0e5, 1.0e7),
            vec![
                "isothermal smooth flat plate",
                "zero pressure gradient",
                "turbulent attached boundary layer",
                "leading-edge correction retained",
            ],
            vec![
                "laminar or transitional plate",
                "separation",
                "finite-width edge effects",
                "surface roughness",
            ],
            0.15,
            DiscrepancyBasis::EngineeringAllowance,
        ),
        CorrelationId::ChurchillBernsteinCylinder => (
            CHURCHILL_BERNSTEIN,
            ValidityDomain::unconstrained()
                .with("Pe", 0.2, 1.0e10)
                .with("Pr", 0.2, 1_000.0)
                .with("Re", 1.0, 1.0e7),
            vec![
                "isolated circular cylinder in crossflow",
                "properties evaluated consistently",
                "crossflow normal to cylinder axis",
            ],
            vec![
                "tube-bank interference",
                "nearby walls",
                "yawed flow",
                "mixed natural convection",
            ],
            0.20,
            DiscrepancyBasis::EngineeringAllowance,
        ),
        CorrelationId::ChurchillChuVerticalPlate => (
            CHURCHILL_CHU,
            ValidityDomain::unconstrained()
                .with("Pr", 0.01, 1.0e5)
                .with("Ra", 0.1, 1.0e12),
            vec![
                "isothermal vertical plate",
                "quiescent ambient",
                "buoyancy-driven external flow",
                "properties evaluated consistently",
            ],
            vec![
                "forced or mixed convection",
                "enclosure confinement",
                "inclined or horizontal surface",
            ],
            0.20,
            DiscrepancyBasis::EngineeringAllowance,
        ),
    };
    CorrelationCard {
        id,
        model: ModelCard::new(
            id.name(),
            "1.0.0",
            Ambition::Solid,
            assumptions.into_iter().map(str::to_string).collect(),
            validity,
            failures.into_iter().map(str::to_string).collect(),
            discrepancy_rel,
        ),
        source,
        discrepancy_basis,
    }
}
