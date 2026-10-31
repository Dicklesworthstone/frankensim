//! Regime-bounded convection correlations whose evidence travels with the
//! predicted heat-transfer coefficient.
//!
//! This crate implements the inexpensive correlation rung below resolved
//! airflow. Every formula has one [`CorrelationCard`], uses the shared
//! [`fs_evidence::ValidityDomain`], refuses missing or out-of-domain groups,
//! and retains its model evidence through the typed coefficient conversion.
//! The optional Robin lowering keeps that evidence beside a downstream
//! conduction row; the row accessor exposes the downstream raw coherent-SI
//! representation and does not give a detached row correlation authority.
//!
//! Formula implementation can be numerically verified; physical prediction
//! quality cannot. Empirical-card discrepancy allowances remain model-form
//! evidence, and analytic limiting rows name their idealizing assumptions.

use core::fmt;
use std::collections::BTreeMap;

mod mixed;
mod natural;
mod sensitivity;

pub use natural::{HorizontalPlateFacing, NaturalConvectionGeometry};

/// Source-bearing compartment-fire correlations and public-data comparisons.
pub mod fire;

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

const RECTANGULAR_DEVELOPING_FIRST_TABLE_GRAETZ: f64 = 10.0;
const RECTANGULAR_DEVELOPING_BRIDGE_DISCLOSURE: &str = "Gz < 10 uses a declared linear-in-Gz engineering bridge from the fully developed CWT limit to the first Table 52 knot; that bridge is not represented as a source-published curve";

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
    /// Rectangular duct, simultaneously developing laminar flow, Pr = 0.72.
    RectangularDuctLaminarCwtDevelopingPr072,
    /// Smooth turbulent tube, Dittus-Boelter relation.
    DittusBoelter,
    /// Smooth turbulent tube, Gnielinski relation.
    Gnielinski,
    /// Isothermal laminar flat plate, average coefficient.
    FlatPlateLaminarAverage,
    /// Isothermal mixed laminar-to-turbulent flat plate with leading-edge correction.
    FlatPlateTurbulentAverage,
    /// Circular cylinder in crossflow, Churchill-Bernstein relation.
    ChurchillBernsteinCylinder,
    /// Isothermal vertical plate in natural convection, Churchill-Chu relation.
    ChurchillChuVerticalPlate,
    /// Isothermal horizontal cylinder in natural convection; length is diameter.
    ChurchillChuHorizontalCylinder,
    /// Buoyantly unstable horizontal plate, laminar range; length is area/perimeter.
    HorizontalPlateUnstableLaminar,
    /// Buoyantly unstable horizontal plate, turbulent range; length is area/perimeter.
    HorizontalPlateUnstableTurbulent,
    /// Buoyantly stable horizontal plate; length is area/perimeter.
    HorizontalPlateStable,
    /// Laminar aiding mixed convection over an isothermal vertical plate.
    ManningQureshiAidingVerticalPlate,
    /// Isolated isothermal sphere in natural convection; length is diameter.
    ChurchillSphere,
}

impl CorrelationId {
    /// Complete catalog order. This is stable audit order, not a selection
    /// preference.
    pub const ALL: [CorrelationId; 18] = [
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
        Self::RectangularDuctLaminarCwtDevelopingPr072,
        Self::ChurchillChuHorizontalCylinder,
        Self::HorizontalPlateUnstableLaminar,
        Self::HorizontalPlateUnstableTurbulent,
        Self::HorizontalPlateStable,
        Self::ManningQureshiAidingVerticalPlate,
        Self::ChurchillSphere,
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
            Self::RectangularDuctLaminarCwtDevelopingPr072 => {
                "convection.rectangular-duct-laminar-cwt-developing-pr072"
            }
            Self::DittusBoelter => "convection.dittus-boelter",
            Self::Gnielinski => "convection.gnielinski",
            Self::FlatPlateLaminarAverage => "convection.flat-plate-laminar-average",
            Self::FlatPlateTurbulentAverage => "convection.flat-plate-turbulent-average",
            Self::ChurchillBernsteinCylinder => "convection.churchill-bernstein-cylinder",
            Self::ChurchillChuVerticalPlate => "convection.churchill-chu-vertical-plate",
            Self::ChurchillChuHorizontalCylinder => "convection.churchill-chu-horizontal-cylinder",
            Self::HorizontalPlateUnstableLaminar => "convection.horizontal-plate-unstable-laminar",
            Self::HorizontalPlateUnstableTurbulent => {
                "convection.horizontal-plate-unstable-turbulent"
            }
            Self::HorizontalPlateStable => "convection.horizontal-plate-stable",
            Self::ManningQureshiAidingVerticalPlate => {
                "convection.manning-qureshi-aiding-vertical-plate"
            }
            Self::ChurchillSphere => "convection.churchill-sphere",
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
            Self::RectangularDuctLaminarCwtDevelopingPr072 => {
                "piecewise-linear Nu_m,T(Gz) from Shah-London Ch. VII Table 52 at a=0.5, Pr=0.72; Gz=Re Pr/(L/Dh), anchored to the CWT limit at Gz=0"
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
            Self::ChurchillChuHorizontalCylinder => {
                "Nu_D = (0.60 + 0.387 Ra_D^(1/6)/(1+(0.559/Pr)^(9/16))^(8/27))^2"
            }
            Self::HorizontalPlateUnstableLaminar => "Nu_L = 0.54 Ra_L^(1/4); L = A/P",
            Self::HorizontalPlateUnstableTurbulent => "Nu_L = 0.15 Ra_L^(1/3); L = A/P",
            Self::HorizontalPlateStable => "Nu_L = 0.27 Ra_L^(1/4); L = A/P",
            Self::ManningQureshiAidingVerticalPlate => {
                "Nu_L = [(2 Ff sqrt(Re))^n + ((4/3) Fn (Ra/Pr)^(1/4))^n]^(1/n); Ff=0.339 Pr^(1/3)(1+0.100 Pr^(-3/4))^(-2/9); Fn=0.503 Pr^(1/4)(1+0.670 Pr^(-9/16))^(-4/9); n=3.5 Pr^0.075"
            }
            Self::ChurchillSphere => {
                "Nu_D = 2 + 0.589 Ra_D^(1/4)/(1+(0.469/Pr)^(9/16))^(4/9)"
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

/// Which source-authority region produced one admitted evaluation.
///
/// This is an execution marker, not a validation colour. In particular,
/// [`Self::TableSliceInterpolation`] includes both cited table knots and the
/// implementation's linear-in-Graetz chords between them; it does not claim
/// that every interpolated value was printed by the source.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EvaluationSourceRegion {
    /// A cited closed-form relation or analytic limiting constant.
    CitedFormula,
    /// A knot or linear-in-Graetz chord inside the cited Table 52 slice.
    TableSliceInterpolation,
    /// A declared engineering bridge outside the source-published table slice.
    DeclaredEngineeringBridge,
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

    /// Aiding mixed-convection point with nonnegative Reynolds and Rayleigh
    /// magnitudes and positive Prandtl number. Both flow drivers cannot vanish.
    ///
    /// The caller declares an isothermal vertical plate with an external stream
    /// parallel to and aiding buoyancy. Dimensionless magnitudes cannot verify
    /// that geometry or alignment; a dimensional adapter must check their signs.
    /// The mixed card admits either pure-convection limit continuously.
    #[must_use]
    pub const fn mixed(reynolds: f64, rayleigh: f64, prandtl: f64) -> Self {
        Self {
            reynolds: Some(reynolds),
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

    fn groups(self, id: CorrelationId) -> Result<BTreeMap<String, f64>, CorrelationError> {
        let mut groups = BTreeMap::new();
        let mixed = id == CorrelationId::ManningQureshiAidingVerticalPlate;
        if let Some(value) = self.reynolds {
            insert_group(&mut groups, "Re", value, mixed)?;
        }
        if let Some(value) = self.prandtl {
            insert_group(&mut groups, "Pr", value, false)?;
        }
        if let Some(value) = self.length_over_hydraulic_diameter {
            insert_group(&mut groups, "L_over_Dh", value, false)?;
        }
        if let Some(value) = self.aspect_ratio {
            insert_group(&mut groups, "aspect_ratio", value, false)?;
        }
        if let Some(value) = self.rayleigh {
            insert_group(
                &mut groups,
                "Ra",
                value,
                mixed || id == CorrelationId::ChurchillSphere,
            )?;
        }
        if let (Some(re), Some(pr)) = (self.reynolds, self.prandtl) {
            insert_group(&mut groups, "Pe", re * pr, mixed)?;
        }
        if let (Some(re), Some(pr), Some(length_ratio)) = (
            self.reynolds,
            self.prandtl,
            self.length_over_hydraulic_diameter,
        ) {
            insert_group(&mut groups, "Gz", re * pr / length_ratio, mixed)?;
        }
        Ok(groups)
    }
}

fn insert_group(
    groups: &mut BTreeMap<String, f64>,
    name: &'static str,
    value: f64,
    allow_zero: bool,
) -> Result<(), CorrelationError> {
    if !(value.is_finite() && (value > 0.0 || (allow_zero && value == 0.0))) {
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
    /// A supplied group is non-finite or violates its sign requirement.
    /// The mixed card admits zero Re/Ra and their derived groups, but
    /// `max(Re,Ra)` must be positive. The sphere admits Ra=0 as its conductive
    /// limit. Other supplied groups must be strictly positive.
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
                "dimensionless group {axis} must be finite and positive (zero is allowed for the mixed card's individual flow drivers or the sphere's Rayleigh number); got bits 0x{value_bits:016x}"
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
    source_region: EvaluationSourceRegion,
}

/// A conduction Robin row paired with the exact correlation evidence that
/// produced its coefficient.
///
/// Fields are private so the pair cannot drift while held as this type.
/// Callers can inspect and clone the downstream boundary row, so a detached
/// [`fs_conduction::ThermalBc`] must not be treated as correlation evidence.
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

    /// Inspect the conduction boundary row created from [`Self::coefficient`].
    ///
    /// The downstream type stores its heat-transfer coefficient as a raw
    /// coherent-SI scalar. A cloned or otherwise detached row does not carry
    /// this crate's source, validity, or model evidence.
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

    /// Machine-readable source-authority region used for this evaluation.
    #[must_use]
    pub const fn source_region(&self) -> EvaluationSourceRegion {
        self.source_region
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
    let groups = inputs.groups(id)?;
    let violations = domain_violations(&card.model.validity, &groups);
    if !violations.is_empty() {
        return Err(CorrelationError::OutOfDomain {
            correlation: id,
            violations,
        });
    }
    if id == CorrelationId::ManningQureshiAidingVerticalPlate
        && groups["Re"] == 0.0
        && groups["Ra"] == 0.0
    {
        return Err(CorrelationError::InvalidGroup {
            axis: "max(Re,Ra)",
            value_bits: 0.0_f64.to_bits(),
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
    let source_region = evaluation_source_region(id, &groups);
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
        source_region,
    })
}

fn evaluation_source_region(
    id: CorrelationId,
    groups: &BTreeMap<String, f64>,
) -> EvaluationSourceRegion {
    match id {
        CorrelationId::RectangularDuctLaminarCwtDevelopingPr072
            if groups["Gz"] < RECTANGULAR_DEVELOPING_FIRST_TABLE_GRAETZ =>
        {
            EvaluationSourceRegion::DeclaredEngineeringBridge
        }
        CorrelationId::RectangularDuctLaminarCwtDevelopingPr072 => {
            EvaluationSourceRegion::TableSliceInterpolation
        }
        _ => EvaluationSourceRegion::CitedFormula,
    }
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
            let graetz = group("Gz");
            3.66 + 0.0668 * graetz / (1.0 + 0.04 * det::pow(graetz, 2.0 / 3.0))
        }
        CorrelationId::RectangularDuctLaminarCwt => rectangular_laminar_cwt(group("aspect_ratio")),
        CorrelationId::RectangularDuctLaminarChf => rectangular_laminar_chf(group("aspect_ratio")),
        CorrelationId::RectangularDuctLaminarCwtDevelopingPr072 => {
            rectangular_laminar_cwt_developing_pr072(group("Gz"))
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
        CorrelationId::ChurchillChuHorizontalCylinder => {
            let denominator = det::pow(1.0 + det::pow(0.559 / group("Pr"), 9.0 / 16.0), 8.0 / 27.0);
            let term = 0.60 + 0.387 * det::pow(group("Ra"), 1.0 / 6.0) / denominator;
            term * term
        }
        CorrelationId::HorizontalPlateUnstableLaminar => 0.54 * det::pow(group("Ra"), 0.25),
        CorrelationId::HorizontalPlateUnstableTurbulent => 0.15 * det::pow(group("Ra"), 1.0 / 3.0),
        CorrelationId::HorizontalPlateStable => 0.27 * det::pow(group("Ra"), 0.25),
        CorrelationId::ManningQureshiAidingVerticalPlate => {
            mixed::AidingVerticalPlate::new(group("Re"), group("Ra"), group("Pr")).nusselt()
        }
        CorrelationId::ChurchillSphere => {
            let denominator = det::pow(1.0 + det::pow(0.469 / group("Pr"), 9.0 / 16.0), 4.0 / 9.0);
            2.0 + 0.589 * det::pow(group("Ra"), 0.25) / denominator
        }
    }
}

fn rectangular_laminar_cwt(aspect_ratio: f64) -> f64 {
    let shape = (-0.548f64)
        .mul_add(aspect_ratio, 2.702)
        .mul_add(aspect_ratio, -5.119)
        .mul_add(aspect_ratio, 4.970)
        .mul_add(aspect_ratio, -2.610)
        .mul_add(aspect_ratio, 1.0);
    7.541 * shape
}

fn rectangular_laminar_chf(aspect_ratio: f64) -> f64 {
    let shape = (-0.1861f64)
        .mul_add(aspect_ratio, 1.0578)
        .mul_add(aspect_ratio, -2.4765)
        .mul_add(aspect_ratio, 3.0853)
        .mul_add(aspect_ratio, -2.0421)
        .mul_add(aspect_ratio, 1.0);
    8.235 * shape
}

fn rectangular_laminar_cwt_developing_pr072(graetz: f64) -> f64 {
    const GRAETZ_KNOTS: [f64; 14] = [
        RECTANGULAR_DEVELOPING_FIRST_TABLE_GRAETZ,
        20.0,
        30.0,
        40.0,
        50.0,
        60.0,
        80.0,
        100.0,
        120.0,
        140.0,
        160.0,
        180.0,
        200.0,
        220.0,
    ];
    const NUSSELT_KNOTS: [f64; 14] = [
        4.20, 4.79, 5.23, 5.61, 5.95, 6.27, 6.88, 7.42, 7.91, 8.37, 8.80, 9.20, 9.60, 10.00,
    ];

    if graetz <= GRAETZ_KNOTS[0] {
        let fully_developed = rectangular_laminar_cwt(0.5);
        return fully_developed + (NUSSELT_KNOTS[0] - fully_developed) * graetz / GRAETZ_KNOTS[0];
    }
    for index in 1..GRAETZ_KNOTS.len() {
        if graetz <= GRAETZ_KNOTS[index] {
            let fraction = (graetz - GRAETZ_KNOTS[index - 1])
                / (GRAETZ_KNOTS[index] - GRAETZ_KNOTS[index - 1]);
            return NUSSELT_KNOTS[index - 1]
                + fraction * (NUSSELT_KNOTS[index] - NUSSELT_KNOTS[index - 1]);
        }
    }
    NUSSELT_KNOTS[NUSSELT_KNOTS.len() - 1]
}

fn card_for(id: CorrelationId) -> CorrelationCard {
    const SHAH_LONDON: SourceProvenance = SourceProvenance {
        citation: "R. K. Shah and A. L. London, Laminar Flow Forced Convection in Ducts, Academic Press, 1978",
        identifier: "ISBN 978-0-12-020051-1",
    };
    const SHAH_LONDON_TABLE_52: SourceProvenance = SourceProvenance {
        citation: "R. K. Shah and A. L. London, Laminar Flow Forced Convection in Ducts, Academic Press, 1978, Chapter VII, Table 52 (Nu_m,T, alpha*=0.5, Pr=0.72)",
        identifier: "doi:10.1016/C2013-0-06152-X; ISBN 978-0-12-020051-1",
    };
    const DITTUS_BOELTER: SourceProvenance = SourceProvenance {
        citation: "F. W. Dittus and L. M. K. Boelter, Heat Transfer in Automobile Radiators of the Tubular Type, University of California Publications in Engineering 2(13), 1930, 443-461",
        identifier: "OCLC 8882310",
    };
    const GNIELINSKI: SourceProvenance = SourceProvenance {
        citation: "V. Gnielinski, New Equations for Heat and Mass Transfer in Turbulent Pipe and Channel Flow, International Chemical Engineering 16(2), 1976, 359-368",
        identifier: "ISSN 0020-6318; volume 16 issue 2 pages 359-368",
    };
    const FLAT_PLATE_LAMINAR: SourceProvenance = SourceProvenance {
        citation: "E. Pohlhausen, Der Waermeaustausch zwischen festen Koerpern und Fluessigkeiten mit kleiner Reibung und kleiner Waermeleitung, ZAMM 1, 1921, 115-121",
        identifier: "ZAMM 1 (1921) 115-121",
    };
    const FLAT_PLATE_MIXED: SourceProvenance = SourceProvenance {
        citation: "HoSung Lee, Thermal Design: Heat Sinks, Thermoelectrics, Heat Pipes, Compact Heat Exchangers, and Solar Cells, 2nd ed., Wiley, 2022, Eq. (1.62)",
        identifier: "doi:10.1002/9781119686040; ISBN 978-1-119-68597-5",
    };
    const CHURCHILL_BERNSTEIN: SourceProvenance = SourceProvenance {
        citation: "S. W. Churchill and M. Bernstein, A Correlating Equation for Forced Convection from Gases and Liquids to a Circular Cylinder in Crossflow, Journal of Heat Transfer 99(2), 1977, 300-306",
        identifier: "doi:10.1115/1.3450685",
    };
    const CHURCHILL_CHU: SourceProvenance = SourceProvenance {
        citation: "S. W. Churchill and H. H. S. Chu, Correlating Equations for Laminar and Turbulent Free Convection from a Vertical Plate, International Journal of Heat and Mass Transfer 18, 1975, 1323-1329",
        identifier: "doi:10.1016/0017-9310(75)90243-4",
    };
    const CHURCHILL_CHU_CYLINDER: SourceProvenance = SourceProvenance {
        citation: "S. W. Churchill and H. H. S. Chu, Correlating Equations for Laminar and Turbulent Free Convection from a Horizontal Cylinder, International Journal of Heat and Mass Transfer 18, 1975, 1049-1053; ASHRAE Handbook Fundamentals 2025, Ch. 4, Table 9 (T9.10)",
        identifier: "doi:10.1016/0017-9310(75)90222-7",
    };
    const HORIZONTAL_PLATE: SourceProvenance = SourceProvenance {
        citation: "ASHRAE Handbook Fundamentals 2025, SI, Ch. 4, Table 9, Eqs. (T9.7)-(T9.9); handbook heat-transfer use of Lloyd-Moran and Goldstein-Sparrow-Jones correlations",
        identifier: "https://handbook.ashrae.org/Handbooks/F25/SI/f25_ch04/f25_ch04_si.aspx#Table9",
    };
    const MANNING_QURESHI: SourceProvenance = SourceProvenance {
        citation: "K. S. Manning and Z. H. Qureshi, Transport Correlations for Laminar Aiding Mixed Convection Over a Vertical Isothermal Surface, Journal of Heat Transfer 116(3), 1994, 777-780, Eqs. (4), (8), (10), (11)",
        identifier: "doi:10.1115/1.2910940",
    };
    const CHURCHILL_SPHERE: SourceProvenance = SourceProvenance {
        citation: "Churchill isolated-sphere relation; COMSOL 6.4 Heat Transfer Theory, External Natural Convection, Sphere, Eq. (4-194), citing Incropera et al., Fundamentals of Heat and Mass Transfer, 6th ed., Eq. (9.35); also ASHRAE Fundamentals 2025, Ch. 4, Table 9 (T9.11), citing Churchill (1990)",
        identifier: "https://doc.comsol.com/6.4/doc/com.comsol.help.heat/heat_ug_theory.07.100.html",
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
        CorrelationId::RectangularDuctLaminarCwtDevelopingPr072 => (
            SHAH_LONDON_TABLE_52,
            ValidityDomain::unconstrained()
                .with("Gz", 0.01, 220.0)
                .with("Pr", 0.72, 0.72)
                .with("Re", 1.0, 2_300.0)
                .with("aspect_ratio", 0.5, 0.5),
            vec![
                "smooth rectangular duct",
                "simultaneously developing velocity and thermal fields",
                "constant wall temperature on all four walls",
                "constant properties with Pr = 0.72",
                "aspect ratio is minor side over major side and equals 0.5",
                "channel-average Nusselt number",
                RECTANGULAR_DEVELOPING_BRIDGE_DISCLOSURE,
            ],
            vec![
                "hydrodynamically developed inlet",
                "axial conduction, buoyancy, or viscous dissipation",
                "strong property variation",
                "rough or interrupted fins",
                "shroud, bypass, or array-interference effects",
            ],
            0.15,
            DiscrepancyBasis::EngineeringAllowance,
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
            FLAT_PLATE_LAMINAR,
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
            FLAT_PLATE_MIXED,
            ValidityDomain::unconstrained()
                .with("Pr", 0.6, 60.0)
                .with("Re", 5.0e5, 1.0e7),
            vec![
                "isothermal smooth flat plate",
                "zero pressure gradient",
                "laminar leading edge followed by attached turbulent flow",
                "transition Reynolds number fixed at Re_x = 5e5 by the leading-edge correction",
            ],
            vec![
                "transition Reynolds number materially different from 5e5",
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
        CorrelationId::ChurchillChuHorizontalCylinder => (
            CHURCHILL_CHU_CYLINDER,
            ValidityDomain::unconstrained()
                .with("Pr", 0.01, 1.0e5)
                .with("Ra", 1.0e-6, 1.0e13),
            vec![
                "isolated long isothermal horizontal circular cylinder; characteristic length is diameter",
                "quiescent single-phase fluid with consistently evaluated film properties",
                "Pr interval is an implementation admission policy, not a source-published validation range",
            ],
            vec![
                "forced or mixed convection",
                "confinement or cylinder arrays",
                "axial end effects",
            ],
            0.20,
            DiscrepancyBasis::EngineeringAllowance,
        ),
        CorrelationId::ManningQureshiAidingVerticalPlate => (
            MANNING_QURESHI,
            ValidityDomain::unconstrained()
                .with("Pr", 0.003, 100.0)
                .with("Ra", 0.0, 1.0e9_f64.next_down())
                .with("Re", 0.0, 5.0e5_f64.next_down()),
            vec![
                "average coefficient over an isothermal vertical plate; characteristic length is height",
                "laminar external stream parallel to and aiding buoyancy; caller must establish alignment",
                "constant-property Boussinesq boundary-layer theory; uniform external stream",
                "strict source laminar bounds Re < 5e5 and Ra < 1e9, with either pure-flow limit admitted",
                "Pr interval [0.003,100] is an implementation policy spanning source comparison points, not a continuously validated experimental range",
                "film-mean properties, if used by an adapter, are a frozen-property approximation and not a source-mandated convention",
                "the two flow drivers cannot both vanish; no zero-Rayleigh conduction intercept is supplied",
            ],
            vec![
                "opposing or oblique flow, separation, transition or turbulence",
                "enclosure confinement, finite-width edge effects or non-isothermal walls",
                "strong property variation or non-Boussinesq density response",
                "creeping-flow or diffusion-dominated regimes outside boundary-layer theory",
            ],
            0.20,
            DiscrepancyBasis::EngineeringAllowance,
        ),
        CorrelationId::ChurchillSphere => (
            CHURCHILL_SPHERE,
            ValidityDomain::unconstrained()
                .with("Pr", 0.7, f64::MAX)
                .with("Ra", 0.0, 1.0e11),
            vec![
                "average coefficient over an isolated isothermal sphere; characteristic length is diameter",
                "quiescent unbounded single-phase Newtonian fluid; transport properties at mean wall/ambient temperature",
                "source bounds Ra <= 1e11 and Pr >= 0.7; no finite upper Pr is published, so its numeric upper bound is the largest finite f64",
                "a beta-based Rayleigh number uses the dimensional adapter's declared Boussinesq property approximation",
                "Nu=2 at Ra=0 is the steady conductive limit to an unbounded fluid and is already included in the coefficient",
            ],
            vec![
                "forced or mixed convection, confinement or particle arrays",
                "non-isothermal sphere or use as a local angular heat-transfer coefficient",
                "strong property variation, transient fluid diffusion or non-Boussinesq density response",
            ],
            0.20,
            DiscrepancyBasis::EngineeringAllowance,
        ),
        CorrelationId::HorizontalPlateUnstableLaminar
        | CorrelationId::HorizontalPlateUnstableTurbulent
        | CorrelationId::HorizontalPlateStable => {
            let (low, high, stability) = match id {
                CorrelationId::HorizontalPlateUnstableLaminar => (
                    2.2e4,
                    8.0e6,
                    "buoyantly unstable face: heated upward or cooled downward for positive expansion",
                ),
                CorrelationId::HorizontalPlateUnstableTurbulent => (
                    8.0e6,
                    1.5e9,
                    "buoyantly unstable face: heated upward or cooled downward for positive expansion",
                ),
                _ => (
                    1.0e5,
                    1.0e10,
                    "buoyantly stable face: heated downward or cooled upward for positive expansion",
                ),
            };
            (
                HORIZONTAL_PLATE,
                ValidityDomain::unconstrained().with("Ra", low, high),
                vec![
                    "isolated isothermal horizontal plate; characteristic length is plate area divided by perimeter",
                    "single-phase quiescent fluid; properties at mean wall/ambient temperature",
                    stability,
                    "handbook empirical fit has no independent Prandtl dependence or published universal Pr range",
                    "inclusive Rayleigh endpoints are the implementation's boundary convention",
                ],
                vec![
                    "forced or mixed convection",
                    "enclosure confinement",
                    "inclined surface",
                    "unvalidated fluid or non-Boussinesq property variation",
                ],
                0.25,
                DiscrepancyBasis::EngineeringAllowance,
            )
        }
    };
    let model_version = match id {
        CorrelationId::RectangularDuctLaminarCwtDevelopingPr072 => "1.1.0",
        _ => "1.0.0",
    };
    CorrelationCard {
        id,
        model: ModelCard::new(
            id.name(),
            model_version,
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
