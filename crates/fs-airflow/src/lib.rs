//! Evidence-bearing fan curves and enclosure-airflow networks.
//!
//! The crate solves the intersection of a monotone piecewise-linear fan
//! pressure curve and a quadratic loss network. Interval Newton certifies the
//! unique root of the declared nominal model. Manufacturer tolerance, loss
//! coefficients, and leakage remain model-form estimates; a numerical root
//! certificate does not promote those physical uncertainties to an enclosure.

pub mod qoi;

use core::fmt;
use fs_convection::CorrelationInputs;
use fs_evidence::{
    Ambition, Evidence, ModelCard, ModelEvidence, NumericalCertificate, ProvenanceHash,
    SensitivitySummary, StatisticalCertificate, ValidityDomain,
};
use fs_ivl::{Interval, RootBox, RootSearchConfig, newton_roots_bounded};
use fs_math::det;
use fs_qty::{
    Area, Density, DynViscosity, Length, Pressure, Qty, QtyAny, Velocity, VolumetricFlowRate,
};
use fs_regime::{RegimeAuditCard, Role, RoleInput, standard_groups};

/// Quadratic pressure-loss resistance in Pa/(m^3/s)^2.
pub type LossResistance = Qty<-7, 1, 0, 0, 0, 0>;

/// Current standard identity for the fan-law scaling relation.
pub const FAN_LAW_SOURCE_IDENTIFIER: &str = "ANSI/AMCA 210-25; ANSI/ASHRAE 51-25";

/// Bibliographic or dataset identity retained by an airflow model.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SourceProvenance {
    /// Human-readable citation.
    pub citation: String,
    /// Stable source, report, or dataset identifier.
    pub identifier: String,
}

impl SourceProvenance {
    /// Construct an explicit source identity.
    #[must_use]
    pub fn new(citation: impl Into<String>, identifier: impl Into<String>) -> Self {
        Self {
            citation: citation.into(),
            identifier: identifier.into(),
        }
    }
}

/// Authority behind a declared relative tolerance.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ToleranceBasis {
    /// A tolerance published by the named manufacturer or retained dataset.
    ManufacturerDeclared,
    /// A caller-declared engineering allowance, not source-published data.
    EngineeringAllowance,
    /// An analytic coefficient with no empirical tolerance contribution.
    Analytic,
}

/// One typed point on a pressure-versus-volume-flow curve.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FanPoint {
    /// Volume flow in m^3/s.
    pub flow: VolumetricFlowRate,
    /// Static pressure rise in Pa.
    pub pressure: Pressure,
}

impl FanPoint {
    /// Construct a curve point from coherent-SI typed values.
    #[must_use]
    pub const fn new(flow: VolumetricFlowRate, pressure: Pressure) -> Self {
        Self { flow, pressure }
    }
}

/// Validated monotone fan data and its physical authority.
#[derive(Debug, Clone, PartialEq)]
pub struct FanCurve {
    name: String,
    points: Vec<FanPoint>,
    source: SourceProvenance,
    pressure_tolerance_rel: f64,
    tolerance_basis: ToleranceBasis,
    admissible_min_flow: VolumetricFlowRate,
    speed_ratio_domain: (f64, f64),
}

impl FanCurve {
    /// Validate typed fan data.
    ///
    /// Flow must increase strictly and pressure must not increase. The interval
    /// below `admissible_min_flow` is an explicit stall/refusal region.
    ///
    /// # Errors
    /// Returns a structured [`AirflowError`] for malformed data, tolerance, or
    /// validity bounds.
    pub fn new(
        name: impl Into<String>,
        points: Vec<FanPoint>,
        source: SourceProvenance,
        pressure_tolerance_rel: f64,
        tolerance_basis: ToleranceBasis,
        admissible_min_flow: VolumetricFlowRate,
        speed_ratio_domain: (f64, f64),
    ) -> Result<Self, AirflowError> {
        let name = name.into();
        if points.len() < 2 {
            return Err(AirflowError::TooFewFanPoints {
                count: points.len(),
            });
        }
        if !(pressure_tolerance_rel.is_finite() && (0.0..1.0).contains(&pressure_tolerance_rel)) {
            return Err(AirflowError::InvalidTolerance {
                context: "fan pressure",
                value_bits: pressure_tolerance_rel.to_bits(),
            });
        }
        let (speed_low, speed_high) = speed_ratio_domain;
        if !(speed_low.is_finite()
            && speed_high.is_finite()
            && speed_low > 0.0
            && speed_low <= 1.0
            && speed_high >= 1.0)
        {
            return Err(AirflowError::InvalidSpeedDomain {
                low_bits: speed_low.to_bits(),
                high_bits: speed_high.to_bits(),
            });
        }
        for (index, point) in points.iter().enumerate() {
            let q = point.flow.value();
            let p = point.pressure.value();
            if !(q.is_finite() && q >= 0.0 && p.is_finite() && p >= 0.0) {
                return Err(AirflowError::InvalidFanPoint {
                    index,
                    flow_bits: q.to_bits(),
                    pressure_bits: p.to_bits(),
                });
            }
        }
        for (index, pair) in points.windows(2).enumerate() {
            if pair[1].flow.value() <= pair[0].flow.value() {
                return Err(AirflowError::NonMonotoneFlow {
                    right_index: index + 1,
                });
            }
            if pair[1].pressure.value() > pair[0].pressure.value() {
                return Err(AirflowError::PressureRise {
                    right_index: index + 1,
                });
            }
        }
        let admissible = admissible_min_flow.value();
        let q_min = points[0].flow.value();
        let q_max = points.last().expect("length checked").flow.value();
        if !(admissible.is_finite() && admissible >= q_min && admissible < q_max) {
            return Err(AirflowError::InvalidStallBoundary {
                value_bits: admissible.to_bits(),
                curve_low_bits: q_min.to_bits(),
                curve_high_bits: q_max.to_bits(),
            });
        }
        Ok(Self {
            name,
            points,
            source,
            pressure_tolerance_rel,
            tolerance_basis,
            admissible_min_flow,
            speed_ratio_domain,
        })
    }

    /// Stable curve name.
    #[must_use]
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Source attached to the pressure-flow data.
    #[must_use]
    pub fn source(&self) -> &SourceProvenance {
        &self.source
    }

    /// Declared relative pressure tolerance.
    #[must_use]
    pub const fn pressure_tolerance_rel(&self) -> f64 {
        self.pressure_tolerance_rel
    }

    /// Authority behind the pressure tolerance.
    #[must_use]
    pub const fn tolerance_basis(&self) -> ToleranceBasis {
        self.tolerance_basis
    }

    fn domain(&self) -> (f64, f64) {
        (
            self.admissible_min_flow.value(),
            self.points.last().expect("length checked").flow.value(),
        )
    }

    fn pressure_interval(&self, flow: Interval) -> Interval {
        let mut output: Option<Interval> = None;
        for pair in self.points.windows(2) {
            let q0 = pair[0].flow.value();
            let q1 = pair[1].flow.value();
            let Some(local) = flow.intersect(Interval::new(q0, q1)) else {
                continue;
            };
            let p0 = Interval::point(pair[0].pressure.value());
            let p1 = Interval::point(pair[1].pressure.value());
            let slope = (p1 - p0) / (Interval::point(q1) - Interval::point(q0));
            let image = p0 + slope * (local - Interval::point(q0));
            output = Some(output.map_or(image, |current| current.hull(image)));
        }
        output.expect("validated solver domain overlaps the fan curve")
    }

    fn slope_interval(&self, flow: Interval) -> Interval {
        let mut output: Option<Interval> = None;
        for pair in self.points.windows(2) {
            let q0 = pair[0].flow.value();
            let q1 = pair[1].flow.value();
            if flow.intersect(Interval::new(q0, q1)).is_none() {
                continue;
            }
            let slope = (Interval::point(pair[1].pressure.value())
                - Interval::point(pair[0].pressure.value()))
                / (Interval::point(q1) - Interval::point(q0));
            output = Some(output.map_or(slope, |current| current.hull(slope)));
        }
        output.expect("validated solver domain overlaps the fan curve")
    }

    /// Shared model card used by final operating-envelope admission.
    ///
    /// Its flow axis is the per-fan base-curve flow before bank arrangement
    /// and fan-law speed scaling; `speed_ratio` is the declared fan-law ratio.
    #[must_use]
    pub fn model_card(&self) -> ModelCard {
        let basis = match self.tolerance_basis {
            ToleranceBasis::ManufacturerDeclared => "manufacturer-declared pressure tolerance",
            ToleranceBasis::EngineeringAllowance => {
                "caller-declared engineering pressure allowance"
            }
            ToleranceBasis::Analytic => "analytic pressure curve",
        };
        ModelCard::new(
            format!("airflow.fan.{}", self.name),
            "1",
            Ambition::Solid,
            vec![
                "monotone piecewise-linear interpolation".to_string(),
                basis.to_string(),
            ],
            ValidityDomain::unconstrained()
                .with("flow_m3_s", self.domain().0, self.domain().1)
                .with(
                    "speed_ratio",
                    self.speed_ratio_domain.0,
                    self.speed_ratio_domain.1,
                ),
            vec![
                "operation below the declared stall boundary".to_string(),
                "system effects not represented by the retained curve".to_string(),
            ],
            self.pressure_tolerance_rel,
        )
    }
}

/// Arrangement for a bank of identical fans.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FanArrangement {
    /// Equal total flow passes through every fan; pressure rises add.
    Series,
    /// Total flow divides equally; every fan sees the same pressure rise.
    Parallel,
}

/// One or more identical fans at a declared speed ratio.
#[derive(Debug, Clone, PartialEq)]
pub struct FanBank {
    curve: FanCurve,
    count: usize,
    arrangement: FanArrangement,
    speed_ratio: f64,
}

impl FanBank {
    /// Construct a series or parallel identical-fan bank.
    ///
    /// # Errors
    /// Refuses zero fans or a speed ratio outside the curve's declared
    /// fan-law scaling domain.
    pub fn new(
        curve: FanCurve,
        count: usize,
        arrangement: FanArrangement,
        speed_ratio: f64,
    ) -> Result<Self, AirflowError> {
        if count == 0 {
            return Err(AirflowError::EmptyFanBank);
        }
        let (low, high) = curve.speed_ratio_domain;
        if !(speed_ratio.is_finite() && speed_ratio >= low && speed_ratio <= high) {
            return Err(AirflowError::SpeedOutOfDomain {
                speed_ratio_bits: speed_ratio.to_bits(),
                low_bits: low.to_bits(),
                high_bits: high.to_bits(),
            });
        }
        Ok(Self {
            curve,
            count,
            arrangement,
            speed_ratio,
        })
    }

    fn flow_factor(&self) -> f64 {
        match self.arrangement {
            FanArrangement::Series => self.speed_ratio,
            FanArrangement::Parallel => self.speed_ratio * self.count as f64,
        }
    }

    fn pressure_factor(&self) -> f64 {
        let speed_squared = self.speed_ratio * self.speed_ratio;
        match self.arrangement {
            FanArrangement::Series => speed_squared * self.count as f64,
            FanArrangement::Parallel => speed_squared,
        }
    }

    fn domain(&self) -> Interval {
        let (low, high) = self.curve.domain();
        Interval::new(low * self.flow_factor(), high * self.flow_factor())
    }

    fn pressure_interval(&self, total_flow: Interval, pressure_scale: f64) -> Interval {
        let base_flow = total_flow / Interval::point(self.flow_factor());
        self.curve.pressure_interval(base_flow)
            * Interval::point(self.pressure_factor() * pressure_scale)
    }

    fn derivative_interval(&self, total_flow: Interval, pressure_scale: f64) -> Interval {
        let base_flow = total_flow / Interval::point(self.flow_factor());
        self.curve.slope_interval(base_flow)
            * Interval::point(self.pressure_factor() * pressure_scale)
            / Interval::point(self.flow_factor())
    }

    /// Interpolated nominal pressure rise at a typed total flow.
    ///
    /// # Errors
    /// Refuses the stall region or flow beyond the retained curve.
    pub fn pressure_at(&self, total_flow: VolumetricFlowRate) -> Result<Pressure, AirflowError> {
        let q = total_flow.value();
        let domain = self.domain();
        if q < domain.lo() {
            return Err(AirflowError::StallRegion {
                intersection_flow_bits: q.to_bits(),
                admissible_min_bits: domain.lo().to_bits(),
            });
        }
        if !(q.is_finite() && q <= domain.hi()) {
            return Err(AirflowError::FlowOutOfDomain {
                flow_bits: q.to_bits(),
                low_bits: domain.lo().to_bits(),
                high_bits: domain.hi().to_bits(),
            });
        }
        Ok(Pressure::new(
            self.pressure_interval(Interval::point(q), 1.0).midpoint(),
        ))
    }
}

/// One quadratic pressure-loss element.
#[derive(Debug, Clone, PartialEq)]
pub struct LossElement {
    /// Stable branch or component name.
    pub name: String,
    /// Nominal resistance in Pa/(m^3/s)^2.
    pub resistance: LossResistance,
    /// Declared relative resistance uncertainty.
    pub uncertainty_rel: f64,
    /// Authority for the coefficient and uncertainty.
    pub source: SourceProvenance,
    /// Authority classification for the uncertainty.
    pub tolerance_basis: ToleranceBasis,
    /// Explicit owner-declared validity for final product admission.
    ///
    /// Absence means that no validated operating domain was supplied. It must
    /// remain cardless rather than being widened to an unconstrained claim.
    regime_validity: Option<ValidityDomain>,
}

impl LossElement {
    /// Construct and validate one loss element.
    ///
    /// # Errors
    /// Refuses an empty name, non-positive resistance, or invalid uncertainty.
    pub fn new(
        name: impl Into<String>,
        resistance: LossResistance,
        uncertainty_rel: f64,
        source: SourceProvenance,
        tolerance_basis: ToleranceBasis,
    ) -> Result<Self, AirflowError> {
        let name = name.into();
        if name.trim().is_empty() {
            return Err(AirflowError::EmptyElementName);
        }
        let value = resistance.value();
        if !(value.is_finite() && value > 0.0) {
            return Err(AirflowError::InvalidResistance {
                element: name,
                value_bits: value.to_bits(),
            });
        }
        if !(uncertainty_rel.is_finite() && (0.0..1.0).contains(&uncertainty_rel)) {
            return Err(AirflowError::InvalidTolerance {
                context: "loss resistance",
                value_bits: uncertainty_rel.to_bits(),
            });
        }
        Ok(Self {
            name,
            resistance,
            uncertainty_rel,
            source,
            tolerance_basis,
            regime_validity: None,
        })
    }

    /// Attach the source-owned operating domain for final product admission.
    ///
    /// The resulting card identity is `airflow.loss.<name>` and its version is
    /// the element's retained [`SourceProvenance::identifier`]. An
    /// unconstrained, unusable, or anonymously sourced domain refuses rather
    /// than creating an apparently validated loss card.
    ///
    /// # Errors
    ///
    /// Returns a structured [`AirflowError`] when the projected identity or
    /// any validity bound cannot support deterministic product admission.
    pub fn with_regime_validity(mut self, validity: ValidityDomain) -> Result<Self, AirflowError> {
        let card_name = format!("airflow.loss.{}", self.name);
        if card_name.trim() != card_name {
            return Err(AirflowError::InvalidLossRegimeIdentity {
                element: self.name,
                field: "card name",
            });
        }
        if self.source.identifier.is_empty()
            || self.source.identifier.trim() != self.source.identifier
        {
            return Err(AirflowError::InvalidLossRegimeIdentity {
                element: self.name,
                field: "source identifier",
            });
        }
        if validity.bounds().is_empty() {
            return Err(AirflowError::EmptyLossRegimeDomain { element: self.name });
        }
        for (axis, &(lo, hi)) in validity.bounds() {
            if axis.is_empty()
                || axis.trim() != axis
                || !lo.is_finite()
                || !hi.is_finite()
                || lo > hi
            {
                return Err(AirflowError::InvalidLossRegimeBound {
                    element: self.name,
                    axis: axis.clone(),
                    low_bits: lo.to_bits(),
                    high_bits: hi.to_bits(),
                });
            }
        }
        self.regime_validity = Some(validity);
        Ok(self)
    }

    /// Project this element's exact owner identity and validity for auditing.
    ///
    /// `None` is an intentional no-claim boundary: the coefficient may still
    /// participate in an estimated solve, but it has no validated-domain card.
    #[must_use]
    pub fn regime_audit_card(&self) -> Option<RegimeAuditCard> {
        self.regime_validity.as_ref().map(|validity| {
            RegimeAuditCard::new(
                format!("airflow.loss.{}", self.name),
                self.source.identifier.clone(),
                validity.clone(),
            )
        })
    }
}

/// Explicit leakage path, kept distinct so a network cannot silently omit it.
#[derive(Debug, Clone, PartialEq)]
pub struct LeakageElement(LossElement);

impl LeakageElement {
    /// Mark a validated loss element as the enclosure leakage path.
    #[must_use]
    pub const fn new(element: LossElement) -> Self {
        Self(element)
    }
}

/// Recursive series/parallel quadratic loss network.
#[derive(Debug, Clone, PartialEq)]
pub enum LossNetwork {
    /// One terminal flow path.
    Element(LossElement),
    /// Elements carrying equal flow, whose pressure losses add.
    Series(Vec<LossNetwork>),
    /// Elements sharing pressure drop, whose flows add.
    Parallel(Vec<LossNetwork>),
}

impl LossNetwork {
    /// Create a non-empty series group.
    ///
    /// # Errors
    /// Refuses an empty group.
    pub fn series(children: Vec<Self>) -> Result<Self, AirflowError> {
        if children.is_empty() {
            Err(AirflowError::EmptyNetworkGroup { kind: "series" })
        } else {
            Ok(Self::Series(children))
        }
    }

    /// Create a non-empty parallel group.
    ///
    /// # Errors
    /// Refuses an empty group.
    pub fn parallel(children: Vec<Self>) -> Result<Self, AirflowError> {
        if children.is_empty() {
            Err(AirflowError::EmptyNetworkGroup { kind: "parallel" })
        } else {
            Ok(Self::Parallel(children))
        }
    }

    /// Nominal equivalent quadratic resistance.
    #[must_use]
    pub fn equivalent_resistance(&self) -> LossResistance {
        LossResistance::new(self.equivalent_scalar(ResistanceBound::Nominal))
    }

    /// Every explicitly validated terminal loss card in traversal order.
    ///
    /// Duplicate identities are retained so the final `fs-regime` audit can
    /// refuse ambiguity instead of silently selecting or deduplicating cards.
    #[must_use]
    pub fn regime_audit_cards(&self) -> Vec<RegimeAuditCard> {
        let mut cards = Vec::new();
        self.collect_regime_audit_cards(&mut cards);
        cards
    }

    fn equivalent_scalar(&self, bound: ResistanceBound) -> f64 {
        match self {
            Self::Element(element) => {
                let factor = match bound {
                    ResistanceBound::Low => 1.0 - element.uncertainty_rel,
                    ResistanceBound::Nominal => 1.0,
                    ResistanceBound::High => 1.0 + element.uncertainty_rel,
                };
                element.resistance.value() * factor
            }
            Self::Series(children) => children
                .iter()
                .map(|child| child.equivalent_scalar(bound))
                .sum(),
            Self::Parallel(children) => {
                let conductance: f64 = children
                    .iter()
                    .map(|child| 1.0 / det::sqrt(child.equivalent_scalar(bound)))
                    .sum();
                1.0 / (conductance * conductance)
            }
        }
    }

    fn allocate(&self, incoming: f64, out: &mut Vec<(String, f64)>) {
        match self {
            Self::Element(element) => out.push((element.name.clone(), incoming)),
            Self::Series(children) => {
                for child in children {
                    child.allocate(incoming, out);
                }
            }
            Self::Parallel(children) => {
                let weights: Vec<f64> = children
                    .iter()
                    .map(|child| 1.0 / det::sqrt(child.equivalent_scalar(ResistanceBound::Nominal)))
                    .collect();
                let total_weight: f64 = weights.iter().sum();
                for (child, weight) in children.iter().zip(weights) {
                    child.allocate(incoming * weight / total_weight, out);
                }
            }
        }
    }

    fn collect_regime_audit_cards(&self, cards: &mut Vec<RegimeAuditCard>) {
        match self {
            Self::Element(element) => {
                if let Some(card) = element.regime_audit_card() {
                    cards.push(card);
                }
            }
            Self::Series(children) | Self::Parallel(children) => {
                for child in children {
                    child.collect_regime_audit_cards(cards);
                }
            }
        }
    }

    fn collect_model(&self, cards: &mut Vec<String>, assumptions: &mut Vec<String>, rel: &mut f64) {
        match self {
            Self::Element(element) => {
                cards.push(format!("airflow.loss.{}", element.name));
                assumptions.push(format!(
                    "{} uses a quadratic pressure-loss coefficient ({:?})",
                    element.name, element.tolerance_basis
                ));
                *rel += element.uncertainty_rel;
            }
            Self::Series(children) | Self::Parallel(children) => {
                for child in children {
                    child.collect_model(cards, assumptions, rel);
                }
            }
        }
    }
}

#[derive(Debug, Clone, Copy)]
enum ResistanceBound {
    Low,
    Nominal,
    High,
}

/// Enclosure network with one mandatory, explicit leakage branch.
#[derive(Debug, Clone, PartialEq)]
pub struct EnclosureNetwork {
    full: LossNetwork,
    leakage_name: String,
}

impl EnclosureNetwork {
    /// Put the explicit leakage element in parallel with the primary network.
    #[must_use]
    pub fn new(primary: LossNetwork, leakage: LeakageElement) -> Self {
        let leakage_name = leakage.0.name.clone();
        Self {
            full: LossNetwork::Parallel(vec![primary, LossNetwork::Element(leakage.0)]),
            leakage_name,
        }
    }

    /// Nominal equivalent resistance including leakage.
    #[must_use]
    pub fn equivalent_resistance(&self) -> LossResistance {
        self.full.equivalent_resistance()
    }

    /// Every explicitly validated loss-element card consumed by this network.
    #[must_use]
    pub fn regime_audit_cards(&self) -> Vec<RegimeAuditCard> {
        self.full.regime_audit_cards()
    }
}

/// Interval-Newton certificate for the unique nominal-model root.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CertifiedOperatingBracket {
    /// Outward-rounded root interval in m^3/s.
    pub flow: Interval,
    /// Number of interval boxes evaluated.
    pub boxes_examined: usize,
}

/// Evidence-bearing terminal branch flow.
#[derive(Debug, Clone, PartialEq)]
pub struct BranchFlow {
    /// Terminal loss-element name.
    pub path: String,
    /// Flow estimate with the same model authority as the operating point.
    pub flow: Evidence<VolumetricFlowRate>,
    /// True for the network's mandatory leakage path.
    pub leakage: bool,
}

/// Solved nominal intersection plus explicitly weaker physical uncertainty.
#[derive(Debug, Clone, PartialEq)]
pub struct OperatingPoint {
    /// Total volume-flow estimate.
    pub flow: Evidence<VolumetricFlowRate>,
    /// System pressure-drop estimate.
    pub pressure: Evidence<Pressure>,
    /// Unique numerical root certificate for the nominal declared model.
    pub nominal_root: CertifiedOperatingBracket,
    /// Per-terminal nominal flow splits with uncertainty attached.
    pub branches: Vec<BranchFlow>,
    /// Nominal leakage fraction of total flow.
    pub leakage_fraction: f64,
}

impl OperatingPoint {
    /// Produce the typed velocity/Reynolds input handed to `fs-convection`.
    ///
    /// # Errors
    /// Refuses an unknown branch or non-positive/non-finite geometry, fluid
    /// properties, or Prandtl number.
    pub fn correlation_handoff(
        &self,
        branch: &str,
        area: Area,
        density: Density,
        dynamic_viscosity: DynViscosity,
        hydraulic_diameter: Length,
        prandtl: f64,
    ) -> Result<CorrelationHandoff, AirflowError> {
        let branch_flow = self
            .branches
            .iter()
            .find(|candidate| candidate.path == branch)
            .ok_or_else(|| AirflowError::UnknownBranch {
                branch: branch.to_string(),
            })?;
        for (field, value) in [
            ("flow area", area.value()),
            ("density", density.value()),
            ("dynamic viscosity", dynamic_viscosity.value()),
            ("hydraulic diameter", hydraulic_diameter.value()),
            ("Prandtl number", prandtl),
        ] {
            if !(value.is_finite() && value > 0.0) {
                return Err(AirflowError::InvalidHandoffInput {
                    field,
                    value_bits: value.to_bits(),
                });
            }
        }
        let velocity_value = branch_flow.flow.value.value() / area.value();
        let inputs = [
            RoleInput {
                role: Role::Density,
                qty: QtyAny::new(density.value(), Density::DIMS),
            },
            RoleInput {
                role: Role::Velocity,
                qty: QtyAny::new(velocity_value, Velocity::DIMS),
            },
            RoleInput {
                role: Role::Length,
                qty: QtyAny::new(hydraulic_diameter.value(), Length::DIMS),
            },
            RoleInput {
                role: Role::DynViscosity,
                qty: QtyAny::new(dynamic_viscosity.value(), DynViscosity::DIMS),
            },
        ];
        let groups = standard_groups(&inputs).map_err(|error| AirflowError::Regime {
            detail: error.to_string(),
        })?;
        let reynolds = groups
            .iter()
            .find(|group| group.name == "Re")
            .expect("all Reynolds roles supplied")
            .value;
        let provenance = ProvenanceHash::chain(
            "airflow-correlation-handoff-v1",
            &[branch_flow.flow.provenance],
        );
        let numerical = NumericalCertificate::estimate(
            branch_flow.flow.numerical.lo / area.value(),
            branch_flow.flow.numerical.hi / area.value(),
        );
        let velocity = Evidence {
            value: Velocity::new(velocity_value),
            qoi: velocity_value,
            numerical,
            statistical: StatisticalCertificate::None,
            model: branch_flow.flow.model.clone(),
            sensitivity: branch_flow.flow.sensitivity.clone(),
            provenance,
            adjoint_ref: None,
        };
        Ok(CorrelationHandoff {
            branch_flow: branch_flow.flow.clone(),
            velocity,
            reynolds,
            correlation_inputs: CorrelationInputs::forced(reynolds, prandtl),
        })
    }
}

/// Typed bridge from a solved branch to forced-convection correlations.
#[derive(Debug, Clone, PartialEq)]
pub struct CorrelationHandoff {
    /// Evidence-bearing branch volume flow.
    pub branch_flow: Evidence<VolumetricFlowRate>,
    /// Evidence-bearing mean branch velocity.
    pub velocity: Evidence<Velocity>,
    /// Reynolds number computed by role-tagged `fs-regime` dimensional logic.
    pub reynolds: f64,
    /// Forced-convection Re/Pr input accepted by `fs-convection`.
    pub correlation_inputs: CorrelationInputs,
}

/// Deterministically solve a fan-bank/network operating point.
///
/// # Errors
/// Refuses a stall-side or absent intersection, incomplete interval search,
/// or a root that interval Newton cannot prove unique.
pub fn solve_operating_point(
    fan: &FanBank,
    network: &EnclosureNetwork,
) -> Result<OperatingPoint, AirflowError> {
    let domain = fan.domain();
    let resistance = network.full.equivalent_scalar(ResistanceBound::Nominal);
    let low_residual = scalar_residual(fan, resistance, domain.lo(), 1.0);
    if low_residual < 0.0 {
        return Err(AirflowError::StallRegion {
            intersection_flow_bits: approximate_unconstrained_root(fan, resistance, domain.lo())
                .to_bits(),
            admissible_min_bits: domain.lo().to_bits(),
        });
    }
    let high_residual = scalar_residual(fan, resistance, domain.hi(), 1.0);
    if high_residual > 0.0 {
        return Err(AirflowError::NoIntersection {
            low_residual_bits: low_residual.to_bits(),
            high_residual_bits: high_residual.to_bits(),
        });
    }
    let nominal = certified_root(fan, resistance, 1.0)?;
    let fan_rel = fan.curve.pressure_tolerance_rel;
    let resistance_low = network.full.equivalent_scalar(ResistanceBound::Low);
    let resistance_high = network.full.equivalent_scalar(ResistanceBound::High);
    let low_flow = certified_root(fan, resistance_high, 1.0 - fan_rel)?
        .flow
        .lo();
    let high_flow = certified_root(fan, resistance_low, 1.0 + fan_rel)?
        .flow
        .hi();
    let nominal_flow = nominal.flow.midpoint();

    let mut cards = vec![fan.curve.model_card().name];
    let mut assumptions = vec![format!(
        "{} identical fans in {:?} at speed ratio {} using {}",
        fan.count, fan.arrangement, fan.speed_ratio, FAN_LAW_SOURCE_IDENTIFIER
    )];
    let mut network_rel = 0.0;
    network
        .full
        .collect_model(&mut cards, &mut assumptions, &mut network_rel);
    cards.sort_unstable();
    cards.dedup();
    assumptions.sort_unstable();
    assumptions.dedup();
    let model = ModelEvidence {
        cards,
        assumptions,
        validity: ValidityDomain::unconstrained()
            .with("flow_m3_s", domain.lo(), domain.hi())
            .with(
                "speed_ratio",
                fan.curve.speed_ratio_domain.0,
                fan.curve.speed_ratio_domain.1,
            ),
        discrepancy_rel: fan_rel + network_rel,
        in_domain: true,
    };
    let provenance = operating_provenance(fan, network);
    let flow = Evidence {
        value: VolumetricFlowRate::new(nominal_flow),
        qoi: nominal_flow,
        numerical: NumericalCertificate::estimate(low_flow, high_flow),
        statistical: StatisticalCertificate::None,
        model: model.clone(),
        sensitivity: SensitivitySummary::default(),
        provenance,
        adjoint_ref: None,
    };
    let pressure_value = resistance * nominal_flow * nominal_flow;
    let pressure_low = resistance_low * low_flow * low_flow;
    let pressure_high = resistance_high * high_flow * high_flow;
    let pressure = Evidence {
        value: Pressure::new(pressure_value),
        qoi: pressure_value,
        numerical: NumericalCertificate::estimate(pressure_low, pressure_high),
        statistical: StatisticalCertificate::None,
        model: model.clone(),
        sensitivity: SensitivitySummary::default(),
        provenance: ProvenanceHash::chain("airflow-operating-pressure-v1", &[provenance]),
        adjoint_ref: None,
    };

    let mut allocated = Vec::new();
    network.full.allocate(nominal_flow, &mut allocated);
    let combined_rel = (fan_rel + network_rel).min(0.99);
    let branches: Vec<BranchFlow> = allocated
        .into_iter()
        .map(|(path, value)| {
            let branch_provenance =
                ProvenanceHash::chain(&format!("airflow-branch-{path}"), &[provenance]);
            BranchFlow {
                leakage: path == network.leakage_name,
                path,
                flow: Evidence {
                    value: VolumetricFlowRate::new(value),
                    qoi: value,
                    numerical: NumericalCertificate::estimate(
                        value * (1.0 - combined_rel),
                        value * (1.0 + combined_rel),
                    ),
                    statistical: StatisticalCertificate::None,
                    model: model.clone(),
                    sensitivity: SensitivitySummary::default(),
                    provenance: branch_provenance,
                    adjoint_ref: None,
                },
            }
        })
        .collect();
    let leakage_flow: f64 = branches
        .iter()
        .filter(|branch| branch.leakage)
        .map(|branch| branch.flow.value.value())
        .sum();
    Ok(OperatingPoint {
        flow,
        pressure,
        nominal_root: nominal,
        branches,
        leakage_fraction: leakage_flow / nominal_flow,
    })
}

fn certified_root(
    fan: &FanBank,
    resistance: f64,
    pressure_scale: f64,
) -> Result<CertifiedOperatingBracket, AirflowError> {
    let domain = fan.domain();
    let f = |flow: Interval| {
        fan.pressure_interval(flow, pressure_scale) - Interval::point(resistance) * flow * flow
    };
    let fp = |flow: Interval| {
        fan.derivative_interval(flow, pressure_scale) - Interval::point(2.0 * resistance) * flow
    };
    let report = newton_roots_bounded(
        &f,
        &fp,
        domain,
        RootSearchConfig {
            min_width: (domain.width() * 1.0e-12).max(1.0e-14),
            max_boxes: 65_536,
        },
    )
    .map_err(|error| AirflowError::RootSearch {
        detail: error.to_string(),
    })?;
    if !report.complete {
        return Err(AirflowError::UncertifiedOperatingPoint {
            boxes_examined: report.boxes_examined,
            roots: report.roots.len(),
            possible: report
                .roots
                .iter()
                .filter(|root| !root.is_certified())
                .count(),
        });
    }
    let certified: Vec<Interval> = report
        .roots
        .iter()
        .filter_map(|root| match root {
            RootBox::Certified(interval) => Some(*interval),
            RootBox::Possible(_) => None,
        })
        .collect();
    let possible = report.roots.len() - certified.len();
    if certified.len() != 1 || possible != 0 {
        return Err(AirflowError::UncertifiedOperatingPoint {
            boxes_examined: report.boxes_examined,
            roots: report.roots.len(),
            possible,
        });
    }
    Ok(CertifiedOperatingBracket {
        flow: certified[0],
        boxes_examined: report.boxes_examined,
    })
}

fn scalar_residual(fan: &FanBank, resistance: f64, q: f64, pressure_scale: f64) -> f64 {
    fan.pressure_interval(Interval::point(q), pressure_scale)
        .midpoint()
        - resistance * q * q
}

fn approximate_unconstrained_root(fan: &FanBank, resistance: f64, lower: f64) -> f64 {
    let pressure = fan
        .pressure_interval(Interval::point(lower), 1.0)
        .midpoint();
    det::sqrt((pressure / resistance).max(0.0))
}

fn operating_provenance(fan: &FanBank, network: &EnclosureNetwork) -> ProvenanceHash {
    let mut bytes = Vec::new();
    bytes.extend_from_slice(b"org.frankensim.fs-airflow.operating-point.v1\0");
    push_string_identity(&mut bytes, &fan.curve.name);
    bytes.extend_from_slice(&(fan.curve.points.len() as u64).to_le_bytes());
    for point in &fan.curve.points {
        bytes.extend_from_slice(&point.flow.value().to_bits().to_le_bytes());
        bytes.extend_from_slice(&point.pressure.value().to_bits().to_le_bytes());
    }
    push_source_identity(&mut bytes, &fan.curve.source);
    bytes.extend_from_slice(&fan.curve.pressure_tolerance_rel.to_bits().to_le_bytes());
    push_tolerance_identity(&mut bytes, fan.curve.tolerance_basis);
    bytes.extend_from_slice(
        &fan.curve
            .admissible_min_flow
            .value()
            .to_bits()
            .to_le_bytes(),
    );
    bytes.extend_from_slice(&fan.curve.speed_ratio_domain.0.to_bits().to_le_bytes());
    bytes.extend_from_slice(&fan.curve.speed_ratio_domain.1.to_bits().to_le_bytes());
    bytes.extend_from_slice(&(fan.count as u64).to_le_bytes());
    bytes.push(match fan.arrangement {
        FanArrangement::Series => 0,
        FanArrangement::Parallel => 1,
    });
    bytes.extend_from_slice(&fan.speed_ratio.to_bits().to_le_bytes());
    push_network_identity(&mut bytes, &network.full);
    push_string_identity(&mut bytes, &network.leakage_name);
    ProvenanceHash::of_bytes(&bytes)
}

fn push_network_identity(bytes: &mut Vec<u8>, network: &LossNetwork) {
    match network {
        LossNetwork::Element(element) => {
            bytes.push(0);
            push_string_identity(bytes, &element.name);
            bytes.extend_from_slice(&element.resistance.value().to_bits().to_le_bytes());
            bytes.extend_from_slice(&element.uncertainty_rel.to_bits().to_le_bytes());
            push_source_identity(bytes, &element.source);
            push_tolerance_identity(bytes, element.tolerance_basis);
        }
        LossNetwork::Series(children) => {
            bytes.push(1);
            bytes.extend_from_slice(&(children.len() as u64).to_le_bytes());
            for child in children {
                push_network_identity(bytes, child);
            }
        }
        LossNetwork::Parallel(children) => {
            bytes.push(2);
            bytes.extend_from_slice(&(children.len() as u64).to_le_bytes());
            for child in children {
                push_network_identity(bytes, child);
            }
        }
    }
}

fn push_source_identity(bytes: &mut Vec<u8>, source: &SourceProvenance) {
    push_string_identity(bytes, &source.citation);
    push_string_identity(bytes, &source.identifier);
}

fn push_tolerance_identity(bytes: &mut Vec<u8>, basis: ToleranceBasis) {
    bytes.push(match basis {
        ToleranceBasis::ManufacturerDeclared => 0,
        ToleranceBasis::EngineeringAllowance => 1,
        ToleranceBasis::Analytic => 2,
    });
}

fn push_string_identity(bytes: &mut Vec<u8>, value: &str) {
    bytes.extend_from_slice(&(value.len() as u64).to_le_bytes());
    bytes.extend_from_slice(value.as_bytes());
}

/// Structured airflow-model refusal.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AirflowError {
    /// At least two fan points are required.
    TooFewFanPoints {
        /// Number of points supplied by the caller.
        count: usize,
    },
    /// A fan point has negative or non-finite data.
    InvalidFanPoint {
        /// Zero-based index of the rejected point.
        index: usize,
        /// Raw IEEE-754 bits of the flow coordinate.
        flow_bits: u64,
        /// Raw IEEE-754 bits of the pressure coordinate.
        pressure_bits: u64,
    },
    /// Fan flow coordinates do not increase strictly.
    NonMonotoneFlow {
        /// Zero-based index of the right-hand point in the rejected pair.
        right_index: usize,
    },
    /// Pressure increases with flow.
    PressureRise {
        /// Zero-based index of the right-hand point in the rejected pair.
        right_index: usize,
    },
    /// A declared relative tolerance is not in `[0, 1)`.
    InvalidTolerance {
        /// Model component whose tolerance was rejected.
        context: &'static str,
        /// Raw IEEE-754 bits of the rejected tolerance.
        value_bits: u64,
    },
    /// Fan-law scaling domain is malformed or excludes the reference speed.
    InvalidSpeedDomain {
        /// Raw IEEE-754 bits of the proposed lower bound.
        low_bits: u64,
        /// Raw IEEE-754 bits of the proposed upper bound.
        high_bits: u64,
    },
    /// Requested speed ratio is outside the declared fan-law domain.
    SpeedOutOfDomain {
        /// Raw IEEE-754 bits of the requested speed ratio.
        speed_ratio_bits: u64,
        /// Raw IEEE-754 bits of the admitted lower bound.
        low_bits: u64,
        /// Raw IEEE-754 bits of the admitted upper bound.
        high_bits: u64,
    },
    /// Stall boundary is outside the retained curve.
    InvalidStallBoundary {
        /// Raw IEEE-754 bits of the proposed stall boundary.
        value_bits: u64,
        /// Raw IEEE-754 bits of the curve's lower flow bound.
        curve_low_bits: u64,
        /// Raw IEEE-754 bits of the curve's upper flow bound.
        curve_high_bits: u64,
    },
    /// Zero fans do not form a bank.
    EmptyFanBank,
    /// Flow falls outside the retained curve.
    FlowOutOfDomain {
        /// Raw IEEE-754 bits of the requested flow.
        flow_bits: u64,
        /// Raw IEEE-754 bits of the curve's lower flow bound.
        low_bits: u64,
        /// Raw IEEE-754 bits of the curve's upper flow bound.
        high_bits: u64,
    },
    /// The operating point lies in the declared non-admissible stall region.
    StallRegion {
        /// Raw IEEE-754 bits of the rejected intersection flow.
        intersection_flow_bits: u64,
        /// Raw IEEE-754 bits of the minimum admissible flow.
        admissible_min_bits: u64,
    },
    /// Loss-element name is empty.
    EmptyElementName,
    /// Quadratic loss resistance is not finite and positive.
    InvalidResistance {
        /// Name of the rejected loss element.
        element: String,
        /// Raw IEEE-754 bits of the rejected resistance.
        value_bits: u64,
    },
    /// Loss-card name or retained source-version identity is malformed.
    InvalidLossRegimeIdentity {
        /// Name of the rejected loss element.
        element: String,
        /// Identity field that failed validation.
        field: &'static str,
    },
    /// A loss validity declaration contains no constrained axes.
    EmptyLossRegimeDomain {
        /// Name of the rejected loss element.
        element: String,
    },
    /// A loss validity axis or its finite interval is unusable.
    InvalidLossRegimeBound {
        /// Name of the rejected loss element.
        element: String,
        /// Axis carrying the unusable identity or interval.
        axis: String,
        /// Raw IEEE-754 bits of the lower bound.
        low_bits: u64,
        /// Raw IEEE-754 bits of the upper bound.
        high_bits: u64,
    },
    /// A series or parallel node has no children.
    EmptyNetworkGroup {
        /// Network composition kind that had no children.
        kind: &'static str,
    },
    /// The admissible fan curve does not cross the system curve.
    NoIntersection {
        /// Raw IEEE-754 bits of the residual at the lower search bound.
        low_residual_bits: u64,
        /// Raw IEEE-754 bits of the residual at the upper search bound.
        high_residual_bits: u64,
    },
    /// Interval root-search configuration failed.
    RootSearch {
        /// Diagnostic returned by the interval root-search layer.
        detail: String,
    },
    /// The bounded search did not prove exactly one root.
    UncertifiedOperatingPoint {
        /// Number of interval boxes examined by the bounded search.
        boxes_examined: usize,
        /// Number of boxes certified to contain a unique root.
        roots: usize,
        /// Number of unresolved possible-root boxes.
        possible: usize,
    },
    /// Requested terminal branch does not exist.
    UnknownBranch {
        /// Requested terminal branch name.
        branch: String,
    },
    /// A handoff quantity is non-finite or non-positive.
    InvalidHandoffInput {
        /// Name of the rejected handoff field.
        field: &'static str,
        /// Raw IEEE-754 bits of the rejected value.
        value_bits: u64,
    },
    /// Role-tagged Reynolds construction refused.
    Regime {
        /// Diagnostic returned by the regime-classification layer.
        detail: String,
    },
}

impl fmt::Display for AirflowError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::TooFewFanPoints { count } => {
                write!(
                    formatter,
                    "fan curve needs at least two points; received {count}"
                )
            }
            Self::InvalidFanPoint { index, .. } => {
                write!(
                    formatter,
                    "fan point {index} must contain finite non-negative SI data"
                )
            }
            Self::NonMonotoneFlow { right_index } => write!(
                formatter,
                "fan flow must increase strictly at point {right_index}"
            ),
            Self::PressureRise { right_index } => write!(
                formatter,
                "fan pressure must not increase at point {right_index}"
            ),
            Self::InvalidTolerance { context, .. } => {
                write!(
                    formatter,
                    "{context} relative tolerance must be finite in [0, 1)"
                )
            }
            Self::InvalidSpeedDomain { .. } => write!(
                formatter,
                "speed-ratio domain must be finite, positive, ordered, and contain 1"
            ),
            Self::SpeedOutOfDomain { .. } => {
                write!(
                    formatter,
                    "requested speed ratio is outside the declared fan-law domain"
                )
            }
            Self::InvalidStallBoundary { .. } => write!(
                formatter,
                "stall boundary must lie inside the retained fan-flow span"
            ),
            Self::EmptyFanBank => write!(formatter, "a fan bank must contain at least one fan"),
            Self::FlowOutOfDomain { .. } => {
                write!(formatter, "flow is outside the retained fan curve")
            }
            Self::StallRegion { .. } => write!(
                formatter,
                "operating point is below the declared stall boundary; extrapolation refused"
            ),
            Self::EmptyElementName => write!(formatter, "loss-element name must not be empty"),
            Self::InvalidResistance { element, .. } => write!(
                formatter,
                "loss resistance for {element:?} must be finite and positive"
            ),
            Self::InvalidLossRegimeIdentity { element, field } => write!(
                formatter,
                "loss element {element:?} cannot publish a regime card with invalid {field}"
            ),
            Self::EmptyLossRegimeDomain { element } => write!(
                formatter,
                "loss element {element:?} cannot publish an unconstrained regime card"
            ),
            Self::InvalidLossRegimeBound { element, axis, .. } => write!(
                formatter,
                "loss element {element:?} has an unusable regime bound for axis {axis:?}"
            ),
            Self::EmptyNetworkGroup { kind } => {
                write!(formatter, "{kind} loss-network group must not be empty")
            }
            Self::NoIntersection { .. } => {
                write!(
                    formatter,
                    "fan and system curves do not cross in the admissible domain"
                )
            }
            Self::RootSearch { detail } => {
                write!(formatter, "interval root search refused: {detail}")
            }
            Self::UncertifiedOperatingPoint {
                boxes_examined,
                roots,
                possible,
            } => write!(
                formatter,
                "operating point not certified: {roots} root boxes ({possible} possible) after {boxes_examined} evaluations"
            ),
            Self::UnknownBranch { branch } => {
                write!(formatter, "loss-network branch {branch:?} does not exist")
            }
            Self::InvalidHandoffInput { field, .. } => {
                write!(
                    formatter,
                    "{field} must be finite and positive for correlation handoff"
                )
            }
            Self::Regime { detail } => write!(formatter, "Reynolds construction refused: {detail}"),
        }
    }
}

impl std::error::Error for AirflowError {}
