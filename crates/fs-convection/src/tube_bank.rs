//! Mean convection for bare circular tubes in crossflow, including array effects.
//!
//! Implements the `1000 <= Re_max <= 200000` rows of Zukauskas (1987),
//! as reproduced in Flownex CS-052, Tables 1--2 and Eqs. (2), (7)--(9):
//! <https://support.flownex.com/hc/en-us/articles/28978533677980>.
//! Bulk transport belongs at the mean fluid temperature; only the Prandtl
//! correction uses fluid at the wall temperature. This is a bank-average
//! coefficient, not an isolated-cylinder or individual-row correlation.

use super::{
    Ambition, BTreeMap, CorrelationError, CorrelationRobinBoundary, DomainViolation, Evidence,
    HeatTransferCoefficient, Length, ModelCard, ModelEvidence, NumericalCertificate,
    ProvenanceHash, SensitivitySummary, SourceProvenance, StatisticalCertificate, Temperature,
    ThermalConductivity, ValidityDomain, domain_violations,
};

mod kernel;
pub use kernel::TubeBankLayout;

/// Actual tube diameter, center-to-center pitches and streamwise row count.
/// The bounded pitch subset is an implementation admission, not a validation claim.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TubeBankGeometry {
    layout: TubeBankLayout,
    rows: u32,
    diameter: Length,
    transverse_over_diameter: f64,
    longitudinal_over_diameter: f64,
    gap_velocity_ratio: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum TubeBankError {
    InvalidInput {
        field: &'static str,
        value_bits: u64,
    },
    OutOfDomain {
        layout: TubeBankLayout,
        violations: Vec<DomainViolation>,
    },
    Correlation(CorrelationError),
}
impl core::fmt::Display for TubeBankError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "tube-bank convection: {self:?}")
    }
}
impl std::error::Error for TubeBankError {}
impl From<CorrelationError> for TubeBankError {
    fn from(error: CorrelationError) -> Self {
        Self::Correlation(error)
    }
}
fn positive(field: &'static str, value: f64) -> Result<f64, TubeBankError> {
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(TubeBankError::InvalidInput {
            field,
            value_bits: value.to_bits(),
        })
    }
}

impl TubeBankGeometry {
    /// Admit bare tubes in a uniform perpendicular approach flow without bypass.
    /// Transverse pitch/D is in [1.25, 2.5]; longitudinal pitch/D is in [1, 2.5].
    /// Touching inline tubes refuse; staggered pitch ratios above two refuse.
    /// Staggered banks require at least two rows: one row has no longitudinal
    /// neighbor spacing to justify the correlation's pitch-ratio dependence.
    pub fn try_new(
        layout: TubeBankLayout,
        rows: u32,
        diameter: Length,
        transverse_pitch: Length,
        longitudinal_pitch: Length,
    ) -> Result<Self, TubeBankError> {
        let d = positive("outside diameter [m]", diameter.value())?;
        let transverse = positive("transverse pitch [m]", transverse_pitch.value())? / d;
        let longitudinal = positive("longitudinal pitch [m]", longitudinal_pitch.value())? / d;
        let groups = BTreeMap::from([
            ("rows".into(), f64::from(rows)),
            ("S_T_over_D".into(), transverse),
            ("S_L_over_D".into(), longitudinal),
            ("S_T_over_S_L".into(), transverse / longitudinal),
        ]);
        let minimum_rows = if layout == TubeBankLayout::Staggered {
            2.0
        } else {
            1.0
        };
        let mut domain = ValidityDomain::unconstrained()
            .with("rows", minimum_rows, f64::from(u32::MAX))
            .with("S_T_over_D", 1.25, 2.5)
            .with("S_L_over_D", 1.0, 2.5);
        if layout == TubeBankLayout::Staggered {
            domain = domain.with("S_T_over_S_L", 0.5, 2.0);
        }
        let violations = domain_violations(&domain, &groups);
        if !violations.is_empty() {
            return Err(TubeBankError::OutOfDomain { layout, violations });
        }
        if layout == TubeBankLayout::Inline && longitudinal <= 1.0 {
            return Err(TubeBankError::InvalidInput {
                field: "inline longitudinal clearance",
                value_bits: (longitudinal - 1.0).to_bits(),
            });
        }
        let gap_velocity_ratio = positive(
            "gap velocity ratio",
            kernel::gap_velocity_ratio(layout, transverse, longitudinal),
        )?;
        Ok(Self {
            layout,
            rows,
            diameter,
            transverse_over_diameter: transverse,
            longitudinal_over_diameter: longitudinal,
            gap_velocity_ratio,
        })
    }
    pub const fn layout(self) -> TubeBankLayout {
        self.layout
    }
    pub const fn rows(self) -> u32 {
        self.rows
    }
    pub const fn diameter(self) -> Length {
        self.diameter
    }
    /// Maximum inter-tube speed divided by upstream unobstructed speed.
    pub const fn gap_velocity_ratio(self) -> f64 {
        self.gap_velocity_ratio
    }
    pub fn maximum_gap_speed(self, approach_speed_m_s: f64) -> Result<f64, TubeBankError> {
        positive(
            "maximum gap speed [m/s]",
            positive("approach speed [m/s]", approach_speed_m_s)? * self.gap_velocity_ratio,
        )
    }
}

pub const SOURCE: SourceProvenance = SourceProvenance {
    citation: "A. Zukauskas (1987), Heat Transfer from Tubes in Cross Flow, Handbook of Single Phase Convective Heat Transfer; Flownex CS-052, Tables 1-2, Eqs. (2), (7)-(9)",
    identifier: "https://support.flownex.com/hc/en-us/articles/28978533677980",
};

/// Source-bearing bank mean. Model discrepancy is explicitly unquantified.
#[derive(Clone, Debug, PartialEq)]
pub struct TubeBankEvaluation {
    geometry: TubeBankGeometry,
    card: ModelCard,
    groups: BTreeMap<String, f64>,
    evidence: Evidence<f64>,
}
impl TubeBankEvaluation {
    pub const fn geometry(&self) -> TubeBankGeometry {
        self.geometry
    }
    pub const fn card(&self) -> &ModelCard {
        &self.card
    }
    pub const fn source(&self) -> &SourceProvenance {
        &SOURCE
    }
    pub const fn groups(&self) -> &BTreeMap<String, f64> {
        &self.groups
    }
    pub const fn evidence(&self) -> &Evidence<f64> {
        &self.evidence
    }

    /// Use bulk-fluid conductivity and the ORIGINAL outside diameter.
    pub fn heat_transfer_coefficient(
        &self,
        conductivity: ThermalConductivity,
    ) -> Result<Evidence<HeatTransferCoefficient>, TubeBankError> {
        let k = positive("bulk conductivity [W/(m K)]", conductivity.value())?;
        let d = self.geometry.diameter.value();
        let h = positive(
            "heat-transfer coefficient [W/(m2 K)]",
            self.evidence.value * k / d,
        )?;
        let mut dimensional = [0_u8; 16];
        dimensional[..8].copy_from_slice(&k.to_bits().to_le_bytes());
        dimensional[8..].copy_from_slice(&d.to_bits().to_le_bytes());
        Ok(Evidence {
            value: HeatTransferCoefficient::new(h),
            qoi: h,
            numerical: NumericalCertificate::estimate(h, h),
            statistical: StatisticalCertificate::None,
            model: self.evidence.model.clone(),
            sensitivity: SensitivitySummary::default(),
            provenance: ProvenanceHash::chain(
                "tube-bank-nu-to-h-v1",
                &[
                    self.evidence.provenance,
                    ProvenanceHash::of_bytes(&dimensional),
                ],
            ),
            adjoint_ref: None,
        })
    }

    /// Lower the SAME bank mean into the existing conduction boundary.
    /// Use an appropriate effective bulk reference (e.g. a log-mean difference
    /// for a warming stream). Applying inlet temperature to the entire bank
    /// is not a finite-capacity-fluid solution and overpredicts heat transfer.
    pub fn robin_boundary(
        &self,
        conductivity: ThermalConductivity,
        reference_temperature: Temperature,
    ) -> Result<CorrelationRobinBoundary, TubeBankError> {
        let temperature = positive(
            "Robin reference temperature [K]",
            reference_temperature.value(),
        )?;
        let coefficient = self.heat_transfer_coefficient(conductivity)?;
        let boundary = fs_conduction::ThermalBc::robin(coefficient.value.value(), temperature)
            .map_err(|_| CorrelationError::NonFiniteResult {
                stage: "tube-bank Robin boundary",
                value_bits: coefficient.value.value().to_bits(),
            })?;
        Ok(CorrelationRobinBoundary {
            coefficient,
            reference_temperature,
            boundary,
        })
    }
}

/// Evaluate the bank from Reynolds number based on APPROACH speed and outside
/// diameter. The geometry computes Re_max; callers must not supply gap Re here.
/// Pr and k use the mean bulk temperature, Pr_surface the same fluid at the wall.
///
/// Re_max in [1e3,2e5] and Pr/Pr_surface in [0.7,500] are the admitted subset;
/// lower/high-Re branches, fins, short-pin end effects, roughness, buoyancy,
/// bypass and pressure drop are absent. Source row interpolation, the finite
/// pitch box and the wall-Pr box are explicit implementation policies, not a
/// jointly validated data envelope. No model error percentage is fabricated.
pub fn evaluate(
    geometry: TubeBankGeometry,
    approach_reynolds: f64,
    prandtl: f64,
    surface_prandtl: f64,
) -> Result<TubeBankEvaluation, TubeBankError> {
    positive("approach Reynolds number", approach_reynolds)?;
    positive("bulk Prandtl number", prandtl)?;
    positive("surface Prandtl number", surface_prandtl)?;
    let reynolds = positive(
        "maximum-gap Reynolds number",
        approach_reynolds * geometry.gap_velocity_ratio,
    )?;
    let groups = BTreeMap::from([
        ("Re_max".into(), reynolds),
        ("Pr".into(), prandtl),
        ("Pr_surface".into(), surface_prandtl),
        ("rows".into(), f64::from(geometry.rows)),
        ("S_T_over_D".into(), geometry.transverse_over_diameter),
        ("S_L_over_D".into(), geometry.longitudinal_over_diameter),
        (
            "S_T_over_S_L".into(),
            geometry.transverse_over_diameter / geometry.longitudinal_over_diameter,
        ),
    ]);
    let name = match geometry.layout {
        TubeBankLayout::Inline => "convection.zukauskas-inline-tube-bank",
        TubeBankLayout::Staggered => "convection.zukauskas-staggered-tube-bank",
    };
    let minimum_rows = if geometry.layout == TubeBankLayout::Staggered {
        2.0
    } else {
        1.0
    };
    let mut validity = ValidityDomain::unconstrained()
        .with("Re_max", 1.0e3, 2.0e5)
        .with("Pr", 0.7, 500.0)
        .with("Pr_surface", 0.7, 500.0)
        .with("rows", minimum_rows, f64::from(u32::MAX))
        .with("S_T_over_D", 1.25, 2.5)
        .with("S_L_over_D", 1.0, 2.5);
    if geometry.layout == TubeBankLayout::Staggered {
        validity = validity.with("S_T_over_S_L", 0.5, 2.0);
    }
    let card = ModelCard::new(name, "1.0.0", Ambition::Solid, [
        "bare smooth circular isothermal tubes in uniform perpendicular single-phase crossflow without bypass",
        "mean bulk properties except Pr_surface, evaluated for the SAME fluid at the wall temperature",
        "outside diameter and maximum inter-tube speed define Re_max; bank mean, not individual-row transfer",
        "finite-row correction linearly interpolates Flownex CS-052 Table 2; it is unity at 16 or more rows",
        "a staggered bank requires at least two rows to give its longitudinal pitch physical meaning",
        "pitch and surface-Pr boxes are implementation admissions, not a jointly validated source envelope",
        "model-form discrepancy is unquantified; numerical formula checks do not establish experimental validation",
    ].into_iter().map(str::to_owned).collect(),
        validity,
        ["finned or short-pin arrays, bypass, roughness, natural/mixed convection, compressibility and phase change",
         "lower or higher Reynolds regimes; local tube/row coefficients; pressure loss; finite-fluid energy balance",
         "using inlet temperature on the whole array despite appreciable fluid heating or cooling"]
            .into_iter().map(str::to_owned).collect(), f64::INFINITY);
    let violations = domain_violations(&card.validity, &groups);
    if !violations.is_empty() {
        return Err(TubeBankError::OutOfDomain {
            layout: geometry.layout,
            violations,
        });
    }
    let nu = positive(
        "bank Nusselt number",
        kernel::mean_nusselt(
            geometry.layout,
            geometry.rows,
            reynolds,
            prandtl,
            surface_prandtl,
            geometry.transverse_over_diameter / geometry.longitudinal_over_diameter,
        ),
    )?;
    let mut bytes = Vec::with_capacity(128);
    bytes.extend_from_slice(name.as_bytes());
    for (axis, value) in &groups {
        bytes.extend_from_slice(axis.as_bytes());
        bytes.extend_from_slice(&value.to_bits().to_le_bytes());
    }
    let evidence = Evidence {
        value: nu,
        qoi: nu,
        numerical: NumericalCertificate::estimate(nu, nu),
        statistical: StatisticalCertificate::None,
        model: ModelEvidence::from_card(&card, &groups),
        sensitivity: SensitivitySummary::default(),
        provenance: ProvenanceHash::chain("tube-bank-mean-v1", &[ProvenanceHash::of_bytes(&bytes)]),
        adjoint_ref: None,
    };
    Ok(TubeBankEvaluation {
        geometry,
        card,
        groups,
        evidence,
    })
}
