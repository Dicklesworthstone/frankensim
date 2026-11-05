//! Specimen-bound JK drop-weight impact response and its energy inverse.
//!
//! `t10 = A * (1 - exp(-b * Ecs))`, where t10 is the mass fraction passing
//! one tenth of the original particle size. A is stored as a fraction (not the
//! percent convention in laboratory reports); b has units kg/J and Ecs J/kg.
//! Source: JKTech / SMC Testing, *An addition to reporting A and b values —
//! The SAG Circuit Specific Energy (SCSE)*, p.1,
//! <https://www.smctesting.com/documents/Reporting_JK_Drop_Weight_and_SMC_Test_Results.pdf>.
//!
//! A and b must come from the same scalar material observation. A*b alone does
//! not identify this curve. The source and its receipts remain attached to the
//! model. Its explicit energy interval is a caller declaration, not inferred
//! laboratory support. No full product-size distribution, abrasion ta, selection
//! rate, liberation, Bond index or AG/SAG circuit power is inferred from t10.

use fs_blake3::ContentHash;
use fs_matdb::{
    MatDbError, MaterialAnswer, NormalizedPack, ObservationDataset, PropertyClaim, PropertyKey,
    PropertyValue, QueryPoint, SelectionPolicy,
};
use fs_math::det;
use fs_qty::{Dimensionless, Qty};

use super::SpecificGrindingEnergy;

/// Reciprocal specific impact energy in coherent SI kg/J.
pub type ImpactEnergySensitivity = Qty<-2, 0, 2, 0, 0, 0>;

/// Independent parameters; their names disambiguate percent and fraction A.
#[must_use]
pub fn parameter_keys() -> [PropertyKey; 2] {
    [
        PropertyKey::new(
            "jk-drop-weight-asymptotic-t10-fraction",
            Dimensionless::DIMS,
        ),
        PropertyKey::new(
            "jk-drop-weight-impact-energy-sensitivity",
            ImpactEnergySensitivity::DIMS,
        ),
    ]
}

/// Source refusal, invalid physical data or an unsupported energy request.
#[derive(Debug, Clone, PartialEq)]
pub enum DropWeightError {
    /// Original MatDB query or receipt failure.
    Material(MatDbError),
    /// Malformed calibration, invalid target or unrepresentable arithmetic.
    Invalid(&'static str),
    /// Requested or inverted energy lies outside the declared interval.
    OutsideEnergyDomain,
}
impl core::fmt::Display for DropWeightError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "drop-weight impact: {self:?}")
    }
}
impl std::error::Error for DropWeightError {}
impl From<MatDbError> for DropWeightError {
    fn from(error: MatDbError) -> Self {
        Self::Material(error)
    }
}
fn require(condition: bool, message: &'static str) -> Result<(), DropWeightError> {
    if condition {
        Ok(())
    } else {
        Err(DropWeightError::Invalid(message))
    }
}

/// Finite inclusive interval for the caller's intended impact-energy calculation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ImpactEnergyDomain {
    bounds: [SpecificGrindingEnergy; 2],
}
impl ImpactEnergyDomain {
    /// Supply J/kg bounds. A nonnegative ordered interval is required.
    pub fn new(bounds: [SpecificGrindingEnergy; 2]) -> Result<Self, DropWeightError> {
        let [low, high] = bounds.map(SpecificGrindingEnergy::value);
        require(
            low.is_finite() && high.is_finite() && low >= 0.0 && high >= low,
            "impact-energy bounds must be finite, nonnegative and ordered",
        )?;
        Ok(Self { bounds })
    }
    /// Original declared lower and upper energies in J/kg.
    #[must_use]
    pub const fn bounds(self) -> [SpecificGrindingEnergy; 2] {
        self.bounds
    }
    fn check(self, energy: SpecificGrindingEnergy) -> Result<(), DropWeightError> {
        let value = energy.value();
        require(
            value.is_finite() && value >= 0.0,
            "impact energy must be finite and nonnegative",
        )?;
        if value < self.bounds[0].value() || value > self.bounds[1].value() {
            return Err(DropWeightError::OutsideEnergyDomain);
        }
        Ok(())
    }
}

/// Original material evidence retained without replacing uncertainty or caveats.
#[derive(Debug, Clone, PartialEq)]
pub struct DropWeightCalibrationSource {
    pack_id: String,
    pack_hash: ContentHash,
    query: QueryPoint,
    observation: ObservationDataset,
    claims: [PropertyClaim; 2],
    answers: [MaterialAnswer; 2],
}
impl DropWeightCalibrationSource {
    /// Original material-pack identity.
    #[must_use]
    pub fn pack_id(&self) -> &str {
        &self.pack_id
    }
    /// Original complete normalized-pack content hash.
    #[must_use]
    pub const fn pack_hash(&self) -> ContentHash {
        self.pack_hash
    }
    /// Original typed source-selection query.
    #[must_use]
    pub const fn query(&self) -> &QueryPoint {
        &self.query
    }
    /// Original specimen, method, artifact, caveats and source terms.
    #[must_use]
    pub const fn observation(&self) -> &ObservationDataset {
        &self.observation
    }
    /// Scalar A-fraction and b claims, in that order.
    #[must_use]
    pub const fn claims(&self) -> &[PropertyClaim; 2] {
        &self.claims
    }
    /// Selected evidence and verified usage receipts, in parameter order.
    #[must_use]
    pub const fn answers(&self) -> &[MaterialAnswer; 2] {
        &self.answers
    }
}

/// Both sides of the one-tenth-size threshold; neither defines a full PSD.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DropWeightResponse {
    passing: Dimensionless,
    retained: Dimensionless,
}
impl DropWeightResponse {
    /// t10 as a mass fraction, not a percent.
    #[must_use]
    pub const fn passing_one_tenth(self) -> Dimensionless {
        self.passing
    }
    /// Mass fraction retained above the one-tenth-size threshold.
    /// Evaluated independently to preserve a small exponential tail when A=1.
    #[must_use]
    pub const fn retained_one_tenth(self) -> Dimensionless {
        self.retained
    }
}

/// Source-bound A,b response for one original test observation.
#[derive(Debug, Clone, PartialEq)]
pub struct DropWeightCard {
    a: Dimensionless,
    b: ImpactEnergySensitivity,
    domain: ImpactEnergyDomain,
    source: DropWeightCalibrationSource,
}
impl DropWeightCard {
    /// Maximum t10 mass fraction as impact energy tends to infinity.
    #[must_use]
    pub const fn asymptotic_t10(&self) -> Dimensionless {
        self.a
    }
    /// Original normalized b in kg/J, inverse specific energy.
    #[must_use]
    pub const fn energy_sensitivity(&self) -> ImpactEnergySensitivity {
        self.b
    }
    /// Explicit calculation domain, not a new experimental validity claim.
    #[must_use]
    pub const fn domain(&self) -> ImpactEnergyDomain {
        self.domain
    }
    /// The original selected laboratory observation and its evidence.
    #[must_use]
    pub const fn calibration_source(&self) -> &DropWeightCalibrationSource {
        &self.source
    }

    /// Calculate the fines statistic at an admitted specific impact energy.
    /// `expm1` preserves the linear small-energy limit; the exponential complement
    /// is evaluated independently so rounding t10 to one does not delete its tail.
    pub fn response(
        &self,
        energy: SpecificGrindingEnergy,
    ) -> Result<DropWeightResponse, DropWeightError> {
        self.domain.check(energy)?;
        let exponent = self.b.value() * energy.value();
        let passing = self.a.value() * -det::expm1(-exponent);
        let retained = (1.0 - self.a.value()) + self.a.value() * det::exp(-exponent);
        require(
            energy.value() == 0.0 || passing > 0.0,
            "positive impact response underflows binary64",
        )?;
        Ok(DropWeightResponse {
            passing: Dimensionless::new(passing),
            retained: Dimensionless::new(retained),
        })
    }

    /// Invert a reachable t10 fraction to J/kg within the declared domain.
    /// A finite energy cannot reach the asymptote A; t10 >= A refuses instead
    /// of returning an infinite energy or silently capping the requested fines.
    pub fn energy_for_t10(
        &self,
        t10: Dimensionless,
    ) -> Result<SpecificGrindingEnergy, DropWeightError> {
        let target = t10.value();
        let a = self.a.value();
        require(
            target.is_finite() && target >= 0.0 && target < a,
            "target t10 must be finite, nonnegative and strictly below A",
        )?;
        let fraction = target / a;
        let log_ratio = if fraction <= 0.5 {
            // Stable -log1p(-fraction) using the deterministic ln kernel.
            // Correct rounding of 1-fraction loses low bits, restored by the ratio.
            let complement = 1.0 - fraction;
            if complement == 1.0 {
                fraction
            } else {
                -det::ln(complement) * (fraction / (1.0 - complement))
            }
        } else {
            // Subtract before dividing to keep targets immediately below A
            // reachable even when target/A would round to one.
            det::ln(a) - det::ln(a - target)
        };
        let energy = SpecificGrindingEnergy::new(log_ratio / self.b.value());
        require(
            energy.value().is_finite() && (target == 0.0 || energy.value() > 0.0),
            "inverse impact energy is unrepresentable",
        )?;
        self.domain.check(energy)?;
        Ok(energy)
    }
}

/// Resolve both independently reported JKDWT parameters from one original row.
/// Percent A must be normalized to a fraction before insertion; b is kg/J.
/// Missing parameters, mixed observations/support, citation-only or interpolated
/// fields, unsupported query coordinates and inadmissible parameters refuse.
/// The caller supplies the intended energy domain explicitly; source test
/// sizes or discrete energy settings are not widened into a validity rectangle.
pub fn resolve_jk_drop_weight(
    pack: &NormalizedPack,
    point: &QueryPoint,
    policy: SelectionPolicy,
    domain: ImpactEnergyDomain,
) -> Result<DropWeightCard, DropWeightError> {
    let resolve = |key: &PropertyKey| -> Result<_, DropWeightError> {
        let answer = pack.claims().query_typed(key, point, policy)?;
        pack.claims().verify_receipt(&answer.receipt)?;
        let claim =
            pack.claims()
                .claim(answer.receipt.selected)
                .ok_or(DropWeightError::Invalid(
                    "selected impact parameter is absent",
                ))?;
        require(
            matches!(claim.value, PropertyValue::Scalar { .. }) && claim.observations.len() == 1,
            "impact parameters require scalars and one original observation",
        )?;
        require(
            point
                .axes()
                .keys()
                .all(|axis| claim.validity.bounds().contains_key(axis)),
            "query includes conditions unsupported by the original impact test",
        )?;
        Ok((claim.clone(), answer))
    };
    let [a_key, b_key] = parameter_keys();
    let (a_claim, a_answer) = resolve(&a_key)?;
    let (b_claim, b_answer) = resolve(&b_key)?;
    require(
        a_claim.observations == b_claim.observations && a_claim.validity == b_claim.validity,
        "A and b must retain the same original specimen observation and source support",
    )?;
    let a = a_answer.evidence.value.value;
    let b = b_answer.evidence.value.value;
    require(
        a.is_finite() && a > 0.0 && a <= 1.0,
        "A must be a finite mass fraction in (0,1]",
    )?;
    require(
        b.is_finite() && b > 0.0,
        "b must be finite and positive in kg/J",
    )?;
    let observation = pack
        .claims()
        .observation(a_claim.observations[0])
        .ok_or(DropWeightError::Invalid(
            "original impact-test observation is absent",
        ))?
        .clone();
    Ok(DropWeightCard {
        a: Dimensionless::new(a),
        b: ImpactEnergySensitivity::new(b),
        domain,
        source: DropWeightCalibrationSource {
            pack_id: pack.pack_id().to_owned(),
            pack_hash: pack.content_hash(),
            query: point.clone(),
            observation,
            claims: [a_claim, b_claim],
            answers: [a_answer, b_answer],
        },
    })
}

#[cfg(test)]
#[path = "drop_weight/tests.rs"]
mod tests;
