//! Experimental Austin selection and daughter fits on their original size grid.
//!
//! Selection: `S_i = a_T (x_i / 1 mm)^alpha / (1 + (x_i/mu)^lambda)`.
//! Daughters: `B_ij = phi (x_(i-1)/x_j)^gamma + (1-phi)(x_(i-1)/x_j)^beta`,
//! `b_ij = B_ij - B_(i+1),j`, for descending source size classes and `i > j`.
//! In particular `B_(j+1),j = 1`. This discrete index convention is not replaced
//! with `x_i/x_j`, which would discard part of each selected parent's mass.
//! Reference: Petrakis, Stamboliadis and Komnitsas (2017), equations 11–15,
//! <https://doi.org/10.14356/kona.2017007>.
//!
//! Seven scalar parameters and an exact class-size curve must share one source
//! observation and support. The last class is an explicitly unresolved sink:
//! zero selection, all remaining cumulative daughter mass, no sub-grid PSD.
//! No clipping, renormalization, rebinning, ore defaults or plant scale-up occurs.
//! The resulting model uses the existing constant-rate population-balance owner;
//! original claims and receipts remain attached through its prepared passes.

use super::breakage::{
    BreakageColumn, BreakageError, BreakageTestContext, BreakageTransfer, MAX_SIZE_CLASSES,
    SizeClass, SizeRedistributionCard,
};
use super::kinetics::{
    CONTACT_TIME_AXIS, GrindingKineticsError, RATE_PREFIX, SELECTION_BREAKAGE_LAW,
    SelectionBreakageModel,
};
use fs_blake3::ContentHash;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSelection, ConstitutiveModelCard, EnvelopeAnswer, InitialStatePolicy, InterpolationPolicy,
    LawId, LawParameter, MatDbError, NormalizedPack, ObservationDataset, PropertyClaim,
    PropertyKey, PropertySupportError, PropertyValue, QueryPoint, SelectionPolicy,
};
use fs_math::det;
use fs_qty::{Dimensionless, Frequency, Length, QuantitySpec, Time};

/// Source curve coordinate. Exact consecutive indices start at zero, coarse first.
pub const CLASS_INDEX_AXIS: &str = "austin-size-class-index";

/// The scalar parameters, in amplitude, alpha, mu, lambda, phi, gamma, beta order.
#[must_use]
pub fn parameter_keys() -> [PropertyKey; 7] {
    [
        ("austin-selection-amplitude", Frequency::DIMS),
        ("austin-selection-size-exponent", Dimensionless::DIMS),
        ("austin-selection-rolloff-size", Length::DIMS),
        ("austin-selection-rolloff-exponent", Dimensionless::DIMS),
        ("austin-daughter-fines-weight", Dimensionless::DIMS),
        ("austin-daughter-fines-exponent", Dimensionless::DIMS),
        ("austin-daughter-coarse-exponent", Dimensionless::DIMS),
    ]
    .map(|(name, dims)| PropertyKey::new(name, dims))
}

/// Exact source size-class ordinates in metres, including the final sink node.
#[must_use]
pub fn class_size_key() -> PropertyKey {
    PropertyKey::new("austin-size-class-upper-size", Length::DIMS)
}

/// Source or physical admission refusal; original underlying errors survive.
#[derive(Clone, Debug, PartialEq)]
pub enum AustinSourceError {
    /// A typed query or receipt could not be admitted.
    Material(MatDbError),
    /// One original claim does not cover the complete requested model interval.
    Support(PropertySupportError),
    /// Existing conservative redistribution owner refused the derived matrix.
    Breakage(BreakageError),
    /// Existing selection/population-balance owner refused the derived model.
    Kinetics(GrindingKineticsError),
    /// Invalid schema, mixed calibration or unrepresentable coefficient.
    Invalid(&'static str),
}
impl core::fmt::Display for AustinSourceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "Austin grinding source: {self:?}")
    }
}
impl std::error::Error for AustinSourceError {}
impl From<MatDbError> for AustinSourceError {
    fn from(value: MatDbError) -> Self {
        Self::Material(value)
    }
}
impl From<PropertySupportError> for AustinSourceError {
    fn from(value: PropertySupportError) -> Self {
        Self::Support(value)
    }
}
impl From<BreakageError> for AustinSourceError {
    fn from(value: BreakageError) -> Self {
        Self::Breakage(value)
    }
}
impl From<GrindingKineticsError> for AustinSourceError {
    fn from(value: GrindingKineticsError) -> Self {
        Self::Kinetics(value)
    }
}
type Result<T> = core::result::Result<T, AustinSourceError>;
fn require(condition: bool, message: &'static str) -> Result<()> {
    if condition {
        Ok(())
    } else {
        Err(AustinSourceError::Invalid(message))
    }
}

/// Immutable original fit, not a claim that its predictions are measurements.
#[derive(Clone, Debug, PartialEq)]
pub struct AustinCalibrationSource {
    pack_id: String,
    pack_hash: ContentHash,
    conditions: QueryPoint,
    maximum_contact_time: Time,
    observation: ObservationDataset,
    claims: Vec<PropertyClaim>,
    answers: Vec<EnvelopeAnswer>,
}
impl AustinCalibrationSource {
    /// Original material-pack identifier.
    #[must_use]
    pub fn pack_id(&self) -> &str {
        &self.pack_id
    }
    /// Original complete material-pack content identity.
    #[must_use]
    pub const fn pack_hash(&self) -> ContentHash {
        self.pack_hash
    }
    /// Fixed operating conditions used to select this original fit.
    #[must_use]
    pub const fn conditions(&self) -> &QueryPoint {
        &self.conditions
    }
    /// Requested interval [0, maximum], admitted against every original claim.
    #[must_use]
    pub const fn maximum_contact_time(&self) -> Time {
        self.maximum_contact_time
    }
    /// Original test method, specimen, uncertainty caveats and source terms.
    #[must_use]
    pub const fn observation(&self) -> &ObservationDataset {
        &self.observation
    }
    /// Original seven scalar claims followed by the exact class-size curve.
    #[must_use]
    pub fn claims(&self) -> &[PropertyClaim] {
        &self.claims
    }
    /// Verified lower/upper evidence: seven parameters, then each source class.
    #[must_use]
    pub fn answers(&self) -> &[EnvelopeAnswer] {
        &self.answers
    }
}

#[derive(Clone, Copy)]
struct Coefficients {
    amplitude: f64,
    alpha: f64,
    mu: f64,
    lambda: f64,
    phi: f64,
    gamma: f64,
    beta: f64,
}
impl Coefficients {
    fn new(values: [f64; 7]) -> Result<Self> {
        let [amplitude, alpha, mu, lambda, phi, gamma, beta] = values;
        require(
            values.iter().all(|v| v.is_finite())
                && amplitude > 0.0
                && alpha > 0.0
                && mu > 0.0
                && lambda >= 0.0
                && (0.0..=1.0).contains(&phi)
                && gamma > 0.0
                && beta > 0.0,
            "finite positive amplitude, alpha, mu and daughter exponents; lambda >= 0; phi in [0,1] required",
        )?;
        Ok(Self {
            amplitude,
            alpha,
            mu,
            lambda,
            phi,
            gamma,
            beta,
        })
    }
    fn selection(self, diameter: Length) -> Result<Frequency> {
        let numerator = self.amplitude * det::pow(diameter.value() / 0.001, self.alpha);
        let denominator = 1.0 + det::pow(diameter.value() / self.mu, self.lambda);
        let rate = numerator / denominator;
        require(
            numerator.is_finite() && denominator.is_finite() && rate.is_finite() && rate > 0.0,
            "positive source selection rate loses arithmetic range",
        )?;
        Ok(Frequency::new(rate))
    }
    fn cumulative(self, ratio: f64) -> Result<f64> {
        require(
            ratio.is_finite() && ratio > 0.0 && ratio <= 1.0,
            "invalid source daughter size ratio",
        )?;
        if ratio == 1.0 {
            return Ok(1.0);
        }
        let value =
            self.phi * det::pow(ratio, self.gamma) + (1.0 - self.phi) * det::pow(ratio, self.beta);
        require(
            value.is_finite() && value > 0.0 && value <= 1.0,
            "positive cumulative daughter fraction loses range",
        )?;
        Ok(value)
    }
    fn columns(self, classes: &[SizeClass]) -> Result<Vec<BreakageColumn>> {
        let mut columns = Vec::with_capacity(classes.len());
        for (j, parent) in classes.iter().enumerate() {
            let mut transfers = Vec::new();
            let mut cumulative = 1.0;
            for i in j + 1..classes.len() {
                // B_(i+1),j uses x_i/x_j. The final unresolved sink has no
                // smaller represented boundary; it receives all remaining mass.
                let next = if i + 1 == classes.len() {
                    0.0
                } else {
                    self.cumulative(classes[i].diameter.value() / parent.diameter.value())?
                };
                let fraction = cumulative - next;
                require(
                    fraction.is_finite() && fraction > 0.0,
                    "daughter difference is negative or loses positive mass",
                )?;
                transfers.push(BreakageTransfer {
                    output: classes[i].id.clone(),
                    mass_fraction: Dimensionless::new(fraction),
                });
                cumulative = next;
            }
            if j + 1 == classes.len() {
                transfers.push(BreakageTransfer {
                    output: parent.id.clone(),
                    mass_fraction: Dimensionless::new(1.0),
                });
            }
            columns.push(BreakageColumn {
                input: parent.id.clone(),
                transfers,
            });
        }
        Ok(columns)
    }
}

fn point(conditions: &QueryPoint, index: usize, elapsed: Time) -> Result<QueryPoint> {
    Ok(conditions
        .clone()
        .with_quantity(
            CLASS_INDEX_AXIS,
            QuantitySpec::dimensional(Dimensionless::DIMS),
            index as f64,
        )?
        .with_quantity(
            CONTACT_TIME_AXIS,
            QuantitySpec::dimensional(Time::DIMS),
            elapsed.value(),
        )?)
}

fn source_claim<'a>(
    pack: &'a NormalizedPack,
    answer: &EnvelopeAnswer,
    query: &QueryPoint,
) -> Result<&'a PropertyClaim> {
    pack.claims().verify_receipt(&answer.lower.receipt)?;
    pack.claims().verify_receipt(&answer.upper.receipt)?;
    let claim =
        pack.claims()
            .claim(answer.lower.receipt.selected)
            .ok_or(AustinSourceError::Invalid(
                "selected original claim is absent",
            ))?;
    require(
        claim.observations.len() == 1,
        "one original fit observation is required",
    )?;
    require(
        query
            .axes()
            .keys()
            .all(|axis| claim.validity.bounds().contains_key(axis)),
        "unsupported operating conditions",
    )?;
    Ok(claim)
}

fn source_grid(claim: &PropertyClaim) -> Result<Vec<SizeClass>> {
    let PropertyValue::Curve {
        abscissa,
        abscissa_dims,
        knots,
        dims,
    } = &claim.value
    else {
        return Err(AustinSourceError::Invalid(
            "class sizes require an exact source curve",
        ));
    };
    require(
        abscissa == CLASS_INDEX_AXIS
            && *abscissa_dims == Dimensionless::DIMS
            && *dims == Length::DIMS
            && claim.interpolation == InterpolationPolicy::TabulatedOnly
            && (2..=MAX_SIZE_CLASSES).contains(&knots.len()),
        "class sizes require 2..128 exact consecutive source indices",
    )?;
    require(
        knots
            .iter()
            .enumerate()
            .all(|(i, &(x, size))| x == i as f64 && size.is_finite() && size > 0.0)
            && knots.windows(2).all(|pair| pair[0].1 > pair[1].1),
        "source class indices/sizes must be consecutive and strictly descending",
    )?;
    Ok(knots
        .iter()
        .enumerate()
        .map(|(i, &(_, size))| SizeClass {
            id: if i + 1 == knots.len() {
                "fines-sink".into()
            } else {
                format!("size-{i:03}")
            },
            diameter: Length::new(size),
        })
        .collect())
}

/// Resolve one complete fitted source into the existing batch population balance.
///
/// Source class sizes are immutable, coarse first; the final node is a sink.
/// The active class size is the upper-size coordinate in the selection law,
/// not a newly inferred geometric-mean diameter. The sink node is only a
/// representative label for unresolved fines and acquires no predicted rate.
/// Every scalar parameter and each exact grid sample must cover the complete
/// [0, maximum_contact_time] interval at the supplied typed conditions.
/// A broad fitted-model domain is not proof of experimental accuracy throughout it.
pub fn resolve_austin_grinding(
    pack: &NormalizedPack,
    conditions: &QueryPoint,
    maximum_contact_time: Time,
    policy: SelectionPolicy,
) -> Result<SelectionBreakageModel> {
    require(
        maximum_contact_time.value().is_finite() && maximum_contact_time.value() > 0.0,
        "maximum contact time must be finite and positive",
    )?;
    require(
        !conditions.axes().contains_key(CLASS_INDEX_AXIS)
            && !conditions.axes().contains_key(CONTACT_TIME_AXIS),
        "class indices and contact interval are owned by the source resolver",
    )?;
    let lower = point(conditions, 0, Time::new(0.0))?;
    let first_upper = point(conditions, 0, maximum_contact_time)?;
    let first_grid = pack.claims().query_envelope_typed(
        &class_size_key(),
        &lower,
        &first_upper,
        ClaimSelection::Policy(policy),
    )?;
    let grid_claim = source_claim(pack, &first_grid, &lower)?;
    let classes = source_grid(grid_claim)?;
    let upper = point(conditions, classes.len() - 1, maximum_contact_time)?;
    let mut claims = Vec::with_capacity(8);
    let mut answers = Vec::with_capacity(7 + classes.len());
    for key in parameter_keys() {
        let answer = pack.claims().query_envelope_typed(
            &key,
            &lower,
            &upper,
            ClaimSelection::Policy(policy),
        )?;
        let claim = source_claim(pack, &answer, &lower)?;
        require(
            matches!(claim.value, PropertyValue::Scalar { .. }),
            "Austin coefficients require original scalar fit parameters",
        )?;
        require(
            claim.observations == grid_claim.observations && claim.validity == grid_claim.validity,
            "all coefficients and class sizes must retain the same original fit and source support",
        )?;
        claims.push(claim.clone());
        answers.push(answer);
    }
    answers.push(first_grid);
    for i in 1..classes.len() {
        let low = point(conditions, i, Time::new(0.0))?;
        let high = point(conditions, i, maximum_contact_time)?;
        let answer = pack.claims().query_envelope_typed(
            &class_size_key(),
            &low,
            &high,
            ClaimSelection::Policy(policy),
        )?;
        require(
            source_claim(pack, &answer, &low)? == grid_claim,
            "source class grid cannot be stitched from different claims",
        )?;
        answers.push(answer);
    }
    claims.push(grid_claim.clone());
    let observation = pack
        .claims()
        .observation(grid_claim.observations[0])
        .ok_or(AustinSourceError::Invalid(
            "original fit observation is absent",
        ))?
        .clone();
    build_model(
        classes,
        AustinCalibrationSource {
            pack_id: pack.pack_id().to_owned(),
            pack_hash: pack.content_hash(),
            conditions: conditions.clone(),
            maximum_contact_time,
            observation,
            claims,
            answers,
        },
    )
}

fn build_model(
    classes: Vec<SizeClass>,
    source: AustinCalibrationSource,
) -> Result<SelectionBreakageModel> {
    let coefficients = Coefficients::new(core::array::from_fn(|i| {
        source.answers[i].lower.evidence.value.value
    }))?;
    let observation = &source.observation;
    let mut validity = ValidityDomain::unconstrained().with_quantity(
        CONTACT_TIME_AXIS,
        QuantitySpec::dimensional(Time::DIMS),
        0.0,
        source.maximum_contact_time.value(),
    );
    for (axis, &value) in source.conditions.axes() {
        let quantity =
            source
                .conditions
                .axis_quantities()
                .get(axis)
                .ok_or(AustinSourceError::Invalid(
                    "operating conditions must be typed",
                ))?;
        validity = validity.with_quantity(axis, quantity.clone(), value, value);
    }
    let context = BreakageTestContext {
        specimen: observation.specimen.clone(),
        method: format!(
            "computed Austin Eq11-15 on exact source class grid; terminal unresolved-fines sink; original fit: {}",
            observation.method
        ),
        provenance: observation.provenance.clone(),
    };
    let parameters = classes
        .iter()
        .enumerate()
        .map(|(i, class)| {
            let rate = if i + 1 == classes.len() {
                Frequency::new(0.0)
            } else {
                coefficients.selection(class.diameter)?
            };
            Ok((
                format!("{RATE_PREFIX}{}", class.id),
                LawParameter {
                    value: rate.value(),
                    dims: Frequency::DIMS,
                },
            ))
        })
        .collect::<Result<_>>()?;
    let daughters = SizeRedistributionCard::new(
        classes.clone(),
        coefficients.columns(&classes)?,
        context,
        validity.clone(),
    )?;
    let mut model = SelectionBreakageModel::from_card(
        daughters,
        ConstitutiveModelCard {
            law: LawId(SELECTION_BREAKAGE_LAW.into()),
            law_version: 1,
            state_schema_version: 0,
            initial_state: InitialStatePolicy::ZeroInternalState,
            parameters,
            validity,
            sources: vec![source.pack_hash, observation.artifact],
            provenance: observation.provenance.clone(),
        },
    )?;
    model.austin_source = Some(Box::new(source));
    Ok(model)
}

#[cfg(test)]
#[path = "austin/tests.rs"]
mod tests;
