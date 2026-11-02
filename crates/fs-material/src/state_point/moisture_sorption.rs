//! Invertible source isotherms for conservative moisture partitioning.
//!
//! Different materials at the same water activity generally retain different
//! dry-basis moisture contents. This chart supplies that local equilibrium
//! relation and its storage capacity; it does not make the specimen equilibrate
//! instantaneously or prescribe a transport coefficient. At fixed temperature,
//! the existing relative-humidity coordinate is the dilute-vapor activity.

use super::*;
use fs_qty::semantic::MoistureBasis;

const ACTIVITY_AXIS: &str = "relative_humidity";
type Result<T> = core::result::Result<T, MaterialStatePointError>;

fn invalid(quantity: &'static str) -> MaterialStatePointError {
    MaterialStatePointError::InvalidDerived { quantity }
}
fn quantity(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}

/// Prepared, strictly increasing dry-basis isotherm M(a) at fixed temperature.
///
/// The original typed source-line query admits every intervening segment and
/// surface cell. Positive finite slopes give a unique inverse and positive
/// storage in an activity-based transport solve. Flat or decreasing branches,
/// hysteresis and source gaps cannot be silently inverted. Source preparation
/// is bounded by the original 4096-knot cap; every evaluation allocates nothing.
#[derive(Clone, Debug, PartialEq)]
pub struct MonotoneMoistureSorption {
    knots: Box<[(f64, f64)]>,
    capacities: Box<[f64]>,
    source_point: QueryPoint,
    card_identity: ContentHash,
    temperature_kelvin: f64,
    selected_claim: ClaimId,
    identity: ContentHash,
}

impl MonotoneMoistureSorption {
    /// Prepare the selected kg-water/kg-dry-solid property on an activity span.
    ///
    /// Both canonical input coordinates must have their exact static semantics:
    /// absolute `temperature` and `relative_humidity`. The latter is constrained
    /// to [0,1]. All other source context stays fixed. `activity_bounds` must
    /// have positive width; a single point cannot establish invertible storage.
    /// The binding and immutable material card must agree. Preparation neither
    /// changes the selected source nor extrapolates its measured/model domain.
    pub fn from_binding(
        card: &MaterialCard,
        binding: &ResolvedMaterialStatePoint,
        property: &str,
        activity_bounds: [f64; 2],
    ) -> Result<Self> {
        let dry = quantity(QuantityKind::Moisture(MoistureBasis::DryMass));
        let activity = quantity(QuantityKind::Moisture(MoistureBasis::RelativeHumidity));
        let absolute = quantity(QuantityKind::AbsoluteTemperature);
        let point = binding.query();
        if card.content_hash() != binding.card_identity()
            || point.axis_quantities().get(ACTIVITY_AXIS) != Some(&activity)
            || point.axis_quantities().get("temperature") != Some(&absolute)
            || !point.axes().get("temperature").is_some_and(|t| *t > 0.0)
            || !point
                .axes()
                .get(ACTIVITY_AXIS)
                .is_some_and(|a| (0.0..=1.0).contains(a))
            || point.axes().contains_key("T")
            || point.axis_quantities().iter().any(|(axis, q)| {
                (axis != ACTIVITY_AXIS && *q == activity)
                    || (axis != "temperature" && *q == absolute)
            })
            || activity_bounds
                .iter()
                .any(|a| !a.is_finite() || !(0.0..=1.0).contains(a))
            || activity_bounds[0] >= activity_bounds[1]
        {
            return Err(invalid(
                "sorption chart needs its original card, fixed absolute temperature and a positive-width activity interval in [0,1]",
            ));
        }
        let bound = binding
            .property(property)
            .ok_or_else(|| invalid("sorption chart needs a uniquely bound dry-basis property"))?;
        if bound.requirement().quantity() != dry
            || !bound.has_only_owned_and_fixed_external_conditions(ACTIVITY_AXIS, activity)
        {
            return Err(invalid(
                "sorption chart requires a static dry-mass moisture property on the owned activity coordinate",
            ));
        }
        let selected_claim = bound.answer().receipt.selected;
        let claim = card
            .claims()
            .claim(selected_claim)
            .ok_or_else(|| invalid("selected sorption source is absent"))?;
        let lower = point
            .clone()
            .with(ACTIVITY_AXIS, activity_bounds[0])
            .map_err(|_| invalid("sorption chart lower activity"))?;
        let upper = point
            .with(ACTIVITY_AXIS, activity_bounds[1])
            .map_err(|_| invalid("sorption chart upper activity"))?;
        let path = card
            .claims()
            .query_line_typed(
                &claim.key,
                ACTIVITY_AXIS,
                &lower,
                &upper,
                fs_matdb::ClaimSelection::Pinned(selected_claim),
            )
            .map_err(|source| MaterialStatePointError::Support {
                property: property.into(),
                source: Box::new(source),
            })?;
        let mut identity =
            DomainHasher::new("org.frankensim.fs-material.monotone-moisture-sorption.v1");
        identity.update(binding.identity().as_bytes());
        identity.update(selected_claim.0.as_bytes());
        let mut knots: Vec<(f64, f64)> = Vec::with_capacity(path.knots.len());
        let mut capacities = Vec::with_capacity(path.knots.len() - 1);
        for (a, answer) in path.knots {
            let moisture = answer.evidence.value.value;
            if !moisture.is_finite() || moisture < 0.0 {
                return Err(MaterialStatePointError::OutsideConsumerDomain {
                    property: property.into(),
                    value: moisture,
                    admissibility: ScalarAdmissibility::NonNegative,
                });
            }
            if let Some(&(previous_a, previous_m)) = knots.last() {
                let capacity = (moisture - previous_m) / (a - previous_a);
                if !capacity.is_finite() || capacity <= 0.0 {
                    return Err(invalid(
                        "sorption chart needs strictly increasing moisture and finite positive storage on every segment",
                    ));
                }
                capacities.push(capacity);
            }
            identity.update(&a.to_bits().to_le_bytes());
            identity.update(answer.receipt.content_hash().as_bytes());
            knots.push((a, moisture));
        }
        Ok(Self {
            knots: knots.into_boxed_slice(),
            capacities: capacities.into_boxed_slice(),
            source_point: binding.query(),
            card_identity: card.content_hash(),
            temperature_kelvin: binding
                .state_coordinate("temperature")
                .expect("admitted temperature"),
            selected_claim,
            identity: identity.finalize(),
        })
    }

    /// Supported water-activity interval, equal to relative humidity at held T.
    pub fn activity_bounds(&self) -> [f64; 2] {
        [self.knots[0].0, self.knots[self.knots.len() - 1].0]
    }
    /// Supported kg-water/kg-dry-solid interval, with no wet-basis conversion.
    pub fn moisture_bounds(&self) -> [f64; 2] {
        [self.knots[0].1, self.knots[self.knots.len() - 1].1]
    }
    pub fn knot_count(&self) -> usize {
        self.knots.len()
    }
    /// Original typed preparation context. A transport consumer must match all
    /// held coordinates, not only the card and temperature, when combining laws.
    /// Only the canonical relative-humidity coordinate evolves in this chart.
    pub const fn source_point(&self) -> &QueryPoint {
        &self.source_point
    }
    pub const fn card_identity(&self) -> ContentHash {
        self.card_identity
    }
    pub const fn temperature_kelvin(&self) -> f64 {
        self.temperature_kelvin
    }
    pub const fn selected_claim(&self) -> ClaimId {
        self.selected_claim
    }
    /// Binds original card, claim, environment and the complete supported path.
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    fn admit_activity(&self, activity: f64) -> Result<()> {
        let [lower, upper] = self.activity_bounds();
        if !activity.is_finite() || activity < lower || activity > upper {
            Err(invalid(
                "activity lies outside the prepared sorption source path",
            ))
        } else {
            Ok(())
        }
    }

    /// Original piecewise-affine equilibrium moisture M(a), in kg/kg dry solid.
    pub fn moisture(&self, activity: f64) -> Result<f64> {
        self.admit_activity(activity)?;
        let i = self.knots.partition_point(|&(a, _)| a < activity);
        if i == 0 || self.knots[i].0 == activity {
            return Ok(self.knots[i].1);
        }
        let (a, m) = self.knots[i - 1];
        let (b, n) = self.knots[i];
        let fraction = (activity - a) / (b - a);
        let value = if fraction <= 0.5 {
            fraction.mul_add(n - m, m)
        } else {
            (1.0 - fraction).mul_add(m - n, n)
        };
        if value.is_finite() {
            Ok(value)
        } else {
            Err(invalid("unrepresentable source sorption interpolation"))
        }
    }

    /// Unique source-supported activity corresponding to a local moisture value.
    pub fn activity(&self, moisture: f64) -> Result<f64> {
        let [lower, upper] = self.moisture_bounds();
        if !moisture.is_finite() || moisture < lower || moisture > upper {
            return Err(invalid(
                "moisture lies outside the prepared sorption source path",
            ));
        }
        let i = self.knots.partition_point(|&(_, m)| m < moisture);
        if i == 0 || self.knots[i].1 == moisture {
            return Ok(self.knots[i].0);
        }
        let (a, m) = self.knots[i - 1];
        let (b, n) = self.knots[i];
        let fraction = (moisture - m) / (n - m);
        let value = if fraction <= 0.5 {
            fraction.mul_add(b - a, a)
        } else {
            (1.0 - fraction).mul_add(a - b, b)
        };
        if value.is_finite() {
            Ok(value)
        } else {
            Err(invalid(
                "unrepresentable inverse source sorption interpolation",
            ))
        }
    }

    /// Positive local dM/da. At a kink this selects the segment to the right;
    /// the upper endpoint uses the final segment. This is a one-sided derivative,
    /// not a claim that a piecewise-affine table is differentiable at its knots.
    pub fn capacity(&self, activity: f64) -> Result<f64> {
        self.admit_activity(activity)?;
        let right = self.knots.partition_point(|&(a, _)| a <= activity);
        Ok(self.capacities[(right - 1).min(self.capacities.len() - 1)])
    }

    /// Positive secant [M(b)-M(a)]/(b-a), symmetric in its endpoints.
    ///
    /// The interval-weighted source slopes are accumulated directly, avoiding
    /// loss of a small storage change by subtracting large moisture datums.
    /// Evaluation visits only crossed segments and allocates nothing.
    pub fn mean_capacity(&self, a: f64, b: f64) -> Result<f64> {
        self.admit_activity(a)?;
        self.admit_activity(b)?;
        if a == b {
            return self.capacity(a);
        }
        let (lower, upper) = if a < b { (a, b) } else { (b, a) };
        let width = upper - lower;
        let mut x = lower;
        let mut right = self.knots.partition_point(|&(a, _)| a <= lower);
        let (mut result, mut correction) = (0.0, 0.0);
        loop {
            let end = self.knots[right].0.min(upper);
            let term = self.capacities[right - 1] * ((end - x) / width);
            let add = term - correction;
            let next = result + add;
            correction = (next - result) - add;
            result = next;
            if end == upper {
                break;
            }
            x = end;
            right += 1;
        }
        if result.is_finite() && result > 0.0 {
            Ok(result)
        } else {
            Err(invalid("unrepresentable source sorption storage capacity"))
        }
    }

    /// Signed local storage change M(b)-M(a), integrated without a large datum.
    pub fn moisture_change(&self, a: f64, b: f64) -> Result<f64> {
        let value = self.mean_capacity(a, b)? * (b - a);
        if value.is_finite() && (a == b || value != 0.0) {
            Ok(value)
        } else {
            Err(invalid("unrepresentable source sorption moisture change"))
        }
    }
}
