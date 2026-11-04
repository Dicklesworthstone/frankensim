//! Loading-dependent moisture sorption energy on the original source path.
//!
//! The admitted property is the differential isosteric heat
//! `q_st = h_vapor - partial_h_bound`, in J/kg WATER, not J/kg dry solid.
//! Positive values denote exothermic uptake. Its integral over dry-basis
//! moisture is an energy per kg dry skeleton. No heat of condensation or
//! material-specific latent heat is supplied implicitly.
//!
//! For a rigid solid with negligible bound-phase pV, ideal vapor at held T
//! gives `partial_u_bound - u_vapor = R_water*T - q_st`. The finite-reservoir
//! consumer owns that conversion and the complete gas/solid energy balance.
//! Definition: IUPAC Gold Book, https://doi.org/10.1351/goldbook.I03313 .
//! The rigid-solid `delta_U - delta_H ~= RT` relation is stated in the original
//! manual: https://old.iupac.org/reports/2001/colloid_2001/manual_of_s_and_t/node31.html .

use super::*;
use fs_math::dd::Dd;
use fs_qty::semantic::MoistureBasis;

/// Differential vapor-to-bound-water enthalpy release [J/kg water].
pub const MOISTURE_ISOSTERIC_HEAT_PROPERTY: &str = "moisture_isosteric_heat";
pub const MOISTURE_ISOSTERIC_HEAT_DIMS: Dims = Dims([2, 0, -2, 0, 0, 0]);
const AXIS: &str = "moisture_content";
type Result<T> = core::result::Result<T, MaterialStatePointError>;

fn invalid(quantity: &'static str) -> MaterialStatePointError {
    MaterialStatePointError::InvalidDerived { quantity }
}
fn quantity(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}

/// Prepared differential isosteric heat q_st(M) at fixed original conditions.
///
/// The complete source path, including intervening surface cells, is admitted
/// once. Integration uses the source's piecewise-affine interpolation exactly
/// up to arithmetic roundoff. No allocation or source re-selection occurs in
/// evaluation; the original source owner limits preparation to 4096 knots.
/// This chart does not establish consistency with a separately supplied
/// temperature-dependent isotherm, hysteresis, or an entropy-production bound.
#[derive(Clone, Debug, PartialEq)]
pub struct MoistureSorptionEnergy {
    knots: Box<[(f64, f64)]>,
    source_point: QueryPoint,
    card_identity: ContentHash,
    selected_claim: ClaimId,
    identity: ContentHash,
}

impl MoistureSorptionEnergy {
    /// Pin the named differential isosteric heat on an original moisture span.
    ///
    /// Temperature is absolute, moisture is kg water/kg dry solid, and all
    /// other source context remains fixed. Signed finite q_st is permitted:
    /// neither exothermic uptake nor a positive source is invented. The caller
    /// explicitly supplies TOTAL vapor-to-bound differential enthalpy, including
    /// any condensation contribution present in its source; a net heat relative
    /// to liquid water cannot be passed as though it had the same meaning.
    pub fn from_binding(
        card: &MaterialCard,
        binding: &ResolvedMaterialStatePoint,
        property: &str,
        moisture_bounds: [f64; 2],
    ) -> Result<Self> {
        let dry = quantity(QuantityKind::Moisture(MoistureBasis::DryMass));
        let absolute = quantity(QuantityKind::AbsoluteTemperature);
        let point = binding.query();
        if card.content_hash() != binding.card_identity()
            || point.axis_quantities().get(AXIS) != Some(&dry)
            || point.axis_quantities().get("temperature") != Some(&absolute)
            || !point.axes().get("temperature").is_some_and(|t| *t > 0.0)
            || !point.axes().get(AXIS).is_some_and(|m| *m >= 0.0)
            || point.axes().contains_key("T")
            || point.axis_quantities().iter().any(|(axis, q)| {
                (axis != AXIS && *q == dry) || (axis != "temperature" && *q == absolute)
            })
            || moisture_bounds.iter().any(|m| !m.is_finite() || *m < 0.0)
            || moisture_bounds[0] > moisture_bounds[1]
        {
            return Err(invalid(
                "sorption energy needs its original card, fixed absolute temperature and ordered dry-basis moisture",
            ));
        }
        let bound = binding
            .property(property)
            .ok_or_else(|| invalid("sorption energy needs a uniquely bound heat property"))?;
        if bound.requirement().dims() != MOISTURE_ISOSTERIC_HEAT_DIMS {
            return Err(MaterialStatePointError::DimensionMismatch {
                property: property.into(),
                expected: MOISTURE_ISOSTERIC_HEAT_DIMS,
                found: bound.requirement().dims(),
            });
        }
        if !bound.has_only_owned_and_fixed_external_conditions(AXIS, dry) {
            return Err(invalid(
                "sorption energy requires a static differential heat on the owned dry-basis moisture coordinate",
            ));
        }
        let selected_claim = bound.answer().receipt.selected;
        let claim = card
            .claims()
            .claim(selected_claim)
            .ok_or_else(|| invalid("selected sorption heat source is absent"))?;
        let lower = point
            .clone()
            .with(AXIS, moisture_bounds[0])
            .map_err(|_| invalid("sorption heat lower moisture"))?;
        let upper = point
            .clone()
            .with(AXIS, moisture_bounds[1])
            .map_err(|_| invalid("sorption heat upper moisture"))?;
        let path = card
            .claims()
            .query_line_typed(
                &claim.key,
                AXIS,
                &lower,
                &upper,
                fs_matdb::ClaimSelection::Pinned(selected_claim),
            )
            .map_err(|source| MaterialStatePointError::Support {
                property: property.into(),
                source: Box::new(source),
            })?;
        let mut identity =
            DomainHasher::new("org.frankensim.fs-material.moisture-sorption-energy.v1");
        identity.update(binding.identity().as_bytes());
        identity.update(selected_claim.0.as_bytes());
        let mut knots = Vec::with_capacity(path.knots.len());
        for (m, answer) in path.knots {
            let value = answer.evidence.value.value;
            if !value.is_finite() {
                return Err(invalid("sorption heat source must remain finite"));
            }
            identity.update(&m.to_bits().to_le_bytes());
            identity.update(answer.receipt.content_hash().as_bytes());
            knots.push((m, value));
        }
        Ok(Self {
            knots: knots.into_boxed_slice(),
            source_point: point,
            card_identity: card.content_hash(),
            selected_claim,
            identity: identity.finalize(),
        })
    }

    pub fn moisture_bounds(&self) -> [f64; 2] {
        [self.knots[0].0, self.knots[self.knots.len() - 1].0]
    }
    pub fn knot_count(&self) -> usize {
        self.knots.len()
    }
    pub const fn source_point(&self) -> &QueryPoint {
        &self.source_point
    }
    pub const fn card_identity(&self) -> ContentHash {
        self.card_identity
    }
    pub const fn selected_claim(&self) -> ClaimId {
        self.selected_claim
    }
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    fn admit(&self, m: f64) -> Result<()> {
        let [lo, hi] = self.moisture_bounds();
        if !m.is_finite() || m < lo || m > hi {
            Err(invalid(
                "moisture lies outside the prepared sorption heat source",
            ))
        } else {
            Ok(())
        }
    }

    /// Differential enthalpy release [J/kg water] at the actual local loading.
    pub fn isosteric_heat(&self, moisture: f64) -> Result<f64> {
        self.admit(moisture)?;
        let i = self.knots.partition_point(|&(m, _)| m < moisture);
        if i == 0 || self.knots[i].0 == moisture {
            return Ok(self.knots[i].1);
        }
        let (a, x) = self.knots[i - 1];
        let (b, y) = self.knots[i];
        let f = (moisture - a) / (b - a);
        let value = (1.0 - f).mul_add(x, f * y);
        if value.is_finite() {
            Ok(value)
        } else {
            Err(invalid("unrepresentable sorption heat interpolation"))
        }
    }

    /// Exact number of source segments visited by this integral, for caller
    /// work admission. An unchanged state still validates support but costs zero.
    pub fn integration_segments(&self, a: f64, b: f64) -> Result<usize> {
        self.admit(a)?;
        self.admit(b)?;
        if a == b {
            return Ok(0);
        }
        let (lo, hi) = (a.min(b), a.max(b));
        Ok(1 + self.knots.partition_point(|&(m, _)| m < hi)
            - self.knots.partition_point(|&(m, _)| m <= lo))
    }

    /// Signed integral of q_st(M) dM [J/kg dry skeleton].
    ///
    /// Every crossed source segment contributes. Local integrals are evaluated
    /// directly, never by subtracting two large cumulative energy potentials.
    /// Reversing an interval reverses its sign; zero increments still require
    /// support. Finite signed heat permits endothermic source data without clipping.
    pub fn integral(&self, a: f64, b: f64) -> Result<f64> {
        self.admit(a)?;
        self.admit(b)?;
        if a == b {
            return Ok(0.0);
        }
        let (lo, hi) = (a.min(b), a.max(b));
        let width = hi - lo;
        let mut x = lo;
        let mut left = self.isosteric_heat(lo)?;
        let mut i = self.knots.partition_point(|&(m, _)| m <= lo);
        let mut mean = Dd::ZERO;
        loop {
            let end = self.knots[i].0.min(hi);
            let right = if end == self.knots[i].0 {
                self.knots[i].1
            } else {
                self.isosteric_heat(end)?
            };
            let midpoint = Dd::from_f64(0.5 * left) + Dd::from_f64(0.5 * right);
            mean = mean + midpoint * Dd::from_f64((end - x) / width);
            if end == hi {
                break;
            }
            x = end;
            left = right;
            i += 1;
        }
        let value = (mean * Dd::from_f64(b - a)).to_f64();
        if value.is_finite() {
            Ok(value)
        } else {
            Err(invalid("unrepresentable integrated sorption energy"))
        }
    }
}
