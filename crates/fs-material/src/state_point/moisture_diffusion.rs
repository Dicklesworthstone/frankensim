//! Exact nonlinear moisture diffusion on a continuously supported source line.
//!
//! For one homogeneous law, `D(M) grad M = grad K(M)` with
//! `K(M) = integral D(m) dm`. The secant of K stays positive across a wet/dry
//! interface even when D vanishes at the dry endpoint. A harmonic average of
//! endpoint D values instead treats that endpoint as an impermeable layer.
//! These are different models: this chart describes a state-dependent law,
//! not a spatial mixture or an interface between distinct materials.

use super::*;
use fs_qty::semantic::MoistureBasis;

const DIFFUSIVITY_DIMS: Dims = Dims([2, 0, -1, 0, 0, 0]);
const MOISTURE_AXIS: &str = "moisture_content";
type Result<T> = core::result::Result<T, MaterialStatePointError>;

fn invalid(quantity: &'static str) -> MaterialStatePointError {
    MaterialStatePointError::InvalidDerived { quantity }
}

/// Prepared nonnegative D(M) [m²/s] at fixed external conditions.
///
/// Preparation checks the complete original source path once, including all
/// intervening table segments and surface cells. Evaluation allocates nothing:
/// a point lookup is logarithmic in the knot count, and an integral visits only
/// the crossed segments. The existing source owner bounds the knot count.
/// Integration is exact for its piecewise-affine interpolation, up to floating
/// point roundoff; this is not an uncertainty or continuum-error certificate.
#[derive(Clone, Debug, PartialEq)]
pub struct KirchhoffMoistureDiffusivity {
    knots: Box<[(f64, f64)]>,
    selected_claim: ClaimId,
    identity: ContentHash,
}

impl KirchhoffMoistureDiffusivity {
    /// Prepare the selected property over an ordered, closed moisture interval.
    ///
    /// `binding` must belong to this immutable card and already select the
    /// diffusivity property. All other coordinates remain exactly as bound.
    /// Moisture is explicitly kg water/kg dry solid; temperature is explicitly
    /// absolute. Zero D is admitted, including at a dry front. Negative D at
    /// ANY intervening knot, an unsupported source gap, and extrapolation all
    /// refuse. Scalar claims and fixed-context cuts through curves or surfaces
    /// use the same canonical `query_line_typed` owner.
    pub fn from_binding(
        card: &MaterialCard,
        binding: &ResolvedMaterialStatePoint,
        property: &str,
        moisture_bounds: [f64; 2],
    ) -> Result<Self> {
        let dry = QuantitySpec::semantic(SemanticType::new(
            QuantityKind::Moisture(MoistureBasis::DryMass),
            ValueForm::Static,
        ));
        let absolute = QuantitySpec::semantic(SemanticType::new(
            QuantityKind::AbsoluteTemperature,
            ValueForm::Static,
        ));
        let point = binding.query();
        if card.content_hash() != binding.card_identity()
            || point.axis_quantities().get(MOISTURE_AXIS) != Some(&dry)
            || point.axis_quantities().get("temperature") != Some(&absolute)
            || !point.axes().get("temperature").is_some_and(|t| *t > 0.0)
            || !point.axes().get(MOISTURE_AXIS).is_some_and(|m| *m >= 0.0)
            || point.axes().contains_key("T")
            || point.axis_quantities().iter().any(|(axis, quantity)| {
                (axis != MOISTURE_AXIS && *quantity == dry)
                    || (axis != "temperature" && *quantity == absolute)
            })
            || moisture_bounds.iter().any(|m| !m.is_finite() || *m < 0.0)
            || moisture_bounds[0] > moisture_bounds[1]
        {
            return Err(invalid(
                "diffusion chart needs its original card, fixed absolute temperature and ordered dry-basis moisture",
            ));
        }
        let bound = binding
            .property(property)
            .ok_or_else(|| invalid("diffusion chart needs a uniquely bound property"))?;
        if bound.requirement().dims() != DIFFUSIVITY_DIMS {
            return Err(MaterialStatePointError::DimensionMismatch {
                property: property.into(),
                expected: DIFFUSIVITY_DIMS,
                found: bound.requirement().dims(),
            });
        }
        if !bound.is_static_at_fixed_external_conditions()
            || bound
                .source_axes()
                .iter()
                .any(|(axis, dims)| axis == MOISTURE_AXIS && *dims != Dims::NONE)
        {
            return Err(invalid(
                "diffusion chart requires a static environmental law on dimensionless dry-basis moisture",
            ));
        }
        let selected_claim = bound.answer().receipt.selected;
        let claim = card
            .claims()
            .claim(selected_claim)
            .ok_or_else(|| invalid("selected diffusion source is absent"))?;
        let lower = point
            .clone()
            .with(MOISTURE_AXIS, moisture_bounds[0])
            .map_err(|_| invalid("diffusion chart lower moisture"))?;
        let upper = point
            .with(MOISTURE_AXIS, moisture_bounds[1])
            .map_err(|_| invalid("diffusion chart upper moisture"))?;
        let path = card
            .claims()
            .query_line_typed(
                &claim.key,
                MOISTURE_AXIS,
                &lower,
                &upper,
                fs_matdb::ClaimSelection::Pinned(selected_claim),
            )
            .map_err(|source| MaterialStatePointError::Support {
                property: property.into(),
                source: Box::new(source),
            })?;
        let mut identity =
            DomainHasher::new("org.frankensim.fs-material.kirchhoff-moisture-diffusivity.v1");
        identity.update(binding.identity().as_bytes());
        identity.update(selected_claim.0.as_bytes());
        let mut knots = Vec::with_capacity(path.knots.len());
        for (moisture, answer) in path.knots {
            let value = answer.evidence.value.value;
            if !value.is_finite() || value < 0.0 {
                return Err(MaterialStatePointError::OutsideConsumerDomain {
                    property: property.into(),
                    value,
                    admissibility: ScalarAdmissibility::NonNegative,
                });
            }
            identity.update(&moisture.to_bits().to_le_bytes());
            identity.update(answer.receipt.content_hash().as_bytes());
            knots.push((moisture, value));
        }
        Ok(Self {
            knots: knots.into_boxed_slice(),
            selected_claim,
            identity: identity.finalize(),
        })
    }

    /// Closed supported interval in kg water/kg dry solid.
    pub fn moisture_bounds(&self) -> [f64; 2] {
        [self.knots[0].0, self.knots[self.knots.len() - 1].0]
    }

    /// Work bound: an interval evaluation visits at most this many knots.
    /// The canonical source-line owner admits at most 4096.
    pub fn knot_count(&self) -> usize {
        self.knots.len()
    }

    /// Exact immutable claim used throughout this chart.
    pub const fn selected_claim(&self) -> ClaimId {
        self.selected_claim
    }

    /// Binds card, selected claim, held environment and complete admitted path.
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    fn admit(&self, moisture: f64) -> Result<()> {
        let [lower, upper] = self.moisture_bounds();
        if !moisture.is_finite() || moisture < lower || moisture > upper {
            Err(invalid(
                "moisture lies outside the prepared diffusion source path",
            ))
        } else {
            Ok(())
        }
    }

    /// Original piecewise-affine diffusivity at an admitted moisture point.
    pub fn diffusivity(&self, moisture: f64) -> Result<f64> {
        self.admit(moisture)?;
        let i = self.knots.partition_point(|&(x, _)| x < moisture);
        if i == 0 || self.knots[i].0 == moisture {
            return Ok(self.knots[i].1);
        }
        let (left, a) = self.knots[i - 1];
        let (right, b) = self.knots[i];
        let fraction = (moisture - left) / (right - left);
        let value = (1.0 - fraction).mul_add(a, fraction * b);
        if value.is_finite() && value >= 0.0 {
            Ok(value)
        } else {
            Err(invalid("unrepresentable source diffusivity interpolation"))
        }
    }

    /// Integral mean D across an interval, with the continuous equal-state limit.
    ///
    /// Contributions are averaged directly over the requested interval, never
    /// obtained by subtracting two large cumulative potentials. This preserves
    /// a small local driving gradient even far from the chart's lower bound.
    /// The result is symmetric in its endpoints and nonnegative without clipping.
    pub fn mean_diffusivity(&self, a: f64, b: f64) -> Result<f64> {
        self.admit(a)?;
        self.admit(b)?;
        if a == b {
            return self.diffusivity(a);
        }
        let (lower, upper) = if a < b { (a, b) } else { (b, a) };
        let width = upper - lower;
        let mut x = lower;
        let mut left_value = self.diffusivity(lower)?;
        let mut i = self.knots.partition_point(|&(x, _)| x <= lower);
        let mut result = 0.0;
        let mut correction = 0.0;
        loop {
            let end = self.knots[i].0.min(upper);
            let right_value = if end == self.knots[i].0 {
                self.knots[i].1
            } else {
                self.diffusivity(end)?
            };
            let low = left_value.min(right_value);
            let mean = low + 0.5 * (left_value.max(right_value) - low);
            let term = mean * ((end - x) / width);
            // Compensated nonnegative accumulation; no cancelling K(b)-K(a).
            let add = term - correction;
            let next = result + add;
            correction = (next - result) - add;
            result = next;
            if end == upper {
                break;
            }
            x = end;
            left_value = right_value;
            i += 1;
        }
        if result.is_finite() && result >= 0.0 {
            Ok(result)
        } else {
            Err(invalid("unrepresentable mean source diffusivity"))
        }
    }

    /// Signed `K(b)-K(a) = integral_a^b D(M) dM`, with no arbitrary datum.
    pub fn integral(&self, a: f64, b: f64) -> Result<f64> {
        let value = self.mean_diffusivity(a, b)? * (b - a);
        if value.is_finite() {
            Ok(value)
        } else {
            Err(invalid("unrepresentable integrated source diffusivity"))
        }
    }
}
