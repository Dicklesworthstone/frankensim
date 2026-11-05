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
/// both point and interval queries take logarithmic work in the knot count.
/// A balanced tree of nonnegative interval means takes linear preparation work
/// and storage. The existing source owner bounds the knot count.
/// Integration is exact for its piecewise-affine interpolation, up to floating
/// point roundoff; this is not an uncertainty or continuum-error certificate.
#[derive(Clone, Debug, PartialEq)]
pub struct KirchhoffMoistureDiffusivity {
    knots: Box<[(f64, f64)]>,
    interval_means: Box<[f64]>,
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
        let segments = knots.len() - 1;
        let mut interval_means = if segments == 0 {
            Vec::new()
        } else {
            vec![0.0; 2 * segments.next_power_of_two()]
        };
        if segments != 0 {
            prepare_means(&knots, &mut interval_means, 1, 0, segments);
        }
        Ok(Self {
            knots: knots.into_boxed_slice(),
            interval_means: interval_means.into_boxed_slice(),
            selected_claim,
            identity: identity.finalize(),
        })
    }

    /// Closed supported interval in kg water/kg dry solid.
    pub fn moisture_bounds(&self) -> [f64; 2] {
        [self.knots[0].0, self.knots[self.knots.len() - 1].0]
    }

    /// Number of admitted source knots; the canonical owner admits at most 4096.
    /// Interval queries reuse prepared means rather than visiting every knot.
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
        let value = interpolate(self.knots[i - 1], self.knots[i], moisture);
        if value.is_finite() && value >= 0.0 {
            Ok(value)
        } else {
            Err(invalid("unrepresentable source diffusivity interpolation"))
        }
    }

    /// Integral mean D across an interval, with the continuous equal-state limit.
    ///
    /// Fully covered source intervals reuse their prepared nonnegative means;
    /// only the two partial end segments require interpolation. Means combine
    /// by physical interval length, never by subtracting cumulative potentials
    /// or forming potentially overflowing absolute integrals. This preserves a
    /// small local driving gradient even far from the chart's lower bound.
    /// The result is symmetric in its endpoints and nonnegative without clipping.
    pub fn mean_diffusivity(&self, a: f64, b: f64) -> Result<f64> {
        self.admit(a)?;
        self.admit(b)?;
        if a == b {
            return self.diffusivity(a);
        }
        let (lower, upper) = if a < b { (a, b) } else { (b, a) };
        let result = self.interval_mean(lower, upper, 1, 0, self.knots.len() - 1);
        if result.is_finite() && result >= 0.0 {
            Ok(result)
        } else {
            Err(invalid("unrepresentable mean source diffusivity"))
        }
    }

    fn interval_mean(&self, lower: f64, upper: f64, node: usize, l: usize, r: usize) -> f64 {
        if lower == self.knots[l].0 && upper == self.knots[r].0 {
            return self.interval_means[node];
        }
        if r - l == 1 {
            return midpoint(
                interpolate(self.knots[l], self.knots[r], lower),
                interpolate(self.knots[l], self.knots[r], upper),
            );
        }
        let m = l + (r - l) / 2;
        let split = self.knots[m].0;
        if upper <= split {
            self.interval_mean(lower, upper, 2 * node, l, m)
        } else if lower >= split {
            self.interval_mean(lower, upper, 2 * node + 1, m, r)
        } else {
            join_means(
                self.interval_mean(lower, split, 2 * node, l, m),
                self.interval_mean(split, upper, 2 * node + 1, m, r),
                [lower, split, upper],
            )
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

fn interpolate((left, a): (f64, f64), (right, b): (f64, f64), x: f64) -> f64 {
    if x == left {
        return a;
    }
    if x == right {
        return b;
    }
    let fraction = (x - left) / (right - left);
    (1.0 - fraction).mul_add(a, fraction * b)
}

fn midpoint(a: f64, b: f64) -> f64 {
    a.min(b) + 0.5 * (a.max(b) - a.min(b))
}

fn prepare_means(knots: &[(f64, f64)], means: &mut [f64], node: usize, l: usize, r: usize) -> f64 {
    let mean = if r - l == 1 {
        midpoint(knots[l].1, knots[r].1)
    } else {
        let m = l + (r - l) / 2;
        join_means(
            prepare_means(knots, means, 2 * node, l, m),
            prepare_means(knots, means, 2 * node + 1, m, r),
            [knots[l].0, knots[m].0, knots[r].0],
        )
    };
    means[node] = mean;
    mean
}

// Convex averaging uses the smaller side's weight and never sums two large
// positive products. Equal laws retain their exact constant limit. Coordinates
// supply the complete width directly, without accumulating rounded lengths.
fn join_means(a: f64, b: f64, [left, split, right]: [f64; 3]) -> f64 {
    if a == b {
        return a;
    }
    let (low, high, low_width, high_width) = if a < b {
        (a, b, split - left, right - split)
    } else {
        (b, a, right - split, split - left)
    };
    if high_width <= low_width {
        low + weighted_difference(high - low, high_width, right - left)
    } else {
        high - weighted_difference(high - low, low_width, right - left)
    }
}

// A tiny geometric weight can underflow even when its product with a large D
// remains significant. Keep the ordinary path cheap; only subnormal weights
// need normalized mantissas and one final gradual-underflow rounding.
fn weighted_difference(value: f64, width: f64, total: f64) -> f64 {
    let weight = width / total;
    if weight >= f64::MIN_POSITIVE || value == 0.0 {
        return value * weight;
    }
    let parts = |x: f64| {
        let (x, shift) = if x < f64::MIN_POSITIVE {
            (x * 18_014_398_509_481_984.0, -54)
        } else {
            (x, 0)
        };
        let bits = x.to_bits();
        (
            f64::from_bits((bits & 0x000f_ffff_ffff_ffff) | 0x3ff0_0000_0000_0000),
            ((bits >> 52) as i32) - 1023 + shift,
        )
    };
    let (v, ve) = parts(value);
    let (w, we) = parts(width);
    let (t, te) = parts(total);
    let mut mantissa = v * w / t;
    let mut exponent = ve + we - te;
    if mantissa >= 2.0 {
        mantissa *= 0.5;
        exponent += 1;
    } else if mantissa < 1.0 {
        mantissa *= 2.0;
        exponent -= 1;
    }
    if exponent >= -1022 {
        mantissa * f64::from_bits(((exponent + 1023) as u64) << 52)
    } else if exponent < -1075 {
        0.0
    } else {
        (mantissa * f64::MIN_POSITIVE) * f64::from_bits(((exponent + 2045) as u64) << 52)
    }
}
