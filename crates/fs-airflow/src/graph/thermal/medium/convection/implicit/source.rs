//! Bounded original scalar/temperature-curve coefficients, prepared once.

use super::{LiquidExteriorConvectionError as Error, finite};
use fs_matdb::{ClaimSelection, InterpolationPolicy, PropertyLineSlice, PropertyValue};
use fs_material::fluid_transport::FluidTransportState;
use fs_material::liquid_temperature::LiquidTemperatureLaw;
use fs_material::state_point::MaterialStatePointError;
use fs_qty::Temperature;

pub(super) struct TransportSample {
    pub state: FluidTransportState,
    // Derivatives in rho, Cp, k, dynamic-mu order, before the film chain factor.
    pub derivatives: [f64; 4],
}

pub(super) struct PreparedScalar {
    pub source: PropertyLineSlice,
    segments: Vec<Segment>,
    positive: bool,
}

struct Segment {
    // Original native coefficients, not a slope between clipped point queries.
    endpoints: [(f64, f64); 2],
    derivative: f64,
}

impl PreparedScalar {
    pub fn new(
        liquid: &LiquidTemperatureLaw<'_>,
        name: &str,
        bounds: [f64; 2],
        positive: bool,
    ) -> Result<Self, Error> {
        let card = liquid.card();
        let original = liquid
            .reference_state()
            .material()
            .property(name)
            .ok_or(Error::Invalid(
                "the original liquid bundle lacks its required source",
            ))?;
        let selected = original.answer().receipt.selected;
        let claim = card.claims().claim(selected).ok_or(Error::Invalid(
            "the original selected liquid claim is absent",
        ))?;
        let axis = liquid.temperature_axis();
        match &claim.value {
            PropertyValue::Scalar { .. }
                if claim.interpolation == InterpolationPolicy::ConstantWithinValidity
                    && claim.validity.bound(axis).is_some_and(|(lo, hi)| hi > lo) => {}
            PropertyValue::Curve {
                abscissa,
                abscissa_dims,
                ..
            } if abscissa == axis
                && *abscissa_dims == Temperature::DIMS
                && claim.interpolation == InterpolationPolicy::LinearInside => {}
            _ => {
                return Err(Error::Invalid(
                    "implicit liquid sources require a constant with thermal validity or an original linear temperature curve",
                ));
            }
        }
        let point = |temperature| {
            liquid
                .reference_point()
                .clone()
                .with(axis, temperature)
                .map_err(|source| MaterialStatePointError::Query {
                    property: name.to_owned(),
                    source,
                })
        };
        let source = card
            .claims()
            .query_line_typed(
                &claim.key,
                axis,
                &point(bounds[0])?,
                &point(bounds[1])?,
                ClaimSelection::Pinned(selected),
            )
            .map_err(|source| MaterialStatePointError::Support {
                property: name.to_owned(),
                source: Box::new(source),
            })?;
        if source.knots.len() < 2
            || source.knots.iter().any(|(_, answer)| {
                let value = answer.evidence.value.value;
                !value.is_finite() || (positive && value <= 0.0)
            })
        {
            return Err(Error::Invalid(
                "liquid source path must have finite physical values and nonzero span",
            ));
        }
        let mut segments = Vec::with_capacity(source.knots.len() - 1);
        for span in source.knots.windows(2) {
            let endpoints = match &claim.value {
                PropertyValue::Scalar { value, .. } => [(span[0].0, *value), (span[1].0, *value)],
                PropertyValue::Curve { knots, .. } => {
                    let index = knots.partition_point(|(t, _)| *t <= span[0].0) - 1;
                    if index + 1 >= knots.len() || span[1].0 > knots[index + 1].0 {
                        return Err(Error::Invalid(
                            "prepared source span crosses an original curve segment",
                        ));
                    }
                    [knots[index], knots[index + 1]]
                }
                _ => unreachable!("source form checked before preparation"),
            };
            let width = endpoints[1].0 - endpoints[0].0;
            let difference = endpoints[1].1 - endpoints[0].1;
            let derivative = finite(difference / width)?;
            if !width.is_finite()
                || width <= 0.0
                || !difference.is_finite()
                || (difference != 0.0 && derivative == 0.0)
            {
                return Err(Error::Invalid(
                    "original liquid source tangent is unrepresentable",
                ));
            }
            segments.push(Segment {
                endpoints,
                derivative,
            });
        }
        Ok(Self {
            source,
            segments,
            positive,
        })
    }

    pub fn sample(&self, temperature: f64, needs_derivative: bool) -> Result<(f64, f64), Error> {
        let knots = &self.source.knots;
        let last = knots.len() - 1;
        if !temperature.is_finite() || temperature < knots[0].0 || temperature > knots[last].0 {
            return Err(Error::Invalid(
                "fluid temperature leaves the original prepared source interval",
            ));
        }
        let upper = knots.partition_point(|(t, _)| *t <= temperature);
        let index = (upper - 1).min(last - 1);
        let segment = &self.segments[index];
        if needs_derivative && index > 0 && temperature == knots[index].0 {
            let left = self.segments[index - 1].derivative;
            let right = segment.derivative;
            // Only rounding-equivalent native slopes can represent one tangent.
            if (left - right).abs() > 128.0 * f64::EPSILON * left.abs().max(right.abs()) {
                return Err(Error::Invalid(
                    "liquid source has no unique derivative at this native knot",
                ));
            }
        }
        let value = if temperature == knots[index].0 {
            knots[index].1.evidence.value.value
        } else if temperature == knots[index + 1].0 {
            knots[index + 1].1.evidence.value.value
        } else {
            let [(a, ya), (b, yb)] = segment.endpoints;
            let fraction = (temperature - a) / (b - a);
            finite((1.0 - fraction) * ya + fraction * yb)?
        };
        if !value.is_finite() || (self.positive && value <= 0.0) {
            return Err(Error::Invalid(
                "original source interpolation produced a nonphysical value",
            ));
        }
        Ok((
            value,
            if needs_derivative {
                segment.derivative
            } else {
                0.0
            },
        ))
    }
}
