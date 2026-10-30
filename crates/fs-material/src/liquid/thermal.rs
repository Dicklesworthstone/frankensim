//! A liquid's thermal expansion from its original density interpolation.
//!
//! At fixed mass, V=m/rho, so (1/V)(dV/dT) = -(1/rho)(d rho/dT).
//! Pressure, composition and every other source coordinate remain fixed.
//! Differentiating an admitted source line adds neither measurements nor an
//! equation of state outside that line. Source uncertainty is not a bound on
//! its derivative, and a table kink is not a differentiable physical law.

use std::sync::Arc;

use fs_matdb::{
    ClaimSelection, EvaluationDecision, MaterialCard, PropertyClaim, PropertyLineSlice,
    PropertyValue, QueryPoint,
};
use fs_qty::Temperature;

use super::ResolvedLiquidState;
use crate::liquid_temperature::validate_temperature;
use crate::state_point::{MaterialStatePointError, ResolvedMaterialStatePoint};

#[cfg(test)]
use crate::state_point::DENSITY_PROPERTY;

#[derive(Clone, Debug, PartialEq)]
struct DensityExpansionPath {
    context: QueryPoint,
    density_property: String,
    source: PropertyLineSlice,
    slopes_kg_m3_k: Vec<f64>,
}

/// One local density tangent and its complete, source-pinned temperature path.
/// The path is shared across updated liquid states, not copied per evaluation.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedLiquidExpansion {
    path: Arc<DensityExpansionPath>,
    temperature_k: f64,
    density_slope_kg_m3_k: f64,
    volumetric_expansion_per_k: f64,
}

impl ResolvedLiquidExpansion {
    /// Absolute temperature of the SAME density and transport state [K].
    #[must_use]
    pub const fn temperature_k(&self) -> f64 {
        self.temperature_k
    }

    /// Signed derivative of the selected density interpolation [kg/(m3 K)].
    #[must_use]
    pub const fn density_slope_kg_m3_k(&self) -> f64 {
        self.density_slope_kg_m3_k
    }

    /// Signed isobaric expansion, `-(d rho/dT)/rho` [1/K].
    #[must_use]
    pub const fn volumetric_expansion_per_k(&self) -> f64 {
        self.volumetric_expansion_per_k
    }

    /// Source-line breakpoints and the original usage receipt at every point.
    #[must_use]
    pub fn source(&self) -> &PropertyLineSlice {
        &self.path.source
    }

    pub(super) fn at_material(
        &self,
        material: &ResolvedMaterialStatePoint,
        density: f64,
    ) -> Result<Self, MaterialStatePointError> {
        sample(Arc::clone(&self.path), material, density)
    }
}

impl ResolvedLiquidState {
    /// Bind the current density source's continuous thermal interval, retaining
    /// its original card, claim, all nonthermal coordinates and source receipts.
    ///
    /// A density scalar needs an explicit temperature-validity interval and
    /// constant interpolation; an isolated density measurement does not imply
    /// zero expansion. Curves and fixed-context surface slices use the native
    /// piecewise-affine interpolation. Gaps and discrete-only spans refuse.
    /// Zero and negative expansion are preserved. Interior slope discontinuities
    /// refuse rather than selecting or averaging a side; interval endpoints
    /// use one-sided interior derivatives. This is a local coefficient, not a
    /// finite-contrast Boussinesq or natural-convection accuracy claim.
    ///
    /// Later `at_query_point` calls update expansion with transport. Changing a
    /// nonthermal coordinate requires an explicit new path and cannot reuse the
    /// old pressure/composition derivative. A retained independent expansion
    /// claim refuses this binding: the caller must choose one source of beta.
    pub fn with_density_expansion(
        mut self,
        card: &MaterialCard,
        temperature_axis: &str,
        lower: Temperature,
        upper: Temperature,
    ) -> Result<Self, MaterialStatePointError> {
        if self.material.card_identity() != card.content_hash() {
            return Err(invalid(
                "liquid density expansion requires the original card",
            ));
        }
        if self.volumetric_expansion_per_k.is_some() {
            return Err(invalid(
                "liquid density expansion cannot replace an independently selected expansion claim",
            ));
        }
        let context = self.material.query();
        validate_temperature(&context, temperature_axis)?;
        let lo = lower.value();
        let hi = upper.value();
        if !lo.is_finite() || !hi.is_finite() || lo <= 0.0 || hi <= lo {
            return Err(invalid(
                "liquid density expansion needs an increasing positive temperature interval",
            ));
        }
        let density_property = self.requirements.properties()[0].name().to_owned();
        let property = self
            .material
            .property(&density_property)
            .expect("liquid resolution retains its density requirement");
        let selected = property.answer().receipt.selected;
        let claim = card
            .claims()
            .claim(selected)
            .ok_or_else(|| invalid("selected liquid density claim is absent"))?;
        if claim
            .value
            .axes()
            .any(|(axis, dims)| axis == temperature_axis && dims != Temperature::DIMS)
        {
            return Err(invalid(
                "liquid density source temperature has incorrect dimensions",
            ));
        }
        let thermal_source = claim.value.axes().any(|(axis, _)| axis == temperature_axis)
            || claim
                .validity
                .bound(temperature_axis)
                .is_some_and(|(a, b)| b > a);
        if !thermal_source {
            return Err(invalid(
                "a density value without temperature support cannot supply expansion",
            ));
        }
        let from = context
            .clone()
            .with(temperature_axis, lo)
            .map_err(|_| invalid("invalid liquid density path start"))?;
        let to = context
            .clone()
            .with(temperature_axis, hi)
            .map_err(|_| invalid("invalid liquid density path end"))?;
        let source = card
            .claims()
            .query_line_typed(
                &claim.key,
                temperature_axis,
                &from,
                &to,
                ClaimSelection::Pinned(selected),
            )
            .map_err(|source| MaterialStatePointError::Support {
                property: density_property.clone(),
                source: Box::new(source),
            })?;
        if source.knots.len() < 2
            || source.knots.iter().any(|(_, answer)| {
                !answer.evidence.value.value.is_finite() || answer.evidence.value.value <= 0.0
            })
        {
            return Err(invalid(
                "liquid density expansion needs a positive continuous density path",
            ));
        }
        let slopes_kg_m3_k = native_slopes(card, claim, &context, &source)?;
        self.density_expansion = Some(sample(
            Arc::new(DensityExpansionPath {
                context,
                density_property,
                source,
                slopes_kg_m3_k,
            }),
            &self.material,
            self.density_kg_m3,
        )?);
        Ok(self)
    }

    /// Source-bound expansion, when a continuous density path was admitted.
    /// Absence remains unknown; it does not mean zero expansion.
    #[must_use]
    pub const fn density_expansion(&self) -> Option<&ResolvedLiquidExpansion> {
        self.density_expansion.as_ref()
    }
}

fn sample(
    path: Arc<DensityExpansionPath>,
    material: &ResolvedMaterialStatePoint,
    density: f64,
) -> Result<ResolvedLiquidExpansion, MaterialStatePointError> {
    let point = material.query();
    let axis = &path.source.axis;
    if point.axis_quantities() != path.context.axis_quantities()
        || point.axes().keys().ne(path.context.axes().keys())
        || point.axes().iter().any(|(name, value)| {
            name != axis && value.to_bits() != path.context.axes()[name].to_bits()
        })
    {
        return Err(invalid(
            "liquid density expansion cannot change its fixed pressure/composition context",
        ));
    }
    let temperature = point.axes()[axis];
    let knots = &path.source.knots;
    let last = knots.len() - 1;
    if temperature < knots[0].0 || temperature > knots[last].0 {
        return Err(invalid(
            "liquid temperature leaves its admitted density-expansion interval",
        ));
    }
    let density_property = material
        .property(&path.density_property)
        .expect("liquid resolution retains its density requirement");
    if density_property.answer().receipt.selected != knots[0].1.receipt.selected {
        return Err(invalid(
            "liquid density expansion cannot switch source claims",
        ));
    }
    let upper = knots.partition_point(|(t, _)| *t <= temperature);
    let segment = (upper - 1).min(last - 1);
    let density_slope = path.slopes_kg_m3_k[segment];
    if segment > 0 && temperature == knots[segment].0 {
        let left = path.slopes_kg_m3_k[segment - 1];
        let scale = left.abs().max(density_slope.abs());
        if (left - density_slope).abs() > 128.0 * f64::EPSILON * scale {
            return Err(invalid(
                "liquid density tangent is discontinuous at this source knot",
            ));
        }
    }
    let expansion = -density_slope / density;
    if !expansion.is_finite() || (density_slope != 0.0 && expansion == 0.0) {
        return Err(invalid("liquid volumetric expansion is unrepresentable"));
    }
    Ok(ResolvedLiquidExpansion {
        path,
        temperature_k: temperature,
        density_slope_kg_m3_k: density_slope,
        volumetric_expansion_per_k: expansion,
    })
}

// Differentiate the native source coefficients, not densities evaluated at
// caller-clipped endpoints. Those can round to the same density on a narrow
// interval even when the original source slope is readily representable.
fn native_slopes(
    card: &MaterialCard,
    claim: &PropertyClaim,
    context: &QueryPoint,
    source: &PropertyLineSlice,
) -> Result<Vec<f64>, MaterialStatePointError> {
    source
        .knots
        .windows(2)
        .map(|span| {
            let (numerator, denominator) = match &claim.value {
                PropertyValue::Curve {
                    abscissa, knots, ..
                } if abscissa == &source.axis => {
                    let upper = knots.partition_point(|(x, _)| *x <= span[0].0);
                    let lower = upper - 1;
                    if lower + 1 >= knots.len() || span[1].0 > knots[lower + 1].0 {
                        return Err(invalid(
                            "liquid density interval crosses a native source segment",
                        ));
                    }
                    (
                        knots[lower + 1].1 - knots[lower].1,
                        knots[lower + 1].0 - knots[lower].0,
                    )
                }
                PropertyValue::Surface2D {
                    axes,
                    samples,
                    triangles,
                    ..
                } if axes.contains(&source.axis) => {
                    let interior = span[0].0 + 0.5 * (span[1].0 - span[0].0);
                    if interior <= span[0].0 || interior >= span[1].0 {
                        return Err(invalid(
                            "liquid density cell needs a representable interior coordinate",
                        ));
                    }
                    let query = context
                        .clone()
                        .with(&source.axis, interior)
                        .map_err(|_| invalid("invalid liquid density cell coordinate"))?;
                    let answer = card
                        .claims()
                        .query_pinned_typed(&claim.key, &query, span[0].1.receipt.selected)
                        .map_err(|source| MaterialStatePointError::Query {
                            property: claim.key.name().to_owned(),
                            source,
                        })?;
                    let EvaluationDecision::LinearInsideTriangle { triangle } =
                        answer.receipt.decision
                    else {
                        return Err(invalid(
                            "liquid density tangent needs an original affine source cell",
                        ));
                    };
                    let [a, b, c] = triangles[triangle as usize].map(|i| samples[i as usize]);
                    // Match the source surface owner's vertex pivot. The
                    // derivative uses the same affine triangle and fixed
                    // remaining coordinate, with no perturbed point queries.
                    let ax = a.0[0] - c.0[0];
                    let ay = a.0[1] - c.0[1];
                    let bx = b.0[0] - c.0[0];
                    let by = b.0[1] - c.0[1];
                    let av = a.1 - c.1;
                    let bv = b.1 - c.1;
                    let numerator = if axes[0] == source.axis {
                        av * by - bv * ay
                    } else {
                        ax * bv - bx * av
                    };
                    (numerator, ax * by - ay * bx)
                }
                // The admitted response is independent of temperature. The
                // path gate still requires explicit thermal validity support.
                _ => (0.0, 1.0),
            };
            let slope = numerator / denominator;
            if !numerator.is_finite()
                || !denominator.is_finite()
                || denominator == 0.0
                || !slope.is_finite()
                || (numerator != 0.0 && slope == 0.0)
            {
                Err(invalid("liquid density tangent is unrepresentable"))
            } else {
                Ok(slope)
            }
        })
        .collect()
}

fn invalid(quantity: &'static str) -> MaterialStatePointError {
    MaterialStatePointError::InvalidDerived { quantity }
}

#[cfg(test)]
mod tests;
