//! Pinned scalar source response along a simultaneous affine state change.

use super::{
    ClaimId, ClaimSelection, ClaimSet, InterpolationPolicy, MAX_PROPERTY_LINE_KNOTS, MatDbError,
    MaterialAnswer, PropertyKey, PropertySupportError, PropertyValue, QueryPoint,
};

/// One continuously supported original scalar law along an affine state path.
///
/// The path is `start + t * (end - start)`, with dimensionless `t` in `[0,1]`.
/// All coordinates retain their original quantity descriptors. Every original
/// crossed source knot/cell edge is retained, with a canonical point answer.
/// This is numerical interpolation support, not a new property claim, mixing
/// law, derivative enclosure, or assertion of support in the surrounding box.
#[derive(Debug, Clone, PartialEq)]
pub struct PropertyPathSlice {
    /// Complete original physical query at affine parameter zero.
    pub start: QueryPoint,
    /// Complete original physical query at affine parameter one.
    pub end: QueryPoint,
    /// Ordered path parameters and original source answers, including endpoints.
    pub knots: Vec<(f64, MaterialAnswer)>,
}

impl ClaimSet {
    /// Follow one explicitly pinned claim along a simultaneous state change.
    ///
    /// The caller owns source selection; this method never switches claims.
    /// Source validity is rectangular, hence endpoint validity covers the
    /// segment's applicability coordinates. Sampled response support is a
    /// separate obligation: original cell connectivity must cover the whole
    /// path, including every interior point. Axis-aligned paths use the existing
    /// source-line owner. Work is bounded by the original source and line caps.
    pub fn query_path_pinned_typed(
        &self,
        property: &PropertyKey,
        start: &QueryPoint,
        end: &QueryPoint,
        selected: ClaimId,
    ) -> Result<PropertyPathSlice, PropertySupportError> {
        let invalid = |reason: &str| PropertySupportError::InvalidEnvelope {
            reason: reason.to_owned(),
        };
        if start.axes().keys().ne(end.axes().keys())
            || start.axis_quantities() != end.axis_quantities()
            || start
                .axes()
                .iter()
                .any(|(axis, &a)| !(end.axes()[axis] - a).is_finite())
        {
            return Err(invalid(
                "an affine path requires identical coordinate names/quantities and finite spans",
            ));
        }
        let query = |point: &QueryPoint| {
            self.query_pinned_typed(property, point, selected)
                .and_then(|answer| {
                    self.verify_receipt(&answer.receipt)?;
                    Ok(answer)
                })
                .map_err(|error| PropertySupportError::Evaluation {
                    point: point.clone(),
                    error: Box::new(error),
                })
        };
        let first = query(start)?;
        let last = query(end)?;
        let varying: Vec<_> = start
            .axes()
            .iter()
            .filter_map(|(axis, &value)| (value != end.axes()[axis]).then_some(axis.as_str()))
            .collect();
        if varying.len() == 1 {
            let axis = varying[0];
            let a = start.axes()[axis];
            let b = end.axes()[axis];
            let (lower, upper) = if a < b { (start, end) } else { (end, start) };
            let line = self.query_line_typed(
                property,
                axis,
                lower,
                upper,
                ClaimSelection::Pinned(selected),
            )?;
            let mut knots: Vec<_> = line
                .knots
                .into_iter()
                .map(|(x, answer)| {
                    let t = if x == a {
                        0.0
                    } else if x == b {
                        1.0
                    } else {
                        (x - a) / (b - a)
                    };
                    (t, answer)
                })
                .collect();
            if a > b {
                knots.reverse();
            }
            if knots.windows(2).any(|pair| pair[0].0 >= pair[1].0) {
                return Err(invalid(
                    "source line knots collapse in affine path coordinates",
                ));
            }
            return Ok(PropertyPathSlice {
                start: start.clone(),
                end: end.clone(),
                knots,
            });
        }
        let claim = self
            .claim(selected)
            .expect("the original query admitted this exact claim");
        let mut parameters = vec![(0.0, None::<(&str, f64)>), (1.0, None)];
        match &claim.value {
            PropertyValue::Scalar { .. } => {}
            PropertyValue::Curve {
                abscissa, knots, ..
            } => {
                let a = start.axes()[abscissa];
                let b = end.axes()[abscissa];
                if a != b {
                    if claim.interpolation == InterpolationPolicy::TabulatedOnly {
                        return Err(PropertySupportError::DiscreteSupport {
                            axis: abscissa.clone(),
                            lower: a.min(b),
                            upper: a.max(b),
                            claim: selected,
                        });
                    }
                    for &(x, _) in knots {
                        if a.min(b) < x && x < a.max(b) {
                            if parameters.len() >= MAX_PROPERTY_LINE_KNOTS {
                                return Err(invalid(
                                    "source path exceeds the original breakpoint budget",
                                ));
                            }
                            let t = (x - a) / (b - a);
                            if !(0.0 < t && t < 1.0) {
                                return Err(invalid(
                                    "source curve knot collapses onto a distinct path endpoint",
                                ));
                            }
                            // Preserve the source abscissa bits in its receipt;
                            // other physical coordinates use the same path t.
                            parameters.push((t, Some((abscissa, x))));
                        }
                    }
                }
            }
            PropertyValue::Surface2D {
                axes,
                samples,
                triangles,
                ..
            } => {
                let a = [start.axes()[&axes[0]], start.axes()[&axes[1]]];
                let b = [end.axes()[&axes[0]], end.axes()[&axes[1]]];
                if a != b {
                    if claim.interpolation != InterpolationPolicy::LinearInside {
                        return Err(invalid(
                            "a varying surface path needs original affine cell interpolation",
                        ));
                    }
                    parameters.extend(
                        crate::surface::segment_knots(samples, triangles, a, b)
                            .map_err(|error| PropertySupportError::Evaluation {
                                point: start.clone(),
                                error: Box::new(error),
                            })?
                            .into_iter()
                            .map(|t| (t, None)),
                    );
                }
            }
        }
        if parameters.len() > MAX_PROPERTY_LINE_KNOTS {
            return Err(invalid(
                "source path exceeds the original breakpoint budget",
            ));
        }
        parameters.sort_by(|a, b| a.0.total_cmp(&b.0));
        if parameters.windows(2).any(|pair| pair[0].0 == pair[1].0) {
            return Err(invalid(
                "distinct source path knots collapse to one parameter",
            ));
        }
        let mut knots = Vec::with_capacity(parameters.len());
        for (t, source_coordinate) in parameters {
            let answer = if t == 0.0 {
                first.clone()
            } else if t == 1.0 {
                last.clone()
            } else {
                let point = affine_point(start, end, t)
                    .and_then(|point| {
                        if let Some((axis, value)) = source_coordinate {
                            point.with(axis, value)
                        } else {
                            Ok(point)
                        }
                    })
                    .map_err(|error| PropertySupportError::Evaluation {
                        point: start.clone(),
                        error: Box::new(error),
                    })?;
                query(&point)?
            };
            knots.push((t, answer));
        }
        Ok(PropertyPathSlice {
            start: start.clone(),
            end: end.clone(),
            knots,
        })
    }
}

fn affine_point(start: &QueryPoint, end: &QueryPoint, t: f64) -> Result<QueryPoint, MatDbError> {
    let mut point = start.clone();
    for (axis, &a) in start.axes() {
        let b = end.axes()[axis];
        if a != b {
            let value = if t <= 0.5 {
                (b - a).mul_add(t, a)
            } else {
                (a - b).mul_add(1.0 - t, b)
            };
            point = point.with(axis, value)?;
        }
    }
    Ok(point)
}
