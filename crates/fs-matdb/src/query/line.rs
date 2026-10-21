//! Fixed-context scalar source paths for existing physical integral owners.

use super::{
    ClaimSelection, ClaimSet, MaterialAnswer, PropertyKey, PropertySupportError, PropertyValue,
    QueryPoint,
};

/// Maximum retained breakpoints in one material source line.
pub const MAX_PROPERTY_LINE_KNOTS: usize = 4_096;

/// One continuously supported, piecewise-affine material law along one axis.
///
/// Each knot is a canonical point answer with its original source receipt.
/// Coordinates are increasing and include both requested endpoints; a zero
/// span contains one knot. This retains source interpolation, not an integral,
/// derivative, uncertainty propagation, or certified floating-point enclosure.
#[derive(Debug, Clone, PartialEq)]
pub struct PropertyLineSlice {
    /// The only coordinate permitted to vary.
    pub axis: String,
    /// Source breakpoints and their canonical material answers.
    pub knots: Vec<(f64, MaterialAnswer)>,
}

impl ClaimSet {
    /// Resolve a source law along one ordered coordinate interval while every
    /// other coordinate and descriptor remains fixed.
    ///
    /// Scalar and curve support/selection use the ordinary envelope owner.
    /// For a sampled surface, declared cells must connect the endpoints along
    /// this exact fixed-coordinate line. Shared source-edge/vertex identities,
    /// never epsilon-merged interval endpoints, establish that continuity.
    /// Every crossed cell edge is retained as a breakpoint. Source holes,
    /// discrete-only spans, source switching, and unresolved cuts refuse.
    /// Ordinary multidimensional envelope semantics remain unchanged.
    ///
    /// # Errors
    /// Returns the existing support and canonical point-query errors. More
    /// than one varying coordinate, an absent axis, or a reversed interval is
    /// an invalid request. A zero interval still requires a valid point query.
    pub fn query_line_typed(
        &self,
        property: &PropertyKey,
        axis: &str,
        lower: &QueryPoint,
        upper: &QueryPoint,
        selection: ClaimSelection,
    ) -> Result<PropertyLineSlice, PropertySupportError> {
        let invalid = |reason: &str| PropertySupportError::InvalidEnvelope {
            reason: reason.to_owned(),
        };
        let (Some(&lo), Some(&hi)) = (lower.axes().get(axis), upper.axes().get(axis)) else {
            return Err(invalid(
                "a source line needs its named coordinate at both endpoints",
            ));
        };
        if lower.axes().keys().ne(upper.axes().keys())
            || lower.axis_quantities() != upper.axis_quantities()
            || lo > hi
            || lower.axes().iter().any(|(name, value)| {
                name != axis && value.to_bits() != upper.axes()[name].to_bits()
            })
        {
            return Err(invalid(
                "a source line requires one ordered axis and identical remaining coordinates and quantities",
            ));
        }
        self.require_quantity(property)
            .map_err(|error| PropertySupportError::Evaluation {
                point: lower.clone(),
                error: Box::new(error),
            })?;
        let (support, mut coordinates) = self.query_envelope_selected_with_line(
            property.name(),
            Some(property),
            lower,
            upper,
            selection,
            Some(axis),
        )?;
        let selected = support.lower.receipt.selected;
        let claim = self
            .claim(selected)
            .expect("support selected an immutable claim");
        if let PropertyValue::Curve {
            abscissa, knots, ..
        } = &claim.value
            && abscissa == axis
        {
            for &(x, _) in knots {
                if lo < x && x < hi {
                    if coordinates.len() >= MAX_PROPERTY_LINE_KNOTS - 2 {
                        return Err(invalid("source line exceeds its breakpoint budget"));
                    }
                    coordinates.push(x);
                }
            }
        }
        coordinates.extend([lo, hi]);
        coordinates.sort_by(f64::total_cmp);
        coordinates.dedup();
        if coordinates.len() > MAX_PROPERTY_LINE_KNOTS {
            return Err(invalid("source line exceeds its breakpoint budget"));
        }
        let mut knots = Vec::with_capacity(coordinates.len());
        for coordinate in coordinates {
            let point = lower.clone().with(axis, coordinate).map_err(|error| {
                PropertySupportError::Evaluation {
                    point: lower.clone(),
                    error: Box::new(error),
                }
            })?;
            let answer = self
                .query_selected(property.name(), Some(property), &point, selection)
                .and_then(|answer| {
                    self.verify_receipt(&answer.receipt)?;
                    Ok(answer)
                })
                .map_err(|error| PropertySupportError::Evaluation {
                    point: point.clone(),
                    error: Box::new(error),
                })?;
            if answer.receipt.selected != selected {
                return Err(PropertySupportError::ClaimChanges {
                    point,
                    selected,
                    other: answer.receipt.selected,
                });
            }
            knots.push((coordinate, answer));
        }
        Ok(PropertyLineSlice {
            axis: axis.to_owned(),
            knots,
        })
    }
}
