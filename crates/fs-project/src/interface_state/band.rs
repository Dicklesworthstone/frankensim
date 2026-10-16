//! Resistance extrema over a declared manufacturing band, not a probability law.
use super::{StatePoint, continuous_coordinate, invalid, query_envelope, query_point};
use crate::{CONTACT_RESISTANCE_DIMS, CONTACT_RESISTANCE_PROPERTY, InterfaceState, ProjectError};
use fs_matdb::{ClaimId, ClaimSelection, ClaimSet, InterpolationPolicy, PropertyValue};

/// Source-curve extrema at fixed geometry, temperature-independent resistance
/// and the original categorical state. These are floating-point source queries,
/// not outward-rounded bounds on physical contact or on a thermal response.
#[derive(Debug, Clone, PartialEq)]
pub struct ResistanceBand {
    pub axis: &'static str,
    pub unit: &'static str,
    pub coordinate_low: f64,
    pub coordinate_high: f64,
    pub nominal_resistance: f64,
    pub minimum_resistance: f64,
    pub maximum_resistance: f64,
    pub claim: ClaimId,
    /// Point declarations: the original band is not applied a second time.
    pub minimum_state: InterfaceState,
    pub maximum_state: InterfaceState,
}

fn point_state(state: &InterfaceState, x: f64) -> InterfaceState {
    let mut point = state.clone();
    let (value, width) = match &mut point {
        InterfaceState::DryContact { pressure, pressure_half_width, .. } => (pressure, pressure_half_width),
        InterfaceState::Tim { thickness, thickness_half_width }
        | InterfaceState::Adhesive { thickness, thickness_half_width } => (thickness, thickness_half_width),
        InterfaceState::GapWithFluid { gap, gap_half_width, .. } => (gap, gap_half_width),
        InterfaceState::BoltedWithPattern { torque, torque_half_width, .. } => (torque, torque_half_width),
    };
    value.value = x;
    width.value = 0.0;
    point
}

/// Query both band ends AND every enclosed piecewise-linear source knot.
/// Unequal slopes are admitted: propagation is not differentiation. Full-box
/// support and claim selection are checked before scanning; every candidate
/// must reproduce that same positive, dimensionally correct resistance claim.
/// A temperature-varying contact law refuses, matching native steady admission.
/// Scalar and fixed-other-axis sources remain constant in the varying state.
///
/// `max_candidates` bounds source knots inspected as well as query count.
/// No collection of query points or changed claims is constructed. A refused
/// query or exhausted allowance returns no partial band. The caller owns
/// cancellation between these bounded source scans and subsequent solves.
pub fn resistance_band(
    claims: &ClaimSet, state: &InterfaceState, temperature_axis: &str,
    low_k: f64, high_k: f64, selection: ClaimSelection, max_candidates: usize,
) -> Result<ResistanceBand, ProjectError> {
    if max_candidates < 3 { return Err(invalid("joint-band source allowance needs at least three candidates")); }
    let answer = query_envelope(claims, CONTACT_RESISTANCE_PROPERTY, state,
        temperature_axis, low_k, high_k, selection).map_err(|e| invalid(e.to_string()))?;
    let claim_id = answer.lower.receipt.selected;
    let claim = claims.claim(claim_id).ok_or_else(|| invalid("selected contact claim is absent"))?;
    let nominal = &answer.lower.evidence.value;
    if nominal.dims != CONTACT_RESISTANCE_DIMS || !nominal.value.is_finite() || nominal.value <= 0.0
        || answer.upper.evidence.value.value.to_bits() != nominal.value.to_bits()
    { return Err(invalid("joint-band propagation requires positive temperature-independent contact resistance")); }
    let (axis, unit) = continuous_coordinate(state);
    let point = query_point(state, temperature_axis, low_k, StatePoint::Nominal)?;
    let low = query_point(state, temperature_axis, low_k, StatePoint::Lower)?.axes()[axis];
    let high = query_point(state, temperature_axis, low_k, StatePoint::Upper)?.axes()[axis];
    let x0 = point.axes()[axis];
    let (abscissa, knots) = match (&claim.value, claim.interpolation) {
        (PropertyValue::Scalar { .. }, InterpolationPolicy::ConstantWithinValidity) => ("", &[][..]),
        (PropertyValue::Curve { abscissa, knots, .. }, InterpolationPolicy::LinearInside) => {
            if knots.len() > max_candidates.saturating_sub(3) {
                return Err(invalid("joint-band source-knot allowance exhausted"));
            }
            if abscissa == temperature_axis && knots.iter().any(|&(t, r)|
                t >= low_k && t <= high_k && r.to_bits() != nominal.value.to_bits())
            { return Err(invalid("joint-band source varies inside the native temperature range")); }
            (abscissa.as_str(), knots.as_slice())
        }
        _ => return Err(invalid("joint-band propagation requires continuous source interpolation, not exact samples")),
    };
    let mut minimum = (nominal.value, x0);
    let mut maximum = minimum;
    let candidates = [low, high].into_iter().chain(knots.iter()
        .filter(|(x, _)| abscissa == axis && *x > low && *x < high).map(|(x, _)| *x));
    for x in candidates {
        let candidate = point.clone().with(axis, x).map_err(|e| invalid(e.to_string()))?;
        let sampled = claims.query_envelope(CONTACT_RESISTANCE_PROPERTY,
            &candidate, &candidate, selection).map_err(|e| invalid(e.to_string()))?;
        let value = &sampled.lower.evidence.value;
        if sampled.lower.receipt.selected != claim_id || value.dims != CONTACT_RESISTANCE_DIMS
            || !value.value.is_finite() || value.value <= 0.0
        { return Err(invalid("joint-band candidate changed its source or has invalid resistance")); }
        // Retain the nominal coordinate on equal values. A flat claim need
        // not move the nominal manufacturing state just to report zero effect.
        if value.value < minimum.0 { minimum = (value.value, x); }
        if value.value > maximum.0 { maximum = (value.value, x); }
    }
    Ok(ResistanceBand {
        axis, unit, coordinate_low: low, coordinate_high: high, nominal_resistance: nominal.value,
        minimum_resistance: minimum.0, maximum_resistance: maximum.0, claim: claim_id,
        minimum_state: point_state(state, minimum.1), maximum_state: point_state(state, maximum.1),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_matdb::{PropertyClaim, PropertyKey, Provenance, SelectionPolicy, UncertaintyModel};
    use fs_qty::QtyAny;
    const SINGLE: ClaimSelection = ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly);
    fn state() -> InterfaceState {
        InterfaceState::DryContact { pressure: QtyAny::new(3.0, crate::spec::dims::PRESSURE),
            pressure_half_width: QtyAny::new(2.0, crate::spec::dims::PRESSURE), finish: "kept-finish".into() }
    }
    fn card(knots: Vec<(f64, f64)>) -> ClaimSet {
        let mut claims = ClaimSet::new();
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(CONTACT_RESISTANCE_PROPERTY, CONTACT_RESISTANCE_DIMS),
            validity: fs_evidence::ValidityDomain::unconstrained().with("T", 200.0, 450.0)
                .with("normal_pressure", 1.0, 5.0),
            value: PropertyValue::Curve { abscissa: "normal_pressure".into(),
                abscissa_dims: crate::spec::dims::PRESSURE, knots, dims: CONTACT_RESISTANCE_DIMS },
            uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::LinearInside,
            observations: Vec::new(), provenance: Provenance { source: "synthetic band, not measurements".into(),
                license: "CC0-1.0".into(), artifact: None },
        }).unwrap();
        claims
    }
    #[test]
    fn interior_kinks_are_real_band_extrema_not_an_endpoint_or_derivative_shortcut() {
        let claims = card(vec![(1.0, 0.2), (2.0, 0.6), (3.0, 0.3), (4.0, 0.1), (5.0, 0.2)]);
        let original = state();
        let band = resistance_band(&claims, &original, "T", 250.0, 400.0, SINGLE, 8).unwrap();
        assert_eq!((band.minimum_resistance, band.maximum_resistance), (0.1, 0.6));
        assert_eq!(band.minimum_state, point_state(&original, 4.0));
        assert_eq!(band.maximum_state, point_state(&original, 2.0));
        assert_eq!((band.axis, band.unit), ("normal_pressure", "Pa"));
        assert_eq!(original, state());
        // Zero half-width at a point avoids imposing [x-width,x+width] AGAIN.
        for s in [&band.minimum_state, &band.maximum_state] {
            query_envelope(&claims, CONTACT_RESISTANCE_PROPERTY, s, "T", 250.0, 400.0, SINGLE).unwrap();
        }
    }
    #[test]
    fn constant_band_keeps_nominal_coordinate_and_source_refusals_return_no_partial_range() {
        let claims = card(vec![(1.0, 0.2), (5.0, 0.2)]);
        let band = resistance_band(&claims, &state(), "T", 250.0, 400.0, SINGLE, 5).unwrap();
        assert_eq!(band.minimum_state, point_state(&state(), 3.0));
        assert_eq!(band.minimum_state, band.maximum_state);
        assert_eq!(band.minimum_resistance, band.nominal_resistance);
        assert!(resistance_band(&claims, &state(), "T", 250.0, 400.0, SINGLE, 4).is_err());
        assert!(resistance_band(&claims, &state(), "T", 190.0, 400.0, SINGLE, 5).is_err());
        let mut wide = state();
        if let InterfaceState::DryContact { pressure_half_width, .. } = &mut wide { pressure_half_width.value = 2.5; }
        assert!(resistance_band(&claims, &wide, "T", 250.0, 400.0, SINGLE, 5).is_err());
    }
    #[test]
    fn nonpositive_interior_values_and_temperature_variation_cannot_be_hidden_by_endpoints() {
        let bad = card(vec![(1.0, 0.2), (3.0, -0.1), (5.0, 0.2)]);
        assert!(resistance_band(&bad, &state(), "T", 250.0, 400.0, SINGLE, 6).is_err());
        let mut claim = card(vec![(1.0, 0.2), (5.0, 0.2)]).claims_for(CONTACT_RESISTANCE_PROPERTY)[0].1.clone();
        claim.value = PropertyValue::Curve { abscissa: "T".into(),
            abscissa_dims: crate::spec::dims::TEMPERATURE,
            knots: vec![(250.0, 0.2), (300.0, 0.5), (400.0, 0.2)], dims: CONTACT_RESISTANCE_DIMS };
        let mut claims = ClaimSet::new(); claims.insert_claim(claim).unwrap();
        assert!(resistance_band(&claims, &state(), "T", 250.0, 400.0, SINGLE, 6).is_err());
    }
}
