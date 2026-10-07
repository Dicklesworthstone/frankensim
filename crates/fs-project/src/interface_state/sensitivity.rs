//! Smooth card response to the continuous coordinate of a manufactured joint.
//! This is a constitutive partial, not a thermal solve or uncertainty bound.
use super::{StatePoint, invalid, query_envelope, query_point};
use crate::{CONTACT_RESISTANCE_PROPERTY, CONTACT_RESISTANCE_DIMS, InterfaceState, ProjectError};
use fs_matdb::{ClaimId, ClaimSelection, ClaimSet, InterpolationPolicy, PropertyValue};

/// Response of the selected resistance claim at the nominal manufactured state.
#[derive(Debug, Clone, PartialEq)]
pub struct StateSensitivity {
    pub axis: &'static str,
    pub unit: &'static str,
    pub nominal: f64,
    /// Area-specific resistance, m² K/W, from the unchanged selected card.
    pub resistance: f64,
    /// dR''/dx in (m² K/W) per coherent manufactured-coordinate unit.
    pub derivative: f64,
    pub claim: ClaimId,
}

/// The one continuous coordinate of an interface declaration. Categorical
/// finish/fluid/pattern and discrete bolt count are deliberately not controls.
#[must_use]
pub fn continuous_coordinate(state: &InterfaceState) -> (&'static str, &'static str) {
    match state {
        InterfaceState::DryContact { .. } => ("normal_pressure", "Pa"),
        InterfaceState::Tim { .. } | InterfaceState::Adhesive { .. } => ("thickness", "m"),
        InterfaceState::GapWithFluid { .. } => ("gap", "m"),
        InterfaceState::BoltedWithPattern { .. } => ("torque", "N*m"),
    }
}

/// Differentiate the same claim admitted over the temperature/manufacturing
/// box. No source mutation, numerical differencing, axis alias or mechanical
/// conversion occurs. The derivative is evaluated at `low_k`; the native
/// steady contact consumer separately requires temperature-independent R''.
///
/// # Errors
/// Existing source/units/support/selection refusals, nonpositive resistance,
/// exact-sample-only data, validity edges, unequal-slope knots or nonfinite
/// arithmetic. A zero-width manufacturing band does not prove smoothness:
/// both adjacent representable query coordinates are admitted separately.
pub fn resistance_sensitivity(
    claims: &ClaimSet, state: &InterfaceState, temperature_axis: &str,
    low_k: f64, high_k: f64, selection: ClaimSelection,
) -> Result<StateSensitivity, ProjectError> {
    let answer = query_envelope(claims, CONTACT_RESISTANCE_PROPERTY, state,
        temperature_axis, low_k, high_k, selection).map_err(|e| invalid(e.to_string()))?;
    let claim_id = answer.lower.receipt.selected;
    let claim = claims.claim(claim_id).ok_or_else(|| invalid("selected contact claim is absent"))?;
    let sample = &answer.lower.evidence.value;
    if sample.dims != CONTACT_RESISTANCE_DIMS || !sample.value.is_finite() || sample.value <= 0.0 {
        return Err(invalid("contact response requires positive finite area-specific resistance"));
    }
    let (axis, unit) = continuous_coordinate(state);
    let point = query_point(state, temperature_axis, low_k, StatePoint::Nominal)?;
    let x = point.axes()[axis];
    if claim.validity.bound(axis).is_some_and(|(lo, hi)| x <= lo || x >= hi) {
        return Err(invalid("a two-sided joint-state derivative is unavailable at a source validity edge"));
    }
    if !(x.next_down() > 0.0 && x.next_up().is_finite()) {
        return Err(invalid("joint-state derivative has no finite positive two-sided neighborhood"));
    }
    let lower = point.clone().with(axis, x.next_down()).map_err(|e| invalid(e.to_string()))?;
    let upper = point.with(axis, x.next_up()).map_err(|e| invalid(e.to_string()))?;
    let neighborhood = claims.query_envelope(CONTACT_RESISTANCE_PROPERTY, &lower, &upper, selection)
        .map_err(|e| invalid(e.to_string()))?;
    if neighborhood.lower.receipt.selected != claim_id {
        return Err(invalid("joint-state derivative changes the admitted source claim"));
    }
    let derivative = match (&claim.value, claim.interpolation) {
        (PropertyValue::Scalar { .. }, InterpolationPolicy::ConstantWithinValidity) => 0.0,
        (PropertyValue::Curve { abscissa, knots, .. }, InterpolationPolicy::LinearInside)
            if abscissa == axis => slope(knots, x)?,
        (PropertyValue::Curve { .. }, InterpolationPolicy::LinearInside) => 0.0,
        _ => return Err(invalid("joint-state derivatives need continuous source interpolation, not exact samples")),
    };
    if !derivative.is_finite() {
        return Err(invalid("joint-state derivative arithmetic is nonfinite"));
    }
    Ok(StateSensitivity { axis, unit, nominal: x, resistance: sample.value,
        derivative, claim: claim_id })
}

fn slope(knots: &[(f64, f64)], x: f64) -> Result<f64, ProjectError> {
    if knots.len() < 2 || x <= knots[0].0 || x >= knots[knots.len() - 1].0 {
        return Err(invalid("joint-state derivative needs an interior source-curve point"));
    }
    let segment = |i: usize| -> Result<f64, ProjectError> {
        let width = knots[i + 1].0 - knots[i].0;
        let delta = knots[i + 1].1 - knots[i].1;
        if !width.is_finite() || width <= 0.0 || !delta.is_finite() {
            return Err(invalid("joint-state source-curve differences are not finite"));
        }
        let value = delta / width;
        if value.is_finite() { Ok(value) } else { Err(invalid("joint-state slope is not representable")) }
    };
    let right = knots.partition_point(|(at, _)| *at < x);
    let left = segment(right - 1)?;
    if knots[right].0 == x && segment(right)? != left {
        return Err(invalid("joint-state derivative is undefined at an unequal-slope source knot"));
    }
    Ok(left)
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_matdb::{PropertyClaim, PropertyKey, Provenance, SelectionPolicy, UncertaintyModel};
    use fs_qty::QtyAny;
    const SINGLE: ClaimSelection = ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly);
    fn state(x: f64) -> InterfaceState {
        InterfaceState::DryContact { pressure: QtyAny::new(x, crate::spec::dims::PRESSURE),
            pressure_half_width: QtyAny::new(0.0, crate::spec::dims::PRESSURE), finish: "fixture".into() }
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
            observations: Vec::new(), provenance: Provenance { source: "synthetic response".into(),
                license: "CC0-1.0".into(), artifact: None },
        }).unwrap();
        claims
    }
    #[test]
    fn response_is_the_original_card_value_and_analytic_coordinate_partial() {
        let claims = card(vec![(1.0, 0.5), (3.0, 0.25), (5.0, 0.125)]);
        for (x, resistance, derivative) in [(2.0, 0.375, -0.125), (4.0, 0.1875, -0.0625)] {
            let actual = resistance_sensitivity(&claims, &state(x), "T", 250.0, 400.0, SINGLE).unwrap();
            assert_eq!(actual.axis, "normal_pressure"); assert_eq!(actual.unit, "Pa");
            assert_eq!(actual.nominal, x); assert_eq!(actual.resistance, resistance);
            assert_eq!(actual.derivative, derivative);
        }
    }
    #[test]
    fn equal_slopes_are_smooth_but_kinks_and_two_sided_edges_refuse() {
        let kink = card(vec![(1.0, 0.5), (3.0, 0.25), (5.0, 0.125)]);
        for x in [1.0, 3.0, 5.0] {
            assert!(resistance_sensitivity(&kink, &state(x), "T", 250.0, 400.0, SINGLE).is_err());
        }
        let smooth = card(vec![(1.0, 0.75), (3.0, 0.5), (5.0, 0.25)]);
        assert_eq!(resistance_sensitivity(&smooth, &state(3.0), "T", 250.0, 400.0, SINGLE)
            .unwrap().derivative, -0.125);
    }
    #[test]
    fn constant_contact_is_zero_not_an_unsupported_or_manufactured_slope() {
        let mut claims = ClaimSet::new();
        let mut claim = card(vec![(1.0, 0.5), (5.0, 0.25)])
            .claims_for(CONTACT_RESISTANCE_PROPERTY)[0].1.clone();
        claim.value = PropertyValue::Scalar { value: 0.13, dims: CONTACT_RESISTANCE_DIMS };
        claim.interpolation = InterpolationPolicy::ConstantWithinValidity;
        claims.insert_claim(claim).unwrap();
        assert_eq!(resistance_sensitivity(&claims, &state(3.0), "T", 250.0, 400.0, SINGLE)
            .unwrap().derivative, 0.0);
        assert!(resistance_sensitivity(&claims, &state(1.0), "T", 250.0, 400.0, SINGLE).is_err());
    }
}
