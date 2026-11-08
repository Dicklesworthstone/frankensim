//! Conservative feedback arithmetic for an original affine emissivity path.

/// Signed slope with its magnitude rounded outward [kg/J]. The sensible
/// argument is (whole-interval delta_T, minimum ORIGINAL Cp), when retained.
pub(super) fn emissivity_enthalpy_slope(
    emissivity: [f64; 2],
    parameter_width: f64,
    enthalpy_span: f64,
    sensible: Option<(f64, f64)>,
) -> Option<f64> {
    if emissivity
        .iter()
        .any(|e| !e.is_finite() || *e <= 0.0 || *e > 1.0)
        || !parameter_width.is_finite()
        || parameter_width <= 0.0
        || parameter_width > 1.0
        || !enthalpy_span.is_finite()
        || enthalpy_span <= 0.0
    {
        return None;
    }
    let minimum_span = if let Some((temperature_span, capacity)) = sensible {
        if !temperature_span.is_finite()
            || temperature_span <= 0.0
            || !capacity.is_finite()
            || capacity <= 0.0
        {
            return None;
        }
        let span = temperature_span.next_down() * capacity;
        if !span.is_finite() {
            return None;
        }
        span.next_down()
    } else {
        enthalpy_span.next_down()
    };
    let change = emissivity[1] - emissivity[0];
    if change == 0.0 {
        return Some(0.0);
    }
    let width = (parameter_width.next_down() * minimum_span).next_down();
    if !width.is_finite() || width <= 0.0 {
        return None;
    }
    let magnitude =
        ((change.abs().next_up() / width).next_up() * (1.0 + 64.0 * f64::EPSILON)).next_up();
    (magnitude.is_finite() && magnitude > 0.0).then_some(magnitude.copysign(change))
}

/// Bound positive inward dP_radiation/dh on one original source interval.
/// The negative -4*epsilon*sigma*A*T^3*dT/dh term is deliberately dropped.
pub(super) fn positive_radiative_feedback(
    emissivity_slope: f64,
    temperatures_k: [f64; 2],
    area_m2: f64,
    reservoir_k: f64,
    sigma: f64,
) -> Option<f64> {
    if !emissivity_slope.is_finite()
        || temperatures_k.iter().any(|t| !t.is_finite() || *t <= 0.0)
        || temperatures_k[0] > temperatures_k[1]
        || !area_m2.is_finite()
        || area_m2 <= 0.0
        || !reservoir_k.is_finite()
        || reservoir_k <= 0.0
        || !sigma.is_finite()
        || sigma <= 0.0
    {
        return None;
    }
    if emissivity_slope == 0.0 {
        return Some(0.0);
    }
    let temperature = if emissivity_slope > 0.0 {
        temperatures_k[0]
    } else {
        temperatures_k[1]
    };
    let contrast = if emissivity_slope > 0.0 {
        reservoir_k - temperature
    } else {
        temperature - reservoir_k
    };
    if contrast <= 0.0 {
        return Some(0.0);
    }
    // Same factored fourth-power difference as the native heat owner. Each
    // positive operation rounds outwards; no heat is evaluated by this bound.
    let sum = (temperature + reservoir_k).next_up();
    let squares =
        ((temperature * temperature).next_up() + (reservoir_k * reservoir_k).next_up()).next_up();
    let mut bound = sigma;
    for factor in [
        area_m2,
        emissivity_slope.abs(),
        contrast.next_up(),
        sum,
        squares,
    ] {
        bound = (bound * factor).next_up();
    }
    (bound.is_finite() && bound > 0.0).then_some(bound)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn bounded(actual: f64, exact: f64) {
        assert!(actual >= exact, "{actual} < {exact}");
        assert!(actual <= exact * (1.0 + 1.0e-12), "{actual} > {exact}");
    }

    #[test]
    fn g1_latent_feedback_keeps_the_inward_heat_sign_and_equilibrium() {
        let decreasing = emissivity_enthalpy_slope([0.8, 0.2], 1.0, 100.0, None).unwrap();
        let increasing = -decreasing;
        bounded(-decreasing, 0.006);
        bounded(
            positive_radiative_feedback(decreasing, [400.0; 2], 2.0, 300.0, 1.0).unwrap(),
            2.0 * 0.006 * (400.0_f64.powi(4) - 300.0_f64.powi(4)),
        );
        bounded(
            positive_radiative_feedback(increasing, [400.0; 2], 2.0, 500.0, 1.0).unwrap(),
            2.0 * 0.006 * (500.0_f64.powi(4) - 400.0_f64.powi(4)),
        );
        for (slope, bath) in [
            (increasing, 300.0),
            (decreasing, 500.0),
            (decreasing, 400.0),
        ] {
            assert_eq!(
                positive_radiative_feedback(slope, [400.0; 2], 2.0, bath, 1.0),
                Some(0.0)
            );
        }
    }

    #[test]
    fn g1_retained_cp_minimum_bounds_the_actual_nonlinear_clock() {
        // Cp=2→8 over 100 K stores 500 J/kg. Its minimum, not the chord
        // capacity 5, sets the maximum path speed and epsilon slope 0.003.
        let slope = emissivity_enthalpy_slope([0.2, 0.8], 1.0, 500.0, Some((100.0, 2.0))).unwrap();
        bounded(slope, 0.003);
        bounded(
            positive_radiative_feedback(slope, [300.0, 400.0], 1.0, 500.0, 1.0).unwrap(),
            0.003 * (500.0_f64.powi(4) - 300.0_f64.powi(4)),
        );
    }

    #[test]
    fn g0_zero_feedback_stays_zero_and_invalid_bounds_refuse() {
        assert_eq!(
            emissivity_enthalpy_slope([0.5; 2], 1.0, 100.0, None),
            Some(0.0)
        );
        assert_eq!(
            positive_radiative_feedback(0.0, [f64::MAX; 2], 1.0, f64::MAX, 1.0),
            Some(0.0)
        );
        assert!(emissivity_enthalpy_slope([0.2, 0.8], 1.0, 100.0, Some((f64::MAX, 2.0))).is_none());
        for invalid in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(positive_radiative_feedback(0.0, [300.0; 2], invalid, 300.0, 1.0).is_none());
            assert!(emissivity_enthalpy_slope([0.2, 0.8], invalid, 100.0, None).is_none());
        }
    }
}
