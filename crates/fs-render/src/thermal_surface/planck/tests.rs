use super::*;
use crate::thermal_surface::{SIGMA_W_M2_K4, planck_radiance_w_m2_sr_nm};

fn relative_close(actual: f64, expected: f64, tolerance: f64) {
    assert!(actual.is_finite());
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs() + 2.0 * f64::from_bits(1),
        "actual={actual:.17e}, expected={expected:.17e}"
    );
}

#[test]
fn adjacent_temperatures_match_independent_high_precision_si_references() {
    // 120-digit decimal SI Planck evaluations, starting from the exact binary64
    // input temperatures/wavelengths, then rounded once for these references.
    // They independently cover the thermal IR, visible/UV Wien, and RJ limits.
    let cases: [(f64, f64, f64, f64, f64); 3] = [
        (
            300.0,
            10_000.0,
            9.093330803086814e-18,
            1.5997156725132198e-4,
            1.5997156725132194e-4,
        ),
        (
            6000.0,
            100.0,
            1.6680314834091934e-15,
            1.8340200115048581e-3,
            1.834020011504855e-3,
        ),
        (
            1e8,
            1e9,
            1.233542434409982e-31,
            8.27816314690484e-24,
            8.27816314690484e-24,
        ),
    ];
    for (cold, wavelength, difference, secant, tangent) in cases {
        let hot = f64::from_bits(cold.to_bits() + 1);
        let forward = planck_radiance_difference_w_m2_sr_nm(hot, cold, wavelength).unwrap();
        let reverse = planck_radiance_difference_w_m2_sr_nm(cold, hot, wavelength).unwrap();
        assert!(forward > 0.0);
        assert_eq!(forward.to_bits(), (-reverse).to_bits());
        relative_close(forward, difference, 3e-12);
        let d = planck_radiance_temperature_secant_w_m2_sr_nm_k(hot, cold, wavelength).unwrap();
        relative_close(d, secant, 3e-12);
        assert_eq!(
            d,
            planck_radiance_temperature_secant_w_m2_sr_nm_k(cold, hot, wavelength).unwrap()
        );
        relative_close(
            planck_radiance_temperature_derivative_w_m2_sr_nm_k(cold, wavelength).unwrap(),
            tangent,
            3e-12,
        );
    }
}

#[test]
fn finite_temperature_changes_and_deep_wien_slopes_match_si_references() {
    for (cold, hot, wavelength, difference) in [
        (299.0, 301.0, 10_000.0, 3.1994367088632835e-4),
        (1000.0, 2000.0, 500.0, 2.150546430117862),
    ] {
        relative_close(
            planck_radiance_difference_w_m2_sr_nm(hot, cold, wavelength).unwrap(),
            difference,
            3e-12,
        );
        let slope = planck_radiance_temperature_secant_w_m2_sr_nm_k(hot, cold, wavelength).unwrap();
        relative_close(slope, difference / (hot - cold), 3e-12);
        let lower = planck_radiance_temperature_derivative_w_m2_sr_nm_k(cold, wavelength).unwrap();
        let upper = planck_radiance_temperature_derivative_w_m2_sr_nm_k(hot, wavelength).unwrap();
        // The Planck source is increasing and convex in kelvin; a secant lies
        // between its endpoint tangents, which bound a passive cooling step.
        assert!(lower <= slope && slope <= upper);
    }
    relative_close(
        planck_radiance_temperature_derivative_w_m2_sr_nm_k(100.0, 1000.0).unwrap(),
        5.605766728980531e-58,
        3e-12,
    );
}

#[test]
fn integrated_difference_and_tangent_recover_bolometric_energy_and_slope() {
    // Independent global blackbody identities test SI/nm and derivative factors.
    // Integrate the wavelength APIs through x=hc/(lambda*k*T), including its
    // Jacobian and pi for Lambertian hemispherical exitance.
    let cold: f64 = 300.0;
    let hot = f64::from_bits(cold.to_bits() + 1);
    let n = 16_384;
    let dx = 50.0 / f64::from(n);
    let mut difference = 0.0;
    let mut derivative = 0.0;
    for i in 0..n {
        let x = (f64::from(i) + 0.5) * dx;
        let wavelength = C2_NM_K / (cold * x);
        let jacobian = core::f64::consts::PI * C2_NM_K / (cold * x * x) * dx;
        difference +=
            planck_radiance_difference_w_m2_sr_nm(hot, cold, wavelength).unwrap() * jacobian;
        derivative += planck_radiance_temperature_derivative_w_m2_sr_nm_k(cold, wavelength)
            .unwrap()
            * jacobian;
    }
    let exact_difference = SIGMA_W_M2_K4 * (hot - cold) * (hot + cold) * (hot * hot + cold * cold);
    relative_close(difference, exact_difference, 2e-9);
    relative_close(derivative, 4.0 * SIGMA_W_M2_K4 * cold.powi(3), 2e-9);
}

#[test]
fn scaled_heat_survives_underflowed_spectra_and_unrepresentable_intermediate_products() {
    let cold: f64 = 1e-65;
    let hot = f64::from_bits(cold.to_bits() + 1);
    let wavelength = 1e71;
    assert_eq!(planck_radiance_w_m2_sr_nm(cold, wavelength), Ok(0.0));
    assert_eq!(
        planck_radiance_difference_w_m2_sr_nm(hot, cold, wavelength),
        Ok(0.0)
    );
    relative_close(
        planck_radiance_temperature_secant_w_m2_sr_nm_k(hot, cold, wavelength).unwrap(),
        9.66924990440166e-276,
        4e-12,
    );
    // K*weight overflows, and unscaled Delta B underflows, but both final rates
    // are ordinary representable numbers; neither intermediate may erase heat.
    let (heat, slope) =
        planck_radiative_exchange_w_and_w_k(hot, Some(cold), wavelength, 1e300, 1e100).unwrap();
    relative_close(heat, 2.0387029290953024e44, 4e-12);
    relative_close(slope, 9.66924990440166e124, 4e-12);
    let reversed =
        planck_radiative_exchange_w_and_w_k(cold, Some(hot), wavelength, 1e300, 1e100).unwrap();
    assert_eq!(heat.to_bits(), (-reversed.0).to_bits());
    assert_eq!(slope, reversed.1);

    let cold: f64 = 1e65;
    let hot = f64::from_bits(cold.to_bits() + 1);
    let wavelength = 3e-59;
    assert_eq!(
        planck_radiance_w_m2_sr_nm(cold, wavelength),
        Err(ThermalSurfaceError::Unrepresentable)
    );
    relative_close(
        planck_radiance_difference_w_m2_sr_nm(hot, cold, wavelength).unwrap(),
        2.309123036683434e295,
        4e-12,
    );
    // Conversely, K*weight underflows while the final physical power does not.
    let (heat, slope) =
        planck_radiative_exchange_w_and_w_k(hot, Some(cold), wavelength, 1e-250, 1e-100).unwrap();
    relative_close(heat, 2.309123036683434e-55, 4e-12);
    relative_close(slope, 1.974957620386691e-104, 4e-12);
    let (vacuum_heat, vacuum_slope) =
        planck_radiative_exchange_w_and_w_k(cold, None, wavelength, 1e-100, 1.0).unwrap();
    relative_close(vacuum_heat, 4.083964333362426e210, 4e-12);
    relative_close(vacuum_slope, 4.083964333362426e145, 4e-12);
}

#[test]
fn equilibrium_vacuum_limits_and_subnormal_tangents_keep_distinct_meanings() {
    for (t, wavelength) in [(300.0, 10_000.0), (6000.0, 500.0), (1e-77, 3e83)] {
        assert_eq!(
            planck_radiance_difference_w_m2_sr_nm(t, t, wavelength),
            Ok(0.0)
        );
        let tangent = planck_radiance_temperature_derivative_w_m2_sr_nm_k(t, wavelength).unwrap();
        assert_eq!(
            planck_radiance_temperature_secant_w_m2_sr_nm_k(t, t, wavelength),
            Ok(tangent)
        );
        let equilibrium =
            planck_radiative_exchange_w_and_w_k(t, Some(t), wavelength, 1.0, 1.0).unwrap();
        assert_eq!(equilibrium, (0.0, tangent));
    }
    let t = 1e-77;
    let wavelength = 3e83;
    let tangent = planck_radiance_temperature_derivative_w_m2_sr_nm_k(t, wavelength).unwrap();
    relative_close(tangent, 1.9749576203866902e-322, 0.0);
    assert!(tangent > 0.0);
    let scaled = planck_radiative_exchange_w_and_w_k(t, Some(t), wavelength, 1e300, 1.0).unwrap();
    assert_eq!(scaled.0, 0.0);
    relative_close(scaled.1, 1.9749576203866902e-22, 4e-12);
    let (heat, slope) =
        planck_radiative_exchange_w_and_w_k(300.0, None, 10_000.0, 1.0, 1.0).unwrap();
    relative_close(heat, 9.924033330070695e-3, 3e-12);
    relative_close(slope, 9.924033330070695e-3 / 300.0, 3e-12);
    assert!(slope < planck_radiance_temperature_derivative_w_m2_sr_nm_k(300.0, 10_000.0).unwrap());
}

#[test]
fn invalid_domains_and_unrepresentable_final_rates_refuse_explicitly() {
    for t in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert_eq!(
            planck_radiance_difference_w_m2_sr_nm(t, t, 500.0),
            Err(ThermalSurfaceError::Temperature)
        );
        assert_eq!(
            planck_radiance_temperature_secant_w_m2_sr_nm_k(300.0, t, 500.0),
            Err(ThermalSurfaceError::Temperature)
        );
        assert_eq!(
            planck_radiance_temperature_derivative_w_m2_sr_nm_k(t, 500.0),
            Err(ThermalSurfaceError::Temperature)
        );
    }
    for wavelength in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert_eq!(
            planck_radiative_exchange_w_and_w_k(300.0, Some(300.0), wavelength, 0.0, 1.0),
            Err(ThermalSurfaceError::Wavelength)
        );
    }
    for invalid in [-1.0, f64::INFINITY, f64::NAN] {
        assert_eq!(
            planck_radiative_exchange_w_and_w_k(300.0, None, 500.0, invalid, 1.0),
            Err(ThermalSurfaceError::ExchangeCoefficient)
        );
        assert_eq!(
            planck_radiative_exchange_w_and_w_k(300.0, None, 500.0, 1.0, invalid),
            Err(ThermalSurfaceError::ExchangeCoefficient)
        );
    }
    assert_eq!(
        planck_radiative_exchange_w_and_w_k(1e65, None, 3e-59, 0.0, 1.0),
        Ok((0.0, 0.0))
    );
    assert_eq!(
        planck_radiative_exchange_w_and_w_k(1e65, None, 3e-59, 1.0, 0.0),
        Ok((0.0, 0.0))
    );
    assert_eq!(
        planck_radiative_exchange_w_and_w_k(1e65, None, 3e-59, 1.0, 1.0),
        Err(ThermalSurfaceError::Unrepresentable)
    );
    // The scaled heat is about 8.278e-288 W, but its secant is about
    // 8.278e-588 W/K. Returning g=0 would bypass the native cooling limiter.
    assert_eq!(
        planck_radiative_exchange_w_and_w_k(1e300, None, 1e300, 1e300, 1e300),
        Err(ThermalSurfaceError::Unrepresentable)
    );
    // Both final rates genuinely underflow without that large scale; no
    // nonzero heat is published and no artificial refusal is necessary.
    assert_eq!(
        planck_radiative_exchange_w_and_w_k(1e300, None, 1e300, 1.0, 1.0),
        Ok((0.0, 0.0))
    );
    // C2/lambda alone overflows, although x=C2/(lambda*T) is ordinary.
    // The new exchange path must refuse truly overflowing watts, not mistake
    // that intermediate infinity for a dark Wien tail and return zero heat.
    assert_eq!(
        planck_radiative_exchange_w_and_w_k(1e308, None, 5e-302, 1.0, 1.0),
        Err(ThermalSurfaceError::Unrepresentable)
    );
}
