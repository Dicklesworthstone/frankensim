use super::*;

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 8.0 * f64::EPSILON * expected.abs().max(f64::MIN_POSITIVE),
        "actual {actual:e}, expected {expected:e}",
    );
}

#[test]
fn g1_exact_si_photon_and_current_reference() {
    // Work function is an explicit synthetic example; this yield is not a
    // measured copper property. References are independently evaluated from SI.
    let law = PhotoelectricLaw::try_new(4.31 * ELEMENTARY_CHARGE_C, 1.0e-4, 0.5).unwrap();
    let response = law.evaluate(250.0e-9, 1.0).unwrap();
    close(response.photon_energy_j, 7.945_783_428_595_714e-19);
    close(
        response.maximum_kinetic_energy_j / ELEMENTARY_CHARGE_C,
        0.649_367_937_328_010_6,
    );
    close(response.stopping_potential_v, 0.649_367_937_328_010_6);
    close(response.absorbed_photon_rate_s, 1.258_529_141_885_677_3e18);
    close(
        response.surface_charging_current_a,
        2.016_385_984_337_303e-5,
    );
    close(
        response.mean_kinetic_energy_j,
        response.maximum_kinetic_energy_j * 0.5,
    );
    assert_eq!(law.work_function_j(), 4.31 * ELEMENTARY_CHARGE_C);
    assert_eq!(law.absorbed_quantum_yield(), 1.0e-4);
    assert_eq!(law.mean_excess_fraction(), 0.5);
}

#[test]
fn g0_threshold_is_strict_and_independent_of_intensity() {
    let threshold_wavelength = 250.0e-9;
    let work = PLANCK_J_S * LIGHT_M_S / threshold_wavelength;
    let law = PhotoelectricLaw::try_new(work, 1.0, 1.0).unwrap();
    for wavelength in [threshold_wavelength, 500.0e-9] {
        for power in [0.0, 1.0, 1.0e100] {
            let response = law.evaluate(wavelength, power).unwrap();
            assert_eq!(response.emitted_electron_rate_s, 0.0);
            assert_eq!(response.maximum_kinetic_energy_j, 0.0);
            assert_eq!(response.stopping_potential_v, 0.0);
            assert_eq!(response.removal_power_w, 0.0);
            assert_eq!(response.kinetic_power_w, 0.0);
            assert_eq!(response.heat_power_w, power);
        }
    }
    let above = law.evaluate(125.0e-9, 1.0).unwrap();
    close(above.maximum_kinetic_energy_j, work);
    assert!(above.emitted_electron_rate_s > 0.0);
}

#[test]
fn g3_power_scales_rates_but_not_electron_energy() {
    let law = PhotoelectricLaw::try_new(2.0 * ELEMENTARY_CHARGE_C, 0.3, 0.4).unwrap();
    let a = law.evaluate(300.0e-9, 0.25).unwrap();
    let b = law.evaluate(300.0e-9, 2.0).unwrap();
    assert_eq!(a.maximum_kinetic_energy_j, b.maximum_kinetic_energy_j);
    assert_eq!(a.mean_kinetic_energy_j, b.mean_kinetic_energy_j);
    assert_eq!(a.stopping_potential_v, b.stopping_potential_v);
    for (first, second) in [
        (a.absorbed_photon_rate_s, b.absorbed_photon_rate_s),
        (a.emitted_electron_rate_s, b.emitted_electron_rate_s),
        (a.surface_charging_current_a, b.surface_charging_current_a),
        (a.kinetic_power_w, b.kinetic_power_w),
        (a.removal_power_w, b.removal_power_w),
        (a.heat_power_w, b.heat_power_w),
    ] {
        close(second, 8.0 * first);
    }
}

#[test]
fn g0_charge_and_energy_balance_at_all_fraction_endpoints() {
    for yield_per_absorbed in [0.0, 0.2, 0.999_999_999_999_999_9, 1.0] {
        for mean_fraction in [0.0, 0.37, 0.999_999_999_999_999_9, 1.0] {
            let law = PhotoelectricLaw::try_new(
                3.0 * ELEMENTARY_CHARGE_C,
                yield_per_absorbed,
                mean_fraction,
            )
            .unwrap();
            let response = law.evaluate(200.0e-9, 2.5).unwrap();
            assert_eq!(
                response.surface_charging_current_a + response.emitted_electron_current_a,
                0.0,
            );
            assert!(response.emitted_electron_rate_s <= response.absorbed_photon_rate_s);
            assert!(response.mean_kinetic_energy_j <= response.maximum_kinetic_energy_j);
            close(
                response.heat_power_w + response.removal_power_w + response.kinetic_power_w,
                response.absorbed_power_w,
            );
            close(
                response.removal_power_w,
                response.emitted_electron_rate_s * law.work_function_j(),
            );
            close(
                response.kinetic_power_w,
                response.emitted_electron_rate_s * response.mean_kinetic_energy_j,
            );
            assert!(response.heat_power_w >= 0.0);
            if yield_per_absorbed == 1.0 && mean_fraction == 1.0 {
                assert_eq!(response.heat_power_w, 0.0);
            } else {
                // Small positive heat survives even next to unit fractions.
                assert!(response.heat_power_w > 0.0);
            }
        }
    }
}

#[test]
fn g0_dark_and_invalid_or_unrepresentable_inputs() {
    let law = PhotoelectricLaw::try_new(4.0 * ELEMENTARY_CHARGE_C, 0.25, 0.5).unwrap();
    let dark = law.evaluate(200.0e-9, 0.0).unwrap();
    assert_eq!(dark.absorbed_photon_rate_s, 0.0);
    assert_eq!(dark.emitted_electron_rate_s, 0.0);
    assert_eq!(dark.heat_power_w, 0.0);
    assert_eq!(dark.removal_power_w, 0.0);
    assert_eq!(dark.kinetic_power_w, 0.0);
    assert!(dark.maximum_kinetic_energy_j > 0.0);
    for work in [-1.0, 0.0, f64::NAN, f64::INFINITY] {
        assert!(PhotoelectricLaw::try_new(work, 0.2, 0.5).is_err());
    }
    for fraction in [-0.1, 1.1, f64::NAN, f64::INFINITY] {
        assert!(PhotoelectricLaw::try_new(1.0e-19, fraction, 0.5).is_err());
        assert!(PhotoelectricLaw::try_new(1.0e-19, 0.2, fraction).is_err());
    }
    for wavelength in [-1.0, 0.0, f64::NAN, f64::INFINITY, f64::MAX] {
        assert!(law.evaluate(wavelength, 1.0).is_err());
    }
    for power in [-1.0, f64::NAN, f64::INFINITY, f64::MAX] {
        assert!(law.evaluate(200.0e-9, power).is_err());
    }
    // Positive charge and heat must not silently disappear through underflow.
    assert!(law.evaluate(200.0e-9, f64::from_bits(1)).is_err());
    assert!(law.evaluate(f64::from_bits(1), 1.0).is_err());
}

#[test]
fn g3_spectral_power_splitting_preserves_charge_and_heat() {
    let law = PhotoelectricLaw::try_new(4.0 * ELEMENTARY_CHARGE_C, 0.2, 0.5).unwrap();
    let uv = law.evaluate(200.0e-9, 1.0).unwrap();
    let uv_a = law.evaluate(200.0e-9, 0.25).unwrap();
    let uv_b = law.evaluate(200.0e-9, 0.75).unwrap();
    let visible = law.evaluate(600.0e-9, 100.0).unwrap();
    close(
        uv_a.emitted_electron_rate_s + uv_b.emitted_electron_rate_s,
        uv.emitted_electron_rate_s,
    );
    close(
        uv_a.kinetic_power_w + uv_b.kinetic_power_w,
        uv.kinetic_power_w,
    );
    close(
        uv_a.removal_power_w + uv_b.removal_power_w,
        uv.removal_power_w,
    );
    close(uv_a.heat_power_w + uv_b.heat_power_w, uv.heat_power_w);
    assert_eq!(visible.emitted_electron_rate_s, 0.0);
    assert_eq!(
        uv.emitted_electron_rate_s + visible.emitted_electron_rate_s,
        uv.emitted_electron_rate_s,
    );
    close(
        uv.heat_power_w + visible.heat_power_w + uv.kinetic_power_w + uv.removal_power_w,
        101.0,
    );
}
