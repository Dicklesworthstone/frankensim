use super::*;

fn state(p: f64, t: f64, region: If97Region) -> If97State {
    If97State::at_pressure_temperature(Pressure::new(p), Temperature::new(t), region).unwrap()
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:.17e}, expected={expected:.17e}, error={:.8e}, allowance={tolerance:.8e}",
        (actual - expected).abs()
    );
}

#[test]
fn original_region1_table5_in_coherent_si() {
    // Original release Table 5. Tolerances are half the final printed digit
    // after converting the published kJ values to coherent SI units.
    for (p, t, expected, tolerance) in [
        (
            3e6,
            300.0,
            [
                0.00100215168,
                115331.273,
                112324.818,
                392.294792,
                4173.01218,
                1507.73921,
            ],
            [5e-12, 5e-4, 5e-4, 5e-7, 5e-6, 5e-6],
        ),
        (
            80e6,
            300.0,
            [
                0.000971180894,
                184142.828,
                106448.356,
                368.563852,
                4010.08987,
                1634.69054,
            ],
            [5e-13, 5e-4, 5e-4, 5e-7, 5e-6, 5e-6],
        ),
        (
            3e6,
            500.0,
            [
                0.00120241800,
                975542.239,
                971934.985,
                2580.41912,
                4655.80682,
                1240.71337,
            ],
            [5e-12, 5e-4, 5e-4, 5e-6, 5e-6, 5e-6],
        ),
    ] {
        let s = state(p, t, If97Region::Region1);
        let actual = [
            s.specific_volume().value(),
            s.specific_enthalpy().value(),
            s.specific_internal_energy().value(),
            s.specific_entropy().value(),
            s.cp().value(),
            s.sound_speed().value(),
        ];
        for k in 0..actual.len() {
            close(actual[k], expected[k], tolerance[k]);
        }
    }
}

#[test]
fn original_region2_table15_in_coherent_si() {
    // Original release Table 15, including the dense admitted Region 2 point.
    for (p, t, expected, tolerance) in [
        (
            3500.0,
            300.0,
            [
                39.4913866, 2549911.45, 2411691.60, 8522.38967, 1913.00162, 427.920172,
            ],
            [5e-8, 5e-3, 5e-3, 5e-6, 5e-6, 5e-7],
        ),
        (
            3500.0,
            700.0,
            [
                92.3015898, 3335683.75, 3012628.19, 10174.9996, 2081.41274, 644.289068,
            ],
            [5e-8, 5e-3, 5e-3, 5e-5, 5e-6, 5e-7],
        ),
        (
            30e6,
            700.0,
            [
                0.00542946619,
                2631494.74,
                2468610.76,
                5175.40298,
                10350.5092,
                480.386523,
            ],
            [5e-12, 5e-3, 5e-3, 5e-6, 5e-5, 5e-7],
        ),
    ] {
        let s = state(p, t, If97Region::Region2);
        let actual = [
            s.specific_volume().value(),
            s.specific_enthalpy().value(),
            s.specific_internal_energy().value(),
            s.specific_entropy().value(),
            s.cp().value(),
            s.sound_speed().value(),
        ];
        for k in 0..actual.len() {
            close(actual[k], expected[k], tolerance[k]);
        }
    }
}

#[test]
fn original_region4_tables35_36_and_b23_reference() {
    for (t, p, allowance) in [
        (300.0, 3536.58941, 5e-6),
        (500.0, 2638897.76, 5e-3),
        (600.0, 12344314.6, 5e-2),
    ] {
        close(
            saturation_pressure(Temperature::new(t)).unwrap().value(),
            p,
            allowance,
        );
    }
    for (p, t) in [(1e5, 372.755919), (1e6, 453.035632), (1e7, 584.149488)] {
        close(
            saturation_temperature(Pressure::new(p)).unwrap().value(),
            t,
            5e-7,
        );
    }
    close(
        boundary23_pressure(Temperature::new(623.15))
            .unwrap()
            .value(),
        16.5291643e6,
        0.05,
    );
}

fn derivative(mut f: impl FnMut(f64) -> f64, x: f64, step: f64) -> f64 {
    (f(x - 2.0 * step) - 8.0 * f(x - step) + 8.0 * f(x + step) - f(x + 2.0 * step)) / (12.0 * step)
}

#[test]
fn analytic_derivatives_obey_independent_caloric_and_maxwell_differences() {
    for (p, t, region) in [
        (3e6, 300.0, If97Region::Region1),
        (80e6, 500.0, If97Region::Region1),
        (1e5, 500.0, If97Region::Region2),
        (30e6, 700.0, If97Region::Region2),
        (90e6, 1000.0, If97Region::Region2),
    ] {
        let s = state(p, t, region);
        let dh_dt = derivative(
            |temperature| state(p, temperature, region).specific_enthalpy().value(),
            t,
            0.04,
        );
        let ds_dt = derivative(
            |temperature| state(p, temperature, region).specific_entropy().value(),
            t,
            0.04,
        );
        let dv_dt = derivative(
            |temperature| state(p, temperature, region).specific_volume().value(),
            t,
            0.04,
        );
        let dp = p * 1e-4;
        let dh_dp = derivative(
            |pressure| state(pressure, t, region).specific_enthalpy().value(),
            p,
            dp,
        );
        let ds_dp = derivative(
            |pressure| state(pressure, t, region).specific_entropy().value(),
            p,
            dp,
        );
        let dv_dp = derivative(
            |pressure| state(pressure, t, region).specific_volume().value(),
            p,
            dp,
        );
        close(dh_dt, s.cp().value(), 3e-7 * s.cp().value());
        close(t * ds_dt, dh_dt, 3e-7 * s.cp().value());
        close(dv_dt, s.dv_dt_at_pressure().value(), 3e-7 * dv_dt.abs());
        close(dh_dp, s.dh_dp_at_temperature().value(), 3e-6 * dh_dp.abs());
        close(ds_dp, s.ds_dp_at_temperature().value(), 3e-6 * ds_dp.abs());
        close(dv_dp, s.dv_dp_at_temperature().value(), 3e-7 * dv_dp.abs());
        // Independent finite derivatives verify T ds = dh - v dp.
        close(
            t * ds_dp,
            dh_dp - s.specific_volume().value(),
            3e-6 * (t * ds_dp).abs(),
        );
        close(
            s.specific_enthalpy().value() - s.specific_internal_energy().value(),
            p * s.specific_volume().value(),
            2e-9,
        );
        close(
            s.specific_gibbs_energy().value(),
            s.specific_enthalpy().value() - t * s.specific_entropy().value(),
            3e-9,
        );
    }
}

fn at_volume_temperature(volume: f64, t: f64, initial_p: f64, region: If97Region) -> If97State {
    let mut p = initial_p;
    for _ in 0..12 {
        let s = state(p, t, region);
        let change = (s.specific_volume().value() - volume) / s.dv_dp_at_temperature().value();
        let next = p - change;
        if next == p {
            return s;
        }
        p = next;
    }
    let s = state(p, t, region);
    close(s.specific_volume().value(), volume, 2e-14 * volume);
    s
}

#[test]
fn cv_is_actual_constant_volume_internal_energy_derivative() {
    for (p, t, region) in [
        (3e6, 300.0, If97Region::Region1),
        (80e6, 500.0, If97Region::Region1),
        (1e5, 500.0, If97Region::Region2),
        (30e6, 700.0, If97Region::Region2),
    ] {
        let s = state(p, t, region);
        let v = s.specific_volume().value();
        let cv = derivative(
            |temperature| {
                at_volume_temperature(v, temperature, p, region)
                    .specific_internal_energy()
                    .value()
            },
            t,
            0.02,
        );
        close(cv, s.cv().value(), 3e-6 * s.cv().value());
        let alpha = s.thermal_expansion().value();
        let kappa = s.isothermal_compressibility().value();
        close(
            s.cp().value() - s.cv().value(),
            t * alpha * alpha * v / kappa,
            5e-10,
        );
    }
}

#[test]
fn density_anomaly_and_dilute_vapor_are_retained() {
    let cold = state(1e5, MINIMUM_TEMPERATURE_K, If97Region::Region1);
    let warm = state(1e5, 300.0, If97Region::Region1);
    assert!(cold.thermal_expansion().value() < 0.0);
    assert!(warm.thermal_expansion().value() > 0.0);
    let dilute = state(1e-100, 700.0, If97Region::Region2);
    close(
        dilute.cp().value() - dilute.cv().value(),
        R_SPECIFIC_J_KG_K,
        4e-13,
    );
    close(
        dilute.pressure().value() * dilute.specific_volume().value(),
        R_SPECIFIC_J_KG_K * dilute.temperature().value(),
        6e-11,
    );
    assert!(dilute.specific_entropy().value().is_finite());
    assert!(dilute.dv_dp_at_temperature().value().is_finite());
    assert_eq!(warm, state(1e5, 300.0, If97Region::Region1));
}

#[test]
fn region4_roundtrip_slope_and_closed_anchors() {
    for i in 0..257 {
        let t = MINIMUM_TEMPERATURE_K
            + (CRITICAL_TEMPERATURE_K - MINIMUM_TEMPERATURE_K) * f64::from(i) / 256.0;
        let p = saturation_pressure(Temperature::new(t)).unwrap();
        close(saturation_temperature(p).unwrap().value(), t, 4e-10);
    }
    for t in [
        MINIMUM_TEMPERATURE_K,
        REGION1_MAXIMUM_TEMPERATURE_K,
        CRITICAL_TEMPERATURE_K,
    ] {
        let p = saturation_pressure(Temperature::new(t)).unwrap();
        assert_eq!(
            saturation_temperature(p).unwrap().value().to_bits(),
            t.to_bits()
        );
    }
    for t in [300.0, 500.0, 600.0, 646.0] {
        let numerical = derivative(
            |temperature| {
                saturation_pressure(Temperature::new(temperature))
                    .unwrap()
                    .value()
            },
            t,
            0.002,
        );
        let analytic = saturation_pressure_derivative(Temperature::new(t))
            .unwrap()
            .value();
        close(analytic, numerical, 3e-7 * numerical);
    }
}

#[test]
fn explicit_saturation_branches_retain_original_finite_gibbs_seam() {
    let mut observed_nonzero_gap = false;
    for t in [
        MINIMUM_TEMPERATURE_K,
        300.0,
        400.0,
        500.0,
        600.0,
        REGION1_MAXIMUM_TEMPERATURE_K,
    ] {
        let pair = If97SaturationPair::at_temperature(Temperature::new(t)).unwrap();
        assert_eq!(pair.temperature().value().to_bits(), t.to_bits());
        assert_eq!(pair.liquid().phase(), If97Phase::SaturatedLiquid);
        assert_eq!(pair.vapor().phase(), If97Phase::SaturatedVapor);
        assert_eq!(
            *pair.liquid(),
            state(pair.pressure().value(), t, If97Region::Region1)
        );
        assert_eq!(
            *pair.vapor(),
            state(pair.pressure().value(), t, If97Region::Region2)
        );
        assert!(pair.latent_enthalpy().value() > 0.0);
        assert!(pair.latent_entropy().value() > 0.0);
        assert!(pair.liquid().density().value() > pair.vapor().density().value());
        let gap = pair.gibbs_energy_gap().value();
        close(
            gap,
            pair.latent_enthalpy().value() - t * pair.latent_entropy().value(),
            3e-9,
        );
        observed_nonzero_gap |= gap.abs() > 0.1;
        let inverse_pair = If97SaturationPair::at_pressure(pair.pressure()).unwrap();
        assert_eq!(
            inverse_pair.pressure().value().to_bits(),
            pair.pressure().value().to_bits()
        );
        assert_eq!(
            *inverse_pair.liquid(),
            state(
                inverse_pair.pressure().value(),
                inverse_pair.temperature().value(),
                If97Region::Region1
            )
        );
        assert_eq!(
            *inverse_pair.vapor(),
            state(
                inverse_pair.pressure().value(),
                inverse_pair.temperature().value(),
                If97Region::Region2
            )
        );
    }
    assert!(
        observed_nonzero_gap,
        "raw source coexistence mismatch must not be forced to zero"
    );
    let pair = If97SaturationPair::at_pressure(Pressure::new(101325.0)).unwrap();
    assert_eq!(pair.pressure().value(), 101325.0);
    // This source's engineering lower bound differs from parent SR1-86.
    assert!(super::super::SaturatedWaterState::at_temperature(Temperature::new(273.15)).is_err());
}

#[test]
fn full_implemented_domains_and_unimplemented_regions_are_explicit() {
    for (p, t, region) in [
        (100e6, MINIMUM_TEMPERATURE_K, If97Region::Region1),
        (100e6, REGION1_MAXIMUM_TEMPERATURE_K, If97Region::Region1),
        (100e6, 863.15, If97Region::Region2),
        (100e6, REGION2_MAXIMUM_TEMPERATURE_K, If97Region::Region2),
    ] {
        state(p, t, region);
    }
    for (p, t, region) in [
        (3e6, MINIMUM_TEMPERATURE_K.next_down(), If97Region::Region1),
        (
            3e6,
            REGION1_MAXIMUM_TEMPERATURE_K.next_up(),
            If97Region::Region1,
        ),
        (100e6_f64.next_up(), 300.0, If97Region::Region1),
        (
            100e6,
            REGION2_MAXIMUM_TEMPERATURE_K.next_up(),
            If97Region::Region2,
        ),
        (25e6, 650.0, If97Region::Region1),
        (25e6, 650.0, If97Region::Region2),
        (1e5, 1200.0, If97Region::Region2),
        (0.0, 300.0, If97Region::Region2),
        (f64::NAN, 300.0, If97Region::Region2),
        (1e5, f64::INFINITY, If97Region::Region2),
    ] {
        assert!(
            If97State::at_pressure_temperature(Pressure::new(p), Temperature::new(t), region)
                .is_err()
        );
    }
    let p = saturation_pressure(Temperature::new(500.0))
        .unwrap()
        .value();
    assert!(
        If97State::at_pressure_temperature(
            Pressure::new(p * (1.0 - 1e-8)),
            Temperature::new(500.0),
            If97Region::Region1
        )
        .is_err()
    );
    assert!(
        If97State::at_pressure_temperature(
            Pressure::new(p * (1.0 + 1e-8)),
            Temperature::new(500.0),
            If97Region::Region2
        )
        .is_err()
    );
    for t in [623.15_f64.next_up(), CRITICAL_TEMPERATURE_K] {
        assert!(matches!(
            If97SaturationPair::at_temperature(Temperature::new(t)),
            Err(If97Error::UnsupportedSaturation { .. })
        ));
    }
    let minimum = saturation_pressure(Temperature::new(MINIMUM_TEMPERATURE_K))
        .unwrap()
        .value();
    let maximum = saturation_pressure(Temperature::new(CRITICAL_TEMPERATURE_K))
        .unwrap()
        .value();
    let pair_maximum = saturation_pressure(Temperature::new(REGION1_MAXIMUM_TEMPERATURE_K))
        .unwrap()
        .value();
    assert!(If97SaturationPair::at_pressure(Pressure::new(pair_maximum)).is_ok());
    assert!(If97SaturationPair::at_pressure(Pressure::new(pair_maximum.next_up())).is_err());
    assert!(saturation_temperature(Pressure::new(minimum.next_down())).is_err());
    assert!(saturation_temperature(Pressure::new(maximum.next_up())).is_err());
    assert!(boundary23_pressure(Temperature::new(623.15_f64.next_down())).is_err());
    assert!(boundary23_pressure(Temperature::new(863.15_f64.next_up())).is_err());
}
