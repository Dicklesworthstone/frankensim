use super::*;

fn state(p: f64, t: f64) -> IceIhState {
    IceIhState::at_pressure_temperature(Pressure::new(p), Temperature::new(t)).unwrap()
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        a.is_finite() && b.is_finite() && (a - b).abs() <= tolerance,
        "actual {a:.17e}, expected {b:.17e}, error {:.4e}, tolerance {tolerance:.4e}",
        a - b
    );
}
fn solve() -> IceMeltingSolve {
    IceMeltingSolve {
        pressure_tolerance_pa: 1e-5,
        temperature_tolerance_k: 1e-10,
        maximum_iterations: 96,
    }
}

#[test]
fn every_original_ice_verification_property_retains_the_iapws95_reference() {
    // IAPWS R10-06(2009), Table 6, all three columns and all sixteen rows.
    let columns = [
        (
            273.16,
            611.657,
            [
                0.611_784_135,
                0.109_085_812_737e-2,
                0.122_069_433_940e4,
                -0.128_495_941_571e-12,
                0.174_387_964_700e-6,
                -0.767_602_985_875e1,
                -0.333_444_253_966e6,
                -0.554_468_75e-1,
                -0.333_444_921_197e6,
                -0.122_069_433_940e4,
                0.209_678_431_622e4,
                0.916_709_492_200e3,
                0.159_863_102_566e-3,
                0.135_714_764_659e7,
                0.117_793_449_348e-9,
                0.114_161_597_779e-9,
            ],
        ),
        (
            273.152_519,
            101_325.0,
            [
                0.101_342_740_69e3,
                0.109_084_388_214e-2,
                0.122_076_932_550e4,
                -0.128_485_364_928e-12,
                0.174_362_219_972e-6,
                -0.767_598_233_365e1,
                -0.333_354_873_637e6,
                -0.918_701_567e1,
                -0.333_465_403_393e6,
                -0.122_076_932_550e4,
                0.209_671_391_024e4,
                0.916_721_463_419e3,
                0.159_841_589_458e-3,
                0.135_705_899_321e7,
                0.117_785_291_765e-9,
                0.114_154_442_556e-9,
            ],
        ),
        (
            100.0,
            100e6,
            [
                -0.222_296_513_088e6,
                0.106_193_389_260e-2,
                0.261_195_122_589e4,
                -0.941_807_981_761e-13,
                0.274_505_162_488e-7,
                -0.866_333_195_517e1,
                -0.483_491_635_676e6,
                -0.328_489_902_347e6,
                -0.589_685_024_936e6,
                -0.261_195_122_589e4,
                0.866_333_195_517e3,
                0.941_678_203_297e3,
                0.258_495_528_207e-4,
                0.291_466_166_994e6,
                0.886_880_048_115e-10,
                0.886_060_982_687e-10,
            ],
        ),
    ];
    for (t, p, expected) in columns {
        let s = state(p, t);
        let actual = [
            s.specific_gibbs_energy().value(),
            s.specific_volume().value(),
            -s.specific_entropy().value(),
            s.volume_pressure_derivative().value(),
            s.volume_temperature_derivative().value(),
            s.gibbs_temperature_second_derivative().value(),
            s.specific_enthalpy().value(),
            s.specific_helmholtz_energy().value(),
            s.specific_internal_energy().value(),
            s.specific_entropy().value(),
            s.cp().value(),
            s.density().value(),
            s.thermal_expansion().value(),
            s.pressure_coefficient().value(),
            s.isothermal_compressibility().value(),
            s.isentropic_compressibility().value(),
        ];
        for (i, (a, b)) in actual.into_iter().zip(expected).enumerate() {
            // Table rounding plus evaluation rounding, not physical uncertainty.
            let absolute = if i == 0 || i == 7 { 6e-9 } else { 0.0 };
            near(a, b, absolute + 6e-12 * b.abs());
        }
        assert_eq!(s.pressure().value(), p);
        assert_eq!(s.temperature().value(), t);
        assert!(s.cv().value() > 0.0 && s.cv().value() < s.cp().value());
    }
}

#[test]
fn analytic_ice_derivatives_close_energy_and_match_independent_differences() {
    for (p, t) in [
        (0.0, 50.0),
        (101_325.0, 253.15),
        (100e6, 100.0),
        (210e6, 260.0),
    ] {
        let s = state(p, t);
        near(
            s.specific_enthalpy().value() - s.specific_internal_energy().value(),
            p * s.specific_volume().value(),
            1e-10,
        );
        near(
            s.specific_gibbs_energy().value() + t * s.specific_entropy().value(),
            s.specific_enthalpy().value(),
            1e-10,
        );
        near(
            s.isentropic_compressibility().value(),
            s.isothermal_compressibility().value() * s.cv().value() / s.cp().value(),
            1e-24,
        );
    }
    for (p, t) in [(101_325.0, 250.0), (10e6, 200.0), (100e6, 100.0)] {
        let s = state(p, t);
        let (a, b) = (state(p, t - 0.01), state(p, t + 0.01));
        near(
            (b.specific_enthalpy().value() - a.specific_enthalpy().value()) / 0.02,
            s.cp().value(),
            2e-6,
        );
        near(
            (b.specific_entropy().value() - a.specific_entropy().value()) / 0.02,
            s.cp().value() / t,
            3e-8,
        );
        near(
            (b.specific_volume().value() - a.specific_volume().value()) / 0.02,
            s.volume_temperature_derivative().value(),
            3e-16,
        );
        let (a, b) = (state(p - 100.0, t), state(p + 100.0, t));
        near(
            (b.specific_gibbs_energy().value() - a.specific_gibbs_energy().value()) / 200.0,
            s.specific_volume().value(),
            3e-12,
        );
        near(
            (b.specific_entropy().value() - a.specific_entropy().value()) / 200.0,
            -s.volume_temperature_derivative().value(),
            3e-14,
        );
        near(
            (b.specific_enthalpy().value() - a.specific_enthalpy().value()) / 200.0,
            s.specific_volume().value() - t * s.volume_temperature_derivative().value(),
            3e-12,
        );
    }
    for (p, t) in [
        (f64::NAN, 250.0),
        (0.0, 49.9),
        (-1.0, 250.0),
        (210e6 + 1.0, 250.0),
        (0.0, 273.1601),
    ] {
        assert!(
            IceIhState::at_pressure_temperature(Pressure::new(p), Temperature::new(t)).is_err()
        );
    }
}

#[test]
fn stable_ice_differences_preserve_tiny_heating_reversal_and_integrals() {
    for p in [0.0, 101_325.0, 100e6, 210e6] {
        for (a, b) in [
            (50.0, 273.16),
            (230.0, 273.15),
            (260.0, 260.0_f64.next_up()),
        ] {
            let delta =
                isobaric_difference(Pressure::new(p), Temperature::new(a), Temperature::new(b))
                    .unwrap();
            let reverse =
                isobaric_difference(Pressure::new(p), Temperature::new(b), Temperature::new(a))
                    .unwrap();
            assert_eq!(
                delta.specific_enthalpy_change().value(),
                -reverse.specific_enthalpy_change().value()
            );
            assert_eq!(
                delta.specific_entropy_change().value(),
                -reverse.specific_entropy_change().value()
            );
            assert!(delta.specific_enthalpy_change().value() > 0.0);
            near(
                delta.specific_internal_energy_change().value()
                    + p * delta.specific_volume_change().value(),
                delta.specific_enthalpy_change().value(),
                1e-10,
            );
            if b - a < 1e-10 {
                near(
                    delta.specific_enthalpy_change().value() / (b - a),
                    state(p, a).cp().value(),
                    3e-10,
                );
                near(
                    delta.specific_entropy_change().value() / (b - a),
                    state(p, a).cp().value() / a,
                    3e-12,
                );
            } else {
                near(
                    delta.specific_enthalpy_change().value(),
                    state(p, b).specific_enthalpy().value()
                        - state(p, a).specific_enthalpy().value(),
                    3e-9,
                );
                near(
                    delta.specific_entropy_change().value(),
                    state(p, b).specific_entropy().value() - state(p, a).specific_entropy().value(),
                    2e-11,
                );
                let middle = (a + b) * 0.5;
                let one = isobaric_difference(
                    Pressure::new(p),
                    Temperature::new(a),
                    Temperature::new(middle),
                )
                .unwrap();
                let two = isobaric_difference(
                    Pressure::new(p),
                    Temperature::new(middle),
                    Temperature::new(b),
                )
                .unwrap();
                near(
                    one.specific_enthalpy_change().value() + two.specific_enthalpy_change().value(),
                    delta.specific_enthalpy_change().value(),
                    3e-9,
                );
            }
        }
        let (a, b, n) = (233.15, 273.15, 512);
        let dx = (b - a) / (n as f64);
        let mut h = 0.0;
        let mut entropy = 0.0;
        for i in 0..=n {
            let t = a + (i as f64) * dx;
            let w = if i == 0 || i == n {
                1.0
            } else if i % 2 == 0 {
                2.0
            } else {
                4.0
            };
            h += w * state(p, t).cp().value();
            entropy += w * state(p, t).cp().value() / t;
        }
        let delta = isobaric_difference(Pressure::new(p), Temperature::new(a), Temperature::new(b))
            .unwrap();
        near(delta.specific_enthalpy_change().value(), h * dx / 3.0, 2e-7);
        near(
            delta.specific_entropy_change().value(),
            entropy * dx / 3.0,
            2e-9,
        );
        let zero = isobaric_difference(Pressure::new(p), Temperature::new(a), Temperature::new(a))
            .unwrap();
        assert_eq!(zero.specific_enthalpy_change().value(), 0.0);
        assert_eq!(zero.specific_internal_energy_change().value(), 0.0);
    }
}

#[test]
fn original_melting_curve_round_trips_and_preserves_source_and_overlap_anchors() {
    // R14-08(2011), Table 3; tabulated pressure precision is 0.001 MPa.
    near(
        melting_pressure(Temperature::new(260.0)).unwrap().value(),
        138.268e6,
        500.0,
    );
    let mut previous = f64::INFINITY;
    for t in [251.165, 260.0, 270.0, 273.15, 273.16] {
        let p = melting_pressure(Temperature::new(t)).unwrap();
        assert!(p.value() < previous);
        previous = p.value();
        let result = melting_temperature(p, solve(), || false).unwrap();
        near(result.temperature().value(), t, 1e-10);
        assert!(result.pressure_residual_pa().abs() <= 1e-5);
        assert_eq!(result.requested_pressure(), p);
        if [251.165, 273.15, 273.16].contains(&t) {
            assert_eq!(result.temperature().value(), t);
            assert_eq!(result.temperature_bounds_k(), [t, t]);
            assert_eq!(result.iterations(), 0);
        }
    }
    let ambient = melting_temperature(Pressure::new(101_325.0), solve(), || false).unwrap();
    // The compact correlation is distinct from solving the two full Gibbs EOS.
    near(ambient.temperature().value(), 273.152_519, 3e-6);
    let t = 270.0;
    let (ta, tb) = (t - 0.0001, t + 0.0001);
    let dp = (melting_pressure(Temperature::new(tb)).unwrap().value()
        - melting_pressure(Temperature::new(ta)).unwrap().value())
        / (tb - ta);
    near(
        dp,
        melting_pressure_derivative(Temperature::new(t))
            .unwrap()
            .value(),
        0.001,
    );
    assert!(melting_pressure(Temperature::new(251.164)).is_err());
    assert!(melting_temperature(Pressure::new(611.656), solve(), || false).is_err());
    assert!(melting_temperature(Pressure::new(209e6), solve(), || false).is_err());
    let mut short = solve();
    short.maximum_iterations = 1;
    assert!(matches!(
        melting_temperature(Pressure::new(101_325.0), short, || false),
        Err(IceIhError::NotConverged { .. })
    ));
    assert!(matches!(
        melting_temperature(Pressure::new(101_325.0), solve(), || true),
        Err(IceIhError::Cancelled)
    ));
    let mut count = 0;
    assert!(matches!(
        melting_temperature(Pressure::new(611.657), solve(), || {
            count += 1;
            count > 1
        }),
        Err(IceIhError::Cancelled)
    ));
}
