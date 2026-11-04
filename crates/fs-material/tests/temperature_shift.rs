//! G0/G1/G3 temperature shifts of the actual tensor spectrum and history owner.
//! Coefficients are authored numerical fixtures, not calibrated material data.

use fs_material::plane_stress_viscosity::MaterialPlane;
use fs_material::visco::temperature::{TemperatureShift, TemperatureShiftModel};
use fs_material::visco::{PronyBranch, TensorProny, ViscoError};
use fs_rand::StreamKey;

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1.0),
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
    );
}

fn shift() -> TemperatureShift {
    TemperatureShift::new(
        TemperatureShiftModel::WilliamsLandelFerry {
            c1: 4.0,
            c2_k: 80.0,
        },
        300.0,
        [270.0, 380.0],
    )
    .unwrap()
}

fn tensor(with_kelvin: bool) -> TensorProny {
    let equilibrium = core::array::from_fn(|i| {
        core::array::from_fn(|j| if i == j { 1e8 * (i + 1) as f64 } else { 0.0 })
    });
    let rows = [
        [2e3, -1e3, 1e3, 0.0, 5e2, -2e2],
        [0.0, 3e3, 2e3, -1e3, 0.0, 0.0],
        [0.0; 6],
        [0.0; 6],
        [0.0; 6],
        [0.0; 6],
    ];
    let branches = vec![
        PronyBranch::from_gram_rows(rows, 0.01).unwrap(),
        PronyBranch::from_gram_rows(rows.map(|row| row.map(|x| 2.0 * x)), 0.5).unwrap(),
    ];
    let law = TensorProny::new_mandel(equilibrium, branches, 0.02)
        .unwrap()
        .with_band((1.0, 1e5))
        .unwrap();
    if with_kelvin {
        let v = [30.0, -20.0, 10.0, 0.0, 5.0, -1.0];
        let kv = core::array::from_fn(|i| core::array::from_fn(|j| v[i] * v[j]));
        law.with_kelvin_viscosity_mandel(kv).unwrap()
    } else {
        law
    }
}

#[test]
fn g0_wlf_and_arrhenius_are_monotone_over_declared_calibrations() {
    let mut rng = StreamKey {
        seed: 0x574c_4641_5252,
        kernel: 0,
        tile: 0,
    }
    .stream();
    for _ in 0..4096 {
        let reference = 290.0 + 100.0 * rng.next_f64();
        let band = [reference - 10.0, reference + 25.0];
        let models = [
            TemperatureShiftModel::WilliamsLandelFerry {
                c1: 0.5 + 5.0 * rng.next_f64(),
                c2_k: 30.0 + 100.0 * rng.next_f64(),
            },
            TemperatureShiftModel::Arrhenius {
                activation_energy_j_mol: 1e3 + 8e4 * rng.next_f64(),
            },
        ];
        for model in models {
            let declaration = TemperatureShift::new(model, reference, band).unwrap();
            let cold = declaration.evaluate(band[0]).unwrap();
            let at_reference = declaration.evaluate(reference).unwrap();
            let hot = declaration.evaluate(band[1]).unwrap();
            assert!(cold.factor() > 1.0 && hot.factor() < 1.0 && hot.factor() > 0.0);
            assert_eq!(at_reference.factor().to_bits(), 1.0f64.to_bits());
            assert_eq!(at_reference.log_factor().to_bits(), 0.0f64.to_bits());
            assert!(cold.log_derivative_per_k() < 0.0 && hot.log_derivative_per_k() < 0.0);
            assert_eq!(*hot.declaration(), declaration);
            assert_eq!(declaration.reference_temperature_k(), reference);
            assert_eq!(declaration.validity_k(), band);
            assert_eq!(declaration.model(), model);
        }
    }
}

#[test]
fn g0_temperature_support_and_wlf_pole_refuse_without_extrapolation() {
    let wlf = TemperatureShiftModel::WilliamsLandelFerry {
        c1: 4.0,
        c2_k: 80.0,
    };
    for band in [
        [220.0, 350.0],
        [219.0, 350.0],
        [0.0, 350.0],
        [301.0, 350.0],
        [300.0, f64::INFINITY],
        [350.0, 280.0],
    ] {
        assert!(TemperatureShift::new(wlf, 300.0, band).is_err());
    }
    for model in [
        TemperatureShiftModel::WilliamsLandelFerry {
            c1: 0.0,
            c2_k: 80.0,
        },
        TemperatureShiftModel::WilliamsLandelFerry {
            c1: 4.0,
            c2_k: -1.0,
        },
        TemperatureShiftModel::Arrhenius {
            activation_energy_j_mol: 0.0,
        },
        TemperatureShiftModel::Arrhenius {
            activation_energy_j_mol: f64::NAN,
        },
    ] {
        assert!(TemperatureShift::new(model, 300.0, [280.0, 320.0]).is_err());
    }
    let declaration = shift();
    for t in [
        0.0,
        -1.0,
        f64::NAN,
        f64::INFINITY,
        270.0f64.next_down(),
        380.0f64.next_up(),
    ] {
        assert!(matches!(
            declaration.evaluate(t),
            Err(ViscoError::Parameters { .. })
        ));
    }
    // A broad, unjustified extrapolation cannot create infinite/zero kinetics.
    assert!(
        TemperatureShift::new(
            TemperatureShiftModel::Arrhenius {
                activation_energy_j_mol: 5e5
            },
            300.0,
            [1.0, 600.0],
        )
        .is_err()
    );
}

#[test]
fn g1_wlf_and_arrhenius_apparent_activation_energy_identity() {
    let wlf = shift();
    let arrhenius = TemperatureShift::new(
        TemperatureShiftModel::Arrhenius {
            activation_energy_j_mol: 50_000.0,
        },
        300.0,
        [270.0, 380.0],
    )
    .unwrap();
    let r = 8.314_462_618_153_24;
    for temperature in [275.0, 300.0, 370.0] {
        let inverse = 1.0 / temperature;
        let h = 1e-8;
        for declaration in [wlf, arrhenius] {
            let derivative = (declaration
                .evaluate(1.0 / (inverse + h))
                .unwrap()
                .log_factor()
                - declaration
                    .evaluate(1.0 / (inverse - h))
                    .unwrap()
                    .log_factor())
                / (2.0 * h);
            let expected = match declaration.model() {
                TemperatureShiftModel::WilliamsLandelFerry { c1, c2_k } => {
                    let ratio = temperature / (c2_k + temperature - 300.0);
                    core::f64::consts::LN_10 * c1 * c2_k * ratio * ratio
                }
                TemperatureShiftModel::Arrhenius {
                    activation_energy_j_mol,
                } => activation_energy_j_mol / r,
            };
            close(derivative, expected, 1e-8);
            close(
                -temperature
                    * temperature
                    * declaration
                        .evaluate(temperature)
                        .unwrap()
                        .log_derivative_per_k(),
                expected,
                5e-15,
            );
        }
    }
}

#[test]
fn g1_shifted_tensor_matches_hand_rescaled_spectrum_and_plane_stress() {
    let reference = tensor(true);
    assert_eq!(
        shift().shift_tensor(&reference, 300.0).unwrap().law(),
        &reference
    );
    for temperature in [280.0, 330.0] {
        let shifted = shift().shift_tensor(&reference, temperature).unwrap();
        let law = shifted.law();
        let a = shifted.state().factor();
        assert_eq!(law.equilibrium_mandel(), reference.equilibrium_mandel());
        assert_eq!(law.factors(), reference.factors());
        assert_eq!(law.maximum_strain_norm(), reference.maximum_strain_norm());
        let hand = TensorProny::new(
            *reference.c_eq(),
            reference.kv_pa_s().map(|c| c.map(|row| row.map(|x| x * a))),
            reference
                .branches()
                .iter()
                .map(|b| PronyBranch::new(b.c_pa, b.tau_s * a).unwrap())
                .collect(),
            (1.0 / a, 1e5 / a),
        )
        .unwrap();
        for reference_omega in [2.0, 10.0, 100.0, 1000.0, 50_000.0] {
            let omega = reference_omega / a;
            // Direct coefficient rescaling is bit-identical in the harmonic owner.
            assert_eq!(
                law.complex_modulus(omega).unwrap(),
                hand.complex_modulus(omega).unwrap()
            );
            let (actual_storage, actual_loss) = law.complex_modulus(omega).unwrap();
            let (expected_storage, expected_loss) =
                reference.complex_modulus(reference_omega).unwrap();
            for (actual, expected) in actual_storage
                .iter()
                .flatten()
                .chain(actual_loss.iter().flatten())
                .zip(
                    expected_storage
                        .iter()
                        .flatten()
                        .chain(expected_loss.iter().flatten()),
                )
            {
                close(*actual, *expected, 3e-14);
            }
            for plane in [
                MaterialPlane::Axes12,
                MaterialPlane::Axes13,
                MaterialPlane::Axes23,
            ] {
                let (ar, ai) = law.plane_stress_complex(plane, omega).unwrap();
                let (er, ei) = reference
                    .plane_stress_complex(plane, reference_omega)
                    .unwrap();
                for (actual, expected) in ar.into_iter().chain(ai).zip(er.into_iter().chain(ei)) {
                    close(actual, expected, 3e-13);
                }
            }
        }
    }
}

#[test]
fn g1_shifted_time_histories_keep_physical_dissipation_and_strain_limits() {
    let reference = tensor(true);
    for temperature in [280.0, 330.0] {
        let shifted = shift().shift_tensor(&reference, temperature).unwrap();
        let a = shifted.state().factor();
        let mut actual = shifted.law().history();
        let mut expected = reference.history();
        let mut mandel = shifted.law().state();
        let mut mandel_reference = reference.state();
        for (strain, dt) in [
            ([0.001, -0.002, 0.0, 0.003, 0.0, -0.001], 0.01),
            ([0.001, -0.002, 0.0, 0.003, 0.0, -0.001], 0.03),
            ([-0.002, 0.001, -0.001, -0.002, 0.001, 0.0], 0.2),
            ([0.0; 6], 0.001),
        ] {
            let observed = actual.advance(strain, dt, || false).unwrap();
            let control = expected.advance(strain, dt / a, || false).unwrap();
            for (left, right) in observed.stress.into_iter().zip(control.stress) {
                close(left, right, 5e-13);
            }
            close(observed.work, control.work, 5e-13);
            close(observed.dissipated, control.dissipated, 5e-13);
            close(observed.stored_energy, control.stored_energy, 5e-13);
            assert!(observed.dissipated >= 0.0);
            assert!(observed.balance_residual.abs() < 1e-9);
            shifted.law().step(&mut mandel, strain, dt).unwrap();
            reference
                .step(&mut mandel_reference, strain, dt / a)
                .unwrap();
            close(mandel.stored, mandel_reference.stored, 5e-13);
            close(mandel.dissipated, mandel_reference.dissipated, 5e-13);
        }
        let saved = actual.clone();
        assert!(actual.advance([0.1; 6], 0.01, || false).is_err());
        assert_eq!(actual, saved);
        assert!(actual.advance([0.0; 6], 0.01, || true).is_err());
        assert_eq!(actual, saved);
    }
}

#[test]
fn g3_heating_accelerates_the_same_material_relaxation() {
    let reference = tensor(false);
    let cold = shift().shift_tensor(&reference, 280.0).unwrap();
    let hot = shift().shift_tensor(&reference, 330.0).unwrap();
    let mut cold_state = cold.law().state();
    let mut hot_state = hot.law().state();
    let strain = [0.001, 0.0, 0.0, 0.0, 0.0, 0.0];
    cold.law().step(&mut cold_state, strain, 0.0).unwrap();
    hot.law().step(&mut hot_state, strain, 0.0).unwrap();
    assert_eq!(cold_state.stored, hot_state.stored);
    cold.law().step(&mut cold_state, strain, 0.1).unwrap();
    hot.law().step(&mut hot_state, strain, 0.1).unwrap();
    assert!(hot_state.stress[0] < cold_state.stress[0]);
    assert!(hot_state.dissipated > cold_state.dissipated);
    close(
        hot_state.work,
        hot_state.stored + hot_state.dissipated,
        1e-12,
    );
    close(
        cold_state.work,
        cold_state.stored + cold_state.dissipated,
        1e-12,
    );
}

#[test]
fn g0_shifted_frequency_support_and_unrepresentable_kinetics_are_honest() {
    let reference = tensor(true);
    let saved = reference.clone();
    for temperature in [280.0, 330.0] {
        let shifted = shift().shift_tensor(&reference, temperature).unwrap();
        let a = shifted.state().factor();
        let band = shifted.law().band_rad_s();
        assert_eq!(band, (1.0 / a, 1e5 / a));
        assert!(shifted.law().complex_modulus(band.0).is_ok());
        assert!(shifted.law().complex_modulus(band.1).is_ok());
        assert!(matches!(
            shifted.law().complex_modulus(band.0.next_down()),
            Err(ViscoError::OutOfBand { .. })
        ));
        assert!(matches!(
            shifted.law().complex_modulus(band.1.next_up()),
            Err(ViscoError::OutOfBand { .. })
        ));
    }
    for (tau, temperature) in [(f64::MAX, 280.0), (f64::from_bits(1), 330.0)] {
        let law = TensorProny::new(
            *reference.c_eq(),
            None,
            vec![PronyBranch::new(reference.branches()[0].c_pa, tau).unwrap()],
            (0.0, f64::MAX),
        )
        .unwrap();
        let saved = law.clone();
        assert!(shift().shift_tensor(&law, temperature).is_err());
        assert_eq!(law, saved);
    }
    let unlimited = reference.clone().with_band((0.0, f64::MAX)).unwrap();
    let hotter = shift().shift_tensor(&unlimited, 330.0).unwrap();
    assert_eq!(hotter.law().band_rad_s(), (0.0, f64::MAX));
    assert_eq!(reference, saved);
}
