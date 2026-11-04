//! Independent analytic finite-sensor controls and a point-quadrature oracle.
use super::*;
use fs_material::porous::layers::{
    AcousticWindowBudget, MAX_ACOUSTIC_WINDOW_PAIRS, RectangularAcousticWindow,
};

fn pairs(n: usize) -> AcousticWindowBudget {
    AcousticWindowBudget {
        max_pairs: n * (n + 1) / 2,
    }
}
fn rectangle(center_m: [f64; 3], extent_m: [f64; 2]) -> RectangularAcousticWindow {
    RectangularAcousticWindow { center_m, extent_m }
}
fn sinc(x: f64) -> f64 {
    if x == 0.0 { 1.0 } else { x.sin() / x }
}

#[test]
fn finite_sensor_distinguishes_mean_pressure_from_mean_energy() {
    let fluid = fluid();
    let stack = TransmittingLayers::new(&fluid, &fluid, OMEGA, &[]).unwrap();
    let direction = [0.36, 0.48, 0.8];
    let amplitude = C64::new(0.7, -0.3);
    let beam = stack
        .prepare_angular_spectrum(
            &[IncidentPlaneWave {
                direction,
                pressure_pa: amplitude,
            }],
            budget(1),
            &mut || false,
        )
        .unwrap();
    let k = OMEGA / 340.0;
    let window = rectangle([0.07, -0.02, 0.11], [0.21, 0.17]);
    let sample = beam
        .transmitted_window(window, pairs(1), &mut || false)
        .unwrap();
    let phase = phase(
        -k * direction
            .into_iter()
            .zip(window.center_m)
            .map(|(a, b)| a * b)
            .sum::<f64>(),
    );
    let factor = sinc(0.5 * k * direction[0] * window.extent_m[0])
        * sinc(0.5 * k * direction[1] * window.extent_m[1]);
    near(
        sample.mean_pressure_pa,
        amplitude * phase.scale(factor),
        3e-15,
    );
    near_real(
        sample.mean_squared_pressure_pa2,
        amplitude.norm_sq() * 0.5,
        2e-16,
    );
    for axis in 0..3 {
        near_real(
            sample.active_intensity_w_m2[axis],
            amplitude.norm_sq() * 0.5 * direction[axis] / 408.0,
            3e-18,
        );
        near_real(sample.reactive_intensity_w_m2[axis], 0.0, 3e-18);
    }
    assert!(sample.mean_squared_pressure_pa2 > sample.mean_pressure_pa.norm_sq() * 0.5);
}

#[test]
fn crossing_beam_window_retains_exact_pressure_and_reactive_interference() {
    let fluid = fluid();
    let stack = TransmittingLayers::new(&fluid, &fluid, OMEGA, &[]).unwrap();
    let waves = [-0.6, 0.6].map(|x| IncidentPlaneWave {
        direction: [x, 0.0, 0.8],
        pressure_pa: C64::ONE,
    });
    let beam = stack
        .prepare_angular_spectrum(&waves, budget(2), &mut || false)
        .unwrap();
    let q = OMEGA / 340.0 * 0.6;
    for width in [1e-8, 0.11, core::f64::consts::PI / q] {
        let window = rectangle([0.023, 0.0, 0.04], [width, 0.2]);
        let sample = beam
            .transmitted_window(window, pairs(2), &mut || false)
            .unwrap();
        let interference = (2.0 * q * window.center_m[0]).cos() * sinc(q * width);
        near(
            sample.mean_pressure_pa,
            phase(-OMEGA / 340.0 * 0.8 * 0.04)
                .scale(2.0 * (q * 0.023).cos() * sinc(q * width / 2.0)),
            5e-15,
        );
        near_real(sample.mean_squared_pressure_pa2, 1.0 + interference, 5e-15);
        near_real(
            sample.mean_squared_particle_velocity_m2_s2[0],
            (0.6 / 408.0_f64).powi(2) * (1.0 - interference),
            1e-20,
        );
        near_real(
            sample.active_intensity_w_m2[2],
            0.8 / 408.0 * (1.0 + interference),
            1e-17,
        );
        near_real(sample.active_intensity_w_m2[0], 0.0, 1e-18);
        near_real(
            sample.reactive_intensity_w_m2[0],
            0.6 / 408.0 * (2.0 * q * 0.023).sin() * sinc(q * width),
            3e-18,
        );
    }
}

#[test]
fn subwavelength_dark_fringe_retains_positive_window_energy_and_heat() {
    let fluid = fluid();
    let lossless = TransmittingLayers::new(&fluid, &fluid, OMEGA, &[]).unwrap();
    let waves = [
        IncidentPlaneWave {
            direction: [0.6, 0.0, 0.8],
            pressure_pa: C64::ONE,
        },
        IncidentPlaneWave {
            direction: [-0.6, 0.0, 0.8],
            pressure_pa: -C64::ONE,
        },
    ];
    let beam = lossless
        .prepare_angular_spectrum(&waves, budget(2), &mut || false)
        .unwrap();
    let law = PorousFluidLaw::try_new(PorousCard::basotect_foam02(), 1.81e-5).unwrap();
    let porous = TransmittingLayers::new(
        &fluid,
        &fluid,
        OMEGA,
        &[Layer::Porous {
            law: &law,
            thickness_m: 0.06,
        }],
    )
    .unwrap();
    let porous_beam = porous
        .prepare_angular_spectrum(&waves, budget(2), &mut || false)
        .unwrap();
    let interior = porous_beam.prepare_layer_field(0, &mut || false).unwrap();
    let one_beam = porous
        .prepare_angular_spectrum(&waves[..1], budget(1), &mut || false)
        .unwrap();
    let one_interior = one_beam.prepare_layer_field(0, &mut || false).unwrap();
    let single_heat = one_interior
        .field_with_heat([0.0, 0.0, 0.023], &mut || false)
        .unwrap()
        .thermal_heating_w_m3;
    let q = OMEGA / 340.0 * 0.6;
    for phase_extent in [1e-6, 1e-8, 1e-10] {
        let width = phase_extent / q;
        let window = rectangle([0.0; 3], [width, 0.1]);
        let sample = beam
            .transmitted_window(window, pairs(2), &mut || false)
            .unwrap();
        // p=-2i sin(q*x): <|p|²>/2=1-sinc(q*width). Evaluate the
        // independent small-window series without subtracting rounded ones.
        let square = phase_extent * phase_extent;
        let expected = square / 6.0 * (1.0 - square / 20.0);
        assert_eq!(sample.mean_pressure_pa, C64::ZERO);
        assert!(sample.mean_squared_pressure_pa2 > 0.0);
        near_real(sample.mean_squared_pressure_pa2, expected, expected * 2e-10);
        near_real(
            sample.active_intensity_w_m2[2],
            expected * 0.8 / 408.0,
            expected * 0.8 / 408.0 * 2e-10,
        );
        let heat = interior
            .window_with_heat(
                rectangle([0.0, 0.0, 0.023], [width, 0.1]),
                pairs(2),
                &mut || false,
            )
            .unwrap();
        // Equal/opposite components have identical local normal propagation;
        // their thermal heat retains the same positive lateral fringe factor.
        let expected_heat = 2.0 * single_heat * expected;
        assert!(heat.thermal_heating_w_m3 > 0.0);
        near_real(
            heat.thermal_heating_w_m3,
            expected_heat,
            expected_heat * 2e-10,
        );
    }
}

#[test]
fn window_observes_refraction_interface_work_and_evanescent_decay() {
    let fluid = fluid();
    let receiver = FluidAcousticState::try_new(0.6, 680.0).unwrap();
    let stack = TransmittingLayers::new(&fluid, &receiver, OMEGA, &[]).unwrap();
    let waves = [
        IncidentPlaneWave {
            direction: [0.8, 0.0, 0.6],
            pressure_pa: C64::new(0.8, 0.2),
        },
        IncidentPlaneWave {
            direction: [0.0, 0.0, 1.0],
            pressure_pa: C64::new(-0.3, 0.7),
        },
    ];
    let beam = stack
        .prepare_angular_spectrum(&waves, budget(2), &mut || false)
        .unwrap();
    let face = rectangle([0.031, -0.01, 0.0], [0.2, 0.1]);
    let left = beam.incident_window(face, pairs(2), &mut || false).unwrap();
    let right = beam
        .transmitted_window(face, pairs(2), &mut || false)
        .unwrap();
    near(left.mean_pressure_pa, right.mean_pressure_pa, 3e-15);
    near(
        left.mean_particle_velocity_m_s[2],
        right.mean_particle_velocity_m_s[2],
        2e-17,
    );
    near_real(
        left.active_intensity_w_m2[2],
        right.active_intensity_w_m2[2],
        1e-17,
    );
    near_real(
        left.reactive_intensity_w_m2[2],
        right.reactive_intensity_w_m2[2],
        1e-17,
    );
    let evanescent = stack
        .prepare_angular_spectrum(&waves[..1], budget(1), &mut || false)
        .unwrap();
    let at_face = evanescent
        .transmitted_window(face, pairs(1), &mut || false)
        .unwrap();
    let depth = 0.13;
    let inside = evanescent
        .transmitted_window(
            rectangle([face.center_m[0], face.center_m[1], depth], face.extent_m),
            pairs(1),
            &mut || false,
        )
        .unwrap();
    let alpha = (OMEGA / 340.0) * f64::sqrt(0.8 * 0.8 - 0.25);
    let decay = (-2.0 * alpha * depth).exp();
    near_real(
        inside.mean_squared_pressure_pa2,
        at_face.mean_squared_pressure_pa2 * decay,
        5e-17,
    );
    near_real(inside.active_intensity_w_m2[2], 0.0, 1e-18);
    near_real(
        inside.reactive_intensity_w_m2[2],
        alpha / (0.6 * OMEGA) * inside.mean_squared_pressure_pa2,
        2e-18,
    );
}

#[test]
fn porous_window_heat_matches_integrated_point_fields_and_keeps_frozen_sources() {
    let incident = CountedFluid {
        state: fluid(),
        samples: Cell::new(0),
    };
    let receiver = MutableReceiver {
        state: Cell::new(FluidAcousticState::try_new(1.8, 420.0).unwrap()),
        samples: Cell::new(0),
    };
    let law = PorousFluidLaw::try_new(PorousCard::basotect_foam02(), 1.81e-5).unwrap();
    let stack = TransmittingLayers::new(
        &incident,
        &receiver,
        OMEGA,
        &[Layer::Porous {
            law: &law,
            thickness_m: 0.06,
        }],
    )
    .unwrap();
    let waves = [
        IncidentPlaneWave {
            direction: [0.36, 0.48, 0.8],
            pressure_pa: C64::new(0.8, 0.2),
        },
        IncidentPlaneWave {
            direction: [-0.6, 0.0, 0.8],
            pressure_pa: C64::new(-0.3, 0.7),
        },
    ];
    let beam = stack
        .prepare_angular_spectrum(&waves, budget(2), &mut || false)
        .unwrap();
    let interior = beam.prepare_layer_field(0, &mut || false).unwrap();
    let window = rectangle([0.031, -0.017, 0.023], [0.08, 0.09]);
    let exact = interior
        .window_with_heat(window, pairs(2), &mut || false)
        .unwrap();
    // Composite Simpson quadrature of the independent public point evaluator.
    let mut reference = [0.0; 5];
    const N: usize = 40;
    for i in 0..=N {
        for j in 0..=N {
            let weight = |n| {
                if n == 0 || n == N {
                    1.0
                } else if n % 2 == 0 {
                    2.0
                } else {
                    4.0
                }
            };
            let position = [
                window.center_m[0] + window.extent_m[0] * (i as f64 / N as f64 - 0.5),
                window.center_m[1] + window.extent_m[1] * (j as f64 / N as f64 - 0.5),
                window.center_m[2],
            ];
            let point = interior.field_with_heat(position, &mut || false).unwrap();
            let values = [
                0.5 * point.field.pressure_pa.norm_sq(),
                point.field.active_intensity_w_m2[2],
                point.field.reactive_intensity_w_m2[0],
                point.viscous_heating_w_m3,
                point.thermal_heating_w_m3,
            ];
            for n in 0..5 {
                reference[n] += values[n] * weight(i) * weight(j) / (9.0 * (N * N) as f64);
            }
        }
    }
    let values = [
        exact.field.mean_squared_pressure_pa2,
        exact.field.active_intensity_w_m2[2],
        exact.field.reactive_intensity_w_m2[0],
        exact.viscous_heating_w_m3,
        exact.thermal_heating_w_m3,
    ];
    for (actual, expected) in values.into_iter().zip(reference) {
        near_real(actual, expected, 2e-8 * expected.abs().max(1e-9));
    }
    assert!(exact.viscous_heating_w_m3 > 0.0 && exact.thermal_heating_w_m3 > 0.0);
    receiver
        .state
        .set(FluidAcousticState::try_new(50.0, 900.0).unwrap());
    assert_eq!(
        interior
            .window_with_heat(window, pairs(2), &mut || false)
            .unwrap(),
        exact
    );
    assert_eq!(incident.samples.get(), 1);
    assert_eq!(receiver.samples.get(), 1);
}

#[test]
fn coherent_window_cancellation_and_lossless_identity_are_exact() {
    let fluid = fluid();
    let stack = TransmittingLayers::new(
        &fluid,
        &fluid,
        OMEGA,
        &[
            Layer::FluidGap { thickness_m: 0.0 },
            Layer::FluidGap { thickness_m: 0.03 },
        ],
    )
    .unwrap();
    let waves = [C64::new(0.7, 0.3), C64::new(-0.7, -0.3)].map(|pressure_pa| IncidentPlaneWave {
        direction: [0.6, 0.0, 0.8],
        pressure_pa,
    });
    let beam = stack
        .prepare_angular_spectrum(&waves, budget(2), &mut || false)
        .unwrap();
    let interior = beam.prepare_layer_field(1, &mut || false).unwrap();
    let window = rectangle([0.01, 0.02, 0.015], [0.03, 0.04]);
    let sample = interior
        .window_with_heat(window, pairs(2), &mut || false)
        .unwrap();
    assert_eq!(sample.field.mean_pressure_pa, C64::ZERO);
    assert_eq!(sample.field.mean_squared_pressure_pa2, 0.0);
    assert_eq!(sample.field.active_intensity_w_m2, [0.0; 3]);
    assert_eq!(sample.field.reactive_intensity_w_m2, [0.0; 3]);
    assert_eq!(
        [sample.viscous_heating_w_m3, sample.thermal_heating_w_m3],
        [0.0; 2]
    );
    let identity = beam.prepare_layer_field(0, &mut || false).unwrap();
    let face = rectangle([0.01, 0.02, 0.0], [0.03, 0.04]);
    assert!(identity.window(face, pairs(2), &mut || false).is_ok());
    assert!(
        identity
            .window_with_heat(face, pairs(2), &mut || false)
            .is_err()
    );
}

#[test]
fn acoustic_windows_bound_pair_work_validate_domains_and_cancel_mid_integral() {
    let fluid = fluid();
    let stack = TransmittingLayers::new(
        &fluid,
        &fluid,
        OMEGA,
        &[Layer::FluidGap { thickness_m: 0.02 }],
    )
    .unwrap();
    let waves = [IncidentPlaneWave {
        direction: [0.0, 0.0, 1.0],
        pressure_pa: C64::ONE,
    }; 3];
    let beam = stack
        .prepare_angular_spectrum(&waves, budget(3), &mut || false)
        .unwrap();
    let valid = rectangle([0.0; 3], [0.1; 2]);
    for max_pairs in [0, 5, MAX_ACOUSTIC_WINDOW_PAIRS + 1] {
        assert!(matches!(
            beam.transmitted_window(valid, AcousticWindowBudget { max_pairs }, &mut || false),
            Err(AcousticSpectrumError::Invalid(_))
        ));
    }
    for extent in [[0.0, 0.1], [-0.1, 0.1], [f64::INFINITY, 0.1]] {
        assert!(
            beam.transmitted_window(rectangle([0.0; 3], extent), pairs(3), &mut || false)
                .is_err()
        );
    }
    assert!(
        beam.incident_window(rectangle([0.0, 0.0, 0.01], [0.1; 2]), pairs(3), &mut || {
            false
        })
        .is_err()
    );
    assert!(
        beam.transmitted_window(
            rectangle([0.0, 0.0, -0.01], [0.1; 2]),
            pairs(3),
            &mut || false
        )
        .is_err()
    );
    assert!(
        beam.transmitted_window(
            rectangle([f64::NAN, 0.0, 0.0], [0.1; 2]),
            pairs(3),
            &mut || false
        )
        .is_err()
    );
    let interior = beam.prepare_layer_field(0, &mut || false).unwrap();
    assert!(
        interior
            .window(
                rectangle([0.0, 0.0, 0.021], [0.1; 2]),
                pairs(3),
                &mut || false
            )
            .is_err()
    );
    let mut calls = 0;
    assert_eq!(
        beam.transmitted_window(valid, pairs(3), &mut || {
            calls += 1;
            calls == 8
        }),
        Err(AcousticSpectrumError::Cancelled)
    );
    assert_eq!(calls, 8);
}
