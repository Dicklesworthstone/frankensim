//! G1/G3/G4: original volume-stress sources observed inside a refracting layer.

use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::quadrupole::{PreparedQuadrupoleObserver, QuadrupolePhasor};
use fs_aeroac::curle3d::refraction::PlanarFluidInterface;
use fs_aeroac::curle3d::refraction::full_wave::layer::interior::{
    PreparedInteriorFullWaveForce, PreparedInteriorFullWaveQuadrupole,
    PreparedLayerInteriorFullWaveMonopole,
};
use fs_aeroac::curle3d::refraction::full_wave::layer::reflection::PreparedLayeredReflectedFullWaveQuadrupole;
use fs_aeroac::curle3d::refraction::full_wave::layer::{
    PlanarFluidLayer, PreparedLayeredFullWaveQuadrupole,
};
use fs_aeroac::curle3d::refraction::full_wave::quadrupole::{
    FullWaveQuadrupoleLimits, PreparedFullWaveQuadrupole,
};
use fs_aeroac::curle3d::refraction::full_wave::{FullWaveForceLimits, FullWaveMonopoleLimits};
use fs_aeroac::observation::HarmonicFieldLine;
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::c64::C64;
use std::f64::consts::TAU;

const AXES: [(usize, usize); 6] = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)];

fn fluid(rho: f64, speed: f64) -> FluidAcousticState {
    FluidAcousticState::try_new(rho, speed).unwrap()
}
fn layer(states: [FluidAcousticState; 3], thickness: f64) -> PlanarFluidLayer {
    PlanarFluidLayer::new(
        states[0],
        states[1],
        states[2],
        [0.0; 3],
        [0.0, 0.0, 1.0],
        thickness,
    )
    .unwrap()
}
fn limits() -> FullWaveQuadrupoleLimits {
    FullWaveQuadrupoleLimits {
        absolute_pressure_tolerance_pa_per_n_m: 2e-6,
        absolute_velocity_tolerance_m_s_per_n_m: 2e-9,
        relative_tolerance: 1e-10,
        maximum_evaluations: 262_144,
        maximum_depth: 40,
    }
}
fn prepare(
    slab: PlanarFluidLayer,
    hz: f64,
    source: [f64; 3],
    observer: [f64; 3],
) -> PreparedInteriorFullWaveQuadrupole {
    PreparedInteriorFullWaveQuadrupole::new(slab, hz, source, observer, limits(), &mut || false)
        .unwrap()
}
fn basis(entry: usize) -> QuadrupolePhasor {
    let mut value = QuadrupolePhasor::default();
    value.stress_moment_n_m[entry] = C64::ONE;
    value
}
fn stress() -> QuadrupolePhasor {
    QuadrupolePhasor {
        stress_moment_n_m: [
            C64::new(0.0007, 0.0001),
            C64::new(-0.0003, 0.0002),
            C64::new(0.0009, -0.0004),
            C64::new(0.00023, -0.00015),
            C64::new(-0.00011, 0.0003),
            C64::new(0.0004, 0.0002),
        ],
    }
}
fn amplitudes(line: HarmonicFieldLine<3>) -> [C64; 4] {
    let [x, y, z] = line.particle_velocity_m_s;
    [line.pressure_pa, x, y, z]
}
fn close(actual: C64, expected: C64, relative: f64, absolute: f64) {
    assert!(
        (actual - expected).abs() <= absolute + relative * expected.abs(),
        "actual={actual:?}, expected={expected:?}, relative={relative:e}, absolute={absolute:e}"
    );
}

#[test]
fn absent_rear_face_retains_original_interface_and_homogeneous_limits() {
    let first = fluid(1.2, 340.0);
    let source = [0.003, -0.005, -0.015];
    let observer = [0.025, 0.008, 0.011];
    let hz = 1100.0;
    for last in [first, fluid(2.4, 340.0), fluid(2.4, 510.0)] {
        let slab = layer([first, last, last], 0.025);
        let original = PreparedFullWaveQuadrupole::new(
            PlanarFluidInterface::new(first, last, [0.0; 3], [0.0, 0.0, 1.0]).unwrap(),
            hz,
            source,
            observer,
            limits(),
            &mut || false,
        )
        .unwrap();
        let inside = prepare(slab, hz, source, observer);
        assert!(inside.diagnostics().reduced_to_single_interface);
        assert_eq!(
            inside.diagnostics().evaluations,
            original.diagnostics().evaluations
        );
        assert_eq!(inside.source_position_m(), source);
        assert_eq!(inside.observer_position_m(), observer);
        assert_eq!(inside.frequency_hz(), hz);
        assert_eq!(inside.layer().layer_state(), last);
        for entry in 0..6 {
            assert_eq!(
                amplitudes(inside.evaluate_harmonic_line(basis(entry)).unwrap()),
                amplitudes(original.evaluate_harmonic_line(basis(entry)).unwrap()),
            );
        }
    }
    let same = PreparedInteriorFullWaveQuadrupole::new(
        layer([first; 3], 0.025),
        hz,
        source,
        observer,
        FullWaveQuadrupoleLimits {
            maximum_evaluations: 1,
            ..limits()
        },
        &mut || false,
    )
    .unwrap();
    assert_eq!(same.diagnostics().evaluations, 0);
    let direct = PreparedQuadrupoleObserver::new(first, hz, observer, &[source], 1).unwrap();
    assert_eq!(
        amplitudes(same.evaluate_harmonic_line(stress()).unwrap()),
        amplitudes(direct.evaluate_harmonic_line(&[stress()]).unwrap()),
    );
}

#[test]
fn all_stresses_match_independent_forward_and_returning_spherical_image_families() {
    let states = [fluid(1.2, 300.0), fluid(5.0, 300.0), fluid(2.4, 300.0)];
    let h = 0.018;
    let hz = 700.0;
    let source = [0.002, -0.003, -0.025];
    let observer = [0.035, 0.01, 0.011];
    let inside = prepare(layer(states, h), hz, source, observer);
    assert!(!inside.diagnostics().reduced_to_single_interface);
    let transmission = 2.0 * 5.0 / (1.2 + 5.0);
    let first_return = (1.2 - 5.0) / (1.2 + 5.0);
    let second_return = (2.4 - 5.0) / (2.4 + 5.0);
    let mut reference = [[C64::ZERO; 4]; 6];
    let mut first_pair = [[C64::ZERO; 4]; 6];
    let mut weight = transmission;
    for image in 0..64 {
        let shift = 2.0 * f64::from(image) * h;
        let forward = [source[0], source[1], source[2] - shift];
        let returning = [source[0], source[1], 2.0 * h - source[2] + shift];
        // Homogeneous reference uses the actual receiving-layer density.
        // An odd image reflection applies S M S, S=diag(1,1,-1); each
        // complete internal round trip preserves the original tensor.
        let native =
            PreparedQuadrupoleObserver::new(states[1], hz, observer, &[forward, returning], 2)
                .unwrap();
        for entry in 0..6 {
            let mut mirrored = basis(entry);
            mirrored.stress_moment_n_m[entry] =
                C64::from_re(second_return * if entry >= 4 { -1.0 } else { 1.0 });
            let values = amplitudes(
                native
                    .evaluate_harmonic_line(&[basis(entry), mirrored])
                    .unwrap(),
            );
            for row in 0..4 {
                let value = values[row].scale(weight);
                reference[entry][row] = reference[entry][row] + value;
                if image == 0 {
                    first_pair[entry][row] = value;
                }
            }
        }
        weight *= first_return * second_return;
    }
    for entry in 0..6 {
        let actual = amplitudes(inside.evaluate_harmonic_line(basis(entry)).unwrap());
        for row in 0..4 {
            close(
                actual[row],
                reference[entry][row],
                4e-8,
                if row == 0 { 3e-7 } else { 3e-10 },
            );
        }
    }
    assert!((reference[2][0] - first_pair[2][0]).abs() > 1e-3 * reference[2][0].abs());
}

#[test]
fn original_source_hessians_match_native_force_derivatives_in_both_speed_orders() {
    let source = [0.004, -0.006, -0.012];
    let observer = [0.028, 0.009, 0.008];
    let hz = 950.0;
    let step = 4e-5;
    for states in [
        [fluid(1.1, 300.0), fluid(2.4, 600.0), fluid(2.7, 400.0)],
        [fluid(1.1, 600.0), fluid(2.4, 450.0), fluid(2.7, 300.0)],
    ] {
        let slab = layer(states, 0.02);
        let inside = prepare(slab, hz, source, observer);
        let mut gradient = [[[C64::ZERO; 3]; 4]; 3];
        for direction in 0..3 {
            for (offset, weight) in [(-2.0, 1.0), (-1.0, -8.0), (1.0, 8.0), (2.0, -1.0)] {
                let mut position = source;
                position[direction] += offset * step;
                let force = PreparedInteriorFullWaveForce::new(
                    slab,
                    hz,
                    position,
                    observer,
                    FullWaveForceLimits {
                        absolute_pressure_tolerance_pa_per_n: 1e-9,
                        absolute_velocity_tolerance_m_s_per_n: 1e-12,
                        relative_tolerance: 1e-11,
                        maximum_evaluations: 262_144,
                        maximum_depth: 40,
                    },
                    &mut || false,
                )
                .unwrap();
                for column in 0..3 {
                    let mut amplitude = [C64::ZERO; 3];
                    amplitude[column] = C64::ONE;
                    let values = amplitudes(force.evaluate_harmonic_line(amplitude).unwrap());
                    for row in 0..4 {
                        gradient[direction][row][column] = gradient[direction][row][column]
                            + values[row].scale(weight / (12.0 * step));
                    }
                }
            }
        }
        for (entry, (i, j)) in AXES.into_iter().enumerate() {
            let actual = amplitudes(inside.evaluate_harmonic_line(basis(entry)).unwrap());
            for row in 0..4 {
                let expected = if i == j {
                    gradient[j][row][i]
                } else {
                    gradient[j][row][i] + gradient[i][row][j]
                };
                close(
                    actual[row],
                    expected,
                    3e-6,
                    if row == 0 { 3e-5 } else { 3e-8 },
                );
            }
        }
    }
}

#[test]
fn both_original_faces_match_exterior_pressure_and_normal_velocity_for_every_tensor_entry() {
    let states = [fluid(1.2, 300.0), fluid(2.4, 600.0), fluid(1.8, 400.0)];
    let h = 0.018;
    let slab = layer(states, h);
    let hz = 870.0;
    let source = [0.003, -0.004, -0.022];
    let epsilon = 1e-9;
    let before = PreparedLayeredReflectedFullWaveQuadrupole::new(
        slab,
        hz,
        source,
        [0.021, 0.013, -epsilon],
        limits(),
        &mut || false,
    )
    .unwrap();
    let near = prepare(slab, hz, source, [0.021, 0.013, epsilon]);
    let far = prepare(slab, hz, source, [0.021, 0.013, h - epsilon]);
    let after = PreparedLayeredFullWaveQuadrupole::new(
        slab,
        hz,
        source,
        [0.021, 0.013, h + epsilon],
        limits(),
        &mut || false,
    )
    .unwrap();
    for entry in 0..6 {
        let a = amplitudes(before.evaluate_total_harmonic_line(basis(entry)).unwrap());
        let b = amplitudes(near.evaluate_harmonic_line(basis(entry)).unwrap());
        let c = amplitudes(far.evaluate_harmonic_line(basis(entry)).unwrap());
        let d = amplitudes(after.evaluate_harmonic_line(basis(entry)).unwrap());
        for (left, right) in [(a, b), (c, d)] {
            close(left[0], right[0], 3e-6, 5e-6);
            close(left[3], right[3], 3e-6, 5e-9);
        }
        // Continuous tangential pressure gradient gives rho_left*u_t_left =
        // rho_right*u_t_right; tangential velocity itself need not agree.
        for row in [1, 2] {
            close(a[row].scale(1.2), b[row].scale(2.4), 3e-6, 5e-9);
            close(c[row].scale(2.4), d[row].scale(1.8), 3e-6, 5e-9);
        }
    }
}

#[test]
fn trace_axis_shear_and_receiving_momentum_retain_original_media_and_clock() {
    let states = [fluid(1.2, 600.0), fluid(2.4, 450.0), fluid(1.8, 300.0)];
    let slab = layer(states, 0.018);
    let source = [0.0, 0.0, -0.014];
    let hz = 870.0;
    let omega = TAU * hz;
    let k1 = omega / states[0].sound_speed_m_s();
    let moment = C64::new(0.001, -0.0002);
    let trace = QuadrupolePhasor {
        stress_moment_n_m: [moment, moment, moment, C64::ZERO, C64::ZERO, C64::ZERO],
    };
    let q = C64::new(0.0, -k1 * k1 / (omega * states[0].density_kg_m3())) * moment;
    let observer = [0.031, 0.015, 0.008];
    for point in [[0.0, 0.0, 0.008], observer] {
        let actual = amplitudes(
            prepare(slab, hz, source, point)
                .evaluate_harmonic_line(trace)
                .unwrap(),
        );
        let native = PreparedLayerInteriorFullWaveMonopole::new(
            slab,
            hz,
            source,
            point,
            FullWaveMonopoleLimits {
                absolute_pressure_tolerance_pa_per_m3_s: 1e-8,
                absolute_velocity_tolerance_per_m2: 1e-11,
                relative_tolerance: 1e-11,
                maximum_evaluations: 262_144,
                maximum_depth: 40,
            },
            &mut || false,
        )
        .unwrap();
        let expected = amplitudes(native.evaluate_harmonic_line(q).unwrap());
        for row in 0..4 {
            close(
                actual[row],
                expected[row],
                1e-6,
                if row == 0 { 3e-8 } else { 3e-11 },
            );
        }
    }
    let axis = prepare(slab, hz, source, [0.0, 0.0, 0.008]);
    let xz = axis.evaluate_harmonic_line(basis(4)).unwrap();
    let yz = axis.evaluate_harmonic_line(basis(5)).unwrap();
    assert_eq!(xz.pressure_pa, C64::ZERO);
    assert_eq!(yz.pressure_pa, C64::ZERO);
    assert!(xz.particle_velocity_m_s[0].abs() > 1e-3);
    assert_eq!(xz.particle_velocity_m_s[1], C64::ZERO);
    assert_eq!(xz.particle_velocity_m_s[2], C64::ZERO);
    assert_eq!(yz.particle_velocity_m_s[0], C64::ZERO);
    assert_eq!(yz.particle_velocity_m_s[2], C64::ZERO);
    assert_eq!(xz.particle_velocity_m_s[0], yz.particle_velocity_m_s[1]);

    let original = prepare(slab, hz, source, observer);
    let field = original.evaluate(stress()).unwrap();
    let line = original.evaluate_harmonic_line(stress()).unwrap();
    assert_eq!(line.frequency_hz, hz);
    assert_eq!(line.pressure_pa, field.pressure_pa);
    assert_eq!(line.particle_velocity_m_s, field.particle_velocity_m_s);
    let rho2 = states[1].density_kg_m3();
    let c2 = states[1].sound_speed_m_s();
    let step = 2e-6;
    let mut divergence = C64::ZERO;
    for direction in 0..3 {
        let mut low = observer;
        let mut high = observer;
        low[direction] -= step;
        high[direction] += step;
        let low = prepare(slab, hz, source, low)
            .evaluate_harmonic_line(stress())
            .unwrap();
        let high = prepare(slab, hz, source, high)
            .evaluate_harmonic_line(stress())
            .unwrap();
        close(
            (high.pressure_pa - low.pressure_pa).scale(0.5 / step),
            C64::new(0.0, omega * rho2) * field.particle_velocity_m_s[direction],
            3e-5,
            4e-5,
        );
        divergence = divergence
            + (high.particle_velocity_m_s[direction] - low.particle_velocity_m_s[direction])
                .scale(0.5 / step);
    }
    close(
        divergence,
        C64::new(0.0, omega / (rho2 * c2 * c2)) * field.pressure_pa,
        3e-5,
        5e-8,
    );
    let scale = 7.0;
    let scaled = prepare(
        layer(
            states.map(|s| fluid(scale * s.density_kg_m3(), s.sound_speed_m_s())),
            0.018,
        ),
        hz,
        source,
        observer,
    )
    .evaluate(stress())
    .unwrap();
    close(scaled.pressure_pa, field.pressure_pa, 3e-8, 3e-8);
    for direction in 0..3 {
        close(
            scaled.particle_velocity_m_s[direction].scale(scale),
            field.particle_velocity_m_s[direction],
            3e-8,
            3e-11,
        );
        let intensity =
            field.pressure_pa.scale(0.5) * field.particle_velocity_m_s[direction].conj();
        assert_eq!(field.active_intensity_w_m2[direction], intensity.re);
        assert_eq!(field.reactive_intensity_w_m2[direction], intensity.im);
    }
}

#[test]
fn original_domains_tensor_tail_norms_and_bounded_work_refuse_partial_fields() {
    let states = [fluid(1.2, 300.0), fluid(2.4, 600.0), fluid(1.8, 400.0)];
    let slab = layer(states, 0.018);
    let source = [0.0, 0.0, -0.025];
    let observer = [0.04, 0.0, 0.011];
    let inside = prepare(slab, 700.0, source, observer);
    let evidence = inside.diagnostics();
    assert!(evidence.evaluations > 0 && evidence.evaluations <= limits().maximum_evaluations);
    let norm = |values: &[f64]| values.iter().map(|v| v * v).sum::<f64>().sqrt();
    assert!(
        norm(&evidence.pressure_tail_bound_pa_per_n_m)
            <= 0.251 * limits().absolute_pressure_tolerance_pa_per_n_m
    );
    assert!(
        norm(&evidence.velocity_tail_bound_m_s_per_n_m)
            <= 0.251 * limits().absolute_velocity_tolerance_m_s_per_n_m
    );
    assert_eq!(
        amplitudes(
            inside
                .evaluate_harmonic_line(QuadrupolePhasor::default())
                .unwrap()
        ),
        [C64::ZERO; 4]
    );
    let mut bad = stress();
    bad.stress_moment_n_m[5] = C64::new(f64::NAN, 0.0);
    assert!(inside.evaluate(bad).is_err());
    for point in [
        [0.0, 0.0, -1.0],
        [0.0; 3],
        [0.0, 0.0, 0.018],
        [0.0, 0.0, 0.019],
    ] {
        assert!(
            PreparedInteriorFullWaveQuadrupole::new(
                slab,
                700.0,
                source,
                point,
                limits(),
                &mut || false
            )
            .is_err()
        );
    }
    for point in [[0.0; 3], [0.0, 0.0, 0.01], [f64::NAN, 0.0, -0.02]] {
        assert!(
            PreparedInteriorFullWaveQuadrupole::new(
                slab,
                700.0,
                point,
                observer,
                limits(),
                &mut || false
            )
            .is_err()
        );
    }
    for hz in [0.0, -1.0, f64::INFINITY] {
        assert!(
            PreparedInteriorFullWaveQuadrupole::new(
                slab,
                hz,
                source,
                observer,
                limits(),
                &mut || false
            )
            .is_err()
        );
    }
    assert!(
        PreparedInteriorFullWaveQuadrupole::new(
            slab,
            700.0,
            source,
            observer,
            FullWaveQuadrupoleLimits {
                maximum_evaluations: 30,
                ..limits()
            },
            &mut || false,
        )
        .is_err()
    );
    assert!(matches!(
        PreparedInteriorFullWaveQuadrupole::new(
            slab,
            700.0,
            source,
            observer,
            limits(),
            &mut || true,
        ),
        Err(AeroacError::Cancelled { .. })
    ));
    let mut polls = 0;
    assert!(matches!(
        PreparedInteriorFullWaveQuadrupole::new(
            slab,
            700.0,
            source,
            observer,
            limits(),
            &mut || {
                polls += 1;
                polls == 9
            },
        ),
        Err(AeroacError::Cancelled { .. })
    ));
    assert!(
        PlanarFluidLayer::new(
            states[0],
            fluid(1.0, 200.0),
            states[2],
            [0.0; 3],
            [0.0, 0.0, 1.0],
            0.018
        )
        .is_err()
    );
}
