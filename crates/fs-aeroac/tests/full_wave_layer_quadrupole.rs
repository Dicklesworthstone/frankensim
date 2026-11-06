//! G1/G3: original volume stresses through a coherent finite fluid layer.

use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::quadrupole::{PreparedQuadrupoleObserver, QuadrupolePhasor};
use fs_aeroac::curle3d::refraction::PlanarFluidInterface;
use fs_aeroac::curle3d::refraction::full_wave::layer::{
    PlanarFluidLayer, PreparedLayeredFullWaveForce, PreparedLayeredFullWaveMonopole,
    PreparedLayeredFullWaveQuadrupole,
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
    receiver: [f64; 3],
) -> PreparedLayeredFullWaveQuadrupole {
    PreparedLayeredFullWaveQuadrupole::new(slab, hz, source, receiver, limits(), &mut || false)
        .unwrap()
}
fn basis(entry: usize) -> QuadrupolePhasor {
    let mut moment = QuadrupolePhasor::default();
    moment.stress_moment_n_m[entry] = C64::ONE;
    moment
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
        "actual={actual:?}, expected={expected:?}, rel={relative:e}, abs={absolute:e}"
    );
}

#[test]
fn absent_and_identical_layers_reduce_to_the_original_quadrupole_owners() {
    let source = [0.003, -0.005, -0.015];
    let receiver = [0.025, 0.008, 0.047];
    let hz = 1100.0;
    let first = fluid(1.2, 343.5);
    let last = fluid(2.8, 470.0);
    let plane = PlanarFluidInterface::new(first, last, [0.0; 3], [0.0, 0.0, 1.0]).unwrap();
    let original =
        PreparedFullWaveQuadrupole::new(plane, hz, source, receiver, limits(), &mut || false)
            .unwrap();
    // An absent slow medium has no physical path or guided-mode pole.
    let absent = prepare(
        layer([first, fluid(1000.0, 100.0), last], 0.0),
        hz,
        source,
        receiver,
    );
    assert!(absent.diagnostics().reduced_to_single_interface);
    for entry in 0..6 {
        assert_eq!(
            amplitudes(absent.evaluate_harmonic_line(basis(entry)).unwrap()),
            amplitudes(original.evaluate_harmonic_line(basis(entry)).unwrap())
        );
    }
    let same = PreparedLayeredFullWaveQuadrupole::new(
        layer([first; 3], 0.02),
        hz,
        source,
        receiver,
        FullWaveQuadrupoleLimits {
            maximum_evaluations: 1,
            ..limits()
        },
        &mut || false,
    )
    .unwrap();
    assert!(same.diagnostics().reduced_to_single_interface);
    assert_eq!(same.diagnostics().evaluations, 0);
    let native = PreparedQuadrupoleObserver::new(first, hz, receiver, &[source], 1).unwrap();
    assert_eq!(
        amplitudes(same.evaluate_harmonic_line(stress()).unwrap()),
        amplitudes(native.evaluate_harmonic_line(&[stress()]).unwrap())
    );
}

#[test]
fn equal_speeds_match_independent_complete_native_spherical_image_series() {
    // Angle-independent R/T gives an independent spherical image reference.
    // Each return contains two reflections: the ORIGINAL stress orientation
    // remains intact, while receiving velocity uses the third fluid density.
    let states = [fluid(1.2, 300.0), fluid(5.0, 300.0), fluid(2.4, 300.0)];
    let h = 0.018;
    let hz = 700.0;
    let source = [0.002, -0.003, -0.025];
    let receiver = [0.035, 0.01, h + 0.018];
    let prepared = prepare(layer(states, h), hz, source, receiver);
    assert!(!prepared.diagnostics().reduced_to_single_interface);
    assert!(prepared.diagnostics().evaluations > 0);
    let transmission = (2.0 * 5.0 / (1.2 + 5.0)) * (2.0 * 2.4 / (5.0 + 2.4));
    let ratio = ((1.2 - 5.0) / (1.2 + 5.0)) * ((2.4 - 5.0) / (2.4 + 5.0));
    let mut reference = [[C64::ZERO; 4]; 6];
    let mut first_path = [[C64::ZERO; 4]; 6];
    let mut weight = transmission;
    for image in 0..64 {
        let position = [source[0], source[1], source[2] - 2.0 * f64::from(image) * h];
        let native =
            PreparedQuadrupoleObserver::new(states[0], hz, receiver, &[position], 1).unwrap();
        for entry in 0..6 {
            let values = amplitudes(native.evaluate_harmonic_line(&[basis(entry)]).unwrap());
            for row in 0..4 {
                let value = values[row].scale(weight * if row == 0 { 1.0 } else { 1.2 / 2.4 });
                reference[entry][row] = reference[entry][row] + value;
                if image == 0 {
                    first_path[entry][row] = value;
                }
            }
        }
        weight *= ratio;
    }
    for entry in 0..6 {
        let actual = amplitudes(prepared.evaluate_harmonic_line(basis(entry)).unwrap());
        for row in 0..4 {
            close(actual[row], reference[entry][row], 4e-8, 2e-8);
        }
    }
    assert!((reference[2][0] - first_path[2][0]).abs() > 1e-3 * reference[2][0].abs());
}

#[test]
fn every_stress_entry_matches_independent_layer_force_source_derivatives() {
    let source = [0.004, -0.006, -0.012];
    let receiver = [0.028, 0.009, 0.038];
    let hz = 950.0;
    let step = 4e-5;
    for states in [
        [fluid(1.1, 300.0), fluid(2.4, 600.0), fluid(2.7, 400.0)],
        [fluid(1.1, 600.0), fluid(2.4, 450.0), fluid(2.7, 300.0)],
    ] {
        let slab = layer(states, 0.02);
        let prepared = prepare(slab, hz, source, receiver);
        let mut gradient = [[[C64::ZERO; 3]; 4]; 3];
        for direction in 0..3 {
            for (offset, weight) in [(-2.0, 1.0), (-1.0, -8.0), (1.0, 8.0), (2.0, -1.0)] {
                let mut position = source;
                position[direction] += offset * step;
                let force = PreparedLayeredFullWaveForce::new(
                    slab,
                    hz,
                    position,
                    receiver,
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
                    let value = amplitudes(force.evaluate_harmonic_line(amplitude).unwrap());
                    for row in 0..4 {
                        gradient[direction][row][column] = gradient[direction][row][column]
                            + value[row].scale(weight / (12.0 * step));
                    }
                }
            }
        }
        for (entry, (i, j)) in AXES.into_iter().enumerate() {
            let actual = amplitudes(prepared.evaluate_harmonic_line(basis(entry)).unwrap());
            for row in 0..4 {
                let mut expected = gradient[j][row][i];
                // Physical symmetric shears contribute both contractions.
                if i != j {
                    expected = expected + gradient[i][row][j];
                }
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
fn original_source_trace_and_axis_shear_retain_near_and_transverse_fields() {
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
    // Source Helmholtz: trace(M d_s d_s H)=-k1² M H. A volume source
    // is -i omega rho1 Q H, so this relation uses the ORIGINAL source medium.
    let q = C64::new(0.0, -k1 * k1 / (omega * states[0].density_kg_m3())) * moment;
    for receiver in [
        [0.0, 0.0, 0.038],
        [1e-9, -2e-9, 0.038],
        [0.047, 0.019, 0.038],
    ] {
        let actual = amplitudes(
            prepare(slab, hz, source, receiver)
                .evaluate_harmonic_line(trace)
                .unwrap(),
        );
        let mono = PreparedLayeredFullWaveMonopole::new(
            slab,
            hz,
            source,
            receiver,
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
        let expected = amplitudes(mono.evaluate_harmonic_line(q).unwrap());
        for row in 0..4 {
            close(
                actual[row],
                expected[row],
                1e-6,
                if row == 0 { 3e-8 } else { 3e-11 },
            );
        }
    }
    let axis = prepare(slab, hz, source, [0.0, 0.0, 0.038]);
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
}

#[test]
fn receiving_momentum_continuity_density_scaling_and_original_clock_are_physical() {
    let states = [fluid(1.2, 300.0), fluid(2.4, 600.0), fluid(1.8, 400.0)];
    let slab = layer(states, 0.018);
    let source = [0.002, -0.004, -0.018];
    let receiver = [0.031, 0.015, 0.04];
    let hz = 1200.0;
    let moment = stress();
    let prepared = prepare(slab, hz, source, receiver);
    let actual = prepared.evaluate(moment).unwrap();
    let line = prepared.evaluate_harmonic_line(moment).unwrap();
    assert_eq!(line.frequency_hz, hz);
    assert_eq!(line.pressure_pa, actual.pressure_pa);
    assert_eq!(line.particle_velocity_m_s, actual.particle_velocity_m_s);
    let rho3 = states[2].density_kg_m3();
    let c3 = states[2].sound_speed_m_s();
    let step = 1e-5;
    let mut divergence = C64::ZERO;
    for axis in 0..3 {
        let mut low = receiver;
        let mut high = receiver;
        low[axis] -= step;
        high[axis] += step;
        let low = prepare(slab, hz, source, low)
            .evaluate_harmonic_line(moment)
            .unwrap();
        let high = prepare(slab, hz, source, high)
            .evaluate_harmonic_line(moment)
            .unwrap();
        let gradient = (high.pressure_pa - low.pressure_pa).scale(0.5 / step);
        close(
            gradient,
            C64::new(0.0, TAU * hz * rho3) * actual.particle_velocity_m_s[axis],
            3e-5,
            4e-5,
        );
        divergence = divergence
            + (high.particle_velocity_m_s[axis] - low.particle_velocity_m_s[axis])
                .scale(0.5 / step);
    }
    close(
        divergence,
        C64::new(0.0, TAU * hz / (rho3 * c3 * c3)) * actual.pressure_pa,
        3e-5,
        5e-8,
    );
    let factor = 7.0;
    let scaled = prepare(
        layer(
            states.map(|s| fluid(factor * s.density_kg_m3(), s.sound_speed_m_s())),
            0.018,
        ),
        hz,
        source,
        receiver,
    )
    .evaluate(moment)
    .unwrap();
    close(scaled.pressure_pa, actual.pressure_pa, 3e-8, 3e-8);
    for axis in 0..3 {
        close(
            scaled.particle_velocity_m_s[axis].scale(factor),
            actual.particle_velocity_m_s[axis],
            3e-8,
            3e-11,
        );
        let intensity = actual.pressure_pa.scale(0.5) * actual.particle_velocity_m_s[axis].conj();
        assert_eq!(actual.active_intensity_w_m2[axis], intensity.re);
        assert_eq!(actual.reactive_intensity_w_m2[axis], intensity.im);
    }
}

#[test]
fn error_norms_work_and_cancellation_remain_bounded_before_publication() {
    let states = [fluid(1.2, 300.0), fluid(2.4, 600.0), fluid(1.8, 400.0)];
    let slab = layer(states, 0.018);
    let source = [0.0, 0.0, -0.025];
    let receiver = [0.04, 0.0, 0.041];
    let prepared = prepare(slab, 700.0, source, receiver);
    assert!(prepared.diagnostics().evaluations <= limits().maximum_evaluations);
    let norm = |values: &[f64]| values.iter().map(|v| v * v).sum::<f64>().sqrt();
    assert!(
        norm(&prepared.diagnostics().pressure_tail_bound_pa_per_n_m)
            <= 0.251 * limits().absolute_pressure_tolerance_pa_per_n_m
    );
    assert!(
        norm(&prepared.diagnostics().velocity_tail_bound_m_s_per_n_m)
            <= 0.251 * limits().absolute_velocity_tolerance_m_s_per_n_m
    );
    assert_eq!(
        amplitudes(
            prepared
                .evaluate_harmonic_line(QuadrupolePhasor::default())
                .unwrap()
        ),
        [C64::ZERO; 4]
    );
    let mut invalid = stress();
    invalid.stress_moment_n_m[5] = C64::new(f64::NAN, 0.0);
    assert!(prepared.evaluate(invalid).is_err());
    assert!(
        PreparedLayeredFullWaveQuadrupole::new(
            slab,
            700.0,
            source,
            receiver,
            FullWaveQuadrupoleLimits {
                maximum_evaluations: 30,
                ..limits()
            },
            &mut || false
        )
        .is_err()
    );
    assert!(matches!(
        PreparedLayeredFullWaveQuadrupole::new(
            slab,
            700.0,
            source,
            receiver,
            limits(),
            &mut || true
        ),
        Err(AeroacError::Cancelled { .. })
    ));
    let mut polls = 0;
    assert!(matches!(
        PreparedLayeredFullWaveQuadrupole::new(
            slab,
            700.0,
            source,
            receiver,
            limits(),
            &mut || {
                polls += 1;
                polls == 9
            }
        ),
        Err(AeroacError::Cancelled { .. })
    ));
    for (hz, position) in [
        (0.0, receiver),
        (700.0, [0.0; 3]),
        (700.0, [0.0, 0.0, 0.01]),
    ] {
        assert!(
            PreparedLayeredFullWaveQuadrupole::new(
                slab,
                hz,
                source,
                position,
                limits(),
                &mut || false
            )
            .is_err()
        );
    }
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
