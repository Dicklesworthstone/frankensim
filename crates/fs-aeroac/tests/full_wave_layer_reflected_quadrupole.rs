//! G1/G3/G4: complete physical stresses through coherent finite-layer returns.

use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::quadrupole::{PreparedQuadrupoleObserver, QuadrupolePhasor};
use fs_aeroac::curle3d::refraction::PlanarFluidInterface;
use fs_aeroac::curle3d::refraction::full_wave::layer::PlanarFluidLayer;
use fs_aeroac::curle3d::refraction::full_wave::layer::reflection::{
    PreparedLayeredReflectedFullWaveForce, PreparedLayeredReflectedFullWaveMonopole,
    PreparedLayeredReflectedFullWaveQuadrupole,
};
use fs_aeroac::curle3d::refraction::full_wave::quadrupole::FullWaveQuadrupoleLimits;
use fs_aeroac::curle3d::refraction::full_wave::reflection::PreparedReflectedFullWaveQuadrupole;
use fs_aeroac::curle3d::refraction::full_wave::{FullWaveForceLimits, FullWaveMonopoleLimits};
use fs_aeroac::observation::HarmonicFieldLine;
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::c64::C64;
use std::f64::consts::TAU;

const ENTRIES: [(usize, usize); 6] = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)];

fn fluid(rho: f64, speed: f64) -> FluidAcousticState {
    FluidAcousticState::try_new(rho, speed).unwrap()
}
fn slab(states: [FluidAcousticState; 3], thickness: f64) -> PlanarFluidLayer {
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
    layer: PlanarFluidLayer,
    hz: f64,
    source: [f64; 3],
    receiver: [f64; 3],
) -> PreparedLayeredReflectedFullWaveQuadrupole {
    PreparedLayeredReflectedFullWaveQuadrupole::new(
        layer,
        hz,
        source,
        receiver,
        limits(),
        &mut || false,
    )
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
fn values(line: HarmonicFieldLine<3>) -> [C64; 4] {
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
fn absent_second_face_and_homogeneous_layers_preserve_the_native_reflected_tensor() {
    let first = fluid(1.2, 300.0);
    let last = fluid(2.4, 450.0);
    let source = [-0.01, 0.02, -0.025];
    let receiver = [0.035, -0.01, -0.041];
    let plane = PlanarFluidInterface::new(first, last, [0.0; 3], [0.0, 0.0, 1.0]).unwrap();
    let original = PreparedReflectedFullWaveQuadrupole::new(
        plane,
        700.0,
        source,
        receiver,
        limits(),
        &mut || false,
    )
    .unwrap();
    for layer in [
        slab([first, fluid(0.8, 200.0), last], 0.0),
        slab([first, last, last], 0.018),
    ] {
        let returned = prepare(layer, 700.0, source, receiver);
        assert!(returned.diagnostics().reduced_to_single_interface);
        assert_eq!(
            returned.diagnostics().evaluations,
            original.diagnostics().evaluations
        );
        for entry in 0..7 {
            let moment = if entry < 6 { basis(entry) } else { stress() };
            assert_eq!(
                values(returned.evaluate_reflected_harmonic_line(moment).unwrap()),
                values(original.evaluate_reflected_harmonic_line(moment).unwrap()),
            );
            assert_eq!(
                values(returned.evaluate_total_harmonic_line(moment).unwrap()),
                values(original.evaluate_total_harmonic_line(moment).unwrap()),
            );
        }
    }
    let homogeneous = PreparedLayeredReflectedFullWaveQuadrupole::new(
        slab([first; 3], 0.02),
        700.0,
        source,
        receiver,
        FullWaveQuadrupoleLimits {
            maximum_evaluations: 1,
            ..limits()
        },
        &mut || false,
    )
    .unwrap();
    assert_eq!(homogeneous.diagnostics().evaluations, 0);
    assert_eq!(
        values(
            homogeneous
                .evaluate_reflected_harmonic_line(stress())
                .unwrap()
        ),
        [C64::ZERO; 4]
    );
    let direct = PreparedQuadrupoleObserver::new(first, 700.0, receiver, &[source], 1).unwrap();
    assert_eq!(
        values(homogeneous.evaluate_total_harmonic_line(stress()).unwrap()),
        values(direct.evaluate_harmonic_line(&[stress()]).unwrap()),
    );
}

fn reflected_stress(moment: QuadrupolePhasor, normal: [f64; 3]) -> QuadrupolePhasor {
    // Independent tensor reflection S M S, with both physical shear entries.
    let mut matrix = [[C64::ZERO; 3]; 3];
    for (value, (i, j)) in moment.stress_moment_n_m.into_iter().zip(ENTRIES) {
        matrix[i][j] = value;
        matrix[j][i] = value;
    }
    let mirror = |i, j| (if i == j { 1.0 } else { 0.0 }) - 2.0 * normal[i] * normal[j];
    QuadrupolePhasor {
        stress_moment_n_m: ENTRIES.map(|(i, j)| {
            let mut value = C64::ZERO;
            for a in 0..3 {
                for b in 0..3 {
                    value = value + matrix[a][b].scale(mirror(i, a) * mirror(j, b));
                }
            }
            value
        }),
    }
}

#[test]
fn rotated_all_component_image_series_retains_internal_returns_and_total_interference() {
    let (r1, r2, r3, speed, h, hz) = (1.2, 5.0, 2.4, 350.0, 0.08, 1100.0);
    let media = [fluid(r1, speed), fluid(r2, speed), fluid(r3, speed)];
    let layer =
        PlanarFluidLayer::new(media[0], media[1], media[2], [0.0; 3], [2.0, 1.0, 2.0], h).unwrap();
    let n = layer.receiving_normal();
    let tangent = [-1.0 / 5.0_f64.sqrt(), 2.0 / 5.0_f64.sqrt(), 0.0];
    let (a, b, range) = (0.03, 0.04, 0.06);
    let source = n.map(|x| -a * x);
    let receiver = core::array::from_fn(|i| range * tangent[i] - b * n[i]);
    let returned = prepare(layer, hz, source, receiver);
    let r12 = (r2 - r1) / (r2 + r1);
    let r23 = (r3 - r2) / (r3 + r2);
    let mut coefficient = (1.0 - r12 * r12) * r23;
    let mut positions = Vec::new();
    let mut coefficients = Vec::new();
    for image in 0..80 {
        positions.push(n.map(|x| (a + 2.0 * f64::from(image) * h) * x));
        coefficients.push(if image == 0 { r12 } else { coefficient });
        if image > 0 {
            coefficient *= -r12 * r23;
        }
    }
    assert!(coefficient.abs() < 1e-50);
    let images =
        PreparedQuadrupoleObserver::new(media[0], hz, receiver, &positions, positions.len())
            .unwrap();
    let direct = PreparedQuadrupoleObserver::new(media[0], hz, receiver, &[source], 1).unwrap();
    for entry in 0..7 {
        let moment = if entry < 6 { basis(entry) } else { stress() };
        let image = reflected_stress(moment, n);
        let image_sources: Vec<_> = coefficients
            .iter()
            .map(|&coefficient| QuadrupolePhasor {
                stress_moment_n_m: image.stress_moment_n_m.map(|v| v.scale(coefficient)),
            })
            .collect();
        let expected_return = values(images.evaluate_harmonic_line(&image_sources).unwrap());
        let expected_direct = values(direct.evaluate_harmonic_line(&[moment]).unwrap());
        for total in [false, true] {
            let expected: [C64; 4] = core::array::from_fn(|row| {
                expected_return[row]
                    + if total {
                        expected_direct[row]
                    } else {
                        C64::ZERO
                    }
            });
            let actual = if total {
                returned.evaluate_total(moment)
            } else {
                returned.evaluate_reflected(moment)
            }
            .unwrap();
            close(actual.pressure_pa, expected[0], 2e-7, 3e-6);
            for axis in 0..3 {
                close(
                    actual.particle_velocity_m_s[axis],
                    expected[axis + 1],
                    2e-7,
                    3e-9,
                );
                let flux = expected[0].scale(0.5) * expected[axis + 1].conj();
                close(
                    C64::new(
                        actual.active_intensity_w_m2[axis],
                        actual.reactive_intensity_w_m2[axis],
                    ),
                    flux,
                    5e-7,
                    3e-8,
                );
            }
        }
    }
    assert!(returned.has_total_field());
    assert_eq!(returned.source_position_m(), source);
    assert_eq!(returned.observer_position_m(), receiver);
    assert_eq!(returned.layer().thickness_m(), h);
    let line = returned.evaluate_total_harmonic_line(stress()).unwrap();
    assert_eq!(line.frequency_hz, hz);
    assert_eq!(
        line.pressure_pa,
        returned.evaluate_total(stress()).unwrap().pressure_pa
    );
}

#[test]
fn every_stress_entry_matches_independent_source_derivatives_of_layer_force_fields() {
    let source = [0.004, -0.006, -0.025];
    let receiver = [0.04, 0.009, -0.023];
    let hz = 700.0;
    let step = 4e-5;
    for states in [
        [fluid(1.2, 300.0), fluid(2.4, 600.0), fluid(1.8, 400.0)],
        [fluid(2.4, 600.0), fluid(1.4, 450.0), fluid(1.2, 300.0)],
    ] {
        let layer = slab(states, 0.018);
        let returned = prepare(layer, hz, source, receiver);
        let mut gradient = [[[C64::ZERO; 3]; 4]; 3];
        for direction in 0..3 {
            for (offset, weight) in [(-2.0, 1.0), (-1.0, -8.0), (1.0, 8.0), (2.0, -1.0)] {
                let mut position = source;
                position[direction] += offset * step;
                let force = PreparedLayeredReflectedFullWaveForce::new(
                    layer,
                    hz,
                    position,
                    receiver,
                    FullWaveForceLimits {
                        absolute_pressure_tolerance_pa_per_n: 1e-10,
                        absolute_velocity_tolerance_m_s_per_n: 1e-13,
                        relative_tolerance: 1e-11,
                        maximum_evaluations: 262_144,
                        maximum_depth: 40,
                    },
                    &mut || false,
                )
                .unwrap();
                for column in 0..3 {
                    let amplitude =
                        core::array::from_fn(|i| if i == column { C64::ONE } else { C64::ZERO });
                    let field = values(force.evaluate_reflected_harmonic_line(amplitude).unwrap());
                    for row in 0..4 {
                        gradient[direction][row][column] = gradient[direction][row][column]
                            + field[row].scale(weight / (12.0 * step));
                    }
                }
            }
        }
        for (entry, (i, j)) in ENTRIES.into_iter().enumerate() {
            let actual = values(
                returned
                    .evaluate_reflected_harmonic_line(basis(entry))
                    .unwrap(),
            );
            for row in 0..4 {
                let expected = gradient[j][row][i]
                    + if i != j {
                        gradient[i][row][j]
                    } else {
                        C64::ZERO
                    };
                close(
                    actual[row],
                    expected,
                    4e-6,
                    if row == 0 { 3e-5 } else { 3e-8 },
                );
            }
        }
    }
}

#[test]
fn original_source_helmholtz_trace_and_axis_shear_survive_layer_reflection() {
    let states = [fluid(1.2, 600.0), fluid(2.4, 450.0), fluid(1.8, 300.0)];
    let layer = slab(states, 0.018);
    let source = [0.0, 0.0, -0.025];
    let hz = 870.0;
    let omega = TAU * hz;
    let k1 = omega / states[0].sound_speed_m_s();
    let moment = C64::new(0.001, -0.0002);
    let trace = QuadrupolePhasor {
        stress_moment_n_m: [moment, moment, moment, C64::ZERO, C64::ZERO, C64::ZERO],
    };
    let q = C64::new(0.0, -k1 * k1 / (omega * states[0].density_kg_m3())) * moment;
    for receiver in [
        [0.0, 0.0, -0.041],
        [1e-9, -2e-9, -0.041],
        [0.047, 0.019, -0.041],
    ] {
        let actual = values(
            prepare(layer, hz, source, receiver)
                .evaluate_reflected_harmonic_line(trace)
                .unwrap(),
        );
        let mono = PreparedLayeredReflectedFullWaveMonopole::new(
            layer,
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
        let expected = values(mono.evaluate_reflected_harmonic_line(q).unwrap());
        for row in 0..4 {
            close(
                actual[row],
                expected[row],
                1e-6,
                if row == 0 { 3e-8 } else { 3e-11 },
            );
        }
    }
    let axis = prepare(layer, hz, source, [0.0, 0.0, -0.041]);
    let xz = axis.evaluate_reflected_harmonic_line(basis(4)).unwrap();
    let yz = axis.evaluate_reflected_harmonic_line(basis(5)).unwrap();
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
fn receiving_momentum_and_continuity_keep_original_density_and_return_direction() {
    let states = [fluid(1.2, 300.0), fluid(2.4, 600.0), fluid(1.8, 400.0)];
    let layer = slab(states, 0.018);
    let source = [0.002, -0.004, -0.025];
    let receiver = [0.04, 0.015, -0.023];
    let hz = 1200.0;
    let moment = stress();
    let returned = prepare(layer, hz, source, receiver);
    let actual = returned.evaluate_reflected(moment).unwrap();
    let rho1 = states[0].density_kg_m3();
    let c1 = states[0].sound_speed_m_s();
    let step = 1e-5;
    let mut divergence = C64::ZERO;
    for axis in 0..3 {
        let mut low = receiver;
        let mut high = receiver;
        low[axis] -= step;
        high[axis] += step;
        let low = prepare(layer, hz, source, low)
            .evaluate_reflected_harmonic_line(moment)
            .unwrap();
        let high = prepare(layer, hz, source, high)
            .evaluate_reflected_harmonic_line(moment)
            .unwrap();
        close(
            (high.pressure_pa - low.pressure_pa).scale(0.5 / step),
            C64::new(0.0, TAU * hz * rho1) * actual.particle_velocity_m_s[axis],
            3e-5,
            4e-5,
        );
        divergence = divergence
            + (high.particle_velocity_m_s[axis] - low.particle_velocity_m_s[axis])
                .scale(0.5 / step);
    }
    close(
        divergence,
        C64::new(0.0, TAU * hz / (rho1 * c1 * c1)) * actual.pressure_pa,
        3e-5,
        5e-8,
    );
    let scale = 7.0;
    let scaled = prepare(
        slab(
            states.map(|s| fluid(scale * s.density_kg_m3(), s.sound_speed_m_s())),
            0.018,
        ),
        hz,
        source,
        receiver,
    )
    .evaluate_reflected(moment)
    .unwrap();
    close(scaled.pressure_pa, actual.pressure_pa, 3e-8, 3e-8);
    for axis in 0..3 {
        close(
            scaled.particle_velocity_m_s[axis].scale(scale),
            actual.particle_velocity_m_s[axis],
            3e-8,
            3e-11,
        );
    }
    // An absent first face gives the same return from the actual second plane.
    let matched = prepare(
        slab([states[0], states[0], states[2]], 0.018),
        hz,
        source,
        receiver,
    );
    let plane = PlanarFluidInterface::new(states[0], states[2], [0.0, 0.0, 0.018], [0.0, 0.0, 1.0])
        .unwrap();
    let shifted = PreparedReflectedFullWaveQuadrupole::new(
        plane,
        hz,
        source,
        receiver,
        limits(),
        &mut || false,
    )
    .unwrap();
    let a = values(matched.evaluate_reflected_harmonic_line(moment).unwrap());
    let b = values(shifted.evaluate_reflected_harmonic_line(moment).unwrap());
    for row in 0..4 {
        close(a[row], b[row], 3e-8, if row == 0 { 3e-8 } else { 3e-11 });
    }
}

#[test]
fn finite_return_self_fields_tail_norms_and_work_refusals_are_explicit() {
    let layer = slab(
        [fluid(1.2, 300.0), fluid(2.4, 600.0), fluid(1.8, 400.0)],
        0.018,
    );
    let source = [0.0, 0.0, -0.025];
    let receiver = [0.04, 0.0, -0.023];
    let returned = prepare(layer, 700.0, source, receiver);
    let norm = |values: &[f64]| values.iter().map(|v| v * v).sum::<f64>().sqrt();
    assert!(
        norm(&returned.diagnostics().pressure_tail_bound_pa_per_n_m)
            <= 0.251 * limits().absolute_pressure_tolerance_pa_per_n_m
    );
    assert!(
        norm(&returned.diagnostics().velocity_tail_bound_m_s_per_n_m)
            <= 0.251 * limits().absolute_velocity_tolerance_m_s_per_n_m
    );
    assert!(returned.diagnostics().evaluations <= limits().maximum_evaluations);
    assert!(!returned.diagnostics().reduced_to_single_interface);
    let self_field = prepare(layer, 700.0, source, source);
    assert!(!self_field.has_total_field());
    for entry in 0..6 {
        assert!(
            values(
                self_field
                    .evaluate_reflected_harmonic_line(basis(entry))
                    .unwrap()
            )
            .iter()
            .all(|v| v.re.is_finite() && v.im.is_finite())
        );
    }
    assert!(
        self_field
            .evaluate_total(QuadrupolePhasor::default())
            .is_err()
    );
    assert!(self_field.evaluate_total_harmonic_line(basis(0)).is_err());
    assert_eq!(
        values(
            returned
                .evaluate_reflected_harmonic_line(QuadrupolePhasor::default())
                .unwrap()
        ),
        [C64::ZERO; 4]
    );
    let mut invalid = stress();
    invalid.stress_moment_n_m[5] = C64::new(f64::NAN, 0.0);
    assert!(returned.evaluate_reflected(invalid).is_err());
    assert!(
        PreparedLayeredReflectedFullWaveQuadrupole::new(
            layer,
            700.0,
            source,
            receiver,
            FullWaveQuadrupoleLimits {
                maximum_evaluations: 1,
                ..limits()
            },
            &mut || false,
        )
        .is_err()
    );
    assert!(matches!(
        PreparedLayeredReflectedFullWaveQuadrupole::new(
            layer,
            700.0,
            source,
            receiver,
            limits(),
            &mut || true,
        ),
        Err(AeroacError::Cancelled { .. })
    ));
    let mut polls = 0;
    assert!(matches!(
        PreparedLayeredReflectedFullWaveQuadrupole::new(
            layer,
            700.0,
            source,
            receiver,
            limits(),
            &mut || {
                polls += 1;
                polls == 9
            },
        ),
        Err(AeroacError::Cancelled { .. })
    ));
    for position in [[0.0; 3], [0.0, 0.0, 0.02], [0.0, 0.0, f64::NAN]] {
        assert!(
            PreparedLayeredReflectedFullWaveQuadrupole::new(
                layer,
                700.0,
                source,
                position,
                limits(),
                &mut || false
            )
            .is_err()
        );
        assert!(
            PreparedLayeredReflectedFullWaveQuadrupole::new(
                layer,
                700.0,
                position,
                receiver,
                limits(),
                &mut || false
            )
            .is_err()
        );
    }
    assert!(
        PreparedLayeredReflectedFullWaveQuadrupole::new(
            layer,
            0.0,
            source,
            receiver,
            limits(),
            &mut || false
        )
        .is_err()
    );
}
