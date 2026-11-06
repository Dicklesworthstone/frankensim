//! G1/G3: original interior stresses, both exterior directions and source-fluid derivatives.

use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::quadrupole::{PreparedQuadrupoleObserver, QuadrupolePhasor};
use fs_aeroac::curle3d::refraction::PlanarFluidInterface;
use fs_aeroac::curle3d::refraction::full_wave::layer::PlanarFluidLayer;
use fs_aeroac::curle3d::refraction::full_wave::layer::interior::source::{
    PreparedInteriorSourceFullWaveForce, PreparedInteriorSourceFullWaveMonopole,
    PreparedInteriorSourceFullWaveQuadrupole,
};
use fs_aeroac::curle3d::refraction::full_wave::quadrupole::{
    FullWaveQuadrupoleLimits, PreparedFullWaveQuadrupole,
};
use fs_aeroac::curle3d::refraction::full_wave::{FullWaveForceLimits, FullWaveMonopoleLimits};
use fs_aeroac::observation::HarmonicFieldLine;
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::c64::C64;
use std::f64::consts::TAU;

const ENTRIES: [(usize, usize); 6] = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)];
fn fluid(rho: f64, c: f64) -> FluidAcousticState {
    FluidAcousticState::try_new(rho, c).unwrap()
}
fn slab(media: [FluidAcousticState; 3], h: f64) -> PlanarFluidLayer {
    PlanarFluidLayer::new(media[0], media[1], media[2], [0.0; 3], [0.0, 0.0, 1.0], h).unwrap()
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
fn force_limits() -> FullWaveForceLimits {
    FullWaveForceLimits {
        absolute_pressure_tolerance_pa_per_n: 1e-10,
        absolute_velocity_tolerance_m_s_per_n: 1e-13,
        relative_tolerance: 1e-11,
        maximum_evaluations: 262_144,
        maximum_depth: 40,
    }
}
fn q_limits() -> FullWaveMonopoleLimits {
    FullWaveMonopoleLimits {
        absolute_pressure_tolerance_pa_per_m3_s: 1e-10,
        absolute_velocity_tolerance_per_m2: 1e-13,
        relative_tolerance: 1e-11,
        maximum_evaluations: 262_144,
        maximum_depth: 40,
    }
}
fn prepare(
    layer: PlanarFluidLayer,
    hz: f64,
    source: [f64; 3],
    observer: [f64; 3],
) -> PreparedInteriorSourceFullWaveQuadrupole {
    PreparedInteriorSourceFullWaveQuadrupole::new(
        layer,
        hz,
        source,
        observer,
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
            C64::new(-0.00019, 0.00011),
            C64::new(0.00031, 0.00017),
        ],
    }
}
fn values(line: HarmonicFieldLine<3>) -> [C64; 4] {
    let [x, y, z] = line.particle_velocity_m_s;
    [line.pressure_pa, x, y, z]
}
fn close(actual: C64, expected: C64, relative: f64, absolute: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs() + absolute,
        "actual={actual:?}, expected={expected:?}, relative={relative:e}, absolute={absolute:e}"
    );
}
fn fields_close(actual: [C64; 4], expected: [C64; 4], relative: f64, pressure: f64, velocity: f64) {
    for i in 0..4 {
        close(
            actual[i],
            expected[i],
            relative,
            if i == 0 { pressure } else { velocity },
        );
    }
}

#[test]
fn absent_rear_face_preserves_the_original_reversed_interface_and_homogeneous_field() {
    let outside = fluid(1.2, 300.0);
    let source = [0.006, -0.003, 0.011];
    let receiver = [0.029, 0.014, -0.025];
    for inside in [outside, fluid(3.0, 300.0), fluid(2.4, 450.0)] {
        let layer = slab([outside, inside, inside], 0.018);
        let prepared = prepare(layer, 700.0, source, receiver);
        let original = PreparedFullWaveQuadrupole::new(
            PlanarFluidInterface::new(inside, outside, [0.0; 3], [0.0, 0.0, -1.0]).unwrap(),
            700.0,
            source,
            receiver,
            limits(),
            &mut || false,
        )
        .unwrap();
        for entry in 0..7 {
            let moment = if entry < 6 { basis(entry) } else { stress() };
            assert_eq!(
                values(prepared.evaluate_harmonic_line(moment).unwrap()),
                values(original.evaluate_harmonic_line(moment).unwrap()),
            );
        }
        assert!(prepared.diagnostics().reduced_to_single_interface);
        assert_eq!(
            prepared.diagnostics().evaluations,
            original.diagnostics().evaluations
        );
        if inside.sound_speed_m_s() == outside.sound_speed_m_s() {
            let mut exact = limits();
            exact.maximum_evaluations = 1;
            assert_eq!(
                PreparedInteriorSourceFullWaveQuadrupole::new(
                    layer,
                    700.0,
                    source,
                    receiver,
                    exact,
                    &mut || false,
                )
                .unwrap()
                .diagnostics()
                .evaluations,
                0,
            );
        }
    }
    let homogeneous = prepare(slab([outside; 3], 0.018), 700.0, source, receiver);
    let direct = PreparedQuadrupoleObserver::new(outside, 700.0, receiver, &[source], 1).unwrap();
    fields_close(
        values(homogeneous.evaluate_harmonic_line(stress()).unwrap()),
        values(direct.evaluate_harmonic_line(&[stress()]).unwrap()),
        3e-13,
        1e-12,
        1e-14,
    );
}

fn reflected_stress(moment: QuadrupolePhasor, normal: [f64; 3]) -> QuadrupolePhasor {
    // Independent physical tensor reflection S M S, with S=I-2nn^T.
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
fn rotated_full_tensor_image_families_radiate_into_both_original_exteriors() {
    let media = [fluid(1.2, 350.0), fluid(5.0, 350.0), fluid(2.4, 350.0)];
    let (h, s, a, r, hz) = (0.08, 0.031, 0.03, 0.06, 1100.0);
    let origin = [0.13, -0.2, 0.07];
    let original =
        PlanarFluidLayer::new(media[0], media[1], media[2], origin, [2.0, 1.0, 2.0], h).unwrap();
    let normal = original.receiving_normal();
    let tangent = [-1.0 / 5.0_f64.sqrt(), 2.0 / 5.0_f64.sqrt(), 0.0];
    let source = core::array::from_fn(|i| origin[i] + s * normal[i]);
    for opposite in [false, true] {
        let layer = if opposite {
            PlanarFluidLayer::new(
                media[2],
                media[1],
                media[0],
                core::array::from_fn(|i| origin[i] + h * normal[i]),
                normal.map(|v| -v),
                h,
            )
            .unwrap()
        } else {
            original
        };
        let n = layer.receiving_normal();
        let face = layer.first_face_point_m();
        let depth = if opposite { h - s } else { s };
        let receiver = core::array::from_fn(|i| face[i] - a * n[i] + r * tangent[i]);
        let prepared = prepare(layer, hz, source, receiver);
        let rho1 = layer.incident_state().density_kg_m3();
        let rho2 = layer.layer_state().density_kg_m3();
        let rho3 = layer.receiving_state().density_kg_m3();
        let r21 = (rho1 - rho2) / (rho1 + rho2);
        let r23 = (rho3 - rho2) / (rho3 + rho2);
        let mut gain = 2.0 * rho1 / (rho1 + rho2);
        let mut positions = Vec::new();
        let mut gains = Vec::new();
        for order in 0..80 {
            let shift = 2.0 * f64::from(order) * h;
            positions.push(core::array::from_fn(|i| source[i] + shift * n[i]));
            gains.push(gain);
            positions.push(core::array::from_fn(|i| {
                source[i] + (2.0 * (h - depth) + shift) * n[i]
            }));
            gains.push(gain * r23);
            gain *= r21 * r23;
        }
        assert!(gain.abs() < 1e-50);
        let images = PreparedQuadrupoleObserver::new(
            layer.incident_state(),
            hz,
            receiver,
            &positions,
            positions.len(),
        )
        .unwrap();
        for entry in 0..7 {
            let moment = if entry < 6 { basis(entry) } else { stress() };
            let mirror = reflected_stress(moment, n);
            let amplitudes: Vec<_> = gains
                .iter()
                .enumerate()
                .map(|(index, &gain)| {
                    let tensor = if index % 2 == 0 { moment } else { mirror };
                    QuadrupolePhasor {
                        stress_moment_n_m: tensor.stress_moment_n_m.map(|v| v.scale(gain)),
                    }
                })
                .collect();
            let expected = values(images.evaluate_harmonic_line(&amplitudes).unwrap());
            let actual = prepared.evaluate(moment).unwrap();
            let line = prepared.evaluate_harmonic_line(moment).unwrap();
            assert_eq!(line.frequency_hz, hz);
            assert_eq!(line.pressure_pa, actual.pressure_pa);
            assert_eq!(line.particle_velocity_m_s, actual.particle_velocity_m_s);
            fields_close(values(line), expected, 5e-8, 2e-5, 2e-8);
            for i in 0..3 {
                let flux = expected[0].scale(0.5) * expected[i + 1].conj();
                close(
                    C64::new(
                        actual.active_intensity_w_m2[i],
                        actual.reactive_intensity_w_m2[i],
                    ),
                    flux,
                    2e-7,
                    2e-6,
                );
            }
        }
        assert_eq!(prepared.source_position_m(), source);
        assert_eq!(prepared.observer_position_m(), receiver);
        assert_eq!(prepared.frequency_hz(), hz);
        assert_eq!(prepared.layer().layer_state(), layer.layer_state());
        assert!(!prepared.diagnostics().reduced_to_single_interface);
        assert!(prepared.diagnostics().evaluations > 0);
    }
}

#[test]
fn every_stress_entry_differentiates_the_native_original_interior_force() {
    // Independent fourth-order source-coordinate derivatives of the existing
    // seven-channel force owner. Physical off-diagonal moments contribute both
    // derivative permutations; no product finite differences are introduced.
    for media in [
        [fluid(1.2, 300.0), fluid(2.4, 600.0), fluid(1.8, 400.0)],
        [fluid(2.4, 600.0), fluid(1.4, 450.0), fluid(1.2, 300.0)],
    ] {
        let layer = slab(media, 0.025);
        let (source, receiver, hz) = ([0.004, -0.003, 0.011], [0.04, 0.015, -0.025], 900.0);
        let prepared = prepare(layer, hz, source, receiver);
        let step = 2e-5;
        let mut derivative = [[[C64::ZERO; 4]; 3]; 3];
        for axis in 0..3 {
            for (offset, coefficient) in [(-2.0, 1.0), (-1.0, -8.0), (1.0, 8.0), (2.0, -1.0)] {
                let mut shifted = source;
                shifted[axis] += offset * step;
                let force = PreparedInteriorSourceFullWaveForce::new(
                    layer,
                    hz,
                    shifted,
                    receiver,
                    force_limits(),
                    &mut || false,
                )
                .unwrap();
                for direction in 0..3 {
                    let unit =
                        core::array::from_fn(|i| if i == direction { C64::ONE } else { C64::ZERO });
                    let response = values(force.evaluate_harmonic_line(unit).unwrap());
                    for row in 0..4 {
                        derivative[axis][direction][row] = derivative[axis][direction][row]
                            + response[row].scale(coefficient / (12.0 * step));
                    }
                }
            }
        }
        for (entry, (i, j)) in ENTRIES.into_iter().enumerate() {
            let expected = core::array::from_fn(|row| {
                derivative[i][j][row]
                    + if i == j {
                        C64::ZERO
                    } else {
                        derivative[j][i][row]
                    }
            });
            fields_close(
                values(prepared.evaluate_harmonic_line(basis(entry)).unwrap()),
                expected,
                3e-6,
                2e-5,
                2e-8,
            );
        }
    }
}

#[test]
fn source_helmholtz_trace_uses_layer_wavenumber_and_axis_shear_remains_physical() {
    let media = [fluid(1.2, 300.0), fluid(2.4, 600.0), fluid(1.8, 400.0)];
    let layer = slab(media, 0.025);
    let source = [0.0, 0.0, 0.011];
    let hz = 700.0;
    let omega = TAU * hz;
    let k2 = omega / media[1].sound_speed_m_s();
    let amplitude = C64::new(0.001, -0.0002);
    let trace = QuadrupolePhasor {
        stress_moment_n_m: [
            amplitude,
            amplitude,
            amplitude,
            C64::ZERO,
            C64::ZERO,
            C64::ZERO,
        ],
    };
    let q = C64::new(0.0, -k2 * k2 / (omega * media[1].density_kg_m3())) * amplitude;
    for receiver in [
        [0.0, 0.0, -0.025],
        [1e-9, -2e-9, -0.025],
        [0.04, 0.015, -0.025],
    ] {
        let quadrupole = prepare(layer, hz, source, receiver);
        let monopole = PreparedInteriorSourceFullWaveMonopole::new(
            layer,
            hz,
            source,
            receiver,
            q_limits(),
            &mut || false,
        )
        .unwrap();
        fields_close(
            values(quadrupole.evaluate_harmonic_line(trace).unwrap()),
            values(monopole.evaluate_harmonic_line(q).unwrap()),
            3e-7,
            1e-7,
            1e-10,
        );
    }
    let axis = prepare(layer, hz, source, [0.0, 0.0, -0.025]);
    let xz = axis.evaluate(basis(4)).unwrap();
    let yz = axis.evaluate(basis(5)).unwrap();
    assert_eq!(xz.pressure_pa, C64::ZERO);
    assert_eq!(xz.particle_velocity_m_s[1], C64::ZERO);
    assert_eq!(xz.particle_velocity_m_s[2], C64::ZERO);
    assert!(xz.particle_velocity_m_s[0].abs() > 0.0);
    assert_eq!(yz.pressure_pa, C64::ZERO);
    close(
        xz.particle_velocity_m_s[0],
        yz.particle_velocity_m_s[1],
        3e-13,
        1e-13,
    );
}

#[test]
fn exterior_momentum_continuity_and_uniform_density_scaling_use_original_fluids() {
    let media = [fluid(1.2, 300.0), fluid(2.4, 600.0), fluid(1.8, 400.0)];
    let layer = slab(media, 0.025);
    let (source, receiver, hz) = ([0.004, -0.003, 0.011], [0.04, 0.015, -0.025], 900.0);
    let moment = stress();
    let actual = prepare(layer, hz, source, receiver)
        .evaluate(moment)
        .unwrap();
    let step = 2e-5;
    let mut divergence = C64::ZERO;
    for axis in 0..3 {
        let mut derivative = [C64::ZERO; 4];
        for (offset, coefficient) in [(-2.0, 1.0), (-1.0, -8.0), (1.0, 8.0), (2.0, -1.0)] {
            let mut shifted = receiver;
            shifted[axis] += offset * step;
            let field = values(
                prepare(layer, hz, source, shifted)
                    .evaluate_harmonic_line(moment)
                    .unwrap(),
            );
            for row in 0..4 {
                derivative[row] = derivative[row] + field[row].scale(coefficient / (12.0 * step));
            }
        }
        close(
            derivative[0],
            C64::new(0.0, TAU * hz * media[0].density_kg_m3()) * actual.particle_velocity_m_s[axis],
            3e-6,
            3e-5,
        );
        divergence = divergence + derivative[axis + 1];
    }
    close(
        divergence,
        C64::new(
            0.0,
            TAU * hz / (media[0].density_kg_m3() * media[0].sound_speed_m_s().powi(2)),
        ) * actual.pressure_pa,
        3e-6,
        5e-8,
    );
    let scaled = slab(
        media.map(|state| fluid(7.0 * state.density_kg_m3(), state.sound_speed_m_s())),
        0.025,
    );
    let scaled = prepare(scaled, hz, source, receiver)
        .evaluate(moment)
        .unwrap();
    close(scaled.pressure_pa, actual.pressure_pa, 2e-7, 1e-7);
    for axis in 0..3 {
        close(
            scaled.particle_velocity_m_s[axis].scale(7.0),
            actual.particle_velocity_m_s[axis],
            2e-7,
            1e-10,
        );
    }
}

#[test]
fn tail_norms_strict_original_domains_and_bounded_work_are_explicit() {
    let media = [fluid(1.2, 300.0), fluid(2.4, 600.0), fluid(1.8, 400.0)];
    let layer = slab(media, 0.025);
    let source = [0.004, -0.003, 0.011];
    let receiver = [0.04, 0.015, -0.025];
    let prepared = prepare(layer, 900.0, source, receiver);
    let diagnostics = prepared.diagnostics();
    assert!(diagnostics.evaluations > 0 && diagnostics.evaluations <= limits().maximum_evaluations);
    let pressure_tail = diagnostics
        .pressure_tail_bound_pa_per_n_m
        .iter()
        .map(|x| x * x)
        .sum::<f64>()
        .sqrt();
    let velocity_tail = diagnostics
        .velocity_tail_bound_m_s_per_n_m
        .iter()
        .map(|x| x * x)
        .sum::<f64>()
        .sqrt();
    assert!(pressure_tail <= 0.251 * limits().absolute_pressure_tolerance_pa_per_n_m);
    assert!(velocity_tail <= 0.251 * limits().absolute_velocity_tolerance_m_s_per_n_m);
    assert_eq!(
        values(
            prepared
                .evaluate_harmonic_line(QuadrupolePhasor::default())
                .unwrap()
        ),
        [C64::ZERO; 4]
    );
    let mut invalid = stress();
    invalid.stress_moment_n_m[5] = C64::new(f64::NAN, 0.0);
    assert!(prepared.evaluate(invalid).is_err());
    let mut small = limits();
    small.maximum_evaluations = 1;
    assert!(
        PreparedInteriorSourceFullWaveQuadrupole::new(
            layer,
            900.0,
            source,
            receiver,
            small,
            &mut || false
        )
        .is_err()
    );
    for z in [-0.01, 0.0, 0.025, 0.03, f64::NAN] {
        assert!(
            PreparedInteriorSourceFullWaveQuadrupole::new(
                layer,
                900.0,
                [0.0, 0.0, z],
                receiver,
                limits(),
                &mut || false
            )
            .is_err()
        );
    }
    for z in [0.0, 0.01, 0.025, 0.03, f64::NAN] {
        assert!(
            PreparedInteriorSourceFullWaveQuadrupole::new(
                layer,
                900.0,
                source,
                [0.0, 0.0, z],
                limits(),
                &mut || false
            )
            .is_err()
        );
    }
    assert!(
        PreparedInteriorSourceFullWaveQuadrupole::new(
            slab(media, 0.0),
            900.0,
            source,
            receiver,
            limits(),
            &mut || false
        )
        .is_err()
    );
    assert!(
        PreparedInteriorSourceFullWaveQuadrupole::new(
            layer,
            0.0,
            source,
            receiver,
            limits(),
            &mut || false
        )
        .is_err()
    );
    assert!(matches!(
        PreparedInteriorSourceFullWaveQuadrupole::new(
            layer,
            900.0,
            source,
            receiver,
            limits(),
            &mut || true
        ),
        Err(AeroacError::Cancelled { .. }),
    ));
    let mut polls = 0;
    assert!(matches!(
        PreparedInteriorSourceFullWaveQuadrupole::new(
            layer,
            900.0,
            source,
            receiver,
            limits(),
            &mut || {
                polls += 1;
                polls == 9
            }
        ),
        Err(AeroacError::Cancelled { .. }),
    ));
    assert_eq!(polls, 9);
}
