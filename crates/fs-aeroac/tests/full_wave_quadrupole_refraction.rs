//! G1/G3/G4: original volume-stress radiation through a refracting fluid plane.

use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::quadrupole::{PreparedQuadrupoleObserver, QuadrupolePhasor};
use fs_aeroac::curle3d::refraction::PlanarFluidInterface;
use fs_aeroac::curle3d::refraction::full_wave::quadrupole::{
    FullWaveQuadrupoleLimits, PreparedFullWaveQuadrupole,
};
use fs_aeroac::curle3d::refraction::full_wave::{
    FullWaveForceLimits, FullWaveMonopoleLimits, PreparedFullWaveForce, PreparedFullWaveMonopole,
};
use fs_aeroac::observation::HarmonicFieldLine;
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::c64::C64;
use std::f64::consts::TAU;

const AXES: [(usize, usize); 6] = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)];

fn fluid(rho: f64, speed: f64) -> FluidAcousticState {
    FluidAcousticState::try_new(rho, speed).unwrap()
}
fn interface(rho: [f64; 2], speed: [f64; 2]) -> PlanarFluidInterface {
    PlanarFluidInterface::new(
        fluid(rho[0], speed[0]),
        fluid(rho[1], speed[1]),
        [0.0; 3],
        [0.0, 0.0, 1.0],
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
fn force_limits() -> FullWaveForceLimits {
    FullWaveForceLimits {
        absolute_pressure_tolerance_pa_per_n: 1e-9,
        absolute_velocity_tolerance_m_s_per_n: 1e-12,
        relative_tolerance: 1e-11,
        maximum_evaluations: 262_144,
        maximum_depth: 40,
    }
}
fn prepare(
    interface: PlanarFluidInterface,
    hz: f64,
    source: [f64; 3],
    observer: [f64; 3],
) -> PreparedFullWaveQuadrupole {
    PreparedFullWaveQuadrupole::new(interface, hz, source, observer, limits(), &mut || false)
        .unwrap()
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
fn basis(entry: usize) -> QuadrupolePhasor {
    let mut source = QuadrupolePhasor::default();
    source.stress_moment_n_m[entry] = C64::ONE;
    source
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
fn original_homogeneous_and_density_jump_tensors_are_exact_without_quadrature() {
    let source = [0.003, -0.005, -0.015];
    let receiver = [0.025, 0.008, 0.027];
    let hz = 1100.0;
    let rho1 = 1.2;
    let original =
        PreparedQuadrupoleObserver::new(fluid(rho1, 343.5), hz, receiver, &[source], 1).unwrap();
    for rho2 in [rho1, 0.7, 8.3] {
        let mut controls = limits();
        controls.maximum_evaluations = 1;
        let prepared = PreparedFullWaveQuadrupole::new(
            interface([rho1, rho2], [343.5; 2]),
            hz,
            source,
            receiver,
            controls,
            &mut || false,
        )
        .unwrap();
        assert!(prepared.diagnostics().exact_equal_speed);
        assert_eq!(prepared.diagnostics().evaluations, 0);
        assert_eq!(
            prepared.diagnostics().pressure_tail_bound_pa_per_n_m,
            [0.0; 4]
        );
        assert_eq!(
            prepared.diagnostics().velocity_tail_bound_m_s_per_n_m,
            [0.0; 10]
        );
        for moment in [stress(), basis(0), basis(3), basis(4), basis(5)] {
            let reference = original.evaluate_harmonic_line(&[moment]).unwrap();
            let actual = prepared.evaluate_harmonic_line(moment).unwrap();
            close(
                actual.pressure_pa,
                reference.pressure_pa.scale(2.0 * rho2 / (rho1 + rho2)),
                4e-15,
                0.0,
            );
            for axis in 0..3 {
                close(
                    actual.particle_velocity_m_s[axis],
                    reference.particle_velocity_m_s[axis].scale(2.0 * rho1 / (rho1 + rho2)),
                    4e-15,
                    0.0,
                );
            }
            if rho2 == rho1 {
                assert_eq!(actual.pressure_pa, reference.pressure_pa);
                assert_eq!(
                    actual.particle_velocity_m_s,
                    reference.particle_velocity_m_s
                );
            }
        }
    }
}

#[test]
fn all_physical_stress_entries_match_independent_original_force_source_derivatives() {
    let source = [0.004, -0.006, -0.012];
    let receiver = [0.028, 0.009, 0.018];
    let hz = 950.0;
    let step = 4e-5;
    for speed in [[330.0, 510.0], [510.0, 330.0]] {
        let plane = interface([1.1, 2.7], speed);
        let prepared = prepare(plane, hz, source, receiver);
        // Differentiate each ORIGINAL Cartesian force column with respect to
        // its original source position. No quadrupole kernel, angular moments
        // or receiver displacement derivative supplies this reference.
        let mut gradient = [[[C64::ZERO; 3]; 4]; 3];
        for direction in 0..3 {
            for (offset, weight) in [(-2.0, 1.0), (-1.0, -8.0), (1.0, 8.0), (2.0, -1.0)] {
                let mut position = source;
                position[direction] += offset * step;
                let force = PreparedFullWaveForce::new(
                    plane,
                    hz,
                    position,
                    receiver,
                    force_limits(),
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
                if i != j {
                    // The input is a physical symmetric shear, so both source
                    // contractions survive; no engineering-shear convention.
                    expected = expected + gradient[i][row][j];
                }
                close(
                    actual[row],
                    expected,
                    2e-6,
                    if row == 0 { 2e-5 } else { 2e-8 },
                );
            }
        }
    }
}

#[test]
fn source_trace_satisfies_helmholtz_and_axis_shear_retains_transverse_velocity() {
    let plane = interface([1.2, 2.8], [510.0, 320.0]);
    let source = [0.0, 0.0, -0.014];
    let hz = 870.0;
    let omega = TAU * hz;
    let k1 = omega / plane.incident_state().sound_speed_m_s();
    let moment = C64::new(0.001, -0.0002);
    let trace = QuadrupolePhasor {
        stress_moment_n_m: [moment, moment, moment, C64::ZERO, C64::ZERO, C64::ZERO],
    };
    // Laplacian_source H=-k1²H; pQ=-i omega rho1 H Q. The actual INCIDENT
    // wavenumber, not a refracted direction or receiving wavenumber, is required.
    let q = C64::new(
        0.0,
        -k1 * k1 / (omega * plane.incident_state().density_kg_m3()),
    ) * moment;
    for receiver in [
        [0.0, 0.0, 0.020],
        [1e-9, -2e-9, 0.020],
        [0.047, 0.019, 0.020],
    ] {
        let prepared = prepare(plane, hz, source, receiver);
        let monopole = PreparedFullWaveMonopole::new(
            plane,
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
        let actual = amplitudes(prepared.evaluate_harmonic_line(trace).unwrap());
        let expected = amplitudes(monopole.evaluate_harmonic_line(q).unwrap());
        for row in 0..4 {
            close(
                actual[row],
                expected[row],
                8e-7,
                if row == 0 { 3e-8 } else { 3e-11 },
            );
        }
    }
    let axis = prepare(plane, hz, source, [0.0, 0.0, 0.020]);
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
    let nearby = prepare(plane, hz, source, [1e-9, -2e-9, 0.020])
        .evaluate_harmonic_line(basis(4))
        .unwrap();
    for axis in 0..3 {
        close(
            nearby.particle_velocity_m_s[axis],
            xz.particle_velocity_m_s[axis],
            1e-5,
            2e-5,
        );
    }
}

fn rotate(v: [f64; 3]) -> [f64; 3] {
    [
        0.36 * v[0] - 0.48 * v[1] + 0.8 * v[2],
        0.8 * v[0] + 0.6 * v[1],
        -0.48 * v[0] + 0.64 * v[1] + 0.6 * v[2],
    ]
}
fn rotated_moment(source: QuadrupolePhasor) -> QuadrupolePhasor {
    let rotation = [
        rotate([1.0, 0.0, 0.0]),
        rotate([0.0, 1.0, 0.0]),
        rotate([0.0, 0.0, 1.0]),
    ];
    let mut matrix = [[C64::ZERO; 3]; 3];
    for (entry, (i, j)) in AXES.into_iter().enumerate() {
        matrix[i][j] = source.stress_moment_n_m[entry];
        matrix[j][i] = source.stress_moment_n_m[entry];
    }
    QuadrupolePhasor {
        stress_moment_n_m: AXES.map(|(i, j)| {
            let mut value = C64::ZERO;
            for a in 0..3 {
                for b in 0..3 {
                    value = value + matrix[a][b].scale(rotation[a][i] * rotation[b][j]);
                }
            }
            value
        }),
    }
}

#[test]
fn receiving_momentum_divergence_full_tensor_rotation_and_density_scaling_are_physical() {
    let plane = interface([1.2, 2.8], [320.0, 470.0]);
    let source = [0.002, -0.004, -0.018];
    let receiver = [0.031, 0.015, 0.022];
    let hz = 1200.0;
    let moment = stress();
    let prepared = prepare(plane, hz, source, receiver);
    let actual = prepared.evaluate(moment).unwrap();
    let line = prepared.evaluate_harmonic_line(moment).unwrap();
    assert_eq!(line.frequency_hz, hz);
    assert_eq!(line.pressure_pa, actual.pressure_pa);
    assert_eq!(line.particle_velocity_m_s, actual.particle_velocity_m_s);
    assert!(
        actual
            .reactive_intensity_w_m2
            .iter()
            .any(|v| v.abs() > 1e-8)
    );
    let rho2 = plane.receiving_state().density_kg_m3();
    let c2 = plane.receiving_state().sound_speed_m_s();
    let step = 1e-5;
    let mut divergence = C64::ZERO;
    for axis in 0..3 {
        let mut low = receiver;
        let mut high = receiver;
        low[axis] -= step;
        high[axis] += step;
        let low = prepare(plane, hz, source, low)
            .evaluate_harmonic_line(moment)
            .unwrap();
        let high = prepare(plane, hz, source, high)
            .evaluate_harmonic_line(moment)
            .unwrap();
        let gradient = (high.pressure_pa - low.pressure_pa).scale(0.5 / step);
        close(
            gradient,
            C64::new(0.0, TAU * hz * rho2) * actual.particle_velocity_m_s[axis],
            3e-5,
            4e-5,
        );
        divergence = divergence
            + (high.particle_velocity_m_s[axis] - low.particle_velocity_m_s[axis])
                .scale(0.5 / step);
    }
    close(
        divergence,
        C64::new(0.0, TAU * hz / (rho2 * c2 * c2)) * actual.pressure_pa,
        1e-4,
        2e-6,
    );

    let shift = [0.25, -0.125, 0.5];
    let moved = |v: [f64; 3]| {
        let v = rotate(v);
        core::array::from_fn(|axis| v[axis] + shift[axis])
    };
    let rotated_plane = PlanarFluidInterface::new(
        plane.incident_state(),
        plane.receiving_state(),
        shift,
        rotate([0.0, 0.0, 1.0]),
    )
    .unwrap();
    let rotated = prepare(rotated_plane, hz, moved(source), moved(receiver))
        .evaluate(rotated_moment(moment))
        .unwrap();
    close(rotated.pressure_pa, actual.pressure_pa, 2e-6, 5e-9);
    let real = rotate(actual.particle_velocity_m_s.map(|v| v.re));
    let imaginary = rotate(actual.particle_velocity_m_s.map(|v| v.im));
    for axis in 0..3 {
        close(
            rotated.particle_velocity_m_s[axis],
            C64::new(real[axis], imaginary[axis]),
            2e-6,
            5e-12,
        );
    }
    // Stress has fixed physical N m units. Common density scaling leaves p
    // unchanged and scales actual velocity inversely; no acoustic gain is fitted.
    let scaled = prepare(
        interface([1.2 * 17.0, 2.8 * 17.0], [320.0, 470.0]),
        hz,
        source,
        receiver,
    )
    .evaluate(moment)
    .unwrap();
    close(scaled.pressure_pa, actual.pressure_pa, 5e-7, 5e-9);
    for axis in 0..3 {
        close(
            scaled.particle_velocity_m_s[axis].scale(17.0),
            actual.particle_velocity_m_s[axis],
            5e-7,
            5e-12,
        );
    }
}

#[test]
fn shared_norm_budgets_geometry_cancellation_and_invalid_moments_refuse_atomically() {
    let plane = interface([1.2, 2.8], [320.0, 470.0]);
    let source = [0.002, -0.004, -0.018];
    let receiver = [0.031, 0.015, 0.022];
    let hz = 1200.0;
    let mut controls = limits();
    controls.relative_tolerance = 0.0;
    let mut polls = 0;
    let prepared =
        PreparedFullWaveQuadrupole::new(plane, hz, source, receiver, controls, &mut || {
            polls += 1;
            false
        })
        .unwrap();
    assert_eq!(prepared.source_position_m(), source);
    assert_eq!(prepared.observer_position_m(), receiver);
    assert_eq!(prepared.frequency_hz(), hz);
    assert_eq!(
        prepared.interface().incident_state(),
        plane.incident_state()
    );
    let d = prepared.diagnostics();
    assert!(!d.exact_equal_speed && d.evaluations > 0);
    let pressure_error = (0..4).fold(0.0_f64, |norm, i| {
        norm.hypot(d.pressure_quadrature_error_pa_per_n_m[i] + d.pressure_tail_bound_pa_per_n_m[i])
    });
    let velocity_error = (0..10).fold(0.0_f64, |norm, i| {
        norm.hypot(
            d.velocity_quadrature_error_m_s_per_n_m[i] + d.velocity_tail_bound_m_s_per_n_m[i],
        )
    });
    assert!(pressure_error <= controls.absolute_pressure_tolerance_pa_per_n_m);
    assert!(velocity_error <= controls.absolute_velocity_tolerance_m_s_per_n_m);
    let mut capped = controls;
    capped.maximum_evaluations = d.evaluations - 1;
    assert!(
        PreparedFullWaveQuadrupole::new(plane, hz, source, receiver, capped, &mut || false)
            .is_err()
    );
    assert!(matches!(
        PreparedFullWaveQuadrupole::new(plane, hz, source, receiver, controls, &mut || true),
        Err(AeroacError::Cancelled { .. })
    ));
    let mut late = 0;
    assert!(matches!(
        PreparedFullWaveQuadrupole::new(plane, hz, source, receiver, controls, &mut || {
            late += 1;
            late >= polls
        }),
        Err(AeroacError::Cancelled { .. })
    ));
    for (bad_source, bad_receiver, bad_hz) in [
        ([0.0; 3], receiver, hz),
        (source, [0.0; 3], hz),
        ([0.0, 0.0, 0.01], receiver, hz),
        (source, [0.0, 0.0, -0.01], hz),
        ([f64::NAN, 0.0, -0.01], receiver, hz),
        (source, receiver, 0.0),
        (source, receiver, f64::INFINITY),
    ] {
        assert!(
            PreparedFullWaveQuadrupole::new(
                plane,
                bad_hz,
                bad_source,
                bad_receiver,
                controls,
                &mut || false
            )
            .is_err()
        );
    }
    let before = prepared.evaluate(stress()).unwrap();
    for value in [C64::new(f64::NAN, 0.0), C64::new(0.0, f64::INFINITY)] {
        let mut invalid = stress();
        invalid.stress_moment_n_m[4] = value;
        assert!(prepared.evaluate(invalid).is_err());
        assert!(prepared.evaluate_harmonic_line(invalid).is_err());
    }
    assert_eq!(
        prepared.evaluate(stress()).unwrap().pressure_pa,
        before.pressure_pa
    );
    let silence = prepared.evaluate(QuadrupolePhasor::default()).unwrap();
    assert_eq!(silence.pressure_pa, C64::ZERO);
    assert_eq!(silence.particle_velocity_m_s, [C64::ZERO; 3]);
    assert_eq!(silence.active_intensity_w_m2, [0.0; 3]);
    assert_eq!(silence.reactive_intensity_w_m2, [0.0; 3]);
}
