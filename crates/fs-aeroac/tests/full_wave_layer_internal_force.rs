//! G1/G3: physical force derivatives and coherent boundary returns inside a layer.

use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::refraction::PlanarFluidInterface;
use fs_aeroac::curle3d::refraction::full_wave::layer::PlanarFluidLayer;
use fs_aeroac::curle3d::refraction::full_wave::layer::interior::internal::{
    FULL_WAVE_LAYER_INTERNAL_REFLECTED_FORCE_SCOPE, PreparedLayerInternalFullWaveForce,
    PreparedLayerInternalFullWaveMonopole,
};
use fs_aeroac::curle3d::refraction::full_wave::reflection::PreparedReflectedFullWaveForce;
use fs_aeroac::curle3d::refraction::full_wave::{FullWaveForceLimits, FullWaveMonopoleLimits};
use fs_aeroac::curle3d::{CompactSourcePhasor, PreparedCompactSourceObserver};
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::c64::C64;

fn fluid(rho: f64, c: f64) -> FluidAcousticState {
    FluidAcousticState::try_new(rho, c).unwrap()
}
fn slab(c: [f64; 3], h: f64) -> PlanarFluidLayer {
    PlanarFluidLayer::new(
        fluid(1.2, c[0]),
        fluid(2.4, c[1]),
        fluid(1.8, c[2]),
        [0.0; 3],
        [0.0, 0.0, 1.0],
        h,
    )
    .unwrap()
}
fn limits() -> FullWaveForceLimits {
    FullWaveForceLimits {
        absolute_pressure_tolerance_pa_per_n: 1e-8,
        absolute_velocity_tolerance_m_s_per_n: 1e-11,
        relative_tolerance: 1e-9,
        maximum_evaluations: 262_144,
        maximum_depth: 36,
    }
}
fn q_limits() -> FullWaveMonopoleLimits {
    FullWaveMonopoleLimits {
        absolute_pressure_tolerance_pa_per_m3_s: 1e-8,
        absolute_velocity_tolerance_per_m2: 1e-12,
        relative_tolerance: 1e-11,
        maximum_evaluations: 262_144,
        maximum_depth: 36,
    }
}
fn unit(axis: usize) -> [C64; 3] {
    core::array::from_fn(|i| if i == axis { C64::ONE } else { C64::ZERO })
}
fn close(a: C64, b: C64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * b.abs().max(1e-3),
        "actual={a:?}, expected={b:?}"
    );
}

#[test]
fn all_seven_channels_match_independent_35_and_50_digit_direct_q_references() {
    // Explicit r21/r23 path series and source/receiver derivatives integrated
    // directly in q; all values agree to 24 digits under precision refinement.
    for (c, h, s, t, r, expected) in [
        (
            [300.0, 600.0, 400.0],
            0.018,
            0.006,
            0.011,
            0.04,
            [
                C64::new(-15.465344208769404, 0.8804956031122443),
                C64::new(2.9021958372425747, -2.0115839513446443),
                C64::new(0.0019254959009253602, -0.06211983689182278),
                C64::new(0.00208534564889871, 0.0366277675214772),
                C64::new(0.0007411934541203309, 0.015522078500736166),
                C64::new(-0.0007497623568098268, -0.017426731315231735),
                C64::new(0.007067837698148903, -0.03308264026150945),
            ],
        ),
        (
            [600.0, 450.0, 300.0],
            0.02,
            0.014,
            0.004,
            0.025,
            [
                C64::new(-29.351532384452723, -0.07831837465054822),
                C64::new(12.715426650511777, 3.8720074724637965),
                C64::new(-0.0003064208695632187, -0.10226656883623798),
                C64::new(-0.00029677991572756525, 0.11122479681684515),
                C64::new(-0.0006043226311285035, 0.0911542798388669),
                C64::new(0.0006798239459062637, -0.07883025592225426),
                C64::new(0.005898879236231193, 0.010728470980039715),
            ],
        ),
    ] {
        let prepared = PreparedLayerInternalFullWaveForce::new(
            slab(c, h),
            700.0,
            [0.0, 0.0, s],
            [r, 0.0, t],
            limits(),
            &mut || false,
        )
        .unwrap();
        let x = prepared.evaluate_reflected_harmonic_line(unit(0)).unwrap();
        let y = prepared.evaluate_reflected_harmonic_line(unit(1)).unwrap();
        let z = prepared.evaluate_reflected_harmonic_line(unit(2)).unwrap();
        let actual = [
            x.pressure_pa,
            z.pressure_pa,
            x.particle_velocity_m_s[0],
            y.particle_velocity_m_s[1],
            z.particle_velocity_m_s[0],
            x.particle_velocity_m_s[2],
            z.particle_velocity_m_s[2],
        ];
        for i in 0..7 {
            close(actual[i], expected[i], 3e-8);
        }
        assert_eq!(y.pressure_pa, C64::ZERO);
        assert!(prepared.diagnostics().evaluations > 0);
    }
}

#[test]
fn rotated_image_forces_reverse_only_the_odd_reflection_normal_component() {
    let point = [0.3, -0.5, 0.2];
    let normal = [0.0, 0.6, 0.8];
    let h = 0.08;
    let layer = PlanarFluidLayer::new(
        fluid(1.2, 350.0),
        fluid(2.4, 350.0),
        fluid(1.8, 350.0),
        point,
        normal,
        h,
    )
    .unwrap();
    let s = 0.025;
    let t = 0.04;
    let source = core::array::from_fn(|i| point[i] + normal[i] * s);
    let observer =
        core::array::from_fn(|i| point[i] + normal[i] * t + if i == 0 { 0.04 } else { 0.0 });
    let prepared = PreparedLayerInternalFullWaveForce::new(
        layer,
        700.0,
        source,
        observer,
        limits(),
        &mut || false,
    )
    .unwrap();
    let left = -1.0 / 3.0;
    let right = -1.0 / 7.0;
    let mut positions = Vec::new();
    let mut weights = Vec::new();
    let mut cycle = 1.0;
    for m in 0..32 {
        let shift = 2.0 * h * m as f64;
        for (depth, gain, odd) in [
            (-s - shift, left * cycle, true),
            (2.0 * h - s + shift, right * cycle, true),
            (s - 2.0 * h - shift, left * right * cycle, false),
            (s + 2.0 * h + shift, left * right * cycle, false),
        ] {
            positions.push(core::array::from_fn(|i| point[i] + normal[i] * depth));
            weights.push((gain, odd));
        }
        cycle *= left * right;
    }
    let images = PreparedCompactSourceObserver::new(
        layer.layer_state(),
        700.0,
        observer,
        &positions,
        positions.len(),
    )
    .unwrap();
    for column in 0..3 {
        let force = unit(column).map(|v| v * C64::new(0.4, -0.3));
        let normal_force = force
            .iter()
            .zip(normal)
            .fold(C64::ZERO, |sum, (v, n)| sum + v.scale(n));
        let amplitudes: Vec<_> = weights
            .iter()
            .map(|&(gain, odd)| {
                CompactSourcePhasor::force_dipole(core::array::from_fn(|i| {
                    (force[i]
                        - if odd {
                            normal_force.scale(2.0 * normal[i])
                        } else {
                            C64::ZERO
                        })
                    .scale(gain)
                }))
            })
            .collect();
        let expected = images.evaluate_harmonic_line(&amplitudes).unwrap();
        let actual = prepared.evaluate_reflected(force).unwrap();
        close(actual.pressure_pa, expected.pressure_pa, 5e-8);
        for axis in 0..3 {
            close(
                actual.particle_velocity_m_s[axis],
                expected.particle_velocity_m_s[axis],
                5e-8,
            );
        }
        let line = prepared.evaluate_reflected_harmonic_line(force).unwrap();
        assert_eq!(line.pressure_pa, actual.pressure_pa);
        assert_eq!(line.particle_velocity_m_s, actual.particle_velocity_m_s);
        assert_eq!(line.frequency_hz, 700.0);
        assert_eq!(actual.scope, FULL_WAVE_LAYER_INTERNAL_REFLECTED_FORCE_SCOPE);
    }
}

#[test]
fn original_source_derivatives_receiver_momentum_and_cartesian_reciprocity() {
    let layer = slab([300.0, 600.0, 400.0], 0.018);
    let source = [0.0, 0.0, 0.006];
    let observer = [0.04, 0.013, 0.011];
    let prepared = PreparedLayerInternalFullWaveForce::new(
        layer,
        700.0,
        source,
        observer,
        limits(),
        &mut || false,
    )
    .unwrap();
    let reverse = PreparedLayerInternalFullWaveForce::new(
        layer,
        700.0,
        observer,
        source,
        limits(),
        &mut || false,
    )
    .unwrap();
    let step = 2e-6;
    let source_factor = C64::new(0.0, 1.0 / (core::f64::consts::TAU * 700.0 * 2.4));
    for axis in 0..3 {
        let mut plus = source;
        let mut minus = source;
        plus[axis] += step;
        minus[axis] -= step;
        let a = PreparedLayerInternalFullWaveMonopole::new(
            layer,
            700.0,
            plus,
            observer,
            q_limits(),
            &mut || false,
        )
        .unwrap()
        .evaluate_reflected_harmonic_line(C64::ONE)
        .unwrap();
        let b = PreparedLayerInternalFullWaveMonopole::new(
            layer,
            700.0,
            minus,
            observer,
            q_limits(),
            &mut || false,
        )
        .unwrap()
        .evaluate_reflected_harmonic_line(C64::ONE)
        .unwrap();
        let actual = prepared
            .evaluate_reflected_harmonic_line(unit(axis))
            .unwrap();
        close(
            actual.pressure_pa,
            (a.pressure_pa - b.pressure_pa).scale(0.5 / step) * source_factor,
            3e-5,
        );
        for row in 0..3 {
            close(
                actual.particle_velocity_m_s[row],
                (a.particle_velocity_m_s[row] - b.particle_velocity_m_s[row]).scale(0.5 / step)
                    * source_factor,
                3e-5,
            );
            let swapped = reverse.evaluate_reflected_harmonic_line(unit(row)).unwrap();
            close(
                actual.particle_velocity_m_s[row],
                swapped.particle_velocity_m_s[axis],
                5e-8,
            );
        }
    }
    let force = [
        C64::new(0.2, -0.1),
        C64::new(-0.15, 0.05),
        C64::new(0.3, 0.2),
    ];
    let actual = prepared.evaluate_reflected_harmonic_line(force).unwrap();
    let direct =
        PreparedCompactSourceObserver::new(layer.layer_state(), 700.0, observer, &[source], 1)
            .unwrap()
            .evaluate_harmonic_line(&[CompactSourcePhasor::force_dipole(force)])
            .unwrap();
    let total = prepared.evaluate_total_harmonic_line(force).unwrap();
    close(
        total.pressure_pa,
        direct.pressure_pa + actual.pressure_pa,
        2e-14,
    );
    for axis in 0..3 {
        close(
            total.particle_velocity_m_s[axis],
            direct.particle_velocity_m_s[axis] + actual.particle_velocity_m_s[axis],
            2e-14,
        );
    }
    for axis in 0..3 {
        let mut plus = observer;
        let mut minus = observer;
        plus[axis] += step;
        minus[axis] -= step;
        let a = PreparedLayerInternalFullWaveForce::new(
            layer,
            700.0,
            source,
            plus,
            limits(),
            &mut || false,
        )
        .unwrap()
        .evaluate_reflected_harmonic_line(force)
        .unwrap();
        let b = PreparedLayerInternalFullWaveForce::new(
            layer,
            700.0,
            source,
            minus,
            limits(),
            &mut || false,
        )
        .unwrap()
        .evaluate_reflected_harmonic_line(force)
        .unwrap();
        close(
            actual.particle_velocity_m_s[axis],
            (a.pressure_pa - b.pressure_pa).scale(0.5 / step) * -source_factor,
            1e-4,
        );
    }
}

#[test]
fn invisible_faces_self_fields_and_bounded_force_admission() {
    let state = fluid(2.4, 350.0);
    let other = fluid(1.2, 350.0);
    let source = [0.0, 0.0, 0.025];
    let observer = [0.04, 0.0, 0.04];
    for (states, image_z, gain, exact) in [
        ([state; 3], -0.025, 0.0, true),
        ([other, state, state], -0.025, -1.0 / 3.0, true),
        ([state, state, other], 0.135, -1.0 / 3.0, false),
    ] {
        let layer = PlanarFluidLayer::new(
            states[0],
            states[1],
            states[2],
            [0.0; 3],
            [0.0, 0.0, 1.0],
            0.08,
        )
        .unwrap();
        let prepared = PreparedLayerInternalFullWaveForce::new(
            layer,
            700.0,
            source,
            observer,
            FullWaveForceLimits {
                maximum_evaluations: if exact {
                    1
                } else {
                    limits().maximum_evaluations
                },
                ..limits()
            },
            &mut || false,
        )
        .unwrap();
        let images = PreparedCompactSourceObserver::new(
            state,
            700.0,
            observer,
            &[source, [0.0, 0.0, image_z]],
            2,
        )
        .unwrap();
        for axis in 0..3 {
            let force = unit(axis);
            let image = core::array::from_fn(|i| force[i].scale(if i == 2 { -gain } else { gain }));
            let expected = images
                .evaluate_harmonic_line(&[
                    CompactSourcePhasor::force_dipole(force),
                    CompactSourcePhasor::force_dipole(image),
                ])
                .unwrap();
            let actual = prepared.evaluate_total_harmonic_line(force).unwrap();
            let tolerance = if exact { 3e-13 } else { 5e-8 };
            close(actual.pressure_pa, expected.pressure_pa, tolerance);
            for row in 0..3 {
                close(
                    actual.particle_velocity_m_s[row],
                    expected.particle_velocity_m_s[row],
                    tolerance,
                );
            }
        }
        assert_eq!(prepared.diagnostics().evaluations == 0, exact);
        assert_eq!(prepared.diagnostics().reduced_to_single_interface, exact);
        assert!(prepared.total_field_available());
    }
    let layer = slab([300.0, 600.0, 400.0], 0.018);
    let source = [0.0, 0.0, 0.006];
    let self_field = PreparedLayerInternalFullWaveForce::new(
        layer,
        700.0,
        source,
        source,
        limits(),
        &mut || false,
    )
    .unwrap();
    let transverse = self_field.evaluate_reflected(unit(0)).unwrap();
    assert_eq!(transverse.pressure_pa, C64::ZERO);
    assert!(transverse.particle_velocity_m_s[0].abs() > 0.0);
    assert!(!self_field.total_field_available());
    assert!(self_field.evaluate_total(unit(0)).is_err());
    assert!(
        self_field
            .evaluate_reflected([C64::new(f64::INFINITY, 0.0); 3])
            .is_err()
    );
    assert!(
        PreparedLayerInternalFullWaveForce::new(
            layer,
            700.0,
            [0.0; 3],
            source,
            limits(),
            &mut || false
        )
        .is_err()
    );
    assert!(
        PreparedLayerInternalFullWaveForce::new(
            layer,
            700.0,
            source,
            source,
            FullWaveForceLimits {
                maximum_evaluations: 1,
                ..limits()
            },
            &mut || false
        )
        .is_err()
    );
    // The budget covers both returning families, not a fresh cap per face.
    let work = self_field.diagnostics().evaluations;
    assert!(work > 1);
    assert!(
        PreparedLayerInternalFullWaveForce::new(
            layer,
            700.0,
            source,
            source,
            FullWaveForceLimits {
                maximum_evaluations: work - 1,
                ..limits()
            },
            &mut || false,
        )
        .is_err()
    );
    assert!(matches!(
        PreparedLayerInternalFullWaveForce::new(
            layer,
            700.0,
            source,
            source,
            limits(),
            &mut || true
        ),
        Err(AeroacError::Cancelled { .. })
    ));
}

#[test]
fn invisible_second_face_uses_original_unequal_speed_reflected_force() {
    let inside = fluid(2.4, 600.0);
    let outside = fluid(1.2, 300.0);
    let point = [0.2, -0.3, 0.1];
    let normal = [0.0, 0.6, 0.8];
    let source = core::array::from_fn(|i| point[i] + normal[i] * 0.006);
    let observer =
        core::array::from_fn(|i| point[i] + normal[i] * 0.011 + if i == 0 { 0.04 } else { 0.0 });
    let layer = PlanarFluidLayer::new(outside, inside, inside, point, normal, 0.018).unwrap();
    let interface = PlanarFluidInterface::new(
        inside,
        outside,
        layer.first_face_point_m(),
        layer.receiving_normal().map(|v| -v),
    )
    .unwrap();
    let prepared = PreparedLayerInternalFullWaveForce::new(
        layer,
        700.0,
        source,
        observer,
        limits(),
        &mut || false,
    )
    .unwrap();
    let original = PreparedReflectedFullWaveForce::new(
        interface,
        700.0,
        source,
        observer,
        limits(),
        &mut || false,
    )
    .unwrap();
    assert!(prepared.diagnostics().reduced_to_single_interface);
    assert!(prepared.diagnostics().evaluations > 0);
    assert_eq!(
        prepared.diagnostics().evaluations,
        original.diagnostics().evaluations
    );
    assert_eq!(
        prepared.diagnostics().pressure_tail_bound_pa_per_n,
        original.diagnostics().pressure_tail_bound_pa_per_n,
    );
    for axis in 0..3 {
        let actual = prepared
            .evaluate_reflected_harmonic_line(unit(axis))
            .unwrap();
        let expected = original
            .evaluate_reflected_harmonic_line(unit(axis))
            .unwrap();
        assert_eq!(actual.pressure_pa, expected.pressure_pa);
        assert_eq!(actual.particle_velocity_m_s, expected.particle_velocity_m_s);
    }
}
