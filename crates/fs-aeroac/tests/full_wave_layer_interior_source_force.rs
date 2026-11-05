//! G1/G3: physical force sources in the layer and actual exterior momentum.
use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::refraction::PlanarFluidInterface;
use fs_aeroac::curle3d::refraction::full_wave::layer::PlanarFluidLayer;
use fs_aeroac::curle3d::refraction::full_wave::layer::interior::source::{
    PreparedInteriorSourceFullWaveForce, PreparedInteriorSourceFullWaveMonopole,
};
use fs_aeroac::curle3d::refraction::full_wave::layer::interior::{
    PreparedInteriorFullWaveForce, PreparedLayerInteriorFullWaveMonopole,
};
use fs_aeroac::curle3d::refraction::full_wave::{
    FullWaveForceLimits, FullWaveMonopoleLimits, PreparedFullWaveForce,
};
use fs_aeroac::curle3d::{CompactSourcePhasor, PreparedCompactSourceObserver};
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::c64::C64;
use std::f64::consts::TAU;

fn fluid(rho: f64, c: f64) -> FluidAcousticState {
    FluidAcousticState::try_new(rho, c).unwrap()
}
fn limits() -> FullWaveForceLimits {
    FullWaveForceLimits {
        absolute_pressure_tolerance_pa_per_n: 1e-12,
        absolute_velocity_tolerance_m_s_per_n: 1e-15,
        relative_tolerance: 1e-10,
        maximum_evaluations: 262_144,
        maximum_depth: 40,
    }
}
fn q_limits() -> FullWaveMonopoleLimits {
    FullWaveMonopoleLimits {
        absolute_pressure_tolerance_pa_per_m3_s: 1e-11,
        absolute_velocity_tolerance_per_m2: 1e-14,
        relative_tolerance: 1e-11,
        maximum_evaluations: 262_144,
        maximum_depth: 40,
    }
}
fn slab(m: [FluidAcousticState; 3], h: f64) -> PlanarFluidLayer {
    PlanarFluidLayer::new(m[0], m[1], m[2], [0.0; 3], [0.0, 0.0, 1.0], h).unwrap()
}
fn prepare(
    l: PlanarFluidLayer,
    hz: f64,
    source: [f64; 3],
    observer: [f64; 3],
) -> PreparedInteriorSourceFullWaveForce {
    PreparedInteriorSourceFullWaveForce::new(l, hz, source, observer, limits(), &mut || false)
        .unwrap()
}
fn unit(axis: usize) -> [C64; 3] {
    core::array::from_fn(|i| if i == axis { C64::ONE } else { C64::ZERO })
}
fn close(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-12),
        "actual={actual:?}, expected={expected:?}, tolerance={tolerance:e}"
    );
}

#[test]
fn seven_original_source_derivatives_match_direct_spectra_and_cartesian_reciprocity() {
    // Independent direct-q tanh-sinh references at 35/50 digits (mpmath 1.3.0),
    // with the analytic interior-source and exterior-receiver derivatives of
    // the same reciprocal Green function. No production branch maps or nodes.
    // Reciprocal conversion changes signs/basis, not their refinement errors:
    // pressure <2e-35, velocity <3e-37 absolute between precision levels.
    let cases = [
        (
            [fluid(1.2, 300.0), fluid(2.4, 600.0), fluid(1.8, 400.0)],
            700.0,
            0.025,
            0.018,
            0.011,
            0.04,
            [
                C64::new(15.331593181211314, 1.4199767956941838),
                C64::new(-17.354082268289794, -3.750606775861318),
                C64::new(0.006326400549341643, 0.03755890072184411),
                C64::new(0.006726081134241659, -0.07262198929356481),
                C64::new(0.0019931314747980494, -0.1194211535381207),
                C64::new(-0.00034905601881881075, -0.09542008046462222),
                C64::new(0.007106339111729671, 0.024017404551634758),
            ],
        ),
        (
            [fluid(2.4, 600.0), fluid(1.4, 450.0), fluid(1.2, 300.0)],
            1500.0,
            0.04,
            0.07,
            0.04,
            0.12,
            [
                C64::new(1.6304738743883518, 8.881854047704289),
                C64::new(-1.8869283199528588, -11.525633575252407),
                C64::new(-0.0016937780012448163, 0.005540770310457754),
                C64::new(0.0032722013322057787, -0.0006006897608590283),
                C64::new(0.0036248095727295176, -0.007213870063770369),
                C64::new(0.003271518095074575, -0.002605107948677715),
                C64::new(0.0005220121609890306, 0.0037403702890234365),
            ],
        ),
    ];
    for (m, hz, a, h, s, r, reference) in cases {
        let l = slab(m, h);
        let source = [0.0, 0.0, s];
        let receiver = [r, 0.0, -a];
        let prepared = prepare(l, hz, source, receiver);
        let fields: [_; 3] = core::array::from_fn(|i| prepared.evaluate(unit(i)).unwrap());
        let actual = [
            fields[0].pressure_pa,
            fields[2].pressure_pa,
            fields[0].particle_velocity_m_s[0],
            fields[1].particle_velocity_m_s[1],
            fields[2].particle_velocity_m_s[0],
            fields[0].particle_velocity_m_s[2],
            fields[2].particle_velocity_m_s[2],
        ];
        for (value, expected) in actual.into_iter().zip(reference) {
            close(value, expected, 4e-8);
        }
        assert_eq!(fields[1].pressure_pa, C64::ZERO);
        assert!((actual[4] - actual[5]).abs() > 1e-3);
        let reciprocal_q = PreparedLayerInteriorFullWaveMonopole::new(
            l,
            hz,
            receiver,
            source,
            q_limits(),
            &mut || false,
        )
        .unwrap()
        .evaluate(C64::ONE)
        .unwrap();
        let reciprocal_f =
            PreparedInteriorFullWaveForce::new(l, hz, receiver, source, limits(), &mut || false)
                .unwrap();
        for column in 0..3 {
            close(
                fields[column].pressure_pa,
                -reciprocal_q.particle_velocity_m_s[column],
                4e-8,
            );
            for row in 0..3 {
                close(
                    fields[column].particle_velocity_m_s[row],
                    reciprocal_f
                        .evaluate(unit(row))
                        .unwrap()
                        .particle_velocity_m_s[column],
                    4e-8,
                );
            }
        }
        let line = prepared.evaluate_harmonic_line(unit(2)).unwrap();
        assert_eq!(line.frequency_hz, hz);
        assert_eq!(line.pressure_pa, fields[2].pressure_pa);
        assert_eq!(line.particle_velocity_m_s, fields[2].particle_velocity_m_s);
        assert!(!prepared.diagnostics().reduced_to_single_interface);
        assert!(prepared.diagnostics().evaluations > 0);
        assert!(prepared.diagnostics().evaluations <= limits().maximum_evaluations);
    }
}

#[test]
fn rotated_internal_force_images_radiate_into_either_original_exterior() {
    let m = [fluid(1.2, 350.0), fluid(5.0, 350.0), fluid(2.4, 350.0)];
    let (h, s, a, r, hz) = (0.08, 0.031, 0.03, 0.06, 1100.0);
    let origin = [0.13, -0.2, 0.07];
    let original = PlanarFluidLayer::new(m[0], m[1], m[2], origin, [2.0, 1.0, 2.0], h).unwrap();
    let normal = original.receiving_normal();
    let tangent = [-1.0 / 5.0_f64.sqrt(), 2.0 / 5.0_f64.sqrt(), 0.0];
    let source = core::array::from_fn(|i| origin[i] + s * normal[i]);
    let force = [
        C64::new(0.2, -0.1),
        C64::new(-0.3, 0.2),
        C64::new(0.1, 0.05),
    ];
    for opposite in [false, true] {
        let layer = if opposite {
            PlanarFluidLayer::new(
                m[2],
                m[1],
                m[0],
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
        let normal_force = force[0].scale(n[0]) + force[1].scale(n[1]) + force[2].scale(n[2]);
        let reflected: [C64; 3] =
            core::array::from_fn(|i| force[i] - normal_force.scale(2.0 * n[i]));
        let mut gain = 2.0 * rho1 / (rho1 + rho2);
        let mut positions = Vec::new();
        let mut amplitudes = Vec::new();
        for m in 0..80 {
            let shift = 2.0 * f64::from(m) * h;
            positions.push(core::array::from_fn(|i| source[i] + shift * n[i]));
            amplitudes.push(CompactSourcePhasor::force_dipole(
                force.map(|f| f.scale(gain)),
            ));
            positions.push(core::array::from_fn(|i| {
                source[i] + (2.0 * (h - depth) + shift) * n[i]
            }));
            amplitudes.push(CompactSourcePhasor::force_dipole(
                reflected.map(|f| f.scale(gain * r23)),
            ));
            gain *= r21 * r23;
        }
        assert!(gain.abs() < 1e-50);
        let expected = PreparedCompactSourceObserver::new(
            layer.incident_state(),
            hz,
            receiver,
            &positions,
            positions.len(),
        )
        .unwrap()
        .evaluate(&amplitudes)
        .unwrap();
        let actual = prepared.evaluate(force).unwrap();
        close(actual.pressure_pa, expected.pressure_pa, 4e-8);
        for i in 0..3 {
            close(
                actual.particle_velocity_m_s[i],
                expected.particle_velocity_m_s[i],
                4e-8,
            );
            close(
                C64::new(
                    actual.active_intensity_w_m2[i],
                    actual.reactive_intensity_w_m2[i],
                ),
                C64::new(
                    expected.active_intensity_w_m2[i],
                    expected.reactive_intensity_w_m2[i],
                ),
                8e-8,
            );
        }
        assert_eq!(prepared.source_position_m(), source);
        assert_eq!(prepared.observer_position_m(), receiver);
        assert_eq!(prepared.frequency_hz(), hz);
        assert_eq!(prepared.layer().layer_state(), layer.layer_state());
    }
}

#[test]
fn original_source_derivatives_and_exterior_momentum_use_their_actual_densities() {
    let m = [fluid(1.2, 300.0), fluid(2.4, 600.0), fluid(1.8, 400.0)];
    let l = slab(m, 0.018);
    let hz = 700.0;
    let source = [0.0, 0.0, 0.011];
    let receiver = [0.04, 0.015, -0.025];
    let prepared = prepare(l, hz, source, receiver);
    let step = 1e-6;
    for axis in 0..3 {
        let mut a = source;
        let mut b = source;
        a[axis] += step;
        b[axis] -= step;
        let q = |s| {
            PreparedInteriorSourceFullWaveMonopole::new(l, hz, s, receiver, q_limits(), &mut || {
                false
            })
            .unwrap()
            .evaluate_harmonic_line(C64::ONE)
            .unwrap()
        };
        let a = q(a);
        let b = q(b);
        let scale = C64::new(0.0, 1.0 / (2.0 * step * TAU * hz * m[1].density_kg_m3()));
        let actual = prepared.evaluate(unit(axis)).unwrap();
        close(
            actual.pressure_pa,
            (a.pressure_pa - b.pressure_pa) * scale,
            2e-5,
        );
        for row in 0..3 {
            close(
                actual.particle_velocity_m_s[row],
                (a.particle_velocity_m_s[row] - b.particle_velocity_m_s[row]) * scale,
                2e-5,
            );
        }
    }
    let force = [
        C64::new(0.2, -0.1),
        C64::new(-0.3, 0.2),
        C64::new(0.1, 0.05),
    ];
    let actual = prepared.evaluate(force).unwrap();
    for axis in 0..3 {
        let mut a = receiver;
        let mut b = receiver;
        a[axis] += step;
        b[axis] -= step;
        let derivative = (prepare(l, hz, source, a)
            .evaluate(force)
            .unwrap()
            .pressure_pa
            - prepare(l, hz, source, b)
                .evaluate(force)
                .unwrap()
                .pressure_pa)
            .scale(0.5 / step);
        close(
            actual.particle_velocity_m_s[axis],
            derivative * C64::new(0.0, -1.0 / (TAU * hz * m[0].density_kg_m3())),
            2e-5,
        );
    }
}

#[test]
fn exact_reductions_axial_motion_and_strict_source_domains_are_retained() {
    let first = fluid(1.2, 300.0);
    let source = [0.006, -0.003, 0.011];
    let receiver = [0.029, 0.014, -0.025];
    let hz = 700.0;
    let force = [
        C64::new(0.2, -0.1),
        C64::new(-0.3, 0.2),
        C64::new(0.1, 0.05),
    ];
    for inside in [first, fluid(3.0, 300.0), fluid(2.4, 450.0)] {
        let layer = slab([first, inside, inside], 0.018);
        let prepared = prepare(layer, hz, source, receiver);
        let original = PreparedFullWaveForce::new(
            PlanarFluidInterface::new(inside, first, [0.0; 3], [0.0, 0.0, -1.0]).unwrap(),
            hz,
            source,
            receiver,
            limits(),
            &mut || false,
        )
        .unwrap();
        let actual = prepared.evaluate(force).unwrap();
        let expected = original.evaluate(force).unwrap();
        assert_eq!(actual.pressure_pa, expected.pressure_pa);
        assert_eq!(actual.particle_velocity_m_s, expected.particle_velocity_m_s);
        assert!(prepared.diagnostics().reduced_to_single_interface);
        if first.sound_speed_m_s() == inside.sound_speed_m_s() {
            assert_eq!(prepared.diagnostics().evaluations, 0);
        }
    }
    let m = [first, fluid(2.4, 600.0), fluid(1.8, 400.0)];
    let layer = slab(m, 0.018);
    let axis = prepare(layer, hz, [0.0, 0.0, 0.011], [0.0, 0.0, -0.025]);
    let transverse = axis.evaluate(unit(0)).unwrap();
    assert_eq!(transverse.pressure_pa, C64::ZERO);
    assert!(transverse.particle_velocity_m_s[0].abs() > 0.0);
    assert_eq!(transverse.particle_velocity_m_s[2], C64::ZERO);
    close(
        axis.evaluate(unit(1)).unwrap().particle_velocity_m_s[1],
        transverse.particle_velocity_m_s[0],
        2e-13,
    );
    assert_eq!(
        axis.evaluate([C64::ZERO; 3]).unwrap().particle_velocity_m_s,
        [C64::ZERO; 3]
    );
    assert!(axis.evaluate([C64::new(f64::NAN, 0.0); 3]).is_err());
    let mut small = limits();
    small.maximum_evaluations = 1;
    assert!(
        PreparedInteriorSourceFullWaveForce::new(layer, hz, source, receiver, small, &mut || false)
            .is_err()
    );
    assert!(matches!(
        PreparedInteriorSourceFullWaveForce::new(
            layer,
            hz,
            source,
            receiver,
            limits(),
            &mut || true
        ),
        Err(AeroacError::Cancelled { .. })
    ));
    for z in [-0.01, 0.0, 0.018, 0.03, f64::NAN] {
        assert!(
            PreparedInteriorSourceFullWaveForce::new(
                layer,
                hz,
                [0.0, 0.0, z],
                receiver,
                limits(),
                &mut || false
            )
            .is_err()
        );
    }
    for z in [0.0, 0.01, 0.018, 0.03, f64::NAN] {
        assert!(
            PreparedInteriorSourceFullWaveForce::new(
                layer,
                hz,
                source,
                [0.0, 0.0, z],
                limits(),
                &mut || false
            )
            .is_err()
        );
    }
    assert!(
        PreparedInteriorSourceFullWaveForce::new(
            slab(m, 0.0),
            hz,
            source,
            receiver,
            limits(),
            &mut || false
        )
        .is_err()
    );
    assert!(
        PreparedInteriorSourceFullWaveForce::new(
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
