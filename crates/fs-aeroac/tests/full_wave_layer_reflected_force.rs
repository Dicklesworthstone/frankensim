//! G1/G3: original force derivatives through every coherent fluid-layer return.

use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::refraction::full_wave::layer::PlanarFluidLayer;
use fs_aeroac::curle3d::refraction::full_wave::layer::reflection::{
    PreparedLayeredReflectedFullWaveForce, PreparedLayeredReflectedFullWaveMonopole,
};
use fs_aeroac::curle3d::refraction::full_wave::{FullWaveForceLimits, FullWaveMonopoleLimits};
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
fn slab(states: [FluidAcousticState; 3], h: f64) -> PlanarFluidLayer {
    PlanarFluidLayer::new(
        states[0],
        states[1],
        states[2],
        [0.0; 3],
        [0.0, 0.0, 1.0],
        h,
    )
    .unwrap()
}
fn prepare(
    layer: PlanarFluidLayer,
    hz: f64,
    source: [f64; 3],
    observer: [f64; 3],
) -> PreparedLayeredReflectedFullWaveForce {
    PreparedLayeredReflectedFullWaveForce::new(layer, hz, source, observer, limits(), &mut || false)
        .unwrap()
}
fn close(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-12),
        "actual={actual:?}, expected={expected:?}, tolerance={tolerance:e}"
    );
}
fn unit(column: usize) -> [C64; 3] {
    core::array::from_fn(|i| if i == column { C64::ONE } else { C64::ZERO })
}

#[test]
fn all_seven_force_channels_match_independently_refined_direct_spectra() {
    // Independent direct-q tanh-sinh at 35/50 digits (mpmath 1.3.0), without
    // production branch maps, normalized denominator or quadrature samples.
    // Refinement differences: force pressure <5e-19, velocity <3e-21 absolute.
    let cases = [
        (
            [fluid(1.2, 300.0), fluid(2.4, 600.0), fluid(1.8, 400.0)],
            700.0,
            0.025,
            0.018,
            0.023,
            0.04,
            [
                C64::new(4.750322168227399, 0.08008661562454795),
                C64::new(5.521433239708101, -1.2992113722244687),
                C64::new(0.0004896283477262852, 0.003851206008633753),
                C64::new(0.0003793506176234353, -0.02250110876048811),
                C64::new(0.001311435651648497, 0.02855944754777732),
                C64::new(-0.001311435651648497, -0.02855944754777732),
                C64::new(-0.0073569746036248525, -0.005570689064708436),
            ],
        ),
        (
            [fluid(2.4, 600.0), fluid(1.4, 450.0), fluid(1.2, 300.0)],
            1500.0,
            0.04,
            0.07,
            0.06,
            0.12,
            [
                C64::new(-0.7431756243168111, -1.6834602495061538),
                C64::new(-0.16135936449573414, -1.6588071640414331),
                C64::new(-0.0003308620229371325, -0.0008288356613157551),
                C64::new(-0.0006202106949250459, 0.0002737964680449673),
                C64::new(0.00036585614653054185, -0.0010580087978838389),
                C64::new(-0.00036585614653054185, 0.0010580087978838389),
                C64::new(-0.00013827600299375264, 0.0006627538627258825),
            ],
        ),
    ];
    for (media, hz, a, h, b, r, expected) in cases {
        let observer = prepare(slab(media, h), hz, [0.0, 0.0, -a], [r, 0.0, -b]);
        let x = observer.evaluate_reflected(unit(0)).unwrap();
        let y = observer.evaluate_reflected(unit(1)).unwrap();
        let z = observer.evaluate_reflected(unit(2)).unwrap();
        let actual = [
            x.pressure_pa,
            z.pressure_pa,
            x.particle_velocity_m_s[0],
            y.particle_velocity_m_s[1],
            z.particle_velocity_m_s[0],
            x.particle_velocity_m_s[2],
            z.particle_velocity_m_s[2],
        ];
        for (a, e) in actual.into_iter().zip(expected) {
            close(a, e, 4e-8);
        }
        assert_eq!(y.pressure_pa, C64::ZERO);
        assert!(observer.diagnostics().evaluations > 0);
        assert!(!observer.diagnostics().reduced_to_single_interface);
        assert!(observer.diagnostics().evaluations <= limits().maximum_evaluations);
        let line = observer.evaluate_reflected_harmonic_line(unit(2)).unwrap();
        assert_eq!(line.frequency_hz, hz);
        assert_eq!(line.pressure_pa, z.pressure_pa);
        assert_eq!(line.particle_velocity_m_s, z.particle_velocity_m_s);
    }
}

#[test]
fn rotated_equal_speed_image_forces_keep_normal_reversal_and_total_intensity() {
    let (r1, r2, r3, c, h, hz) = (1.2, 5.0, 2.4, 350.0, 0.08, 1100.0);
    let media = [fluid(r1, c), fluid(r2, c), fluid(r3, c)];
    let layer =
        PlanarFluidLayer::new(media[0], media[1], media[2], [0.0; 3], [2.0, 1.0, 2.0], h).unwrap();
    let n = layer.receiving_normal();
    let t = [-1.0 / 5.0_f64.sqrt(), 2.0 / 5.0_f64.sqrt(), 0.0];
    let (a, b, r) = (0.03, 0.04, 0.06);
    let source = n.map(|x| -a * x);
    let receiver = core::array::from_fn(|i| r * t[i] - b * n[i]);
    let force = [
        C64::new(0.2, -0.1),
        C64::new(-0.3, 0.2),
        C64::new(0.1, 0.05),
    ];
    let dot = force[0].scale(n[0]) + force[1].scale(n[1]) + force[2].scale(n[2]);
    let image_force: [C64; 3] = core::array::from_fn(|i| force[i] - dot.scale(2.0 * n[i]));
    let observer = prepare(layer, hz, source, receiver);
    let r12 = (r2 - r1) / (r2 + r1);
    let r23 = (r3 - r2) / (r3 + r2);
    let mut coefficient = (1.0 - r12 * r12) * r23;
    let mut positions = Vec::new();
    let mut sources = Vec::new();
    for m in 0..80 {
        let scale = if m == 0 { r12 } else { coefficient };
        positions.push(n.map(|x| (a + 2.0 * (m as f64) * h) * x));
        sources.push(CompactSourcePhasor::force_dipole(
            image_force.map(|f| f.scale(scale)),
        ));
        if m > 0 {
            coefficient *= -r12 * r23;
        }
    }
    assert!(coefficient.abs() < 1e-50);
    for total in [false, true] {
        if total {
            positions.push(source);
            sources.push(CompactSourcePhasor::force_dipole(force));
        }
        let expected =
            PreparedCompactSourceObserver::new(media[0], hz, receiver, &positions, positions.len())
                .unwrap()
                .evaluate(&sources)
                .unwrap();
        let actual = if total {
            observer.evaluate_total(force)
        } else {
            observer.evaluate_reflected(force)
        }
        .unwrap();
        close(actual.pressure_pa, expected.pressure_pa, 4e-8);
        for axis in 0..3 {
            close(
                actual.particle_velocity_m_s[axis],
                expected.particle_velocity_m_s[axis],
                4e-8,
            );
            close(
                C64::new(
                    actual.active_intensity_w_m2[axis],
                    actual.reactive_intensity_w_m2[axis],
                ),
                C64::new(
                    expected.active_intensity_w_m2[axis],
                    expected.reactive_intensity_w_m2[axis],
                ),
                8e-8,
            );
        }
    }
    assert!(observer.has_total_field());
    assert_eq!(observer.frequency_hz(), hz);
    assert_eq!(observer.source_position_m(), source);
    assert_eq!(observer.observer_position_m(), receiver);
    assert_eq!(observer.layer().thickness_m(), h);
}

#[test]
fn force_is_the_original_source_derivative_of_the_reflected_monopole() {
    let media = [fluid(1.2, 300.0), fluid(2.4, 600.0), fluid(1.8, 400.0)];
    let layer = slab(media, 0.018);
    let hz = 700.0;
    let source = [0.0, 0.0, -0.025];
    let receiver = [0.04, 0.0, -0.023];
    let prepared = prepare(layer, hz, source, receiver);
    let monopole_limits = FullWaveMonopoleLimits {
        absolute_pressure_tolerance_pa_per_m3_s: 1e-10,
        absolute_velocity_tolerance_per_m2: 1e-13,
        relative_tolerance: 1e-10,
        maximum_evaluations: 262_144,
        maximum_depth: 40,
    };
    for axis in [0, 2] {
        let step = 1e-6;
        let mut plus = source;
        let mut minus = source;
        plus[axis] += step;
        minus[axis] -= step;
        let monopole = |s| {
            PreparedLayeredReflectedFullWaveMonopole::new(
                layer,
                hz,
                s,
                receiver,
                monopole_limits,
                &mut || false,
            )
            .unwrap()
            .evaluate_reflected(C64::ONE)
            .unwrap()
        };
        let plus = monopole(plus);
        let minus = monopole(minus);
        let scale = C64::new(
            0.0,
            1.0 / (2.0 * step * TAU * hz * media[0].density_kg_m3()),
        );
        let expected_pressure = (plus.pressure_pa - minus.pressure_pa) * scale;
        let expected_velocity: [C64; 3] = core::array::from_fn(|i| {
            (plus.particle_velocity_m_s[i] - minus.particle_velocity_m_s[i]) * scale
        });
        let actual = prepared.evaluate_reflected(unit(axis)).unwrap();
        close(actual.pressure_pa, expected_pressure, 3e-6);
        for i in 0..3 {
            close(actual.particle_velocity_m_s[i], expected_velocity[i], 3e-6);
        }
    }
}

#[test]
fn exact_reductions_axis_motion_self_fields_and_work_refusals_are_explicit() {
    let media = [fluid(1.2, 300.0), fluid(2.4, 600.0), fluid(1.8, 400.0)];
    let layer = slab(media, 0.018);
    let source = [0.0, 0.0, -0.025];
    let receiver = [0.0, 0.0, -0.05];
    let hz = 700.0;
    let axis = prepare(layer, hz, source, receiver);
    let transverse = axis.evaluate_reflected(unit(0)).unwrap();
    assert_eq!(transverse.pressure_pa, C64::ZERO);
    assert!(transverse.particle_velocity_m_s[0].abs() > 0.0);
    assert_eq!(transverse.particle_velocity_m_s[1], C64::ZERO);
    assert_eq!(transverse.particle_velocity_m_s[2], C64::ZERO);
    let self_field = prepare(layer, hz, source, source);
    assert!(!self_field.has_total_field());
    assert!(
        self_field
            .evaluate_reflected(unit(2))
            .unwrap()
            .pressure_pa
            .abs()
            .is_finite()
    );
    assert!(self_field.evaluate_total([C64::ZERO; 3]).is_err());
    assert!(self_field.evaluate_total_harmonic_line(unit(0)).is_err());
    assert!(
        axis.evaluate_reflected([C64::new(f64::NAN, 0.0); 3])
            .is_err()
    );
    assert_eq!(
        axis.evaluate_reflected([C64::ZERO; 3])
            .unwrap()
            .particle_velocity_m_s,
        [C64::ZERO; 3]
    );
    let line = axis.evaluate_total_harmonic_line(unit(2)).unwrap();
    assert_eq!(line.frequency_hz, hz);
    assert_eq!(
        line.pressure_pa,
        axis.evaluate_total(unit(2)).unwrap().pressure_pa
    );
    let same = slab([media[0]; 3], 0.02);
    let homogeneous = prepare(same, hz, source, receiver);
    assert!(homogeneous.diagnostics().reduced_to_single_interface);
    assert_eq!(homogeneous.diagnostics().evaluations, 0);
    let force = [
        C64::new(0.2, 0.1),
        C64::new(-0.1, 0.03),
        C64::new(0.1, -0.2),
    ];
    assert_eq!(
        homogeneous.evaluate_reflected(force).unwrap().pressure_pa,
        C64::ZERO
    );
    let direct = PreparedCompactSourceObserver::new(media[0], hz, receiver, &[source], 1)
        .unwrap()
        .evaluate(&[CompactSourcePhasor::force_dipole(force)])
        .unwrap();
    let actual = homogeneous.evaluate_total(force).unwrap();
    close(actual.pressure_pa, direct.pressure_pa, 2e-13);
    for i in 0..3 {
        close(
            actual.particle_velocity_m_s[i],
            direct.particle_velocity_m_s[i],
            2e-13,
        );
    }
    let mut budget = limits();
    budget.maximum_evaluations = 1;
    assert!(
        PreparedLayeredReflectedFullWaveForce::new(
            layer,
            hz,
            source,
            receiver,
            budget,
            &mut || false
        )
        .is_err()
    );
    assert!(matches!(
        PreparedLayeredReflectedFullWaveForce::new(
            layer,
            hz,
            source,
            receiver,
            limits(),
            &mut || true
        ),
        Err(AeroacError::Cancelled { .. })
    ));
    assert!(
        PreparedLayeredReflectedFullWaveForce::new(
            layer,
            0.0,
            source,
            receiver,
            limits(),
            &mut || false
        )
        .is_err()
    );
    for p in [[0.0; 3], [0.0, 0.0, 0.02], [0.0, 0.0, f64::NAN]] {
        assert!(
            PreparedLayeredReflectedFullWaveForce::new(layer, hz, source, p, limits(), &mut || {
                false
            })
            .is_err()
        );
    }
}
