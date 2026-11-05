//! G1/G3: original force derivatives measured inside the actual fluid layer.

use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::refraction::PlanarFluidInterface;
use fs_aeroac::curle3d::refraction::full_wave::layer::PlanarFluidLayer;
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
) -> PreparedInteriorFullWaveForce {
    PreparedInteriorFullWaveForce::new(layer, hz, source, observer, limits(), &mut || false)
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
fn all_seven_interior_force_channels_match_independently_refined_direct_spectra() {
    // Independent direct-q tanh-sinh at 35/50 digits (mpmath 1.3.0), without
    // production branch maps, normalized weights or Gauss--Kronrod samples.
    // Refinement differences: force pressure <4e-35, velocity <3e-37 absolute.
    let cases = [
        (
            [fluid(1.2, 300.0), fluid(2.4, 600.0), fluid(1.8, 400.0)],
            700.0,
            0.025,
            0.018,
            0.011,
            0.04,
            [
                C64::new(30.66318636242263, 2.8399535913883676),
                C64::new(26.909100026390597, 0.6331377127688922),
                C64::new(0.006326400549341643, 0.03755890072184411),
                C64::new(0.006726081134241659, -0.07262198929356481),
                C64::new(0.00034905601881881075, 0.09542008046462222),
                C64::new(-0.0019931314747980494, 0.1194211535381207),
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
                C64::new(0.9511097600598719, 5.1810815278275015),
                C64::new(0.5530247980535719, 2.1472466372896144),
                C64::new(-0.0016937780012448163, 0.005540770310457754),
                C64::new(0.0032722013322057787, -0.0006006897608590283),
                C64::new(-0.003271518095074575, 0.002605107948677715),
                C64::new(-0.0036248095727295176, 0.007213870063770369),
                C64::new(0.0005220121609890306, 0.0037403702890234365),
            ],
        ),
    ];
    for (media, hz, a, h, s, r, expected) in cases {
        let observer = prepare(slab(media, h), hz, [0.0, 0.0, -a], [r, 0.0, s]);
        let x = observer.evaluate(unit(0)).unwrap();
        let y = observer.evaluate(unit(1)).unwrap();
        let z = observer.evaluate(unit(2)).unwrap();
        let actual = [
            x.pressure_pa,
            z.pressure_pa,
            x.particle_velocity_m_s[0],
            y.particle_velocity_m_s[1],
            z.particle_velocity_m_s[0],
            x.particle_velocity_m_s[2],
            z.particle_velocity_m_s[2],
        ];
        for (actual, expected) in actual.into_iter().zip(expected) {
            close(actual, expected, 4e-8);
        }
        // Source and receiving normal derivatives retain distinct internal
        // standing-wave factors; the mixed entries need not agree or oppose.
        assert!((actual[4] - actual[5]).abs() > 1e-3);
        assert_eq!(y.pressure_pa, C64::ZERO);
        assert!(y.particle_velocity_m_s[1].abs() > 0.0);
        assert!(!observer.diagnostics().reduced_to_single_interface);
        assert!(observer.diagnostics().evaluations > 0);
        assert!(observer.diagnostics().evaluations <= limits().maximum_evaluations);
    }
}

#[test]
fn rotated_internal_image_forces_retain_density_gain_and_backward_normal_reversal() {
    let (r1, r2, r3, c, h, hz) = (1.2, 5.0, 2.4, 350.0, 0.08, 1100.0);
    let media = [fluid(r1, c), fluid(r2, c), fluid(r3, c)];
    let origin = [0.13, -0.2, 0.07];
    let layer =
        PlanarFluidLayer::new(media[0], media[1], media[2], origin, [2.0, 1.0, 2.0], h).unwrap();
    let n = layer.receiving_normal();
    let t = [-1.0 / 5.0_f64.sqrt(), 2.0 / 5.0_f64.sqrt(), 0.0];
    let (a, s, r) = (0.03, 0.031, 0.06);
    let source = core::array::from_fn(|i| origin[i] - a * n[i]);
    let receiver = core::array::from_fn(|i| origin[i] + r * t[i] + s * n[i]);
    let force = [
        C64::new(0.2, -0.1),
        C64::new(-0.3, 0.2),
        C64::new(0.1, 0.05),
    ];
    let normal_force = force[0].scale(n[0]) + force[1].scale(n[1]) + force[2].scale(n[2]);
    let image_force: [C64; 3] = core::array::from_fn(|i| force[i] - normal_force.scale(2.0 * n[i]));
    let observer = prepare(layer, hz, source, receiver);
    let r12 = (r2 - r1) / (r2 + r1);
    let r23 = (r3 - r2) / (r3 + r2);
    // The force derivative contributes rho_layer/rho_incident relative to
    // the equivalent volume-source images in the original layer fluid.
    let mut gain = 2.0 * r2 / (r1 + r2);
    let mut positions = Vec::new();
    let mut amplitudes = Vec::new();
    for m in 0..80 {
        let shift = 2.0 * f64::from(m) * h;
        positions.push(core::array::from_fn(|i| origin[i] + (-a - shift) * n[i]));
        amplitudes.push(CompactSourcePhasor::force_dipole(
            force.map(|f| f.scale(gain)),
        ));
        positions.push(core::array::from_fn(|i| {
            origin[i] + (a + 2.0 * h + shift) * n[i]
        }));
        amplitudes.push(CompactSourcePhasor::force_dipole(
            image_force.map(|f| f.scale(gain * r23)),
        ));
        gain *= -r12 * r23;
    }
    assert!(gain.abs() < 1e-50);
    let expected =
        PreparedCompactSourceObserver::new(media[1], hz, receiver, &positions, positions.len())
            .unwrap()
            .evaluate(&amplitudes)
            .unwrap();
    let actual = observer.evaluate(force).unwrap();
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
    assert_eq!(observer.frequency_hz(), hz);
    assert_eq!(observer.source_position_m(), source);
    assert_eq!(observer.observer_position_m(), receiver);
    assert_eq!(observer.layer().layer_state(), media[1]);
    assert_eq!(observer.layer().thickness_m(), h);
}

#[test]
fn each_original_force_axis_is_the_source_derivative_of_the_interior_monopole() {
    let media = [fluid(1.2, 300.0), fluid(2.4, 600.0), fluid(1.8, 400.0)];
    let layer = slab(media, 0.018);
    let hz = 700.0;
    let source = [0.0, 0.0, -0.025];
    let receiver = [0.04, 0.015, 0.011];
    let observer = prepare(layer, hz, source, receiver);
    let monopole_limits = FullWaveMonopoleLimits {
        absolute_pressure_tolerance_pa_per_m3_s: 1e-11,
        absolute_velocity_tolerance_per_m2: 1e-14,
        relative_tolerance: 1e-11,
        maximum_evaluations: 262_144,
        maximum_depth: 40,
    };
    let step = 1e-6;
    for axis in 0..3 {
        let mut plus = source;
        let mut minus = source;
        plus[axis] += step;
        minus[axis] -= step;
        let monopole = |position| {
            PreparedLayerInteriorFullWaveMonopole::new(
                layer,
                hz,
                position,
                receiver,
                monopole_limits,
                &mut || false,
            )
            .unwrap()
            .evaluate_harmonic_line(C64::ONE)
            .unwrap()
        };
        let plus = monopole(plus);
        let minus = monopole(minus);
        let conversion = C64::new(
            0.0,
            1.0 / (2.0 * step * TAU * hz * media[0].density_kg_m3()),
        );
        let actual = observer.evaluate(unit(axis)).unwrap();
        close(
            actual.pressure_pa,
            (plus.pressure_pa - minus.pressure_pa) * conversion,
            2e-5,
        );
        for row in 0..3 {
            close(
                actual.particle_velocity_m_s[row],
                (plus.particle_velocity_m_s[row] - minus.particle_velocity_m_s[row]) * conversion,
                2e-5,
            );
        }
    }
}

#[test]
fn exact_reductions_axis_motion_original_clock_and_bounded_domains_are_retained() {
    let source = [-0.01, 0.02, -0.025];
    let receiver = [0.04, -0.01, 0.011];
    let hz = 700.0;
    let first = fluid(1.2, 300.0);
    let last = fluid(2.4, 450.0);
    let force = [
        C64::new(0.2, 0.1),
        C64::new(-0.1, 0.03),
        C64::new(0.1, -0.2),
    ];
    let invisible = prepare(slab([first, last, last], 0.018), hz, source, receiver);
    let interface = PlanarFluidInterface::new(first, last, [0.0; 3], [0.0, 0.0, 1.0]).unwrap();
    let original =
        PreparedFullWaveForce::new(interface, hz, source, receiver, limits(), &mut || false)
            .unwrap();
    let actual = invisible.evaluate(force).unwrap();
    let expected = original.evaluate(force).unwrap();
    assert_eq!(actual.pressure_pa, expected.pressure_pa);
    assert_eq!(actual.particle_velocity_m_s, expected.particle_velocity_m_s);
    assert_eq!(
        invisible.diagnostics().evaluations,
        original.diagnostics().evaluations
    );
    assert!(invisible.diagnostics().reduced_to_single_interface);
    let mut tiny = limits();
    tiny.maximum_evaluations = 1;
    let homogeneous = PreparedInteriorFullWaveForce::new(
        slab([first; 3], 0.018),
        hz,
        source,
        receiver,
        tiny,
        &mut || false,
    )
    .unwrap();
    let free = PreparedCompactSourceObserver::new(first, hz, receiver, &[source], 1)
        .unwrap()
        .evaluate(&[CompactSourcePhasor::force_dipole(force)])
        .unwrap();
    let actual = homogeneous.evaluate(force).unwrap();
    close(actual.pressure_pa, free.pressure_pa, 2e-13);
    for i in 0..3 {
        close(
            actual.particle_velocity_m_s[i],
            free.particle_velocity_m_s[i],
            2e-13,
        );
    }
    assert_eq!(homogeneous.diagnostics().evaluations, 0);
    let media = [first, fluid(2.4, 600.0), fluid(1.8, 400.0)];
    let layer = slab(media, 0.018);
    let axis_source = [0.0, 0.0, -0.025];
    let axis_receiver = [0.0, 0.0, 0.011];
    let axis = prepare(layer, hz, axis_source, axis_receiver);
    let x = axis.evaluate(unit(0)).unwrap();
    assert_eq!(x.pressure_pa, C64::ZERO);
    assert!(x.particle_velocity_m_s[0].abs() > 0.0);
    assert_eq!(x.particle_velocity_m_s[1], C64::ZERO);
    assert_eq!(x.particle_velocity_m_s[2], C64::ZERO);
    close(
        axis.evaluate(unit(1)).unwrap().particle_velocity_m_s[1],
        x.particle_velocity_m_s[0],
        2e-13,
    );
    let line = axis.evaluate_harmonic_line(force).unwrap();
    let field = axis.evaluate(force).unwrap();
    assert_eq!(line.frequency_hz, hz);
    assert_eq!(line.pressure_pa, field.pressure_pa);
    assert_eq!(line.particle_velocity_m_s, field.particle_velocity_m_s);
    assert_eq!(
        axis.evaluate([C64::ZERO; 3]).unwrap().particle_velocity_m_s,
        [C64::ZERO; 3]
    );
    assert!(axis.evaluate([C64::new(f64::NAN, 0.0); 3]).is_err());
    assert!(
        PreparedInteriorFullWaveForce::new(
            layer,
            hz,
            axis_source,
            axis_receiver,
            tiny,
            &mut || false
        )
        .is_err()
    );
    assert!(matches!(
        PreparedInteriorFullWaveForce::new(
            layer,
            hz,
            axis_source,
            axis_receiver,
            limits(),
            &mut || true,
        ),
        Err(AeroacError::Cancelled { .. })
    ));
    assert!(
        PreparedInteriorFullWaveForce::new(layer, 0.0, source, receiver, limits(), &mut || false)
            .is_err()
    );
    assert!(
        PreparedInteriorFullWaveForce::new(
            slab(media, 0.0),
            hz,
            source,
            receiver,
            limits(),
            &mut || false
        )
        .is_err()
    );
    for observer in [
        [0.0; 3],
        [0.0, 0.0, 0.018],
        [0.0, 0.0, -0.01],
        [0.0, 0.0, 0.03],
        [0.0, 0.0, f64::NAN],
    ] {
        assert!(
            PreparedInteriorFullWaveForce::new(layer, hz, source, observer, limits(), &mut || {
                false
            })
            .is_err()
        );
    }
    for source in [[0.0; 3], receiver, [0.0, 0.0, f64::NAN]] {
        assert!(
            PreparedInteriorFullWaveForce::new(layer, hz, source, receiver, limits(), &mut || {
                false
            })
            .is_err()
        );
    }
}
