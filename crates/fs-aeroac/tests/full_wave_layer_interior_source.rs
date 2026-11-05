//! G1/G3: physical layer sources, reciprocal pressure and exterior momentum.

use fs_aeroac::curle3d::refraction::PlanarFluidInterface;
use fs_aeroac::curle3d::refraction::full_wave::layer::PlanarFluidLayer;
use fs_aeroac::curle3d::refraction::full_wave::layer::interior::PreparedLayerInteriorFullWaveMonopole;
use fs_aeroac::curle3d::refraction::full_wave::layer::interior::source::PreparedInteriorSourceFullWaveMonopole;
use fs_aeroac::curle3d::refraction::full_wave::{FullWaveMonopoleLimits, PreparedFullWaveMonopole};
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::c64::C64;
use std::f64::consts::{PI, TAU};

fn fluid(rho: f64, c: f64) -> FluidAcousticState {
    FluidAcousticState::try_new(rho, c).unwrap()
}
fn limits() -> FullWaveMonopoleLimits {
    FullWaveMonopoleLimits {
        absolute_pressure_tolerance_pa_per_m3_s: 1e-9,
        absolute_velocity_tolerance_per_m2: 1e-12,
        relative_tolerance: 1e-10,
        maximum_evaluations: 262_144,
        maximum_depth: 40,
    }
}
fn slab(media: [FluidAcousticState; 3], h: f64) -> PlanarFluidLayer {
    PlanarFluidLayer::new(media[0], media[1], media[2], [0.0; 3], [0.0, 0.0, 1.0], h).unwrap()
}
fn prepare(
    layer: PlanarFluidLayer,
    hz: f64,
    source: [f64; 3],
    receiver: [f64; 3],
) -> PreparedInteriorSourceFullWaveMonopole {
    PreparedInteriorSourceFullWaveMonopole::new(layer, hz, source, receiver, limits(), &mut || {
        false
    })
    .unwrap()
}
fn close(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-12),
        "actual={actual:?} expected={expected:?} relative tolerance={tolerance:e}"
    );
}

#[test]
fn interior_source_recovers_original_single_interface_when_second_face_is_absent() {
    let source = [0.006, -0.003, 0.011];
    let receiver = [0.029, 0.014, -0.017];
    let q = C64::new(0.002, -0.001);
    for (outside, inside) in [
        (fluid(1.2, 340.0), fluid(1.2, 340.0)),
        (fluid(1.2, 340.0), fluid(4.0, 490.0)),
        (fluid(1.2, 340.0), fluid(4.0, 340.0)),
    ] {
        let layer = slab([outside, inside, inside], 0.03);
        let original = PreparedFullWaveMonopole::new(
            PlanarFluidInterface::new(inside, outside, [0.0; 3], [0.0, 0.0, -1.0]).unwrap(),
            800.0,
            source,
            receiver,
            limits(),
            &mut || false,
        )
        .unwrap();
        let prepared = prepare(layer, 800.0, source, receiver);
        let actual = prepared.evaluate(q).unwrap();
        let expected = original.evaluate(q).unwrap();
        assert_eq!(actual.pressure_pa, expected.pressure_pa);
        assert_eq!(actual.particle_velocity_m_s, expected.particle_velocity_m_s);
        assert!(prepared.diagnostics().reduced_to_single_interface);
        assert_eq!(
            prepared.diagnostics().evaluations,
            original.diagnostics().evaluations
        );
        let line = prepared.evaluate_harmonic_line(q).unwrap();
        assert_eq!(line.pressure_pa, actual.pressure_pa);
        assert_eq!(line.particle_velocity_m_s, actual.particle_velocity_m_s);
    }
}

#[test]
fn source_inside_retains_both_coherent_spherical_image_families() {
    let [rho1, rho2, rho3] = [1.2, 3.0, 0.8];
    let c = 340.0;
    let h = 0.024;
    let hz = 1700.0;
    let source = [0.004, -0.002, 0.009];
    let receiver = [0.027, 0.012, -0.014];
    let layer = slab([fluid(rho1, c), fluid(rho2, c), fluid(rho3, c)], h);
    let q = C64::new(0.001, 0.0003);
    let actual = prepare(layer, hz, source, receiver).evaluate(q).unwrap();
    let omega = TAU * hz;
    let k = omega / c;
    let r21 = (rho1 - rho2) / (rho1 + rho2);
    let r23 = (rho3 - rho2) / (rho3 + rho2);
    let transmission = 2.0 * rho1 / (rho1 + rho2);
    let mut expected = [C64::ZERO; 4];
    let mut return_weight = transmission;
    for m in 0..100 {
        for (image_z, reflection) in [
            (source[2] + 2.0 * m as f64 * h, 1.0),
            (2.0 * h - source[2] + 2.0 * m as f64 * h, r23),
        ] {
            let displacement = [
                receiver[0] - source[0],
                receiver[1] - source[1],
                receiver[2] - image_z,
            ];
            let radius = displacement.iter().map(|x| x * x).sum::<f64>().sqrt();
            let (sin, cos) = (k * radius).sin_cos();
            let p = q
                * C64::new(0.0, -omega * rho2)
                * C64::new(cos, sin).scale(return_weight * reflection / (4.0 * PI * radius));
            expected[0] = expected[0] + p;
            let radial_velocity = p * C64::new(1.0 / (rho1 * c), 1.0 / (omega * rho1 * radius));
            for axis in 0..3 {
                expected[axis + 1] =
                    expected[axis + 1] + radial_velocity.scale(displacement[axis] / radius);
            }
        }
        return_weight *= r21 * r23;
    }
    close(actual.pressure_pa, expected[0], 2e-8);
    for axis in 0..3 {
        close(actual.particle_velocity_m_s[axis], expected[axis + 1], 2e-8);
    }
    assert!(
        actual
            .reactive_intensity_w_m2
            .iter()
            .any(|v| v.abs() > 1e-9)
    );
}

#[test]
fn pressure_is_reciprocal_but_velocity_uses_the_actual_receiving_fluid() {
    let layer = slab(
        [fluid(1.2, 300.0), fluid(2.8, 510.0), fluid(0.9, 380.0)],
        0.02,
    );
    let source = [0.003, -0.004, 0.007];
    let receiver = [0.031, 0.009, -0.018];
    let hz = 900.0;
    let q = C64::new(0.001, -0.0002);
    let prepared = prepare(layer, hz, source, receiver);
    let actual = prepared.evaluate(q).unwrap();
    let reciprocal = PreparedLayerInteriorFullWaveMonopole::new(
        layer,
        hz,
        receiver,
        source,
        limits(),
        &mut || false,
    )
    .unwrap()
    .evaluate(q)
    .unwrap();
    close(actual.pressure_pa, reciprocal.pressure_pa, 2e-9);
    let step = 2e-7;
    for axis in 0..3 {
        let mut lower = receiver;
        let mut upper = receiver;
        lower[axis] -= step;
        upper[axis] += step;
        let derivative = (prepare(layer, hz, source, upper)
            .evaluate(q)
            .unwrap()
            .pressure_pa
            - prepare(layer, hz, source, lower)
                .evaluate(q)
                .unwrap()
                .pressure_pa)
            .scale(0.5 / step);
        let momentum = C64::new(0.0, TAU * hz * layer.incident_state().density_kg_m3())
            * actual.particle_velocity_m_s[axis];
        close(momentum, derivative, 3e-6);
    }
    assert_eq!(prepared.source_position_m(), source);
    assert_eq!(prepared.observer_position_m(), receiver);
}

#[test]
fn rotated_source_radiation_and_native_observation_clock_keep_original_field() {
    let layer = slab(
        [fluid(1.1, 320.0), fluid(2.5, 460.0), fluid(0.8, 370.0)],
        0.03,
    );
    let source = [0.004, 0.002, 0.011];
    let receiver = [0.032, -0.014, -0.02];
    // Proper cyclic rotation, followed by a finite translation.
    let rotate = |v: [f64; 3]| [v[2], v[0], v[1]];
    let shift = [0.5, -0.25, 0.75];
    let point = |v| {
        let v = rotate(v);
        core::array::from_fn(|i| v[i] + shift[i])
    };
    let transformed = PlanarFluidLayer::new(
        layer.incident_state(),
        layer.layer_state(),
        layer.receiving_state(),
        shift,
        rotate([0.0, 0.0, 1.0]),
        layer.thickness_m(),
    )
    .unwrap();
    let q = C64::new(0.0004, -0.0003);
    let base = prepare(layer, 600.0, source, receiver).evaluate(q).unwrap();
    let prepared = prepare(transformed, 600.0, point(source), point(receiver));
    let moved = prepared.evaluate(q).unwrap();
    close(moved.pressure_pa, base.pressure_pa, 2e-8);
    for (destination, original) in [(0, 2), (1, 0), (2, 1)] {
        close(
            moved.particle_velocity_m_s[destination],
            base.particle_velocity_m_s[original],
            2e-8,
        );
    }
    let line = prepared.evaluate_harmonic_line(q).unwrap();
    assert_eq!(line.pressure_pa, moved.pressure_pa);
    assert_eq!(line.particle_velocity_m_s, moved.particle_velocity_m_s);
    assert_eq!(line.frequency_hz, 600.0);
}

#[test]
fn interior_source_refuses_wrong_sides_missing_layer_exhaustion_and_cancellation() {
    let layer = slab(
        [fluid(1.2, 300.0), fluid(3.0, 500.0), fluid(0.8, 360.0)],
        0.02,
    );
    let source = [0.0, 0.0, 0.007];
    let receiver = [0.02, 0.0, -0.01];
    for z in [-0.01, 0.0, 0.02, 0.03, f64::NAN] {
        assert!(
            PreparedInteriorSourceFullWaveMonopole::new(
                layer,
                800.0,
                [0.0, 0.0, z],
                receiver,
                limits(),
                &mut || false
            )
            .is_err()
        );
    }
    for z in [0.0, 0.01, 0.03] {
        assert!(
            PreparedInteriorSourceFullWaveMonopole::new(
                layer,
                800.0,
                source,
                [0.0, 0.0, z],
                limits(),
                &mut || false
            )
            .is_err()
        );
    }
    assert!(
        PreparedInteriorSourceFullWaveMonopole::new(
            layer,
            800.0,
            source,
            receiver,
            limits(),
            &mut || true
        )
        .is_err()
    );
    let mut tiny = limits();
    tiny.maximum_evaluations = 1;
    assert!(
        PreparedInteriorSourceFullWaveMonopole::new(
            layer,
            800.0,
            source,
            receiver,
            tiny,
            &mut || false
        )
        .is_err()
    );
    let absent = slab(
        [
            layer.incident_state(),
            layer.layer_state(),
            layer.receiving_state(),
        ],
        0.0,
    );
    assert!(
        PreparedInteriorSourceFullWaveMonopole::new(
            absent,
            800.0,
            source,
            receiver,
            limits(),
            &mut || false
        )
        .is_err()
    );
    assert!(
        prepare(layer, 800.0, source, receiver)
            .evaluate(C64::new(f64::NAN, 0.0))
            .is_err()
    );
}
