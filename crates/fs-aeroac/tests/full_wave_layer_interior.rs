//! G1/G3: physical standing fields, original interface continuity and momentum.

use fs_aeroac::curle3d::refraction::PlanarFluidInterface;
use fs_aeroac::curle3d::refraction::full_wave::layer::interior::PreparedLayerInteriorFullWaveMonopole;
use fs_aeroac::curle3d::refraction::full_wave::layer::reflection::PreparedLayeredReflectedFullWaveMonopole;
use fs_aeroac::curle3d::refraction::full_wave::layer::{
    PlanarFluidLayer, PreparedLayeredFullWaveMonopole,
};
use fs_aeroac::curle3d::refraction::full_wave::{FullWaveMonopoleLimits, PreparedFullWaveMonopole};
use fs_aeroac::curle3d::{
    CoherentCompactField, CompactSourcePhasor, PreparedCompactSourceObserver,
};
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
    frequency: f64,
    source: [f64; 3],
    receiver: [f64; 3],
) -> PreparedLayerInteriorFullWaveMonopole {
    PreparedLayerInteriorFullWaveMonopole::new(
        layer,
        frequency,
        source,
        receiver,
        limits(),
        &mut || false,
    )
    .unwrap()
}
fn close(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-11),
        "actual={actual:?}, expected={expected:?}, tolerance={tolerance:e}"
    );
}
fn compare(actual: &CoherentCompactField, expected: &CoherentCompactField, tolerance: f64) {
    close(actual.pressure_pa, expected.pressure_pa, tolerance);
    for axis in 0..3 {
        close(
            actual.particle_velocity_m_s[axis],
            expected.particle_velocity_m_s[axis],
            tolerance,
        );
    }
}

#[test]
fn matched_exit_keeps_exact_original_spherical_and_single_interface_fields() {
    let receiving = fluid(2.4, 450.0);
    let source = [-0.02, 0.01, -0.03];
    let receiver = [0.035, -0.025, 0.012];
    let q = C64::new(0.002, -0.001);
    for incident in [receiving, fluid(1.2, 300.0)] {
        let layer = slab([incident, receiving, receiving], 0.04);
        let interior = prepare(layer, 900.0, source, receiver);
        let interface =
            PlanarFluidInterface::new(incident, receiving, [0.0; 3], [0.0, 0.0, 1.0]).unwrap();
        let original = PreparedFullWaveMonopole::new(
            interface,
            900.0,
            source,
            receiver,
            limits(),
            &mut || false,
        )
        .unwrap();
        let expected = original.evaluate(q).unwrap();
        let actual = interior.evaluate(q).unwrap();
        assert_eq!(actual.pressure_pa, expected.pressure_pa);
        assert_eq!(actual.particle_velocity_m_s, expected.particle_velocity_m_s);
        assert!(interior.diagnostics().reduced_to_single_interface);
        if incident == receiving {
            let spherical =
                PreparedCompactSourceObserver::new(receiving, 900.0, receiver, &[source], 1)
                    .unwrap()
                    .evaluate(&[CompactSourcePhasor::monopole(q)])
                    .unwrap();
            compare(&actual, &spherical, 8e-14);
            assert_eq!(interior.diagnostics().evaluations, 0);
            assert!(actual.reactive_intensity_w_m2.iter().any(|v| v.abs() > 0.0));
        }
        let line = interior.evaluate_harmonic_line(q).unwrap();
        assert_eq!(line.frequency_hz, 900.0);
        assert_eq!(line.pressure_pa, actual.pressure_pa);
        assert_eq!(line.particle_velocity_m_s, actual.particle_velocity_m_s);
        assert_eq!(interior.source_position_m(), source);
        assert_eq!(interior.observer_position_m(), receiver);
        assert_eq!(interior.layer().thickness_m(), 0.04);
        let zero = interior.evaluate(C64::ZERO).unwrap();
        assert_eq!(zero.pressure_pa, C64::ZERO);
        assert_eq!(zero.particle_velocity_m_s, [C64::ZERO; 3]);
        assert!(interior.evaluate(C64::new(f64::NAN, 0.0)).is_err());
    }
}

// Independent closed spherical-image series when the three speeds agree.
// Each pass contributes a source behind the entrance and its far-face image.
// Pressure uses the original rho1*Q; velocity uses the receiving rho2 and the
// actual signed image direction, retaining both near-field induction terms.
fn image_series(rho: [f64; 3], c: f64, hz: f64, a: f64, h: f64, receiver: [f64; 3]) -> [C64; 4] {
    let omega = TAU * hz;
    let k = omega / c;
    let t12 = 2.0 * rho[1] / (rho[0] + rho[1]);
    let r21 = (rho[0] - rho[1]) / (rho[0] + rho[1]);
    let r23 = (rho[2] - rho[1]) / (rho[2] + rho[1]);
    let mut total = [C64::ZERO; 4];
    let mut return_weight = 1.0;
    for bounce in 0..160 {
        for (source_z, weight) in [
            (-a - 2.0 * bounce as f64 * h, return_weight),
            (a + 2.0 * (bounce + 1) as f64 * h, return_weight * r23),
        ] {
            let delta = [receiver[0], receiver[1], receiver[2] - source_z];
            let distance = delta.iter().map(|v| v * v).sum::<f64>().sqrt();
            let (sin, cos) = (k * distance).sin_cos();
            let pressure =
                C64::new(sin, -cos).scale(omega * rho[0] * t12 * weight / (4.0 * PI * distance));
            total[0] = total[0] + pressure;
            let radial = pressure * C64::new(1.0, 1.0 / (k * distance)).scale(1.0 / (rho[1] * c));
            for axis in 0..3 {
                total[axis + 1] = total[axis + 1] + radial.scale(delta[axis] / distance);
            }
        }
        return_weight *= r21 * r23;
    }
    total
}

#[test]
fn independent_spherical_images_recover_both_internal_directions_and_near_fields() {
    for (rho, h, z) in [
        ([1.2, 3.4, 1.7], 0.04, 0.019),
        ([1.2, 15.0, 2.4], 0.06, 0.052),
        // At this near-entrance receiver the evanescent tail has alpha*h>700.
        // An exterior transfer times a growing interior propagator underflows.
        ([1.2, 3.4, 1.7], 0.6, 0.001),
    ] {
        let media = rho.map(|v| fluid(v, 300.0));
        let receiver = [0.024, 0.013, z];
        let prepared = prepare(slab(media, h), 700.0, [0.0, 0.0, -0.015], receiver);
        let actual = prepared.evaluate(C64::ONE).unwrap();
        let expected = image_series(rho, 300.0, 700.0, 0.015, h, receiver);
        close(actual.pressure_pa, expected[0], 5e-8);
        for axis in 0..3 {
            close(actual.particle_velocity_m_s[axis], expected[axis + 1], 5e-8);
        }
        assert!(prepared.diagnostics().evaluations > 0);
        assert!(
            prepared.diagnostics().pressure_tail_bound_pa_per_m3_s
                <= limits().absolute_pressure_tolerance_pa_per_m3_s / 4.0
        );
        assert!(
            prepared
                .diagnostics()
                .velocity_tail_bound_per_m2
                .iter()
                .all(|v| *v <= limits().absolute_velocity_tolerance_per_m2 / 4.0)
        );
    }
}

#[test]
fn thick_layer_near_field_resolves_without_reintegrating_the_canceled_evanescent_exponent() {
    let rho = [1.2, 3.4, 1.7];
    let h = 100_000.0;
    let hz = 0.00001;
    let receiver = [0.024, 0.013, 0.01];
    let mut accuracy = limits();
    accuracy.absolute_pressure_tolerance_pa_per_m3_s = 1e-14;
    let prepared = PreparedLayerInteriorFullWaveMonopole::new(
        slab(rho.map(|v| fluid(v, 300.0)), h),
        hz,
        [0.0, 0.0, -0.015],
        receiver,
        accuracy,
        &mut || false,
    )
    .unwrap();
    let actual = prepared.evaluate(C64::ONE).unwrap();
    let expected = image_series(rho, 300.0, hz, 0.015, h, receiver);
    close(actual.pressure_pa, expected[0], 5e-8);
    for axis in 0..3 {
        close(actual.particle_velocity_m_s[axis], expected[axis + 1], 5e-8);
    }
    assert!(prepared.diagnostics().evaluations < accuracy.maximum_evaluations);
}

#[test]
fn both_faces_preserve_pressure_normal_velocity_and_tangential_momentum() {
    let media = [fluid(1.2, 300.0), fluid(2.4, 600.0), fluid(1.8, 400.0)];
    let h = 0.018;
    let layer = slab(media, h);
    let source = [0.0, 0.0, -0.025];
    let eps = 1e-8;
    let first = prepare(layer, 700.0, source, [0.04, 0.01, eps])
        .evaluate(C64::ONE)
        .unwrap();
    let before = PreparedLayeredReflectedFullWaveMonopole::new(
        layer,
        700.0,
        source,
        [0.04, 0.01, -eps],
        limits(),
        &mut || false,
    )
    .unwrap()
    .evaluate_total(C64::ONE)
    .unwrap();
    let last = prepare(layer, 700.0, source, [0.04, 0.01, h - eps])
        .evaluate(C64::ONE)
        .unwrap();
    let after = PreparedLayeredFullWaveMonopole::new(
        layer,
        700.0,
        source,
        [0.04, 0.01, h + eps],
        limits(),
        &mut || false,
    )
    .unwrap()
    .evaluate(C64::ONE)
    .unwrap();
    for (inside, outside, exterior) in [(first, before, media[0]), (last, after, media[2])] {
        close(inside.pressure_pa, outside.pressure_pa, 4e-6);
        close(
            inside.particle_velocity_m_s[2],
            outside.particle_velocity_m_s[2],
            4e-6,
        );
        for axis in 0..2 {
            close(
                inside.particle_velocity_m_s[axis].scale(media[1].density_kg_m3()),
                outside.particle_velocity_m_s[axis].scale(exterior.density_kg_m3()),
                4e-6,
            );
        }
    }
}

#[test]
fn interior_velocity_obeys_receiving_momentum_and_rotates_with_original_geometry() {
    let media = [fluid(1.2, 300.0), fluid(2.4, 600.0), fluid(1.8, 400.0)];
    let layer = slab(media, 0.018);
    let source = [0.0, 0.0, -0.025];
    let receiver = [0.04, 0.013, 0.007];
    let hz = 700.0;
    let original = prepare(layer, hz, source, receiver)
        .evaluate(C64::ONE)
        .unwrap();
    let step = 1e-6;
    let mut divergence = C64::ZERO;
    for axis in 0..3 {
        let mut before = receiver;
        let mut after = receiver;
        before[axis] -= step;
        after[axis] += step;
        let after = prepare(layer, hz, source, after)
            .evaluate(C64::ONE)
            .unwrap();
        let before = prepare(layer, hz, source, before)
            .evaluate(C64::ONE)
            .unwrap();
        let derivative = (after.pressure_pa - before.pressure_pa).scale(0.5 / step);
        divergence = divergence
            + (after.particle_velocity_m_s[axis] - before.particle_velocity_m_s[axis])
                .scale(0.5 / step);
        close(
            derivative,
            original.particle_velocity_m_s[axis]
                * C64::new(0.0, TAU * hz * media[1].density_kg_m3()),
            3e-7,
        );
    }
    // Continuity tests the local compressibility as well as the momentum
    // density: a receiving-speed substitution would violate this identity.
    let bulk = media[1].density_kg_m3() * media[1].sound_speed_m_s().powi(2);
    close(
        divergence,
        original.pressure_pa * C64::new(0.0, TAU * hz / bulk),
        2e-6,
    );
    let rotate = |v: [f64; 3]| [v[2], v[0], v[1]];
    let offset = [2.0, -3.0, 5.0];
    let transform = |v: [f64; 3]| {
        let v = rotate(v);
        core::array::from_fn(|i| v[i] + offset[i])
    };
    let rotated = PlanarFluidLayer::new(
        media[0],
        media[1],
        media[2],
        offset,
        rotate([0.0, 0.0, 1.0]),
        0.018,
    )
    .unwrap();
    let moved = prepare(rotated, hz, transform(source), transform(receiver))
        .evaluate(C64::ONE)
        .unwrap();
    close(moved.pressure_pa, original.pressure_pa, 5e-9);
    for (actual, index) in moved.particle_velocity_m_s.iter().zip([2, 0, 1]) {
        close(*actual, original.particle_velocity_m_s[index], 5e-9);
    }
    let axial = prepare(layer, hz, source, [0.0, 0.0, 0.007])
        .evaluate(C64::ONE)
        .unwrap();
    assert_eq!(axial.particle_velocity_m_s[0], C64::ZERO);
    assert_eq!(axial.particle_velocity_m_s[1], C64::ZERO);
}

#[test]
fn original_work_limits_cancellation_and_interior_geometry_refuse_without_an_owner() {
    let layer = slab(
        [fluid(1.2, 300.0), fluid(2.4, 600.0), fluid(1.8, 400.0)],
        0.018,
    );
    let source = [0.0, 0.0, -0.025];
    let receiver = [0.04, 0.01, 0.007];
    let mut tiny = limits();
    tiny.maximum_evaluations = 1;
    assert!(
        PreparedLayerInteriorFullWaveMonopole::new(
            layer,
            700.0,
            source,
            receiver,
            tiny,
            &mut || false
        )
        .is_err()
    );
    assert!(
        PreparedLayerInteriorFullWaveMonopole::new(
            layer,
            700.0,
            source,
            receiver,
            limits(),
            &mut || true
        )
        .is_err()
    );
    let mut calls = 0;
    assert!(
        PreparedLayerInteriorFullWaveMonopole::new(
            layer,
            700.0,
            source,
            receiver,
            limits(),
            &mut || {
                calls += 1;
                calls >= 80
            }
        )
        .is_err()
    );
    assert_eq!(calls, 80);
    for position in [
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.018],
        [0.0, 0.0, -0.01],
        [0.0, 0.0, 0.02],
        [f64::INFINITY, 0.0, 0.007],
    ] {
        assert!(
            PreparedLayerInteriorFullWaveMonopole::new(
                layer,
                700.0,
                source,
                position,
                limits(),
                &mut || false
            )
            .is_err()
        );
    }
    assert!(
        PreparedLayerInteriorFullWaveMonopole::new(
            layer,
            700.0,
            [0.0; 3],
            receiver,
            limits(),
            &mut || false
        )
        .is_err()
    );
}
