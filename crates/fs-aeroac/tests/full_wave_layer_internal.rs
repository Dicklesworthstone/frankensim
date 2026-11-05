//! G1/G3: original sources and receivers inside a refracting finite layer.

use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::refraction::PlanarFluidInterface;
use fs_aeroac::curle3d::refraction::full_wave::FullWaveMonopoleLimits;
use fs_aeroac::curle3d::refraction::full_wave::layer::PlanarFluidLayer;
use fs_aeroac::curle3d::refraction::full_wave::layer::interior::internal::PreparedLayerInternalFullWaveMonopole;
use fs_aeroac::curle3d::refraction::full_wave::layer::interior::source::PreparedInteriorSourceFullWaveMonopole;
use fs_aeroac::curle3d::refraction::full_wave::reflection::PreparedReflectedFullWaveMonopole;
use fs_aeroac::curle3d::{CompactSourcePhasor, PreparedCompactSourceObserver};
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::c64::C64;
use std::f64::consts::{PI, TAU};

fn fluid(rho: f64, c: f64) -> FluidAcousticState {
    FluidAcousticState::try_new(rho, c).unwrap()
}
fn limits() -> FullWaveMonopoleLimits {
    FullWaveMonopoleLimits {
        absolute_pressure_tolerance_pa_per_m3_s: 1e-6,
        absolute_velocity_tolerance_per_m2: 1e-9,
        relative_tolerance: 1e-9,
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
) -> PreparedLayerInternalFullWaveMonopole {
    PreparedLayerInternalFullWaveMonopole::new(layer, hz, source, receiver, limits(), &mut || false)
        .unwrap()
}
fn close(actual: C64, expected: C64, relative: f64, absolute: f64) {
    assert!(
        (actual - expected).abs() <= absolute + relative * expected.abs(),
        "actual={actual:?}, expected={expected:?}, rel={relative:e}, abs={absolute:e}"
    );
}

#[test]
fn absent_faces_recover_original_interface_and_exact_homogeneous_sphere() {
    let source = [0.004, -0.003, 0.009];
    let receiver = [0.023, 0.014, 0.021];
    let q = C64::new(0.001, -0.0003);
    let inside = fluid(3.0, 340.0);
    for outside in [inside, fluid(1.2, 340.0), fluid(1.2, 470.0)] {
        let layer = slab([outside, inside, inside], 0.03);
        let prepared = prepare(layer, 800.0, source, receiver);
        let original = PreparedReflectedFullWaveMonopole::new(
            PlanarFluidInterface::new(inside, outside, [0.0; 3], [0.0, 0.0, -1.0]).unwrap(),
            800.0,
            source,
            receiver,
            limits(),
            &mut || false,
        )
        .unwrap();
        let actual = prepared.evaluate_total(q).unwrap();
        let expected = original.evaluate_total(q).unwrap();
        assert_eq!(actual.pressure_pa, expected.pressure_pa);
        assert_eq!(actual.particle_velocity_m_s, expected.particle_velocity_m_s);
        let reflected = prepared.evaluate_reflected(q).unwrap();
        let original_reflected = original.evaluate_reflected(q).unwrap();
        assert_eq!(reflected.pressure_pa, original_reflected.pressure_pa);
        assert_eq!(
            reflected.particle_velocity_m_s,
            original_reflected.particle_velocity_m_s
        );
        assert!(prepared.diagnostics().reduced_to_single_interface);
        assert_eq!(
            prepared.diagnostics().evaluations,
            original.diagnostics().evaluations
        );
        if outside == inside {
            let free = PreparedCompactSourceObserver::new(inside, 800.0, receiver, &[source], 1)
                .unwrap()
                .evaluate(&[CompactSourcePhasor::monopole(q)])
                .unwrap();
            close(actual.pressure_pa, free.pressure_pa, 2e-15, 0.0);
            for axis in 0..3 {
                close(
                    actual.particle_velocity_m_s[axis],
                    free.particle_velocity_m_s[axis],
                    2e-15,
                    0.0,
                );
            }
            assert_eq!(reflected.pressure_pa, C64::ZERO);
            assert_eq!(reflected.particle_velocity_m_s, [C64::ZERO; 3]);
            assert_eq!(prepared.diagnostics().evaluations, 0);
        }
    }
    // Only the upper face remains. Keep its original location; do not infer a
    // first-face coordinate or silently exchange the two source-fluid states.
    let outside = fluid(0.9, 390.0);
    let layer = slab([inside, inside, outside], 0.03);
    let prepared = prepare(layer, 800.0, source, receiver);
    let original = PreparedReflectedFullWaveMonopole::new(
        PlanarFluidInterface::new(inside, outside, [0.0, 0.0, 0.03], [0.0, 0.0, 1.0]).unwrap(),
        800.0,
        source,
        receiver,
        limits(),
        &mut || false,
    )
    .unwrap();
    let actual = prepared.evaluate_total(q).unwrap();
    let expected = original.evaluate_total(q).unwrap();
    close(actual.pressure_pa, expected.pressure_pa, 2e-7, 1e-9);
    for axis in 0..3 {
        close(
            actual.particle_velocity_m_s[axis],
            expected.particle_velocity_m_s[axis],
            2e-7,
            1e-11,
        );
    }
}

fn spherical(rho: f64, c: f64, hz: f64, source: [f64; 3], receiver: [f64; 3], q: C64) -> [C64; 4] {
    let displacement = core::array::from_fn::<_, 3, _>(|axis| receiver[axis] - source[axis]);
    let radius = displacement.iter().map(|v| v * v).sum::<f64>().sqrt();
    let omega = TAU * hz;
    let (sin, cos) = (omega / c * radius).sin_cos();
    let pressure =
        q * C64::new(0.0, -omega * rho) * C64::new(cos, sin).scale(1.0 / (4.0 * PI * radius));
    let velocity = pressure * C64::new(1.0 / (rho * c), 1.0 / (omega * rho * radius));
    [
        pressure,
        velocity.scale(displacement[0] / radius),
        velocity.scale(displacement[1] / radius),
        velocity.scale(displacement[2] / radius),
    ]
}

fn reflected_images(
    rho: [f64; 3],
    c: f64,
    h: f64,
    hz: f64,
    source: [f64; 3],
    receiver: [f64; 3],
    q: C64,
) -> [C64; 4] {
    let left = (rho[0] - rho[1]) / (rho[0] + rho[1]);
    let right = (rho[2] - rho[1]) / (rho[2] + rho[1]);
    let round_trip = left * right;
    let mut total = [C64::ZERO; 4];
    let mut weight = 1.0;
    // Four independent original spherical image families; no spectral kernel,
    // native boundary coefficient or numerical derivative is used here.
    for m in 0..160 {
        for (z, reflection) in [
            (-source[2] - 2.0 * m as f64 * h, left),
            (2.0 * h - source[2] + 2.0 * m as f64 * h, right),
            (source[2] - 2.0 * (m + 1) as f64 * h, round_trip),
            (source[2] + 2.0 * (m + 1) as f64 * h, round_trip),
        ] {
            let value = spherical(
                rho[1],
                c,
                hz,
                [source[0], source[1], z],
                receiver,
                q.scale(weight * reflection),
            );
            for row in 0..4 {
                total[row] = total[row] + value[row];
            }
        }
        weight *= round_trip;
    }
    total
}

#[test]
fn all_four_spherical_image_families_survive_equal_depth_and_deep_evanescence() {
    let rho = [1.2, 3.0, 0.8];
    let c = 340.0;
    let hz = 1700.0;
    let q = C64::new(0.0005, 0.0002);
    for (h, source, receiver) in [
        (0.03, [0.004, -0.002, 0.009], [0.027, 0.012, 0.018]),
        (0.03, [0.004, -0.002, 0.009], [0.027, 0.012, 0.009]),
        (0.6, [0.0, 0.0, 0.0002], [0.0007, 0.0003, 0.0003]),
        (0.6, [0.0, 0.0, 0.5998], [0.0007, 0.0003, 0.5997]),
    ] {
        let prepared = prepare(slab(rho.map(|r| fluid(r, c)), h), hz, source, receiver);
        let expected_reflected = reflected_images(rho, c, h, hz, source, receiver, q);
        let direct = spherical(rho[1], c, hz, source, receiver, q);
        let reflected = prepared.evaluate_reflected(q).unwrap();
        let total = prepared.evaluate_total(q).unwrap();
        close(reflected.pressure_pa, expected_reflected[0], 2e-7, 1e-9);
        close(
            total.pressure_pa,
            expected_reflected[0] + direct[0],
            2e-7,
            1e-9,
        );
        for axis in 0..3 {
            close(
                reflected.particle_velocity_m_s[axis],
                expected_reflected[axis + 1],
                2e-7,
                1e-11,
            );
            close(
                total.particle_velocity_m_s[axis],
                expected_reflected[axis + 1] + direct[axis + 1],
                2e-7,
                1e-11,
            );
        }
        if h > 0.1 {
            let cutoff = prepared.diagnostics().transverse_wavenumber_cutoff_rad_m;
            assert!((cutoff * cutoff - (TAU * hz / c).powi(2)).sqrt() * h > 700.0);
        }
    }
}

#[test]
fn reciprocal_pressure_and_all_receiving_momentum_components_use_actual_layer_fluid() {
    let layer = slab(
        [fluid(1.2, 300.0), fluid(2.8, 510.0), fluid(0.9, 380.0)],
        0.02,
    );
    let source = [0.003, -0.004, 0.007];
    let receiver = [0.031, 0.009, 0.013];
    let hz = 900.0;
    let q = C64::new(0.001, -0.0002);
    let prepared = prepare(layer, hz, source, receiver);
    let actual = prepared.evaluate_total(q).unwrap();
    let reciprocal = prepare(layer, hz, receiver, source);
    close(
        actual.pressure_pa,
        reciprocal.evaluate_total(q).unwrap().pressure_pa,
        2e-7,
        1e-9,
    );
    close(
        prepared.evaluate_reflected(q).unwrap().pressure_pa,
        reciprocal.evaluate_reflected(q).unwrap().pressure_pa,
        2e-7,
        1e-9,
    );
    let step = 1e-6;
    let mut divergence = C64::ZERO;
    for axis in 0..3 {
        let mut lower = receiver;
        let mut upper = receiver;
        lower[axis] -= step;
        upper[axis] += step;
        let lower = prepare(layer, hz, source, lower).evaluate_total(q).unwrap();
        let upper = prepare(layer, hz, source, upper).evaluate_total(q).unwrap();
        let gradient = (upper.pressure_pa - lower.pressure_pa).scale(0.5 / step);
        let momentum = C64::new(0.0, TAU * hz * 2.8) * actual.particle_velocity_m_s[axis];
        close(momentum, gradient, 2e-5, 0.002);
        divergence = divergence
            + (upper.particle_velocity_m_s[axis] - lower.particle_velocity_m_s[axis])
                .scale(0.5 / step);
    }
    close(
        divergence,
        C64::new(0.0, TAU * hz / (2.8 * 510.0 * 510.0)) * actual.pressure_pa,
        2e-4,
        2e-6,
    );
    let line = prepared.evaluate_total_harmonic_line(q).unwrap();
    assert_eq!(line.frequency_hz, hz);
    assert_eq!(line.pressure_pa, actual.pressure_pa);
    assert_eq!(line.particle_velocity_m_s, actual.particle_velocity_m_s);
    assert!(
        actual
            .reactive_intensity_w_m2
            .iter()
            .any(|v| v.abs() > 1e-10)
    );

    let rotate = |v: [f64; 3]| [v[2], v[0], v[1]];
    let shift = [0.5, -0.25, 0.75];
    let point = |v| {
        let v = rotate(v);
        core::array::from_fn(|i| v[i] + shift[i])
    };
    let moved_layer = PlanarFluidLayer::new(
        layer.incident_state(),
        layer.layer_state(),
        layer.receiving_state(),
        shift,
        rotate([0.0, 0.0, 1.0]),
        layer.thickness_m(),
    )
    .unwrap();
    let moved = prepare(moved_layer, hz, point(source), point(receiver))
        .evaluate_total(q)
        .unwrap();
    close(moved.pressure_pa, actual.pressure_pa, 2e-7, 1e-9);
    for (destination, original) in [(0, 2), (1, 0), (2, 1)] {
        close(
            moved.particle_velocity_m_s[destination],
            actual.particle_velocity_m_s[original],
            2e-7,
            1e-11,
        );
    }
}

#[test]
fn total_field_meets_both_original_fluid_faces() {
    let layer = slab(
        [fluid(1.1, 320.0), fluid(2.5, 480.0), fluid(0.8, 365.0)],
        0.026,
    );
    let reversed = PlanarFluidLayer::new(
        layer.receiving_state(),
        layer.layer_state(),
        layer.incident_state(),
        [0.0, 0.0, layer.thickness_m()],
        [0.0, 0.0, -1.0],
        layer.thickness_m(),
    )
    .unwrap();
    let source = [0.001, -0.002, 0.009];
    let q = C64::new(0.001, 0.0002);
    let epsilon = 1e-8;
    for (external_layer, z, sign, rho) in [
        (layer, 0.0, 1.0, layer.incident_state().density_kg_m3()),
        (
            reversed,
            layer.thickness_m(),
            -1.0,
            layer.receiving_state().density_kg_m3(),
        ),
    ] {
        let inside = prepare(layer, 800.0, source, [0.014, 0.006, z + sign * epsilon])
            .evaluate_total(q)
            .unwrap();
        let outside = PreparedInteriorSourceFullWaveMonopole::new(
            external_layer,
            800.0,
            source,
            [0.014, 0.006, z - sign * epsilon],
            limits(),
            &mut || false,
        )
        .unwrap()
        .evaluate(q)
        .unwrap();
        close(inside.pressure_pa, outside.pressure_pa, 5e-6, 1e-9);
        close(
            inside.particle_velocity_m_s[2],
            outside.particle_velocity_m_s[2],
            5e-6,
            1e-11,
        );
        for axis in 0..2 {
            close(
                inside.particle_velocity_m_s[axis].scale(2.5),
                outside.particle_velocity_m_s[axis].scale(rho),
                5e-6,
                1e-11,
            );
        }
    }
}

#[test]
fn symmetric_returns_cancel_normally_but_the_reflected_self_field_is_finite() {
    let exterior = fluid(1.2, 340.0);
    let layer = slab([exterior, fluid(2.5, 470.0), exterior], 0.02);
    let source = [0.0, 0.0, 0.01];
    let mut controls = limits();
    controls.relative_tolerance = 0.5;
    let prepared = PreparedLayerInternalFullWaveMonopole::new(
        layer,
        750.0,
        source,
        source,
        controls,
        &mut || false,
    )
    .unwrap();
    assert!(!prepared.total_field_available());
    assert!(prepared.evaluate_total(C64::ONE).is_err());
    assert!(prepared.evaluate_total(C64::ZERO).is_err());
    assert!(prepared.evaluate_total_harmonic_line(C64::ONE).is_err());
    let reflected = prepared.evaluate_reflected(C64::ONE).unwrap();
    assert!(reflected.pressure_pa.abs().is_finite() && reflected.pressure_pa.abs() > 1.0);
    assert_eq!(reflected.particle_velocity_m_s, [C64::ZERO; 3]);
    let d = prepared.diagnostics();
    assert!(
        d.pressure_quadrature_error_pa_per_m3_s + d.pressure_tail_bound_pa_per_m3_s
            <= controls.absolute_pressure_tolerance_pa_per_m3_s
    );
    let velocity_error = core::array::from_fn::<_, 2, _>(|i| {
        d.velocity_quadrature_error_per_m2[i] + d.velocity_tail_bound_per_m2[i]
    });
    assert!(
        velocity_error[0].hypot(velocity_error[1]) <= controls.absolute_velocity_tolerance_per_m2
    );
    let line = prepared.evaluate_reflected_harmonic_line(C64::ONE).unwrap();
    assert_eq!(line.pressure_pa, reflected.pressure_pa);
    assert_eq!(line.particle_velocity_m_s, reflected.particle_velocity_m_s);
    let lateral = prepare(layer, 750.0, source, [0.025, 0.0, 0.01])
        .evaluate_total(C64::ONE)
        .unwrap();
    assert_eq!(lateral.particle_velocity_m_s[2], C64::ZERO);
    assert_eq!(lateral.particle_velocity_m_s[1], C64::ZERO);
    assert!(lateral.particle_velocity_m_s[0].abs() > 0.0);
}

#[test]
fn geometry_cancellation_shared_budget_and_phasor_refusals_are_atomic() {
    let layer = slab(
        [fluid(1.2, 300.0), fluid(2.8, 510.0), fluid(0.9, 380.0)],
        0.02,
    );
    let source = [0.003, -0.004, 0.007];
    let receiver = [0.031, 0.009, 0.013];
    let mut polls = 0;
    let prepared = PreparedLayerInternalFullWaveMonopole::new(
        layer,
        900.0,
        source,
        receiver,
        limits(),
        &mut || {
            polls += 1;
            false
        },
    )
    .unwrap();
    assert_eq!(prepared.layer().layer_state(), layer.layer_state());
    assert_eq!(prepared.layer().thickness_m(), layer.thickness_m());
    assert_eq!(prepared.frequency_hz(), 900.0);
    assert_eq!(prepared.source_position_m(), source);
    assert_eq!(prepared.observer_position_m(), receiver);
    let mut capped = limits();
    capped.maximum_evaluations = prepared.diagnostics().evaluations - 1;
    assert!(
        PreparedLayerInternalFullWaveMonopole::new(
            layer,
            900.0,
            source,
            receiver,
            capped,
            &mut || false
        )
        .is_err()
    );
    assert!(matches!(
        PreparedLayerInternalFullWaveMonopole::new(
            layer,
            900.0,
            source,
            receiver,
            limits(),
            &mut || true
        ),
        Err(AeroacError::Cancelled { .. })
    ));
    let mut late = 0;
    assert!(matches!(
        PreparedLayerInternalFullWaveMonopole::new(
            layer,
            900.0,
            source,
            receiver,
            limits(),
            &mut || {
                late += 1;
                late >= polls
            }
        ),
        Err(AeroacError::Cancelled { .. })
    ));
    for (bad_source, bad_receiver) in [
        ([0.0; 3], receiver),
        (source, [0.0, 0.0, 0.02]),
        ([f64::NAN, 0.0, 0.01], receiver),
        (source, [0.0, 0.0, -0.01]),
    ] {
        assert!(
            PreparedLayerInternalFullWaveMonopole::new(
                layer,
                900.0,
                bad_source,
                bad_receiver,
                limits(),
                &mut || false
            )
            .is_err()
        );
    }
    assert!(
        PreparedLayerInternalFullWaveMonopole::new(
            layer,
            0.0,
            source,
            receiver,
            limits(),
            &mut || false
        )
        .is_err()
    );
    let zero = slab(
        [
            layer.incident_state(),
            layer.layer_state(),
            layer.receiving_state(),
        ],
        0.0,
    );
    assert!(
        PreparedLayerInternalFullWaveMonopole::new(
            zero,
            900.0,
            source,
            receiver,
            limits(),
            &mut || false
        )
        .is_err()
    );
    let q = C64::new(0.001, -0.0002);
    let before = prepared.evaluate_total(q).unwrap();
    assert!(
        prepared
            .evaluate_reflected(C64::new(f64::NAN, 0.0))
            .is_err()
    );
    assert!(
        prepared
            .evaluate_total_harmonic_line(C64::new(0.0, f64::INFINITY))
            .is_err()
    );
    assert_eq!(
        prepared.evaluate_total(q).unwrap().pressure_pa,
        before.pressure_pa
    );
    let silence = prepared.evaluate_total(C64::ZERO).unwrap();
    assert_eq!(silence.pressure_pa, C64::ZERO);
    assert_eq!(silence.particle_velocity_m_s, [C64::ZERO; 3]);
    assert_eq!(silence.active_intensity_w_m2, [0.0; 3]);
    assert_eq!(silence.reactive_intensity_w_m2, [0.0; 3]);
}
