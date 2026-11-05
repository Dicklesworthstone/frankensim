//! G1/G3: coherent finite-layer returns, source clocks and full particle motion.

use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::refraction::PlanarFluidInterface;
use fs_aeroac::curle3d::refraction::full_wave::FullWaveMonopoleLimits;
use fs_aeroac::curle3d::refraction::full_wave::layer::PlanarFluidLayer;
use fs_aeroac::curle3d::refraction::full_wave::layer::reflection::PreparedLayeredReflectedFullWaveMonopole;
use fs_aeroac::curle3d::refraction::full_wave::reflection::PreparedReflectedFullWaveMonopole;
use fs_aeroac::curle3d::{CompactSourcePhasor, PreparedCompactSourceObserver};
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::c64::C64;
use std::f64::consts::TAU;

fn fluid(rho: f64, c: f64) -> FluidAcousticState {
    FluidAcousticState::try_new(rho, c).unwrap()
}
fn limits() -> FullWaveMonopoleLimits {
    FullWaveMonopoleLimits {
        absolute_pressure_tolerance_pa_per_m3_s: 1e-10,
        absolute_velocity_tolerance_per_m2: 1e-13,
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
    receiver: [f64; 3],
) -> PreparedLayeredReflectedFullWaveMonopole {
    PreparedLayeredReflectedFullWaveMonopole::new(
        layer,
        hz,
        source,
        receiver,
        limits(),
        &mut || false,
    )
    .unwrap()
}
fn close(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-10),
        "actual={actual:?}, expected={expected:?}, tolerance={tolerance:e}"
    );
}

#[test]
fn independent_direct_spectra_include_all_three_critical_circles() {
    // Independent mpmath 1.3.0 direct-q tanh-sinh at 35 and 50 digits.
    // Ordinary complex C/S, separately split critical points and lateral
    // oscillations: no production transforms, scaling or GK samples reused.
    // Refinement differences: pressure <5e-16, velocity <5e-19 absolute.
    let cases = [
        (
            [fluid(1.2, 300.0), fluid(2.4, 600.0), fluid(1.8, 400.0)],
            700.0,
            0.025,
            0.018,
            0.023,
            0.04,
            [
                C64::new(1066.0835825391047, -1695.0660938541113),
                C64::new(4.750322168227399, 0.08008661562454795),
                C64::new(-5.521433239708101, 1.2992113722244687),
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
                C64::new(-1685.4152679511676, -2525.2156809146957),
                C64::new(-0.7431756243168111, -1.6834602495061538),
                C64::new(0.16135936449573414, 1.6588071640414331),
            ],
        ),
    ];
    for (media, hz, a, h, b, r, expected) in cases {
        let observer = prepare(slab(media, h), hz, [0.0, 0.0, -a], [r, 0.0, -b]);
        let actual = observer.evaluate_reflected(C64::ONE).unwrap();
        close(actual.pressure_pa, expected[0], 4e-8);
        close(actual.particle_velocity_m_s[0], expected[1], 4e-8);
        close(actual.particle_velocity_m_s[2], expected[2], 4e-8);
        assert_eq!(actual.particle_velocity_m_s[1], C64::ZERO);
        assert!(!observer.diagnostics().reduced_to_single_interface);
        assert!(observer.diagnostics().evaluations > 0);
        assert!(observer.diagnostics().evaluations <= limits().maximum_evaluations);
        let line = observer.evaluate_reflected_harmonic_line(C64::ONE).unwrap();
        assert_eq!(line.frequency_hz, hz);
        assert_eq!(line.pressure_pa, actual.pressure_pa);
        assert_eq!(line.particle_velocity_m_s, actual.particle_velocity_m_s);
        let velocity = actual.particle_velocity_m_s;
        for axis in 0..3 {
            let intensity = actual.pressure_pa.scale(0.5) * velocity[axis].conj();
            assert_eq!(intensity.re, actual.active_intensity_w_m2[axis]);
            assert_eq!(intensity.im, actual.reactive_intensity_w_m2[axis]);
        }
    }
}

#[test]
fn absent_matched_and_translated_original_interfaces_keep_their_fields() {
    let m1 = fluid(1.2, 300.0);
    let m3 = fluid(2.4, 450.0);
    let source = [-0.01, 0.02, -0.025];
    let receiver = [0.035, -0.01, -0.041];
    let q = C64::new(0.003, -0.001);
    let interface = PlanarFluidInterface::new(m1, m3, [0.0; 3], [0.0, 0.0, 1.0]).unwrap();
    let original = PreparedReflectedFullWaveMonopole::new(
        interface,
        700.0,
        source,
        receiver,
        limits(),
        &mut || false,
    )
    .unwrap();
    for layer in [
        slab([m1, fluid(0.8, 200.0), m3], 0.0),
        slab([m1, m3, m3], 0.018),
    ] {
        let observer = prepare(layer, 700.0, source, receiver);
        for total in [false, true] {
            let actual = if total {
                observer.evaluate_total(q)
            } else {
                observer.evaluate_reflected(q)
            }
            .unwrap();
            let expected = if total {
                original.evaluate_total(q)
            } else {
                original.evaluate_reflected(q)
            }
            .unwrap();
            assert_eq!(actual.pressure_pa, expected.pressure_pa);
            assert_eq!(actual.particle_velocity_m_s, expected.particle_velocity_m_s);
        }
        assert!(observer.diagnostics().reduced_to_single_interface);
        assert_eq!(
            observer.diagnostics().evaluations,
            original.diagnostics().evaluations
        );
    }
    let h = 0.018;
    let matched_incident = prepare(slab([m1, m1, m3], h), 700.0, source, receiver);
    let shifted = PlanarFluidInterface::new(m1, m3, [0.0, 0.0, h], [0.0, 0.0, 1.0]).unwrap();
    let shifted = PreparedReflectedFullWaveMonopole::new(
        shifted,
        700.0,
        source,
        receiver,
        limits(),
        &mut || false,
    )
    .unwrap();
    let actual = matched_incident.evaluate_reflected(q).unwrap();
    let expected = shifted.evaluate_reflected(q).unwrap();
    close(actual.pressure_pa, expected.pressure_pa, 4e-8);
    for axis in 0..3 {
        close(
            actual.particle_velocity_m_s[axis],
            expected.particle_velocity_m_s[axis],
            4e-8,
        );
    }
    let homogeneous = prepare(slab([m1; 3], h), 700.0, source, receiver);
    assert_eq!(
        homogeneous.evaluate_reflected(q).unwrap().pressure_pa,
        C64::ZERO
    );
    assert_eq!(
        homogeneous
            .evaluate_reflected(q)
            .unwrap()
            .particle_velocity_m_s,
        [C64::ZERO; 3]
    );
    assert_eq!(homogeneous.diagnostics().evaluations, 0);
    let direct = PreparedCompactSourceObserver::new(m1, 700.0, receiver, &[source], 1)
        .unwrap()
        .evaluate(&[CompactSourcePhasor::monopole(q)])
        .unwrap();
    let actual = homogeneous.evaluate_total(q).unwrap();
    assert_eq!(actual.pressure_pa, direct.pressure_pa);
    assert_eq!(actual.particle_velocity_m_s, direct.particle_velocity_m_s);
}

#[test]
fn common_speed_density_layer_matches_coherent_image_train_and_total_interference() {
    let (r1, r2, r3, c) = (1.2, 5.0, 2.4, 350.0);
    let (a, h, b, r, hz) = (0.03, 0.08, 0.04, 0.06, 1100.0);
    let media = [fluid(r1, c), fluid(r2, c), fluid(r3, c)];
    let source = [0.0, 0.0, -a];
    let receiver = [r, 0.0, -b];
    let layer = prepare(slab(media, h), hz, source, receiver);
    let r12 = (r2 - r1) / (r2 + r1);
    let r23 = (r3 - r2) / (r3 + r2);
    let mut positions = Vec::new();
    let mut amplitudes = Vec::new();
    let mut returning = (1.0 - r12 * r12) * r23;
    for m in 0..80 {
        let coefficient = if m == 0 { r12 } else { returning };
        positions.push([0.0, 0.0, a + 2.0 * (m as f64) * h]);
        amplitudes.push(CompactSourcePhasor::monopole(C64::from_re(coefficient)));
        if m > 0 {
            returning *= -r12 * r23;
        }
    }
    assert!(returning.abs() < 1e-50);
    let expected = PreparedCompactSourceObserver::new(media[0], hz, receiver, &positions, 80)
        .unwrap()
        .evaluate(&amplitudes)
        .unwrap();
    let reflected = layer.evaluate_reflected(C64::ONE).unwrap();
    close(reflected.pressure_pa, expected.pressure_pa, 4e-8);
    for axis in 0..3 {
        close(
            reflected.particle_velocity_m_s[axis],
            expected.particle_velocity_m_s[axis],
            4e-8,
        );
    }
    positions.push(source);
    amplitudes.push(CompactSourcePhasor::monopole(C64::ONE));
    let expected = PreparedCompactSourceObserver::new(media[0], hz, receiver, &positions, 81)
        .unwrap()
        .evaluate(&amplitudes)
        .unwrap();
    let total = layer.evaluate_total(C64::ONE).unwrap();
    close(total.pressure_pa, expected.pressure_pa, 4e-8);
    for axis in 0..3 {
        close(
            total.particle_velocity_m_s[axis],
            expected.particle_velocity_m_s[axis],
            4e-8,
        );
        close(
            C64::new(
                total.active_intensity_w_m2[axis],
                total.reactive_intensity_w_m2[axis],
            ),
            C64::new(
                expected.active_intensity_w_m2[axis],
                expected.reactive_intensity_w_m2[axis],
            ),
            8e-8,
        );
    }
    assert!(layer.diagnostics().evaluations > 0);
}

#[test]
fn momentum_reciprocity_finite_self_field_and_bounded_refusals() {
    let media = [fluid(1.2, 300.0), fluid(2.4, 600.0), fluid(1.8, 400.0)];
    let model = slab(media, 0.018);
    let hz = 700.0;
    let source = [0.0, 0.0, -0.025];
    let receiver = [0.04, 0.0, -0.023];
    let prepared = prepare(model, hz, source, receiver);
    let actual = prepared.evaluate_reflected(C64::ONE).unwrap();
    let reversed = prepare(model, hz, receiver, source)
        .evaluate_reflected(C64::ONE)
        .unwrap();
    close(actual.pressure_pa, reversed.pressure_pa, 2e-12);
    for axis in [0, 2] {
        let step = 1e-6;
        let mut plus = receiver;
        let mut minus = receiver;
        plus[axis] += step;
        minus[axis] -= step;
        let dp = (prepare(model, hz, source, plus)
            .evaluate_reflected(C64::ONE)
            .unwrap()
            .pressure_pa
            - prepare(model, hz, source, minus)
                .evaluate_reflected(C64::ONE)
                .unwrap()
                .pressure_pa)
            .scale(0.5 / step);
        let velocity = dp * C64::new(0.0, -1.0 / (TAU * hz * media[0].density_kg_m3()));
        close(actual.particle_velocity_m_s[axis], velocity, 2e-6);
    }
    let at_source = prepare(model, hz, source, source);
    assert!(!at_source.total_field_available());
    assert!(
        at_source
            .evaluate_reflected(C64::ONE)
            .unwrap()
            .pressure_pa
            .abs()
            .is_finite()
    );
    assert!(at_source.evaluate_total(C64::ZERO).is_err());
    assert!(at_source.evaluate_total_harmonic_line(C64::ONE).is_err());
    assert!(prepared.total_field_available());
    assert_eq!(
        prepared.evaluate_reflected(C64::ZERO).unwrap().pressure_pa,
        C64::ZERO
    );
    assert!(
        prepared
            .evaluate_reflected(C64::new(f64::NAN, 0.0))
            .is_err()
    );
    assert_eq!(prepared.source_position_m(), source);
    assert_eq!(prepared.observer_position_m(), receiver);
    assert_eq!(prepared.frequency_hz(), hz);
    assert_eq!(prepared.layer().thickness_m(), 0.018);
    let total_line = prepared.evaluate_total_harmonic_line(C64::ONE).unwrap();
    assert_eq!(total_line.frequency_hz, hz);
    assert_eq!(
        total_line.pressure_pa,
        prepared.evaluate_total(C64::ONE).unwrap().pressure_pa
    );
    for invalid in [[0.0, 0.0, 0.0], [0.0, 0.0, 0.02], [f64::NAN, 0.0, -0.02]] {
        assert!(
            PreparedLayeredReflectedFullWaveMonopole::new(
                model,
                hz,
                source,
                invalid,
                limits(),
                &mut || false
            )
            .is_err()
        );
    }
    let mut budget = limits();
    budget.maximum_evaluations = 1;
    assert!(
        PreparedLayeredReflectedFullWaveMonopole::new(
            model,
            hz,
            source,
            receiver,
            budget,
            &mut || false
        )
        .is_err()
    );
    assert!(matches!(
        PreparedLayeredReflectedFullWaveMonopole::new(
            model,
            hz,
            source,
            receiver,
            limits(),
            &mut || true
        ),
        Err(AeroacError::Cancelled { .. })
    ));
    assert!(
        PlanarFluidLayer::new(
            media[0],
            fluid(1.0, 200.0),
            media[2],
            [0.0; 3],
            [0.0, 0.0, 1.0],
            0.1
        )
        .is_err()
    );
}
