//! G1/G3: coherent incident-side fields, critical reflection and original clocks.

use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::refraction::PlanarFluidInterface;
use fs_aeroac::curle3d::refraction::full_wave::reflection::PreparedReflectedFullWaveMonopole;
use fs_aeroac::curle3d::refraction::full_wave::{FullWaveMonopoleLimits, PreparedFullWaveMonopole};
use fs_aeroac::curle3d::{CompactSourcePhasor, PreparedCompactSourceObserver};
use fs_aeroac::observation::{HarmonicObservationBudget, PreparedHarmonicField};
use fs_fft::harmonic_window::HarmonicWindow;
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::c64::C64;
use std::f64::consts::TAU;

fn fluid(rho: f64, c: f64) -> FluidAcousticState {
    FluidAcousticState::try_new(rho, c).unwrap()
}
fn plane(media: [FluidAcousticState; 2]) -> PlanarFluidInterface {
    PlanarFluidInterface::new(media[0], media[1], [0.0; 3], [0.0, 0.0, 1.0]).unwrap()
}
fn limits() -> FullWaveMonopoleLimits {
    FullWaveMonopoleLimits {
        absolute_pressure_tolerance_pa_per_m3_s: 1e-8,
        absolute_velocity_tolerance_per_m2: 1e-11,
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
) -> PreparedReflectedFullWaveMonopole {
    PreparedReflectedFullWaveMonopole::new(interface, hz, source, observer, limits(), &mut || false)
        .unwrap()
}
fn close(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-12),
        "actual={actual:?}, expected={expected:?}, tolerance={tolerance:e}"
    );
}

#[test]
fn equal_speed_images_preserve_induction_and_coherent_direct_reflected_interference() {
    let gas = fluid(1.2, 343.0);
    let source = [-0.01, 0.005, -0.015];
    let image = [source[0], source[1], -source[2]];
    let receiver = [0.02, 0.025, -0.01];
    let hz = 700.0;
    let q = C64::new(0.002, -0.001);
    let original =
        PreparedCompactSourceObserver::new(gas, hz, receiver, &[source, image], 2).unwrap();
    let direct = original
        .evaluate(&[
            CompactSourcePhasor::monopole(q),
            CompactSourcePhasor::default(),
        ])
        .unwrap();
    for rho2 in [0.3, 1.2, 4.8] {
        let r = (rho2 - 1.2) / (rho2 + 1.2);
        let interface = plane([gas, fluid(rho2, 343.0)]);
        let prepared = prepare(interface, hz, source, receiver);
        let reflected = prepared.evaluate_reflected(q).unwrap();
        let total = prepared.evaluate_total(q).unwrap();
        let expected_reflected = original
            .evaluate(&[
                CompactSourcePhasor::default(),
                CompactSourcePhasor::monopole(q.scale(r)),
            ])
            .unwrap();
        let expected_total = original
            .evaluate(&[
                CompactSourcePhasor::monopole(q),
                CompactSourcePhasor::monopole(q.scale(r)),
            ])
            .unwrap();
        close(reflected.pressure_pa, expected_reflected.pressure_pa, 3e-13);
        close(total.pressure_pa, expected_total.pressure_pa, 3e-13);
        for axis in 0..3 {
            close(
                reflected.particle_velocity_m_s[axis],
                expected_reflected.particle_velocity_m_s[axis],
                3e-13,
            );
            close(
                total.particle_velocity_m_s[axis],
                expected_total.particle_velocity_m_s[axis],
                3e-13,
            );
            let expected_flux = expected_total.pressure_pa.scale(0.5)
                * expected_total.particle_velocity_m_s[axis].conj();
            close(
                C64::new(
                    total.active_intensity_w_m2[axis],
                    total.reactive_intensity_w_m2[axis],
                ),
                expected_flux,
                6e-13,
            );
        }
        if r != 0.0 {
            assert!(
                (total.active_intensity_w_m2[0]
                    - direct.active_intensity_w_m2[0]
                    - reflected.active_intensity_w_m2[0])
                    .abs()
                    > 1e-9
            );
        } else {
            assert_eq!(reflected.pressure_pa, C64::ZERO);
            assert_eq!(reflected.particle_velocity_m_s, [C64::ZERO; 3]);
        }
        assert!(prepared.diagnostics().exact_equal_speed);
        assert_eq!(prepared.diagnostics().evaluations, 0);
        assert_eq!(
            prepared.interface().receiving_state(),
            interface.receiving_state()
        );
        assert_eq!(prepared.source_position_m(), source);
        assert_eq!(prepared.observer_position_m(), receiver);
        assert!(prepared.total_field_available());
    }
}

#[test]
fn original_reflected_spectra_match_independent_critical_and_evanescent_references() {
    // Direct-q mpmath tanh-sinh at 35/50 digits, split at both original critical
    // points. No production branch map, log normalizer or GK nodes are reused.
    let cases = [
        (
            [fluid(1.2, 300.0), fluid(2.4, 600.0)],
            700.0,
            0.025,
            0.023,
            0.04,
            [
                C64::new(1642.7391157234595, -2112.2708350825457),
                C64::new(5.49280389483676, 0.18714327914220203),
                C64::new(-7.5010642986821575, 1.5210671685188152),
            ],
        ),
        (
            [fluid(2.4, 600.0), fluid(1.2, 300.0)],
            1500.0,
            0.04,
            0.06,
            0.12,
            [
                C64::new(-6057.50764907886, -4821.342708846944),
                C64::new(-2.7693821678200504, -3.2719312014191795),
                C64::new(1.9180828231501152, 3.9989741181713915),
            ],
        ),
        (
            [fluid(1.2, 343.0), fluid(2.4, 343.01)],
            1000.0,
            0.03,
            0.05,
            0.07,
            [
                C64::new(1749.7409903669735, 691.6633916926491),
                C64::new(2.2310390405505446, 2.544095296107232),
                C64::new(-2.549850443244649, -2.9075820070279734),
            ],
        ),
    ];
    for (media, hz, a, b, r, expected) in cases {
        let prepared = prepare(plane(media), hz, [0.0, 0.0, -a], [r, 0.0, -b]);
        let field = prepared.evaluate_reflected(C64::ONE).unwrap();
        close(field.pressure_pa, expected[0], 3e-8);
        close(field.particle_velocity_m_s[0], expected[1], 3e-8);
        close(field.particle_velocity_m_s[2], expected[2], 3e-8);
        assert_eq!(field.particle_velocity_m_s[1], C64::ZERO);
        assert!(!prepared.diagnostics().exact_equal_speed);
        assert!(prepared.diagnostics().evaluations > 0);
        assert!(
            prepared.diagnostics().pressure_tail_bound_pa_per_m3_s
                <= limits().absolute_pressure_tolerance_pa_per_m3_s * 0.25
        );
        let tail = prepared.diagnostics().velocity_tail_bound_per_m2;
        assert!(tail[0].hypot(tail[1]) <= limits().absolute_velocity_tolerance_per_m2 * 0.25);
    }
}

#[test]
fn continuity_momentum_reciprocity_and_rotation_join_the_original_transmission_owner() {
    let hz = 700.0;
    let source = [0.0, 0.0, -0.04];
    for media in [
        [fluid(1.2, 300.0), fluid(2.4, 600.0)],
        [fluid(2.4, 600.0), fluid(1.2, 300.0)],
    ] {
        let interface = plane(media);
        let epsilon = 1e-7;
        let total = prepare(interface, hz, source, [0.06, 0.0, -epsilon])
            .evaluate_total(C64::ONE)
            .unwrap();
        let through = PreparedFullWaveMonopole::new(
            interface,
            hz,
            source,
            [0.06, 0.0, epsilon],
            limits(),
            &mut || false,
        )
        .unwrap()
        .evaluate(C64::ONE)
        .unwrap();
        close(total.pressure_pa, through.pressure_pa, 3e-5);
        close(
            total.particle_velocity_m_s[2],
            through.particle_velocity_m_s[2],
            3e-5,
        );
        close(
            total.particle_velocity_m_s[0].scale(media[0].density_kg_m3()),
            through.particle_velocity_m_s[0].scale(media[1].density_kg_m3()),
            3e-5,
        );
    }
    let interface = plane([fluid(1.2, 300.0), fluid(2.4, 600.0)]);
    let receiver = [0.06, 0.0, -0.03];
    let prepared = prepare(interface, hz, source, receiver);
    let reflected = prepared.evaluate_reflected(C64::ONE).unwrap();
    let total = prepared.evaluate_total(C64::ONE).unwrap();
    let reverse = prepare(interface, hz, receiver, source);
    close(
        reverse.evaluate_reflected(C64::ONE).unwrap().pressure_pa,
        reflected.pressure_pa,
        3e-13,
    );
    close(
        reverse.evaluate_total(C64::ONE).unwrap().pressure_pa,
        total.pressure_pa,
        3e-13,
    );
    for axis in [0, 2] {
        let step = 1e-6;
        let mut plus = receiver;
        let mut minus = receiver;
        plus[axis] += step;
        minus[axis] -= step;
        let p = prepare(interface, hz, source, plus);
        let m = prepare(interface, hz, source, minus);
        let momentum = C64::new(0.0, TAU * hz * 1.2 * 2.0 * step);
        close(
            (p.evaluate_reflected(C64::ONE).unwrap().pressure_pa
                - m.evaluate_reflected(C64::ONE).unwrap().pressure_pa)
                / momentum,
            reflected.particle_velocity_m_s[axis],
            2e-6,
        );
        close(
            (p.evaluate_total(C64::ONE).unwrap().pressure_pa
                - m.evaluate_total(C64::ONE).unwrap().pressure_pa)
                / momentum,
            total.particle_velocity_m_s[axis],
            2e-6,
        );
    }
    let rotate = |v: [f64; 3]| [v[2], v[0], v[1]];
    let rotated_interface = PlanarFluidInterface::new(
        interface.incident_state(),
        interface.receiving_state(),
        [0.0; 3],
        rotate([0.0, 0.0, 1.0]),
    )
    .unwrap();
    let rotated = prepare(rotated_interface, hz, rotate(source), rotate(receiver))
        .evaluate_total(C64::ONE)
        .unwrap();
    close(rotated.pressure_pa, total.pressure_pa, 2e-12);
    for axis in 0..3 {
        close(
            rotated.particle_velocity_m_s[axis],
            total.particle_velocity_m_s[(axis + 2) % 3],
            2e-12,
        );
    }
}

#[test]
fn finite_reflected_self_field_native_clock_and_explicit_refusals() {
    let interface = plane([fluid(1.2, 300.0), fluid(2.4, 600.0)]);
    let hz = 700.0;
    let source = [0.0, 0.0, -0.025];
    let receiver = [0.04, 0.0, -0.023];
    let self_field = prepare(interface, hz, source, source);
    assert!(!self_field.total_field_available());
    let reflection = self_field.evaluate_reflected(C64::ONE).unwrap();
    assert!(reflection.pressure_pa.abs().is_finite());
    assert!(reflection.pressure_pa.abs() > 0.0);
    assert_eq!(reflection.particle_velocity_m_s[0], C64::ZERO);
    assert_eq!(reflection.particle_velocity_m_s[1], C64::ZERO);
    assert!(self_field.evaluate_total(C64::ONE).is_err());
    assert!(self_field.evaluate_total_harmonic_line(C64::ZERO).is_err());
    let prepared = prepare(interface, hz, source, receiver);
    let q = C64::new(0.002, -0.001);
    let total = prepared.evaluate_total(q).unwrap();
    let line = prepared.evaluate_total_harmonic_line(q).unwrap();
    assert_eq!(line.pressure_pa, total.pressure_pa);
    assert_eq!(line.particle_velocity_m_s, total.particle_velocity_m_s);
    let reflected = prepared.evaluate_reflected(q).unwrap();
    assert_eq!(
        prepared
            .evaluate_reflected_harmonic_line(q)
            .unwrap()
            .pressure_pa,
        reflected.pressure_pa
    );
    let observed = PreparedHarmonicField::new(
        &[line],
        HarmonicObservationBudget {
            max_lines: 1,
            max_line_pairs: 1,
        },
    )
    .unwrap()
    .observe_window(HarmonicWindow::new(0.137, 1.0 / prepared.frequency_hz()).unwrap())
    .unwrap();
    let expected = 0.5 * total.pressure_pa.norm_sq();
    assert!((observed.mean_square_pressure_pa2 - expected).abs() < 3e-12 * expected);
    for axis in 0..3 {
        assert!(
            (observed.mean_intensity_w_m2[axis] - total.active_intensity_w_m2[axis]).abs()
                <= 3e-12 * total.active_intensity_w_m2[axis].abs().max(1e-12)
        );
    }
    assert_eq!(
        prepared.evaluate_total(C64::ZERO).unwrap().pressure_pa,
        C64::ZERO
    );
    assert!(
        prepared
            .evaluate_reflected(C64::new(f64::NAN, 0.0))
            .is_err()
    );
    for bad in [[0.0; 3], [0.0, 0.0, 0.01], [f64::INFINITY, 0.0, -0.01]] {
        assert!(
            PreparedReflectedFullWaveMonopole::new(
                interface,
                hz,
                source,
                bad,
                limits(),
                &mut || false
            )
            .is_err()
        );
    }
    for limit in [
        FullWaveMonopoleLimits {
            maximum_evaluations: 44,
            ..limits()
        },
        FullWaveMonopoleLimits {
            relative_tolerance: 1.0,
            ..limits()
        },
        FullWaveMonopoleLimits {
            maximum_depth: 0,
            ..limits()
        },
    ] {
        assert!(
            PreparedReflectedFullWaveMonopole::new(
                interface,
                hz,
                source,
                receiver,
                limit,
                &mut || false
            )
            .is_err()
        );
    }
    assert!(
        PreparedReflectedFullWaveMonopole::new(
            interface,
            0.0,
            source,
            receiver,
            limits(),
            &mut || false
        )
        .is_err()
    );
    let mut polls = 0;
    assert!(matches!(
        PreparedReflectedFullWaveMonopole::new(
            interface,
            hz,
            source,
            receiver,
            limits(),
            &mut || {
                polls += 1;
                polls > 12
            }
        ),
        Err(AeroacError::Cancelled { .. })
    ));
}
