//! G1/G3: exact near fields, independent Sommerfeld references and observations.

use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::refraction::full_wave::{FullWaveMonopoleLimits, PreparedFullWaveMonopole};
use fs_aeroac::curle3d::refraction::{
    PlanarFluidInterface, PreparedRefractedMonopole, RefractedMonopoleLimits,
};
use fs_aeroac::curle3d::{CompactSourcePhasor, PreparedCompactSourceObserver};
use fs_aeroac::observation::{HarmonicObservationBudget, PreparedHarmonicField};
use fs_fft::harmonic_window::HarmonicWindow;
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::c64::C64;
use std::f64::consts::TAU;

fn fluid(rho: f64, speed: f64) -> FluidAcousticState {
    FluidAcousticState::try_new(rho, speed).unwrap()
}
fn limits() -> FullWaveMonopoleLimits {
    FullWaveMonopoleLimits {
        absolute_pressure_tolerance_pa_per_m3_s: 1e-7,
        absolute_velocity_tolerance_per_m2: 1e-10,
        relative_tolerance: 1e-10,
        maximum_evaluations: 65_536,
        maximum_depth: 32,
    }
}
fn plane(media: [FluidAcousticState; 2]) -> PlanarFluidInterface {
    PlanarFluidInterface::new(media[0], media[1], [0.0; 3], [0.0, 0.0, 1.0]).unwrap()
}
fn prepare(
    media: [FluidAcousticState; 2],
    source: [f64; 3],
    receiver: [f64; 3],
    hz: f64,
) -> PreparedFullWaveMonopole {
    PreparedFullWaveMonopole::new(plane(media), hz, source, receiver, limits(), &mut || false)
        .unwrap()
}
fn close(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-12),
        "actual={actual:?}, expected={expected:?}, relative tolerance={tolerance:e}"
    );
}

#[test]
fn equal_speed_is_exact_including_near_field_and_receiving_density() {
    let source = [-0.01, 0.005, -0.015];
    let receiver = [0.02, 0.025, 0.01];
    let gas = fluid(1.2, 320.0);
    let frequency = 500.0;
    let q = C64::new(0.002, -0.001);
    let reference = PreparedCompactSourceObserver::new(gas, frequency, receiver, &[source], 1)
        .unwrap()
        .evaluate(&[CompactSourcePhasor::monopole(q)])
        .unwrap();
    for rho2 in [1.2, 4.8] {
        let observer = prepare([gas, fluid(rho2, 320.0)], source, receiver, frequency);
        let field = observer.evaluate(q).unwrap();
        let pressure_factor = 2.0 * rho2 / (1.2 + rho2);
        let velocity_factor = 2.0 * 1.2 / (1.2 + rho2);
        close(
            field.pressure_pa,
            reference.pressure_pa.scale(pressure_factor),
            3e-14,
        );
        for axis in 0..3 {
            close(
                field.particle_velocity_m_s[axis],
                reference.particle_velocity_m_s[axis].scale(velocity_factor),
                3e-14,
            );
        }
        assert!(field.reactive_intensity_w_m2[2].abs() > field.active_intensity_w_m2[2].abs());
        assert!(observer.diagnostics().exact_equal_speed);
        assert_eq!(observer.diagnostics().evaluations, 0);
        assert_eq!(observer.interface().receiving_state(), fluid(rho2, 320.0));
        assert_eq!(observer.source_position_m(), source);
        assert_eq!(observer.observer_position_m(), receiver);
    }
    // Near-equal original speeds must still use their distinct spectra.
    let near = prepare([gas, fluid(4.8, 320.000_01)], source, receiver, frequency);
    assert!(!near.diagnostics().exact_equal_speed);
    assert!(near.diagnostics().evaluations > 0);
    close(
        near.evaluate(q).unwrap().pressure_pa,
        reference.pressure_pa.scale(1.6),
        1e-7,
    );
}

#[test]
fn full_wave_matches_independent_high_precision_near_and_critical_spectra() {
    // Independent mpmath 1.3.0 direct-q tanh-sinh quadrature at 35 and 50
    // decimal digits, split at k1,k2 and oscillation-resolved evanescent tails.
    // It does not use the production branch substitutions or GK panel rule.
    // Pressure changes between precisions were <1e-32, velocities <3e-35.
    let cases = [
        (
            [fluid(1.2, 300.0), fluid(2.4, 450.0)],
            700.0,
            0.025,
            0.018,
            0.035,
            C64::new(6598.877_704_006_922, -7519.932_291_762_797),
            C64::new(13.341_760_979_269_208, 1.330_183_142_541_795),
            C64::new(16.187_025_631_170_134, 3.599_986_748_220_998_6),
            -6778.373_423_398_142,
        ),
        (
            [fluid(2.4, 600.0), fluid(1.2, 300.0)],
            1500.0,
            0.04,
            0.06,
            0.12,
            C64::new(-4918.510_262_256_88, 4111.726_119_948_162),
            C64::new(-11.169_383_935_675_89, 4.104_322_477_105_612),
            C64::new(-10.123_606_344_340_23, 7.873_620_695_991_724),
            883.429_804_583_029_7,
        ),
    ];
    for (media, hz, a, b, r, p, ur, un, evanescent_imaginary_pressure) in cases {
        let observer = prepare(media, [0.0, 0.0, -a], [r, 0.0, b], hz);
        let field = observer.evaluate(C64::ONE).unwrap();
        close(field.pressure_pa, p, 3e-10);
        close(field.particle_velocity_m_s[0], ur, 3e-10);
        close(field.particle_velocity_m_s[2], un, 3e-10);
        assert_eq!(field.particle_velocity_m_s[1], C64::ZERO);
        let without_both_evanescent = p - C64::new(0.0, evanescent_imaginary_pressure);
        assert!(
            (field.pressure_pa - without_both_evanescent).abs() > 0.1 * p.abs(),
            "the physically important high-q spherical source tail must survive"
        );
        assert!(observer.diagnostics().maximum_bessel_argument < 3000.0);
    }
}

#[test]
fn reciprocal_transmission_and_rotated_geometry_preserve_the_original_field() {
    let media = [fluid(1.2, 300.0), fluid(4.3, 510.0)];
    let source = [-0.03, 0.04, -0.02];
    let receiver = [0.07, -0.05, 0.03];
    let hz = 1900.0;
    let forward = prepare(media, source, receiver, hz);
    let reverse_plane =
        PlanarFluidInterface::new(media[1], media[0], [0.0; 3], [0.0, 0.0, -1.0]).unwrap();
    let reverse =
        PreparedFullWaveMonopole::new(reverse_plane, hz, receiver, source, limits(), &mut || false)
            .unwrap();
    let field = forward.evaluate(C64::ONE).unwrap();
    close(
        reverse.evaluate(C64::ONE).unwrap().pressure_pa,
        field.pressure_pa,
        4e-10,
    );
    let transform = |p: [f64; 3]| [1.0 + p[2], -2.0 + p[0], 3.0 + p[1]];
    let rotated_plane =
        PlanarFluidInterface::new(media[0], media[1], transform([0.0; 3]), [7.0, 0.0, 0.0])
            .unwrap();
    let rotated = PreparedFullWaveMonopole::new(
        rotated_plane,
        hz,
        transform(source),
        transform(receiver),
        limits(),
        &mut || false,
    )
    .unwrap();
    let other = rotated.evaluate(C64::ONE).unwrap();
    close(other.pressure_pa, field.pressure_pa, 4e-10);
    for (new, old) in [(0, 2), (1, 0), (2, 1)] {
        close(
            other.particle_velocity_m_s[new],
            field.particle_velocity_m_s[old],
            4e-10,
        );
    }
}

#[test]
fn velocity_is_the_spatial_gradient_of_the_same_transmitted_pressure() {
    let media = [fluid(1.2, 300.0), fluid(2.4, 450.0)];
    let source = [-0.01, 0.005, -0.025];
    let receiver = [0.035, -0.007, 0.018];
    let hz = 700.0;
    let tight = FullWaveMonopoleLimits {
        absolute_pressure_tolerance_pa_per_m3_s: 1e-9,
        absolute_velocity_tolerance_per_m2: 1e-12,
        relative_tolerance: 1e-12,
        ..limits()
    };
    let observe = |position| {
        PreparedFullWaveMonopole::new(plane(media), hz, source, position, tight, &mut || false)
            .unwrap()
            .evaluate(C64::ONE)
            .unwrap()
    };
    let field = observe(receiver);
    let step = 1e-5;
    for axis in 0..3 {
        let mut left = receiver;
        let mut right = receiver;
        left[axis] -= step;
        right[axis] += step;
        let derivative = (observe(right).pressure_pa - observe(left).pressure_pa).scale(0.5 / step);
        let from_pressure = derivative / C64::new(0.0, TAU * hz * media[1].density_kg_m3());
        close(field.particle_velocity_m_s[axis], from_pressure, 5e-7);
    }
}

#[test]
fn on_axis_symmetry_and_far_field_approach_the_refracted_ray_tube() {
    let media = [fluid(1.2, 300.0), fluid(2.4, 450.0)];
    let source = [0.0, 0.0, -10.0];
    let receiver = [0.0, 0.0, 15.0];
    let hz = 1000.0;
    let wave = prepare(media, source, receiver, hz);
    let ray = PreparedRefractedMonopole::new(
        plane(media),
        hz,
        source,
        receiver,
        RefractedMonopoleLimits {
            minimum_kh: 100.0,
            minimum_direction_cosine: 0.5,
            path_tolerance_m: 1e-10,
            maximum_iterations: 64,
        },
        &mut || false,
    )
    .unwrap();
    let full = wave.evaluate(C64::ONE).unwrap();
    let leading = ray.evaluate(C64::ONE).unwrap();
    assert_eq!(full.particle_velocity_m_s[0], C64::ZERO);
    assert_eq!(full.particle_velocity_m_s[1], C64::ZERO);
    assert_eq!(wave.diagnostics().velocity_tail_bound_per_m2[0], 0.0);
    close(full.pressure_pa, leading.pressure_pa, 0.01);
    close(
        full.particle_velocity_m_s[2],
        leading.particle_velocity_m_s[2],
        0.01,
    );
}

#[test]
fn refinement_controls_tail_and_finite_quadrature_separately() {
    let media = [fluid(1.2, 300.0), fluid(2.4, 450.0)];
    let source = [0.0, 0.0, -0.025];
    let receiver = [0.035, 0.0, 0.018];
    let coarse_limits = FullWaveMonopoleLimits {
        absolute_pressure_tolerance_pa_per_m3_s: 1e-2,
        absolute_velocity_tolerance_per_m2: 1e-5,
        relative_tolerance: 0.0,
        ..limits()
    };
    let coarse = PreparedFullWaveMonopole::new(
        plane(media),
        700.0,
        source,
        receiver,
        coarse_limits,
        &mut || false,
    )
    .unwrap();
    let fine = prepare(media, source, receiver, 700.0);
    let c = coarse.diagnostics();
    let f = fine.diagnostics();
    assert!(f.transverse_wavenumber_cutoff_rad_m > c.transverse_wavenumber_cutoff_rad_m);
    assert!(f.pressure_tail_bound_pa_per_m3_s < c.pressure_tail_bound_pa_per_m3_s * 1e-3);
    assert!(c.pressure_tail_bound_pa_per_m3_s > 0.0);
    assert!(c.pressure_quadrature_error_pa_per_m3_s > 0.0);
    let fc = coarse.evaluate(C64::ONE).unwrap();
    let ff = fine.evaluate(C64::ONE).unwrap();
    assert!(
        (fc.pressure_pa - ff.pressure_pa).abs()
            < coarse_limits.absolute_pressure_tolerance_pa_per_m3_s
    );
    for axis in 0..3 {
        assert!(
            (fc.particle_velocity_m_s[axis] - ff.particle_velocity_m_s[axis]).abs()
                < coarse_limits.absolute_velocity_tolerance_per_m2
        );
    }
}

#[test]
fn prepared_source_clock_feeds_existing_physical_microphone_windows() {
    let observer = prepare(
        [fluid(1.2, 300.0), fluid(2.4, 450.0)],
        [0.0, 0.0, -0.025],
        [0.035, 0.0, 0.018],
        700.0,
    );
    let q = C64::new(0.002, -0.001);
    let line = observer.evaluate_harmonic_line(q).unwrap();
    let field = observer.evaluate(q).unwrap();
    assert_eq!(line.pressure_pa, field.pressure_pa);
    assert_eq!(line.particle_velocity_m_s, field.particle_velocity_m_s);
    let observed = PreparedHarmonicField::new(
        &[line],
        HarmonicObservationBudget {
            max_lines: 1,
            max_line_pairs: 1,
        },
    )
    .unwrap()
    .observe_window(HarmonicWindow::new(0.137, 1.0 / observer.frequency_hz()).unwrap())
    .unwrap();
    let expected = 0.5 * field.pressure_pa.norm_sq();
    assert!((observed.mean_square_pressure_pa2 - expected).abs() < 3e-12 * expected);
    for axis in [0, 2] {
        assert!(
            (observed.mean_intensity_w_m2[axis] - field.active_intensity_w_m2[axis]).abs()
                < 3e-12 * field.active_intensity_w_m2[axis].abs()
        );
    }
}

#[test]
fn invalid_geometry_work_exhaustion_and_mid_quadrature_cancellation_refuse() {
    let interface = plane([fluid(1.2, 300.0), fluid(2.4, 450.0)]);
    let source = [0.0, 0.0, -0.025];
    let receiver = [0.035, 0.0, 0.018];
    let attempt = |hz, s, r, budget| {
        PreparedFullWaveMonopole::new(interface, hz, s, r, budget, &mut || false)
    };
    for bad in [
        FullWaveMonopoleLimits {
            absolute_pressure_tolerance_pa_per_m3_s: 0.0,
            ..limits()
        },
        FullWaveMonopoleLimits {
            absolute_velocity_tolerance_per_m2: f64::NAN,
            ..limits()
        },
        FullWaveMonopoleLimits {
            relative_tolerance: 1.0,
            ..limits()
        },
        FullWaveMonopoleLimits {
            maximum_depth: 49,
            ..limits()
        },
        FullWaveMonopoleLimits {
            maximum_evaluations: 0,
            ..limits()
        },
        FullWaveMonopoleLimits {
            maximum_evaluations: 30,
            ..limits()
        },
    ] {
        assert!(attempt(700.0, source, receiver, bad).is_err());
    }
    assert!(attempt(0.0, source, receiver, limits()).is_err());
    assert!(attempt(700.0, receiver, source, limits()).is_err());
    assert!(attempt(700.0, source, [0.0; 3], limits()).is_err());
    assert!(matches!(
        attempt(
            700.0,
            source,
            receiver,
            FullWaveMonopoleLimits {
                maximum_depth: 1,
                ..limits()
            }
        ),
        Err(AeroacError::NotConverged { .. })
    ));
    let mut polls = 0;
    assert!(matches!(
        PreparedFullWaveMonopole::new(interface, 700.0, source, receiver, limits(), &mut || {
            polls += 1;
            polls == 100
        }),
        Err(AeroacError::Cancelled { .. })
    ));
    let observer = attempt(700.0, source, receiver, limits()).unwrap();
    assert!(observer.evaluate(C64::new(f64::INFINITY, 0.0)).is_err());
    assert_eq!(observer.evaluate(C64::ZERO).unwrap().pressure_pa, C64::ZERO);
    let replay = attempt(700.0, source, receiver, limits()).unwrap();
    assert_eq!(
        observer.evaluate(C64::ONE).unwrap().pressure_pa,
        replay.evaluate(C64::ONE).unwrap().pressure_pa
    );
}
