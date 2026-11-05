//! Physical force-to-field refraction: source derivatives and vector near fields.

use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::refraction::PlanarFluidInterface;
use fs_aeroac::curle3d::refraction::full_wave::{
    FullWaveForceLimits, FullWaveMonopoleLimits, PreparedFullWaveForce, PreparedFullWaveMonopole,
};
use fs_aeroac::curle3d::{
    CoherentCompactField, CompactSourcePhasor, PreparedCompactSourceObserver,
};
use fs_aeroac::observation::{HarmonicObservationBudget, PreparedHarmonicField};
use fs_fft::harmonic_window::HarmonicWindow;
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::c64::C64;
use std::f64::consts::TAU;

fn fluid(rho: f64, speed: f64) -> FluidAcousticState {
    FluidAcousticState::try_new(rho, speed).unwrap()
}
fn plane(media: [FluidAcousticState; 2]) -> PlanarFluidInterface {
    PlanarFluidInterface::new(media[0], media[1], [0.0; 3], [0.0, 0.0, 1.0]).unwrap()
}
fn limits() -> FullWaveForceLimits {
    FullWaveForceLimits {
        absolute_pressure_tolerance_pa_per_n: 1e-10,
        absolute_velocity_tolerance_m_s_per_n: 1e-12,
        relative_tolerance: 1e-11,
        maximum_evaluations: 65_536,
        maximum_depth: 32,
    }
}
fn prepare(
    media: [FluidAcousticState; 2],
    source: [f64; 3],
    receiver: [f64; 3],
    hz: f64,
) -> PreparedFullWaveForce {
    PreparedFullWaveForce::new(plane(media), hz, source, receiver, limits(), &mut || false).unwrap()
}
fn close(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-12),
        "actual={actual:?}, expected={expected:?}, relative tolerance={tolerance:e}"
    );
}
fn components(field: &CoherentCompactField) -> [C64; 4] {
    let [ux, uy, uz] = field.particle_velocity_m_s;
    [field.pressure_pa, ux, uy, uz]
}
fn unit_force(axis: usize) -> [C64; 3] {
    core::array::from_fn(|i| if i == axis { C64::ONE } else { C64::ZERO })
}

#[test]
fn equal_speed_matches_complete_free_space_dipole_with_original_density_jump() {
    let source = [-0.01, 0.005, -0.015];
    let receiver = [0.02, 0.025, 0.01];
    let gas = fluid(1.2, 320.0);
    let hz = 500.0;
    let force = [C64::new(0.3, -0.1), C64::new(-0.2, 0.4), C64::new(0.1, 0.7)];
    let reference = PreparedCompactSourceObserver::new(gas, hz, receiver, &[source], 1)
        .unwrap()
        .evaluate(&[CompactSourcePhasor::force_dipole(force)])
        .unwrap();
    for rho2 in [1.2, 4.8] {
        let observer = PreparedFullWaveForce::new(
            plane([gas, fluid(rho2, 320.0)]),
            hz,
            source,
            receiver,
            FullWaveForceLimits {
                maximum_evaluations: 1,
                ..limits()
            },
            &mut || false,
        )
        .unwrap();
        let field = observer.evaluate(force).unwrap();
        close(
            field.pressure_pa,
            reference.pressure_pa.scale(2.0 * rho2 / (1.2 + rho2)),
            5e-14,
        );
        for axis in 0..3 {
            close(
                field.particle_velocity_m_s[axis],
                reference.particle_velocity_m_s[axis].scale(2.0 * 1.2 / (1.2 + rho2)),
                5e-14,
            );
        }
        assert!(observer.diagnostics().exact_equal_speed);
        assert_eq!(observer.diagnostics().evaluations, 0);
        assert_eq!(observer.source_position_m(), source);
        assert_eq!(observer.observer_position_m(), receiver);
        assert_eq!(observer.interface().receiving_state(), fluid(rho2, 320.0));
    }
    let near = prepare([gas, fluid(4.8, 320.000_01)], source, receiver, hz);
    assert!(!near.diagnostics().exact_equal_speed);
    assert!(near.diagnostics().evaluations > 0);
    close(
        near.evaluate(force).unwrap().pressure_pa,
        reference.pressure_pa.scale(1.6),
        1e-7,
    );
}

#[test]
fn independent_high_precision_spectrum_pins_every_force_tensor_channel() {
    // Independent mpmath 1.3.0 direct-q tanh-sinh integration at 35 and 50
    // decimal digits, split at both k values and oscillation-resolved tails.
    // No production branch transformation/GK rule is reused. The worst
    // precision-refinement change is <5.4e-35 in these seven SI coefficients.
    // Ordered p/Fx,p/Fz,ux/Fx,uy/Fy,ux/Fz,uz/Fx,uz/Fz; x radial and z normal.
    let cases = [
        (
            [fluid(1.2, 300.0), fluid(2.4, 450.0)],
            700.0,
            0.025,
            0.018,
            0.035,
            [
                C64::new(26.683_521_958_538_415, 2.660_366_285_083_59),
                C64::new(32.746_057_353_409_584, 1.048_789_043_349_652_7),
                C64::new(0.006_874_401_615_394_208, 0.005_102_761_451_954_21),
                C64::new(0.007_200_857_903_164_554, -0.072_224_735_013_648_32),
                C64::new(0.000_285_518_358_253_676_4, 0.093_314_650_306_768_4),
                C64::new(-0.001_739_838_799_484_127_7, 0.097_847_030_440_593_11),
                C64::new(0.004_645_536_719_890_781, 0.047_045_613_905_070_83),
            ],
        ),
        (
            [fluid(2.4, 600.0), fluid(1.2, 300.0)],
            1500.0,
            0.04,
            0.06,
            0.12,
            [
                C64::new(-5.584_691_967_837_945, 2.052_161_238_552_806),
                C64::new(-2.512_515_560_131_056_7, -1.195_184_003_693_876_3),
                C64::new(-0.010_196_139_852_918_536, 0.001_672_729_953_711_517_7),
                C64::new(0.001_512_090_764_524_376_5, 0.004_114_959_871_884_32),
                C64::new(-0.005_654_370_183_657_802, -0.006_618_914_865_884_216),
                C64::new(-0.010_959_929_516_859_454, 0.001_313_634_963_544_348_2),
                C64::new(-0.002_895_443_199_858_246, -0.000_288_436_771_985_338_85),
            ],
        ),
    ];
    for (media, hz, a, b, range, expected) in cases {
        let observer = prepare(media, [0.0, 0.0, -a], [range, 0.0, b], hz);
        let fx = observer.evaluate(unit_force(0)).unwrap();
        let fy = observer.evaluate(unit_force(1)).unwrap();
        let fz = observer.evaluate(unit_force(2)).unwrap();
        let actual = [
            fx.pressure_pa,
            fz.pressure_pa,
            fx.particle_velocity_m_s[0],
            fy.particle_velocity_m_s[1],
            fz.particle_velocity_m_s[0],
            fx.particle_velocity_m_s[2],
            fz.particle_velocity_m_s[2],
        ];
        for (value, reference) in actual.into_iter().zip(expected) {
            close(value, reference, 3e-10);
        }
        // A force perpendicular to the source-receiver plane has no pressure,
        // yet carries real transverse particle motion. Scalar impedance fails.
        assert_eq!(fy.pressure_pa, C64::ZERO);
        assert_eq!(fy.particle_velocity_m_s[0], C64::ZERO);
        assert_eq!(fy.particle_velocity_m_s[2], C64::ZERO);
        assert!(fy.particle_velocity_m_s[1].abs() > 0.001);
        assert!((actual[4] - actual[5]).abs() > 0.001);
        assert!(observer.diagnostics().maximum_bessel_argument < 3000.0);
    }
}

#[test]
fn physical_force_is_the_analytic_source_derivative_of_the_same_volume_field() {
    let media = [fluid(1.2, 300.0), fluid(2.4, 450.0)];
    let source = [-0.01, 0.005, -0.025];
    let receiver = [0.035, -0.007, 0.018];
    let hz = 700.0;
    let force = prepare(media, source, receiver, hz);
    let source_field = |position| {
        PreparedFullWaveMonopole::new(
            plane(media),
            hz,
            position,
            receiver,
            FullWaveMonopoleLimits {
                absolute_pressure_tolerance_pa_per_m3_s: 1e-10,
                absolute_velocity_tolerance_per_m2: 1e-13,
                relative_tolerance: 1e-12,
                maximum_evaluations: 65_536,
                maximum_depth: 32,
            },
            &mut || false,
        )
        .unwrap()
        .evaluate(C64::ONE)
        .unwrap()
    };
    let step = 2e-6;
    let factor = C64::new(0.0, 1.0 / (TAU * hz * media[0].density_kg_m3()));
    for axis in 0..3 {
        let mut left = source;
        let mut right = source;
        left[axis] -= step;
        right[axis] += step;
        let lhs = components(&source_field(left));
        let rhs = components(&source_field(right));
        let exact = components(&force.evaluate(unit_force(axis)).unwrap());
        for component in 0..4 {
            let difference = (rhs[component] - lhs[component]).scale(0.5 / step) * factor;
            close(exact[component], difference, 3e-7);
        }
    }
}

#[test]
fn every_force_velocity_component_satisfies_original_fluid_momentum() {
    let media = [fluid(1.2, 300.0), fluid(2.4, 450.0)];
    let source = [-0.01, 0.005, -0.025];
    let receiver = [0.035, -0.007, 0.018];
    let hz = 700.0;
    let force = [C64::new(0.3, -0.1), C64::new(-0.2, 0.4), C64::new(0.1, 0.7)];
    let field = prepare(media, source, receiver, hz)
        .evaluate(force)
        .unwrap();
    let step = 2e-6;
    for axis in 0..3 {
        let mut left = receiver;
        let mut right = receiver;
        left[axis] -= step;
        right[axis] += step;
        let lhs = prepare(media, source, left, hz)
            .evaluate(force)
            .unwrap()
            .pressure_pa;
        let rhs = prepare(media, source, right, hz)
            .evaluate(force)
            .unwrap()
            .pressure_pa;
        let gradient = (rhs - lhs).scale(0.5 / step);
        close(
            field.particle_velocity_m_s[axis],
            gradient / C64::new(0.0, TAU * hz * media[1].density_kg_m3()),
            3e-7,
        );
    }
}

#[test]
fn axis_retains_two_transverse_velocity_channels_and_rotation_is_covariant() {
    let media = [fluid(1.2, 300.0), fluid(2.4, 450.0)];
    let source = [0.0, 0.0, -0.025];
    let receiver = [0.0, 0.0, 0.018];
    let hz = 700.0;
    let axis = prepare(media, source, receiver, hz);
    let fx = axis.evaluate(unit_force(0)).unwrap();
    let fy = axis.evaluate(unit_force(1)).unwrap();
    assert_eq!(fx.pressure_pa, C64::ZERO);
    assert_eq!(fy.pressure_pa, C64::ZERO);
    assert_eq!(fx.particle_velocity_m_s[0], fy.particle_velocity_m_s[1]);
    assert!(fx.particle_velocity_m_s[0].abs() > 0.001);
    assert_eq!(fx.particle_velocity_m_s[1], C64::ZERO);
    assert_eq!(fx.particle_velocity_m_s[2], C64::ZERO);
    assert_eq!(axis.diagnostics().pressure_tail_bound_pa_per_n[0], 0.0);
    assert_eq!(axis.diagnostics().velocity_tail_bound_m_s_per_n[2], 0.0);
    assert_eq!(axis.diagnostics().velocity_tail_bound_m_s_per_n[3], 0.0);
    assert!(axis.diagnostics().velocity_tail_bound_m_s_per_n[1] > 0.0);

    let receiver = [0.035, -0.007, 0.018];
    let force = [C64::new(0.3, -0.1), C64::new(-0.2, 0.4), C64::new(0.1, 0.7)];
    let reference = prepare(media, source, receiver, hz)
        .evaluate(force)
        .unwrap();
    let transform = |p: [f64; 3]| [1.0 + p[2], -2.0 + p[0], 3.0 + p[1]];
    let rotated_plane =
        PlanarFluidInterface::new(media[0], media[1], transform([0.0; 3]), [7.0, 0.0, 0.0])
            .unwrap();
    let rotated = PreparedFullWaveForce::new(
        rotated_plane,
        hz,
        transform(source),
        transform(receiver),
        limits(),
        &mut || false,
    )
    .unwrap()
    .evaluate([force[2], force[0], force[1]])
    .unwrap();
    close(rotated.pressure_pa, reference.pressure_pa, 3e-10);
    for (new, old) in [(0, 2), (1, 0), (2, 1)] {
        close(
            rotated.particle_velocity_m_s[new],
            reference.particle_velocity_m_s[old],
            3e-10,
        );
    }
}

#[test]
fn tensor_norm_refinement_reduces_tail_and_controls_arbitrary_force() {
    let media = [fluid(1.2, 300.0), fluid(2.4, 450.0)];
    let source = [0.0, 0.0, -0.025];
    let receiver = [0.035, 0.0, 0.018];
    let coarse_limits = FullWaveForceLimits {
        absolute_pressure_tolerance_pa_per_n: 1e-4,
        absolute_velocity_tolerance_m_s_per_n: 1e-6,
        relative_tolerance: 0.0,
        ..limits()
    };
    let coarse = PreparedFullWaveForce::new(
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
    assert!(f.pressure_tail_bound_pa_per_n[1] < c.pressure_tail_bound_pa_per_n[1] * 1e-3);
    let p_error = core::array::from_fn::<_, 2, _>(|i| {
        c.pressure_quadrature_error_pa_per_n[i] + c.pressure_tail_bound_pa_per_n[i]
    });
    let v_error = core::array::from_fn::<_, 5, _>(|i| {
        c.velocity_quadrature_error_m_s_per_n[i] + c.velocity_tail_bound_m_s_per_n[i]
    });
    assert!(
        p_error.iter().map(|v| v * v).sum::<f64>().sqrt()
            <= coarse_limits.absolute_pressure_tolerance_pa_per_n
    );
    assert!(
        v_error.iter().map(|v| v * v).sum::<f64>().sqrt()
            <= coarse_limits.absolute_velocity_tolerance_m_s_per_n
    );
    let force = [C64::new(0.3, -0.1), C64::new(-0.2, 0.4), C64::new(0.1, 0.7)];
    let force_norm = force.iter().map(|v| v.norm_sq()).sum::<f64>().sqrt();
    let fc = coarse.evaluate(force).unwrap();
    let ff = fine.evaluate(force).unwrap();
    assert!(
        (fc.pressure_pa - ff.pressure_pa).abs()
            < force_norm * coarse_limits.absolute_pressure_tolerance_pa_per_n
    );
    let velocity_error = (0..3)
        .map(|i| (fc.particle_velocity_m_s[i] - ff.particle_velocity_m_s[i]).norm_sq())
        .sum::<f64>()
        .sqrt();
    assert!(velocity_error < force_norm * coarse_limits.absolute_velocity_tolerance_m_s_per_n);
}

#[test]
fn harmonic_observer_retains_physical_force_clock_and_active_intensity() {
    let observer = prepare(
        [fluid(1.2, 300.0), fluid(2.4, 450.0)],
        [0.0, 0.0, -0.025],
        [0.035, 0.0, 0.018],
        700.0,
    );
    let force = [C64::new(0.3, -0.1), C64::new(-0.2, 0.4), C64::new(0.1, 0.7)];
    let line = observer.evaluate_harmonic_line(force).unwrap();
    let field = observer.evaluate(force).unwrap();
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
    for axis in 0..3 {
        assert!(
            (observed.mean_intensity_w_m2[axis] - field.active_intensity_w_m2[axis]).abs()
                < 3e-12 * field.active_intensity_w_m2[axis].abs()
        );
    }
}

#[test]
fn invalid_inputs_work_exhaustion_and_cancellation_return_no_prepared_force() {
    let interface = plane([fluid(1.2, 300.0), fluid(2.4, 450.0)]);
    let source = [0.0, 0.0, -0.025];
    let receiver = [0.035, 0.0, 0.018];
    let attempt =
        |hz, s, r, budget| PreparedFullWaveForce::new(interface, hz, s, r, budget, &mut || false);
    for bad in [
        FullWaveForceLimits {
            absolute_pressure_tolerance_pa_per_n: 0.0,
            ..limits()
        },
        FullWaveForceLimits {
            absolute_velocity_tolerance_m_s_per_n: f64::NAN,
            ..limits()
        },
        FullWaveForceLimits {
            relative_tolerance: 1.0,
            ..limits()
        },
        FullWaveForceLimits {
            maximum_depth: 49,
            ..limits()
        },
        FullWaveForceLimits {
            maximum_evaluations: 0,
            ..limits()
        },
        FullWaveForceLimits {
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
            FullWaveForceLimits {
                maximum_depth: 1,
                ..limits()
            }
        ),
        Err(AeroacError::NotConverged { .. })
    ));
    let mut polls = 0;
    assert!(matches!(
        PreparedFullWaveForce::new(interface, 700.0, source, receiver, limits(), &mut || {
            polls += 1;
            polls == 100
        }),
        Err(AeroacError::Cancelled { .. })
    ));
    let observer = attempt(700.0, source, receiver, limits()).unwrap();
    assert!(
        observer
            .evaluate([C64::new(f64::INFINITY, 0.0); 3])
            .is_err()
    );
    let silence = observer.evaluate([C64::ZERO; 3]).unwrap();
    assert_eq!(silence.pressure_pa, C64::ZERO);
    assert_eq!(silence.particle_velocity_m_s, [C64::ZERO; 3]);
    let replay = attempt(700.0, source, receiver, limits()).unwrap();
    assert_eq!(
        components(&observer.evaluate(unit_force(2)).unwrap()),
        components(&replay.evaluate(unit_force(2)).unwrap())
    );
}
