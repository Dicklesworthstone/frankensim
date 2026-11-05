//! G1/G3: physical force fields on the original incident side of a fluid plane.

use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::refraction::PlanarFluidInterface;
use fs_aeroac::curle3d::refraction::full_wave::reflection::{
    PreparedReflectedFullWaveForce, PreparedReflectedFullWaveMonopole,
};
use fs_aeroac::curle3d::refraction::full_wave::{FullWaveForceLimits, FullWaveMonopoleLimits};
use fs_aeroac::curle3d::{CompactSourcePhasor, PreparedCompactSourceObserver};
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::c64::C64;
use std::f64::consts::TAU;

fn fluid(rho: f64, c: f64) -> FluidAcousticState {
    FluidAcousticState::try_new(rho, c).unwrap()
}
fn interface(states: [FluidAcousticState; 2]) -> PlanarFluidInterface {
    PlanarFluidInterface::new(states[0], states[1], [0.0; 3], [0.0, 0.0, 1.0]).unwrap()
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
fn all_seven_reflected_force_channels_match_independent_direct_spectrum_references() {
    // Independent mpmath 1.3.0 direct-q tanh-sinh, split at both critical q.
    // Refinement from 35 to 50 decimal digits; no production substitutions,
    // logarithmic coefficients or GK nodes. Both speed orderings and a small
    // original speed contrast retain their complete evanescent spectra.
    let cases = [
        (
            [fluid(1.2, 300.0), fluid(2.4, 600.0)],
            700.0,
            0.025,
            0.023,
            0.04,
            [
                C64::new(5.49280389483676, 0.18714327914220203),
                C64::new(7.501064298682157, -1.5210671685188153),
                C64::new(0.0009941298080441253, 0.0019337064234342956),
                C64::new(0.0008864517244615431, -0.026018062241002547),
                C64::new(0.0016484727167333123, 0.031643179976080235),
                C64::new(-0.0016484727167333123, -0.031643179976080235),
                C64::new(-0.01079487460733584, -0.007785969744400459),
            ],
        ),
        (
            [fluid(2.4, 600.0), fluid(1.2, 300.0)],
            1500.0,
            0.04,
            0.06,
            0.12,
            [
                C64::new(-2.7693821678200505, -3.2719312014191795),
                C64::new(-1.9180828231501153, -3.9989741181713915),
                C64::new(-0.001613709698032413, -0.0014310251505219166),
                C64::new(-0.0012054259818575608, 0.0010202797715720136),
                C64::new(0.00008130281066098722, -0.0028789613590562453),
                C64::new(-0.00008130281066098722, 0.0028789613590562453),
                C64::new(0.00010211607988956937, 0.0019143620230789085),
            ],
        ),
        (
            [fluid(1.2, 343.0), fluid(2.4, 343.01)],
            1000.0,
            0.03,
            0.05,
            0.07,
            [
                C64::new(2.2310390405505447, 2.544095296107232),
                C64::new(2.549850443244649, 2.9075820070279734),
                C64::new(0.003028086404592746, 0.0030420959228500817),
                C64::new(0.004820301691336966, -0.0042271534706502225),
                C64::new(-0.002048122262367726, 0.008307776539815938),
                C64::new(0.002048122262367726, -0.008307776539815938),
                C64::new(-0.002479773724205308, -0.005267725481965003),
            ],
        ),
    ];
    for (states, hz, a, b, r, expected) in cases {
        let prepared = PreparedReflectedFullWaveForce::new(
            interface(states),
            hz,
            [0.0, 0.0, -a],
            [r, 0.0, -b],
            limits(),
            &mut || false,
        )
        .unwrap();
        let x = prepared.evaluate_reflected(unit(0)).unwrap();
        let y = prepared.evaluate_reflected(unit(1)).unwrap();
        let z = prepared.evaluate_reflected(unit(2)).unwrap();
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
        close(actual[4], -actual[5], 2e-12);
        assert!(prepared.diagnostics().evaluations > 0);
        assert!(prepared.diagnostics().evaluations <= limits().maximum_evaluations);
        assert!(!prepared.diagnostics().exact_equal_speed);
    }
}

#[test]
fn original_equal_speed_image_force_reverses_only_its_normal_component() {
    let incident = fluid(1.2, 343.0);
    let receiving = fluid(2.4, 343.0);
    let normal = [1.0 / 3.0, 2.0 / 3.0, -2.0 / 3.0];
    let tangent = [2.0 / 5.0_f64.sqrt(), -1.0 / 5.0_f64.sqrt(), 0.0];
    let origin = [0.2, -0.1, 0.4];
    let source = core::array::from_fn(|i| origin[i] - 0.08 * normal[i] + 0.03 * tangent[i]);
    let observer = core::array::from_fn(|i| origin[i] - 0.13 * normal[i] + 0.2 * tangent[i]);
    let image = core::array::from_fn(|i| source[i] + 0.16 * normal[i]);
    let force = [C64::new(1.0, 0.3), C64::new(-0.2, 0.8), C64::new(0.7, -0.1)];
    let normal_force = force
        .into_iter()
        .zip(normal)
        .fold(C64::ZERO, |s, (f, n)| s + f.scale(n));
    let image_force =
        core::array::from_fn(|i| (force[i] - normal_force.scale(2.0 * normal[i])).scale(1.0 / 3.0));
    let plane = PlanarFluidInterface::new(incident, receiving, origin, normal).unwrap();
    let prepared = PreparedReflectedFullWaveForce::new(
        plane,
        913.0,
        source,
        observer,
        FullWaveForceLimits {
            maximum_evaluations: 1,
            ..limits()
        },
        &mut || false,
    )
    .unwrap();
    let original =
        PreparedCompactSourceObserver::new(incident, 913.0, observer, &[source, image], 2).unwrap();
    let expected_reflected = original
        .evaluate(&[
            CompactSourcePhasor::default(),
            CompactSourcePhasor::force_dipole(image_force),
        ])
        .unwrap();
    let expected_total = original
        .evaluate(&[
            CompactSourcePhasor::force_dipole(force),
            CompactSourcePhasor::force_dipole(image_force),
        ])
        .unwrap();
    for (actual, expected) in [
        (
            prepared.evaluate_reflected(force).unwrap(),
            expected_reflected,
        ),
        (prepared.evaluate_total(force).unwrap(), expected_total),
    ] {
        close(actual.pressure_pa, expected.pressure_pa, 2e-12);
        for i in 0..3 {
            close(
                actual.particle_velocity_m_s[i],
                expected.particle_velocity_m_s[i],
                2e-12,
            );
            close(
                C64::new(
                    actual.active_intensity_w_m2[i],
                    actual.reactive_intensity_w_m2[i],
                ),
                expected.pressure_pa.scale(0.5) * expected.particle_velocity_m_s[i].conj(),
                4e-12,
            );
        }
    }
    assert!(prepared.diagnostics().exact_equal_speed);
    assert_eq!(prepared.diagnostics().evaluations, 0);
    let matched = PlanarFluidInterface::new(incident, incident, origin, normal).unwrap();
    let no_boundary = PreparedReflectedFullWaveForce::new(
        matched,
        913.0,
        source,
        observer,
        limits(),
        &mut || false,
    )
    .unwrap();
    let reflected = no_boundary.evaluate_reflected(force).unwrap();
    assert_eq!(reflected.pressure_pa, C64::ZERO);
    assert_eq!(reflected.particle_velocity_m_s, [C64::ZERO; 3]);
    let direct = original
        .evaluate(&[
            CompactSourcePhasor::force_dipole(force),
            CompactSourcePhasor::default(),
        ])
        .unwrap();
    close(
        no_boundary.evaluate_total(force).unwrap().pressure_pa,
        direct.pressure_pa,
        2e-13,
    );
}

#[test]
fn reflected_force_is_the_source_derivative_and_obeys_incident_fluid_momentum() {
    let plane = interface([fluid(1.2, 300.0), fluid(2.4, 600.0)]);
    let hz = 700.0;
    let source = [0.0, 0.0, -0.025];
    let observer = [0.04, 0.0, -0.023];
    let prepared =
        PreparedReflectedFullWaveForce::new(plane, hz, source, observer, limits(), &mut || false)
            .unwrap();
    let force = [C64::new(0.3, 0.1), C64::new(-0.2, 0.4), C64::new(0.7, -0.2)];
    let actual = prepared.evaluate_reflected(force).unwrap();
    let delta = 2e-6;
    let mono_limits = FullWaveMonopoleLimits {
        absolute_pressure_tolerance_pa_per_m3_s: 1e-9,
        absolute_velocity_tolerance_per_m2: 1e-12,
        relative_tolerance: 1e-12,
        maximum_evaluations: 262_144,
        maximum_depth: 40,
    };
    for axis in 0..3 {
        let mut plus = source;
        let mut minus = source;
        plus[axis] += delta;
        minus[axis] -= delta;
        let mono = |s| {
            PreparedReflectedFullWaveMonopole::new(plane, hz, s, observer, mono_limits, &mut || {
                false
            })
            .unwrap()
            .evaluate_reflected(C64::ONE)
            .unwrap()
            .pressure_pa
        };
        let derivative =
            (mono(plus) - mono(minus)).scale(0.5 / delta) * C64::new(0.0, 1.0 / (TAU * hz * 1.2));
        close(
            prepared.evaluate_reflected(unit(axis)).unwrap().pressure_pa,
            derivative,
            2e-5,
        );
        let mut plus = observer;
        let mut minus = observer;
        plus[axis] += delta;
        minus[axis] -= delta;
        let pressure = |x| {
            PreparedReflectedFullWaveForce::new(plane, hz, source, x, limits(), &mut || false)
                .unwrap()
                .evaluate_reflected(force)
                .unwrap()
                .pressure_pa
        };
        let velocity = (pressure(plus) - pressure(minus)).scale(0.5 / delta)
            * C64::new(0.0, -1.0 / (TAU * hz * 1.2));
        close(actual.particle_velocity_m_s[axis], velocity, 2e-5);
    }
}

#[test]
fn reflected_self_field_axis_motion_clock_and_bounded_refusals_are_explicit() {
    let plane = interface([fluid(1.2, 300.0), fluid(2.4, 600.0)]);
    let source = [0.0, 0.0, -0.04];
    let observer = [0.0, 0.0, -0.07];
    let hz = 500.0;
    let prepare = |s, o, l| PreparedReflectedFullWaveForce::new(plane, hz, s, o, l, &mut || false);
    let prepared = prepare(source, observer, limits()).unwrap();
    let x = prepared.evaluate_reflected(unit(0)).unwrap();
    let y = prepared.evaluate_reflected(unit(1)).unwrap();
    assert_eq!(x.pressure_pa, C64::ZERO);
    assert!(x.particle_velocity_m_s[0].abs() > 1e-6);
    close(
        x.particle_velocity_m_s[0],
        y.particle_velocity_m_s[1],
        2e-13,
    );
    let force = [C64::new(0.2, -0.4), C64::new(0.1, 0.3), C64::new(-0.5, 0.2)];
    for (field, line) in [
        (
            prepared.evaluate_reflected(force).unwrap(),
            prepared.evaluate_reflected_harmonic_line(force).unwrap(),
        ),
        (
            prepared.evaluate_total(force).unwrap(),
            prepared.evaluate_total_harmonic_line(force).unwrap(),
        ),
    ] {
        assert_eq!(line.frequency_hz, hz);
        assert_eq!(line.pressure_pa, field.pressure_pa);
        assert_eq!(line.particle_velocity_m_s, field.particle_velocity_m_s);
    }
    let at_source = prepare(source, source, limits()).unwrap();
    assert!(!at_source.has_total_field());
    assert!(
        at_source
            .evaluate_reflected(force)
            .unwrap()
            .pressure_pa
            .abs()
            .is_finite()
    );
    assert!(at_source.evaluate_total(force).is_err());
    assert!(at_source.evaluate_total_harmonic_line(force).is_err());
    assert!(prepared.has_total_field());
    let silent = prepared.evaluate_total([C64::ZERO; 3]).unwrap();
    assert_eq!(silent.pressure_pa, C64::ZERO);
    assert_eq!(silent.particle_velocity_m_s, [C64::ZERO; 3]);
    assert!(
        prepared
            .evaluate_reflected([C64::new(f64::NAN, 0.0); 3])
            .is_err()
    );
    assert!(
        prepare(
            source,
            observer,
            FullWaveForceLimits {
                maximum_evaluations: 1,
                ..limits()
            }
        )
        .is_err()
    );
    for point in [
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.04],
        [f64::INFINITY, 0.0, -0.04],
    ] {
        assert!(prepare(point, observer, limits()).is_err());
        assert!(prepare(source, point, limits()).is_err());
    }
    assert!(matches!(
        PreparedReflectedFullWaveForce::new(plane, hz, source, observer, limits(), &mut || true),
        Err(AeroacError::Cancelled { .. })
    ));
    assert!(
        PreparedReflectedFullWaveForce::new(plane, 0.0, source, observer, limits(), &mut || false)
            .is_err()
    );
}
