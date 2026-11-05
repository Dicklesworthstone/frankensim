//! Refraction preserves physical Q/F source coherence and global work bounds.

use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::refraction::PlanarFluidInterface;
use fs_aeroac::curle3d::refraction::full_wave::array::{
    FullWaveCompactSourceLimits, PreparedFullWaveCompactSourceObserver,
};
use fs_aeroac::curle3d::refraction::full_wave::{
    FullWaveForceLimits, FullWaveMonopoleLimits, PreparedFullWaveForce, PreparedFullWaveMonopole,
};
use fs_aeroac::curle3d::{CoherenceBudget, CompactSourcePhasor, PreparedCompactSourceObserver};
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::c64::C64;

fn limits() -> FullWaveCompactSourceLimits {
    FullWaveCompactSourceLimits {
        max_sources: 4,
        maximum_total_evaluations: 524_288,
        monopole: FullWaveMonopoleLimits {
            absolute_pressure_tolerance_pa_per_m3_s: 1e-7,
            absolute_velocity_tolerance_per_m2: 1e-10,
            relative_tolerance: 1e-9,
            maximum_evaluations: 65_536,
            maximum_depth: 32,
        },
        force: FullWaveForceLimits {
            absolute_pressure_tolerance_pa_per_n: 1e-8,
            absolute_velocity_tolerance_m_s_per_n: 1e-11,
            relative_tolerance: 1e-9,
            maximum_evaluations: 65_536,
            maximum_depth: 32,
        },
    }
}

fn plane(c2: f64) -> PlanarFluidInterface {
    PlanarFluidInterface::new(
        FluidAcousticState::try_new(1.2, 300.0).unwrap(),
        FluidAcousticState::try_new(2.4, c2).unwrap(),
        [0.0; 3],
        [0.0, 0.0, 1.0],
    )
    .unwrap()
}

fn sources() -> [CompactSourcePhasor; 2] {
    [
        CompactSourcePhasor {
            volume_flow_m3_s: C64::new(0.0001, -0.0003),
            force_on_fluid_n: [
                C64::new(0.2, 0.1),
                C64::new(-0.1, 0.05),
                C64::new(0.1, -0.2),
            ],
        },
        CompactSourcePhasor {
            volume_flow_m3_s: C64::new(-0.0002, 0.0001),
            force_on_fluid_n: [
                C64::new(-0.04, 0.03),
                C64::new(0.06, -0.02),
                C64::new(0.08, 0.09),
            ],
        },
    ]
}

fn close(actual: C64, expected: C64) {
    assert!(
        (actual - expected).abs() <= 2e-12 * expected.abs().max(1e-12),
        "actual={actual:?}, expected={expected:?}"
    );
}

#[test]
fn exact_equal_speed_array_preserves_q_force_cross_terms_and_receiving_density() {
    let interface = plane(300.0);
    let positions = [[0.0, 0.0, -0.025], [-0.01, 0.005, -0.04]];
    let receiver = [0.035, 0.0, 0.018];
    let mut budget = limits();
    budget.maximum_total_evaluations = 1;
    budget.monopole.maximum_evaluations = 1;
    budget.force.maximum_evaluations = 1;
    let refracted = PreparedFullWaveCompactSourceObserver::new(
        interface,
        700.0,
        receiver,
        &positions,
        budget,
        &mut || false,
    )
    .unwrap();
    let homogeneous = PreparedCompactSourceObserver::new(
        interface.incident_state(),
        700.0,
        receiver,
        &positions,
        2,
    )
    .unwrap();
    let actual = refracted.evaluate(&sources()).unwrap();
    let expected = homogeneous.evaluate(&sources()).unwrap();
    close(actual.pressure_pa, expected.pressure_pa.scale(4.0 / 3.0));
    for axis in 0..3 {
        close(
            actual.particle_velocity_m_s[axis],
            expected.particle_velocity_m_s[axis].scale(2.0 / 3.0),
        );
    }
    assert_eq!(refracted.evaluations(), 0);
    assert_eq!(refracted.source_count(), 2);
    assert_eq!(refracted.source_position_m(1), Some(positions[1]));
    assert_eq!(refracted.source_position_m(2), None);
    for index in 0..2 {
        let diagnostics = refracted.source_diagnostics(index).unwrap();
        assert!(diagnostics.monopole.exact_equal_speed && diagnostics.force.exact_equal_speed);
    }
    let line = refracted.evaluate_harmonic_line(&sources()).unwrap();
    assert_eq!(line.frequency_hz, 700.0);
    assert_eq!(line.pressure_pa, actual.pressure_pa);
    assert_eq!(line.particle_velocity_m_s, actual.particle_velocity_m_s);
    let one_mode = refracted
        .evaluate_coherence_modes(
            &[&sources()],
            CoherenceBudget {
                max_modes: 1,
                max_source_mode_pairs: 2,
            },
        )
        .unwrap();
    close(
        C64::from_re(one_mode.mean_square_pressure_pa2),
        C64::from_re(0.5 * actual.pressure_pa.abs().powi(2)),
    );
    for axis in 0..3 {
        close(
            C64::new(
                one_mode.active_intensity_w_m2[axis],
                one_mode.reactive_intensity_w_m2[axis],
            ),
            actual.pressure_pa.scale(0.5) * actual.particle_velocity_m_s[axis].conj(),
        );
    }
}

#[test]
fn refracted_array_matches_original_kernels_and_keeps_explicit_source_statistics() {
    let interface = plane(450.0);
    let positions = [[0.0, 0.0, -0.025], [-0.01, 0.005, -0.04]];
    let receiver = [0.035, 0.0, 0.018];
    let source = sources();
    let observer = PreparedFullWaveCompactSourceObserver::new(
        interface,
        700.0,
        receiver,
        &positions,
        limits(),
        &mut || false,
    )
    .unwrap();
    let mut expected = [C64::ZERO; 4];
    let mut independent_pressure_square = 0.0;
    let mut expected_work = 0;
    for (position, amplitude) in positions.into_iter().zip(source) {
        let q = PreparedFullWaveMonopole::new(
            interface,
            700.0,
            position,
            receiver,
            limits().monopole,
            &mut || false,
        )
        .unwrap();
        let f = PreparedFullWaveForce::new(
            interface,
            700.0,
            position,
            receiver,
            limits().force,
            &mut || false,
        )
        .unwrap();
        expected_work += q.diagnostics().evaluations + f.diagnostics().evaluations;
        let q = q
            .evaluate_harmonic_line(amplitude.volume_flow_m3_s)
            .unwrap();
        let f = f
            .evaluate_harmonic_line(amplitude.force_on_fluid_n)
            .unwrap();
        let p = q.pressure_pa + f.pressure_pa;
        independent_pressure_square += 0.5 * p.abs().powi(2);
        expected[0] = expected[0] + p;
        for axis in 0..3 {
            expected[axis + 1] =
                expected[axis + 1] + q.particle_velocity_m_s[axis] + f.particle_velocity_m_s[axis];
        }
    }
    let actual = observer.evaluate(&source).unwrap();
    close(actual.pressure_pa, expected[0]);
    for axis in 0..3 {
        close(actual.particle_velocity_m_s[axis], expected[axis + 1]);
    }
    assert_eq!(observer.evaluations(), expected_work);
    let independent = observer.evaluate_independent_sources(&source).unwrap();
    close(
        C64::from_re(independent.mean_square_pressure_pa2),
        C64::from_re(independent_pressure_square),
    );
    let first = [source[0], CompactSourcePhasor::default()];
    let second = [CompactSourcePhasor::default(), source[1]];
    let factored = observer
        .evaluate_coherence_modes(
            &[&first, &second],
            CoherenceBudget {
                max_modes: 2,
                max_source_mode_pairs: 4,
            },
        )
        .unwrap();
    close(
        C64::from_re(independent.mean_square_pressure_pa2),
        C64::from_re(factored.mean_square_pressure_pa2),
    );
    assert_eq!(
        independent.velocity_covariance_m2_s2,
        factored.velocity_covariance_m2_s2
    );
    assert_eq!(
        independent.active_intensity_w_m2,
        factored.active_intensity_w_m2
    );
    assert!(
        (independent.mean_square_pressure_pa2 - 0.5 * actual.pressure_pa.abs().powi(2)).abs()
            > 1e-3
    );
    // A tighter shared cap must fail even though both per-kernel caps still fit.
    let tight = FullWaveCompactSourceLimits {
        maximum_total_evaluations: expected_work - 1,
        ..limits()
    };
    assert!(
        PreparedFullWaveCompactSourceObserver::new(
            interface,
            700.0,
            receiver,
            &positions,
            tight,
            &mut || false,
        )
        .is_err()
    );
}

#[test]
fn coherent_nulls_independent_energy_empty_arrays_and_refusals_remain_distinct() {
    let interface = plane(300.0);
    let receiver = [0.035, 0.0, 0.018];
    let position = [0.0, 0.0, -0.025];
    let observer = PreparedFullWaveCompactSourceObserver::new(
        interface,
        700.0,
        receiver,
        &[position, position],
        limits(),
        &mut || false,
    )
    .unwrap();
    let first = sources()[0];
    let opposite = CompactSourcePhasor {
        volume_flow_m3_s: first.volume_flow_m3_s.scale(-1.0),
        force_on_fluid_n: first.force_on_fluid_n.map(|v| v.scale(-1.0)),
    };
    let source = [first, opposite];
    let field = observer.evaluate(&source).unwrap();
    assert_eq!(field.pressure_pa, C64::ZERO);
    assert_eq!(field.particle_velocity_m_s, [C64::ZERO; 3]);
    assert!(
        observer
            .evaluate_independent_sources(&source)
            .unwrap()
            .mean_square_pressure_pa2
            > 0.0
    );
    assert!(observer.evaluate(&source[..1]).is_err());
    let bad = CompactSourcePhasor::monopole(C64::new(f64::NAN, 0.0));
    assert!(observer.evaluate(&[bad, first]).is_err());
    assert!(
        observer
            .evaluate_coherence_modes(
                &[&source],
                CoherenceBudget {
                    max_modes: 1,
                    max_source_mode_pairs: 1
                },
            )
            .is_err()
    );
    assert!(matches!(
        PreparedFullWaveCompactSourceObserver::new(
            interface,
            700.0,
            receiver,
            &[position],
            limits(),
            &mut || true,
        ),
        Err(AeroacError::Cancelled { .. })
    ));
    assert!(
        PreparedFullWaveCompactSourceObserver::new(
            interface,
            700.0,
            receiver,
            &[position],
            FullWaveCompactSourceLimits {
                max_sources: 0,
                ..limits()
            },
            &mut || false,
        )
        .is_err()
    );
    let empty = PreparedFullWaveCompactSourceObserver::new(
        interface,
        700.0,
        receiver,
        &[],
        FullWaveCompactSourceLimits {
            max_sources: 0,
            ..limits()
        },
        &mut || false,
    )
    .unwrap();
    assert_eq!(empty.evaluate(&[]).unwrap().pressure_pa, C64::ZERO);
    assert_eq!(empty.evaluations(), 0);
    assert!(
        PreparedFullWaveCompactSourceObserver::new(
            interface,
            f64::NAN,
            receiver,
            &[],
            limits(),
            &mut || false,
        )
        .is_err()
    );
    assert!(
        PreparedFullWaveCompactSourceObserver::new(
            interface,
            700.0,
            [0.0; 3],
            &[],
            limits(),
            &mut || false,
        )
        .is_err()
    );
}
