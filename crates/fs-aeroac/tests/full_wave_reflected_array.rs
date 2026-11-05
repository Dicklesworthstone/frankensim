//! G1/G3: original source coherence survives incident-side full-wave observation.

use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::refraction::PlanarFluidInterface;
use fs_aeroac::curle3d::refraction::full_wave::array::{
    FULL_WAVE_REFLECTED_ARRAY_SCOPE, FULL_WAVE_TOTAL_ARRAY_SCOPE, FullWaveCompactFieldKind,
    FullWaveCompactSourceLimits, PreparedFullWaveCompactSourceObserver,
};
use fs_aeroac::curle3d::refraction::full_wave::reflection::{
    PreparedReflectedFullWaveForce, PreparedReflectedFullWaveMonopole,
};
use fs_aeroac::curle3d::refraction::full_wave::{FullWaveForceLimits, FullWaveMonopoleLimits};
use fs_aeroac::curle3d::{
    CoherenceBudget, CoherentCompactField, CompactFieldStatistics, CompactSourcePhasor,
    PreparedCompactSourceObserver,
};
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::c64::C64;

fn plane(c2: f64) -> PlanarFluidInterface {
    PlanarFluidInterface::new(
        FluidAcousticState::try_new(1.2, 300.0).unwrap(),
        FluidAcousticState::try_new(2.4, c2).unwrap(),
        [0.0; 3],
        [0.0, 0.0, 1.0],
    )
    .unwrap()
}

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

fn image_source(source: CompactSourcePhasor) -> CompactSourcePhasor {
    CompactSourcePhasor {
        volume_flow_m3_s: source.volume_flow_m3_s.scale(1.0 / 3.0),
        force_on_fluid_n: core::array::from_fn(|i| {
            source.force_on_fluid_n[i].scale(if i == 2 { -1.0 / 3.0 } else { 1.0 / 3.0 })
        }),
    }
}

fn close(actual: C64, expected: C64) {
    assert!(
        (actual - expected).abs() <= 3e-12 * expected.abs().max(1e-12),
        "actual={actual:?}, expected={expected:?}"
    );
}

fn field_close(actual: &CoherentCompactField, expected: &CoherentCompactField) {
    close(actual.pressure_pa, expected.pressure_pa);
    for i in 0..3 {
        close(
            actual.particle_velocity_m_s[i],
            expected.particle_velocity_m_s[i],
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
        );
    }
}

fn statistics_close(actual: &CompactFieldStatistics, expected: &CompactFieldStatistics) {
    close(
        C64::from_re(actual.mean_square_pressure_pa2),
        C64::from_re(expected.mean_square_pressure_pa2),
    );
    for i in 0..3 {
        close(
            C64::new(
                actual.active_intensity_w_m2[i],
                actual.reactive_intensity_w_m2[i],
            ),
            C64::new(
                expected.active_intensity_w_m2[i],
                expected.reactive_intensity_w_m2[i],
            ),
        );
        for j in 0..3 {
            close(
                actual.velocity_covariance_m2_s2[i][j],
                expected.velocity_covariance_m2_s2[i][j],
            );
        }
    }
}

#[test]
fn density_jump_images_retain_q_force_and_direct_reflected_statistical_cross_terms() {
    let interface = plane(300.0);
    let positions = [[0.0, 0.0, -0.025], [-0.01, 0.005, -0.04]];
    let receiver = [0.035, -0.012, -0.018];
    let amplitudes = sources();
    let image_positions = positions.map(|p| [p[0], p[1], -p[2]]);
    let image_amplitudes = amplitudes.map(image_source);
    let mut budget = limits();
    budget.maximum_total_evaluations = 1;
    budget.monopole.maximum_evaluations = 1;
    budget.force.maximum_evaluations = 1;
    let reflected = PreparedFullWaveCompactSourceObserver::new_reflected(
        interface,
        700.0,
        receiver,
        &positions,
        budget,
        &mut || false,
    )
    .unwrap();
    let total = PreparedFullWaveCompactSourceObserver::new_total(
        interface,
        700.0,
        receiver,
        &positions,
        budget,
        &mut || false,
    )
    .unwrap();
    let image = PreparedCompactSourceObserver::new(
        interface.incident_state(),
        700.0,
        receiver,
        &image_positions,
        2,
    )
    .unwrap();
    let combined = PreparedCompactSourceObserver::new(
        interface.incident_state(),
        700.0,
        receiver,
        &[
            positions[0],
            positions[1],
            image_positions[0],
            image_positions[1],
        ],
        4,
    )
    .unwrap();
    let actual = total.evaluate(&amplitudes).unwrap();
    field_close(
        &reflected.evaluate(&amplitudes).unwrap(),
        &image.evaluate(&image_amplitudes).unwrap(),
    );
    field_close(
        &actual,
        &combined
            .evaluate(&[
                amplitudes[0],
                amplitudes[1],
                image_amplitudes[0],
                image_amplitudes[1],
            ])
            .unwrap(),
    );
    assert_eq!(reflected.field_kind(), FullWaveCompactFieldKind::Reflected);
    assert_eq!(total.field_kind(), FullWaveCompactFieldKind::IncidentTotal);
    assert_eq!(reflected.scope(), FULL_WAVE_REFLECTED_ARRAY_SCOPE);
    assert_eq!(actual.scope, FULL_WAVE_TOTAL_ARRAY_SCOPE);
    assert_eq!(total.evaluations(), 0);
    assert_eq!(reflected.evaluations(), 0);
    assert_eq!(total.source_count(), 2);
    assert_eq!(total.source_position_m(1), Some(positions[1]));
    for i in 0..2 {
        let d = total.source_diagnostics(i).unwrap();
        assert!(d.monopole.exact_equal_speed && d.force.exact_equal_speed);
    }
    let zero = CompactSourcePhasor::default();
    let expected = combined
        .evaluate_coherence_modes(
            &[
                &[amplitudes[0], zero, image_amplitudes[0], zero],
                &[zero, amplitudes[1], zero, image_amplitudes[1]],
            ],
            CoherenceBudget {
                max_modes: 2,
                max_source_mode_pairs: 8,
            },
        )
        .unwrap();
    let independent = total.evaluate_independent_sources(&amplitudes).unwrap();
    statistics_close(&independent, &expected);
    assert_eq!(independent.scope, FULL_WAVE_TOTAL_ARRAY_SCOPE);
    let factored = total
        .evaluate_coherence_modes(
            &[&[amplitudes[0], zero], &[zero, amplitudes[1]]],
            CoherenceBudget {
                max_modes: 2,
                max_source_mode_pairs: 4,
            },
        )
        .unwrap();
    statistics_close(&factored, &expected);
    let original = PreparedCompactSourceObserver::new(
        interface.incident_state(),
        700.0,
        receiver,
        &positions,
        2,
    )
    .unwrap();
    let incorrectly_independent_paths = original
        .evaluate_independent_sources(&amplitudes)
        .unwrap()
        .mean_square_pressure_pa2
        + reflected
            .evaluate_independent_sources(&amplitudes)
            .unwrap()
            .mean_square_pressure_pa2;
    assert!(
        (independent.mean_square_pressure_pa2 - incorrectly_independent_paths).abs()
            > 1e-5 * independent.mean_square_pressure_pa2
    );
    let line = total.evaluate_harmonic_line(&amplitudes).unwrap();
    assert_eq!(line.frequency_hz, 700.0);
    assert_eq!(line.pressure_pa, actual.pressure_pa);
    assert_eq!(line.particle_velocity_m_s, actual.particle_velocity_m_s);
}

#[test]
fn unequal_speed_arrays_reuse_original_spherical_kernels_and_one_global_work_budget() {
    let interface = plane(600.0);
    let positions = [[0.0, 0.0, -0.025], [-0.01, 0.005, -0.04]];
    let receiver = [0.04, -0.012, -0.023];
    let amplitudes = sources();
    let reflected = PreparedFullWaveCompactSourceObserver::new_reflected(
        interface,
        700.0,
        receiver,
        &positions,
        limits(),
        &mut || false,
    )
    .unwrap();
    let total = PreparedFullWaveCompactSourceObserver::new_total(
        interface,
        700.0,
        receiver,
        &positions,
        limits(),
        &mut || false,
    )
    .unwrap();
    let mut expected_reflected = [C64::ZERO; 4];
    let mut expected_total = [C64::ZERO; 4];
    let mut work = 0;
    for (position, amplitude) in positions.into_iter().zip(amplitudes) {
        let q = PreparedReflectedFullWaveMonopole::new(
            interface,
            700.0,
            position,
            receiver,
            limits().monopole,
            &mut || false,
        )
        .unwrap();
        let f = PreparedReflectedFullWaveForce::new(
            interface,
            700.0,
            position,
            receiver,
            limits().force,
            &mut || false,
        )
        .unwrap();
        work += q.diagnostics().evaluations + f.diagnostics().evaluations;
        for (expected, q, f) in [
            (
                &mut expected_reflected,
                q.evaluate_reflected_harmonic_line(amplitude.volume_flow_m3_s)
                    .unwrap(),
                f.evaluate_reflected_harmonic_line(amplitude.force_on_fluid_n)
                    .unwrap(),
            ),
            (
                &mut expected_total,
                q.evaluate_total_harmonic_line(amplitude.volume_flow_m3_s)
                    .unwrap(),
                f.evaluate_total_harmonic_line(amplitude.force_on_fluid_n)
                    .unwrap(),
            ),
        ] {
            expected[0] = expected[0] + q.pressure_pa + f.pressure_pa;
            for i in 0..3 {
                expected[i + 1] =
                    expected[i + 1] + q.particle_velocity_m_s[i] + f.particle_velocity_m_s[i];
            }
        }
    }
    for (observer, expected) in [(&reflected, expected_reflected), (&total, expected_total)] {
        let actual = observer.evaluate(&amplitudes).unwrap();
        close(actual.pressure_pa, expected[0]);
        for i in 0..3 {
            close(actual.particle_velocity_m_s[i], expected[i + 1]);
        }
        assert_eq!(observer.evaluations(), work);
        assert!(work > 0 && work <= limits().maximum_total_evaluations);
    }
    assert!(
        PreparedFullWaveCompactSourceObserver::new_total(
            interface,
            700.0,
            receiver,
            &positions,
            FullWaveCompactSourceLimits {
                maximum_total_evaluations: work - 1,
                ..limits()
            },
            &mut || false,
        )
        .is_err()
    );
}

#[test]
fn finite_reflected_self_fields_coherent_nulls_and_singular_totals_remain_distinct() {
    let interface = plane(300.0);
    let position = [0.01, -0.02, -0.03];
    let first = sources()[0];
    let opposite = CompactSourcePhasor {
        volume_flow_m3_s: -first.volume_flow_m3_s,
        force_on_fluid_n: first.force_on_fluid_n.map(|v| -v),
    };
    let reflected = PreparedFullWaveCompactSourceObserver::new_reflected(
        interface,
        700.0,
        position,
        &[position, position],
        limits(),
        &mut || false,
    )
    .unwrap();
    let cancelled = reflected.evaluate(&[first, opposite]).unwrap();
    assert_eq!(cancelled.pressure_pa, C64::ZERO);
    assert_eq!(cancelled.particle_velocity_m_s, [C64::ZERO; 3]);
    assert!(
        reflected
            .evaluate_independent_sources(&[first, opposite])
            .unwrap()
            .mean_square_pressure_pa2
            > 0.0
    );
    assert!(
        PreparedFullWaveCompactSourceObserver::new_total(
            interface,
            700.0,
            position,
            &[position],
            limits(),
            &mut || false,
        )
        .is_err()
    );
    assert!(reflected.evaluate(&[first]).is_err());
    assert!(
        reflected
            .evaluate(&[
                first,
                CompactSourcePhasor::monopole(C64::new(f64::NAN, 0.0)),
            ])
            .is_err()
    );
    assert!(
        reflected
            .evaluate_coherence_modes(
                &[&[first, opposite]],
                CoherenceBudget {
                    max_modes: 1,
                    max_source_mode_pairs: 1
                },
            )
            .is_err()
    );
    assert!(matches!(
        PreparedFullWaveCompactSourceObserver::new_reflected(
            interface,
            700.0,
            position,
            &[position],
            limits(),
            &mut || true,
        ),
        Err(AeroacError::Cancelled { .. })
    ));
    assert!(
        PreparedFullWaveCompactSourceObserver::new_reflected(
            interface,
            700.0,
            position,
            &[position],
            FullWaveCompactSourceLimits {
                max_sources: 0,
                ..limits()
            },
            &mut || false,
        )
        .is_err()
    );
    for receiver in [[0.0; 3], [0.0, 0.0, 0.01], [f64::NAN, 0.0, -0.01]] {
        assert!(
            PreparedFullWaveCompactSourceObserver::new_reflected(
                interface,
                700.0,
                receiver,
                &[],
                limits(),
                &mut || false,
            )
            .is_err()
        );
        assert!(
            PreparedFullWaveCompactSourceObserver::new_total(
                interface,
                700.0,
                receiver,
                &[],
                limits(),
                &mut || false,
            )
            .is_err()
        );
    }
    for source in [[0.0; 3], [0.0, 0.0, 0.01], [f64::NAN, 0.0, -0.01]] {
        assert!(
            PreparedFullWaveCompactSourceObserver::new_total(
                interface,
                700.0,
                position,
                &[source],
                limits(),
                &mut || false,
            )
            .is_err()
        );
    }
    assert!(
        PreparedFullWaveCompactSourceObserver::new_total(
            interface,
            f64::NAN,
            position,
            &[],
            limits(),
            &mut || false,
        )
        .is_err()
    );
    let empty = PreparedFullWaveCompactSourceObserver::new_total(
        interface,
        700.0,
        position,
        &[],
        limits(),
        &mut || false,
    )
    .unwrap();
    assert_eq!(empty.evaluate(&[]).unwrap().pressure_pa, C64::ZERO);
    assert_eq!(
        empty.evaluate_independent_sources(&[]).unwrap().scope,
        FULL_WAVE_TOTAL_ARRAY_SCOPE
    );
    assert_eq!(empty.evaluations(), 0);
}
