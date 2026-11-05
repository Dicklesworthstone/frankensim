//! G1/G3: physical source coherence and bounded work through the native fluid layer.

use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::refraction::PlanarFluidInterface;
use fs_aeroac::curle3d::refraction::full_wave::array::{
    FullWaveCompactFieldKind, FullWaveCompactSourceLimits, PreparedFullWaveCompactSourceObserver,
};
use fs_aeroac::curle3d::refraction::full_wave::layer::array::{
    FULL_WAVE_LAYER_ARRAY_SCOPE, FULL_WAVE_LAYER_REFLECTED_ARRAY_SCOPE,
    FULL_WAVE_LAYER_TOTAL_ARRAY_SCOPE, PreparedLayeredFullWaveCompactSourceObserver,
};
use fs_aeroac::curle3d::refraction::full_wave::layer::reflection::{
    PreparedLayeredReflectedFullWaveForce, PreparedLayeredReflectedFullWaveMonopole,
};
use fs_aeroac::curle3d::refraction::full_wave::layer::{
    PlanarFluidLayer, PreparedLayeredFullWaveForce, PreparedLayeredFullWaveMonopole,
};
use fs_aeroac::curle3d::refraction::full_wave::{FullWaveForceLimits, FullWaveMonopoleLimits};
use fs_aeroac::curle3d::{CoherenceBudget, CompactSourcePhasor};
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::c64::C64;

fn fluid(rho: f64, c: f64) -> FluidAcousticState {
    FluidAcousticState::try_new(rho, c).unwrap()
}
fn layer(states: [FluidAcousticState; 3], thickness: f64) -> PlanarFluidLayer {
    PlanarFluidLayer::new(
        states[0],
        states[1],
        states[2],
        [0.0; 3],
        [0.0, 0.0, 1.0],
        thickness,
    )
    .unwrap()
}
fn limits() -> FullWaveCompactSourceLimits {
    FullWaveCompactSourceLimits {
        max_sources: 4,
        maximum_total_evaluations: 1_048_576,
        monopole: FullWaveMonopoleLimits {
            absolute_pressure_tolerance_pa_per_m3_s: 1e-7,
            absolute_velocity_tolerance_per_m2: 1e-10,
            relative_tolerance: 1e-9,
            maximum_evaluations: 262_144,
            maximum_depth: 36,
        },
        force: FullWaveForceLimits {
            absolute_pressure_tolerance_pa_per_n: 1e-8,
            absolute_velocity_tolerance_m_s_per_n: 1e-11,
            relative_tolerance: 1e-9,
            maximum_evaluations: 262_144,
            maximum_depth: 36,
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
fn close(actual: C64, expected: C64) {
    assert!(
        (actual - expected).abs() <= 4e-12 * expected.abs().max(1e-12),
        "actual={actual:?}, expected={expected:?}"
    );
}

#[test]
fn absent_and_identical_layers_retain_original_array_values_and_exact_work() {
    let positions = [[0.0, 0.0, -0.025], [-0.01, 0.005, -0.04]];
    let receiver = [0.04, -0.01, 0.06];
    let first = fluid(1.2, 300.0);
    for slab in [
        layer([first, fluid(1000.0, 100.0), fluid(2.4, 300.0)], 0.0),
        layer([first; 3], 0.02),
    ] {
        let mut budget = limits();
        budget.maximum_total_evaluations = 1;
        budget.monopole.maximum_evaluations = 1;
        budget.force.maximum_evaluations = 1;
        let array = PreparedLayeredFullWaveCompactSourceObserver::new(
            slab,
            700.0,
            receiver,
            &positions,
            budget,
            &mut || false,
        )
        .unwrap();
        let interface = PlanarFluidInterface::new(
            slab.incident_state(),
            slab.receiving_state(),
            slab.first_face_point_m(),
            slab.receiving_normal(),
        )
        .unwrap();
        let original = PreparedFullWaveCompactSourceObserver::new(
            interface,
            700.0,
            receiver,
            &positions,
            budget,
            &mut || false,
        )
        .unwrap();
        let actual = array.evaluate(&sources()).unwrap();
        let expected = original.evaluate(&sources()).unwrap();
        assert_eq!(actual.pressure_pa, expected.pressure_pa);
        assert_eq!(actual.particle_velocity_m_s, expected.particle_velocity_m_s);
        assert_eq!(actual.active_intensity_w_m2, expected.active_intensity_w_m2);
        assert_eq!(actual.scope, FULL_WAVE_LAYER_ARRAY_SCOPE);
        let a = array.evaluate_independent_sources(&sources()).unwrap();
        let b = original.evaluate_independent_sources(&sources()).unwrap();
        assert_eq!(a.mean_square_pressure_pa2, b.mean_square_pressure_pa2);
        assert_eq!(a.velocity_covariance_m2_s2, b.velocity_covariance_m2_s2);
        assert_eq!(a.active_intensity_w_m2, b.active_intensity_w_m2);
        assert_eq!(array.evaluations(), 0);
        assert_eq!(array.layer().layer_state(), slab.layer_state());
        assert_eq!(array.layer().thickness_m(), slab.thickness_m());
        assert_eq!(array.source_count(), 2);
        assert_eq!(array.source_position_m(1), Some(positions[1]));
        assert_eq!(array.source_position_m(2), None);
        for i in 0..2 {
            let d = array.source_diagnostics(i).unwrap();
            assert!(d.monopole.reduced_to_single_interface && d.force.reduced_to_single_interface);
        }
    }
}

#[test]
fn finite_layer_transfers_keep_cross_mechanism_phase_all_moments_and_global_work() {
    let slab = layer(
        [fluid(1.2, 300.0), fluid(2.4, 600.0), fluid(1.8, 400.0)],
        0.018,
    );
    let positions = [[0.0, 0.0, -0.025], [-0.01, 0.005, -0.04]];
    let receiver = [0.04, -0.012, 0.041];
    let amplitudes = sources();
    let array = PreparedLayeredFullWaveCompactSourceObserver::new(
        slab,
        700.0,
        receiver,
        &positions,
        limits(),
        &mut || false,
    )
    .unwrap();
    let mut fields = [[C64::ZERO; 4]; 2];
    let mut work = 0;
    for (index, (&position, &amplitude)) in positions.iter().zip(&amplitudes).enumerate() {
        let q = PreparedLayeredFullWaveMonopole::new(
            slab,
            700.0,
            position,
            receiver,
            limits().monopole,
            &mut || false,
        )
        .unwrap();
        let f = PreparedLayeredFullWaveForce::new(
            slab,
            700.0,
            position,
            receiver,
            limits().force,
            &mut || false,
        )
        .unwrap();
        work += q.diagnostics().evaluations + f.diagnostics().evaluations;
        let d = array.source_diagnostics(index).unwrap();
        assert_eq!(d.monopole.evaluations, q.diagnostics().evaluations);
        assert_eq!(d.force.evaluations, f.diagnostics().evaluations);
        assert!(!d.monopole.reduced_to_single_interface && !d.force.reduced_to_single_interface);
        let q = q
            .evaluate_harmonic_line(amplitude.volume_flow_m3_s)
            .unwrap();
        let f = f
            .evaluate_harmonic_line(amplitude.force_on_fluid_n)
            .unwrap();
        fields[index][0] = q.pressure_pa + f.pressure_pa;
        for i in 0..3 {
            fields[index][i + 1] = q.particle_velocity_m_s[i] + f.particle_velocity_m_s[i];
        }
    }
    let actual = array.evaluate(&amplitudes).unwrap();
    close(actual.pressure_pa, fields[0][0] + fields[1][0]);
    for i in 0..3 {
        close(
            actual.particle_velocity_m_s[i],
            fields[0][i + 1] + fields[1][i + 1],
        );
    }
    let independent = array.evaluate_independent_sources(&amplitudes).unwrap();
    close(
        C64::from_re(independent.mean_square_pressure_pa2),
        (fields[0][0] * fields[0][0].conj() + fields[1][0] * fields[1][0].conj()).scale(0.5),
    );
    for i in 0..3 {
        close(
            C64::new(
                independent.active_intensity_w_m2[i],
                independent.reactive_intensity_w_m2[i],
            ),
            (fields[0][0] * fields[0][i + 1].conj() + fields[1][0] * fields[1][i + 1].conj())
                .scale(0.5),
        );
        for j in 0..3 {
            close(
                independent.velocity_covariance_m2_s2[i][j],
                (fields[0][i + 1] * fields[0][j + 1].conj()
                    + fields[1][i + 1] * fields[1][j + 1].conj())
                .scale(0.5),
            );
        }
    }
    let factor = array
        .evaluate_coherence_modes(
            &[&amplitudes],
            CoherenceBudget {
                max_modes: 1,
                max_source_mode_pairs: 2,
            },
        )
        .unwrap();
    close(
        C64::from_re(factor.mean_square_pressure_pa2),
        (actual.pressure_pa * actual.pressure_pa.conj()).scale(0.5),
    );
    assert_eq!(factor.scope, FULL_WAVE_LAYER_ARRAY_SCOPE);
    let line = array.evaluate_harmonic_line(&amplitudes).unwrap();
    assert_eq!(line.frequency_hz, 700.0);
    assert_eq!(line.pressure_pa, actual.pressure_pa);
    assert_eq!(line.particle_velocity_m_s, actual.particle_velocity_m_s);
    assert_eq!(array.evaluations(), work);
    assert!(work > 0 && work <= limits().maximum_total_evaluations);
    assert!(
        PreparedLayeredFullWaveCompactSourceObserver::new(
            slab,
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
fn equal_speeds_do_not_erase_density_layer_work_and_empty_arrays_still_validate_domain() {
    let slab = layer(
        [fluid(1.2, 300.0), fluid(5.0, 300.0), fluid(2.4, 300.0)],
        0.08,
    );
    let position = [0.0, 0.0, -0.025];
    let receiver = [0.035, 0.01, 0.098];
    assert!(
        PreparedLayeredFullWaveCompactSourceObserver::new(
            slab,
            700.0,
            receiver,
            &[position],
            FullWaveCompactSourceLimits {
                maximum_total_evaluations: 1,
                ..limits()
            },
            &mut || false,
        )
        .is_err()
    );
    let empty = PreparedLayeredFullWaveCompactSourceObserver::new(
        slab,
        700.0,
        receiver,
        &[],
        limits(),
        &mut || false,
    )
    .unwrap();
    assert_eq!(empty.evaluate(&[]).unwrap().pressure_pa, C64::ZERO);
    assert_eq!(empty.evaluations(), 0);
    assert!(empty.evaluate(&sources()).is_err());
    assert!(
        empty
            .evaluate_coherence_modes(
                &[&[]],
                CoherenceBudget {
                    max_modes: 0,
                    max_source_mode_pairs: 0
                }
            )
            .is_err()
    );
    for observer in [
        [0.0, 0.0, 0.08],
        [0.0, 0.0, 0.04],
        [0.0, 0.0, -0.01],
        [f64::NAN, 0.0, 0.1],
    ] {
        assert!(
            PreparedLayeredFullWaveCompactSourceObserver::new(
                slab,
                700.0,
                observer,
                &[],
                limits(),
                &mut || false,
            )
            .is_err()
        );
    }
    for source in [[0.0; 3], [0.0, 0.0, 0.04], [f64::NAN, 0.0, -0.01]] {
        assert!(
            PreparedLayeredFullWaveCompactSourceObserver::new(
                slab,
                700.0,
                receiver,
                &[source],
                limits(),
                &mut || false,
            )
            .is_err()
        );
    }
    assert!(
        PreparedLayeredFullWaveCompactSourceObserver::new(
            slab,
            0.0,
            receiver,
            &[],
            limits(),
            &mut || false,
        )
        .is_err()
    );
    assert!(
        PreparedLayeredFullWaveCompactSourceObserver::new(
            slab,
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
    assert!(matches!(
        PreparedLayeredFullWaveCompactSourceObserver::new(
            slab,
            700.0,
            receiver,
            &[position],
            limits(),
            &mut || true,
        ),
        Err(AeroacError::Cancelled { .. })
    ));
}

#[test]
fn incident_side_collapsed_or_invisible_second_faces_keep_exact_original_arrays() {
    let first = fluid(1.2, 300.0);
    let last = fluid(2.4, 300.0);
    let positions = [[0.0, 0.0, -0.025], [-0.01, 0.005, -0.04]];
    let receiver = [0.035, -0.012, -0.018];
    let interface = PlanarFluidInterface::new(first, last, [0.0; 3], [0.0, 0.0, 1.0]).unwrap();
    let mut budget = limits();
    budget.maximum_total_evaluations = 1;
    budget.monopole.maximum_evaluations = 1;
    budget.force.maximum_evaluations = 1;
    for slab in [
        layer([first, fluid(1000.0, 100.0), last], 0.0),
        layer([first, last, last], 0.08),
    ] {
        let reflected = PreparedLayeredFullWaveCompactSourceObserver::new_reflected(
            slab,
            700.0,
            receiver,
            &positions,
            budget,
            &mut || false,
        )
        .unwrap();
        let total = PreparedLayeredFullWaveCompactSourceObserver::new_total(
            slab,
            700.0,
            receiver,
            &positions,
            budget,
            &mut || false,
        )
        .unwrap();
        let original_reflected = PreparedFullWaveCompactSourceObserver::new_reflected(
            interface,
            700.0,
            receiver,
            &positions,
            budget,
            &mut || false,
        )
        .unwrap();
        let original_total = PreparedFullWaveCompactSourceObserver::new_total(
            interface,
            700.0,
            receiver,
            &positions,
            budget,
            &mut || false,
        )
        .unwrap();
        assert_eq!(reflected.field_kind(), FullWaveCompactFieldKind::Reflected);
        assert_eq!(total.field_kind(), FullWaveCompactFieldKind::IncidentTotal);
        for (array, original, scope) in [
            (
                &reflected,
                &original_reflected,
                FULL_WAVE_LAYER_REFLECTED_ARRAY_SCOPE,
            ),
            (&total, &original_total, FULL_WAVE_LAYER_TOTAL_ARRAY_SCOPE),
        ] {
            let actual = array.evaluate(&sources()).unwrap();
            let expected = original.evaluate(&sources()).unwrap();
            assert_eq!(actual.pressure_pa, expected.pressure_pa);
            assert_eq!(actual.particle_velocity_m_s, expected.particle_velocity_m_s);
            assert_eq!(actual.active_intensity_w_m2, expected.active_intensity_w_m2);
            assert_eq!(actual.scope, scope);
            let a = array.evaluate_independent_sources(&sources()).unwrap();
            let b = original.evaluate_independent_sources(&sources()).unwrap();
            assert_eq!(a.mean_square_pressure_pa2, b.mean_square_pressure_pa2);
            assert_eq!(a.velocity_covariance_m2_s2, b.velocity_covariance_m2_s2);
            assert_eq!(a.scope, scope);
            assert_eq!(array.evaluations(), 0);
            assert!(
                array
                    .source_diagnostics(0)
                    .unwrap()
                    .force
                    .reduced_to_single_interface
            );
        }
        let self_field = PreparedLayeredFullWaveCompactSourceObserver::new_reflected(
            slab,
            700.0,
            receiver,
            &[receiver],
            budget,
            &mut || false,
        )
        .unwrap();
        assert!(
            self_field
                .evaluate(&sources()[..1])
                .unwrap()
                .pressure_pa
                .abs()
                > 0.0
        );
        assert!(
            PreparedLayeredFullWaveCompactSourceObserver::new_total(
                slab,
                700.0,
                receiver,
                &[receiver],
                budget,
                &mut || false,
            )
            .is_err()
        );
        for observer in [[0.0; 3], [0.0, 0.0, 0.01]] {
            assert!(
                PreparedLayeredFullWaveCompactSourceObserver::new_total(
                    slab,
                    700.0,
                    observer,
                    &[],
                    budget,
                    &mut || false,
                )
                .is_err()
            );
        }
    }
}

#[test]
fn complete_layer_reflection_and_original_direct_fields_stay_coherent_per_source() {
    let slab = layer(
        [fluid(1.2, 300.0), fluid(2.4, 600.0), fluid(1.8, 400.0)],
        0.018,
    );
    let positions = [[0.0, 0.0, -0.025], [-0.01, 0.005, -0.04]];
    let receiver = [0.04, -0.012, -0.023];
    let amplitudes = sources();
    let reflected = PreparedLayeredFullWaveCompactSourceObserver::new_reflected(
        slab,
        700.0,
        receiver,
        &positions,
        limits(),
        &mut || false,
    )
    .unwrap();
    let total = PreparedLayeredFullWaveCompactSourceObserver::new_total(
        slab,
        700.0,
        receiver,
        &positions,
        limits(),
        &mut || false,
    )
    .unwrap();
    let mut reflected_fields = [[C64::ZERO; 4]; 2];
    let mut total_fields = [[C64::ZERO; 4]; 2];
    let mut work = 0;
    for (index, (&position, &amplitude)) in positions.iter().zip(&amplitudes).enumerate() {
        let q = PreparedLayeredReflectedFullWaveMonopole::new(
            slab,
            700.0,
            position,
            receiver,
            limits().monopole,
            &mut || false,
        )
        .unwrap();
        let f = PreparedLayeredReflectedFullWaveForce::new(
            slab,
            700.0,
            position,
            receiver,
            limits().force,
            &mut || false,
        )
        .unwrap();
        work += q.diagnostics().evaluations + f.diagnostics().evaluations;
        for (fields, q, f) in [
            (
                &mut reflected_fields,
                q.evaluate_reflected_harmonic_line(amplitude.volume_flow_m3_s)
                    .unwrap(),
                f.evaluate_reflected_harmonic_line(amplitude.force_on_fluid_n)
                    .unwrap(),
            ),
            (
                &mut total_fields,
                q.evaluate_total_harmonic_line(amplitude.volume_flow_m3_s)
                    .unwrap(),
                f.evaluate_total_harmonic_line(amplitude.force_on_fluid_n)
                    .unwrap(),
            ),
        ] {
            fields[index][0] = q.pressure_pa + f.pressure_pa;
            for i in 0..3 {
                fields[index][i + 1] = q.particle_velocity_m_s[i] + f.particle_velocity_m_s[i];
            }
        }
    }
    for (array, fields) in [(&reflected, reflected_fields), (&total, total_fields)] {
        let field = array.evaluate(&amplitudes).unwrap();
        close(field.pressure_pa, fields[0][0] + fields[1][0]);
        for i in 0..3 {
            close(
                field.particle_velocity_m_s[i],
                fields[0][i + 1] + fields[1][i + 1],
            );
        }
        let independent = array.evaluate_independent_sources(&amplitudes).unwrap();
        close(
            C64::from_re(independent.mean_square_pressure_pa2),
            (fields[0][0] * fields[0][0].conj() + fields[1][0] * fields[1][0].conj()).scale(0.5),
        );
        let zero = CompactSourcePhasor::default();
        let factored = array
            .evaluate_coherence_modes(
                &[&[amplitudes[0], zero], &[zero, amplitudes[1]]],
                CoherenceBudget {
                    max_modes: 2,
                    max_source_mode_pairs: 4,
                },
            )
            .unwrap();
        assert_eq!(
            factored.mean_square_pressure_pa2,
            independent.mean_square_pressure_pa2
        );
        assert_eq!(
            factored.velocity_covariance_m2_s2,
            independent.velocity_covariance_m2_s2
        );
        assert_eq!(
            factored.active_intensity_w_m2,
            independent.active_intensity_w_m2
        );
        let line = array.evaluate_harmonic_line(&amplitudes).unwrap();
        assert_eq!(line.frequency_hz, 700.0);
        assert_eq!(line.pressure_pa, field.pressure_pa);
        assert_eq!(line.particle_velocity_m_s, field.particle_velocity_m_s);
        assert_eq!(array.evaluations(), work);
    }
    assert!(work > 0);
    assert!(
        PreparedLayeredFullWaveCompactSourceObserver::new_total(
            slab,
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
