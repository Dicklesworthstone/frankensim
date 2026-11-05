//! G1/G3/G4: original Q/F arrays retain coherent fields and statistics inside a layer.

use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::refraction::PlanarFluidInterface;
use fs_aeroac::curle3d::refraction::full_wave::array::{
    FullWaveCompactSourceLimits, PreparedFullWaveCompactSourceObserver,
};
use fs_aeroac::curle3d::refraction::full_wave::layer::PlanarFluidLayer;
use fs_aeroac::curle3d::refraction::full_wave::layer::interior::internal::array::{
    FULL_WAVE_LAYER_INTERNAL_REFLECTED_ARRAY_SCOPE, FULL_WAVE_LAYER_INTERNAL_TOTAL_ARRAY_SCOPE,
    PreparedLayerInternalFullWaveCompactSourceObserver,
};
use fs_aeroac::curle3d::refraction::full_wave::layer::interior::internal::{
    PreparedLayerInternalFullWaveForce, PreparedLayerInternalFullWaveMonopole,
};
use fs_aeroac::curle3d::refraction::full_wave::{FullWaveForceLimits, FullWaveMonopoleLimits};
use fs_aeroac::curle3d::{
    CoherenceBudget, CoherentCompactField, CompactFieldStatistics, CompactSourcePhasor,
    PreparedCompactSourceObserver,
};
use fs_aeroac::observation::{
    HarmonicObservationBudget, HarmonicSampleClock, PreparedHarmonicField,
};
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::{c64::C64, dd::Dd};

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
fn prepare(
    slab: PlanarFluidLayer,
    receiver: [f64; 3],
    positions: &[[f64; 3]],
    total: bool,
    budget: FullWaveCompactSourceLimits,
) -> Result<PreparedLayerInternalFullWaveCompactSourceObserver, AeroacError> {
    if total {
        PreparedLayerInternalFullWaveCompactSourceObserver::new(
            slab,
            700.0,
            receiver,
            positions,
            budget,
            &mut || false,
        )
    } else {
        PreparedLayerInternalFullWaveCompactSourceObserver::new_reflected(
            slab,
            700.0,
            receiver,
            positions,
            budget,
            &mut || false,
        )
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
fn scaled(source: CompactSourcePhasor, weight: C64) -> CompactSourcePhasor {
    CompactSourcePhasor {
        volume_flow_m3_s: source.volume_flow_m3_s * weight,
        force_on_fluid_n: source.force_on_fluid_n.map(|v| v * weight),
    }
}
fn close(actual: C64, expected: C64) {
    assert!(
        (actual - expected).abs() <= 4e-12 * expected.abs().max(1e-10),
        "actual={actual:?}, expected={expected:?}"
    );
}
fn values(field: &CoherentCompactField) -> [C64; 4] {
    [
        field.pressure_pa,
        field.particle_velocity_m_s[0],
        field.particle_velocity_m_s[1],
        field.particle_velocity_m_s[2],
    ]
}
fn field_close(actual: &CoherentCompactField, expected: &CoherentCompactField) {
    for (a, b) in values(actual).into_iter().zip(values(expected)) {
        close(a, b);
    }
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
    }
}
fn check_moments(actual: &CompactFieldStatistics, fields: &[[C64; 4]]) {
    let p2 = fields
        .iter()
        .fold(C64::ZERO, |sum, f| sum + f[0] * f[0].conj())
        .scale(0.5);
    close(C64::from_re(actual.mean_square_pressure_pa2), p2);
    for i in 0..3 {
        let intensity = fields
            .iter()
            .fold(C64::ZERO, |sum, f| sum + f[0] * f[i + 1].conj())
            .scale(0.5);
        close(
            C64::new(
                actual.active_intensity_w_m2[i],
                actual.reactive_intensity_w_m2[i],
            ),
            intensity,
        );
        for j in 0..3 {
            let covariance = fields
                .iter()
                .fold(C64::ZERO, |sum, f| sum + f[i + 1] * f[j + 1].conj())
                .scale(0.5);
            close(actual.velocity_covariance_m2_s2[i][j], covariance);
        }
    }
}

#[test]
fn original_homogeneous_and_single_interface_arrays_keep_source_units_and_clock() {
    let state = fluid(2.4, 300.0);
    let positions = [[0.0, 0.0, 0.025], [-0.01, 0.005, 0.04]];
    let receiver = [0.04, -0.012, 0.055];
    let mut budget = limits();
    budget.maximum_total_evaluations = 1;
    budget.monopole.maximum_evaluations = 1;
    budget.force.maximum_evaluations = 1;
    for exterior in [state, fluid(1.2, 300.0)] {
        let slab = layer([exterior, state, state], 0.08);
        let reversed =
            PlanarFluidInterface::new(state, exterior, [0.0; 3], [0.0, 0.0, -1.0]).unwrap();
        for total in [true, false] {
            let array = prepare(slab, receiver, &positions, total, budget).unwrap();
            let original = if total {
                PreparedFullWaveCompactSourceObserver::new_total(
                    reversed,
                    700.0,
                    receiver,
                    &positions,
                    budget,
                    &mut || false,
                )
            } else {
                PreparedFullWaveCompactSourceObserver::new_reflected(
                    reversed,
                    700.0,
                    receiver,
                    &positions,
                    budget,
                    &mut || false,
                )
            }
            .unwrap();
            let actual = array.evaluate(&sources()).unwrap();
            field_close(&actual, &original.evaluate(&sources()).unwrap());
            assert_eq!(array.evaluations(), 0);
            assert_eq!(array.includes_direct_field(), total);
            assert_eq!(array.layer().layer_state(), state);
            assert_eq!(array.observer_position_m(), receiver);
            assert_eq!(array.frequency_hz(), 700.0);
            assert_eq!(array.source_count(), 2);
            assert_eq!(array.source_position_m(1), Some(positions[1]));
            assert_eq!(array.source_position_m(2), None);
            assert!(array.source_diagnostics(2).is_none());
            for i in 0..2 {
                let d = array.source_diagnostics(i).unwrap();
                assert!(
                    d.monopole.reduced_to_single_interface && d.force.reduced_to_single_interface
                );
                assert_eq!(d.monopole.evaluations + d.force.evaluations, 0);
            }
            let expected_scope = if total {
                FULL_WAVE_LAYER_INTERNAL_TOTAL_ARRAY_SCOPE
            } else {
                FULL_WAVE_LAYER_INTERNAL_REFLECTED_ARRAY_SCOPE
            };
            assert_eq!(array.scope(), expected_scope);
            assert_eq!(actual.scope, expected_scope);
            let line = array.evaluate_harmonic_line(&sources()).unwrap();
            let original_line = original.evaluate_harmonic_line(&sources()).unwrap();
            assert_eq!(line.frequency_hz, 700.0);
            close(line.pressure_pa, original_line.pressure_pa);
            for i in 0..3 {
                close(
                    line.particle_velocity_m_s[i],
                    original_line.particle_velocity_m_s[i],
                );
            }
            // At one quarter period after an exact large epoch, exp(-i wt)
            // selects the original phasor's imaginary part. Propagation has
            // already occurred in the native field and must not be delayed again.
            let realization = PreparedHarmonicField::new(
                &[line],
                HarmonicObservationBudget {
                    max_lines: 1,
                    max_line_pairs: 1,
                },
            )
            .unwrap();
            let clock = HarmonicSampleClock::new(Dd::from_f64(1_073_741_824.0), 2800).unwrap();
            let sample = realization
                .field_at_clock_dd(clock.time_at(1).unwrap())
                .unwrap();
            close(
                C64::from_re(sample.pressure_pa),
                C64::from_re(original_line.pressure_pa.im),
            );
            for i in 0..3 {
                close(
                    C64::from_re(sample.particle_velocity_m_s[i]),
                    C64::from_re(original_line.particle_velocity_m_s[i].im),
                );
            }
            if exterior == state {
                if total {
                    let direct =
                        PreparedCompactSourceObserver::new(state, 700.0, receiver, &positions, 2)
                            .unwrap();
                    field_close(&actual, &direct.evaluate(&sources()).unwrap());
                } else {
                    assert_eq!(values(&actual), [C64::ZERO; 4]);
                }
            }
        }
    }
}

#[test]
fn two_face_arrays_retain_q_force_path_coherence_and_one_aggregate_budget() {
    let slab = layer(
        [fluid(1.2, 300.0), fluid(2.4, 600.0), fluid(1.8, 400.0)],
        0.018,
    );
    let positions = [[0.0, 0.0, 0.006], [-0.01, 0.005, 0.014]];
    let receiver = [0.04, -0.012, 0.011];
    let amplitudes = sources();
    let total = prepare(slab, receiver, &positions, true, limits()).unwrap();
    let reflected = prepare(slab, receiver, &positions, false, limits()).unwrap();
    let direct =
        PreparedCompactSourceObserver::new(slab.layer_state(), 700.0, receiver, &positions, 2)
            .unwrap();
    let mut fields = [[C64::ZERO; 4]; 2];
    let mut work = 0;
    for (i, (&position, &amplitude)) in positions.iter().zip(&amplitudes).enumerate() {
        let q = PreparedLayerInternalFullWaveMonopole::new(
            slab,
            700.0,
            position,
            receiver,
            limits().monopole,
            &mut || false,
        )
        .unwrap();
        let f = PreparedLayerInternalFullWaveForce::new(
            slab,
            700.0,
            position,
            receiver,
            limits().force,
            &mut || false,
        )
        .unwrap();
        let d = total.source_diagnostics(i).unwrap();
        assert_eq!(d.monopole.evaluations, q.diagnostics().evaluations);
        assert_eq!(d.force.evaluations, f.diagnostics().evaluations);
        assert!(!d.monopole.reduced_to_single_interface && !d.force.reduced_to_single_interface);
        work += d.monopole.evaluations + d.force.evaluations;
        let q = q
            .evaluate_total_harmonic_line(amplitude.volume_flow_m3_s)
            .unwrap();
        let f = f
            .evaluate_total_harmonic_line(amplitude.force_on_fluid_n)
            .unwrap();
        fields[i][0] = q.pressure_pa + f.pressure_pa;
        for axis in 0..3 {
            fields[i][axis + 1] = q.particle_velocity_m_s[axis] + f.particle_velocity_m_s[axis];
        }
    }
    let coherent: [C64; 4] = core::array::from_fn(|i| fields[0][i] + fields[1][i]);
    let actual = total.evaluate(&amplitudes).unwrap();
    for (a, b) in values(&actual).into_iter().zip(coherent) {
        close(a, b);
    }
    for axis in 0..3 {
        close(
            C64::new(
                actual.active_intensity_w_m2[axis],
                actual.reactive_intensity_w_m2[axis],
            ),
            (coherent[0] * coherent[axis + 1].conj()).scale(0.5),
        );
    }
    let returned = values(&reflected.evaluate(&amplitudes).unwrap());
    let free = values(&direct.evaluate(&amplitudes).unwrap());
    for i in 0..4 {
        close(coherent[i], free[i] + returned[i]);
    }
    let independent = total.evaluate_independent_sources(&amplitudes).unwrap();
    check_moments(&independent, &fields);
    let shared_weight = C64::new(0.3, 0.4);
    let other_weight = C64::new(-0.2, 0.1);
    let first_mode = [amplitudes[0], scaled(amplitudes[1], shared_weight)];
    let second_mode = [
        CompactSourcePhasor::default(),
        scaled(amplitudes[1], other_weight),
    ];
    let moments = total
        .evaluate_coherence_modes(
            &[&first_mode, &second_mode],
            CoherenceBudget {
                max_modes: 2,
                max_source_mode_pairs: 4,
            },
        )
        .unwrap();
    check_moments(
        &moments,
        &[
            core::array::from_fn(|i| fields[0][i] + fields[1][i] * shared_weight),
            fields[1].map(|v| v * other_weight),
        ],
    );
    assert_eq!(moments.scope, FULL_WAVE_LAYER_INTERNAL_TOTAL_ARRAY_SCOPE);
    // Independent physical sources do not make each source's direct and
    // returning histories independent; those interference terms must survive.
    let wrong_path_power = direct
        .evaluate_independent_sources(&amplitudes)
        .unwrap()
        .mean_square_pressure_pa2
        + reflected
            .evaluate_independent_sources(&amplitudes)
            .unwrap()
            .mean_square_pressure_pa2;
    assert!(
        (independent.mean_square_pressure_pa2 - wrong_path_power).abs()
            > 1e-5 * independent.mean_square_pressure_pa2
    );
    assert_eq!(total.evaluations(), work);
    assert_eq!(reflected.evaluations(), work);
    assert!(work > 0 && work <= limits().maximum_total_evaluations);
    for total in [true, false] {
        assert!(
            prepare(
                slab,
                receiver,
                &positions,
                total,
                FullWaveCompactSourceLimits {
                    maximum_total_evaluations: work - 1,
                    ..limits()
                },
            )
            .is_err()
        );
    }
}

#[test]
fn reflected_self_fields_and_coherent_nulls_remain_distinct_from_independence() {
    let state = fluid(2.4, 300.0);
    let slab = layer([fluid(1.2, 300.0), state, state], 0.08);
    let point = [0.0, 0.0, 0.025];
    let amplitudes = sources();
    assert!(prepare(slab, point, &[point], true, limits()).is_err());
    let self_field = prepare(slab, point, &[point], false, limits()).unwrap();
    let response = self_field.evaluate(&amplitudes[..1]).unwrap();
    assert!(
        values(&response)
            .iter()
            .all(|v| v.re.is_finite() && v.im.is_finite())
    );
    assert!(response.pressure_pa.abs() > 0.0);
    assert_eq!(
        response.scope,
        FULL_WAVE_LAYER_INTERNAL_REFLECTED_ARRAY_SCOPE
    );
    for total in [true, false] {
        let array = prepare(slab, [0.035, 0.01, 0.05], &[point, point], total, limits()).unwrap();
        let first = amplitudes[0];
        let opposite = scaled(first, C64::from_re(-1.0));
        let null = array.evaluate(&[first, opposite]).unwrap();
        assert_eq!(values(&null), [C64::ZERO; 4]);
        assert_eq!(null.active_intensity_w_m2, [0.0; 3]);
        assert!(
            array
                .evaluate_independent_sources(&[first, opposite])
                .unwrap()
                .mean_square_pressure_pa2
                > 0.0
        );
        assert!(array.evaluate(&[first]).is_err());
        assert!(
            array
                .evaluate(&[
                    first,
                    CompactSourcePhasor::monopole(C64::new(f64::NAN, 0.0))
                ])
                .is_err()
        );
        assert!(
            array
                .evaluate_coherence_modes(
                    &[&[first, opposite]],
                    CoherenceBudget {
                        max_modes: 1,
                        max_source_mode_pairs: 1
                    },
                )
                .is_err()
        );
        assert!(
            array
                .evaluate_coherence_modes(
                    &[&[first]],
                    CoherenceBudget {
                        max_modes: 1,
                        max_source_mode_pairs: 2
                    },
                )
                .is_err()
        );
    }
}

#[test]
fn strict_dd_domains_empty_arrays_source_limits_and_cancellation_are_admitted() {
    let state = fluid(1.2, 300.0);
    let slab = layer([state; 3], 0.08);
    let source = [0.0, 0.0, 0.025];
    let receiver = [0.035, 0.01, 0.05];
    let outside = [
        [0.0; 3],
        [0.0, 0.0, -0.01],
        [0.0, 0.0, 0.08],
        [0.0, 0.0, 0.09],
        [f64::NAN, 0.0, 0.01],
    ];
    for total in [true, false] {
        let empty = prepare(slab, receiver, &[], total, limits()).unwrap();
        assert_eq!(values(&empty.evaluate(&[]).unwrap()), [C64::ZERO; 4]);
        assert_eq!(empty.evaluations(), 0);
        for bad in outside {
            assert!(prepare(slab, bad, &[], total, limits()).is_err());
            assert!(prepare(slab, receiver, &[bad], total, limits()).is_err());
        }
        assert!(prepare(layer([state; 3], 0.0), receiver, &[], total, limits()).is_err());
        for budget in [
            FullWaveCompactSourceLimits {
                max_sources: 0,
                ..limits()
            },
            FullWaveCompactSourceLimits {
                maximum_total_evaluations: 0,
                ..limits()
            },
        ] {
            assert!(prepare(slab, receiver, &[source], total, budget).is_err());
        }
        // Two unequal-density faces still need quadrature at equal speeds.
        let density_layer = layer([state, fluid(5.0, 300.0), fluid(2.4, 300.0)], 0.08);
        assert!(
            prepare(
                density_layer,
                receiver,
                &[source],
                total,
                FullWaveCompactSourceLimits {
                    maximum_total_evaluations: 1,
                    ..limits()
                },
            )
            .is_err()
        );
        // DD admission is necessary even with no scalar/force source owners.
        let diagonal = PlanarFluidLayer::new(state, state, state, [0.0; 3], [1.0; 3], 1.0).unwrap();
        let interior = [1e16, 1.0, -1e16];
        assert!(prepare(diagonal, interior, &[], total, limits()).is_ok());
        assert!(prepare(diagonal, interior.map(|v| -v), &[], total, limits()).is_err());
    }
    for frequency in [0.0, f64::NAN] {
        assert!(
            PreparedLayerInternalFullWaveCompactSourceObserver::new(
                slab,
                frequency,
                receiver,
                &[],
                limits(),
                &mut || false,
            )
            .is_err()
        );
    }
    let mut polls = 0;
    assert!(matches!(
        PreparedLayerInternalFullWaveCompactSourceObserver::new(
            slab,
            700.0,
            receiver,
            &[source],
            limits(),
            &mut || {
                polls += 1;
                polls >= 4
            },
        ),
        Err(AeroacError::Cancelled { .. })
    ));
    assert!(matches!(
        PreparedLayerInternalFullWaveCompactSourceObserver::new_reflected(
            slab,
            700.0,
            receiver,
            &[source],
            limits(),
            &mut || true,
        ),
        Err(AeroacError::Cancelled { .. })
    ));
}
