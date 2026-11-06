//! G1/G3/G4: physical volume cells and their original joint recorded spectra.

use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::quadrupole::PreparedQuadrupoleObserver;
use fs_aeroac::curle3d::quadrupole::recorded::{
    QuadrupoleRecordingBudget, QuadrupoleSampleClock, QuadrupoleWelchStream,
};
use fs_aeroac::curle3d::quadrupole::{QuadrupolePhasor, resolved_reynolds_moment_n_m};
use fs_aeroac::curle3d::refraction::PlanarFluidInterface;
use fs_aeroac::curle3d::refraction::full_wave::quadrupole::array::{
    FullWaveQuadrupoleField, FullWaveQuadrupoleSourceLimits, PreparedFullWaveQuadrupoleObserver,
};
use fs_aeroac::curle3d::refraction::full_wave::quadrupole::{
    FullWaveQuadrupoleLimits, PreparedFullWaveQuadrupole,
};
use fs_aeroac::curle3d::refraction::full_wave::reflection::PreparedReflectedFullWaveQuadrupole;
use fs_aeroac::curle3d::{CoherenceBudget, CompactFieldStatistics};
use fs_fft::welch::{WelchBudget, WelchConfig, WelchDetrend, WelchWindow};
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::c64::C64;

const POSITIONS: [[f64; 3]; 2] = [[-0.02, 0.01, -0.08], [0.035, -0.02, -0.06]];
const RECEIVER: [f64; 3] = [0.03, 0.02, 0.13];

fn interface(rho2: f64, c2: f64) -> PlanarFluidInterface {
    PlanarFluidInterface::new(
        FluidAcousticState::try_new(1.2, 340.0).unwrap(),
        FluidAcousticState::try_new(rho2, c2).unwrap(),
        [0.0; 3],
        [0.0, 0.0, 1.0],
    )
    .unwrap()
}
fn limits() -> FullWaveQuadrupoleSourceLimits {
    FullWaveQuadrupoleSourceLimits {
        max_sources: 8,
        maximum_total_evaluations: 1_048_576,
        quadrupole: FullWaveQuadrupoleLimits {
            absolute_pressure_tolerance_pa_per_n_m: 1e-6,
            absolute_velocity_tolerance_m_s_per_n_m: 1e-9,
            relative_tolerance: 1e-9,
            maximum_evaluations: 262_144,
            maximum_depth: 40,
        },
    }
}
fn coherence() -> CoherenceBudget {
    CoherenceBudget {
        max_modes: 8,
        max_source_mode_pairs: 64,
    }
}
fn moment(index: usize) -> QuadrupolePhasor {
    let a = (index + 1) as f64 * 1e-4;
    QuadrupolePhasor {
        stress_moment_n_m: [
            C64::new(a, 0.3 * a),
            C64::new(-0.2 * a, a),
            C64::from_re(0.7 * a),
            C64::new(0.4 * a, -0.5 * a),
            C64::from_re(-0.6 * a),
            C64::new(a, 0.1 * a),
        ],
    }
}
fn near(actual: C64, expected: C64) {
    assert!(
        (actual - expected).abs() <= 1e-12 * expected.abs() + 1e-14,
        "actual={actual:?}, expected={expected:?}"
    );
}
fn stats_scale(actual: &CompactFieldStatistics, original: &CompactFieldStatistics, p: f64, u: f64) {
    near(
        C64::from_re(actual.mean_square_pressure_pa2),
        C64::from_re(original.mean_square_pressure_pa2 * p * p),
    );
    for i in 0..3 {
        near(
            C64::new(
                actual.active_intensity_w_m2[i],
                actual.reactive_intensity_w_m2[i],
            ),
            C64::new(
                original.active_intensity_w_m2[i],
                original.reactive_intensity_w_m2[i],
            )
            .scale(p * u),
        );
        for j in 0..3 {
            near(
                actual.velocity_covariance_m2_s2[i][j],
                original.velocity_covariance_m2_s2[i][j].scale(u * u),
            );
        }
    }
}

#[test]
fn coherent_cells_match_native_tensor_fields_before_forming_flux() {
    for c2 in [340.0, 510.0] {
        let boundary = interface(2.4, c2);
        let prepared = PreparedFullWaveQuadrupoleObserver::new(
            boundary,
            700.0,
            RECEIVER,
            &POSITIONS,
            limits(),
            &mut || false,
        )
        .unwrap();
        assert_eq!(prepared.source_count(), 2);
        assert_eq!(prepared.source_position_m(1), Some(POSITIONS[1]));
        assert_eq!(prepared.source_position_m(2), None);
        assert_eq!(prepared.observer_position_m(), RECEIVER);
        let sources = [moment(0), moment(1)];
        let mut expected = [C64::ZERO; 4];
        let mut evaluations = 0;
        for (i, (&position, &source)) in POSITIONS.iter().zip(&sources).enumerate() {
            let original = PreparedFullWaveQuadrupole::new(
                boundary,
                700.0,
                position,
                RECEIVER,
                limits().quadrupole,
                &mut || false,
            )
            .unwrap();
            evaluations += original.diagnostics().evaluations;
            assert_eq!(
                prepared.source_diagnostics(i).unwrap().evaluations,
                original.diagnostics().evaluations
            );
            let field = original.evaluate(source).unwrap();
            expected[0] = expected[0] + field.pressure_pa;
            for j in 0..3 {
                expected[j + 1] = expected[j + 1] + field.particle_velocity_m_s[j];
            }
        }
        assert_eq!(prepared.evaluations(), evaluations);
        let field = prepared.evaluate(&sources).unwrap();
        let line = prepared.evaluate_harmonic_line(&sources).unwrap();
        assert_eq!(line.frequency_hz, 700.0);
        assert_eq!(line.pressure_pa, field.pressure_pa);
        near(field.pressure_pa, expected[0]);
        for j in 0..3 {
            near(field.particle_velocity_m_s[j], expected[j + 1]);
            let flux = expected[0].scale(0.5) * expected[j + 1].conj();
            near(
                C64::new(
                    field.active_intensity_w_m2[j],
                    field.reactive_intensity_w_m2[j],
                ),
                flux,
            );
        }
        let stats = prepared
            .evaluate_coherence_modes(&[&sources], coherence(), || false)
            .unwrap();
        near(
            C64::from_re(stats.mean_square_pressure_pa2),
            C64::from_re(0.5 * field.pressure_pa.abs().powi(2)),
        );
        for j in 0..3 {
            near(
                C64::new(
                    stats.active_intensity_w_m2[j],
                    stats.reactive_intensity_w_m2[j],
                ),
                C64::new(
                    field.active_intensity_w_m2[j],
                    field.reactive_intensity_w_m2[j],
                ),
            );
        }
    }
}

#[test]
fn tensor_and_cell_coherence_survive_pressure_nulls() {
    let prepared = PreparedFullWaveQuadrupoleObserver::new(
        interface(2.4, 510.0),
        700.0,
        RECEIVER,
        &[POSITIONS[0]; 2],
        limits(),
        &mut || false,
    )
    .unwrap();
    let first = moment(0);
    let opposite = QuadrupolePhasor {
        stress_moment_n_m: first.stress_moment_n_m.map(|q| q.scale(-1.0)),
    };
    let sources = [first, opposite];
    let coherent = prepared.evaluate(&sources).unwrap();
    assert_eq!(coherent.pressure_pa, C64::ZERO);
    assert_eq!(coherent.particle_velocity_m_s, [C64::ZERO; 3]);
    let mode = prepared
        .evaluate_coherence_modes(&[&sources], coherence(), || false)
        .unwrap();
    assert_eq!(mode.mean_square_pressure_pa2, 0.0);
    let independent = prepared.evaluate_independent(&sources, || false).unwrap();
    assert!(independent.mean_square_pressure_pa2 > 0.0);
    let single_sources = [first, QuadrupolePhasor::default()];
    let single = prepared
        .evaluate_coherence_modes(&[&single_sources], coherence(), || false)
        .unwrap();
    stats_scale(&independent, &single, 2.0_f64.sqrt(), 2.0_f64.sqrt());
    let second_sources = [QuadrupolePhasor::default(), opposite];
    let factors = prepared
        .evaluate_coherence_modes(&[&single_sources, &second_sources], coherence(), || false)
        .unwrap();
    stats_scale(&factors, &independent, 1.0, 1.0);
    let partial = prepared
        .evaluate_coherence_modes(&[&sources, &single_sources], coherence(), || false)
        .unwrap();
    stats_scale(&partial, &single, 1.0, 1.0);
}

#[test]
fn recorded_reynolds_stresses_keep_native_clock_geometry_and_welch_normalization() {
    let clock = QuadrupoleSampleClock::new(0.125, 1.0 / 8192.0, 0.0).unwrap();
    let config = WelchConfig {
        segment_samples: 64,
        hop_samples: 32,
        window: WelchWindow::PeriodicHann,
        detrend: WelchDetrend::SegmentMean,
    };
    let budget = QuadrupoleRecordingBudget {
        maximum_sources: 2,
        welch: WelchBudget {
            maximum_channels: 12,
            maximum_samples_per_channel: 96,
            maximum_transforms: 128,
            maximum_workspace_bytes: 1_000_000,
        },
    };
    let mut acquisition =
        QuadrupoleWelchStream::new(&POSITIONS, clock, 96, config, budget, || false).unwrap();
    for sample in 0..96 {
        let phase = core::f64::consts::TAU * sample as f64 / 64.0;
        let moments = core::array::from_fn::<_, 2, _>(|cell| {
            let phase = phase + cell as f64 * 0.17;
            resolved_reynolds_moment_n_m(
                1.2,
                0.002,
                [
                    0.4 + (3.0 * phase).cos(),
                    (2.0 * phase).sin(),
                    0.2 * (7.0 * phase).cos(),
                ],
            )
            .unwrap()
        });
        acquisition
            .push_sample(&POSITIONS, clock.time_at(sample).unwrap(), &moments, || {
                false
            })
            .unwrap();
    }
    let recorded = acquisition.finish().unwrap();
    let boundary = interface(2.4, 340.0);
    let refracted = recorded
        .prepare_refracted_observer(boundary, RECEIVER, 6, limits(), || false)
        .unwrap();
    assert_eq!(refracted.evaluations(), 0);
    assert_eq!(recorded.clock(), clock);
    assert_eq!(refracted.frequency_hz(), 768.0);
    let original = recorded
        .prepare_observer(
            FluidAcousticState::try_new(1.2, 340.0).unwrap(),
            RECEIVER,
            6,
            || false,
        )
        .unwrap();
    let expected = recorded
        .evaluate_bin(&original, 6, coherence(), || false)
        .unwrap();
    let actual = recorded
        .evaluate_refracted_bin(&refracted, 6, coherence(), || false)
        .unwrap();
    assert!(actual.mean_square_pressure_pa2 > 0.0);
    // Independent equal-speed transmission: p*=2rho2/(rho1+rho2),
    // u*=2rho1/(rho1+rho2). The source stress already contains rho and dV.
    stats_scale(&actual, &expected, 4.0 / 3.0, 2.0 / 3.0);
    let reordered = PreparedFullWaveQuadrupoleObserver::new(
        boundary,
        768.0,
        RECEIVER,
        &[POSITIONS[1], POSITIONS[0]],
        limits(),
        &mut || false,
    )
    .unwrap();
    assert!(
        recorded
            .evaluate_refracted_bin(&reordered, 6, coherence(), || false)
            .is_err()
    );
    assert!(
        recorded
            .evaluate_refracted_bin(&refracted, 7, coherence(), || false)
            .is_err()
    );
    assert!(
        recorded
            .prepare_refracted_observer(boundary, RECEIVER, 0, limits(), || false)
            .is_err()
    );
    let insufficient = CoherenceBudget {
        max_modes: 1,
        max_source_mode_pairs: 2,
    };
    assert!(
        recorded
            .evaluate_refracted_bin(&refracted, 6, insufficient, || false)
            .is_err()
    );
    assert!(matches!(
        recorded.evaluate_refracted_bin(&refracted, 6, coherence(), || true),
        Err(AeroacError::Cancelled { .. })
    ));
}

#[test]
fn aggregate_work_empty_geometry_and_late_cancellation_do_not_publish_partial_arrays() {
    let boundary = interface(2.4, 510.0);
    let one = PreparedFullWaveQuadrupoleObserver::new(
        boundary,
        700.0,
        RECEIVER,
        &POSITIONS[..1],
        limits(),
        &mut || false,
    )
    .unwrap();
    assert!(one.evaluations() > 0);
    let mut capped = limits();
    capped.maximum_total_evaluations = 2 * one.evaluations() - 1;
    assert!(
        PreparedFullWaveQuadrupoleObserver::new(
            boundary,
            700.0,
            RECEIVER,
            &[POSITIONS[0]; 2],
            capped,
            &mut || false,
        )
        .is_err()
    );
    capped = limits();
    capped.max_sources = 0;
    assert!(
        PreparedFullWaveQuadrupoleObserver::new(
            boundary,
            700.0,
            RECEIVER,
            &POSITIONS[..1],
            capped,
            &mut || false,
        )
        .is_err()
    );
    let empty = PreparedFullWaveQuadrupoleObserver::new(
        boundary,
        700.0,
        RECEIVER,
        &[],
        capped,
        &mut || false,
    )
    .unwrap();
    assert_eq!(empty.evaluations(), 0);
    assert_eq!(empty.evaluate(&[]).unwrap().pressure_pa, C64::ZERO);
    assert!(
        PreparedFullWaveQuadrupoleObserver::new(boundary, 0.0, RECEIVER, &[], capped, &mut || {
            false
        },)
        .is_err()
    );
    assert!(
        PreparedFullWaveQuadrupoleObserver::new(
            boundary,
            700.0,
            POSITIONS[0],
            &[],
            capped,
            &mut || false,
        )
        .is_err()
    );
    assert!(matches!(
        PreparedFullWaveQuadrupoleObserver::new(
            boundary,
            700.0,
            RECEIVER,
            &POSITIONS,
            limits(),
            &mut || true,
        ),
        Err(AeroacError::Cancelled { .. })
    ));
    let before = one.evaluate(&[moment(0)]).unwrap();
    let mut polls = 0;
    assert!(matches!(
        one.evaluate_cancellable(&[moment(0)], || {
            polls += 1;
            polls >= 3
        }),
        Err(AeroacError::Cancelled { .. })
    ));
    let mut bad = moment(0);
    bad.stress_moment_n_m[5] = C64::new(f64::NAN, 0.0);
    assert!(one.evaluate(&[bad]).is_err());
    assert!(one.evaluate(&[]).is_err());
    assert_eq!(
        one.evaluate(&[moment(0)]).unwrap().pressure_pa,
        before.pressure_pa
    );
}

#[test]
fn incident_arrays_keep_native_tensor_returns_and_direct_reflected_interference() {
    let receiver = [0.03, 0.02, -0.13];
    let sources = [moment(0), moment(1)];
    for c2 in [340.0, 510.0] {
        let boundary = interface(2.4, c2);
        let reflected = PreparedFullWaveQuadrupoleObserver::new_reflected(
            boundary,
            700.0,
            receiver,
            &POSITIONS,
            limits(),
            &mut || false,
        )
        .unwrap();
        let total = PreparedFullWaveQuadrupoleObserver::new_total(
            boundary,
            700.0,
            receiver,
            &POSITIONS,
            limits(),
            &mut || false,
        )
        .unwrap();
        assert_eq!(reflected.field_kind(), FullWaveQuadrupoleField::Reflected);
        assert_eq!(total.field_kind(), FullWaveQuadrupoleField::IncidentTotal);
        assert_ne!(reflected.scope(), total.scope());
        let direct = PreparedQuadrupoleObserver::new(
            boundary.incident_state(),
            700.0,
            receiver,
            &POSITIONS,
            2,
        )
        .unwrap()
        .evaluate(&sources)
        .unwrap();
        let mut expected_return = [C64::ZERO; 4];
        let mut evaluations = 0;
        for (&position, &source) in POSITIONS.iter().zip(&sources) {
            let native = PreparedReflectedFullWaveQuadrupole::new(
                boundary,
                700.0,
                position,
                receiver,
                limits().quadrupole,
                &mut || false,
            )
            .unwrap();
            evaluations += native.diagnostics().evaluations;
            let field = native.evaluate_reflected(source).unwrap();
            expected_return[0] = expected_return[0] + field.pressure_pa;
            for axis in 0..3 {
                expected_return[axis + 1] =
                    expected_return[axis + 1] + field.particle_velocity_m_s[axis];
            }
        }
        assert_eq!(reflected.evaluations(), evaluations);
        assert_eq!(total.evaluations(), evaluations);
        let returned = reflected.evaluate(&sources).unwrap();
        let combined = total.evaluate(&sources).unwrap();
        near(returned.pressure_pa, expected_return[0]);
        near(
            combined.pressure_pa,
            direct.pressure_pa + expected_return[0],
        );
        let statistics = total
            .evaluate_coherence_modes(&[&sources], coherence(), || false)
            .unwrap();
        near(
            C64::from_re(statistics.mean_square_pressure_pa2),
            C64::from_re(0.5 * combined.pressure_pa.abs().powi(2)),
        );
        let mut interference = 0.0_f64;
        for axis in 0..3 {
            let velocity = direct.particle_velocity_m_s[axis] + expected_return[axis + 1];
            near(
                returned.particle_velocity_m_s[axis],
                expected_return[axis + 1],
            );
            near(combined.particle_velocity_m_s[axis], velocity);
            let flux = (direct.pressure_pa + expected_return[0]).scale(0.5) * velocity.conj();
            near(
                C64::new(
                    combined.active_intensity_w_m2[axis],
                    combined.reactive_intensity_w_m2[axis],
                ),
                flux,
            );
            near(
                C64::new(
                    statistics.active_intensity_w_m2[axis],
                    statistics.reactive_intensity_w_m2[axis],
                ),
                flux,
            );
            interference = interference.max(
                (combined.active_intensity_w_m2[axis]
                    - direct.active_intensity_w_m2[axis]
                    - returned.active_intensity_w_m2[axis])
                    .abs(),
            );
        }
        assert!(interference > 1e-12);
        let line = total.evaluate_harmonic_line(&sources).unwrap();
        assert_eq!(line.frequency_hz, 700.0);
        assert_eq!(line.pressure_pa, combined.pressure_pa);
    }
}

#[test]
fn recorded_total_stresses_match_an_independent_direct_and_mirrored_cell_recording() {
    let receiver = [0.03, 0.02, -0.13];
    let image_positions = POSITIONS.map(|p| [p[0], p[1], -p[2]]);
    let all_positions = [
        POSITIONS[0],
        POSITIONS[1],
        image_positions[0],
        image_positions[1],
    ];
    let clock = QuadrupoleSampleClock::new(0.125, 1.0 / 8192.0, 0.0).unwrap();
    let config = WelchConfig {
        segment_samples: 64,
        hop_samples: 32,
        window: WelchWindow::PeriodicHann,
        detrend: WelchDetrend::SegmentMean,
    };
    let budget = QuadrupoleRecordingBudget {
        maximum_sources: 4,
        welch: WelchBudget {
            maximum_channels: 24,
            maximum_samples_per_channel: 96,
            maximum_transforms: 128,
            maximum_workspace_bytes: 1_000_000,
        },
    };
    let mut original =
        QuadrupoleWelchStream::new(&POSITIONS, clock, 96, config, budget, || false).unwrap();
    let mut image_reference =
        QuadrupoleWelchStream::new(&all_positions, clock, 96, config, budget, || false).unwrap();
    for sample in 0..96 {
        let phase = core::f64::consts::TAU * sample as f64 / 64.0;
        let moments: [[f64; 6]; 2] = core::array::from_fn(|cell| {
            core::array::from_fn(|entry| {
                1e-4 * ((entry + 1) as f64 * phase + cell as f64 * 0.23).cos()
            })
        });
        // Independent equal-speed plane image: S M S, S=diag(1,1,-1),
        // and R=(2.4-1.2)/(2.4+1.2)=1/3. Density and dV are already in M.
        let images = moments.map(|m| {
            core::array::from_fn(|entry| m[entry] * if entry >= 4 { -1.0 / 3.0 } else { 1.0 / 3.0 })
        });
        let all_moments = [moments[0], moments[1], images[0], images[1]];
        let time = clock.time_at(sample).unwrap();
        original
            .push_sample(&POSITIONS, time, &moments, || false)
            .unwrap();
        image_reference
            .push_sample(&all_positions, time, &all_moments, || false)
            .unwrap();
    }
    let original = original.finish().unwrap();
    let image_reference = image_reference.finish().unwrap();
    let boundary = interface(2.4, 340.0);
    let total = original
        .prepare_total_observer(boundary, receiver, 3, limits(), || false)
        .unwrap();
    assert_eq!(total.field_kind(), FullWaveQuadrupoleField::IncidentTotal);
    assert_eq!(total.evaluations(), 0);
    let expected_observer = image_reference
        .prepare_observer(boundary.incident_state(), receiver, 3, || false)
        .unwrap();
    let expected = image_reference
        .evaluate_bin(&expected_observer, 3, coherence(), || false)
        .unwrap();
    let measured = original
        .evaluate_refracted_bin(&total, 3, coherence(), || false)
        .unwrap();
    assert!(expected.mean_square_pressure_pa2 > 0.0);
    stats_scale(&measured, &expected, 1.0, 1.0);
    let absent = original
        .prepare_reflected_observer(interface(1.2, 340.0), receiver, 3, limits(), || false)
        .unwrap();
    let silent = original
        .evaluate_refracted_bin(&absent, 3, coherence(), || false)
        .unwrap();
    assert_eq!(silent.mean_square_pressure_pa2, 0.0);
    assert_eq!(silent.velocity_covariance_m2_s2, [[C64::ZERO; 3]; 3]);
    let reordered = PreparedFullWaveQuadrupoleObserver::new_total(
        boundary,
        total.frequency_hz(),
        receiver,
        &[POSITIONS[1], POSITIONS[0]],
        limits(),
        &mut || false,
    )
    .unwrap();
    assert!(
        original
            .evaluate_refracted_bin(&reordered, 3, coherence(), || false)
            .is_err()
    );
    assert!(
        original
            .evaluate_refracted_bin(&total, 4, coherence(), || false)
            .is_err()
    );
    assert!(matches!(
        original.evaluate_refracted_bin(&total, 3, coherence(), || true),
        Err(AeroacError::Cancelled { .. })
    ));
}

#[test]
fn incident_arrays_admit_finite_return_self_fields_and_enforce_original_work_and_sides() {
    let boundary = interface(2.4, 510.0);
    let receiver = POSITIONS[0];
    let reflected = PreparedFullWaveQuadrupoleObserver::new_reflected(
        boundary,
        700.0,
        receiver,
        &[receiver],
        limits(),
        &mut || false,
    )
    .unwrap();
    assert!(
        reflected
            .evaluate(&[moment(0)])
            .unwrap()
            .pressure_pa
            .abs()
            .is_finite()
    );
    assert!(
        PreparedFullWaveQuadrupoleObserver::new_total(
            boundary,
            700.0,
            receiver,
            &[receiver],
            limits(),
            &mut || false,
        )
        .is_err()
    );
    let mut capped = limits();
    capped.maximum_total_evaluations = 2 * reflected.evaluations() - 1;
    assert!(
        PreparedFullWaveQuadrupoleObserver::new_reflected(
            boundary,
            700.0,
            receiver,
            &[receiver; 2],
            capped,
            &mut || false,
        )
        .is_err()
    );
    for total in [false, true] {
        let prepare = |observer, positions: &[[f64; 3]]| {
            if total {
                PreparedFullWaveQuadrupoleObserver::new_total(
                    boundary,
                    700.0,
                    observer,
                    positions,
                    limits(),
                    &mut || false,
                )
            } else {
                PreparedFullWaveQuadrupoleObserver::new_reflected(
                    boundary,
                    700.0,
                    observer,
                    positions,
                    limits(),
                    &mut || false,
                )
            }
        };
        let empty = prepare(receiver, &[]).unwrap();
        assert_eq!(empty.evaluate(&[]).unwrap().pressure_pa, C64::ZERO);
        assert!(prepare(RECEIVER, &[]).is_err());
        assert!(prepare([0.0; 3], &[]).is_err());
        assert!(prepare(receiver, &[RECEIVER]).is_err());
        assert!(prepare(receiver, &[[f64::NAN, 0.0, -1.0]]).is_err());
    }
    let opposite = QuadrupolePhasor {
        stress_moment_n_m: moment(0).stress_moment_n_m.map(|m| m.scale(-1.0)),
    };
    let pair = PreparedFullWaveQuadrupoleObserver::new_reflected(
        boundary,
        700.0,
        receiver,
        &[receiver; 2],
        limits(),
        &mut || false,
    )
    .unwrap();
    assert_eq!(
        pair.evaluate(&[moment(0), opposite]).unwrap().pressure_pa,
        C64::ZERO
    );
    assert!(
        pair.evaluate_independent(&[moment(0), opposite], || false)
            .unwrap()
            .mean_square_pressure_pa2
            > 0.0
    );
    assert!(matches!(
        PreparedFullWaveQuadrupoleObserver::new_total(
            boundary,
            700.0,
            receiver,
            &POSITIONS[1..],
            limits(),
            &mut || true,
        ),
        Err(AeroacError::Cancelled { .. })
    ));
}
