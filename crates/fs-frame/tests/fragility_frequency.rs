//! G1/G3/G4: one physical spectrum and horizon across fragility time levels.

use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_frame::{
    FragilityControls, FragilityError, GroundMotion, HistoryError, HistoryLimits, StoryFrame,
    StoryParams, e_stopped_fragility, e_stopped_fragility_with_spectrum, try_losses,
};
use fs_qty::{Dims, QtyAny};
use fs_scenario::{RealizationBudget, SpectralDiscretization, SpectrumModel, StochasticEnsemble};

const TIME: Dims = Dims([0, 0, 1, 0, 0, 0]);
const RATE: Dims = Dims([0, 0, -1, 0, 0, 0]);

fn ensemble() -> StochasticEnsemble {
    StochasticEnsemble {
        name: "physical-fragility".into(),
        seed: 29,
        members: 12,
        duration: QtyAny::new(0.5, TIME),
        dt: QtyAny::new(0.03125, TIME),
        model: SpectrumModel::KanaiTajimi {
            s0: 1e-6,
            omega_g: QtyAny::new(8.0, RATE),
            zeta_g: 0.4,
        },
    }
}

fn grid() -> SpectralDiscretization {
    SpectralDiscretization {
        omega_max_rad_s: 16.0,
        frequency_lines: 8,
    }
}

fn controls() -> FragilityControls {
    FragilityControls {
        drift_limit: 1e-5,
        alpha: 0.05,
        margin: 0.0,
        realization_budget: RealizationBudget {
            max_samples: 128,
            max_work: 8_192,
        },
        history_limits: HistoryLimits::new(128, 30, 1e-12, 1.0),
    }
}

fn with_cx<R>(cancelled: bool, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    if cancelled {
        gate.request();
    }
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 29,
                kernel_id: 1,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-12 * expected.abs() + 1e-22,
        "{actual} != {expected}"
    );
}

#[test]
fn g3_spectral_fragility_refinement_preserves_members_and_horizon() {
    let ensemble = ensemble();
    let unchanged = ensemble.clone();
    let params = StoryParams::default();
    let mut controls = controls();
    let spectrum = grid();
    let peaks: Vec<f64> = (0..ensemble.members)
        .map(|member| {
            let record = ensemble
                .realize_with_spectrum(member, spectrum, controls.realization_budget)
                .unwrap();
            StoryFrame::new(params)
                .run_from_time_zero(
                    GroundMotion::new(&record.values, ensemble.dt.value),
                    controls.history_limits,
                )
                .unwrap()
                .peak_abs_displacement_m
                / params.h
        })
        .collect();
    let mut ordered = peaks.clone();
    ordered.sort_by(f64::total_cmp);
    controls.drift_limit = f64::midpoint(ordered[5], ordered[6]);
    assert!(
        ordered[5] < ordered[6],
        "fixture must discriminate exceedances"
    );
    let report = with_cx(false, |cx| {
        e_stopped_fragility_with_spectrum(&ensemble, params, spectrum, controls, cx).unwrap()
    });
    assert_eq!(ensemble, unchanged);
    assert_eq!(report.members_used, ensemble.members);
    assert_eq!(report.exceedances, 6);
    assert_eq!(report.p_hat, 0.5);
    assert!(!report.stopped_early);
    assert_eq!(report.spectrum, Some(spectrum));
    assert_eq!(report.end_time_s, 15.0 * ensemble.dt.value);
    assert_eq!(report.mlmc.levels.len(), 2);
    assert!(report.mlmc.levels.iter().all(|level| level.samples == 12));
    // Fifteen fine advances; eight coarse advances, including the final dt.
    assert_eq!(report.mlmc.total_cost, 12.0 * (15.0 + 8.0));
    let mean = peaks.iter().sum::<f64>() / 12.0;
    let variance = peaks.iter().map(|peak| (peak - mean).powi(2)).sum::<f64>() / 11.0;
    near(report.mlmc.estimate, mean);
    near(report.mlmc.estimator_variance, variance / 12.0);
    near(
        report.mlmc.levels.iter().map(|level| level.mean).sum(),
        mean,
    );

    // Keep the same physical endpoint when refining the observation clock.
    // The new coarse history is then bit-for-bit the previous fine history:
    // neither its spectral lines nor its initial forcing have changed.
    let mut refined = ensemble.clone();
    refined.dt.value *= 0.5;
    refined.duration.value = 31.0 * refined.dt.value;
    let finer = with_cx(false, |cx| {
        e_stopped_fragility_with_spectrum(&refined, params, spectrum, controls, cx).unwrap()
    });
    assert_eq!(finer.end_time_s, report.end_time_s);
    assert_eq!(
        finer.mlmc.levels[0].mean.to_bits(),
        report.mlmc.estimate.to_bits()
    );
    near(
        finer.mlmc.levels[0].variance / 12.0,
        report.mlmc.estimator_variance,
    );
}

#[test]
fn g4_spectral_fragility_refuses_aliasing_and_preserves_native_limits() {
    let ensemble = ensemble();
    let params = StoryParams::default();
    let original = ensemble.clone();
    let spectrum = grid();
    let controls = controls();
    let aliased = SpectralDiscretization {
        // Fits the fine clock, but not the declared 2dt comparison.
        omega_max_rad_s: 0.75 * core::f64::consts::PI / ensemble.dt.value,
        ..spectrum
    };
    assert!(matches!(
        with_cx(false, |cx| e_stopped_fragility_with_spectrum(
            &ensemble, params, aliased, controls, cx
        )),
        Err(FragilityError::CoarseBandAliased { .. })
    ));
    let no_lines = SpectralDiscretization {
        frequency_lines: 0,
        ..spectrum
    };
    assert!(matches!(
        with_cx(false, |cx| e_stopped_fragility_with_spectrum(
            &ensemble, params, no_lines, controls, cx
        )),
        Err(FragilityError::Realization { member: 0, .. })
    ));
    let no_work = FragilityControls {
        realization_budget: RealizationBudget {
            max_work: 0,
            ..controls.realization_budget
        },
        ..controls
    };
    assert!(matches!(
        with_cx(false, |cx| e_stopped_fragility_with_spectrum(
            &ensemble, params, spectrum, no_work, cx
        )),
        Err(FragilityError::Realization { member: 0, .. })
    ));
    let short_history = FragilityControls {
        history_limits: HistoryLimits {
            max_samples: 15,
            ..controls.history_limits
        },
        ..controls
    };
    assert!(matches!(
        with_cx(false, |cx| e_stopped_fragility_with_spectrum(
            &ensemble,
            params,
            spectrum,
            short_history,
            cx
        )),
        Err(FragilityError::History {
            member: 0,
            coarse: false,
            source: HistoryError::SampleLimitExceeded {
                samples: 16,
                max_samples: 15
            },
        })
    ));
    assert!(matches!(
        with_cx(false, |cx| e_stopped_fragility_with_spectrum(
            &ensemble,
            StoryParams {
                mass: f64::NAN,
                ..params
            },
            spectrum,
            controls,
            cx,
        )),
        Err(FragilityError::History {
            member: 0,
            source: HistoryError::InvalidStoryParameter { .. },
            ..
        })
    ));
    assert!(matches!(
        with_cx(true, |cx| e_stopped_fragility_with_spectrum(
            &ensemble, params, spectrum, controls, cx
        )),
        Err(FragilityError::Cancelled)
    ));
    assert_eq!(ensemble, original);
}

#[test]
fn g1_legacy_risk_uses_time_zero_without_decimating_its_spectrum() {
    let ensemble = ensemble();
    let params = StoryParams::default();
    let expected: Vec<f64> = (0..ensemble.members)
        .map(|member| {
            let record = ensemble.realize(member).unwrap();
            StoryFrame::new(params)
                .run_from_time_zero(
                    GroundMotion::new(&record.values, ensemble.dt.value),
                    controls().history_limits,
                )
                .unwrap()
                .peak_abs_displacement_m
                / params.h
        })
        .collect();
    let limit = expected.iter().sum::<f64>() / 12.0;
    let report = e_stopped_fragility(&ensemble, params, limit, 0.05, 0.0);
    assert_eq!(report.spectrum, None);
    assert_eq!(report.mlmc.levels.len(), 1);
    assert_eq!(report.mlmc.levels[0].samples, 12);
    assert_eq!(report.mlmc.total_cost, 12.0 * 15.0);
    assert_eq!(report.end_time_s, 15.0 * ensemble.dt.value);
    let exceedances = u32::try_from(expected.iter().filter(|&&peak| peak > limit).count()).unwrap();
    assert_eq!(report.exceedances, exceedances);
    near(report.mlmc.estimate, limit);
    assert_eq!(
        try_losses(&ensemble, params, params.scale).unwrap(),
        expected
    );
}
