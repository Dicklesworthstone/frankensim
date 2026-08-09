//! G0/G3/G4 robust-observation admission and disposition audit.

use fs_assimilate::{
    AssimError, BatchOutcome, Belief, CensorDirection, Observation, ObservationBatch,
    ObservationDisposition, RobustObservation, assimilate_all, assimilate_observation_batch,
    point_sensor,
};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey, VirtualClock};

const TEST_STREAM: StreamKey = StreamKey {
    seed: 0x524F_4255_5354,
    kernel_id: 0x524F,
    tile: 0,
    iteration: 0,
};

fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    let clock = VirtualClock::new();
    let result = pool.scope(|arena| {
        let cx = Cx::new(
            gate,
            arena,
            TEST_STREAM,
            Budget::INFINITE,
            ExecMode::Deterministic,
        )
        .with_time_source(&clock);
        f(&cx)
    });
    assert!(pool.stats().quiescent());
    result
}

fn sensor(component: usize, dim: usize, value: f64, noise: f64, id: &str) -> Observation {
    point_sensor(component, dim, value, noise, id).expect("valid sensor")
}

#[test]
fn diagonal_batch_replays_the_existing_core_exactly() {
    let gate = CancelGate::new();
    with_cx(&gate, |cx| {
        let prior = Belief::new(vec![3.0, -2.0], vec![vec![2.0, 0.25], vec![0.25, 1.5]], cx)
            .expect("valid prior");
        let observations = vec![
            sensor(0, 2, 2.5, 0.2, "thermocouple-a"),
            sensor(1, 2, -1.25, 0.4, "pressure-tap-b"),
        ];
        let expected = assimilate_all(&prior, &observations, cx).expect("legacy update");
        let records = observations
            .iter()
            .cloned()
            .map(RobustObservation::available)
            .collect();
        let batch = ObservationBatch::new(records, vec![vec![0.2, 0.0], vec![0.0, 0.4]], cx)
            .expect("checked diagonal batch");
        let result =
            assimilate_observation_batch(&prior, &batch, cx).expect("robust diagonal update");

        assert_eq!(result.posterior(), Some(&expected));
        assert_eq!(result.audit().outcome(), BatchOutcome::Updated);
        assert_eq!(result.audit().effective_dof(), 2);
        assert!(
            result
                .audit()
                .entries()
                .iter()
                .all(|(_, disposition)| *disposition == ObservationDisposition::Accepted)
        );
        assert!(
            result
                .audit()
                .receipt_id()
                .starts_with("robust-observation-audit:v1:")
        );
        let replay =
            assimilate_observation_batch(&prior, &batch, cx).expect("deterministic replay");
        assert_eq!(replay, result);
    });
}

#[test]
fn missing_record_is_excluded_and_reduces_effective_dof() {
    let gate = CancelGate::new();
    with_cx(&gate, |cx| {
        let prior = Belief::scalar(300.0, 4.0).expect("valid prior");
        let observation = sensor(0, 1, 302.0, 0.25, "thermocouple-live");
        let expected =
            assimilate_all(&prior, core::slice::from_ref(&observation), cx).expect("legacy update");
        let batch = ObservationBatch::new(
            vec![
                RobustObservation::available(observation),
                RobustObservation::missing("thermocouple-missing").expect("missing record"),
            ],
            vec![vec![0.25]],
            cx,
        )
        .expect("checked batch");
        let result = assimilate_observation_batch(&prior, &batch, cx).expect("missing exclusion");

        assert_eq!(result.posterior(), Some(&expected));
        assert_eq!(result.audit().effective_dof(), 1);
        assert_eq!(
            result.audit().entries(),
            &[
                (
                    "thermocouple-live".to_owned(),
                    ObservationDisposition::Accepted
                ),
                (
                    "thermocouple-missing".to_owned(),
                    ObservationDisposition::ExcludedMissing
                ),
            ]
        );
    });
}

#[test]
fn censored_saturated_and_delayed_records_refuse_without_partial_posterior() {
    let gate = CancelGate::new();
    with_cx(&gate, |cx| {
        let prior = Belief::scalar(0.0, 1.0).expect("valid prior");
        let ordinary = RobustObservation::available(sensor(0, 1, 0.5, 0.2, "ordinary"));
        let censored = RobustObservation::censored(
            sensor(0, 1, 1.0, 0.2, "censored"),
            1.0,
            CensorDirection::AtOrAbove,
        )
        .expect("censored record");
        let saturated =
            RobustObservation::saturated(sensor(0, 1, 5.0, 0.2, "saturated"), -5.0, 5.0)
                .expect("saturated record");
        let delayed = RobustObservation::delayed(sensor(0, 1, 0.4, 0.2, "delayed"), 0.5, 0.1)
            .expect("delayed record");
        let batch = ObservationBatch::new(
            vec![ordinary, censored, saturated, delayed],
            vec![vec![0.2]],
            cx,
        )
        .expect("checked refusal batch");
        let result =
            assimilate_observation_batch(&prior, &batch, cx).expect("audited refusal result");

        assert_eq!(result.posterior(), None);
        assert_eq!(result.audit().outcome(), BatchOutcome::RefusedPathology);
        assert_eq!(result.audit().effective_dof(), 0);
        let dispositions: Vec<_> = result
            .audit()
            .entries()
            .iter()
            .map(|(_, disposition)| *disposition)
            .collect();
        assert_eq!(
            dispositions,
            vec![
                ObservationDisposition::WithheldByBatchRefusal,
                ObservationDisposition::RefusedCensored,
                ObservationDisposition::RefusedSaturated,
                ObservationDisposition::RefusedDelayed,
            ]
        );
    });
}

#[test]
fn all_missing_batch_returns_an_audited_no_data_refusal() {
    let gate = CancelGate::new();
    with_cx(&gate, |cx| {
        let prior = Belief::scalar(0.0, 1.0).expect("valid prior");
        let batch = ObservationBatch::new(
            vec![RobustObservation::missing("missing-only").expect("missing record")],
            vec![],
            cx,
        )
        .expect("checked no-data batch");
        let result = assimilate_observation_batch(&prior, &batch, cx).expect("audited refusal");

        assert_eq!(result.posterior(), None);
        assert_eq!(
            result.audit().outcome(),
            BatchOutcome::RefusedNoUsableObservations
        );
        assert_eq!(
            result.audit().entries()[0].1,
            ObservationDisposition::ExcludedMissing
        );
    });
}

#[test]
fn correlated_covariance_changes_the_posterior_from_naive_diagonal_noise() {
    let gate = CancelGate::new();
    with_cx(&gate, |cx| {
        let prior = Belief::diagonal(vec![0.0, 0.0], &[1.0, 1.0], cx).expect("valid prior");
        let observations = vec![
            sensor(0, 2, 1.0, 1.0, "correlated-a"),
            sensor(1, 2, 2.0, 1.0, "correlated-b"),
        ];
        let naive = assimilate_all(&prior, &observations, cx).expect("naive diagonal update");
        let batch = ObservationBatch::new(
            observations
                .into_iter()
                .map(RobustObservation::available)
                .collect(),
            vec![vec![1.0, 0.75], vec![0.75, 1.0]],
            cx,
        )
        .expect("checked correlated batch");
        let result = assimilate_observation_batch(&prior, &batch, cx).expect("correlated update");
        let posterior = result.posterior().expect("published posterior");

        assert_ne!(posterior, &naive);
        assert!((posterior.mean()[0] - (0.5 / 3.4375)).abs() < 1.0e-12);
        assert!((posterior.mean()[1] - (3.25 / 3.4375)).abs() < 1.0e-12);
        assert!(posterior.covariance()[0][1] > 0.0);
        posterior.validate(cx).expect("posterior invariants");
    });
}

#[test]
fn covariance_gate_refuses_ambiguous_or_conflicting_noise_authority() {
    let gate = CancelGate::new();
    with_cx(&gate, |cx| {
        let indefinite = ObservationBatch::new(
            vec![
                RobustObservation::available(sensor(0, 2, 0.0, 1.0, "indefinite-a")),
                RobustObservation::available(sensor(1, 2, 0.0, 1.0, "indefinite-b")),
            ],
            vec![vec![1.0, 2.0], vec![2.0, 1.0]],
            cx,
        );
        assert_eq!(
            indefinite,
            Err(AssimError::ObservationCovarianceNotPositiveSemidefinite)
        );

        let singular = ObservationBatch::new(
            vec![
                RobustObservation::available(sensor(0, 2, 0.0, 1.0, "singular-a")),
                RobustObservation::available(sensor(1, 2, 0.0, 1.0, "singular-b")),
            ],
            vec![vec![1.0, 1.0], vec![1.0, 1.0]],
            cx,
        );
        assert_eq!(
            singular,
            Err(AssimError::ObservationCovarianceNotPositiveDefinite { pivot: 1 })
        );

        let mismatch = ObservationBatch::new(
            vec![RobustObservation::available(sensor(
                0, 1, 0.0, 0.25, "mismatch",
            ))],
            vec![vec![0.5]],
            cx,
        );
        assert_eq!(
            mismatch,
            Err(AssimError::ObservationCovarianceNoiseMismatch { index: 0 })
        );

        let duplicate = ObservationBatch::new(
            vec![
                RobustObservation::available(sensor(0, 1, 0.0, 1.0, "duplicate")),
                RobustObservation::missing("duplicate").expect("missing record"),
            ],
            vec![vec![1.0]],
            cx,
        );
        assert_eq!(
            duplicate,
            Err(AssimError::DuplicateObservationInstrument {
                instrument: "duplicate".to_owned()
            })
        );
    });
}

#[test]
fn invalid_pathology_metadata_and_precancel_are_refused() {
    let observation = sensor(0, 1, 0.0, 1.0, "pathology");
    assert_eq!(
        RobustObservation::saturated(observation.clone(), -1.0, 1.0),
        Err(AssimError::InvalidPathologyParameter {
            parameter: "saturated value endpoint"
        })
    );
    assert_eq!(
        RobustObservation::delayed(observation.clone(), 0.0, 1.0),
        Err(AssimError::InvalidPathologyParameter {
            parameter: "delay time constant"
        })
    );

    let gate = CancelGate::new();
    gate.request();
    with_cx(&gate, |cx| {
        let error = ObservationBatch::new(
            vec![RobustObservation::available(observation)],
            vec![vec![1.0]],
            cx,
        )
        .expect_err("pre-cancelled batch must refuse");
        assert!(matches!(
            error,
            AssimError::Cancelled {
                phase: "batch preflight",
                ..
            }
        ));
    });
}

// sj31i.7.3: the dimensional core's batch algebra applied to a two-channel
// robust batch — covariance diagonal dims are the squared reading dims per
// channel, mechanically.
#[test]
fn robust_batch_covariance_dims_follow_the_channel_algebra() {
    use fs_qty::Dims;
    use fs_qty::inference::{ObservationSchema, SlotSchema, StateSchema};
    use fs_qty::semantic::QuantitySpec;

    let length = Dims([1, 0, 0, 0, 0, 0]);
    let velocity = Dims([1, 0, -1, 0, 0, 0]);
    let gate = CancelGate::new();
    with_cx(&gate, |cx| {
        // Two channels with different reading dimensions: a length gauge
        // and a speedometer over a [length, velocity] state.
        let state = StateSchema::try_new(vec![
            SlotSchema::new(QuantitySpec::dimensional(length)),
            SlotSchema::new(QuantitySpec::dimensional(velocity)),
        ])
        .expect("state schema");
        let length_reading = ObservationSchema::new(QuantitySpec::dimensional(length));
        let velocity_reading = ObservationSchema::new(QuantitySpec::dimensional(velocity));

        let batch = ObservationBatch::new(
            vec![
                RobustObservation::available(sensor(0, 2, 10.5, 0.2, "gauge-length")),
                RobustObservation::available(sensor(1, 2, 1.25, 0.1, "gauge-velocity")),
            ],
            vec![vec![0.2, 0.0], vec![0.0, 0.1]],
            cx,
        )
        .expect("batch admits");
        assert_eq!(batch.records().len(), 2);

        // The batch covariance diagonal for channel k carries exactly the
        // squared reading dimensions of channel k, mechanically derived.
        assert_eq!(
            length_reading
                .noise_variance_dims()
                .expect("length variance dims"),
            Dims([2, 0, 0, 0, 0, 0])
        );
        assert_eq!(
            velocity_reading
                .noise_variance_dims()
                .expect("velocity variance dims"),
            Dims([2, 0, -2, 0, 0, 0])
        );
        assert_eq!(state.len(), 2);
        // Cross-channel covariance entries carry the product of the two
        // reading dimensions.
        let cross = length.checked_plus(velocity).expect("cross dims");
        assert_eq!(cross, Dims([2, 0, -1, 0, 0, 0]));
    });
}

// ---------------------------------------------------------------------------
// sj31i.16 shared-calibration covariance battery: common-mode floors,
// independent/correlated limits, duplicate-row replay refusal, and the
// naive-update exposure comparison.
// ---------------------------------------------------------------------------

mod shared_calibration_groups {
    use super::*;
    use fs_assimilate::groups::{
        GROUPED_RECEIPT_PREFIX, GroupedBatch, GroupedRecord, RowId, SharedSource,
        assimilate_grouped,
    };

    fn reading(row: &str, value: f64, noise: f64, source: Option<&str>) -> GroupedRecord {
        GroupedRecord::new(
            point_sensor(0, 1, value, noise, "shared-cal").expect("observation"),
            RowId::try_new(row).expect("row id"),
            source.map(str::to_string),
        )
        .expect("grouped record")
    }

    fn replicated_batch(n: usize, common_variance: f64, noise: f64) -> GroupedBatch {
        let records = (0..n)
            .map(|i| {
                reading(
                    &format!("row-{i}"),
                    1.0 + 0.001 * i as f64,
                    noise,
                    Some("cal"),
                )
            })
            .collect();
        GroupedBatch::try_new(
            records,
            vec![SharedSource::try_new("cal", common_variance).expect("source")],
        )
        .expect("batch")
    }

    #[test]
    fn common_mode_floor_bounds_replicated_readings() {
        let gate = CancelGate::new();
        with_cx(&gate, |cx| {
            let prior = Belief::scalar(0.0, 2.0).expect("prior");
            let batch = replicated_batch(8, 0.5, 0.1);
            let result = assimilate_grouped(&prior, &batch, cx).expect("grouped update");
            let grouped_variance = result
                .robust()
                .posterior()
                .expect("posterior")
                .variance(0)
                .expect("variance");
            let naive_variance = result
                .independent_posterior()
                .variance(0)
                .expect("naive variance");
            let floor = batch.common_mode_floor("cal", 2.0).expect("floor computes");

            // The exact grouped value matches the closed form, the naive
            // independent update dives below the common-mode floor (the
            // spurious N-fold claim this bead rejects), and the grouped
            // update respects it.
            assert!(
                (grouped_variance - floor).abs() <= 1e-9 * floor.max(1.0),
                "grouped {grouped_variance} vs exact floor {floor}"
            );
            assert!(
                naive_variance < floor,
                "naive independent claim {naive_variance} must dive below the floor {floor} (that is the defect)"
            );
            assert!(
                grouped_variance > naive_variance,
                "grouped posterior must be more conservative than the naive one"
            );
            assert!(result.identity().starts_with(GROUPED_RECEIPT_PREFIX));
        });
    }

    #[test]
    fn independent_and_fully_correlated_limits_are_exact() {
        let gate = CancelGate::new();
        with_cx(&gate, |cx| {
            let prior = Belief::scalar(0.0, 2.0).expect("prior");
            // Zero common variance: grouped path reduces to the independent
            // update bit-for-information.
            let independent_batch = replicated_batch(4, 0.0, 0.1);
            let grouped = assimilate_grouped(&prior, &independent_batch, cx).expect("grouped");
            let naive = assimilate_all(
                &prior,
                &(0..4_usize)
                    .map(|i| {
                        point_sensor(0, 1, 1.0 + 0.001 * i as f64, 0.1, format!("plain-{i}"))
                            .expect("obs")
                    })
                    .collect::<Vec<_>>(),
                cx,
            )
            .expect("naive");
            let grouped_var = grouped
                .robust()
                .posterior()
                .expect("p")
                .variance(0)
                .expect("v");
            let naive_var = naive.variance(0).expect("v");
            assert!(
                (grouped_var - naive_var).abs() <= 1e-12,
                "zero common variance must recover the independent limit: {grouped_var} vs {naive_var}"
            );

            // Very large common variance: the readings carry almost no
            // information about the state.
            let dominating = replicated_batch(4, 1.0e8, 0.1);
            let dominated = assimilate_grouped(&prior, &dominating, cx).expect("grouped");
            let dominated_var = dominated
                .robust()
                .posterior()
                .expect("p")
                .variance(0)
                .expect("v");
            assert!(
                dominated_var > 1.9,
                "dominating common mode must leave the prior nearly intact: {dominated_var}"
            );
        });
    }

    #[test]
    fn duplicate_rows_and_source_declaration_refusals_are_typed() {
        let duplicate = GroupedBatch::try_new(
            vec![
                reading("row-a", 1.0, 0.1, None),
                reading("row-a", 1.1, 0.1, None),
            ],
            Vec::new(),
        );
        assert!(matches!(
            duplicate,
            Err(AssimError::DuplicateDatasetRow { .. })
        ));

        let undeclared =
            GroupedBatch::try_new(vec![reading("row-a", 1.0, 0.1, Some("ghost"))], Vec::new());
        assert!(matches!(
            undeclared,
            Err(AssimError::SharedSourceDeclaration { .. })
        ));

        let unused = GroupedBatch::try_new(
            vec![reading("row-a", 1.0, 0.1, None)],
            vec![SharedSource::try_new("idle", 0.5).expect("source")],
        );
        assert!(matches!(
            unused,
            Err(AssimError::SharedSourceDeclaration { .. })
        ));
    }

    #[test]
    fn built_covariance_is_symmetric_psd_and_permutation_stable() {
        let gate = CancelGate::new();
        with_cx(&gate, |cx| {
            let records = vec![
                reading("row-a", 1.0, 0.10, Some("cal")),
                reading("row-b", 1.1, 0.20, Some("cal")),
                reading("row-c", 0.9, 0.05, None),
            ];
            let batch = GroupedBatch::try_new(
                records.clone(),
                vec![SharedSource::try_new("cal", 0.3).expect("source")],
            )
            .expect("batch");
            let covariance = batch.covariance();
            assert!((covariance[0][0] - 0.4).abs() <= 1e-15);
            assert!((covariance[0][1] - 0.3).abs() <= 1e-15);
            assert_eq!(covariance[0][1].to_bits(), covariance[1][0].to_bits());
            assert_eq!(covariance[0][2].to_bits(), 0.0_f64.to_bits());
            assert!((covariance[1][1] - 0.5).abs() <= 1e-15);
            assert!((covariance[2][2] - 0.05).abs() <= 1e-15);

            let prior = Belief::scalar(0.0, 2.0).expect("prior");
            let forward = assimilate_grouped(&prior, &batch, cx).expect("forward");
            let mut reversed_records = records;
            reversed_records.reverse();
            let reversed_batch = GroupedBatch::try_new(
                reversed_records,
                vec![SharedSource::try_new("cal", 0.3).expect("source")],
            )
            .expect("reversed batch");
            let reversed = assimilate_grouped(&prior, &reversed_batch, cx).expect("reversed");
            let forward_mean = forward.robust().posterior().expect("p").mean()[0];
            let reversed_mean = reversed.robust().posterior().expect("p").mean()[0];
            assert!(
                (forward_mean - reversed_mean).abs() <= 1e-12,
                "declaration order must not move the posterior: {forward_mean} vs {reversed_mean}"
            );
        });
    }

    #[test]
    fn replay_is_idempotent_and_cancellation_is_clean() {
        let gate = CancelGate::new();
        with_cx(&gate, |cx| {
            let prior = Belief::scalar(0.0, 2.0).expect("prior");
            let batch = replicated_batch(3, 0.25, 0.1);
            let first = assimilate_grouped(&prior, &batch, cx).expect("first");
            let second = assimilate_grouped(&prior, &batch, cx).expect("second");
            assert_eq!(first.identity(), second.identity());
        });
        let gate = CancelGate::new();
        gate.request();
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        let clock = VirtualClock::new();
        pool.scope(|arena| {
            let cx = Cx::new(
                &gate,
                arena,
                TEST_STREAM,
                Budget::INFINITE,
                ExecMode::Deterministic,
            )
            .with_time_source(&clock);
            let prior = Belief::scalar(0.0, 2.0).expect("prior");
            let batch = replicated_batch(3, 0.25, 0.1);
            let refusal = assimilate_grouped(&prior, &batch, &cx).expect_err("pre-cancel refuses");
            assert!(matches!(
                refusal,
                AssimError::Cancelled { .. } | AssimError::BudgetRefused(_)
            ));
        });
    }
}
