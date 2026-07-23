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
