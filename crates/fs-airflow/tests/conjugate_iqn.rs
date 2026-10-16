//! Opt-in acceleration against a closed-form, two-segment lumped-solid case.
//! These are numerical model tests, not a whole cooling-product validation.
use fs_airflow::AirflowError;
use fs_airflow::conjugate::{
    AirPath, AirSegment, ConjugateConfig, IqnIlsConfig, SolidRegionState,
    solve_conjugate, solve_conjugate_iqn, solve_conjugate_iqn_from,
};
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 17,
                kernel_id: 71,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

fn path() -> AirPath {
    AirPath::new(
        300.0,
        0.001,
        1_000.0,
        vec![
            AirSegment::new("first", 1.0, 50.0).unwrap(),
            AirSegment::new("second", 1.0, 25.0).unwrap(),
        ],
    )
    .unwrap()
}

fn config() -> ConjugateConfig {
    ConjugateConfig {
        max_iterations: 12,
        ..ConjugateConfig::default()
    }
}

fn lumped_solid(_: &Cx<'_>, references: &[f64]) -> Result<Vec<SolidRegionState>, AirflowError> {
    // Each solid obeys its actual steady balance P=hA(T_wall-T_ref).
    Ok([("first", 50.0, 5.0), ("second", 25.0, 3.0)]
        .into_iter()
        .zip(references)
        .map(|((region, conductance, power), &reference)| SolidRegionState {
            region: region.to_owned(),
            area_m2: 1.0,
            mean_wall_temperature_k: reference + power / conductance,
            heat_rate_w: power,
            mean_reference_temperature_k: Some(reference),
        })
        .collect())
}

#[test]
fn iqn_closes_high_ntu_coupling_within_the_declared_budget() {
    with_cx(|cx| {
        let path = path();
        let config = config();
        assert!(matches!(
            solve_conjugate(cx, &path, &config, lumped_solid),
            Err(AirflowError::ConjugateNotConverged { .. })
        ));
        let solution = solve_conjugate_iqn(
            cx, &path, &config, IqnIlsConfig::default(), lumped_solid,
        )
        .unwrap();
        // Independent global balance and closed-form segment effectiveness.
        // Capacity rate is one W/K, so the two air rises are 5 K and 3 K.
        let expected_walls = [
            300.0 + 5.0 / (1.0 - (-50.0_f64).exp()),
            305.0 + 3.0 / (1.0 - (-25.0_f64).exp()),
        ];
        for (state, expected) in solution.solid.iter().zip(expected_walls) {
            assert!((state.mean_wall_temperature_k - expected).abs() < 1.0e-7);
        }
        assert!((solution.march.outlet_temperature_k - 308.0).abs() < 1.0e-7);
        assert!((solution.balance.air_total_w - 8.0).abs() < 1.0e-7);
        assert!(solution.balance.max_region_imbalance_w < 1.0e-7);
        assert!(solution.iterations <= config.max_iterations);
    });
}

#[test]
fn vector_acceleration_does_not_bypass_the_watt_balance_gate() {
    with_cx(|cx| {
        let result = solve_conjugate_iqn(
            cx,
            &path(),
            &config(),
            IqnIlsConfig::default(),
            |cx, references| {
                let mut states = lumped_solid(cx, references)?;
                states[0].heat_rate_w += 1.0;
                Ok(states)
            },
        );
        assert!(matches!(
            result,
            Err(AirflowError::ConjugateBalanceUnclosed { .. })
        ));
    });
}

#[test]
fn invalid_history_policy_refuses_before_invoking_the_solid() {
    with_cx(|cx| {
        let mut called = false;
        let result = solve_conjugate_iqn(
            cx,
            &path(),
            &config(),
            IqnIlsConfig {
                max_history: 0,
                ..IqnIlsConfig::default()
            },
            |cx, references| {
                called = true;
                lumped_solid(cx, references)
            },
        );
        assert!(!called);
        assert!(matches!(
            result,
            Err(AirflowError::InvalidConjugateInput {
                field: "IQN history limit",
                ..
            })
        ));
    });
}

#[test]
fn a_reference_only_restart_uses_fresh_history_and_the_same_physics() {
    with_cx(|cx| {
        let solution = solve_conjugate_iqn_from(
            cx,
            &path(),
            &config(),
            &[302.0, 304.0],
            IqnIlsConfig::default(),
            lumped_solid,
        )
        .unwrap();
        assert!((solution.march.outlet_temperature_k - 308.0).abs() < 1.0e-7);
        assert!(solution.balance.max_region_imbalance_w < 1.0e-7);
    });
}
