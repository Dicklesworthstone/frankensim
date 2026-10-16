//! Two air supplies coupled through a conducting two-node solid. The oracle
//! eliminates the Robin references and solves the physical heat balance,
//! independently of the partitioned iteration and its accelerator.
use fs_airflow::AirflowError;
use fs_airflow::conjugate::{
    AirPath, AirSegment, ConjugateConfig, IqnIlsConfig, SolidRegionState,
    solve_conjugate_iqn,
};
use fs_airflow::graph::thermal::{
    solve_conjugate_branches, solve_conjugate_branches_iqn,
    solve_conjugate_branches_iqn_from,
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
            StreamKey { seed: 17, kernel_id: 72, tile: 0, iteration: 0 },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

fn paths() -> Vec<AirPath> {
    vec![
        AirPath::new(300.0, 0.001, 1_000.0,
            vec![AirSegment::new("cold", 1.0, 50.0).unwrap()]).unwrap(),
        AirPath::new(310.0, 0.002, 1_000.0,
            vec![AirSegment::new("warm", 1.0, 25.0).unwrap()]).unwrap(),
    ]
}

fn config() -> ConjugateConfig {
    ConjugateConfig { max_iterations: 12, ..ConjugateConfig::default() }
}

fn shared_solid(_: &Cx<'_>, refs: &[f64]) -> Result<Vec<SolidRegionState>, AirflowError> {
    // Conductance of 1 W/K between the two nodes, sources 5 and 3 W:
    // [51 -1; -1 26] walls = [5 + 50*ref0; 3 + 25*ref1].
    let b0 = 5.0 + 50.0 * refs[0];
    let b1 = 3.0 + 25.0 * refs[1];
    let determinant = 51.0 * 26.0 - 1.0;
    let walls = [(26.0 * b0 + b1) / determinant, (b0 + 51.0 * b1) / determinant];
    Ok([("cold", 50.0), ("warm", 25.0)]
        .into_iter().enumerate()
        .map(|(index, (region, conductance))| SolidRegionState {
            region: region.to_owned(),
            area_m2: 1.0,
            mean_wall_temperature_k: walls[index],
            heat_rate_w: conductance * (walls[index] - refs[index]),
            mean_reference_temperature_k: Some(refs[index]),
        }).collect())
}

fn physical_solution() -> [f64; 2] {
    // Eliminate the fluid analytically: Q_i = C_i*(1-exp(-NTU_i))*(wall-inlet).
    let g0 = 1.0 - (-50.0_f64).exp();
    let g1 = 2.0 * (1.0 - (-12.5_f64).exp());
    let b0 = 5.0 + g0 * 300.0;
    let b1 = 3.0 + g1 * 310.0;
    let determinant = (g0 + 1.0) * (g1 + 1.0) - 1.0;
    [((g1 + 1.0) * b0 + b1) / determinant,
     (b0 + (g0 + 1.0) * b1) / determinant]
}

#[test]
fn global_secants_resolve_common_solid_coupling_with_one_solve_per_iteration() {
    with_cx(|cx| {
        let paths = paths();
        assert!(matches!(
            solve_conjugate_branches(cx, &paths, &config(), shared_solid),
            Err(AirflowError::ConjugateNotConverged { .. })
        ));
        let mut solves = 0;
        let solution = solve_conjugate_branches_iqn(
            cx, &paths, &config(), IqnIlsConfig::default(), |cx, refs| {
                solves += 1;
                shared_solid(cx, refs)
            },
        ).unwrap();
        assert_eq!(solves, solution.iterations);
        assert!(solves <= config().max_iterations);
        for (branch, expected) in solution.branches.iter().zip(physical_solution()) {
            assert!((branch.solid[0].mean_wall_temperature_k - expected).abs() < 1.0e-7);
            assert!(branch.balance.max_region_imbalance_w < 1.0e-7);
            assert_eq!(branch.history.len(), solves);
        }
        let heat: f64 = solution.branches.iter().map(|b| b.balance.air_total_w).sum();
        assert!((heat - 8.0).abs() < 1.0e-7);
    });
}

#[test]
fn cross_branch_acceleration_keeps_the_smaller_branch_watt_gate() {
    with_cx(|cx| {
        let result = solve_conjugate_branches_iqn(
            cx, &paths(), &config(), IqnIlsConfig::default(), |cx, refs| {
                let mut states = shared_solid(cx, refs)?;
                states[1].heat_rate_w += 1.0e-4;
                Ok(states)
            },
        );
        assert!(matches!(result, Err(AirflowError::ConjugateBalanceUnclosed { .. })));
    });
}

#[test]
fn warm_start_preserves_the_same_cross_branch_physics() {
    with_cx(|cx| {
        let solution = solve_conjugate_branches_iqn_from(
            cx, &paths(), &config(), &[310.0, 300.0],
            IqnIlsConfig::default(), shared_solid,
        ).unwrap();
        for (branch, expected) in solution.branches.iter().zip(physical_solution()) {
            assert!((branch.solid[0].mean_wall_temperature_k - expected).abs() < 1.0e-7);
        }
    });
}

#[test]
fn one_branch_delegates_to_the_single_path_iqn_driver() {
    fn solid(_: &Cx<'_>, refs: &[f64]) -> Result<Vec<SolidRegionState>, AirflowError> {
        Ok(vec![SolidRegionState {
            region: "cold".to_owned(),
            area_m2: 1.0,
            mean_wall_temperature_k: refs[0] + 0.1,
            heat_rate_w: 5.0,
            mean_reference_temperature_k: Some(refs[0]),
        }])
    }
    with_cx(|cx| {
        let paths = paths();
        let single = solve_conjugate_iqn(
            cx, &paths[0], &config(), IqnIlsConfig::default(), solid,
        ).unwrap();
        let wrapped = solve_conjugate_branches_iqn(
            cx, &paths[..1], &config(), IqnIlsConfig::default(), solid,
        ).unwrap();
        assert_eq!(wrapped.branches, vec![single]);
    });
}

#[test]
fn malformed_acceleration_policy_never_runs_the_common_solid() {
    with_cx(|cx| {
        let mut called = false;
        let result = solve_conjugate_branches_iqn(
            cx, &paths(), &config(),
            IqnIlsConfig { max_history: 0, ..IqnIlsConfig::default() },
            |cx, refs| { called = true; shared_solid(cx, refs) },
        );
        assert!(!called);
        assert!(matches!(result, Err(AirflowError::InvalidConjugateInput {
            field: "IQN history limit", ..
        })));
    });
}
