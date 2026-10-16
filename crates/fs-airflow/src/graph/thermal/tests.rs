use super::*;
use crate::conjugate::{AirSegment, solve_conjugate};
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};

fn with_gate<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(gate, arena, StreamKey {
            seed: 73, kernel_id: 17, tile: 0, iteration: 0,
        }, Budget::INFINITE, ExecMode::Deterministic))
    })
}

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    with_gate(&CancelGate::new(), f)
}

fn path(name: &str, inlet: f64, capacity: f64, conductance: f64) -> AirPath {
    AirPath::new(inlet, capacity / 1000.0, 1000.0,
        vec![AirSegment::new(name, 1.0, conductance).expect("segment")])
        .expect("air path")
}

fn response(paths: &[AirPath], references: &[f64], walls: &[f64]) -> Vec<SolidRegionState> {
    paths.iter().flat_map(|p| p.segments()).zip(references).zip(walls)
        .map(|((segment, &reference), &wall)| SolidRegionState {
            region: segment.region().to_string(),
            area_m2: segment.area_m2(),
            mean_wall_temperature_k: wall,
            heat_rate_w: segment.area_m2() * segment.htc_w_per_m2_k() * (wall - reference),
            mean_reference_temperature_k: Some(reference),
        }).collect()
}

fn coupled_paths() -> Vec<AirPath> {
    vec![path("left", 295.0, 2.0, 1.0), path("right", 310.0, 5.0, 3.0)]
}

fn shared_solid(paths: &[AirPath], references: &[f64]) -> Vec<SolidRegionState> {
    let wall = (10.0 + references[0] + 3.0 * references[1]) / 4.0;
    response(paths, references, &[wall, wall])
}

fn close(actual: f64, expected: f64) {
    assert!((actual - expected).abs() < 2.0e-8,
        "actual {actual:.16e}, expected {expected:.16e}");
}

#[test]
fn independent_inlets_and_streamwise_states_match_closed_form() {
    let first = AirPath::new(290.0, 0.002, 1000.0, vec![
        AirSegment::new("upstream", 1.0, 1.0).unwrap(),
        AirSegment::new("downstream", 1.0, 2.0).unwrap(),
    ]).unwrap();
    let paths = vec![first, path("separate", 310.0, 4.0, 2.0)];
    let result = with_cx(|cx| solve_conjugate_branches(cx, &paths,
        &ConjugateConfig::default(), |_, refs| Ok(response(&paths, refs, &[350.0; 3])))).unwrap();
    close(result.branches[0].march.outlet_temperature_k, 350.0 - 60.0 * (-1.5_f64).exp());
    close(result.branches[1].march.outlet_temperature_k, 350.0 - 40.0 * (-0.5_f64).exp());
    assert_eq!(result.branches[1].march.segments[0].inlet_temperature_k, 310.0);
    assert_eq!(result.reference_temperatures_k.len(), 3);
}

#[test]
fn shared_solid_heat_split_matches_independent_energy_solution() {
    let paths = coupled_paths();
    let mut calls = 0;
    let result = with_cx(|cx| solve_conjugate_branches(cx, &paths,
        &ConjugateConfig::default(), |_, refs| {
            calls += 1;
            Ok(shared_solid(&paths, refs))
        })).unwrap();
    let k0 = 2.0 * (1.0 - (-0.5_f64).exp());
    let k1 = 5.0 * (1.0 - (-0.6_f64).exp());
    let wall = (10.0 + k0 * 295.0 + k1 * 310.0) / (k0 + k1);
    close(result.branches[0].solid[0].mean_wall_temperature_k, wall);
    close(result.branches[1].solid[0].mean_wall_temperature_k, wall);
    close(result.branches[0].march.total_heat_rate_w, k0 * (wall - 295.0));
    close(result.branches[1].march.total_heat_rate_w, k1 * (wall - 310.0));
    close(result.branches.iter().map(|b| b.march.total_heat_rate_w).sum(), 10.0);
    assert_eq!(calls, result.iterations);
    for branch in &result.branches {
        assert_eq!(branch.iterations, calls);
        assert_eq!(branch.history.len(), calls);
        assert_eq!(branch.history.last().unwrap().iteration, calls - 1);
    }
}

#[test]
fn two_conductively_connected_solids_match_reduced_two_by_two_system() {
    let paths = coupled_paths();
    let coupling = 2.0;
    let result = with_cx(|cx| solve_conjugate_branches(cx, &paths,
        &ConjugateConfig { max_iterations: 200, ..ConjugateConfig::default() }, |_, refs| {
            let a = 1.0 + coupling;
            let d = 3.0 + coupling;
            let b0 = 10.0 + refs[0];
            let b1 = 4.0 + 3.0 * refs[1];
            let determinant = a * d - coupling * coupling;
            let walls = [(d * b0 + coupling * b1) / determinant,
                         (coupling * b0 + a * b1) / determinant];
            Ok(response(&paths, refs, &walls))
        })).unwrap();
    let k0 = 2.0 * (1.0 - (-0.5_f64).exp());
    let k1 = 5.0 * (1.0 - (-0.6_f64).exp());
    let a = k0 + coupling;
    let d = k1 + coupling;
    let b0 = 10.0 + k0 * 295.0;
    let b1 = 4.0 + k1 * 310.0;
    let determinant = a * d - coupling * coupling;
    close(result.branches[0].solid[0].mean_wall_temperature_k,
        (d * b0 + coupling * b1) / determinant);
    close(result.branches[1].solid[0].mean_wall_temperature_k,
        (coupling * b0 + a * b1) / determinant);
    close(result.branches.iter().map(|b| b.march.total_heat_rate_w).sum(), 14.0);
}

#[test]
fn single_branch_is_exactly_the_existing_driver() {
    let paths = vec![path("one", 300.0, 2.0, 1.0)];
    let config = ConjugateConfig::default();
    let multi = with_cx(|cx| solve_conjugate_branches(cx, &paths, &config,
        |_, refs| Ok(response(&paths, refs, &[350.0])))).unwrap();
    let single = with_cx(|cx| solve_conjugate(cx, &paths[0], &config,
        |_, refs| Ok(response(&paths, refs, &[350.0])))).unwrap();
    assert_eq!(multi.branches, vec![single]);
}

#[test]
fn cancellation_retains_all_branches_and_fixed_relaxation_resumes_exactly() {
    let paths = coupled_paths();
    let config = ConjugateConfig::default();
    let gate = CancelGate::new();
    let mut calls = 0;
    let error = with_gate(&gate, |cx| solve_conjugate_branches(cx, &paths, &config, |_, refs| {
        calls += 1;
        if calls == 3 { gate.request(); }
        Ok(shared_solid(&paths, refs))
    })).unwrap_err();
    let AirflowError::Cancelled { iteration, references_k } = error else {
        panic!("expected cancellation");
    };
    assert_eq!(iteration, 2);
    assert_eq!(references_k.len(), 2);
    let resumed = with_cx(|cx| solve_conjugate_branches_from(cx, &paths, &config, &references_k,
        |_, refs| Ok(shared_solid(&paths, refs)))).unwrap();
    let full = with_cx(|cx| solve_conjugate_branches(cx, &paths, &config,
        |_, refs| Ok(shared_solid(&paths, refs)))).unwrap();
    assert_eq!(resumed.reference_temperatures_k, full.reference_temperatures_k);
    for (a, b) in resumed.branches.iter().zip(&full.branches) {
        assert_eq!(a.march, b.march);
        assert_eq!(a.solid, b.solid);
    }
}

#[test]
fn precancelled_invocation_does_not_run_the_solid() {
    let gate = CancelGate::new();
    gate.request();
    let paths = coupled_paths();
    let error = with_gate(&gate, |cx| solve_conjugate_branches(cx, &paths,
        &ConjugateConfig::default(), |_, _| panic!("cancelled before solid"))).unwrap_err();
    assert!(matches!(error, AirflowError::Cancelled { iteration: 0, references_k } if references_k.len() == 2));
}

#[test]
fn surfaces_cannot_belong_to_two_independent_branches() {
    let paths = vec![path("same", 300.0, 1.0, 1.0), path("same", 310.0, 2.0, 1.0)];
    let error = with_cx(|cx| solve_conjugate_branches(cx, &paths,
        &ConjugateConfig::default(), |_, _| panic!("invalid ownership"))).unwrap_err();
    assert!(matches!(error, AirflowError::DuplicateAirSegment { region } if region == "same"));
}

#[test]
fn shared_iteration_budget_is_not_multiplied_by_branch_count() {
    let paths = coupled_paths();
    let mut calls = 0;
    let config = ConjugateConfig { max_iterations: 1, ..ConjugateConfig::default() };
    let error = with_cx(|cx| solve_conjugate_branches(cx, &paths, &config, |_, refs| {
        calls += 1;
        Ok(shared_solid(&paths, refs))
    })).unwrap_err();
    assert_eq!(calls, 1);
    assert!(matches!(error, AirflowError::ConjugateNotConverged { iterations: 1, .. }));
}

#[test]
fn large_branch_cannot_mask_small_branch_heat_balance_fault() {
    let paths = vec![path("small", 300.0, 1.0, 1.0), path("large", 300.0, 1.0, 1.0)];
    let error = with_cx(|cx| solve_conjugate_branches(cx, &paths,
        &ConjugateConfig::default(), |_, refs| {
            let mut states = response(&paths, refs, &[300.0001, 1000.0]);
            states[0].heat_rate_w += 1.0e-5;
            Ok(states)
        })).unwrap_err();
    assert!(matches!(error, AirflowError::ConjugateBalanceUnclosed { .. }));
}

#[test]
fn malformed_shared_response_uses_existing_wiring_refusals() {
    let paths = coupled_paths();
    let config = ConjugateConfig::default();
    let error = with_cx(|cx| solve_conjugate_branches(cx, &paths, &config, |_, refs| {
        let mut states = shared_solid(&paths, refs);
        states.swap(0, 1);
        Ok(states)
    })).unwrap_err();
    assert!(matches!(error, AirflowError::SegmentRegionMismatch { .. }));
    let error = with_cx(|cx| solve_conjugate_branches(cx, &paths, &config, |_, _| Ok(Vec::new()))).unwrap_err();
    assert!(matches!(error, AirflowError::SolidResponseArity { expected: 2, found: 0 }));
}

#[test]
fn invalid_configuration_and_resume_shape_refuse_before_solid_work() {
    let paths = coupled_paths();
    for config in [
        ConjugateConfig { max_iterations: 0, ..ConjugateConfig::default() },
        ConjugateConfig { temperature_tolerance_k: f64::NAN, ..ConjugateConfig::default() },
        ConjugateConfig { balance_relative_tolerance: 1.0, ..ConjugateConfig::default() },
        ConjugateConfig { relaxation: Relaxation::Fixed { omega: 0.0 }, ..ConjugateConfig::default() },
    ] {
        assert!(with_cx(|cx| solve_conjugate_branches(cx, &paths, &config,
            |_, _| panic!("invalid config"))).is_err());
    }
    assert!(with_cx(|cx| solve_conjugate_branches_from(cx, &paths,
        &ConjugateConfig::default(), &[300.0], |_, _| panic!("invalid resume"))).is_err());
    assert!(with_cx(|cx| solve_conjugate_branches(cx, &[],
        &ConjugateConfig::default(), |_, _| panic!("no branches"))).is_err());
}

#[test]
fn branch_aitken_closes_the_same_physical_fixed_point() {
    let paths = coupled_paths();
    let config = ConjugateConfig {
        relaxation: Relaxation::Aitken { omega_init: 0.5, omega_max: 1.0 },
        max_iterations: 200,
        ..ConjugateConfig::default()
    };
    let result = with_cx(|cx| solve_conjugate_branches(cx, &paths, &config,
        |_, refs| Ok(shared_solid(&paths, refs)))).unwrap();
    let plain = with_cx(|cx| solve_conjugate_branches(cx, &paths, &ConjugateConfig::default(),
        |_, refs| Ok(shared_solid(&paths, refs)))).unwrap();
    for (a, b) in result.branches.iter().zip(&plain.branches) {
        close(a.march.outlet_temperature_k, b.march.outlet_temperature_k);
    }
}
