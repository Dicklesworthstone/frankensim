//! Real tetrahedral Poisson states and coefficient adjoints, not algebraic stubs.
use fs_adjoint::DensityPoisson;
use fs_adjoint::ift::goal::{DensityGoalError, DensityGoalLimits};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};

fn limits() -> DensityGoalLimits {
    DensityGoalLimits { max_dofs: 256, max_cells: 4096,
        max_operator_applications: 2000, restart: 30, max_cycles: 20,
        primal_tolerance: 1e-12, adjoint_tolerance: 1e-12 }
}
fn relative(a: &[f64], b: &[f64]) -> f64 {
    let error: Vec<_> = a.iter().zip(b).map(|(x,y)| x-y).collect();
    fs_solver::norm2(&error) / fs_solver::norm2(b).max(f64::MIN_POSITIVE)
}

#[test]
fn uniform_coefficient_scaling_changes_the_actual_state_and_adjoint() {
    let (mesh, positions) = fs_feec::kuhn_cube(3);
    let model = DensityPoisson::new(&mesh, &positions, vec![1.0; mesh.tets.len()]);
    let load: Vec<_> = (0..model.n()).map(|i| 1.0 + i as f64 / 8.0).collect();
    let weights: Vec<_> = (0..model.n()).map(|i| 0.5 + (i%3) as f64 / 4.0).collect();
    let rho = vec![1.0; mesh.tets.len()];
    let a = model.solve_goal(&rho, &load, &weights, 0.0, limits(), None).unwrap();
    let b = model.solve_goal(&vec![2.0; rho.len()], &load, &weights, 0.25, limits(), None).unwrap();
    assert!(relative(&model.apply_density(&rho, &a.state), &load) < 1e-12);
    assert!(relative(&model.apply_density(&rho, &a.adjoint), &weights) < 1e-12);
    assert!(relative(&b.state, &a.state.iter().map(|v| v*0.5).collect::<Vec<_>>()) < 1e-11);
    assert!(relative(&b.adjoint, &a.adjoint.iter().map(|v| v*0.5).collect::<Vec<_>>()) < 1e-11);
    assert!((b.residual + 0.25 - 0.5*a.residual).abs() < 1e-10);
    assert!((a.gradient.iter().sum::<f64>() + a.residual).abs() < 1e-10);
    assert!(a.primal_relative_residual < limits().primal_tolerance);
    assert!(a.adjoint_relative_residual < limits().adjoint_tolerance);
    assert!(a.operator_applications <= limits().max_operator_applications);
    assert_eq!(a, model.solve_goal(&rho, &load, &weights, 0.0, limits(), None).unwrap());
}

#[test]
fn per_tetrahedron_adjoint_matches_independent_resolved_differences() {
    let (mesh, positions) = fs_feec::kuhn_cube(3);
    let rho: Vec<_> = (0..mesh.tets.len()).map(|i| 0.6 + (i%7) as f64 * 0.2).collect();
    let model = DensityPoisson::new(&mesh, &positions, rho.clone());
    let load: Vec<_> = (0..model.n()).map(|i| 1.0 + (i%3) as f64).collect();
    let weights: Vec<_> = (0..model.n()).map(|i| if i%2==0 {1.0} else {-0.3}).collect();
    let base = model.solve_goal(&rho, &load, &weights, 0.7, limits(), None).unwrap();
    let h = 1e-5;
    for i in 0..rho.len() {
        let mut plus = rho.clone(); plus[i] += h;
        let mut minus = rho.clone(); minus[i] -= h;
        let a = model.solve_goal(&plus, &load, &weights, 0.7, limits(), None).unwrap();
        let b = model.solve_goal(&minus, &load, &weights, 0.7, limits(), None).unwrap();
        let observed = (a.residual - b.residual)/(2.0*h);
        assert!((observed-base.gradient[i]).abs() < 2e-7 * (1.0+observed.abs()),
            "tet {i}: {observed} vs {}", base.gradient[i]);
    }
}

#[test]
fn zero_load_or_goal_has_the_exact_zero_solve_semantics() {
    let (mesh, positions) = fs_feec::kuhn_cube(2);
    let rho = vec![1.0; mesh.tets.len()];
    let model = DensityPoisson::new(&mesh, &positions, rho.clone());
    assert_eq!(model.n(), 1);
    let no_goal = model.solve_goal(&rho, &[1.0], &[0.0], 0.5, limits(), None).unwrap();
    assert_eq!(no_goal.residual, -0.5);
    assert!(no_goal.gradient.iter().all(|g| *g==0.0));
    assert_eq!(no_goal.adjoint_iterations, 0);
    let no_load = model.solve_goal(&rho, &[0.0], &[1.0], 0.5, limits(), None).unwrap();
    assert_eq!(no_load.state, vec![0.0]);
    assert_eq!(no_load.primal_iterations, 0);
    assert!(no_load.gradient.iter().all(|g| *g==0.0));
    let zero = model.solve_goal(&rho, &[0.0], &[0.0], 0.5,
        DensityGoalLimits { max_operator_applications: 0, ..limits() }, None).unwrap();
    assert_eq!(zero.operator_applications, 0);
}

#[test]
fn the_same_hard_operator_cap_funds_both_solves_and_fresh_residuals() {
    let (mesh, positions) = fs_feec::kuhn_cube(2);
    let rho = vec![1.0; mesh.tets.len()];
    let model = DensityPoisson::new(&mesh, &positions, rho.clone());
    for cap in 0..8 {
        let error = model.solve_goal(&rho, &[1.0], &[1.0], 0.0,
            DensityGoalLimits { max_operator_applications: cap, ..limits() }, None).unwrap_err();
        assert!(matches!(error, DensityGoalError::Budget { applications, limit, .. }
            if applications <= cap && limit == cap));
        if cap >= 4 {
            assert!(matches!(error, DensityGoalError::Budget { phase: "adjoint", applications: 4, .. }));
        }
    }
    let exact = model.solve_goal(&rho, &[1.0], &[1.0], 0.0,
        DensityGoalLimits { max_operator_applications: 8, ..limits() }, None).unwrap();
    assert_eq!(exact.operator_applications, 8);
}

#[test]
fn insufficient_krylov_work_never_returns_an_unchecked_gradient() {
    let (mesh, positions) = fs_feec::kuhn_cube(4);
    let rho: Vec<_> = (0..mesh.tets.len()).map(|i| 0.3+(i%11) as f64*0.2).collect();
    let model = DensityPoisson::new(&mesh, &positions, rho.clone());
    let load: Vec<_> = (0..model.n()).map(|i| 0.5+(i%5) as f64).collect();
    let weights = vec![1.0; model.n()];
    let error = model.solve_goal(&rho, &load, &weights, 0.0,
        DensityGoalLimits { restart: 1, max_cycles: 1, ..limits() }, None).unwrap_err();
    assert!(matches!(error, DensityGoalError::NotConverged { phase: "primal", applications: 4, .. }));
    let funded = model.solve_goal(&rho, &load, &weights, 0.0, limits(), None).unwrap();
    assert!(funded.primal_relative_residual < 1e-12);
}

#[test]
fn invalid_inputs_and_cancelled_context_do_not_publish_partial_samples() {
    let (mesh, positions) = fs_feec::kuhn_cube(2);
    let rho = vec![1.0; mesh.tets.len()];
    let model = DensityPoisson::new(&mesh, &positions, rho.clone());
    assert!(model.solve_goal(&rho[..rho.len()-1], &[1.0], &[1.0], 0.0, limits(), None).is_err());
    for value in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        let mut invalid = rho.clone(); invalid[0] = value;
        assert!(model.solve_goal(&invalid, &[1.0], &[1.0], 0.0, limits(), None).is_err());
    }
    assert!(model.solve_goal(&rho, &[1.0], &[1.0], 0.0,
        DensityGoalLimits { max_dofs: 0, ..limits() }, None).is_err());
    let gate = CancelGate::new(); gate.request();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(&gate, arena, StreamKey {seed:0,kernel_id:1,tile:0,iteration:0},
            Budget::INFINITE, ExecMode::Deterministic);
        assert_eq!(model.solve_goal(&rho, &[1.0], &[1.0], 0.0, limits(), Some(&cx)).unwrap_err(),
            DensityGoalError::Cancelled { applications: 0 });
    });
}
