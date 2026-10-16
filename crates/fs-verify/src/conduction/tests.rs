use fs_alloc::{ArenaConfig, ArenaPool};
use fs_conduction::fixtures::{on_box_face, unit_cube};
use fs_conduction::verification::{
    ConductionBoundError, FluxBudget, MeanSolveConfig, MeanTemperatureSolution, TetError,
    solve_with_mean_bound,
};
use fs_conduction::{
    ConductionMesh, ConductionProblem, ConductivityModel, ElementMaterials, InitialGuess,
    MaterialId, MaterialTable, Nonlinearity, ScalarField, SolveConfig, ThermalBc, ThermalBoundary,
    ThermalBoundaryBuilder,
};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};

fn with_cx<T>(f: impl FnOnce(&Cx<'_>) -> T) -> T {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 71,
                kernel_id: 11,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}
fn mesh(n: usize) -> ConductionMesh {
    let (complex, positions) = unit_cube(n);
    ConductionMesh::new(complex, positions).unwrap()
}
fn ends(mesh: &ConductionMesh, left: f64, right: f64) -> ThermalBoundary {
    ThermalBoundaryBuilder::new(mesh)
        .region(
            "left",
            |f| on_box_face(f.centroid[0], 0.0),
            ThermalBc::dirichlet(left).unwrap(),
        )
        .unwrap()
        .region(
            "right",
            |f| on_box_face(f.centroid[0], 1.0),
            ThermalBc::dirichlet(right).unwrap(),
        )
        .unwrap()
        .adiabatic_remainder()
        .finish()
        .unwrap()
}
fn config() -> MeanSolveConfig {
    let mut solve = SolveConfig::default();
    solve.nonlinearity = Nonlinearity::FixedPoint {
        relaxation: 1.0,
        max_backtracks: 4,
    };
    solve.stop.residual_rtol = 1e-11;
    solve.stop.residual_atol = 1e-12;
    solve.linear.tolerance = 1e-13;
    let mut dual = solve.clone();
    dual.initial = InitialGuess::Uniform(0.0);
    MeanSolveConfig {
        primal: solve,
        dual,
        flux: FluxBudget::default(),
    }
}
fn contains(solution: &MeanTemperatureSolution, truth: f64) {
    let interval = solution.bound.enclosure;
    assert!(
        interval.lo <= truth && truth <= interval.hi,
        "{interval:?} excludes {truth}"
    );
}

#[test]
fn actual_primal_and_dual_enclose_mean_with_a_refinement_ladder() {
    with_cx(|cx| {
        let mut previous = f64::INFINITY;
        for n in [2, 4] {
            let mesh = mesh(n);
            let boundary = ends(&mesh, 300.0, 300.0);
            let material = ConductivityModel::isotropic_declared(1.0).unwrap();
            let source = ScalarField::Uniform(2.0);
            let result = solve_with_mean_bound(
                cx,
                ConductionProblem {
                    mesh: &mesh,
                    boundary: &boundary,
                    material: &material,
                    element_materials: None,
                    source: &source,
                },
                config(),
            )
            .unwrap();
            contains(&result, 300.0 + 1.0 / 6.0);
            let width = result.bound.enclosure.hi - result.bound.enclosure.lo;
            assert!(
                width < 0.5 * previous,
                "real-solve mean width n={n}: {width} vs {previous}"
            );
            previous = width;
            assert_eq!(result.primal.temperature.len(), mesh.vertex_count());
            assert_eq!(result.dual.temperature.len(), mesh.vertex_count());
            assert!(result.primal.temperature.iter().all(|t| *t >= 300.0));
            assert!(result.dual.temperature.iter().any(|z| *z > 0.0));
            for &(vertex, _) in boundary.dirichlet() {
                assert_eq!(result.dual.temperature[vertex], 0.0);
            }
        }
    });
}

#[test]
fn load_changes_reach_the_actual_field_and_its_bound() {
    with_cx(|cx| {
        let mesh = mesh(4);
        let boundary = ends(&mesh, 300.0, 300.0);
        let material = ConductivityModel::isotropic_declared(1.0).unwrap();
        let mut means = Vec::new();
        for f in [0.0, 2.0, 4.0] {
            let source = ScalarField::Uniform(f);
            let result = solve_with_mean_bound(
                cx,
                ConductionProblem {
                    mesh: &mesh,
                    boundary: &boundary,
                    material: &material,
                    element_materials: None,
                    source: &source,
                },
                config(),
            )
            .unwrap();
            contains(&result, 300.0 + f / 12.0);
            means.push(0.5 * result.bound.candidate_mean.lo + 0.5 * result.bound.candidate_mean.hi);
        }
        assert!(means[1] > means[0] + 0.1);
        assert!((means[2] - means[0] - 2.0 * (means[1] - means[0])).abs() < 1e-8);
    });
}

#[test]
fn heterogeneous_series_slab_uses_assigned_materials_not_the_fallback() {
    with_cx(|cx| {
        let mesh = mesh(2);
        let boundary = ends(&mesh, 300.0, 301.0);
        let fallback = ConductivityModel::isotropic_declared(99.0).unwrap();
        let table = MaterialTable::new([
            (
                MaterialId(1),
                ConductivityModel::isotropic_declared(1.0).unwrap(),
            ),
            (
                MaterialId(2),
                ConductivityModel::isotropic_declared(4.0).unwrap(),
            ),
        ])
        .unwrap();
        let ids = mesh
            .complex()
            .tets
            .iter()
            .map(|tet| {
                let x = tet
                    .iter()
                    .map(|&v| mesh.positions()[v as usize][0])
                    .sum::<f64>()
                    / 4.0;
                if x < 0.5 {
                    MaterialId(1)
                } else {
                    MaterialId(2)
                }
            })
            .collect();
        let assigned = ElementMaterials::new(table, ids).unwrap();
        let source = ScalarField::Uniform(0.0);
        let result = solve_with_mean_bound(
            cx,
            ConductionProblem {
                mesh: &mesh,
                boundary: &boundary,
                material: &fallback,
                element_materials: Some(&assigned),
                source: &source,
            },
            config(),
        )
        .unwrap();
        // R_left=.5, R_right=.125: interface=300.8; piecewise-linear mean=300.65.
        contains(&result, 300.65);
        assert!(result.bound.enclosure.hi - result.bound.enclosure.lo < 1e-6);
        for (p, &t) in mesh.positions().iter().zip(&result.primal.temperature) {
            if p[0] == 0.5 {
                assert!((t - 300.8).abs() < 1e-8);
            }
        }
    });
}

#[test]
fn robin_only_problem_preserves_reference_and_dual_trace_terms() {
    with_cx(|cx| {
        let mesh = mesh(2);
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .region("ambient", |_| true, ThermalBc::robin(2.0, 300.0).unwrap())
            .unwrap()
            .finish()
            .unwrap();
        let material = ConductivityModel::isotropic_declared(1.0).unwrap();
        let source = ScalarField::Uniform(0.0);
        let result = solve_with_mean_bound(
            cx,
            ConductionProblem {
                mesh: &mesh,
                boundary: &boundary,
                material: &material,
                element_materials: None,
                source: &source,
            },
            config(),
        )
        .unwrap();
        contains(&result, 300.0);
        assert!(result.bound.enclosure.hi - result.bound.enclosure.lo < 1e-6);
        assert!(
            result
                .primal
                .temperature
                .iter()
                .all(|t| (*t - 300.0).abs() < 1e-8)
        );
        assert!(result.dual.temperature.iter().all(|z| *z > 0.0));
    });
}

#[test]
fn nodal_source_is_bound_and_work_budget_still_refuses_before_solving() {
    with_cx(|cx| {
        let mesh = mesh(2);
        let boundary = ends(&mesh, 300.0, 300.0);
        let material = ConductivityModel::isotropic_declared(1.0).unwrap();
        let source = ScalarField::Nodal(mesh.positions().iter().map(|p| 1.0 + p[0]).collect());
        let problem = ConductionProblem {
            mesh: &mesh,
            boundary: &boundary,
            material: &material,
            element_materials: None,
            source: &source,
        };
        let result = solve_with_mean_bound(cx, problem, config()).unwrap();
        contains(&result, 300.125); // exact source 1+x gives a rise of 1/8
        let original = fs_conduction::solve(cx, problem, config().primal).unwrap();
        assert_eq!(result.primal.temperature, original.temperature);
        // Verification has not replaced or projected the input source.
        for (i, p) in mesh.positions().iter().enumerate() {
            assert_eq!(source.at(i), 1.0 + p[0]);
        }
        let mut budget = config();
        budget.flux.max_cells = 1;
        assert!(matches!(
            solve_with_mean_bound(cx, problem, budget),
            Err(ConductionBoundError::Verification(TetError::Budget))
        ));
    });
}

#[test]
fn actual_affine_source_primal_and_dual_enclose_cubic_mean_with_mixed_boundaries() {
    with_cx(|cx| {
        for robin in [false, true] {
            let mesh = mesh(2);
            let boundary = ThermalBoundaryBuilder::new(&mesh)
                .region("left", |f| on_box_face(f.centroid[0], 0.0),
                    ThermalBc::dirichlet(300.0).unwrap()).unwrap()
                .region("right", |f| on_box_face(f.centroid[0], 1.0),
                    if robin { ThermalBc::robin(4.0, 298.0).unwrap() }
                    else { ThermalBc::dirichlet(300.0).unwrap() }).unwrap()
                .adiabatic_remainder().finish().unwrap();
            let material = ConductivityModel::constant_tensor(
                [[4.0,0.0,0.0],[0.0,2.0,0.5],[0.0,0.5,1.0]]).unwrap();
            let source = ScalarField::Nodal(mesh.positions().iter().map(|p| 24.0*p[0]).collect());
            let before = source.clone();
            let problem = ConductionProblem { mesh: &mesh, boundary: &boundary,
                material: &material, element_materials: None, source: &source };
            // T=300+x-x^3, -div(K grad T)=24*x, true mean=300.25 K.
            let result = solve_with_mean_bound(cx, problem, config()).unwrap();
            contains(&result, 300.25);
            let direct = fs_conduction::solve(cx, problem, config().primal).unwrap();
            assert_eq!(result.primal.temperature, direct.temperature);
            let reused = fs_conduction::verification::bound_temperature_mean(cx, problem,
                &direct.temperature, config().dual, FluxBudget::default()).unwrap();
            assert_eq!(reused.bound.enclosure, result.bound.enclosure);
            assert_eq!(reused.dual.temperature, result.dual.temperature);
            assert_eq!(source, before);
            assert!(result.primal.temperature.iter().any(|t| *t > 300.0));
            for &(v, _) in boundary.dirichlet() { assert_eq!(result.dual.temperature[v], 0.0); }
        }
    });
}

#[test]
fn unfinished_affine_source_field_retains_its_error_and_actual_residual_correction() {
    with_cx(|cx| {
        let mesh = mesh(2);
        let boundary = ends(&mesh, 300.0, 300.0);
        let material = ConductivityModel::isotropic_declared(1.0).unwrap();
        let source = ScalarField::Nodal(mesh.positions().iter().map(|p| 6.0*p[0]).collect());
        let problem = ConductionProblem { mesh: &mesh, boundary: &boundary,
            material: &material, element_materials: None, source: &source };
        let candidate = vec![300.0; mesh.vertex_count()];
        let result = fs_conduction::verification::bound_temperature_mean(cx, problem,
            &candidate, config().dual, FluxBudget::default()).unwrap();
        let bounds = result.bound.enclosure;
        assert!(bounds.lo <= 300.25 && 300.25 <= bounds.hi);
        assert!(result.bound.candidate_mean.lo <= 300.0 && result.bound.candidate_mean.hi >= 300.0);
        assert!(result.bound.integral.residual_correction.lo > 0.0);
        assert!(result.bound.integral.primal.energy_error_upper > 0.0);
        assert_eq!(candidate, vec![300.0; mesh.vertex_count()]);
        // A nodal representation of a uniform source retains the old result.
        let uniform = ScalarField::Uniform(2.0);
        let constant_nodes = ScalarField::Nodal(vec![2.0; mesh.vertex_count()]);
        let a = solve_with_mean_bound(cx, ConductionProblem { source: &uniform, ..problem }, config()).unwrap();
        let b = solve_with_mean_bound(cx, ConductionProblem { source: &constant_nodes, ..problem }, config()).unwrap();
        assert_eq!(a.primal.temperature, b.primal.temperature);
        assert_eq!(a.dual.temperature, b.dual.temperature);
        assert_eq!(a.bound.enclosure, b.bound.enclosure);
    });
}

#[test]
fn affine_source_support_does_not_relax_nonconstant_boundary_or_bad_field_refusals() {
    with_cx(|cx| {
        let mesh = mesh(2);
        let material = ConductivityModel::isotropic_declared(1.0).unwrap();
        let source = ScalarField::Nodal(mesh.positions().iter().map(|p| 6.0*p[0]).collect());
        let variable_h = ScalarField::Nodal(mesh.positions().iter().map(|p| 1.0+p[1]).collect());
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .region("ambient", |_| true, ThermalBc::Robin {
                htc: variable_h, t_ref: ScalarField::Uniform(300.0) }).unwrap().finish().unwrap();
        let problem = ConductionProblem { mesh: &mesh, boundary: &boundary,
            material: &material, element_materials: None, source: &source };
        assert!(matches!(solve_with_mean_bound(cx, problem, config()),
            Err(ConductionBoundError::Verification(TetError::Unsupported(_)))));
        let fixed = ends(&mesh, 300.0, 300.0);
        let malformed = ScalarField::Nodal(vec![1.0]);
        assert!(matches!(solve_with_mean_bound(cx, ConductionProblem {
            boundary: &fixed, source: &malformed, ..problem }, config()),
            Err(ConductionBoundError::Conduction(_))));
        assert!(fs_conduction::verification::bound_temperature_mean(cx,
            ConductionProblem { boundary: &fixed, ..problem }, &vec![299.0; mesh.vertex_count()],
            config().dual, FluxBudget::default()).is_err());
    });
}
