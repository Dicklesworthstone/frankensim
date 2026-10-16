//! A small outer update cannot turn a step-only solid stop into equilibrium.
use super::*;

fn step_only_config() -> SolveConfig {
    let mut cfg = config();
    cfg.nonlinearity = Nonlinearity::FixedPoint { relaxation: 1e-6, max_backtracks: 8 };
    cfg.stop.step_atol = 1.0;
    cfg.initial = InitialGuess::Uniform(350.0);
    cfg
}

fn assert_solid_refusal(result: Result<fs_conduction::CoupledRadiationSolution, ConductionError>) {
    assert!(matches!(result, Err(ConductionError::Radiation { surface, what, .. })
        if surface == "enclosure-solid" && what.contains("residual gate")));
}

#[test]
fn an_incomplete_initial_solid_is_not_an_enclosure_solution() {
    let mesh = two_disconnected_slabs();
    let boundary = two_slab_boundary(&mesh);
    let enclosure = two_slab_enclosure(&mesh);
    let material = ConductivityModel::isotropic_declared(10.0).unwrap();
    let source = ScalarField::Uniform(0.0);
    let problem = ConductionProblem { mesh: &mesh, boundary: &boundary, material: &material,
        element_materials: None, source: &source };
    let cfg = step_only_config();
    with_cx(|cx| {
        // The general solver legitimately supports a declared step-only stop.
        // A coupled equilibrium consumer must not mistake that for R = 0.
        let raw = fs_conduction::solve(cx, problem, cfg.clone()).unwrap();
        assert!(raw.report.final_residual > raw.report.residual_threshold,
            "negative control must actually reach the incomplete-solid boundary");
        assert_solid_refusal(solve_with_gray_diffuse_enclosure(cx, problem, None,
            &enclosure, cfg, CoupledRadiationConfig::default()));
    });
}

#[test]
fn the_residual_gate_also_applies_after_radiation_changes_the_boundary() {
    let mesh = two_disconnected_slabs();
    let boundary = two_slab_boundary(&mesh);
    let enclosure = two_slab_enclosure(&mesh);
    let material = ConductivityModel::isotropic_declared(10.0).unwrap();
    let source = ScalarField::Uniform(0.0);
    let problem = ConductionProblem { mesh: &mesh, boundary: &boundary, material: &material,
        element_materials: None, source: &source };
    let exact_initial: Vec<_> = mesh.positions().iter()
        .map(|p| if p[0] < 1.0 { 400.0 } else { 300.0 }).collect();
    let mut cfg = step_only_config();
    cfg.initial = InitialGuess::Free(fs_conduction::DofMap::new(&boundary, mesh.vertex_count())
        .unwrap().gather(&exact_initial));
    with_cx(|cx| {
        let raw = fs_conduction::solve(cx, problem, cfg.clone()).unwrap();
        assert!(raw.report.final_residual <= raw.report.residual_threshold,
            "the initial nonradiating state must pass; refusal belongs to a later solid solve");
        // A very loose outer stop would otherwise accept the first tiny update.
        let outer = CoupledRadiationConfig { surface_temperature_rtol: 0.0,
            surface_temperature_atol_k: 1.0, relaxation: 0.5, max_iterations: 100 };
        assert_solid_refusal(solve_with_gray_diffuse_enclosure(cx, problem, None,
            &enclosure, cfg, outer));
    });
}
