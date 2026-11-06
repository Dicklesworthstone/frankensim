//! Analytic nonlinear equilibrium and conservation with no reservoir.
#![allow(clippy::float_cmp)]
mod support;

use fs_conduction::bc::{
    SurfaceHeatConduction, ThermalBc, ThermalBoundary, ThermalBoundaryBuilder,
};
use fs_conduction::fixtures::{on_box_face, structured_tets};
use fs_conduction::material::{ElementMaterials, MaterialId, MaterialTable};
use fs_conduction::solve::{MeanTemperatureNewtonConfig, StopRule, solve_mean_temperature_newton};
use fs_conduction::{
    ConductionError, ConductionMesh, ConductionProblem, ConductivityModel, ConductivityTable,
    ScalarField,
};
use support::{with_cancelled_cx, with_cx};

fn mesh() -> ConductionMesh {
    let x = [0.0, 0.1, 0.5, 0.8, 1.0];
    let (complex, positions) =
        structured_tets([4, 2, 2], &|i, j, k| [x[i], j as f64 / 2.0, k as f64 / 2.0]);
    ConductionMesh::new(complex, positions).unwrap()
}

fn material() -> ConductivityModel {
    // k(T) = 1 + 0.02*(T - 300), with no extrapolation.
    ConductivityModel::isotropic(
        ConductivityTable::declared_curve(vec![(280.0, 0.6), (430.0, 3.6)]).unwrap(),
    )
}

fn cooling_fluxes(mesh: &ConductionMesh) -> ThermalBoundary {
    ThermalBoundaryBuilder::new(mesh)
        .region(
            "left",
            |f| on_box_face(f.centroid[0], 0.0),
            ThermalBc::neumann(100.0).unwrap(),
        )
        .unwrap()
        .region(
            "right",
            |f| on_box_face(f.centroid[0], 1.0),
            ThermalBc::neumann(-300.0).unwrap(),
        )
        .unwrap()
        .adiabatic_remainder()
        .finish()
        .unwrap()
}

fn volume_mean(mesh: &ConductionMesh, temperature: &[f64]) -> f64 {
    mesh.complex()
        .tets
        .iter()
        .enumerate()
        .map(|(e, tet)| {
            mesh.element_volume(e) * tet.iter().map(|&v| temperature[v as usize]).sum::<f64>() / 4.0
        })
        .sum::<f64>()
        / mesh.total_volume()
}

#[test]
fn nonlinear_conductivity_preserves_the_exact_manufactured_field_and_every_heat_row() {
    let mesh = mesh();
    let boundary = cooling_fluxes(&mesh);
    let material = material();
    // T=300+100*x: -div(k(T)*grad(T)) = -0.02*100^2 = -200 W/m^3.
    // Outward flux is +100 at x=0 and -300 at x=1. The volume mean is 350 K.
    let source = ScalarField::Uniform(-200.0);
    let problem = ConductionProblem {
        mesh: &mesh,
        boundary: &boundary,
        material: &material,
        element_materials: None,
        source: &source,
    };
    let solved = with_cx(|cx| {
        solve_mean_temperature_newton(cx, problem, 350.0, MeanTemperatureNewtonConfig::default())
            .unwrap()
    });
    for (&actual, p) in solved.temperature.iter().zip(mesh.positions()) {
        assert!(
            (actual - (300.0 + 100.0 * p[0])).abs() < 2e-7,
            "{actual} at {p:?}"
        );
    }
    assert!((volume_mean(&mesh, &solved.temperature) - 350.0).abs() < 1e-10);
    assert_eq!(solved.report.free_dofs, mesh.vertex_count());
    assert_eq!(solved.report.energy.dirichlet_in_w, 0.0);
    assert!(solved.report.energy.closure_w.abs() < 1e-10);
    assert!(solved.report.final_residual <= solved.report.residual_threshold);
    assert!(solved.report.iterations > 1);
    assert_eq!(solved.report.linear.len(), 2 * solved.report.iterations);
    assert!(
        solved
            .report
            .linear
            .iter()
            .all(|r| r.converged_true && r.method == "fgmres")
    );
    let nodal_mean = solved.temperature.iter().sum::<f64>() / mesh.vertex_count() as f64;
    assert!((nodal_mean - 350.0).abs() > 1.0);
    let replay = with_cx(|cx| {
        solve_mean_temperature_newton(cx, problem, 350.0, MeanTemperatureNewtonConfig::default())
            .unwrap()
    });
    assert_eq!(solved, replay);
}

#[test]
fn heterogeneous_nonlinear_equilibrium_has_continuous_interface_heat_flux() {
    let mesh = mesh();
    let boundary = cooling_fluxes(&mesh);
    // Left: k=1+0.02*(T-300), T'=100.
    // Right: k=4+0.08*(T-350), T'=50.
    // Both have -div(k*T')=-200; at x=.5 their conductive flux is -200 W/m^2.
    let table = MaterialTable::new([
        (MaterialId(1), material()),
        (
            MaterialId(2),
            ConductivityModel::isotropic(
                ConductivityTable::declared_curve(vec![(320.0, 1.6), (420.0, 9.6)]).unwrap(),
            ),
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
            MaterialId(if x < 0.5 { 1 } else { 2 })
        })
        .collect();
    let assignment = ElementMaterials::new(table, ids).unwrap();
    // The uniform fallback is deliberately inapplicable; the assignment owns the law.
    let fallback = ConductivityModel::isotropic(
        ConductivityTable::declared_curve(vec![(800.0, 8.0), (900.0, 9.0)]).unwrap(),
    );
    let source = ScalarField::Uniform(-200.0);
    let problem = ConductionProblem {
        mesh: &mesh,
        boundary: &boundary,
        material: &fallback,
        element_materials: Some(&assignment),
        source: &source,
    };
    let solved = with_cx(|cx| {
        solve_mean_temperature_newton(cx, problem, 343.75, MeanTemperatureNewtonConfig::default())
            .unwrap()
    });
    for (&actual, p) in solved.temperature.iter().zip(mesh.positions()) {
        let expected = 300.0 + 100.0 * p[0].min(0.5) + 50.0 * (p[0] - 0.5).max(0.0);
        assert!((actual - expected).abs() < 2e-7);
    }
    assert!((volume_mean(&mesh, &solved.temperature) - 343.75).abs() < 1e-10);
    assert_eq!(
        solved.report.element_material_identity,
        Some(assignment.identity())
    );
    assert_eq!(solved.report.energy.dirichlet_in_w, 0.0);
    assert!(solved.report.final_residual <= solved.report.residual_threshold);
}

#[test]
fn insulated_body_with_nonlinear_sheet_relaxes_to_its_declared_mean() {
    let mesh = mesh();
    let material = material();
    let source = ScalarField::Uniform(0.0);
    let sheet_face = mesh
        .boundary()
        .iter()
        .position(|face| on_box_face(face.centroid[1], 0.0))
        .unwrap();
    let sheet = SurfaceHeatConduction::from_isotropic_conductivity(
        sheet_face,
        ConductivityTable::declared_curve(vec![(320.0, 5.0), (380.0, 12.0)]).unwrap(),
        0.01,
    )
    .unwrap();
    let boundary = with_cx(|cx| {
        ThermalBoundaryBuilder::new(&mesh)
            .adiabatic_remainder()
            .finish()
            .unwrap()
            .with_surface_heat_conduction(cx, &mesh, &[sheet])
            .unwrap()
    });
    let problem = ConductionProblem {
        mesh: &mesh,
        boundary: &boundary,
        material: &material,
        element_materials: None,
        source: &source,
    };
    let config = MeanTemperatureNewtonConfig {
        initial_temperature: Some(
            mesh.positions()
                .iter()
                .map(|p| 350.0 + 40.0 * (p[0] - 0.5))
                .collect(),
        ),
        ..MeanTemperatureNewtonConfig::default()
    };
    let solved = with_cx(|cx| solve_mean_temperature_newton(cx, problem, 350.0, config).unwrap());
    assert!(solved.temperature.iter().all(|t| (t - 350.0).abs() < 2e-8));
    assert!((volume_mean(&mesh, &solved.temperature) - 350.0).abs() < 1e-10);
    assert_eq!(solved.report.energy.dirichlet_in_w, 0.0);
    assert!(solved.report.final_residual <= solved.report.residual_threshold);
    with_cancelled_cx(|cx| {
        assert!(matches!(
            solve_mean_temperature_newton(
                cx,
                problem,
                350.0,
                MeanTemperatureNewtonConfig::default()
            ),
            Err(ConductionError::Cancelled { .. })
        ))
    });
}

#[test]
fn incompatible_heat_and_the_omitted_reference_row_cannot_be_hidden() {
    let mesh = mesh();
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .adiabatic_remainder()
        .finish()
        .unwrap();
    let material = material();
    let source = ScalarField::Uniform(1.0);
    let problem = ConductionProblem {
        mesh: &mesh,
        boundary: &boundary,
        material: &material,
        element_materials: None,
        source: &source,
    };
    let config = MeanTemperatureNewtonConfig {
        stop: StopRule {
            residual_atol: 1e9,
            ..StopRule::default()
        },
        ..MeanTemperatureNewtonConfig::default()
    };
    let error =
        with_cx(|cx| solve_mean_temperature_newton(cx, problem, 350.0, config).unwrap_err());
    assert!(error.to_string().contains("no steady heat balance"));
    let small_source = ScalarField::Uniform(1e-4);
    let config = MeanTemperatureNewtonConfig {
        net_heat_tolerance_w: 1e-3,
        ..MeanTemperatureNewtonConfig::default()
    };
    // Compatibility deliberately permits this imbalance. All original rows
    // still cannot balance; dropping row 0 from the final test would hide it.
    assert!(
        with_cx(|cx| solve_mean_temperature_newton(
            cx,
            ConductionProblem {
                source: &small_source,
                ..problem
            },
            350.0,
            config
        ))
        .is_err()
    );
}

#[test]
fn material_support_and_nonlinear_iteration_budget_remain_authoritative() {
    let mesh = mesh();
    let boundary = cooling_fluxes(&mesh);
    let material = material();
    let source = ScalarField::Uniform(-200.0);
    let problem = ConductionProblem {
        mesh: &mesh,
        boundary: &boundary,
        material: &material,
        element_materials: None,
        source: &source,
    };
    assert!(matches!(
        with_cx(|cx| solve_mean_temperature_newton(
            cx,
            problem,
            450.0,
            MeanTemperatureNewtonConfig::default()
        )),
        Err(ConductionError::OutsideTemperatureSpan { .. })
    ));
    let config = MeanTemperatureNewtonConfig {
        stop: StopRule {
            max_iterations: 1,
            ..StopRule::default()
        },
        ..MeanTemperatureNewtonConfig::default()
    };
    assert!(matches!(
        with_cx(|cx| solve_mean_temperature_newton(cx, problem, 350.0, config)),
        Err(ConductionError::NotConverged { iterations: 1, .. })
    ));
}

#[test]
fn inadmissible_full_newton_step_backtracks_inside_the_original_material_domain() {
    let mesh = mesh();
    let boundary = cooling_fluxes(&mesh);
    let material = material();
    let source = ScalarField::Uniform(-200.0);
    let problem = ConductionProblem {
        mesh: &mesh,
        boundary: &boundary,
        material: &material,
        element_materials: None,
        source: &source,
    };
    // The reversed gradient starts in [300,400] K but its full Newton step
    // overshoots the original 430 K upper support. Backtracking is essential.
    let initial: Vec<_> = mesh
        .positions()
        .iter()
        .map(|p| 350.0 - 100.0 * (p[0] - 0.5))
        .collect();
    let config = MeanTemperatureNewtonConfig {
        initial_temperature: Some(initial),
        ..MeanTemperatureNewtonConfig::default()
    };
    let mut no_backtracking = config.clone();
    no_backtracking.line_search.max_backtracks = 0;
    assert!(matches!(
        with_cx(|cx| solve_mean_temperature_newton(cx, problem, 350.0, no_backtracking)),
        Err(ConductionError::OutsideTemperatureSpan { .. })
    ));
    let solved = with_cx(|cx| solve_mean_temperature_newton(cx, problem, 350.0, config).unwrap());
    for (&actual, p) in solved.temperature.iter().zip(mesh.positions()) {
        assert!((actual - (300.0 + 100.0 * p[0])).abs() < 2e-7);
    }
    assert!((volume_mean(&mesh, &solved.temperature) - 350.0).abs() < 1e-10);
}
