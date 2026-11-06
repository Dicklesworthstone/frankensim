//! G1 analytic Neumann fields and G0 heat compatibility/cancellation checks.
#![allow(clippy::float_cmp)]
mod support;

use fs_conduction::bc::{SurfaceHeatConduction, ThermalBc, ThermalBoundary, ThermalBoundaryBuilder};
use fs_conduction::fixtures::{box_grid, on_box_face, structured_tets};
use fs_conduction::material::{ElementMaterials, MaterialId, MaterialTable};
use fs_conduction::solve::{MeanTemperatureConfig, solve_mean_temperature};
use fs_conduction::{ConductionError, ConductionMesh, ConductionProblem, ConductivityModel,
    ConductivityTable, ScalarField};
use fs_rep_mesh::TetComplex;
use support::{with_cancelled_cx, with_cx};

fn mesh() -> ConductionMesh {
    let x = [0.0, 0.1, 0.5, 0.8, 1.0];
    let (complex, positions) = structured_tets([4, 2, 2],
        &|i, j, k| [x[i], j as f64 / 2.0, k as f64 / 2.0]);
    ConductionMesh::new(complex, positions).unwrap()
}

fn fluxes(mesh: &ConductionMesh, heat: f64) -> ThermalBoundary {
    ThermalBoundaryBuilder::new(mesh)
        .region("left", |f| on_box_face(f.centroid[0], 0.0), ThermalBc::neumann(heat).unwrap())
        .unwrap()
        .region("right", |f| on_box_face(f.centroid[0], 1.0), ThermalBc::neumann(-heat).unwrap())
        .unwrap().adiabatic_remainder().finish().unwrap()
}

fn volume_mean(mesh: &ConductionMesh, values: &[f64]) -> f64 {
    mesh.complex().tets.iter().enumerate().map(|(e, tet)| {
        mesh.element_volume(e) * tet.iter().map(|&v| values[v as usize]).sum::<f64>() / 4.0
    }).sum::<f64>() / mesh.total_volume()
}

#[test]
fn anisotropic_linear_flux_sets_volume_mean_and_retains_every_heat_equation() {
    let mesh = mesh();
    let boundary = fluxes(&mesh, 8.0);
    let material = ConductivityModel::constant_tensor([
        [2.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 5.0],
    ]).unwrap();
    let source = ScalarField::Uniform(0.0);
    let problem = ConductionProblem {
        mesh: &mesh, boundary: &boundary, material: &material, element_materials: None, source: &source,
    };
    for mean in [310.0, 810.0] {
        let solved = with_cx(|cx| solve_mean_temperature(cx, problem, mean,
            MeanTemperatureConfig::default()).unwrap());
        for (&actual, p) in solved.temperature.iter().zip(mesh.positions()) {
            assert!((actual - (mean + 4.0 * (p[0] - 0.5))).abs() < 2e-9);
        }
        assert!((volume_mean(&mesh, &solved.temperature) - mean).abs() < 1e-10);
        // Uneven x spacing makes the nodal arithmetic mean the wrong datum.
        let nodal_mean = solved.temperature.iter().sum::<f64>() / mesh.vertex_count() as f64;
        assert!((nodal_mean - mean).abs() > 0.05);
        assert_eq!(solved.report.free_dofs, mesh.vertex_count());
        assert_eq!(solved.report.energy.dirichlet_in_w, 0.0);
        assert!(solved.report.energy.closure_w.abs() < 1e-12);
        assert!(solved.report.final_residual <= solved.report.residual_threshold);
        assert!(solved.report.linear[0].converged_true);
        let replay = with_cx(|cx| solve_mean_temperature(cx, problem, mean,
            MeanTemperatureConfig::default()).unwrap());
        assert_eq!(solved, replay);
    }
}

#[test]
fn assigned_materials_preserve_piecewise_flux_and_ignore_unused_fallback() {
    let mesh = mesh();
    let boundary = fluxes(&mesh, 8.0);
    let table = MaterialTable::new([
        (MaterialId(1), ConductivityModel::isotropic_declared(2.0).unwrap()),
        (MaterialId(2), ConductivityModel::isotropic_declared(4.0).unwrap()),
    ]).unwrap();
    let ids = mesh.complex().tets.iter().map(|tet| {
        let x = tet.iter().map(|&v| mesh.positions()[v as usize][0]).sum::<f64>() / 4.0;
        MaterialId(if x < 0.5 { 1 } else { 2 })
    }).collect();
    let assignment = ElementMaterials::new(table, ids).unwrap();
    // This nonlinear fallback is not part of the assigned physical material.
    let fallback = ConductivityModel::isotropic(
        ConductivityTable::declared_curve(vec![(300.0, 1.0), (400.0, 2.0)]).unwrap());
    let source = ScalarField::Uniform(0.0);
    let problem = ConductionProblem {
        mesh: &mesh, boundary: &boundary, material: &fallback,
        element_materials: Some(&assignment), source: &source,
    };
    let solved = with_cx(|cx| solve_mean_temperature(cx, problem, 320.0,
        MeanTemperatureConfig::default()).unwrap());
    for (&actual, p) in solved.temperature.iter().zip(mesh.positions()) {
        let shape = if p[0] <= 0.5 { 4.0 * p[0] } else { 2.0 * p[0] + 1.0 };
        assert!((actual - (320.0 - 1.75 + shape)).abs() < 2e-9);
    }
    assert_eq!(solved.report.element_material_identity, Some(assignment.identity()));
    assert!((volume_mean(&mesh, &solved.temperature) - 320.0).abs() < 1e-10);
}

#[test]
fn insulated_equilibrium_with_sampled_sheet_needs_no_time_marching_or_reservoir() {
    let mesh = mesh();
    let material = ConductivityModel::isotropic_declared(2.0).unwrap();
    let sheet = SurfaceHeatConduction::from_isotropic_conductivity(0,
        ConductivityTable::declared_curve(vec![(300.0, 5.0), (400.0, 5.0)]).unwrap(), 0.01).unwrap();
    let boundary = with_cx(|cx| ThermalBoundaryBuilder::new(&mesh)
        .adiabatic_remainder().finish().unwrap()
        .with_surface_heat_conduction(cx, &mesh, &[sheet]).unwrap());
    let source = ScalarField::Uniform(0.0);
    let problem = ConductionProblem {
        mesh: &mesh, boundary: &boundary, material: &material, element_materials: None, source: &source,
    };
    let solved = with_cx(|cx| solve_mean_temperature(cx, problem, 321.0,
        MeanTemperatureConfig::default()).unwrap());
    assert_eq!(solved.temperature, vec![321.0; mesh.vertex_count()]);
    assert_eq!(solved.report.linear[0].iterations, 0);
    assert_eq!(solved.report.energy.dirichlet_in_w, 0.0);
    with_cancelled_cx(|cx| assert!(matches!(
        solve_mean_temperature(cx, problem, 321.0, MeanTemperatureConfig::default()),
        Err(ConductionError::Cancelled { .. })
    )));
}

#[test]
fn incompatible_heating_is_refused_instead_of_being_projected_away() {
    let mesh = mesh();
    let boundary = ThermalBoundaryBuilder::new(&mesh).adiabatic_remainder().finish().unwrap();
    let material = ConductivityModel::isotropic_declared(2.0).unwrap();
    let source = ScalarField::Uniform(1.0);
    let problem = ConductionProblem {
        mesh: &mesh, boundary: &boundary, material: &material, element_materials: None, source: &source,
    };
    let config = MeanTemperatureConfig { residual_atol_w: 1e9, ..MeanTemperatureConfig::default() };
    let error = with_cx(|cx| solve_mean_temperature(cx, problem, 300.0, config).unwrap_err());
    assert!(error.to_string().contains("no steady heat balance"));
}

#[test]
fn temporarily_eliminated_vertex_still_has_to_satisfy_the_original_heat_equation() {
    let mesh = mesh();
    let boundary = ThermalBoundaryBuilder::new(&mesh).adiabatic_remainder().finish().unwrap();
    let material = ConductivityModel::isotropic_declared(2.0).unwrap();
    let source = ScalarField::Uniform(1e-4);
    let problem = ConductionProblem {
        mesh: &mesh, boundary: &boundary, material: &material, element_materials: None, source: &source,
    };
    // Compatibility alone allows this imbalance, but the independent full
    // equilibrium residual must catch the artificial algebraic-pin reaction.
    let config = MeanTemperatureConfig { net_heat_tolerance_w: 1e-3, ..MeanTemperatureConfig::default() };
    assert!(matches!(with_cx(|cx| solve_mean_temperature(cx, problem, 300.0, config)),
        Err(ConductionError::NotConverged { .. })));
}

#[test]
fn conductivity_domains_and_nonlinear_laws_remain_authoritative() {
    let mesh = mesh();
    let boundary = ThermalBoundaryBuilder::new(&mesh).adiabatic_remainder().finish().unwrap();
    let source = ScalarField::Uniform(0.0);
    let limited = ConductivityModel::isotropic(
        ConductivityTable::declared_curve(vec![(300.0, 2.0), (310.0, 2.0)]).unwrap());
    let problem = ConductionProblem {
        mesh: &mesh, boundary: &boundary, material: &limited, element_materials: None, source: &source,
    };
    assert!(matches!(with_cx(|cx| solve_mean_temperature(cx, problem, 320.0,
        MeanTemperatureConfig::default())), Err(ConductionError::OutsideTemperatureSpan { .. })));
    let nonlinear = ConductivityModel::isotropic(
        ConductivityTable::declared_curve(vec![(300.0, 2.0), (400.0, 3.0)]).unwrap());
    assert!(matches!(with_cx(|cx| solve_mean_temperature(cx,
        ConductionProblem { material: &nonlinear, ..problem }, 320.0,
        MeanTemperatureConfig::default())), Err(ConductionError::Config { .. })));
    let sheet = SurfaceHeatConduction::from_isotropic_conductivity(0,
        ConductivityTable::declared_curve(vec![(300.0, 2.0), (400.0, 3.0)]).unwrap(), 0.01).unwrap();
    let coated = with_cx(|cx| boundary.clone().with_surface_heat_conduction(cx, &mesh, &[sheet]).unwrap());
    assert!(matches!(with_cx(|cx| solve_mean_temperature(cx,
        ConductionProblem { boundary: &coated, ..problem }, 305.0,
        MeanTemperatureConfig::default())), Err(ConductionError::Config { .. })));
}

#[test]
fn one_mean_does_not_fabricate_equilibria_for_two_disconnected_bodies() {
    let (first, mut positions) = box_grid([1, 1, 1], [1.0, 1.0, 1.0]);
    let offset = positions.len() as u32;
    let extra = positions.iter().map(|p| [p[0] + 3.0, p[1], p[2]]).collect::<Vec<_>>();
    positions.extend(extra);
    let mut tets = first.tets.clone();
    tets.extend(first.tets.iter().map(|tet| tet.map(|v| v + offset)));
    let mesh = ConductionMesh::new(TetComplex::from_tets(positions.len(), tets), positions).unwrap();
    let boundary = ThermalBoundaryBuilder::new(&mesh).adiabatic_remainder().finish().unwrap();
    let material = ConductivityModel::isotropic_declared(2.0).unwrap();
    let source = ScalarField::Uniform(0.0);
    let problem = ConductionProblem {
        mesh: &mesh, boundary: &boundary, material: &material, element_materials: None, source: &source,
    };
    let error = with_cx(|cx| solve_mean_temperature(cx, problem, 300.0,
        MeanTemperatureConfig::default()).unwrap_err());
    assert!(error.to_string().contains("one connected body"));
}
