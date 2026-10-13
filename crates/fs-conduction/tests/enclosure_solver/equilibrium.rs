//! Physical reassembly at actual, rather than lagged, enclosure temperatures.
use super::*;
use fs_conduction::radiation::EnclosureEquilibriumReport;

/// Public assembly route for the test's no-contact fixtures. Both heterogeneous
/// materials and the exact physical Neumann loads are independently rebound.
pub(super) fn actual_equilibrium(
    cx: &fs_exec::Cx<'_>, problem: ConductionProblem<'_>, enclosure: &GrayDiffuseEnclosure,
    solution: &fs_conduction::ConductionSolution, applied: &fs_conduction::RadiosityReport,
) -> EnclosureEquilibriumReport {
    let means = enclosure.surfaces().iter().map(|surface|
        surface.mean_temperature(problem.mesh, &solution.temperature).unwrap()).collect::<Vec<_>>();
    let nonlinear = enclosure.solve(cx, &means).unwrap();
    let rows = enclosure.surfaces().iter().zip(&nonlinear.net_outward_flux_w_m2)
        .map(|(s, &flux)| (s.name().to_owned(), s.face_slots().to_vec(), flux)).collect::<Vec<_>>();
    let boundary = problem.boundary.with_uniform_outward_flux_overlays(problem.mesh, &rows).unwrap();
    let system = match problem.element_materials {
        Some(materials) => fs_conduction::assemble::assemble_operator_with_element_materials(
            cx, problem.mesh, &boundary, problem.material, problem.source, &solution.temperature, materials),
        None => fs_conduction::assemble::assemble_operator(
            cx, problem.mesh, &boundary, problem.material, problem.source, &solution.temperature),
    }.unwrap();
    let dofs = fs_conduction::DofMap::new(&boundary, problem.mesh.vertex_count()).unwrap();
    let residual = fs_conduction::assemble::residual(&system, &dofs, &solution.temperature);
    let mismatch = nonlinear.net_outward_heat_w.iter().zip(&applied.net_outward_heat_w)
        .map(|(new, old)| (new-old).abs()).fold(0.0f64, f64::max);
    EnclosureEquilibriumReport {
        nonlinear_radiosity: nonlinear, residual_w: fs_solver::norm2(&residual),
        residual_threshold_w: solution.report.residual_threshold,
        max_patch_heat_mismatch_w: mismatch,
    }
}

fn loose_outer(iterations: usize) -> CoupledRadiationConfig {
    CoupledRadiationConfig { surface_temperature_rtol: 0.0,
        surface_temperature_atol_k: 1000.0, relaxation: 0.25, max_iterations: iterations }
}

#[test]
fn small_updates_and_closed_total_heat_cannot_accept_wrong_surface_loads() {
    let mesh = two_disconnected_slabs();
    let boundary = two_slab_boundary(&mesh);
    let enclosure = two_slab_enclosure(&mesh);
    let material = ConductivityModel::isotropic_declared(10.0).unwrap();
    let source = ScalarField::Uniform(0.0);
    let problem = ConductionProblem { mesh: &mesh, boundary: &boundary,
        material: &material, element_materials: None, source: &source };
    with_cx(|cx| {
        // Reproduce the single candidate that the temperature-only gate accepts.
        let initial = fs_conduction::solve(cx, problem, config()).unwrap();
        let old_means = enclosure.surfaces().iter().map(|s|
            s.mean_temperature(&mesh, &initial.temperature).unwrap()).collect::<Vec<_>>();
        let applied = enclosure.solve(cx, &old_means).unwrap();
        let rows = enclosure.surfaces().iter().zip(&applied.net_outward_flux_w_m2)
            .map(|(s, &q)| (s.name().to_owned(), s.face_slots().to_vec(), q)).collect::<Vec<_>>();
        let overlaid = boundary.with_uniform_outward_flux_overlays(&mesh, &rows).unwrap();
        let candidate = fs_conduction::solve(cx, ConductionProblem { boundary: &overlaid, ..problem }, config()).unwrap();
        assert!(candidate.report.final_residual <= candidate.report.residual_threshold);
        for (surface, &old) in enclosure.surfaces().iter().zip(&old_means) {
            assert!((surface.mean_temperature(&mesh, &candidate.temperature).unwrap()-old).abs() < 1000.0);
        }
        let check = actual_equilibrium(cx, problem, &enclosure, &candidate, &applied);
        assert!(applied.enclosure_energy_closure_w.abs() < 1e-8);
        assert!(check.nonlinear_radiosity.enclosure_energy_closure_w.abs() < 1e-8);
        assert!(check.max_patch_heat_mismatch_w > 100.0);
        assert!(check.residual_w > 1000.0*check.residual_threshold_w);
        let error = solve_with_gray_diffuse_enclosure(cx, problem, None, &enclosure,
            config(), loose_outer(1)).unwrap_err();
        assert!(matches!(error, ConductionError::Radiation { what, .. }
            if what.contains("did not converge in 1 iterations") && what.contains("nonlinear residual")));
    });
}

#[test]
fn actual_driver_continues_until_the_original_nonlinear_equilibrium_gate_passes() {
    let mesh = two_disconnected_slabs();
    let boundary = two_slab_boundary(&mesh);
    let enclosure = two_slab_enclosure(&mesh);
    let material = ConductivityModel::isotropic_declared(10.0).unwrap();
    let source = ScalarField::Uniform(0.0);
    let problem = ConductionProblem { mesh: &mesh, boundary: &boundary,
        material: &material, element_materials: None, source: &source };
    with_cx(|cx| {
        let solution = solve_with_gray_diffuse_enclosure(cx, problem, None, &enclosure,
            config(), loose_outer(256)).unwrap();
        assert!(solution.radiation.iterations > 1);
        assert_eq!(solution.radiation.equilibrium_checks, solution.radiation.iterations);
        let expected = actual_equilibrium(cx, problem, &enclosure, &solution.conduction,
            &solution.radiation.radiosity);
        assert_eq!(solution.radiation.equilibrium, expected);
        assert!(expected.residual_w <= expected.residual_threshold_w);
        let means = &expected.nonlinear_radiosity.surface_temperatures_k;
        let q = 5.670_374_419e-8*(means[0]-means[1])*(means[0]+means[1])
            *(means[0]*means[0]+means[1]*means[1])/1.5;
        assert!((20.0*(400.0-means[0])-q).abs() < 2e-5);
        assert!((20.0*(means[1]-300.0)-q).abs() < 2e-5);
        let replay = solve_with_gray_diffuse_enclosure(cx, problem, None, &enclosure,
            config(), loose_outer(256)).unwrap();
        assert_eq!(solution, replay);
    });
}

#[test]
fn final_state_check_keeps_nonlinear_assignments_and_ignores_unused_fallback() {
    use fs_conduction::{ConductivityTable, ElementMaterials, MaterialId, MaterialTable};
    let mesh = two_disconnected_slabs();
    let boundary = two_slab_boundary(&mesh);
    let enclosure = two_slab_enclosure(&mesh);
    let materials = ElementMaterials::new(MaterialTable::new([
        (MaterialId(0), ConductivityModel::isotropic(ConductivityTable::declared_curve(
            vec![(250.0,7.0),(500.0,12.0)]).unwrap())),
        (MaterialId(1), ConductivityModel::isotropic(ConductivityTable::declared_curve(
            vec![(250.0,17.0),(500.0,22.0)]).unwrap())),
    ]).unwrap(), mesh.complex().tets.iter().map(|tet| MaterialId(
        if mesh.positions()[tet[0] as usize][0] < 1.0 { 0 } else { 1 })).collect()).unwrap();
    let source = ScalarField::Uniform(0.0);
    with_cx(|cx| {
        let mut previous = None;
        for fallback in [1.0,30.0] {
            let fallback = ConductivityModel::isotropic_declared(fallback).unwrap();
            let problem = ConductionProblem { mesh: &mesh, boundary: &boundary, material: &fallback,
                element_materials: Some(&materials), source: &source };
            let solution = solve_with_gray_diffuse_enclosure(cx, problem, None, &enclosure,
                config(), loose_outer(256)).unwrap();
            let expected = actual_equilibrium(cx, problem, &enclosure, &solution.conduction,
                &solution.radiation.radiosity);
            assert_eq!(solution.radiation.equilibrium, expected);
            assert!(expected.residual_w <= expected.residual_threshold_w);
            if let Some(prior) = previous { assert_eq!(solution, prior); }
            previous = Some(solution);
        }
    });
}
