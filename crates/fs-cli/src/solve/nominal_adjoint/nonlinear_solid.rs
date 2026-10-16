//! Native selected-temperature gradients with the full heterogeneous k(T)
//! Jacobian. Radiation and air consumers select their own COMPLETE feedback
//! branches before reaching this one. Nothing here freezes another state law.

use fs_conduction::{LinearConfig, ThermalInterfaces};
use fs_conduction::adjoint::{RobinGradient, RobinResponse};
use super::{ConductionProblem, Cx, SolveRefusal, lower, poll};

pub(super) fn needed(cx: &Cx<'_>, problem: ConductionProblem<'_>)
    -> Result<bool, SolveRefusal>
{
    poll(cx)?;
    if let Some(materials) = problem.element_materials {
        materials.validate_for(problem.mesh).map_err(lower)?;
        for element in 0..problem.mesh.element_count() {
            if element % 512 == 0 { poll(cx)?; }
            if materials.model_for(element).map_err(lower)?.is_temperature_dependent() {
                return Ok(true);
            }
        }
        // An unused fallback must not select the algorithm or its validity.
        Ok(false)
    } else {
        Ok(problem.material.is_temperature_dependent())
    }
}

pub(super) fn pullback(
    cx: &Cx<'_>, problem: ConductionProblem<'_>, interfaces: Option<&ThermalInterfaces>,
    linear: LinearConfig, temperature: &[f64], weights: &[f64],
) -> Result<RobinGradient, SolveRefusal> {
    // Zero mean-feedback ports: the production owner still assembles the
    // actual material/contact Jacobian, verifies the unmodified primal and
    // solves its true transpose with exact inner-iteration/cancellation gates.
    // A small residual of a frozen conductivity operator is never substituted.
    // No inverse/stability certificate is claimed for this nonlinear mode.
    RobinResponse::pullback_mean_htc(cx, problem, interfaces, linear, temperature,
        &[], &[], weights, 0).map_err(lower)
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_conduction::{ConductionMesh, ConductivityModel, ConductivityTable,
        ElementMaterials, InitialGuess, MaterialId, MaterialTable, ScalarField,
        SolveConfig, ThermalBc, ThermalBoundary, ThermalBoundaryBuilder};
    use fs_conduction::fixtures::{box_grid, on_box_face};
    use fs_alloc::{ArenaConfig, ArenaPool};
    use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};

    fn with_cx<T>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> T) -> T {
        ArenaPool::new(ArenaConfig::default()).scope(|arena| f(&Cx::new(gate, arena,
            StreamKey { seed: 57, kernel_id: 825, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic)))
    }
    fn mesh() -> ConductionMesh {
        let (complex, positions) = box_grid([3, 2, 2], [0.1, 0.04, 0.03]);
        ConductionMesh::new(complex, positions).unwrap()
    }
    fn config() -> SolveConfig {
        let mut cfg = SolveConfig::default();
        cfg.initial = InitialGuess::Uniform(320.0);
        cfg.stop.step_atol = 0.0;
        cfg.stop.residual_rtol = 1e-12;
        cfg.linear.tolerance = 1e-11;
        cfg
    }
    fn material() -> ConductivityModel {
        ConductivityModel::isotropic(ConductivityTable::declared_curve(
            vec![(250.0, 2.0), (450.0, 22.0)]).unwrap())
    }
    fn boundary(mesh: &ConductionMesh, h: f64, reference: f64) -> ThermalBoundary {
        ThermalBoundaryBuilder::new(mesh)
            .region("hot", |f| on_box_face(f.centroid[0], 0.0), ThermalBc::dirichlet(340.0).unwrap()).unwrap()
            .region("wall", |f| on_box_face(f.centroid[0], 0.1), ThermalBc::robin(h, reference).unwrap()).unwrap()
            .adiabatic_remainder().finish().unwrap()
    }
    fn solve(cx: &Cx<'_>, density: f64, h: f64, reference: f64) -> Vec<f64> {
        let mesh = mesh();
        let material = material();
        let source = ScalarField::Nodal(mesh.positions().iter()
            .map(|p| density + 40000.0*p[1]).collect());
        let boundary = boundary(&mesh, h, reference);
        fs_conduction::solve(cx, ConductionProblem { mesh: &mesh, boundary: &boundary,
            material: &material, element_materials: None, source: &source }, config()).unwrap().temperature
    }

    #[test]
    fn nonlinear_native_route_matches_source_and_boundary_physical_resolves() {
        with_cx(&CancelGate::new_clock_free(), |cx| {
            let mesh = mesh();
            let material = material();
            let source = ScalarField::Nodal(mesh.positions().iter()
                .map(|p| 3000.0 + 40000.0*p[1]).collect());
            let boundary = boundary(&mesh, 40.0, 293.0);
            let problem = ConductionProblem { mesh: &mesh, boundary: &boundary,
                material: &material, element_materials: None, source: &source };
            assert!(needed(cx, problem).unwrap());
            let temperature = solve(cx, 3000.0, 40.0, 293.0);
            let retained = temperature.clone();
            let n = mesh.vertex_count();
            let mut weights = vec![0.0; n]; weights[n-1] = 1.0;
            let gradient = pullback(cx, problem, None, config().linear, &temperature, &weights).unwrap();
            assert_eq!(temperature, retained);
            assert!(gradient.relative_residual < config().linear.tolerance);
            let lambda = &gradient.nodal_load;
            let density = mesh.complex().tets.iter().enumerate().map(|(e, tet)|
                mesh.element_volume(e)/4.0 * tet.iter().map(|&v| lambda[v as usize]).sum::<f64>()).sum::<f64>();
            let mut h_bar = 0.0;
            let mut ref_bar = 0.0;
            for face in mesh.boundary().iter().filter(|f| on_box_face(f.centroid[0], 0.1)) {
                let vertices = face.vertices.map(|v| v as usize);
                let bars = super::super::contractions::boundary(face.area,
                    vertices.map(|v| lambda[v]), vertices.map(|v| temperature[v]), 40.0, 293.0);
                h_bar += bars[0]; ref_bar += bars[1];
            }
            for (actual, high, low, step) in [
                (density, solve(cx,3001.0,40.0,293.0), solve(cx,2999.0,40.0,293.0), 1.0),
                (h_bar, solve(cx,3000.0,40.001,293.0), solve(cx,3000.0,39.999,293.0), 0.001),
                (ref_bar, solve(cx,3000.0,40.0,293.01), solve(cx,3000.0,40.0,292.99), 0.01),
            ] {
                let expected = (high[n-1]-low[n-1])/(2.0*step);
                assert!((actual-expected).abs() < 5e-5*expected.abs().max(1e-6),
                    "complete nonlinear adjoint {actual:e} versus physical resolve {expected:e}");
            }
            for &v in fs_conduction::DofMap::new(&boundary, n).unwrap().fixed() {
                assert_eq!(lambda[v], 0.0);
            }
            let assigned = ElementMaterials::new(MaterialTable::new([
                (MaterialId(7), material.clone())]).unwrap(), vec![MaterialId(7); mesh.element_count()]).unwrap();
            let unused = ConductivityModel::isotropic(ConductivityTable::declared_curve(
                vec![(0.0, 1.0), (1.0, 2.0)]).unwrap());
            let assigned_problem = ConductionProblem { material: &unused, element_materials: Some(&assigned), ..problem };
            assert!(needed(cx, assigned_problem).unwrap());
            let assigned_gradient = pullback(cx, assigned_problem, None, config().linear, &temperature, &weights).unwrap();
            for (a,b) in assigned_gradient.nodal_load.iter().zip(lambda) {
                assert!((a-b).abs() <= 1e-9*b.abs().max(1.0));
            }
            let constant = ElementMaterials::new(MaterialTable::new([
                (MaterialId(7), ConductivityModel::isotropic_declared(12.0).unwrap())]).unwrap(),
                vec![MaterialId(7); mesh.element_count()]).unwrap();
            assert!(!needed(cx, ConductionProblem { element_materials: Some(&constant), ..problem }).unwrap());
        });
    }

    #[test]
    fn nonlinear_native_route_refuses_kinks_and_honors_zero_goals_and_cancellation() {
        with_cx(&CancelGate::new_clock_free(), |cx| {
            let mesh = mesh();
            let source = ScalarField::Uniform(0.0);
            let boundary = ThermalBoundaryBuilder::new(&mesh)
                .remainder("wall", ThermalBc::robin(40.0,300.0).unwrap()).unwrap().finish().unwrap();
            let kink = ConductivityModel::isotropic(ConductivityTable::declared_curve(
                vec![(250.0, 2.0), (300.0, 4.0), (450.0, 22.0)]).unwrap());
            let smooth = material();
            let n = mesh.vertex_count();
            let temperature = vec![300.0; n];
            let weights = vec![0.0; n];
            let problem = ConductionProblem { mesh: &mesh, boundary: &boundary,
                material: &kink, element_materials: None, source: &source };
            assert!(pullback(cx, problem, None, config().linear, &temperature, &weights).is_err());
            let problem = ConductionProblem { material: &smooth, ..problem };
            let gradient = pullback(cx, problem, None, config().linear, &temperature, &weights).unwrap();
            assert_eq!(gradient.iterations, 0);
            assert!(gradient.nodal_load.iter().all(|&v| v == 0.0));
            let gate = CancelGate::new_clock_free(); gate.request();
            with_cx(&gate, |cancelled| {
                assert_eq!(pullback(cancelled, problem, None, config().linear,
                    &temperature, &weights).unwrap_err().code, "cli-solve-cancelled");
            });
        });
    }
}
