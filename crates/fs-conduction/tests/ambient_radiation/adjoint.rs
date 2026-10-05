//! Finite differences of actual non-isothermal, nonlinear-material FEM solves.
use super::*;
use fs_conduction::adjoint::RobinResponse;
use fs_conduction::radiation::{AmbientRadiationGradient, pullback_ambient_radiation};
use fs_conduction::{ConductivityTable, LinearConfig};

#[derive(Clone, Copy)]
struct Controls { density: f64, h: f64, reference: f64, reservoir: f64, epsilon: f64 }
fn nominal() -> Controls {
    Controls { density: 100.0, h: 5.0, reference: 300.0, reservoir: 290.0, epsilon: 0.8 }
}
struct State {
    mesh: ConductionMesh,
    boundary: ThermalBoundary,
    material: ConductivityModel,
    source: ScalarField,
    patch: AmbientRadiationPatch,
    solved: AmbientRadiationSolution,
}
impl State {
    fn problem(&self) -> ConductionProblem<'_> {
        ConductionProblem { mesh: &self.mesh, boundary: &self.boundary,
            material: &self.material, element_materials: None, source: &self.source }
    }
    fn weights(&self) -> Vec<f64> {
        let mut w = vec![0.0; self.mesh.vertex_count()];
        *w.last_mut().unwrap() = 1.0;
        w
    }
    fn goal(&self) -> f64 { *self.solved.conduction.temperature.last().unwrap() }
    fn gradient(&self, cx: &fs_exec::Cx<'_>) -> AmbientRadiationGradient {
        pullback_ambient_radiation(cx, self.problem(), None, linear(),
            &self.solved.conduction.temperature, &["cooler"], &[0.0],
            std::slice::from_ref(&self.patch), &self.weights(), 2*self.mesh.vertex_count()).unwrap()
    }
}
fn linear() -> LinearConfig {
    LinearConfig { tolerance: 1e-9, ..controls().linear }
}
fn physical(cx: &fs_exec::Cx<'_>, c: Controls) -> State {
    let (complex, positions) = box_grid([3, 2, 2], [1.0; 3]);
    let mesh = ConductionMesh::new(complex, positions).unwrap();
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region("hot", |f| on_box_face(f.centroid[0], 0.0), ThermalBc::dirichlet(400.0).unwrap()).unwrap()
        .region("cooler", |f| on_box_face(f.centroid[0], 1.0), ThermalBc::robin(c.h, c.reference).unwrap()).unwrap()
        .adiabatic_remainder().finish().unwrap();
    let material = ConductivityModel::isotropic(ConductivityTable::declared_curve(
        vec![(250.0, 6.0), (500.0, 16.0)]).unwrap());
    let source = ScalarField::Nodal(mesh.positions().iter().map(|p| c.density + 500.0*p[1]).collect());
    let patch = patch(c.epsilon, c.reservoir);
    let mut cfg = controls();
    cfg.stop.residual_rtol = 1e-12;
    let solved = solve_with_ambient_radiation(cx, ConductionProblem { mesh: &mesh,
        boundary: &boundary, material: &material, element_materials: None, source: &source },
        None, std::slice::from_ref(&patch), cfg, AmbientRadiationConfig {
            max_iterations: 150, temperature_tolerance_k: 1e-10,
            balance_tolerance_w: 1e-9, balance_relative_tolerance: 1e-12,
            relaxation: 0.5,
        }).unwrap();
    State { mesh, boundary, material, source, patch, solved }
}
fn close(actual: f64, expected: f64) {
    assert!((actual-expected).abs() < 8e-5*expected.abs().max(1e-4),
        "adjoint {actual:e} versus fresh native difference {expected:e}");
}

#[test]
fn radiation_adjoint_matches_power_convection_reservoir_and_emissivity_resolves() {
    with_cx(|cx| {
        let base = physical(cx, nominal());
        let before = base.solved.clone();
        let gradient = base.gradient(cx);
        assert_eq!(base.solved, before, "derivatives must not mutate or re-solve the primal");
        assert_eq!(gradient.regions, ["cooler"]);
        assert_eq!(gradient.radiation_regions, ["cooler"]);
        let g = &gradient.convection;
        let density = base.mesh.complex().tets.iter().enumerate().map(|(e, tet)|
            base.mesh.element_volume(e)/4.0 * tet.iter()
                .map(|&v| g.nodal_load[v as usize]).sum::<f64>()).sum::<f64>();
        let analytic = [density, g.log_htc[0]/nominal().h, g.references[0],
            gradient.reservoir_temperatures[0], gradient.emissivities[0]];
        for (i, step) in [0.25, 0.001, 0.01, 0.01, 0.0001].into_iter().enumerate() {
            let perturb = |sign: f64| {
                let mut c = nominal();
                let value = match i { 0 => &mut c.density, 1 => &mut c.h,
                    2 => &mut c.reference, 3 => &mut c.reservoir, _ => &mut c.epsilon };
                *value += sign*step;
                physical(cx, c).goal()
            };
            close(analytic[i], (perturb(1.0)-perturb(-1.0))/(2.0*step));
        }
        // Spatial variation is essential: an isothermal fixture would hide
        // replacing the secant trace with a different, uniform-flux model.
        let wall: Vec<_> = base.mesh.boundary().iter()
            .filter(|face| on_box_face(face.centroid[0], 1.0))
            .flat_map(|face| face.vertices.map(|v| base.solved.conduction.temperature[v as usize])).collect();
        assert!(wall.iter().copied().fold(f64::NEG_INFINITY, f64::max)
            - wall.iter().copied().fold(f64::INFINITY, f64::min) > 0.1);
        let frozen = RobinResponse::pullback_mean_htc(cx, ConductionProblem {
            boundary: &base.solved.combined_boundary, ..base.problem() }, None, linear(),
            &base.solved.conduction.temperature, &["cooler"], &[0.0], &base.weights(),
            2*base.mesh.vertex_count()).unwrap();
        let wrong_density = base.mesh.complex().tets.iter().enumerate().map(|(e, tet)|
            base.mesh.element_volume(e)/4.0 * tet.iter()
                .map(|&v| frozen.nodal_load[v as usize]).sum::<f64>()).sum::<f64>();
        assert!((density-wrong_density).abs() > 1e-3*density.abs(), "frozen radiation must lose");
        assert!(g.relative_residual < linear().tolerance);
        for &v in fs_conduction::DofMap::new(&base.boundary, base.mesh.vertex_count()).unwrap().fixed() {
            assert_eq!(g.nodal_load[v], 0.0);
        }
    });
}

#[test]
fn radiation_secant_partials_cover_equal_temperatures_and_hot_reservoirs() {
    for (wall, reservoir) in [(350.0, 280.0), (350.0, 350.0), (350.0, 450.0)] {
        let p = patch(0.8, reservoir);
        let slope = p.secant_partials_w_m2_k2(wall).unwrap();
        let step = 0.001;
        close(slope[0], (p.secant_coefficient_w_m2_k(wall+step).unwrap()
            - p.secant_coefficient_w_m2_k(wall-step).unwrap())/(2.0*step));
        close(slope[1], (patch(0.8, reservoir+step).secant_coefficient_w_m2_k(wall).unwrap()
            - patch(0.8, reservoir-step).secant_coefficient_w_m2_k(wall).unwrap())/(2.0*step));
    }
    assert!(patch(0.8, 300.0).secant_partials_w_m2_k2(250.0).is_err());
    assert!(patch(0.8, 300.0).secant_partials_w_m2_k2(500.0).is_err());
}

#[test]
fn radiation_adjoint_rechecks_the_physical_field_and_respects_limits() {
    with_cx(|cx| {
        let base = physical(cx, nominal());
        let n = base.mesh.vertex_count();
        let w = base.weights();
        let run = |cx: &fs_exec::Cx<'_>, t: &[f64], cfg, entries, weights: &[f64]| {
            pullback_ambient_radiation(cx, base.problem(), None, cfg, t, &["cooler"],
                &[0.0], std::slice::from_ref(&base.patch), weights, entries)
        };
        let mut changed = base.solved.conduction.temperature.clone();
        *changed.last_mut().unwrap() += 1.0;
        assert!(run(cx, &changed, linear(), 2*n, &w).is_err(), "public reports cannot substitute a failed physical field");
        assert!(run(cx, &base.solved.conduction.temperature, linear(), 2*n-1, &w).is_err());
        let cfg = LinearConfig { max_iterations: 1, restart: 1, ..linear() };
        assert!(matches!(run(cx, &base.solved.conduction.temperature, cfg, 2*n, &w),
            Err(ConductionError::LinearSolveFailed { krylov_iterations: 1, .. })));
        let zero = run(cx, &base.solved.conduction.temperature, linear(), 2*n, &vec![0.0; n]).unwrap();
        assert_eq!(zero.convection.iterations, 0);
        assert!(zero.convection.nodal_load.iter().all(|&x| x == 0.0));
        assert_eq!(zero.reservoir_temperatures, [0.0]);
        assert_eq!(zero.emissivities, [0.0]);
        with_cancelled_cx(|cancelled| assert!(matches!(run(cancelled,
            &base.solved.conduction.temperature, linear(), 2*n, &w),
            Err(ConductionError::Cancelled { .. }))));
    });
}

const TRANSFER: [f64; 4] = [0.15, 0.05, 0.30, 0.20];
struct CoupledState {
    mesh: ConductionMesh, boundary: ThermalBoundary, material: ConductivityModel,
    source: ScalarField, patches: Vec<AmbientRadiationPatch>, temperature: Vec<f64>,
}
impl CoupledState {
    fn problem(&self) -> ConductionProblem<'_> {
        ConductionProblem { mesh: &self.mesh, boundary: &self.boundary,
            material: &self.material, element_materials: None, source: &self.source }
    }
    fn weights(&self) -> Vec<f64> {
        let mut weights = vec![0.0; self.temperature.len()];
        *weights.last_mut().unwrap() = 1.0;
        weights
    }
}

// An independently iterated two-trace physical FEM with local gray radiation
// and cross-trace reference feedback. The nonuniform wall, full k(T), and
// prescribed nodes distinguish this from a lumped or frozen-solid comparison.
fn coupled_physical(cx: &fs_exec::Cx<'_>, c: Controls, radiation: bool) -> CoupledState {
    let (complex, positions) = box_grid([3, 2, 2], [1.0; 3]);
    let mesh = ConductionMesh::new(complex, positions).unwrap();
    let material = ConductivityModel::isotropic(ConductivityTable::declared_curve(
        vec![(250.0, 6.0), (500.0, 16.0)]).unwrap());
    let source = ScalarField::Nodal(mesh.positions().iter().map(|p| c.density+500.0*p[1]).collect());
    let patches = if radiation { vec![patch(c.epsilon, c.reservoir),
        AmbientRadiationPatch::new("second", patch(0.5,430.0).emissivity().clone(),430.0).unwrap()]
    } else { Vec::new() };
    let mut means = [350.0; 2];
    for _ in 0..160 {
        let references = [c.reference+TRANSFER[0]*(means[0]-300.0)+TRANSFER[1]*(means[1]-300.0),
            305.0+TRANSFER[2]*(means[0]-300.0)+TRANSFER[3]*(means[1]-300.0)];
        let h = [c.h, 7.0];
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .region("hot", |f| on_box_face(f.centroid[0],0.0), ThermalBc::dirichlet(400.0).unwrap()).unwrap()
            .region("cooler", |f| on_box_face(f.centroid[0],1.0) && f.centroid[1]<0.5,
                ThermalBc::robin(h[0],references[0]).unwrap()).unwrap()
            .region("second", |f| on_box_face(f.centroid[0],1.0) && f.centroid[1]>=0.5,
                ThermalBc::robin(h[1],references[1]).unwrap()).unwrap()
            .adiabatic_remainder().finish().unwrap();
        let mut replacements = Vec::new();
        for (i, patch) in patches.iter().enumerate() {
            let hr = patch.secant_coefficient_w_m2_k(means[i]).unwrap();
            let total = h[i]+hr;
            let reference = h[i]/total*references[i]+hr/total*patch.ambient_temperature_k();
            let index = boundary.region_names().iter().position(|r| r==patch.region()).unwrap();
            replacements.push((index,total,reference));
        }
        let combined = boundary.with_uniform_robin_replacements(&replacements).unwrap();
        let mut cfg = controls(); cfg.stop.residual_rtol = 1e-12;
        let solved = fs_conduction::solve(cx, ConductionProblem { mesh: &mesh, boundary: &combined,
            material: &material, element_materials: None, source: &source },cfg).unwrap();
        let mut next = [0.0; 2];
        for (i, name) in ["cooler","second"].into_iter().enumerate() {
            next[i] = solved.report.robin_fluxes.iter().find(|r| r.region==name).unwrap().mean_wall_temperature_k;
        }
        if (next[0]-means[0]).abs().max((next[1]-means[1]).abs()) < 1e-10 {
            return CoupledState { mesh, boundary, material, source, patches, temperature: solved.temperature };
        }
        for i in 0..2 { means[i] = 0.5*(means[i]+next[i]); }
    }
    panic!("coupled physical fixture did not converge");
}

#[test]
fn coupled_radiation_adjoint_keeps_cross_reference_and_material_feedback() {
    with_cx(|cx| {
        for radiation in [false,true] {
            let base = coupled_physical(cx, nominal(), radiation);
            let before = base.temperature.clone();
            let gradient = AmbientRadiationGradient::pullback_with_reference_feedback(cx,
                base.problem(), None, linear(), &base.temperature, &["cooler","second"],
                &[0.0;2], &TRANSFER, &base.patches, &base.weights(), 8*base.mesh.vertex_count()).unwrap();
            assert_eq!(base.temperature,before);
            let g = &gradient.convection;
            let density = base.mesh.complex().tets.iter().enumerate().map(|(e,tet)|
                base.mesh.element_volume(e)/4.0*tet.iter().map(|&v| g.nodal_load[v as usize]).sum::<f64>()).sum::<f64>();
            let mut analytic = vec![density,g.log_htc[0]/nominal().h,g.references[0]];
            if radiation { analytic.extend([gradient.reservoir_temperatures[0],gradient.emissivities[0]]); }
            for (i, (&actual, step)) in analytic.iter().zip([0.25,0.001,0.01,0.01,0.0001]).enumerate() {
                let perturb = |sign: f64| {
                    let mut c = nominal();
                    let value = match i { 0=>&mut c.density,1=>&mut c.h,2=>&mut c.reference,
                        3=>&mut c.reservoir,_=>&mut c.epsilon };
                    *value += sign*step;
                    *coupled_physical(cx,c,radiation).temperature.last().unwrap()
                };
                close(actual,(perturb(1.0)-perturb(-1.0))/(2.0*step));
            }
            let frozen = AmbientRadiationGradient::pullback_with_reference_feedback(cx,
                base.problem(), None, linear(), &base.temperature, &["cooler","second"],
                &[0.0;2], &[0.0;4], &base.patches, &base.weights(), 8*base.mesh.vertex_count()).unwrap();
            assert!((g.references[0]-frozen.convection.references[0]).abs() > 1e-3*g.references[0].abs(),
                "freezing cross-reference feedback must give a detectably wrong derivative");
            assert!(AmbientRadiationGradient::pullback_with_reference_feedback(cx,
                base.problem(), None, linear(), &base.temperature, &["cooler","second"],
                &[0.0;2], &TRANSFER[..3], &base.patches, &base.weights(), 8*base.mesh.vertex_count()).is_err());
        }
    });
}
