//! The source-density pullback must transpose the CONSISTENT load integration,
//! not multiply lambda by lumped volumes or divide by the current watts.
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_conduction::{ConductionMesh, ConductionProblem, ConductivityModel,
    ComponentPower, PowerMap, PowerUncertainty, LinearConfig, ScalarField,
    ThermalBc, ThermalBoundary, ThermalBoundaryBuilder};
use fs_conduction::transient::VolumetricHeatCapacity;
use fs_conduction::transient::backward_euler::{BackwardEuler, StepConfig};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_rep_mesh::TetComplex;

fn context<R>(gate: &CancelGate, run: impl FnOnce(&Cx<'_>) -> R) -> R {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| run(&Cx::new(
        gate, arena, StreamKey {seed: 19, kernel_id: 823, tile: 0, iteration: 0},
        Budget::INFINITE, ExecMode::Deterministic)))
}
fn mesh() -> ConductionMesh {
    ConductionMesh::new(TetComplex::from_tets(4, vec![[0,1,2,3]]),
        vec![[0.0,0.0,0.0],[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]]).unwrap()
}
fn boundary(mesh: &ConductionMesh, fixed: bool) -> ThermalBoundary {
    let builder = ThermalBoundaryBuilder::new(mesh);
    let builder = if fixed {
        builder.region("base", |face| face.centroid[2] == 0.0,
            ThermalBc::dirichlet(300.0).unwrap()).unwrap()
    } else { builder };
    builder.adiabatic_remainder().finish().unwrap()
}
fn config() -> StepConfig {
    StepConfig {linear: LinearConfig {tolerance: 1e-12, max_iterations: 2000, restart: 20},
        energy_tolerance_j: 1e-8}
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a-b).abs() < tolerance, "{a:e} versus {b:e}");
}
fn map(power: f64) -> PowerMap {
    PowerMap::new(vec![ComponentPower::new("dormant-chip", power,
        PowerUncertainty::Unstated, vec![0]).unwrap()], power).unwrap()
}

#[test]
fn consistent_source_mass_preserves_mean_zero_modes_and_has_no_extra_dt() {
    context(&CancelGate::new_clock_free(), |cx| {
        let mesh = mesh(); let boundary = boundary(&mesh, false);
        let material = ConductivityModel::isotropic_declared(2.0).unwrap();
        let source = ScalarField::Nodal(vec![2.0,-3.0,5.0,1.0]);
        let engine = BackwardEuler::uniform(cx, &mesh,
            VolumetricHeatCapacity::declared(5.0).unwrap()).unwrap();
        let lambda = [1.0,-1.0,0.0,0.0];
        for dt in [0.2,2.0] {
            let step = engine.linearize_step(cx, ConductionProblem {mesh: &mesh,
                boundary: &boundary, material: &material, element_materials: None, source: &source},
                None, &[300.0;4], dt, config(), None, &[]).unwrap();
            let g = step.source_density_pullback(cx, &lambda).unwrap();
            for (actual, expected) in g.iter().zip([1.0/120.0,-1.0/120.0,0.0,0.0]) {
                near(*actual, expected, 1e-14);
            }
            let contraction: f64 = g.iter().enumerate().map(|(i,g)|g*source.at(i)).sum();
            near(contraction, step.source_multiplier_pullback(cx,&lambda).unwrap(), 1e-14);
        }
    });
}

#[test]
fn a_zero_watt_component_has_a_nonzero_absolute_power_derivative() {
    context(&CancelGate::new_clock_free(), |cx| {
        let mesh = mesh(); let boundary = boundary(&mesh, false);
        let material = ConductivityModel::isotropic_declared(2.0).unwrap();
        let engine = BackwardEuler::uniform(cx, &mesh,
            VolumetricHeatCapacity::declared(5.0).unwrap()).unwrap();
        let (source,audit) = map(0.0).volumetric_source(&mesh,0.0).unwrap();
        let step = engine.linearize_step(cx, ConductionProblem {mesh: &mesh,
            boundary: &boundary, material: &material, element_materials: None, source: &source},
            None, &[300.0;4], 2.0, config(), None, &[]).unwrap();
        let g = step.pullback(cx, &[0.25;4], &[], &[]).unwrap();
        let density = step.source_density_pullback(cx,&g.nodal_load).unwrap();
        let derivative = density[0] / audit.rows()[0].bound_volume_m3();
        // Insulated tet: average temperature rise = P*dt/(c*V).
        near(derivative, 2.0/(5.0/6.0), 1e-10);
        near(step.source_multiplier_pullback(cx,&g.nodal_load).unwrap(),0.0,1e-14);
        let h = 1e-3;
        let (perturbed,_) = map(h).volumetric_source(&mesh,0.0).unwrap();
        let plus = engine.advance(cx, ConductionProblem {mesh: &mesh,
            boundary: &boundary, material: &material, element_materials: None, source: &perturbed},
            None, &[300.0;4], 2.0, config()).unwrap();
        let nominal: f64 = step.primal().temperature.iter().sum::<f64>()/4.0;
        near((plus.temperature.iter().sum::<f64>()/4.0-nominal)/h,derivative,1e-6);
    });
}

#[test]
fn fixed_equations_do_not_remove_source_support_on_fixed_vertices() {
    context(&CancelGate::new_clock_free(), |cx| {
        let mesh = mesh(); let boundary = boundary(&mesh, true);
        let material = ConductivityModel::isotropic_declared(2.0).unwrap();
        let source = ScalarField::Uniform(0.0);
        let engine = BackwardEuler::uniform(cx, &mesh,
            VolumetricHeatCapacity::declared(5.0).unwrap()).unwrap();
        let step = engine.linearize_step(cx, ConductionProblem {mesh: &mesh,
            boundary: &boundary, material: &material, element_materials: None, source: &source},
            None, &[300.0;4], 2.0, config(), None, &[]).unwrap();
        let actual = step.source_density_pullback(cx,&[7.0,8.0,9.0,1.0]).unwrap();
        assert_eq!(actual, step.source_density_pullback(cx,&[0.0,0.0,0.0,1.0]).unwrap());
        for (a,e) in actual.iter().zip([1.0/120.0,1.0/120.0,1.0/120.0,1.0/60.0]) {
            near(*a,e,1e-14);
        }
        assert!(step.source_density_pullback(cx,&[1.0]).is_err());
        assert!(step.source_density_pullback(cx,&[f64::NAN;4]).is_err());
        let cancelled = CancelGate::new_clock_free(); cancelled.request();
        context(&cancelled, |stopped| {
            assert!(step.source_density_pullback(stopped,&[0.0;4]).is_err());
        });
    });
}
