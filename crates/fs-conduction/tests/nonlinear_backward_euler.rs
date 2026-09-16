//! Independent one-tetrahedron P1 oracle for the nonlinear timestep.
//! The manufactured load is derived from exact element mass/stiffness formulas,
//! not by calling production assembly on a desired answer.
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_conduction::{ConductionError, ConductionMesh, ConductionProblem, ConductivityModel,
    ElementMaterials, LinearConfig, MaterialId, MaterialTable, ScalarField,
    ThermalBc, ThermalBoundaryBuilder};
use fs_conduction::material::ConductivityTable;
use fs_conduction::transient::VolumetricHeatCapacity;
use fs_conduction::transient::backward_euler::{BackwardEuler, NonlinearStepConfig, StepConfig};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_rep_mesh::TetComplex;

fn with_cx<T>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> T) -> T {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| f(&Cx::new(gate, arena,
        StreamKey { seed: 713, kernel_id: 723, tile: 0, iteration: 0 },
        Budget::INFINITE, ExecMode::Deterministic)))
}
fn mesh() -> ConductionMesh {
    ConductionMesh::new(TetComplex::from_tets(4, vec![[0,1,2,3]]),
        vec![[0.,0.,0.], [1.,0.,0.], [0.,1.,0.], [0.,0.,1.]]).unwrap()
}
fn curve() -> ConductivityModel {
    ConductivityModel::isotropic(ConductivityTable::declared_curve(
        vec![(280.0,0.6),(340.0,1.8)]).unwrap())
}
fn config() -> StepConfig {
    StepConfig { linear: LinearConfig { tolerance: 1e-11, max_iterations: 128, restart: 4 },
        energy_tolerance_j: 1e-9 }
}
fn close(a: f64, b: f64, tol: f64) { assert!((a-b).abs() <= tol, "{a:e} != {b:e}"); }
fn manufactured_source() -> ScalarField {
    let old = [299.0,299.5,300.2,300.4];
    let target = [300.0,301.0,302.0,303.0];
    // Unit tetra: V=1/6, M=V/20*(I+11^T), grad(T)=[1,2,3].
    // K*T=V*k(Tmean)*[-6,1,2,3], with k(301.5)=1.03.
    // M^-1 K*T=20*k*[-6,1,2,3] (the conductive load sums to zero).
    ScalarField::Nodal((0..4).map(|i|
        100.0/0.4*(target[i]-old[i])+20.0*1.03*[-6.,1.,2.,3.][i]).collect())
}

#[test]
fn known_nonuniform_endpoint_uses_new_temperature_conductivity_and_replays() {
    let mesh = mesh();
    let boundary = ThermalBoundaryBuilder::new(&mesh).adiabatic_remainder().finish().unwrap();
    let material = curve();
    let source = manufactured_source();
    let old = [299.0,299.5,300.2,300.4];
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let engine = BackwardEuler::uniform(cx,&mesh,VolumetricHeatCapacity::declared(100.0).unwrap()).unwrap();
        let p = ConductionProblem { mesh:&mesh, boundary:&boundary, material:&material,
            element_materials:None, source:&source };
        let a = engine.advance_nonlinear(cx,p,None,&old,0.4,config(),NonlinearStepConfig::default()).unwrap();
        let b = engine.advance_nonlinear(cx,p,None,&old,0.4,config(),NonlinearStepConfig::default()).unwrap();
        assert_eq!(a.step.temperature,b.step.temperature);
        assert_eq!(old,[299.0,299.5,300.2,300.4]);
        for (&t, exact) in a.step.temperature.iter().zip([300.,301.,302.,303.]) { close(t,exact,1e-8); }
        assert!(a.nonlinear_iterations > 1);
        assert!(a.residual_j <= a.threshold_j);
        assert!(a.step.krylov_iterations <= config().linear.max_iterations);
        close(a.step.stored_energy_change_j,0.4*a.step.source_w,1e-9);
        // This deliberately wrong model is a falsifier, not another oracle.
        let frozen = ConductivityModel::isotropic_declared(1.0+0.02*(old.iter().sum::<f64>()/4.0-300.0)).unwrap();
        let wrong = engine.advance(cx,ConductionProblem {material:&frozen,..p},None,&old,0.4,config()).unwrap();
        assert!(wrong.temperature.iter().zip(&a.step.temperature).any(|(a,b)|(a-b).abs()>1e-3));
    });
}

#[test]
fn constant_material_matches_the_existing_linear_robin_step() {
    let mesh = mesh();
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .remainder("air",ThermalBc::robin(2.0,300.0).unwrap()).unwrap().finish().unwrap();
    let material = ConductivityModel::isotropic_declared(3.0).unwrap();
    let source = ScalarField::Uniform(120.0);
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let engine = BackwardEuler::uniform(cx,&mesh,VolumetricHeatCapacity::declared(100.0).unwrap()).unwrap();
        let p = ConductionProblem {mesh:&mesh,boundary:&boundary,material:&material,source:&source,element_materials:None};
        let old = [330.,310.,320.,300.];
        let linear = engine.advance(cx,p,None,&old,0.2,config()).unwrap();
        let nonlinear = engine.advance_nonlinear(cx,p,None,&old,0.2,config(),NonlinearStepConfig::default()).unwrap();
        for (a,b) in linear.temperature.iter().zip(&nonlinear.step.temperature) {close(*a,*b,1e-8);}
        close(linear.robin_out_w,nonlinear.step.robin_out_w,1e-7);
        assert_eq!(nonlinear.nonlinear_iterations,1);
    });
}

#[test]
fn prescribed_temperature_reaction_includes_consistent_storage() {
    let mesh = mesh();
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region("fixed",|face|face.centroid[0]==0.0,ThermalBc::dirichlet(300.0).unwrap()).unwrap()
        .adiabatic_remainder().finish().unwrap();
    let material = curve(); let source = ScalarField::Uniform(100.0);
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let engine = BackwardEuler::uniform(cx,&mesh,VolumetricHeatCapacity::declared(10.0).unwrap()).unwrap();
        let p = ConductionProblem {mesh:&mesh,boundary:&boundary,material:&material,source:&source,element_materials:None};
        let result = engine.advance_nonlinear(cx,p,None,&[300.;4],0.2,config(),NonlinearStepConfig::default()).unwrap().step;
        for &(v,t) in boundary.dirichlet() {assert_eq!(result.temperature[v],t);}
        close(result.stored_energy_change_j,0.2*(result.source_w+result.dirichlet_in_w),1e-9);
        let mut bad = [300.;4]; bad[boundary.dirichlet()[0].0]=301.;
        assert!(engine.advance_nonlinear(cx,p,None,&bad,0.2,config(),NonlinearStepConfig::default()).is_err());
    });
}

#[test]
fn heterogeneous_laws_and_capacities_do_not_collapse_to_the_fallback() {
    let points = vec![[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.],
        [2.,0.,0.],[3.,0.,0.],[2.,1.,0.],[2.,0.,1.]];
    let mesh = ConductionMesh::new(TetComplex::from_tets(8,vec![[0,1,2,3],[4,5,6,7]]),points).unwrap();
    let boundary = ThermalBoundaryBuilder::new(&mesh).adiabatic_remainder().finish().unwrap();
    // An unusable fallback must not replace assigned material laws.
    let fallback = ConductivityModel::isotropic(ConductivityTable::declared_curve(vec![(1.,1.),(2.,2.)]).unwrap());
    let materials = ElementMaterials::new(MaterialTable::new([
        (MaterialId(0),curve()),(MaterialId(1),ConductivityModel::isotropic_declared(2.0).unwrap()),
    ]).unwrap(),vec![MaterialId(0),MaterialId(1)]).unwrap();
    let source = ScalarField::Uniform(2.0);
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let capacities = [10.,20.].map(|v|VolumetricHeatCapacity::declared(v).unwrap());
        let engine = BackwardEuler::per_element(cx,&mesh,&capacities).unwrap();
        let p = ConductionProblem {mesh:&mesh,boundary:&boundary,material:&fallback,source:&source,element_materials:Some(&materials)};
        let result = engine.advance_nonlinear(cx,p,None,&[300.;8],0.5,config(),NonlinearStepConfig::default()).unwrap();
        for &t in &result.step.temperature[..4] {close(t,300.1,1e-9);}
        for &t in &result.step.temperature[4..] {close(t,300.05,1e-9);}
        close(result.step.stored_energy_change_j,1.0/3.0,1e-9);
    });
}

#[test]
fn a_loose_newton_tolerance_cannot_bypass_energy_closure() {
    let mesh = mesh();
    let boundary = ThermalBoundaryBuilder::new(&mesh).adiabatic_remainder().finish().unwrap();
    let material = curve(); let source = ScalarField::Uniform(1.0);
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let engine = BackwardEuler::uniform(cx,&mesh,VolumetricHeatCapacity::declared(100.0).unwrap()).unwrap();
        let p = ConductionProblem {mesh:&mesh,boundary:&boundary,material:&material,source:&source,element_materials:None};
        let lax = NonlinearStepConfig {residual_atol_j:1e6,..NonlinearStepConfig::default()};
        let error = engine.advance_nonlinear(cx,p,None,&[300.;4],1.,config(),lax).unwrap_err();
        assert!(error.to_string().contains("energy residual"));
    });
}

#[test]
fn exhausted_invalid_and_cancelled_trials_leave_history_untouched() {
    let mesh = mesh();
    let boundary = ThermalBoundaryBuilder::new(&mesh).adiabatic_remainder().finish().unwrap();
    let material = curve(); let source = manufactured_source(); let old = [299.,299.5,300.2,300.4];
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let engine = BackwardEuler::uniform(cx,&mesh,VolumetricHeatCapacity::declared(100.0).unwrap()).unwrap();
        let p = ConductionProblem {mesh:&mesh,boundary:&boundary,material:&material,source:&source,element_materials:None};
        let one = NonlinearStepConfig {max_iterations:1,..NonlinearStepConfig::default()};
        assert!(matches!(engine.advance_nonlinear(cx,p,None,&old,0.4,config(),one),Err(ConductionError::NotConverged{..})));
        let mut short = config(); short.linear.max_iterations=1;
        assert!(engine.advance_nonlinear(cx,p,None,&old,0.4,short,NonlinearStepConfig::default()).is_err());
        short=config(); short.linear.restart=0;
        assert!(engine.advance_nonlinear(cx,p,None,&old,0.4,short,NonlinearStepConfig::default()).is_err());
        for dt in [0.,-1.,f64::NAN,f64::INFINITY] {
            assert!(engine.advance_nonlinear(cx,p,None,&old,dt,config(),NonlinearStepConfig::default()).is_err());
        }
        gate.request();
        assert!(matches!(engine.advance_nonlinear(cx,p,None,&old,0.4,config(),NonlinearStepConfig::default()),Err(ConductionError::Cancelled{..})));
        assert_eq!(old,[299.,299.5,300.2,300.4]);
    });
}

#[test]
fn no_material_extrapolation_when_the_endpoint_cannot_fit_its_declared_span() {
    let mesh = mesh();
    let boundary = ThermalBoundaryBuilder::new(&mesh).adiabatic_remainder().finish().unwrap();
    let material = ConductivityModel::isotropic(ConductivityTable::declared_curve(vec![(299.,1.),(301.,2.)]).unwrap());
    let source = ScalarField::Uniform(100.0);
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let engine = BackwardEuler::uniform(cx,&mesh,VolumetricHeatCapacity::declared(100.0).unwrap()).unwrap();
        let p = ConductionProblem {mesh:&mesh,boundary:&boundary,material:&material,source:&source,element_materials:None};
        // Exact mean endpoint would be 305 K; no admissible solution exists.
        assert!(engine.advance_nonlinear(cx,p,None,&[300.;4],5.,config(),NonlinearStepConfig::default()).is_err());
        assert!(matches!(engine.advance_nonlinear(cx,p,None,&[302.;4],0.1,config(),NonlinearStepConfig::default()),
            Err(ConductionError::OutsideTemperatureSpan{..})));
    });
}
