use super::*;
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_conduction::{ConductivityModel, ConductivityTable, LinearConfig, ScalarField,
    ThermalBc, ThermalBoundaryBuilder};
use fs_conduction::fixtures::{box_grid, on_box_face};
use fs_conduction::transient::VolumetricHeatCapacity;
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use crate::transient::variational::{WeakConstraintWindow, WindowControl, WindowObjective};
use crate::transient::variational::joint::{JointWindow, JointWindowStudy, ParameterFamily};
use crate::transient::variational::study::StudySettings;
use crate::StopReason;
use std::cell::Cell;

fn with_cx(f: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(&gate, arena, StreamKey { seed: 61, kernel_id: 820, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic);
        f(&cx, &gate);
    });
}
fn config() -> ConductionWindowConfig {
    ConductionWindowConfig { step: StepConfig { linear: LinearConfig {
        tolerance: 1e-11, max_iterations: 2000, restart: 20 }, energy_tolerance_j: 1e-8 },
        nonlinear: None, max_vertices: 100, max_elements: 100, max_intervals: 8, max_parameters: 1 }
}
struct Fixture { mesh: ConductionMesh, boundary: ThermalBoundary, material: ConductivityModel }
impl Fixture {
    fn new(pinned: bool, nonlinear: bool) -> Self {
        let (complex, positions) = box_grid([2,1,1], [0.1,0.04,0.03]);
        let mesh = ConductionMesh::new(complex, positions).unwrap();
        let builder = ThermalBoundaryBuilder::new(&mesh);
        let builder = if pinned { builder.region("fixed", |face| on_box_face(face.centroid[0],0.0),
            ThermalBc::dirichlet(300.0).unwrap()).unwrap() } else { builder };
        let boundary = builder.adiabatic_remainder().finish().unwrap();
        let material = if nonlinear { ConductivityModel::isotropic(
            ConductivityTable::declared_curve(vec![(250.0,1.0),(450.0,21.0)]).unwrap())
        } else { ConductivityModel::isotropic_declared(2.0).unwrap() };
        Self { mesh, boundary, material }
    }
}
struct Thermal<'a> {
    engine: &'a BackwardEuler<'a>, fixture: &'a Fixture, source: ScalarField,
    missing_partial: Cell<bool>, targets: Vec<f64>,
}
impl ConductionWindowModel for Thermal<'_> {
    fn engine(&self) -> &BackwardEuler<'_> { self.engine }
    fn problem(&self, _: usize) -> Result<ConductionProblem<'_>, ConductionError> {
        Ok(ConductionProblem { mesh: &self.fixture.mesh, boundary: &self.fixture.boundary,
            material: &self.fixture.material, element_materials: None, source: &self.source })
    }
    fn parameter_count(&self) -> usize { 1 }
    fn parameter_pullback(&self, _: usize, cx: &Cx<'_>, step: &StepLinearization<'_>,
        lambda: &[f64], parameters: &mut [f64]) -> Result<(), ConductionError>
    {
        if !self.missing_partial.get() {
            // Uniform source density, INCLUDING the zero-source parameter point.
            parameters[0] = step.source_density_pullback(cx,lambda)?.iter().sum();
        }
        Ok(())
    }
}
fn model<'a>(fixture: &'a Fixture, engine: &'a BackwardEuler<'a>, source: f64) -> Thermal<'a> {
    Thermal { engine, fixture, source: ScalarField::Uniform(source), missing_partial: Cell::new(false), targets: vec![] }
}
fn dot(a: &[f64], b: &[f64]) -> f64 { a.iter().zip(b).map(|(a,b)| a*b).sum() }
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a-b).abs() < tolerance*b.abs().max(1.0), "{a:e} != {b:e}");
}

#[test]
fn native_endpoint_and_history_source_cotangents_match_uniform_heating() {
    with_cx(|cx,_| {
        let f = Fixture::new(false,false);
        let engine = BackwardEuler::uniform(cx,&f.mesh,VolumetricHeatCapacity::declared(5.0).unwrap()).unwrap();
        let m = model(&f,&engine,10.0);
        let policy = ConductionWindowPolicy::new(cx,&f.mesh,&f.boundary,&[0.0,2.0],config(),&mut||false).unwrap();
        let old = vec![300.0;policy.dimension()];
        let tape = policy.record(&m,0,0.0,2.0,&old,&mut||false).unwrap();
        let direct = engine.linearize_step(cx,m.problem(0).unwrap(),None,&old,2.0,config().step,None,&[]).unwrap();
        assert_eq!(tape.endpoint(),direct.primal().temperature);
        for &t in tape.endpoint() { close(t,304.0,1e-11); }
        let seed = vec![1.0/old.len() as f64;old.len()];
        let bar = tape.pullback(&seed,&[0.7],&mut||false).unwrap();
        close(bar.initial.iter().sum(),1.0,1e-10);
        close(bar.parameters[0],1.1,1e-10); // dt/c + direct = 2/5 + 0.7
        assert_eq!(bar.replayed_steps,0);
        let zero = model(&f,&engine,0.0);
        close(policy.record(&zero,0,0.0,2.0,&old,&mut||false).unwrap()
            .pullback(&seed,&[0.0],&mut||false).unwrap().parameters[0],0.4,1e-10);
    });
}

#[test]
fn prescribed_nodes_are_not_controls_and_nonlinear_pullback_is_correct() {
    with_cx(|cx,_| {
        let f = Fixture::new(true,true);
        let engine = BackwardEuler::uniform(cx,&f.mesh,VolumetricHeatCapacity::declared(1000.0).unwrap()).unwrap();
        let m = model(&f,&engine,2000.0);
        let mut cfg = config(); cfg.nonlinear = Some(NonlinearStepConfig::default());
        let policy = ConductionWindowPolicy::new(cx,&f.mesh,&f.boundary,&[0.0,0.3],cfg,&mut||false).unwrap();
        assert!(policy.dimension()<f.mesh.vertex_count());
        let old = vec![305.0;policy.dimension()];
        let expanded = policy.expand_field(&old).unwrap();
        for &(v,t) in f.boundary.dirichlet() { assert_eq!(expanded[v],t); assert_eq!(policy.slot_of(v),None); }
        assert_eq!(policy.gather_field(&expanded).unwrap(),old);
        let mut invalid = expanded.clone(); invalid[f.boundary.dirichlet()[0].0] += 1.0;
        assert!(policy.gather_field(&invalid).is_err());
        let tape = policy.record(&m,0,0.0,0.3,&old,&mut||false).unwrap();
        let seed = (0..old.len()).map(|i|(i+1) as f64/old.len() as f64).collect::<Vec<_>>();
        let bar = tape.pullback(&seed,&[0.0],&mut||false).unwrap();
        for i in 0..old.len() {
            let (mut plus,mut minus) = (old.clone(),old.clone()); plus[i]+=1e-4; minus[i]-=1e-4;
            let a = policy.record(&m,0,0.0,0.3,&plus,&mut||false).unwrap();
            let b = policy.record(&m,0,0.0,0.3,&minus,&mut||false).unwrap();
            close(bar.initial[i],(dot(&seed,a.endpoint())-dot(&seed,b.endpoint()))/2e-4,3e-5);
        }
        let plus = model(&f,&engine,2000.01); let minus = model(&f,&engine,1999.99);
        let a = policy.record(&plus,0,0.0,0.3,&old,&mut||false).unwrap();
        let b = policy.record(&minus,0,0.0,0.3,&old,&mut||false).unwrap();
        close(bar.parameters[0],(dot(&seed,a.endpoint())-dot(&seed,b.endpoint()))/0.02,2e-7);
    });
}

#[test]
fn mesh_and_boundary_identity_are_checked_before_physical_work() {
    with_cx(|cx,_| {
        let f = Fixture::new(true,false); let other = Fixture::new(true,false);
        let engine = BackwardEuler::uniform(cx,&other.mesh,VolumetricHeatCapacity::declared(1000.0).unwrap()).unwrap();
        let m = model(&other,&engine,0.0);
        let p = ConductionWindowPolicy::new(cx,&f.mesh,&f.boundary,&[0.0,1.0],config(),&mut||false).unwrap();
        assert!(p.record(&m,0,0.0,1.0,&vec![300.0;p.dimension()],&mut||false).is_err());
        let changed = ThermalBoundaryBuilder::new(&other.mesh)
            .region("fixed",|face|on_box_face(face.centroid[0],0.0),ThermalBc::dirichlet(301.0).unwrap()).unwrap()
            .adiabatic_remainder().finish().unwrap();
        let p = ConductionWindowPolicy::new(cx,&other.mesh,&changed,&[0.0,1.0],config(),&mut||false).unwrap();
        assert!(p.record(&m,0,0.0,1.0,&vec![300.0;p.dimension()],&mut||false).is_err());
        assert!(p.record(&m,0,0.0,0.9,&vec![300.0;p.dimension()],&mut||false).is_err());
        let mut small = config(); small.max_vertices=1;
        assert!(ConductionWindowPolicy::new(cx,&f.mesh,&f.boundary,&[0.0,1.0],small,&mut||false).is_err());
    });
}

#[test]
fn missing_derivatives_and_both_cancellation_sources_publish_no_gradient() {
    with_cx(|cx,gate| {
        let f = Fixture::new(false,false);
        let engine = BackwardEuler::uniform(cx,&f.mesh,VolumetricHeatCapacity::declared(5.0).unwrap()).unwrap();
        let m = model(&f,&engine,10.0);
        let p = ConductionWindowPolicy::new(cx,&f.mesh,&f.boundary,&[0.0,1.0],config(),&mut||false).unwrap();
        let old = vec![300.0;p.dimension()];
        assert!(matches!(p.record(&m,0,0.0,1.0,&old,&mut||true),Err(WindowError::Cancelled)));
        let tape = p.record(&m,0,0.0,1.0,&old,&mut||false).unwrap();
        let seed = vec![1.0/old.len() as f64;old.len()];
        m.missing_partial.set(true);
        assert!(matches!(tape.pullback(&seed,&[0.0],&mut||false),Err(WindowError::NonFinite(_))));
        m.missing_partial.set(false);
        assert!(tape.pullback(&seed,&[0.0],&mut||false).is_ok());
        gate.request();
        assert!(matches!(tape.pullback(&seed,&[0.0],&mut||false),Err(WindowError::Cancelled)));
        assert!(matches!(p.record(&m,0,0.0,1.0,&old,&mut||false),Err(WindowError::Cancelled)));
    });
}

impl WindowObjective for Thermal<'_> {
    fn evaluate(&self,_:&[f64],_:usize,x:&[f64],out:&mut[f64],_:&mut dyn FnMut()->bool)->Result<f64,String> {
        let mut cost=0.0;
        for ((x,y),g) in x.iter().zip(&self.targets).zip(out) {let r=(x-y)/0.05; *g=r/0.05; cost+=0.5*r*r;}
        Ok(cost)
    }
}
struct Family<'a> { engine: &'a BackwardEuler<'a>, fixture: &'a Fixture, target: Vec<f64> }
impl<'a> ParameterFamily for Family<'a> {
    type Model = Thermal<'a>;
    fn instantiate(&self,p:&[f64],_:&mut dyn FnMut()->bool)->Result<Self::Model,String> {
        let mut m=model(self.fixture,self.engine,p[0]); m.targets=self.target.clone(); Ok(m)
    }
}
// Independent four-variable normal equations for a uniform-field solution.
fn uniform_oracle(nodes:usize)->Vec<f64> {
    let n=nodes as f64; let mut a=vec![vec![0.0;4];4]; let mut b=vec![0.0;4];
    let mut term=|row:[f64;4],target:f64,weight:f64| {
        for i in 0..4 {b[i]+=weight*row[i]*target;for j in 0..4 {a[i][j]+=weight*row[i]*row[j];}}
    };
    for k in 0..3 {let mut row=[0.0;4];row[k]=1.0;term(row,0.5*k as f64,n/0.05_f64.powi(2));}
    term([1.0,0.0,0.0,0.0],0.0,n/4.0);
    term([-1.0,1.0,0.0,-0.05],0.0,n/0.2_f64.powi(2));
    term([0.0,-1.0,1.0,-0.05],0.0,n/0.2_f64.powi(2));
    term([0.0,0.0,0.0,1.0],0.0,0.01);
    for i in 0..4 {let pivot=a[i][i];for j in i..4 {a[i][j]/=pivot;} b[i]/=pivot;
        for k in 0..4 {if k!=i {let q=a[k][i];for j in i..4 {a[k][j]-=q*a[i][j];} b[k]-=q*b[i];}}}
    b
}
#[test]
fn joint_source_and_spatial_states_match_independent_map_and_resume() {
    with_cx(|cx,_| {
        let f=Fixture::new(false,false); let n=f.mesh.vertex_count();
        let engine=BackwardEuler::uniform(cx,&f.mesh,VolumetricHeatCapacity::declared(5.0).unwrap()).unwrap();
        let times=[0.0,0.25,0.5];
        let p=ConductionWindowPolicy::new(cx,&f.mesh,&f.boundary,&times,config(),&mut||false).unwrap();
        let w=WeakConstraintWindow::new(&times,&vec![300.0;3*n],&vec![0.1;n],&vec![2.0;n],&vec![0.2;2*n],1000).unwrap();
        let family=Family{engine:&engine,fixture:&f,target:(0..3).flat_map(|k|vec![300.0+0.5*k as f64;n]).collect()};
        let joint=JointWindow::new(&w,&family,&[0.0],&[1.0],&[10.0],3).unwrap();
        let settings=StudySettings{memory:12,gradient_tolerance:2e-6,max_evaluations:2000,max_optimizer_components:100000};
        let mut control=WindowControl::new(5000,10000,100000);
        let mut straight=JointWindowStudy::new(&joint,&vec![0.0;joint.control_dimension()],p,settings,&mut control,&mut||false).unwrap();
        let mut split=straight.clone();
        let stop=straight.run(400,&mut control,&mut||false).unwrap();
        assert_eq!(stop.reason,StopReason::GradNorm,"{stop:?}");
        for _ in 0..400 {if split.run(1,&mut control,&mut||false).unwrap().reason!=StopReason::IterationCap {break;}}
        assert_eq!(straight.accepted(),split.accepted());
        let oracle=uniform_oracle(n); let got=straight.accepted();
        for k in 0..3 {for &x in &got.window.states[k*n..(k+1)*n] {close(x-300.0,oracle[k],2e-5);}}
        close(got.parameters[0],oracle[3],2e-5);
        assert_eq!(got.controls,straight.optimizer().x);
    });
}
