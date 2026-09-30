use super::*;
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_conduction::{ConductionError, ConductionMesh, ConductionProblem, ConductivityModel,
    ConductivityTable, LinearConfig, ScalarField, ThermalBc, ThermalBoundary, ThermalBoundaryBuilder};
use fs_conduction::fixtures::{box_grid, on_box_face};
use fs_conduction::transient::{VolumetricHeatCapacity,
    backward_euler::{BackwardEuler, NonlinearStepConfig, StepConfig, StepLinearization}};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use crate::conduction_assimilation::ConductionWindowConfig;
use std::cell::{Cell, RefCell};

fn with_cx(f: impl FnOnce(&Cx<'_>)) {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| f(&Cx::new(&gate, arena,
        StreamKey { seed: 61, kernel_id: 820, tile: 0, iteration: 0 }, Budget::INFINITE, ExecMode::Deterministic)));
}
struct Domain { mesh: ConductionMesh, boundary: ThermalBoundary, material: ConductivityModel }
impl Domain {
    fn new(nonlinear: bool) -> Self {
        let (complex, positions) = box_grid([1,1,1], [1.0,1.0,1.0]);
        let mesh = ConductionMesh::new(complex, positions).unwrap();
        let builder = ThermalBoundaryBuilder::new(&mesh);
        let builder = if nonlinear {
            builder.region("fixed", |f| on_box_face(f.centroid[0], 0.0), ThermalBc::dirichlet(300.0).unwrap()).unwrap()
        } else { builder };
        let boundary = builder.adiabatic_remainder().finish().unwrap();
        let material = if nonlinear {
            ConductivityModel::isotropic(ConductivityTable::declared_curve(vec![(250.0,1.0),(450.0,21.0)]).unwrap())
        } else { ConductivityModel::isotropic_declared(2.0).unwrap() };
        Self { mesh, boundary, material }
    }
}
struct Model<'a> {
    engine: &'a BackwardEuler<'a>, domain: &'a Domain, source: ScalarField, double: ScalarField,
    scheduled: bool, changed: Cell<bool>, broken: Cell<bool>, calls: Cell<usize>,
    forward_times: RefCell<Vec<ConductionSubstep>>, reverse_times: RefCell<Vec<ConductionSubstep>>,
}
impl<'a> Model<'a> {
    fn new(engine: &'a BackwardEuler<'a>, domain: &'a Domain, amplitude: f64, scheduled: bool) -> Self {
        Self { engine, domain, source: ScalarField::Uniform(amplitude), double: ScalarField::Uniform(2.0*amplitude),
            scheduled, changed: Cell::new(false), broken: Cell::new(false), calls: Cell::new(0),
            forward_times: RefCell::new(Vec::new()), reverse_times: RefCell::new(Vec::new()) }
    }
    fn weight(&self, time: ConductionSubstep) -> f64 {
        if self.scheduled && time.index % 2 == 1 { 2.0 } else { 1.0 }
    }
}
impl ConductionWindowModel for Model<'_> {
    fn engine(&self) -> &BackwardEuler<'_> { self.engine }
    fn problem(&self, _: usize) -> Result<ConductionProblem<'_>, ConductionError> {
        self.calls.set(self.calls.get()+1);
        Ok(ConductionProblem { mesh: &self.domain.mesh, boundary: &self.domain.boundary,
            material: &self.domain.material, element_materials: None, source: &self.source })
    }
    fn problem_at(&self, time: ConductionSubstep) -> Result<ConductionProblem<'_>, ConductionError> {
        self.forward_times.borrow_mut().push(time);
        let mut problem = self.problem(time.interval)?;
        if self.weight(time) == 2.0 || self.changed.get() { problem.source = &self.double; }
        Ok(problem)
    }
    fn parameter_count(&self) -> usize { 1 }
    fn parameter_pullback(&self, _: usize, cx: &Cx<'_>, step: &StepLinearization<'_>, seed: &[f64], out: &mut [f64])
        -> Result<(), ConductionError> {
        out[0] = step.source_density_pullback(cx, seed)?.iter().sum(); Ok(())
    }
    fn parameter_pullback_at(&self, time: ConductionSubstep, cx: &Cx<'_>, step: &StepLinearization<'_>, seed: &[f64], out: &mut [f64])
        -> Result<(), ConductionError> {
        self.reverse_times.borrow_mut().push(time);
        self.parameter_pullback(time.interval, cx, step, seed, out)?;
        out[0] *= self.weight(time);
        if self.broken.get() { out[0] = f64::NAN; }
        Ok(())
    }
}
fn config(nonlinear: bool) -> ConductionWindowConfig {
    ConductionWindowConfig { step: StepConfig { linear: LinearConfig {
        tolerance: 1e-12, max_iterations: 2000, restart: 16 }, energy_tolerance_j: 1e-8 },
        nonlinear: nonlinear.then_some(NonlinearStepConfig::default()), max_vertices: 32,
        max_elements: 128, max_intervals: 8, max_parameters: 4 }
}
fn limits() -> SubstepLimits {
    SubstepLimits { max_steps: 64, max_record_components: 1024, checkpoints: 8, replayed_steps: 1024 }
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a-b).abs() <= tolerance*b.abs().max(1.0), "{a:e} != {b:e}");
}
fn dot(a: &[f64], b: &[f64]) -> f64 { a.iter().zip(b).map(|(a,b)| a*b).sum() }

#[test]
fn one_substep_retains_native_endpoint_pullback_and_zero_replay_behavior() {
    with_cx(|cx| {
        let d = Domain::new(false);
        let engine = BackwardEuler::uniform(cx, &d.mesh, VolumetricHeatCapacity::declared(5.0).unwrap()).unwrap();
        let m = Model::new(&engine, &d, 10.0, false);
        let base = ConductionWindowPolicy::new(cx, &d.mesh, &d.boundary, &[0.0,2.0], config(false), &mut || false).unwrap();
        let policy = ConductionSubsteps::new(base.clone(), &[1], limits(), &mut || false).unwrap();
        let initial = vec![300.0; base.dimension()]; let seed = vec![1.0/base.dimension() as f64; base.dimension()];
        let ordinary = base.record(&m, 0, 0.0, 2.0, &initial, &mut || false).unwrap();
        let refined = policy.record(&m, 0, 0.0, 2.0, &initial, &mut || false).unwrap();
        assert_eq!(ordinary.endpoint(), refined.endpoint());
        assert_eq!(ordinary.pullback(&seed, &[0.3], &mut || false).unwrap(), refined.pullback(&seed, &[0.3], &mut || false).unwrap());
        assert_eq!(refined.accepted_steps(), 1);
    });
}

#[test]
fn scheduled_uniform_heat_has_exact_chain_rule_and_identical_replay_contexts() {
    with_cx(|cx| {
        let d = Domain::new(false);
        let engine = BackwardEuler::uniform(cx, &d.mesh, VolumetricHeatCapacity::declared(5.0).unwrap()).unwrap();
        let m = Model::new(&engine, &d, 10.0, true);
        let base = ConductionWindowPolicy::new(cx, &d.mesh, &d.boundary, &[0.0,2.0], config(false), &mut || false).unwrap();
        let policy = ConductionSubsteps::new(base, &[4], limits(), &mut || false).unwrap();
        let n = policy.base().dimension(); let initial = vec![300.0; n]; let seed = vec![1.0/n as f64; n];
        let tape = policy.record(&m, 0, 0.0, 2.0, &initial, &mut || false).unwrap();
        for &x in tape.endpoint() { close(x, 306.0, 1e-10); }
        let recorded = m.forward_times.borrow().clone();
        let gradient = tape.pullback(&seed, &[0.3], &mut || false).unwrap();
        close(gradient.initial.iter().sum(), 1.0, 1e-10);
        close(gradient.parameters[0], 0.9, 1e-10);
        assert_eq!(gradient.replayed_steps, 8); assert_eq!(gradient.peak_checkpoints, 3);
        assert_eq!(m.reverse_times.borrow().as_slice(), recorded.iter().rev().copied().collect::<Vec<_>>());
        for time in m.forward_times.borrow().iter() { assert_eq!(*time, recorded[time.index]); }
        assert_eq!(m.calls.get(), 4+gradient.replayed_steps);
    });
}

#[test]
fn checkpointed_nonlinear_constrained_gradients_match_perturbed_production_solutions() {
    with_cx(|cx| {
        let d = Domain::new(true);
        let engine = BackwardEuler::uniform(cx, &d.mesh, VolumetricHeatCapacity::declared(5.0).unwrap()).unwrap();
        let base = ConductionWindowPolicy::new(cx, &d.mesh, &d.boundary, &[0.0,0.6], config(true), &mut || false).unwrap();
        let policy = ConductionSubsteps::new(base, &[3], limits(), &mut || false).unwrap();
        let n = policy.base().dimension(); let x: Vec<_> = (0..n).map(|i| 301.0+0.1*i as f64).collect();
        let seed: Vec<_> = (0..n).map(|i| 0.1+0.02*i as f64).collect();
        let m = Model::new(&engine, &d, 2.0, true);
        let tape = policy.record(&m, 0, 0.0, 0.6, &x, &mut || false).unwrap();
        let g = tape.pullback(&seed, &[0.0], &mut || false).unwrap();
        let evaluate = |point: &[f64], p: f64| {
            let model = Model::new(&engine, &d, p, true);
            let t = policy.record(&model, 0, 0.0, 0.6, point, &mut || false).unwrap();
            dot(&seed, t.endpoint())
        };
        let eps = 1e-4;
        for i in 0..n {
            let (mut a, mut b) = (x.clone(), x.clone()); a[i] += eps; b[i] -= eps;
            close(g.initial[i], (evaluate(&a,2.0)-evaluate(&b,2.0))/(2.0*eps), 3e-5);
        }
        close(g.parameters[0], (evaluate(&x,2.0+eps)-evaluate(&x,2.0-eps))/(2.0*eps), 3e-5);
        let full = policy.base().expand_field(tape.endpoint()).unwrap();
        for &(node, value) in d.boundary.dirichlet() { assert_eq!(full[node], value); }
    });
}

#[test]
fn replay_detects_changed_forward_and_rejects_missing_parameter_derivatives() {
    with_cx(|cx| {
        let d = Domain::new(false);
        let engine = BackwardEuler::uniform(cx, &d.mesh, VolumetricHeatCapacity::declared(5.0).unwrap()).unwrap();
        let m = Model::new(&engine, &d, 1.0, false);
        let base = ConductionWindowPolicy::new(cx, &d.mesh, &d.boundary, &[0.0,1.0], config(false), &mut || false).unwrap();
        let policy = ConductionSubsteps::new(base, &[4], limits(), &mut || false).unwrap();
        let n = policy.base().dimension(); let initial = vec![300.0;n]; let seed = vec![1.0/n as f64;n];
        let tape = policy.record(&m,0,0.0,1.0,&initial,&mut || false).unwrap(); let before = tape.endpoint().to_vec();
        let expected = tape.pullback(&seed,&[0.0],&mut || false).unwrap();
        m.changed.set(true);
        assert!(matches!(tape.pullback(&seed,&[0.0],&mut || false), Err(WindowError::Integrator { phase: "substep replay", .. })));
        m.changed.set(false); m.broken.set(true);
        assert!(matches!(tape.pullback(&seed,&[0.0],&mut || false), Err(WindowError::NonFinite(_))));
        m.broken.set(false);
        assert_eq!(tape.endpoint(), before); assert_eq!(tape.pullback(&seed,&[0.0],&mut || false).unwrap(), expected);
    });
}

#[test]
fn cancelled_reverse_is_retryable_and_replay_limit_is_hard() {
    with_cx(|cx| {
        let d = Domain::new(false);
        let engine = BackwardEuler::uniform(cx, &d.mesh, VolumetricHeatCapacity::declared(5.0).unwrap()).unwrap();
        let m = Model::new(&engine, &d, 1.0, false);
        let base = ConductionWindowPolicy::new(cx, &d.mesh, &d.boundary, &[0.0,1.0], config(false), &mut || false).unwrap();
        let policy = ConductionSubsteps::new(base.clone(), &[4], limits(), &mut || false).unwrap();
        let n = policy.base().dimension(); let initial = vec![300.0;n]; let seed = vec![1.0/n as f64;n];
        let tape = policy.record(&m,0,0.0,1.0,&initial,&mut || false).unwrap();
        let before = tape.endpoint().to_vec(); let calls = m.calls.get();
        assert!(matches!(tape.pullback(&seed,&[0.0],&mut || m.calls.get()>calls), Err(WindowError::Cancelled)));
        assert_eq!(tape.endpoint(), before);
        let expected = tape.pullback(&seed,&[0.0],&mut || false).unwrap();
        let short = ConductionSubsteps::new(base, &[4], SubstepLimits {
            replayed_steps: expected.replayed_steps-1, ..limits()
        }, &mut || false).unwrap();
        let t = short.record(&m,0,0.0,1.0,&initial,&mut || false).unwrap(); let calls = m.calls.get();
        assert!(matches!(t.pullback(&seed,&[0.0],&mut || false), Err(WindowError::Integrator { phase: "substep replay", .. })));
        assert_eq!(m.calls.get()-calls, expected.replayed_steps-1);
        assert_eq!(t.endpoint(), before);
        assert_eq!(tape.pullback(&seed,&[0.0],&mut || false).unwrap(), expected);
    });
}

#[test]
fn grid_and_storage_admission_happens_before_physical_forecasting() {
    with_cx(|cx| {
        let d = Domain::new(false);
        let base = ConductionWindowPolicy::new(cx,&d.mesh,&d.boundary,&[0.0,1.0,1.3],config(false),&mut || false).unwrap();
        for counts in [vec![0,1],vec![1],vec![usize::MAX,1]] {
            assert!(ConductionSubsteps::new(base.clone(),&counts,limits(),&mut || false).is_err());
        }
        for cap in [SubstepLimits { max_steps: 3, ..limits() }, SubstepLimits { checkpoints: 1, ..limits() },
            SubstepLimits { max_record_components: 1, ..limits() }] {
            assert!(ConductionSubsteps::new(base.clone(),&[4,4],cap,&mut || false).is_err());
        }
        let policy = ConductionSubsteps::new(base,&[3,7],limits(),&mut || false).unwrap();
        assert_eq!(policy.times(), &[0.0,1.0,1.3]);
        for k in 0..2 {
            let times = policy.substep_times(k).unwrap();
            assert_eq!(times[0].to_bits(), policy.times()[k].to_bits());
            assert_eq!(times[times.len()-1].to_bits(), policy.times()[k+1].to_bits());
            assert!(times.windows(2).all(|t| t[0]<t[1]));
        }
        assert!(policy.substep_times(2).is_none());
        let large = 9_007_199_254_740_992.0;
        let base = ConductionWindowPolicy::new(cx,&d.mesh,&d.boundary,&[large,large+2.0],config(false),&mut || false).unwrap();
        assert!(ConductionSubsteps::new(base,&[4],limits(),&mut || false).is_err());
    });
}

#[test]
fn refinement_adds_no_model_error_controls_and_study_continuation_is_exact() {
    use crate::transient::variational::{WeakConstraintWindow, WindowControl, WindowObjective};
    use crate::transient::variational::study::{WeakConstraintStudy, StudySettings};
    struct Loss;
    impl WindowObjective for Loss {
        fn evaluate(&self, _: &[f64], n: usize, states: &[f64], bar: &mut [f64], _: &mut dyn FnMut()->bool) -> Result<f64,String> {
            let mut cost = 0.0;
            for (i,(x,b)) in states.iter().zip(bar).enumerate() {
                let target = 300.2+0.3*(i/n) as f64; *b = (x-target)/0.04;
                cost += 0.5*(x-target)*(*b);
            }
            Ok(cost)
        }
    }
    with_cx(|cx| {
        let d = Domain::new(false);
        let engine = BackwardEuler::uniform(cx,&d.mesh,VolumetricHeatCapacity::declared(5.0).unwrap()).unwrap();
        let m = Model::new(&engine,&d,1.0,false);
        let base = ConductionWindowPolicy::new(cx,&d.mesh,&d.boundary,&[0.0,1.0,2.0],config(false),&mut || false).unwrap();
        let n = base.dimension();
        let refined = ConductionSubsteps::new(base,&[4,8],limits(),&mut || false).unwrap();
        let w = WeakConstraintWindow::new(refined.times(),&vec![300.0;3*n],&vec![1.0;n],&vec![0.5;n],&vec![0.2;2*n],1000).unwrap();
        assert_eq!(w.control_dimension(),3*n);
        let settings = StudySettings { memory: 8, gradient_tolerance: 1e-6, max_evaluations: 400, max_optimizer_components: 10000 };
        let (mut a,mut b) = (WindowControl::new(400,800,10000),WindowControl::new(400,800,10000));
        let mut one = WeakConstraintStudy::new(&w,&m,&Loss,&vec![0.0;3*n],refined.clone(),settings,&mut a,&mut || false).unwrap();
        let mut split = WeakConstraintStudy::new(&w,&m,&Loss,&vec![0.0;3*n],refined,settings,&mut b,&mut || false).unwrap();
        one.run(100,&mut a,&mut || false).unwrap();
        split.run(2,&mut b,&mut || false).unwrap(); let mut split = split.clone();
        split.run(98,&mut b,&mut || false).unwrap();
        assert_eq!(one.accepted(),split.accepted());
        assert_eq!(one.optimizer().x,split.optimizer().x);
        assert_eq!(one.optimizer().history,split.optimizer().history);
        assert_eq!(one.accepted().states.len(),3*n);
        assert_eq!(one.accepted().defects.len(),2*n);
        assert_eq!(one.accepted().accepted_steps,12);
        assert!(one.accepted().observation_value<1.0);
    });
}
