use super::*;
use std::cell::Cell;

struct Model {
    p: Cell<f64>, calls: Cell<usize>, live: Cell<usize>, peak: Cell<usize>, broken: Cell<bool>,
}
impl Model {
    fn new(p: f64) -> Self {
        Self { p: Cell::new(p), calls: Cell::new(0), live: Cell::new(0), peak: Cell::new(0), broken: Cell::new(false) }
    }
}
#[derive(Clone)]
struct Scheme;
// Nonlinear, nonautonomous map deliberately NOT implementing OdeVjp. The global
// substep index changes the forcing and must survive coarse-window partitioning.
fn map(x: &[f64], p: f64, i: usize, h: f64) -> Vec<f64> {
    vec![x[0]+h*(p*x[1]+0.2*x[0]*x[0]), x[1]+h*(-0.3*x[0]+p*(i+1) as f64)]
}
struct MapTape<'a> { m: &'a Model, old: Vec<f64>, value: Vec<f64>, end: f64, h: f64, i: usize }
impl Drop for MapTape<'_> { fn drop(&mut self) { self.m.live.set(self.m.live.get()-1); } }
impl IntervalScheme<Model> for Scheme {
    type Tape<'a> = MapTape<'a> where Self: 'a, Model: 'a;
    fn dimension(&self, _: &Model) -> usize { 2 }
    fn parameter_count(&self, _: &Model) -> usize { 1 }
    fn validate(&self, _: &[f64]) -> Result<(), WindowError> { Ok(()) }
    fn record<'a>(&'a self, model: &'a Model, i: usize, start: f64, end: f64, x: &[f64],
        check: &mut dyn FnMut() -> bool) -> Result<Self::Tape<'a>, WindowError>
    {
        poll(check)?; model.calls.set(model.calls.get()+1);
        model.live.set(model.live.get()+1); model.peak.set(model.peak.get().max(model.live.get()));
        Ok(MapTape { m: model, old: x.to_vec(), value: map(x, model.p.get(), i, end-start), end, h: end-start, i })
    }
}
impl IntervalTape for MapTape<'_> {
    fn endpoint(&self) -> &[f64] { &self.value }
    fn end_time(&self) -> f64 { self.end }
    fn accepted_steps(&self) -> usize { 1 }
    fn pullback(&self, b: &[f64], direct: &[f64], check: &mut dyn FnMut() -> bool) -> Result<TrajectoryGradient, WindowError> {
        poll(check)?;
        Ok(TrajectoryGradient {
            initial: if self.m.broken.get() { vec![f64::NAN] } else {
                vec![(1.0+0.4*self.h*self.old[0])*b[0]-0.3*self.h*b[1], self.h*self.m.p.get()*b[0]+b[1]]
            },
            parameters: vec![direct[0]+self.h*(self.old[1]*b[0]+(self.i+1) as f64*b[1])],
            replayed_steps: 0, peak_checkpoints: 0,
        })
    }
}
fn budget() -> SubstepBudget { SubstepBudget { max_state_components: 2, max_parameters: 1, checkpoints: 16, replayed_substeps: 1000 } }
fn policy(count: usize) -> CheckpointedIntervals<Scheme> {
    CheckpointedIntervals::new(Scheme, SubstepGrid::uniform(&[0.0, 0.7], count, 128).unwrap(), budget())
}
fn close(a: f64, b: f64, tol: f64) { assert!((a-b).abs() < tol*b.abs().max(1.0), "{a} != {b}"); }
fn replay_count(n: usize) -> usize {
    if n == 1 { 1 } else { let m = n/2+n%2; m+replay_count(m)+replay_count(n-m) }
}

#[test]
fn nonlinear_chain_matches_full_storage_and_independent_differences() {
    for count in [1, 2, 3, 8, 17] {
        let model = Model::new(0.2); let policy = policy(count); let x = [0.5,-0.2]; let seed = [0.7,-0.3];
        let tape = policy.record(&model, 0, 0.0, 0.7, &x, &mut || false).unwrap();
        let got = tape.pullback(&seed, &[0.13], &mut || false).unwrap();
        let mut states = vec![x.to_vec()];
        for (i, t) in policy.grid.fine.windows(2).enumerate() { states.push(map(states.last().unwrap(), model.p.get(), i, t[1]-t[0])); }
        assert_eq!(tape.endpoint(), states.last().unwrap());
        let mut expected = seed.to_vec(); let mut pbar = vec![0.13];
        for i in (0..count).rev() {
            let t = &policy.grid.fine;
            let step = Scheme.record(&model, i, t[i], t[i+1], &states[i], &mut || false).unwrap();
            let g = step.pullback(&expected, &pbar, &mut || false).unwrap(); expected = g.initial; pbar = g.parameters;
        }
        assert_eq!(got.initial, expected); assert_eq!(got.parameters, pbar);
        assert_eq!(got.replayed_steps, replay_count(count));
        assert!(got.peak_checkpoints <= budget().checkpoints);
        assert_eq!(model.peak.get(), 1, "never retain more than one inner tape");
        let value = |x: &[f64], p: f64| {
            let mut y = x.to_vec();
            for (i,t) in policy.grid.fine.windows(2).enumerate() { y = map(&y,p,i,t[1]-t[0]); }
            y.iter().zip(seed).map(|(a,b)| a*b).sum::<f64>()+0.13*p
        };
        let eps=1e-6;
        for i in 0..2 { let (mut plus,mut minus)=(x,x);plus[i]+=eps;minus[i]-=eps;
            close(got.initial[i],(value(&plus,0.2)-value(&minus,0.2))/(2.0*eps),2e-8); }
        close(got.parameters[0],(value(&x,0.2+eps)-value(&x,0.2-eps))/(2.0*eps),2e-8);
    }
}

#[test]
fn partitions_keep_global_forcing_addresses_and_add_direct_partials_once() {
    let grid=SubstepGrid::new(&[0.0,0.1,0.3,0.6,0.7],&[0,2,4],8).unwrap();
    let p=CheckpointedIntervals::new(Scheme,grid,budget()); let m=Model::new(0.3); let x=[0.4,-0.1];
    let first=p.record(&m,0,0.0,0.3,&x,&mut||false).unwrap();
    let second=p.record(&m,1,0.3,0.7,first.endpoint(),&mut||false).unwrap();
    let b=second.pullback(&[1.0,-0.2],&[7.0],&mut||false).unwrap();
    let got=first.pullback(&b.initial,&b.parameters,&mut||false).unwrap();
    let all=CheckpointedIntervals::new(Scheme,SubstepGrid::new(p.grid.fine_times(),&[0,4],8).unwrap(),budget());
    let t=all.record(&m,0,0.0,0.7,&x,&mut||false).unwrap();
    assert_eq!(t.endpoint(),second.endpoint());
    let expected=t.pullback(&[1.0,-0.2],&[7.0],&mut||false).unwrap();
    assert_eq!(got.initial,expected.initial);assert_eq!(got.parameters,expected.parameters);
}

#[test]
fn validates_clocks_and_limits_before_model_work() {
    for (fine,knots) in [(vec![0.0,0.1,0.1],vec![0,2]),(vec![0.0,f64::NAN,1.0],vec![0,2]),
        (vec![0.0,0.5,1.0],vec![1,2]),(vec![0.0,0.5,1.0],vec![0,3]),(vec![0.0,0.5,1.0],vec![0,0,2])] {
        assert!(SubstepGrid::new(&fine,&knots,8).is_err());
    }
    assert!(SubstepGrid::uniform(&[0.0,1.0],0,16).is_err());
    assert!(SubstepGrid::uniform(&[0.0,1.0],usize::MAX,16).is_err());
    assert!(SubstepGrid::uniform(&[1.0,1.0+f64::EPSILON],2,16).is_err());
    assert!(SubstepGrid::uniform(&[0.0,1.0,2.0],9,16).is_err());
    let model=Model::new(0.2);let mut p=policy(8);p.budget.checkpoints=3;
    assert!(p.record(&model,0,0.0,0.7,&[0.5,-0.2],&mut||false).is_err());
    assert_eq!(model.calls.get(),0);
    p.budget=budget();p.budget.max_state_components=1;
    assert!(p.record(&model,0,0.0,0.7,&[0.5,-0.2],&mut||false).is_err());
    assert_eq!(model.calls.get(),0);
}

#[test]
fn replay_budget_is_exact_and_failed_sweeps_leave_recording_retryable() {
    let model=Model::new(0.2);let mut p=policy(9);p.budget.replayed_substeps=replay_count(9)-1;
    let t=p.record(&model,0,0.0,0.7,&[0.5,-0.2],&mut||false).unwrap();let endpoint=t.endpoint().to_vec();
    assert!(t.pullback(&[1.0,0.0],&[0.0],&mut||false).is_err());assert_eq!(endpoint,t.endpoint());
    let mut p=policy(9);p.budget.replayed_substeps=replay_count(9);
    let t=p.record(&model,0,0.0,0.7,&[0.5,-0.2],&mut||false).unwrap();
    assert_eq!(t.pullback(&[1.0,0.0],&[0.0],&mut||false).unwrap().replayed_steps,replay_count(9));
}

#[test]
fn changed_replay_and_bad_pullbacks_are_refused_without_partial_results() {
    let model=Model::new(0.2);let p=policy(5);
    let t=p.record(&model,0,0.0,0.7,&[0.5,-0.2],&mut||false).unwrap();
    let expected=t.pullback(&[1.0,0.0],&[0.0],&mut||false).unwrap();
    model.p.set(0.3);
    assert!(matches!(t.pullback(&[1.0,0.0],&[0.0],&mut||false),Err(WindowError::Integrator{phase:"substep replay",..})));
    model.p.set(0.2);model.broken.set(true);
    assert!(matches!(t.pullback(&[1.0,0.0],&[0.0],&mut||false),Err(WindowError::IntervalOutput{..})));
    model.broken.set(false);
    assert_eq!(expected,t.pullback(&[1.0,0.0],&[0.0],&mut||false).unwrap());
}

#[test]
fn cancellation_during_forward_and_reverse_retains_no_inner_tapes() {
    let model=Model::new(0.2);let p=policy(7);let x=[0.5,-0.2];
    assert!(matches!(p.record(&model,0,0.0,0.7,&x,&mut||model.calls.get()>=3),Err(WindowError::Cancelled)));
    assert_eq!(model.live.get(),0);
    let t=p.record(&model,0,0.0,0.7,&x,&mut||false).unwrap();let before=t.endpoint().to_vec();let calls=model.calls.get();
    assert!(matches!(t.pullback(&[1.0,0.0],&[0.0],&mut||model.calls.get()>calls+2),Err(WindowError::Cancelled)));
    assert_eq!(model.live.get(),0);assert_eq!(t.endpoint(),before);
    assert!(t.pullback(&[1.0,0.0],&[0.0],&mut||false).is_ok());
}

#[test]
fn refined_solver_does_not_add_window_controls_or_model_error_terms() {
    use crate::transient::variational::{WeakConstraintWindow,WindowControl,WindowObjective};
    struct Empty;
    impl WindowObjective for Empty {
        fn evaluate(&self,_:&[f64],_:usize,_:&[f64],bar:&mut[f64],_:&mut dyn FnMut()->bool)->Result<f64,String> {bar.fill(0.0);Ok(0.0)}
    }
    let p=policy(8);let m=Model::new(0.2);let x=[0.5,-0.2];
    let t=p.record(&m,0,0.0,0.7,&x,&mut||false).unwrap();
    let reference=[x[0],x[1],t.endpoint()[0]+0.1,t.endpoint()[1]-0.2];
    let w=WeakConstraintWindow::new(&[0.0,0.7],&reference,&[1.0,1.0],&[1.0,1.0],&[0.5,0.5],32).unwrap();
    let result=w.evaluate_using(&m,&Empty,&[0.0;4],&p,&mut WindowControl::new(1,1,128),&mut||false).unwrap();
    assert_eq!(w.control_dimension(),4);assert_eq!(result.defects.len(),2);assert_eq!(result.accepted_steps,8);
    close(result.model_value,0.1,1e-13);
    let expected=t.pullback(&[-0.4,0.8],&[0.0],&mut||false).unwrap();
    for (a,b) in result.gradient[..2].iter().zip(&expected.initial) {close(*a,*b,1e-12);}
    close(result.parameter_gradient[0],expected.parameters[0],1e-12);
}

#[test]
fn production_rk45_substeps_match_retained_native_tapes_bitwise() {
    use crate::transient::variational::IntervalPolicy;
    use fs_time::{AdaptiveState, PiController};
    use fs_time::adaptive::adjoint::OdeVjp;
    use fs_time::adaptive::adjoint::trajectory::{RecordedRk45, RecordingConfig, ReplayBudget};
    struct Decay;
    impl OdeVjp for Decay {
        fn dimension(&self)->usize {1} fn parameter_count(&self)->usize {1}
        fn rhs(&self,t:f64,x:&[f64],out:&mut[f64]) {out[0]=-0.7*x[0]+0.03*t;}
        fn rhs_vjp(&self,_:f64,x:&[f64],seed:&[f64],xb:&mut[f64],pb:&mut[f64])->Result<(),String> {
            xb[0]=-0.7*seed[0];pb[0]=-x[0]*seed[0];Ok(())
        }
    }
    let base=IntervalPolicy {recording:RecordingConfig {end:1.0,rtol:1e-10,atol:1e-12,
        controller:PiController::default(),max_workspace_components:4096},initial_step:0.11,
        max_attempts:1000,max_records:1000,replay:ReplayBudget {checkpoints:32,replayed_steps:10000}};
    let grid=SubstepGrid::uniform(&[0.0,0.7],5,16).unwrap();
    let p=CheckpointedIntervals::new(base.clone(),grid.clone(),budget());let model=Decay;
    let t=p.record(&model,0,0.0,0.7,&[1.2],&mut||false).unwrap();
    let mut states=vec![1.2];let mut tapes=Vec::new();
    for times in grid.fine_times().windows(2) {
        let mut config=base.recording.clone();config.end=times[1];
        let mut native=RecordedRk45::new(&model,AdaptiveState::new(times[0],&states,base.initial_step),config).unwrap();
        native.advance(1000,1000,&mut||false).unwrap();states=native.state().u.clone();tapes.push(native);
    }
    assert_eq!(t.endpoint(),states);
    let got=t.pullback(&[0.4],&[0.19],&mut||false).unwrap();
    let mut bar=vec![0.4];let mut pb=vec![0.19];
    for native in tapes.iter().rev() {
        let b=native.pullback(&bar,&pb,base.replay,&mut||false).unwrap();bar=b.initial;pb=b.parameters;
    }
    assert_eq!(got.initial,bar);assert_eq!(got.parameters,pb);
    assert!(got.replayed_steps>replay_count(5));
}
