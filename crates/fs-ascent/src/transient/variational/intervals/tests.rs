use super::*;
use super::super::{WeakConstraintWindow, WindowControl, WindowObjective};
use super::super::study::{StudySettings, WeakConstraintStudy};
use std::cell::Cell;

// Independent closed affine-map fixture, deliberately NOT an OdeVjp. It checks
// that a discrete solver's state/parameter derivatives enter the same objective.
struct Affine { parameter: f64, calls: Cell<usize> }
#[derive(Clone, Copy)]
struct Scheme { malformed: u8 }
struct Tape { value: Vec<f64>, time: f64, a: f64, b: f64, malformed: u8 }
impl IntervalScheme<Affine> for Scheme {
    type Tape<'a> = Tape;
    fn dimension(&self, _: &Affine) -> usize { 1 }
    fn parameter_count(&self, _: &Affine) -> usize { 1 }
    fn validate(&self, _: &[f64]) -> Result<(), WindowError> { Ok(()) }
    fn record<'a>(&'a self, model: &'a Affine, _: usize, start: f64, end: f64,
        initial: &[f64], _: &mut dyn FnMut() -> bool) -> Result<Tape, WindowError>
    {
        model.calls.set(model.calls.get()+1);
        let dt = end-start; let a = 1.0/(1.0+2.0*dt); let b = dt*a;
        let value = if self.malformed == 1 { vec![] } else { vec![a*initial[0]+b*model.parameter] };
        Ok(Tape { value, time: if self.malformed == 2 { start } else { end }, a, b, malformed: self.malformed })
    }
}
impl IntervalTape for Tape {
    fn endpoint(&self) -> &[f64] { &self.value }
    fn end_time(&self) -> f64 { self.time }
    fn accepted_steps(&self) -> usize { 1 }
    fn pullback(&self, seed: &[f64], direct: &[f64], _: &mut dyn FnMut() -> bool)
        -> Result<TrajectoryGradient, WindowError>
    {
        Ok(TrajectoryGradient {
            initial: if self.malformed == 3 { vec![] } else { vec![self.a*seed[0]] },
            parameters: vec![if self.malformed == 4 { f64::NAN } else { self.b*seed[0]+direct[0] }],
            replayed_steps: 1, peak_checkpoints: 1,
        })
    }
}
struct Loss;
impl WindowObjective for Loss {
    fn evaluate(&self, _: &[f64], _: usize, x: &[f64], bar: &mut [f64], _: &mut dyn FnMut() -> bool)
        -> Result<f64, String>
    {
        for (g,x) in bar.iter_mut().zip(x) { *g = x-0.7; }
        Ok(bar.iter().map(|g| 0.5*g*g).sum())
    }
}
fn window() -> WeakConstraintWindow {
    WeakConstraintWindow::new(&[0.0,0.25,0.75], &[1.2,1.0,0.9], &[2.0], &[0.5], &[0.3,0.4], 32).unwrap()
}
fn control() -> WindowControl { WindowControl::new(2000,4000,1000) }
fn settings() -> StudySettings { StudySettings {
    memory: 5, gradient_tolerance: 1e-8, max_evaluations: 1000, max_optimizer_components: 10000,
} }

#[test]
fn discrete_only_model_includes_both_endpoint_terms_and_parameter_gradient() {
    let w=window(); let m=Affine {parameter:0.2,calls:Cell::new(0)};
    let z=[0.1,-0.05,0.12]; let mut c=control();
    let result=w.evaluate_using(&m,&Loss,&z,&Scheme{malformed:0},&mut c,&mut||false).unwrap();
    let x=&result.states; let mut gradient=x.iter().map(|x|x-0.7).collect::<Vec<_>>();
    gradient[0]+=(x[0]-1.2)/0.25;
    let mut gp=0.0;
    for (k,(dt,sigma)) in [(0.25,0.3),(0.5,0.4)].into_iter().enumerate() {
        let a=1.0/(1.0+2.0*dt);let b=dt*a;let defect=x[k+1]-a*x[k]-b*m.parameter;
        let r=defect/(sigma*sigma);gradient[k]-=a*r;gradient[k+1]+=r;gp-=b*r;
    }
    for (a,b) in result.gradient.iter().zip(gradient) {assert!((a-2.0*b).abs()<1e-12);}
    assert!((result.parameter_gradient[0]-gp).abs()<1e-12);
    assert_eq!(result.accepted_steps,2);assert_eq!(m.calls.get(),2);assert_eq!(c.interval_attempts(),2);
}

#[test]
fn existing_optimizer_can_solve_and_resume_a_discrete_only_window() {
    let w=window();let m=Affine{parameter:0.2,calls:Cell::new(0)};
    let mut c=control();let mut whole=WeakConstraintStudy::new(&w,&m,&Loss,&[0.0;3],Scheme{malformed:0},settings(),&mut c,&mut||false).unwrap();
    let before=whole.accepted().value;let mut split=whole.clone();let mut split_c=control();
    let end=whole.run(100,&mut c,&mut||false).unwrap();
    for _ in 0..100 {if split.run(1,&mut split_c,&mut||false).unwrap().reason!=crate::StopReason::IterationCap {break;}}
    assert_eq!(end.reason,crate::StopReason::GradNorm);
    assert!(whole.accepted().value<before);assert_eq!(whole.accepted(),split.accepted());
    assert_eq!(whole.optimizer().history,split.optimizer().history);
}

#[test]
fn malformed_backend_outputs_never_escape_as_partial_gradients() {
    let w=window();let m=Affine{parameter:0.2,calls:Cell::new(0)};
    for malformed in 1..=4 {
        let mut c=control();let result=w.evaluate_using(&m,&Loss,&[0.0;3],&Scheme{malformed},&mut c,&mut||false);
        assert!(matches!(result,Err(WindowError::IntervalOutput{..})|Err(WindowError::NonFinite(_))));
        assert_eq!(c.evaluations(),1);assert_eq!(c.interval_attempts(),1);
    }
}

#[test]
fn cancellation_after_forecast_is_charged_and_retryable() {
    let w=window();let m=Affine{parameter:0.2,calls:Cell::new(0)};let mut c=control();
    assert_eq!(w.evaluate_using(&m,&Loss,&[0.0;3],&Scheme{malformed:0},&mut c,&mut||m.calls.get()>0),Err(WindowError::Cancelled));
    assert_eq!(c.interval_attempts(),1);
    assert!(w.evaluate_using(&m,&Loss,&[0.0;3],&Scheme{malformed:0},&mut c,&mut||false).is_ok());
    assert_eq!(c.evaluations(),2);assert_eq!(c.interval_attempts(),3);
}

#[test]
fn rk45_adapter_reuses_the_exact_recorded_endpoint_and_pullback() {
    struct Decay;
    impl OdeVjp for Decay {
        fn dimension(&self)->usize {1} fn parameter_count(&self)->usize {0}
        fn rhs(&self,_:f64,x:&[f64],out:&mut[f64]) {out[0]=-0.7*x[0];}
        fn rhs_vjp(&self,_:f64,_:&[f64],b:&[f64],x:&mut[f64],_:&mut[f64])->Result<(),String> {x[0]=-0.7*b[0];Ok(())}
    }
    let p=IntervalPolicy { recording:fs_time::adaptive::adjoint::trajectory::RecordingConfig {
        end:0.75,rtol:1e-10,atol:1e-12,controller:fs_time::PiController::default(),max_workspace_components:256,
    },initial_step:0.1,max_attempts:10000,max_records:10000,
        replay:fs_time::adaptive::adjoint::trajectory::ReplayBudget{checkpoints:32,replayed_steps:10000} };
    let tape=p.record(&Decay,0,0.0,0.75,&[1.2],&mut||false).unwrap();
    let mut direct=RecordedRk45::new(&Decay,AdaptiveState::new(0.0,&[1.2],0.1),p.recording.clone()).unwrap();
    direct.advance(p.max_attempts,p.max_records,&mut||false).unwrap();
    assert_eq!(tape.endpoint(),direct.state().u);
    assert_eq!(tape.pullback(&[0.8],&[],&mut||false).unwrap(),direct.pullback(&[0.8],&[],p.replay,&mut||false).unwrap());
}
