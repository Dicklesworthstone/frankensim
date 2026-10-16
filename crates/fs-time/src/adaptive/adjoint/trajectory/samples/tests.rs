use super::*;
use std::cell::{Cell, RefCell};

struct Decay { rate: Cell<f64>, calls: Cell<usize> }
impl Decay { fn new() -> Self { Self { rate: Cell::new(-0.7), calls: Cell::new(0) } } }
impl OdeVjp for Decay {
    fn dimension(&self) -> usize { 1 }
    fn parameter_count(&self) -> usize { 1 }
    fn rhs(&self, _: f64, u: &[f64], out: &mut [f64]) {
        self.calls.set(self.calls.get() + 1); out[0] = self.rate.get() * u[0];
    }
    fn rhs_vjp(&self, _: f64, u: &[f64], seed: &[f64], x: &mut [f64], p: &mut [f64]) -> Result<(), String> {
        x[0] = self.rate.get() * seed[0]; p[0] = u[0] * seed[0]; Ok(())
    }
}
fn config(end: f64) -> RecordingConfig { RecordingConfig {
    end, rtol: 1e-10, atol: 1e-12, controller: PiController::default(), max_workspace_components: 4096,
} }
fn budget() -> ReplayBudget { ReplayBudget { checkpoints: 64, replayed_steps: 100_000 } }
fn record<'a>(m: &'a Decay, times: &[f64]) -> RecordedRk45<'a, Decay> {
    RecordedRk45::new_sampled(m, AdaptiveState::new(0.0, &[1.2], 1.0), config(2.0), times, 100).unwrap()
}
struct Loss { seen: RefCell<Vec<usize>>, fail: Cell<bool> }
impl Loss { fn new() -> Self { Self { seen: RefCell::new(Vec::new()), fail: Cell::new(false) } } }
impl SampleObjective for Loss {
    fn evaluate(&self, i: usize, t: f64, u: &[f64], x: &mut [f64], p: &mut [f64]) -> Result<f64, String> {
        self.seen.borrow_mut().push(i);
        if self.fail.get() { return Err("sensor unavailable".into()); }
        let residual = u[0] - (0.5 + 0.1 * t);
        x[0] = residual;
        // A direct parameter term equal to 0.2*p at the fixed p=-0.7.
        p[0] = 0.2;
        Ok(0.5 * residual * residual - 0.14)
    }
}

#[test]
fn all_observations_match_exact_decay_with_initial_and_repeated_times() {
    let model = Decay::new(); let times = [0.0, 0.13, 0.13, 0.71, 1.4, 2.0];
    let mut tape = record(&model, &times);
    assert_eq!(tape.advance(10000, 10000, &mut || false).unwrap().status, RecordingStatus::ReachedEnd);
    for &time in &times[1..] { assert!(tape.step_schedule().any(|(_, end, _)| end.to_bits() == time.to_bits())); }
    let objective = Loss::new(); let got = tape.pullback_samples(&objective, budget(), &mut || false).unwrap();
    let (mut value, mut gx, mut gp) = (0.0, 0.0, 0.0);
    for &t in times.iter().rev() {
        let e = (-0.7 * t).exp(); let y = 1.2 * e; let residual = y - (0.5 + 0.1 * t);
        value += 0.5 * residual * residual - 0.14; gx += residual * e; gp += residual * t * y + 0.2;
    }
    assert!((got.value - value).abs() < 2e-9);
    assert!((got.gradient.initial[0] - gx).abs() < 2e-9);
    assert!((got.gradient.parameters[0] - gp).abs() < 2e-9);
    assert_eq!(got.observations, times.len());
    assert_eq!(*objective.seen.borrow(), (0..times.len()).rev().collect::<Vec<_>>());
}

#[test]
fn discrete_sample_gradient_matches_full_storage_and_frozen_mesh_difference() {
    let model = Decay::new(); let times = [0.0, 0.23, 0.87, 1.23, 2.0];
    let mut tape = record(&model, &times); tape.advance(10000, 10000, &mut || false).unwrap();
    let got = tape.pullback_samples(&Loss::new(), budget(), &mut || false).unwrap();
    let mut states = vec![tape.initial.clone()];
    for r in &tape.records { states.push(replay(&model, states.last().unwrap(), r.start, r.end, r.h, &mut || false).unwrap().next); }
    let mut bar = Cotangent { initial: vec![0.0], parameters: vec![0.0] };
    let mut cost = 0.0; let mut count = 0;
    for i in (0..tape.records.len()).rev() {
        let r = &tape.records[i]; let lo = times.partition_point(|t| *t < r.end); let hi = times.partition_point(|t| *t <= r.end);
        tape.observe_range(&Loss::new(), lo..hi, &states[i+1], &mut bar, &mut cost, &mut count, &mut || false).unwrap();
        let w = replay(&model, &states[i], r.start, r.end, r.h, &mut || false).unwrap();
        let v = reverse(&model, &states[i], r.start, r.end, r.h, &bar.initial, w, &mut || false).unwrap();
        bar.initial = v.initial; bar.parameters[0] += v.parameters[0];
    }
    tape.observe_range(&Loss::new(), 0..1, &states[0], &mut bar, &mut cost, &mut count, &mut || false).unwrap();
    assert_eq!(got.value.to_bits(), cost.to_bits()); assert_eq!(got.gradient.initial, bar.initial); assert_eq!(got.gradient.parameters, bar.parameters);
    let frozen = |x0: f64, rate: f64| {
        let perturbed = Decay { rate: Cell::new(rate), calls: Cell::new(0) };
        let mut u = vec![x0]; let residual = x0 - 0.5; let mut cost = 0.5 * residual * residual + 0.2 * rate;
        for r in &tape.records {
            u = replay(&perturbed, &u, r.start, r.end, r.h, &mut || false).unwrap().next;
            for &t in &times[1..] { if t == r.end { let e = u[0] - (0.5 + 0.1*t); cost += 0.5*e*e + 0.2*rate; } }
        }
        cost
    };
    let h = 1e-6;
    assert!((got.gradient.initial[0] - (frozen(1.2+h,-0.7)-frozen(1.2-h,-0.7))/(2.0*h)).abs() < 2e-8);
    assert!((got.gradient.parameters[0] - (frozen(1.2,-0.7+h)-frozen(1.2,-0.7-h))/(2.0*h)).abs() < 2e-8);
}

#[test]
fn observation_aligned_recording_resumes_at_every_attempt_boundary() {
    let model = Decay::new(); let times = [0.0, 0.05, 0.05, 0.18, 1.0, 2.0];
    let mut straight = record(&model, &times); straight.advance(10000,10000,&mut||false).unwrap();
    let mut split = record(&model, &times);
    assert_eq!(split.advance(10000,1,&mut||false).unwrap().status,RecordingStatus::RecordLimit);
    let mut split = split.clone();
    for _ in 0..10000 { if split.advance(1,10000,&mut||false).unwrap().status==RecordingStatus::ReachedEnd { break; } }
    assert_eq!(straight.records,split.records); assert_eq!(straight.state.u,split.state.u);
    assert_eq!(straight.state.h.to_bits(),split.state.h.to_bits()); assert_eq!(straight.state.err_prev.to_bits(),split.state.err_prev.to_bits());
    assert_eq!(straight.pullback_samples(&Loss::new(),budget(),&mut||false).unwrap(),split.pullback_samples(&Loss::new(),budget(),&mut||false).unwrap());
}

#[test]
fn observation_failure_cancellation_and_budget_leave_no_partial_gradient() {
    let model=Decay::new();let mut tape=record(&model,&[0.0,0.5,2.0]);
    assert!(matches!(tape.pullback_samples(&Loss::new(),budget(),&mut||false),Err(TrajectoryError::Incomplete)));
    tape.advance(10000,10000,&mut||false).unwrap();let before=tape.records.clone();let objective=Loss::new();
    objective.fail.set(true);
    assert_eq!(tape.pullback_samples(&objective,budget(),&mut||false),Err(TrajectoryError::Observation("sensor unavailable".into())));
    objective.fail.set(false);let expected=tape.pullback_samples(&objective,budget(),&mut||false).unwrap();
    let short=ReplayBudget{replayed_steps:expected.gradient.replayed_steps-1,..budget()};
    assert_eq!(tape.pullback_samples(&objective,short,&mut||false),Err(TrajectoryError::ReplayLimit));
    objective.seen.borrow_mut().clear();
    assert!(matches!(tape.pullback_samples(&objective,budget(),&mut||!objective.seen.borrow().is_empty()),Err(TrajectoryError::Step(AdjointError::Cancelled))));
    assert_eq!(tape.records,before);assert_eq!(tape.pullback_samples(&objective,budget(),&mut||false).unwrap(),expected);
}

#[test]
fn rejects_missing_partials_and_nonfinite_objective_values() {
    struct Broken(bool);
    impl SampleObjective for Broken {
        fn evaluate(&self,_:usize,_:f64,_:&[f64],x:&mut[f64],p:&mut[f64])->Result<f64,String> {
            x.fill(0.0); if self.0 {p.fill(0.0);Ok(f64::INFINITY)} else {Ok(0.0)}
        }
    }
    let model=Decay::new();let mut tape=record(&model,&[0.5,2.0]);tape.advance(10000,10000,&mut||false).unwrap();
    for flag in [false,true] {assert!(matches!(tape.pullback_samples(&Broken(flag),budget(),&mut||false),Err(TrajectoryError::Observation(_))));}
}

#[test]
fn validates_timetable_before_rhs_and_handles_zero_length_interval() {
    let model=Decay::new();
    for times in [vec![f64::NAN],vec![-0.1],vec![2.1],vec![0.5,0.4]] {
        assert!(RecordedRk45::new_sampled(&model,AdaptiveState::new(0.0,&[1.2],0.1),config(2.0),&times,10).is_err());
    }
    assert!(RecordedRk45::new_sampled(&model,AdaptiveState::new(0.0,&[1.2],0.1),config(2.0),&[0.5,1.0],1).is_err());
    assert_eq!(model.calls.get(),0);
    let tape=RecordedRk45::new_sampled(&model,AdaptiveState::new(2.0,&[1.2],0.1),config(2.0),&[2.0,2.0],2).unwrap();
    let got=tape.pullback_samples(&Loss::new(),ReplayBudget{checkpoints:0,replayed_steps:0},&mut||false).unwrap();
    assert_eq!(got.observations,2);assert!((got.gradient.initial[0]-1.0).abs()<1e-15);assert_eq!(got.gradient.parameters,vec![0.4]);
    assert_eq!(model.calls.get(),0);
}

#[test]
fn empty_timetable_preserves_the_existing_terminal_path() {
    let model=Decay::new();let mut sampled=record(&model,&[]);
    let mut ordinary=RecordedRk45::new(&model,AdaptiveState::new(0.0,&[1.2],1.0),config(2.0)).unwrap();
    sampled.advance(10000,10000,&mut||false).unwrap();ordinary.advance(10000,10000,&mut||false).unwrap();
    assert_eq!(sampled.records,ordinary.records);
    assert_eq!(sampled.pullback(&[1.0],&[0.0],budget(),&mut||false).unwrap(),ordinary.pullback(&[1.0],&[0.0],budget(),&mut||false).unwrap());
    assert!(sampled.pullback_samples(&Loss::new(),budget(),&mut||false).is_err());
}
