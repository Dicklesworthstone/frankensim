use super::*;
use std::cell::{Cell, RefCell};

struct Decay {
    rate: f64,
    order: RefCell<Vec<usize>>,
    fault: Cell<u8>,
    changed: Cell<bool>,
    loss_calls: Cell<usize>,
}
impl Decay {
    fn new(rate: f64) -> Self {
        Self { rate, order: RefCell::new(Vec::new()), fault: Cell::new(0),
            changed: Cell::new(false), loss_calls: Cell::new(0) }
    }
}
impl OdeVjp for Decay {
    fn dimension(&self) -> usize { 1 }
    fn parameter_count(&self) -> usize { 1 }
    fn rhs(&self, _: f64, state: &[f64], out: &mut [f64]) { out[0] = self.rate * state[0]; }
    fn rhs_vjp(&self, _: f64, state: &[f64], seed: &[f64], x: &mut [f64], p: &mut [f64]) -> Result<(), String> {
        x[0] = self.rate * seed[0]; p[0] = state[0] * seed[0]; Ok(())
    }
}
impl JointSampleObjective for Decay {
    fn observe(&self, sample: usize, _: f64, state: &[f64]) -> Result<f64, String> {
        self.order.borrow_mut().push(sample);
        if self.fault.get() == 1 { return Ok(f64::NAN); }
        Ok(state[0] + 0.5 * self.rate * self.rate + if self.changed.get() { 0.01 } else { 0.0 })
    }
    fn loss<C: FnMut() -> bool>(&self, y: &[f64], seeds: &mut [f64], direct: &mut [f64], _: &mut C)
        -> Result<f64, TrajectoryError>
    {
        self.loss_calls.set(self.loss_calls.get() + 1);
        // C = I + 2*11^T. This independent oracle uses its explicit inverse.
        let sum: f64 = y.iter().map(|v| v - 0.4).sum();
        let shift = 2.0 * sum / (1.0 + 2.0 * y.len() as f64);
        let mut value = 0.25 * self.rate * self.rate;
        for (i, (seed, value_i)) in seeds.iter_mut().zip(y).enumerate() {
            *seed = value_i - 0.4 - shift;
            value += 0.5 * (value_i - 0.4) * *seed;
            if i == 0 && self.fault.get() == 2 { *seed = f64::NAN; }
        }
        direct[0] = 0.5 * self.rate;
        if self.fault.get() == 3 { return Err(TrajectoryError::Observation("joint loss refused".into())); }
        if self.fault.get() == 4 { self.changed.set(true); }
        Ok(value)
    }
    fn observation_vjp(&self, _: usize, _: f64, _: &[f64], seed: f64, x: &mut [f64], p: &mut [f64])
        -> Result<(), String>
    {
        if self.fault.get() == 5 { return Err("sensor derivative refused".into()); }
        x[0] = seed; p[0] = self.rate * seed; Ok(())
    }
}
fn budget() -> ReplayBudget { ReplayBudget { checkpoints: 64, replayed_steps: 100_000 } }
fn tape<'a>(model: &'a Decay, times: &[f64], initial: f64) -> RecordedRk45<'a, Decay> {
    let cfg = RecordingConfig { end: *times.last().unwrap(), rtol: 1e-11, atol: 1e-13,
        controller: PiController::default(), max_workspace_components: 1024 };
    let mut result = RecordedRk45::new_sampled(model, AdaptiveState::new(0.0, &[initial], 0.4), cfg, times, 64).unwrap();
    result.advance(10_000, 10_000, &mut || false).unwrap(); result
}
fn analytic(times: &[f64], initial: f64, rate: f64) -> (f64, f64, f64) {
    let y: Vec<_> = times.iter().map(|t| initial * (rate*t).exp() + 0.5*rate*rate).collect();
    let sum: f64 = y.iter().map(|v| v-0.4).sum(); let shift = 2.0*sum/(1.0+2.0*y.len() as f64);
    let (mut value, mut x, mut p) = (0.25*rate*rate, 0.0, 0.5*rate);
    for (&t, &prediction) in times.iter().zip(&y) {
        let seed = prediction-0.4-shift;
        value += 0.5*(prediction-0.4)*seed;
        x += seed*(rate*t).exp(); p += seed*(t*initial*(rate*t).exp()+rate);
    }
    (value, x, p)
}

#[test]
fn cross_time_joint_gradient_includes_sensor_and_direct_parameter_terms() {
    let times = [0.0, 0.17, 0.17, 0.8, 1.5]; let model = Decay::new(-0.7);
    let recording = tape(&model, &times, 1.2);
    let got = recording.pullback_joint(&model, budget(), 64, &mut || false).unwrap();
    let expected = analytic(&times, 1.2, -0.7);
    assert!((got.value-expected.0).abs()<2e-10);
    assert!((got.gradient.initial[0]-expected.1).abs()<2e-9);
    assert!((got.gradient.parameters[0]-expected.2).abs()<2e-9);
    assert_eq!(got.observations, times.len()); assert_eq!(model.loss_calls.get(), 1);
    let expected_order: Vec<_> = (0..times.len()).chain((0..times.len()).rev()).collect();
    assert_eq!(*model.order.borrow(), expected_order);
}

#[test]
fn finite_difference_of_full_joint_simulation_matches_pullback() {
    let times = [0.0, 0.19, 0.53, 0.95]; let model = Decay::new(-0.4);
    let got = tape(&model, &times, 0.9).pullback_joint(&model, budget(), 64, &mut || false).unwrap();
    let value = |initial: f64, rate: f64| {
        let m = Decay::new(rate);
        tape(&m, &times, initial).pullback_joint(&m, budget(), 64, &mut || false).unwrap().value
    };
    let h=1e-5;
    assert!((got.gradient.initial[0]-(value(0.9+h,-0.4)-value(0.9-h,-0.4))/(2.0*h)).abs()<2e-8);
    assert!((got.gradient.parameters[0]-(value(0.9,-0.4+h)-value(0.9,-0.4-h))/(2.0*h)).abs()<2e-8);
}

#[test]
fn collection_and_reverse_share_one_replay_allowance() {
    let model=Decay::new(-0.7); let recording=tape(&model,&[0.0,0.4,1.0],1.2);
    let expected=recording.pullback_joint(&model,budget(),64,&mut||false).unwrap();
    // Obtain the old sweep's work count with a zero objective, not the joint callback.
    struct Zero;
    impl SampleObjective for Zero {
        fn evaluate(&self,_:usize,_:f64,_:&[f64],x:&mut[f64],p:&mut[f64])->Result<f64,String> {
            x.fill(0.0);p.fill(0.0);Ok(0.0)
        }
    }
    let old=recording.pullback_samples(&Zero,budget(),&mut||false).unwrap();
    assert_eq!(expected.gradient.replayed_steps,old.gradient.replayed_steps+recording.accepted_steps());
    let exact=ReplayBudget {replayed_steps:expected.gradient.replayed_steps,..budget()};
    assert_eq!(recording.pullback_joint(&model,exact,64,&mut||false).unwrap(),expected);
    let short=ReplayBudget {replayed_steps:exact.replayed_steps-1,..exact};
    assert_eq!(recording.pullback_joint(&model,short,64,&mut||false),Err(TrajectoryError::ReplayLimit));
}

#[test]
fn failures_in_either_pass_never_publish_a_partial_gradient() {
    let model=Decay::new(-0.7);let recording=tape(&model,&[0.0,0.4,1.0],1.2);
    let before=recording.records.clone();
    for fault in 1..=5 {
        model.fault.set(fault);model.changed.set(false);
        assert!(matches!(recording.pullback_joint(&model,budget(),64,&mut||false),Err(TrajectoryError::Observation(_))));
        assert_eq!(recording.records,before);
    }
    model.fault.set(0);model.changed.set(false);
    assert!(recording.pullback_joint(&model,budget(),64,&mut||false).is_ok());
}

#[test]
fn cancellation_and_invalid_caps_preserve_recording() {
    let model=Decay::new(-0.7);let recording=tape(&model,&[0.0,0.4,1.0],1.2);
    assert!(recording.pullback_joint(&model,budget(),2,&mut||false).is_err());
    assert!(model.order.borrow().is_empty());
    let mut polls=0;
    let expected=recording.pullback_joint(&model,budget(),64,&mut||{polls+=1;false}).unwrap();
    for at in [1,7,polls/2,polls] {
        let mut current=0;
        assert_eq!(recording.pullback_joint(&model,budget(),64,&mut||{current+=1;current==at}),
            Err(TrajectoryError::Step(AdjointError::Cancelled)));
    }
    assert_eq!(recording.pullback_joint(&model,budget(),64,&mut||false).unwrap(),expected);
}

#[test]
fn initial_only_joint_objective_needs_no_replays() {
    let model=Decay::new(-0.7);let times=[0.0,0.0];let recording=tape(&model,&times,1.2);
    let got=recording.pullback_joint(&model,ReplayBudget {checkpoints:0,replayed_steps:0},2,&mut||false).unwrap();
    let expected=analytic(&times,1.2,-0.7);
    assert!((got.value-expected.0).abs()<1e-14);
    assert!((got.gradient.initial[0]-expected.1).abs()<1e-14);
    assert!((got.gradient.parameters[0]-expected.2).abs()<1e-14);
    assert_eq!(got.gradient.replayed_steps,0);assert_eq!(got.gradient.peak_checkpoints,0);
}
