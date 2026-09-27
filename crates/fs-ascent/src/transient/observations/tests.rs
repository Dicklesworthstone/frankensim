use super::*;
use crate::transient::{evaluate_transient, TransientConfig, TransientStudy};
use crate::sqp::SqpStop;
use fs_time::PiController;
use fs_time::adaptive::adjoint::trajectory::{RecordingConfig, ReplayBudget};

pub(crate) const BOUNDS: [[f64; 2]; 2] = [[0.1, 2.0], [-1.0, 1.0]];
// Different experiments use different initial amplitudes. Rate and sensor bias
// are shared decision coordinates; a second channel observes twice the state.
pub(crate) struct DecayFamily { pub amplitude: f64 }
pub(crate) struct DecayModel { rate: f64, bias: f64, initial: [f64; 1] }
impl SensorFamily for DecayFamily {
    type Model = DecayModel;
    fn bounds(&self) -> &[[f64; 2]] { &BOUNDS }
    fn instantiate(&self, p: &[f64]) -> Result<DecayModel, String> {
        Ok(DecayModel { rate: p[0], bias: p[1], initial: [self.amplitude] })
    }
}
impl OdeVjp for DecayModel {
    fn dimension(&self) -> usize { 1 }
    fn parameter_count(&self) -> usize { 2 }
    fn rhs(&self, _: f64, state: &[f64], out: &mut [f64]) { out[0] = -self.rate * state[0]; }
    fn rhs_vjp(&self, _: f64, state: &[f64], seed: &[f64], x: &mut [f64], p: &mut [f64]) -> Result<(), String> {
        x[0] = -self.rate * seed[0]; p[0] = -state[0] * seed[0]; p[1] = 0.0; Ok(())
    }
}
impl SensorModel for DecayModel {
    fn initial_values(&self) -> &[f64] { &self.initial }
    fn initial_vjp(&self, _: &[f64], p: &mut [f64]) -> Result<(), String> { p.fill(0.0); Ok(()) }
    fn predict(&self, channel: u64, _: f64, x: &[f64]) -> Result<f64, String> {
        if channel > 1 { return Err("unknown channel".into()); }
        Ok((channel + 1) as f64 * x[0] + self.bias)
    }
    fn prediction_vjp(&self, channel: u64, _: f64, _: &[f64], seed: f64, x: &mut [f64], p: &mut [f64]) -> Result<(), String> {
        x[0] = (channel + 1) as f64 * seed; p[0] = 0.0; p[1] = seed; Ok(())
    }
}
pub(crate) fn config(end: f64) -> TransientConfig {
    TransientConfig { start: 0.0, initial_step: 0.1,
        recording: RecordingConfig { end, rtol: 1e-10, atol: 1e-12, controller: PiController::default(), max_workspace_components: 4096 },
        max_state_components: 1, max_samples: 64, max_attempts: 10000, max_records: 10000,
        replay: ReplayBudget { checkpoints: 64, replayed_steps: 100000 }, max_kkt_dimension: 6 }
}

#[test]
fn quadratic_noise_scaling_and_huber_branches_have_matched_derivatives() {
    let q = SensorReading::new(0.0, 0, 2.0, 0.5, SensorLoss::Quadratic).unwrap();
    assert_eq!(q.score(3.0).unwrap(), (2.0, 4.0));
    let h = SensorReading::new(0.0, 0, 2.0, 0.5, SensorLoss::Huber { threshold: 1.5 }).unwrap();
    assert_eq!(h.score(3.0).unwrap(), (1.875, 3.0));
    assert_eq!(h.score(1.0).unwrap(), (1.875, -3.0));
    for reading in [q, h] {
        for pred in [1.0, 1.6, 2.0, 2.4, 3.0] {
            let step = 1e-6;
            let fd = (reading.score(pred + step).unwrap().0 - reading.score(pred - step).unwrap().0) / (2.0 * step);
            assert!((fd - reading.score(pred).unwrap().1).abs() < 2e-9);
        }
    }
}

#[test]
fn noisy_multi_channel_gradients_include_direct_sensor_bias() {
    let rows = [(0.0, 0, 1.0, 0.4), (0.2, 0, 0.8, 0.3), (0.2, 1, 1.8, 0.7), (1.0, 1, 1.5, 0.2)]
        .into_iter().map(|(t,c,y,s)| SensorReading::new(t,c,y,s,SensorLoss::Huber{threshold:1.5}).unwrap()).collect::<Vec<_>>();
    let f = ObservedFamily::new(DecayFamily { amplitude: 1.2 }, SensorData::new(&rows, 4).unwrap());
    let p = [0.7, 0.15]; let cfg = config(1.0);
    let got = evaluate_transient(&f, &cfg, &p, &mut || false).unwrap().unwrap();
    let (mut value, mut rate, mut bias) = (0.0, 0.0, 0.0);
    for r in &rows {
        let x = 1.2 * (-p[0] * r.time()).exp(); let gain = (r.channel()+1) as f64;
        let (v, b) = r.score(gain*x+p[1]).unwrap(); value += v; rate += b*gain*(-r.time())*x; bias += b;
    }
    assert!((got.value-value).abs() < 1e-8);
    assert!((got.gradient[0]-rate).abs() < 1e-8);
    assert!((got.gradient[1]-bias).abs() < 1e-8);
    assert_eq!(got.observations, 4);
}

#[test]
fn robust_loss_limits_the_effect_of_one_bad_reading_in_a_fit() {
    let fit = |loss| {
        let rows = (0..25).map(|i| {
            let t = i as f64 / 8.0;
            let y = 1.2 * (-0.7*t).exp() + 0.15 + if i == 12 { 3.0 } else { 0.0 };
            SensorReading::new(t,0,y,0.1,loss).unwrap()
        }).collect::<Vec<_>>();
        let f = ObservedFamily::new(DecayFamily{amplitude:1.2},SensorData::new(&rows,25).unwrap());
        let mut study = TransientStudy::new(&f,&[0.9,0.0],config(3.0),&mut||false).unwrap();
        let report = study.run(1e-6,100,1500,&mut||false).unwrap();
        assert_eq!(report.stop,SqpStop::Converged,"{report:?}");
        let p = study.optimizer().point(); (p[0]-0.7).abs()+(p[1]-0.15).abs()
    };
    assert!(fit(SensorLoss::Huber{threshold:1.5}) < 0.2*fit(SensorLoss::Quadratic));
}

#[test]
fn invalid_data_and_sensor_bindings_refuse_instead_of_reusing_outputs() {
    for sigma in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        assert!(SensorReading::new(0.0,0,1.0,sigma,SensorLoss::Quadratic).is_err());
    }
    assert!(SensorReading::new(0.0,0,1.0,1.0,SensorLoss::Huber{threshold:0.0}).is_err());
    let row = SensorReading::new(0.2,0,1.0,1.0,SensorLoss::Quadratic).unwrap();
    assert!(SensorData::new(&[row.clone()],0).is_err());
    assert!(SensorData::new(&[],1).is_err());
    let early = SensorReading::new(0.1,0,1.0,1.0,SensorLoss::Quadratic).unwrap();
    assert!(SensorData::new(&[row.clone(),early],2).is_err());
    let f = ObservedFamily::new(DecayFamily{amplitude:1.2},SensorData::new(&[row],1).unwrap());
    let m = f.instantiate(&[0.7,0.15]).unwrap();
    assert!(m.evaluate(0,0.3,&[1.0],&mut[0.0],&mut[0.0;2]).is_err());
    assert!(m.evaluate(1,0.2,&[1.0],&mut[0.0],&mut[0.0;2]).is_err());
    assert!(m.evaluate(0,0.2,&[1.0],&mut[],&mut[0.0;2]).is_err());
    assert!(m.data.readings[0].score(f64::NAN).is_err());
}
