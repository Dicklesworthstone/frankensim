//! End-to-end pressure against analytic acceleration, not against the old
//! discrete difference. Legacy propagation tests remain separate and unchanged.
use super::*;
use std::f64::consts::TAU;

fn medium() -> RayleighMedium { RayleighMedium { density: 1.2, sound_speed: 320.0 } }
fn patch(x: f64, shape: Vec<f64>) -> SurfaceSample {
    SurfaceSample { position_m: [x,0.0,0.0], area_m2: 0.1, mode_shape: shape }
}

#[test]
fn pressure_matches_continuous_acceleration_amplitude_and_phase_through_the_audio_band() {
    let rate = 48_000;
    let travel = 75.5;
    let position = [0.0,0.0,travel*320.0/f64::from(rate)];
    let gain = 1.2*0.1/(TAU*position[2]);
    for ratio in [1,3,4,16] {
        for frequency in [440.0,4800.0,15000.0,21600.0] {
            let mut mic = BaffledPressure::new_multirate_differentiating(
                &[patch(0.0,vec![1.0])],1,rate*ratio,rate,position,medium()).unwrap();
            let latency = travel+mic.filter_delay_samples();
            let mut trace = vec![0.0;ratio as usize];
            let (mut error,mut signal) = (0.0,0.0);
            for frame in 0..2400 {
                for (sub,value) in trace.iter_mut().enumerate() {
                    let t = (frame*ratio as usize+sub) as f64/f64::from(rate*ratio);
                    *value = (TAU*frequency*t).sin();
                }
                let pressure = mic.step_trace(&trace).unwrap();
                if frame>500 {
                    let t = ((frame+1)*ratio as usize-1) as f64/f64::from(rate*ratio);
                    let exact = gain*TAU*frequency*(TAU*frequency*(t-latency/f64::from(rate))).cos();
                    error += (pressure-exact).powi(2);
                    signal += exact.powi(2);
                }
            }
            assert!((error/signal).sqrt()<5e-4,
                "ratio {ratio}, frequency {frequency}: relative RMS {}",(error/signal).sqrt());
        }
    }
}

#[test]
fn derivative_keeps_the_highest_mode_si_gain_and_all_prepared_histories() {
    let count = MAX_RAYLEIGH_INPUTS;
    let mut shape = vec![0.0;count];shape[count-1]=1.0;
    let surface = [patch(0.0,shape)];
    let mut a = BaffledPressure::new_multirate_differentiating(
        &surface,count,128_000,32_000,[0.0,0.0,1.005],medium()).unwrap();
    let mut b = BaffledPressure::new_multirate_differentiating(
        &surface,count,128_000,32_000,[0.0,0.0,1.005],medium()).unwrap();
    assert_eq!(a.filter_delay_samples(),84.0);
    assert_eq!(a.interpolation_delay_samples(),32);
    let pointers = (a.history.as_ptr(),a.acceleration.as_ptr(),a.previous_velocity.as_ptr());
    let gain = 1.2*0.1/(TAU*1.005);
    let mut trace = vec![0.0;4*count];
    let mut components = vec![0.0;count];
    for frame in 0..500 {
        assert!(a.step_trace(&trace[..trace.len()-1]).is_err());
        trace[count-1]=f64::NAN;
        assert!(a.step_trace(&trace).is_err());
        for sub in 0..4 { trace[sub*count+count-1]=(frame*4+sub+1) as f64/128_000.0; }
        let pressure = a.step_trace(&trace).unwrap();
        assert_eq!(pressure,b.step_trace(&trace).unwrap());
        a.mode_pressures(&mut components).unwrap();
        assert!((pressure-components.iter().sum::<f64>()).abs()<1e-12);
        assert!(components[..count-1].iter().all(|v| *v==0.0));
        if frame>350 { assert!((pressure/gain-1.0).abs()<1e-9); }
    }
    assert_eq!(pointers,(a.history.as_ptr(),a.acceleration.as_ptr(),a.previous_velocity.as_ptr()));
}

#[test]
fn the_derivative_preserves_intrasample_attacks_and_signed_surface_cancellation() {
    let surface = [patch(0.0,vec![1.0])];
    let mut mic = BaffledPressure::new_multirate_differentiating(
        &surface,1,128_000,32_000,[0.0,0.0,1.0],medium()).unwrap();
    let mut peak = 0.0_f64;
    for frame in 0..500 {
        let trace = if frame==0 { [0.01,0.0,0.0,0.0] } else { [0.0;4] };
        let pressure = mic.step_trace(&trace).unwrap();
        if frame<mic.delay_samples.0 { assert_eq!(pressure,0.0); }
        peak=peak.max(pressure.abs());
    }
    assert!(peak>1e-5,"endpoint-zero substep attack was erased");
    let surface = [patch(-0.3,vec![1.0]),patch(0.3,vec![-1.0])];
    let mut axis = BaffledPressure::new_differentiating(
        &surface,1,48_000,[0.0,0.0,1.0],medium()).unwrap();
    let mut side = BaffledPressure::new_differentiating(
        &surface,1,48_000,[0.8,0.0,1.0],medium()).unwrap();
    assert_eq!(axis.filter_delay_samples(),80.0);
    let mut side_peak = 0.0_f64;
    for frame in 0..1000 {
        let velocity = [(frame as f64*0.07).sin()];
        assert!(axis.step(&velocity).unwrap().abs()<1e-12);
        side_peak=side_peak.max(side.step(&velocity).unwrap().abs());
    }
    assert!(side_peak>1.0);
}

#[test]
fn differentiating_frontdoor_preserves_complete_clock_and_geometry_admission() {
    let surface = [patch(0.0,vec![1.0])];
    for (input,output) in [(0,48_000),(48_000,0),(48_001,48_000),(816_000,48_000)] {
        assert!(BaffledPressure::new_multirate_differentiating(
            &surface,1,input,output,[0.0,0.0,1.0],medium()).is_err());
    }
    assert!(BaffledPressure::new_differentiating(
        &surface,2,48_000,[0.0,0.0,1.0],medium()).is_err());
    assert!(BaffledPressure::new_differentiating(
        &surface,1,48_000,[0.0,0.0,0.0],medium()).is_err());
}
