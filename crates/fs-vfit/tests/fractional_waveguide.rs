use fs_vfit::waveguide::{PassiveWaveguide, WaveguideSpec};
use fs_vfit::waveguide::fractional::FractionalTransit;

fn spec(n: usize, reflection: f64) -> WaveguideSpec {
    WaveguideSpec { one_way_samples: n, impedance_pa_s_m3: 1e6,
        time_step_s: 1.0/48000.0, reflection, max_memory_bytes: 1<<20 }
}

#[test]
fn zero_fraction_retains_every_original_integer_frame() {
    for reflection in [-1.0, -0.7, 0.0, 1.0] {
        let mut old=PassiveWaveguide::new(spec(7,reflection)).unwrap();
        let mut new=PassiveWaveguide::with_fractional_transit(spec(7,reflection),
            FractionalTransit::new(0.0).unwrap()).unwrap();
        for n in 0..1000 {
            let x=if n<900 { f64::from(n%31)-15.0 } else {0.0};
            assert_eq!(old.step(x).unwrap(),new.step(x).unwrap());
        }
        assert!(new.fractional_transit().is_none());
    }
}

#[test]
fn fractional_impulse_matches_independent_allpass_difference_equation() {
    for delay in [0.001,0.05,0.25,0.5,0.9,0.999] {
        let transit=FractionalTransit::new(delay).unwrap();
        let mut line=PassiveWaveguide::with_fractional_transit(spec(7,0.0),transit).unwrap();
        let a=(1.0-delay)/(1.0+delay);
        let (mut x1,mut y1)=(0.0,0.0);
        for n in 0..2048 {
            let x=if n==0 {1.0} else {0.0};
            let delayed_x=if n==7 {1.0} else {0.0};
            let expected=a*delayed_x+x1-a*y1;
            x1=delayed_x;y1=expected;
            let frame=line.step(x).unwrap();
            assert!((frame.terminal_pressure_pa-expected).abs()<2e-13);
            if n<7 {assert_eq!(frame.terminal_pressure_pa,0.0);}
        }
    }
}

#[test]
fn both_fractional_memories_store_energy_and_do_not_become_fake_damping() {
    for delay in [0.01,0.2,0.5,0.99] {
        for reflection in [-1.0,-0.7,0.0,1.0] {
            let s=spec(3,reflection);
            let mut line=PassiveWaveguide::with_fractional_transit(s,FractionalTransit::new(delay).unwrap()).unwrap();
            let mut delivered=0.0;let mut escaped=0.0;
            for n in 0..4096 {
                let old=line.stored_energy_j();
                let x=if n<128 {(n as f64*0.7).sin()} else {0.0};
                let frame=line.step(x).unwrap();
                delivered+=frame.inlet_work_j;escaped+=frame.terminal_loss_j;
                let scale=old+frame.stored_energy_j+frame.inlet_work_j.abs()+frame.terminal_loss_j;
                assert!(frame.balance_residual_j().abs()<=3e-13*scale.max(1e-300));
                assert!(frame.terminal_loss_j>=0.0);
                let memory=line.fractional_state_pa().iter().map(|s|s*s).sum::<f64>()*s.time_step_s/s.impedance_pa_s_m3;
                assert!(memory<=frame.stored_energy_j+1e-25);
                assert_eq!(line.stored_energy_j(),frame.stored_energy_j);
            }
            assert!((delivered-escaped-line.stored_energy_j()).abs()<1e-22);
        }
    }
}

#[test]
fn measured_passband_has_unit_gain_and_the_declared_phase() {
    let delay=0.37;
    let transit=FractionalTransit::new(delay).unwrap();
    for bin in [1,17,61,101] {
        let omega=core::f64::consts::TAU*f64::from(bin)/512.0;
        let phase=7.0*omega+2.0*(transit.realized_dc_samples()*(omega/2.0).tan()).atan();
        let mut line=PassiveWaveguide::with_fractional_transit(spec(7,0.0),transit).unwrap();
        let (mut re,mut im)=(0.0,0.0);
        for n in 0..2560 {
            let angle=f64::from(n)*omega;
            let y=line.step(angle.sin()).unwrap().terminal_pressure_pa;
            if n>=2048 {re+=y*angle.sin();im+=y*angle.cos();}
        }
        re/=256.0;im/=256.0;
        assert!((re-phase.cos()).abs()<1e-11);
        assert!((im+phase.sin()).abs()<1e-11);
        assert!((re.hypot(im)-1.0).abs()<1e-11);
        let discrepancy=(phase-(7.0+delay)*omega).abs();
        assert!(discrepancy<=transit.phase_error_bound_rad(omega).unwrap()+2e-14);
    }
    // A linear fractional interpolator would quietly remove high-frequency energy.
    let linear_gain=((1.0-delay).powi(2)+delay*delay).sqrt();
    assert!(linear_gain<0.8);
}

#[test]
fn fractional_refusal_preserves_return_waves_and_memory() {
    let mut line=PassiveWaveguide::with_fractional_transit(spec(3,-0.8),FractionalTransit::new(0.3).unwrap()).unwrap();
    for n in 0..30 {line.step(f64::from(n%7)-3.0).unwrap();}
    let mut reference=line.clone();
    for bad in [f64::NAN,f64::INFINITY,f64::MAX] {
        assert!(line.step(bad).is_err());
        assert_eq!(line.fractional_state_pa(),reference.fractional_state_pa());
        assert_eq!(line.stored_energy_j(),reference.stored_energy_j());
    }
    for n in 0..100 {assert_eq!(line.step(f64::from(n%3)).unwrap(),reference.step(f64::from(n%3)).unwrap());}
    for bad in [-0.1,1.0,f64::NAN,f64::INFINITY,f64::MIN_POSITIVE] {
        assert!(FractionalTransit::new(bad).is_err());
    }
    assert!(FractionalTransit::new(0.3).unwrap().phase_error_bound_rad(core::f64::consts::PI).is_err());
}
