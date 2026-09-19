use super::*;
use crate::broadband_radiation::{ComplexShTrainingSample, DirectFarFieldHeldOutSample,
    HarmonicTimeConvention, RadiationSampleDiagnostics, WeightPreset};
use super::super::{ImpactBody, ImpactConfig};

const RATE: u32 = 48_000;
fn samples() -> SampledRadiationData {
    // Analytic pulsating sphere, a=0.08m: F/v=-i*w*rho*a^2
    // exp(-i*k*a)/(1-i*k*a). This is NOT an authored drum transfer.
    let transfer = |f: f64| {
        let w = core::f64::consts::TAU*f;
        C64::new(0.0,-w*1.2*0.08*0.08)
            * C64::new(1.0,-w*0.08/343.0).recip()
            * C64::new(det::cos(w*0.08/343.0),-det::sin(w*0.08/343.0))
    };
    let diag = RadiationSampleDiagnostics { captured_fraction:1.0,
        panels_per_wavelength:100.0, condition_lower_bound:1.0 };
    let root = (4.0*core::f64::consts::PI).sqrt();
    SampledRadiationData { source_id:"analytic pulsating sphere; no mesh certificate".into(),
        harmonic_time_convention:HarmonicTimeConvention::ExpNegativeIOmegaT,
        l_max:0,input_ids:vec!["surface_velocity_m_s".into()],
        training:(0..24).map(|i| {let hz=30.0+f64::from(i)*80.0;
            ComplexShTrainingSample{omega_rad_s:core::f64::consts::TAU*hz,
                coefficients_by_input:vec![vec![transfer(hz).scale(root)]],diagnostics:diag}}).collect(),
        held_out:(0..23).map(|i| {let hz=70.0+f64::from(i)*80.0;
            DirectFarFieldHeldOutSample{omega_rad_s:core::f64::consts::TAU*hz,
                directions:vec![[0.0,0.0,1.0]],far_field_by_input:vec![vec![transfer(hz)]],diagnostics:diag}}).collect() }
}
fn controls() -> BroadbandRadiationControls {
    BroadbandRadiationControls { sample_rate_hz:f64::from(RATE),minimum_captured_fraction:0.99,
        fit_order:1,fit_iterations:8,fit_weights:WeightPreset::Uniform,fit_d:true,
        far_field_signal_floor:1e-10,maximum_normalized_error:0.03,rms_normalized_error:0.02 }
}
fn artifact() -> ImpactRadiation {
    ImpactRadiation::from_velocity_samples(&samples(),controls(),0.08,343.0).unwrap()
}
fn system(steps:u64) -> ImpactSystem {
    let (body,_) = ImpactBody::free_mass(1.0,0.0,0.0).unwrap();
    ImpactSystem::new(vec![body],vec![],vec![],vec![],ImpactConfig{dt_s:1.0/f64::from(RATE),
        max_steps:steps,maximum_energy_j:10.0,energy_absolute_tolerance_j:1e-12,
        energy_relative_tolerance:1e-8,maximum_generalized_force:1000.0}).unwrap()
}
fn listener(range:f64) -> ImpactListener {
    ImpactListener{position_m:[0.0,0.0,range],maximum_delay_error_s:1.01/f64::from(RATE),
        maximum_delay_samples:4096,maximum_abs_pressure_pa:1000.0}
}
fn renderer(a:&ImpactRadiation,range:f64,steps:u64)->ImpactPressureRenderer<'_> {
    ImpactPressureRenderer::new(system(steps),a,vec![VelocityProjection{
        input_id:"surface_velocity_m_s".into(),weights:vec![1.0]}],vec![0.0],listener(range),RATE,97).unwrap()
}
#[test]
fn sphere_phase_shift_and_acceleration_units_match_an_independent_impedance() {
    let a=artifact();
    assert!((a.source_time_shift_s()-0.08/343.0).abs()<1e-15);
    for row in samples().held_out {
        let fitted=a.bank.inputs[0].filters[0].eval(row.omega_rad_s).unwrap().conj()
            .scale(1.0/(4.0*core::f64::consts::PI).sqrt());
        let expected=C64::from_re(1.2*0.08*0.08)
            * C64::new(1.0,-row.omega_rad_s*0.08/343.0).recip();
        assert!((fitted-expected).abs()<0.03*expected.abs());
        // Negative control: retaining the sphere's -a/c phase fails after
        // the conversion even though that wrong transfer has equal magnitude.
        if row.omega_rad_s>3000.0 {
            let angle=row.omega_rad_s*0.08/343.0;
            let wrong=expected*C64::new(angle.cos(),-angle.sin());
            assert!((wrong-expected).abs()>0.5*expected.abs());
        }
    }
    let mut corrupt=samples();
    for row in &mut corrupt.held_out {row.far_field_by_input[0][0]=row.far_field_by_input[0][0].scale(-1.0);}
    assert!(ImpactRadiation::from_velocity_samples(&corrupt,controls(),0.08,343.0).is_err());
}
#[test]
fn pressure_has_real_travel_time_range_scaling_and_no_steady_velocity_tone() {
    let a=artifact();let near_range=0.08+343.0*32.0/f64::from(RATE);
    let far_range=near_range+343.0*32.0/f64::from(RATE);
    let mut near=renderer(&a,near_range,1000);let mut far=renderer(&a,far_range,1000);
    let delta=far.timing().propagation_samples-near.timing().propagation_samples;
    assert!((31..=33).contains(&delta));
    let mut np=Vec::new();let mut fp=Vec::new();
    for i in 0..1000 {
        let force=if i<16 {100.0}else{0.0};near.set_forces(&[force]).unwrap();far.set_forces(&[force]).unwrap();
        let(mut x,mut y)=([0.0],[0.0]);near.block(&mut x).unwrap();far.block(&mut y).unwrap();
        np.push(x[0]);fp.push(y[0]);
    }
    assert!(np[..near.timing().propagation_samples].iter().all(|p|*p==0.0));
    assert!(np.iter().any(|p|p.abs()>0.01));
    for i in delta..1000 {assert!((fp[i]*far_range-np[i-delta]*near_range).abs()<1e-12);}
    assert!(np[900..].iter().all(|p|p.abs()<1e-10));
    assert!(near.mechanics().state()[1]>0.0); // still coasting: no synthetic velocity hum
}
#[test]
fn callback_partitions_and_pause_resume_retain_all_physical_histories() {
    let a=artifact();let mut baseline=renderer(&a,1.0,600);baseline.set_forces(&[40.0]).unwrap();
    let mut expected=vec![0.0;600];for c in expected.chunks_mut(97) {baseline.block(c).unwrap();}
    for size in [1,7,32,97] {
        let mut actual=renderer(&a,1.0,600);actual.set_forces(&[40.0]).unwrap();
        let mut out=vec![0.0;600];for c in out.chunks_mut(size) {actual.block(c).unwrap();}
        assert_eq!(out.iter().map(|p|p.to_bits()).collect::<Vec<_>>(),
            expected.iter().map(|p|p.to_bits()).collect::<Vec<_>>());
        assert_eq!(actual.mechanics().state(),baseline.mechanics().state());
    }
}
#[test]
fn admission_and_atomic_controls_never_erase_accepted_sound() {
    let a=artifact();let mut r=renderer(&a,1.0,10);
    assert!(r.validate_sample_rate(44100).is_err());assert!(r.validate_sample_count(11).is_err());
    r.set_forces(&[10.0]).unwrap();assert!(r.set_forces(&[f64::NAN]).is_err());assert!(r.set_forces(&[]).is_err());
    assert_eq!(r.forces,[10.0]);assert_eq!(r.samples_rendered(),0);
    assert!(r.block(&mut []).is_err());r.block(&mut [0.0;10]).unwrap();
    assert!(r.block(&mut [0.0]).is_err());assert_eq!(r.samples_rendered(),10);
    let p=vec![VelocityProjection{input_id:"wrong_basis".into(),weights:vec![1.0]}];
    assert!(ImpactPressureRenderer::new(system(1),&a,p,vec![0.0],listener(1.0),RATE,1).is_err());
    let p=||vec![VelocityProjection{input_id:"surface_velocity_m_s".into(),weights:vec![1.0]}];
    assert!(ImpactPressureRenderer::new(system(1),&a,p(),vec![0.0],listener(0.07),RATE,1).is_err());
    let mut capped=listener(1.0);capped.maximum_delay_samples=1;
    assert!(ImpactPressureRenderer::new(system(1),&a,p(),vec![0.0],capped,RATE,1).is_err());
}
#[test]
fn failed_observation_poisons_the_host_instead_of_replaying_a_partial_callback() {
    let a=artifact();let mut r=renderer(&a,1.0,10);
    r.pressure_limit=1e-15;r.set_forces(&[100.0]).unwrap();
    assert!(r.block(&mut [0.0;4]).is_err());assert!(matches!(r.block(&mut [0.0]),Err(RenderError::Poisoned)));
    assert!(r.set_forces(&[0.0]).is_err());
}
