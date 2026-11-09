//! Original Maxwell startup, retained memory and transactional endpoint handoff.

use fs_fft::harmonic_window::{HarmonicWindow, HarmonicWindowProfile, phase_angular_interval_from_clock_dd};
use fs_material::visco::{Mandel6, PronyBranch, TensorProny, TensorPronyState, ViscoError};
use fs_math::{c64::C64, dd::Dd};

fn matrix(diagonal: f64, row: [f64; 6]) -> Mandel6 {
    core::array::from_fn(|i| core::array::from_fn(|j| row[i] * row[j] + if i == j { diagonal } else { 0.0 }))
}
fn law() -> TensorProny {
    TensorProny::new_mandel(matrix(7.0, [1.0,-0.5,0.2,0.4,-0.3,0.6]), vec![
        PronyBranch::new(matrix(0.0, [1.0,2.0,-0.5,0.8,-0.4,0.3]),0.2).unwrap(),
        PronyBranch::new(matrix(0.4, [-0.2,0.3,1.0,-0.6,0.7,0.8]),1.4).unwrap(),
    ],0.05).unwrap().with_kelvin_viscosity_mandel(matrix(0.01,[0.1,-0.2,0.05,0.08,-0.03,0.04])).unwrap()
        .with_band((0.0,100.0)).unwrap()
}
fn fields() -> [[C64; 6]; 2] {
    core::array::from_fn(|k| core::array::from_fn(|i| C64::new(
        (i as f64-2.0)*0.001/(k+1) as f64, (3.0-i as f64)*0.0007*(k+1) as f64,
    )))
}
fn apply(c: &Mandel6, x: &[f64; 6]) -> [f64; 6] {
    core::array::from_fn(|i| (0..6).map(|j| c[i][j]*x[j]).sum())
}
fn dot(x: &[f64; 6], y: &[f64; 6]) -> f64 { x.iter().zip(y).map(|(a,b)|a*b).sum() }
fn near(a: f64, b: f64, scale: f64, tolerance: f64) {
    assert!((a-b).abs()<=tolerance*scale, "actual={a:e} expected={b:e} scale={scale:e}");
}
fn waveform(frequencies: &[f64], fields: &[[C64; 6]], bias: &[f64; 6], t: f64) -> ([f64; 6],[f64; 6]) {
    let (mut strain,mut rate)=(*bias,[0.0;6]);
    for (&omega,field) in frequencies.iter().zip(fields) {
        let (s,c)=(omega*t).sin_cos();
        for i in 0..6 { strain[i]+=field[i].re*c+field[i].im*s; rate[i]+=omega*(field[i].im*c-field[i].re*s); }
    }
    (strain,rate)
}
fn stored(law:&TensorProny,strain:&[f64;6],q:&[[f64;6]]) -> f64 {
    0.5*dot(strain,&apply(law.c_eq(),strain)) + q.iter().zip(law.branches())
        .map(|(q,b)|0.5*dot(q,&apply(&b.c_pa,q))).sum::<f64>()
}
fn state(law:&TensorProny,strain:[f64;6],q:Vec<[f64;6]>) -> TensorPronyState {
    let energy=stored(law,&strain,&q);
    let mut stress=apply(law.c_eq(),&strain);
    for (branch,memory) in law.branches().iter().zip(&q) {
        let value=apply(&branch.c_pa,memory);
        for i in 0..6 { stress[i]+=value[i]; }
    }
    TensorPronyState {q,strain,stress,stored:energy,work:energy,dissipated:0.0}
}

// Independent real-time RK4 integration of the original branch ODE and the
// physical stress work/heat rates; no production temporal or Gram owner enters.
fn oracle(law:&TensorProny,frequencies:&[f64],fields:&[[C64;6]],bias:&[f64;6],initial:&TensorPronyState,t0:f64,h:f64)
    -> (Vec<[f64;6]>,[f64;3])
{
    let n=6*law.branches().len(); let mut y=vec![0.0;n+3];
    for (chunk,q) in y[..n].chunks_exact_mut(6).zip(&initial.q) { chunk.copy_from_slice(q); }
    let derivative=|t:f64,y:&[f64]| {
        let (eps,rate)=waveform(frequencies,fields,bias,t);
        let mut dy=vec![0.0;n+3]; let mut stress=apply(law.c_eq(),&eps);
        if let Some(kv)=law.kv_pa_s() {
            let viscous=apply(kv,&rate); dy[n]=dot(&viscous,&rate);
            for i in 0..6 { stress[i]+=viscous[i]; }
        }
        for (b,branch) in law.branches().iter().enumerate() {
            let q: [f64;6]=y[6*b..6*b+6].try_into().unwrap();
            let s=apply(&branch.c_pa,&q); dy[n+1]+=dot(&q,&s)/branch.tau_s;
            for i in 0..6 { stress[i]+=s[i]; dy[6*b+i]=rate[i]-q[i]/branch.tau_s; }
        }
        dy[n+2]=dot(&stress,&rate); dy
    };
    let steps=8192; let dt=h/steps as f64;
    for k in 0..steps {
        let t=t0+k as f64*dt; let a=derivative(t,&y);
        let z:Vec<_>=y.iter().zip(&a).map(|(x,d)|x+0.5*dt*d).collect(); let b=derivative(t+0.5*dt,&z);
        let z:Vec<_>=y.iter().zip(&b).map(|(x,d)|x+0.5*dt*d).collect(); let c=derivative(t+0.5*dt,&z);
        let z:Vec<_>=y.iter().zip(&c).map(|(x,d)|x+dt*d).collect(); let d=derivative(t+dt,&z);
        for i in 0..y.len() { y[i]+=dt*(a[i]+2.0*b[i]+2.0*c[i]+d[i])/6.0; }
    }
    (y[..n].chunks_exact(6).map(|q|q.try_into().unwrap()).collect(),y[n..].try_into().unwrap())
}

#[test]
fn retained_anisotropic_memory_matches_real_ode_work_heat_and_endpoint_handoff() {
    let law=law(); let w=[2.3,3.9]; let fields=fields(); let bias=[0.001,-0.002,0.0005,0.0,0.0,0.0];
    let t0=0.37; let h=0.29; let mut initial=law.state();
    law.step(&mut initial,waveform(&w,&fields,&bias,t0).0,0.43).unwrap();
    let before=initial.clone(); let source=law.prepare_harmonic_spectrum(&w).unwrap();
    let temporal=source.prepare_causal_window(HarmonicWindow::new(t0,h).unwrap()).unwrap();
    assert!(temporal.numerical_payload_bytes()<=source.causal_window_payload_bound_bytes().unwrap());
    assert!(source.causal_window_work_bound().unwrap()>0);
    let result=temporal.evaluate_with_bias(&fields,&bias,&initial).unwrap(); assert_eq!(initial,before);
    let (q,integrals)=oracle(&law,&w,&fields,&bias,&initial,t0,h);
    for (actual,expected) in result.state.q.iter().flatten().zip(q.iter().flatten()) { near(*actual,*expected,0.01,3e-12); }
    let point=&result.material;
    near(point.kelvin_heat_j_m3,integrals[0],integrals[0],3e-12);
    near(point.maxwell_heat_j_m3,integrals[1],integrals[1],3e-12);
    near(point.stress_work_j_m3,integrals[2],integrals[2].abs(),3e-12);
    let (end_strain,end_rate)=waveform(&w,&fields,&bias,t0+h);
    near(point.end_stored_energy_j_m3,stored(&law,&end_strain,&q),point.end_stored_energy_j_m3,3e-12);
    let mut stress=apply(law.c_eq(),&end_strain); let viscous=apply(law.kv_pa_s().unwrap(),&end_rate);
    for (branch,q) in law.branches().iter().zip(&q) {
        let s=apply(&branch.c_pa,q); for i in 0..6 { stress[i]+=s[i]; }
    }
    for i in 0..6 { near(result.state.stress[i],stress[i]+viscous[i],0.1,3e-12); }
    near(point.energy_balance_residual_j_m3,0.0,point.material_heat_j_m3+point.stress_work_j_m3.abs(),3e-13);
    let settled=source.prepare_window(HarmonicWindow::new(t0,h).unwrap()).unwrap().evaluate_with_bias(&fields,&bias).unwrap();
    assert!((point.maxwell_heat_j_m3-settled.maxwell_heat_j_m3).abs()>0.01*point.maxwell_heat_j_m3);
    let mut next=result.state.clone(); let old=next.clone();
    let target=core::array::from_fn(|i|next.strain[i]+0.0001*(i as f64-2.0));
    let step=law.step_with_mean(&mut next,target,0.07).unwrap();
    near(step.step_work-(next.stored-old.stored)-step.step_dissipated,0.0,old.stored+next.stored,2e-13);
    near(result.state.work-initial.work,point.stress_work_j_m3,point.stress_work_j_m3.abs(),2e-14);
}

#[test]
fn cosine_maximum_startup_keeps_true_quintic_heat_in_a_tiny_window() {
    let law=TensorProny::new_mandel(matrix(1e-40,[0.0;6]),vec![PronyBranch::new(matrix(1.0,[0.0;6]),0.8).unwrap()],0.1).unwrap();
    let amplitude=0.002; let omega=3.0; let h=1e-10; let mut strain=[0.0;6]; strain[0]=amplitude;
    let initial=state(&law,strain,vec![[0.0;6]]); let mut field=[C64::ZERO;6]; field[0]=C64::from_re(amplitude);
    let source=law.prepare_harmonic_spectrum(&[omega]).unwrap();
    let result=source.prepare_causal_window(HarmonicWindow::new(0.0,h).unwrap()).unwrap().evaluate(&[field],&initial).unwrap();
    let expected_q=-0.5*amplitude*omega*omega*h*h;
    let expected_heat=amplitude*amplitude*omega.powi(4)*h.powi(5)/(20.0*0.8);
    near(result.state.q[0][0],expected_q,expected_q.abs(),3e-10);
    near(result.material.maxwell_heat_j_m3,expected_heat,expected_heat,3e-10);
    assert!(result.material.maxwell_heat_j_m3>0.0);
    near(result.material.energy_balance_residual_j_m3,0.0,result.material.stress_work_j_m3.abs(),3e-14);
}

#[test]
fn retained_zero_drive_heat_survives_equal_rounded_endpoint_energies() {
    let law=TensorProny::new_mandel(matrix(1.0,[0.0;6]),vec![PronyBranch::new(matrix(2.0,[0.0;6]),0.8).unwrap()],0.1).unwrap()
        .with_band((1.0,10.0)).unwrap();
    let q=[0.01,-0.005,0.0,0.0,0.0,0.0]; let initial=state(&law,[0.0;6],vec![q]); let h=1e-20;
    let source=law.prepare_harmonic_spectrum(&[]).unwrap();
    let result=source.prepare_causal_window(HarmonicWindow::new(4.0,h).unwrap()).unwrap().evaluate(&[],&initial).unwrap();
    let expected=2.0*initial.stored*h/0.8;
    near(result.material.maxwell_heat_j_m3,expected,expected,4e-15);
    near(result.material.stored_energy_change_j_m3,-expected,expected,4e-15);
    assert_eq!(result.material.start_stored_energy_j_m3,result.material.end_stored_energy_j_m3);
    assert_eq!(result.material.stress_work_j_m3,0.0);
    assert_eq!(result.state.q,initial.q);
}

#[test]
fn already_settled_history_matches_the_original_coherent_observer() {
    let law=law(); let w=[2.3,3.9]; let fields=fields(); let bias=[0.001;6]; let t0=0.37;
    let q=law.branches().iter().map(|branch| {
        let mut q=[0.0;6];
        for (&omega,e) in w.iter().zip(&fields) {
            let x=omega*branch.tau_s; let den=1.0+x*x; let (s,c)=(omega*t0).sin_cos();
            for i in 0..6 { q[i]+=((x*x*e[i].re+x*e[i].im)*c+(x*x*e[i].im-x*e[i].re)*s)/den; }
        } q
    }).collect();
    let initial=state(&law,waveform(&w,&fields,&bias,t0).0,q);
    let source=law.prepare_harmonic_spectrum(&w).unwrap(); let window=HarmonicWindow::new(t0,0.43).unwrap();
    let actual=source.prepare_causal_window(window).unwrap().evaluate_with_bias(&fields,&bias,&initial).unwrap().material;
    let expected=source.prepare_window(window).unwrap().evaluate_with_bias(&fields,&bias).unwrap();
    for (a,b) in [
        (actual.start_stored_energy_j_m3,expected.start_stored_energy_j_m3),
        (actual.end_stored_energy_j_m3,expected.end_stored_energy_j_m3),
        (actual.stored_energy_change_j_m3,expected.stored_energy_change_j_m3),
        (actual.material_heat_j_m3,expected.material_heat_j_m3),
        (actual.stress_work_j_m3,expected.stress_work_j_m3),
    ] { near(a,b,b.abs(),4e-13); }
}

#[test]
fn accepted_intervals_keep_dd_clock_and_compose_original_memory_and_ledgers() {
    let law=law(); let w=[2.3,3.9]; let fields=fields(); let start=Dd::from_pair(1e16,0.3);
    let mut strain=[0.0;6];
    for (&omega,e) in w.iter().zip(&fields) {
        let phase=phase_angular_interval_from_clock_dd(omega,start,0.125).unwrap().start;
        for i in 0..6 { strain[i]+=(Dd::from_f64(e[i].re)*phase[0]-Dd::from_f64(e[i].im)*phase[1]).to_f64(); }
    }
    let initial=state(&law,strain,vec![[0.0007;6],[-0.0003;6]]); let source=law.prepare_harmonic_spectrum(&w).unwrap();
    let full=source.prepare_causal_window_from_clock_dd(start,0.375).unwrap().evaluate(&fields,&initial).unwrap();
    let first=source.prepare_causal_window_from_clock_dd(start,0.125).unwrap().evaluate(&fields,&initial).unwrap();
    let second_start=start+Dd::from_f64(0.125); assert_eq!(second_start.hi,start.hi); assert_ne!(second_start.lo,start.lo);
    let second=source.prepare_causal_window_from_clock_dd(second_start,0.25).unwrap().evaluate(&fields,&first.state).unwrap();
    for (a,b) in full.state.q.iter().flatten().zip(second.state.q.iter().flatten()) { near(*a,*b,0.01,3e-13); }
    near(full.state.work,second.state.work,full.state.work.abs(),4e-13);
    near(full.state.dissipated,second.state.dissipated,full.state.dissipated,4e-13);
    near(full.material.stored_energy_change_j_m3,first.material.stored_energy_change_j_m3+second.material.stored_energy_change_j_m3,full.state.stored+initial.stored,4e-13);
}

#[test]
fn cancellation_source_domain_and_continuity_refusals_preserve_actual_history() {
    let law=law().with_band((1.0,10.0)).unwrap(); let source=law.prepare_harmonic_spectrum(&[2.0,2.0]).unwrap();
    let initial=state(&law,[0.0;6],vec![[0.001;6],[-0.002;6]]); let before=initial.clone();
    let field=[C64::new(1e100,-1e100);6]; let fields=[field,field.map(|z|-z)];
    let window=source.prepare_causal_window(HarmonicWindow::new(0.0,0.2).unwrap()).unwrap();
    let actual=window.evaluate(&fields,&initial).unwrap(); assert!(actual.material.material_heat_j_m3>0.0);
    let empty=law.prepare_harmonic_spectrum(&[]).unwrap();
    let held=empty.prepare_causal_window(HarmonicWindow::new(0.0,0.2).unwrap()).unwrap().evaluate(&[],&initial).unwrap();
    near(actual.material.material_heat_j_m3,held.material.material_heat_j_m3,held.material.material_heat_j_m3,3e-14);
    let mut polls=0;
    assert!(matches!(window.evaluate_with_checkpoint(&fields,&initial,|| {
        polls+=1; if polls==7 {Err(ViscoError::Parameters{what:"cancelled"})}else{Ok(())}
    }),Err(ViscoError::Parameters{what:"cancelled"})));
    let mut wrong=initial.clone(); wrong.strain[0]=0.001; assert!(window.evaluate(&fields,&wrong).is_err());
    wrong=initial.clone(); wrong.stored*=2.0; assert!(window.evaluate(&fields,&wrong).is_err());
    wrong=initial.clone(); wrong.q.pop(); assert!(window.evaluate(&fields,&wrong).is_err());
    assert!(window.evaluate_with_bias(&fields,&[0.001;6],&initial).is_err());
    assert!(window.evaluate(&[[C64::from_re(0.1);6],[C64::ZERO;6]],&initial).is_err());
    assert!(source.prepare_causal_window(HarmonicWindow::with_profile(0.0,0.2,HarmonicWindowProfile::Triangular).unwrap()).is_err());
    assert_eq!(initial,before);
}

#[test]
fn original_stress_work_does_not_require_a_representable_stress_squared_moment() {
    let law=TensorProny::new_mandel(matrix(1.0,[0.0;6]),vec![
        PronyBranch::new(matrix(1e200,[0.0;6]),0.8).unwrap(),
    ],0.1).unwrap();
    let mut memory=[0.0;6]; memory[0]=0.001;
    let initial=state(&law,[0.0;6],vec![memory]);
    let mut field=[C64::ZERO;6]; field[0]=C64::new(0.0,0.002);
    let fields=[field]; let frequencies=[3.0]; let h=0.1;
    // Stress is physical and finite, although its unneeded square overflows.
    assert!(initial.stress[0].is_finite());
    assert!(!(initial.stress[0]*initial.stress[0]).is_finite());
    let source=law.prepare_harmonic_spectrum(&frequencies).unwrap();
    let result=source.prepare_causal_window(HarmonicWindow::new(0.0,h).unwrap()).unwrap()
        .evaluate(&fields,&initial).unwrap();
    let (_,integrals)=oracle(&law,&frequencies,&fields,&[0.0;6],&initial,0.0,h);
    near(result.material.maxwell_heat_j_m3,integrals[1],integrals[1],3e-12);
    near(result.material.stress_work_j_m3,integrals[2],integrals[2].abs(),3e-12);
    assert!(result.state.stress.iter().all(|value|value.is_finite()));
}
