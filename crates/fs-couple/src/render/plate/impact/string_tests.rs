use super::*;
use super::super::{ImpactConfig, ImpactSystem, ImpactSubstepConfig};
use fs_dcontact::Obstacle;
use fs_exec::CancelGate;

fn wire(n: usize) -> WireSpan {
    WireSpan { endpoints_m: [[0.0,0.0],[0.3,0.0]], linear_density_kg_m: 0.01,
        tension_n: 1.0, bending_stiffness_n_m2: 1e-6, damping_per_s: vec![0.0;n] }
}
fn law(ea: f64) -> StringStretching { StringStretching { axial_rigidity_n: ea, maximum_slope: 0.2 } }
fn model(span: &WireSpan, ea: f64) -> StringPotential {
    let b=span.stretching_body(vec![ModalAcousticState::default();span.damping_per_s.len()],law(ea)).unwrap();
    let BodyPotential::String(s)=b.potential else {panic!()};s
}
fn config(dt_s: f64) -> ImpactConfig {
    ImpactConfig { dt_s, max_steps: 10_000, maximum_energy_j: 2.0,
        energy_absolute_tolerance_j: 1e-10, energy_relative_tolerance: 1e-8,
        maximum_generalized_force: 1e4 }
}

#[test]
fn existing_stress_channel_retains_bending_and_matches_spatial_extension_energy() {
    let span=wire(3);let ea=2000.0;let s=model(&span,ea);
    let linear=span.body(vec![ModalAcousticState::default();3]).unwrap();
    let BodyPotential::Linear(w)=linear.potential else {panic!()};
    assert_eq!(s.omegas(),w);
    let mass=0.5*span.linear_density_kg_m*span.length_m();
    let physical=[0.001,-0.0003,0.0001];let q=physical.map(|a|a*det::sqrt(mass));
    let mut slope2=0.0;let cells=1024;
    for i in 0..cells {
        let x=span.length_m()*(i as f64+0.5)/cells as f64;
        let slope: f64=physical.iter().enumerate().map(|(k,a)| {
            let k=(k+1) as f64*std::f64::consts::PI/span.length_m();a*k*(k*x).cos()
        }).sum();
        slope2+=slope*slope*span.length_m()/cells as f64;
    }
    let expected_tension=span.tension_n+ea*slope2/(2.0*span.length_m());
    let expected_energy=ea*slope2*slope2/(8.0*span.length_m());
    let observation=s.observe(&q).unwrap();
    assert!((observation.tension_n-expected_tension).abs()<1e-12);
    assert!((observation.stretching_energy_j-expected_energy).abs()<1e-16);
    let linear_energy: f64=w.iter().zip(q).map(|(w,q)|0.5*(w*q).powi(2)).sum();
    assert!((s.potential(&q)-linear_energy-expected_energy).abs()<1e-16);
    assert!(observation.slope_bound>0.0 && observation.additional_strain>0.0);
    let zero=model(&span,0.0);assert_eq!(zero.observe(&q).unwrap().tension_n,span.tension_n);
    assert!((zero.potential(&q)-linear_energy).abs()<1e-16);
}

#[test]
fn analytic_string_force_and_tangent_use_the_same_nonlinear_storage() {
    let s=model(&wire(3),5000.0);let q=[0.00008,-0.00002,0.00001];let d=[0.2,-0.1,0.3];
    let (mut g,mut hd)=([0.0;3],[0.0;3]);s.gradient(&q,&mut g);s.hessian_vector(&q,&d,&mut hd);
    for i in 0..3 {
        let h=1e-9;let mut hi=q;let mut lo=q;hi[i]+=h;lo[i]-=h;
        let difference=(s.potential(&hi)-s.potential(&lo))/(2.0*h);
        assert!((g[i]-difference).abs()<2e-7*g[i].abs().max(1.0));
    }
    let h=1e-9;let mut hi=q;let mut lo=q;
    for i in 0..3 {hi[i]+=h*d[i];lo[i]-=h*d[i];}
    let (mut gp,mut gm)=([0.0;3],[0.0;3]);s.gradient(&hi,&mut gp);s.gradient(&lo,&mut gm);
    for i in 0..3 {assert!((hd[i]-(gp[i]-gm[i])/(2.0*h)).abs()<2e-7*hd[i].abs().max(1.0));}
    // An excited partial changes another partial's incremental restoring force.
    let (mut at_rest,mut excited)=([0.0;3],[0.0;3]);
    s.hessian_vector(&[0.0;3],&[0.0,1.0,0.0],&mut at_rest);
    s.hessian_vector(&q,&[0.0,1.0,0.0],&mut excited);
    assert!(excited[1]>at_rest[1]);assert!(excited[0].abs()>0.0);
}

fn period(dt: f64, amplitude_m: f64, ea: f64) -> f64 {
    let mut span=wire(1);span.bending_stiffness_n_m2=0.0;
    let a=amplitude_m*det::sqrt(0.5*span.linear_density_kg_m*span.length_m());
    let body=span.stretching_body(vec![ModalAcousticState {
        displacement_m_sqrt_kg:a,velocity_m_sqrt_kg_per_s:0.0 }],law(ea)).unwrap();
    let mut system=ImpactSystem::new(vec![body],vec![],vec![],vec![],config(dt)).unwrap().prepare_analytic().unwrap();
    let initial=system.stored_energy_j();let gate=CancelGate::new_clock_free();let mut previous=a;
    for i in 1..=1000 {
        let f=system.step(&[0.0],&gate).unwrap();let q=system.state()[0];
        assert!((f.stored_energy_j-initial).abs()<1e-8);
        if q<=0.0 {return 4.0*((i-1) as f64+previous/(previous-q))*dt;}
        previous=q;
    }
    panic!("no first quarter-cycle");
}
#[test]
fn finite_amplitude_period_matches_continuous_duffing_quadrature_under_refinement() {
    let amplitude=0.003;let ea=2000.0;let s=model(&wire(1),ea);
    let a=amplitude*det::sqrt(0.5*0.01*0.3);
    let beta=fs_nlmodal::single_mode_beta(&s.storage,0);
    let omega=std::f64::consts::PI/0.3*det::sqrt(1.0/0.01); // zero-EI fixture above
    let f=|theta:f64| (1.0+beta*a*a/(2.0*omega*omega)*(1.0+theta.sin().powi(2))).sqrt().recip();
    let n=1024;let h=std::f64::consts::FRAC_PI_2/n as f64;let mut integral=f(0.0)+f(std::f64::consts::FRAC_PI_2);
    for i in 1..n {integral+=if i%2==0 {2.0*f(i as f64*h)}else{4.0*f(i as f64*h)};}
    let exact=4.0/omega*h*integral/3.0;
    let coarse=period(0.0001,amplitude,ea);let fine=period(0.00005,amplitude,ea);
    assert!((fine/exact-1.0).abs()<1e-4);
    assert!((fine-exact).abs()<0.6*(coarse-exact).abs());
    assert!(fine<0.95*period(0.00005,amplitude,0.0));
}

#[test]
fn nonlinear_wire_and_receiver_share_contact_work_and_keep_exact_retry() {
    let make=|| {
        let span=wire(2);let (receiver,w)=ImpactBody::free_mass(0.02,-0.00001,0.4).unwrap();
        let body=span.stretching_body(vec![ModalAcousticState::default();2],law(5000.0)).unwrap();
        let shapes=fs_dcontact::string_collocation(span.length_m(),span.linear_density_kg_m,&[0.11],2).unwrap();
        let contact=Obstacle::new(vec![w,-shapes[0],-shapes[1]],1,3,vec![0.0],vec![1.0],2e7,1.5,
            "authored impact on a stretching filament".into()).unwrap().with_internal_loss(0.05).unwrap();
        ImpactSystem::new(vec![receiver,body],vec![contact],vec![],vec![],config(2e-6)).unwrap()
    };
    let mut fd=make().prepare().unwrap();let mut exact=make().prepare_analytic().unwrap();
    let mut retry=make().prepare_analytic().unwrap();let initial=exact.stored_energy_j();let mut net=0.0;
    let mut peak_tension=1.0_f64;let gate=CancelGate::new_clock_free();
    for i in 0..512 {
        let force=if i<64 {[0.01,0.0,0.0]}else{[0.0;3]};
        if i==128 {
            let before=retry.state().to_vec();let samples=retry.samples();
            retry.set_iteration_limit(0).unwrap();assert!(retry.step(&force,&gate).is_err());
            assert_eq!(retry.state(),before);assert_eq!(retry.samples(),samples);retry.set_iteration_limit(50).unwrap();
        }
        fd.step(&force,&gate).unwrap();let f=exact.step(&force,&gate).unwrap();retry.step(&force,&gate).unwrap();
        assert_eq!(exact.state(),retry.state());
        for (a,b) in fd.state().iter().zip(exact.state()) {assert!((a-b).abs()<1e-7);}
        net+=f.supplied_work_j-f.dissipated_energy_j;
        assert!((f.stored_energy_j-initial-net).abs()<1e-7);
        peak_tension=peak_tension.max(exact.string_observation(1).unwrap().tension_n);
        assert!(exact.string_observation(0).is_none());
    }
    assert!(peak_tension>1.000001);assert!(exact.state()[2..].iter().any(|x|x.abs()>1e-8));
}

#[test]
fn continuous_slope_limit_refuses_initial_and_endpoint_motion_without_state_change() {
    let span=wire(1);let tiny=StringStretching {maximum_slope:0.0001,..law(100.0)};
    let state=ModalAcousticState {displacement_m_sqrt_kg:0.001,velocity_m_sqrt_kg_per_s:0.0};
    assert!(span.stretching_body(vec![state],tiny).is_err());
    for substeps in [false,true] {
        let body=span.stretching_body(vec![ModalAcousticState {
            displacement_m_sqrt_kg:0.0,velocity_m_sqrt_kg_per_s:0.02 }],tiny).unwrap();
        let mut prepared=ImpactSystem::new(vec![body],vec![],vec![],vec![],config(0.001)).unwrap().prepare_analytic().unwrap();
        let gate=CancelGate::new_clock_free();let before=prepared.state().to_vec();
        if substeps {
            let mut s=prepared.with_substeps(ImpactSubstepConfig {max_depth:4,max_attempts:31}).unwrap();
            assert!(s.step(&[0.0],&gate).is_err());assert_eq!(s.state(),before);assert_eq!(s.samples(),0);
        } else {
            assert!(prepared.step(&[0.0],&gate).is_err());assert_eq!(prepared.state(),before);assert_eq!(prepared.samples(),0);
        }
    }
    for ea in [-1.0,f64::INFINITY,f64::NAN] {assert!(law(ea).validate().is_err());}
    assert!(StringStretching {maximum_slope:0.31,..law(10.0)}.validate().is_err());
    assert!(wire(MAX_IMPACT_MODES+1).stretching_body(vec![ModalAcousticState::default();MAX_IMPACT_MODES+1],law(1.0)).is_err());
}
