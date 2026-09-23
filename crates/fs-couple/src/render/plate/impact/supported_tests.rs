use super::*;
use super::super::{ImpactConfig, ImpactSystem, ImpactSubstepConfig};
use fs_exec::CancelGate;

fn support()->TranslatingSupport {TranslatingSupport {mass_kg:0.01,stiffness_n_m:12.0,
    damping_n_s_m:0.02,initial_position_m:0.0,initial_velocity_m_s:0.0,
    maximum_travel_m:0.02,maximum_slope:0.2}}
fn wire(n:usize)->WireSpan {WireSpan {endpoints_m:[[-0.1,0.0],[0.1,0.0]],linear_density_kg_m:0.01,
    tension_n:0.8,bending_stiffness_n_m2:1e-7,damping_per_s:vec![2.0;n]}}
fn model()->SupportedStrings {SupportedStrings::new(&[wire(3),wire(2)],
    Some(StringStretching {axial_rigidity_n:100.0,maximum_slope:0.2}),support()).unwrap()}
fn config()->ImpactConfig {ImpactConfig {dt_s:2e-5,max_steps:1200,maximum_energy_j:1.0,
    energy_absolute_tolerance_j:1e-10,energy_relative_tolerance:1e-7,maximum_generalized_force:1e4}}

#[test]
fn rigid_translation_keeps_all_wire_mass_and_relative_rest() {
    for n in [1,2,8,16] {
        let s=SupportedStrings::new(&[wire(n),wire(n)],None,
            TranslatingSupport {initial_position_m:0.001,initial_velocity_m_s:0.07,..support()}).unwrap();
        let body=s.body();let state=s.initial_state();
        let kinetic=body.initial.iter().map(|v|0.5*v.velocity_m_sqrt_kg_per_s.powi(2)).sum::<f64>();
        assert!((kinetic-0.5*(0.01+2.0*0.002)*0.07_f64.powi(2)).abs()<1e-18);
        let z:Vec<_>=state.chunks_exact(2).map(|s|s[0]).collect();
        assert!((s.potential(&z)-0.5*12.0*0.001_f64.powi(2)).abs()<1e-18);
        let observation=s.observe_interleaved(&state,0).unwrap();
        assert!((observation.position_m-0.001).abs()<1e-18);
        assert!(observation.maximum_wire_slope<1e-15);
        assert!(s.force_weight().recip().powi(2)>support().mass_kg,
            "unresolved sine translation mass must not disappear");
    }
}

#[test]
fn kinetic_normalization_matches_independent_spatial_quadrature_and_force_work() {
    let s=model();let n=s.support_coordinate();let mut dz=vec![0.0;s.mode_count()];
    for (i,v) in dz.iter_mut().enumerate() {*v=0.017*(i as f64+0.4);}
    let v=dz[n]*s.force_weight();let mut kinetic=0.5*support().mass_kg*v*v;
    for (wire,r) in s.0.wires.iter().zip(&s.0.ranges) {
        let cells=4096;let dx=wire.length_m()/cells as f64;
        for cell in 0..cells {
            let x=(cell as f64+0.5)*dx;let mut speed=v;
            for (mode,i) in (r.clone()).enumerate() {
                let phi=(2.0/(wire.linear_density_kg_m*wire.length_m())).sqrt()
                    *((mode+1) as f64*core::f64::consts::PI*x/wire.length_m()).sin();
                speed+=phi*(dz[i]-s.0.translation[i]*v);
            }
            kinetic+=0.5*wire.linear_density_kg_m*speed*speed*dx;
        }
    }
    let expected=0.5*dz.iter().map(|v|v*v).sum::<f64>();
    assert!((kinetic-expected).abs()<2e-8*expected);
    let f=0.37;assert!(((f*s.force_weight())*dz[n]-f*v).abs()<4.0*f64::EPSILON*(f*v).abs());
}

#[test]
fn supported_energy_tangent_and_relative_damping_retain_the_transpose_reaction() {
    let s=model();let n=s.mode_count();let z:Vec<_>=(0..n).map(|i|1e-5*(i as f64-1.2)).collect();
    let d:Vec<_>=(0..n).map(|i|0.003*(i as f64+0.7)).collect();
    let (mut g,mut hv)=(vec![0.0;n],vec![0.0;n]);s.gradient(&z,&mut g);s.hessian_vector(&z,&d,&mut hv);
    for i in 0..n {
        let h=1e-8;let mut a=z.clone();let mut b=z.clone();a[i]+=h;b[i]-=h;
        let fd=(s.potential(&a)-s.potential(&b))/(2.0*h);
        assert!((g[i]-fd).abs()<1e-6*g[i].abs().max(1e-5));
    }
    let h=1e-6;let a:Vec<_>=z.iter().zip(&d).map(|(z,d)|z+h*d).collect();
    let b:Vec<_>=z.iter().zip(&d).map(|(z,d)|z-h*d).collect();
    let (mut ga,mut gb)=(vec![0.0;n],vec![0.0;n]);s.gradient(&a,&mut ga);s.gradient(&b,&mut gb);
    for i in 0..n {assert!((hv[i]-(ga[i]-gb[i])/(2.0*h)).abs()<1e-7*hv[i].abs().max(1.0));}
    let dim=2*n;let mut r=vec![0.0;dim*dim];s.add_resistance(0,dim,&mut r).unwrap();
    let mut power=0.0;
    for i in 0..n {for j in 0..n {power+=d[i]*r[(2*i+1)*dim+2*j+1]*d[j];}}
    let velocity=d[n-1]*s.force_weight();
    let expected=0.02*velocity*velocity+(0..n-1).map(|i|s.0.damping[i]*(d[i]-s.0.translation[i]*velocity).powi(2)).sum::<f64>();
    assert!((power-expected).abs()<1e-14 && power>0.0);
}

#[test]
fn absolute_wire_contact_preserves_gaps_units_and_relative_kinematics() {
    let s=model();let line=LineContact::uniform(0.2,8,0.00002,1e6,1.5,0.05,"fixture".into()).unwrap();
    let shapes=vec![vec![2.0];8];let n=1+s.mode_count();
    let ob=s.contact(0,&line,&shapes,0..1,1,n).unwrap();
    let local=s.0.wires[0].contact(&line,&shapes,0..1,1..4,n).unwrap();
    let x=0.001;let mut state=vec![0.0;n];state[0]=0.0001;
    for i in 0..s.support_coordinate() {state[1+i]=s.0.translation[i]*x;}
    state[n-1]=x/s.force_weight();
    for (row,original) in ob.collocation().chunks_exact(n).zip(local.collocation().chunks_exact(n)) {
        assert_eq!(&row[..n-1],&original[..n-1]);
        let relative: f64=row.iter().zip(&state).map(|(b,q)|b*q).sum();
        assert!((relative-(2.0*state[0]-x)).abs()<1e-17);
    }
    assert_eq!(ob.gaps(),line.gaps_m);assert_eq!(ob.weights(),line.measures_m);assert_eq!(ob.internal_loss(),0.05);
    assert!(s.contact(0,&line,&shapes,1..2,1,n).is_err());
}

fn experiment()->(ImpactSystem,usize,f64) {
    let s=SupportedStrings::new(&[wire(2)],None,TranslatingSupport {
        initial_position_m:0.0001,initial_velocity_m_s:-0.12,stiffness_n_m:0.0,damping_n_s_m:0.0,..support()}).unwrap();
    let target=ImpactBody {potential:BodyPotential::Linear(vec![1000.0]),
        initial:vec![ModalAcousticState::default()],damping_per_s:vec![0.0]};
    let line=LineContact::uniform(0.2,5,0.00002,1e6,1.5,0.2,"moving support impact".into()).unwrap();
    let n=1+s.mode_count();let ob=s.contact(0,&line,&vec![vec![2.0];5],0..1,1,n).unwrap();
    let system=ImpactSystem::new(vec![target,s.body()],vec![ob],vec![],vec![],config()).unwrap();
    (system,n-1,s.force_weight())
}

#[test]
fn moving_support_releases_real_contact_without_resets_and_preserves_retry_work() {
    let (system,port,weight)=experiment();let initial=system.stored_energy_j();
    let mut moved=system.prepare_analytic().unwrap().with_substeps(ImpactSubstepConfig {max_depth:4,max_attempts:31}).unwrap();
    let (system,_,_)=experiment();let mut clean=system.prepare_analytic().unwrap().with_substeps(ImpactSubstepConfig {max_depth:4,max_attempts:31}).unwrap();
    let gate=CancelGate::new_clock_free();let mut force=vec![0.0;port+1];let mut net=0.0;
    let (mut head,mut loss,mut minimum)=(0.0_f64,0.0,0.0001_f64);
    for tick in 0..600 {
        // A real force reverses the rail and withdraws the wires after contact.
        force[port]=if tick>=150 {0.3*weight}else{0.0};
        if tick==175 {
            let before=moved.state().to_vec();let cancel=CancelGate::new_clock_free();cancel.request();
            assert!(moved.step(&force,&cancel).is_err());assert_eq!(moved.state(),before);
        }
        let f=moved.step(&force,&gate).unwrap();clean.step(&force,&gate).unwrap();assert_eq!(moved.state(),clean.state());
        net+=f.supplied_work_j-f.dissipated_energy_j;loss+=f.dissipated_energy_j;
        assert!((f.stored_energy_j-initial-net).abs()<1e-7);
        head=head.max(moved.state()[0].abs());minimum=minimum.min(moved.support_observation(1).unwrap().position_m);
    }
    assert!(head>1e-8 && loss>0.0 && minimum<0.0);
    assert!(moved.support_observation(1).unwrap().position_m>0.0002);
    assert_eq!(moved.samples(),600);
}

#[test]
fn support_travel_slope_and_complete_state_budgets_are_enforced() {
    assert!(SupportedStrings::new(&[wire(256)],None,support()).is_err());
    assert!(SupportedStrings::new(&[wire(2)],None,TranslatingSupport {mass_kg:0.0,..support()}).is_err());
    let s=model();let mut state=s.initial_state();let n=s.support_coordinate();
    state[2*n]=0.021/s.force_weight();assert!(s.observe_interleaved(&state,0).is_err());
    state=s.initial_state();state[0]=0.1;assert!(s.observe_interleaved(&state,0).is_err());
    let s=SupportedStrings::new(&[wire(1)],None,TranslatingSupport {maximum_travel_m:1e-7,..support()}).unwrap();
    let mut system=ImpactSystem::new(vec![s.body()],vec![],vec![],vec![],config()).unwrap().prepare_analytic().unwrap();
    let before=system.state().to_vec();let gate=CancelGate::new_clock_free();
    assert!(system.step(&[0.0,100.0],&gate).is_err());assert_eq!(system.state(),before);assert_eq!(system.samples(),0);
    system.step(&[0.0;2],&gate).unwrap();assert_eq!(system.samples(),1);
}
