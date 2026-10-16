use super::*;
use super::super::KelvinBranch;
use super::super::super::{ImpactConfig, ImpactSystem};
use fs_exec::CancelGate;
use fs_material::fiber::WoolFelt;
use fs_phs::Storage;

fn carrier()->Carrier {Carrier {mass_kg:0.4,stiffness_n_m:0.0,resistance_n_s_m:0.0,
    initial_position_m:0.0,initial_velocity_m_s:0.0}}
fn sites(weight:f64,creep:bool)->Vec<MountSite> {
    [-1.0,1.0].into_iter().map(|sign|MountSite {
        compression_per_carrier_m:-sign,
        pad:FeltPad {area_m2:0.0003,thickness_m:0.006,precompression_m:0.0003,
            weights:vec![sign*weight],law:WoolFelt::new(30000.0,0.2,2.2,3.0,0.15,0.7).unwrap(),
            prior_maximum_strain:0.05,creep:if creep {vec![KelvinBranch {
                stiffness_n_m:4500.0,viscosity_n_s_m:24.0}]}else{vec![]}},
    }).collect()
}
fn config(steps:u64)->ImpactConfig {ImpactConfig {dt_s:1e-5,max_steps:steps,maximum_energy_j:20.0,
    energy_absolute_tolerance_j:1e-9,energy_relative_tolerance:1e-6,maximum_generalized_force:1e5}}
fn system(c:Carrier,creep:bool,steps:u64)->(ImpactSystem,CarrierPort) {
    let (source,w)=ImpactBody::free_mass(0.1,0.0,0.03).unwrap();
    let mount=MountedFelt::new(1,c,&sites(w,creep)).unwrap();
    let s=ImpactSystem::new(vec![source,mount.body],vec![],mount.pads,vec![],config(steps)).unwrap();
    (s,mount.port)
}

#[test]
fn one_carrier_retains_every_local_row_and_the_explicit_ground_energy() {
    let original=sites(2.0,true);let mut c=carrier();
    c.stiffness_n_m=3200.0;c.resistance_n_s_m=8.0;c.initial_position_m=0.002;c.initial_velocity_m_s=-0.04;
    let m=MountedFelt::new(1,c,&original).unwrap();
    assert_eq!(m.port.coordinate,1);assert_eq!(m.pads.len(),2);
    assert_eq!(m.body.damping_per_s,vec![20.0]);
    for (pad,site) in m.pads.iter().zip(&original) {
        assert_eq!(&pad.weights[..1],&site.pad.weights);
        assert_eq!(pad.weights[1],site.compression_per_carrier_m/mass_root());
        assert_eq!(pad.precompression_m,site.pad.precompression_m);
        assert_eq!(pad.prior_maximum_strain,site.pad.prior_maximum_strain);
        assert_eq!(pad.creep[0].stiffness_n_m,4500.0);
        assert_eq!(pad.creep[0].viscosity_n_s_m,24.0);
    }
    let q=m.body.initial[0].displacement_m_sqrt_kg;
    let p=m.body.initial[0].velocity_m_sqrt_kg_per_s;
    let energy=m.body.potential.potential(&[q])+0.5*p*p;
    assert!((energy-0.5*c.stiffness_n_m*c.initial_position_m.powi(2)
        -0.5*c.mass_kg*c.initial_velocity_m_s.powi(2)).abs()<1e-16);
    let motion=m.port.motion(&[0.0,0.0,q,p,9.0]).unwrap();
    assert!((motion[0]-c.initial_position_m).abs()<1e-17);
    assert!((motion[1]-c.initial_velocity_m_s).abs()<1e-16);
}
fn mass_root()->f64 {0.4_f64.sqrt()}

#[test]
fn both_backings_translate_together_and_contact_reactions_balance() {
    let (s,port)=system(carrier(),false,1);let mut x=s.state().to_vec();
    x[0]=0.0001/0.1_f64.sqrt().recip();x[1]=0.0;
    let mut g=vec![0.0;x.len()];s.mechanical.gradient(&x,&mut g);
    assert!(g[0]>0.0 && g[2]<0.0);
    assert!((g[0]*0.1_f64.sqrt()+g[2]/port.inverse_sqrt_mass).abs()<1e-12);
    let before=s.mechanical.hamiltonian(&x);
    x[0]+=0.001*0.1_f64.sqrt();x[2]+=0.001/port.inverse_sqrt_mass;
    assert!((s.mechanical.hamiltonian(&x)-before).abs()<1e-14);
    // A displaced source loads opposite washers differently. They are NOT
    // replaced by a mean-shape spring or separate carriers with doubled mass.
    assert_ne!(s.pads[0].strain(&x,2),s.pads[1].strain(&x,2));
}

#[test]
fn the_free_carrier_receives_real_motion_and_preserves_total_momentum() {
    let (mut s,port)=system(carrier(),false,200);let initial=s.stored_energy_j();
    let momentum=0.1*0.03;let mut loss=0.0;
    let gate=CancelGate::new_clock_free();
    for _ in 0..200 {
        let frame=s.step(&[0.0,0.0],&gate).unwrap();loss+=frame.dissipated_energy_j;
        assert!(frame.balance_residual_j.abs()<1e-8);
        let actual=0.1_f64.sqrt()*s.state()[1]+mass_root()*s.state()[3];
        assert!((actual-momentum).abs()<1e-11);
        assert!((frame.stored_energy_j+loss-initial).abs()<1e-7);
    }
    assert!(port.motion(s.state()).unwrap()[1]>1e-6);
}

#[test]
fn carrier_felt_and_creep_share_prepared_rollback_and_accepted_time() {
    let make=||system(carrier(),true,80).0.prepare_analytic().unwrap();
    let mut trial=make();let mut clean=make();let gate=CancelGate::new_clock_free();
    for tick in 0..80 {
        if tick==30 {
            let before=trial.state().to_vec();let count=trial.samples();let energy=trial.stored_energy_j();
            let history=trial.felt_history(0).unwrap();
            let cancel=CancelGate::new_clock_free();cancel.request();
            assert!(trial.step(&[0.0,0.0],&cancel).is_err());
            assert!(trial.step(&[0.0,1e6],&gate).is_err());
            assert_eq!(trial.state(),before);assert_eq!(trial.samples(),count);
            assert_eq!(trial.stored_energy_j(),energy);
            assert_eq!(trial.felt_history(0).unwrap(),history);
        }
        let a=trial.step(&[0.0,0.2],&gate).unwrap();
        let b=clean.step(&[0.0,0.2],&gate).unwrap();
        assert_eq!(trial.state(),clean.state());assert_eq!(a.felt_crush_loss_j,b.felt_crush_loss_j);
        assert_eq!(a.stored_energy_j,b.stored_energy_j);
    }
    assert!(trial.state()[4..].iter().any(|v|v.abs()>1e-12));
    let before=trial.state().to_vec();assert!(trial.step(&[0.0,0.0],&gate).is_err());
    assert_eq!(trial.state(),before);assert_eq!(trial.samples(),80);
}

#[test]
fn carrier_mass_changes_the_coupled_response_not_just_a_readout() {
    let mut c=carrier();c.mass_kg=0.04;
    let (mut light,_)=system(c,false,300);let (mut heavy,_)=system(carrier(),false,300);
    let gate=CancelGate::new_clock_free();
    for _ in 0..300 {light.step(&[0.0,0.0],&gate).unwrap();heavy.step(&[0.0,0.0],&gate).unwrap();}
    assert!((light.state()[0]-heavy.state()[0]).abs()>1e-8);
    assert!((light.state()[1]-heavy.state()[1]).abs()>1e-7);
}

#[test]
fn invalid_mounts_refuse_without_mutating_the_original_pads() {
    let original=sites(1.0,true);let before=original[0].pad.weights.clone();
    for c in [Carrier {mass_kg:0.0,..carrier()},Carrier {mass_kg:f64::NAN,..carrier()},
        Carrier {stiffness_n_m:-1.0,..carrier()},Carrier {resistance_n_s_m:-1.0,..carrier()},
        Carrier {initial_velocity_m_s:f64::INFINITY,..carrier()},
        Carrier {mass_kg:f64::MAX,stiffness_n_m:f64::MIN_POSITIVE,..carrier()}] {
        assert!(MountedFelt::new(1,c,&original).is_err());
    }
    for n in [0,MAX_IMPACT_MODES,usize::MAX] {assert!(MountedFelt::new(n,carrier(),&original).is_err());}
    for projection in [1.01,f64::NAN] {
        let mut s=original.clone();s[0].compression_per_carrier_m=projection;
        assert!(MountedFelt::new(1,carrier(),&s).is_err());
    }
    let mut uncoupled=original.clone();for s in &mut uncoupled {s.compression_per_carrier_m=0.0;}
    assert!(MountedFelt::new(1,carrier(),&uncoupled).is_err());
    assert!(MountedFelt::new(1,carrier(),&[]).is_err());
    assert!(MountedFelt::new(1,carrier(),&vec![original[0].clone();17]).is_err());
    assert_eq!(original[0].pad.weights,before);
    let port=CarrierPort {coordinate:usize::MAX,inverse_sqrt_mass:1.0};
    assert!(port.motion(&[]).is_err());
}
