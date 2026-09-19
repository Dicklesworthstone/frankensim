//! Same physical input through the prepared linear image and independent checks.
use fs_couple::modal_acoustic_time::{ModalAcousticState, ModalAcousticTimeBudget};
use fs_couple::render::plate::impact::{BodyPotential, ImpactBody, ImpactConfig, ImpactSystem, VolumeSpring};
use fs_couple::render::plate::impact::linear::{LinearImpactConfig, LinearImpactSystem, VolumeConnection};
use fs_couple::render::schedule::force::coupled::{ModalCouplingConfig,
    contact::ModalContactConfig, contact::multiple::MultiContactConfig};
use fs_dcontact::Obstacle;
use fs_exec::CancelGate;

fn config(rate: u32, steps: u64) -> LinearImpactConfig {
    LinearImpactConfig {
        sample_rate_hz: rate, max_steps: steps, maximum_generalized_force: 1e6,
        component: ModalAcousticTimeBudget { maximum_total_energy_j: 10.0,
            ..ModalAcousticTimeBudget::audible_reference() },
        coupling: ModalCouplingConfig { max_modes: 64, max_connections: 8, max_setup_terms: 100000,
            nyquist_guard_fraction: 0.9, maximum_total_energy_j: 10.0, maximum_abs_pressure_pa: 1e6,
            maximum_abs_connection_force_n: 1e5, solve_relative_tolerance: 1e-10,
            energy_absolute_tolerance_j: 1e-10, energy_relative_tolerance: 1e-8 },
        contact: ModalContactConfig { max_iterations: 100, maximum_force_n: 1e4,
            maximum_penetration_m: 0.01, force_absolute_tolerance_n: 1e-10,
            force_relative_tolerance: 1e-9 },
        multiple: MultiContactConfig { max_contacts: 32, max_sweeps: 100, max_setup_terms: 100000 },
    }
}
fn body(omega: f64, q: f64, v: f64) -> ImpactBody {
    ImpactBody { potential: BodyPotential::Linear(vec![omega]),
        initial: vec![ModalAcousticState { displacement_m_sqrt_kg: q, velocity_m_sqrt_kg_per_s: v }],
        damping_per_s: vec![if omega==0.0 {0.0} else {2.0}] }
}
fn physical_input() -> (Vec<ImpactBody>, Obstacle, VolumeSpring) {
    let (striker, b) = ImpactBody::free_mass(0.04, -0.0002, 0.8).unwrap();
    let bodies = vec![striker, body(1200.0,0.0,0.0), body(1700.0,0.0,0.0)];
    let contact = Obstacle::new(vec![b,-10.0,0.0],1,3,vec![0.0],vec![1.0],
        1e7,1.5,"synthetic tip test, not identified wood".into()).unwrap();
    let volume = VolumeSpring {bulk_modulus_pa: 1.4e5, volume_m3: 0.02, areas: vec![0.0,0.01,-0.01]};
    (bodies,contact,volume)
}
fn model(rate:u32, steps:u64) -> LinearImpactSystem {
    let (bodies,contact,volume) = physical_input();
    LinearImpactSystem::new(bodies,vec![contact],vec![VolumeConnection {spring:volume,reference_area_m2:0.02}],
        config(rate,steps),&CancelGate::new_clock_free()).unwrap()
}
fn reference(rate:u32, steps:u64) -> ImpactSystem {
    let (bodies,contact,volume) = physical_input();
    ImpactSystem::new(bodies,vec![contact],vec![],vec![volume],ImpactConfig {
        dt_s:1.0/f64::from(rate),max_steps:steps,maximum_energy_j:10.0,
        energy_absolute_tolerance_j:1e-10,energy_relative_tolerance:1e-8,maximum_generalized_force:1e6,
    }).unwrap()
}

#[test]
fn volume_coordinate_scale_preserves_physical_energy_and_two_way_motion() {
    let mut a = Vec::new();
    for reference_area_m2 in [0.005,0.02,0.08] {
        let (_,_,volume)=physical_input();
        let bodies=vec![body(0.0,0.0,0.0),body(1200.0,2e-4,0.0),body(1700.0,0.0,0.0)];
        let m=LinearImpactSystem::new(bodies,vec![],vec![VolumeConnection {spring:volume,reference_area_m2}],
            config(192000,512),&CancelGate::new_clock_free()).unwrap();
        let expected=0.5*(1200.0_f64*2e-4).powi(2)+0.5*(1.4e5/0.02)*(0.01_f64*2e-4).powi(2);
        assert!((m.frame().stored_energy_j-expected).abs()<1e-14);
        a.push(m);
    }
    let gate=CancelGate::new_clock_free();
    let mut receiver_peak=0.0_f64;
    for _ in 0..512 {
        for m in &mut a {m.step(&[0.0;3],&gate).unwrap();}
        receiver_peak=receiver_peak.max(a[0].state()[5].abs());
        for other in &a[1..] {for (x,y) in a[0].state().iter().zip(other.state()) {assert!((x-y).abs()<1e-10);}}
    }
    assert!(receiver_peak>1e-5,"unstruck second head must receive cavity work");
}

#[test]
fn prepared_and_gonzalez_images_converge_for_the_same_struck_two_head_system() {
    let gate=CancelGate::new_clock_free();
    let mut errors=Vec::new();
    for (rate,steps) in [(96000,384),(192000,768),(384000,1536)] {
        let mut modal=model(rate,steps); let mut reference=reference(rate,steps);
        assert!((modal.frame().stored_energy_j-0.5*0.04*0.8*0.8).abs()<1e-14);
        for _ in 0..steps {
            let frame=modal.step(&[0.0;3],&gate).unwrap();
            reference.step(&[0.0;3],&gate).unwrap();
            assert!(frame.balance_residual_j.abs()<1e-9);
        }
        let error=modal.state().iter().zip(reference.state()).enumerate().map(|(i,(a,b))|
            ((a-b)*if i%2==0 {1700.0}else{1.0}).powi(2)).sum::<f64>().sqrt();
        errors.push(error);
        assert!(modal.state()[1]<0.0,"physical striker must rebound");
        assert!(modal.state()[5].abs()>1e-6,"unstruck head must not stay silent");
    }
    // This is agreement under refinement, NOT bit identity of two integrators.
    assert!(errors[2]<0.4*errors[0],"image discrepancy did not refine: {errors:?}");
}

#[test]
fn distributed_contacts_are_solved_together_without_direct_receiver_drive() {
    let gate=CancelGate::new_clock_free();
    let (striker,b)=ImpactBody::free_mass(0.04,-0.0002,0.8).unwrap();
    let make=|stiffness| {
        let contact=Obstacle::new(vec![b,-10.0,0.0,b,0.0,-10.0],2,3,vec![0.0;2],vec![0.5;2],
            stiffness,1.5,"synthetic distributed impact".into()).unwrap();
        LinearImpactSystem::new(vec![striker.clone(),body(1200.,0.,0.),body(1200.,0.,0.)],vec![contact],
            vec![],config(192000,512),&gate).unwrap()
    };
    let mut active=make(1e7); let mut off=make(0.0); let mut peak=0.0_f64;
    assert_eq!(active.contact_count(),2);
    for _ in 0..512 {
        active.step(&[0.;3],&gate).unwrap(); off.step(&[0.;3],&gate).unwrap();
        peak=peak.max(active.state()[3].abs());
        assert!((active.state()[3]-active.state()[5]).abs()<1e-8);
        assert!(off.state()[2..].iter().all(|x|*x==0.0));
    }
    assert!(peak>1e-4);
}

#[test]
fn input_refusal_cancellation_and_budget_extension_keep_the_complete_state() {
    let gate=CancelGate::new_clock_free(); let mut actual=model(192000,5); let mut baseline=model(192000,9);
    for _ in 0..3 {actual.step(&[0.;3],&gate).unwrap();baseline.step(&[0.;3],&gate).unwrap();}
    let state=actual.state().to_vec(); let frame=*actual.frame();
    assert!(actual.step(&[f64::NAN,0.,0.],&gate).is_err());
    assert!(actual.step(&[1e7,0.,0.],&gate).is_err());
    assert!(actual.step(&[0.;2],&gate).is_err());
    let cancel=CancelGate::new_clock_free();cancel.request();
    assert!(actual.step(&[0.;3],&cancel).is_err());
    assert_eq!(actual.state(),state);assert_eq!(*actual.frame(),frame);
    for _ in 0..2 {actual.step(&[0.;3],&gate).unwrap();baseline.step(&[0.;3],&gate).unwrap();}
    assert!(actual.step(&[0.;3],&gate).is_err());
    assert_eq!(actual.remaining_steps(),0);
    actual.extend_step_budget(9).unwrap();
    for _ in 0..4 {actual.step(&[0.;3],&gate).unwrap();baseline.step(&[0.;3],&gate).unwrap();}
    assert_eq!(actual.frame(),baseline.frame());
    assert_eq!(actual.state().iter().map(|x|x.to_bits()).collect::<Vec<_>>(),
        baseline.state().iter().map(|x|x.to_bits()).collect::<Vec<_>>());
    assert!(actual.extend_step_budget(9).is_err());
}

#[test]
fn malformed_or_unsupported_input_is_not_repaired_or_silently_dropped() {
    let gate=CancelGate::new_clock_free();
    let mut drag=body(0.,0.,0.);drag.damping_per_s[0]=1.0;
    assert!(LinearImpactSystem::new(vec![drag],vec![],vec![],config(192000,1),&gate).is_err());
    let (bodies,_,mut volume)=physical_input(); volume.areas[0]=0.01;
    assert!(LinearImpactSystem::new(bodies.clone(),vec![],vec![VolumeConnection {spring:volume,
        reference_area_m2:0.02}],config(192000,1),&gate).is_err(),"a three-body volume is unsupported, not reduced to two");
    let raw=Obstacle::from_raw_parts(vec![1.0],1,vec![],vec![],1e7,1.5,"malformed".into());
    assert!(LinearImpactSystem::new(bodies.clone(),vec![raw],vec![],config(192000,1),&gate).is_err());
    let all=Obstacle::new(vec![1.,-1.,-1.],1,3,vec![0.],vec![1.],1e7,1.5,"three bodies".into()).unwrap();
    assert!(LinearImpactSystem::new(bodies,vec![all],vec![],config(192000,1),&gate).is_err());
}
