//! Same physical input through the prepared linear image and independent checks.
use fs_couple::modal_acoustic_time::{ModalAcousticState, ModalAcousticTimeBudget};
use fs_couple::render::plate::impact::{BodyPotential, ImpactBody, ImpactConfig, ImpactSystem, VolumeSpring};
use fs_couple::render::plate::impact::linear::{LinearImpactConfig, LinearImpactSystem, VolumeConnection};
use fs_couple::render::plate::impact::damping::ViscousDamper;
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
fn shared_volume_and_damper_keep_all_body_cross_terms_and_separate_budgets() {
    let gate = CancelGate::new_clock_free();
    let weights = [1.0, -2.0, 0.5];
    let mut q = [2e-4, -1e-4, 3e-4];
    let mut v = [0.08, -0.03, 0.05];
    let bodies = q.iter().zip(&v).map(|(&q,&v)| body(0.0,q,v)).collect();
    let mut limits = config(96000,512);
    // The individual energies are 0.0032, 0.00045 and 0.00125 J. Flattening
    // these bodies would wrongly reject their combined 0.0049 J at admission.
    limits.component.maximum_total_energy_j = 0.0035;
    let area = 0.02;
    let stiffness = 3200.0;
    let damping = 3.0;
    let volume = VolumeConnection { reference_area_m2: area,
        spring: VolumeSpring { bulk_modulus_pa: 80000.0, volume_m3: 0.01,
            areas: weights.iter().map(|b| area*b).collect() } };
    let damper = ViscousDamper { weights: weights.to_vec(), damping_n_s_m: damping };
    let mut system = LinearImpactSystem::new_with_dampers(bodies,vec![],vec![volume],
        vec![damper],limits,&gate).unwrap();
    let initial = system.frame().stored_energy_j;
    let dt = system.sample_period_s();
    let norm = weights.iter().map(|b| b*b).sum::<f64>();
    let denominator = 1.0 + 0.5*damping*norm*dt + 0.25*stiffness*norm*dt*dt;
    let mut loss = 0.0;
    for _ in 0..512 {
        // Independent scalar midpoint solution: x=b^T q has inverse mass
        // b^T b. All modal momentum changes lie along b; its orthogonal
        // complement is unchanged. Pairwise springs or diagonalized damping
        // would both violate this exact three-body solution.
        let x = weights.iter().zip(&q).map(|(b,q)| b*q).sum::<f64>();
        let speed = weights.iter().zip(&v).map(|(b,v)| b*v).sum::<f64>();
        let next_speed = ((1.0-0.5*damping*norm*dt-0.25*stiffness*norm*dt*dt)*speed
            - stiffness*norm*dt*x)/denominator;
        for i in 0..3 {
            let next_v = v[i] + weights[i]*(next_speed-speed)/norm;
            q[i] += 0.5*dt*(v[i]+next_v);
            v[i] = next_v;
        }
        let frame = system.step(&[0.0;3],&gate).unwrap();
        assert_eq!(frame.supplied_work_j,0.0);
        loss += frame.dissipated_energy_j;
        for i in 0..3 {
            assert!((system.state()[2*i]-q[i]).abs()<1e-13);
            assert!((system.state()[2*i+1]-v[i]).abs()<1e-12);
        }
        let x = weights.iter().zip(&q).map(|(b,q)| b*q).sum::<f64>();
        let energy = 0.5*v.iter().map(|v| v*v).sum::<f64>() + 0.5*stiffness*x*x;
        assert!((frame.stored_energy_j-energy).abs()<1e-13);
        assert!((frame.stored_energy_j+loss-initial).abs()<1e-12);
    }
    assert!(loss>1e-5,"the complete physical damper must dissipate energy");
    assert!((2.0*system.state()[1]+system.state()[3]-0.13).abs()<1e-12);
    assert!((-0.5*system.state()[1]+system.state()[5]-0.01).abs()<1e-12);
}

#[test]
fn multi_body_volume_and_damper_share_the_joint_flexible_contact_solve() {
    let gate = CancelGate::new_clock_free();
    let parts = vec![body(0.0,-0.0002,0.8),body(1200.0,0.0,0.0),
        body(2300.0,0.0,0.0),body(3100.0,0.0,0.0)];
    let grouped = ImpactBody {
        potential: BodyPotential::Linear(vec![0.0,1200.0,2300.0,3100.0]),
        initial: parts.iter().flat_map(|b| b.initial.iter().copied()).collect(),
        damping_per_s: parts.iter().flat_map(|b| b.damping_per_s.iter().copied()).collect(),
    };
    let make = |bodies| LinearImpactSystem::new_with_dampers(bodies,
        vec![Obstacle::new(vec![1.0,-1.0,0.4,0.0, 1.0,0.3,0.0,-0.7],2,4,
            vec![0.0;2],vec![1.0;2],1e7,1.5,"manufactured flexible cavity assembly".into()).unwrap()],
        vec![VolumeConnection { reference_area_m2:0.02,
            spring:VolumeSpring {bulk_modulus_pa:1.4e5,volume_m3:0.02,
                areas:vec![0.0002,0.01,-0.008,0.004]} }],
        vec![ViscousDamper {weights:vec![0.2,1.0,-0.5,0.7],damping_n_s_m:20.0}],
        config(192000,512),&gate).unwrap();
    let mut split = make(parts);
    let mut grouped = make(vec![grouped]);
    let initial = split.frame().stored_energy_j;
    let mut loss = 0.0;
    let mut peaks = [0.0_f64;3];
    for _ in 0..512 {
        let frame = split.step(&[0.0;4],&gate).unwrap();
        grouped.step(&[0.0;4],&gate).unwrap();
        for (a,b) in split.state().iter().zip(grouped.state()) { assert!((a-b).abs()<1e-12); }
        for (i,peak) in peaks.iter_mut().enumerate() { *peak=peak.max(split.state()[2*i+3].abs()); }
        loss += frame.dissipated_energy_j;
        assert!(frame.balance_residual_j.abs()<1e-9);
        assert!((frame.stored_energy_j+loss-initial).abs()<1e-9);
    }
    assert_eq!(split.contact_count(),2);
    assert!(peaks.iter().all(|p| *p>1e-4),"every elastic body must receive work: {peaks:?}");
    assert!(loss>0.0);
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

/// G1/G3: a contact row is physical geometry, independent of how the same
/// diagonal modes are partitioned into rigid, shaft and resonator bodies.
#[test]
fn flexible_contact_rows_preserve_all_bodies_and_joint_reactions() {
    let gate = CancelGate::new_clock_free();
    for points in [1, 2] {
        let parts = vec![body(0.0,-0.0002,0.8),body(1200.0,0.0,0.0),body(2300.0,0.0,0.0)];
        let grouped = ImpactBody {
            potential: BodyPotential::Linear(vec![0.0,1200.0,2300.0]),
            initial: parts.iter().flat_map(|b| b.initial.iter().copied()).collect(),
            damping_per_s: parts.iter().flat_map(|b| b.damping_per_s.iter().copied()).collect(),
        };
        let rows = [1.0,-1.0,0.4, 1.0,0.3,-0.7];
        let contact = || Obstacle::new(rows[..3*points].to_vec(),points,3,vec![0.0;points],
            vec![1.0;points],1e7,1.5,"manufactured rigid/shaft/head contact".into()).unwrap();
        let mut split = LinearImpactSystem::new(parts,vec![contact()],vec![],config(192000,512),&gate).unwrap();
        let mut single = LinearImpactSystem::new(vec![grouped],vec![contact()],vec![],config(192000,512),&gate).unwrap();
        assert_eq!(split.contact_count(),points);
        let initial = split.frame().stored_energy_j;
        let mut loss = 0.0;
        let mut peaks = [0.0_f64;2];
        for tick in 0..512 {
            if tick == 100 {
                let state = split.state().to_vec(); let frame = *split.frame();
                let stopped = CancelGate::new_clock_free(); stopped.request();
                assert!(split.step(&[0.0;3],&stopped).is_err());
                assert!(split.step(&[0.0,f64::NAN,0.0],&gate).is_err());
                assert_eq!(split.state(),state); assert_eq!(*split.frame(),frame);
            }
            let a = split.step(&[0.0;3],&gate).unwrap();
            let b = single.step(&[0.0;3],&gate).unwrap();
            assert_eq!(a.sample,b.sample);
            assert_eq!(a.supplied_work_j,0.0);
            assert!(a.balance_residual_j.abs()<1e-9);
            loss += a.dissipated_energy_j;
            for (x,y) in split.state().iter().zip(single.state()) {assert!((x-y).abs()<1e-12);}
            for (i,peak) in peaks.iter_mut().enumerate() {*peak=peak.max(split.state()[2*i+3].abs());}
        }
        assert!(peaks.iter().all(|v| *v>1e-4),"contact must excite both elastic bodies: {peaks:?}");
        assert!((split.frame().stored_energy_j+loss-initial).abs()<1e-9);
    }
}

#[test]
fn multi_body_contact_keeps_the_original_component_energy_budgets() {
    let gate=CancelGate::new_clock_free();
    let mut limits=config(192000,1);
    limits.component.maximum_total_energy_j=0.006;
    let contact=Obstacle::new(vec![1.0,-1.0,0.4],1,3,vec![0.0],vec![1.0],
        1e7,1.5,"separately budgeted rigid/shaft/head contact".into()).unwrap();
    let mut system=LinearImpactSystem::new(
        vec![body(0.0,0.0,0.1),body(1200.0,0.0,0.1),body(2300.0,0.0,0.1)],
        vec![contact],vec![],limits,&gate).unwrap();
    // Each body has 0.005 J. Combining their modes into a single component
    // would wrongly reject the same admitted physical state at 0.015 J.
    assert!((system.frame().stored_energy_j-0.015).abs()<1e-16);
    system.step(&[0.0;3],&gate).unwrap();
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
    let mut drag=body(0.,0.,0.);drag.damping_per_s[0]=-1.0;
    assert!(LinearImpactSystem::new(vec![drag],vec![],vec![],config(192000,1),&gate).is_err());
    let (bodies,_,mut volume)=physical_input(); volume.areas.pop();
    assert!(LinearImpactSystem::new(bodies.clone(),vec![],vec![VolumeConnection {spring:volume,
        reference_area_m2:0.02}],config(192000,1),&gate).is_err(),"a volume must declare every original coordinate");
    let raw=Obstacle::from_raw_parts(vec![1.0],1,vec![],vec![],1e7,1.5,"malformed".into());
    assert!(LinearImpactSystem::new(bodies.clone(),vec![raw],vec![],config(192000,1),&gate).is_err());
    let zero=Obstacle::new(vec![0.;3],1,3,vec![0.],vec![1.],1e7,1.5,"no moving body".into()).unwrap();
    assert!(LinearImpactSystem::new(bodies,vec![zero],vec![],config(192000,1),&gate).is_err());
}
