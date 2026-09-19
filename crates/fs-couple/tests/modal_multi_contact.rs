use fs_couple::modal_acoustic_time::{ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeModel};
use fs_couple::render::schedule::force::coupled::{CoupledModalSystem, ModalAttachment, ModalCouplingConfig, ModalCouplingError};
use fs_couple::render::schedule::force::coupled::contact::{ContactModalSystem, ModalContact, ModalContactConfig};
use fs_couple::render::schedule::force::coupled::contact::multiple::{MultiContactConfig, MultiContactModalSystem};
use fs_dcontact::Obstacle;
use fs_exec::CancelGate;
use fs_math::c64::C64;

fn coupling() -> ModalCouplingConfig {
    ModalCouplingConfig { max_modes: 16, max_connections: 4, max_setup_terms: 4096,
        nyquist_guard_fraction: 0.9, maximum_total_energy_j: 1000.0,
        maximum_abs_pressure_pa: 1e6, maximum_abs_connection_force_n: 1e6,
        solve_relative_tolerance: 1e-11, energy_absolute_tolerance_j: 1e-10, energy_relative_tolerance: 1e-8 }
}
fn limits() -> ModalContactConfig {
    ModalContactConfig { max_iterations: 96, maximum_force_n: 1e6, maximum_penetration_m: 0.1,
        force_absolute_tolerance_n: 1e-8, force_relative_tolerance: 1e-9 }
}
fn multiple() -> MultiContactConfig { MultiContactConfig { max_contacts: 8, max_sweeps: 128, max_setup_terms: 4096 } }
fn mode(rate: u32, w: f64, q: f64, v: f64, gain: f64) -> ModalAcousticTimeModel {
    let mut model = ModalAcousticTimeModel::try_new(rate, vec![ModalAcousticMode {
        angular_frequency_rad_s: w, damping_ratio: 0.0, pressure_per_modal_velocity: C64::new(gain,0.0),
    }], ModalAcousticTimeBudget::audible_reference()).unwrap();
    model.restore_states(&[ModalAcousticState { displacement_m_sqrt_kg: q, velocity_m_sqrt_kg_per_s: v }]).unwrap();
    model
}
fn network(models: Vec<ModalAcousticTimeModel>) -> CoupledModalSystem {
    CoupledModalSystem::new(models, vec![], coupling(), &CancelGate::new()).unwrap()
}
fn contact(a: usize, b: usize, gap: f64, k: f64, alpha: f64, chi: f64) -> (ModalContact, ModalContactConfig) {
    (ModalContact {
        left: ModalAttachment { component:a, shapes:vec![1.0] },
        right: ModalAttachment { component:b, shapes:vec![1.0] },
        law: Obstacle::new(vec![-1.0],1,1,vec![gap],vec![1.0],k,alpha,"authored multi-contact test".into())
            .unwrap().with_internal_loss(chi).unwrap(),
    },limits())
}
fn launch(k: f64) -> MultiContactModalSystem {
    let models=vec![mode(48000,300.0,0.0,1.0,0.0),mode(48000,800.0,0.0,0.0,1.0),mode(48000,450.0,0.0,-0.7,0.0)];
    MultiContactModalSystem::new(network(models),vec![contact(0,1,0.0001,k,1.5,0.12),
        contact(1,2,0.00012,k*0.8,1.5,0.2)],multiple(),&CancelGate::new()).unwrap()
}
fn states(system: &MultiContactModalSystem) -> Vec<ModalAcousticState> {
    system.components().iter().flat_map(|m|m.states()).copied().collect()
}

#[test]
fn two_simultaneous_contacts_drive_the_shared_unforced_receiver_and_close_total_energy() {
    let mut system=launch(3e6);let initial=system.total_energy_j().unwrap();
    let mut overlap=false;let mut separation=false;let mut peak=0.0_f64;let mut losses=0.0;
    for _ in 0..1201 {
        let frame=system.step(&[0.0;3]).unwrap();
        overlap |= frame.contacts.iter().all(|p|p.normal_force_n>1e-6);
        separation |= overlap && frame.contacts.iter().any(|p|p.normal_force_n==0.0);
        for point in &frame.contacts {
            assert!(point.normal_force_n>=0.0);
            assert!(point.dissipation_j>=0.0);
            assert!(point.constitutive_residual_n.abs()<=point.force_tolerance_n);
        }
        assert_eq!(frame.external_work_j,0.0);
        assert!(frame.energy_residual_j.abs()<=frame.energy_tolerance_j);
        peak=peak.max(frame.observer_pressure_pa.abs());
        losses+=frame.network_dissipation_j+frame.contact_dissipation_j;
    }
    assert!(overlap && separation);
    assert!(peak>0.05);
    assert!((system.total_energy_j().unwrap()+losses-initial).abs()<1e-8);
    let mut disabled=launch(0.0);
    for _ in 0..1201 { assert_eq!(disabled.step(&[0.0;3]).unwrap().observer_pressure_pa,0.0); }
}

#[test]
fn each_reaction_uses_the_same_actual_shared_body_endpoints() {
    let mut system=MultiContactModalSystem::new(network(vec![mode(48000,500.0,0.001,0.0,0.0),
        mode(48000,700.0,0.0,0.0,1.0),mode(48000,600.0,-0.002,0.0,0.0)]),
        vec![contact(0,1,0.0,1e8,1.0,0.0),contact(1,2,0.0,1e8,1.0,0.0)],multiple(),&CancelGate::new()).unwrap();
    let before=states(&system);
    let frame=system.step(&[0.0;3]).unwrap().clone();let after=states(&system);
    assert!(frame.sweeps>1,"one sequential collision pass is not a joint solution");
    for i in 0..2 {
        let x0=before[i].displacement_m_sqrt_kg-before[i+1].displacement_m_sqrt_kg;
        let x1=after[i].displacement_m_sqrt_kg-after[i+1].displacement_m_sqrt_kg;
        assert!(x0>0.0 && x1>0.0);
        let independent=1e8*(x0+x1)/2.0;
        assert!((frame.contacts[i].normal_force_n-independent).abs()<=frame.contacts[i].force_tolerance_n);
    }
}

#[test]
fn redundant_contact_maps_share_reactions_instead_of_requiring_invertible_contact_compliance() {
    let make=||network(vec![mode(48000,800.0,0.001,0.0,1.0),mode(48000,800.0,0.0,0.0,1.0)]);
    let (one,c)=contact(0,1,0.0,3e5,1.0,0.0);
    let mut single=ContactModalSystem::new(make(),one,c,&CancelGate::new()).unwrap();
    let mut duplicate=MultiContactModalSystem::new(make(),vec![contact(0,1,0.0,1.5e5,1.0,0.0),
        contact(0,1,0.0,1.5e5,1.0,0.0)],multiple(),&CancelGate::new()).unwrap();
    for _ in 0..128 {
        single.step(&[0.0;2]).unwrap();duplicate.step(&[0.0;2]).unwrap();
        for (a,b) in single.components().iter().zip(duplicate.components()) {
            assert!((a.states()[0].displacement_m_sqrt_kg-b.states()[0].displacement_m_sqrt_kg).abs()<1e-10);
            assert!((a.states()[0].velocity_m_sqrt_kg_per_s-b.states()[0].velocity_m_sqrt_kg_per_s).abs()<1e-7);
        }
    }
}

#[test]
fn fully_active_linear_contacts_converge_to_independent_chain_normal_mode() {
    let mut errors=Vec::new();let shifted=(800.0_f64.powi(2)+3e5).sqrt();let t=1.0/3000.0;
    for rate in [24000,48000,96000] {
        let models=vec![mode(rate,800.0,0.001,0.0,0.0),mode(rate,800.0,0.0,0.0,1.0),mode(rate,800.0,-0.001,0.0,0.0)];
        let mut s=MultiContactModalSystem::new(network(models),vec![contact(0,1,0.0,3e5,1.0,0.0),
            contact(1,2,0.0,3e5,1.0,0.0)],multiple(),&CancelGate::new()).unwrap();
        for _ in 0..rate/3000 { s.step(&[0.0;3]).unwrap(); }
        let expected_q=0.001*(shifted*t).cos();let expected_v=-0.001*shifted*(shifted*t).sin();
        let state=s.components()[0].states()[0];
        errors.push(((state.displacement_m_sqrt_kg-expected_q)*800.0).hypot(state.velocity_m_sqrt_kg_per_s-expected_v));
    }
    assert!(errors[0]/errors[1]>3.7 && errors[1]/errors[2]>3.7,"{errors:?}");
}

#[test]
fn nonconverged_and_cancelled_trials_do_not_publish_partial_contacts() {
    let mut system=launch(3e6);let mut reference=launch(3e6);
    system.step(&[0.0;3]).unwrap();reference.step(&[0.0;3]).unwrap();
    let before=states(&system);let last=system.last_frame().unwrap().clone();
    let gate=CancelGate::new();gate.request();
    assert!(matches!(system.step_under_gate(&[0.0;3],&gate),Err(ModalCouplingError::Cancelled)));
    for f in [vec![f64::NAN,0.0,0.0],vec![],vec![1e100,0.0,0.0]] {
        assert!(system.step(&f).is_err());assert_eq!(states(&system),before);
        assert_eq!(system.last_frame(),Some(&last));
    }
    system.step(&[0.0;3]).unwrap();reference.step(&[0.0;3]).unwrap();
    assert_eq!(states(&system),states(&reference));assert_eq!(system.last_frame(),reference.last_frame());
    let mut budget=multiple();budget.max_sweeps=1;
    let models=vec![mode(48000,500.0,0.001,0.0,0.0),mode(48000,700.0,0.0,0.0,1.0),mode(48000,600.0,-0.002,0.0,0.0)];
    let mut limited=MultiContactModalSystem::new(network(models),vec![contact(0,1,0.0,1e8,1.0,0.0),
        contact(1,2,0.0,1e8,1.0,0.0)],budget,&CancelGate::new()).unwrap();
    let old=states(&limited);
    assert!(matches!(limited.step(&[0.0;3]),Err(ModalCouplingError::ContactSolve {..})));
    assert_eq!(states(&limited),old);assert_eq!(limited.samples_rendered(),0);assert!(limited.last_frame().is_none());
}

#[test]
fn singleton_and_contact_order_permutations_agree_within_declared_solver_accuracy() {
    let make=||network(vec![mode(48000,300.0,0.0,1.0,0.0),mode(48000,800.0,0.0,0.0,1.0)]);
    let (point,c)=contact(0,1,0.0001,3e6,1.5,0.12);
    let mut single=ContactModalSystem::new(make(),point.clone(),c,&CancelGate::new()).unwrap();
    let mut multi=MultiContactModalSystem::new(make(),vec![(point,c)],multiple(),&CancelGate::new()).unwrap();
    for _ in 0..200 {
        let a=single.step(&[0.0;2]).unwrap().observer_pressure_pa;
        let b=multi.step(&[0.0;2]).unwrap().observer_pressure_pa;assert!((a-b).abs()<1e-7);
    }
    let make=||network(vec![mode(48000,300.0,0.0,1.0,0.0),mode(48000,800.0,0.0,0.0,1.0),mode(48000,450.0,0.0,-0.7,0.0)]);
    let points=vec![contact(0,1,0.0001,3e6,1.5,0.12),contact(1,2,0.00012,2.4e6,1.5,0.2)];
    let mut reversed=points.clone();reversed.reverse();
    let mut a=MultiContactModalSystem::new(make(),points,multiple(),&CancelGate::new()).unwrap();
    let mut b=MultiContactModalSystem::new(make(),reversed,multiple(),&CancelGate::new()).unwrap();
    for _ in 0..200 {
        let pa=a.step(&[0.0;3]).unwrap().observer_pressure_pa;
        let pb=b.step(&[0.0;3]).unwrap().observer_pressure_pa;assert!((pa-pb).abs()<1e-7);
    }
}

#[test]
fn contact_count_setup_and_physical_admission_precede_allocation_or_motion() {
    let make=||network(vec![mode(48000,300.0,0.0,0.0,0.0),mode(48000,800.0,0.0,0.0,1.0)]);
    for c in [MultiContactConfig {max_contacts:0,..multiple()},MultiContactConfig {max_contacts:33,..multiple()},
        MultiContactConfig {max_sweeps:0,..multiple()},MultiContactConfig {max_setup_terms:6,..multiple()}] {
        assert!(MultiContactModalSystem::new(make(),vec![contact(0,1,0.0,3e6,1.5,0.0)],c,&CancelGate::new()).is_err());
    }
    // n=2, p=1, k=0: 1*(2*2+1)+2*1=7 setup terms.
    let c=MultiContactConfig {max_setup_terms:7,..multiple()};
    assert!(MultiContactModalSystem::new(make(),vec![contact(0,1,0.0,3e6,1.5,0.0)],c,&CancelGate::new()).is_ok());
    let (mut malformed,c)=contact(0,1,0.0,1.0,1.0,0.0);malformed.right.shapes.push(1.0);
    assert!(MultiContactModalSystem::new(make(),vec![(malformed,c)],multiple(),&CancelGate::new()).is_err());
    assert!(MultiContactModalSystem::new(make(),vec![],multiple(),&CancelGate::new()).is_err());
}
