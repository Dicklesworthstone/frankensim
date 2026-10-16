//! Real two-body engagement/separation through the existing modal/contact owners.
use fs_couple::modal_acoustic_time::{ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeModel};
use fs_couple::render::schedule::force::coupled::{CoupledModalSystem, ModalAttachment, ModalConnection, ModalCouplingConfig, ModalCouplingError};
use fs_couple::render::schedule::force::coupled::contact::{ContactModalSystem, ModalContact, ModalContactConfig};
use fs_dcontact::Obstacle;
use fs_exec::CancelGate;
use fs_math::c64::C64;

fn budgets() -> ModalCouplingConfig {
    ModalCouplingConfig { max_modes:8, max_connections:4, max_setup_terms:4096,
        nyquist_guard_fraction:0.9, maximum_total_energy_j:100.0, maximum_abs_pressure_pa:1e6,
        maximum_abs_connection_force_n:1e6, solve_relative_tolerance:1e-11,
        energy_absolute_tolerance_j:1e-10, energy_relative_tolerance:1e-8 }
}
fn contact_config() -> ModalContactConfig {
    ModalContactConfig { max_iterations:80, maximum_force_n:1e5, maximum_penetration_m:0.002,
        force_absolute_tolerance_n:1e-10, force_relative_tolerance:1e-10 }
}
fn image(rate:u32, omega:f64, damping:f64, q:f64, v:f64, gain:f64) -> ModalAcousticTimeModel {
    let mut m=ModalAcousticTimeModel::try_new(rate,vec![ModalAcousticMode {
        angular_frequency_rad_s:omega,damping_ratio:damping,pressure_per_modal_velocity:C64::new(gain,0.0)
    }],ModalAcousticTimeBudget::audible_reference()).unwrap();
    m.restore_states(&[ModalAcousticState {displacement_m_sqrt_kg:q,velocity_m_sqrt_kg_per_s:v}]).unwrap();m
}
fn models(rate:u32) -> Vec<ModalAcousticTimeModel> {
    vec![image(rate,800.0,0.0,0.0,1.0,0.0),image(rate,1100.0,0.0,0.0,0.0,1.0)]
}
fn attachment(component:usize) -> ModalAttachment { ModalAttachment {component,shapes:vec![1.0]} }
fn contact(k:f64, alpha:f64, chi:f64, gap:f64) -> ModalContact {
    ModalContact { left:attachment(0),right:attachment(1),law:Obstacle::new(vec![-1.0],1,1,
        vec![gap],vec![1.0],k,alpha,"authored:two-flexible-bodies".into()).unwrap().with_internal_loss(chi).unwrap() }
}
fn system(chi:f64) -> ContactModalSystem {
    ContactModalSystem::new(CoupledModalSystem::new(models(48000),vec![],budgets(),&CancelGate::new()).unwrap(),
        contact(3e8,1.5,chi,0.0005),contact_config(),&CancelGate::new()).unwrap()
}
fn states(s:&ContactModalSystem) -> Vec<ModalAcousticState> {
    s.components().iter().flat_map(|m|m.states().iter().copied()).collect()
}

#[test]
fn collision_drives_an_unforced_receiver_then_separates_without_attraction() {
    for chi in [0.0,0.2,50.0] {
        let mut s=system(chi);let initial=s.total_energy_j().unwrap();
        let (mut active,mut separated,mut peak,mut loss)=(false,false,0.0_f64,0.0);
        for _ in 0..1200 {
            let frame=*s.step(&[0.0,0.0]).unwrap();
            assert!(frame.normal_force_n>=0.0);
            assert!(frame.contact_dissipation_j>=0.0);
            assert!(frame.energy_residual_j.abs()<=frame.energy_tolerance_j);
            assert!(frame.constitutive_residual_n.abs()<=frame.force_tolerance_n);
            if frame.normal_force_n>0.0 {active=true;}
            else if active && frame.penetration_after_m==0.0 {separated=true;}
            else if !active {assert_eq!(frame.observer_pressure_pa,0.0);}
            loss+=frame.network_dissipation_j+frame.contact_dissipation_j;
            peak=peak.max(frame.observer_pressure_pa.abs());
        }
        assert!(active && separated && peak>0.1, "must touch, excite the receiver and release");
        assert!((s.total_energy_j().unwrap()-initial+loss).abs()<1e-7);
        if chi==0.0 {assert!((s.total_energy_j().unwrap()-initial).abs()<1e-7);}
        else {assert!(loss>0.01);}
    }
}

#[test]
fn disabling_contact_is_exactly_the_original_network_in_all_damping_regimes() {
    for zeta in [0.0,0.1,1.0,2.0] {
        let build=||vec![image(48000,800.0,zeta,0.0,1.0,0.5),image(48000,1100.0,zeta,0.0,0.0,1.0)];
        let connections=||vec![ModalConnection {left:attachment(0),right:attachment(1),
            stiffness_n_m:1e5,damping_n_s_m:3.0,rest_extension_m:0.0}];
        let mut reference=CoupledModalSystem::new(build(),connections(),budgets(),&CancelGate::new()).unwrap();
        let base=CoupledModalSystem::new(build(),connections(),budgets(),&CancelGate::new()).unwrap();
        let mut s=ContactModalSystem::new(base,contact(0.0,1.5,0.2,0.0005),contact_config(),&CancelGate::new()).unwrap();
        for i in 0..128 {
            let forces=if i<37 {[0.5,-0.25]}else{[0.0;2]};
            let expected=reference.step(&forces).unwrap().observer_pressure_pa;
            let actual=s.step(&forces).unwrap().observer_pressure_pa;
            assert_eq!(actual.to_bits(),expected.to_bits());
            for (a,b) in s.components().iter().zip(reference.components()) {assert_eq!(a.states(),b.states());}
        }
    }
}

#[test]
fn nonlinear_contact_is_solved_with_the_bilateral_network_not_against_free_bodies() {
    let build=||CoupledModalSystem::new(models(48000),vec![ModalConnection {
        left:attachment(0),right:attachment(1),stiffness_n_m:2e5,damping_n_s_m:5.0,rest_extension_m:0.0,
    }],budgets(),&CancelGate::new()).unwrap();
    let mut s=ContactModalSystem::new(build(),contact(3e8,1.5,0.5,0.0001),contact_config(),&CancelGate::new()).unwrap();
    let initial=s.total_energy_j().unwrap();let (mut work,mut loss,mut contacts)=(0.0,0.0,0);
    for i in 0..300 {
        let f=*s.step(&[if i<60 {3.0}else{0.0},-0.5]).unwrap();
        work+=f.external_work_j;loss+=f.network_dissipation_j+f.contact_dissipation_j;
        contacts+=usize::from(f.normal_force_n>0.0);
    }
    assert!(contacts>0);
    assert!((s.total_energy_j().unwrap()-initial+loss-work).abs()<1e-7);
}

#[test]
fn smooth_always_closed_contact_refines_to_the_independent_linear_solution() {
    let omega=800.0_f64;let stiffness=1e5;let gap=-0.001;
    let frequency=(omega*omega+2.0*stiffness).sqrt();let equilibrium=2.0*stiffness*gap/(frequency*frequency);
    let duration=0.002;let exact_q=0.5*equilibrium*(1.0-(frequency*duration).cos());
    let exact_v=0.5*equilibrium*frequency*(frequency*duration).sin();
    let mut errors=Vec::new();
    for rate in [24000,48000,96000] {
        let m=vec![image(rate,omega,0.0,0.0,0.0,0.0),image(rate,omega,0.0,0.0,0.0,1.0)];
        let base=CoupledModalSystem::new(m,vec![],budgets(),&CancelGate::new()).unwrap();
        let mut s=ContactModalSystem::new(base,contact(stiffness,1.0,0.0,gap),contact_config(),&CancelGate::new()).unwrap();
        for _ in 0..rate/500 {assert!(s.step(&[0.0;2]).unwrap().penetration_after_m>0.0);}
        let a=s.components()[0].states()[0];
        errors.push(((a.displacement_m_sqrt_kg-exact_q)*omega).hypot(a.velocity_m_sqrt_kg_per_s-exact_v));
    }
    assert!(errors[0]/errors[1]>3.7 && errors[1]/errors[2]>3.7,"{errors:?}");
}

#[test]
fn refused_and_cancelled_trials_leave_the_entire_contact_network_retryable() {
    let mut s=system(0.2);let mut reference=system(0.2);
    for _ in 0..32 {s.step(&[0.0;2]).unwrap();reference.step(&[0.0;2]).unwrap();}
    let before=states(&s);let report=s.last_frame().copied();
    let gate=CancelGate::new();gate.request();
    assert!(matches!(s.step_under_gate(&[0.0;2],&gate),Err(ModalCouplingError::Cancelled)));
    for f in [vec![],vec![f64::NAN,0.0],vec![1e100,0.0]] {
        assert!(s.step(&f).is_err());assert_eq!(states(&s),before);assert_eq!(s.last_frame().copied(),report);
    }
    s.step_under_gate(&[0.0;2],&CancelGate::new()).unwrap();reference.step(&[0.0;2]).unwrap();
    assert_eq!(states(&s),states(&reference));assert_eq!(s.last_frame(),reference.last_frame());
}

#[test]
fn exhausted_root_and_final_penetration_gates_do_not_publish_candidates() {
    for (iterations,max_pen) in [(1,0.002),(80,1e-8)] {
        let m=vec![image(48000,800.0,0.0,0.000499,1.0,0.0),image(48000,1100.0,0.0,0.0,0.0,1.0)];
        let base=CoupledModalSystem::new(m,vec![],budgets(),&CancelGate::new()).unwrap();
        let mut c=contact_config();c.max_iterations=iterations;c.maximum_penetration_m=max_pen;
        let mut s=ContactModalSystem::new(base,contact(3e8,1.5,0.2,0.0005),c,&CancelGate::new()).unwrap();
        let before=states(&s);assert!(s.step(&[0.0;2]).is_err());
        assert_eq!(states(&s),before);assert_eq!(s.samples_rendered(),0);assert!(s.last_frame().is_none());
    }
}

#[test]
fn untrusted_contact_shape_law_and_admission_budgets_refuse_without_panics() {
    let base=||CoupledModalSystem::new(models(48000),vec![],budgets(),&CancelGate::new()).unwrap();
    for law in [
        Obstacle::from_raw_parts(vec![-1.0],1,vec![],vec![1.0],1.0,1.5,"raw".into()),
        Obstacle::from_raw_parts(vec![1.0],1,vec![0.0],vec![1.0],1.0,1.5,"raw".into()),
        Obstacle::from_raw_parts(vec![-1.0],1,vec![0.0],vec![-1.0],1.0,1.5,"raw".into()),
        Obstacle::from_raw_parts(vec![-1.0],1,vec![0.0],vec![1.0],f64::NAN,1.5,"raw".into()),
        Obstacle::from_raw_parts(vec![-1.0],1,vec![0.0],vec![1.0],1.0,1.5,"".into()),
    ] {
        let mut spec=contact(1.0,1.5,0.0,0.0);spec.law=law;
        assert!(ContactModalSystem::new(base(),spec,contact_config(),&CancelGate::new()).is_err());
    }
    let mut spec=contact(1.0,1.5,0.0,0.0);spec.left.component=10;
    assert!(ContactModalSystem::new(base(),spec,contact_config(),&CancelGate::new()).is_err());
    let mut cfg=contact_config();cfg.max_iterations=129;
    assert!(ContactModalSystem::new(base(),contact(1.0,1.5,0.0,0.0),cfg,&CancelGate::new()).is_err());
    let mut advanced=base();advanced.step(&[0.0;2]).unwrap();
    assert!(ContactModalSystem::new(advanced,contact(1.0,1.5,0.0,0.0),contact_config(),&CancelGate::new()).is_err());
    assert_eq!(system(0.2).contact_law().provenance(),"authored:two-flexible-bodies");
}
