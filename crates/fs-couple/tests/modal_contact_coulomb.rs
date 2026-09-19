//! Real shared modal bodies under the 2-D set-valued Coulomb graph.
use fs_couple::modal_acoustic_time::{ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeModel};
use fs_couple::render::schedule::force::coupled::{CoupledModalSystem, ModalAttachment, ModalConnection, ModalCouplingConfig, ModalCouplingError};
use fs_couple::render::schedule::force::coupled::contact::{ModalContact, ModalContactConfig};
use fs_couple::render::schedule::force::coupled::contact::multiple::{MultiContactConfig, MultiContactModalSystem};
use fs_couple::render::schedule::force::coupled::contact::multiple::coulomb::{ModalCoulombFriction, CoulombContactRegime};
use fs_couple::render::schedule::force::coupled::contact::multiple::friction::ModalFriction;
use fs_exec::CancelGate;
use fs_dcontact::Obstacle;
use fs_math::c64::C64;
use fs_tribo::{ContactFrame, FrictionLaw, InputAuthority, InterfaceMedium, InterfaceSystemRef, TangentialSlip};

fn close(a: f64, b: f64, tolerance: f64) { assert!((a-b).abs() <= tolerance, "{a:e} != {b:e} within {tolerance:e}"); }
fn magnitude(v: [f64; 2]) -> f64 { v[0].hypot(v[1]) }
fn interface() -> InterfaceSystemRef {
    InterfaceSystemRef::new("synthetic:ordered-pair", "history:stateless", "authored-coulomb-test",
        InputAuthority::SyntheticFixture, InterfaceMedium::Dry).unwrap()
}
fn spec(mu: f64, mixed: f64) -> ModalCoulombFriction {
    ModalCoulombFriction { left_shapes: [vec![mixed,1.0,0.0],vec![0.0,0.0,1.0]],
        right_shapes: [vec![mixed,1.0,0.0],vec![0.0,0.0,1.0]], coefficient:mu,
        maximum_force_n:1e4, velocity_tolerance_m_s:1e-8, interface:interface() }
}
fn limits() -> MultiContactConfig { MultiContactConfig { max_contacts:4, max_sweeps:128, max_setup_terms:4096 } }
fn normal_config() -> ModalContactConfig {
    ModalContactConfig { max_iterations:128, maximum_force_n:1e5, maximum_penetration_m:0.01,
        force_absolute_tolerance_n:1e-10, force_relative_tolerance:1e-10 }
}
fn model(normal_q: f64, velocity: [f64;2], last_omega: f64) -> ModalAcousticTimeModel {
    let mut m=ModalAcousticTimeModel::try_new(48000,[800.0,600.0,last_omega].into_iter().map(|w| ModalAcousticMode {
        angular_frequency_rad_s:w,damping_ratio:0.0,pressure_per_modal_velocity:C64::new(1.0,0.0),
    }).collect(), ModalAcousticTimeBudget::audible_reference()).unwrap();
    m.restore_states(&[ModalAcousticState { displacement_m_sqrt_kg:normal_q,velocity_m_sqrt_kg_per_s:0.0 },
        ModalAcousticState { displacement_m_sqrt_kg:0.0,velocity_m_sqrt_kg_per_s:velocity[0] },
        ModalAcousticState { displacement_m_sqrt_kg:0.0,velocity_m_sqrt_kg_per_s:velocity[1] }]).unwrap(); m
}
fn attachment(component:usize, shapes:[f64;3]) -> ModalAttachment { ModalAttachment {component,shapes:shapes.to_vec()} }
fn normal(a:usize,b:usize,gap:f64) -> ModalContact {
    ModalContact {left:attachment(a,[1.0,0.0,0.0]),right:attachment(b,[1.0,0.0,0.0]),
        law:Obstacle::new(vec![-1.0],1,1,vec![gap],vec![1.0],1e5,1.0,"synthetic-coulomb-normal".into()).unwrap()}
}
fn network(models:Vec<ModalAcousticTimeModel>, contacts:Vec<ModalContact>, links:Vec<ModalConnection>,
    limits:MultiContactConfig, config:ModalContactConfig) -> MultiContactModalSystem
{
    let base=CoupledModalSystem::new(models,links,ModalCouplingConfig {max_modes:16,max_connections:4,max_setup_terms:4096,
        nyquist_guard_fraction:0.9,maximum_total_energy_j:100.0,maximum_abs_pressure_pa:1e6,
        maximum_abs_connection_force_n:1e6,solve_relative_tolerance:1e-11,
        energy_absolute_tolerance_j:1e-10,energy_relative_tolerance:1e-8},&CancelGate::new()).unwrap();
    MultiContactModalSystem::new(base,contacts.into_iter().map(|c|(c,config)).collect(),limits,&CancelGate::new()).unwrap()
}
fn pair(velocity:[f64;2],equilibrium:bool,gap:f64,limits:MultiContactConfig,config:ModalContactConfig) -> MultiContactModalSystem {
    let force=if equilibrium {100.0/(1.0+2e5/(800.0*800.0))} else {0.0};
    network(vec![model(-force/(800.0*800.0),velocity,600.0),model(force/(800.0*800.0),[0.0;2],600.0)],
        vec![normal(0,1,gap)],vec![],limits,config)
}
fn shared(reverse:bool) -> MultiContactModalSystem {
    let mut contacts=vec![normal(0,1,-0.0006),normal(1,2,-0.0006)];if reverse {contacts.reverse();}
    network(vec![model(0.0,[1.0,0.2],900.0),model(0.0,[0.2,-0.1],900.0),model(0.0,[-0.3,0.4],900.0)],contacts,
        vec![ModalConnection {left:attachment(0,[0.5,1.0,0.2]),right:attachment(2,[0.5,1.0,0.2]),
            stiffness_n_m:8e4,damping_n_s_m:0.4,rest_extension_m:0.0}],limits(),normal_config())
}
fn shared_specs() -> Vec<Option<ModalCoulombFriction>> {
    [(0.3,0.4),(0.5,-0.2)].into_iter().map(|(mu,mixed)| {
        let mut s=spec(mu,mixed);s.left_shapes[1]=vec![0.15,0.1,0.8];s.right_shapes[1]=s.left_shapes[1].clone();Some(s)
    }).collect()
}
fn states(s:&MultiContactModalSystem) -> Vec<ModalAcousticState> {s.components().iter().flat_map(|m|m.states().iter().copied()).collect()}
fn q(x:&[ModalAcousticState],i:usize) -> f64 {x[i].displacement_m_sqrt_kg}

#[test]
fn subthreshold_two_direction_load_is_held_without_regularized_creep() {
    let mut system=pair([0.0;2],true,-0.001,limits(),normal_config())
        .with_coulomb_friction(vec![Some(spec(0.3,0.0))],&CancelGate::new()).unwrap();
    let initial=states(&system);
    for _ in 0..200 {
        let f=system.step(&[0.0,3.0,4.0,0.0,-3.0,-4.0]).unwrap();let t=f.coulomb_friction[0].unwrap();
        assert_eq!(t.regime,CoulombContactRegime::Sticking);
        close(t.reaction_n[0],3.0,1e-8);close(t.reaction_n[1],4.0,1e-8);
        assert!(magnitude(t.reaction_n)<t.capacity_n);
        assert!(magnitude(t.slip_velocity_m_s)<1e-10);
        assert!(t.dissipation_j<1e-12 && t.work_residual_j.abs()<=t.work_tolerance_j);
        assert!(f.friction.is_empty());
    }
    for (a,b) in states(&system).iter().zip(initial) {
        close(a.displacement_m_sqrt_kg,b.displacement_m_sqrt_kg,1e-12);
        close(a.velocity_m_sqrt_kg_per_s,b.velocity_m_sqrt_kg_per_s,1e-9);
    }
}

#[test]
fn diagonal_sliding_matches_the_circular_capacity_not_a_friction_pyramid() {
    for sign in [-1.0,1.0] {
        let mut s=pair([0.6*sign,0.8*sign],false,-0.001,limits(),normal_config())
            .with_coulomb_friction(vec![Some(spec(0.3,0.0))],&CancelGate::new()).unwrap();
        let dt=s.sample_period_s();let dn=(1.0-(800.0*dt).cos())/(800.0*800.0);
        let expected_normal=100.0/(1.0+1e5*dn);
        let f=s.step(&[0.0;6]).unwrap().clone();let t=f.coulomb_friction[0].unwrap();
        close(f.contacts[0].normal_force_n,expected_normal,2e-8);
        close(t.reaction_n[0],0.6*sign*0.3*expected_normal,2e-8);
        close(t.reaction_n[1],0.8*sign*0.3*expected_normal,2e-8);
        close(magnitude(t.reaction_n),t.capacity_n,1e-9);
        assert_eq!(t.regime,CoulombContactRegime::Sliding);
        let frame=ContactFrame::new([0.0,0.0,1.0]).unwrap();
        let slip=TangentialSlip::new(&frame,[t.slip_velocity_m_s[0],t.slip_velocity_m_s[1],0.0]).unwrap();
        let owner=FrictionLaw::Coulomb {static_mu:0.3,kinetic_mu:0.3}.evaluate(&interface(),expected_normal,slip).unwrap();
        close(t.reaction_n[0],-owner.traction_n()[0],2e-8);close(t.reaction_n[1],-owner.traction_n()[1],2e-8);
        close(t.dissipation_j,owner.dissipated_power_w()*dt,1e-11);
    }
}

#[test]
fn overload_transitions_from_sticking_to_sliding_with_resolved_two_axis_work() {
    let mut s=pair([0.0;2],true,-0.001,limits(),normal_config())
        .with_coulomb_friction(vec![Some(spec(0.3,0.0))],&CancelGate::new()).unwrap();
    for _ in 0..8 {assert_eq!(s.step(&[0.0,3.0,4.0,0.0,-3.0,-4.0]).unwrap().coulomb_friction[0].unwrap().regime,CoulombContactRegime::Sticking);}
    let t=s.step(&[0.0,30.0,40.0,0.0,-30.0,-40.0]).unwrap().coulomb_friction[0].unwrap();
    assert_eq!(t.regime,CoulombContactRegime::Sliding);assert!(t.dissipation_j>0.0);
    assert!(t.slip_velocity_m_s.iter().all(|v|*v>0.0));
    close(magnitude(t.reaction_n),t.capacity_n,1e-8);
}

#[test]
fn shared_body_tangents_feed_back_into_normals_and_close_independently_reconstructed_work() {
    let specifications=shared_specs();
    let mut s=shared(false).with_coulomb_friction(specifications.clone(),&CancelGate::new()).unwrap();
    let mut reference=shared(false);let initial=s.total_energy_j().unwrap();
    let(mut work,mut loss,mut feedback,mut both)=(0.0,0.0,0.0_f64,false);
    for sample in 0..240 {
        let old=states(&s);let mut forces=[0.0;9];forces[1]=if sample<40 {0.5}else{0.0};forces[8]=-0.25;
        let f=s.step(&forces).unwrap().clone();let base=reference.step(&forces).unwrap();let new=states(&s);
        let external: f64=(0..9).map(|i|forces[i]*(q(&new,i)-q(&old,i))).sum();
        close(f.external_work_j,external,1e-12);
        let mut friction_loss=0.0;
        for i in 0..2 {
            let t=f.coulomb_friction[i].unwrap();let spec=specifications[i].as_ref().unwrap();let mut delta=[0.0;2];
            for axis in 0..2 {for mode in 0..3 {
                delta[axis]+=spec.left_shapes[axis][mode]*(q(&new,3*i+mode)-q(&old,3*i+mode))
                    -spec.right_shapes[axis][mode]*(q(&new,3*i+3+mode)-q(&old,3*i+3+mode));
            }}
            close(t.slip_velocity_m_s[0],delta[0]/s.sample_period_s(),1e-12);
            close(t.slip_velocity_m_s[1],delta[1]/s.sample_period_s(),1e-12);
            close(t.capacity_n,spec.coefficient*f.contacts[i].normal_force_n,1e-12);
            close(t.removed_work_j,t.reaction_n[0]*delta[0]+t.reaction_n[1]*delta[1],1e-12);
            close(t.dissipation_j,t.capacity_n*magnitude(delta),1e-12);
            assert!(t.graph_residual_n<=t.force_tolerance_n && t.work_residual_j.abs()<=t.work_tolerance_j);
            assert!(t.dissipation_j>=0.0);friction_loss+=t.dissipation_j;
            feedback=feedback.max((f.contacts[i].normal_force_n-base.contacts[i].normal_force_n).abs());
        }
        both|=f.contacts.iter().all(|p|p.normal_force_n>1.0)&&f.coulomb_friction.iter().all(|p|magnitude(p.unwrap().reaction_n)>0.1);
        close(f.friction_dissipation_j,friction_loss,1e-12);
        assert!(f.energy_residual_j.abs()<=f.energy_tolerance_j);
        work+=external;loss+=f.network_dissipation_j+f.contact_dissipation_j+friction_loss;
    }
    assert!(both&&feedback>1e-5&&loss>1e-3);close(s.total_energy_j().unwrap()-initial+loss-work,0.0,1e-7);
}

#[test]
fn arbitrary_tangent_rotation_and_contact_permutation_preserve_physics() {
    let original=shared_specs();let mut rotated=original.clone();let angle=0.731_f64;
    for s in rotated.iter_mut().flatten() {for shapes in [&mut s.left_shapes,&mut s.right_shapes] {
        let old=shapes.clone();for mode in 0..3 {
            shapes[0][mode]=angle.cos()*old[0][mode]+angle.sin()*old[1][mode];
            shapes[1][mode]=-angle.sin()*old[0][mode]+angle.cos()*old[1][mode];
        }
    }}
    let mut reversed=original.clone();reversed.reverse();
    let mut a=shared(false).with_coulomb_friction(original,&CancelGate::new()).unwrap();
    let mut b=shared(false).with_coulomb_friction(rotated,&CancelGate::new()).unwrap();
    let mut c=shared(true).with_coulomb_friction(reversed,&CancelGate::new()).unwrap();
    for _ in 0..100 {
        let fa=a.step(&[0.0;9]).unwrap().clone();let fb=b.step(&[0.0;9]).unwrap().clone();c.step(&[0.0;9]).unwrap();
        for other in [&b,&c] {for (x,y) in states(&a).iter().zip(states(other)) {
            close(x.displacement_m_sqrt_kg,y.displacement_m_sqrt_kg,1e-10);close(x.velocity_m_sqrt_kg_per_s,y.velocity_m_sqrt_kg_per_s,1e-8);
        }}
        for i in 0..2 {
            let x=fa.coulomb_friction[i].unwrap();let y=fb.coulomb_friction[i].unwrap();
            close(angle.cos()*x.reaction_n[0]+angle.sin()*x.reaction_n[1],y.reaction_n[0],1e-7);
            close(-angle.sin()*x.reaction_n[0]+angle.cos()*x.reaction_n[1],y.reaction_n[1],1e-7);
            close(x.dissipation_j,y.dissipation_j,1e-10);
        }
    }
}

#[test]
fn zero_capacity_and_explicit_absence_keep_original_normal_states() {
    for (gap,coefficient) in [(1.0,0.3),(-0.001,0.0)] {
        let make=||pair([0.6,0.8],false,gap,limits(),normal_config());let mut baseline=make();
        let mut s=make().with_coulomb_friction(vec![Some(spec(coefficient,0.2))],&CancelGate::new()).unwrap();
        let mut absent=make().with_coulomb_friction(vec![None],&CancelGate::new()).unwrap();
        for _ in 0..80 {
            let reference=baseline.step(&[0.0;6]).unwrap().clone();let f=s.step(&[0.0;6]).unwrap().clone();absent.step(&[0.0;6]).unwrap();
            assert_eq!(states(&s),states(&baseline));assert_eq!(states(&absent),states(&baseline));assert_eq!(reference.contacts,f.contacts);
            let t=f.coulomb_friction[0].unwrap();assert_eq!(t.regime,CoulombContactRegime::Inactive);assert_eq!(t.reaction_n,[0.0;2]);
            assert_eq!(f.friction_dissipation_j,0.0);assert!(absent.last_frame().unwrap().coulomb_friction[0].is_none());
        }
    }
}

#[test]
fn refused_trials_and_cancellation_preserve_clock_diagnostics_and_retry() {
    let make=||pair([0.6,0.8],false,-0.001,limits(),normal_config()).with_coulomb_friction(vec![Some(spec(0.3,0.2))],&CancelGate::new()).unwrap();
    let(mut s,mut reference)=(make(),make());for _ in 0..8{s.step(&[0.0;6]).unwrap();reference.step(&[0.0;6]).unwrap();}
    let old=states(&s);let report=s.last_frame().cloned();let count=s.samples_rendered();let gate=CancelGate::new();gate.request();
    assert!(matches!(s.step_under_gate(&[0.0;6],&gate),Err(ModalCouplingError::Cancelled)));
    for f in [vec![],vec![f64::NAN;6],vec![1e100;6]] {
        assert!(s.step(&f).is_err());assert_eq!(states(&s),old);assert_eq!(s.last_frame(),report.as_ref());assert_eq!(s.samples_rendered(),count);
    }
    s.step_under_gate(&[0.0;6],&CancelGate::new()).unwrap();reference.step(&[0.0;6]).unwrap();
    assert_eq!(states(&s),states(&reference));assert_eq!(s.last_frame(),reference.last_frame());
    let mut low=spec(0.3,0.0);low.maximum_force_n=1e-12;
    let mut capped=pair([0.6,0.8],false,-0.001,limits(),normal_config()).with_coulomb_friction(vec![Some(low)],&CancelGate::new()).unwrap();
    let before=states(&capped);assert!(capped.step(&[0.0;6]).is_err());assert_eq!(states(&capped),before);assert!(capped.last_frame().is_none());
    // A starved joint root budget must not publish any component or report.
    let mut c=normal_config();c.max_iterations=1;
    let mut starved=pair([0.0;2],true,-0.001,limits(),c).with_coulomb_friction(vec![Some(spec(0.3,0.0))],&CancelGate::new()).unwrap();
    let before=states(&starved);assert!(starved.step(&[0.0,3.0,4.0,0.0,-3.0,-4.0]).is_err());
    assert_eq!(states(&starved),before);assert_eq!(starved.samples_rendered(),0);
}

#[test]
fn redundant_contacts_and_mixed_absence_do_not_require_global_invertibility() {
    let base=network(vec![model(0.0,[0.6,0.8],600.0),model(0.0,[0.0;2],600.0)],
        vec![normal(0,1,-0.001),normal(0,1,-0.001)],vec![],limits(),normal_config());
    let mut s=base.with_coulomb_friction(vec![Some(spec(0.3,0.0)),Some(spec(0.5,0.0))],&CancelGate::new()).unwrap();
    let f=s.step(&[0.0;6]).unwrap();assert!(f.contacts.iter().all(|c|c.normal_force_n>0.0));
    assert!(f.coulomb_friction.iter().all(|t|t.unwrap().dissipation_j>0.0));
    let mut s=shared(false).with_coulomb_friction(vec![None,Some(spec(0.3,0.0))],&CancelGate::new()).unwrap();
    let f=s.step(&[0.0;9]).unwrap();assert!(f.coulomb_friction[0].is_none()&&f.coulomb_friction[1].is_some());
}

#[test]
fn admission_checks_plane_rank_complete_maps_source_and_exact_setup_budget() {
    let valid=spec(0.3,0.0);let make=||pair([0.0;2],false,-0.001,limits(),normal_config());
    for bad in [ModalCoulombFriction {coefficient:f64::NAN,..valid.clone()},
        ModalCoulombFriction {coefficient:-0.1,..valid.clone()},
        ModalCoulombFriction {maximum_force_n:0.0,..valid.clone()},
        ModalCoulombFriction {velocity_tolerance_m_s:0.0,..valid.clone()},
        ModalCoulombFriction {velocity_tolerance_m_s:f64::INFINITY,..valid.clone()},
        ModalCoulombFriction {left_shapes:[vec![],vec![0.0,0.0,1.0]],..valid.clone()},
        ModalCoulombFriction {left_shapes:[vec![0.0,f64::NAN,0.0],vec![0.0,0.0,1.0]],..valid.clone()},
        ModalCoulombFriction {left_shapes:[vec![0.0,1.0,0.0],vec![0.0,1.0,0.0]],
            right_shapes:[vec![0.0,1.0,0.0],vec![0.0,1.0,0.0]],..valid.clone()}] {
        assert!(make().with_coulomb_friction(vec![Some(bad)],&CancelGate::new()).is_err());
    }
    assert!(InterfaceSystemRef::new("x","y","",InputAuthority::CallerDeclared,InterfaceMedium::Dry).is_err());
    assert!(make().with_coulomb_friction(vec![],&CancelGate::new()).is_err());
    // n=6,p=1,k=0: 5*(6*2+1)+9*6=119, combined original + new setup.
    for cap in [118,119] {
        let l=MultiContactConfig {max_setup_terms:cap,..limits()};
        assert_eq!(pair([0.0;2],false,-0.001,l,normal_config()).with_coulomb_friction(vec![Some(valid.clone())],&CancelGate::new()).is_ok(),cap==119);
    }
    let mut advanced=make();advanced.step(&[0.0;6]).unwrap();assert!(advanced.with_coulomb_friction(vec![Some(valid.clone())],&CancelGate::new()).is_err());
    let attached=make().with_coulomb_friction(vec![Some(valid.clone())],&CancelGate::new()).unwrap();
    assert_eq!(attached.coulomb_friction_law(0),Some(&valid));assert!(attached.friction_law(0).is_none());
    assert!(attached.with_friction(vec![Some(ModalFriction {left_shapes:vec![0.0,1.0,0.0],right_shapes:vec![0.0,1.0,0.0],
        coefficient:0.3,regularization_speed_m_s:0.01,maximum_force_n:1e4,source:"old-rung".into()})],&CancelGate::new()).is_err());
    let gate=CancelGate::new();gate.request();
    assert!(matches!(make().with_coulomb_friction(vec![Some(valid)],&gate),Err(ModalCouplingError::Cancelled)));
}
