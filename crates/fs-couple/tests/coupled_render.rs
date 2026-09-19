use fs_couple::modal_acoustic_time::{ModalAcousticMode, ModalAcousticTimeBudget, ModalAcousticTimeModel};
use fs_couple::render::{ControlDelta, RenderContext, RenderError, RenderVoice, GatedRenderOutcome};
use fs_couple::render::schedule::ScheduledRenderer;
use fs_couple::render::schedule::force::{ForceInitialization, ForceRenderConfig, ModalForceVoice, ModalForceEvent};
use fs_couple::render::schedule::force::coupled::{CoupledModalSystem, ModalAttachment, ModalConnection, ModalCouplingConfig};
use fs_couple::render::schedule::force::coupled::render::CoupledModalVoice;
use fs_exec::CancelGate;
use fs_math::c64::C64;

fn config() -> ModalCouplingConfig {
    ModalCouplingConfig { max_modes: 8, max_connections: 4, max_setup_terms: 4096,
        nyquist_guard_fraction: 0.9, maximum_total_energy_j: 1000.0,
        maximum_abs_pressure_pa: 1e6, maximum_abs_connection_force_n: 1e8,
        solve_relative_tolerance: 1e-11, energy_absolute_tolerance_j: 1e-12,
        energy_relative_tolerance: 1e-10 }
}
fn models() -> Vec<ModalAcousticTimeModel> {
    [800.0,1100.0].into_iter().enumerate().map(|(i,w)| ModalAcousticTimeModel::try_new(48_000,
        vec![ModalAcousticMode { angular_frequency_rad_s: w, damping_ratio: 0.03,
            pressure_per_modal_velocity: C64::new(i as f64,0.0) }],
        ModalAcousticTimeBudget::audible_reference()).unwrap()).collect()
}
fn connections() -> Vec<ModalConnection> {
    vec![ModalConnection { left: ModalAttachment { component: 0, shapes: vec![1.0] },
        right: ModalAttachment { component: 1, shapes: vec![1.0] },
        stiffness_n_m: 3e5, damping_n_s_m: 10.0, rest_extension_m: 0.0 }]
}
fn system() -> CoupledModalSystem {
    CoupledModalSystem::new(models(),connections(),config(),&CancelGate::new()).unwrap()
}
fn renderer() -> ScheduledRenderer {
    let voices = models().into_iter().enumerate().map(|(i,m)|
        ModalForceVoice::new(m,vec![vec![1.0]],vec![if i==0 { 100.0 } else { 0.0 }],
            ForceInitialization::StaticPreload).unwrap()).collect();
    ScheduledRenderer::from_coupled_modal_forces(voices,vec![
        ModalForceEvent { sample:37, voice:0, port:0, force_n:0.0 },
        ModalForceEvent { sample:71, voice:1, port:0, force_n:20.0 },
    ],ForceRenderConfig { sample_rate_hz:48000,max_block:512,max_events:2,max_controls:2,max_projection_terms:4 },
        connections(),config(),&CancelGate::new()).unwrap()
}

#[test]
fn static_preload_solves_the_full_connected_stiffness_not_independent_compliance() {
    let mut s = system();
    s.initialize_static_equilibrium(&[100.0,0.0],&CancelGate::new()).unwrap();
    let (a,b,k) = (800.0_f64.powi(2)+3e5,1100.0_f64.powi(2)+3e5,3e5_f64);
    let determinant = a*b-k*k;
    let expected = [100.0*b/determinant,100.0*k/determinant];
    for (m,q) in s.components().iter().zip(expected) {
        assert!((m.states()[0].displacement_m_sqrt_kg-q).abs() < 1e-18);
    }
    assert!(expected[1]>0.0, "the unforced receiver must be statically displaced");
    assert!((expected[0]-100.0/800.0_f64.powi(2)).abs()>1e-5);
    for _ in 0..100 {
        s.step(&[100.0,0.0]).unwrap();
        assert!(s.components().iter().all(|m| m.states()[0].velocity_m_sqrt_kg_per_s.abs()<1e-10));
    }
    let old:Vec<_> = s.components().iter().map(|m| m.states()[0]).collect();
    assert!(s.initialize_static_equilibrium(&[1.0,0.0],&CancelGate::new()).is_err());
    assert_eq!(s.components().iter().map(|m|m.states()[0]).collect::<Vec<_>>(),old);
}

#[test]
fn preload_refusal_preserves_every_component_and_allows_retry() {
    let mut c=config();c.maximum_total_energy_j=1e-4;
    let mut s=CoupledModalSystem::new(models(),connections(),c,&CancelGate::new()).unwrap();
    assert!(s.initialize_static_equilibrium(&[100.0,0.0],&CancelGate::new()).is_err());
    assert!(s.components().iter().all(|m|m.states()[0].displacement_m_sqrt_kg==0.0));
    assert_eq!(s.samples_rendered(),0);
    let gate=CancelGate::new();gate.request();
    assert!(s.initialize_static_equilibrium(&[1.0,0.0],&gate).is_err());
    s.initialize_static_equilibrium(&[1.0,0.0],&CancelGate::new()).unwrap();
    assert!(s.components()[1].states()[0].displacement_m_sqrt_kg>0.0);
}

#[test]
fn physical_schedule_preserves_coupled_mechanics_across_every_callback_partition() {
    let mut direct=system();
    direct.initialize_static_equilibrium(&[100.0,0.0],&CancelGate::new()).unwrap();
    let expected:Vec<_> = (0..257).map(|i|direct.step(&[if i<37 {100.0}else{0.0},if i>=71 {20.0}else{0.0}])
        .unwrap().observer_pressure_pa.to_bits()).collect();
    for block in [1,7,37,64,128,512] {
        let mut render=renderer();let mut out=vec![0.0;257];
        for chunk in out.chunks_mut(block) { render.block(chunk).unwrap(); }
        assert_eq!(out.iter().map(|v|v.to_bits()).collect::<Vec<_>>(),expected);
        assert_eq!(render.samples_rendered(),257);
        assert!(render.pending_controls().is_empty());
        assert!(render.validate_sample_rate(44_100).is_err());
    }
}

#[test]
fn cancellation_resumes_the_same_network_and_pending_controls() {
    let mut reference=renderer();let mut expected=[0.0;192];reference.block(&mut expected).unwrap();
    let mut render=renderer();let mut actual=[0.0;192];render.block(&mut actual[..64]).unwrap();
    let gate=CancelGate::new();gate.request();let mut untouched=[-7.0;64];
    assert_eq!(render.render_under_gate(&gate,&mut untouched,64,1).unwrap(),GatedRenderOutcome::Cancelled {blocks:0});
    assert_eq!(untouched,[-7.0;64]);assert_eq!(render.samples_rendered(),64);
    render.render_under_gate(&CancelGate::new(),&mut actual[64..],64,2).unwrap();
    assert_eq!(actual.map(f64::to_bits),expected.map(f64::to_bits));
}

#[test]
fn coupled_control_batches_are_atomic_and_physics_errors_poison_only_the_host() {
    let make=||RenderContext::new(vec![RenderVoice::CoupledModal(Box::new(
        CoupledModalVoice::new(system(),vec![0.0;2]).unwrap()))],32);
    let valid=ControlDelta::SetModalForce {voice:0,mode:0,force_n_per_sqrt_kg:100.0};
    for invalid in [ControlDelta::SetModalForce {voice:0,mode:2,force_n_per_sqrt_kg:1.0},
        ControlDelta::SetBlowingPressure {voice:0,pressure_pa:100.0},
        ControlDelta::SetPlateForce {voice:0,force_n:100.0}] {
        let mut context=make();assert!(context.apply_controls(&[valid,invalid]).is_err());
        assert!(context.control_log().is_empty());let mut out=[-1.0;32];context.block(&mut out).unwrap();
        assert!(out.iter().all(|x|*x==0.0));
    }
    let mut context=make();context.apply_controls(&[ControlDelta::SetModalForce {voice:0,mode:0,force_n_per_sqrt_kg:1e10}]).unwrap();
    assert!(matches!(context.block(&mut [0.0;32]),Err(RenderError::Coupled(_))));
    let mut untouched=[7.0;32];assert!(matches!(context.block(&mut untouched),Err(RenderError::Poisoned)));
    assert_eq!(untouched,[7.0;32]);
}

#[test]
fn mixed_preload_and_retained_vibration_are_not_silently_reinterpreted() {
    let voices=models().into_iter().enumerate().map(|(i,m)|ModalForceVoice::new(m,vec![vec![1.0]],vec![1.0],
        if i==0 {ForceInitialization::StaticPreload}else{ForceInitialization::RetainState}).unwrap()).collect();
    assert!(ScheduledRenderer::from_coupled_modal_forces(voices,vec![],ForceRenderConfig {
        sample_rate_hz:48000,max_block:64,max_events:0,max_controls:0,max_projection_terms:2,
    },connections(),config(),&CancelGate::new()).is_err());
}
