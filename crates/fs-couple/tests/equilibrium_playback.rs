//! Fitted state must survive the transition from static design to audio.
use fs_couple::modal_acoustic_time::{ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeModel};
use fs_couple::render::schedule::force::ForceRenderConfig;
use fs_couple::render::schedule::force::file::{ModalPerformance, design::EquilibriumDesignFile};
use fs_couple::render::schedule::force::coupled::equilibrium::sensitivity::objective::design::{DesignControl, DesignWork};
use fs_couple::render::schedule::force::coupled::equilibrium::sensitivity::objective::design::forward::playback::*;
use fs_couple::render::GatedRenderOutcome;
use fs_exec::CancelGate;
use fs_math::c64::C64;

const MODEL: &str = "frankensim-modal-performance-v2\nsample_rate_hz 48000\nsamples 1\nfull_scale_pa 1\nlimits 0.9 10 1000 1000 10000\ncompile_limits 0 1\nvoices 1\nvoice retain-state 1 1\nmode 2 0 1 0 0 0\nport 0 1\ncoupling_limits 0 16384 0.9 1000 10000 10000 1e-10 1e-11 1e-9\nconnections 0\nevents 0\n";
const DESIGN: &str = "frankensim-equilibrium-design-v1\npreload_limits 0 1 16384\nsensitivity_limits 0 16384 16384 0\ndesign_limits 1 1 1 4\ncases 1\ncase two-loads 2 1\nload 0 0 1\nload 0 0 3\ntarget 0 0 1 1 1\nvariables 1\nvariable force-N 1 1 0 8 1\nbind actuator-force 0 0\n";
fn config(block: usize) -> DesignPlaybackConfig {
    DesignPlaybackConfig { samples:257, force:ForceRenderConfig { sample_rate_hz:48000,
        max_block:block, max_events:8, max_controls:16, max_projection_terms:1024 } }
}
fn events() -> Vec<CaseForceEvent> {
    vec![CaseForceEvent {sample:37,load:0,force_n:0.0}, CaseForceEvent {sample:71,load:1,force_n:0.0}]
}
fn loaded() -> EquilibriumDesignFile {
    EquilibriumDesignFile::from_bytes(MODEL.as_bytes(),DESIGN.as_bytes(),&CancelGate::new()).unwrap()
}
fn actual(x: f64, block: usize) -> Vec<f64> {
    let p=loaded(); let gate=CancelGate::new();
    let playback=p.problem().playback_case(&[x],0,events(),config(block),&mut DesignControl::new(1,1),&gate).unwrap();
    assert_eq!(playback.info().physical_parameters,[1.0+x]);
    let mut r=playback.into_renderer(); let mut out=vec![0.0;257];
    for part in out.chunks_mut(block) {r.block(part).unwrap();}
    out
}

#[test]
fn parameterized_initial_loads_and_independent_releases_match_the_original_stepper() {
    for x in [0.0,2.0] {
        let mut model=ModalAcousticTimeModel::try_new(48000,vec![ModalAcousticMode {
            angular_frequency_rad_s:2.0,damping_ratio:0.0,pressure_per_modal_velocity:C64::new(1.0,0.0),
        }],ModalAcousticTimeBudget {nyquist_guard_fraction:0.9,maximum_abs_displacement_m_sqrt_kg:10.0,
            maximum_abs_velocity_m_sqrt_kg_per_s:1000.0,maximum_total_energy_j:1000.0,maximum_abs_pressure_pa:10000.0}).unwrap();
        model.restore_states(&[ModalAcousticState {displacement_m_sqrt_kg:(4.0+x)/4.0,
            velocity_m_sqrt_kg_per_s:0.0}]).unwrap();
        let expected:Vec<f64>=(0..257).map(|i|model.step(&[if i<37 {4.0+x}else if i<71 {3.0}else{0.0}])
            .unwrap().observer_pressure_pa).collect();
        for block in [1,7,37,64,257] { assert_eq!(actual(x,block),expected); }
        assert!(expected[..37].iter().all(|p|p.abs()<1e-12));
        assert!(expected[71..].iter().any(|p|p.abs()>1e-3));
    }
    assert_ne!(actual(0.0,64),actual(2.0,64));
}

#[test]
fn fitted_contact_parameters_continue_the_same_loaded_network_as_an_authored_performance() {
    let model=include_str!("../examples/equilibrium-design.model").replace("mode 100 0.02 0 0", "mode 100 0.02 1 0");
    let design=include_str!("../examples/equilibrium-design.fit");
    let source=EquilibriumDesignFile::from_bytes(model.as_bytes(),design.as_bytes(),&CancelGate::new()).unwrap();
    let playback=source.problem().playback_case(&[0.5,0.8,0.5],0,
        vec![CaseForceEvent {sample:37,load:0,force_n:0.0}],config(64),
        &mut DesignControl::new(1,3),&CancelGate::new()).unwrap();
    assert_eq!(playback.info().physical_parameters,[600.0,1.8e8,0.0002]);
    assert!(playback.info().initial_energy_j>0.0);
    let authored=model.replace("performance-v3","performance-v4").replace("samples 1\n","samples 257\n")
        .replace("compile_limits 0 2","compile_limits 16 1024")
        .replace("voice free-mass 1 1","voice free-mass-preload 1 1")
        .replace("voice retain-state 1 1","voice static-preload 1 1")
        .replace("port 0 5","port 0.8 5").replace("connection 400 0 0","connection 600 0 0")
        .replace("contact_limits 128", "multi_contact_limits 4 128 16384\ncontacts 1\ncontact_limits 128")
        .replace("contact 100000000 2 0 0.0001", "contact 180000000 2 0 0.0002")
        .replace("events 0", "events 1\nforce 37 0 0 0");
    let mut direct=ModalPerformance::from_bytes(authored.as_bytes(),64).unwrap().into_renderer();
    let mut rendered=playback.into_renderer();
    let (mut a,mut b)=(vec![0.0;257],vec![0.0;257]);
    for (a,b) in a.chunks_mut(64).zip(b.chunks_mut(64)) {direct.block(a).unwrap();rendered.block(b).unwrap();}
    assert_eq!(a,b);
    assert!(a[..37].iter().all(|p|p.abs()<1e-9));
    assert!(a[37..].iter().any(|p|p.abs()>1e-4));
}

#[test]
fn cancellation_at_release_keeps_pending_controls_and_resumes_without_reinitializing() {
    let p=loaded(); let gate=CancelGate::new(); let mut work=DesignControl::new(1,1);
    let playback=p.problem().playback_case(&[2.0],0,events(),config(64),&mut work,&gate).unwrap();
    assert_eq!(work.work(),DesignWork {evaluations:1,case_solves:1});
    let mut r=playback.into_renderer();let mut out=vec![0.0;257];r.block(&mut out[..37]).unwrap();
    let pending=r.pending_controls().to_vec();let cancelled=CancelGate::new();cancelled.request();
    let mut untouched=[9.0;64];
    assert_eq!(r.render_under_gate(&cancelled,&mut untouched,64,1).unwrap(),GatedRenderOutcome::Cancelled {blocks:0});
    assert_eq!(untouched,[9.0;64]);assert_eq!(r.pending_controls(),pending.as_slice());
    for part in out[37..].chunks_mut(64) {r.block(part).unwrap();}
    assert_eq!(out,actual(2.0,64));
}

#[test]
fn bad_requests_refuse_before_preload_and_do_not_mutate_the_design() {
    let p=loaded();let gate=CancelGate::new();let mut work=DesignControl::new(16,16);
    for (case,ev) in [(1,events()),(0,vec![CaseForceEvent {sample:257,load:0,force_n:0.0}]),
        (0,vec![CaseForceEvent {sample:0,load:2,force_n:0.0}]),
        (0,vec![CaseForceEvent {sample:0,load:0,force_n:f64::NAN}])] {
        assert!(p.problem().playback_case(&[0.0],case,ev,config(64),&mut work,&gate).is_err());
        assert_eq!(work.work(),DesignWork::default());
    }
    let mut wrong=config(64);wrong.force.sample_rate_hz=44100;
    assert!(p.problem().playback_case(&[0.0],0,events(),wrong,&mut work,&gate).is_err());
    assert_eq!(work.work().case_solves,0);
    assert!(p.problem().playback_case(&[9.0],0,events(),config(64),&mut work,&gate).is_err());
    assert_eq!(work.work().case_solves,0);
    let cancelled=CancelGate::new();cancelled.request();let before=work.work();
    assert!(p.problem().playback_case(&[0.0],0,events(),config(64),&mut work,&cancelled).is_err());
    assert_eq!(work.work(),before);
    assert_eq!(p.problem().load_cases()[0].loads[0].force_n,1.0);
    let baseline=p.problem().evaluate(&[0.0],&mut DesignControl::new(1,1),&gate).unwrap();
    assert_eq!(baseline.value,0.0);
}
