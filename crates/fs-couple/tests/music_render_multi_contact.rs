//! Public runtime, file decoder and actual command must drive the same contacts.
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::sync::atomic::{AtomicUsize, Ordering};
use fs_blake3::hash_domain;
use fs_couple::modal_acoustic_time::{ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeModel};
use fs_couple::pcm_wav::encode_pcm16_wav;
use fs_couple::render::{ControlDelta, GatedRenderOutcome, RenderError};
use fs_couple::render::schedule::ScheduledRenderer;
use fs_couple::render::schedule::force::coupled::{CoupledModalSystem, ModalAttachment, ModalCouplingConfig};
use fs_couple::render::schedule::force::coupled::contact::{ModalContact, ModalContactConfig};
use fs_couple::render::schedule::force::coupled::contact::multiple::{MultiContactConfig, MultiContactModalSystem};
use fs_couple::render::schedule::force::file::{ModalPerformance, MODAL_MULTI_CONTACT_PERFORMANCE_SCHEMA, MODAL_MULTI_CONTACT_PERFORMANCE_HASH_DOMAIN};
use fs_dcontact::Obstacle;
use fs_exec::CancelGate;
use fs_math::c64::C64;

const INPUT: &str = include_str!("../examples/modal-multi-contact.performance");
const WAV_DOMAIN: &str = "org.frankensim.fs-couple.music-render-wav.v1";
static NEXT: AtomicUsize = AtomicUsize::new(0);

fn direct_system() -> MultiContactModalSystem {
    let models = [(300.0,0.005,0.0,1.0),(800.0,0.01,1.0,0.0),(450.0,0.008,0.0,-0.7)]
        .into_iter().map(|(omega,zeta,gain,velocity)| {
            let mut m = ModalAcousticTimeModel::try_new(48_000, vec![ModalAcousticMode {
                angular_frequency_rad_s: omega, damping_ratio: zeta, pressure_per_modal_velocity: C64::new(gain,0.0),
            }], ModalAcousticTimeBudget { nyquist_guard_fraction:0.9, maximum_abs_displacement_m_sqrt_kg:1.0,
                maximum_abs_velocity_m_sqrt_kg_per_s:1000.0, maximum_total_energy_j:1000.0, maximum_abs_pressure_pa:1000.0 }).unwrap();
            m.restore_states(&[ModalAcousticState { displacement_m_sqrt_kg:0.0, velocity_m_sqrt_kg_per_s:velocity }]).unwrap();
            m
        }).collect();
    let network = CoupledModalSystem::new(models,vec![],ModalCouplingConfig {
        max_modes:4096,max_connections:4,max_setup_terms:4096,nyquist_guard_fraction:0.9,
        maximum_total_energy_j:1000.0,maximum_abs_pressure_pa:1000.0,maximum_abs_connection_force_n:1e6,
        solve_relative_tolerance:1e-11,energy_absolute_tolerance_j:1e-10,energy_relative_tolerance:1e-8,
    },&CancelGate::new()).unwrap();
    let contacts = [(0,1,3e6,0.12,0.0001),(1,2,2.4e6,0.2,0.00012)].into_iter().map(|(a,b,k,chi,gap)| {
        (ModalContact { left:ModalAttachment {component:a,shapes:vec![1.0]},right:ModalAttachment {component:b,shapes:vec![1.0]},
            law:Obstacle::new(vec![-1.0],1,1,vec![gap],vec![1.0],k,1.5,"authored integration comparison".into())
                .unwrap().with_internal_loss(chi).unwrap() },
         ModalContactConfig {max_iterations:96,maximum_force_n:1e6,maximum_penetration_m:0.1,
            force_absolute_tolerance_n:1e-8,force_relative_tolerance:1e-9})
    }).collect();
    MultiContactModalSystem::new(network,contacts,MultiContactConfig {
        max_contacts:4,max_sweeps:64,max_setup_terms:4096,
    },&CancelGate::new()).unwrap()
}
fn renderer(text: &str, block: usize) -> ScheduledRenderer {
    ModalPerformance::from_bytes(text.as_bytes(),block).unwrap().into_renderer()
}
fn pressure(text: &str, block: usize) -> Vec<f64> {
    let mut r=renderer(text,block);let mut out=vec![0.0;1201];
    for chunk in out.chunks_mut(block) { r.block(chunk).unwrap(); }
    out
}
fn directory() -> PathBuf {
    let serial=NEXT.fetch_add(1,Ordering::Relaxed);
    let stamp=std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos();
    let path=std::env::temp_dir().join(format!("fs-multi-contact-{}-{stamp}-{serial}",std::process::id()));
    std::fs::create_dir(&path).unwrap();path
}
fn run(input:&Path,out:&Path,block:usize) -> Output {
    Command::new(env!("CARGO_BIN_EXE_music_render")).arg("modal").arg(input).arg(out)
        .arg("--block").arg(block.to_string()).output().unwrap()
}

#[test]
fn file_schedule_matches_direct_joint_physics_at_all_callback_partitions() {
    let mut direct=direct_system();let mut expected=Vec::new();let mut overlap=false;
    for sample in 0..1201 {
        let frame=direct.step(&[if (300..340).contains(&sample) {1.0}else{0.0},0.0,0.0]).unwrap();
        overlap |= frame.contacts.iter().all(|p|p.normal_force_n>1e-6);
        expected.push(frame.observer_pressure_pa.to_bits());
    }
    assert!(overlap);
    for block in [1,7,37,64,512,1201] {
        assert_eq!(pressure(INPUT,block).iter().map(|p|p.to_bits()).collect::<Vec<_>>(),expected);
    }
    let p=ModalPerformance::from_bytes(INPUT.as_bytes(),37).unwrap();let info=p.info();
    assert_eq!(info.schema,MODAL_MULTI_CONTACT_PERFORMANCE_SCHEMA);
    assert_eq!((info.voices,info.modes,info.contacts,info.connections,info.force_events),(3,3,2,0,2));
    assert_eq!(info.input_hash,hash_domain(MODAL_MULTI_CONTACT_PERFORMANCE_HASH_DOMAIN,INPUT.as_bytes()));
}

#[test]
fn cancelled_callback_retains_joint_state_pending_controls_and_unwritten_output() {
    let expected=pressure(INPUT,64);let mut r=renderer(INPUT,64);let mut out=vec![0.0;1201];
    r.block(&mut out[..64]).unwrap();let pending=r.pending_controls().to_vec();
    let gate=CancelGate::new();gate.request();let mut untouched=[7.0;64];
    assert_eq!(r.render_under_gate(&gate,&mut untouched,64,1).unwrap(),GatedRenderOutcome::Cancelled {blocks:0});
    assert_eq!(untouched,[7.0;64]);assert_eq!(r.samples_rendered(),64);assert_eq!(r.pending_controls(),pending.as_slice());
    for chunk in out[64..].chunks_mut(64) { r.block(chunk).unwrap(); }
    assert_eq!(out,expected);assert!(r.pending_controls().is_empty());
    assert!(r.validate_sample_rate(44_100).is_err());
}

#[test]
fn multi_contact_controls_validate_atomically_and_failed_callbacks_cannot_resume() {
    let make=||renderer(INPUT,32).into_context();
    let mut baseline=make();let mut expected=[0.0;32];baseline.block(&mut expected).unwrap();
    let valid=ControlDelta::SetModalForce {voice:0,mode:0,force_n_per_sqrt_kg:100.0};
    for invalid in [ControlDelta::SetModalForce {voice:0,mode:3,force_n_per_sqrt_kg:1.0},
        ControlDelta::SetBlowingPressure {voice:0,pressure_pa:1.0},ControlDelta::SetPlateForce {voice:0,force_n:1.0}] {
        let mut context=make();assert!(context.apply_controls(&[valid,invalid]).is_err());
        assert!(context.control_log().is_empty());let mut got=[0.0;32];context.block(&mut got).unwrap();assert_eq!(got,expected);
    }
    let mut context=make();context.apply_controls(&[ControlDelta::SetModalForce {voice:0,mode:0,force_n_per_sqrt_kg:1e100}]).unwrap();
    assert!(context.block(&mut [0.0;32]).is_err());let mut untouched=[7.0;32];
    assert!(matches!(context.block(&mut untouched),Err(RenderError::Poisoned)));assert_eq!(untouched,[7.0;32]);
}

#[test]
fn complete_contact_sets_and_budgets_are_required_before_any_render() {
    for text in [INPUT.replace("contacts 2","contacts 0"),INPUT.replace("contacts 2","contacts 33"),
        INPUT.replace("multi_contact_limits 4 64 4096","multi_contact_limits 4 0 4096"),
        INPUT.replace("multi_contact_limits 4 64 4096","multi_contact_limits 4 64 1"),
        INPUT.replace("contact_right 2 1","contact_right 999 1"),
        INPUT.replace("contact_right 2 1","contact_right 2 1 1"),
        INPUT.replace("2400000 1.5 0.2","2400000 1.5 NaN"),
        INPUT.replacen("voice retain-state","voice static-preload",1),format!("{INPUT}ignored\n")] {
        assert!(ModalPerformance::from_bytes(text.as_bytes(),37).is_err(),"{text}");
    }
    for (end,_) in INPUT.match_indices('\n') {
        if end+1<INPUT.len() { assert!(ModalPerformance::from_bytes(&INPUT.as_bytes()[..end+1],37).is_err()); }
    }
}

#[test]
fn actual_command_streams_joint_contact_audio_with_relocatable_identity_and_short_tail() {
    let dir=directory();let input=dir.join("source.performance");std::fs::write(&input,INPUT).unwrap();
    let wave=pressure(INPUT,37);assert!(wave.iter().any(|p|p.abs()>0.05));
    let (expected,clips)=encode_pcm16_wav(&wave,48000,2.0).unwrap();assert_eq!(clips,0);
    let wav_hash=hash_domain(WAV_DOMAIN,&expected).to_hex();
    let input_hash=hash_domain(MODAL_MULTI_CONTACT_PERFORMANCE_HASH_DOMAIN,INPUT.as_bytes()).to_hex();
    for block in [37,512] {
        let out=dir.join(format!("block-{block}.wav"));let result=run(&input,&out,block);
        assert!(result.status.success(),"{}",String::from_utf8_lossy(&result.stdout));
        let bytes=std::fs::read(&out).unwrap();assert_eq!(bytes,expected);assert_eq!(bytes.len(),44+2*1201);
        assert_eq!(u32::from_le_bytes(bytes[40..44].try_into().unwrap()),2*1201);
        let sidecar=std::fs::read_to_string(out.with_extension("provenance.json")).unwrap();
        for token in [format!("\"schema\":\"{MODAL_MULTI_CONTACT_PERFORMANCE_SCHEMA}\""),"\"contacts\":2".to_string(),
            "\"contact\":\"simultaneous-implicit-nonadhesive-power-law\"".to_string(),format!("\"wav_blake3\":\"{wav_hash}\""),
            format!("\"blake3\":\"{input_hash}\"")] { assert!(sidecar.contains(&token),"missing {token}: {sidecar}"); }
    }
    let relocated=dir.join("relocated.performance");std::fs::write(&relocated,INPUT).unwrap();let out=dir.join("replay.wav");
    assert!(run(&relocated,&out,37).status.success());assert_eq!(std::fs::read(&out).unwrap(),expected);
    let sidecar=out.with_extension("provenance.json");let metadata=std::fs::read(&sidecar).unwrap();
    assert_eq!(metadata,std::fs::read(dir.join("block-37.provenance.json")).unwrap());
    assert!(!run(&input,&out,37).status.success());assert_eq!(std::fs::read(&out).unwrap(),expected);
    assert_eq!(std::fs::read(sidecar).unwrap(),metadata);
}

#[test]
fn disabling_all_contacts_silences_the_unforced_receiver_and_bad_sets_create_no_artifacts() {
    let disabled=INPUT.replace("contact 3000000 ","contact 0 ").replace("contact 2400000 ","contact 0 ");
    assert!(pressure(&disabled,37).iter().all(|p|*p==0.0));
    let dir=directory();let input=dir.join("disabled.performance");let out=dir.join("silent.wav");
    std::fs::write(&input,disabled).unwrap();let result=run(&input,&out,37);
    assert!(result.status.success(),"{}",String::from_utf8_lossy(&result.stdout));
    assert!(std::fs::read(&out).unwrap()[44..].iter().all(|b|*b==0));
    std::fs::write(&input,INPUT.replace("contact_right 2 1","contact_right 9 1")).unwrap();
    let out=dir.join("refused.wav");assert!(!run(&input,&out,37).status.success());
    assert!(!out.exists());assert!(!out.with_extension("provenance.json").exists());
}
