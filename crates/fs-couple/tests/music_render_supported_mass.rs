//! Physical supported-mass preload through native force compilation and files.
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::sync::atomic::{AtomicUsize, Ordering};
use fs_blake3::hash_domain;
use fs_couple::modal_acoustic_time::{ModalAcousticMode, ModalAcousticTimeBudget, ModalAcousticTimeModel};
use fs_couple::pcm_wav::encode_pcm16_wav;
use fs_couple::render::GatedRenderOutcome;
use fs_couple::render::schedule::ScheduledRenderer;
use fs_couple::render::schedule::force::{ForceInitialization, ForceRenderConfig, ModalForceEvent, ModalForceVoice};
use fs_couple::render::schedule::force::coupled::{ModalAttachment, ModalConnection, ModalCouplingConfig};
use fs_couple::render::schedule::force::coupled::contact::{ModalContact, ModalContactConfig};
use fs_couple::render::schedule::force::file::{ModalPerformance, MODAL_CONTACT_PERFORMANCE_HASH_DOMAIN};
use fs_dcontact::Obstacle;
use fs_exec::CancelGate;
use fs_math::c64::C64;

const INPUT: &str = include_str!("../examples/supported-mass-preload.performance");
const WAV_DOMAIN: &str = "org.frankensim.fs-couple.music-render-wav.v1";
static NEXT: AtomicUsize = AtomicUsize::new(0);

fn direct() -> ScheduledRenderer {
    let budget=ModalAcousticTimeBudget {nyquist_guard_fraction:0.9,maximum_abs_displacement_m_sqrt_kg:1.0,
        maximum_abs_velocity_m_sqrt_kg_per_s:1000.0,maximum_total_energy_j:1000.0,maximum_abs_pressure_pa:1000.0};
    let mass=ModalAcousticTimeModel::try_free_mass(48000,0.04,0.0,0.0,budget).unwrap();
    let receiver=ModalAcousticTimeModel::try_new(48000,vec![
        ModalAcousticMode {angular_frequency_rad_s:800.0,damping_ratio:0.01,pressure_per_modal_velocity:C64::new(1.0,0.0)},
        ModalAcousticMode {angular_frequency_rad_s:1800.0,damping_ratio:0.02,pressure_per_modal_velocity:C64::new(0.4,0.0)},
    ],budget).unwrap();
    let voices=vec![ModalForceVoice::new(mass,vec![vec![5.0]],vec![2.0],ForceInitialization::StaticPreload).unwrap(),
        ModalForceVoice::new(receiver,vec![vec![1.0,0.5]],vec![0.0],ForceInitialization::StaticPreload).unwrap()];
    let support=ModalConnection {left:ModalAttachment {component:0,shapes:vec![5.0]},
        right:ModalAttachment {component:0,shapes:vec![0.0]},stiffness_n_m:400.0,damping_n_s_m:0.05,rest_extension_m:0.0001};
    let contact=ModalContact {left:ModalAttachment {component:0,shapes:vec![5.0]},
        right:ModalAttachment {component:1,shapes:vec![1.0,0.5]},
        law:Obstacle::new(vec![-1.0],1,1,vec![0.0002],vec![1.0],3e6,1.5,"authored-supported-mass".into())
            .unwrap().with_internal_loss(0.12).unwrap()};
    let events=[(37,0.0),(1201,1.0),(1800,0.0)].into_iter().map(|(sample,force_n)|
        ModalForceEvent {sample,voice:0,port:0,force_n}).collect();
    ScheduledRenderer::from_contact_modal_forces(voices,events,
        ForceRenderConfig {sample_rate_hz:48000,max_block:512,max_events:3,max_controls:16,max_projection_terms:64},
        vec![support],ModalCouplingConfig {max_modes:4096,max_connections:4,max_setup_terms:4096,
            nyquist_guard_fraction:0.9,maximum_total_energy_j:1000.0,maximum_abs_pressure_pa:1000.0,
            maximum_abs_connection_force_n:1e6,solve_relative_tolerance:1e-11,
            energy_absolute_tolerance_j:1e-10,energy_relative_tolerance:1e-8},contact,
        ModalContactConfig {max_iterations:96,maximum_force_n:1e6,maximum_penetration_m:0.1,
            force_absolute_tolerance_n:1e-8,force_relative_tolerance:1e-9},&CancelGate::new()).unwrap()
}
fn pressure(text:&str,block:usize) -> Vec<f64> {
    let p=ModalPerformance::from_bytes(text.as_bytes(),block).unwrap();
    let mut out=vec![0.0;p.info().samples as usize];let mut r=p.into_renderer();
    for chunk in out.chunks_mut(block) {r.block(chunk).unwrap();}
    out
}
fn directory() -> PathBuf {
    let serial=NEXT.fetch_add(1,Ordering::Relaxed);
    let stamp=std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos();
    let path=std::env::temp_dir().join(format!("fs-supported-mass-{}-{stamp}-{serial}",std::process::id()));
    std::fs::create_dir(&path).unwrap();path
}
fn run(input:&Path,out:&Path,block:usize) -> Output {
    Command::new(env!("CARGO_BIN_EXE_music_render")).arg("modal").arg(input).arg(out)
        .arg("--block").arg(block.to_string()).output().unwrap()
}

#[test]
fn supported_mass_file_matches_the_native_physics_and_is_quiet_until_release() {
    let mut r=direct();let mut expected=vec![0.0;4801];
    for chunk in expected.chunks_mut(37) {r.block(chunk).unwrap();}
    assert!(expected[..37].iter().all(|p|p.abs()<1e-8));
    assert!(expected[37..].iter().any(|p|p.abs()>0.001));
    for block in [1,7,37,64,512] {assert_eq!(pressure(INPUT,block),expected);}
    let info=ModalPerformance::from_bytes(INPUT.as_bytes(),37).unwrap().info();
    assert_eq!((info.voices,info.modes,info.connections,info.contacts,info.force_events),(2,3,1,1,3));
    assert_eq!(info.input_hash,hash_domain(MODAL_CONTACT_PERFORMANCE_HASH_DOMAIN,INPUT.as_bytes()));
    let disabled=INPUT.replace("contact 3000000 ","contact 0 ");
    assert!(pressure(&disabled,37).iter().all(|p|*p==0.0));
    let heavier=INPUT.replace("mass 0.04 0 0","mass 0.16 0 0")
        .replace("port 2 5","port 2 2.5").replace("left 0 5","left 0 2.5");
    assert_ne!(pressure(&heavier,37),expected);
}

#[test]
fn cancellation_at_the_release_boundary_keeps_the_settled_state_and_pending_force() {
    let expected=pressure(INPUT,64);let mut out=vec![0.0;4801];
    let mut r=ModalPerformance::from_bytes(INPUT.as_bytes(),64).unwrap().into_renderer();
    r.block(&mut out[..37]).unwrap();let pending=r.pending_controls().to_vec();
    assert_eq!(pending[0].sample,37);
    let gate=CancelGate::new();gate.request();let mut untouched=[7.0;64];
    assert_eq!(r.render_under_gate(&gate,&mut untouched,64,1).unwrap(),GatedRenderOutcome::Cancelled {blocks:0});
    assert_eq!(untouched,[7.0;64]);assert_eq!(r.samples_rendered(),37);assert_eq!(r.pending_controls(),pending.as_slice());
    for chunk in out[37..].chunks_mut(64) {r.block(chunk).unwrap();}
    assert_eq!(out,expected);assert!(r.pending_controls().is_empty());
}

#[test]
fn malformed_or_unsupported_preloads_refuse_before_any_output_is_created() {
    let dir=directory();let input=dir.join("bad.performance");
    let cases=[INPUT.replace("mass 0.04 0 0","mass 0.04 0.001 0"),
        INPUT.replace("mass 0.04 0 0","mass 0.04 0 1"),INPUT.replace("mass 0.04 0 0","mass 0 0 0"),
        INPUT.replace("free-mass-preload 1 1","free-mass-preload 2 1"),
        INPUT.replace("voice static-preload","voice retain-state"),
        INPUT.replace("connection 400 0.05 0.0001","connection 0 0.05 0.0001"),
        INPUT.replace("coupling_limits 4 4096","coupling_limits 4 14")];
    for (i,text) in cases.iter().enumerate() {
        assert!(ModalPerformance::from_bytes(text.as_bytes(),37).is_err(),"{text}");
        std::fs::write(&input,text).unwrap();let out=dir.join(format!("refused-{i}.wav"));
        assert!(!run(&input,&out,37).status.success());
        assert!(!out.exists());assert!(!out.with_extension("provenance.json").exists());
    }
    // Ordinary retained free-mass input still carries supplied state, not preload.
    let retained=INPUT.replace("free-mass-preload","free-mass").replace("voice static-preload","voice retain-state");
    assert!(ModalPerformance::from_bytes(retained.as_bytes(),37).is_ok());
    let cut=retained.find("coupling_limits ").unwrap();
    let end=retained.find("events ").unwrap();
    let independent=format!("{}{}",&retained[..cut],&retained[end..])
        .replace("frankensim-modal-performance-v3","frankensim-modal-performance-v1");
    assert!(ModalPerformance::from_bytes(independent.as_bytes(),37).is_ok());
    assert!(ModalPerformance::from_bytes(independent.replace("voice free-mass ","voice free-mass-preload ").as_bytes(),37).is_err());
}

#[test]
fn actual_command_preserves_waveform_hash_short_tail_relocation_and_existing_files() {
    let dir=directory();let input=dir.join("supported.performance");std::fs::write(&input,INPUT).unwrap();
    let (expected,clips)=encode_pcm16_wav(&pressure(INPUT,37),48000,0.1).unwrap();assert_eq!(clips,0);
    let wav_hash=hash_domain(WAV_DOMAIN,&expected).to_hex();
    let input_hash=hash_domain(MODAL_CONTACT_PERFORMANCE_HASH_DOMAIN,INPUT.as_bytes()).to_hex();
    for block in [37,512] {
        let out=dir.join(format!("block-{block}.wav"));let result=run(&input,&out,block);
        assert!(result.status.success(),"{}",String::from_utf8_lossy(&result.stdout));
        let bytes=std::fs::read(&out).unwrap();assert_eq!(bytes,expected);assert_eq!(bytes.len(),44+2*4801);
        assert_eq!(u32::from_le_bytes(bytes[40..44].try_into().unwrap()),2*4801);
        let sidecar=std::fs::read_to_string(out.with_extension("provenance.json")).unwrap();
        assert!(sidecar.contains(&format!("\"wav_blake3\":\"{wav_hash}\"")));
        assert!(sidecar.contains(&format!("\"blake3\":\"{input_hash}\"")));
        assert!(sidecar.contains("\"contacts\":1"));
    }
    let moved=dir.join("relocated.performance");std::fs::write(&moved,INPUT).unwrap();let out=dir.join("replay.wav");
    assert!(run(&moved,&out,37).status.success());assert_eq!(std::fs::read(&out).unwrap(),expected);
    let metadata=std::fs::read(out.with_extension("provenance.json")).unwrap();
    assert_eq!(metadata,std::fs::read(dir.join("block-37.provenance.json")).unwrap());
    assert!(!run(&input,&out,37).status.success());assert_eq!(std::fs::read(&out).unwrap(),expected);
    assert_eq!(metadata,std::fs::read(out.with_extension("provenance.json")).unwrap());
}
