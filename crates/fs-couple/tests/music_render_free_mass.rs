//! Physical mass declarations feed the existing contact/scheduler/WAV path.
use std::path::Path;
use std::process::{Command, Output};
use fs_blake3::hash_domain;
use fs_couple::modal_acoustic_time::{ModalAcousticMode, ModalAcousticTimeBudget, ModalAcousticTimeModel};
use fs_couple::pcm_wav::encode_pcm16_wav;
use fs_couple::render::GatedRenderOutcome;
use fs_couple::render::schedule::force::file::{ModalPerformance, MODAL_CONTACT_PERFORMANCE_HASH_DOMAIN};
use fs_couple::render::schedule::force::coupled::{CoupledModalSystem, ModalAttachment, ModalCouplingConfig};
use fs_couple::render::schedule::force::coupled::contact::{ContactModalSystem, ModalContact, ModalContactConfig};
use fs_dcontact::Obstacle;
use fs_exec::CancelGate;
use fs_math::c64::C64;

const INPUT: &str = include_str!("../examples/free-striker.performance");
const WAV_DOMAIN: &str = "org.frankensim.fs-couple.music-render-wav.v1";

fn direct() -> ContactModalSystem {
    let budget = ModalAcousticTimeBudget { nyquist_guard_fraction:0.9,
        maximum_abs_displacement_m_sqrt_kg:1.0, maximum_abs_velocity_m_sqrt_kg_per_s:1000.0,
        maximum_total_energy_j:1000.0, maximum_abs_pressure_pa:1000.0 };
    let striker = ModalAcousticTimeModel::try_free_mass(48_000,0.04,0.0,1.0,budget).unwrap();
    let receiver = ModalAcousticTimeModel::try_new(48_000,vec![
        ModalAcousticMode {angular_frequency_rad_s:800.0,damping_ratio:0.01,pressure_per_modal_velocity:C64::new(1.0,0.0)},
        ModalAcousticMode {angular_frequency_rad_s:1800.0,damping_ratio:0.02,pressure_per_modal_velocity:C64::new(0.4,0.0)},
    ],budget).unwrap();
    let network = CoupledModalSystem::new(vec![striker,receiver],vec![],ModalCouplingConfig {
        max_modes:4096,max_connections:4,max_setup_terms:4096,nyquist_guard_fraction:0.9,
        maximum_total_energy_j:1000.0,maximum_abs_pressure_pa:1000.0,maximum_abs_connection_force_n:1e6,
        solve_relative_tolerance:1e-11,energy_absolute_tolerance_j:1e-10,energy_relative_tolerance:1e-8,
    },&CancelGate::new()).unwrap();
    ContactModalSystem::new(network,ModalContact {
        left:ModalAttachment {component:0,shapes:vec![5.0]},
        right:ModalAttachment {component:1,shapes:vec![1.0,0.5]},
        law:Obstacle::new(vec![-1.0],1,1,vec![0.0005],vec![1.0],3e6,1.5,
            "authored-free-striker".into()).unwrap().with_internal_loss(0.12).unwrap(),
    },ModalContactConfig {max_iterations:96,maximum_force_n:1e6,maximum_penetration_m:0.1,
        force_absolute_tolerance_n:1e-8,force_relative_tolerance:1e-9},&CancelGate::new()).unwrap()
}
fn waveform(text:&str,block:usize) -> Vec<f64> {
    let performance=ModalPerformance::from_bytes(text.as_bytes(),block).unwrap();
    let mut out=vec![0.0;performance.info().samples as usize];let mut r=performance.into_renderer();
    for chunk in out.chunks_mut(block) {r.block(chunk).unwrap();}
    out
}
fn command(input:&Path,out:&Path,block:usize) -> Output {
    Command::new(env!("CARGO_BIN_EXE_music_render")).arg("modal").arg(input).arg(out)
        .arg("--block").arg(block.to_string()).output().unwrap()
}

#[test]
fn mass_file_matches_physical_contact_api_and_rebound_is_not_spring_driven() {
    let mut system=direct();let mut expected=Vec::new();let mut engaged=false;let mut coast=None;
    for sample in 0..4801 {
        let f=system.step(&[if (300..340).contains(&sample) {5.0}else{0.0},0.0,0.0]).unwrap();
        engaged |= f.normal_force_n>1e-6;
        expected.push(f.observer_pressure_pa);
        if sample>=512 {
            assert_eq!(f.normal_force_n,0.0);
            let v=system.components()[0].states()[0].velocity_m_sqrt_kg_per_s;
            assert!(v<0.0,"striker must rebound");
            if let Some(bits)=coast {assert_eq!(v.to_bits(),bits);} else {coast=Some(v.to_bits());}
        }
    }
    assert!(engaged);assert!(expected[..20].iter().all(|p|*p==0.0));
    assert!(expected.iter().any(|p|p.abs()>0.01));
    for block in [1,7,37,64,512,4801] {assert_eq!(waveform(INPUT,block),expected);}
    let p=ModalPerformance::from_bytes(INPUT.as_bytes(),64).unwrap();
    assert_eq!((p.info().voices,p.info().modes,p.info().contacts),(2,3,1));
    assert_eq!(p.info().input_hash,hash_domain(MODAL_CONTACT_PERFORMANCE_HASH_DOMAIN,INPUT.as_bytes()));
    let mut r=p.into_renderer();let mut out=vec![0.0;4801];r.block(&mut out[..64]).unwrap();
    let gate=CancelGate::new();gate.request();let mut untouched=[7.0;64];
    assert_eq!(r.render_under_gate(&gate,&mut untouched,64,1).unwrap(),GatedRenderOutcome::Cancelled {blocks:0});
    assert_eq!(untouched,[7.0;64]);assert_eq!(r.samples_rendered(),64);
    for chunk in out[64..].chunks_mut(64) {r.block(chunk).unwrap();}
    assert_eq!(out,expected);
}

#[test]
fn free_mass_input_is_explicit_complete_and_cannot_request_undefined_preload() {
    for text in [
        INPUT.replace("voice free-mass 1 1","voice free-mass 2 1"),
        INPUT.replace("voice free-mass 1 1","voice free-mass 1 0"),
        INPUT.replace("voice free-mass 1 1","voice free-mass 1 65537"),
        INPUT.replace("mass 0.04 0 1","mass -0.04 0 1"),
        INPUT.replace("mass 0.04 0 1","mass 0 0 1"),
        INPUT.replace("mass 0.04 0 1","mass NaN 0 1"),
        INPUT.replace("mass 0.04 0 1","mass 0.04 0 inf"),
        INPUT.replace("mass 0.04 0 1","mass 0.04 0"),
        INPUT.replace("mass 0.04 0 1","mass 0.04 0 1 99"),
        INPUT.replace("mass 0.04 0 1","mode 0 0 0 0 0 0"),
        INPUT.replace("mode 800 0.01","mode 0 0.01"),
        INPUT.replace("voice retain-state 2 1","voice static-preload 2 1"),
    ] {assert!(ModalPerformance::from_bytes(text.as_bytes(),37).is_err(),"{text}");}
    for (end,_) in INPUT.match_indices('\n') {
        if end+1<INPUT.len() {assert!(ModalPerformance::from_bytes(&INPUT.as_bytes()[..end+1],37).is_err());}
    }
    // A free mass is legal without any contact too, but its motion alone is
    // not reinterpreted as a narrow-band audio signal.
    let drift="frankensim-modal-performance-v1\nsample_rate_hz 48000\nsamples 65\nfull_scale_pa 1\nlimits 0.9 1 1000 1000 1000\ncompile_limits 2 2\nvoices 1\nvoice free-mass 1 1\nmass 0.04 0 1\nport 1 5\nevents 1\nforce 37 0 0 0\n";
    assert!(waveform(drift,37).iter().all(|p|*p==0.0));
}

#[test]
fn receiver_sound_depends_on_collision_mass_and_launch_velocity() {
    let base=waveform(INPUT,37);
    for text in [INPUT.replace("contact 3000000 ","contact 0 "),
        INPUT.replace("0.0005 1 authored-free-striker","1 1 authored-free-striker")] {
        assert!(waveform(&text,37).iter().all(|p|*p==0.0));
    }
    let slow=INPUT.replace("mass 0.04 0 1","mass 0.04 0 0.5");
    // A physical mass change must also change its conjugate attachment/port
    // normalization. Shapes are never silently reinterpreted by the reader.
    let heavy=INPUT.replace("mass 0.04 0 1","mass 0.16 0 1")
        .replace("port 0 5","port 0 2.5").replace("contact_left 0 5","contact_left 0 2.5");
    for text in [slow,heavy] {
        let changed=waveform(&text,37);
        assert!(base.iter().zip(changed).any(|(a,b)|(a-b).abs()>0.001));
    }
}

#[test]
fn actual_command_streams_the_strike_preserves_tail_and_refuses_invalid_or_existing_output() {
    let stamp=std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos();
    let dir=std::env::temp_dir().join(format!("fs-free-striker-{}-{stamp}",std::process::id()));
    std::fs::create_dir(&dir).unwrap();let input=dir.join("striker.performance");
    std::fs::write(&input,INPUT).unwrap();
    let (expected,clipped)=encode_pcm16_wav(&waveform(INPUT,37),48_000,1.0).unwrap();
    assert_eq!(clipped,0);let wav_hash=hash_domain(WAV_DOMAIN,&expected).to_hex();
    for block in [37,512] {
        let out=dir.join(format!("strike-{block}.wav"));let result=command(&input,&out,block);
        assert!(result.status.success(),"stdout={} stderr={}",String::from_utf8_lossy(&result.stdout),String::from_utf8_lossy(&result.stderr));
        let bytes=std::fs::read(&out).unwrap();assert_eq!(bytes,expected);assert_eq!(bytes.len(),44+2*4801);
        let metadata=std::fs::read_to_string(out.with_extension("provenance.json")).unwrap();
        assert!(metadata.contains(&format!("\"wav_blake3\":\"{wav_hash}\"")));
        assert!(!command(&input,&out,block).status.success());assert_eq!(std::fs::read(&out).unwrap(),bytes);
        assert_eq!(std::fs::read_to_string(out.with_extension("provenance.json")).unwrap(),metadata);
    }
    std::fs::write(&input,INPUT.replace("mass 0.04 0 1","mass -1 0 1")).unwrap();
    let out=dir.join("invalid.wav");assert!(!command(&input,&out,37).status.success());
    assert!(!out.exists());assert!(!out.with_extension("provenance.json").exists());
}
