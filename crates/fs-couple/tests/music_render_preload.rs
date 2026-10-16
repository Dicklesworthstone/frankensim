//! Complete preload -> release -> scheduler -> actual command regression.
use std::path::{Path,PathBuf};
use std::process::{Command,Output};
use std::sync::atomic::{AtomicUsize,Ordering};
use fs_blake3::hash_domain;
use fs_couple::modal_acoustic_time::{ModalAcousticMode,ModalAcousticTimeBudget,ModalAcousticTimeModel};
use fs_couple::pcm_wav::encode_pcm16_wav;
use fs_couple::render::schedule::force::coupled::{CoupledModalSystem,ModalAttachment,ModalCouplingConfig};
use fs_couple::render::schedule::force::coupled::contact::{ModalContact,ModalContactConfig};
use fs_couple::render::schedule::force::coupled::contact::multiple::{MultiContactConfig,MultiContactModalSystem};
use fs_couple::render::schedule::force::file::ModalPerformance;
use fs_dcontact::Obstacle;
use fs_exec::CancelGate;
use fs_math::c64::C64;
const INPUT:&str=include_str!("../examples/modal-contact-preload.performance");
static NEXT:AtomicUsize=AtomicUsize::new(0);
fn directory()->PathBuf {
    let stamp=std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos();
    let dir=std::env::temp_dir().join(format!("fs-contact-preload-{}-{stamp}-{}",std::process::id(),NEXT.fetch_add(1,Ordering::Relaxed)));
    std::fs::create_dir(&dir).unwrap();dir
}
fn run(input:&Path,out:&Path,block:usize)->Output {
    Command::new(env!("CARGO_BIN_EXE_music_render")).arg("modal").arg(input).arg(out)
        .arg("--block").arg(block.to_string()).output().unwrap()
}
fn pressure(text:&str,block:usize)->Vec<f64> {
    let performance=ModalPerformance::from_bytes(text.as_bytes(),block).unwrap();
    let mut out=vec![0.0;performance.info().samples as usize];let mut render=performance.into_renderer();
    for chunk in out.chunks_mut(block){render.block(chunk).unwrap();}out
}
fn direct()->MultiContactModalSystem {
    let modes=[800.0,1000.0,1200.0].into_iter().enumerate().map(|(i,w)|ModalAcousticTimeModel::try_new(48000,
        vec![ModalAcousticMode {angular_frequency_rad_s:w,damping_ratio:0.01,
            pressure_per_modal_velocity:C64::new(if i==1{1.0}else{0.0},0.0)}],ModalAcousticTimeBudget {
            nyquist_guard_fraction:0.9,maximum_abs_displacement_m_sqrt_kg:1.0,maximum_abs_velocity_m_sqrt_kg_per_s:1000.0,
            maximum_total_energy_j:1000.0,maximum_abs_pressure_pa:1000.0}).unwrap()).collect();
    let gate=CancelGate::new();let mut network=CoupledModalSystem::new(modes,vec![],ModalCouplingConfig {
        max_modes:4096,max_connections:4,max_setup_terms:4096,nyquist_guard_fraction:0.9,maximum_total_energy_j:1000.0,
        maximum_abs_pressure_pa:1000.0,maximum_abs_connection_force_n:1e6,solve_relative_tolerance:1e-11,
        energy_absolute_tolerance_j:1e-10,energy_relative_tolerance:1e-8},&gate).unwrap();
    let contacts=[(0,1,1e8,0.1,0.0001),(1,2,8e7,0.2,0.00015)].into_iter().map(|(a,b,k,chi,gap)|(
        ModalContact {left:ModalAttachment{component:a,shapes:vec![1.0]},right:ModalAttachment{component:b,shapes:vec![1.0]},
            law:Obstacle::new(vec![-1.0],1,1,vec![gap],vec![1.0],k,2.0,"authored preload comparison".into())
                .unwrap().with_internal_loss(chi).unwrap()},
        ModalContactConfig{max_iterations:96,maximum_force_n:1e6,maximum_penetration_m:0.1,
            force_absolute_tolerance_n:1e-10,force_relative_tolerance:1e-11})).collect::<Vec<_>>();
    let config=MultiContactConfig {max_contacts:4,max_sweeps:64,max_setup_terms:4096};
    network.initialize_contact_equilibrium(&[400.0,0.0,-300.0],&contacts,config,&gate).unwrap();
    MultiContactModalSystem::new(network,contacts,config,&gate).unwrap()
}

#[test]
fn file_preload_retains_direct_physics_and_is_quiet_until_the_scheduled_release() {
    let mut system=direct();let mut expected=Vec::new();
    for i in 0..4801 {
        let forces=if i<37{[400.0,0.0,-300.0]}else if (300..340).contains(&i){[-60.0,0.0,0.0]}else{[0.0;3]};
        expected.push(system.step(&forces).unwrap().observer_pressure_pa);
    }
    assert!(expected[..37].iter().all(|p|p.abs()<1e-7));
    assert!(expected[37..].iter().any(|p|p.abs()>1e-3));
    for block in [1,37,64,512] {
        let actual=pressure(INPUT,block);assert_eq!(actual,expected,"block {block}");
    }
}

#[test]
fn existing_v3_single_contact_schema_also_accepts_explicit_preload() {
    let input="frankensim-modal-performance-v3\nsample_rate_hz 48000\nsamples 257\nfull_scale_pa 1\nlimits 0.9 1 1000 1000 1000\ncompile_limits 16 64\nvoices 2\nvoice static-preload 1 1\nmode 800 0.01 0 0 0 0\nport 400 1\nvoice static-preload 1 1\nmode 800 0.01 1 0 0 0\nport 0 1\ncoupling_limits 4 4096 0.9 1000 1000 1000000 1e-11 1e-10 1e-8\nconnections 0\ncontact_limits 96 1000000 0.1 1e-10 1e-11\ncontact 100000000 2 0.3 0.0001 1 authored-preload\ncontact_left 0 1\ncontact_right 1 1\nevents 1\nforce 37 0 0 0\n";
    let a=pressure(input,37);let b=pressure(input,64);assert_eq!(a,b);
    assert!(a[..37].iter().all(|x|x.abs()<1e-7));assert!(a[37..].iter().any(|x|x.abs()>1e-3));
}

#[test]
fn actual_command_has_no_startup_impact_keeps_short_tail_and_never_overwrites() {
    let dir=directory();let input=dir.join("preload.performance");std::fs::write(&input,INPUT).unwrap();
    let expected_pressure=pressure(INPUT,37);let(expected,clips)=encode_pcm16_wav(&expected_pressure,48000,0.2).unwrap();
    assert_eq!(clips,0);assert!(expected[44..44+37*2].iter().all(|b|*b==0));
    let hash=hash_domain("org.frankensim.fs-couple.music-render-wav.v1",&expected).to_hex();
    for block in [37,512] {
        let out=dir.join(format!("block-{block}.wav"));let result=run(&input,&out,block);
        assert!(result.status.success(),"{}",String::from_utf8_lossy(&result.stdout));
        assert_eq!(std::fs::read(&out).unwrap(),expected);assert_eq!(expected.len(),44+4801*2);
        let metadata=std::fs::read(out.with_extension("provenance.json")).unwrap();
        assert!(String::from_utf8_lossy(&metadata).contains(&hash));
        assert!(!run(&input,&out,block).status.success());assert_eq!(std::fs::read(&out).unwrap(),expected);
        assert_eq!(std::fs::read(out.with_extension("provenance.json")).unwrap(),metadata);
    }
    let moved=dir.join("same-bytes.performance");std::fs::write(&moved,INPUT).unwrap();let out=dir.join("replay.wav");
    assert!(run(&moved,&out,37).status.success());assert_eq!(std::fs::read(&out).unwrap(),expected);
    assert_eq!(std::fs::read(out.with_extension("provenance.json")).unwrap(),std::fs::read(dir.join("block-37.provenance.json")).unwrap());
}

#[test]
fn mixed_nonzero_or_unfunded_preloads_refuse_before_creating_output() {
    let dir=directory();let input=dir.join("bad.performance");let out=dir.join("absent.wav");
    for bad in [INPUT.replacen("voice static-preload","voice retain-state",1),
        INPUT.replace("mode 800 0.01 0 0 0 0","mode 800 0.01 0 0 0.001 0"),
        INPUT.replace("multi_contact_limits 4 64 4096","multi_contact_limits 4 1 4096"),
        INPUT.replace("contact_limits 96 1000000","contact_limits 96 1"),
        INPUT.replace("multi_contact_limits 4 64 4096","multi_contact_limits 4 64 1")] {
        assert!(ModalPerformance::from_bytes(bad.as_bytes(),37).is_err());
        std::fs::write(&input,bad).unwrap();assert!(!run(&input,&out,37).status.success());
        assert!(!out.exists());assert!(!out.with_extension("provenance.json").exists());
    }
}
