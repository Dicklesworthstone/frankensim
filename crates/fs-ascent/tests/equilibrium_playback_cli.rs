//! The actual fitting binary must render its accepted parameters, not its template.
use std::path::{Path,PathBuf};
use std::process::{Command,Output};
use std::sync::atomic::{AtomicUsize,Ordering};
use fs_couple::pcm_wav::encode_pcm16_wav;
use fs_couple::render::schedule::force::ForceRenderConfig;
use fs_couple::render::schedule::force::file::design::EquilibriumDesignFile;
use fs_couple::render::schedule::force::coupled::equilibrium::sensitivity::objective::design::DesignControl;
use fs_couple::render::schedule::force::coupled::equilibrium::sensitivity::objective::design::forward::playback::{CaseForceEvent,DesignPlaybackConfig};
use fs_exec::CancelGate;
const MODEL:&str=include_str!("../../fs-couple/examples/equilibrium-playback.model");
const DESIGN:&str=include_str!("../../fs-couple/examples/equilibrium-playback.fit");
static NEXT:AtomicUsize=AtomicUsize::new(0);
fn files(model:&str)->(PathBuf,PathBuf,PathBuf) {
    let stamp=std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos();
    let dir=std::env::temp_dir().join(format!("fit-playback-{}-{stamp}-{}",std::process::id(),NEXT.fetch_add(1,Ordering::Relaxed)));
    std::fs::create_dir(&dir).unwrap();
    let m=dir.join("model");let d=dir.join("design");let w=dir.join("result.wav");
    std::fs::write(&m,model).unwrap();std::fs::write(&d,DESIGN).unwrap();(m,d,w)
}
fn run(m:&Path,d:&Path,w:&Path,block:usize,extra:&[&str])->Output {
    Command::new(env!("CARGO_BIN_EXE_equilibrium_fit")).arg(m).arg(d)
        .arg("--playback-wav").arg(w).args(["--playback-case","tuned-load","--playback-samples","4801",
            "--playback-release","37","--playback-full-scale-pa","0.004","--playback-block"])
        .arg(block.to_string()).args(extra).output().unwrap()
}
fn success(out:Output)->String {
    assert!(out.status.success(),"{}",String::from_utf8_lossy(&out.stderr));String::from_utf8(out.stdout).unwrap()
}
fn number(s:&str,key:&str)->f64 {s.split_once(key).unwrap().1.split([',','}',']']).next().unwrap().parse().unwrap()}
fn expected(x:f64)->Vec<u8> {
    let gate=CancelGate::new();let p=EquilibriumDesignFile::from_bytes(MODEL.as_bytes(),DESIGN.as_bytes(),&gate).unwrap();
    let mut r=p.problem().playback_case(&[x],0,vec![CaseForceEvent {sample:37,load:0,force_n:0.0}],
        DesignPlaybackConfig {samples:4801,force:ForceRenderConfig {sample_rate_hz:48000,max_block:64,
            max_events:1,max_controls:262144,max_projection_terms:16777216}},
        &mut DesignControl::new(1,1),&gate).unwrap().into_renderer();
    let mut wave=vec![0.0;4801];for chunk in wave.chunks_mut(64) {r.block(chunk).unwrap();}
    let (bytes,clips)=encode_pcm16_wav(&wave,48000,0.004).unwrap();assert_eq!(clips,0);bytes
}
#[test]
fn fitted_parameter_is_rendered_with_the_original_physics_at_multiple_block_sizes() {
    let mut prior=None;
    for block in [37,512] {
        let (m,d,w)=files(MODEL);let text=success(run(&m,&d,&w,block,&[]));
        assert!(text.contains("\"converged\":true"));assert!(number(&text,"\"objective\":")<1e-12);
        let x=number(&text,"\"decision\":");assert!((x-1.0).abs()<1e-6);
        let bytes=std::fs::read(&w).unwrap();assert_eq!(bytes,expected(x));
        assert_eq!(bytes.len(),44+2*4801);
        assert!(bytes[44..44+2*37].iter().all(|x|*x==0));
        assert!(bytes[44+2*37..].iter().any(|x|*x!=0));
        assert_eq!(number(&text,"\"preload_evaluations\":"),1.0);
        assert!(number(&text,"\"evaluations_including_audit\":")<=256.0);
        if let Some(previous)=prior {assert_eq!(bytes,previous);} prior=Some(bytes.clone());
        let refused=run(&m,&d,&w,block,&[]);assert!(!refused.status.success());assert!(refused.stdout.is_empty());
        assert_eq!(std::fs::read(&w).unwrap(),bytes);
        assert_eq!(std::fs::read_to_string(m).unwrap(),MODEL);assert_eq!(std::fs::read_to_string(d).unwrap(),DESIGN);
    }
}
#[test]
fn minimum_budget_keeps_the_initial_design_and_accounts_for_audit_and_playback() {
    let (m,d,w)=files(MODEL);let text=success(run(&m,&d,&w,37,&["--evaluations","3"]));
    assert!(text.contains("\"converged\":false"));assert!(text.contains("\"stop\":\"EvaluationLimit\""));
    assert_eq!(number(&text,"\"decision\":"),0.0);
    assert_eq!(number(&text,"\"evaluations_including_audit\":"),3.0);
    assert_eq!(number(&text,"\"case_solves\":"),3.0);
    assert_eq!(std::fs::read(w).unwrap(),expected(0.0));
}
#[test]
fn no_acoustic_transfer_remains_silent_instead_of_inventing_an_observer_gain() {
    let (m,d,w)=files(&MODEL.replace("mode 1000 0.01 1 0","mode 1000 0.01 0 0"));
    success(run(&m,&d,&w,37,&[]));assert!(std::fs::read(w).unwrap()[44..].iter().all(|b|*b==0));
}
#[test]
fn bad_or_ambiguous_playback_requests_create_no_output_or_partial_result() {
    for extra in [vec!["--evaluations","2"],vec!["--scenarios","unselected-realization"],
        vec!["--playback-case","duplicate"],vec!["--playback-release","4801"]] {
        let (m,d,w)=files(MODEL);let result=run(&m,&d,&w,37,&extra);
        assert!(!result.status.success());assert!(result.stdout.is_empty());assert!(!w.exists());
    }
    let (m,d,w)=files(MODEL);
    let result=Command::new(env!("CARGO_BIN_EXE_equilibrium_fit")).arg(m).arg(d)
        .arg("--playback-wav").arg(&w).output().unwrap();
    assert!(!result.status.success());assert!(!w.exists());
}
