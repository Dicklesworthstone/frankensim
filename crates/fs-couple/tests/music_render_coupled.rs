//! File-to-command regressions for genuine two-way modal coupling.
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::sync::atomic::{AtomicUsize, Ordering};
use fs_blake3::hash_domain;
use fs_couple::pcm_wav::encode_pcm16_wav;
use fs_couple::render::schedule::force::file::{
    ModalPerformance, MODAL_COUPLED_PERFORMANCE_SCHEMA, MODAL_COUPLED_PERFORMANCE_HASH_DOMAIN,
    MODAL_PERFORMANCE_SCHEMA, MODAL_PERFORMANCE_HASH_DOMAIN,
};

const EXAMPLE: &str = include_str!("../examples/coupled-modal.performance");
const WAV_DOMAIN: &str = "org.frankensim.fs-couple.music-render-wav.v1";
static NEXT: AtomicUsize = AtomicUsize::new(0);

fn directory() -> PathBuf {
    let serial = NEXT.fetch_add(1, Ordering::Relaxed);
    let stamp = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos();
    let path = std::env::temp_dir().join(format!("fs-coupled-render-{}-{stamp}-{serial}", std::process::id()));
    std::fs::create_dir(&path).unwrap();
    path
}
fn run(input: &Path, output: &Path, block: usize) -> Output {
    Command::new(env!("CARGO_BIN_EXE_music_render"))
        .arg("modal").arg(input).arg(output).arg("--block").arg(block.to_string())
        .output().unwrap()
}
fn waveform(text: &str, block: usize) -> Vec<f64> {
    let performance = ModalPerformance::from_bytes(text.as_bytes(), block).unwrap();
    let mut pressure = vec![0.0; performance.info().samples as usize];
    let mut renderer = performance.into_renderer();
    for chunk in pressure.chunks_mut(block) { renderer.block(chunk).unwrap(); }
    pressure
}
fn independent_v1() -> String {
    let start = EXAMPLE.find("coupling_limits ").unwrap();
    let end = EXAMPLE.find("events ").unwrap();
    format!("{}{}", &EXAMPLE[..start], &EXAMPLE[end..])
        .replacen(MODAL_COUPLED_PERFORMANCE_SCHEMA, MODAL_PERFORMANCE_SCHEMA, 1)
}

#[test]
fn versioned_connection_records_are_strict_and_v1_keeps_independent_semantics() {
    let p = ModalPerformance::from_bytes(EXAMPLE.as_bytes(), 37).unwrap();
    let info = p.info();
    assert_eq!(info.schema, MODAL_COUPLED_PERFORMANCE_SCHEMA);
    assert_eq!((info.voices, info.modes, info.connections, info.force_events), (2,2,1,3));
    assert_eq!(info.input_hash, hash_domain(MODAL_COUPLED_PERFORMANCE_HASH_DOMAIN, EXAMPLE.as_bytes()));
    assert_ne!(info.input_hash, hash_domain(MODAL_PERFORMANCE_HASH_DOMAIN, EXAMPLE.as_bytes()));
    let v1 = independent_v1();
    let old = ModalPerformance::from_bytes(v1.as_bytes(), 37).unwrap();
    assert_eq!(old.info().schema, MODAL_PERFORMANCE_SCHEMA);
    assert_eq!(old.info().connections, 0);
    assert_eq!(old.info().input_hash, hash_domain(MODAL_PERFORMANCE_HASH_DOMAIN, v1.as_bytes()));
    assert!(waveform(&v1,37).iter().all(|x| *x == 0.0));

    // Removing any complete record suffix must not create a partial network.
    for (end,_) in EXAMPLE.match_indices('\n') {
        if end + 1 < EXAMPLE.len() { assert!(ModalPerformance::from_bytes(&EXAMPLE.as_bytes()[..end+1], 37).is_err()); }
    }
    for text in [
        EXAMPLE.replace("connections 1", "connections 65"),
        EXAMPLE.replace("left 0 1", "left 7 1"),
        EXAMPLE.replace("left 0 1", "left 0 1 2"),
        EXAMPLE.replace("right 1 1", "right 1"),
        EXAMPLE.replace("connection 300000 10 0", "connection -1 10 0"),
        EXAMPLE.replace("connection 300000 10 0", "connection 300000 NaN 0"),
        EXAMPLE.replace("coupling_limits 4 1024", "coupling_limits 4 7"),
        EXAMPLE.replacen("voice static-preload", "voice retain-state", 1),
        format!("{EXAMPLE}ignored\n"),
    ] { assert!(ModalPerformance::from_bytes(text.as_bytes(),37).is_err(), "{text}"); }
}

#[test]
fn actual_command_uses_the_connected_runtime_and_preserves_short_tails_and_replay() {
    let dir = directory();
    let input = dir.join("network.performance");
    std::fs::write(&input, EXAMPLE).unwrap();
    let pressure = waveform(EXAMPLE,37);
    assert!(pressure.iter().any(|p| p.abs()>1e-3), "only the mechanically driven receiver is observed");
    let (expected, clips) = encode_pcm16_wav(&pressure,48_000,0.2).unwrap();
    let expected_hash = hash_domain(WAV_DOMAIN, &expected).to_hex();
    let input_hash = hash_domain(MODAL_COUPLED_PERFORMANCE_HASH_DOMAIN,EXAMPLE.as_bytes()).to_hex();
    for block in [1,37,512] {
        let output = dir.join(format!("block-{block}.wav"));
        let result = run(&input,&output,block);
        assert!(result.status.success(), "{}", String::from_utf8_lossy(&result.stdout));
        let wav = std::fs::read(&output).unwrap();
        assert_eq!(wav,expected);
        assert_eq!(wav.len(),44+2*4801);
        assert_eq!(u32::from_le_bytes(wav[40..44].try_into().unwrap()),2*4801);
        let sidecar = std::fs::read_to_string(output.with_extension("provenance.json")).unwrap();
        for marker in [
            format!("\"schema\":\"{MODAL_COUPLED_PERFORMANCE_SCHEMA}\""),
            "\"connections\":1".to_string(), "\"coupling\":\"implicit-bilateral-spring-damper\"".to_string(),
            format!("\"wav_blake3\":\"{expected_hash}\""), format!("\"blake3\":\"{input_hash}\""),
            format!("\"clipped_samples\":{clips}"),
        ] { assert!(sidecar.contains(&marker), "missing {marker}: {sidecar}"); }
    }
    let relocated = dir.join("relocated.performance");
    std::fs::write(&relocated, EXAMPLE).unwrap();
    let replay = dir.join("replayed.wav");
    assert!(run(&relocated,&replay,37).status.success());
    assert_eq!(std::fs::read(&replay).unwrap(),expected);
    assert_eq!(std::fs::read(replay.with_extension("provenance.json")).unwrap(),
        std::fs::read(dir.join("block-37.provenance.json")).unwrap());
}

#[test]
fn removing_the_mechanical_connection_silences_the_unforced_receiver() {
    let disconnected = EXAMPLE.replace("connection 300000 10 0", "connection 0 0 0");
    let dir = directory();
    let input = dir.join("disconnected.performance");
    let output = dir.join("disconnected.wav");
    std::fs::write(&input,&disconnected).unwrap();
    assert!(waveform(&disconnected,37).iter().all(|x| *x == 0.0));
    let result = run(&input,&output,37);
    assert!(result.status.success(), "{}",String::from_utf8_lossy(&result.stdout));
    let bytes = std::fs::read(output).unwrap();
    assert!(bytes[44..].iter().all(|b| *b == 0));
    assert!(waveform(EXAMPLE,64).iter().any(|x| x.abs()>1e-3));
}

#[test]
fn invalid_connections_and_existing_outputs_refuse_without_replacing_evidence() {
    let dir = directory();
    let input = dir.join("invalid.performance");
    let output = dir.join("must-not-exist.wav");
    std::fs::write(&input,EXAMPLE.replace("right 1 1", "right 999 1")).unwrap();
    assert!(!run(&input,&output,37).status.success());
    assert!(!output.exists());
    assert!(!output.with_extension("provenance.json").exists());
    std::fs::write(&input,EXAMPLE).unwrap();
    std::fs::write(&output,b"existing evidence").unwrap();
    assert!(!run(&input,&output,37).status.success());
    assert_eq!(std::fs::read(&output).unwrap(),b"existing evidence");
    assert!(!output.with_extension("provenance.json").exists());
}
