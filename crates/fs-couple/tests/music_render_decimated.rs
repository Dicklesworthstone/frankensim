use std::path::{Path,PathBuf};
use std::process::{Command,Output};
use std::sync::atomic::{AtomicUsize,Ordering};
use fs_blake3::hash_domain;
use fs_couple::pcm_wav::{decimate::Decimator,encode_pcm16_wav};
use fs_couple::render::schedule::force::file::ModalPerformance;

const INPUT:&str=include_str!("../examples/free-striker-192k.performance");
const LEGACY:&str=include_str!("../examples/free-striker.performance");
const WAV_DOMAIN:&str="org.frankensim.fs-couple.music-render-wav.v1";
static NEXT:AtomicUsize=AtomicUsize::new(0);
fn directory()->PathBuf {
    let tick=std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos();
    let dir=std::env::temp_dir().join(format!("fs-decimated-{}-{tick}-{}",std::process::id(),NEXT.fetch_add(1,Ordering::Relaxed)));
    std::fs::create_dir(&dir).unwrap();dir
}
fn run(input:&Path,output:&Path,block:usize,decimate:bool)->Output {
    let mut command=Command::new(env!("CARGO_BIN_EXE_music_render"));
    command.arg("modal").arg(input).arg(output).arg("--block").arg(block.to_string());
    if decimate {command.arg("--decimate");}
    command.output().unwrap()
}
fn raw_pressure(text:&str)->Vec<f64> {
    let p=ModalPerformance::from_bytes(text.as_bytes(),37).unwrap();
    let mut out=vec![0.0;p.info().samples as usize];let mut r=p.into_renderer();
    for block in out.chunks_mut(37) {r.block(block).unwrap();}
    out
}
fn observed(text:&str,ratio:usize)->Vec<f64> {
    let raw=raw_pressure(text);let mut filter=Decimator::new(ratio,1).unwrap();
    assert_eq!(raw.len()%ratio,0);
    raw.chunks_exact(ratio).map(|group| {let v=filter.preview(group).unwrap()[0];filter.commit();v}).collect()
}

#[test]
fn actual_striker_command_matches_high_rate_mechanics_then_shared_filter_and_encoder() {
    let dir=directory();let input=dir.join("high-rate.performance");std::fs::write(&input,INPUT).unwrap();
    let pressure=observed(INPUT,4);assert_eq!(pressure.len(),4801);
    assert!(pressure.iter().any(|v|v.abs()>0.01));
    let (expected,clips)=encode_pcm16_wav(&pressure,48_000,1.0).unwrap();assert_eq!(clips,0);
    let hash=hash_domain(WAV_DOMAIN,&expected).to_hex();
    let info=ModalPerformance::from_bytes(INPUT.as_bytes(),37).unwrap().info();
    for block in [1,37,512] {
        let out=dir.join(format!("block-{block}.wav"));let result=run(&input,&out,block,true);
        assert!(result.status.success(),"{}",String::from_utf8_lossy(&result.stdout));
        let wav=std::fs::read(&out).unwrap();assert_eq!(wav,expected);assert_eq!(wav.len(),44+4801*2);
        assert_eq!(u32::from_le_bytes(wav[24..28].try_into().unwrap()),48_000);
        assert_eq!(u32::from_le_bytes(wav[40..44].try_into().unwrap()),4801*2);
        let sidecar=std::fs::read_to_string(out.with_extension("provenance.json")).unwrap();
        for token in ["\"mechanics_sample_rate_hz\":192000".to_string(),"\"mechanics_samples\":19204".to_string(),
            "\"samples\":4801".to_string(),"\"ratio\":4".to_string(),"\"delay_output_samples\":44".to_string(),
            "\"first_output_source_index\":3".to_string(),"\"delay_compensated\":false".to_string(),
            "\"tail\":\"no-flush-declared-window\"".to_string(),format!("\"wav_blake3\":\"{hash}\""),
            format!("\"blake3\":\"{}\"",info.input_hash.to_hex())] {
            assert!(sidecar.contains(&token),"missing {token}: {sidecar}");
        }
    }
    let relocated=dir.join("relocated.performance");std::fs::write(&relocated,INPUT).unwrap();
    let out=dir.join("replay.wav");assert!(run(&relocated,&out,37,true).status.success());
    assert_eq!(std::fs::read(&out).unwrap(),expected);
    let sidecar=out.with_extension("provenance.json");let metadata=std::fs::read(&sidecar).unwrap();
    assert_eq!(metadata,std::fs::read(dir.join("block-37.provenance.json")).unwrap());
    assert!(!run(&input,&out,37,true).status.success());assert_eq!(std::fs::read(&out).unwrap(),expected);
    assert_eq!(std::fs::read(&sidecar).unwrap(),metadata);
}

#[test]
fn explicit_rate_choice_and_complete_intervals_are_required_before_artifact_creation() {
    let dir=directory();let input=dir.join("input.performance");
    for (index,(text,flag)) in [
        (INPUT.to_string(),false),
        (INPUT.replace("samples 19204","samples 19205"),true),
        (INPUT.replace("sample_rate_hz 192000","sample_rate_hz 100000"),true),
        (INPUT.replace("sample_rate_hz 192000","sample_rate_hz 44100"),true),
        (LEGACY.to_string(),true),
    ].into_iter().enumerate() {
        std::fs::write(&input,text).unwrap();let out=dir.join(format!("refused-{index}.wav"));
        assert!(!run(&input,&out,37,flag).status.success());assert!(!out.exists());
        assert!(!out.with_extension("provenance.json").exists());
    }
}

#[test]
fn ordinary_48khz_path_keeps_its_original_samples_and_unfiltered_provenance() {
    let dir=directory();let input=dir.join("legacy.performance");let out=dir.join("legacy.wav");
    std::fs::write(&input,LEGACY).unwrap();
    let pressure=raw_pressure(LEGACY);let (expected,_)=encode_pcm16_wav(&pressure,48_000,1.0).unwrap();
    let result=run(&input,&out,37,false);
    assert!(result.status.success(),"{}",String::from_utf8_lossy(&result.stdout));
    assert_eq!(std::fs::read(&out).unwrap(),expected);
    let sidecar=std::fs::read_to_string(out.with_extension("provenance.json")).unwrap();
    assert!(!sidecar.contains("\"observation\""));assert!(!sidecar.contains("delay_output_samples"));
    assert!(sidecar.contains("\"samples\":4801"));
    assert!(sidecar.contains(&format!("\"wav_blake3\":\"{}\"",hash_domain(WAV_DOMAIN,&expected).to_hex())));
}

#[test]
fn high_frequency_source_is_not_merely_relabelled_or_unfiltered_subsampled_by_the_command() {
    let dir=directory();let input=dir.join("tone.performance");
    let mut levels=Vec::new();
    for frequency in [6000.0,36000.0] {
        let text=format!("frankensim-modal-performance-v1\nsample_rate_hz 96000\nsamples 3072\nfull_scale_pa 2\n\
            limits 0.9 1 1000 1000 1000\ncompile_limits 0 1\nvoices 1\nvoice retain-state 1 1\n\
            mode {} 0 1 0 0 1\nport 0 1\nevents 0\n",core::f64::consts::TAU*frequency);
        std::fs::write(&input,&text).unwrap();let out=dir.join(format!("tone-{frequency}.wav"));
        let result=run(&input,&out,37,true);assert!(result.status.success(),"{}",String::from_utf8_lossy(&result.stdout));
        let wav=std::fs::read(&out).unwrap();assert_eq!(wav.len(),44+2*1536);
        let pressure:Vec<_>=wav[44..].chunks_exact(2).map(|b|f64::from(i16::from_le_bytes(b.try_into().unwrap()))*2.0/32767.0).collect();
        levels.push((pressure[512..].iter().map(|v|v*v).sum::<f64>()/1024.0).sqrt());
    }
    assert!((levels[0]-std::f64::consts::FRAC_1_SQRT_2).abs()<3e-4,"{levels:?}");
    assert!(levels[1]<1e-5,"{levels:?}");
}
