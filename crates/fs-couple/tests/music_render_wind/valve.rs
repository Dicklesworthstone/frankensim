use super::{directory, success, run};
use std::process::Command;
use fs_couple::bernoulli_aperture::performance::file::PlateValvePerformance;
use fs_couple::pcm_wav::{encode_pcm16_wav,decimate::Decimator};
use fs_couple::pcm_wav::observation::PressureRenderer;
use fs_exec::CancelGate;
const INPUT:&str=include_str!("../../examples/plate-valve.performance");
fn reference(text:&str)->(Vec<f64>,usize) {
    let p=PlateValvePerformance::from_bytes(text.as_bytes(),37,&CancelGate::new()).unwrap();
    let info=p.info();let mut source=p.into_renderer();let mut physical=vec![0.0;info.samples as usize];
    for out in physical.chunks_mut(37) {source.block(out).unwrap();}
    let ratio=info.sample_rate_hz as usize/48000;
    let mut filter=Decimator::new(ratio,1).unwrap();let delay=filter.delay_output_frames() as usize;
    let output=physical.chunks_exact(ratio).map(|samples| {
        let p=filter.preview(samples).unwrap()[0];filter.commit();p
    }).collect();(output,delay)
}
fn native()->String {INPUT.replace("audio 96000 9602","audio 48000 4801")}

#[test]
fn supplied_plate_valve_reaches_exact_native_and_decimated_pcm_without_fixture_substitution() {
    let dir=directory();
    for (index,text) in [INPUT.to_string(),native()].iter().enumerate() {
        let input=dir.join(format!("source-{index}.performance"));std::fs::write(&input,text).unwrap();
        let (pressure,delay)=reference(text);
        let (expected,clips)=encode_pcm16_wav(&pressure,48000,20.0).unwrap();
        let hash=fs_blake3::hash_domain("org.frankensim.fs-couple.music-render-wav.v1",&expected).to_hex();
        for block in [1,37,512] {
            let output=dir.join(format!("valve-{index}-{block}.wav"));
            success(run(&input,&output,block,index==0));
            assert_eq!(std::fs::read(&output).unwrap(),expected);
            let metadata=std::fs::read_to_string(output.with_extension("provenance.json")).unwrap();
            for token in ["\"fixture\":\"plate-valve-input\"","frankensim-plate-valve-performance-v1",
                "\"nodes\":21","\"triangles\":24","\"memory_branches\":2","\"pressure_point\":\"inlet\"",
                "not an exterior microphone","\"samples\":4801"] {assert!(metadata.contains(token),"{token}: {metadata}");}
            assert!(metadata.contains(&hash));assert!(metadata.contains(&format!("\"clipped_samples\":{clips}")));
            assert!(metadata.contains(&format!("\"delay_output_samples\":{delay}")));
        }
        let relocated=dir.join(format!("relocated-{index}.data"));std::fs::write(&relocated,text).unwrap();
        let replay=dir.join(format!("replay-{index}.wav"));success(run(&relocated,&replay,37,index==0));
        assert_eq!(std::fs::read(replay.with_extension("provenance.json")).unwrap(),
            std::fs::read(dir.join(format!("valve-{index}-37.provenance.json"))).unwrap());
        assert!(!run(&input,&replay,37,index==0).status.success());assert_eq!(std::fs::read(&replay).unwrap(),expected);
    }
}

#[test]
fn mixed_rate_valve_ensemble_retains_physical_pressure_and_ignores_source_pcm_scales() {
    let dir=directory();let high=dir.join("high.performance");let low=dir.join("low.performance");
    std::fs::write(&low,native()).unwrap();let (a,da)=reference(INPUT);let (b,db)=reference(&native());
    assert!(da>db);let delay=da-db;let mut expected=vec![0.0;4801];
    for (n,out) in expected.iter_mut().enumerate() {*out+=a[n];if n>=delay {*out+=b[n-delay];}}
    let (wav,_)=encode_pcm16_wav(&expected,48000,40.0).unwrap();
    for (index,text) in [INPUT.to_string(),INPUT.replace("9602 20","9602 0.001")].iter().enumerate() {
        std::fs::write(&high,text).unwrap();let output=dir.join(format!("ensemble-{index}.wav"));
        success(Command::new(env!("CARGO_BIN_EXE_music_render")).arg("ensemble").arg(&output)
            .arg("--valve").arg(&high).arg("--valve").arg(&low)
            .args(["--decimate","--block","37","--full-scale-pa","40"]).output().unwrap());
        assert_eq!(std::fs::read(&output).unwrap(),wav);
        let metadata=std::fs::read_to_string(output.with_extension("provenance.json")).unwrap();
        for token in ["\"kind\":\"valve\"","\"source_pcm_scale_applied\":false",
            "internal coupled tube pressure","\"common_delay_output_samples\":40"] {assert!(metadata.contains(token));}
    }
}

#[test]
fn source_physics_and_clock_refusals_precede_output_creation() {
    let dir=directory();let input=dir.join("input.performance");
    for (index,(text,decimate)) in [
        (INPUT.to_string(),false),
        (native(),true),
        (INPUT.replace("audio 96000 9602","audio 100000 10002"),true),
        (INPUT.replace("9602 20","9603 20"),true),
        (INPUT.replace("region 0 4000000000","region 0 3000000000"),true),
        (INPUT.replace("edge 6 13","edge 0 8"),true),
        (INPUT.replace("contact 1000000000000","contact NaN"),true),
        (INPUT.replace("observation inlet","observation exterior"),true),
    ].iter().enumerate() {
        std::fs::write(&input,text).unwrap();let output=dir.join(format!("refused-{index}.wav"));
        assert!(!run(&input,&output,37,*decimate).status.success());
        assert!(!output.exists());assert!(!output.with_extension("provenance.json").exists());
    }
}

const EXTERIOR:&str=include_str!("../../examples/plate-valve-exterior.performance");
#[test]
fn supplied_exterior_receiver_reaches_native_and_decimated_wav_with_physical_delays() {
    let dir=directory();let input=dir.join("baffled.performance");
    for (i,text) in [EXTERIOR.to_string(),EXTERIOR.replace("audio 96000 9602","audio 48000 4801")].iter().enumerate() {
        std::fs::write(&input,text).unwrap();let (pressure,_)=reference(text);
        let (expected,_)=encode_pcm16_wav(&pressure,48000,0.01).unwrap();
        assert!(pressure.iter().any(|p|p.abs()>1e-10));
        assert!(expected[44..].iter().any(|b|*b!=0),"exterior pressure must reach PCM at its declared scale");
        for block in [1,37,512] {
            let output=dir.join(format!("exterior-{i}-{block}.wav"));success(run(&input,&output,block,i==0));
            assert_eq!(std::fs::read(&output).unwrap(),expected);
            let meta=std::fs::read_to_string(output.with_extension("provenance.json")).unwrap();
            for key in ["\"pressure_point\":\"baffled-outlet\"","one-way exterior Rayleigh pressure",
                "\"propagation_delay_mechanical_samples\":[","\"radiation_feedback_added\":false",
                "\"radial_rings\":8","\"angular_points\":32"] {assert!(meta.contains(key),"{key}: {meta}");}
        }
    }
}

#[test]
fn exterior_ensemble_preserves_physical_travel_delays_not_just_filter_alignment() {
    let dir=directory();let high=dir.join("close.performance");let low=dir.join("far.performance");
    let distant=EXTERIOR.replace("audio 96000 9602","audio 48000 4801")
        .replace("baffled-outlet 0 0 0.2","baffled-outlet 0.1 0 0.4");
    std::fs::write(&low,&distant).unwrap();let (a,da)=reference(EXTERIOR);let (b,db)=reference(&distant);
    let delay=da-db;let mut mixed=a.clone();
    for (out,p) in mixed[delay..].iter_mut().zip(&b) {*out+=p;}
    let (expected,_)=encode_pcm16_wav(&mixed,48000,0.02).unwrap();
    for (i,text) in [EXTERIOR.to_string(),EXTERIOR.replace("9602 0.01","9602 100")].iter().enumerate() {
        std::fs::write(&high,text).unwrap();let output=dir.join(format!("receivers-{i}.wav"));
        success(Command::new(env!("CARGO_BIN_EXE_music_render")).arg("ensemble").arg(&output)
            .arg("--valve").arg(&high).arg("--valve").arg(&low)
            .args(["--decimate","--block","37","--full-scale-pa","0.02"]).output().unwrap());
        assert_eq!(std::fs::read(&output).unwrap(),expected);
        let meta=std::fs::read_to_string(output.with_extension("provenance.json")).unwrap();
        assert_eq!(meta.matches("\"radiation_feedback_added\":false").count(),2);
        assert_eq!(meta.matches("\"propagation_delay_mechanical_samples\":").count(),2);
        assert!(!meta.contains("internal coupled tube pressure"));
    }
}

#[test]
fn incomplete_or_underresolved_exterior_observations_refuse_before_wav_creation() {
    let dir=directory();let input=dir.join("receiver.performance");
    for (i,record) in [
        "observation baffled-outlet", "observation baffled-outlet 0 0 0 8 32 4000",
        "observation baffled-outlet 0 0 NaN 8 32 4000", "observation baffled-outlet 0 0 0.2 0 32 4000",
        "observation baffled-outlet 0 0 0.2 8 513 4000", "observation baffled-outlet 0 0 0.2 8 32 9601",
        "observation baffled-outlet 0 0 0.2 1 8 4000", "observation baffled-outlet 0 0 0.2 8 32 4000 ignored",
    ].iter().enumerate() {
        std::fs::write(&input,INPUT.replace("observation inlet",record)).unwrap();
        let output=dir.join(format!("bad-receiver-{i}.wav"));assert!(!run(&input,&output,37,true).status.success());
        assert!(!output.exists());assert!(!output.with_extension("provenance.json").exists());
    }
}

#[path = "radiation.rs"]
mod radiation;
