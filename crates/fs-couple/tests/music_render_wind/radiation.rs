//! Actual command-line consumer of the geometry-derived terminal feedback.
use super::*;
const RADIATING:&str=include_str!("../../examples/plate-valve-radiating.performance");

#[test]
fn radiation_loaded_native_and_decimated_wav_preserve_source_clocks_and_explicit_scale() {
    let dir=directory();let input=dir.join("radiating.performance");
    for (i,text) in [RADIATING.to_string(),RADIATING.replace("audio 96000 9602","audio 48000 4801")].iter().enumerate() {
        std::fs::write(&input,text).unwrap();let (pressure,_)=reference(text);
        let (expected,clips)=encode_pcm16_wav(&pressure,48000,0.01).unwrap();
        assert!(expected[44..].iter().any(|b|*b!=0));
        for block in [1,37,512] {
            let output=dir.join(format!("radiating-{i}-{block}.wav"));success(run(&input,&output,block,i==0));
            assert_eq!(std::fs::read(&output).unwrap(),expected);
            let metadata=std::fs::read_to_string(output.with_extension("provenance.json")).unwrap();
            for key in ["\"pressure_point\":\"baffled-outlet\"","\"radiation_feedback_added\":true",
                "compact-baffled-piston-positive-real-v1","\"extra_end_correction_m\":0",
                "\"replaces_memoryless_reflection\":true","\"propagation_delay_mechanical_samples\":[",
                "checked_max_complex_relative_error","checked_max_resistance_relative_error"] {
                assert!(metadata.contains(key),"{key}: {metadata}");
            }
            assert!(!metadata.contains("one-way exterior Rayleigh pressure"));
            assert!(metadata.contains(&format!("\"clipped_samples\":{clips}")));
        }
        let moved=dir.join(format!("relocated-{i}.data"));std::fs::write(&moved,text).unwrap();
        let output=dir.join(format!("relocated-{i}.wav"));success(run(&moved,&output,37,i==0));
        assert_eq!(std::fs::read(output.with_extension("provenance.json")).unwrap(),
            std::fs::read(dir.join(format!("radiating-{i}-37.provenance.json"))).unwrap());
        assert!(!run(&moved,&output,37,i==0).status.success());assert_eq!(std::fs::read(&output).unwrap(),expected);
    }
}

#[test]
fn loaded_ensemble_uses_accepted_flow_and_preserves_propagation_without_source_gains() {
    let dir=directory();let high=dir.join("near.performance");let low=dir.join("far.performance");
    let far=RADIATING.replace("audio 96000 9602","audio 48000 4801")
        .replace("baffled-outlet 0 0 0.2","baffled-outlet 0.1 0 0.4");
    std::fs::write(&low,&far).unwrap();let (a,da)=reference(RADIATING);let (b,db)=reference(&far);
    let delay=da-db;let mut mixed=a.clone();
    for (sample,other) in mixed[delay..].iter_mut().zip(&b) {*sample+=other;}
    let expected=encode_pcm16_wav(&mixed,48000,0.02).unwrap().0;
    for (i,text) in [RADIATING.to_string(),RADIATING.replace("9602 0.01","9602 100")].iter().enumerate() {
        std::fs::write(&high,text).unwrap();let output=dir.join(format!("loaded-mix-{i}.wav"));
        success(Command::new(env!("CARGO_BIN_EXE_music_render")).arg("ensemble").arg(&output)
            .arg("--valve").arg(&high).arg("--valve").arg(&low)
            .args(["--decimate","--block","37","--full-scale-pa","0.02"]).output().unwrap());
        assert_eq!(std::fs::read(&output).unwrap(),expected);
        let metadata=std::fs::read_to_string(output.with_extension("provenance.json")).unwrap();
        assert_eq!(metadata.matches("\"radiation_feedback_added\":true").count(),2);
        assert_eq!(metadata.matches("\"source_pcm_scale_applied\":false").count(),2);
        assert_eq!(metadata.matches("\"propagation_delay_mechanical_samples\":").count(),2);
        assert_eq!(metadata.matches("\"radiation_load\":").count(),2);
    }
}

#[test]
fn invalid_radiation_cannot_create_wav_or_silently_fall_back_to_fixed_reflection() {
    let dir=directory();let input=dir.join("radiation.performance");
    for (i,source) in [
        RADIATING.replace("baffled-low-ka 1500","baffled-low-ka 1500 -0.8"),
        RADIATING.replace("baffled-low-ka 1500","baffled-low-ka 0"),
        RADIATING.replace("baffled-low-ka 1500","baffled-low-ka 10000"),
        RADIATING.replace("baffled-low-ka 1500","baffled-low-ka NaN"),
        RADIATING.replace("8 32 1500","8 32 2000"),
        RADIATING.replace("audio 96000 9602","audio 96000 9603"),
    ].iter().enumerate() {
        std::fs::write(&input,source).unwrap();let output=dir.join(format!("refused-radiation-{i}.wav"));
        assert!(!run(&input,&output,37,true).status.success());
        assert!(!output.exists());assert!(!output.with_extension("provenance.json").exists());
    }
    // Internal pressure is still a useful diagnostic of the same loaded source.
    let source=RADIATING.replace("observation baffled-outlet 0 0 0.2 8 32 1500","observation terminal");
    std::fs::write(&input,&source).unwrap();let output=dir.join("loaded-internal.wav");
    success(run(&input,&output,37,true));
    let (pressure,_)=reference(&source);
    assert_eq!(std::fs::read(&output).unwrap(),encode_pcm16_wav(&pressure,48000,0.01).unwrap().0);
    let metadata=std::fs::read_to_string(output.with_extension("provenance.json")).unwrap();
    assert!(metadata.contains("\"pressure_point\":\"terminal\""));
    assert!(metadata.contains("\"radiation_load\":"));assert!(!metadata.contains("outlet_receiver"));
    assert!(metadata.contains("not an exterior microphone"));
}
