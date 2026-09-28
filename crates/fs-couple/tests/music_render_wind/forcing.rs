use super::*;
const FORCED:&str=include_str!("../../examples/plate-valve-forced.performance");

#[test]
fn physical_force_file_reaches_exact_pcm_and_mixed_ensembles_without_per_source_gains() {
    let dir=directory();let input=dir.join("force.performance");std::fs::write(&input,FORCED).unwrap();
    let (pressure,delay)=reference(FORCED);let (wav,_)=encode_pcm16_wav(&pressure,48000,1.0).unwrap();
    assert!(wav[44..].iter().any(|v|*v!=0));
    for block in [1,37,512] {
        let output=dir.join(format!("forced-{block}.wav"));success(run(&input,&output,block,true));
        assert_eq!(std::fs::read(&output).unwrap(),wav);
        let meta=std::fs::read_to_string(output.with_extension("provenance.json")).unwrap();
        for token in ["\"mechanical_forces\":", "\"node\":13", "\"triangles\":[10, 11]", "\"events\":7", "not lip-body"] {
            assert!(meta.contains(token),"{token}: {meta}");
        }
    }
    let low=dir.join("native.performance");std::fs::write(&low,native()).unwrap();
    let (b,db)=reference(&native());let mut mix=pressure;
    for (a,b) in mix[(delay-db)..].iter_mut().zip(&b){*a+=b;}
    let (mixed,_)=encode_pcm16_wav(&mix,48000,40.0).unwrap();
    for (i,text) in [FORCED.to_string(),FORCED.replace("9602 1","9602 100")].iter().enumerate() {
        std::fs::write(&input,text).unwrap();let output=dir.join(format!("ensemble-{i}.wav"));
        success(Command::new(env!("CARGO_BIN_EXE_music_render")).arg("ensemble").arg(&output)
            .arg("--valve").arg(&input).arg("--valve").arg(&low)
            .args(["--decimate","--block","37","--full-scale-pa","40"]).output().unwrap());
        assert_eq!(std::fs::read(&output).unwrap(),mixed);
        let meta=std::fs::read_to_string(output.with_extension("provenance.json")).unwrap();
        assert!(meta.contains("\"mechanical_forces\":"));assert!(meta.contains("\"source_pcm_scale_applied\":false"));
    }
}

#[test]
fn invalid_force_sources_refuse_before_creating_waveforms_or_sidecars() {
    let dir=directory();let input=dir.join("force.performance");
    for (i,text) in [FORCED.replace("force_port node 13","force_port node 21"),
        FORCED.replace("force_event 73 0 -0.001","force_event 9602 0 -0.001"),
        FORCED.replace("force_event 73 0 -0.001","force_event 73 0 NaN"),
        FORCED.replace("force_port patch 2 10 11","force_port patch 2 10 10"),
        FORCED.replace("force_ports 2","force_ports 18446744073709551615")].iter().enumerate() {
        std::fs::write(&input,text).unwrap();let output=dir.join(format!("refused-{i}.wav"));
        assert!(!run(&input,&output,37,true).status.success());
        assert!(!output.exists());assert!(!output.with_extension("provenance.json").exists());
    }
}
