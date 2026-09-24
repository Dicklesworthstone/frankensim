use super::*;
use fs_couple::pcm_wav::observation::{DecimatedRenderer,PressureRenderer};
use fs_couple::pcm_wav::encode_pcm16_wav;
const LOSSY:&str=include_str!("../../examples/plate-valve-viscothermal.performance");
fn expected(text:&str)->Vec<u8> {
    let p=PlateValvePerformance::from_bytes(text.as_bytes(),37,&CancelGate::new()).unwrap();let i=p.info();
    let mut r=DecimatedRenderer::new(p.into_renderer(),i.sample_rate_hz,48000,37).unwrap();
    let samples=r.output_samples_for(i.samples).unwrap() as usize;let mut out=vec![0.0;samples];
    for block in out.chunks_mut(37) {r.block(block).unwrap();}
    encode_pcm16_wav(&out,48000,i.full_scale_pa).unwrap().0
}

#[test]
fn wind_cli_keeps_gas_losses_through_callback_partitions_and_explicit_decimation() {
    let dir=directory();let input=dir.join("lossy.performance");std::fs::write(&input,LOSSY).unwrap();
    let wav=expected(LOSSY);assert!(wav[44..].iter().any(|&v|v!=0));
    let lossless=LOSSY.replace(" viscothermal 100 1000 12 8","");
    assert_ne!(wav[44..],expected(&lossless)[44..],"gas-wall physics must reach PCM");
    for block in [1,37,512] {
        let output=dir.join(format!("lossy-{block}.wav"));success(run(&input,&output,block,true));
        assert_eq!(std::fs::read(&output).unwrap(),wav);
        let sidecar=std::fs::read_to_string(output.with_extension("provenance.json")).unwrap();
        for token in ["wide-tube-zk-passive-rl-rc-v1","\"cells\":12","\"arms_per_load\":8",
            "\"loss_node_range\":[","\"propagation_section_range\":[","\"checked_max_complex_relative_error\":",
            "\"extra_inviscid_inertia_compliance\":false","\"ratio\":2"] {
            assert!(sidecar.contains(token),"{token}: {sidecar}");
        }
    }
    let output=dir.join("lossy-37.wav");assert!(!run(&input,&output,37,true).status.success());
    assert_eq!(std::fs::read(&output).unwrap(),wav);
}

#[test]
fn viscothermal_valve_ensemble_preserves_physical_pressure_and_has_no_source_scale_gain() {
    let dir=directory();let input=dir.join("lossy.performance");
    let mut previous=None;
    for (i,text) in [LOSSY.to_string(),LOSSY.replace("9602 0.01","9602 100")].iter().enumerate() {
        std::fs::write(&input,text).unwrap();let output=dir.join(format!("ensemble-{i}.wav"));
        success(Command::new(env!("CARGO_BIN_EXE_music_render")).arg("ensemble").arg(&output)
            .arg("--valve").arg(&input).args(["--full-scale-pa","0.01","--decimate","--block","37"]).output().unwrap());
        let actual=std::fs::read(&output).unwrap();assert_eq!(actual,expected(LOSSY));
        if let Some(ref p)=previous {assert_eq!(&actual,p);}previous=Some(actual);
        let sidecar=std::fs::read_to_string(output.with_extension("provenance.json")).unwrap();
        assert!(sidecar.contains("wide-tube-zk-passive-rl-rc-v1"));
        assert!(sidecar.contains("\"source_pcm_scale_applied\":false"));
    }
}

#[test]
fn invalid_gas_loss_inputs_fail_before_the_cli_creates_any_output() {
    let dir=directory();let input=dir.join("bad.performance");
    for (i,text) in [LOSSY.replace("viscothermal 100 1000 12 8","viscothermal 100 1000 1 8"),
        LOSSY.replace("viscothermal 100 1000 12 8","viscothermal 100 1000 12 9"),
        LOSSY.replace("0.002 1048576 viscothermal","0.002 1024 viscothermal"),
        LOSSY.replace("8 32 1000","8 32 1500"),LOSSY.replace("viscothermal 100","viscothermal 0")].iter().enumerate() {
        std::fs::write(&input,text).unwrap();let output=dir.join(format!("bad-{i}.wav"));
        assert!(!run(&input,&output,37,true).status.success());assert!(!output.exists());
        assert!(!output.with_extension("provenance.json").exists());
    }
    std::fs::write(&input,LOSSY).unwrap();let output=dir.join("implicit-decimation.wav");
    assert!(!run(&input,&output,37,false).status.success());assert!(!output.exists());
}
