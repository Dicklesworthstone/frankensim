use super::*;
use fs_material::gas::GasState;
const REGIONAL:&str=include_str!("../../examples/plate-valve-regional-gas.performance");

#[test]
fn regional_wind_and_ensemble_write_the_actual_pcm_and_outlet_medium() {
    let dir=directory();let input=dir.join("regional.performance");std::fs::write(&input,REGIONAL).unwrap();
    let (pressure,_)=reference(REGIONAL);let expected=encode_pcm16_wav(&pressure,48000,0.01).unwrap().0;
    assert!(expected[44..].iter().any(|&b|b!=0));
    let outlet=GasState::try_new_moist_air(313.15,101325.0,0.5).unwrap();
    for block in [37,512] {
        let out=dir.join(format!("regional-{block}.wav"));success(run(&input,&out,block,true));
        assert_eq!(std::fs::read(&out).unwrap(),expected);
        let meta=std::fs::read_to_string(out.with_extension("provenance.json")).unwrap();
        assert!(meta.contains("\"gas\":{"));assert!(meta.contains("wide-tube-zk-passive-rl-rc-v1"));
        let receiver=meta.split("\"outlet_receiver\":").nth(1).unwrap().split("\"radial_rings\"").next().unwrap();
        assert!(receiver.contains(&format!("\"density_kg_m3\":{:e}",outlet.density)));
        assert!(receiver.contains(&format!("\"sound_speed_m_s\":{:e}",outlet.sound_speed)));
    }
    let changed_scale=REGIONAL.replace("9602 0.01","9602 100");std::fs::write(&input,changed_scale).unwrap();
    let out=dir.join("ensemble.wav");
    success(Command::new(env!("CARGO_BIN_EXE_music_render")).arg("ensemble").arg(&out).arg("--valve").arg(&input)
        .args(["--decimate","--full-scale-pa","0.01","--block","37"]).output().unwrap());
    assert_eq!(std::fs::read(&out).unwrap(),expected,"source scale is not a per-region gain");
}

#[test]
fn regional_source_refusals_precede_output_creation() {
    let dir=directory();let input=dir.join("bad.performance");
    for (i,field) in ["gas 313.15 2","gas NaN 0.5","gas 313.15 0.5 gas 293.15 0"].iter().enumerate() {
        std::fs::write(&input,REGIONAL.replace("gas 313.15 0.5",field)).unwrap();
        let out=dir.join(format!("bad-{i}.wav"));assert!(!run(&input,&out,37,true).status.success());
        assert!(!out.exists());assert!(!out.with_extension("provenance.json").exists());
    }
}
