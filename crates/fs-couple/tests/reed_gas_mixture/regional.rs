use super::*;
const REGIONAL: &str = include_str!("../../examples/reed-regional-mixture.performance");
const OVERRIDE: &str = " gas mixture 293.15 100000 mole 3 N2 0.4 O2 0.2 CO2 0.4";

#[test]
fn regional_wind_and_ensemble_reach_the_same_native_and_decimated_pcm_owner() {
    let dir = directory();
    for high in [false, true] {
        let text = if high { high_rate(REGIONAL) } else { REGIONAL.to_string() };
        let uniform = text.replace(OVERRIDE, "");
        let input = dir.join(format!("regional-{high}.performance"));
        let other = dir.join(format!("uniform-{high}.performance"));
        std::fs::write(&input, &text).unwrap();
        std::fs::write(&other, &uniform).unwrap();
        let a = pressure(&text); let b = pressure(&uniform);
        let (expected, _) = encode_pcm16_wav(&a, 48000, 20000.0).unwrap();
        let (unchanged, _) = encode_pcm16_wav(&b, 48000, 20000.0).unwrap();
        assert_ne!(&expected[44..], &unchanged[44..], "downstream gas must reach PCM");
        for block in [37, 113] {
            let output = dir.join(format!("regional-{high}-{block}.wav"));
            success(wind(&input, &output, block, high));
            assert_eq!(std::fs::read(&output).unwrap(), expected);
        }
        let mixed: Vec<f64> = a.iter().zip(&b).map(|(p,q)| p+q).collect();
        let (expected, _) = encode_pcm16_wav(&mixed, 48000, 40000.0).unwrap();
        let output = dir.join(format!("ensemble-{high}.wav"));
        let mut command = Command::new(env!("CARGO_BIN_EXE_music_render"));
        command.arg("ensemble").arg(&output).arg("--reed").arg(&input).arg("--reed").arg(&other)
            .args(["--block", "37", "--full-scale-pa", "40000"]);
        if high { command.arg("--decimate"); }
        success(command.output().unwrap());
        assert_eq!(std::fs::read(output).unwrap(), expected);
    }
}

#[test]
fn invalid_regional_media_refuse_before_creating_audio_or_sidecar() {
    let dir = directory();
    let input = dir.join("regional.performance");
    for (index, record) in [" gas mixture 293.15 50000 mole 1 N2 1",
        " gas mixture 293.15 100000 mole 1 He 1", " gas 293.15 0 extra",
        " gas mixture 293.15 100000 mole 1 H2O 1"].into_iter().enumerate() {
        std::fs::write(&input, REGIONAL.replace(OVERRIDE, record)).unwrap();
        let output = dir.join(format!("refused-regional-{index}.wav"));
        assert!(!wind(&input, &output, 37, false).status.success());
        assert!(!output.exists());
        assert!(!output.with_extension("provenance.json").exists());
    }
}
