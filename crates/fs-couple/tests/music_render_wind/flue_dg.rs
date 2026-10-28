use super::{directory, run, success};
use fs_couple::flue_loop::dg::performance::file::FlueDgFilePerformance;
use fs_couple::pcm_wav::encode_pcm16_wav;
use fs_couple::pcm_wav::observation::{DecimatedRenderer, PressureRenderer};
use fs_exec::CancelGate;
use std::process::Command;

fn input() -> String {
    let example = include_str!("../../examples/flue-dg.performance");
    let header = example.split_once("\ngesture\n").unwrap().0
        .replace("output 2048", "output 512");
    format!("{header}\ngesture\nfrankensim-gesture-schedule-v1\ncontrol_rate_hz\t1000\ntracks\t1\ntrack\tblow\tblowing-pressure\ninitial\tpressure-pa\t2e1\nevents\t0\n")
}

#[test]
fn original_dg_flue_reaches_the_existing_wind_decimator_and_pcm_writer() {
    let data = input();
    let source = FlueDgFilePerformance::from_bytes(data.as_bytes(), 37,
        &CancelGate::new()).unwrap();
    let info = source.info();
    let mut renderer = DecimatedRenderer::new(source.into_renderer(),
        info.sample_rate_hz, 48000, 37).unwrap();
    let mut pressure = vec![0.0; renderer.output_samples_for(info.samples).unwrap() as usize];
    for chunk in pressure.chunks_mut(37) { renderer.block(chunk).unwrap(); }
    let expected = encode_pcm16_wav(&pressure, 48000, info.full_scale_pa).unwrap().0;
    assert!(expected[44..].iter().any(|byte| *byte != 0));
    let directory = directory();
    let path = directory.join("flue-dg.performance");
    std::fs::write(&path, &data).unwrap();
    for block in [1, 37] {
        let output = directory.join(format!("flue-dg-{block}.wav"));
        success(run(&path, &output, block, true));
        assert_eq!(std::fs::read(&output).unwrap(), expected);
        let metadata = std::fs::read_to_string(output.with_extension("provenance.json")).unwrap();
        for token in [
            "frankensim-flue-dg-performance-v1",
            "original-exact-dg",
            "\"pressure_point\":\"mean-inlet\"",
            "not an exterior microphone",
        ] {
            assert!(metadata.contains(token), "{token}: {metadata}");
        }
    }
    let output = directory.join("acceptance.wav");
    let refused = Command::new(env!("CARGO_BIN_EXE_music_render"))
        .arg("wind").arg(&path).arg(&output).args(["--decimate", "--acceptance"])
        .output().unwrap();
    assert!(!refused.status.success());
    assert!(!output.exists());
    assert!(!output.with_extension("provenance.json").exists());
}

#[test]
fn native_dg_flue_observer_gas_reaches_mouth_terminal_and_pcm() {
    let data = include_str!("../../examples/flue-dg-baffled.performance")
        .replace("output 2048 20", "output 512 0.02")
        .replace("50000000\n", "50000000\nmouth_radiation -0.03 0 4 2\n")
        .replace("events\t1\nevent\t1e-2\t0e0\tpressure-pa\t0e0\n", "events\t0\n");
    let source = FlueDgFilePerformance::from_bytes(data.as_bytes(), 37,
        &CancelGate::new_clock_free()).unwrap();
    assert!(source.renderer().observer_gas().is_some());
    assert!(source.renderer().mouth_observer_gas().is_some());
    let info = source.info();
    let mut renderer = DecimatedRenderer::new(source.into_renderer(),
        info.sample_rate_hz, 48000, 37).unwrap();
    let mut pressure = vec![0.0; 256];
    for chunk in pressure.chunks_mut(37) { renderer.block(chunk).unwrap(); }
    let expected = encode_pcm16_wav(&pressure, 48000, info.full_scale_pa).unwrap().0;
    let dir = directory();
    let input = dir.join("gas-mouth.performance");
    std::fs::write(&input, data).unwrap();
    let output = dir.join("gas-mouth.wav");
    success(run(&input, &output, 37, true));
    assert_eq!(std::fs::read(&output).unwrap(), expected);
    let metadata = std::fs::read_to_string(output.with_extension("provenance.json")).unwrap();
    for token in ["\"observer_gas\":{\"law\":\"stokes-kirchhoff\"",
        "\"mouth_observer_gas\":{\"law\":\"stokes-kirchhoff\"",
        "\"molecular_relaxation_omitted\":true", "\"additional_latency_samples\":"] {
        assert!(metadata.contains(token), "{token}");
    }
}
