use super::{directory, run, success};
use fs_couple::bernoulli_aperture::performance::file::PlateValvePerformance;
use fs_couple::pcm_wav::encode_pcm16_wav;
use fs_couple::pcm_wav::observation::{DecimatedRenderer, PressureRenderer};
use fs_exec::CancelGate;
use std::process::Command;

fn input() -> String {
    let example = include_str!("../../examples/plate-valve-dg.performance");
    let header = example
        .split_once("\nschedule\n")
        .unwrap()
        .0
        .replace("audio 96000 9602", "audio 96000 512");
    format!(
        "{header}\nschedule\nfrankensim-gesture-schedule-v1\ncontrol_rate_hz\t700\ntracks\t1\ntrack\tmouth\tblowing-pressure\ninitial\tpressure-pa\t5e0\nevents\t0\n"
    )
}

#[test]
fn exact_dg_valve_reaches_wind_and_ensemble_pcm_through_the_existing_renderer() {
    let data = input();
    let source =
        PlateValvePerformance::from_bytes(data.as_bytes(), 37, &CancelGate::new()).unwrap();
    let info = source.info();
    let mut renderer =
        DecimatedRenderer::new(source.into_renderer(), info.sample_rate_hz, 48000, 37).unwrap();
    let mut pressure = vec![0.0; renderer.output_samples_for(info.samples).unwrap() as usize];
    for chunk in pressure.chunks_mut(37) {
        renderer.block(chunk).unwrap();
    }
    let expected = encode_pcm16_wav(&pressure, 48000, info.full_scale_pa)
        .unwrap()
        .0;
    assert!(expected[44..].iter().any(|b| *b != 0));
    let directory = directory();
    let input = directory.join("dg.performance");
    std::fs::write(&input, &data).unwrap();
    for block in [1, 37] {
        let output = directory.join(format!("dg-{block}.wav"));
        success(run(&input, &output, block, true));
        assert_eq!(std::fs::read(&output).unwrap(), expected);
        let metadata = std::fs::read_to_string(output.with_extension("provenance.json")).unwrap();
        for token in [
            "exact-zoh-bessel-dg-v1",
            "\"total_dofs\":108",
            "\"series_poles\":8",
            "\"thermal_poles\":8",
            "\"pressure_point\":\"inlet\"",
            "\"memory_branches\":2",
            "not an exterior microphone",
            // faba3eb39 renamed the disclaimer from "wall/terminal loss" to "loss ledger".
            "loss ledger inferred, not an independent heat integral",
        ] {
            assert!(metadata.contains(token), "{token}: {metadata}");
        }
    }
    let output = directory.join("dg-ensemble.wav");
    success(
        Command::new(env!("CARGO_BIN_EXE_music_render"))
            .arg("ensemble")
            .arg(&output)
            .arg("--valve")
            .arg(&input)
            .args(["--full-scale-pa", "20", "--decimate", "--block", "37"])
            .output()
            .unwrap(),
    );
    assert_eq!(std::fs::read(output).unwrap(), expected);
}

#[test]
fn unbound_dg_exterior_observation_refuses_before_output_creation() {
    let directory = directory();
    let path = directory.join("bad-dg.performance");
    std::fs::write(
        &path,
        input().replace("observation inlet", "observation terminal"),
    )
    .unwrap();
    let output = directory.join("bad-dg.wav");
    assert!(!run(&path, &output, 37, true).status.success());
    assert!(!output.exists());
    assert!(!output.with_extension("provenance.json").exists());
}
