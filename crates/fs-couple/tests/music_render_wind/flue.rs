//! Supplied jet/card/gas data through the actual physical renderer and CLI.
use super::*;
use fs_aeroac::jetcard::mint_tonal_interim_card;
use fs_couple::flue_loop::network::performance::file::FlueFilePerformance;
use fs_exec::CancelGate;

const FLUE: &str = include_str!("../../examples/flue.performance");
mod network;
fn short() -> String {
    let head = FLUE.split_once("\ngesture\n").unwrap().0
        .replace("output 24000 5", "output 1024 1");
    format!("{head}\ngesture\nfrankensim-gesture-schedule-v1\ncontrol_rate_hz\t1000\ntracks\t1\ntrack\tblow\tblowing-pressure\ninitial\tpressure-pa\t3e2\nevents\t0\n")
}
fn physical(text: &str) -> (Vec<f64>, u64) {
    let p = FlueFilePerformance::from_bytes(text.as_bytes(), 37, &CancelGate::new_clock_free()).unwrap();
    let info = p.info();
    let mut renderer = p.into_renderer();
    let mut pressure = vec![0.0; info.samples as usize];
    for block in pressure.chunks_mut(37) { renderer.block(block).unwrap(); }
    assert_eq!(renderer.voice().samples_rendered(), 2*info.samples);
    assert!(pressure.iter().all(|p| p.is_finite()));
    assert!(pressure.iter().any(|p| p.abs() > 1e-10));
    (pressure, info.jet_card_hash)
}

#[test]
fn flue_file_reaches_the_actual_cli_receiver_and_sole_pcm_owner() {
    let dir = directory(); let input = dir.join("flue.performance");
    let source = short(); std::fs::write(&input, &source).unwrap();
    let (pressure, _) = physical(&source);
    let (expected, _) = encode_pcm16_wav(&pressure, 48_000, 1.0).unwrap();
    for block in [1, 37, 512] {
        let output = dir.join(format!("flue-{block}.wav"));
        success(run(&input, &output, block, false));
        assert_eq!(std::fs::read(&output).unwrap(), expected);
        let sidecar = std::fs::read_to_string(output.with_extension("provenance.json")).unwrap();
        for token in ["frankensim-flue-performance-v1", "\"fixture\":\"flue-input\"",
            "\"mechanics_sample_rate_hz\":96000", "\"mechanics_samples\":2048",
            "\"output_sample_rate_hz\":48000", "\"composition_authority\":\"EstimateOnly\"",
            "\"viscothermal_sections\":1", "one-way infinite-baffle outlet pressure"] {
            assert!(sidecar.contains(token), "{token}: {sidecar}");
        }
    }
}

#[test]
fn inline_card_and_data_only_gas_and_geometry_changes_reach_the_physical_source() {
    let source = short(); let reference = physical(&source);
    let card = mint_tonal_interim_card().unwrap().to_json();
    assert_eq!(physical(&source.replace("jet_card tonal-interim", &format!("jet_card json {card}"))), reference);
    let energy: f64 = reference.0.iter().map(|p| p*p).sum();
    for changed in [source.replace("gas 293.15 101325 0", "gas mixture 293.15 101325 mole 2 N2 0.6 O2 0.4"),
        source.replace("exterior same", "exterior 323.15 101325 0"),
        source.replace("tube 0.3", "tube 0.32"),
        source.replace("windway none", "windway slit 0.015")] {
        let pressure = physical(&changed).0;
        let difference: f64 = pressure.iter().zip(&reference.0).map(|(a,b)| (a-b).powi(2)).sum();
        assert!(difference > 1e-10*energy, "declared physical input must alter pressure");
    }
}

#[test]
fn invalid_flue_sources_and_duplicate_observation_refuse_before_output_creation() {
    let dir = directory(); let input = dir.join("flue.performance"); let source = short();
    for (i, text) in [source.replace("tonal-interim", "none"),
        source.replace("96000 48000", "95000 48000"), source.replace("gas 293.15", "gas NaN"),
        source.replace("exterior same", "exterior 293.15 80000 0"),
        source.replace("windway none", "windway slit -0.01"),
        source.replace("windway none", "windway guessed"),
        source.replace("gas 293.15 101325 0", "gas thermochemical 293.15 101325 mole 1 N2 1"),
        format!("{source}ignored\n")].iter().enumerate() {
        std::fs::write(&input, text).unwrap(); let output = dir.join(format!("refused-{i}.wav"));
        assert!(!run(&input, &output, 37, false).status.success());
        assert!(!output.exists()); assert!(!output.with_extension("provenance.json").exists());
    }
    std::fs::write(&input, source).unwrap();
    for flag in ["--decimate", "--acceptance"] {
        let output = dir.join(format!("{flag}.wav"));
        let result = Command::new(env!("CARGO_BIN_EXE_music_render")).arg("wind").arg(&input).arg(&output).arg(flag).output().unwrap();
        assert!(!result.status.success()); assert!(!output.exists());
        assert!(!output.with_extension("provenance.json").exists());
    }
}

#[test]
fn declared_windway_reaches_the_source_and_cli_pressure_without_a_second_model() {
    let source = short();
    // Old V1 files retain exactly the original ideal jet when the row is absent.
    assert_eq!(physical(&source), physical(&source.replace("windway none\n", "")));
    let source = source.replace("windway none", "windway slit 0.015");
    let p = FlueFilePerformance::from_bytes(source.as_bytes(), 37, &CancelGate::new_clock_free()).unwrap();
    let gas = p.renderer().voice().network().unwrap().source().gases[0];
    let law = p.renderer().voice().viscous_windway().unwrap();
    assert_eq!(law.length_m(), 0.015);
    assert_eq!(law.gas(), Some(&gas));
    let dir = directory();
    let input = dir.join("windway.performance");
    let output = dir.join("windway.wav");
    std::fs::write(&input, &source).unwrap();
    let expected = encode_pcm16_wav(&physical(&source).0, 48_000, 1.0).unwrap().0;
    success(run(&input, &output, 37, false));
    assert_eq!(std::fs::read(&output).unwrap(), expected);
    let sidecar = std::fs::read_to_string(output.with_extension("provenance.json")).unwrap();
    assert!(sidecar.contains("\"windway_length_m\":1.5e-2"));
    assert!(sidecar.contains("\"windway_dynamic_viscosity_pa_s\":"));
}

#[path = "flue/material.rs"]
mod material;
