//! File, gas, observer, schedule and actual CLI composition over existing owners.
use super::*;
use fs_couple::render::schedule::reed::ReedPerformanceError;

const SOURCE: &str = include_str!("../../examples/reed-baffled-mixture.performance");
const RECEIVER: &str = "baffled_outlet 0 0 0.5 4 16 1000";

#[test]
fn outlet_file_reaches_native_and_decimated_wav_without_the_proxy() {
    let dir = directory();
    for high in [false, true] {
        let text = if high { high_rate(SOURCE) } else { SOURCE.to_string() };
        let input = dir.join(format!("outlet-{high}.performance"));
        std::fs::write(&input, &text).unwrap();
        let source = ReedPerformance::from_bytes(text.as_bytes(), 37).unwrap();
        let info = source.outlet_receiver().unwrap();
        assert_eq!(info.gas, source.segment_gases()[1]);
        assert_eq!(info.sample_rate_hz, source.info().sample_rate_hz);
        assert_eq!(info.receiver.position_m, [0.0, 0.0, 0.5]);
        let raw = pressure(&text);
        let (expected, _) = encode_pcm16_wav(&raw, 48000, 2.0).unwrap();
        assert!(expected[44..].iter().any(|byte| *byte != 0));
        let proxy = text.replace(RECEIVER, "listener_m 0.5");
        assert!(ReedPerformance::from_bytes(proxy.as_bytes(), 37).unwrap().outlet_receiver().is_none());
        assert_ne!(raw, pressure(&proxy));
        for block in [37, 113] {
            let output = dir.join(format!("outlet-{high}-{block}.wav"));
            success(wind(&input, &output, block, high));
            assert_eq!(std::fs::read(&output).unwrap(), expected);
            let meta = std::fs::read_to_string(output.with_extension("provenance.json")).unwrap();
            assert!(meta.contains("\"pressure_point\":\"baffled-outlet\""));
            assert!(meta.contains("terminal-only Rayleigh pressure"));
            assert!(meta.contains("\"priming_source_sample\":-1"));
            assert!(meta.contains("\"axial_delay_mechanical_samples\":"));
            assert!(meta.contains("\"exterior_delay_mechanical_samples\":"));
            assert!(!meta.contains("bore-pressure plus compact-jet proxy"));
        }
    }
}

#[test]
fn downstream_composition_and_common_pressure_reach_terminal_sound() {
    let unchanged = SOURCE.replace(
        "gas mixture 293.15 100000 mole 3 N2 0.4 O2 0.2 CO2 0.4",
        "gas mixture 293.15 100000 mole 2 N2 0.7 O2 0.3",
    );
    let low = SOURCE.replace("mixture 293.15 100000", "mixture 293.15 50000");
    let base = ReedPerformance::from_bytes(SOURCE.as_bytes(), 37).unwrap();
    let thin = ReedPerformance::from_bytes(low.as_bytes(), 37).unwrap();
    let a = base.outlet_receiver().unwrap();
    let b = thin.outlet_receiver().unwrap();
    assert_eq!(a.gas.sound_speed, b.gas.sound_speed);
    assert_eq!(a.gas.density, 2.0 * b.gas.density);
    assert_eq!(a.radius_m, b.radius_m);
    assert_eq!(a.receiver, b.receiver);
    let actual = pressure(SOURCE);
    for other in [&unchanged, &low] {
        let compared = pressure(other);
        let signal: f64 = actual.iter().map(|p| p * p).sum();
        let difference: f64 = actual.iter().zip(&compared).map(|(a, b)| (a - b).powi(2)).sum();
        assert!(signal > 0.0 && difference > 1e-10 * signal);
    }
}

#[test]
fn invalid_or_mismatched_observers_refuse_before_either_output_exists() {
    let dir = directory();
    let input = dir.join("bad-outlet.performance");
    let mut cases: Vec<String> = ["closed", "ideal-open", "unflanged"].iter()
        .map(|term| SOURCE.replace("termination flanged", &format!("termination {term}"))).collect();
    for record in [
        "baffled_outlet 0 0 0.5 4 16", "baffled_outlet 0 0 0.5 4 16 1000 extra",
        "baffled_outlet NaN 0 0.5 4 16 1000", "baffled_outlet 0 0 0.01 4 16 1000",
        "baffled_outlet 0 0 0.5 0 16 1000", "baffled_outlet 0 0 0.5 129 16 1000",
        "baffled_outlet 0 0 0.5 4 7 1000", "baffled_outlet 0 0 0.5 4 513 1000",
        "baffled_outlet 0 0 0.5 4 16 10000", "baffled_outlet 0 0 0.5 4 16 0",
        "baffled_outlet 0 0 0.5 18446744073709551615 16 1000",
    ] { cases.push(SOURCE.replace(RECEIVER, record)); }
    for (index, text) in cases.iter().enumerate() {
        assert!(ReedPerformance::from_bytes(text.as_bytes(), 37).is_err());
        std::fs::write(&input, text).unwrap();
        let output = dir.join(format!("bad-{index}.wav"));
        assert!(!wind(&input, &output, 37, false).status.success());
        assert!(!output.exists());
        assert!(!output.with_extension("provenance.json").exists());
    }
    let invalid = SOURCE.replace(RECEIVER, "baffled_outlet 0 0 0.01 4 16 1000");
    let error = match ReedPerformance::from_bytes(invalid.as_bytes(), 37) {
        Ok(_) => panic!("invalid receiver admitted"), Err(error) => error,
    };
    assert!(matches!(&error, ReedPerformanceError::Outlet(_)));
    assert!(std::error::Error::source(&error).is_some());
}

#[test]
fn outlet_ensemble_sums_actual_pressure_not_source_encoding_scales() {
    let dir = directory();
    let first = dir.join("first.performance");
    let second = dir.join("second.performance");
    // One listener gas (ensemble guard): the second part differs by its receiver.
    let far = SOURCE.replace(RECEIVER, "baffled_outlet 0 0 0.7 4 16 1000");
    std::fs::write(&second, &far).unwrap();
    let mut mixed = pressure(SOURCE);
    for (a, b) in mixed.iter_mut().zip(pressure(&far)) { *a += b; }
    let (expected, _) = encode_pcm16_wav(&mixed, 48000, 4.0).unwrap();
    for (index, scale) in ["2", "2000000"].iter().enumerate() {
        let text = SOURCE.replace("audio 48000 4801 2", &format!("audio 48000 4801 {scale}"));
        std::fs::write(&first, text).unwrap();
        let output = dir.join(format!("ensemble-{index}.wav"));
        success(Command::new(env!("CARGO_BIN_EXE_music_render")).arg("ensemble").arg(&output)
            .arg("--reed").arg(&first).arg("--reed").arg(&second)
            .args(["--block", "37", "--full-scale-pa", "4"]).output().unwrap());
        assert_eq!(std::fs::read(&output).unwrap(), expected);
        let meta = std::fs::read_to_string(output.with_extension("provenance.json")).unwrap();
        assert_eq!(meta.matches("\"pressure_point\":\"baffled-outlet\"").count(), 2);
        assert_eq!(meta.matches("\"priming_source_sample\":-1").count(), 2);
    }
}

#[test]
fn separate_exterior_material_reaches_actual_cli_pcm_and_invalid_records_create_nothing() {
    let dir = directory();
    let input = dir.join("external.performance");
    let text = SOURCE.replace(RECEIVER,
        &format!("{RECEIVER}\nexterior mixture 300 100000 mass 2 N2 0.7 O2 0.3"));
    let admitted = ReedPerformance::from_bytes(text.as_bytes(), 37).unwrap();
    assert_ne!(admitted.outlet_receiver().unwrap().gas, admitted.segment_gases()[1]);
    let (expected, _) = encode_pcm16_wav(&pressure(&text), 48000, 2.0).unwrap();
    std::fs::write(&input, &text).unwrap();
    let output = dir.join("external.wav");
    success(wind(&input, &output, 113, false));
    assert_eq!(std::fs::read(&output).unwrap(), expected);
    let metadata = std::fs::read_to_string(output.with_extension("provenance.json")).unwrap();
    assert!(metadata.contains("terminal-only Rayleigh pressure in the bound exterior gas"));
    for (index, record) in ["exterior 300 90000 0",
        "exterior mixture 300 100000 mole 1 Xe 1", "exterior 300 100000 0 extra"].iter().enumerate() {
        std::fs::write(&input, SOURCE.replace(RECEIVER, &format!("{RECEIVER}\n{record}"))).unwrap();
        let output = dir.join(format!("bad-exterior-{index}.wav"));
        assert!(!wind(&input, &output, 37, false).status.success());
        assert!(!output.exists());
        assert!(!output.with_extension("provenance.json").exists());
    }
}

#[test]
fn data_only_two_medium_tone_holes_reach_cli_pcm_at_both_block_sizes() {
    let dir = directory();
    let fixture = include_str!("../../examples/reed-two-medium-tone-hole.performance");
    for sigma in [1.0, 0.3] {
        let exterior = fixture.replace("hole 0.0005 0.002 0.0015 0.3",
            &format!("hole 0.0005 0.002 0.0015 {sigma}"));
        let text = exterior.replace("exterior mixture 303.15 100000 mole 2 N2 0.8 O2 0.2\n", "");
        let source = ReedPerformance::from_bytes(exterior.as_bytes(), 37).unwrap();
        assert_eq!(source.segment_gases().len(), 4);
        assert_ne!(source.outlet_receiver().unwrap().gas, source.segment_gases()[1]);
        let pcm = pressure(&exterior);
        assert!(pcm.iter().all(|x| x.is_finite()));
        assert_ne!(pcm, pressure(&text), "only the exterior material record changed");
        let (expected, _) = encode_pcm16_wav(&pcm, 48000, 2.0).unwrap();
        assert!(expected[44..].iter().any(|byte| *byte != 0));
        let input = dir.join(format!("two-medium-{sigma}.performance"));
        std::fs::write(&input, exterior).unwrap();
        for block in [37, 113] {
            let output = dir.join(format!("two-medium-{sigma}-{block}.wav"));
            success(wind(&input, &output, block, false));
            assert_eq!(std::fs::read(&output).unwrap(), expected);
            let metadata = std::fs::read_to_string(output.with_extension("provenance.json")).unwrap();
            assert!(metadata.contains("\"pressure_point\":\"baffled-outlet\""));
            assert!(metadata.contains("bound exterior gas"));
        }
    }
}
