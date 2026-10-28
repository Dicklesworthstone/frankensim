//! The native entry point must use the same original physical owners and PCM.
use super::{directory, success};
use fs_couple::nsf_dg::performance::file::ResolvedDgFilePerformance;
use fs_couple::pcm_wav::encode_pcm16_wav;
use fs_couple::pcm_wav::observation::{DecimatedRenderer, PressureRenderer};
use fs_exec::CancelGate;
use std::process::Command;
use std::sync::Arc;

fn input() -> String {
    include_str!("../../examples/resolved-dg.performance")
        .replace("output 2048 20", "output 256 0.1")
        .replace("boundary 3\nhold 128 -0.1\nhold 256 0.1\nhold 2048 0", "boundary 2\nhold 128 -0.1\nhold 256 0")
}

#[test]
fn resolved_cli_keeps_original_pressure_decimation_and_physical_source_receipts() {
    let dir = directory();
    // Temperature changes the ORIGINAL lattice/DG clock at the same mesh and
    // microphone timing. No field dimensions, source duration or output gain
    // change in the second description.
    for (case, text) in [input(), input().replace("gas dry-air 300 100000", "gas dry-air 310 100000")]
        .into_iter().enumerate()
    {
        let source = ResolvedDgFilePerformance::from_bytes(text.as_bytes(), 37,
            Arc::new(CancelGate::new_clock_free())).unwrap();
        let info = source.info();
        let dt = source.renderer().system().load().propagator().dt();
        let mut renderer = DecimatedRenderer::new(source.into_renderer(), info.sample_rate_hz, 48000, 37).unwrap();
        let mut pressure = vec![0.0; renderer.output_samples_for(info.samples).unwrap() as usize];
        for chunk in pressure.chunks_mut(37) { renderer.block(chunk).unwrap(); }
        let expected = encode_pcm16_wav(&pressure, 48000, info.full_scale_pa).unwrap().0;
        let input = dir.join(format!("resolved-{case}.performance"));
        std::fs::write(&input, &text).unwrap();
        for block in [1, 37] {
            let output = dir.join(format!("resolved-{case}-{block}.wav"));
            let result = Command::new(env!("CARGO_BIN_EXE_music_render"))
                .arg("resolved").arg(&input).arg(&output)
                .args(["--block", &block.to_string(), "--decimate", "--self-contained"])
                .output().unwrap();
            success(result);
            assert_eq!(std::fs::read(&output).unwrap(), expected);
            let metadata = std::fs::read_to_string(output.with_extension("provenance.json")).unwrap();
            for token in ["\"fixture\":\"resolved-dg-input\"", "\"resolved_input\":{",
                "original-NsfGrid3-D3Q27", "\"solid_cells\":4", "\"pressure_point\":\"baffled-outlet\"",
                "\"observer_gas\":{", "numerical-fixture-Cramer2012", "\"sample_rate_hz\":48000",
                "\"observation_sample_rate_hz\":96000", "\"kernel_support_observation_samples\":",
                "\"declared_physical_intervals\":", "\"bandlimited_antialias_filter\":false"]
            {
                assert!(metadata.contains(token), "{token}: {metadata}");
            }
            assert!(metadata.contains(&format!("\"nsf_dt_s\":{dt:e}")));
        }
    }
}

#[test]
fn resolved_route_refuses_a_different_source_or_unqualified_acceptance_before_audio() {
    let dir = directory();
    for (case, text, flags, expected) in [
        (0, include_str!("../../examples/flue-dg.performance").to_owned(), vec![], "requires frankensim-resolved-dg"),
        (1, input(), vec!["--acceptance"], "requires source geometry"),
        (2, input().replace("gas dry-air 300 100000", "gas thermochemical 300 100000 mole 1 N2 1"),
            vec![], "original thermochemical source packs"),
    ] {
        let input = dir.join(format!("refused-{case}.performance"));
        let output = dir.join(format!("refused-{case}.wav"));
        std::fs::write(&input, text).unwrap();
        let result = Command::new(env!("CARGO_BIN_EXE_music_render"))
            .arg("resolved").arg(&input).arg(&output).arg("--decimate").args(flags).output().unwrap();
        assert!(!result.status.success());
        let message = format!("{}{}", String::from_utf8_lossy(&result.stdout), String::from_utf8_lossy(&result.stderr));
        assert!(message.contains(expected), "{message}");
        assert!(!output.exists());
        assert!(!output.with_extension("provenance.json").exists());
    }
}
