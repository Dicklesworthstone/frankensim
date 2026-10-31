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

fn side_aperture_input() -> String {
    include_str!("../../examples/resolved-dg-side-aperture.performance")
        .replace("output 2048 20", "output 256 0.1")
        .replace("boundary pressure 3\nhold 128 20\nhold 256 -20\nhold 2048 0",
            "boundary pressure 2\nhold 128 20\nhold 256 0")
}

#[test]
fn resolved_cli_keeps_original_pressure_decimation_and_physical_source_receipts() {
    let dir = directory();
    // Temperature changes the ORIGINAL lattice/DG clock at the same mesh and
    // microphone timing. No field dimensions, source duration or output gain
    // change in the second description.
    for (case, text) in [input(), input().replace("gas dry-air 300 100000", "gas dry-air 310 100000"), side_aperture_input()]
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
            if case == 2 {
                assert!(metadata.contains("\"reservoir_observation\":{\"kind\":\"original-tile-continuum-flow-Rayleigh\""));
                assert!(metadata.contains("\"grid_cells\":[[5,2,1]]"));
                assert!(metadata.contains("\"bore_observation_scope\":"));
                assert!(metadata.contains("\"reservoir_aperture\":{\"kind\":\"prescribed-exterior-ghost-state\""));
                assert!(metadata.contains("\"cell_count\":1"));
                assert!(metadata.contains("\"mass_account\":{\"field\":\"reservoir_mass_out_kg\",\"unit\":\"kg\",\"positive\":\"out of resolved fluid\"}"));
                assert!(metadata.contains("\"energy_account\":{\"field\":\"flow.open_energy_j\",\"unit\":\"J\",\"positive\":\"into resolved fluid\"}"));
            }
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

#[test]
fn both_aperture_roles_survive_native_render_and_sidecar_on_either_face() {
    let dir = directory();
    for high in [false, true] {
        let text = include_str!("../../examples/resolved-dg-split-aperture.performance")
            .replace("output 2048 20", "output 4 20")
            .replace("characteristic 0 high", if high { "characteristic 0 high" } else { "characteristic 0 low" })
            .replace("coupled_patch 1 2 1 2 pressure 0",
                "coupled_patch 1 2 1 2 characteristic 0\nsource_patch 0 2 0 1 pressure 0")
            .replace("boundary pressure 3\nhold 128 20\nhold 256 -20\nhold 2048 0",
                "boundary pressure 1\nhold 4 0");
        let source = ResolvedDgFilePerformance::from_bytes(text.as_bytes(), 4,
            Arc::new(CancelGate::new_clock_free())).unwrap();
        assert_eq!(source.renderer().system().coupled_face_is_high(), high);
        let mut renderer = source.into_renderer();
        let mut pressure = [0.0; 4];
        renderer.block(&mut pressure).unwrap();
        let expected = encode_pcm16_wav(&pressure, 96000, 20.0).unwrap().0;
        let input = dir.join(format!("two-apertures-{high}.performance"));
        let output = dir.join(format!("two-apertures-{high}.wav"));
        std::fs::write(&input, text).unwrap();
        success(Command::new(env!("CARGO_BIN_EXE_music_render"))
            .arg("resolved").arg(&input).arg(&output)
            .args(["--block", "3", "--self-contained"]).output().unwrap());
        assert_eq!(std::fs::read(&output).unwrap(), expected);
        let metadata = std::fs::read_to_string(output.with_extension("provenance.json")).unwrap();
        let coupled = metadata.split("\"coupled_patch\":{").nth(1).unwrap()
            .split("\"source_patch\":{").next().unwrap();
        let source = metadata.split("\"source_patch\":{").nth(1).unwrap()
            .split("\"external_boundary\":{").next().unwrap();
        assert!(coupled.contains(&format!("\"high_face\":{high}")));
        assert!(source.contains(&format!("\"high_face\":{}", !high)));
        assert!(coupled.contains("\"start\":[1,1],\"end\":[2,2]"));
        assert!(source.contains("\"start\":[0,0],\"end\":[2,1]"));
        assert!(coupled.contains("\"kind\":\"incoming-characteristic\""));
        assert!(source.contains("\"kind\":\"prescribed-pressure-offset\""));
    }
}

#[test]
fn resolved_command_accepts_native_thermoviscous_field_schema() {
    let dir = directory();
    let text = include_str!("../../examples/resolved-thermoviscous.performance")
        .replace("output 1920 1", "output 16 1")
        .replace("ramp 384 0 293.15 20 295\nhold 768 20 295\nramp 1152 20 295 0 293.15\nhold 1920 0 293.15",
            "hold 16 0 293.15");
    let input = dir.join("resolved-thermoviscous.performance");
    let output = dir.join("resolved-thermoviscous.wav");
    std::fs::write(&input, text).unwrap();
    success(Command::new(env!("CARGO_BIN_EXE_music_render"))
        .arg("resolved").arg(&input).arg(&output)
        .args(["--block", "4", "--decimate", "--self-contained"]).output().unwrap());
    assert!(output.exists());
    let metadata = std::fs::read_to_string(output.with_extension("provenance.json")).unwrap();
    for token in ["\\\"fixture\\\":\\\"resolved-thermoviscous-input\\\"",
        "\\\"resolved_thermoviscous_input\\\":{",
        "frankensim-resolved-thermoviscous-dg-performance-v1",
        "\\\"vertices\\\":8", "\\\"tetrahedra\\\":6",
        "original nonlinear NSF source into complete 3-D thermoviscous field"]
    {
        assert!(metadata.contains(token), "{token}: {metadata}");
    }
}
