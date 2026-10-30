//! Generated gas reaches the existing command, exterior receiver and PCM owner.
use fs_couple::flue_loop::dg::performance::file::FlueDgFilePerformance;
use fs_couple::pcm_wav::{encode_pcm16_wav, observation::{DecimatedRenderer, PressureRenderer}};
use fs_exec::CancelGate;
use fs_material::gas_species::input::{ThermochemicalPhase, ThermochemicalSources};
use std::process::Command;

#[path = "../src/render/schedule/reed/field_input/tests/support.rs"]
mod support;

#[test]
fn authored_steady_dg_gas_reaches_original_microphone_decimation_and_pcm() {
    let a = support::packs("N2", 0.028, 3.5, 200.0);
    let b = support::packs("CO2", 0.044, 4.5, 200.0);
    let bindings = [a.binding(), b.binding()];
    let sources = ThermochemicalSources { species: &bindings, water_species: None, phase: ThermochemicalPhase::FrozenGas };
    let text = include_str!("../examples/flue-dg-steady-gas.performance");
    let source = FlueDgFilePerformance::from_bytes_with_gas_sources(
        text.as_bytes(), 37, &CancelGate::new_clock_free(), Some(&sources),
    ).unwrap();
    let load = source.renderer().voice().dg().unwrap().load();
    let field = load.duct().gas_field().unwrap();
    assert_eq!(field.stations_m(), &[0.0, 0.025, 0.05]);
    assert_ne!(load.duct().inlet_gas(), load.duct().outlet_gas());
    assert_ne!(load.duct().outlet_gas(), load.radiation_record().unwrap().gas_record().state());
    let info = source.info();
    let mut renderer = DecimatedRenderer::new(source.into_renderer(), info.sample_rate_hz, 48000, 37).unwrap();
    let mut pressure = vec![0.0; renderer.output_samples_for(info.samples).unwrap() as usize];
    for chunk in pressure.chunks_mut(37) { renderer.block(chunk).unwrap(); }
    assert!(pressure.iter().any(|p| p.abs() > 1e-12));
    let expected = encode_pcm16_wav(&pressure, 48000, info.full_scale_pa).unwrap().0;
    let stamp = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos();
    let dir = std::env::temp_dir().join(format!("fs-dg-steady-{}-{stamp}", std::process::id()));
    std::fs::create_dir(&dir).unwrap();
    let input = dir.join("source.performance");
    let output = dir.join("pressure.wav");
    std::fs::write(&input, text).unwrap();
    let mut command = Command::new(env!("CARGO_BIN_EXE_music_render"));
    command.arg("wind").arg(&input).arg(&output)
        .args(["--block", "37", "--decimate", "--gas-water", "none", "--gas-phase", "frozen"]);
    for (index, packs) in [a, b].iter().enumerate() {
        let paths = [dir.join(format!("species-{index}.pack")), dir.join(format!("nasa9-{index}.pack")),
            dir.join(format!("transport-{index}.pack"))];
        for (path, bytes) in paths.iter().zip([packs.species.to_bytes(), packs.nasa9.to_bytes(), packs.transport.to_bytes()]) {
            std::fs::write(path, bytes).unwrap();
        }
        command.arg("--gas-component").args(&paths);
    }
    let result = command.output().unwrap();
    assert!(result.status.success(), "{}{}", String::from_utf8_lossy(&result.stdout), String::from_utf8_lossy(&result.stderr));
    assert_eq!(std::fs::read(&output).unwrap(), expected);
    let metadata = std::fs::read_to_string(output.with_extension("provenance.json")).unwrap();
    for token in ["frankensim-flue-dg-performance-v1", "original-exact-dg", "baffled-outlet", "stokes-kirchhoff"] {
        assert!(metadata.contains(token), "{token}: {metadata}");
    }
}
