//! A modal file's observer transfers come from its declared geometry and gas.
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::sync::atomic::{AtomicU64, Ordering};
use fs_couple::render::schedule::force::file::ModalPerformance;
use fs_couple::render::schedule::force::file::radiation::compact_baffled_transfer;
use fs_material::gas::GasState;
use fs_material::gas_species::input::read_gas_state;

const INPUT: &str = include_str!("../examples/modal-radiating.performance");
const AIR: &str = "ambient 293.15 101325 0";
static NEXT: AtomicU64 = AtomicU64::new(0);
fn directory() -> PathBuf {
    let time = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos();
    let path = std::env::temp_dir().join(format!("fs-modal-gas-{}-{time}-{}",
        std::process::id(), NEXT.fetch_add(1, Ordering::Relaxed)));
    std::fs::create_dir(&path).unwrap();
    path
}
fn invoke(input: &Path, output: &Path) -> Output {
    Command::new(env!("CARGO_BIN_EXE_music_render")).arg("modal").arg(input).arg(output)
        .args(["--block", "37"]).output().unwrap()
}
fn gas(record: &str) -> GasState {
    let record = record.strip_prefix("ambient ").unwrap();
    read_gas_state(&mut record.split_ascii_whitespace(), None).unwrap()
}
fn render(text: &str, samples: usize) -> Vec<f64> {
    let mut renderer = ModalPerformance::from_bytes(text.as_bytes(), 571).unwrap().into_renderer();
    let mut out = vec![0.0; samples];
    for chunk in out.chunks_mut(571) { renderer.block(chunk).unwrap(); }
    out
}
/// The same mechanics with each mode's transfer written into its H fields.
fn authored(text: &str, gas: &GasState) -> String {
    let (mechanics, _) = text.split_once("\nradiation ").unwrap();
    let areas = [[0.02, -0.01].as_slice(), [0.015].as_slice()];
    let omegas = [[1000.0, 2000.0].as_slice(), [1500.0].as_slice()];
    let mut out = String::new();
    let mut voice = 0;
    let mut mode = 0;
    for line in mechanics.lines() {
        if line.starts_with("voice ") {
            if line != "voice static-preload 2 2" { voice += 1; }
            mode = 0;
        }
        if let Some(rest) = line.strip_prefix("mode ") {
            let fields: Vec<&str> = rest.split_ascii_whitespace().collect();
            let omega: f64 = fields[0].parse().unwrap();
            assert_eq!(omega, omegas[voice][mode]);
            let h = compact_baffled_transfer(omega, areas[voice][mode], gas, 1.0);
            out.push_str(&format!("mode {} {} {:e} {:e} {} {}\n", fields[0], fields[1], h.re, h.im, fields[4], fields[5]));
            mode += 1;
        } else {
            out.push_str(line);
            out.push('\n');
        }
    }
    out
}

#[test]
fn the_radiation_block_is_exactly_the_computed_transfer_and_scales_with_density() {
    let air = gas(AIR);
    let derived = render(INPUT, 4801);
    let written = render(&authored(INPUT, &air), 4801);
    assert_eq!(derived.iter().map(|p| p.to_bits()).collect::<Vec<_>>(),
        written.iter().map(|p| p.to_bits()).collect::<Vec<_>>());
    let peak = derived.iter().fold(0.0_f64, |m, p| m.max(p.abs()));
    assert!(peak > 1e-4, "{peak}");
    // Only the gas changes; the mechanics ignore it, so p scales with rho at equal c.
    let thin_record = "ambient 293.15 70927.5 0";
    let thin = gas(thin_record);
    assert_eq!(thin.sound_speed.to_bits(), air.sound_speed.to_bits());
    let ratio = thin.density / air.density;
    let scaled = render(&INPUT.replace(AIR, thin_record), 4801);
    for (a, b) in derived.iter().zip(&scaled) { assert!((b - ratio * a).abs() <= 1e-12 * peak); }
    let info = ModalPerformance::from_bytes(INPUT.as_bytes(), 64).unwrap().info();
    assert_eq!(info.radiation.unwrap().gas, air);
    assert_eq!(info.radiation.unwrap().listener_m, 1.0);
}

#[test]
fn a_ringing_mode_is_heard_at_its_retarded_time_in_each_declared_gas() {
    // One undamped mode ringing from q = 0, qdot = V: q = (V/w) sin(w t), so the
    // compact baffled monopole is p(t) = -rho A w V sin(w (t - r/c)) / (2 pi r).
    let (w, v, area, r) = (2000.0, 0.01, 0.02, 3.0);
    for record in [AIR, "ambient 293.15 70927.5 0", "ambient mixture 293.15 101325 mole 1 CO2 1"] {
        let text = format!("frankensim-modal-performance-v1\nsample_rate_hz 48000\nsamples 4800\n\
            full_scale_pa 1\nlimits 0.9 1 100000 1e9 1e7\ncompile_limits 64 256\nvoices 1\n\
            voice retain-state 1 1\nmode {w} 0 0 0 0 {v}\nport 0 1\nevents 0\n\
            radiation compact-baffled\n{record}\nobserver {r}\nareas 1 {area}\n");
        let gas = gas(record);
        let actual = render(&text, 4800);
        let amplitude = gas.density * area * w * v / (core::f64::consts::TAU * r);
        for (n, p) in actual.iter().enumerate() {
            let t = (n + 1) as f64 / 48_000.0;
            let expected = -amplitude * (w * (t - r / gas.sound_speed)).sin();
            assert!((p - expected).abs() <= 1e-9 * amplitude, "{record} n={n}: {p} vs {expected}");
        }
    }
    // The CO2 lag differs from air's by w r (1/c_co2 - 1/c_air): about 4.7 rad here.
    let lag = 2000.0 * 3.0 * (1.0 / gas("ambient mixture 293.15 101325 mole 1 CO2 1").sound_speed
        - 1.0 / gas(AIR).sound_speed);
    assert!(lag > 1.0, "{lag}");
}

#[test]
fn a_radiation_block_admits_only_derived_transfers_and_a_resolved_gas() {
    for (old, new) in [
        (AIR, "ambient 293.15 101325"),
        (AIR, "ambient mixture 293.15 101325 mole 1 Xe 1"),
        (AIR, "ambient thermochemical 293.15 101325 mole 1 N2 1"),
        (AIR, "ambient 293.15 0 0"),
        ("mode 2000 0.04 0 0 0 0", "mode 2000 0.04 0.75 0 0 0"),
        ("areas 1 0.015", "areas 1 0.015\nareas 1 0.015"),
        ("\nareas 1 0.015", ""),
        ("areas 2 0.02 -0.01", "areas 1 0.02"),
        ("areas 1 0.015", "areas 1 NaN"),
        ("observer 1.0", "observer 0"),
        ("radiation compact-baffled", "radiation guessed"),
    ] {
        assert!(INPUT.contains(old));
        let text = INPUT.replace(old, new);
        assert!(ModalPerformance::from_bytes(text.as_bytes(), 64).is_err(), "{old} -> {new}");
    }
    let error = ModalPerformance::from_bytes(INPUT.replace(AIR, "ambient mixture 293.15 101325 mole 1 Xe 1").as_bytes(), 64)
        .err().unwrap();
    assert!(error.to_string().contains("line 22: ambient gas refused"), "{error}");
}

#[test]
fn the_command_reports_the_gas_and_refuses_an_unresolved_one_without_artifacts() {
    let dir = directory();
    let input = dir.join("radiating.performance");
    std::fs::write(&input, INPUT).unwrap();
    let output = dir.join("radiating.wav");
    let result = invoke(&input, &output);
    assert!(result.status.success(), "{}", String::from_utf8_lossy(&result.stdout));
    let sidecar = std::fs::read_to_string(output.with_extension("provenance.json")).unwrap();
    assert!(sidecar.contains(&format!("\"sound_speed_m_s\":{:e}", gas(AIR).sound_speed)), "{sidecar}");
    assert!(sidecar.contains("\"model\":\"compact-baffled narrow-band\""));
    let bad = dir.join("bad.performance");
    std::fs::write(&bad, INPUT.replace(AIR, "ambient mixture 293.15 101325 mole 1 Xe 1")).unwrap();
    let refused = dir.join("bad.wav");
    let result = invoke(&bad, &refused);
    assert!(!result.status.success());
    assert!(String::from_utf8_lossy(&result.stdout).contains("ambient gas refused"));
    assert!(!refused.exists() && !refused.with_extension("provenance.json").exists());
}
