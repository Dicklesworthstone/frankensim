//! The bowed listening command takes its surrounding gas from a record file.
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};
use fs_material::gas::GasState;
use fs_material::gas_species::input::read_gas_state;

static NEXT: AtomicU64 = AtomicU64::new(0);
fn scratch() -> PathBuf {
    let stamp = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos();
    let path = std::env::temp_dir().join(format!("fs-bowed-listen-{}-{stamp}-{}",
        std::process::id(), NEXT.fetch_add(1, Ordering::Relaxed)));
    std::fs::create_dir(&path).unwrap();
    path // Retained for diagnostics; existing paths are never removed.
}
fn listen(dir: &Path, name: &str, record: &str) -> std::process::Output {
    let gas = dir.join(format!("{name}.gas"));
    std::fs::write(&gas, record).unwrap();
    Command::new(env!("CARGO_BIN_EXE_bowed_listen"))
        .arg("--gas-record").arg(&gas)
        .arg("--samples").arg("4800")
        .arg("--wav").arg(dir.join(format!("{name}.wav")))
        .arg("--receipt").arg(dir.join(format!("{name}.jsonl")))
        .output().unwrap()
}
fn field(receipt: &str, key: &str) -> f64 {
    let start = receipt.find(&format!("\"{key}\":")).unwrap_or_else(|| panic!("{key} in {receipt}")) + key.len() + 3;
    let end = start + receipt[start..].find([',', '}']).unwrap();
    receipt[start..end].parse().unwrap()
}
fn state(record: &str) -> GasState {
    read_gas_state(&mut record.split_ascii_whitespace(), None).unwrap()
}

#[test]
fn only_the_gas_record_changes_the_radiated_pressure_scale() {
    let dir = scratch();
    let air = "293.15 101325 0\n";
    let shipped = std::fs::read_to_string(concat!(env!("CARGO_MANIFEST_DIR"), "/examples/dry-air-20c.gas")).unwrap();
    assert_eq!(shipped, air, "the shipped record is the one measured here");
    let mut full_scale = Vec::new();
    for (name, record) in [("air", air), ("thin", "293.15 70927.5 0\n"),
        ("co2", "mixture 293.15 101325 mole 1 CO2 1\n"),
        ("heliox", "mixture 293.15 101325 mole 2 He 0.8 O2 0.2\n")] {
        let out = listen(&dir, name, record);
        assert!(out.status.success(), "{name}: {}", String::from_utf8_lossy(&out.stderr));
        let receipt = std::fs::read_to_string(dir.join(format!("{name}.jsonl"))).unwrap();
        let gas = state(record);
        assert_eq!(field(&receipt, "density_kg_m3"), gas.density);
        assert_eq!(field(&receipt, "sound_speed_m_s"), gas.sound_speed);
        full_scale.push((field(&receipt, "full_scale_pa"), gas.density));
    }
    // The one-way body motion does not depend on the gas, so the radiated
    // pressure, and hence its peak, scales with density alone.
    let (p_air, rho_air) = full_scale[0];
    for &(p, rho) in &full_scale[1..] {
        assert!((p / p_air - rho / rho_air).abs() < 1e-12 * (rho / rho_air), "{p} {p_air} {rho} {rho_air}");
    }
    assert!((full_scale[1].1 / rho_air - 0.7).abs() < 1e-12);
    assert!(full_scale[2].1 > 1.4 * rho_air);
    assert!(full_scale[3].1 < 0.35 * rho_air);
}

#[test]
fn a_missing_or_unresolvable_gas_record_refuses_and_writes_nothing() {
    let dir = scratch();
    for (name, record) in [("short", "293.15 101325\n"), ("xenon", "mixture 293.15 101325 mole 1 Xe 1\n"),
        ("trailing", "293.15 101325 0 7\n"), ("thermo", "thermochemical 293.15 101325 mole 1 N2 1\n")] {
        let out = listen(&dir, name, record);
        assert!(!out.status.success(), "{name}");
        assert!(String::from_utf8_lossy(&out.stderr).contains("gas record"), "{name}");
        assert!(!dir.join(format!("{name}.wav")).exists());
        assert!(!dir.join(format!("{name}.jsonl")).exists());
    }
    let out = Command::new(env!("CARGO_BIN_EXE_bowed_listen"))
        .arg("--wav").arg(dir.join("none.wav")).arg("--receipt").arg(dir.join("none.jsonl"))
        .output().unwrap();
    assert!(!out.status.success());
    assert!(String::from_utf8_lossy(&out.stderr).contains("--gas-record"));
    assert!(!dir.join("none.wav").exists());
}
