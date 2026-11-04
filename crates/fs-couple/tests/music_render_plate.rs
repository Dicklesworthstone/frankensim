//! The actual command consumes geometry/materials, renders and hashes a WAV.
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};
use fs_couple::pcm_wav::encode_pcm16_wav;
use fs_couple::render::plate::file::PlatePerformance;

const SOURCE: &[u8] = include_bytes!("../examples/plate-mesh.performance");
static NEXT: AtomicU64 = AtomicU64::new(0);
fn scratch() -> PathBuf {
    let stamp = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos();
    let path = std::env::temp_dir().join(format!("fs-music-plate-{}-{stamp}-{}",
        std::process::id(), NEXT.fetch_add(1, Ordering::Relaxed)));
    std::fs::create_dir(&path).unwrap();
    path // Retained on failure/success for diagnostics; existing paths never removed.
}
fn run(input: &Path, output: &Path, block: usize) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_music_render"))
        .arg("plate").arg(input).arg(output).arg("--block").arg(block.to_string())
        .output().unwrap()
}
fn rendered(input: &Path, output: &Path, block: usize) -> Vec<u8> {
    let result = run(input, output, block);
    assert!(result.status.success(), "stdout={} stderr={}",
        String::from_utf8_lossy(&result.stdout), String::from_utf8_lossy(&result.stderr));
    assert!(String::from_utf8_lossy(&result.stdout).contains("\"verdict\":\"rendered\""));
    std::fs::read(output).unwrap()
}

#[test]
fn command_matches_public_physics_and_pcm_across_partitions_and_relocation() {
    let root = scratch();
    let source = root.join("mesh.performance");
    let moved = root.join("same-mesh-elsewhere.performance");
    std::fs::write(&source, SOURCE).unwrap();
    std::fs::write(&moved, SOURCE).unwrap();
    let mut performance = PlatePerformance::from_bytes(SOURCE, 4801).unwrap().into_renderer();
    let mut pressure = vec![0.0; 4801];
    performance.block(&mut pressure).unwrap();
    let (expected, clips) = encode_pcm16_wav(&pressure, 48_000, 1.0).unwrap();
    assert!(expected[44..].iter().any(|&b| b!=0));
    let hash = fs_blake3::hash_domain("org.frankensim.fs-couple.music-render-wav.v1", &expected).to_hex();
    let a = root.join("a.wav");
    let b = root.join("b.wav");
    assert_eq!(rendered(&source, &a, 37), expected);
    assert_eq!(rendered(&moved, &b, 571), expected);
    assert_eq!(expected.len(), 44+4801*2);
    for output in [&a,&b] {
        let sidecar = std::fs::read_to_string(output.with_extension("provenance.json")).unwrap();
        for fragment in ["\"nodes\":30", "\"triangles\":40", "\"sections\":2",
            "\"requested_modes\":3", "\"force_events\":3", "\"samples\":4801"] {
            assert!(sidecar.contains(fragment));
        }
        assert!(sidecar.contains(&format!("\"clipped_samples\":{clips}")));
        assert!(sidecar.contains(&hash));
        assert!(!sidecar.contains(&source.display().to_string()));
    }
    let previous = std::fs::read(&a).unwrap();
    assert!(!run(&source, &a, 64).status.success());
    assert_eq!(std::fs::read(&a).unwrap(), previous);
}

#[test]
fn modifying_a_constitutive_region_changes_the_produced_sound() {
    let root = scratch();
    let a = root.join("first.performance");
    let b = root.join("stiffer.performance");
    std::fs::write(&a, SOURCE).unwrap();
    let changed = std::str::from_utf8(SOURCE).unwrap().replace("11000000000.0", "22000000000.0");
    std::fs::write(&b, changed).unwrap();
    let first = rendered(&a, &root.join("first.wav"), 64);
    let second = rendered(&b, &root.join("stiffer.wav"), 64);
    assert_ne!(&first[44..], &second[44..], "material assignment must reach the audio, not only metadata");
}

#[test]
fn explicit_material_gas_reaction_reaches_the_command_and_sidecar() {
    let root = scratch();
    let text = std::str::from_utf8(SOURCE).unwrap().replace("observer 1.0\n",
        "observer 1.0\nradiation-load baffled-piston\n");
    let input = root.join("loaded.performance");
    let output = root.join("loaded.wav");
    std::fs::write(&input, &text).unwrap();
    let loaded = PlatePerformance::from_bytes(text.as_bytes(), 4801).unwrap();
    assert!(loaded.info().baffled_piston_loading);
    let mut performance = loaded.into_renderer();
    let mut pressure = vec![0.0; 4801];
    performance.block(&mut pressure).unwrap();
    let mut unloaded = vec![0.0; 4801];
    PlatePerformance::from_bytes(SOURCE, 4801).unwrap().into_renderer().block(&mut unloaded).unwrap();
    assert!(pressure.iter().zip(&unloaded).any(|(a,b)| (a-b).abs() > 1e-7),
        "the declared reaction must change the physical signal");
    let (expected, _) = encode_pcm16_wav(&pressure, 48_000, 1.0).unwrap();
    assert_eq!(rendered(&input, &output, 37), expected);
    let sidecar = std::fs::read_to_string(output.with_extension("provenance.json")).unwrap();
    assert!(sidecar.contains("\"radiation_load\":\"baffled-piston\""));
    assert!(!sidecar.contains("no radiation loading"));
}

#[test]
fn data_only_pressure_change_alters_loaded_motion_beyond_pressure_scaling() {
    let text = std::str::from_utf8(SOURCE).unwrap();
    let render = |loaded: bool, doubled_pressure: bool| {
        let mut source = text.to_string();
        if loaded { source = source.replace("observer 1.0\n", "observer 1.0\nradiation-load baffled-piston\n"); }
        if doubled_pressure { source = source.replace("ambient 293.15 101325 0", "ambient 293.15 202650 0"); }
        pressure(&source)
    };
    for loaded in [false, true] {
        let (one, air) = render(loaded, false);
        let (two, dense) = render(loaded, true);
        assert_eq!(air.sound_speed, dense.sound_speed);
        assert_eq!(dense.density, 2.0*air.density);
        let scale = one.iter().map(|p| p*p).sum::<f64>();
        let defect = one.iter().zip(two).map(|(a,b)| (b-2.0*a).powi(2)).sum::<f64>();
        assert!(scale > 1e-12);
        if loaded { assert!(defect/scale > 1e-8, "gas did not reach mechanical reaction: {defect}/{scale}"); }
        else { assert!(defect/scale < 1e-24, "unloaded same-c reference must scale only with density"); }
    }
}

#[test]
fn radiation_selection_is_explicit_and_unknown_or_duplicate_loads_refuse() {
    let text = std::str::from_utf8(SOURCE).unwrap();
    let source = |row: &str| text.replace("observer 1.0\n", &format!("observer 1.0\n{row}\n"));
    assert!(!PlatePerformance::from_bytes(source("radiation-load none").as_bytes(), 64)
        .unwrap().info().baffled_piston_loading);
    for row in ["radiation-load", "radiation-load automatic", "radiation-load baffled-piston extra",
        "radiation-load none\nradiation-load baffled-piston"] {
        assert!(PlatePerformance::from_bytes(source(row).as_bytes(), 64).is_err(), "{row}");
    }
}

#[test]
fn invalid_geometry_and_conflicting_options_leave_no_output_artifact() {
    let root = scratch();
    let input = root.join("invalid.performance");
    let out = root.join("invalid.wav");
    let invalid = std::str::from_utf8(SOURCE).unwrap().replace("triangle 0 1 7 0", "triangle 0 1 7 50");
    std::fs::write(&input, invalid).unwrap();
    assert!(!run(&input, &out, 37).status.success());
    assert!(!out.exists());
    assert!(!out.with_extension("provenance.json").exists());
    std::fs::write(&input, SOURCE).unwrap();
    let result = Command::new(env!("CARGO_BIN_EXE_music_render"))
        .arg("plate").arg(&input).arg(&out).args(["--seconds", "1"]).output().unwrap();
    assert!(!result.status.success());
    assert!(!out.exists());
}

fn with_ambient(record: &str) -> String {
    std::str::from_utf8(SOURCE).unwrap().replace("ambient 293.15 101325 0", record)
}
fn pressure(text: &str) -> (Vec<f64>, fs_material::gas::GasState) {
    let performance = PlatePerformance::from_bytes(text.as_bytes(), 4801).unwrap();
    let gas = performance.info().gas;
    let mut renderer = performance.into_renderer();
    let mut out = vec![0.0; 4801];
    renderer.block(&mut out).unwrap();
    (out, gas)
}

#[test]
fn changing_only_the_gas_record_rescales_and_retards_the_sound_physically() {
    let air = with_ambient("ambient 293.15 101325 0");
    let thin = with_ambient("ambient 293.15 70927.5 0");
    let co2 = with_ambient("ambient mixture 293.15 101325 mole 1 CO2 1");
    let heliox = with_ambient("ambient mixture 293.15 101325 mole 2 He 0.8 O2 0.2");
    let (pa, ga) = pressure(&air);
    let (pt, gt) = pressure(&thin);
    let (pc, gc) = pressure(&co2);
    let (ph, gh) = pressure(&heliox);
    // Force-driven plate, no reaction load: at equal c the far field scales with rho.
    assert_eq!(ga.sound_speed.to_bits(), gt.sound_speed.to_bits());
    let ratio = gt.density / ga.density;
    assert!((ratio - 0.7).abs() < 1e-12, "{ratio}");
    let peak = pa.iter().fold(0.0_f64, |m, p| m.max(p.abs()));
    assert!(peak > 1e-6);
    for (a, t) in pa.iter().zip(&pt) { assert!((t - ratio * a).abs() <= 1e-12 * peak); }
    // The first wavefront reaches the 1 m listener r/c after the drive starts.
    for (out, gas) in [(&pa, ga), (&pt, gt), (&pc, gc), (&ph, gh)] {
        let onset = (48_000.0 / gas.sound_speed).floor() as usize;
        assert!(out[..onset].iter().all(|p| *p == 0.0));
        assert!(out[onset] != 0.0);
    }
    assert!(gc.sound_speed < 0.8 * ga.sound_speed && gc.density > 1.4 * ga.density);
    // 80/20 He-O2 is light and fast: the same motion, heard earlier and quieter
    // by rho. Over equal windows after each arrival, energy scales with rho^2.
    assert!(gh.sound_speed > 1.7 * ga.sound_speed && gh.density < 0.35 * ga.density);
    let window = |out: &[f64], gas: fs_material::gas::GasState| {
        let onset = (48_000.0 / gas.sound_speed).floor() as usize;
        out[onset..onset + 4000].iter().map(|p| p * p).sum::<f64>()
    };
    let energy = window(&ph, gh) / window(&pa, ga);
    let expected = (gh.density / ga.density).powi(2);
    assert!((energy / expected - 1.0).abs() < 1e-2, "{energy} vs {expected}");
    // The command consumes the same edited file.
    let root = scratch();
    for (name, text, gas) in [("air", &air, ga), ("co2", &co2, gc)] {
        let input = root.join(format!("{name}.performance"));
        std::fs::write(&input, text).unwrap();
        let output = root.join(format!("{name}.wav"));
        let wav = rendered(&input, &output, 64);
        let onset = (48_000.0 / gas.sound_speed).floor() as usize;
        let pcm = |i: usize| i16::from_le_bytes([wav[44 + 2 * i], wav[45 + 2 * i]]);
        assert!((0..onset).all(|i| pcm(i) == 0), "{name}: silence before r/c in the WAV");
        let sidecar = std::fs::read_to_string(output.with_extension("provenance.json")).unwrap();
        assert!(sidecar.contains(&format!("\"sound_speed_m_s\":{:e}", gas.sound_speed)), "{sidecar}");
        assert!(sidecar.contains(&format!("\"propagation_delay_s\":{:e}", 1.0 / gas.sound_speed)));
    }
}

#[test]
fn a_gas_record_that_does_not_resolve_refuses_by_name_and_writes_nothing() {
    let root = scratch();
    for record in ["ambient 293.15 101325", "ambient mixture 293.15 101325 mole 1 Xe 1",
        "ambient thermochemical 293.15 101325 mole 1 N2 1"] {
        let input = root.join("bad-gas.performance");
        std::fs::write(&input, with_ambient(record)).unwrap();
        let out = root.join("bad-gas.wav");
        let result = run(&input, &out, 37);
        assert!(!result.status.success(), "{record}");
        let stdout = String::from_utf8_lossy(&result.stdout);
        assert!(stdout.contains("\"verdict\":\"refused\"") && stdout.contains("line 3"), "{record}: {stdout}");
        assert!(!out.exists());
        assert!(!out.with_extension("provenance.json").exists());
    }
    let input = root.join("record.performance");
    std::fs::write(&input, SOURCE).unwrap();
    let out = root.join("record.wav");
    let result = Command::new(env!("CARGO_BIN_EXE_music_render"))
        .arg("plate").arg(&input).arg(&out).arg("--gas-record").arg(&input).output().unwrap();
    assert!(!result.status.success());
    assert!(!out.exists());
}
