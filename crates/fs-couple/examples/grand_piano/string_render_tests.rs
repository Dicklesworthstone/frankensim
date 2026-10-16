//! Public render preparation, retained material admission, stereo and score time.
use super::*;
use std::{io::Write, path::PathBuf, sync::atomic::{AtomicU64, Ordering}};

fn options(args: &[&str]) -> Result<Options, String> {
    Options::parse(&args.iter().map(|s| (*s).to_owned()).collect::<Vec<_>>())
}
fn fresh_path(extension: &str) -> PathBuf {
    static SERIAL: AtomicU64 = AtomicU64::new(0);
    let time = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos();
    std::env::temp_dir().join(format!("fs-piano-string-{}-{time}-{}.{}",
        std::process::id(), SERIAL.fetch_add(1, Ordering::Relaxed), extension))
}
fn material_file(slope: f64) -> String {
    let path = fresh_path("fsps");
    let mut file = std::fs::OpenOptions::new().write(true).create_new(true).open(&path).unwrap();
    write!(file, "{}\nstretch,69,150000,{slope}\n", linear::string_stretching::HEADER).unwrap();
    path.to_str().unwrap().to_owned()
}

#[test]
fn admitted_material_survives_input_changes_without_a_second_read_or_linear_fallback() {
    let path = material_file(0.2);
    let mut o = options(&["--preset", "steinway-d", "--render", "p.wav",
        "--string-stretching", &path, "--dampers", "estimated"]).unwrap();
    o.modes = 12;
    let c = selected_scale(None, &o).unwrap()[48];
    let admitted = load_string_stretching(&[c], &o).unwrap();
    // Mutate only this test's own fresh input. A prepared render owns the value,
    // not a path whose contents can replace physical parameters during baking.
    std::fs::write(&path, "invalid replacement material").unwrap();
    let modes = board::demonstration();
    let mut piano = prepare_instrument_with_string_material(vec![c], &modes, &o, admitted.as_ref()).unwrap();
    assert!(piano.bank.has_string_stretching());
    assert!(piano.damper_resolution().is_some());
    assert!(prepare_instrument(vec![c], &modes, &o).is_err());
    assert!(prepare_instrument_with_string_material(vec![c], &modes, &o, None).is_err());
    piano.jack_on(69, 70., 0.007).unwrap();
    let mut extension = 0.0_f64;
    for _ in 0..2400 {
        piano.step().unwrap();
        extension = extension.max(piano.bank.string_stretching_observation(0).unwrap().stretching_energy_j);
    }
    assert!(extension > 0. && piano.accounting.felt_loss_j > 0. && piano.accounting.shank_loss_j > 0.);
    assert!((piano.accounting.input_work_j - piano.energy_j() - piano.accounting.dissipated_j()).abs() < 1e-7);
}

#[test]
fn physical_slope_refusal_does_not_publish_a_partial_or_linear_substitute_wav() {
    let path = material_file(1e-10); let output = fresh_path("wav");
    let o = options(&["--render", output.to_str().unwrap(), "--string-stretching", &path,
        "--note", "69", "--velocity", "2", "--duration", "0.03", "--modes", "12"]).unwrap();
    let c = geometry::demonstration_scale().unwrap()[48];
    let failure = render(output.to_str().unwrap(), vec![c], &board::demonstration(), None, &o).unwrap_err();
    assert!(failure.contains("slope"), "unexpected refusal: {failure}");
    assert!(!output.exists(), "a rejected physical trajectory must not publish a partial WAV");
}
