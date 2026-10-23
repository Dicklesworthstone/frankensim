//! Only the declared air changes; the radiated piano pressure follows rho and c.
use super::*;

fn pressure(gas_record: Option<&str>, directory: &std::path::Path, name: &str) -> Vec<f64> {
    let course = geometry::demonstration_scale().unwrap()[48];
    let modes = board::demonstration();
    let surface: Vec<_> = [-0.25, 0.25].into_iter().enumerate().map(|(side, x)| {
        board_geometry::SurfaceSample {
            position_m: [x, 0.0, 0.0],
            area_m2: 0.1,
            mode_shape: modes.iter().enumerate().map(|(i, mode)| {
                mode.volume / 0.2 * if (i + side) % 2 == 0 { 1.0 } else { -0.5 }
            }).collect(),
        }
    }).collect();
    let record = gas_record.map(|text| {
        let path = directory.join(format!("{name}.gas"));
        std::fs::write(&path, text).unwrap();
        path.to_str().unwrap().to_owned()
    });
    let wav = directory.join(format!("{name}.wav"));
    let csv = directory.join(format!("{name}.csv"));
    let options = Options {
        render: Some(wav.to_str().unwrap().to_owned()),
        note: Some(69), velocity: Some(2.0), duration: 1400.0 / 48_000.0,
        sample_rate: 48_000, modes: 8, pcm_full_scale_pa: 64.0,
        microphone: Some([-0.4, 0.1, 0.8]),
        receiver_pressure_csv: Some(csv.to_str().unwrap().to_owned()),
        gas_record: record,
        ..Options::default()
    };
    render(wav.to_str().unwrap(), vec![course], &modes, Some(&surface), &options).unwrap();
    std::fs::read_to_string(&csv).unwrap().lines().skip(1)
        .map(|row| row.split(',').nth(2).unwrap().parse::<f64>().unwrap()).collect()
}

#[test]
fn the_radiated_pressure_scales_with_density_and_arrives_with_the_sound_speed() {
    let stamp = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos();
    let directory = std::env::temp_dir().join(format!("frankensim-piano-ambient-{}-{stamp}", std::process::id()));
    std::fs::create_dir(&directory).unwrap();
    let air = pressure(None, &directory, "reference");
    let dry = pressure(Some("293.15 101325 0\n"), &directory, "dry");
    let thin = pressure(Some("293.15 70927.5 0\n"), &directory, "thin");
    let co2 = pressure(Some("mixture 293.15 101325 mole 1 CO2 1\n"), &directory, "co2");
    let heliox = pressure(Some("mixture 293.15 101325 mole 2 He 0.8 O2 0.2\n"), &directory, "heliox");
    let reference = ambient::Ambient::reference().gas;
    let dense = ambient::parse("mixture 293.15 101325 mole 1 CO2 1").unwrap();
    let still = ambient::parse("293.15 101325 0").unwrap();
    let ratio = ambient::parse("293.15 70927.5 0").unwrap().density / still.density;
    let peak = dry.iter().fold(0.0_f64, |m, p| m.max(p.abs()));
    assert!(peak > 0.0);
    // No radiation backreaction: the motion ignores the air, and the Rayleigh
    // pressure is rho-linear with identical delays at equal c. Equal c needs
    // equal composition, so the pair is explicitly dry at both pressures.
    for (a, t) in dry.iter().zip(&thin) { assert!((t - ratio * a).abs() <= 1e-12 * peak, "{a} {t}"); }
    // The built-in reference is humid (65 % RH): lighter and faster than dry air.
    assert!(reference.water_mole_fraction > 0.0 && reference.sound_speed > still.sound_speed);
    assert!(reference.density < still.density);
    assert!(air.iter().zip(&dry).any(|(h, d)| h != d), "the humid reference rendered as dry air");
    // The CO2 trace lags by the path times (1/c_co2 - 1/c_air); the two surface
    // samples lie 0.82 m and 1.04 m from the microphone.
    let per_metre = (1.0 / dense.sound_speed - 1.0 / reference.sound_speed) * 48_000.0;
    let correlation = |lag: usize| air.iter().zip(&co2[lag..]).map(|(a, c)| a * c).sum::<f64>();
    let lag = (0..200).max_by(|a, b| correlation(*a).total_cmp(&correlation(*b))).unwrap() as f64;
    let near = (0.15_f64.powi(2) + 0.1_f64.powi(2) + 0.8_f64.powi(2)).sqrt();
    let far = (0.65_f64.powi(2) + 0.1_f64.powi(2) + 0.8_f64.powi(2)).sqrt();
    assert!(lag >= near * per_metre - 1.0 && lag <= far * per_metre + 1.0,
        "lag {lag} outside [{}, {}]", near * per_metre, far * per_metre);
    // 80/20 He-O2 is faster: the air trace now lags the heliox one.
    let light = ambient::parse("mixture 293.15 101325 mole 2 He 0.8 O2 0.2").unwrap();
    let ahead = (1.0 / reference.sound_speed - 1.0 / light.sound_speed) * 48_000.0;
    let correlation = |lag: usize| heliox.iter().zip(&air[lag..]).map(|(h, a)| h * a).sum::<f64>();
    let lead = (0..200).max_by(|a, b| correlation(*a).total_cmp(&correlation(*b))).unwrap() as f64;
    assert!(lead >= near * ahead - 1.0 && lead <= far * ahead + 1.0,
        "heliox lead {lead} outside [{}, {}]", near * ahead, far * ahead);
    // A record that does not resolve refuses before any output is written.
    let bad = directory.join("bad.gas");
    std::fs::write(&bad, "293.15 101325\n").unwrap();
    assert!(ambient::Ambient::selected(bad.to_str()).is_err());
}

#[test]
fn a_thermochemical_record_resolves_only_through_its_supplied_source_packs() {
    // Synthetic packs (the music_render CLI's own fixture), not a real gas.
    let directory = std::env::temp_dir().join(format!("fs-piano-thermochemical-{}-{}",
        std::process::id(), std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos()));
    std::fs::create_dir(&directory).unwrap();
    let packs = thermochemical_support::packs("specimen-gas", 0.028, 3.5, 200.0);
    let [species, nasa9, transport] = ["species", "nasa9", "transport"].map(|name| directory.join(format!("{name}.pack")));
    std::fs::write(&species, packs.species.to_bytes()).unwrap();
    std::fs::write(&nasa9, packs.nasa9.to_bytes()).unwrap();
    std::fs::write(&transport, packs.transport.claims_pack().to_bytes()).unwrap();
    let record = directory.join("specimen.gas");
    std::fs::write(&record, "thermochemical 293.15 101325 mole 1 specimen-gas 1\n").unwrap();
    let args = |extra: &[&str]| -> Vec<String> {
        ["--gas-record", record.to_str().unwrap()].iter().chain(extra).map(|s| s.to_string()).collect()
    };
    // Without its source packs the record refuses by name; nothing is defaulted.
    let bare = Options::parse(&args(&[])).unwrap();
    let refused = selected_air(&bare).unwrap_err();
    assert!(refused.contains("explicit source registry"), "{refused}");
    // With them it is that gas: ideal-gas density p M / (R T) of the pack's species.
    let supplied = Options::parse(&args(&["--gas-component", species.to_str().unwrap(), nasa9.to_str().unwrap(),
        transport.to_str().unwrap(), "--gas-water", "none", "--gas-phase", "frozen"])).unwrap();
    let air = selected_air(&supplied).unwrap();
    assert_eq!(air.source, ambient::Source::CommandLine);
    let density = 101_325.0 * 0.028 / (fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * 293.15);
    assert!((air.gas.density / density - 1.0).abs() < 1e-12, "{} vs {density}", air.gas.density);
    assert!(air.receipt().contains("--gas-record"));
    // The packs without a water identity or phase refuse at the front door.
    assert!(Options::parse(&args(&["--gas-component", species.to_str().unwrap(), nasa9.to_str().unwrap(),
        transport.to_str().unwrap()])).is_err());
}
