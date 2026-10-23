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
    let thin = pressure(Some("293.15 70927.5 0\n"), &directory, "thin");
    let co2 = pressure(Some("mixture 293.15 101325 mole 1 CO2 1\n"), &directory, "co2");
    let reference = ambient::Ambient::reference().gas;
    let dense = ambient::parse("mixture 293.15 101325 mole 1 CO2 1").unwrap();
    let ratio = ambient::parse("293.15 70927.5 0").unwrap().density / reference.density;
    let peak = air.iter().fold(0.0_f64, |m, p| m.max(p.abs()));
    assert!(peak > 0.0);
    // No radiation backreaction: the motion ignores the air, and the Rayleigh
    // pressure is rho-linear with identical delays at equal c.
    for (a, t) in air.iter().zip(&thin) { assert!((t - ratio * a).abs() <= 1e-12 * peak, "{a} {t}"); }
    // The CO2 trace lags by the path times (1/c_co2 - 1/c_air); the two surface
    // samples lie 0.82 m and 1.04 m from the microphone.
    let per_metre = (1.0 / dense.sound_speed - 1.0 / reference.sound_speed) * 48_000.0;
    let correlation = |lag: usize| air.iter().zip(&co2[lag..]).map(|(a, c)| a * c).sum::<f64>();
    let lag = (0..200).max_by(|a, b| correlation(*a).total_cmp(&correlation(*b))).unwrap() as f64;
    let near = (0.15_f64.powi(2) + 0.1_f64.powi(2) + 0.8_f64.powi(2)).sqrt();
    let far = (0.65_f64.powi(2) + 0.1_f64.powi(2) + 0.8_f64.powi(2)).sqrt();
    assert!(lag >= near * per_metre - 1.0 && lag <= far * per_metre + 1.0,
        "lag {lag} outside [{}, {}]", near * per_metre, far * per_metre);
    // A record that does not resolve refuses before any output is written.
    let bad = directory.join("bad.gas");
    std::fs::write(&bad, "293.15 101325\n").unwrap();
    assert!(ambient::Ambient::selected(bad.to_str()).is_err());
}
