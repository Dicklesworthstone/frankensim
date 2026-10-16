//! A TMM driving-point reflectance is a [`DelayedFilter`]: the same
//! characteristic port as a muffler, an HVAC run, or a pulse tube.

use fs_duct::{Duct, LossModel, Segment, Termination, impedance_sweep};
use fs_material::gas::{GasSpec, GasState};
use fs_math::c64::C64;
use fs_vfit::FitOptions;
use fs_vfit::discretize::{DelayedFilter, reflectance};

fn air() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).expect("air")
}

#[test]
fn open_cylinder_reflectance_returns_an_inverted_pulse() {
    let length = 0.34;
    let radius = 7.0e-3;
    let duct = Duct {
        segments: vec![Segment::Cylinder { radius, length }],
    };
    let gas = air();
    let dt = 1.0 / 16_000.0;
    let delay_samples = 2.0 * length / gas.sound_speed / dt;
    let area = core::f64::consts::PI * radius * radius;
    let zc = gas.density * gas.sound_speed / area;
    let omega0 = core::f64::consts::PI * gas.sound_speed / (2.0 * length);
    let sweep = impedance_sweep(
        &duct,
        &gas,
        0.25 * omega0,
        8.0 * omega0,
        80,
        LossModel::AllRegime,
        Termination::UnflangedOpen,
    )
    .expect("sweep");
    let omega: Vec<f64> = sweep.iter().map(|r| r.omega).collect();
    let h: Vec<C64> = sweep
        .iter()
        .map(|r| {
            let rac = reflectance(r.impedance, zc);
            rac.conj()
        })
        .collect();
    let mut opts = FitOptions::new(4);
    opts.fit_e = false;
    opts.iterations = 10;
    let mut line = DelayedFilter::from_tabulated(&omega, &h, delay_samples, dt, &opts, omega0)
        .expect("realize");
    let mut hist = Vec::new();
    for k in 0..80 {
        hist.push(line.push(if k == 0 { 1.0 } else { 0.0 }).expect("step"));
    }
    let k_ret = delay_samples.round() as usize;
    let peak = hist[k_ret.saturating_sub(2)..=(k_ret + 2).min(hist.len() - 1)]
        .iter()
        .copied()
        .fold(0.0_f64, |a, v| if v.abs() > a.abs() { v } else { a });
    assert!(
        peak < -0.15,
        "open-end characteristic return should invert, got {peak} near sample {k_ret}"
    );
    assert!(peak > -1.05, "passive |R| <= 1, got {peak}");
}

#[test]
fn scattering_passivity_caps_an_active_residual() {
    let filter = fs_vfit::discretize::DigitalFilter {
        sections: Vec::new(),
        direct: 1.25,
        t_s: 1.0 / 16_000.0,
        prewarp: 0.0,
    };
    let mut line = DelayedFilter::new(8.0, filter).expect("line");
    line.enforce_scattering_passivity(&[200.0, 800.0, 2_000.0]);
    let mut y = 0.0;
    for _ in 0..12 {
        y = line.push(1.0).expect("step");
    }
    assert!(
        y.abs() <= 1.0 + 1.0e-12,
        "passive residual cannot amplify, got {y}"
    );
}
