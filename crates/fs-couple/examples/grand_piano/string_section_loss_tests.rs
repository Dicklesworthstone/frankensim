//! Material substitution must reach the existing reciprocal time-domain owner.
use super::*;
use super::super::linear::{Bank, BoardMode};
use std::f64::consts::{PI, TAU};

fn template() -> Course {
    Course { length_m: 1.0, tension_n: 400.0, unison: 1,
        duplex_length_m: 0.0, detune_cents: 0.0,
        ..super::super::geometry::demonstration_scale().unwrap()[48] }
}

fn text(eta: f64, lo: f64, hi: f64) -> String {
    format!("{LOSS_HEADER}\n{LOSS_ROW}\nlayer,69,0.001,8000,80000000000,{eta},{lo},{hi},synthetic loss fixture,CC0-1.0\n")
}

fn course(eta: f64) -> Course {
    apply(&text(eta, 0.0, 100_000.0), &[template()]).unwrap()[0]
}

fn board(bridge: f64) -> [BoardMode; 1] {
    [BoardMode { frequency_hz: 100.0, damping_ratio: 0.0,
        bridge: [bridge; 88], volume: 0.0 }]
}

fn bank(course: Course, bridge: f64, modes: usize) -> Bank {
    Bank::new(&[course], &board(bridge), 48_000, 3_000.0, modes, true).unwrap()
}

fn near(actual: f64, expected: f64) {
    assert!((actual - expected).abs() <= 2e-12 * expected.abs().max(1e-12),
        "{actual} != {expected}");
}

#[test]
fn layered_viscosity_drives_every_member_span_and_transverse_direction() {
    let base = Course { unison: 3, duplex_length_m: 0.25, detune_cents: 2.0, ..template() };
    let authored = format!("{LOSS_HEADER}\n{LOSS_ROW}\n\
        layer,69,0.0006,8000,80000000000,100000000,0,100000,core fixture,CC0-1.0\n\
        layer,69,0.001,10000,120000000000,400000000,10,80000,outer fixture,CC0-1.0\n");
    let bound = apply(&authored, &[base]).unwrap()[0];
    let i_core = PI * 0.0006_f64.powi(4) / 4.0;
    let i_outer = PI * (0.001_f64.powi(4) - 0.0006_f64.powi(4)) / 4.0;
    let cb = 1e8 * i_core + 4e8 * i_outer;
    let loss = bound.bonded_section.unwrap().bending_loss().unwrap();
    near(loss.viscous_rigidity_nm2s(), cb);
    assert_eq!(loss.omega_band_rad_s(), (10.0, 80_000.0));
    near(bound.bonded_section.unwrap().axial_rigidity_n(),
        PI * (80e9 * 0.0006_f64.powi(2)
            + 120e9 * (0.001_f64.powi(2) - 0.0006_f64.powi(2))));
    let mechanics = Bank::new_with_transverse_bridge(&[bound], &board(0.03),
        48_000, 3_000.0, 4, true, Some(&[vec![0.02]])).unwrap();
    assert_eq!(mechanics.strings.len(), 12);
    let image = mechanics.gpu_image().unwrap();
    let rates = image.string_loss.unwrap();
    for string in &mechanics.strings {
        let length = if string.duplex { bound.duplex_length_m } else { bound.length_m };
        assert_eq!(string.modes.len(), 4);
        for (i, index) in string.modes.clone().enumerate() {
            let k = (i + 1) as f64 * PI / length;
            near(rates.strain_rate[index], cb * k.powi(4) / bound.linear_density_kg_m);
            // Existing separately authored velocity drag is not represented
            // as an internal material law, nor counted twice.
            assert_eq!(rates.mass_rate[index], 0.60);
        }
    }
}

#[test]
fn viscosity_only_substitution_changes_actual_decay_and_closes_energy() {
    let soft_loss = course(1e8);
    let strong_loss = course(8e8);
    assert_eq!(soft_loss.linear_density_kg_m, strong_loss.linear_density_kg_m);
    assert_eq!(soft_loss.flexural_rigidity_nm2, strong_loss.flexural_rigidity_nm2);
    assert_eq!(soft_loss.tension_n, strong_loss.tension_n);
    assert_eq!(soft_loss.partial_hz(1, soft_loss.tension_n),
        strong_loss.partial_hz(1, strong_loss.tension_n));
    let mut endpoints = Vec::new();
    for material in [soft_loss, strong_loss] {
        let mut mechanics = bank(material, 0.0, 3);
        mechanics.q[0] = 1e-4 * material.modal_mass_kg().sqrt();
        let initial = mechanics.energy();
        let mut deposited = 0.0;
        for _ in 0..1024 {
            mechanics.predict();
            mechanics.finish(&[0.0]);
            assert!(mechanics.last_modal_loss_j >= -1e-15 * initial);
            deposited += mechanics.last_modal_loss_j;
            mechanics.commit();
        }
        assert!((mechanics.energy() + deposited - initial).abs() < 1e-10 * initial);
        assert!(deposited > 0.0);
        endpoints.push((mechanics.energy(), mechanics.q.clone(), mechanics.v.clone()));
    }
    assert!(endpoints[1].0 < endpoints[0].0);
    assert_ne!(endpoints[0].1, endpoints[1].1);
    assert_ne!(endpoints[0].2, endpoints[1].2);
}

#[test]
fn material_bending_loss_reacts_on_the_same_moving_bridge() {
    let material = course(4e8);
    let mut mechanics = bank(material, 0.05, 3);
    let n = mechanics.modes.len();
    mechanics.v[0] = 0.01;
    mechanics.v[n] = -0.02;
    let initial = mechanics.energy();
    let mut deposited = 0.0;
    for _ in 0..64 {
        mechanics.predict();
        mechanics.finish(&[0.0]);
        deposited += mechanics.last_modal_loss_j;
        mechanics.commit();
    }
    assert!(deposited > 0.0);
    assert!((mechanics.energy() + deposited - initial).abs() < 1e-10 * initial);
    assert_ne!(mechanics.v[n], -0.02);
}

#[test]
fn actual_retained_frequencies_are_checked_including_duplex_and_detuning() {
    let base = template();
    let fundamental = TAU * course(1e8).partial_hz(1, base.tension_n);
    let limited = apply(&text(1e8, 0.0, fundamental * 1.1), &[base]).unwrap()[0];
    assert!(Bank::new(&[limited], &board(0.0), 48_000, 3_000.0, 1, true).is_ok());
    assert!(Bank::new(&[limited], &board(0.0), 48_000, 3_000.0, 2, true).is_err());
    let duplex = Course { duplex_length_m: 0.25, ..limited };
    assert!(Bank::new(&[duplex], &board(0.0), 48_000, 3_000.0, 1, true).is_err());

    let detuned_base = Course { unison: 3, detune_cents: 20.0, ..base };
    let tight = apply(&text(1e8, fundamental * 0.9, fundamental * 1.005),
        &[detuned_base]).unwrap()[0];
    assert!(Bank::new(&[tight], &board(0.0), 48_000, 3_000.0, 1, true).is_err());

    // A disabled dissipative operator does not evaluate a material loss law.
    assert!(Bank::new(&[limited], &board(0.0), 48_000, 3_000.0, 2, false).is_ok());
    let above = apply(&text(0.0, fundamental * 1.1, 100_000.0), &[base]).unwrap()[0];
    assert!(Bank::new(&[above], &board(0.0), 48_000, 3_000.0, 1, true).is_err());
}

#[test]
fn zero_is_explicit_missing_loss_refuses_and_replacement_clears_old_binding() {
    let source = [template()];
    for invalid in [text(-1.0, 0.0, 100_000.0), text(f64::NAN, 0.0, 100_000.0),
        text(1e8, 10.0, 1.0), text(1e8, 0.0, f64::INFINITY),
        format!("{LOSS_HEADER}\n{LOSS_ROW}\nlayer,69,0.001,8000,80000000000,synthetic,CC0-1.0\n"),
        format!("{}layer,69,0.002,8000,80000000000,1,200000,300000,synthetic,CC0-1.0\n",
            text(1e8, 0.0, 100_000.0))] {
        assert!(apply(&invalid, &source).is_err());
    }
    let zero = course(0.0);
    let image = bank(zero, 0.0, 3).gpu_image().unwrap().string_loss.unwrap();
    assert!(image.strain_rate.iter().all(|rate| *rate == 0.0));
    assert!(image.mass_rate.iter().all(|rate| *rate == 0.60));
    let legacy = format!("{HEADER}\n{ROW}\nlayer,69,0.001,8000,80000000000,synthetic,CC0-1.0\n");
    let replaced = apply(&legacy, &[zero]).unwrap()[0];
    assert!(replaced.bonded_section.unwrap().bending_loss().is_none());
    let image = bank(replaced, 0.0, 3).gpu_image().unwrap().string_loss.unwrap();
    assert!(image.strain_rate.iter().all(|rate| *rate > 0.0));
    assert!(Bank::new_with_string_damping(&[zero], &board(0.0),
        48_000, 3_000.0, 3, true, None, true).is_err());
}
