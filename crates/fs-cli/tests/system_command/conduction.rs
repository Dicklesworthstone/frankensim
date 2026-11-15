//! `conduction_body`: fs-conduction's backward-Euler finite-element body,
//! lowered from a scene.
//!
//! The check is the one-dimensional closed form of fs-conduction
//! tests/conduction_pair.rs, from a scene. Two equal slabs,
//! examples/system/two_slabs.json, start at 400 K and 300 K with insulated
//! outer ends and exchange heat through a 0.5 W/K conductance (a contact of
//! 8e-4 m^2 K/W over the 20 x 20 mm end, Bi = 5). The claim, declared before
//! the first run: the energy that has crossed the wall by t* = 0.1 L^2 /
//! alpha converges to the closed form at order at least
//! [`CLAIMED_ORDER`] between 4 and 8 cells along each slab, with the step
//! t* / (2 n^2); the pair's stored energy and the ledger close; the wall's
//! last rate is its conductance times the two faces' reported temperatures.
//! Then the refusals.

use std::f64::consts::{FRAC_PI_2, PI};

use super::*;
use fs_cli::SystemScene;

const LENGTH_M: f64 = 0.03;
const CROSS_M: f64 = 0.02;
const RHO_C_J_M3_K: f64 = 3.6e6;
const CONDUCTANCE_W_K: f64 = 0.5;
/// Bi = 2 L / (k R''), R'' = A / G, for the half-problem about the wall.
const BIOT: f64 = 2.0 * LENGTH_M * CONDUCTANCE_W_K / (15.0 * CROSS_M * CROSS_M);
const HOT_K: f64 = 400.0;
const COLD_K: f64 = 300.0;
/// The comparison time in units of L^2 / alpha.
const FOURIER: f64 = 0.1;
/// The least order between the two rungs that meets the claim.
const CLAIMED_ORDER: f64 = 1.8;

/// The excess heat the hot slab holds over the mean at the start [J].
fn excess_j() -> f64 {
    RHO_C_J_M3_K * LENGTH_M * CROSS_M * CROSS_M * 0.5 * (HOT_K - COLD_K)
}

/// The energy that has crossed the wall by t*, from the closed form: each
/// slab is the plane wall insulated at its outer end with a Robin end of
/// `h = 2 G / A`, still holding `sum C_n sin(b_n) / b_n exp(-b_n^2 Fo)` of its
/// excess, `b_n tan b_n = Bi`, `C_n = 4 sin b_n / (2 b_n + sin 2 b_n)`.
fn closed_form_crossed_j() -> f64 {
    let mut held = 0.0;
    for n in 0..64 {
        // b tan b - Bi rises from -Bi to +infinity on (n pi, n pi + pi/2).
        let start = f64::from(n) * PI;
        let (mut lo, mut hi) = (start, start + FRAC_PI_2);
        loop {
            let mid = 0.5 * (lo + hi);
            if mid <= lo || mid >= hi {
                break;
            }
            if mid * mid.tan() < BIOT {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        let b = 0.5 * (lo + hi);
        let term =
            4.0 * b.sin() / (2.0 * b + (2.0 * b).sin()) * b.sin() / b * (-b * b * FOURIER).exp();
        held += term;
        if term.abs() < 1.0e-17 {
            break;
        }
    }
    excess_j() * (1.0 - held)
}

fn probe(scene: &SystemScene, path: &str) -> f64 {
    scene
        .probe(path)
        .unwrap()
        .unwrap_or_else(|| panic!("{path} is not reported"))
}

/// The example with `cells` along each slab and its step t* / (2 n^2).
fn two_slabs(cells: u32) -> String {
    let text = std::fs::read_to_string(example("two_slabs.json")).unwrap();
    assert_eq!(text.matches("\"cells\": [4, 2, 2]").count(), 2);
    assert!(text.contains("\"step_s\": 0.675,") && text.contains("\"duration_s\": 21.6,"));
    let step_s = 21.6 / f64::from(2 * cells * cells);
    text.replace(
        "\"cells\": [4, 2, 2]",
        &format!("\"cells\": [{cells}, 2, 2]"),
    )
    .replace("\"step_s\": 0.675,", &format!("\"step_s\": {step_s},"))
}

#[test]
fn two_conduction_bodies_from_a_scene_converge_to_the_closed_form_at_second_order() {
    let exact_j = closed_form_crossed_j();
    // The pair's stored energy, rho c V (T_a + T_b).
    let held_j = RHO_C_J_M3_K * LENGTH_M * CROSS_M * CROSS_M * (HOT_K + COLD_K);
    eprintln!("closed form: {exact_j:.15e} J crossed (Bi {BIOT})");
    let mut errors = Vec::new();
    for cells in [4, 8] {
        let text = two_slabs(cells);
        scene_context(|cx, _| {
            let mut scene = SystemScene::from_text(&text, "two_slabs.json").unwrap();
            let mut crossed_j = 0.0;
            let mut time_s = scene.time_s();
            while scene.advance(f64::MAX, cx).unwrap() {
                let now_s = scene.time_s();
                crossed_j += probe(&scene, "wall.heat_rate_w") * (now_s - time_s);
                time_s = now_s;
            }
            let (face_a, face_b, rate) = (
                probe(&scene, "a.x_max.temperature_k"),
                probe(&scene, "b.x_min.temperature_k"),
                probe(&scene, "wall.heat_rate_w"),
            );
            let wall = rate / (CONDUCTANCE_W_K * (face_a - face_b)) - 1.0;
            let ledger = scene.ledger();
            let (stored, residual) = (ledger.stored_change_j(), ledger.energy_residual_j());
            eprintln!(
                "n {cells}: {} steps to {time_s} s; crossed {crossed_j:.15e} J; faces {face_a} K, {face_b} K, wall rate / G dT - 1 = {wall:e}; a {} .. {} K, b {} .. {} K; stored change {stored:e} J, ledger residual {residual:e} J",
                scene.steps(),
                probe(&scene, "a.min_temperature_k"),
                probe(&scene, "a.max_temperature_k"),
                probe(&scene, "b.min_temperature_k"),
                probe(&scene, "b.max_temperature_k"),
            );
            assert!(
                (time_s / 21.6 - 1.0).abs() <= 1e-12,
                "the run ends at t*: {time_s}"
            );
            // The wall reads the faces the bodies return on the settled load
            // line.
            assert!(wall.abs() <= 1e-9, "wall rate off its faces by {wall:e}");
            // The pair is closed: nothing enters, so what one slab gives the
            // other takes, up to the owner's linear solves.
            assert!(
                stored.abs() <= 1e-9 * held_j,
                "the pair changed by {stored:e} J"
            );
            assert!(
                residual.abs() <= 1e-9 * held_j,
                "ledger residual {residual:e} J"
            );
            errors.push((crossed_j - exact_j).abs() / excess_j());
        });
    }
    let order = (errors[0] / errors[1]).log2();
    eprintln!("relative errors {errors:?}; order {order:.3}");
    assert!(
        order >= CLAIMED_ORDER,
        "order {order:.3} below {CLAIMED_ORDER}"
    );
}

#[test]
fn the_two_slab_example_runs_through_the_binary_and_reports_its_faces() {
    let path = example("two_slabs.json");
    let (code, json, stderr) = run_binary(&["--json", "system", &path]);
    assert_eq!(code, 0, "{stderr}");
    let (face_a, face_b) = (
        value(&json, &element("a"), "x_max.temperature_k"),
        value(&json, &element("b"), "x_min.temperature_k"),
    );
    let rate = value(&json, &element("wall"), "heat_rate_w");
    assert!(
        COLD_K < face_b && face_b < face_a && face_a < HOT_K && rate > 0.0,
        "faces {face_a} K, {face_b} K, rate {rate} W"
    );
    let (code, text, stderr) = run_binary(&["system", &path]);
    assert_eq!(code, 0, "{stderr}");
    assert!(
        text.contains("conduction bodies step fs-conduction"),
        "{text}"
    );
}

#[test]
fn conduction_body_scenes_refuse_bad_boxes_ports_materials_and_solves_by_name() {
    let base = two_slabs(4);
    for (name, from, to, reason) in [
        (
            "face",
            "\"ports\": [\"x_max\"]",
            "\"ports\": [\"w_max\"]",
            "ports name box faces: x_min, x_max",
        ),
        (
            "twice",
            "\"ports\": [\"x_max\"]",
            "\"ports\": [\"x_max\", \"x_max\"]",
            "port x_max is named twice",
        ),
        (
            "no ports",
            "\"ports\": [\"x_max\"]",
            "\"ports\": []",
            "ports is a nonempty array of box faces",
        ),
        (
            "fraction",
            "\"cells\": [4, 2, 2]",
            "\"cells\": [4.5, 2, 2]",
            "cells are whole numbers from 1 to 64",
        ),
        (
            "too many",
            "\"cells\": [4, 2, 2]",
            "\"cells\": [64, 64, 2]",
            "at most 4096 cells",
        ),
        (
            "box",
            "\"box_m\": [0.03, 0.02, 0.02]",
            "\"box_m\": [0.03, 0.02]",
            "box_m is an array of three numbers",
        ),
        (
            "capacity",
            "\"volumetric_heat_capacity_j_m3_k\": 3600000.0",
            "\"volumetric_heat_capacity_j_m3_k\": -3600000.0",
            "is not positive",
        ),
        (
            "solve",
            ",\n      \"solve\": {\"energy_tolerance_j\": 1e-8, \"tolerance\": 1e-13, \"max_iterations\": 4000, \"restart\": 20}",
            "",
            "missing solve (energy_tolerance_j)",
        ),
        (
            "energy",
            "\"energy_tolerance_j\": 1e-8, ",
            "",
            "missing energy_tolerance_j",
        ),
        (
            "tolerance",
            "\"tolerance\": 1e-13",
            "\"tolerance\": 1.5",
            "tolerance must lie in (0, 1)",
        ),
        (
            "restart",
            "\"restart\": 20",
            "\"restart\": 0",
            "restart is a whole number from 1 to 1000",
        ),
        (
            "field",
            "\"temperature_k\": 400.0,",
            "\"temperature_k\": 400.0, \"emissivity\": 0.5,",
            "unsupported field \"emissivity\"",
        ),
    ] {
        assert!(base.contains(from), "{name}: test fixture changed");
        let (code, stderr) = refused(name, &base.replacen(from, to, 1));
        assert_eq!(code, 3, "{name}: {stderr}");
        assert!(stderr.contains(reason), "{name}: {stderr}");
    }
}
