//! `thermal_timeline`: fs-thermal enthalpy cells in a scene, against the
//! owner's explicit rule written out here, the continuous solutions it
//! approximates, and the fs-mbd damped pair's closed-form heat.

use super::*;

/// h = 1000 (T - 290) J/kg: a sensible chart of 1000 J/(kg K).
const CHART: &str = r#"{"phase": "solid", "enthalpy_curve": [[0.0, 290.0, 0.0, 1000.0], [200000.0, 490.0, 0.0, 1000.0]]}"#;

/// A cell of the 1000 J/(kg K) chart at `kelvin`, with its geometry.
fn cell(id: u32, mass_kg: f64, kelvin: f64, extra: &str) -> String {
    let h = 1000.0 * (kelvin - 290.0);
    let volume = mass_kg / 1000.0;
    let x = f64::from(id) * 0.1;
    let end = x + 0.1;
    format!(
        r#"{{"id": {id}, "law": "metal", "mass_kg": {mass_kg:?}, "specific_enthalpy_j_kg": {h:?},
          "volume_m3": {volume:?}, "surface_area_m2": 0.01,
          "start_m": [{x:?}, 0.0, 0.0], "end_m": [{end:?}, 0.0, 0.0]{extra}}}"#
    )
}

/// Two cells of 1000 J/K and 2000 J/K joined inside the timeline by
/// 10 W/K, with a 50 W heater on the first cell's port.
fn heated_pair() -> String {
    format!(
        r#"{{
  "schema": "frankensim.system.v0", "duration_s": 20.0, "step_s": 0.5,
  "elements": {{
    "bar": {{"kind": "thermal_timeline", "max_step_s": 1.0,
      "laws": {{"metal": {CHART}}},
      "cells": [{}, {}],
      "links": [{{"a": 0, "b": 1, "conductance_w_k": 10.0, "area_m2": 0.001}}],
      "ports": [0]}},
    "heater": {{"kind": "heater", "power_w": 50.0}}
  }},
  "connections": [{{"join": ["bar.cell0", "heater.heat"]}}]
}}"#,
        cell(0, 1.0, 300.0, ""),
        cell(1, 2.0, 290.0, "")
    )
}

#[test]
fn a_heated_cell_and_its_linked_neighbour_follow_the_owner_s_explicit_rule() {
    let path = scratch("heated_pair.json", &heated_pair());
    let (code, json, stderr) = run_binary(&["--json", "system", path.to_str().unwrap()]);
    assert_eq!(code, 0, "{stderr}");
    let (c0, c1, g, p, dt) = (1000.0, 2000.0, 10.0, 50.0, 0.5);
    // The owner's rule: the heat each cell gets over a step is computed at
    // the temperatures the step starts from; the port's heat is P dt.
    let (mut t0, mut t1) = (300.0_f64, 290.0_f64);
    for _ in 0..40 {
        let q = g * (t0 - t1);
        t0 += dt * (p - q) / c0;
        t1 += dt * q / c1;
    }
    let found0 = value(&json, &element("bar"), "cell0.temperature_k");
    let found1 = value(&json, &element("bar"), "cell1.temperature_k");
    assert!((found0 - t0).abs() < 1.0e-9, "{found0} vs {t0}");
    assert!((found1 - t1).abs() < 1.0e-9, "{found1} vs {t1}");
    assert!(
        (value(&json, &element("bar"), "cell0.specific_enthalpy_j_kg") - 1000.0 * (t0 - 290.0))
            .abs()
            < 1.0e-6
    );
    // The continuous solution: the mean rises at P / (C0 + C1) and the
    // difference relaxes at G (1/C0 + 1/C1) towards P / (C0 rate).
    let rate = g * (1.0 / c0 + 1.0 / c1);
    let limit = p / (c0 * rate);
    let difference = limit + (10.0 - limit) * (-rate * 20.0_f64).exp();
    let mean = (c0 * 300.0 + c1 * 290.0) / (c0 + c1) + p * 20.0 / (c0 + c1);
    let continuous0 = mean + difference * c1 / (c0 + c1);
    assert!(
        (found0 - continuous0).abs() < 0.01,
        "{found0} vs {continuous0}"
    );
    // Every joule the heater gave is stored, and the ledger closes.
    let stored = c0 * (found0 - 300.0) + c1 * (found1 - 290.0);
    assert!((stored - 1000.0).abs() < 1.0e-7, "{stored}");
    assert!((value(&json, &element("bar"), "stored_change_j") - 1000.0).abs() < 1.0e-7);
    assert!((value(&json, &element("heater"), "energy_in_j") + 1000.0).abs() < 1.0e-9);
    assert!(
        value(&json, "\"totals\"", "energy_residual_j").abs() < 1.0e-9,
        "{json}"
    );
    let steps = value(&json, "\"schema\"", "steps");
    let shortenings = value(&json, "\"schema\"", "shortenings");
    assert!(
        (steps - 40.0).abs() <= 0.0 && shortenings.abs() <= 0.0,
        "{steps} steps, {shortenings} shortenings"
    );
    assert!(json.contains("fs-thermal owner's explicit rule"), "{json}");
}

/// One cell of 1000 J/K at 350 K with a film of 50 W/(m2 K) over 0.01 m2
/// to 290 K, carried on the environment wire when `environment` is set
/// (otherwise the cell has a port, joined to a heater of no power).
fn cooling_cup(environment: bool) -> String {
    let (wire, heater, connection) = if environment {
        (
            r#""environment": true"#,
            "",
            r#"{"prescribe": "cup.environment", "temperature_k": 290.0}"#,
        )
    } else {
        (
            r#""ports": [0]"#,
            r#", "heater": {"kind": "heater", "power_w": 0.0}"#,
            r#"{"join": ["cup.cell0", "heater.heat"]}"#,
        )
    };
    format!(
        r#"{{
  "schema": "frankensim.system.v0", "duration_s": 200.0, "step_s": 1.0,
  "elements": {{
    "cup": {{"kind": "thermal_timeline", "max_step_s": 2.0,
      "laws": {{"metal": {CHART}}},
      "cells": [{}],
      {wire}}}{heater}
  }},
  "connections": [{connection}]
}}"#,
        cell(
            0,
            1.0,
            350.0,
            r#", "boundary": {"ambient_k": 290.0, "htc_w_m2_k": 50.0}"#
        )
    )
}

#[test]
fn a_cell_cools_through_its_film_and_the_heat_leaves_on_the_environment_wire() {
    let path = scratch("cooling_cup.json", &cooling_cup(true));
    let (code, json, stderr) = run_binary(&["--json", "system", path.to_str().unwrap()]);
    assert_eq!(code, 0, "{stderr}");
    let (capacity, film, dt) = (1000.0, 0.5, 1.0);
    let mut excess = 60.0_f64;
    let mut last_rate = f64::NAN;
    for _ in 0..200 {
        last_rate = -film * excess;
        excess += dt * last_rate / capacity;
    }
    let found = value(&json, &element("cup"), "cell0.temperature_k");
    assert!(
        (found - (290.0 + excess)).abs() < 1.0e-9,
        "{found} vs {}",
        290.0 + excess
    );
    assert!(
        (found - (290.0 + 60.0 * (-0.1_f64).exp())).abs() < 0.01,
        "{found}"
    );
    // The last step's mean boundary heat, and all of it on the ledger's
    // prescribed environment port.
    let rate = value(&json, &element("cup"), "environment_heat_w");
    assert!((rate - last_rate).abs() < 1.0e-9, "{rate} vs {last_rate}");
    let crossed = value(&json, "\"port\":\"cup.environment\"", "energy_in_j");
    assert!(
        (crossed - capacity * (found - 350.0)).abs() < 1.0e-7,
        "{crossed}"
    );
    assert!(
        value(&json, "\"totals\"", "energy_residual_j").abs() < 1.0e-9,
        "{json}"
    );

    // Without the environment port the boundary heat has no wire: refused.
    let (code, stderr) = refused("no_environment.json", &cooling_cup(false));
    assert_eq!(code, 3, "{stderr}");
    assert!(
        stderr.contains("cell 0 has a boundary that exchanges heat and no environment port"),
        "{stderr}"
    );
}

#[test]
fn a_contact_made_and_broken_and_a_film_switched_on_take_effect_at_their_times() {
    let scene = format!(
        r#"{{
  "schema": "frankensim.system.v0", "duration_s": 10.0, "step_s": 0.5,
  "sample_every_s": 0.5,
  "elements": {{
    "pair": {{"kind": "thermal_timeline", "max_step_s": 1.0,
      "laws": {{"metal": {CHART}}},
      "cells": [{}, {}],
      "events": [
        {{"time_s": 2.0, "contact": {{"id": 7, "a": 0, "b": 1, "conductance_w_k": 10.0, "area_m2": 0.001}}}},
        {{"time_s": 6.0, "separation": 7}},
        {{"time_s": 8.0, "boundary": {{"cell": 1, "ambient_k": 280.0, "htc_w_m2_k": 50.0}}}}
      ],
      "environment": true, "publish": [0]}}
  }},
  "connections": [{{"prescribe": "pair.environment", "temperature_k": 280.0}}]
}}"#,
        cell(0, 1.0, 300.0, ""),
        cell(1, 2.0, 290.0, "")
    );
    let path = scratch("events.json", &scene);
    let (code, json, stderr) = run_binary(&["--json", "system", path.to_str().unwrap()]);
    assert_eq!(code, 0, "{stderr}");
    let (c0, c1, g, film, dt) = (1000.0, 2000.0, 10.0, 0.5, 0.5);
    let (mut t0, mut t1) = (300.0_f64, 290.0_f64);
    let mut expected = vec![(0.0, t0, t1)];
    for step in 0..20 {
        let start = dt * f64::from(step);
        let q = if (2.0..6.0).contains(&start) {
            g * (t0 - t1)
        } else {
            0.0
        };
        let lost = if start >= 8.0 {
            film * (t1 - 280.0)
        } else {
            0.0
        };
        t0 -= dt * q / c0;
        t1 += dt * (q - lost) / c1;
        expected.push((start + dt, t0, t1));
    }
    let samples: Vec<&str> = json.split("{\"time_s\":").skip(1).collect();
    assert_eq!(samples.len(), expected.len(), "{json}");
    for (sample, (time, kelvin0, kelvin1)) in samples.iter().zip(&expected) {
        let at: f64 = sample.split_once(',').unwrap().0.parse().unwrap();
        assert!((at - time).abs() <= 0.0, "{at} vs {time}");
        let found0 = value(sample, "", "pair.cell0.temperature_k");
        let found1 = value(sample, "", "pair.cell1.temperature_k");
        assert!(
            (found0 - kelvin0).abs() < 1.0e-9,
            "t = {at}: {found0} vs {kelvin0}"
        );
        assert!(
            (found1 - kelvin1).abs() < 1.0e-9,
            "t = {at}: {found1} vs {kelvin1}"
        );
    }
    // Nothing moved before the contact, and the first cell kept what it
    // had once the link was broken.
    let first = |index: usize| value(samples[index], "", "pair.cell0.temperature_k");
    assert!((first(4) - 300.0).abs() <= 0.0, "{}", first(4));
    assert!(first(12) < 300.0 && (first(20) - first(12)).abs() <= 0.0);
    // The film's heat left on the environment wire; the ledger closes.
    let crossed = value(&json, "\"port\":\"pair.environment\"", "energy_in_j");
    let stored = c0 * (t0 - 300.0) + c1 * (t1 - 290.0);
    assert!((crossed - stored).abs() < 1.0e-7, "{crossed} vs {stored}");
    assert!(
        value(&json, "\"totals\"", "energy_residual_j").abs() < 1.0e-9,
        "{json}"
    );
    // Saved at 4 s, where the link the event made is active and its
    // separation and the film are still ahead, and resumed by another
    // process, the run reports what the run with no pause reports.
    let (scene, checkpoint) = (path.to_str().unwrap(), path.with_extension("fssnap"));
    let checkpoint = checkpoint.to_str().unwrap();
    let save = ["--until", "4", "--whole-steps", "--checkpoint", checkpoint];
    let (code, prefix, stderr) = run_binary(&[&["--json", "system", scene], &save[..]].concat());
    assert_eq!(code, 0, "{stderr}");
    let root = prefix.split_once("\"resume_root\":\"").unwrap().1;
    let resume = ["--resume", checkpoint, "--resume-root", &root[..64]];
    let (code, resumed, stderr) = run_binary(&[&["--json", "system", scene], &resume[..]].concat());
    assert_eq!(code, 0, "{stderr}");
    assert_eq!(resumed, json);
}

#[test]
fn the_damped_pair_heats_a_lining_of_enthalpy_cells_with_its_closed_form_heat() {
    let (code, json, stderr) =
        run_binary(&["--json", "system", &example("damped_pair_lining.json")]);
    assert_eq!(code, 0, "{stderr}");
    // Two 2 kg bodies closing at 2 m/s through 0.7 N s/m: the relative
    // velocity decays as exp(-0.7 t) and the kinetic energy lost by t is
    // 2 (1 - exp(-1.4 t)) J, all of it routed into the first cell.
    let heat = 2.0 * (1.0 - (-1.4_f64 * 4.0).exp());
    let capacity = 0.01 * 500.0;
    let kelvin: Vec<f64> = (0..3)
        .map(|id| {
            value(
                &json,
                &element("lining"),
                &format!("cell{id}.temperature_k"),
            )
        })
        .collect();
    let stored: f64 = kelvin.iter().map(|t| capacity * (t - 300.0)).sum();
    assert!((stored - heat).abs() < 1.0e-9, "{stored} vs {heat}");
    assert!((value(&json, &element("lining"), "stored_change_j") - heat).abs() < 1.0e-9);
    assert!((value(&json, &element("pair"), "stored_change_j") + heat).abs() < 1.0e-9);
    // The heat enters at the first cell and spreads along the links.
    assert!(
        kelvin[0] > kelvin[1] && kelvin[1] > kelvin[2] && kelvin[2] > 300.0,
        "{kelvin:?}"
    );
    assert!(
        value(&json, "\"totals\"", "energy_residual_j").abs() < 1.0e-9,
        "{json}"
    );
}

#[test]
fn thermal_timelines_refuse_bad_declarations_and_adaptive_steps_by_name() {
    let scene = heated_pair();
    let first = cell(0, 1.0, 300.0, "");
    let cells = format!("{}, {}", first, cell(1, 2.0, 290.0, ""));
    let duplicate = format!("{}, {}", first, cell(0, 2.0, 290.0, ""));
    let flat = first.replace("\"end_m\": [0.1, 0.0, 0.0]", "\"end_m\": [0.0, 0.0, 0.0]");
    // Every one is an input refusal, exit 3.
    for (name, from, to, reason) in [
        (
            "law.json",
            "\"laws\": {\"metal\"",
            "\"laws\": {\"brass\"",
            "no law named \"metal\"",
        ),
        (
            "twice.json",
            cells.as_str(),
            duplicate.as_str(),
            "cell 0 is declared twice",
        ),
        (
            "port.json",
            "\"ports\": [0]",
            "\"ports\": [5]",
            "cell 5 has a port and is not active",
        ),
        (
            "step.json",
            "\"max_step_s\": 1.0,",
            "",
            "missing max_step_s",
        ),
        (
            "flat.json",
            first.as_str(),
            flat.as_str(),
            "invalid birth geometry",
        ),
        (
            "birth.json",
            "\"ports\": [0]",
            r#""ports": [0], "events": [{"time_s": 1.0, "birth": {"id": 2}}]"#,
            "an event is {time_s, contact | separation | boundary}",
        ),
        (
            "contact.json",
            "\"ports\": [0]",
            r#""ports": [0], "events": [{"time_s": 1.0, "contact": {"a": 0, "b": 1, "conductance_w_k": 1.0, "area_m2": 0.001}}]"#,
            "a link made by an event names its id",
        ),
        (
            "field.json",
            "\"ports\": [0]",
            "\"ports\": [0], \"radiation\": true",
            "unsupported field \"radiation\"",
        ),
        (
            "adaptive.json",
            "\"step_s\": 0.5,",
            r#""step_s": 0.5, "adaptive": {"tolerances": {"heater.energy_in_j": {"absolute": 1e-6, "relative": 0}}},"#,
            "bar is a thermal_timeline",
        ),
    ] {
        assert_eq!(scene.matches(from).count(), 1, "{name}: {from}");
        let (code, stderr) = refused(name, &scene.replacen(from, to, 1));
        assert_eq!(code, 3, "{name}: {stderr}");
        assert!(stderr.contains(reason), "{name}: {stderr}");
    }
}
