//! `contact_interface`: fs-dcontact's one-point normal contact, lowered from
//! a scene onto the kernel's load line.
//!
//! The checks are the owner's own closed form and the energy the contact
//! must hand back. A mass on a Hertz wall bounces to the closed-form depth
//! and leaves at the speed it came with (the same case and bounds as
//! fs-dcontact tests/component.rs). A stick strikes a mode of
//! examples/system/struck_mode.json through a contact that reads the mode's
//! end, and the stick, the mode and the contact end holding the energy the
//! stick brought. Then the refusals.

use super::*;
use fs_cli::SystemScene;

/// The live drum's Hertz tip, as in fs-dcontact's tests: `K = 4/3 E*
/// sqrt(R)`, `E* = 0.8 GPa`, `R = 3 mm` [N/m^1.5].
const HERTZ_K: f64 = 4.0 / 3.0 * 0.8e9 * 0.054_772_255_750_516_61;
const ALPHA: f64 = 1.5;
const MASS_KG: f64 = 0.01;
const SPEED_M_S: f64 = 2.0;
/// The mass starts this far from the wall [m].
const APPROACH_M: f64 = 1.0e-5;
/// `contact_interface`'s default search tolerances.
const ABSOLUTE_N: f64 = 1.0e-9;
const RELATIVE: f64 = 1.0e-12;
/// `I = integral from 0 to 1 of (1 - s^(alpha + 1))^(-1/2) ds` for alpha =
/// 3/2, (2/5) B(2/5, 1/2), which fs-dcontact's test computes by Simpson's
/// rule and checks to 1e-5. The contact lasts `2 I p_max / v`.
const DURATION_INTEGRAL: f64 = 1.471_639;

/// A mass at `SPEED_M_S` towards a wall held at rest, the contact between
/// them, in steps of `step_s` for `duration_s`.
fn bounce_scene(step_s: f64, duration_s: f64) -> String {
    format!(
        r#"{{
  "schema": "frankensim.system.v0", "duration_s": {duration_s}, "step_s": {step_s},
  "elements": {{
    "ball": {{"kind": "mass", "mass_kg": {MASS_KG}, "velocity_m_s": -{SPEED_M_S}}},
    "tip": {{"kind": "contact_interface", "stiffness": {HERTZ_K}, "exponent": {ALPHA},
            "closure_m": -{APPROACH_M},
            "solve": {{"maximum_force_n": 1000000.0, "maximum_penetration_m": 0.01}}}}
  }},
  "connections": [
    {{"join": ["ball.p", "tip.a"]}},
    {{"prescribe": "tip.b", "velocity_m_s": 0}}
  ]
}}"#
    )
}

fn probe(scene: &SystemScene, path: &str) -> f64 {
    scene
        .probe(path)
        .unwrap()
        .unwrap_or_else(|| panic!("{path} is not reported"))
}

/// The most work the settled reactions can do beyond the secants over a
/// contact that closes to `peak_m` and opens again: the residual tolerance
/// at the deepest force over the travel in and out.
fn defect_j(peak_m: f64) -> f64 {
    (ABSOLUTE_N + RELATIVE * HERTZ_K * peak_m.powf(ALPHA)) * 2.0 * peak_m
}

#[test]
fn a_mass_on_a_hertz_wall_from_a_scene_bounces_at_the_closed_form_depth_and_speed() {
    // Energy: m v^2 / 2 = K p^(alpha + 1) / (alpha + 1) at the deepest point.
    let peak = ((ALPHA + 1.0) * MASS_KG * SPEED_M_S * SPEED_M_S / (2.0 * HERTZ_K))
        .powf(1.0 / (ALPHA + 1.0));
    let duration = 2.0 * DURATION_INTEGRAL * peak / SPEED_M_S;
    let kinetic_j = MASS_KG * SPEED_M_S * SPEED_M_S / 2.0;
    let steps = 400.0;
    let step_s = duration / steps;
    let text = bounce_scene(step_s, 2.0 * duration);
    scene_context(|cx, _| {
        let mut scene = SystemScene::from_text(&text, "bounce.json").unwrap();
        let (mut entered, mut left, mut deepest) = (false, false, 0.0_f64);
        let (mut after_exit, mut run) = (0, 0_u32);
        while after_exit < 10 && scene.advance(f64::MAX, cx).unwrap() {
            run += 1;
            let closure = probe(&scene, "tip.closure_m");
            deepest = deepest.max(closure);
            entered |= closure > 0.0;
            left |= entered && closure <= 0.0;
            if left {
                after_exit += 1;
            }
        }
        assert!(left, "the ball leaves the wall within the run");
        let speed = probe(&scene, "ball.velocity_m_s") / SPEED_M_S - 1.0;
        let depth = deepest / peak - 1.0;
        let sampling =
            (ALPHA + 1.0) * DURATION_INTEGRAL * DURATION_INTEGRAL / (4.0 * steps * steps);
        let residual = scene.ledger().energy_residual_j();
        let rounding = 8.0 * f64::from(run) * f64::EPSILON * kinetic_j;
        eprintln!(
            "bounce: {run} steps; exit speed {speed:e} (bound {:e}); depth {depth:e} (bound -{sampling:e}); ledger residual {residual:e} J (bound {:e})",
            defect_j(peak) / (2.0 * kinetic_j) + 1e-12,
            defect_j(peak) + rounding
        );
        // The secant and the mass's midpoint rule conserve the energy, so
        // the ball leaves at the speed it came with, to the reaction's
        // tolerance; the 1e-12 is fs-dcontact's allowance for rounding.
        assert!(
            speed.abs() <= defect_j(peak) / (2.0 * kinetic_j) + 1e-12,
            "exit speed off by {speed:e}"
        );
        // A step lands within half a step of the turning point, where the
        // closure lies below the peak by at most (alpha + 1) I^2 / (4 n^2) of
        // it for n steps in the contact; never above it.
        assert!(
            depth <= 1e-12 && depth >= -sampling,
            "deepest closure {deepest} against {peak}: {depth:e}"
        );
        // The ledger books each step's work and stored change; besides the
        // settled reactions' defect it may round at a few ulps a step.
        assert!(
            residual.abs() <= defect_j(peak) + rounding,
            "ledger residual {residual:e} J"
        );
    });
}

#[test]
fn a_stick_strikes_a_mode_through_a_contact_that_reads_its_end_and_the_energy_stays() {
    let path = example("struck_mode.json");
    let text = std::fs::read_to_string(&path).unwrap();
    // The example's stick, mode and attachment weight (1/sqrt(0.05 kg)).
    let (stick_kg, speed, omega) = (0.01, 2.0, core::f64::consts::TAU * 200.0);
    let brought_j = stick_kg * speed * speed / 2.0;
    // No contact can close deeper than the stick's whole energy allows.
    let ceiling_m = ((ALPHA + 1.0) * brought_j / HERTZ_K).powf(1.0 / (ALPHA + 1.0));
    scene_context(|cx, _| {
        let mut scene = SystemScene::from_text(&text, &path).unwrap();
        let (mut contacts, mut touching, mut steps) = (0_u32, false, 0_u32);
        while scene.advance(f64::MAX, cx).unwrap() {
            steps += 1;
            let closing = probe(&scene, "tip.closure_m") > 0.0;
            if closing && !touching {
                contacts += 1;
            }
            touching = closing;
        }
        assert!(contacts >= 1, "the stick reaches the mode");
        let stick = probe(&scene, "stick.velocity_m_s");
        assert!(stick > 0.0, "the stick comes back: {stick} m/s");
        let (q, v) = (
            probe(&scene, "head.mode0.displacement"),
            probe(&scene, "head.mode0.velocity"),
        );
        let mode_j = 0.5 * (v * v + omega * omega * q * q);
        assert!(mode_j > 0.0, "the mode rings");
        let closure = probe(&scene, "tip.closure_m").max(0.0);
        let contact_j = HERTZ_K / (ALPHA + 1.0) * closure.powf(ALPHA + 1.0);
        let held_j = 0.5 * stick_kg * stick * stick + mode_j + contact_j;
        // Each contact can miss by its settled reactions' defect at most;
        // the sums of the three energies round at a few ulps a step.
        let allowed = f64::from(contacts) * defect_j(ceiling_m)
            + 8.0 * f64::from(steps) * f64::EPSILON * brought_j;
        let residual = scene.ledger().energy_residual_j();
        eprintln!(
            "struck mode: {steps} steps, {contacts} contacts; stick {stick} m/s, mode {mode_j:e} J, contact {contact_j:e} J; held - brought {:e} J, ledger residual {residual:e} J (bound {allowed:e})",
            held_j - brought_j
        );
        assert!(
            (held_j - brought_j).abs() <= allowed,
            "held {held_j} J of {brought_j} J ({contacts} contacts, allowed {allowed:e})"
        );
        assert!(residual.abs() <= allowed, "ledger residual {residual:e} J");
    });
    let (code, json, stderr) = run_binary(&["--json", "system", &path]);
    assert_eq!(code, 0, "{stderr}");
    assert!(
        value(&json, &element("tip"), "reaction_n").is_finite(),
        "{json}"
    );
}

#[test]
fn contact_interface_scenes_refuse_bad_laws_budgets_and_ends_by_name() {
    let base = bounce_scene(1.0e-6, 1.0e-4);
    for (name, from, to, reason) in [
        (
            "exponent",
            "\"exponent\": 1.5",
            "\"exponent\": 0.5",
            "contact exponent must be at least 1",
        ),
        (
            "stiffness",
            "\"stiffness\": ",
            "\"stiffness\": -",
            "contact stiffness must be non-negative and finite",
        ),
        (
            "solve",
            ",\n            \"solve\": {\"maximum_force_n\": 1000000.0, \"maximum_penetration_m\": 0.01}",
            "",
            "missing solve",
        ),
        (
            "ceiling",
            "\"maximum_force_n\": 1000000.0, ",
            "",
            "missing maximum_force_n",
        ),
        (
            "iterations",
            "\"maximum_force_n\": 1000000.0",
            "\"maximum_force_n\": 1000000.0, \"max_iterations\": 200",
            "max_iterations is a whole number from 1 to 128",
        ),
        (
            "tolerance",
            "\"maximum_force_n\": 1000000.0",
            "\"maximum_force_n\": 1000000.0, \"force_relative_tolerance\": 1.5",
            "the relative one below 1",
        ),
        (
            "ends",
            "\"closure_m\": ",
            "\"reads_ends\": [\"c\"], \"closure_m\": ",
            "reads_ends names \"a\" or \"b\"",
        ),
        (
            "field",
            "\"closure_m\": ",
            "\"loss\": 0.1, \"closure_m\": ",
            "unsupported field \"loss\"",
        ),
    ] {
        assert!(base.contains(from), "{name}: test fixture changed");
        let (code, stderr) = refused(name, &base.replace(from, to));
        assert_eq!(code, 3, "{name}: {stderr}");
        assert!(stderr.contains(reason), "{name}: {stderr}");
    }
}

#[test]
fn lossy_contact_scene_turns_missing_rebound_energy_into_heat_and_replays() {
    let path = example("lossy_impact.json");
    let text = std::fs::read_to_string(&path).unwrap();
    scene_context(|cx, _| {
        let mut scene = SystemScene::from_text(&text, &path).unwrap();
        while scene.time_s() < 0.05 {
            assert!(scene.advance(f64::MAX, cx).unwrap());
        }
        assert!(probe(&scene, "tip.closure_m") > 0.0);
        let (bytes, root) = scene.checkpoint().unwrap();
        let mut resumed = SystemScene::from_text(&text, &path).unwrap();
        resumed.resume_checkpoint(&bytes, root).unwrap();
        assert_eq!(resumed.render_json(), scene.render_json());

        scene.run(cx).unwrap();
        resumed.run(cx).unwrap();
        assert_eq!(resumed.render_json(), scene.render_json());
        let rebound = probe(&scene, "ball.velocity_m_s");
        let heat_j = 2.0 * (probe(&scene, "body.temperature_k") - 300.0);
        assert!(rebound > 0.0 && rebound < 1.0, "{rebound}");
        assert!(probe(&scene, "tip.closure_m") < 0.0);
        assert!(heat_j > 0.0, "{heat_j}");
        // The wall is stationary and supplies no work. Contact potential
        // vanishes after separation; only kinetic and thermal energy remain.
        assert!((0.5 * rebound * rebound + heat_j - 0.5).abs() < 5e-9);
        assert!(scene.ledger().energy_residual_j().abs() < 5e-9);

        let changed = text.replace(
            "\"internal_loss_s_per_m\": 0.3",
            "\"internal_loss_s_per_m\": 0.4",
        );
        assert_ne!(changed, text);
        let mut incompatible = SystemScene::from_text(&changed, &path).unwrap();
        assert!(incompatible.resume_checkpoint(&bytes, root).is_err());
    });
    let (code, json, stderr) = run_binary(&["--json", "system", &path]);
    assert_eq!(code, 0, "{stderr}");
    let rebound = value(&json, &element("ball"), "velocity_m_s");
    let heat_j = 2.0 * (value(&json, &element("body"), "temperature_k") - 300.0);
    assert!(rebound > 0.0 && rebound < 1.0);
    assert!(heat_j > 0.0);
    assert!((0.5 * rebound * rebound + heat_j - 0.5).abs() < 5e-9);
}

#[test]
fn zero_contact_loss_keeps_elastic_scene_bytes_and_negative_loss_refuses() {
    let base = bounce_scene(1.0e-6, 1.0e-4);
    let explicit = base.replace(
        "\"closure_m\": ",
        "\"internal_loss_s_per_m\": 0.0, \"closure_m\": ",
    );
    assert_ne!(base, explicit);
    scene_context(|cx, _| {
        let mut original = SystemScene::from_text(&base, "bounce.json").unwrap();
        let mut zero = SystemScene::from_text(&explicit, "bounce.json").unwrap();
        original.run(cx).unwrap();
        zero.run(cx).unwrap();
        assert_eq!(original.render_json(), zero.render_json());
    });
    let negative = base.replace(
        "\"closure_m\": ",
        "\"internal_loss_s_per_m\": -0.1, \"closure_m\": ",
    );
    let (code, error) = refused("negative_contact_loss", &negative);
    assert_eq!(code, 3, "{error}");
    assert!(error.contains("Hunt-Crossley chi"), "{error}");
    let unconnected = base.replace(
        "\"closure_m\": ",
        "\"internal_loss_s_per_m\": 0.1, \"closure_m\": ",
    );
    let error = SystemScene::from_text(&unconnected, "unconnected-contact-heat.json")
        .err()
        .expect("dissipation needs a declared heat destination");
    assert!(error.message.contains("heat"), "{}", error.message);
}
