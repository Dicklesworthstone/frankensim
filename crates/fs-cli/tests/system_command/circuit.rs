//! G1/G3: `kirchhoff_circuit` scenes step fs-phs netlists on the midpoint
//! rule's own closed forms, heating one body per resistor.

use super::scene_context;

fn close(actual: f64, expected: f64, what: &str) {
    assert!(
        (actual - expected).abs() <= 1e-9 * (1.0 + expected.abs()),
        "{what}: {actual} != {expected}"
    );
}

fn probe(scene: &fs_cli::SystemScene, path: &str) -> f64 {
    scene
        .probe(path)
        .unwrap()
        .unwrap_or_else(|| panic!("{path} is null"))
}

fn run(text: &str) -> fs_cli::SystemScene {
    let mut scene = fs_cli::SystemScene::from_text(text, "circuit.json").unwrap();
    scene_context(|cx, _| scene.run(cx).unwrap());
    close(scene.time_s(), 10.0, "time");
    assert!(scene.ledger().energy_residual_j().abs() < 1e-9);
    scene
}

/// A scene of 100 steps of 0.1 s: the circuit `k` with the `fields` given,
/// one heat capacity per entry of `bodies` (name, heat port, J/K) from 300 K
/// joined to that port of `k`, and port 0 held at `volts` if given.
fn scene(fields: &str, bodies: &[(&str, &str, f64)], volts: Option<f64>) -> String {
    let elements = bodies
        .iter()
        .map(|(body, _, capacity)| {
            format!(
                r#",
    "{body}": {{"kind": "heat_capacity", "capacity_j_k": {capacity}, "temperature_k": 300.0}}"#
            )
        })
        .collect::<Vec<_>>()
        .concat();
    let connections = volts
        .map(|volts| format!(r#"{{"prescribe": "k.source0", "voltage_v": {volts}}}"#))
        .into_iter()
        .chain(
            bodies
                .iter()
                .map(|(body, port, _)| format!(r#"{{"join": ["k.{port}", "{body}.heat"]}}"#)),
        )
        .collect::<Vec<_>>()
        .join(", ");
    format!(
        r#"{{
  "schema": "frankensim.system.v0", "duration_s": 10.0, "step_s": 0.1,
  "elements": {{
    "k": {{"kind": "kirchhoff_circuit", {fields},
      "energy_tolerance_j": 1e-9}}{elements}
  }},
  "connections": [{connections}]
}}"#
    )
}

/// A charge `Q0` on `C` across node 1 and ground, discharged through `R1`
/// from node 1 to node 2 and `R2` from node 2 to ground.
const C: f64 = 1.0;
const Q0: f64 = 2.0;
const R1: f64 = 2.0;
const R2: f64 = 3.0;

#[test]
fn a_charged_capacitor_heats_its_two_resistors_in_proportion() {
    let text = scene(
        &format!(
            r#""node_count": 3, "branches": [
        {{"nodes": [1, 0], "capacitance_f": {C}, "charge_c": {Q0}}},
        {{"nodes": [1, 2], "resistance_ohm": {R1}}},
        {{"nodes": [2, 0], "resistance_ohm": {R2}}}]"#
        ),
        &[("one", "heat1", 10.0), ("two", "heat2", 20.0)],
        None,
    );
    let scene = run(&text);
    // The midpoint rule on q' = -q / ((R1 + R2) C).
    let a = 0.1 / (2.0 * (R1 + R2) * C);
    let charge = |steps: i32| Q0 * ((1.0 - a) / (1.0 + a)).powi(steps);
    let q = charge(100);
    close(probe(&scene, "k.charge0_c"), q, "charge");
    // The released field splits R1 : R2; the last step's current is the
    // midpoint charge's over (R1 + R2) C.
    let released = (Q0 * Q0 - q * q) / (2.0 * C);
    let current = 0.5 * (charge(99) + q) / ((R1 + R2) * C);
    for (body, port, ohms, capacity) in [("one", "heat1", R1, 10.0), ("two", "heat2", R2, 20.0)] {
        let share = released * ohms / (R1 + R2);
        close(
            probe(&scene, &format!("{body}.temperature_k")),
            300.0 + share / capacity,
            body,
        );
        close(
            probe(&scene, &format!("k.{port}_w")),
            ohms * current * current,
            port,
        );
    }
    let report = scene.render_json();
    assert!(
        report.contains("a kirchhoff_circuit is fs-phs's descriptor circuit"),
        "{report}"
    );
    assert_eq!(
        probe(&scene, "k.energy_tolerance_j").to_bits(),
        f64::to_bits(1e-9)
    );
}

/// The capacitor first in branch order, so the owner's state (fluxes, then
/// charges) is not in branch order: `CM` with `QM` across node 1 and ground
/// through `RA`, and `V` on port 0 across node 2 and ground driving `RB`
/// from node 2 to node 3 and `L` from node 3 to ground, from flux `PHI0`.
const CM: f64 = 0.5;
const QM: f64 = 1.0;
const RA: f64 = 4.0;
const V: f64 = 6.0;
const RB: f64 = 2.0;
const L: f64 = 0.5;
const PHI0: f64 = 0.25;

fn mixed() -> String {
    scene(
        &format!(
            r#""node_count": 4, "branches": [
        {{"nodes": [1, 0], "capacitance_f": {CM}, "charge_c": {QM}}},
        {{"nodes": [1, 0], "resistance_ohm": {RA}}},
        {{"nodes": [2, 0], "voltage_port": 0}},
        {{"nodes": [2, 3], "resistance_ohm": {RB}}},
        {{"nodes": [3, 0], "inductance_h": {L}, "flux_wb": {PHI0}}}]"#
        ),
        &[("rc_body", "heat1", 50.0), ("rl_body", "heat3", 100.0)],
        Some(V),
    )
}

#[test]
fn a_netlist_reports_its_state_by_name_and_its_source_current() {
    let scene = run(&mixed());
    let a = 0.1 / (2.0 * RA * CM);
    let charge = |steps: i32| QM * ((1.0 - a) / (1.0 + a)).powi(steps);
    // The midpoint rule on L i' = V - RB i from i = PHI0 / L.
    let b = RB * 0.1 / (2.0 * L);
    let (mut current, mut mean, mut supplied) = (PHI0 / L, 0.0, 0.0);
    for _ in 0..100 {
        let next = (current * (1.0 - b) + 0.1 * V / L) / (1.0 + b);
        mean = 0.5 * (current + next);
        supplied += V * mean * 0.1;
        current = next;
    }
    close(probe(&scene, "k.charge0_c"), charge(100), "charge");
    close(probe(&scene, "k.flux4_wb"), L * current, "flux");
    // The source delivers the step's mean current into the circuit.
    close(probe(&scene, "k.source0.current_a"), mean, "source current");
    close(probe(&scene, "k.heat3_w"), RB * mean * mean, "RL heat rate");
    let field = 0.5 * L * (current * current - (PHI0 / L) * (PHI0 / L));
    close(
        probe(&scene, "rl_body.temperature_k"),
        300.0 + (supplied - field) / 100.0,
        "RL body",
    );
    let q = charge(100);
    close(
        probe(&scene, "rc_body.temperature_k"),
        300.0 + (QM * QM - q * q) / (2.0 * CM) / 50.0,
        "RC body",
    );
}

/// `V` on port 0 across the primary of a transformer of ratio `N` (v1 = N
/// v2); the secondary drives `RS` from node 2 to node 3 and `L1` from node 3
/// to ground. `L2` from node 4 to ground, coupled to `L1` by `M`, closes
/// through `RL` from node 4 to ground. The time constants, 1.55 s and
/// 6.45 s, keep the run in its transient.
const N: f64 = 2.0;
const VT: f64 = 10.0;
const RS: f64 = 1.0;
const L1: f64 = 2.0;
const L2: f64 = 3.0;
const M: f64 = 1.0;
const RL: f64 = 0.5;

#[test]
fn coupled_inductors_behind_a_transformer_follow_the_midpoint_rule() {
    // `examples/system/transformer_coupled_loops.json`: the circuit on these
    // constants, `VT` held on port 0, and 40 J/K bodies `series` on heat1 and
    // `loop` on heat4.
    let text = std::fs::read_to_string(super::example("transformer_coupled_loops.json")).unwrap();
    let scene = run(&text);
    // (L / dt + R / 2) i' = (L / dt - R / 2) i + (VT / N, 0) for the
    // inductance matrix L = [[L1, M], [M, L2]] and R = diag(RS, RL): the
    // secondary loop sees VT / N, the second loop only the coupling.
    let dt = 0.1;
    let (a11, a12, a22) = (L1 / dt + RS / 2.0, M / dt, L2 / dt + RL / 2.0);
    let (b11, b22) = (L1 / dt - RS / 2.0, L2 / dt - RL / 2.0);
    let (mut i, mut mean, mut heats) = ([0.0_f64; 2], [0.0_f64; 2], [0.0_f64; 2]);
    for _ in 0..100 {
        let rhs = [b11 * i[0] + a12 * i[1] + VT / N, a12 * i[0] + b22 * i[1]];
        let det = a11 * a22 - a12 * a12;
        let next = [
            (rhs[0] * a22 - a12 * rhs[1]) / det,
            (a11 * rhs[1] - a12 * rhs[0]) / det,
        ];
        mean = [0.5 * (i[0] + next[0]), 0.5 * (i[1] + next[1])];
        heats[0] += RS * mean[0] * mean[0] * dt;
        heats[1] += RL * mean[1] * mean[1] * dt;
        i = next;
    }
    close(probe(&scene, "k.flux2_wb"), L1 * i[0] + M * i[1], "flux 1");
    close(probe(&scene, "k.flux3_wb"), M * i[0] + L2 * i[1], "flux 2");
    // The source's power VT I is the secondary's VT / N i1.
    close(
        probe(&scene, "k.source0.current_a"),
        mean[0] / N,
        "source current",
    );
    close(
        probe(&scene, "k.heat4_w"),
        RL * mean[1] * mean[1],
        "loop heat rate",
    );
    for (body, heat) in ["series", "loop"].into_iter().zip(heats) {
        close(
            probe(&scene, &format!("{body}.temperature_k")),
            300.0 + heat / 40.0,
            body,
        );
    }
}

/// The mixed scene stepped adaptively at `order`, with a tolerance for every
/// account except `without`.
fn adaptive(order: u32, without: &str) -> String {
    let tolerances = [
        "k.flux4_wb",
        "k.charge0_c",
        "rc_body.temperature_k",
        "rl_body.temperature_k",
    ]
    .into_iter()
    .filter(|name| *name != without)
    .map(|name| format!(r#""{name}": {{"absolute": 1e-10, "relative": 1e-10}}"#))
    .collect::<Vec<_>>()
    .join(", ");
    mixed().replace(
        r#""step_s": 0.1,"#,
        &format!(
            r#""step_s": 0.1, "adaptive": {{"order": {order}, "tolerances": {{{tolerances}}}}},"#
        ),
    )
}

#[test]
fn an_adaptive_netlist_accounts_its_state_at_order_two() {
    let scene = run(&adaptive(2, ""));
    // The exact solutions, which the controlled steps approach to 1e-6 of
    // each coordinate's scale (its start, or its end if larger).
    let exact = [
        ("k.charge0_c", QM * (-10.0 / (RA * CM)).exp(), QM),
        (
            "k.flux4_wb",
            L * (V / RB + (PHI0 / L - V / RB) * (-RB * 10.0 / L).exp()),
            L * V / RB,
        ),
    ];
    for (name, value, scale) in exact {
        let reported = probe(&scene, name);
        assert!(
            (reported - value).abs() <= 1e-6 * scale,
            "{name}: {reported} vs {value}"
        );
    }
    for (bad, expected) in [
        (
            adaptive(3, ""),
            "is a kirchhoff_circuit, whose step is of order 2",
        ),
        (
            adaptive(2, "k.charge0_c"),
            "k.charge0_c: missing explicit tolerance",
        ),
    ] {
        let error = fs_cli::SystemScene::from_text(&bad, "circuit.json")
            .err()
            .unwrap_or_else(|| panic!("admitted: {expected}"))
            .to_string();
        assert!(error.contains(expected), "{expected}: {error}");
    }
}

#[test]
fn a_netlist_reports_null_before_a_step_and_bad_netlists_are_refused() {
    let good = mixed();
    let admitted = fs_cli::SystemScene::from_text(&good, "circuit.json").unwrap();
    assert_eq!(admitted.probe("k.source0.current_a").unwrap(), None);
    assert_eq!(admitted.probe("k.heat1_w").unwrap(), None);
    close(probe(&admitted, "k.flux4_wb"), PHI0, "start flux");
    for (bad, expected) in [
        (
            good.replace(
                r#""nodes": [3, 0], "inductance_h""#,
                r#""nodes": [3, 0], "voltage_port": 0, "inductance_h""#,
            ),
            "exactly one of resistance_ohm",
        ),
        (
            good.replace(
                r#""resistance_ohm": 4}"#,
                r#""resistance_ohm": 4}, {"nodes": [1, 0], "voltage_port": 0}"#,
            ),
            "an external port drives more than one source",
        ),
        (
            good.replace(r#""flux_wb": 0.25"#, r#""charge_c": 0.25"#),
            "charge_c belongs to a branch with capacitance_f",
        ),
        (
            good.replace(r#""nodes": [2, 3]"#, r#""nodes": [2, 4]"#),
            "nodes: two whole numbers below 4",
        ),
        (
            good.replace(r#""node_count": 4"#, r#""node_count": 5"#),
            "FloatingNode { node: 4 }",
        ),
        (
            good.replace(
                r#""energy_tolerance_j": 1e-9"#,
                r#""energy_tolerance_j": 1e-9, "mutual": [{"branches": [1, 4], "inductance_h": 0.1}]"#,
            ),
            "circuit owner refused",
        ),
        (
            good.replace(",\n      \"energy_tolerance_j\": 1e-9", ""),
            "missing energy_tolerance_j",
        ),
        (
            good.replace(r#""node_count": 4"#, r#""node_count": 4, "laws": []"#),
            "unsupported field \"laws\"",
        ),
    ] {
        assert_ne!(bad, good, "the replacement must change the scene");
        let error = fs_cli::SystemScene::from_text(&bad, "circuit.json")
            .err()
            .unwrap_or_else(|| panic!("admitted: {expected}"))
            .to_string();
        assert!(error.contains("elements.k"), "{error}");
        assert!(error.contains(expected), "{expected}: {error}");
    }
}
