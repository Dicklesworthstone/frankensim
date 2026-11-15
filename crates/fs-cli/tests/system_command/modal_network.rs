//! `modal_structure` with a `network`: fs-nlmodal solves the modes of a
//! lumped stiffness and mass, against the chain wall - 400 N/m - 1 kg -
//! 200 N/m - 0.5 kg solved by hand here, and against the same chain built
//! from the kernel's mass and spring elements.
//!
//! The chain's pencil, K = [[600, -200], [-200, 200]] N/m and
//! M = diag(1, 0.5) kg: det(K - l M) = 0.5 l^2 - 500 l + 80000, so
//! w^2 = 200 and 800 rad^2/s^2, with mass-normalized shapes (1, 2) / sqrt(3)
//! and (1, -1) sqrt(2 / 3). The tip's squared shape entries are 4/3 and 2/3.

use super::*;

const W2: [f64; 2] = [200.0, 800.0];

/// A network chain with the tip attached (and the root, when `root`),
/// both prescribed free, damped by `rayleigh`, started at `start`.
fn chain(rayleigh: &str, start: &str, duration_s: f64) -> String {
    format!(
        r#"{{
  "schema": "frankensim.system.v0", "duration_s": {duration_s}, "step_s": 0.001,
  "elements": {{
    "chain": {{"kind": "modal_structure", "sample_rate_hz": 1000,
              "network": {{"stiffness": [[600.0, -200.0], [-200.0, 200.0]],
                           "mass": [[1.0, 0.0], [0.0, 0.5]]{rayleigh}}},
              {start}
              "attachments": {{"tip": {{"domain": "translation", "weights": [[0.0, 1.0]]}},
                               "root": {{"domain": "translation", "weights": [[1.0, 0.0]]}}}}}},
    "case": {{"kind": "heat_capacity", "capacity_j_k": 1.0, "temperature_k": 300.0}}
  }},
  "connections": [
    {{"prescribe": "chain.tip", "force_n": 0}},
    {{"prescribe": "chain.root", "force_n": 0}},
    {{"join": ["case.heat", "chain.heat"]}}
  ]
}}"#
    )
}

/// Each mode's share of a start of 0.01 m at the tip, in the tip's and the
/// root's positions: the tip's are 0.02/3 and 0.01/3, the root's 0.01/3
/// and -0.01/3 (the shapes times the modal start Phi^T M x).
const TIP_SHARE: [f64; 2] = [0.02 / 3.0, 0.01 / 3.0];
const ROOT_SHARE: [f64; 2] = [0.01 / 3.0, -0.01 / 3.0];

#[test]
fn a_chain_released_from_a_stretch_rings_as_its_hand_solved_modes() {
    let scene = chain("", r#""start": {"displacement": [0.0, 0.01]},"#, 1.0);
    let path = scratch("chain_ring.json", &scene);
    let (code, json, stderr) = run_binary(&["--json", "system", path.to_str().unwrap()]);
    assert_eq!(code, 0, "{stderr}");
    // Undamped from rest at a stretch: x = sum share cos(w t), so the
    // velocity is -sum share w sin(w t).
    let t = 1.0_f64;
    let velocity = |shares: [f64; 2]| -> f64 {
        shares
            .iter()
            .zip(W2)
            .map(|(share, square)| -share * square.sqrt() * (square.sqrt() * t).sin())
            .sum()
    };
    for (row, shares) in [("tip", TIP_SHARE), ("root", ROOT_SHARE)] {
        let found = value(&json, &element("chain"), &format!("{row}.velocity_m_s"));
        let expected = velocity(shares);
        assert!(
            (found - expected).abs() <= 1.0e-11,
            "{row}: {found} vs {expected}"
        );
    }
    // The modes hold the stretch's energy, 1/2 x^T K x = 0.01 J, undamped.
    let mut energy = 0.0;
    for (mode, square) in W2.iter().enumerate() {
        let q = value(
            &json,
            &element("chain"),
            &format!("mode{mode}.displacement"),
        );
        let v = value(&json, &element("chain"), &format!("mode{mode}.velocity"));
        energy += 0.5 * square.mul_add(q * q, v * v);
    }
    assert!((energy - 0.01).abs() <= 1.0e-14, "{energy}");
    assert!(value(&json, &element("chain"), "heat_w").abs() <= 0.0);
    assert!(
        json.contains("solves from an authored stiffness and mass"),
        "{json}"
    );
}

#[test]
fn rayleigh_damping_rings_each_mode_down_at_its_own_ratio_and_heats_the_case() {
    let (a, b) = (0.5_f64, 1.0e-3_f64);
    let scene = chain(
        r#", "rayleigh": [0.5, 0.001]"#,
        r#""start": {"displacement": [0.0, 0.01]},"#,
        1.0,
    );
    let path = scratch("chain_rayleigh.json", &scene);
    let (code, json, stderr) = run_binary(&["--json", "system", path.to_str().unwrap()]);
    assert_eq!(code, 0, "{stderr}");
    // C = a M + b K is diagonal in the modes: c_k = a + b w_k^2 = 2 zeta w.
    // From rest at q0: q' = -q0 w^2 / wd e^(-zeta w t) sin(wd t).
    let t = 1.0_f64;
    let (mut tip, mut energy) = (0.0, 0.0);
    // The modal starts' squares, (Phi^T M x)^2: 1e-4 / 3 and 1e-4 / 6.
    for ((share, square), start_squared) in
        TIP_SHARE.iter().zip(W2).zip([1.0e-4 / 3.0, 1.0e-4 / 6.0])
    {
        let w = square.sqrt();
        let zeta = b.mul_add(square, a) / (2.0 * w);
        let wd = w * (1.0 - zeta * zeta).sqrt();
        let decay = (-zeta * w * t).exp();
        tip += -share * square / wd * decay * (wd * t).sin();
        let shape = (wd * t).cos() + zeta * w / wd * (wd * t).sin();
        let q2 = start_squared * (decay * shape).powi(2);
        let v2 = start_squared * (square / wd * decay * (wd * t).sin()).powi(2);
        energy += 0.5 * square.mul_add(q2, v2);
    }
    let found = value(&json, &element("chain"), "tip.velocity_m_s");
    assert!((found - tip).abs() <= 1.0e-11, "{found} vs {tip}");
    // What the modes lost is the case's heat.
    let heat = 0.01 - energy;
    assert!(heat > 0.001, "{heat}");
    let dissipated = value(&json, &element("chain"), "dissipated_j");
    assert!(
        (dissipated - heat).abs() <= 1.0e-12,
        "{dissipated} vs {heat}"
    );
    let kelvin = value(&json, &element("case"), "temperature_k");
    assert!((kelvin - (300.0 + heat)).abs() <= 1.0e-10, "{kelvin}");
    assert!(
        value(&json, "\"totals\"", "energy_residual_j").abs() < 1.0e-9,
        "{json}"
    );
}

#[test]
fn a_held_tip_force_moves_the_network_exactly_and_the_element_chain_at_second_order() {
    // From rest under F = 1 N at the tip: each mode q = s F / w^2 (1 - cos w t)
    // with s its tip entry, so the tip's velocity is F sum s^2 sin(w t) / w.
    let t = 1.0_f64;
    let squared_shapes = [4.0 / 3.0, 2.0 / 3.0];
    let exact: f64 = squared_shapes
        .iter()
        .zip(W2)
        .map(|(s2, square)| s2 * (square.sqrt() * t).sin() / square.sqrt())
        .sum();
    let scale: f64 = squared_shapes
        .iter()
        .zip(W2)
        .map(|(s2, square)| s2 / square.sqrt())
        .sum();
    let network = r#"{
  "schema": "frankensim.system.v0", "duration_s": 1.0, "step_s": 0.001,
  "elements": {
    "chain": {"kind": "modal_structure", "sample_rate_hz": 1000,
              "network": {"stiffness": [[600.0, -200.0], [-200.0, 200.0]],
                          "mass": [[1.0, 0.0], [0.0, 0.5]]},
              "attachments": {"tip": {"domain": "translation", "weights": [[0.0, 1.0]]}}},
    "push": {"kind": "force", "force_n": 1.0},
    "case": {"kind": "heat_capacity", "capacity_j_k": 1.0, "temperature_k": 300.0}
  },
  "connections": [
    {"join": ["chain.tip", "push.p"]},
    {"join": ["case.heat", "chain.heat"]}
  ]
}"#;
    let path = scratch("chain_pushed.json", network);
    let (code, json, stderr) = run_binary(&["--json", "system", path.to_str().unwrap()]);
    assert_eq!(code, 0, "{stderr}");
    let found = value(&json, &element("chain"), "tip.velocity_m_s");
    assert!(
        (found - exact).abs() <= 1.0e-11 * scale,
        "{found} vs {exact}"
    );
    // The force's work is stored in the modes, undamped.
    let work = -value(&json, &element("push"), "energy_in_j");
    let stored = value(&json, &element("chain"), "stored_change_j");
    assert!((work - stored).abs() <= 1.0e-12, "{work} vs {stored}");

    // The same chain from the kernel's mass and spring elements takes the
    // midpoint step: the same motion, but not exact at this step.
    let elements = r#"{
  "schema": "frankensim.system.v0", "duration_s": 1.0, "step_s": 0.001,
  "elements": {
    "m1": {"kind": "mass", "mass_kg": 1.0},
    "m2": {"kind": "mass", "mass_kg": 0.5},
    "k1": {"kind": "spring", "stiffness_n_m": 400.0},
    "k2": {"kind": "spring", "stiffness_n_m": 200.0},
    "push": {"kind": "force", "force_n": 1.0}
  },
  "connections": [
    {"join": ["m1.p", "k1.a", "k2.a"]},
    {"join": ["m2.p", "k2.b", "push.p"]},
    {"prescribe": "k1.b", "velocity_m_s": 0}
  ]
}"#;
    let path = scratch("element_chain.json", elements);
    let (code, json, stderr) = run_binary(&["--json", "system", path.to_str().unwrap()]);
    assert_eq!(code, 0, "{stderr}");
    let midpoint = value(&json, &element("m2"), "velocity_m_s");
    let error = (midpoint - exact).abs();
    assert!(
        error > 1.0e-7 * scale && error < 0.1 * scale,
        "{midpoint} vs {exact}"
    );
}

#[test]
fn modal_networks_refuse_what_they_cannot_solve_by_name() {
    let scene = chain("", "", 0.01);
    for (name, from, to, reason) in [
        (
            "both.json",
            "\"network\":",
            r#""modes": [{"frequency_hz": 5.0, "damping_ratio": 0.0}], "network":"#,
            "give either modes or network",
        ),
        (
            "square.json",
            "[[600.0, -200.0], [-200.0, 200.0]]",
            "[[600.0, -200.0, 0.0], [-200.0, 200.0, 0.0]]",
            "column count must be 2",
        ),
        (
            "mass.json",
            "[[1.0, 0.0], [0.0, 0.5]]",
            "[[1.0, 0.0], [0.0, -0.5]]",
            "the loaded mass was refused",
        ),
        (
            "asymmetric.json",
            "[[1.0, 0.0], [0.0, 0.5]]",
            "[[1.0, 0.1], [0.0, 0.5]]",
            "differs from its mirror",
        ),
        (
            "free.json",
            "[[600.0, -200.0], [-200.0, 200.0]]",
            "[[200.0, -200.0], [-200.0, 200.0]]",
            "mode 0 has",
        ),
        (
            "weights.json",
            "[[0.0, 1.0]]",
            "[[1.0]]",
            "one weight per degree of freedom, 2 here",
        ),
        (
            "start.json",
            "\"attachments\"",
            r#""start": {"velocity": [1.0]}, "attachments""#,
            "start.velocity: 2 finite numbers, one per degree of freedom",
        ),
        (
            "rayleigh.json",
            "\"mass\": [[1.0, 0.0], [0.0, 0.5]]",
            r#""mass": [[1.0, 0.0], [0.0, 0.5]], "rayleigh": [0.5]"#,
            "rayleigh: [a, b] of finite numbers",
        ),
        (
            "negative.json",
            "\"mass\": [[1.0, 0.0], [0.0, 0.5]]",
            r#""mass": [[1.0, 0.0], [0.0, 0.5]], "rayleigh": [-100.0, 0.0]"#,
            "non-negative finite damping",
        ),
    ] {
        assert_eq!(scene.matches(from).count(), 1, "{name}: {from}");
        let (code, stderr) = refused(name, &scene.replacen(from, to, 1));
        assert_eq!(code, 3, "{name}: {stderr}");
        assert!(stderr.contains(reason), "{name}: {stderr}");
    }
    // A start written in a mode row belongs to the modes form only.
    let authored = r#"{
  "schema": "frankensim.system.v0", "duration_s": 0.01, "step_s": 0.001,
  "elements": {
    "rod": {"kind": "modal_structure", "sample_rate_hz": 1000,
            "modes": [{"frequency_hz": 5.0, "damping_ratio": 0.0}],
            "start": {"displacement": [0.01]},
            "attachments": {"tip": {"domain": "translation", "weights": [[1.0]]}}}
  },
  "connections": [{"prescribe": "rod.tip", "force_n": 0}]
}"#;
    let (code, stderr) = refused("authored_start.json", authored);
    assert_eq!(code, 3, "{stderr}");
    assert!(stderr.contains("start belongs to a network"), "{stderr}");
}
