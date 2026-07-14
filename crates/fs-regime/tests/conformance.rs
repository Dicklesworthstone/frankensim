//! fs-regime conformance suite (the tfz.3 bead). Acceptance: Pi-group
//! extraction algebraically correct on textbook batteries; G3 unit-
//! rescaling invariance of every group (the definitive test); validity
//! predicates catch seeded misuse with ranked alternatives; scaling
//! recommendations measurably improve conditioning; flagship fixtures
//! match hand calculations; reports reproducible.

use fs_math::det;
use fs_qty::{Dims, QtyAny};
use fs_regime::{
    Input, Role, RoleInput, assess, condition_number, flux_model_cards, pi_groups, standard_groups,
};
use fs_regime::{ScalingMap, cards::admit};
use std::collections::BTreeMap;

fn verdict(case: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-regime/conformance\",\"case\":\"{case}\",\"verdict\":\"pass\",\
         \"detail\":\"{detail}\"}}"
    );
}

fn role(role: Role, value: f64, dims: [i8; 6]) -> RoleInput {
    RoleInput {
        role,
        qty: QtyAny::new(value, Dims(dims)),
    }
}

/// The spout-pour fluid state (viscous suspension, mid-tilt).
fn spout_inputs() -> Vec<RoleInput> {
    vec![
        role(Role::Density, 1200.0, [-3, 1, 0, 0, 0, 0]),
        role(Role::Velocity, 0.3, [1, 0, -1, 0, 0, 0]),
        role(Role::Length, 0.02, [1, 0, 0, 0, 0, 0]),
        role(Role::DynViscosity, 8.0, [-1, 1, -1, 0, 0, 0]),
        role(Role::SurfaceTension, 0.06, [0, 1, -2, 0, 0, 0]),
        role(Role::Gravity, 9.806_65, [1, 0, -2, 0, 0, 0]),
    ]
}

/// Battery of textbook sets: (inputs, expected rank, expected count).
fn textbook_batteries() -> Vec<(Vec<Input>, usize, usize)> {
    let q = |v: f64, d: [i8; 6]| QtyAny::new(v, Dims(d));
    vec![
        (
            vec![
                Input {
                    name: "rho".into(),
                    qty: q(998.0, [-3, 1, 0, 0, 0, 0]),
                },
                Input {
                    name: "v".into(),
                    qty: q(1.5, [1, 0, -1, 0, 0, 0]),
                },
                Input {
                    name: "d".into(),
                    qty: q(0.1, [1, 0, 0, 0, 0, 0]),
                },
                Input {
                    name: "mu".into(),
                    qty: q(1e-3, [-1, 1, -1, 0, 0, 0]),
                },
            ],
            3,
            1,
        ),
        (
            // Drag: (F, rho, V, L, mu) -> Cd and Re.
            vec![
                Input {
                    name: "f".into(),
                    qty: q(12.0, [1, 1, -2, 0, 0, 0]),
                },
                Input {
                    name: "rho".into(),
                    qty: q(1.225, [-3, 1, 0, 0, 0, 0]),
                },
                Input {
                    name: "v".into(),
                    qty: q(8.0, [1, 0, -1, 0, 0, 0]),
                },
                Input {
                    name: "l".into(),
                    qty: q(0.12, [1, 0, 0, 0, 0, 0]),
                },
                Input {
                    name: "mu".into(),
                    qty: q(1.81e-5, [-1, 1, -1, 0, 0, 0]),
                },
            ],
            3,
            2,
        ),
        (
            // Heat convection: (h, k, L, rho, mu, cp, V) -> Nu, Re, Pr.
            vec![
                Input {
                    name: "h".into(),
                    qty: q(25.0, [0, 1, -3, -1, 0, 0]),
                },
                Input {
                    name: "k".into(),
                    qty: q(0.6, [1, 1, -3, -1, 0, 0]),
                },
                Input {
                    name: "l".into(),
                    qty: q(0.05, [1, 0, 0, 0, 0, 0]),
                },
                Input {
                    name: "rho".into(),
                    qty: q(998.0, [-3, 1, 0, 0, 0, 0]),
                },
                Input {
                    name: "mu".into(),
                    qty: q(1e-3, [-1, 1, -1, 0, 0, 0]),
                },
                Input {
                    name: "cp".into(),
                    qty: q(4184.0, [2, 0, -2, -1, 0, 0]),
                },
                Input {
                    name: "v".into(),
                    qty: q(0.8, [1, 0, -1, 0, 0, 0]),
                },
            ],
            4,
            3,
        ),
        (
            // Amount-bearing battery: (n, L, c) with c in mol/m³.
            vec![
                Input {
                    name: "n".into(),
                    qty: q(2.0, [0, 0, 0, 0, 0, 1]),
                },
                Input {
                    name: "l".into(),
                    qty: q(3.0, [1, 0, 0, 0, 0, 0]),
                },
                Input {
                    name: "c".into(),
                    qty: q(2.0 / 27.0, [-3, 0, 0, 0, 0, 1]),
                },
            ],
            2,
            1,
        ),
    ]
}

#[test]
fn rg_001_pi_batteries_are_algebraically_correct() {
    for (i, (inputs, want_rank, want_count)) in textbook_batteries().iter().enumerate() {
        let basis = pi_groups(inputs).expect("pi extraction");
        assert_eq!(basis.rank, *want_rank, "battery {i} rank");
        assert_eq!(
            basis.groups.len(),
            *want_count,
            "battery {i} Buckingham count"
        );
        // Every group is exactly dimensionless by integer arithmetic.
        for g in &basis.groups {
            let mut residual = [0i64; 6];
            for (input, &e) in inputs.iter().zip(&g.exponents) {
                for (slot, &d) in residual.iter_mut().zip(&input.qty.dims.0) {
                    *slot += i64::from(d) * e;
                }
            }
            assert_eq!(residual, [0; 6], "battery {i}: group not dimensionless");
            assert!(g.value.is_finite() && g.value > 0.0);
        }
    }
    verdict(
        "rg-001",
        "pipe/drag/convection/amount batteries: rank, Buckingham count, integer dimensionlessness",
    );
}

#[test]
fn rg_002_unit_rescaling_invariance_is_exact() {
    // THE definitive G3 test: re-express every input in a different
    // coherent unit system (scale factors per SI base dim) and recompute
    // the groups from the RAW rescaled numbers. Dimensionless products
    // must be invariant because the scale factors cancel exactly.
    // The amount scale is deliberately nontrivial; these legacy mol-free
    // fixtures must remain invariant because their sixth exponent is zero.
    let scales = [1000.0f64, 0.001, 60.0, 1.8, 3.0, 7.0];
    let rescale = |q: QtyAny| -> QtyAny {
        let mut v = q.value;
        for (&s, d) in scales.iter().zip(q.dims.0) {
            v *= det::powi(s, i32::from(d));
        }
        QtyAny::new(v, q.dims)
    };
    let base = spout_inputs();
    let rescaled: Vec<RoleInput> = base
        .iter()
        .map(|r| RoleInput {
            role: r.role,
            qty: rescale(r.qty),
        })
        .collect();
    let g_base = standard_groups(&base).expect("base groups");
    let g_resc = standard_groups(&rescaled).expect("rescaled groups");
    assert_eq!(g_base.len(), g_resc.len());
    for (a, b) in g_base.iter().zip(&g_resc) {
        assert_eq!(a.name, b.name);
        assert!(
            (a.value - b.value).abs() / a.value < 1e-12,
            "{}: {} vs {} under rescaling",
            a.name,
            a.value,
            b.value
        );
    }
    // Same law for the raw Pi machinery.
    let pi_a = pi_groups(
        &base
            .iter()
            .map(|r| Input {
                name: r.role.tag().to_string(),
                qty: r.qty,
            })
            .collect::<Vec<_>>(),
    )
    .expect("pi base");
    let pi_b = pi_groups(
        &rescaled
            .iter()
            .map(|r| Input {
                name: r.role.tag().to_string(),
                qty: r.qty,
            })
            .collect::<Vec<_>>(),
    )
    .expect("pi rescaled");
    for (a, b) in pi_a.groups.iter().zip(&pi_b.groups) {
        assert_eq!(a.exponents, b.exponents, "basis must be unit-independent");
        assert!(
            (a.value - b.value).abs() / a.value.abs().max(1e-300) < 1e-12,
            "pi value must be unit-invariant"
        );
    }
    verdict(
        "rg-002",
        "groups + pi basis exactly invariant under base-unit rescaling",
    );
}

#[test]
fn rg_003_validity_predicates_catch_seeded_misuse() {
    let registry = flux_model_cards();
    // Seeded misuse: a creeping-flow solver at Re = 10^4.
    let mut groups = BTreeMap::new();
    groups.insert("Re".to_string(), 1.0e4);
    groups.insert("Ma".to_string(), 0.05);
    let verdict_bad = admit(&registry, &groups, "flux.stokes-creeping").expect("known model");
    assert!(!verdict_bad.allowed, "Re=1e4 must refuse a creeping solver");
    assert!(
        verdict_bad.reasons.iter().any(|r| r.contains("Re")),
        "refusal must name the violated group"
    );
    // Ranked alternatives: LES (valid here) must rank strictly ahead of
    // potential flow and the refused Stokes card is not listed.
    let names: Vec<&str> = verdict_bad
        .alternatives
        .iter()
        .map(|(n, _)| n.as_str())
        .collect();
    let les_pos = names
        .iter()
        .position(|&n| n == "flux.les-ns")
        .expect("les listed");
    assert!(
        verdict_bad.alternatives[les_pos].1.abs() < 1e-15,
        "LES is valid at Re=1e4/Ma=0.05 (distance 0)"
    );
    // A valid choice admits cleanly.
    let mut ok_groups = BTreeMap::new();
    ok_groups.insert("Re".to_string(), 120.0);
    ok_groups.insert("Ma".to_string(), 0.02);
    let verdict_ok = admit(&registry, &ok_groups, "flux.free-surface-lbm").expect("model");
    assert!(verdict_ok.allowed);
    assert!(verdict_ok.reasons.is_empty());
    // Unknown models are structured errors, not panics.
    assert!(admit(&registry, &ok_groups, "flux.warp-drive").is_err());
    verdict(
        "rg-003",
        "creeping solver at Re=1e4 refused with named bound + LES ranked as distance-0 \
         alternative; valid LBM admitted",
    );
}

#[test]
fn rg_003b_unbounded_above_cards_admit_in_range_points() {
    // Regression: cards with a one-sided "no upper limit" domain
    // (potential-flow Re ≥ 1e4, euler-bernoulli slenderness ≥ 20, timoshenko
    // ≥ 5) previously used f64::INFINITY, which ValidityDomain treats as an
    // unusable/empty domain — so they could NEVER be admitted and rejected
    // valid points with an EMPTY reason list (violating Invariant 4, and
    // self-contradicting distance_to_validity = 0). They now use f64::MAX.
    let registry = flux_model_cards();

    // A slender beam well inside Euler-Bernoulli's [20, ∞) admits cleanly.
    let mut slender = BTreeMap::new();
    slender.insert("slenderness".to_string(), 30.0);
    let v = admit(&registry, &slender, "solid.euler-bernoulli").expect("known model");
    assert!(v.allowed, "slenderness=30 is valid for Euler-Bernoulli");
    assert!(v.reasons.is_empty(), "an admitted card has no reasons: {:?}", v.reasons);

    // "Unbounded above" really is unbounded: an astronomically slender member
    // must still admit (the pre-fix INFINITY made even this fail).
    let mut very_slender = BTreeMap::new();
    very_slender.insert("slenderness".to_string(), 1.0e6);
    assert!(
        admit(&registry, &very_slender, "solid.euler-bernoulli").expect("model").allowed,
        "there is no upper slenderness limit for Euler-Bernoulli"
    );

    // High-Re potential flow admits.
    let mut hi_re = BTreeMap::new();
    hi_re.insert("Re".to_string(), 5.0e4);
    assert!(
        admit(&registry, &hi_re, "flux.potential-flow").expect("model").allowed,
        "Re=5e4 is valid for potential flow"
    );

    // The lower bound is still enforced, and its violation is NAMED
    // (Invariant 4 holds): a stubby beam (slenderness 10 < 20) is refused.
    let mut stubby = BTreeMap::new();
    stubby.insert("slenderness".to_string(), 10.0);
    let refused = admit(&registry, &stubby, "solid.euler-bernoulli").expect("model");
    assert!(!refused.allowed, "slenderness=10 is below Euler-Bernoulli's floor");
    assert!(
        refused.reasons.iter().any(|r| r.contains("slenderness")),
        "refusal must name the violated bound, not be empty: {:?}",
        refused.reasons
    );

    verdict(
        "rg-003b",
        "one-sided (unbounded-above) cards admit in-range points and still name a \
         violated lower bound",
    );
}

#[test]
fn rg_004_recommended_scaling_improves_conditioning() {
    // Fixture: a 3-DOF elastic system assembled in raw SI with mixed
    // magnitudes (GPa stiffness, mm displacements, kN loads) vs the same
    // physics assembled in regime-recommended nondimensional variables.
    let inputs = vec![
        role(Role::Density, 7850.0, [-3, 1, 0, 0, 0, 0]),
        role(Role::Velocity, 2.0, [1, 0, -1, 0, 0, 0]),
        role(Role::Length, 0.004, [1, 0, 0, 0, 0, 0]),
    ];
    let map = ScalingMap::recommend(&inputs).expect("scales");
    // Raw system: rows mix stiffness (Pa·m ~ 1e9), geometry (1e-3) and a
    // dimensionless coupling — a classic badly-scaled assembly.
    let k = 2.0e11 * 0.004; // E·L, N/m-ish row scale
    let raw = [
        k,
        -k * 0.5,
        0.0, //
        -k * 0.5,
        2e-3,
        1.0e-3, //
        0.0,
        1.0e-3,
        2.5,
    ];
    let cond_raw = condition_number(&raw, 3).expect("cond raw");
    // Nondimensionalized: divide force-like rows by the pressure*area
    // scale and length-like columns by L* — here rows were built from
    // (stiffness, length, unity) quantities, so scale them by their
    // dimension factors.
    let f_force = map.factor([0, 1, -2, 0, 0, 0]); // force scale
    let f_len = map.factor([1, 0, 0, 0, 0, 0]);
    let scaled = [
        raw[0] * f_len / f_force,
        raw[1] * f_len / f_force,
        raw[2] * f_len / f_force,
        raw[3] * f_len / f_force,
        raw[4],
        raw[5],
        raw[6],
        raw[7],
        raw[8],
    ];
    let cond_scaled = condition_number(&scaled, 3).expect("cond scaled");
    println!(
        "{{\"suite\":\"fs-regime/conformance\",\"metric\":\"conditioning\",\
         \"before\":{cond_raw:.4e},\"after\":{cond_scaled:.4e}}}"
    );
    assert!(
        cond_scaled * 100.0 < cond_raw,
        "nondimensionalization must improve conditioning >=100x: {cond_raw:.3e} -> \
         {cond_scaled:.3e}"
    );
    verdict(
        "rg-004",
        "fixture system condition number improved and ledgered",
    );
}

#[test]
fn rg_005_similarity_engine_finds_the_neighbor() {
    let inputs = vec![
        role(Role::Density, 1.225, [-3, 1, 0, 0, 0, 0]),
        role(Role::Velocity, 1.47, [1, 0, -1, 0, 0, 0]),
        role(Role::Length, 0.001, [1, 0, 0, 0, 0, 0]),
        role(Role::DynViscosity, 1.834e-5, [-1, 1, -1, 0, 0, 0]),
    ];
    // Re = 1.225*1.47*0.001/1.834e-5 ~ 98.2 — near the Re=100 benchmark.
    let ev = assess(&inputs).expect("assess");
    let report = &ev.value;
    let re = report.groups["Re"];
    assert!((re - 98.2).abs() < 1.0, "Re fixture ~98 (got {re})");
    let near = report.nearest_benchmark.as_ref().expect("neighbor found");
    assert_eq!(near.name, "cylinder-crossflow-Re100");
    assert!(near.expectation.contains("Cd ~ 1.33"));
    assert_eq!(
        near.grade, "info",
        "Re=98 vs 100 is close (distance {})",
        near.distance
    );
    // Determinism: same inputs, same provenance hash.
    let ev2 = assess(&inputs).expect("assess again");
    assert_eq!(ev.provenance, ev2.provenance);
    assert_eq!(ev.value, ev2.value);
    verdict(
        "rg-005",
        "Re=98 matched to the Re=100 cylinder benchmark with Cd expectation; reports \
         reproducible",
    );
}

#[test]
fn rg_006_flagship_fixtures_match_hand_calculations() {
    // Spout (mid-tilt): hand values from the definitions.
    let spout = standard_groups(&spout_inputs()).expect("spout groups");
    let get = |name: &str| {
        spout
            .iter()
            .find(|g| g.name == name)
            .unwrap_or_else(|| panic!("{name} missing"))
            .value
    };
    let re = 1200.0 * 0.3 * 0.02 / 8.0; // 0.9
    let we = 1200.0 * 0.09 * 0.02 / 0.06; // 36
    let ca = 8.0 * 0.3 / 0.06; // 40
    let oh = 8.0 / (1200.0f64 * 0.06 * 0.02).sqrt(); // ~6.66
    let bo = 1200.0 * 9.806_65 * 0.0004 / 0.06; // ~78.45
    assert!((get("Re") - re).abs() / re < 1e-12);
    assert!((get("We") - we).abs() / we < 1e-12);
    assert!((get("Ca") - ca).abs() / ca < 1e-12);
    assert!((get("Oh") - oh).abs() / oh < 1e-9);
    assert!((get("Bo") - bo).abs() / bo < 1e-12);
    // The pour is viscous-dominated: the report says so and refuses LES.
    let ev = assess(&spout_inputs()).expect("assess spout");
    assert!(ev.value.dominant_balance.contains("viscous"));
    assert!(
        ev.value
            .valid_models
            .contains(&"flux.stokes-creeping".to_string()),
        "Re=0.9 admits creeping flow"
    );
    assert!(
        ev.value
            .invalid_models
            .iter()
            .any(|(n, _)| n == "flux.les-ns"),
        "Re=0.9 refuses LES"
    );

    // Ornithoid gait: Re, Ma, St from (air, V=8, chord=0.12, f=4Hz).
    let ornithoid = vec![
        role(Role::Density, 1.225, [-3, 1, 0, 0, 0, 0]),
        role(Role::Velocity, 8.0, [1, 0, -1, 0, 0, 0]),
        role(Role::Length, 0.12, [1, 0, 0, 0, 0, 0]),
        role(Role::DynViscosity, 1.81e-5, [-1, 1, -1, 0, 0, 0]),
        role(Role::SoundSpeed, 343.0, [1, 0, -1, 0, 0, 0]),
        role(Role::Frequency, 4.0, [0, 0, -1, 0, 0, 0]),
    ];
    let og = standard_groups(&ornithoid).expect("ornithoid");
    let find = |n: &str| og.iter().find(|g| g.name == n).expect("group").value;
    assert!((find("Re") - 1.225 * 8.0 * 0.12 / 1.81e-5).abs() < 1.0);
    assert!((find("Ma") - 8.0 / 343.0).abs() < 1e-12);
    assert!((find("St") - 4.0 * 0.12 / 8.0).abs() < 1e-12, "St = 0.06");

    // Frame member: slenderness + damping ratio.
    let frame = vec![
        role(Role::Length, 3.2, [1, 0, 0, 0, 0, 0]),
        role(Role::GyrationRadius, 0.04, [1, 0, 0, 0, 0, 0]),
        role(Role::Damping, 800.0, [0, 1, -1, 0, 0, 0]),
        role(Role::Stiffness, 2.0e6, [0, 1, -2, 0, 0, 0]),
        role(Role::Mass, 500.0, [0, 1, 0, 0, 0, 0]),
    ];
    let fg = standard_groups(&frame).expect("frame");
    let findf = |n: &str| fg.iter().find(|g| g.name == n).expect("group").value;
    assert!((findf("slenderness") - 80.0).abs() < 1e-12);
    let zeta = 800.0 / (2.0 * (2.0e6f64 * 500.0).sqrt()); // 0.0126...
    assert!((findf("zeta") - zeta).abs() / zeta < 1e-9);
    verdict(
        "rg-006",
        "spout Re/We/Ca/Oh/Bo, ornithoid Re/Ma/St, frame slenderness/zeta all match hand \
         calculations; spout admits creeping and refuses LES",
    );
}
